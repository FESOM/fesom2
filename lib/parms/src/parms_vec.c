/*--------------------------------------------------------------------
  parms_VecAXPY          : y = y + alpha x
  parms_VecAYPX          : y = scalar y + x
  parms_VecDOT           : inner product of two vectors.
  parms_VecDotArray      : inner products for the vector to an array
                           of vector object.

  parms_VecGetNorm2      : get the 2-norm of a vector.
  parms_VecPerm          : permute a vector object.

  parms_VecScale         : scale the components of a vector.
  A code fragment for using vector functions:

  parms_Vec vec;
  parms_Map map;

  // assume PE0 contains global variables (0, 2, 4, 6, 8), PE1
  // contains global variables (1, 3, 5, 7, 9). Both PE0 and PE1 call
  // the following parms_VecSetValues function. Each PE will pick and
  // insert variables which belong to the local PE according to the
  // partionting information stored in map.

  // all vector computation functions can be called
  parms_VecGetNorm2(vec, &value, map);
  parms_VecDot(vec, vec, &value, map);
  ...

  $Id: parms_vec.c,v 1.3 2007-05-10 14:11:21 zzli Exp $
  ------------------------------------------------------------------*/
#include "parms_sys.h"
#include "parms_viewer.h"
#include "parms_vec.h"
#include "parms_map_impl.h"
#if defined(__ICC)
#include <mathimf.h>
#else
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#endif  

typedef struct ext_vec_data {
  /*! \var im The external variables in global indexing
   */
  int     *im;
 /*!  \var values The external variable contribution.
  */
  FLOAT    *values;
 /*!  \var n Number if external variable contribution.
  */
  int    n; 
 /*!  \var space Size of work space.
  */
  int    space;    
} *ext_vec_data;


/** 
 * Scale a vector.
 *
 * All components of parms_Vec object self times scalar. 
 * \f$self = scalar \times self\f$.
 *  
 * @param self   A vector object.
 * @param scalar A scalar.      
 * 
 * @return 0 on success.
 */
int parms_VecScale(FLOAT *self, FLOAT scalar, parms_Map map) 
{
  int lsize, i;

  lsize = parms_MapGetLocalSize(map);
  for (i = 0; i < lsize; i++) 
  {
    self[i] *= scalar;
  }

  return 0;
}


/** 
 * Perform \f$self := scalar \times x + self\f$.
 * 
 * @param self   A vector object.      
 * @param x      Another vector object.
 * @param scalar A scalar.
 * 
 * @return 0 on success.
 */
int parms_VecAXPY(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map)   
{
  int i, lsize;

  lsize = parms_MapGetLocalSize(map);


  for (i = 0; i < lsize; i++) 
 {
    self[i] += scalar * x[i];
  }


  return 0;
}

/** 
 * Perform \f$self = scalar \times self + x\f$.
 * 
 * @param self    A vector object.      
 * @param x 	  Another vector object.
 * @param scalar  A scalar.
 * 
 * @return 0 on success.
 */
int parms_VecAYPX(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map) 
{
  int lsize, i;
  lsize = parms_MapGetLocalSize(map);

  for (i = 0; i < lsize; i++) {
    self[i] = scalar * self[i] + x[i];
  }

  return 0;
}


/* Dot product on local PE 
*/
static int vec_LDot(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map map)    
{
  int lsize, i;

  lsize = parms_MapGetLocalSize(map);

#if defined(DBL_CMPLX)
  FLOAT dot = 0.0 + 0.0*I;
  for (i = 0; i < lsize; i++)
  	dot  +=  self[i] * x[i];  
#else
  FLOAT dot = 0.0;  
  for (i = 0; i < lsize; i++) {
      dot  +=  self[i] * x[i];
  }
#endif

  *value = dot;


  return 0;
}

/* Dot product on local PE - uses complex conjugate
*/
static int vec_LDotc(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map map)    
{
  int lsize, i;

  lsize = parms_MapGetLocalSize(map);

#if defined(DBL_CMPLX)
  FLOAT dot = 0.0 + 0.0*I;
  for (i = 0; i < lsize; i++)
  	dot  +=  self[i] * conj(x[i]);
#else
  FLOAT dot = 0.0;  
  for (i = 0; i < lsize; i++) {
      dot  +=  self[i] * x[i];
  }
#endif

  *value = dot;

  return 0;
}

/** 
 * Perform the (global) inner product of two vectors.
 * 
 *  value = self x^{T}. If self
 * 
 *
 * @param self   A vector object.                  
 * @param x 	 Another vector object.            
 * @param value  The inner product returned.       
 * 
 * @return 0 on success.
 */
int parms_VecDOT(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map is)
{
  MPI_Comm comm;

  comm = is->comm;

  if (is->isserial) {
    vec_LDot(self, x, value, is);
  }
  else {
    FLOAT ldot;
    vec_LDot(self, x, &ldot, is);
    
#if defined(DBL_CMPLX)
    MPI_Allreduce(&ldot, value, 1, MPI_CMPLX, MPI_CMPLX_SUM, comm);
#else
    MPI_Allreduce(&ldot, value, 1, MPI_DOUBLE, MPI_SUM, comm);
#endif    
  }

  return 0;
}

/** 
 * Perform the (global) inner product of two vectors.
 * 
 * If self and x are real vectors, value = self x^{T}. If self
 * and x are complex vectors, value = self \overline{x}^{T}.
 *
 * @param self   A vector object.                  
 * @param x 	 Another vector object.            
 * @param value  The inner product returned.       
 * 
 * @return 0 on success.
 */
int parms_VecDOTC(FLOAT *self, FLOAT *x, REAL *value, parms_Map is)
{
  MPI_Comm comm;

  comm = is->comm;
  
  if (is->isserial) {
#if defined(DBL_CMPLX)  
    FLOAT val;    
    vec_LDotc(self, x, &val, is);
    *value = creal(val);    
#else
    vec_LDot(self, x, value, is);
#endif   
  }
  else {
#if defined(DBL_CMPLX)  
    FLOAT ldot;
    REAL lval;
    vec_LDotc(self, x, &ldot, is);
    lval = creal(ldot);

    MPI_Allreduce(&lval, value, 1, MPI_DOUBLE, MPI_SUM, comm);
#else
    FLOAT ldot;
    vec_LDot(self, x, &ldot, is);
    MPI_Allreduce(&ldot, value, 1, MPI_DOUBLE, MPI_SUM, comm);
#endif    
  }

  return 0;
}

/** 
 * Return the 2-norm of the vector.
 * 
 * @param self  A vector object.    
 * @param value The 2-norm returned.
 * 
 * @return 0 on success.
 */
int parms_VecGetNorm2(FLOAT *self, REAL *value, parms_Map is)
{ 

  if (is->isserial) {

    FLOAT dot;
    vec_LDotc(self, self, &dot, is);
    *value = ABS_VALUE(dot);
    *value = sqrt(*value);

  }
  else {
    REAL dot;
#if defined(DBL_CMPLX)    
    parms_VecDOTC(self, self, &dot, is);
#else
    parms_VecDOT(self, self, &dot, is);
#endif 
    *value = fabs(dot);
    *value = sqrt(*value);
  }
  return 0;
}

/** 
 * Perform the inner product between self and an array of vectors.  
 *
 * The pseudo code:
 *
 *  \f{verbatim}
 *  for (i = 0; i < n; i++) {
 *    result[i] = self * vecarray[i];
 *  }
 *  \f}
 *
 * @param self       A vector object.                              
 * @param n 	     The size of vecarray.                         
 * @param vecarray   An array of vector objects.                   
 * @param aux 	     An auxiliary array.                           
 * @param result     An array of size n to store inner products.   
 * 
 * @return 0 on success.
 */
int parms_VecDotArray(FLOAT *self, int n, FLOAT
			**vecarray, FLOAT *result, parms_Map is)
{
  int i;
  FLOAT *aux;
  MPI_Comm comm;

  comm = is->comm;
  if (is->isserial) {
    for (i = 0; i < n; i++)
      vec_LDotc(self, vecarray[i], &result[i], is);
  }
  else {
    PARMS_NEWARRAY(aux, n);
    for (i = 0; i < n; i++) 
      vec_LDotc(self, vecarray[i], &aux[i], is);
#if defined(DBL_CMPLX)
      MPI_Allreduce(aux, result, n, MPI_CMPLX, MPI_CMPLX_SUM, comm);
#else
      MPI_Allreduce(aux, result, n, MPI_DOUBLE, MPI_SUM, comm);
#endif    
    PARMS_FREE(aux); 
  }
  return 0;
}

int parms_VecPerm(FLOAT *self, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  FLOAT        *newvec;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  size = parms_MapGetLocalSize(map);
  if (map->isperm && (map->isvecperm == false)) { 
    perm = map->perm;
    PARMS_NEWARRAY(newvec, size);
    for (i = 0; i < size; i++) {
      k = perm[i];
      newvec[k] = self[i];
    }
    memcpy(self, newvec, size*sizeof(FLOAT));

    PARMS_FREE(newvec);
  }
  return 0;
}

int parms_VecPermAux(FLOAT *self, FLOAT *aux, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  size = parms_MapGetLocalSize(map);
  if (map->isperm) { 
    perm = map->perm;
    for (i = 0; i < size; i++) {
      k = perm[i];
      aux[k] = self[i];
    }
  }
  else
  {
    for (i = 0; i < size; i++) {
      aux[i] = self[i];
    }
   }     
  return 0;
}

int parms_VecInvPerm(FLOAT *self, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  FLOAT        *newvec;
  
  size = parms_MapGetLocalSize(map);
  if (map->isperm && (map->isvecperm == false)) { 
    perm = map->iperm;
    PARMS_NEWARRAY(newvec, size);
    for (i = 0; i < size; i++) {
      k = perm[i];
      newvec[k] = self[i];
    }
    memcpy(self, newvec, size*sizeof(FLOAT));

    PARMS_FREE(newvec);
  }
  return 0;
}

int parms_VecInvPermAux(FLOAT *self, FLOAT *aux, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  size = parms_MapGetLocalSize(map);
  if (map->isperm) { 
    perm = map->iperm;
    for (i = 0; i < size; i++) {
      k = perm[i];
      aux[k] = self[i];
    }
  }
  else
  {
    for (i = 0; i < size; i++) {
      aux[i] = self[i];
    }
   }  
  return 0;
}

/** 
 * Insert values to parms_Vec object self.
 * 
 *  A pseudo code from the global point of view:
 *
 *  \f{verbatim}
 *  for (i = 0; i < m; i++) {
 *    self[im[i]] = values[i]; 
 *  }
 *  \f}
 *  
 * @param self   A vector object.                          
 * @param m      The number of variables to be inserted.     
 * @param im     An array of global variable indices.     
 * @param value  An array of values to be inserted to self.s 
 * @param mode   The style of set values:
 *               -ADD    add values to parms_Vec object self. 
 *               -INSERT assign values to parms_Vec object self.
 *               
 * @return 0 on success.
 */
int parms_VecSetValues(FLOAT *self, int m, int *im, FLOAT
			       *values, INSERTMODE mode, parms_Map map) 
{
  int lsize, rindex, index, lrindex, *rp, i;
  int offset = map->start;

  lsize = parms_MapGetLocalSize(map);

  if (!map->isserial) {
    if (map->isvecperm) {
      for (i = 0; i < m; i++) {
       rindex = im[i] - offset;
	rp = parms_MapGlobalToLocal(map, rindex);
	if ((rp != NULL) && ((index = *rp) < lsize)) {
	  lrindex = map->perm[index];
	  if (mode == ADD) {
	    self[lrindex] += values[i];
	  }
	  else if (mode == INSERT) {
	    self[lrindex] = values[i];
	  }
	}
      }
    }
    else {
      for (i = 0; i < m; i++) {
       rindex = im[i] - offset;
	rp = parms_MapGlobalToLocal(map, rindex);
	if ((rp != NULL) && ((lrindex = *rp) < lsize)) {
	  if (mode == ADD) {
	    self[lrindex] += values[i];
	  }
	  else if(mode == INSERT) {
	    self[lrindex] = values[i];
	  }
	}
      }
    }
  }
  else {
    for (i = 0; i < m; i++) {
      rindex = im[i] - offset;
      if (map->isvecperm) {
	rindex = map->perm[rindex];
      }
      if (mode == ADD) {
	self[rindex] += values[i];
      }
      else if(mode == INSERT) {
	self[rindex] = values[i];
      }
    }
  }
  return 0;
}

/** 
 * Insert values to parms_Vec object self. This assumes the vector 
 * values are being set element-by-element. A call to parms_VecAssembleElementVector
 * is required to complete the vector once all entries have been added.
 * 
 *  A pseudo code from the global point of view:
 *
 *  \f{verbatim}
 *  for (i = 0; i < m; i++) {
 *    self[im[i]] = values[i]; 
 *  }
 *  \f}
 *  
 * @param self   A vector object.                          
 * @param m      The number of variables to be inserted.     
 * @param im     An array of global variable indices.     
 * @param value  An array of values to be inserted to self.s 
 * @param mode   The style of set values:
 *               -ADD    add values to parms_Vec object self. 
 *               -INSERT assign values to parms_Vec object self.
 *               
 * @return 0 on success.
 */
int parms_VecSetElementVector(FLOAT *self, int m, int *im, FLOAT
			       *values, INSERTMODE mode, parms_Map map) 
{
  int lsize, rindex, index, lrindex, *rp, i,pos,k;
  int offset = map->start, space;
  ext_vec_data vdata;

  lsize = parms_MapGetLocalSize(map);
  
  if(map->isserial) {
/* call vecsetvalues for serial version */
    return parms_VecSetValues(self, m, im, values, mode, map);  
  }  

/* Parallel call */
  if(map->isdatalloc)
  {
    vdata = (ext_vec_data)map->data;
  }
  else
  {
  /* Allocate some memory for vdata */
    map->isdatalloc = true;
    PARMS_NEW(vdata);
    vdata->space = lsize;
    PARMS_NEWARRAY(vdata->im, vdata->space);
    PARMS_NEWARRAY0(vdata->values, vdata->space);
    vdata->n = 0;
    map->data = vdata;
  }

  if (map->isvecperm) {
    for (i = 0; i < m; i++) {
     rindex = im[i] - offset;
     rp = parms_MapGlobalToLocal(map, rindex);
     if ((rp != NULL) && ((index = *rp) < lsize)) {
	  lrindex = map->perm[index];
	  if (mode == ADD) {
	    self[lrindex] += values[i];
	  }
	  else if (mode == INSERT) {
	    self[lrindex] = values[i];
	  }
	}
      else{
	/* external contribution - first check if there is enough memory*/
	  space = vdata->space;
        if(vdata->n == space){
	  /* reallocate memory for holding new entry */
	  space += 10;
        PARMS_RESIZE(vdata->im,   space);
        PARMS_RESIZE(vdata->values, space);
        vdata->space = space;
        }	  

        pos = 0;
        for(k = 0; k<vdata->n; k++) /* check if row has been added already */
        {
          if(vdata->im[k] == im[i])
              break;
        }
        if(k < vdata->n) // already added row
        {
          pos = k;
          vdata->values[pos] += values[i];
        }
        else          // new row contribution
        {
          pos = vdata->n;
          vdata->im[vdata->n] = im[i];
          vdata->values[vdata->n++] = values[i];
        }
	}
     }
    }
    else {
      for (i = 0; i < m; i++) {
       rindex = im[i] - offset;
	rp = parms_MapGlobalToLocal(map, rindex);
	if ((rp != NULL) && ((lrindex = *rp) < lsize)) {
	  if (mode == ADD) {
	    self[lrindex] += values[i];
	  }
	  else if(mode == INSERT) {
	    self[lrindex] = values[i];
	  }
	}
      else{
	/* external contribution - first check if there is enough memory*/
	  space = vdata->space;
        if(vdata->n == space){
	  /* reallocate memory for holding new entry */
	  space += 10;
        PARMS_RESIZE(vdata->im,   space);
        PARMS_RESIZE(vdata->values, space);
        vdata->space = space;
        }	  

        pos = 0;
        for(k = 0; k<vdata->n; k++) /* check if row has been added already */
        {
          if(vdata->im[k] == im[i])
              break;
        }
        if(k < vdata->n) // already added row
        {
          pos = k;
          vdata->values[pos] += values[i];
        }
        else          // new row contribution
        {
          pos = vdata->n;
          vdata->im[pos] = im[i];
          vdata->values[pos] = values[i];
          vdata->n++;
        }
	}
    }
   }
 
  return 0;
}

/** 
 * Completes setting up values for the distributed vector
 *
 * @param self   A vector object.                          
 * @param map    A pARMS map object
 *               
 * @return 0 on success.
 */
int parms_VecAssembleElementVector(FLOAT *self, parms_Map map) 
{
  FLOAT *gvec;
  int i, k, gnodv,lsize, mypid, npro, lrindex, index;
  int *info, *disp, *pos, *gextvars;
  MPI_Comm comm;
  ext_vec_data vdata;
  
  lsize = map->lsize;
  mypid = map->pid;
  npro = map->npro;
  comm = map->comm;

  if(map->isdatalloc){
    vdata = (ext_vec_data)map->data;
/* Allocate some memory */
    PARMS_NEWARRAY(info, npro);
    PARMS_NEWARRAY(disp, npro+1);

/* collect contributions from other processors */
    MPI_Allgather(&vdata->n, 1, MPI_INT, info, 1, MPI_INT, comm);
    disp[0] = 0;
    for (i = 0; i < npro; i++) {
      disp[i+1] = disp[i] + info[i];
    }
  
    gnodv = disp[npro];
    PARMS_NEWARRAY(gextvars, gnodv);
    PARMS_NEWARRAY0(gvec, gnodv);
    
    MPI_Allgatherv(vdata->im, vdata->n, MPI_INT, gextvars,
		   info, disp, MPI_INT, comm);
#if defined(DBL_CMPLX)
    MPI_Allgatherv(vdata->values, vdata->n, MPI_CMPLX, gvec,
		   info, disp, MPI_CMPLX, comm);   
#else
    MPI_Allgatherv(vdata->values, vdata->n, MPI_DOUBLE, gvec,
		   info, disp, MPI_DOUBLE, comm);  
#endif    
    for(i=0; i<npro; i++)
    {
      if(mypid != i)
      {
        for(k = disp[i]; k<disp[i+1]; k++)
        {
           index = gextvars[k] - map->start;
           pos = parms_MapGlobalToLocal(map,index);
           if ((pos != NULL) && ((lrindex = *pos) < lsize)) {
             self[lrindex] += gvec[k];
           }
         }
       }
     }    
/* Free memory */
    PARMS_FREE(info);
    PARMS_FREE(disp);
    PARMS_FREE(gvec);
    PARMS_FREE(gextvars);
    PARMS_FREE(vdata->im);
    PARMS_FREE(vdata->values);
    PARMS_FREE(vdata);
    map->isdatalloc = false;
  } 
   
  return 0;
}

/** 
 * Gather distributed vector to a global array.
 * 
 * @param self The distributed vector.
 * @param ga A global vector.
 * 
 * @return 0 on success.
 */
int parms_VecGather(FLOAT *self, FLOAT *ga, parms_Map map)
{
  BOOL isserial;
  int i, j, incx, lsize, index, npro, pid;
  int *num_recv, *ind_array, *ind_snd, maxnum;
  FLOAT *wk;

  isserial = map->isserial;
  pid      = map->pid;
  npro     = map->npro;

  lsize = parms_MapGetLocalSize(map);
  if (isserial) {
      for (i = 0; i < lsize; i++) ga[i] = self[i];
  }
  else {
    PARMS_NEWARRAY0(num_recv, npro);
    PARMS_NEWARRAY(ind_snd,   lsize);
    MPI_Allgather(&lsize, 1, MPI_INT, num_recv, 1, MPI_INT, map->comm);
    maxnum = 0;
    for (i = 0; i < npro; i++) {
      if (maxnum < num_recv[i]) {
	maxnum = num_recv[i];
      }
    }
    PARMS_NEWARRAY(wk, maxnum);
    PARMS_NEWARRAY(ind_array, maxnum);

    for (i = 0; i < lsize; i++) {
	index = map->lvars[i];
	ga[index] = self[i];
	ind_snd[i] = map->lvars[i];
    }

#if defined(DBL_CMPLX)    
    for (i = 0; i < npro; i++) {
      if (pid == i) {
	MPI_Bcast(ind_snd, lsize, MPI_INT, i, map->comm);
	MPI_Bcast(self, lsize, MPI_CMPLX, i, map->comm);
      }
      else {
	MPI_Bcast(ind_array, num_recv[i], MPI_INT, i, map->comm);
	MPI_Bcast(wk, num_recv[i], MPI_CMPLX, i, map->comm);
	for (j = 0; j < num_recv[i]; j++) {
	  index = ind_array[j];
	  ga[index] = wk[j];
	}
      }
    }
#else
    for (i = 0; i < npro; i++) {
      if (pid == i) {
	MPI_Bcast(ind_snd, lsize, MPI_INT, i, map->comm);
	MPI_Bcast(self, lsize, MPI_DOUBLE, i, map->comm);
      }
      else {
	MPI_Bcast(ind_array, num_recv[i], MPI_INT, i, map->comm);
	MPI_Bcast(wk, num_recv[i], MPI_DOUBLE, i, map->comm);
	for (j = 0; j < num_recv[i]; j++) {
	  index = ind_array[j];
	  ga[index] = wk[j];
	}
      }
    }
#endif    
    PARMS_FREE(ind_snd);
    PARMS_FREE(ind_array);
    PARMS_FREE(wk);
  }

  return 0;
}

