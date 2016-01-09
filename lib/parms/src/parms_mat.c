/*--------------------------------------------------------------------
  parms_MatVec          : perform the multiplication of
                            matrix-vector product. 
  parms_MatCreate         : create a matrix object.
  parms_MatExtend         : extend submatrix by including equations
                            correspond to the immediate neighbouring
                            variables. 
  parms_MatFree           : free the memory for a matrix object.
  parms_MatFreeSubMat     : free the memory for the submatrix object.
  parms_MatGetCommHandler : get the communication handler.
  parms_MatGetDiag        : get the diagonal part of the local matrix. 
  parms_MatGetSubMat      : get the local matrix.
  parms_MatMVPY           : perform \f$z = \alpha \times self \times x
                                           + beta \times y\f$. 
  parms_MatSetCommType    : set the communication type.
  parms_MatSetValues      : insert/add values to the matrix object.
  parms_MatSetup          : set up a matrix object.
  parms_MatVecOffDiag     : perform the multiplication of the
                            off-diagonal matrix and the external vars. 
  parms_MatView           : dump a matrix object.

  A code fragment for using matrix functions:

  parms_Mat mat;
  parms_Map map;
  parms_PC  pc;

  ...

  // create a matrix object.
  parms_MatCreate(&mat, map);

  // insert/add values into the matrix
  parms_MatSetValues(mat, ...);

  // setup the matrix. 
  // this function divides local variables into two categories:
  // interior variables and interface variables. interior variables
  // listed first. when the matrix-vector product is performed, the
  // vector is permuted automatically if the vector is created
  parms_MatSetup(mat);

  // matrix-vector product
  parms_MatVec(mat, vec, y);

  // create a preconditioner object.
  parms_PCCreate(&pc, mat);

  $Id: parms_mat.c,v 1.4 2006-12-18 22:23:40 zzli Exp $
  ------------------------------------------------------------------*/
#include "include/parms_mat_impl.h"

/*
int parms_MatCreate_vcsr(parms_Mat self);
int parms_MatCreate_dvcsr(parms_Mat self);
int parms_MatFree_dvcsr(parms_Mat *self);
int parms_MatView_vcsr(parms_Mat self, parms_Viewer v);
int parms_MatView_dvcsr(parms_Mat self, parms_Viewer v);
*/
/** 
 * Free the parms_Mat object pointed to by self.
 * 
 * @param self A pointer to a parms_Mat object.
 * 
 * @return 0 on success.
 */
int parms_MatFree(parms_Mat *self) 
{
  int         m, i, nnz;

  (*self)->ref--;
  if ((*self)->ref == 0 ) {
    if(!(*self)->isserial)
	parms_MatFree_dvcsr(self);
    parms_MapFree(&(*self)->is);
    if ((*self)->isalloc) {
      m = (*self)->m;
      for (i = 0; i < m; i++) {
	  nnz = (*self)->aux_data->nnzrow[i];
	  if (nnz) {
	    PARMS_FREE((*self)->aux_data->pj[i]);
	    PARMS_FREE((*self)->aux_data->pa[i]);
	  }
      }
      PARMS_FREE((*self)->aux_data->pa);
      PARMS_FREE((*self)->aux_data->pj);
      PARMS_FREE((*self)->aux_data->nnzrow);
      PARMS_FREE((*self)->aux_data);
      PARMS_FREE((*self)->ops);
      PARMS_FREE(*self);
    }
  }
  return 0;
}

/** 
 * Dump parms_Mat object.
 * 
 * @param self A pointer to a parms_Mat object.
 * @param v    A parms_Viewer object.   
 * 
 * @return 0 on success.
 */
int parms_MatView(parms_Mat self, parms_Viewer v)
{
  if((self)->isserial)
  	return parms_MatView_vcsr(self, v);
  else
	return parms_MatView_dvcsr(self, v);
}

/** 
 * Dump parms_Mat object in Coordinate format.
 * 
 * @param self A pointer to a parms_Mat object.
 * @param v    A parms_Viewer object.   
 * 
 * @return 0 on success.
 */
int parms_MatViewCOO(parms_Mat self, parms_Viewer v)
{
  if((self)->isserial)
  	return parms_MatViewCOO_vcsr(self, v);
  else
	return parms_MatViewCOO_dvcsr(self, v);
}

/** 
 * Set the communication type.
 *
 * Set the communication style across processors.
 * communication style:
 *  -P2P       point-to-point (data copied to/from auxilliary buffers).
 *  -DERIVED   derived datatype.
 *
 * @param self  A matrix object.
 * @param ctype Communication style:
 *              - P2P     point-to-point (data copied to/from
 *                        auxilliary buffers).
 *              - DERIVED derived datatype.
 *              
 * @return 0 on success.
 */
int parms_MatSetCommType(parms_Mat self, COMMTYPE ctype)
{
  return self->ops->setcommtype(self, ctype);
}

/** 
 * Create a parms_Mat object.
 *
 * Create a parms_Mat object based on data distribution layout map.
 * 
 * @param self  A pointer to the parms_Mat object created. 
 * @param map   A parms_Map object, which describes the data
 *              distribution among processors.  
 * 
 * @return 0 on success.
 */
int parms_MatCreate(parms_Mat *self, parms_Map map)
{
  parms_Mat  new_mat;

  PARMS_NEW0((new_mat));
  new_mat->ref = 1;
  PARMS_NEW0((new_mat)->ops);
  new_mat->is               = map;
  map->ref++;
  new_mat->type             = MAT_NULL;
  new_mat->ilutype          =PCILU0;     /*-- default local precon is ilu0. See parms_pc.c ---*/
  new_mat->issetup          = false;
  new_mat->isperm           = map->isperm;
  new_mat->isalloc          = false;
  new_mat->isreset          = false; 
  new_mat->isassembled        = false;   

  new_mat->m = parms_MapGetLocalSize(map);
  new_mat->n = new_mat->m;
  new_mat->M = parms_MapGetGlobalSize(map);
  new_mat->N = new_mat->M;
  new_mat->isserial = map->isserial;
  *self = new_mat;

  return 0;
}

/** 
 * Perform \f$y = self \times x\f$.
 * 
 * @param self A parms_Mat object.      
 * @param x    A vector.      
 * @param y    Another vector object.
 * 
 * @return 0 on success.
 */
int parms_MatVec(parms_Mat self, FLOAT *x, FLOAT *y)
{
  return self->ops->apply(self, x, y);
}

/** 
 * Set up parms_Mat object self.
 *
 * This is the most important function for the parms_Mat object. This
 * function combines the function bdry and setup in the old version 
 * of pARMS. The function sets up the data structure needed by the
 * distributed matrix-vector multiplication, divides the variables on
 * the local processors into two categories: interior and interface
 * variables.
 * 
 * @param self A parms_Mat object. 
 * 
 * @return 0 on success.
 */
int parms_MatSetup(parms_Mat self)
{      
  /* free and create new table for offdiagonal column count */
  if(self->odtable)
    parms_TableFree(&self->odtable);

  if ((!self->issetup) && (self->type == MAT_NULL)) {
    self->type = MAT_VCSR;
    if(self->isserial){
	parms_MatCreate_vcsr(self);
    }
    else{
      if(!self->isreset) 
	  parms_MatCreate_dvcsr(self);
    }
  }
  
  return self->ops->setup(self);
}

/** 
 * Insert/add values to the parms_Mat object self.
 * 
 * @param self    A parms_Mat object.                                
 * @param m 	The number of rows inserted.      
 * @param im 	An array of row global indices.                                                 
 * @param ia 	An array of pointer to the beginning of each row in
 *                array ja.    
 * @param ja      An array of column global indices.     
 * @param values  An array of values.                    
 * @param mode 	Insert value mode:                     
 * 		      - INSERT  insert values to parms_Mat self.
 * 		      - ADD     add values to parms_Mat self.
 *
 * NOTE: New entries are always inserted first, so mode does not
 *       matter if this is a new entry. Subsequent calls will either 
 *       replace (mode = INSERT) or add (mode = ADD) to the existing 
 *       entries in a particular position 
 * 		  
 * @return 0 on success.
 */
int parms_MatSetValues(parms_Mat self, int m, int *im, int *ia,
		       int *ja, FLOAT *values, INSERTMODE mode)  
{
  parms_Map  is;
  parms_vcsr aux_data;
  BOOL isserial, isalloc;
  int i, j, k, offset, start, end, rindex, cindex;
  int space, numinf, index, lrindex, lcindex, *rja;
  int *perm, *rp, *cp, *odp,inc, lsize, size;
  FLOAT *ra;
  BOOL found;

  isalloc  = self->isalloc;
  is       = self->is;
  offset    = is->start;
  lsize    = parms_MapGetLocalSize(is);
  isserial = is->isserial;

  if (isalloc) {
    aux_data = self->aux_data;     
  }
  else {
    self->isalloc = true;
    PARMS_NEW(aux_data);
    aux_data->n = lsize;
    PARMS_NEWARRAY(aux_data->space,   aux_data->n);
    PARMS_NEWARRAY0(aux_data->nnzrow, aux_data->n);
    PARMS_NEWARRAY(aux_data->pj,      aux_data->n);
    PARMS_NEWARRAY(aux_data->pa,      aux_data->n);
    for (i = 0; i < aux_data->n; i++) {
      aux_data->space[i] = 8;
      PARMS_NEWARRAY(aux_data->pj[i], aux_data->space[i]); 
      PARMS_NEWARRAY(aux_data->pa[i], aux_data->space[i]);      
    }
    self->aux_data = aux_data;
    /* create table for tracking off-diagonal column count */
    parms_TableCreate(&self->odtable, NULL, lsize);
  }

  if (isserial) {
    for (i = 0; i < m; i++) {
      rindex = im[i] - offset;
      ra     = aux_data->pa[rindex];
      rja    = aux_data->pj[rindex];
      start  = ia[i] - offset;
      end    = ia[i+1] - offset;
      for (j = start; j < end; j++) {
	cindex = ja[j] - offset;
	space = aux_data->space[rindex];
	found = false;
	for (k = 0; k < aux_data->nnzrow[rindex]; k++) {
	  if (rja[k] == cindex) {
	    found = true;
	    index = k;
	    break;
	  }
	}
	if (found) {
	  if (mode == INSERT) {
	    ra[index] = values[j];
	  }
	  else if (mode == ADD) {
	    ra[index] += values[j];
	  }
	}
	else { /* insert a new entry */
	  if (space == aux_data->nnzrow[rindex]) {
	    /* reallocate memory for holding new entry */
	    space += 8;
	    PARMS_RESIZE(aux_data->pa[rindex], space);
	    PARMS_RESIZE(aux_data->pj[rindex], space);
	    aux_data->space[rindex] = space;	    
	    ra  = aux_data->pa[rindex];
	    rja = aux_data->pj[rindex];
	  }
	  rja[aux_data->nnzrow[rindex]]  = cindex;
	  ra[aux_data->nnzrow[rindex]++] = values[j];
	}
      }
    }
  }
  else {
    perm    = is->perm;
    numinf  = is->ninf;
/* Check if matrix is to be re-used with same non-zero pattern */
    if(self->isreset && (self->resetpattern == SAME_NONZERO_STRUCTURE))
    {
      for (i = 0; i < m; i++) {
      rindex = im[i] - offset;
      rp = parms_MapGlobalToLocal(is, rindex);
      /* this row resides the local processor */
      if (rp != NULL && *rp < lsize) {
	/* calculate the local row index */
	lrindex = perm[*rp];
	ra      = aux_data->pa[lrindex];
	rja     = aux_data->pj[lrindex];

	start = ia[i]-offset;
	end   = ia[i+1]-offset;
	for (j = start; j < end; j++) {
	  cindex = ja[j] - offset;
	  /* the local processor contains this column index */   
	  cp = parms_MapGlobalToLocal(is, cindex);
	  if (cp != NULL && *cp < lsize) {
	    /* local column index */
	    lcindex = perm[*cp];

	    found = false;
	    for (k = 0; k < aux_data->nnzrow[lrindex]; k++) {
	      if (rja[k] == lcindex) {
	        found = true;
	        index = k;
	        break;
	      }
	     }
	     if (found) {
	       if (mode == INSERT) {
	         ra[index] = values[j];
	       }
	       else if (mode == ADD) {
	         ra[index] += values[j];
	       }
	     }
	  }
	  else {
	    if(cp != NULL){
	      lcindex = *cp;
	    
	      found = false;
	      for (k = 0; k < aux_data->nnzrow[lrindex]; k++) {
	        if (rja[k] == lcindex) {
	          found = true;
	          index = k;
	          break;
	        }
	       }
	       if (found) {
	         if (mode == INSERT) {
	           ra[index] = values[j];
	         }
	         else if (mode == ADD) {
	           ra[index] += values[j];
	         }
	       }
	     }
	  }
	}
     }
    }
   }
   else{        
    for (i = 0; i < m; i++) {
      rindex = im[i] - offset;
      rp = parms_MapGlobalToLocal(is, rindex);
      /* this row resides the local processor */
      if (rp != NULL && *rp < lsize) {
	/* calculate the local row index */
	lrindex = *rp;
	ra      = aux_data->pa[lrindex];
	rja     = aux_data->pj[lrindex];

	start = ia[i]-offset;
	end   = ia[i+1]-offset;
	for (j = start; j < end; j++) {
	  cindex = ja[j] - offset;
	  /* the local processor contains this column index */   
	  cp = parms_MapGlobalToLocal(is, cindex);
	  if (cp != NULL && *cp < lsize) {
	    /* local column index */
	    lcindex = *cp;
	  }
	  else {
	    if (perm[lrindex] == -1) { /* not marked yet */
	      perm[lrindex] = lsize-1-numinf;
	      numinf++;
	    }
	    lcindex = -cindex-1;
	  }
	  space   = aux_data->space[lrindex];
	  found = false;
	  for (k = 0; k < aux_data->nnzrow[lrindex]; k++) {
	    if (rja[k] == lcindex) {
	      found = true;
	      index = k;
	      break;
	    }
	  }
	  if (found) {
	    if (mode == INSERT) {
	      ra[index] = values[j];
	    }
	    else if (mode == ADD) {
	      ra[index] += values[j];
	    }
	  }
	  else { /* insert the new entry */
	    if (space == aux_data->nnzrow[lrindex]) {
	      /* reallocate memory for holding new entry */
	      space += 8;
	      PARMS_RESIZE(aux_data->pa[lrindex], space);
	      PARMS_RESIZE(aux_data->pj[lrindex], space);
	      aux_data->space[lrindex] = space;
	      ra  = aux_data->pa[lrindex];
	      rja = aux_data->pj[lrindex];
	    }
	    if (cp == NULL) {
	      size = parms_TableGetSize(is->table);
	      parms_TablePut(is->table, cindex, size);
           /* Initialize column count */
	      inc = 0;
	      parms_TablePut(self->odtable,cindex,inc);
	    }
	    else if((cp !=NULL) && (lcindex < 0)){
	      /* update off-diagonal column count */	    
	       odp = parms_TableGet(self->odtable, cindex);
             inc = *odp + 1;
	       parms_TablePut(self->odtable,cindex,inc);
	    }
	    
	    rja[aux_data->nnzrow[lrindex]]  = lcindex;
	    ra[aux_data->nnzrow[lrindex]++] = values[j];
	   }
	  }
       }
      }
     is->ninf = numinf;
    }
  }

  return 0;
}



/** 
 * Insert/add values to the parms_Mat object self. This assumes matrix 
 * values are being added element-by-element
 * 
 * @param self    A parms_Mat object.                                
 * @param m 	The number of rows inserted.      
 * @param im 	An array of row global indices.                                                 
 * @param ia 	An array of pointer to the beginning of each row in
 *                array ja.    
 * @param ja      An array of column global indices.     
 * @param values  An array of values.                    
 * @param mode 	Insert value mode:                     
 * 		      - INSERT  insert values to parms_Mat self.
 * 		      - ADD     add values to parms_Mat self.
 *
 * NOTE: New entries are always inserted first, so mode does not
 *       matter if this is a new entry. Subsequent calls will either 
 *       replace (mode = INSERT) or add (mode = ADD) to the existing 
 *       entries in a particular position 
 * 		  
 * @return 0 on success.
 */
int parms_MatSetElementMatrix(parms_Mat self, int m, int *im, int *ia,
		       int *ja, FLOAT *values, INSERTMODE mode)  
{
  parms_Map  is;
  parms_vcsr aux_data, ext_data;
  BOOL isserial, isalloc, isassembled;
  int i, j, k, offset, start, end, rindex, cindex;
  int space, numinf, index, lrindex, lcindex, *rja, *pj;
  int *perm, *rp, *cp,*odp,inc, *ext_im, lsize, size, ext_nnz;
  FLOAT *ra, *pa;
  BOOL found;

  isalloc  = self->isalloc;
  is       = self->is;
  offset    = is->start;
  lsize    = parms_MapGetLocalSize(is);
  isserial = is->isserial;
  isassembled = self->isassembled;
  
  if(isserial)
    parms_MatSetValues(self, m, im, ia, ja, values, mode);
  else{
    if (isalloc) {
      aux_data = self->aux_data;
      if(isassembled)
      {
      /* reallocate some data for the external data (contributions to 
       * other processors) from the current element matrix 
       */
        PARMS_NEW(ext_data);
        is->n_ext = is->n_ext > 0 ? is->n_ext : lsize;
        ext_data->n = is->n_ext;
        PARMS_NEWARRAY(ext_data->space,   ext_data->n);
        PARMS_NEWARRAY0(ext_data->nnzrow, ext_data->n);
        PARMS_NEWARRAY(ext_data->pj,      ext_data->n);
        PARMS_NEWARRAY(ext_data->pa,      ext_data->n);
        PARMS_NEWARRAY(ext_im,      ext_data->n);
        for (i = 0; i < ext_data->n; i++) {
          ext_data->space[i] = 8;
          PARMS_NEWARRAY(ext_data->pj[i], ext_data->space[i]); 
          PARMS_NEWARRAY(ext_data->pa[i], ext_data->space[i]);      
        }
        self->ext_data = ext_data; 
        is->ext_im = ext_im;
        /* reinitialize number of external contribution to zero */    
        is->n_ext = 0;
        /* reset */
        self->isassembled = false;
      }      
      else{    
        ext_data = self->ext_data; 
        ext_im = is->ext_im;
      }                 
    }
    else {
      self->isalloc = true;
      PARMS_NEW(aux_data);
      aux_data->n = lsize;
      PARMS_NEWARRAY(aux_data->space,   aux_data->n);
      PARMS_NEWARRAY0(aux_data->nnzrow, aux_data->n);
      PARMS_NEWARRAY(aux_data->pj,      aux_data->n);
      PARMS_NEWARRAY(aux_data->pa,      aux_data->n);
      for (i = 0; i < aux_data->n; i++) {
        aux_data->space[i] = 8;
        PARMS_NEWARRAY(aux_data->pj[i], aux_data->space[i]); 
        PARMS_NEWARRAY(aux_data->pa[i], aux_data->space[i]);      
      }
      self->aux_data = aux_data;
      /* create table to track off-diag column count */
      parms_TableCreate(&self->odtable,NULL,lsize);
/* Now allocate some data for the external data (contributions to 
 * other processors) from the current element matrix 
*/
      PARMS_NEW(ext_data);
      ext_data->n = lsize;
      PARMS_NEWARRAY(ext_data->space,   ext_data->n);
      PARMS_NEWARRAY0(ext_data->nnzrow, ext_data->n);
      PARMS_NEWARRAY(ext_data->pj,      ext_data->n);
      PARMS_NEWARRAY(ext_data->pa,      ext_data->n);
      PARMS_NEWARRAY(ext_im,      ext_data->n);
      for (i = 0; i < ext_data->n; i++) {
        ext_data->space[i] = 8;
        PARMS_NEWARRAY(ext_data->pj[i], ext_data->space[i]); 
        PARMS_NEWARRAY(ext_data->pa[i], ext_data->space[i]);      
      }
      self->ext_data = ext_data; 
      is->ext_im = ext_im;    
      
    }
  
    perm    = is->perm;
    numinf  = is->ninf;
/* Check if matrix is to be re-used with same non-zero pattern */
    if(self->isreset && (self->resetpattern == SAME_NONZERO_STRUCTURE))
    {
      for (i = 0; i < m; i++) {
      rindex = im[i] - offset;
      rp = parms_MapGlobalToLocal(is, rindex);
      /* this row resides the local processor */
      if (rp != NULL && *rp < lsize) {
	/* calculate the local row index */
	lrindex = perm[*rp];
	ra      = aux_data->pa[lrindex];
	rja     = aux_data->pj[lrindex];

	start = ia[i]-offset;
	end   = ia[i+1]-offset;
	for (j = start; j < end; j++) {
	  cindex = ja[j] - offset;
	  /* the local processor contains this column index */   
	  cp = parms_MapGlobalToLocal(is, cindex);
	  if (cp != NULL && *cp < lsize) {
	    /* local column index */
	    lcindex = perm[*cp];

	    for (k = 0; k < aux_data->nnzrow[lrindex]; k++) {
	      if (rja[k] == lcindex) {
		ra[k] = values[j];
	        break;
	      }
	     }
	  }
	  else {
	    if(cp != NULL){
	      lcindex = *cp;
	    
	      for (k = 0; k < aux_data->nnzrow[lrindex]; k++) {
	        if (rja[k] == lcindex) {
		  ra[k] = values[j];
		  break;
	        }
	       }
	      
	    }
	  }
	}
     }
    else{

/* This is an external contribution */
/* first check if there is enough memory allocated */
      space = ext_data->n;
      if(is->n_ext == space){
	  /* reallocate memory for holding new entry */
	  space += 8;
        PARMS_RESIZE(ext_data->space,   space);
        PARMS_RESIZE(ext_data->nnzrow, space);
        PARMS_RESIZE(ext_data->pj,      space);
        PARMS_RESIZE(ext_data->pa,      space);
        PARMS_RESIZE(is->ext_im,      space);
        for (k = ext_data->n; k < space; k++) {
          ext_data->space[k] = 8;
          ext_data->nnzrow[k] = 0;
          PARMS_NEWARRAY(ext_data->pj[k], ext_data->space[k]); 
          PARMS_NEWARRAY(ext_data->pa[k], ext_data->space[k]);      
        }
        ext_data->n = space;
        ext_im = is->ext_im;
      }

      rindex = 0;
      for(k = 0; k<is->n_ext; k++) /* check if row has been added already */
      {
        if(ext_im[k] == im[i])
            break;
      }
      if(k < is->n_ext) // already added row
        rindex = k;
      else          // new row contribution
      {
        rindex = is->n_ext;
        ext_im[rindex] = im[i];
      }
/* now get row info */            
	start = ia[i]-offset;
	end   = ia[i+1]-offset;
	pj = ext_data->pj[rindex];
      pa = ext_data->pa[rindex]; 
	ext_nnz = ext_data->nnzrow[rindex];
	cindex = 0;

	for (j = start; j < end; j++) {
        for(k = 0; k < ext_nnz; k++) /* check if column is available */
        { 
            if(pj[k] == ja[j])
              break;
        }
        if(k < ext_nnz){
          cindex = k;
          pa[cindex] += values[j];
        }
        else{
/* check if space is enough */     
	    space   = ext_data->space[rindex];
	    if (space == ext_data->nnzrow[rindex]) {
	      /* reallocate memory for holding new entry */
	      space += 8;
	      PARMS_RESIZE(ext_data->pa[rindex], space);
	      PARMS_RESIZE(ext_data->pj[rindex], space);
	      ext_data->space[rindex] = space;
	      pj = ext_data->pj[rindex];
            pa = ext_data->pa[rindex]; 
	    }
/* enter new column */
          pj[ext_data->nnzrow[rindex]] = ja[j];
          pa[ext_data->nnzrow[rindex]++] = values[j];
        }
      }
/* now update n_ext */
      if(rindex == is->n_ext)
        is->n_ext++;
    }      
   }
  }  
  else{            
    
     for (i = 0; i < m; i++) {
       rindex = im[i] - offset;
       rp = parms_MapGlobalToLocal(is, rindex);
       /* this row resides the local processor */
       if (rp != NULL && *rp < lsize) {
	 /* calculate the local row index */
	 lrindex = *rp;
	 ra      = aux_data->pa[lrindex];
	 rja     = aux_data->pj[lrindex];

	 start = ia[i]-offset;
	 end   = ia[i+1]-offset;
	 for (j = start; j < end; j++) {
	   cindex = ja[j] - offset;
	   /* the local processor contains this column index */   
	   cp = parms_MapGlobalToLocal(is, cindex);
	   if (cp != NULL && *cp < lsize) {
	     /* local column index */
	     lcindex = *cp;
	   }
	   else {
	     if (perm[lrindex] == -1) { /* nor marked yet */
	       perm[lrindex] = lsize-1-numinf;
	       numinf++;
	     }
	     lcindex = -cindex-1;
	    }
	    space   = aux_data->space[lrindex];
	    found = false;
	    for (k = 0; k < aux_data->nnzrow[lrindex]; k++) {
	      if (rja[k] == lcindex) {
	        found = true;
		ra[k] = values[j];
	        break;
	      }
	    }
	    if (found == false) { /* insert the new entry */
	    if (space == aux_data->nnzrow[lrindex]) {
	      /* reallocate memory for holding new entry */
	      space += 8;
	      PARMS_RESIZE(aux_data->pa[lrindex], space);
	      PARMS_RESIZE(aux_data->pj[lrindex], space);
	      aux_data->space[lrindex] = space;
	      ra  = aux_data->pa[lrindex];
	      rja = aux_data->pj[lrindex];
	    }
	    if (cp == NULL) {
	      size = parms_TableGetSize(is->table);
	      parms_TablePut(is->table, cindex, size);
           /* Initialize column count */
	      inc = 0;
	      parms_TablePut(self->odtable,cindex,inc);
	    }
	    else if((cp !=NULL) && (lcindex < 0)){
	      /* update off-diagonal column count */	    
	       odp = parms_TableGet(self->odtable, cindex);
             inc = *odp + 1;
	       parms_TablePut(self->odtable,cindex,inc);
	    }
	     
	    rja[aux_data->nnzrow[lrindex]]  = lcindex;
	    ra[aux_data->nnzrow[lrindex]++] = values[j];
	   }
	 }
      }
      else{

/* This is an external contribution */
/* first check if there is enough memory allocated */
      space = ext_data->n;
      if(is->n_ext == space){
	  /* reallocate memory for holding new entry */
	  space += 8;
        PARMS_RESIZE(ext_data->space,   space);
        PARMS_RESIZE(ext_data->nnzrow, space);
        PARMS_RESIZE(ext_data->pj,      space);
        PARMS_RESIZE(ext_data->pa,      space);
        PARMS_RESIZE(is->ext_im,      space);
        for (k = ext_data->n; k < space; k++) {
          ext_data->space[k] = 8;
          ext_data->nnzrow[k] = 0;
          PARMS_NEWARRAY(ext_data->pj[k], ext_data->space[k]); 
          PARMS_NEWARRAY(ext_data->pa[k], ext_data->space[k]);      
        }
        ext_data->n = space;
        ext_im = is->ext_im;
      }

      rindex = 0;
      for(k = 0; k<is->n_ext; k++) /* check if row has been added already */
      {
        if(ext_im[k] == im[i])
            break;
      }
      if(k < is->n_ext) // already added row
        rindex = k;
      else          // new row contribution
      {
        rindex = is->n_ext;
        ext_im[rindex] = im[i];
      }
/* now get row info */            
	start = ia[i]-offset;
	end   = ia[i+1]-offset;
	pj = ext_data->pj[rindex];
      pa = ext_data->pa[rindex]; 
	ext_nnz = ext_data->nnzrow[rindex];
	cindex = 0;

	for (j = start; j < end; j++) {
        for(k = 0; k < ext_nnz; k++) /* check if column is available */
        { 
            if(pj[k] == ja[j])
              break;
        }
        if(k < ext_nnz){
          cindex = k;
          pa[cindex] += values[j];
        }
        else{
/* check if space is enough */     
	    space   = ext_data->space[rindex];
	    if (space == ext_data->nnzrow[rindex]) {
	      /* reallocate memory for holding new entry */
	      space += 8;
	      PARMS_RESIZE(ext_data->pa[rindex], space);
	      PARMS_RESIZE(ext_data->pj[rindex], space);
	      ext_data->space[rindex] = space;
	      pj = ext_data->pj[rindex];
            pa = ext_data->pa[rindex]; 
	    }
/* enter new column */
          pj[ext_data->nnzrow[rindex]] = ja[j];
          pa[ext_data->nnzrow[rindex]++] = values[j];
        }
      }
/* now update n_ext */
      if(rindex == is->n_ext)
        is->n_ext++;
    }      
   }
   is->ninf = numinf;
  }
 }
  return 0;
}

/** 
 * Assembles the finite element matrix by updating off-processor 
 * contributions.
 * 
 * @param self    A parms_Mat object.                                           
 * @return 0 on success.
 */
int parms_MatAssembleElementMatrix(parms_Mat self)
{
  int i, j,mypid, npro, lsize, k, k1, k2, lnnz, gextnnz, n_ext; 
  int gnodv,index, nrow;
  int *clp, *info,*disp;     
  MPI_Comm comm;

  parms_Map is;
  parms_vcsr ext_data;
  int *nnzstart=NULL, *ptrvrecv=NULL,*ia=NULL, *ja=NULL;
  int *gextvars=NULL, *gextnnzrow=NULL;
  int *gextpja=NULL, *ext_nnzrow, *nnzptr, *ext_im;
  FLOAT *aa=NULL, *gextpa=NULL;

  /* extract information arrays */
  is                = self->is;
  npro              = is->npro;
  mypid             = is->pid;
  comm              = is->comm;
  lsize             = parms_MapGetLocalSize(is);

/* Check to see if there is any external data. If so, then update 
 * The local matrix before setting up the matrix.
*/ 
  ext_data = self->ext_data;
  n_ext = is->n_ext;
  ext_im = is->ext_im;

  PARMS_NEWARRAY(info, npro);
  PARMS_NEWARRAY(disp, npro+1);

  MPI_Allgather(&n_ext, 1, MPI_INT, info, 1, MPI_INT, comm);
  disp[0] = 0;
  for (i = 0; i < npro; i++) {
    disp[i+1] = disp[i] + info[i];
  }  
  /* total number of external variables */
  gnodv = disp[npro];

  if(gnodv) /* complete update of matrix values before matrix setup */
  {
    ext_nnzrow = ext_data->nnzrow;
    /* First gather data for all processors */
    PARMS_NEWARRAY(gextvars,  gnodv);
    PARMS_NEWARRAY0(gextnnzrow,  gnodv);
    MPI_Allgatherv(ext_im, n_ext, MPI_INT, gextvars,
		   info, disp, MPI_INT, comm);
    MPI_Allgatherv(ext_nnzrow, n_ext, MPI_INT, gextnnzrow,
		   info, disp, MPI_INT, comm);		   

    /* reshape pja and pa structs into 1D arrays */
    lnnz = 0;
    for(i = 0; i<n_ext; i++)
      lnnz += ext_nnzrow[i];
    
    PARMS_NEWARRAY(ja, lnnz);
    PARMS_NEWARRAY(aa, lnnz);
    k = 0;
    for(i = 0; i<n_ext; i++)
    {
      for(j = 0; j<ext_nnzrow[i]; j++)
      {
         ja[k] = ext_data->pj[i][j];
         aa[k++] = ext_data->pa[i][j];
       }
    }
    /* get global nnz count for external data */
    PARMS_NEWARRAY(nnzptr, npro+1);  
    MPI_Allgather(&lnnz, 1, MPI_INT, info, 1, MPI_INT, comm);
    nnzptr[0] = 0;
    for (i = 0; i < npro; i++) {
      nnzptr[i+1] = nnzptr[i] + info[i];
    }
    gextnnz = nnzptr[npro];           
    PARMS_NEWARRAY(gextpja, gextnnz);
    PARMS_NEWARRAY(gextpa, gextnnz);
 
    MPI_Allgatherv(ja, lnnz, MPI_INT, gextpja,
		   info, nnzptr, MPI_INT, comm);    
#if defined(DBL_CMPLX)       
    MPI_Allgatherv(aa, lnnz, MPI_CMPLX, gextpa,
		   info, nnzptr, MPI_CMPLX, comm); 
#else   
    MPI_Allgatherv(aa, lnnz, MPI_DOUBLE, gextpa,
		   info, nnzptr, MPI_DOUBLE, comm); 
#endif
    /* Now loop over external variables to get matrix contributions */
    PARMS_NEWARRAY0(nnzstart, gnodv);
    PARMS_NEWARRAY(ptrvrecv, gnodv);
    k = 0;
    k1 = 0;
    k2 = 0;
    for(i=0; i<npro; i++)
    {
      if(i != mypid)
      { 
        k2 = nnzptr[i];
	  for (j = disp[i]; j < disp[i+1]; j++) {
	    index = gextvars[j] - is->start;
	    clp = parms_MapGlobalToLocal(is, index);
          if (clp != NULL && *clp < lsize) {
            ptrvrecv[k] = j;
            nnzstart[k++] = k2;
            k1 += gextnnzrow[j];
          }
          k2 += gextnnzrow[j];
        }
       }
     }
     /* Allocate memory */
     nrow = k;
     PARMS_RESIZE(nnzstart, nrow);
     PARMS_NEWARRAY(ia, nrow+1);
     PARMS_FREE(ja);
     PARMS_FREE(aa); 
     
     ia[0] = is->start;
     /* Now loop again to actually set up csr matrix */
    for(i=0; i<nrow; i++)
    {
        ia[1] = ia[0] + gextnnzrow[ptrvrecv[i]];
        parms_MatSetValues(self, 1, &gextvars[ptrvrecv[i]], ia, &gextpja[nnzstart[i]], &gextpa[nnzstart[i]], ADD);
    }
/* free allocated memory */
    PARMS_FREE(ptrvrecv);
    PARMS_FREE(nnzptr);
    PARMS_FREE(nnzstart);
    PARMS_FREE(gextvars);
    PARMS_FREE(gextnnzrow);
    PARMS_FREE(ia);
    PARMS_FREE(gextpja);
    PARMS_FREE(gextpa);    
  }
    PARMS_FREE(disp);
    PARMS_FREE(info);
/* free allocated memory for external data */
    if(ext_data){
     PARMS_FREE(ext_data->nnzrow);
     PARMS_FREE(ext_data->space);      
     for(i=0;i<ext_data->n; i++)
     {
        PARMS_FREE(ext_data->pa[i]);
        PARMS_FREE(ext_data->pj[i]);
     }
     PARMS_FREE(ext_data->pa);
     PARMS_FREE(ext_data->pj);
     PARMS_FREE(ext_data);     
     PARMS_FREE(ext_im);
    }
/* Done with matrix update */
  self->isassembled = true;
  
  return 0;
}		       

/** 
 * Insert values to the parms_Mat object self. This assumes matrix 
 * values for this row have already been set, and are to be replaced 
 * by the new ones provided as input.
 * 
 * @param self    A parms_Mat object.                                
 * @param m 	The number of rows inserted.      
 * @param im 	An array of row global indices.                                                 
 * @param ia 	An array of pointer to the beginning of each row in
 *                array ja.    
 * @param ja      An array of column global indices.     
 * @param values  An array of values.                      
 * @return 0 on success.
 */
int parms_MatResetRowValues(parms_Mat self, int m, int *im, int *ia,
		       int *ja, FLOAT *values)
{
  int i, j, k, rindex,lrindex, lsize, offset, n_ext, gnodv;
  int npro,myid,gm, cindex;
  int *rp,*ext_im,*disp=NULL,*info=NULL,*gvars=NULL,*rowptr=NULL;
  int *cp, *odp, *pj, nnz, dec;
  parms_Map is;
  parms_vcsr aux_data, ext_data;
  BOOL isserial, isalloc, isassembled;  
  MPI_Comm comm;  
  
  isalloc  = self->isalloc;
  is       = self->is;
  offset    = is->start;
  lsize    = parms_MapGetLocalSize(is);
  isserial = is->isserial;
  n_ext    = is->n_ext; 
  ext_im   = is->ext_im; 
  npro     = is->npro;
  myid     = is->pid;
  comm     = is->comm;
  isassembled = self->isassembled;
  
/* Trivial case */  
  if((isserial)||(!isalloc)){
    return parms_MatSetValues(self, m, im, ia, ja, values, INSERT);
  }    
  
/* Take care of external contributions if any. Check to see 
 * if the external variables have been assembled prior to this 
 * function call
*/
  PARMS_NEWARRAY(disp, npro+1);
  PARMS_NEWARRAY(info, npro);

  MPI_Allgather(&n_ext, 1, MPI_INT, info, 1, MPI_INT, comm);
  disp[0] = 0;
  for (i = 0; i < npro; i++) {
    disp[i+1] = disp[i] + info[i];
  }  
  /* total number of external variables */
  gnodv = disp[npro];
  /* Check ... */
  if((gnodv) && (!isassembled))
  {
     ext_data = self->ext_data;
    /* First gather data for all processors */
    PARMS_NEWARRAY(rowptr, npro+1);  
    MPI_Allgather(&m, 1, MPI_INT, info, 1, MPI_INT, comm);
    rowptr[0] = 0;
    for (i = 0; i < npro; i++) {
      rowptr[i+1] = rowptr[i] + info[i];
    }
    gm = rowptr[npro];          
        
    PARMS_NEWARRAY(gvars,  gm);
    MPI_Allgatherv(im, m, MPI_INT, gvars,
		   info, rowptr, MPI_INT, comm);
    
  /* now loop over external variables to remove their contribution */
    for(i=0; i<npro; i++)
    {
      if(i != myid)
      {
        for(j = rowptr[i]; j<rowptr[i+1]; j++)
        {
          for(k=0; k<n_ext; k++)
          {
             if(ext_im[k] == gvars[j]){
               ext_data->nnzrow[k] = 0;
               break;
             }
          }
        }
      }
    } 
    PARMS_FREE(gvars);   
    PARMS_FREE(rowptr);   
  }
  PARMS_FREE(disp);
  PARMS_FREE(info);
/* Done with external contributions. Now reset the row */

/* Check if matrix is to be re-used with same non-zero pattern */
    if(self->isreset && (self->resetpattern == SAME_NONZERO_STRUCTURE))
    {
      parms_MatSetValues(self, m, im, ia, ja, values, INSERT);    
    }
    else  
    {
      aux_data = self->aux_data;
     /* begin main loop over rows to be reset */
     for(i = 0; i<m; i++)
     {
       rindex = im[i] - offset;
       rp = parms_MapGlobalToLocal(is, rindex);    
       /* this row resides the local processor */
       if (rp != NULL && *rp < lsize) {
       /* calculate the local row index */
	   lrindex = *rp;

       /* Assume row is an internal variable - not an interface variable */
         if(is->perm[lrindex] != -1){
          is->perm[lrindex] = -1;
          is->ninf = is->ninf - 1;
          
       /* loop over off-diagonal nodes to see if any needs to be removed from table */
          nnz = aux_data->nnzrow[lrindex];
          pj = aux_data->pj[lrindex];
          for(j = 0; j<nnz; j++)
          {
             /* consider off-diagonal columns only */
             if(pj[j] >= 0) continue;
               
             cindex = -pj[j] - 1;
             cp = parms_TableGet(is->table, cindex);
           
             if(cp != NULL && *cp >= lsize)
             {
               odp = parms_TableGet(self->odtable, cindex);
               if(*odp == 0)
               {
                    parms_TableRemoveFromLast(is->table, cindex);
               }
               else
               {
             /* Decrement column count and continue */
                    dec = *odp - 1;
                    parms_TablePut(self->odtable,cindex,dec);             
               }
             }
           }      
         }
       /* Reset nnz for current row to overwrite data */
         aux_data->nnzrow[lrindex] = 0;
        /* Now call matsetvalues to insert values into row */
//         parms_MatSetValues(self, 1, &im[i], &ia[i], ja, values, INSERT);   
        }
      }
      /* Now call matsetvalues to insert values into rows */
      parms_MatSetValues(self, m, im, ia, ja, values, INSERT);            
     }
     
  return 0;
}

/** 
 * Reset the matrix to be re-used. 
 * @param self    A parms_Mat object.                                
 * @param nonzerostructure  The nonzero structure:
 *                        SAME_NONZERO_STRUCTURE
 *                        DIFFERENT_NONZERO_STRUCTURE     
 *                     
 * @return 0 on success.
 */
int parms_MatReset(parms_Mat self, NNZSTRUCT nonzerostructure)
{
  parms_Map  is;
  parms_vcsr aux_data;
  int i, j;
  int index;
  int lsize, nnz;
  FLOAT *ra;

  is       = self->is;
  lsize    = parms_MapGetLocalSize(is);
  aux_data = self->aux_data;
  
/*********************/ 
    if(self->issetup){
      
      if(nonzerostructure == SAME_NONZERO_STRUCTURE)
      {
        self->resetpattern = SAME_NONZERO_STRUCTURE;
        /* zero out previous entries */
        for(i=0; i<lsize; i++)
        {
            ra = aux_data->pa[i];
            nnz = aux_data->nnzrow[i];
            for(j=0; j<nnz; j++)
            {
              ra[j] = 0.0;
            }
        }
      }
      else
      {  
         self->resetpattern = DIFFERENT_NONZERO_STRUCTURE;      
        /* Matrix structure is being re-used */
         PARMS_NEWARRAY(aux_data->space,   aux_data->n);
        /* Initially same non-zero structure */      
        for (i = 0; i < aux_data->n; i++) {
          aux_data->space[i] = aux_data->nnzrow[i];
        /* Now zero out entries by simply setting nnzrow[i] = 0 */
          aux_data->nnzrow[i] = 0;      
        }
        /* reset permutation array - perm */   
        for (i = 0; i < lsize; i++) {
           is->perm[i] = -1;
        } 
        /* reset number of interface variables to zero */
        is->ninf = 0;
        /* reset number of external variables to zero (if any) */
        is->n_ext = 0; 
        
        /* create new table for offdiagonal column count */
        parms_TableCreate(&self->odtable, NULL, lsize);        
        /* free hash table associated with previous data, and create new one */
        parms_TableFree(&is->table);
        /* create new hash table for new data */
        parms_TableCreate(&is->table, NULL, lsize);
        for(index = 0; index < lsize; index++)
        {
	    parms_TablePut(is->table, is->lvars[index], index);
	  }
       }
	 self->issetup = false;
       self->isreset = true;
     }

  return 0;
}


/** 
 * Get the diagonal part of the local matrix. 
 * 
 * @param self A parms_Mat object.
 * @param mat  The diagonal part of the local matrix.    
 * 
 * @return 0 on success.
 */
int parms_MatGetDiag(parms_Mat self, void **mat)
{
  return self->ops->getdiag(self, mat);
}
/** 
 * Get the diagonal part of the local matrix. 
 * 
 * @param self A parms_Mat object.
 * @param mat  The diagonal part of the local matrix.    
 * 
 * @return 0 on success.
 */
int parms_MatGetOffDiag(parms_Mat self, void **mat)
{
  return self->ops->getoffdiag(self, mat);
}

/** 
 * Perform \f$z = alpha*self*x + beta*y\f$.
 * 
 * @param self   A matrix object.           
 * @param alpha  A scalar.                  
 * @param x      A vector object.           
 * @param beta   A scalar.                  
 * @param y      A vector object.           
 * @param z      A vector stores the result.
 * 
 * @return 0 on success.
 */
int parms_MatMVPY(parms_Mat self, FLOAT alpha, FLOAT *x, FLOAT beta,
		  FLOAT *y, FLOAT *z)
{
  return self->ops->mvpy(self, alpha, x, beta, y, z);
}

/** 
 * Perform the multiplication of the off-diagonal matrix and the
 * external vars. 
 *
 * The local matrix can be written as follows:
 * 
 *  \f[
 *  \left(
 *  \begin{array}{ccc}
 *    B   &   E & 0\\
 *    F   &   C & M_{ext}
 *  \end{array}
 *  \right)
 *  \f],
 *  where \f$\left(\begin{array}{cc}
 *    B  &   E \\
 *    F  &   C 
 *  \end{array}
 *  \right)\f$ corresponds to the variables on the local PE. This
 *  function performs 
 *  \f[
 *    y[pos..n] = M_{ext} \times x_{ext}
 *  \f]
 *
 * @param self A matrix object.                                      
 * @param x    A vector object.                                      
 * @param y    A vector object.                                      
 * @param pos  The offset of x from the beginning of the local
 *             vector. 
 * 
 * @return 0 on success.
 */
int parms_MatVecOffDiag(parms_Mat self, FLOAT *x, FLOAT *y, int pos)
{
  return self->ops->mvoffd(self, x, y, pos);
}

/** 
 * Free the memory for the submatrix.
 * 
 * @param self A parms_Mat object.             
 * @param mat  The submatrix to be freed.      
 * 
 * @return 0 on success.
 */
int parms_MatFreeSubMat(parms_Mat self, void *mat)
{
  return self->ops->matfree(self, mat);
}

/** 
 * Extend submatrix by including equations correspond to the
 * immediate neighbouring variables.
 * 
 * @param self     A matrix object.                                         
 * @param handler  A communication handler.                                 
 * @param start    The beginning location of mat in the local matrix.       
 * @param mat      The submatrix to be extended.        
 * @param n 	 The size of extended matrix returned.
 * @param ext_mat  The extended matrix created.            
 * 
 * @return 0 on success.
 */
int parms_MatExtend(parms_Mat self, parms_Comm handler, int start,
		    void *mat, int *n, void **ext_mat)
{

  return self->ops->extend(self, handler, start, mat, n, ext_mat);
}

/** 
 * Get the local matrix. 
 * 
 * @param self A matrix object.                               
 * @param mat  The submatrix returned in a specific format.   
 * 
 * @return 0 on success.
 */
int parms_MatGetSubMat(parms_Mat self, void **mat)
{
  return self->ops->getlmat(self, mat);
}

/** 
 * Get the communication handler.
 * 
 * @param self    A matrix object.                     
 * @param handler The communication handler returned.  
 * 
 * @return 0 on success.
 */
int parms_MatGetCommHandler(parms_Mat self, parms_Comm *handler)
{
  return self->ops->gethandler(self, handler);
}
