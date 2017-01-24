/*--------------------------------------------------------------------
  parms_MapCreateFromLocal  : create a local parms_Map object.
  parms_MapCreateFromGlobal : create a parms_Map object based on the 
                              output of Metis.
  parms_MapCreateFromDist   : create a parms_Map object based on the
                              output of ParMetis.
  parms_MapCreateFromPetsc  : create a parms_Map object based on the
                              default partitoning strategy in PETSc
                              (contiguous chunks of rows are
                              distributed among PEs)   
  parms_MapCreateFromPtr    : create a parms_Map object based on a
                              list of all vars and a pointer array to
                              the beginning location of each PE in the
                              list.   
  parms_MapGetGlobalSize    : get global size of data partitioned.
  parms_MapGetLocalSize     : get the size of data on the local PE.
  parms_MapGetNumProcs      : get the number of PEs.
  parms_MapGetPid           : get id of PE.
  parms_MapGlobalToLocal    : get the local index of a var with its
                              global index as the input. 
  parms_MapView             : dump the parms_Map object.
  parms_MapFree             : free the memory for the parms_Map object.

  A code fragment for using parms_Map functions:

  parms_Map map;
  parms_Mat mat;
  
  // partition the graph or mesh using Metis 
  METIS_PartGraphVKway(&n, ..., riord);

  // create a parms_Map object based on the output of Metis
  parms_MapCreateFromGlobal(&map, n, riord, ...);

  // create vector and matrix object based on map
  parms_MatCreate(&mat, map);

  $Id: parms_map.c,v 1.1.1.1 2006-11-27 22:28:01 zzli Exp $
--------------------------------------------------------------------*/

#include "include/parms_map_impl.h"

/** 
 * Free the memory for the parms_Map object pointed to by self.
 * 
 * @param self A pointer to the parms_Map object.
 * 
 * @return 0 on success.
 */
int parms_MapFree(parms_Map *self)
{
   (*self)->ref--;
  /* if the reference number is zero, then remove the object */
  if ((*self)->ref == 0 ) {
    if (!(*self)->isserial) {
      if((*self)->lvars) PARMS_FREE((*self)->lvars);
      if((*self)->vsend) PARMS_FREE((*self)->vsend);
      if((*self)->vstable)parms_TableFree(&(*self)->vstable);             
      parms_TableFree(&(*self)->table);
      MPI_Comm_free(&(*self)->comm);
    }

    if ((*self)->ispermalloc) {
      PARMS_FREE((*self)->perm);
      PARMS_FREE((*self)->iperm);
    }
    PARMS_FREE(*self);
  }
  return 0;
}

/** 
 * Dump the parms_Map object.
 * 
 * @param self A parms_Map object.       
 * @param v    A parms_Viewer object.    
 * 
 * @return 0 on success.
 */
int parms_MapView(parms_Map self, parms_Viewer v)
{
  FILE *fp;
  int i, *val;

  parms_ViewerGetFP(v, &fp);
  fprintf(fp, "=====Dump the content of parms_Map=====\n");
  fprintf(fp, "pid = %d, npro = %d\n", self->pid, self->npro);
  fprintf(fp, "lsize = %d, gsize = %d\n", self->lsize, self->gsize);
  fprintf(fp, "isserial = %d, isperm = %d, isalloc = %d\n",
	  self->isserial, self->isperm, self->ispermalloc);
  if (!self->isserial) {
    fprintf(fp, "size in the table is %d\n",
	    parms_TableGetSize(self->table)); 
    fprintf(fp, "global vertices ==> local vertices on proc%d\n", self->pid);
    for (i = 0; i < self->lsize; i++) {
      val = parms_TableGet(self->table, self->lvars[i]);
      fprintf(fp, "%d==>%d\n", self->lvars[i], *val);
    }
    if (self->isperm) {
      fprintf(fp, "The permutation array is:\n");
      for (i = 0; i < self->lsize; i++) {
	fprintf(fp, "perm[%d]=%d\n", i, self->perm[i]);
      }
    }
  }

  fprintf(fp, "nint = %d\n", self->nint);
  fprintf(fp, "ninf = %d\n", self->ninf);

  fprintf(fp, "ninf_send = %d\n", self->ninf_send);
  for (i = 0; i < self->ninf_send; i++) {
    fprintf(fp, "send(%d)=%d  ", i, self->vsend[i]);
    val = parms_TableGet(self->vstable, self->vsend[i]);
    fprintf(fp, "vsend[%d]=>%d\n", self->vsend[i], *val);
  }

  
  fprintf(fp, "=======================================\n");
  parms_ViewerStoreFP(v, fp);
  return 0;
}

/** 
 * Get the size of variables on the local PE.
 *
 * Return the size of variables not vertices on the local PE.
 * 
 * @param self A parms_Map object.
 * 
 * @return The size of variables on the local PE.
 */
int parms_MapGetLocalSize(parms_Map self)
{
  return self->lsize;
}

/** 
 * Get the global size of variables.
 *
 * Return the global size of variables rather than vertices.
 * 
 * @param self A parms_Map object.
 * 
 * @return The global size of variables rather than vertices.
 */
int parms_MapGetGlobalSize(parms_Map self)
{
  return self->gsize;
}

/** 
 * Get PE's ID.
 * 
 * @param self A parms_Map object.
 * 
 * @return The PE's ID.
 */
int parms_MapGetPid(parms_Map self)
{
  /*-------------------- obtains processor IDs */
  return self->pid;
}

/** 
 * Get the number of PEs.
 * 
 * @param self A parms_Map object.
 * 
 * @return The number of PEs.
 */
int parms_MapGetNumProcs(parms_Map self)
{
  /*-------------------- Get number of PEs */

  return self->npro;
}

/** 
 * Create a local parms_Map object.
 * 
 * @param self   A parms_Map object created.         
 * @param gsize  The global size of variables.       
 * @param offset The start index. 
 *               - 1 FORTRAN
 *               - 0 C
 *
 * @return 0 on success.
 */
int parms_MapCreateFromLocal(parms_Map *self, int gsize, int offset)  
{
  parms_Map newMap;

  PARMS_NEW0((newMap));
  newMap->ref = 1;
  newMap->gsize = gsize;
  newMap->lsize = gsize;
  newMap->start = offset;
  newMap->isserial  = true;
  newMap->isperm    = false;
  newMap->isvecperm    = false;
  MPI_Comm_rank(MPI_COMM_WORLD, &newMap->pid);
  newMap->npro = 1;
  newMap->comm = 0;
  newMap->nint = gsize;
  newMap->ninf = 0;
  newMap->n_ext = 0;  

  newMap->ispermalloc = false;
  newMap->isdatalloc = false;  
  
  *self = newMap;
  
/* Define complex data type for MPI if complex code is compiled */
#if defined(DBL_CMPLX)
  parms_InitComplex();
#endif 
  return 0;
}

/** 
 * Create a parms_Map object based on the Metis partitioning.
 * 
 * @param self     A pointer to the parms_Map object created.       
 * @param gsize    The total number of vertices.                    
 * @param npar 	 An integer array of size gsize. node \f$i\f$     
 *       	       resides on PE npar[i].                           
 * @param comm     MPI communicator.                                
 * @param offset   The start index.                                 
 *                 - 1 FORTRAN
 *                 - 0 C
 * @param dof      The number of variables associated with each
 *                 vertex.   
 * @param VARSTYPE Assuming the variables \f$u_i, v_i\f$ are
 *                 associated with vertex \f$i\f$, two styles of
 *                 numbering variables are as follows:
 *                 - INTERLACED. Variables are numbered in the
 *                   order of \f$u_1, v_1, u_2, v_2, \cdots\f$; 
 *                 - NONINTERLACED. Variables are numbered in the
 *                   order of \f$u_1, u_2, u_3,...,v_1, v_2,...\f$. 
 *
 * @return 0 on success.
 */
int parms_MapCreateFromGlobal(parms_Map *self, int gsize, int *npar,
			      MPI_Comm comm, int offset, int dof,
			      VARSTYPE vtype)
{
  parms_Map newMap;
  int pid, i, j, gindex;
  int hint, gv_size, lsize;
  MPI_Comm newComm;

  MPI_Comm_dup(comm, &newComm);
  PARMS_NEW0((newMap));
  newMap->ref = 1;
  MPI_Comm_rank(newComm, &newMap->pid);
  MPI_Comm_size(newComm, &newMap->npro);
  newMap->comm = newComm;
  pid = newMap->pid;
  gv_size = newMap->gsize = gsize * dof;
  newMap->start = offset;
  newMap->dof = dof;
  newMap->vtype = vtype;
  newMap->isserial =  false;
  if (newMap->npro == 1) {
    newMap->isserial = true;
  }
  newMap->isperm   =  false;
  newMap->isvecperm    = false;
  newMap->ispermalloc = false;
  newMap->isdatalloc = false;  
  
  if (!newMap->isserial) {
    hint = gv_size / newMap->npro;
    /* create a hash table */
    parms_TableCreate(&newMap->table, NULL, hint);
    PARMS_NEWARRAY(newMap->lvars, hint);
    lsize = 0;
    if (vtype == INTERLACED) {
      for (i = 0; i < gsize; i++) {
	if (npar[i] - offset == pid) {
	  for (j = 0; j < dof; j++) {
	    gindex = dof*i+j;
	    /* put the pair (gindex, lsize) into the table */
	    parms_TablePut(newMap->table, gindex, lsize);
	    if (lsize >= hint) {
	      hint += 50;
	      PARMS_RESIZE(newMap->lvars, hint);
	    }
	    newMap->lvars[lsize++] = gindex;
	  }
	}
      }
    }
    else if (vtype == NONINTERLACED) {
      for (i = 0; i < gsize; i++) {
	if (npar[i] - offset == pid) {
	  for (j = 0; j < dof; j++) {
	    gindex = gsize*j+i;
	    /* put the pair (gindex, lsize) into the table */
	    parms_TablePut(newMap->table, gindex, lsize);
	    if (lsize >= hint) {
	      hint += 50;
	      PARMS_RESIZE(newMap->lvars, hint);
	    }
	    newMap->lvars[lsize++] = gindex;
	  }
	}
      }
    }

    newMap->lsize = lsize;
    PARMS_RESIZE(newMap->lvars, lsize);
    newMap->ispermalloc = true;
    PARMS_NEWARRAY0(newMap->perm,  lsize);
    PARMS_NEWARRAY0(newMap->iperm, lsize);
    for (i = 0; i < lsize; i++) {
      newMap->perm[i] = -1;
    }
  }
  else {
    lsize = gv_size;
    newMap->lsize = gv_size;
  }
  newMap->nint     =  lsize;
  newMap->ninf     =  0;
  newMap->n_ext = 0;  
  
  *self = newMap;
  
/* Define complex data type for MPI if complex code is compiled */
#if defined(DBL_CMPLX)
  parms_InitComplex();
#endif 
  
  return 0;
}

/** 
 * Create a parms_Map object based on the output of ParMetis.
 * 
 * @param self     A parms_Map object created.
 * @param vtxdist  An integer array of size np+1, where np is the
 *                 number of PEs. This array indicates the range of
 *                 vertices that are local to each processor.  PE i
 *                 stores vertices in the range of [vtxdist[i],
 *                 vtxdist[i+1]).
 * @param part     An array of size equal to the number of
 *                 locally-stored vertices. part[j] indicates the ID
 *                 of the PE to which the vertex with local index j
 *                 and global index vtxdist[pid]+j belongs (pid is ID
 *                 of local PE). 
 * @param comm     MPI communicator.
 * @param offset   The start index. 
 *                 - 1 FORTRAN
 *                 - 0 C
 * @param dof      The number of variables associated with each
 *                  vertex.  
 * @param vtype    Assuming the variables u_i, v_i are associated
 *                 with vertex i, two styles of numbering variables
 *                 are as follows: 
 *                 - INTERLACED. Variables are numbered in the
 *                   order of \f$u_1, v_1, u_2, v_2, \cdots\f$; 
 *                 - NONINTERLACED. Variables are numbered in the
 *                   order of \f$u_1, u_2, u_3,...,v_1, v_2,...\f$. 
 * 
 * @return 0 on success.
 */
int parms_MapCreateFromDist(parms_Map *self, int *vtxdist, int *part,
			    MPI_Comm comm, int offset, int dof,
			    VARSTYPE vtype)  
{
  parms_Map newMap;
  int npro, pid, i, j, l, gsize;
  int gv_size, lsize, nl, ind, gindex;
  int *num, *nums, *num_rcv, *disp, *snd_buf, *rcv_buf;
  MPI_Comm newComm;

  MPI_Comm_dup(comm, &newComm);
  PARMS_NEW0((newMap));
  newMap->ref = 1;
  MPI_Comm_rank(newComm, &newMap->pid);
  MPI_Comm_size(newComm, &newMap->npro);
  newMap->comm = newComm;
  npro = newMap->npro;
  pid  = newMap->pid;
  /* get the number of local vertices */
  nl =  vtxdist[pid+1] - vtxdist[pid];
  /* calculate the total number of vertices */
  gsize = vtxdist[npro] - vtxdist[0];
  /* total number of variables  */
  gv_size = newMap->gsize = gsize * dof;
  newMap->start = offset;
  newMap->dof = dof;
  newMap->vtype = vtype;
  newMap->isserial =  false;
  if (newMap->npro == 1) {
    newMap->isserial = true;
  }
  newMap->isperm =  false;
  newMap->isvecperm    = false;
  newMap->ispermalloc = false;
  newMap->isdatalloc = false;  
  
  if (!newMap->isserial) {
    PARMS_NEWARRAY(snd_buf, nl);
    /* create a hash table */
    parms_TableCreate(&newMap->table, NULL, nl);
    PARMS_NEWARRAY0(num, npro);
    PARMS_NEWARRAY(nums, npro);
    /* num[i] stores the number of locally-stored variables being  
       distributed to PE i */
    for (i = 0; i < nl; i++) {
      num[part[i]-offset]++;
    }
    MPI_Allreduce(num, nums, npro, MPI_INT, MPI_SUM, newComm); 
    /* nums[i] stores the number of variables on PE i */
    lsize = newMap->lsize = nums[pid]*dof;
    PARMS_NEWARRAY(newMap->lvars, lsize);
    PARMS_FREE(nums);

    PARMS_NEWARRAY(disp,     npro+1);
    PARMS_NEWARRAY(num_rcv,  npro);
    /* num_rcv stores the number of data received from other processors */
    for (i = 0; i < npro; i++) {
      MPI_Gather(&num[i], 1, MPI_INT, num_rcv, 1, MPI_INT, i,
		 newComm); 
      /* snd_buf stores the data sent to PE i */
      ind = 0;
      for (j = 0; j < nl; j++) {
	if (part[j]-offset == i) {
	  snd_buf[ind++] = vtxdist[pid] + j - offset;
	}
      }
      if (pid == i) {
	/* disp is an integer array.  disp[i]  specifies  the
	   displacement relative to rcv_buf  at which to place the  
	   incoming data from PE i  */
	disp[0] = 0;
	for (j = 0; j < npro; j++) {
	  disp[j+1] = disp[j] + num_rcv[j];
	}

	/* variables in rcv_buf are stored in C-style */
	PARMS_NEWARRAY(rcv_buf, disp[npro]);
	MPI_Gatherv(snd_buf, num[i], MPI_INT, newMap->lvars, num_rcv,
		    disp, MPI_INT, i, newComm);
	if (vtype == INTERLACED) {
	  ind = 0;
	  for (j = 0; j < disp[npro]; j++) {
	    for (l = 0; l < dof; l++) {
	      gindex = dof*rcv_buf[j]+l;
	      parms_TablePut(newMap->table, gindex, ind);
	      newMap->lvars[ind++] = gindex;
	    }
	  }
	}
	else if (vtype == NONINTERLACED) {
	  for (j = 0; j < disp[npro]; j++) {
	    parms_TablePut(newMap->table, newMap->lvars[j], j);
	  }
	  ind = disp[npro];
	  for (j = 0; j < disp[npro]; j++) {
	    for (l = 1; l < dof; l++) {
	      gindex = gsize*l+newMap->lvars[j];
	      parms_TablePut(newMap->table, gindex, ind);
	      newMap->lvars[ind++] = gindex;
	    }
	  }
	}
	PARMS_FREE(rcv_buf);
      }
      else {
	MPI_Gatherv(snd_buf, num[i], MPI_INT, rcv_buf, num_rcv,
		    disp, MPI_INT, i, newComm); 
      }
    }

    PARMS_FREE(snd_buf);
    PARMS_FREE(num);
    PARMS_FREE(num_rcv);
    PARMS_FREE(disp);
    newMap->ispermalloc = true;
    PARMS_NEWARRAY0(newMap->perm,  lsize);
    PARMS_NEWARRAY0(newMap->iperm, lsize);
    for (i = 0; i < lsize; i++) {
      newMap->perm[i] = -1;
    }
  }
  else {
    lsize = gv_size;
    newMap->lsize = gv_size;
  }
  newMap->nint    =  lsize;
  newMap->ninf    =  0;
  newMap->n_ext = 0;  

  *self = newMap;

/* Define complex data type for MPI if complex code is compiled */
#if defined(DBL_CMPLX)
  parms_InitComplex();
#endif 

  return 0;
}

/** 
 * Create a parms_Map object.
 * 
 * @param self     A parms_Map object created.
 * @param gsize    The total number of vertices.       
 * @param nodes    A list of all vertices stored PE by PE.        
 * @param p2nodes  An integer array of size np+1, np is the number of
 *                 PEs. If k1 = p2nodes[i], k2 = p2nodes[i+1]
 *                 , then PE i contains the vertices in the 
 *                 range of [[nodes[k1], [nodes[k2-1]]. 
 * @param comm     MPI communication.
 * @param dof      The number of variables associated with each
 *                 node. 
 * @param vtype    Assuming the variables \f$u_i, v_i\f$ are
 *                 associated with vertex \f$i\f$, two styles of
 *                 numbering variables are as follows:
 *                 - INTERLACED. Variables are numbered in the
 *                   order of \f$u_1, v_1, u_2, v_2, \cdots\f$; 
 *                 - NONINTERLACED. Variables are numbered in the
 *                   order of \f$u_1, u_2, u_3,...,v_1, v_2,...\f$. 
 * 
 * @return 0 on success.
 */
int parms_MapCreateFromPtr(parms_Map *self, int gsize, int *nodes, int
			   *p2nodes, MPI_Comm comm, int dof, VARSTYPE
			   vtype)    
{
  int pid, npro, index, i, j, l;
  int start_index, end_index, gindex;
  int offset, start_offset;
  parms_Map newMap;
  MPI_Comm newComm;

  MPI_Comm_dup(comm, &newComm);

  PARMS_NEW0((newMap));
  newMap->ref = 1;
  newMap->gsize = gsize * dof;
  start_offset = newMap->start = p2nodes[0];
  newMap->dof = dof;
  newMap->vtype = vtype;
  newMap->isserial = false;
  newMap->isperm   = false;
  newMap->isvecperm    = false;
  newMap->ispermalloc = false;
  newMap->isdatalloc = false;  
  newMap->nint      = 0;
  MPI_Comm_rank(newComm, &newMap->pid);
  MPI_Comm_size(newComm, &newMap->npro);
  newMap->comm = newComm;
  npro = newMap->npro;
  pid  = newMap->pid;

  if (npro == 1) {
    newMap->isserial = true;
  }

  newMap->lsize  = p2nodes[pid+1] - p2nodes[pid];
  newMap->lsize *= dof;
  

  if (!newMap->isserial) {
    PARMS_NEWARRAY(newMap->lvars, newMap->lsize);

    index = 0;
    /* create a table */
    parms_TableCreate(&newMap->table, NULL, newMap->lsize);

    if (vtype == INTERLACED) {
      /* calculate the range [start_index, end_index) */
      start_index = p2nodes[pid]   - start_offset;
      end_index   = p2nodes[pid+1] - start_offset;
      for (j = start_index; j < end_index; j++) {
	offset = j - start_index;
	for (l = 0; l < dof; l++) {
	  /* global index of variable l on vertex nodes[j] */
	  gindex = dof*(nodes[j]-start_offset)+l;
	  parms_TablePut(newMap->table, gindex, index);
	  /* local variables in global indices */
	  newMap->lvars[dof*offset+l] = gindex;

	  index++;
	}
      }
    }
    else if (vtype == NONINTERLACED) {
      /* calculate the range [start_index, end_index) */
      start_index = p2nodes[pid]   - start_offset;
      end_index   = p2nodes[pid+1] - start_offset;
          
      for (j = start_index; j < end_index; j++) {
	offset = j - start_index;
	for (l = 0; l < dof; l++) {
	  /* global index of variable l on vertex nodes[j] */
	  gindex = gsize*l+nodes[j]-start_offset;
	  parms_TablePut(newMap->table, gindex, index);
	  /* local variables in global indices */
	  newMap->lvars[dof*offset+l] = gindex;
  
	  index++;
	}
      }
    }
    newMap->ispermalloc = true;
    PARMS_NEWARRAY0(newMap->perm,  newMap->lsize);
    PARMS_NEWARRAY0(newMap->iperm, newMap->lsize);
    for (i = 0; i < newMap->lsize; i++) {
      newMap->perm[i] = -1;
    }
  }
  newMap->nint    =  newMap->lsize;
  newMap->ninf    =  0;
  newMap->n_ext = 0;  

  *self = newMap;

/* Define complex data type and ops for MPI if complex code is compiled */
#if defined(DBL_CMPLX)
  parms_InitComplex();
#endif 

  return 0;  
}

/** 
 * Create a parms_Map object with the default partioning strategy in
 * PETSc. 
 * 
 * @param self  A parms_Map object created.  
 * @param m     The local size of variables. 
 * @param M     The global size of variables.
 * @param comm  MPI communicatior.           
 * 
 * @return 0 on success.
 */
int parms_MapCreateFromPetsc(parms_Map *self, int m, int M, MPI_Comm
			     comm) 
{
  int pid, npro, index, i, j;
  int start_index, end_index;
  int *nptr;
  parms_Map newMap;
  MPI_Comm newComm;

  MPI_Comm_dup(comm, &newComm);
  /* create a parms_Map object */

  PARMS_NEW0((newMap));
  newMap->ref = 1;
  newMap->gsize = M;
  newMap->start = 0;
  newMap->dof   = 1;
  newMap->vtype = NONINTERLACED;
  newMap->isserial = false;
  newMap->isperm   = false;
  newMap->isvecperm    = false;
  newMap->ispermalloc = false;
  newMap->isdatalloc = false;    
  newMap->nint     = m;
  newMap->ninf     = 0;
  newMap->n_ext = 0;    
  MPI_Comm_rank(newComm, &newMap->pid);
  MPI_Comm_size(newComm, &newMap->npro);
  newMap->comm = newComm;
  npro = newMap->npro;
  pid  = newMap->pid;

  if (npro == 1) {
    newMap->isserial = true;
  }

  newMap->lsize = m;
  
  /* if the number of PEs is more than 1 */
  if (!newMap->isserial) {
    PARMS_NEWARRAY(newMap->lvars, newMap->lsize);
    PARMS_NEWARRAY(nptr, npro);

    /* create a table */
    parms_TableCreate(&newMap->table, NULL, newMap->lsize);

    MPI_Allgather(&m, 1, MPI_INT, nptr, 1, MPI_INT, newComm);

    /* chunks for contiguous data are distributed across PEs. find
       the range [start_index, end_index) contained in the local PE. */
    start_index = 0;
    for (i = 0; i < pid; i++) {
      start_index += nptr[i];
    }

    end_index   = start_index + m;
    index = 0;
    for (j = start_index; j < end_index; j++) {
      parms_TablePut(newMap->table, j, index);
      newMap->lvars[index++] = j;
    }
    newMap->ispermalloc = true;
    PARMS_NEWARRAY0(newMap->perm,  newMap->lsize);
    PARMS_NEWARRAY0(newMap->iperm, newMap->lsize);
    for (i = 0; i < newMap->lsize; i++) {
      newMap->perm[i] = -1;
    }
  }

  *self = newMap;

/* Define complex data type for MPI if complex code is compiled */
#if defined(DBL_CMPLX)
  parms_InitComplex();
#endif 

  return 0;  
}

/** 
 * Get local index for a given global index.
 *
 * Return a pointer to an integer. If it is NULL, then the variable 
 * with global index gindex doesn't reside on the local
 * PE. Otherwise, it points to an address of a variable 
 * whose value is local index. 
 *
 * @param self    A parms_Map object.
 * @param gindex  A global index.
 * 
 * @return A pointer to an integer whose value is the corresponding
 *         local index. 
 */
int *parms_MapGlobalToLocal(parms_Map self, int gindex)
{
  int *lindex;
  
  if(self->isserial)
     lindex = &gindex;
  else
     lindex = parms_TableGet(self->table, gindex);
  
  return lindex;
}

/**
 * Get global index array for the local PE (see *lvars in parms_map_impl.h)
 *
*/
int parms_MapGetGlobalIndices(parms_Map self, int *im)
{      
	int i, offset = self->start;

	for(i=0; i<self->lsize; i++)
		im[i] = self->lvars[i]+offset;     
	
	return 0;
}
