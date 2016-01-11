/* 
 * fort_part.c
 *
 */

/* 

 a Fortran interface to METIS-5.1 derived from a metis-4.0-interface
  as part of FoSSI, the Family of Simplified Solver Interfaces by
  Stephan Frickenhaus, Computing-Center, AWI-Bremerhaven, Germany
  www.awi-bremerhaven.de
 
*/

/* SF 06/2001: The partitioner Interface to METIS-4.0 */
/* NR 12/2015: The partitioner Interface to METIS-5.1.0 */
/* n: dimension of system (number of equations=number of variables) */
/*    compressed matrix on input: */
/* ptr : rowptr or colptr */
/* adj : colind or rowind */
/* wgt : one weight, e.g., number of levels on each node */
/* np  : partitioning for np processors */
/* part : partitioning vector */
#include <stdio.h>
#include <stdlib.h>
#include "metis.h"

void partit(idx_t *n, idx_t *ptr, idx_t *adj, idx_t *wgt, idx_t *np, idx_t *part)
{


  int i, wgt_type ;
  idx_t opt[METIS_NOPTIONS];  
  idx_t ncon, ec;
  idx_t *wgt_2d3d;
  int ierr;

  if (*np==1) { for(i=0;i<*n;i++) part[i]=0; return;}

#ifdef PART_WEIGHTED
  wgt_type = 2;  /* 3D and 2D equally distributed */
    /* wgt_type = 1; /* only 3D equally distributed, like old partit2  */
#else
  wgt_type = 0; /* only 2D nodes equally distributed */
#endif

  METIS_SetDefaultOptions(opt);

  /* The following values are the result of some tests with mesh_aguv, wgt_type=2. */
  /* METIS_PartGraphRecursive resulted in a far better partition than Kway. According to */
  /* forum entries, one should always compare. There is no rule which one works best. */ 

  opt[METIS_OPTION_CONTIG]    = 0;  /* 1: contiguous partitions, please */
                                   /* With more weights, this makes the partitions uglier... */
                                   /* Ignored by METIS_PartGraphRecursive */

  opt[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_CUT; /* Minimize edge cut. */

  /* Alternative: METIS_OBJTYPE_VOL, minimize number of nodes on the interface. */ 
  /* It is nearly the same, as the nodes all have weight 1 in this regard */
  /* (one of the NULL fields in our case). But it does not matter if a node is connected */ 
  /* to the other partition by e.g., 1 or 5 edges, thus, _VOL is what we want. */ 
   
    /* But: _VOL only works with METIS_PartGraphKway*/

  opt[METIS_OPTION_NUMBERING] = 1; /* Fortran Numbering */
  opt[METIS_OPTION_NCUTS]     = 2; /* Build NCUTS partitions, choose the best */
  opt[METIS_OPTION_NITER]     = 25; /* higher => better quality, more compute time. Default: 10 */

  opt[METIS_OPTION_UFACTOR]   = 1; /* max. allowed load inbalance in promille */
  /* Well, somehow relative. Real life inbalance ends up at about 20%, comparing  */
  /* max(#3D-nodes)/min(#3D-nodes), with METIS_PartGraphKway. */
  /* But: lower value -> lower maximum number of nodes in all partitions. This it what counts in the end! */
  /* Far better load balancing with METIS_PartGraphRecursive. */

/*   opt[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; */
/*   opt[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM; */
/*   opt[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP2SIDED; */
  
  /* Available Options: */
  /* METIS PartGraphRecursive   METIS PartGraphKway */
  /* METIS_OPTION_OBJTYPE,  -           x */
  /* METIS_OPTION_CTYPE,    x           x */
  /* METIS_OPTION_IPTYPE,   x           x */
  /* METIS_OPTION_RTYPE,    x           x */
  /* METIS_OPTION_NO2HOP,   x           x */
  /* METIS_OPTION_NCUTS,    x           x */
  /* METIS_OPTION_NITER,    x           x */
  /* METIS_OPTION_SEED,     x           x */
  /* METIS_OPTION_UFACTOR,  x           x */
  /* METIS_OPTION_MINCONN,  -           x */
  /* METIS_OPTION_CONTIG,   -           x */
  /* METIS_OPTION_NUMBERING x           x */
  /* METIS_OPTION_DBGLVL    x           x */

  if (wgt_type == 0) {
    ncon = 1;   
    printf("Distribution weight: 2D nodes\n"); 
    /* ierr = METIS_PartGraphKway(n,&ncon,ptr,adj,NULL,NULL,NULL,np,NULL,NULL,opt,&ec,part); */
    ierr = METIS_PartGraphRecursive(n,&ncon,ptr,adj,NULL,NULL,NULL,np,NULL,NULL,opt,&ec,part);
    
  } else if  (wgt_type == 1) {
    ncon = 1;
    printf("Distribution weight: 3D nodes\n");
    /* ierr = METIS_PartGraphKway(n,&ncon,ptr,adj,wgt, NULL,NULL,np,NULL,NULL,opt,&ec,part); */
    ierr = METIS_PartGraphRecursive(n,&ncon,ptr,adj,wgt, NULL,NULL,np,NULL,NULL,opt,&ec,part);
    
  } else  {
    /*    quasi Default */
    ncon=2;
    wgt_2d3d = (idx_t*) malloc(2* *n *sizeof(idx_t));
    for (i=0; i<*n; i++){
      wgt_2d3d[2*i]   = 1;
      wgt_2d3d[2*i+1] = wgt[i]+100; /* soften the 3D-criteria to allow for better general quality */
    }
    
    printf("Distribution weight: 2D and 3D nodes\n");
    /* ierr = METIS_PartGraphKway(n,&ncon,ptr,adj,wgt_2d3d,NULL,NULL,np,NULL,NULL,opt,&ec,part);  */
    ierr = METIS_PartGraphRecursive(n,&ncon,ptr,adj,wgt_2d3d,NULL,NULL,np,NULL,NULL,opt,&ec,part);

    free(wgt_2d3d);
  }
   
  if (ierr != METIS_OK) {
    printf("MEITS finished with error, code=%d\n", ierr);
  } else { 
    printf("METIS edgecut %d\n", ec);
    
    for(i=0;i<*n;i++) part[i]--;
  }
}
