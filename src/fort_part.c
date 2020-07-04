/* 
 * fort_part.c
 *
 */

#ifndef METIS_VERSION
#define METIS_VERSION 5
#endif

/* #define GRAPH_OUTPUT          /\* Output graph into a file in addition to partitioning *\/ */
#define USE_EDGE_WEIGHTS      /* Adds edge weights to METIS partitioning */
#define MAX_HIER_LEVELS 10    /* Maximum number of hierarchical partitioning levels 
                                 (should be set equal to one for compatibility with old versions) */

/* 

 a Fortran interface to METIS-5.1 derived from a metis-4.0-interface
  as part of FoSSI, the Family of Simplified Solver Interfaces by
  Stephan Frickenhaus, Computing-Center, AWI-Bremerhaven, Germany
  www.awi-bremerhaven.de

  Hierarchic partitioning added by 
  Vadym Aizinger, Computing-Center, AWI-Bremerhaven, Germany
  www.awi-bremerhaven.de
 
*/

/* SF 06/2001: The partitioner Interface to METIS-4.0 */
/* NR 12/2015: The partitioner Interface to METIS-5.1.0 */
/* NR 11/2016: Combined Interface to Metis 4.0 and 5.1.0, */
/*             switch at compile time with e.g. -DMETISVERISON=4 */
/* VA 05/2017: Hierarchic partitioning added for METIS-5.1.0 interface */
/* n: dimension of system (number of equations=number of variables) */
/*    compressed matrix on input: */
/* ptr : rowptr or colptr */
/* adj : colind or rowind */
/* wgt : one weight, e.g., number of levels on each node */
/* np  : partitioning for np processors */
/* part : partitioning vector */

#include <stdio.h>
#include <stdlib.h>

#if METIS_VERSION == 5 /* ---------------- METIS 5 part ------------------------ */
#include "metis.h"

void partit(idx_t *n, idx_t *ptr, idx_t *adj, idx_t *wgt, idx_t *np, idx_t *part)
{
  int i, j, wgt_type;
  idx_t opt[METIS_NOPTIONS];
  idx_t ncon, ec;
  idx_t *wgt_2d3d = NULL;
  idx_t *wgt_edge = NULL;
  int ierr;
  int n_levels, n_part_low = 1;     /* Number of lower hierarchy levels and partitions in them */

  static int current_level = -1; /* Current level in hierarchical partitioning */

  current_level++;
  if(np[current_level]<1)
  {
    printf("Partitioning on level %d called with less than one partition!\n", current_level+1);
    exit(1);
  }

  if (np[current_level]==1) { for(i=0;i<*n;i++) part[i]=0; return;}

  for (n_levels=current_level+1; n_levels<MAX_HIER_LEVELS; n_levels++)
  {
    if (np[n_levels]==0)
      break;
    n_part_low*=np[n_levels];
  }

  if (current_level==0)
  {
    printf("Starting hierarchic partitioning with %d level(s)\n", n_levels);
    printf("Level:");
    for (i=0; i<n_levels; i++)
      printf("\t%d", i+1);
    printf("   Total\nNpart:");
    for (i=0; i<n_levels; i++)
      printf("\t%d", np[i]);
    printf("   %d\n", n_part_low*np[current_level]);
  }

  printf("Partitioning into %d subdomains on level %d with Metis version %d.%d.%d\n",
         np[current_level], current_level+1, METIS_VER_MAJOR,METIS_VER_MINOR,METIS_VER_SUBMINOR);

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
  opt[METIS_OPTION_NCUTS]     = 10; /* Build NCUTS partitions, choose the best */
  opt[METIS_OPTION_NITER]     = 15; /* higher => better quality, more compute time. Default: 10 */

  opt[METIS_OPTION_UFACTOR]   = 1; /* max. allowed load inbalance in promille */
  /* Well, somehow relative. Real life inbalance ends up at about 20%, comparing  */
  /* max(#3D-nodes)/min(#3D-nodes), with METIS_PartGraphKway. */
  /* But: lower value -> lower maximum number of nodes in all partitions. This it what counts in the end! */
  /* Far better load balancing with METIS_PartGraphRecursive. */

#ifdef METISRANDOMSEED
  opt[METIS_OPTION_SEED] = METISRANDOMSEED; /* Set in the Makefile */
  /* Useful, if the first partition did not work in FESOM */
  /* Try different seeds until the partition is ok. */ 
#endif

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

  switch(wgt_type) 
  {
  case 0:
    ncon = 1;   
    printf("Distribution weight: 2D nodes\n"); 
    wgt_2d3d = NULL;
    /* ierr = METIS_PartGraphKway(n,&ncon,ptr,adj,NULL,NULL,NULL,np,NULL,NULL,opt,&ec,part); */
/*     ierr = METIS_PartGraphRecursive(n,&ncon,ptr,adj,NULL,NULL,wgt_edge,np+current_level,NULL,NULL,opt,&ec,part); */
    break;
  case 1:
    ncon = 1;
    wgt_2d3d = wgt;
    printf("Distribution weight: 3D nodes\n");
    /* ierr = METIS_PartGraphKway(n,&ncon,ptr,adj,wgt, NULL,NULL,np,NULL,NULL,opt,&ec,part); */
/*     ierr = METIS_PartGraphRecursive(n,&ncon,ptr,adj,wgt, NULL,wgt_edge,np+current_level,NULL,NULL,opt,&ec,part); */
    break;
  case 2:
    ncon=2;
    wgt_2d3d = (idx_t*) malloc(*n*2*sizeof(idx_t));
    if (wgt_2d3d == NULL)
    {
      printf("Error allocating dynamic memory for hierarchical partitioning!\n");
      exit(1);
    }
    for (i=0; i<*n; i++){
      wgt_2d3d[2*i]   = 1;
      wgt_2d3d[2*i+1] = wgt!=NULL ? wgt[i]+100 : 100; /* soften the 3D-criteria to allow for better general quality */
    }
    
    printf("Distribution weight: 2D and 3D nodes\n");
    /* ierr = METIS_PartGraphKway(n,&ncon,ptr,adj,wgt_2d3d,NULL,NULL,np,NULL,NULL,opt,&ec,part);  */
/*     ierr = METIS_PartGraphRecursive(n,&ncon,ptr,adj,wgt_2d3d,NULL,wgt_edge,np+current_level,NULL,NULL,opt,&ec,part); */
/*     free(wgt_2d3d); */
    break;
  default:
    printf("Wrong weighting option for partitioning\n");
    exit(1);
  }
   
#ifdef USE_EDGE_WEIGHTS /* At the coarsest level, set edge weights to the number of the 3D nodes */
  if (current_level==0 && wgt_type!=0)
  {
    wgt_edge = calloc(ptr[*n], sizeof(*wgt_edge));
    if (wgt_edge == NULL)
    {
      printf("Error allocating dynamic memory for hierarchical partitioning!\n");
      exit(1);
    }

    for (i=0; i<*n; i++)
      for (j=ptr[i]-1; j<ptr[i+1]-1; j++)
        wgt_edge[j] = wgt!=NULL ? wgt[i]+wgt[adj[j]-1] : 1;
  }
#endif

#ifdef GRAPH_OUTPUT /* Output the graph to file and exit */
  FILE *coord_file=fopen("fesom.graph", "w");
  if (!coord_file)
  {
    printf("Cannot open graph file!\n");
    exit(1);
  }

  /* Output header */
  fprintf(coord_file, "%d %d %s %d\n", *n, (ptr[*n])/2, "11", ncon);

  /* Output the data for each vertex of the graph */
  for (i=0; i<*n; i++)
  {
    for (j=0; j<ncon; j++)
      fprintf(coord_file, " %d", wgt_2d3d!=NULL? wgt_2d3d[i*ncon+j] : 0);
    for (j=ptr[i]-1; j<ptr[i+1]-1; j++)
    {
      fprintf(coord_file, " %d %d", adj[j], wgt_edge!=NULL ? wgt_edge[j] : 1);
    }
    fprintf(coord_file, "\n");
  }
    
  fclose(coord_file);
#endif

  ierr = METIS_PartGraphRecursive(n,&ncon,ptr,adj,wgt_2d3d,NULL,wgt_edge,np+current_level,NULL,NULL,opt,&ec,part);
  if (ierr != METIS_OK) {
    printf("MEITS finished with error, code=%d\n", ierr);
    exit(1);
  } else { 
    printf("METIS edgecut %d\n", ec);
    
    for(i=0;i<*n;i++) part[i]--;
  }

  if (current_level < n_levels-1)
  {
    idx_t *vertex_glob_loc = calloc(*n, sizeof(*vertex_glob_loc)); /* Local (C-)index of a graph vertex */
    idx_t *ptr_loc = calloc(*n+1, sizeof(*ptr_loc));
    idx_t *adj_loc = calloc(ptr[*n], sizeof(*adj_loc));
    idx_t *wgt_loc = calloc(*n, sizeof(*wgt_loc));
    idx_t *part_loc = calloc(*n, sizeof(*part_loc));
    if (vertex_glob_loc == NULL || ptr_loc == NULL || adj_loc == NULL || 
        wgt_loc == NULL || part_loc == NULL)
    {
      printf("Error allocating dynamic memory for hierarchical partitioning!\n");
      exit(1);
    }

    for(i=np[current_level]-1; i>=0; i--)
    {
      idx_t n_loc = 0;

      /* Search for graph vertices belonging to the i-th partition */
      for (j=0; j<*n; j++)
        if (part[j]==i)
          vertex_glob_loc[j] = n_loc++;

      /* Convert graph to lower level indexing */
      ptr_loc[0]=1;
      for (j=0; j<*n; j++)
        if (part[j]==i)
        {
          idx_t k, vertex_loc = vertex_glob_loc[j], index_loc = ptr_loc[vertex_loc]-1;

          for (k=ptr[j]-1; k<ptr[j+1]-1; k++)
            if (part[adj[k]-1]==i)
              adj_loc[index_loc++]=vertex_glob_loc[adj[k]-1]+1;
          ptr_loc[vertex_loc+1]=index_loc+1;
          if (wgt!=NULL)
            wgt_loc[vertex_loc]=wgt[j];
        }

      partit(&n_loc, ptr_loc, adj_loc, wgt_loc, np, part_loc);

      /* Convert the partitioned graph back to the current level indexing */
      for (j=0; j<*n; j++)
        if (part[j]==i)
          part[j] = i*n_part_low+part_loc[vertex_glob_loc[j]];
    }
    free(vertex_glob_loc);
    free(ptr_loc);
    free(adj_loc);
    free(wgt_loc);
    free(part_loc);
  }

  if (wgt_2d3d != wgt_edge && wgt_2d3d != NULL)
    free(wgt_2d3d);
  if (wgt_edge != NULL)
    free(wgt_edge);
  /* Decrement current_level variable before returing to the previous level */
  current_level--;
}

#elif METIS_VERSION == 4 /* ---------------- METIS 4 part ------------------------ */

void partit(int *n, int *ptr, int *adj, int *wgt, int *np, int *part)
{
  int opt[5];
  int numfl=1; /* 0: C-numbering ; 1: F-numbering*/
  int vflg=2; /* weights on vertices */
  int i, wgt_type ;
  int ncon, ec;
  int *wgt_2d3d;


  if (*np==1) { for(i=0;i<*n;i++) part[i]=0; return;}

  printf("Start partitioning with Metis version 4 \n");

#ifdef PART_WEIGHTED
  wgt_type = 2;  /* 3D and 2D equally distributed */
  /* wgt_type = 1; /* only 3D equally distributed, like old partit2  */
#else
  wgt_type = 0; /* only 2D nodes equally distributed */
#endif

opt[0]=0; /* default options */



  if (wgt_type == 0) {
    printf("Distribution weight: 2D nodes\n"); 
    vflg=0; /* no weights */
        METIS_PartGraphKway(n,ptr,adj,NULL,NULL,&vflg,&numfl,np,opt,&ec,part); 
/* METIS_PartGraphRecursive(n,ptr,adj,NULL,NULL,&vflg,&numfl,np,opt,&ec,part); */
    
  } else if  (wgt_type == 1) {
    vflg=2; /* weights on vertices*/
    printf("Distribution weight: 3D nodes\n");
        METIS_PartGraphKway(n,&ncon,ptr,adj,wgt,NULL,&vflg,&numfl,np,opt,&ec,part); 
/* METIS_PartGraphRecursive(n,&ncon,ptr,adj,wgt,NULL,&vflg,&numfl,np,opt,&ec,part); */
    
  } else  {
    /*    quasi Default */
    ncon=2;
    vflg=2; /* weights on vertices*/
    wgt_2d3d = (int*) malloc(2* *n *sizeof(int));
    for (i=0; i<*n; i++){
      wgt_2d3d[2*i]   = 1;
      wgt_2d3d[2*i+1] = wgt[i]+100; /* soften the 3D-criteria to allow for better general quality */
    }
    
    printf("Distribution weight: 2D and 3D nodes\n");
        METIS_mCPartGraphKway(n,&ncon,ptr,adj,wgt_2d3d,NULL,&vflg,&numfl,np,opt,&ec,part);  
/* METIS_mCPartGraphRecursive(n,&ncon,ptr,adj,wgt_2d3d,NULL,&vflg,&numfl,np,opt,&ec,part); */

    free(wgt_2d3d);
  }
   
  if (ec < 0) {
    printf("METIS finished with error, edgecut=%d\n", ec);
  } else { 
    printf("METIS edgecut %d\n", ec);
    for(i=0;i<*n;i++) part[i]--;    
  }

}
#endif
