/* 
 * fort_part.c
 * 3.4 sf
 *
 */

/* 

 a Fortran interface to METIS-4.0 and simple block/range partitioners
  as part of FoSSI, the Family of Simplified Solver Interfaces by
  Stephan Frickenhaus, Computing-Center, AWI-Bremerhaven, Germany
  www.awi-bremerhaven.de
 
*/

/* SF 06/2001 */
/* The partinioner Interface to METIS-4.0 */
/* n: dimension of system (number of equations=number of variables) */
/*    compressed matrix on input: */
/* ptr : rowptr or colptr */
/* adj : colind or rowind */
/* np  : partitioning for np processors */
/* part : partitioning vector */

#include <stdio.h>


#ifdef T3E
#define IndxInt short
#endif
#if SGI | SUN | IBM | LINUX
#define IndxInt int
#endif

#ifdef UNDER_
void partit_block_
#endif
#ifdef UPCASE
void PARTIT_BLOCK
#endif
#ifdef MATCH
void partit_block
#endif
(int *n, IndxInt *ptr, IndxInt *adj, int *np, IndxInt *part)
{
int i ;

if (*np==1) { for(i=0;i<*n;i++) part[i]=0; return;}
for (i=0;i<*n;i++) {
 part[i]=*np * i / *n; 
 if (part[i]==*np) part[i]=*np-1;
 }
}

#ifdef UNDER_
void partit_range_
#endif
#ifdef UPCASE
void PARTIT_RANGE
#endif
#ifdef MATCH
void partit_range
#endif
(int *n, int *np, IndxInt *part)
{
int i,s ;
 if (*np==1) {part[0]=1;part[1]=*n+1;return;}
 s=0;
 part[0]=1;
 for (i=1;i<=*np;i++) 
  { part[i]= part[i-1]+*n/(*np);
    if (i==*np) part[i]=*n+1;
    /*    printf("PART_RANGE %d %d\n",i,part[i]);  */
  }
}

#ifdef UNDER_
void partit_
#endif
#ifdef UPCASE
void PARTIT
#endif
#ifdef MATCH
void partit
#endif
(int *n, IndxInt *ptr, IndxInt *adj, int *np, IndxInt *part)
{
int i ;
int *null=0;
int opt[5];
int numfl=1;
int vflg=0;
int ec;
if (*np==1) { for(i=0;i<*n;i++) part[i]=0; return;}
opt[0]=0;
#ifdef DBGPART
printf("Partit %d %d %d %d %d %d\n",*n,*np,ptr[0],ptr[1],adj[0],adj[1]);
#endif
METIS_PartGraphKway(n,ptr,adj,null,null,&vflg,&numfl,np,opt,&ec,part);
#ifdef DBGPART
printf("Edgecut %d\n", ec);
#endif
if (numfl) for(i=0;i<*n;i++) part[i]--;
}

#ifdef UNDER_
void partit2_
#endif
#ifdef UPCASE
void PARTIT2
#endif
#ifdef MATCH
void partit2
#endif
(int *n, IndxInt *ptr, IndxInt *adj, IndxInt *wgt, int *np, IndxInt *part)


{

int i ;
int *null;
int opt[5];
int numfl=1; /* 0: C-numbering ; 1: F-numbering*/
int vflg=2; /* weights on vertices */
int ec;

if (*np==1) { for(i=0;i<*n;i++) part[i]=0; return;}

opt[0]=0; /* default options */

printf("Partit %d %d %d %d %d %d\n",*n,*np,ptr[0],ptr[1],adj[0],adj[1]);

null=0; /* eso es mui importante */

METIS_PartGraphVKway(n,ptr,adj,wgt,null,&vflg,&numfl,np,opt,&ec,part);

printf("Edgecut %d\n", ec);
for(i=0;i<*n;i++) part[i]--;
}
