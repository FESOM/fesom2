#ifndef __ITSOL_INCLUDED_PROTOS_H__
#define __ITSOL_INCLUDED_PROTOS_H__

#include "globheads.h"

#if defined(FORTRAN_CAPS)
#define qsplit  QSPLIT
#define readmtc READMTC
#define csrcsc  CSRCSC
#define roscal  ROSCAL
#define coscal  COSCAL
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define qsplit  qsplit__
#define readmtc readmtc__
#define csrcsc  csrcsc__
#define roscal  roscal__
#define coscal  coscal__
#elif defined(FORTRAN_UNDERSCORE)
#define qsplit  qsplit_
#define readmtc readmtc_
#define csrcsc  csrcsc_
#define roscal  roscal_
#define coscal  coscal_
#endif

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

/* sets */
extern void errexit(char *f_str, ...);
extern void *Malloc(int nbytes, char *msg); 
extern int setupP4 (p4ptr amat, int Bn, int Cn,  csptr F,  csptr E);
extern int cleanP4(p4ptr amat);
extern int setupILUT(ilutptr amat, int len);
extern int cleanILUT(ilutptr amat, int indic);
extern int setupCS(csptr amat, int len, int job); 
extern int cleanCS(csptr amat);
extern int nnzCS(csptr amat); 
extern int nnz_arms (arms PreSt);
extern int cs_nnz (csptr A) ;
extern int cscpy(csptr amat, csptr bmat);
extern int cleanILU( iluptr lu );
extern int mallocRow( iluptr lu, int nrow );
extern int CSRcs(int n, FLOAT *a, int *ja, int *ia, csptr mat, int rsa);
extern int lev4_nnz(p4ptr levmat, int *lev); 
extern int setupILU( iluptr lu, int n );
extern void setup_arms (arms Levmat);
extern int cleanARMS(arms ArmsPre);
extern int csSplit4(csptr amat, int bsize, int csize, csptr B, csptr F,
		    csptr E, csptr C); 
/* MatOps */
extern void matvec(csptr mata, FLOAT *x, FLOAT *y); 
extern void matvecz(csptr mata, FLOAT *x, FLOAT *y, FLOAT *z);
extern void luinv(int n, FLOAT *a, FLOAT *x, FLOAT *y); 
extern int lusolC( FLOAT *y, FLOAT *x, iluptr lu ); 
extern int rpermC(csptr mat, int *perm); 
extern int cpermC(csptr mat, int *perm) ; 
extern int dpermC(csptr mat, int *perm) ; 
extern int CSparTran(csptr amat, csptr bmat, CompressType *compress);
extern void invsp(int start, ilutptr ilusch, FLOAT *y, FLOAT *x);
extern int armsol2(FLOAT *x,  arms Prec);
extern int ascend (p4ptr levmat, FLOAT *x, FLOAT *wk);
extern int descend(p4ptr levmat, FLOAT *x, FLOAT *wk);
extern p4ptr Lvsol2(FLOAT *x, int nlev, p4ptr levmat, ilutptr ilusch,
		    int flag); 
extern int   Uvsol2(FLOAT *x, int nlev, int n, p4ptr levmat, ilutptr
		    ilusch); 
extern void  SchLsol(ilutptr ilusch, FLOAT *y) ;
extern void  SchUsol(ilutptr ilusch, FLOAT *y) ;
extern void Lsolp(int start, csptr mata, FLOAT *b, FLOAT *y);
extern void Usolp(int start, csptr mata, FLOAT *y, FLOAT *x);
extern int invGauss(int nn, FLOAT *A); 
extern int invSVD(int nn, FLOAT *A) ;
extern void arms_Usol(csptr mata, FLOAT *b, FLOAT *x);
extern void arms_Lsol(csptr mata, FLOAT *b, FLOAT *x);
extern int condestArms(arms armspre, FLOAT *y, FILE *fp );

/* misc.c */
extern int SparTran(csptr amat, csptr bmat, int job, int flag); 
extern int coscalC(csptr mata, double *diag, int nrm);
extern void dscale(int n, double *dd, FLOAT *x, FLOAT * y);
extern void hilosort(csptr mat, int abval, int hilo);
extern void printmat(FILE *ft, csptr A, int i0, int i1);
extern void qqsort(int *ja, FLOAT *ma, int left, int right);
extern void qsort2C(int *ja, FLOAT *ma, int left, int right, int
		    abval); 
extern void qsort3i(int *wa, int *cor1, int *cor2, int left, int
		    right); 
extern void qsortC(int *ja, FLOAT *ma, int left, int right, int
		   abval); 
extern void qsortR2I(double *wa, int *cor1, int *cor2, int left, int
		     right); 
extern int qsplitC(FLOAT *a, int *ind, int n, int ncut);
extern int roscalC(csptr mata, double *diag, int nrm);
extern void swapj(int v[], int i, int j);
extern void swapm(FLOAT v[], int i, int j);
extern void dswapm(double v[], int i, int j);

/* piluNEW.c */
extern int pilu(p4ptr amat, csptr B, csptr C, double *droptol, int
		*lfil, csptr schur);

/* ilutpC.c */
extern int ilutD(csptr amat, double *droptol, int *lfil, ilutptr
		 ilusch);
extern int ilutpC(csptr amat, double *droptol, int *lfil, double
		  permtol, int mband, ilutptr ilusch);

/* PQ.c */
extern int PQperm(csptr mat, int *Pord, int *Qord, int
		  *nnod, double tol, int nbnd);
extern int preSel(csptr mat, int *icor, int *jcor, int job, double
		  tol, int *count, int nbnd);
extern int weightsC(csptr mat, double *w);
extern int add2com(int *nback, int nod, int *iord, int *riord);
extern int add2is(int *last, int nod, int *iord, int *riord);
extern int indsetC(csptr mat, int bsize, int *iord, int *nnod, double
		   tol,int nbnd); 

/* setblks.c */
extern int KeyComp( const void *vfst, const void *vsnd );
extern int init_blocks( csptr csmat, int *pnBlock, int **pnB, int
			**pperm, double eps, double *t_hash, double
			*t_angle );

/* systimer.c */
extern double sys_timer();

#endif 
