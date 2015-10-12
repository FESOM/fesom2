/**
 * @file   parms_sys_cmplx.h
 * @author Zhongze Li
 * @date   Tue Oct 17 11:59:09 2006
 * 
 * @brief  Macros and typedef needed by all other header files. 
 * This is used when the complex version of the code is compiled 
 * (with the -DDBL_CMPLX flag)
 *
 * DOK
 * 
 */

#include <complex.h>

typedef struct{
   double real, imag;
}complex_type;

/* Define complex data-type for MPI */
#ifdef USE_MPI
MPI_Datatype MPI_CMPLX; 
MPI_Op MPI_CMPLX_SUM;   
#endif  

/* external BLAS */
#if defined(DBL_CMPLX)
#define FLOAT   complex double
#define ABS_VALUE(a) (cabs(a)) /* definition of absolute value */
#define GDOTC    ZDOTC
#define GDOTU    ZDOTU
#define GCOPY   ZCOPY
#define GSCAL   ZSCAL
#define GAXPY   ZAXPY
#define GNRM2   DZNRM2
#define GGEMV   ZGEMV
#define GGEMM   ZGEMM
#define GGETRF  ZGETRF
#define GGETRS  ZGETRS
#define GGETRI  ZGETRI
#define GGESVD  ZGESVD
#if defined(FORTRAN_CAPS)
#define zdotc_     ZDOTC
#define zdotu_     ZDOTU
#define zcopy_    ZCOPY
#define zscal_    ZSCAL
#define zaxpy_    ZAXPY
#define dznrm2_    DZNRM2
#define zgemv_    ZGEMV
#define zgemm_    ZGEMM
#define zgetrf_   ZGETRF
#define zgetrs_   ZGETRS
#define zgetri_   ZGETRI
#define zgesvd_   ZGESVD
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define zdotc_     zdotc__
#define zdotu_     zdotu__
#define zcopy_    zcopy__
#define zscal_    zscal__
#define zaxpy_    zaxpy__
#define dznrm2_    dznrm2__
#define zgemv_    zgemv__
#define zgemm_    zgemm__
#define zgetrf_   zgetrf__
#define zgetrs_   zgetrs__
#define zgetri_   zgetri__
#define zgesvd_   zgesvd__
#elif !defined(FORTRAN_UNDERSCORE)
#define BLAS_C_INTERFACE
#include <essl.h>
#define zdotc_     zdotc
#define zdotu_     zdotu
#define zcopy_    zcopy
#define zscal_    zscal
#define zaxpy_    zaxpy
#define dznrm2_    dznrm2
#define zgemv_    zgemv
#define zgemm_    zgemm
#define zgetrf_   zgetrf
#define zgetrs_   zgetrs
#define zgetri_   zgetri
#define zgesvd_   zgesvd
#endif

/* FORTRAN subroutines */
#ifndef BLAS_C_INTERFACE
#define ZDOTC(n,x,incx,y,incy)        zdotc_(&(n),(x),&(incx),(y),&(incy))
#define ZDOTU(n,x,incx,y,incy)        zdotu_(&(n),(x),&(incx),(y),&(incy))
#define ZCOPY(n,x,incx,y,incy)       zcopy_(&(n),(x),&(incx),(y),&(incy))
#define ZSCAL(n,alpha,x,incx)        zscal_(&(n),&(alpha),(x), &(incx))
#define ZAXPY(n,alpha,x,incx,y,incy) zaxpy_(&(n), &(alpha), (x), \
					 &(incx), y, &(incy))
#define DZNRM2(n, x, incx)            dznrm2_(&(n), (x), &(incx))
#define ZGEMV(transa, m, n, alpha, a, lda, x, incx, beta, y, incy)  \
  zgemv_((transa), &(m), &(n), &(alpha), (a), &(lda), (x), &(incx), \
	 &(beta), (y), &(incy))
#define ZGEMM(transa,transb,l,n,m,alpha,a,lda,b,ldb,beta,c,ldc)	\
  zgemm_((transa), (transb), &(l), &(n), &(m), &(alpha), (a),	\
	 &(lda), b, &(ldb), &(beta), (c), &(ldc)) 
#define ZGETRF(m, n, a, lda, ipvt, info)		\
  zgetrf_(&(m), &(n), (a), &(lda), (ipvt), &(info))
#define ZGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info) \
  zgetrs_((trans), &(n), &(nrhs), (a), &(lda), (ipiv), (b), &(ldb), &(info))
#define ZGETRI(n, a, lda, ipvt, work, lwork, info)			\
  zgetri_(&(n), (a), &(lda), (ipvt), (work), &(lwork), &(info))

extern complex double zdotc_(int *n, complex double *x, int *incx, complex double *y, int
		    *incy);  
extern complex double zdotu_(int *n, complex double *x, int *incx, complex double *y, int
		    *incy);		    
extern void   zcopy_(int *n, complex double *x, int *incx, complex double *y, int
		    *incy); 
extern void   zscal_(int *n, complex double *alpha, complex double *x, int *incx);
extern void   zaxpy_(int *n, complex double *alpha, complex double *x, int *incx,
		    complex double *y, int *incy);
extern double dznrm2_(int *n, complex double *x, int *incx);
extern void   zgemv_(char *transa, int *m, int *n, complex double *alpha,
		     complex double *a, int *lda, complex double *x, int *incx, complex double
		     *beta, complex double *y, int *incy);
extern void   zgemm_(char *transa, char *transb, int *l, int *m, int
		     *n, complex double *alpha, complex double *a, int *lda, complex double
		     *b, int *ldb, complex double *beta, complex double *c, int *ldc);       
extern void   zgetrf_(int *m, int *n, complex double *a, int *lda, int *ipvt,
		      int *info); 
extern void   zgetrs_(char *trans, int *n, int *nrhs, complex double *a, int
		      *lda, int *ipiv, complex double *b, int *ldb, int
		      *info);   
extern void   zgetri_(int *n, complex double *a, int *lda, int *ipvt, complex double
		      *work, int *lwork, int *info);
#else
#define ZDOTC(n,x,incx,y,incy)        zdotc_((n), (x), (incx), (y), (incy)) 
#define ZDOTU(n,x,incx,y,incy)        zdotu_((n), (x), (incx), (y), (incy)) 
#define ZCOPY(n,x,incx,y,incy)       zcopy_((n), (x), (incx), (y), \
					   (incy)) 
#define ZSCAL(n,alpha,x,incx)        zscal_((n), (alpha), (x), (incx)) 
#define ZAXPY(n,alpha,x,incx,y,incy) zaxpy_((n), (alpha), (x), (incx), \
					   (y), (incy)) 
#define DZNRM2(n,x,incx)              dznrm2_((n), (x), (incx))

#define ZGEMV(transa,m,n,alpha,a,lda,x,incx,beta,y,incy)		\
  zgemv_((transa), (m), (n),						\
	 (alpha), (a), (lda), (x), (incx),				\
	 (beta), (y), (incy))

#define ZGEMM(transa,transb,l,n,m,alpha,a,lda,b,ldb,beta,c,ldc)		\
  zgemm_((transa),(transb),						\
	 (l),(n),(m),(alpha),(a),					\
	 (lda),(b),(ldb),(beta),(c),(ldc))
#define ZGETRF(m, n, a, lda, ipvt, info)  \
  zgetrf_((m), (n), (a), (lda), (ipvt), &(info))
#define ZGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info) \
  zgetrs_((trans), (n), (nrhs), (a), (lda), (ipiv), (b), (ldb), &(info))
#define ZGETRI(n, a, lda, ipvt, work, lwork, info)		\
  zgetri_((n), (a), (lda), (ipvt), (work), (lwork), &(info))
#endif

extern void zgesvd_(char*, char*, int*, int*, complex double*, int*, double*,
		    complex double *, int*, complex double*, int*, complex double*, int*,
		    double*, int*);  
		    
/* givens rotations*/
double sgn(double x, double y);
double sign(double x);
double abs1(complex double x);
double abssq(complex double x);
extern void zclartg(complex double f, complex double g, double *cs, complex double *sn, complex double *rot);	

/*---- initialize complex datatypes and ops --- */
extern void parms_InitComplex();	
void complex_sum(complex_type *in, complex_type *inout, int *len, MPI_Datatype *data_ptr);       

#endif 
