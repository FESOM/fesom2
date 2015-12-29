/**
 * @file   parms_sys_dbl.h
 * @author Zhongze Li
 * @date   Tue Oct 17 11:59:09 2006
 * 
 * @brief  Macros and typedef needed by all other header files. 
 * This is used when the real version of the code is compiled 
 * (with the -DDBL flag)
 *
 * DOK
 * 
 */

/* external BLAS */
#if defined(DBL)
#define FLOAT   double
#define ABS_VALUE(a) (fabs(a)) /* definition of absolute value */
#define GDOT    DDOT
#define GCOPY   DCOPY
#define GSCAL   DSCAL
#define GAXPY   DAXPY
#define GNRM2   DNRM2
#define GDMIN   IDMIN
#define GGEMV   DGEMV
#define GGEMM   DGEMM
#define GGETRF  DGETRF
#define GGETRS  DGETRS
#define GGETRI  DGETRI
#define GGESVD  DGESVD
#if defined(FORTRAN_CAPS)
#define ddot_     DDOT
#define dcopy_    DCOPY
#define dscal_    DSCAL
#define daxpy_    DAXPY
#define dnrm2_    DNRM2
#define idmin_    IDMIN
#define dgemv_    DGEMV
#define dgemm_    DGEMM
#define dgetrf_   DGETRF
#define dgetrs_   DGETRS
#define dgetri_   DGETRI
#define dgesvd_   DGESVD
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define ddot_     ddot__
#define dcopy_    dcopy__
#define dscal_    dscal__
#define daxpy_    daxpy__
#define dnrm2_    dnrm2__
#define idmin_    idmin__
#define dgemv_    dgemv__
#define dgemm_    dgemm__
#define dgetrf_   dgetrf__
#define dgetrs_   dgetrs__
#define dgetri_   dgetri__
#define dgesvd_   dgesvd__
#elif !defined(FORTRAN_UNDERSCORE)
#define BLAS_C_INTERFACE
#define ddot_     ddot
#define dcopy_    dcopy
#define dscal_    dscal
#define daxpy_    daxpy
#define dnrm2_    dnrm2
#define idmin_    idmin
#define dgemv_    dgemv
#define dgemm_    dgemm
#define dgetrf_   dgetrf
#define dgetrs_   dgetrs
#define dgetri_   dgetri
#define dgesvd_   dgesvd
#endif

/* FORTRAN subroutines */
#ifndef BLAS_C_INTERFACE
#define DDOT(n,x,incx,y,incy)        ddot_(&(n),(x),&(incx),(y),&(incy))
#define DCOPY(n,x,incx,y,incy)       dcopy_(&(n),(x),&(incx),(y),&(incy))
#define DSCAL(n,alpha,x,incx)        dscal_(&(n),&(alpha),(x), &(incx))
#define DAXPY(n,alpha,x,incx,y,incy) daxpy_(&(n), &(alpha), (x), \
					 &(incx), y, &(incy))
#define DNRM2(n, x, incx)            dnrm2_(&(n), (x), &(incx))
#define IDMIN(n, sx, incx)           idmin_((&(n), (sx), &(incx))
#define DGEMV(transa, m, n, alpha, a, lda, x, incx, beta, y, incy)  \
  dgemv_((transa), &(m), &(n), &(alpha), (a), &(lda), (x), &(incx), \
	 &(beta), (y), &(incy))
#define DGEMM(transa,transb,l,n,m,alpha,a,lda,b,ldb,beta,c,ldc)	\
  dgemm_((transa), (transb), &(l), &(n), &(m), &(alpha), (a),	\
	 &(lda), b, &(ldb), &(beta), (c), &(ldc)) 
#define DGETRF(m, n, a, lda, ipvt, info)		\
  dgetrf_(&(m), &(n), (a), &(lda), (ipvt), &(info))
#define DGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info) \
  dgetrs_((trans), &(n), &(nrhs), (a), &(lda), (ipiv), (b), &(ldb), &(info))
#define DGETRI(n, a, lda, ipvt, work, lwork, info)			\
  dgetri_(&(n), (a), &(lda), (ipvt), (work), &(lwork), &(info))

extern double ddot_(int *n, double *x, int *incx, double *y, int
		    *incy);  
extern void   dcopy_(int *n, double *x, int *incx, double *y, int
		    *incy); 
extern void   dscal_(int *n, double *alpha, double *x, int *incx);
extern void   daxpy_(int *n, double *alpha, double *x, int *incx,
		    double *y, int *incy);
extern double dnrm2_(int *n, double *x, int *incx);
extern void   idmin_(int *n, double *sx, int *incx);
extern void   dgemv_(char *transa, int *m, int *n, double *alpha,
		     double *a, int *lda, double *x, int *incx, double
		     *beta, double *y, int *incy);
extern void   dgemm_(char *transa, char *transb, int *l, int *m, int
		     *n, double *alpha, double *a, int *lda, double
		     *b, int *ldb, double *beta, double *c, int *ldc);       
extern void   dgetrf_(int *m, int *n, double *a, int *lda, int *ipvt,
		      int *info); 
extern void   dgetrs_(char *trans, int *n, int *nrhs, double *a, int
		      *lda, int *ipiv, double *b, int *ldb, int
		      *info);   
extern void   dgetri_(int *n, double *a, int *lda, int *ipvt, double
		      *work, int *lwork, int *info);
#else
#define DDOT(n,x,incx,y,incy)        ddot_((n), (x), (incx), (y), (incy)) 
#define DCOPY(n,x,incx,y,incy)       dcopy_((n), (x), (incx), (y), \
					   (incy)) 
#define DSCAL(n,alpha,x,incx)        dscal_((n), (alpha), (x), (incx)) 
#define DAXPY(n,alpha,x,incx,y,incy) daxpy_((n), (alpha), (x), (incx), \
					   (y), (incy)) 
#define DNRM2(n,x,incx)              dnrm2_((n), (x), (incx))

#define IDMIN(n,sx,incx)             idmin_((n), (sx), (incx))
#define DGEMV(transa,m,n,alpha,a,lda,x,incx,beta,y,incy)		\
  dgemv_((transa), (m), (n),						\
	 (alpha), (a), (lda), (x), (incx),				\
	 (beta), (y), (incy))

#define DGEMM(transa,transb,l,n,m,alpha,a,lda,b,ldb,beta,c,ldc)		\
  dgemm_((transa),(transb),						\
	 (l),(n),(m),(alpha),(a),					\
	 (lda),(b),(ldb),(beta),(c),(ldc))
#define DGETRF(m, n, a, lda, ipvt, info)  \
  dgetrf_((m), (n), (a), (lda), (ipvt), &(info))
#define DGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info) \
  dgetrs_((trans), (n), (nrhs), (a), (lda), (ipiv), (b), (ldb), &(info))
#define DGETRI(n, a, lda, ipvt, work, lwork, info)		\
  dgetri_((n), (a), (lda), (ipvt), (work), (lwork), &(info))
#endif

extern void dgesvd_(char*, char*, int*, int*, double*, int*, double*,
		    double *, int*, double*, int*, double*, int*,
		    int*);     

#endif 

