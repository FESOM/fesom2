/**
 * @file   parms_sys.h
 * @author Zhongze Li
 * @date   Tue Oct 17 11:59:09 2006
 * 
 * @brief  Macros and typedef needed by all other header files. 
 * 
 */

#ifndef _PARMS_SYS_H_
#define _PARMS_SYS_H_

#ifdef __STDC__
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >=199901L
#define C99
#elif defined(__STDC_VERSION__) && __STDC_VERSION__>=199409L
#define C89
#else
#define C89
#endif
#endif

#ifdef USE_MPI
#define PARMS_ABORT(error_number)			\
  if (true) {						\
    int flag;						\
    MPI_Initialized(&flag);				\
    if (flag) {						\
      MPI_Abort(MPI_COMM_WORLD, (error_number));	\
    }							\
    else {						\
      exit((error_number));				\
    }							\
  } else
#else
#define PARMS_ABORT(error_number)		\
  if (true) {					\
    exit((error_number));			\
  } else
#endif

#include <stdio.h>
#ifdef C99
#include <stdbool.h>
#endif 
#include "mpi.h"

#if defined(__cplusplus)
#define PARMS_CXX_BEGIN extern "C" {
#define PARMS_CXX_END   }
#else
#define PARMS_CXX_BEGIN
#define PARMS_CXX_END
#endif 

PARMS_CXX_BEGIN

typedef enum VARSTYPE   {INTERLACED, NONINTERLACED} VARSTYPE;
typedef enum INSERTMODE {INSERT, ADD} INSERTMODE;
typedef enum NNZSTRUCT {SAME_NONZERO_STRUCTURE, DIFFERENT_NONZERO_STRUCTURE} NNZSTRUCT;
typedef enum COMMTYPE   {P2P, DERIVED} COMMTYPE;
typedef enum SOLVERTYPE {SOLFGMRES, SOLGMRES, SOLBICGS, SOLCG, SOLPBICGS, SOLPBICGS_RAS, SOLBICGS_RAS} SOLVERTYPE;
typedef enum PARAMTYPE  {MAXITS, KSIZE, DTOL, NEIG} PARAMTYPE;
typedef enum MATTYPE    {MAT_NULL=-1, MAT_VCSR=0, MAT_CSR=1} MATTYPE;
typedef enum PCTYPE     {PCBJ, PCSCHUR, PCRAS, PCSCHURRAS} PCTYPE;
typedef enum PCILUTYPE  {PCILU0, PCILUK, PCILUT, PCARMS} PCILUTYPE;

#ifdef C99
typedef _Bool BOOL;
#else
typedef  int  BOOL;
#define true  1
#define false 0
#endif 

#define PARMS_COUT "parms_cout"
#define PARMS_CERR "parms_cerr"

typedef struct parms_FactParam {
  int    schur_start;		/* start position of the Schur complement */
  int    start;			/* start row index of the submatrix in
				   the whole matrix */
  int    n;			/* the size of the local matrix */
  /*
    mc = multicoloring or not (0-not, 1-yes).
    ipar[0:17]  = integer array to store parameters for both
       arms construction (arms2) and iteration (armsol2).
    
       ipar[0]:=nlev.  number of levels (reduction processes).  see
       also "on return" below.

       ipar[1]:= level-reordering option to be used.  
                 if ipar[1]==0 ARMS uses a block independent set ordering
                  with a sort of reverse cutill Mc Kee ordering to build 
                  the blocks. This yields a symmetric ordering. 
                  in this case the reordering routine called is indsetC
                 if ipar[1] == 1, then a nonsymmetric ordering is used.
                  In this case, the B matrix is constructed to be as
                  diagonally dominant as possible and as sparse as possble.
                  in this case the reordering routine called is ddPQ.
    
       ipar[2]:=bsize. Dimension of the blocks. In this version, bsize
       is only a target block size. The size of each block can vary
       and is >= bsize.
    
       ipar[3]:=iout if (iout > 0) statistics on the run are printed
       to FILE *ft

       The following are not used by arms2 -- but should set before
       calling the preconditionning operation armsol2:

       ipar[4]:= Krylov subspace dimension for last level 
         ipar[4] == 0 means only backward/forward solve is performed
            on last level.                   
       ipar[5]:=  maximum # iterations on last level
    
       ipar[6-9] NOT used [reserved for later use] - must be 
       set to zero. [see however a special use of ipar[5] in fgmresC.] 
    
       The following set method options for arms2. Their default
       values can all be set to zero if desired.  
    
       ipar[10-13] == meth[0:3] = method flags for interlevel blocks
       ipar[14-17] == meth[0:3] = method flags for last level block -
         with the following meaning  
           meth[0] permutations of rows  0:no 1: yes. affects rperm
             NOT USED IN THIS VERSION ** enter 0.. Data: rperm
           meth[1] permutations of columns 0:no 1: yes. So far this is 
             USED ONLY FOR LAST BLOCK [ILUTP instead of ILUT]. (so
             ipar[11] does no matter - enter zero). If ipar[15] is one
             then ILUTP will be used instead of ILUT. Permutation data
             stored in: perm2.  
           meth[2] diag. row scaling. 0:no 1:yes. Data: D1
           meth[3] diag. column scaling. 0:no 1:yes. Data: D2
             similarly for meth[14], ..., meth[17] all transformations
             related to parametres in meth[*] (permutation,
             scaling,..) are applied to the matrix before processing
             it   
    
    droptol  = Threshold parameters for dropping elements in ILU
       factorization. 
      droptol[0:4] = related to the multilevel  block factorization
      droptol[5:5] = related to ILU factorization of last block.
        This flexibility is more than is really needed. one can use a
        single parameter for all. it is preferable to use one value
        for droptol[0:4] and another (smaller) for droptol[5:6] 
      droptol[0] = threshold for dropping  in L [B]. See piluNEW.c:
      droptol[1] = threshold for dropping  in U [B].
      droptol[2] = threshold for dropping  in L^{-1} F 
      droptol[3] = threshold for dropping  in E U^{-1} 
      droptol[4] = threshold for dropping  in Schur complement
      droptol[5] = threshold for dropping  in L in last block [see
        ilutpC.c] 
      droptol[6] = threshold for dropping  in U in last block [see
        ilutpC.c] 
    
    lfil     = lfil[0:6] is an array containing the fill-in parameters.
      similar explanations as above, namely:
      lfil[0] = amount of fill-in kept  in L [B]. 
      lfil[1] = amount of fill-in kept  in U [B].
      etc..
    
    tolind   = tolerance parameter used by the indset function. 
      a row is not accepted into the independent set if the *relative*
      diagonal tolerance is below tolind. see indset function for
      details. Good values are between 0.05 and 0.5 -- larger values
      tend to be better for harder problems.
    
    Note:   The comments above are for arms. The first element for every 
    array is for other preconditioner: ilut, iluk
  */
  int    mc;		
  int    lfil[7];	
  int    ipar[18];	
  double droptol[7];
  double pgfpar[2];
  double tolind;
  BOOL   isalloc;
} *parms_FactParam;

#define ZERO      0.0
#define EPSILON   1.0e-20
#define EPSMAC    1.0e-16  

/* Compile real or complex version of code */
#if defined(DBL_CMPLX)
#include "parms_sys_cmplx.h"
#elif defined(DBL)
#include "parms_sys_dbl.h"
#else
#include "parms_sys_dbl.h"
#endif 

PARMS_CXX_END

#endif 
