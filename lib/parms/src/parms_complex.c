#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parms_sys.h"

#if defined(DBL_CMPLX)
#define epsmac 1.0e-16

/*---------------------------------------------------------------------------*
 * A C implementation of the LAPACK auxilliary routine clartg.f.
 * Computes the plane givens rotation so that:
 *
 *     [  CS         SN  ] . [ F ]  =  [ rot ]
 *     [ -conjg(SN)  CS  ]   [ G ]     [ 0 ]
 *
 *  where CS**2 + |SN|**2 = 1 and Re(rot) >= 0.
 *
 * The algorithm implements both parts of the LAPACK3E clartg routine
 * and algorithm 3 in "On computing Givens Rotations Reliably and Efficiently"
 * by Bindel et al.
 *---------------------------------------------------------------------------*/

/* Auxilliary functions */
/* -------------------------------------------*/
double sgn(double x, double y){
  /* Sign transfer function */
  if (y >= 0.0) return fabs(x);
  else return -fabs(x);
} 
double sign(double x){
  return x/fabs(x);
} 
double abs1(complex double x){
  return fabs(creal(x)) + fabs(cimag(x));
}
double abssq(complex double x){
  return (pow(creal(x), 2) + pow(cimag(x),2));
  //  return cpow(cabs(x),2);
}
/* ----------------------------------------------*/

/* Givens Rotation */
void zclartg(complex double f, complex double g, double *cs, complex double *sn, complex double *rot)
{
	double one=1.0, zero=0.0;
   complex double czero = 0.0 + 0.0*I; 
   double D, F2, G2, FG2;
   if (cabs(g) <= epsmac)   {
     *cs = sgn(one, creal(f));
     *sn = czero;
     *rot = f*(*cs);
   }
    else if (cabs(f) <= epsmac){
      *cs = zero;
      *sn = conj(g)/cabs(g); 
      *rot = cabs( g );
    } 
   else{
     F2 = abssq(f); 
     G2 = abssq( g );
     FG2 = F2 + G2;
     if(fabs(FG2) <= epsmac) FG2 = epsmac;
     D = 1/sqrt(F2*FG2); 
     *cs = F2 * D;
     FG2 = FG2*D;  
     *rot = f*FG2;
     *sn = f*D;
     *sn = conj(g)*(*sn); 
   }
}

/*-------------- complex addition operation for MPI -----*/
/* This function defines the addition operation for an array 
 * of complex numbers, for use in MPI reduction operations.
 * It is defined here, but not use, since everything seems to 
 * work well with MPI_SUM when the complex number is defined 
 * as a contiguous array of two doubles (see definition for MPI_CMPLX).
 */
void complex_sum(complex_type *in, complex_type *inout, int *len, MPI_Datatype *data_ptr)
{ 
  int i;						
  complex_type sum; 					
  for(i = 0; i<(*len); i++)				
  {							
  	sum.real = in->real + inout->real;		
  	sum.imag = in->imag + inout->imag;		
  	*inout = sum;					
  	in++;						
  	inout++;					
  }	
}

/* ----------------- Initialize complex data type and ops for MPI ----*/
void parms_InitComplex()
{
  MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_CMPLX);
  MPI_Type_commit( &MPI_CMPLX );
  
  MPI_Op_create((MPI_User_function *)complex_sum, true, &MPI_CMPLX_SUM);
}
   

#endif
