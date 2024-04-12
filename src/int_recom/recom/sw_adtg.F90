!> \file sw_adtg.f90
!! \BRIEF 
!> Module with sw_adtg function - compute adiabatic temp. gradient from S,T,P
MODULE msw_adtg
CONTAINS
!>  Function to calculate adiabatic temperature gradient as per UNESCO 1983 routines.
FUNCTION sw_adtg  (s,t,p)

  !     ==================================================================
  !     Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================
  USE msingledouble
  IMPLICIT NONE
  !> salinity [psu (PSU-78)]
  REAL(kind=r8) :: s
  !> temperature [degree C (IPTS-68)]
  REAL(kind=r8) :: t
  !> pressure [db]
  REAL(kind=r8) :: p

  REAL(kind=r8) :: a0,a1,a2,a3,b0,b1,c0,c1,c2,c3,d0,d1,e0,e1,e2
  REAL(kind=r8) :: sref

  REAL(kind=r8) :: sw_adtg

  sref = 35.d0
  a0 =  3.5803d-5
  a1 = +8.5258d-6
  a2 = -6.836d-8
  a3 =  6.6228d-10

  b0 = +1.8932d-6
  b1 = -4.2393d-8

  c0 = +1.8741d-8
  c1 = -6.7795d-10
  c2 = +8.733d-12
  c3 = -5.4481d-14

  d0 = -1.1351d-10
  d1 =  2.7759d-12

  e0 = -4.6206d-13
  e1 = +1.8676d-14
  e2 = -2.1687d-16

  sw_adtg =  a0 + (a1 + (a2 + a3*T)*T)*T &
       + (b0 + b1*T)*(S-sref) &
       + ( (c0 + (c1 + (c2 + c3*T)*T)*T) + (d0 + d1*T)*(S-sref) )*P &
       + (  e0 + (e1 + e2*T)*T )*P*P

END FUNCTION sw_adtg


!>  Function to calculate adiabatic temperature gradient as per UNESCO 1983 routines.
!! and derivative with respect to temperature and salinity
FUNCTION sw_adtg_DNAD  (s,t,p)

  !     ==================================================================
  !     Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================
  USE msingledouble
  USE Dual_Num_Auto_Diff
  IMPLICIT NONE
  !> salinity [psu (PSU-78)]
  TYPE(DUAL_NUM) :: s
  !> temperature [degree C (IPTS-68)]
  TYPE(DUAL_NUM) :: t
  !> pressure [db]
  TYPE(DUAL_NUM) :: p

  REAL(kind=r8) :: a0,a1,a2,a3,b0,b1,c0,c1,c2,c3,d0,d1,e0,e1,e2
  REAL(kind=r8) :: sref

  TYPE(DUAL_NUM) :: sw_adtg_DNAD

  sref = 35.d0
  a0 =  3.5803d-5
  a1 = +8.5258d-6
  a2 = -6.836d-8
  a3 =  6.6228d-10

  b0 = +1.8932d-6
  b1 = -4.2393d-8

  c0 = +1.8741d-8
  c1 = -6.7795d-10
  c2 = +8.733d-12
  c3 = -5.4481d-14

  d0 = -1.1351d-10
  d1 =  2.7759d-12

  e0 = -4.6206d-13
  e1 = +1.8676d-14
  e2 = -2.1687d-16

  sw_adtg_DNAD =  a0 + (a1 + (a2 + a3*T)*T)*T &
       + (b0 + b1*T)*(S-sref) &
       + ( (c0 + (c1 + (c2 + c3*T)*T)*T) + (d0 + d1*T)*(S-sref) )*P &
       + (  e0 + (e1 + e2*T)*T )*P*P

END FUNCTION sw_adtg_DNAD
END MODULE msw_adtg
