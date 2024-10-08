!> \file sw_ptmp.f90
!! \BRIEF 
!> Module with sw_ptmp function - compute potential T from in-situ T
MODULE msw_ptmp
CONTAINS
!> Function to calculate potential temperature [C] from in-situ temperature
FUNCTION sw_ptmp  (s,t,p,pr)

  !     ==================================================================
  !     Calculates potential temperature [C] from in-situ Temperature [C]
  !     From UNESCO 1983 report.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = temperature         [degree C (IPTS-68)]
  !     p  = pressure            [db]
  !     pr = reference pressure  [db]

  USE msingledouble
  USE msw_adtg
  IMPLICIT NONE

! Input arguments
  !> salinity [psu (PSS-78)]
  REAL(kind=r8) :: s
  !> temperature [degree C (IPTS-68)]
  REAL(kind=r8) :: t
  !> pressure [db]
  REAL(kind=r8) :: p
  !> reference pressure  [db]  
  REAL(kind=r8) :: pr

! local arguments
  REAL(kind=r8) :: del_P ,del_th, th, q
  REAL(kind=r8) :: onehalf, two, three
  PARAMETER (onehalf = 0.5d0, two = 2.d0, three = 3.d0 )

! REAL(kind=r8) :: sw_adtg
! EXTERNAL sw_adtg

! Output 
  REAL(kind=r8) :: sw_ptmp

  ! theta1
  del_P  = PR - P
  del_th = del_P*sw_adtg(S,T,P)
  th     = T + onehalf*del_th
  q      = del_th

  ! theta2
  del_th = del_P*sw_adtg(S,th,P+onehalf*del_P)
  th     = th + (1.d0 - 1.d0/SQRT(two))*(del_th - q)
  q      = (two-SQRT(two))*del_th + (-two+three/SQRT(two))*q

  ! theta3
  del_th = del_P*sw_adtg(S,th,P+onehalf*del_P)
  th     = th + (1.d0 + 1.d0/SQRT(two))*(del_th - q)
  q      = (two + SQRT(two))*del_th + (-two-three/SQRT(two))*q

  ! theta4
  del_th = del_P*sw_adtg(S,th,P+del_P)
  sw_ptmp     = th + (del_th - two*q)/(two*three)

  RETURN
END FUNCTION sw_ptmp


!> Function to calculate potential temperature [C] from in-situ temperature
!! and derivative with respect to insitu temperature and salinity
FUNCTION sw_ptmp_DNAD  (s,t,p,pr)

  !     It is similar to subroutine 'sw_ptmp' above except that it also computes
  !     partial derivative of potential temperature
  !     with respect to insitu temperature and salinity.
  !     ==================================================================
  !     Calculates potential temperature [C] from in-situ Temperature [C]
  !     From UNESCO 1983 report.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = temperature         [degree C (IPTS-68)]
  !     p  = pressure            [db]
  !     pr = reference pressure  [db]

  USE msingledouble
  USE msw_adtg
  USE Dual_Num_Auto_Diff
  IMPLICIT NONE

! Input arguments
  !> salinity [psu (PSS-78)]
  TYPE(DUAL_NUM) :: s
  !> temperature [degree C (IPTS-68)]
  TYPE(DUAL_NUM) :: t
  !> pressure [db]
  TYPE(DUAL_NUM) :: p
  !> reference pressure  [db]  
  TYPE(DUAL_NUM) :: pr

! local arguments
  TYPE(DUAL_NUM) :: del_P ,del_th, th, q
  REAL(kind=r8) :: onehalf, two, three
  PARAMETER (onehalf = 0.5d0, two = 2.d0, three = 3.d0 )

! Output 
  TYPE(DUAL_NUM) :: sw_ptmp_DNAD

  ! theta1
  del_P  = PR - P
  del_th = del_P*sw_adtg_DNAD(S,T,P)
  th     = T + onehalf*del_th
  q      = del_th

  ! theta2
  del_th = del_P*sw_adtg_DNAD(S,th,P+onehalf*del_P)
  th     = th + (1.d0 - 1.d0/SQRT(two))*(del_th - q)
  q      = (two-SQRT(two))*del_th + (-two+three/SQRT(two))*q

  ! theta3
  del_th = del_P*sw_adtg_DNAD(S,th,P+onehalf*del_P)
  th     = th + (1.d0 + 1.d0/SQRT(two))*(del_th - q)
  q      = (two + SQRT(two))*del_th + (-two-three/SQRT(two))*q

  ! theta4
  del_th = del_P*sw_adtg_DNAD(S,th,P+del_P)
  sw_ptmp_DNAD  = th + (del_th - two*q)/(two*three)

  RETURN
END FUNCTION sw_ptmp_DNAD
END MODULE msw_ptmp
