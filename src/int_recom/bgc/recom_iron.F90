!===============================================================================
! MODULE: recom_iron
!
! Purpose:
!   Iron speciation chemistry for the REcoM biogeochemical model.
!   Contains three pure functions extracted from recom_sms.F90 (Step 3
!   restructuring).  All functions are side-effect free and depend only on
!   their arguments plus the recom_config namelist flag REcoM_Geider_limiter.
!
! Functions:
!   recom_limiter          — nutrient/quota limiter (Geider or exponential form)
!   iron_chemistry         — free-iron solver, single organic ligand (quadratic)
!   iron_chemistry_2ligands — free-iron solver, two organic ligands (cubic)
!
! Usage:
!   use recom_iron
!
! Origin:
!   Extracted verbatim from src/int_recom/recom_sms.F90 lines 7569-7661.
!   Original implementation: REcoM development team, AWI Bremerhaven.
!===============================================================================

module recom_iron

  implicit none
  private

  public :: recom_limiter
  public :: iron_chemistry
  public :: iron_chemistry_2ligands

contains

!-------------------------------------------------------------------------------
! Function for calculating limiter
!-------------------------------------------------------------------------------

function recom_limiter(slope, qa, qb) result(lim)
  use recom_config, only: REcoM_Geider_limiter
  implicit none

  real(kind=8), intent(in) :: slope, qa, qb
  real(kind=8)             :: lim
  real(kind=8)             :: dq

  dq = qa - qb
  if (REcoM_Geider_limiter) then
    lim = max(min(-slope * dq, 1.d0), 0.d0)
  else
    lim = 1.d0 - exp(-slope * (abs(dq) - dq)**2)
  end if

end function recom_limiter

!-------------------------------------------------------------------------------
! Function for iron chemistry — single organic ligand (quadratic solution)
!
! Solves for free (inorganic) iron given:
!   Fe             : Total dissolved iron          [µmol m-3]
!   totalLigand    : Total organic ligand conc.    [µmol m-3]
!   ligandStabConst: Conditional stability constant [M-1]
!
! Returns: free (inorganic) iron concentration [µmol m-3]
!-------------------------------------------------------------------------------

function iron_chemistry(Fe, totalLigand, ligandStabConst) result(freeFe)
  implicit none

  real(kind=8), intent(in) :: Fe, totalLigand, ligandStabConst
  real(kind=8)             :: freeFe
  real(kind=8)             :: ligand, FeL, a, b, c, discrim

  a = ligandStabConst
  b = ligandStabConst * (Fe - totalLigand) + 1.d0
  c = -totalLigand
  discrim = b*b - 4.d0 * a * c

  if (a .ne. 0.d0 .and. discrim .ge. 0.d0) then
    ligand = (-b + sqrt(discrim)) / (2.d0 * a)
    FeL    = totalLigand - ligand
    freeFe = Fe - FeL
  else
    freeFe = 0.d0   ! No free iron
  end if

end function iron_chemistry

!-------------------------------------------------------------------------------
! Function for iron chemistry — two organic ligands (cubic solution)
!
! Solves the 3rd-order polynomial for free iron when two ligand classes
! (L1, L2) compete for complexation.  Returns the physically meaningful
! (largest real) root via Cardano / trigonometric method.
!
! Arguments:
!   fet : Total dissolved iron     [µmol m-3]
!   l1t : Total ligand-1 conc.     [µmol m-3]
!   l2t : Total ligand-2 conc.     [µmol m-3]
!   k1  : Stability constant L1    [M-1]
!   k2  : Stability constant L2    [M-1]
!
! Returns: free iron [µmol m-3]
!-------------------------------------------------------------------------------

function iron_chemistry_2ligands(fet, l1t, l2t, k1, k2) result(free_fe)
  implicit none

  real(kind=8), intent(in) :: fet, l1t, l2t, k1, k2
  real(kind=8)             :: free_fe

  real(kind=8) :: a3, a2, a1, a0, a, b, c
  real(kind=8) :: p, q, discr, rho, phi, amp
  real(kind=8) :: one3rd, one27th
  real(kind=8) :: fe1, fe2, fe3
  real(kind=8), parameter :: pi = 3.1415926535897931d0

  ! Coefficients of the 4th-order polynomial
  a3 = k1 * k2
  a2 = k1*k2*(l1t + l2t - fet) + k1 + k2
  a1 = 1.d0 - (k1 + k2)*fet + k1*l1t + k2*l2t
  a0 = -fet

  ! Normalised polynomial coefficients
  a = a2 / a3
  b = a1 / a3
  c = a0 / a3

  one3rd  = 1.d0 / 3.d0
  one27th = 1.d0 / 27.d0

  ! Cardano / trigonometric solution
  p    = b - a*a * one3rd
  q    = c - a*b * one3rd + 2.d0*a*a*a * one27th
  discr = q*q / 4.d0 + p*p*p * one27th

  rho  = sqrt(-(p*p*p * one27th))
  phi  = acos(-q / (2.d0 * rho))
  amp  = 2.d0 * rho**one3rd

  ! Three real roots — take the largest (physically meaningful free-Fe value)
  fe1 = amp * cos(phi * one3rd)             - a * one3rd
  fe2 = amp * cos((phi + 2.d0*pi) * one3rd) - a * one3rd
  fe3 = amp * cos((phi + 4.d0*pi) * one3rd) - a * one3rd

  free_fe = max(fe1, fe2, fe3)

end function iron_chemistry_2ligands

end module recom_iron

