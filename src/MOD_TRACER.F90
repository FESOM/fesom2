!==========================================================
MODULE MOD_TRACER
USE O_PARAM
USE, intrinsic :: ISO_FORTRAN_ENV
IMPLICIT NONE
SAVE

TYPE T_TRACER
real(kind=WP), allocatable, dimension(:,:)  :: values, valuesAB  ! instant values & Adams-Bashfort interpolation
logical                                     :: smooth_bh_tra=.false.
real(kind=WP)                               :: gamma0_tra, gamma1_tra, gamma2_tra
logical                                     :: i_vert_diff =.false.
character(20)                               :: tra_adv_hor, tra_adv_ver, tra_adv_lim ! type of the advection scheme for this tracer
real(kind=WP)                               :: tra_adv_ph  = 1.  ! a parameter to be used in horizontal advection (for MUSCL it is the fraction of fourth-order contribution in the solution)
real(kind=WP)                               :: tra_adv_pv  = 1.  ! a parameter to be used in horizontal advection (for QR4C  it is the fraction of fourth-order contribution in the solution)
integer                                     :: ID
END TYPE T_TRACER

!auxuary arrays to work with tracers:
real(kind=WP), allocatable         :: del_ttf(:,:)
real(kind=WP), allocatable         :: del_ttf_advhoriz(:,:),del_ttf_advvert(:,:)
!_______________________________________________________________________________
! in case ldiag_DVD=.true. --> calculate discrete variance decay (DVD)
real(kind=WP), allocatable                    :: tr_dvd_horiz(:,:,:), tr_dvd_vert(:,:,:)
! The fct part
real(kind=WP),allocatable,dimension(:,:)      :: fct_LO          ! Low-order solution
real(kind=WP),allocatable,dimension(:,:)      :: adv_flux_hor    ! Antidif. horiz. contrib. from edges / backup for iterafive fct scheme
real(kind=WP),allocatable,dimension(:,:)      :: adv_flux_ver    ! Antidif. vert. fluxes from nodes    / backup for iterafive fct scheme

real(kind=WP),allocatable,dimension(:,:)      :: fct_ttf_max,fct_ttf_min
real(kind=WP),allocatable,dimension(:,:)      :: fct_plus,fct_minus
! MUSCL type reconstruction
integer,allocatable,dimension(:)              :: nn_num, nboundary_lay
integer,allocatable,dimension(:,:)            :: nn_pos
integer,allocatable,dimension(:,:)            :: edge_up_dn_tri
real(kind=WP),allocatable,dimension(:,:,:)    :: edge_up_dn_grad
end module MOD_TRACER
!==========================================================

