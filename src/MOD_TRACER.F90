!==========================================================
MODULE MOD_TRACER
USE O_PARAM
IMPLICIT NONE
SAVE

TYPE T_TRACER
real(kind=WP), allocatable, dimension(:,:)  :: values, valuesAB !instant values & Adams-Bashfort interpolation
logical                                     :: smooth_bh_tra=.false.
real(kind=WP)                               :: gamma0_tra, gamma1_tra, gamma2_tra
logical                                     :: i_vert_diff =.false.
character(20)                               :: tra_adv_hor, tra_adv_ver, tra_adv_lim ! type of the advection scheme for this tracer
real(kind=WP)                               :: tra_adv_ph  = 1.  ! a parameter to be used in horizontal advection (for MUSCL it is the fraction of fourth-order contribution in the solution)
real(kind=WP)                               :: tra_adv_pv  = 1.  ! a parameter to be used in horizontal advection (for QR4C  it is the fraction of fourth-order contribution in the solution)
integer                                     :: ID
END TYPE T_TRACER
end module MOD_TRACER
!==========================================================

