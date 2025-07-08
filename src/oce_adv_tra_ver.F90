!===============================================================================================================================
!**************** routines for vertical tracer advection ***********************
module oce_adv_tra_ver_module
  use MOD_MESH
  use MOD_TRACER
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_comm_auto
  implicit none
  private
  public :: adv_tra_vert_impl, adv_tra_ver_upw1, adv_tra_ver_qr4c, adv_tra_vert_ppm, adv_tra_ver_cdiff
contains
! implicit 1st order upwind vertical advection with to solve for fct_LO
! updates the input tracer ttf
  subroutine adv_tra_vert_impl(dt, w, ttf, partit, mesh)
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    real(kind=WP), intent(in), target  :: dt
    type(t_partit),intent(in), target  :: partit
    type(t_mesh),  intent(in), target  :: mesh
    real(kind=WP), intent(inout)       :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)          :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
  end subroutine adv_tra_vert_impl
!===============================================================================
! 1st order upwind (explicit)
! returns flux given at vertical interfaces of scalar volumes
! IF o_init_zero=.TRUE.  : flux will be set to zero before computation
! IF o_init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
  subroutine adv_tra_ver_upw1(w, ttf, partit, mesh, flux, o_init_zero)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    type(t_partit),intent(in), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
    logical, optional                 :: o_init_zero
  end subroutine adv_tra_ver_upw1
!===============================================================================
! QR (4th order centerd)
! returns flux given at vertical interfaces of scalar volumes
! IF o_init_zero=.TRUE.  : flux will be set to zero before computation
! IF o_init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
  subroutine adv_tra_ver_qr4c(w, ttf, partit, mesh, num_ord, flux, o_init_zero)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    type(t_partit),intent(in), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
    logical, optional                 :: o_init_zero
  end subroutine adv_tra_ver_qr4c
!===============================================================================
! Vertical advection with PPM reconstruction (5th order)
! returns flux given at vertical interfaces of scalar volumes
! IF o_init_zero=.TRUE.  : flux will be set to zero before computation
! IF o_init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
 subroutine adv_tra_vert_ppm(dt, w, ttf, partit, mesh, flux, o_init_zero)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    real(kind=WP), intent(in), target :: dt
    type(t_partit),intent(in), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP)                     :: tvert(mesh%nl), tv
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
    logical, optional                 :: o_init_zero
  end subroutine adv_tra_vert_ppm
! central difference reconstruction (2nd order, use only with FCT)
! returns flux given at vertical interfaces of scalar volumes
! IF o_init_zero=.TRUE.  : flux will be set to zero before computation
! IF o_init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
  subroutine adv_tra_ver_cdiff(w, ttf, partit, mesh, flux, o_init_zero)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    type(t_partit),intent(in), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    integer                           :: n, nz, nl1
    real(kind=WP)                     :: tvert(mesh%nl), tv
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: W  (mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl,  partit%myDim_nod2D)
    logical, optional                 :: o_init_zero
  end subroutine adv_tra_ver_cdiff
end module oce_adv_tra_ver_module
