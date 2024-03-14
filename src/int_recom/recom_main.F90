! 24.03.2023
! OG
!===============================================================================
! Main REcoM 
module recom_interface
    interface
        subroutine recom(ice, dynamics, tracers, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        use MOD_DYN
        use MOD_ICE
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_ice)   , intent(inout), target :: ice
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine
    end interface
end module

subroutine recom(ice, dynamics, tracers, partit, mesh)
    use g_config
    use MOD_MESH
    use MOD_TRACER
    use MOD_DYN
    USE MOD_ICE
    use o_ARRAYS
    use o_PARAM
    USE MOD_PARTIT
    USE MOD_PARSUP

    use recom_declarations
    use recom_locvar
    use recom_glovar
    use recom_config
    use recom_ciso
    use g_clock
    use g_forcing_arrays, only: press_air, u_wind, v_wind, shortwave
    use g_comm_auto
    IMPLICIT NONE

    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    type(t_ice)   , intent(inout), target :: ice
    !___________________________________________________________________________

    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice

! ======================================================================================
!! Depth information

!! zbar(nl) allocate the array for storing the standard depths, it is negativ
!! Z(nl-1)  mid-depths of cells

!! max. number of levels at node n
!! nzmax = nlevels_nod2D(n)
!! u_ice and v_ice are at nodes
!! u_w, v_w are at nodes (interpolated from elements)
!! u_wind and v_wind are always at nodes
! ======================================================================================


end subroutine recom

! ======================================================================================
! Alkalinity restoring to climatology                                 	     
! ======================================================================================
subroutine bio_fluxes(tracers, partit, mesh)

    use recom_declarations
    use recom_locvar
    use recom_glovar
    use recom_config

    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use mod_tracer

    use g_config
    use o_arrays
    use g_comm_auto
    use g_forcing_arrays
    use g_support

    implicit none
    integer                               :: n, elem, elnodes(3),n1
    real(kind=WP)                         :: ralk, net

    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh

    !___________________________________________________________________________
    real(kind=WP), dimension(:,:), pointer :: alkalinity

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    alkalinity => tracers%data(2+ialk)%values(:,:) ! 1 temp, 2 salt, 3 din, 4 dic, 5 alk
!___________________________________________________________________
! on freshwater inflow/outflow or virtual alkalinity:
  ! 1. In zlevel & zstar the freshwater flux is applied in the update of the 
  ! ssh matrix when solving the continuity equation of vertically 
  ! integrated flow. The alkalinity concentration in the first layer will 
  ! be then adjusted according to the change in volume.

  ! In this case ralk is forced to be zero by setting ref_alk=0. and ref_alk_local=.false.

  ! 2. In cases where the volume of the upper layer is fixed (i.e. linfs)  the freshwater flux 
  ! 'ralk*water_flux(n)' is applied as a virtual alkalinity boundary condition via the vertical 
  ! diffusion operator.

  ! --> ralk*water_flux(n) : virtual alkalinity flux 
  ! virtual alkalinity flux

!  if (use_virt_alk) then ! OG in case of virtual alkalinity flux
!     ralk=ref_alk
!     do n=1, myDim_nod2D+eDim_nod2D
!        if (ref_alk_local) ralk = tracers%data(2+ialk)%values(1, n)
!        virtual_alk(n)=ralk*water_flux(n) 
!     end do
!  end if

    !___________________________________________________________________________
    ! Balance alkalinity restoring to climatology
!    if (.not. restore_alkalinity) return
    do n=1, myDim_nod2D+eDim_nod2D
!     relax_alk(n)=surf_relax_Alk*(Alk_surf(n)-tr_arr(1,n,2+ialk)) ! 1 temp, 2 salt
!        relax_alk(n)=surf_relax_Alk * (Alk_surf(n) - tracers%data(2+ialk)%values(1, n)) 
!        relax_alk(n)=surf_relax_Alk * (Alk_surf(n) - alkalinity(ulevels_nod2d(n),n)
        relax_alk(n)=surf_relax_Alk * (Alk_surf(n) - alkalinity(1, n))
    end do

!  if (mype==0) then
!      write(*,*) 'Alk_surf  = ', Alk_surf
!      write(*,*) 'tr_arr    = ', tracers%data(2+ialk)%values(1, :)
!      write(*,*) 'relax_alk = ', relax_alk
!  endif


  ! 2. virtual alkalinity flux
!  if (use_virt_alk) then ! is already zero otherwise
!     call integrate_nod(virtual_alk, net, partit, mesh)
!     virtual_alk=virtual_alk-net/ocean_area
!  end if


  ! 3. restoring to Alkalinity climatology
!  call integrate_nod(relax_alk, net, mesh)
    call integrate_nod(relax_alk, net, partit, mesh)
!  if (mype==0) then
!      write(*,*) 'ocean_area     = ', ocean_area
!      write(*,*) 'net            = ', net
!      write(*,*) 'net/ocean_area = ', net/ocean_area
!  endif

    relax_alk=relax_alk-net/ocean_area  ! at ocean surface layer


!  if (mype==0) then
!     write(*,*) '____________________________________________________________'
!     write(*,*) ' --> relax_alk,  = ', relax_alk
!  endif

end subroutine bio_fluxes
