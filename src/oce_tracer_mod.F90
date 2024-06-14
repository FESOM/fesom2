!============================================================================================
MODULE o_tracers
USE MOD_MESH
USE MOD_TRACER
USE MOD_PARTIT
USE MOD_PARSUP
IMPLICIT NONE

CONTAINS
!
!
!===============================================================================
SUBROUTINE init_tracers_AB(tr_num, tracers, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    use g_config, only: flag_debug
    use o_arrays
    use g_comm_auto
    IMPLICIT NONE
    integer,        intent(in)            :: tr_num
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(inout), target :: tracers
    integer                               :: n,nz 

!$ACC parallel loop collapse(2) default(present) !!!async(1)
do n=1, partit%myDim_nod2D+partit%eDim_nod2D
       do nz=1, mesh%nl-1
       ! del_ttf will contain all advection / diffusion contributions for this tracer. Set it to 0 at the beginning!
       tracers%work%del_ttf          (nz, n) = 0.0_WP
       tracers%work%del_ttf_advhoriz (nz, n) = 0.0_WP
       tracers%work%del_ttf_advvert  (nz, n) = 0.0_WP
       end do
end do
!$ACC end parallel loop
!$OMP PARALLEL DO
    do n=1, partit%myDim_nod2D+partit%eDim_nod2D
       ! AB interpolation
       if (tracers%data(tr_num)%AB_order==2) then
           tracers%data(tr_num)%valuesAB(:, n)  =-(0.5_WP+epsilon)*tracers%data(tr_num)%valuesold(1, :, n)+(1.5_WP+epsilon)*tracers%data(tr_num)%values(:, n)
       elseif (tracers%data(tr_num)%AB_order==3) then
           tracers%data(tr_num)%valuesAB(:, n)  =5.0_WP*tracers%data(tr_num)%valuesold(2, :, n)-16.0_WP*tracers%data(tr_num)%valuesold(1, :, n)+23.0_WP*tracers%data(tr_num)%values(:, n)
           tracers%data(tr_num)%valuesAB(:, n)  =tracers%data(tr_num)%valuesAB(:, n)/12.0_WP
       end if
    end do
!$OMP END PARALLEL DO

    if (tracers%data(tr_num)%AB_order==2) then
!$OMP PARALLEL DO
       do n=1, partit%myDim_nod2d+partit%eDim_nod2D
          tracers%data(tr_num)%valuesold(1, :, n)=tracers%data(tr_num)%values(:, n)
       end do
!$OMP END PARALLEL DO
    elseif (tracers%data(tr_num)%AB_order==3) then
!$OMP PARALLEL DO
       do n=1, partit%myDim_nod2d+partit%eDim_nod2D
          tracers%data(tr_num)%valuesold(2, :, n)=tracers%data(tr_num)%valuesold(1, :, n)
          tracers%data(tr_num)%valuesold(1, :, n)=tracers%data(tr_num)%values(:, n)
       end do
!$OMP END PARALLEL DO
    end if

    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[38m'//'             --> call tracer_gradient_elements'//achar(27)//'[0m'
    call tracer_gradient_elements(tracers%data(tr_num)%valuesAB, partit, mesh)
    call exchange_elem_begin(tr_xy, partit)

    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[38m'//'             --> call tracer_gradient_z'//achar(27)//'[0m'
    call tracer_gradient_z(tracers%data(tr_num)%values, partit, mesh)    !WHY NOT AB HERE? DSIDOREN!
    call exchange_elem_end(partit)      ! tr_xy used in fill_up_dn_grad
!$OMP BARRIER

    call exchange_nod_begin(tr_z, partit) ! not used in fill_up_dn_grad 

    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[38m'//'             --> call fill_up_dn_grad'//achar(27)//'[0m'
    call fill_up_dn_grad(tracers%work, partit, mesh)
    call exchange_nod_end(partit)       ! tr_z halos should have arrived by now.

    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[38m'//'             --> call tracer_gradient_elements'//achar(27)//'[0m'
    call tracer_gradient_elements(tracers%data(tr_num)%values, partit, mesh) !redefine tr_arr to the current timestep
    call exchange_elem(tr_xy, partit)

END SUBROUTINE init_tracers_AB
!
!
!=======================================================================
SUBROUTINE tracer_gradient_elements(ttf, partit, mesh)
    !computes elemental gradient of tracer
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    USE o_PARAM
    USE o_ARRAYS
    IMPLICIT NONE

    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP)                         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    integer                               :: elem,  elnodes(3)
    integer                               :: nz, nzmin, nzmax

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, elnodes, nz, nzmin, nzmax)
    DO elem=1, myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        nzmin = ulevels(elem)
        nzmax = nlevels(elem)
        !!PS DO nz=1, nlevels(elem)-1
        DO nz=nzmin, nzmax-1   
            tr_xy(1,nz, elem)=sum(gradient_sca(1:3,elem)*ttf(nz,elnodes))
            tr_xy(2,nz, elem)=sum(gradient_sca(4:6,elem)*ttf(nz,elnodes))
        END DO
    END DO
!$OMP END PARALLEL DO
END SUBROUTINE tracer_gradient_elements
!
!
!========================================================================================
SUBROUTINE tracer_gradient_z(ttf, partit, mesh)
    !computes vertical gradient of tracer
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    USE o_PARAM
    USE o_ARRAYS
    USE g_CONFIG
    IMPLICIT NONE
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP)            :: ttf(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP)            :: dz
    integer                  :: n, nz, nzmin, nzmax

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax, dz)
    DO n=1, myDim_nod2D+eDim_nod2D
    !!PS nlev=nlevels_nod2D(n)
    nzmax=nlevels_nod2D(n)
    nzmin=ulevels_nod2D(n)
    !!PS DO nz=2,  nlev-1
    DO nz=nzmin+1,  nzmax-1
        dz=0.5_WP*(hnode_new(nz-1,n)+hnode_new(nz,n))
        tr_z(nz, n)=(ttf(nz-1,n)-ttf(nz,n))/dz
    END DO
    !!PS tr_z(1,    n)=0.0_WP
    !!PS tr_z(nlev, n)=0.0_WP
    tr_z(nzmin, n)=0.0_WP
    tr_z(nzmax, n)=0.0_WP
    END DO
!$OMP END PARALLEL DO
END SUBROUTINE tracer_gradient_z
!
!
!========================================================================================
SUBROUTINE relax_to_clim(tr_num, tracers, partit, mesh)
    use g_config,only: dt
    use o_arrays
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    IMPLICIT NONE

    integer,        intent(in)             :: tr_num
    type(t_mesh),   intent(in),    target  :: mesh
    type(t_partit), intent(inout), target  :: partit
    type(t_tracer), intent(inout), target  :: tracers
    integer                                :: n, nzmin, nzmax    
    real(kind=WP), dimension(:,:), pointer :: trarr

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    trarr=>tracers%data(tr_num)%values(:,:)

    if ((clim_relax>1.0e-8_WP).and.(tracers%data(tr_num)%ID==1)) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nzmin, nzmax)
        DO n=1, myDim_nod2D
            nzmin = ulevels_nod2D(n)    
            nzmax = nlevels_nod2D(n)    
            !!PS tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)=tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)+&
            !!PS         relax2clim(n)*dt*(Tclim(1:nlevels_nod2D(n)-1,n)-tr_arr(1:nlevels_nod2D(n)-1,n,tr_num))
            trarr(nzmin:nzmax-1,n)=trarr(nzmin:nzmax-1,n)+&
                    relax2clim(n)*dt*(Tclim(nzmin:nzmax-1,n)-trarr(nzmin:nzmax-1,n))
        END DO
!$OMP END PARALLEL DO
    END if
    if ((clim_relax>1.0e-8_WP).and.(tracers%data(tr_num)%ID==2)) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nzmin, nzmax)
        DO n=1, myDim_nod2D
            nzmin = ulevels_nod2D(n)    
            nzmax = nlevels_nod2D(n)  
            trarr(nzmin:nzmax-1,n)=trarr(nzmin:nzmax-1,n)+&
                    relax2clim(n)*dt*(Sclim(nzmin:nzmax-1,n)-trarr(nzmin:nzmax-1,n))
        END DO
!$OMP END PARALLEL DO
    END IF 
END SUBROUTINE relax_to_clim
END MODULE o_tracers
