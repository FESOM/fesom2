!============================================================================================
MODULE o_tracers
USE MOD_MESH
IMPLICIT NONE

interface
  subroutine tracer_gradient_z(ttf, mesh)
    use g_PARSUP, only: myDim_nod2D, eDim_nod2D
    use mod_mesh
    type(t_mesh), intent(in)  , target :: mesh
    real(kind=WP)                      :: ttf(mesh%nl-1,myDim_nod2D+eDim_nod2D)
  end subroutine
end interface

CONTAINS
!
!
!=======================================================================
SUBROUTINE tracer_gradient_elements(ttf, mesh)
    !computes elemental gradient of tracer 

    USE o_PARAM
    USE o_ARRAYS
    USE g_PARSUP
    IMPLICIT NONE

    type(t_mesh), intent(in) , target :: mesh
    real(kind=WP)            :: ttf(mesh%nl-1,myDim_nod2D+eDim_nod2D)
    integer                  :: elem,  elnodes(3)
    integer                  :: n, nz, nzmin, nzmax


#include  "associate_mesh.h"

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
END SUBROUTINE tracer_gradient_elements
!
!
!========================================================================================
SUBROUTINE init_tracers_AB(tr_num, mesh)
    use g_config, only: flag_debug
    use g_parsup
    use o_arrays
    use g_comm_auto
    use mod_mesh

    IMPLICIT NONE
    integer                    :: tr_num,n,nz 
    type(t_mesh), intent(in)   , target :: mesh

    !filling work arrays
    del_ttf=0.0_WP

    !AB interpolation
    tr_arr_old(:,:,tr_num)=-(0.5_WP+epsilon)*tr_arr_old(:,:,tr_num)+(1.5_WP+epsilon)*tr_arr(:,:,tr_num)

    if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call tracer_gradient_elements'//achar(27)//'[0m'
    call tracer_gradient_elements(tr_arr_old(:,:,tr_num), mesh)
    call exchange_elem_begin(tr_xy)

    if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call tracer_gradient_z'//achar(27)//'[0m'
    call tracer_gradient_z(tr_arr(:,:,tr_num), mesh)
    call exchange_elem_end()      ! tr_xy used in fill_up_dn_grad
    call exchange_nod_begin(tr_z) ! not used in fill_up_dn_grad 

    if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call fill_up_dn_grad'//achar(27)//'[0m'
    call fill_up_dn_grad(mesh)
    call exchange_nod_end()       ! tr_z halos should have arrived by now.

    if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call tracer_gradient_elements'//achar(27)//'[0m'
    call tracer_gradient_elements(tr_arr(:,:,tr_num), mesh) !redefine tr_arr to the current timestep
    call exchange_elem(tr_xy)

END SUBROUTINE init_tracers_AB
!
!
!========================================================================================
SUBROUTINE relax_to_clim(tr_num, mesh)

    use g_config,only: dt
    USE g_PARSUP
    use o_arrays
    IMPLICIT NONE

    type(t_mesh), intent(in) , target :: mesh
    integer                  :: tr_num,n,nz, nzmin, nzmax

#include  "associate_mesh.h"

    if ((clim_relax>1.0e-8_WP).and.(tr_num==1)) then
        DO n=1, myDim_nod2D
            nzmin = ulevels_nod2D(n)    
            nzmax = nlevels_nod2D(n)    
            !!PS tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)=tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)+&
            !!PS         relax2clim(n)*dt*(Tclim(1:nlevels_nod2D(n)-1,n)-tr_arr(1:nlevels_nod2D(n)-1,n,tr_num))
            tr_arr(nzmin:nzmax-1,n,tr_num)=tr_arr(nzmin:nzmax-1,n,tr_num)+&
                    relax2clim(n)*dt*(Tclim(nzmin:nzmax-1,n)-tr_arr(nzmin:nzmax-1,n,tr_num))
        END DO
    END if
    if ((clim_relax>1.0e-8_WP).and.(tr_num==2)) then
        DO n=1, myDim_nod2D
            nzmin = ulevels_nod2D(n)    
            nzmax = nlevels_nod2D(n)  
            !!PS tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)=tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)+&
            !!PS         relax2clim(n)*dt*(Sclim(1:nlevels_nod2D(n)-1,n)-tr_arr(1:nlevels_nod2D(n)-1,n,tr_num))
            tr_arr(nzmin:nzmax-1,n,tr_num)=tr_arr(nzmin:nzmax-1,n,tr_num)+&
                    relax2clim(n)*dt*(Sclim(nzmin:nzmax-1,n)-tr_arr(nzmin:nzmax-1,n,tr_num))
        END DO
    END IF 
END SUBROUTINE relax_to_clim
END MODULE o_tracers
!
!
!========================================================================================
SUBROUTINE tracer_gradient_z(ttf, mesh)
    !computes vertical gradient of tracer
    USE o_PARAM
    USE MOD_MESH
    USE o_ARRAYS
    USE g_PARSUP
    USE g_CONFIG
    IMPLICIT NONE
    type(t_mesh), intent(in) , target :: mesh
    real(kind=WP)            :: ttf(mesh%nl-1,myDim_nod2D+eDim_nod2D)
    real(kind=WP)            :: dz
    integer                  :: n, nz, nzmin, nzmax

#include  "associate_mesh.h"

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
END SUBROUTINE tracer_gradient_z
