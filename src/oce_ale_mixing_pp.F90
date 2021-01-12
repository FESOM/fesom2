!=======================================================================
subroutine oce_mixing_pp(mesh)
    !  Compute Richardson number dependent Av and Kv following
    !  Pacanowski and Philander, 1981
    !  Av = Avmax * factor**2 + Av0,
    !  Kv = Kvmax * factor**3 + Kv0, 
    !  factor=1/(1+5Ri),  Ri=N**2/(dU/dz)**2 is the Richardson number
    !                     N is the buoyancy frequency
    !                     dU/dz is the vertical velocity shear
    !  Avmax, Kvmax are tunable
    
    ! Output: Av(2:nlevels(elem)-1,:)  == vert. visc. coeff. at zbar(2:...)
    !         Kv(2:nlevels_nod2D(node),:) == vert. diff. coeff. at zbar(2:...)   
    ! NR external if for mo
    ! SD no if in Kv computations (only minor differences are introduced)
    !    
    !      
use MOD_MESH
USE o_PARAM
USE o_ARRAYS
USE g_PARSUP
USE g_config
use i_arrays
IMPLICIT NONE

type(t_mesh), intent(in) , target :: mesh
real(kind=WP)            :: dz_inv, bv, shear, a, rho_up, rho_dn, t, s, Kv0_b
integer                  :: node, nz, nzmax, nzmin, elem, elnodes(3), i

#include "associate_mesh.h"
    !___________________________________________________________________________
    do node=1, myDim_nod2D+eDim_nod2D
        nzmin = ulevels_nod2d(node)
        nzmax = nlevels_nod2d(node)
        
        !_______________________________________________________________________
        ! implement changing ALE zlevel at every node in PP mxing 
        !!PS do nz=2,nlevels_nod2d(node)-1
        do nz=nzmin+1,nzmax-1
            dz_inv=1.0_WP/(Z_3d_n(nz-1,node)-Z_3d_n(nz,node))
            shear = (Unode(1,nz-1,node)-Unode(1,nz,node))**2 +&
                    (Unode(2,nz-1,node)-Unode(2,nz,node))**2 
            shear = shear*dz_inv*dz_inv
            Kv(nz,node) = shear/(shear+5._WP*max(bvfreq(nz,node),0.0_WP)+1.0e-14)  ! To avoid NaNs at start
        end do
    end do
    
    !_______________________________________________________________________
    ! viscosity
    DO elem=1, myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        nzmin = ulevels(elem)
        nzmax = nlevels(elem)
        !!PS DO nz=2,nlevels(elem)-1
        DO nz=nzmin+1,nzmax-1
            Av(nz,elem)= mix_coeff_PP*sum(Kv(nz,elnodes)**2)/3.0_WP+A_ver
        END DO
    END DO
    
    !_______________________________________________________________________
    ! diffusivity
    do node=1,myDim_nod2D+eDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)
        !___________________________________________________________________
        ! set constant background diffusivity with namelist value K_ver
        if(Kv0_const) then
            !!PS do nz=2,nlevels_nod2d(node)-1
            do nz=nzmin+1,nzmax-1 
                Kv(nz,node) = mix_coeff_PP*Kv(nz,node)**3+K_ver
            end do
        !___________________________________________________________________________
        ! set latitudinal and depth dependent background diffusivity after 
        ! Q. Wang FESOM1.4 approach
        else
            !!PS do nz=2,nlevels_nod2d(node)-1
            do nz=nzmin+1,nzmax-1 
                call Kv0_background_qiang(Kv0_b,geo_coord_nod2D(2,node)/rad,abs(zbar_3d_n(nz,node)))
                Kv(nz,node) = mix_coeff_PP*Kv(nz,node)**3+Kv0_b
            end do
        end if
    end do

end subroutine oce_mixing_pp
!
!
!_______________________________________________________________________________
! calculate non constant background diffusion coefficient as it is also used in 
! KPP mixing scheme (see. oce_ale_mxing_kpp.F90-->subroutine ri_iwmix(viscA, diffK))
! after suggestions from Qiang Wang in FESOM1.4
subroutine Kv0_background_qiang(Kv0_b,lat,dep)
    ! Kv0_b ... backround mixing coefficient
    use o_PARAM, only: WP
    implicit none
    real(kind=WP), intent(out) :: Kv0_b
    real(kind=WP), intent(in)  :: lat, dep
    real(kind=WP)              :: aux, ratio
    
    !___________________________________________________________________________
    ! set latitudinal and depth dependent background diffusivity after 
    ! Q. Wang FESOM1.4 approach
    !_______________________________________________________________________
    aux = (0.6_WP + 1.0598_WP / 3.1415926_WP * ATAN( 4.5e-3_WP * (dep - 2500.0_WP))) * 1.0e-5_WP
    
    !_______________________________________________________________________
    ! latitudinal equatorial scaling
    if (abs(lat) < 5.0_WP) then
        ratio = 1.0_WP
    else
        ratio = MIN( 1.0_WP + 9.0_WP * (abs(lat) - 5.0_WP) / 10.0_WP, 10.0_WP )
    end if 
    
    !_______________________________________________________________________
    ! latitudinal arctic scaling
    if (lat > 70.0_WP) then
        if (dep <= 50.0_WP)     then
            ratio=4.0_WP + 6.0_WP * (50.0_WP - dep) / 50.0_WP
        else
            ratio=4.0_WP  
        endif   
    end if
    !_______________________________________________________________________
    Kv0_b = aux*ratio
    
end subroutine Kv0_background_qiang
!
!
!_______________________________________________________________________________
! calculate non constant background diffusion coefficient as it is also used in 
! KPP mixing scheme (see. oce_ale_mxing_kpp.F90-->subroutine ri_iwmix(viscA, diffK))
! first implemented my Q.Wang in FESOM1.4
subroutine Kv0_background(Kv0_b,lat,dep)
    ! Kv0_b ... backround mixing coefficient
    use o_PARAM, only: WP
    implicit none
    real(kind=WP), intent(out) :: Kv0_b
    real(kind=WP), intent(in)  :: lat, dep
    real(kind=WP)              :: aux, ratio
    
    !___________________________________________________________________________
    ! set latitudinal and depth dependent background diffusivity after 
    ! Q. Wang FESOM1.4 approach
    !_______________________________________________________________________
    ! latitudinal equatorial scaling
    if (abs(lat) < 5.0_WP) then
        ratio = 1.0_WP
    else
        ratio = MIN( 1.0_WP + 9.0_WP * (abs(lat) - 5.0_WP) / 10.0_WP, 10.0_WP )
    end if 
    
    !_______________________________________________________________________
    ! latitudinal <70° scaling
    if (lat < 70.0_WP) then
        aux = (0.6_WP + 1.0598_WP / 3.1415926_WP * ATAN( 4.5e-3_WP * (dep - 2500.0_WP))) * 1e-5_WP
    ! latitudinal arctic scaling
    else
        aux = (0.6_WP + 1.0598_WP / 3.1415926_WP * ATAN( 4.5e-3_WP * (dep - 2500.0_WP))) * 1.e-6_WP
        ratio=3.0_WP
        if (dep < 80.0_WP)     then
            ratio=1.0_WP
        elseif(dep < 100.0_WP) then
            ratio=2.0_WP  
        endif   
    end if
    !_______________________________________________________________________
    Kv0_b = aux*ratio
    
end subroutine Kv0_background
