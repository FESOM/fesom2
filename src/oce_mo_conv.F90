!
!
!_______________________________________________________________________________
subroutine mo_convect(mesh)
    USE o_PARAM
    USE MOD_MESH
    USE o_ARRAYS
    USE g_PARSUP
    USE g_config
    use i_arrays
    use g_comm_auto
    IMPLICIT NONE

    integer                           :: node, elem, nz, elnodes(3), nzmin, nzmax
    type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"
    !___________________________________________________________________________
    ! add vertical mixing scheme of Timmermann and Beckmann, 2004,"Parameterization 
    ! of vertical mixing in the Weddell Sea!
    ! Computes the mixing length derived from the Monin
    if (use_momix) then
        mo = 0._WP
        do node=1, myDim_nod2D+eDim_nod2D
            nzmax = nlevels_nod2d(node)
            nzmin = ulevels_nod2d(node)
            !___________________________________________________________________
            ! apply Monin Obukov only below certain latitude momix_lat
            if (geo_coord_nod2D(2,node)>momix_lat*rad) cycle
            
            ! apply Monin Obukov only when there is no cavity
            if (nzmin>1) cycle
            
            !___________________________________________________________________
            ! calcualte monin obukhov length
            call mo_length(water_flux(node),heat_flux(node), &         
                    stress_atmoce_x(node),stress_atmoce_y(node), &    
                    u_ice(node),v_ice(node),a_ice(node), &                             
                    dt, mixlength(node))
            
            !___________________________________________________________________
            ! increase vertical diffusion within monin obukov length to namelist
            ! parameter momix_kv. momix_kv in moment set to 0.01 --> that means 
            ! very strong vertical mixing within mixlength
            !!PS do nz = 2,nlevels_nod2d(node)-1
            do nz = nzmin+1,nzmax-1
                if(abs(zbar_3d_n(nz,node)) <= mixlength(node)) then
                    ! write in vertice array to also apply to viscosity
                    mo(nz,node)=momix_kv    ! Potentialy bad place 
                    
                    ! apply monin-obukov mixing to diffusivity
                    Kv(nz,node)=Kv(nz,node)+mo(nz,node)
                end if    
            end do 
        end do
    end if
    
    !
    !___________________________________________________________________________
    ! apply mixing enhancements to vertical diffusivity
    do node=1, myDim_nod2D+eDim_nod2D
        nzmax = nlevels_nod2d(node)
        nzmin = ulevels_nod2d(node)
        !!PS do nz=2, nlevels_nod2d(node)-1
        do nz=nzmin+1, nzmax-1
            ! force on convection if unstable 
            if (use_instabmix .and. bvfreq(nz, node) < 0._WP)  Kv(nz,node)=max(Kv(nz,node),instabmix_kv)
            !!PS if (bvfreq(nz, node) < 0._WP)  Kv(nz,node)=kv_conv ! fesom1.4 style
            
            if (nzmin>1) cycle
                ! force enhanced wind mixing --> from original fesom1.4 pp mixing
                if (use_windmix .and. nz<=windmix_nl+1) Kv(nz,node)=max(Kv(nz,node), windmix_kv)
            
        end do
    end do
    
    !
    !___________________________________________________________________________
    ! apply mixing enhancements to vertical viscosity
    ! elem2D_nodes has no dimension until +eDim_elem2D
    do elem=1, myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        nzmax = nlevels(elem)
        nzmin = ulevels(elem)
        !!PS do nz=2,nlevels(elem)-1
        do nz=nzmin+1,nzmax-1
            
            ! force on convection if unstable 
            if (use_instabmix .and. any(bvfreq(nz, elnodes) < 0._WP)) Av(nz,elem)=max(Av(nz,elem),instabmix_kv)
            !!PS if (any(bvfreq(nz, elnodes) < 0._WP)) Av(nz,elem)=av_conv ! fesom1.4 style
            
            ! if cavity dont apply any monin obukov or wind mixing 
            if (nzmin>1) cycle
                ! apply monin-obukov mixing to viscosity
                if (use_momix .and. sum(geo_coord_nod2D(2,elnodes))/3.0<=momix_lat*rad) Av(nz,elem)=Av(nz,elem)+sum(mo(nz,elnodes))/3.0_WP
                
                ! force enhanced wind mixing --> from original fesom1.4 pp mixing
                if (use_windmix .and. nz<=windmix_nl+1) Av(nz,elem)=max(Av(nz,elem), windmix_kv)
        end do
    end do
    !!PS call exchange_elem(Av)
end subroutine mo_convect
!
!
!_______________________________________________________________________________
subroutine mo_length(water_flux,heat_flux,stress_x,stress_y,  &
       u_ice,v_ice,a_ice,dt,mixlength)
    ! vertical mixing scheme of Timmermann and Beckmann, 2004.
    ! computes the mixing length derived from the Monin-Obukhov length
    ! Ralph Timmermann, 14.06.2006
  
  USE o_PARAM, only: WP
    implicit none

    real(kind=WP)             :: water_flux, heat_flux, stress_x, stress_y
    real(kind=WP)             :: u_ice, v_ice, a_ice, uabs
    real(kind=WP)             :: dt, ret, rtc, mixlength
    real(kind=WP)             :: qfm, qtm, qw
    real(kind=WP)             :: ustar,tau, obuk
    real(kind=WP), parameter :: cosgam = 0.913632  ! cos(24.*3.14/180.)

    qfm            = water_flux * 34._WP       ! note that water_flux>0
    ![psu * m/s]   [m/s]   [psu]    ! f. upward fresh water flux

    qtm            = - 2.38e-7_WP * heat_flux  ! heat_flux>0 f. upward heat flux
    ![K * m/s]

    tau = sqrt(stress_x**2+stress_y**2)
    ustar = sqrt(tau/1030._WP)
    uabs = sqrt(u_ice**2+v_ice**2)

    qw = 1.25_WP * ustar**3 * (1._WP-a_ice) + 0.005_WP * uabs**3 * cosgam * a_ice  !Eq. 8 of TB04

    call pmlktmo(qfm,qtm,qw,obuk)

    rtc=dt/(10.0_WP * 86400.0_WP)     ! (NR: inverse of) time constant of mixed layer retreat

    if (obuk.lt.mixlength) then
       ret = (obuk-mixlength)*rtc
       mixlength = mixlength+ret
    else
       mixlength=obuk
    endif

end subroutine mo_length
!
!
!_______________________________________________________________________________
subroutine pmlktmo(qfm,qtm,qw,obuk)
    ! gives the Monin-Obukhov length
    ! qtm  = Heat Flux into ML                                 [K m/s]
    ! qfm  = salinity flux into ML                             [psu m/s]
    ! qw   = production of turbulent kinetic energy 
    !-----------------------------------------------------------------------   
  USE o_PARAM, only: WP
    implicit none

    integer           :: iter
    real(kind=WP)            :: qtm, qfm, qw, obuk
    real(kind=WP), parameter :: qhw   = 1/7.0                 ! [1/m] !NR inverse
    real(kind=WP), parameter :: betas = 0.0008
    real(kind=WP), parameter :: betat = 0.00004
    real(kind=WP)            :: a1, f0, f1, ttmp, qrho

    qrho=betas*qfm-betat*qtm
    ttmp=60._WP

    if(qrho>0.) then
       ttmp=0._WP
    else
       do iter=1,5
          a1 = exp(-ttmp*qhw) 
          f0 = 2._WP* qw * a1 + 9.81_WP * qrho * ttmp
          f1 =-2._WP* qw * a1 * qhw + 9.81_WP * qrho
          ttmp = ttmp - f0 / f1
          ttmp = max(ttmp,10._WP) 
       enddo
    end if

    obuk=max(ttmp,10.)

end subroutine pmlktmo
