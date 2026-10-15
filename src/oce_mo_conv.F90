!
!
!_______________________________________________________________________________
subroutine mo_convect(ice, partit, mesh)
    USE o_PARAM
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_ICE
    USE o_ARRAYS
    USE g_config
    use g_comm_auto
    IMPLICIT NONE
    type(t_ice), intent(in), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    !___________________________________________________________________________
    integer                               :: node, elem, nz, elnodes(3), nzmin, nzmax
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice, a_ice
    u_ice => ice%uice(:)
    v_ice => ice%vice(:)
    a_ice => ice%data(1)%values(:)
    
    !___________________________________________________________________________
    ! add vertical mixing scheme of Timmermann and Beckmann, 2004,"Parameterization 
    ! of vertical mixing in the Weddell Sea!
    ! Computes the mixing length derived from the Monin
    if (use_momix) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax)
        do node=1, partit%myDim_nod2D+partit%eDim_nod2D
            mo(:, node) = 0._WP
            nzmax = mesh%nlevels_nod2d(node)
            nzmin = mesh%ulevels_nod2d(node)
            !___________________________________________________________________
            ! apply Monin Obukov only below certain latitude momix_lat
            if (mesh%geo_coord_nod2D(2,node)>momix_lat*rad) cycle
            
            ! apply Monin Obukov only when there is no cavity
            if (nzmin>1) cycle
            
            !___________________________________________________________________
            ! calcualte monin obukhov length
!!PS             call mo_length(water_flux(node),heat_flux(node), &         
!!PS                     stress_atmoce_x(node),stress_atmoce_y(node), &    
!!PS                     u_ice(node),v_ice(node),a_ice(node), &                             
!!PS                     dt, mixlength(node))
            call mo_length(water_flux(node),heat_flux(node), &         
                    stress_node_surf(1,node),stress_node_surf(2,node), &    
                    u_ice(node),v_ice(node),a_ice(node), &                             
                    dt, mixlength(node))        
            
            !___________________________________________________________________
            ! increase vertical diffusion within monin obukov length to namelist
            ! parameter momix_kv. momix_kv in moment set to 0.01 --> that means 
            ! very strong vertical mixing within mixlength
            !!PS do nz = 2,nlevels_nod2d(node)-1
            do nz = nzmin+1,nzmax-1
                if(abs(mesh%zbar_3d_n(nz,node)) <= mixlength(node)) then
                    ! write in vertice array to also apply to viscosity
                    mo(nz,node)=momix_kv    ! Potentialy bad place 
                    
                    ! apply monin-obukov mixing to diffusivity
                    Kv(nz,node)=Kv(nz,node)+mo(nz,node)
                end if    
            end do 
        end do
!$OMP END PARALLEL DO
    end if
    !
    !___________________________________________________________________________
    ! apply mixing enhancements to vertical diffusivity
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax)
    do node=1, partit%myDim_nod2D+partit%eDim_nod2D
        nzmax = mesh%nlevels_nod2d(node)
        nzmin = mesh%ulevels_nod2d(node)
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
!$OMP END PARALLEL DO
    !
    !___________________________________________________________________________
    ! apply mixing enhancements to vertical viscosity
    ! elem2D_nodes has no dimension until +eDim_elem2D
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, elnodes, nz, nzmin, nzmax)
    do elem=1, partit%myDim_elem2D
        elnodes=mesh%elem2D_nodes(:,elem)
        nzmax = mesh%nlevels(elem)
        nzmin = mesh%ulevels(elem)
        !!PS do nz=2,nlevels(elem)-1
        do nz=nzmin+1,nzmax-1
            
            ! force on convection if unstable 
            if (use_instabmix .and. any(bvfreq(nz, elnodes) < 0._WP)) Av(nz,elem)=max(Av(nz,elem),instabmix_kv)
            !!PS if (any(bvfreq(nz, elnodes) < 0._WP)) Av(nz,elem)=av_conv ! fesom1.4 style
            
            ! if cavity dont apply any monin obukov or wind mixing 
            if (nzmin>1) cycle
                ! apply monin-obukov mixing to viscosity
                if (use_momix .and. sum(mesh%geo_coord_nod2D(2,elnodes))/3.0<=momix_lat*rad) Av(nz,elem)=Av(nz,elem)+sum(mo(nz,elnodes))/3.0_WP
                
                ! force enhanced wind mixing --> from original fesom1.4 pp mixing
                if (use_windmix .and. nz<=windmix_nl+1) Av(nz,elem)=max(Av(nz,elem), windmix_kv)
        end do
    end do
!$OMP END PARALLEL DO
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
#ifdef ENABLE_OPENACC
!
!
!===============================================================================
! The Kv copy after KPP and mo_convect on the device (Kv_double and Av = viscAE current there); mo_length
! and pmlktmo are inlined with the same arithmetic. Kv and Av are copied back for the host code.
subroutine mo_convect_acc(ice, partit, mesh)
    USE o_PARAM
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_ICE
    USE o_ARRAYS
    USE g_config
    use g_comm_auto
    IMPLICIT NONE
    type(t_ice), intent(in), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    integer                               :: node, elem, nz, n1, n2, n3, nzmin, nzmax, iter
    integer                               :: myDim_nod2D, nod_size, myDim_elem2D, nl
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice, a_ice
    logical, save                         :: data_created = .false.
    real(kind=WP)             :: uabs, ret, rtc, qfm, qtm, qw, ustar, tau, obuk, mlen
    real(kind=WP), parameter  :: cosgam = 0.913632  ! cos(24.*3.14/180.), as in mo_length
    real(kind=WP), parameter  :: qhw   = 1/7.0      ! as in pmlktmo
    real(kind=WP), parameter  :: betas = 0.0008
    real(kind=WP), parameter  :: betat = 0.00004
    real(kind=WP)             :: a1, f0, f1, ttmp, qrho, l_dt

    u_ice => ice%uice(:)
    v_ice => ice%vice(:)
    a_ice => ice%data(1)%values(:)
    myDim_nod2D  = partit%myDim_nod2D
    nod_size     = partit%myDim_nod2D+partit%eDim_nod2D
    myDim_elem2D = partit%myDim_elem2D
    nl           = mesh%nl
    l_dt         = dt
    if (.not. data_created) then
        !$ACC ENTER DATA CREATE(Kv) COPYIN(mesh%geo_coord_nod2D)
        if (use_momix) then
            !$ACC ENTER DATA COPYIN(mo, mixlength)
        end if
        data_created = .true.
    end if
    if (use_momix) then
        ! computed on the host (ice step)
        !$ACC UPDATE DEVICE(u_ice, v_ice, a_ice)
    end if

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(nz)
    do node=1, myDim_nod2D
        do nz=1, nl
            Kv(nz, node)=Kv_double(nz, node, 1)
        end do
    end do
    !$ACC END PARALLEL LOOP
    call exchange_nod(Kv, partit, luse_g2g=.true.)

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) &
    !$ACC PRIVATE(nz, nzmin, nzmax, uabs, ret, rtc, qfm, qtm, qw, ustar, tau, obuk, mlen, a1, f0, f1, ttmp, qrho, iter)
    do node=1, nod_size
        nzmax = mesh%nlevels_nod2d(node)
        nzmin = mesh%ulevels_nod2d(node)
        if (use_momix) then
            do nz=1, nl
                mo(nz, node) = 0._WP
            end do
            if (mesh%geo_coord_nod2D(2,node)<=momix_lat*rad .and. nzmin<=1) then
                ! mo_length
                qfm   = water_flux(node) * 34._WP
                qtm   = - 2.38e-7_WP * heat_flux(node)
                tau   = sqrt(stress_node_surf(1,node)**2+stress_node_surf(2,node)**2)
                ustar = sqrt(tau/1030._WP)
                uabs  = sqrt(u_ice(node)**2+v_ice(node)**2)
                qw    = 1.25_WP * ustar**3 * (1._WP-a_ice(node)) + 0.005_WP * uabs**3 * cosgam * a_ice(node)
                ! pmlktmo
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
                rtc=l_dt/(10.0_WP * 86400.0_WP)
                mlen = mixlength(node)
                if (obuk.lt.mlen) then
                   ret = (obuk-mlen)*rtc
                   mlen = mlen+ret
                else
                   mlen=obuk
                endif
                mixlength(node) = mlen
                do nz = nzmin+1,nzmax-1
                    if(abs(mesh%zbar_3d_n(nz,node)) <= mixlength(node)) then
                        mo(nz,node)=momix_kv
                        Kv(nz,node)=Kv(nz,node)+mo(nz,node)
                    end if
                end do
            end if
        end if
        do nz=nzmin+1, nzmax-1
            if (use_instabmix .and. bvfreq(nz, node) < 0._WP)  Kv(nz,node)=max(Kv(nz,node),instabmix_kv)
            if (nzmin>1) cycle
            if (use_windmix .and. nz<=windmix_nl+1) Kv(nz,node)=max(Kv(nz,node), windmix_kv)
        end do
    end do
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(nz, nzmin, nzmax, n1, n2, n3)
    do elem=1, myDim_elem2D
        n1 = mesh%elem2D_nodes(1,elem)
        n2 = mesh%elem2D_nodes(2,elem)
        n3 = mesh%elem2D_nodes(3,elem)
        nzmax = mesh%nlevels(elem)
        nzmin = mesh%ulevels(elem)
        do nz=nzmin+1,nzmax-1
            if (use_instabmix .and. (bvfreq(nz,n1) < 0._WP .or. bvfreq(nz,n2) < 0._WP .or. bvfreq(nz,n3) < 0._WP)) &
                Av(nz,elem)=max(Av(nz,elem),instabmix_kv)
            if (nzmin>1) cycle
            if (use_momix) then
                if ((mesh%geo_coord_nod2D(2,n1)+mesh%geo_coord_nod2D(2,n2)+mesh%geo_coord_nod2D(2,n3))/3.0<=momix_lat*rad) &
                    Av(nz,elem)=Av(nz,elem)+(mo(nz,n1)+mo(nz,n2)+mo(nz,n3))/3.0_WP
            end if
            if (use_windmix .and. nz<=windmix_nl+1) Av(nz,elem)=max(Av(nz,elem), windmix_kv)
        end do
    end do
    !$ACC END PARALLEL LOOP
    !$ACC UPDATE HOST(Kv, Av)
end subroutine mo_convect_acc
#endif
