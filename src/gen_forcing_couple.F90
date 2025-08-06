module force_flux_consv_interface
  interface
    subroutine force_flux_consv(field2d, mask, n, h, do_stats, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_mesh),   intent(in),    target :: mesh
      type(t_partit), intent(inout), target :: partit
      real(kind=WP), intent (inout) :: field2d(partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent (in)    :: mask(partit%myDim_nod2D+partit%eDim_nod2D)
      integer, intent (in)          :: n, h
      logical, intent (in)          :: do_stats
    end subroutine force_flux_consv
  end interface
end module force_flux_consv_interface
module compute_residual_interface
  interface
    subroutine compute_residual(field2d, mask, n, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_mesh),   intent(in),    target :: mesh
      type(t_partit), intent(inout), target :: partit
      real(kind=WP),  intent (in) :: field2d(partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP),  intent (in) :: mask(partit%myDim_nod2D+partit%eDim_nod2D)
      integer, intent (in)       :: n
    end subroutine compute_residual
  end interface
end module compute_residual_interface
module integrate_2D_interface
  interface
    subroutine integrate_2D(flux_global, flux_local, eff_vol, field2d, mask, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_mesh),   intent(in),    target :: mesh
      type(t_partit), intent(inout), target :: partit
      real(kind=WP), intent (out) :: flux_global(2), flux_local(2)
      real(kind=WP), intent (out) :: eff_vol(2)
      real(kind=WP), intent (in)  :: field2d(partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent (in)  :: mask(partit%myDim_nod2D   +partit%eDim_nod2D)
    end subroutine integrate_2D
  end interface
end module integrate_2D_interface

module update_atm_forcing_interface
    interface
        subroutine update_atm_forcing(istep, ice, tracers, dynamics, partit, mesh)
        USE MOD_TRACER
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        USE MOD_DYN
        integer,        intent(in)            :: istep
        type(t_ice),    intent(inout), target :: ice
        type(t_tracer), intent(in),    target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        type(t_dyn)   , intent(in),    target :: dynamics
        end subroutine update_atm_forcing
    end interface
end module update_atm_forcing_interface

module net_rec_from_atm_interface
  interface
    subroutine net_rec_from_atm(action, partit)
      USE MOD_PARTIT
      USE MOD_PARSUP
      logical,        intent(in)             :: action
      type(t_partit), intent(inout), target  :: partit
    end subroutine net_rec_from_atm
  end interface
end module net_rec_from_atm_interface
! Routines for updating ocean surface forcing fields
!-------------------------------------------------------------------------
subroutine update_atm_forcing(istep, ice, tracers, dynamics, partit, mesh)
  use o_PARAM
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  use MOD_TRACER
  use MOD_ICE
  use MOD_DYN
  use o_arrays
  use g_forcing_param
  use g_forcing_arrays
  use g_clock
  use g_config
  use g_comm_auto
  use g_rotate_grid
  use net_rec_from_atm_interface
  use g_sbf, only: sbc_do
  use g_sbf, only: atmdata, i_totfl, i_xwind, i_ywind, i_xstre, i_ystre, i_humi, i_qsr, i_qlw, i_tair, i_prec, i_mslp, i_cloud, i_snow, &
                                     l_xwind, l_ywind, l_xstre, l_ystre, l_humi, l_qsr, l_qlw, l_tair, l_prec, l_mslp, l_cloud, l_snow
#if defined (__oasis)
  use cpl_driver
#endif
  use gen_bulk
  use force_flux_consv_interface
#if defined (__recom)
  use REcoM_GloVar, only: x_co2atm, GloCO2flux_seaicemask
#endif

  implicit none
  integer,        intent(in)            :: istep
  type(t_ice)   , intent(inout), target :: ice
  type(t_tracer), intent(in),    target :: tracers
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  type(t_dyn)   , intent(in), target :: dynamics
  !_____________________________________________________________________________
  integer		   :: i, itime,n2,n,nz,k,elem
  real(kind=WP)            :: i_coef, aux
  real(kind=WP)	           :: dux, dvy,tx,ty,tvol
  real(kind=WP)            :: t1, t2
  !---wiso-code
  integer                  :: nt1, nt2
  real(kind=WP), parameter :: zwisomin = 1.e-6_WP
  !---wiso-code-end
#ifdef __oasis
  real(kind=WP)        				   :: flux_global(2), flux_local(2), eff_vol(2)
  real(kind=WP), dimension(:), allocatable , save  :: exchange
  real(kind=WP), dimension(:), allocatable , save  :: mask !, weight
  logical, save                                    :: firstcall=.true.
  logical                                          :: action
  logical                                          :: do_rotate_oce_wind=.false.
  logical                                          :: do_rotate_ice_wind=.false.
  INTEGER                                          :: my_global_rank, ierror
  INTEGER 					   :: status(MPI_STATUS_SIZE)
#endif
  !character(15)                         :: vari, filevari
  !character(4)                          :: fileyear
  !integer, parameter                    :: nci=192, ncj=94 ! T62 grid
  !real(kind=WP), dimension(nci,ncj)     :: array_nc, array_nc2,array_nc3,x
  !character(500)                        :: file
  !_____________________________________________________________________________
  ! pointer on necessary derived types
  real(kind=WP), dimension(:), pointer  :: u_ice, v_ice, u_w, v_w
  real(kind=WP), dimension(:), pointer  :: stress_atmice_x, stress_atmice_y
  real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
#if defined (__oasis) || defined (__ifsinterface)
  real(kind=WP), dimension(:), pointer  ::  oce_heat_flux, ice_heat_flux 
  real(kind=WP), dimension(:), pointer  ::  tmp_oce_heat_flux, tmp_ice_heat_flux 
#endif 
#if defined (__oifs) || defined (__ifsinterface)
  real(kind=WP), dimension(:), pointer  :: ice_temp, ice_alb, enthalpyoffuse
  real(kind=WP),               pointer  :: tmelt
  real(kind=WP), dimension(:,:,:), pointer :: UVnode
#endif
  real(kind=WP)              , pointer  :: rhoair
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  u_ice            => ice%uice(:)
  v_ice            => ice%vice(:)
  u_w              => ice%srfoce_u(:)
  v_w              => ice%srfoce_v(:)
  stress_atmice_x  => ice%stress_atmice_x(:)
  stress_atmice_y  => ice%stress_atmice_y(:)
  a_ice            => ice%data(1)%values(:)
  m_ice            => ice%data(2)%values(:)
  m_snow           => ice%data(3)%values(:)
#if defined (__oifs) || defined (__ifsinterface)
  ice_temp         => ice%data(4)%values(:)
  ice_alb          => ice%atmcoupl%ice_alb(:)
  enthalpyoffuse   => ice%atmcoupl%enthalpyoffuse(:)
  tmelt            => ice%thermo%tmelt
  UVnode           => dynamics%uvnode(:,:,:)
#endif      
#if defined (__oasis) || defined (__ifsinterface)
  oce_heat_flux    => ice%atmcoupl%oce_flx_h(:)
  ice_heat_flux    => ice%atmcoupl%ice_flx_h(:)
  tmp_oce_heat_flux=> ice%atmcoupl%tmpoce_flx_h(:)
  tmp_ice_heat_flux=> ice%atmcoupl%tmpice_flx_h(:)
#endif 
  rhoair           => ice%thermo%rhoair
  
  !_____________________________________________________________________________
  t1=MPI_Wtime()
#if defined (__oasis)
     if (firstcall) then
        allocate(exchange(myDim_nod2D+eDim_nod2D), mask(myDim_nod2D+eDim_nod2D))
        allocate(a2o_fcorr_stat(nrecv,6))
        a2o_fcorr_stat=0.
        exchange      =0.
        mask          =0.
        firstcall=.false.
     end if
     do i=1,nsend
         exchange  =0.
         if (i.eq.1) then
#if defined (__oifs)
            ! AWI-CM3 outgoing state vectors
              do n=1,myDim_nod2D+eDim_nod2D
              exchange(n)=tracers%data(1)%values(1, n)+tmelt	           ! sea surface temperature [K]
              end do
            elseif (i.eq.2) then
              exchange(:) = a_ice(:)                                  ! ice concentation [%]
            elseif (i.eq.3) then
              exchange(:) = m_snow(:)                                 ! snow thickness
            elseif (i.eq.4) then
              exchange(:) = ice_temp(:)                               ! ice surface temperature
            elseif (i.eq.5) then
              exchange(:) = ice_alb(:)                                ! ice albedo
            elseif (i.eq.6) then
              do n=1,myDim_nod2D+eDim_nod2D
                exchange(n) = UVnode(1,1,n)
              end do
            elseif (i.eq.7) then
              do n=1,myDim_nod2D+eDim_nod2D
                exchange(n) = UVnode(2,1,n)
              end do
#if defined (__recom)
            elseif (i.eq.8) then
              ! GloCO2flux_seaicemask is in [mmol/m2/s], need [kg/m2/s]
              ! Conversion: 1.0e-3_WP -> mol/s -> kg/s
              ! 1 mol CO2 = 44.0095 g/mol = 0.0440095 kg/mol (NIST 2018)
              exchange(:) = GloCO2flux_seaicemask(:) * 1.0e-3_WP * 0.0440095_WP  ! [kg m⁻² s⁻¹]
#endif
            else    
            print *, 'not installed yet or error in cpl_oasis3mct_send', mype
#else   ! oifs
            ! AWI-CM2 outgoing state vectors
            do n=1,myDim_nod2D+eDim_nod2D
            exchange(n)=tracers%data(1)%values(1, n)                     ! sea surface temperature [°C]
            end do
            elseif (i.eq.2) then
            exchange(:) = m_ice(:)                                  ! ice thickness [m]
            elseif (i.eq.3) then
            exchange(:) = a_ice(:)                                  ! ice concentation [%]
            elseif (i.eq.4) then
            exchange(:) = m_snow(:)                                 ! snow thickness
!---wiso-code
            elseif (i.eq.5) then
             exchange(:) = 0.0_WP
             if (lwiso) then
              nt1 = index_wiso_tracers(1)
              nt2 = index_wiso_tracers(3)
              where (tracers%data(nt2)%values(1,:) > zwisomin)
                exchange(:) = (tracers%data(nt1)%values(1,:)/tracers%data(nt2)%values(1,:)/wiso_smow(1) - 1._WP)*1000._WP ! delta 18O of surface water
              end where
             end if

            elseif (i.eq.6) then
             exchange(:) = 0.0_WP
             if (lwiso) then
              nt1 = index_wiso_tracers(2)
              nt2 = index_wiso_tracers(3)
              where (tracers%data(nt2)%values(1,:) > zwisomin)
                exchange(:) = (tracers%data(nt1)%values(1,:)/tracers%data(nt2)%values(1,:)/wiso_smow(2) - 1._WP)*1000._WP ! delta D of surface water
              end where
             end if

            elseif (i.eq.7) then
            exchange(:) = 0.0_WP                          ! delta H216O of surface water is set to zero permill

            elseif (i.eq.8) then
             exchange(:) = 0.0_WP
             if (lwiso) then
              where (tr_arr_ice(:, 3) > zwisomin)
                exchange(:) = (tr_arr_ice(:, 1)/tr_arr_ice(:, 3)/wiso_smow(1) - 1._WP)*1000._WP ! delta 18O of sea ice
              end where
             end if

            elseif (i.eq.9) then
             exchange(:) = 0.0_WP
             if (lwiso) then
              where (tr_arr_ice(:, 3) > zwisomin)
                exchange(:) = (tr_arr_ice(:, 2)/tr_arr_ice(:, 3)/wiso_smow(2) - 1._WP)*1000._WP ! delta 18O of sea ice
              end where
             end if

            elseif (i.eq.10) then
            exchange(:) = 0.0_WP                          ! delta H216O of sea ice is set to zero permill
!---wiso-code-end
            else	    
            print *, 'not installed yet or error in cpl_oasis3mct_send', mype
#endif  ! oifs
         endif
         call cpl_oasis3mct_send(i, exchange, action, partit)
      end do
#ifdef VERBOSE
      do i=1, nsend 
        if (mype==0) write(*,*) 'SEND: field ', i, ' max val:', maxval(exchange), ' . ACTION? ', action 
      end do
#endif
      mask=1.
      do i=1,nrecv
         exchange =0.0
         call cpl_oasis3mct_recv (i, exchange, action, partit)
         !if (.not. action) cycle
	 !Do not apply a correction at first time step!
    if (i==1 .and. action .and. istep/=1) call net_rec_from_atm(action, partit)
        if (i.eq.1) then
            if (.not. action) cycle
            stress_atmoce_x(:) =  exchange(:)                    ! taux_oce
            do_rotate_oce_wind=.true.
        elseif (i.eq.2) then
            if (.not. action) cycle
            stress_atmoce_y(:) =  exchange(:)                    ! tauy_oce
            do_rotate_oce_wind=.true.
        elseif (i.eq.3) then
            if (.not. action) cycle	
            stress_atmice_x(:) =  exchange(:)                    ! taux_ice
            do_rotate_ice_wind=.true.
        elseif (i.eq.4) then
            if (.not. action) cycle	
            stress_atmice_y(:) =  exchange(:)                    ! tauy_ice
            do_rotate_ice_wind=.true.	     
        elseif (i.eq.5) then
            if (action) then 
                prec_rain(:)    =  exchange(:)	                  ! tot_prec
                mask=1.
                call force_flux_consv(prec_rain, mask, i, 0,action, partit, mesh)
            end if
        elseif (i.eq.6) then 
            if (action) then
                prec_snow(:)    =  exchange(:)                    ! snowfall
                mask=1.
                call force_flux_consv(prec_snow, mask,i,1,action, partit, mesh) ! Northern hemisphere
                call force_flux_consv(prec_snow, mask,i,2,action, partit, mesh) ! Southern Hemisphere
            end if
         elseif (i.eq.7) then
            if (action) then
                evap_no_ifrac(:)     =  exchange(:)        	  ! tot_evap
                tmp_evap_no_ifrac(:) =  exchange(:) 		  ! to reset for flux 
                                                              ! correction
            end if
            mask=1.-a_ice
            evap_no_ifrac(:)     =  tmp_evap_no_ifrac(:)
            call force_flux_consv(evap_no_ifrac,mask,i,0,action, partit, mesh)
        elseif (i.eq.8) then
            if (action) then
                sublimation(:)       =  exchange(:)        	  ! tot_subl
                tmp_sublimation(:)   =  exchange(:) 		  ! to reset for flux 
                                                              ! correction
            end if
            mask=a_ice 
            sublimation(:)       =  tmp_sublimation(:)
            call force_flux_consv(sublimation,mask,i,1,action, partit, mesh) ! Northern hemisphere
            call force_flux_consv(sublimation,mask,i,2,action, partit, mesh) ! Southern Hemisphere
        elseif (i.eq.9) then
            if (action) then
                oce_heat_flux(:)     =  exchange(:)        	  ! heat_oce
                tmp_oce_heat_flux(:) =  exchange(:) 		  ! to reset for flux 
                                                              ! correction
            end if
            mask=1.-a_ice
            oce_heat_flux(:)     =  tmp_oce_heat_flux(:)
            call force_flux_consv(oce_heat_flux, mask, i, 0,action, partit, mesh)
        elseif (i.eq.10) then
            if (action) then
                ice_heat_flux(:)     =  exchange(:)        	  ! heat_ice
                tmp_ice_heat_flux(:) =  exchange(:) 		  ! to reset for flux 
                                                              ! correction
            end if
            mask=a_ice
            ice_heat_flux(:)     =  tmp_ice_heat_flux(:)
            call force_flux_consv(ice_heat_flux, mask, i, 1,action, partit, mesh) ! Northern hemisphere
            call force_flux_consv(ice_heat_flux, mask, i, 2,action, partit, mesh) ! Southern Hemisphere	     
        elseif (i.eq.11) then
            if (action) then
                shortwave(:)         =  exchange(:)		        ! heat_swr
                tmp_shortwave(:)     =  exchange(:) 		    ! to reset for flux 
                                                                ! correction
            end if
            mask=1.-a_ice
            shortwave(:)   =  tmp_shortwave(:)
            call force_flux_consv(shortwave, mask, i, 0,action, partit, mesh)
        elseif (i.eq.12) then
            if (action) then
                runoff(:)            =  exchange(:)        ! AWI-CM2: runoff, AWI-CM3: runoff + excess snow on glaciers
                mask=1.
                call force_flux_consv(runoff, mask, i, 0,action, partit, mesh)
            end if
#if defined (__oifs)

         elseif (i.eq.13) then
             if (action) then
                enthalpyoffuse(:)            =  exchange(:)*333.55*1000000.0        ! enthalpy of fusion via solid water discharge from glaciers
                enthalpyoffuse = -min(enthalpyoffuse, 1000.0)
                runoff(:)            = runoff(:) + exchange(:)                      ! Add calving massflux to the liquid runoff. Heatflux goes into enthalpyoffuse.

                mask=1.
                call force_flux_consv(enthalpyoffuse, mask, i, 0, action, partit, mesh)
             end if
         elseif (i.eq.14) then
             if (action) then
                u_wind(:)                     = exchange(:)        ! zonal wind
             end if
         elseif (i.eq.15) then
             if (action) then
                v_wind(:)                     = exchange(:)        ! meridional wind
             end if
#if defined (__recom)
         elseif (i.eq.16) then
             if (action) then
                ! Convert mass mixing ratio (kg/kg) to mole fraction (dimensionless)
                ! MW_CO2 = 44.0095 g/mol (NIST 2018)
                ! MW_dry_air = 28.9647 g/mol (standard atmosphere composition)
                x_co2atm(:) = exchange(:) * ((28.9647_WP/44.0095_WP)*1e6_WP)  ! [mole fraction]
             end if
#endif
#else ! oifs
         elseif (i.eq.13) then
            if (action) then
                 if (lwiso) then         
                     www3(:)         =  exchange(:)
                 else if (use_icebergs) then    
                     u_wind(:)                     = exchange(:)        ! zonal wind
                 end if
             end if
             mask=1.
             if (lwiso) then         
                 call force_flux_consv(www3, mask, i, 0,action, partit, mesh)
             else if (use_icebergs) then    
                 call force_flux_consv(u_wind, mask, i, 0, action, partit, mesh)
             end if
         elseif (i.eq.14) then
             if (action) then
                 if (lwiso) then         
                     www1(:)         =  exchange(:)
                 else if (use_icebergs) then    
                     v_wind(:)                     = exchange(:)        ! meridional wind
                 end if
             end if
             mask=1.
             if (lwiso) then         
                 call force_flux_consv(www1, mask, i, 0,action, partit, mesh)
             else if (use_icebergs) then    
                 call force_flux_consv(v_wind, mask, i, 0, action, partit, mesh)
             end if
         elseif (i.eq.15) then
             if (action) then
             ! tot_prec_hdo over water: this variable includes (i) rain over open water and sea ice, (ii) snow and evap over open water,  (iii) river runoff
             www2(:)         =  exchange(:)
             end if
             mask=1.
             if (lwiso) then
             call force_flux_consv(www2, mask, i, 0,action, partit, mesh)
             end if
         elseif (i.eq.17) then
             if (action) then
             ! snowfall_o18 over seaice: this variable includes snow and sublimation over sea ice areas
             iii1(:)         =  exchange(:)
             tmp_iii1(:)     =  exchange(:)                   ! to reset for flux correction
             end if
             mask=a_ice
             iii1(:)         =  tmp_iii1(:)
             if (lwiso) then
                 call force_flux_consv(iii1,mask,i,1,action, partit, mesh) ! Northern hemisphere
                 call force_flux_consv(iii1,mask,i,2,action, partit, mesh) ! Southern Hemisphere
             end if
         elseif (i.eq.18) then
             if (action) then
             ! snowfall_hdo over seaice: this variable includes snow and sublimation over sea ice areas
             iii2(:)         =  exchange(:)
             tmp_iii2(:)     =  exchange(:)                   ! to reset for flux correction
             end if
             mask=a_ice
             iii2(:)         =  tmp_iii2(:)
             if (lwiso) then
                 call force_flux_consv(iii2,mask,i,1,action, partit, mesh) ! Northern hemisphere
                 call force_flux_consv(iii2,mask,i,2,action, partit, mesh) ! Southern Hemisphere
             end if
         elseif (i.eq.16) then
             if (action) then
             ! snowfall_o16 over seaice: this variable includes snow and sublimation over sea ice areas
             iii3(:)         =  exchange(:)
             tmp_iii3(:)     =  exchange(:)                   ! to reset for flux correction
             end if
             mask=a_ice
             iii3(:)         =  tmp_iii3(:)
             if (lwiso) then
                 call force_flux_consv(iii3,mask,i,1,action, partit, mesh) ! Northern hemisphere
                 call force_flux_consv(iii3,mask,i,2,action, partit, mesh) ! Southern Hemisphere
             end if
         elseif (i.eq.19) then
             if (action) then
                u_wind(:)                     = exchange(:)        ! meridional wind
             end if
             mask=1
             if (use_icebergs.and.lwiso) then    
                 call force_flux_consv(u_wind, mask, i, 0, action, partit, mesh)
             end if
         elseif (i.eq.20) then
             if (action) then
                v_wind(:)                     = exchange(:)        ! meridional wind
             end if
             mask=1
             if (use_icebergs.and.lwiso) then    
                 call force_flux_consv(v_wind, mask, i, 0, action, partit, mesh)
             end if
#endif !   oifs
         end if

#ifdef VERBOSE
	  if (mype==0) then
		write(*,*) 'FESOM RECV: flux ', i, ', max val: ', maxval(exchange)
	  end if
#endif
      end do

    if ((do_rotate_oce_wind .AND. do_rotate_ice_wind) .AND. rotated_grid) then
        do n=1, myDim_nod2D+eDim_nod2D
            call vector_g2r(stress_atmoce_x(n), stress_atmoce_y(n), coord_nod2D(1, n), coord_nod2D(2, n), 0)
            call vector_g2r(stress_atmice_x(n), stress_atmice_y(n), coord_nod2D(1, n), coord_nod2D(2, n), 0)
        end do
        do_rotate_oce_wind=.false.
        do_rotate_ice_wind=.false.
    end if
#else
#ifndef __ifsinterface
  call sbc_do(partit, mesh)
!$OMP PARALLEL DO
  DO n=1, myDim_nod2D+eDim_nod2D
     u_wind(n)    = atmdata(i_xwind,n)
     v_wind(n)    = atmdata(i_ywind,n)
     if (l_xstre) then
        stress_atmoce_x(n)=atmdata(i_xstre, n)
     end if
     if (l_ystre) then
        stress_atmoce_y(n)=atmdata(i_ystre, n)
     end if
     shum(n)      = atmdata(i_humi ,n)
     longwave(n)  = atmdata(i_qlw  ,n)
     shortwave(n) = atmdata(i_qsr  ,n)
     Tair(n)      = atmdata(i_tair ,n)-273.15_WP
     prec_rain(n) = atmdata(i_prec ,n)/1000._WP
     prec_snow(n) = atmdata(i_snow ,n)/1000._WP
     if (l_mslp) then
        press_air(n) = atmdata(i_mslp ,n) ! unit should be Pa
     end if
  END DO
!$OMP END PARALLEL DO

  if (use_cavity) then 
!$OMP PARALLEL DO
    do i=1,myDim_nod2d+eDim_nod2d
        if (ulevels_nod2d(i)>1) then
            u_wind(i)   = 0.0_WP
            v_wind(i)   = 0.0_WP
            shum(i)     = 0.0_WP
            longwave(i) = 0.0_WP
            Tair(i)     = 0.0_WP
            prec_rain(i)= 0.0_WP
            prec_snow(i)= 0.0_WP
            if (l_mslp) then
               press_air(i)= 0.0_WP 
            end if
            runoff(i)   = 0.0_WP
        end if
    end do
!$OMP END PARALLEL DO
  endif 
  ! second, compute exchange coefficients
  ! 1) drag coefficient 
  if(AOMIP_drag_coeff) then
     call cal_wind_drag_coeff(partit)
  end if
  ! 2) drag coeff. and heat exchange coeff. over ocean in case using ncar formulae
  if(ncar_bulk_formulae) then
!     cd_atm_oce_arr=0.0_WP
!     ch_atm_oce_arr=0.0_WP
!     ce_atm_oce_arr=0.0_WP
     call ncar_ocean_fluxes_mode(ice, partit, mesh)
  elseif(AOMIP_drag_coeff) then
     cd_atm_oce_arr=cd_atm_ice_arr
  end if
  ! third, compute wind stress
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dux, dvy, aux)

  do i=1,myDim_nod2d+eDim_nod2d   
     !__________________________________________________________________________
     if (ulevels_nod2d(i)>1) then
        stress_atmoce_x(i)=0.0_WP
        stress_atmoce_y(i)=0.0_WP
        stress_atmice_x(i)=0.0_WP
        stress_atmice_y(i)=0.0_WP
        cycle
     end if 
     
     !__________________________________________________________________________
     dux=u_wind(i)-(1.0_WP-Swind)*u_w(i) 
     dvy=v_wind(i)-(1.0_WP-Swind)*v_w(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     if ( .not. l_xstre) then
        stress_atmoce_x(i) = Cd_atm_oce_arr(i)*aux*dux
     end if
     if ( .not. l_ystre) then
        stress_atmoce_y(i) = Cd_atm_oce_arr(i)*aux*dvy
     end if
     !__________________________________________________________________________
     dux=u_wind(i)-u_ice(i) 
     dvy=v_wind(i)-v_ice(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     stress_atmice_x(i) = Cd_atm_ice_arr(i)*aux*dux
     stress_atmice_y(i) = Cd_atm_ice_arr(i)*aux*dvy
  end do
!$OMP END PARALLEL DO
  ! heat and fresh water fluxes are treated in i_therm and ice2ocean
#endif /* skip all in case of __ifsinterface */
#endif /* (__oasis) */

  t2=MPI_Wtime()

#ifdef VERBOSE
  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
     write(*,*) 'update forcing data took', t2-t1
  end if
#endif

end subroutine update_atm_forcing
!
!------------------------------------------------------------------------------------
!
#if defined (__oasis)
!
!=================================================================
!
! Modifies the flux 'field2d' in order to conserve 
! the net fluxes.
! We distinguish between NH and SH fluxes (integer 'h') 
! where possible.
!
!================================================================= 
! History :
!  07-11  (D.Sidorenko,	AWI Germany) first routine
!  10-12  (T.Rackow, 	AWI Germany) code reordering and cleanup  
!-----------------------------------------------------------------
!
SUBROUTINE force_flux_consv(field2d, mask, n, h, do_stats, partit, mesh)

  use g_forcing_arrays,	only : 	atm_net_fluxes_north, atm_net_fluxes_south, 	&
  				oce_net_fluxes_north, oce_net_fluxes_south, 	&
				flux_correction_north, flux_correction_south,	&
				flux_correction_total
  use mod_mesh
  USE MOD_PARTIT
  USE MOD_PARSUP
  use cpl_driver,	 only : nrecv, cpl_recv, a2o_fcorr_stat
  use o_PARAM,           only : mstep, WP
  use compute_residual_interface
  use integrate_2D_interface
  IMPLICIT NONE
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit  
  real(kind=WP),  INTENT (INOUT) 	:: field2d(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP),  INTENT (IN)	        :: mask(partit%myDim_nod2D+partit%eDim_nod2D)
  INTEGER,        INTENT (IN)	        :: n
  INTEGER,        INTENT (IN)           :: h !hemisphere: 0=GL, 1=NH, 2=SH
  logical,        INTENT (IN)	        :: do_stats
  real(kind=WP)			        :: rmask(partit%myDim_nod2D+partit%eDim_nod2D)   
  real(kind=WP)		                :: weight(partit%myDim_nod2D+partit%eDim_nod2D)    
  real(kind=WP)                         :: flux_global(2), flux_local(2)
  real(kind=WP)                         :: eff_vol(2)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

#if defined (__oifs)
  return !OIFS-FESOM2 coupling uses OASIS3MCT conservative remapping instead
#endif

  if (mstep==1) then
	if (mype == 0) write(*,*) 'Do not apply a correction at first time step for ', trim(cpl_recv(n))
	return
  end if

  
  !just keep NH or SH part in mask
  rmask=mask
  SELECT CASE(h)
	   CASE( 1 ) ; where (geo_coord_nod2D(2, :)< 0) !just NH 
	                     rmask=0. 
		       end where
	   CASE( 2 ) ; where (geo_coord_nod2D(2, :)>=0) !just SH 
	                     rmask=0. 
		       end where
	   CASE DEFAULT 
	   !keep global mask
  END SELECT
 
  !residual (net) fluxes; computes also oce_net_fluxes_*
  call compute_residual(field2d, rmask, n, partit, mesh)
  
#ifdef VERBOSE
  if (mype == 0) then
  !atm net fluxes, oce net fluxes before modification
  write(*,'(3A,3e15.7)') 'atm NH SH GL / ', trim(cpl_recv(n)), ': ', 		&
  			    atm_net_fluxes_north(n), atm_net_fluxes_south(n), 	&
			    atm_net_fluxes_north(n)+ atm_net_fluxes_south(n)
			    
  write(*,'(3A,3e15.7)') 'oce NH SH GL / ', trim(cpl_recv(n)), ': ', 		&
  			    oce_net_fluxes_north(n), oce_net_fluxes_south(n), 	&
			    oce_net_fluxes_north(n)+ oce_net_fluxes_south(n)
  end if			    
#endif
  if (do_stats)	then
     if (h==0 .or. h==1) then
        a2o_fcorr_stat(n, 1) =a2o_fcorr_stat(n, 1)+atm_net_fluxes_north(n)
        a2o_fcorr_stat(n, 4) =a2o_fcorr_stat(n, 4)+oce_net_fluxes_north(n)
     end if
     if (h==0 .or. h==2) then
        a2o_fcorr_stat(n, 2) =a2o_fcorr_stat(n, 2)+atm_net_fluxes_south(n)
        a2o_fcorr_stat(n, 5) =a2o_fcorr_stat(n, 5)+oce_net_fluxes_south(n)
     end if     
  end if

  !integrate (masked) abs(field2d) to get positive weights
  call integrate_2D(flux_global, flux_local, eff_vol, abs(field2d), rmask, partit, mesh)
  
  !get weight pattern with integral 1
  if (abs(sum(flux_global))>1.e-10) then
    weight=abs(field2d/sum(flux_global))
  else
    !should rarely happen
    weight=1.0_WP / sum(eff_vol)
    write(*,*) 'Warning: Constant redistribution for flux ', trim(cpl_recv(n))
  end if
  
  !weight is still global 2D field, just keep NH or SH part
  where (rmask<1.e-10)
    weight=0.
  end where
  
  !redistribute the residual according to the mask
  SELECT CASE(h)
	   CASE( 0 ) ; field2d=field2d+weight*flux_correction_total(n) ! GL	       
	   CASE( 1 ) ; field2d=field2d+weight*flux_correction_north(n) ! NH
	   CASE( 2 ) ; field2d=field2d+weight*flux_correction_south(n) ! SH
  END SELECT
  
  !check conservation
  call integrate_2D(flux_global, flux_local, eff_vol, field2d, rmask, partit, mesh)
#ifdef VERBOSE
  if (mype == 0) then
  write(*,'(3A,3e15.7)') 'oce NH SH GL / ', trim(cpl_recv(n)), ': ', 		&
  			   flux_global(1), flux_global(2), sum(flux_global)
  end if
#endif
  
  !last flux			   
  if (n==nrecv .AND. mype==0) write(*,*) 'Fluxes have been modified.'  
END SUBROUTINE force_flux_consv

!
! Compute the difference between the net fluxes seen by the atmosphere
! and ocean component (residual flux) for flux n.
!
SUBROUTINE compute_residual(field2d, mask, n, partit, mesh)

  use g_forcing_arrays,	only : 	atm_net_fluxes_north, atm_net_fluxes_south, 	&
  				oce_net_fluxes_north, oce_net_fluxes_south, 	&
				flux_correction_north, flux_correction_south,	&
				flux_correction_total
  use o_PARAM, only : WP 
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  use integrate_2D_interface
 
  IMPLICIT NONE
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit  
  real(kind=WP), INTENT(IN)   :: field2d(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP), INTENT(IN)   :: mask(partit%myDim_nod2D+partit%eDim_nod2D)  
  INTEGER,       INTENT(IN)   :: n
  
  real(kind=WP)               :: flux_global(2), flux_local(2)
  real(kind=WP)               :: eff_vol(2)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  !compute net flux (for flux n) on ocean side
  call integrate_2D(flux_global, flux_local, eff_vol, field2d, mask, partit, mesh)
  oce_net_fluxes_north(n)=flux_global(1)
  oce_net_fluxes_south(n)=flux_global(2)
  
  !compute the residual fluxes for NH, SH and Global
  flux_correction_north(n)= atm_net_fluxes_north(n) - oce_net_fluxes_north(n)
  flux_correction_south(n)= atm_net_fluxes_south(n) - oce_net_fluxes_south(n)
  flux_correction_total(n)= flux_correction_north(n) + flux_correction_south(n)
END SUBROUTINE compute_residual

!
! -field_2d (input) is any (partitioned) 2D field
! -flux_local  (returned) is the net local flux (for current pc)
! -flux_global (returned) is the communicated and summarized flux_local  
!
SUBROUTINE integrate_2D(flux_global, flux_local, eff_vol, field2d, mask, partit, mesh)
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  use o_PARAM, only: WP 
  IMPLICIT NONE
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(in),    target :: partit
  real(kind=WP), INTENT(OUT)  :: flux_global(2), flux_local(2)
  real(kind=WP), INTENT(OUT)  :: eff_vol(2)
  real(kind=WP), INTENT(IN)   :: field2d(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP), INTENT(IN)   :: mask(partit%myDim_nod2D   +partit%eDim_nod2D) 
   
  real(kind=WP)               :: eff_vol_local(2)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  flux_local(1)=sum(lump2d_north*field2d(1:myDim_nod2D)*mask(1:myDim_nod2D))
  flux_local(2)=sum(lump2d_south*field2d(1:myDim_nod2D)*mask(1:myDim_nod2D))
  call MPI_AllREDUCE(flux_local, flux_global, 2, &
  		     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
		     
		     
  eff_vol_local(1)=sum(lump2d_north*mask(1:myDim_nod2D))
  eff_vol_local(2)=sum(lump2d_south*mask(1:myDim_nod2D))
  call MPI_AllREDUCE(eff_vol_local, eff_vol,  2, & 
  		     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
		     
END SUBROUTINE integrate_2D

!
!---------------------------------------------------------------------------------------------------
!  Receive atmospheric net fluxes (atm_net_fluxes_north and atm_net_fluxes_south)
!
!SUBROUTINE net_rec_from_atm(action)
!  
!
!  use g_forcing_arrays,	only : 	atm_net_fluxes_north, atm_net_fluxes_south, 	&
!  				oce_net_fluxes_north, oce_net_fluxes_south, 	&
!				flux_correction_north, flux_correction_south,	&
!				flux_correction_total
!  use cpl_driver,	only : 	nrecv
!  
!  IMPLICIT NONE
!
!  LOGICAL, INTENT(IN)   :: action
!
!  if (action) then
!  atm_net_fluxes_north=0.
!  atm_net_fluxes_south=0.
!
!  atm_net_fluxes_north=(/ 1.,2.,3.,4. , 0.5e7, -0.1e10, -0.3e8, 0.5e-13, -0.7e17, -0.4e-6, 0.5e5, 0.0/)
!  atm_net_fluxes_south=(/ 1.,2.,3.,4. , 0.4e7, 0.3e06, -0.3e8, -0.3,-0.7e17, -0.4e-6, 0.2e6, 0.0/)
!  end if
!END SUBROUTINE net_rec_from_atm
!
!
!---------------------------------------------------------------------------------------------------
!  Receieve atmospheric net fluxes (atm_net_fluxes_north and atm_net_fluxes_south)
!
SUBROUTINE net_rec_from_atm(action, partit)
!
  use g_forcing_arrays
  use cpl_driver
  use o_PARAM, only: WP
  USE MOD_PARTIT
  USE MOD_PARSUP
  IMPLICIT NONE

  LOGICAL,        INTENT (IN)   		  :: action
  type(t_partit), intent(inout), target           :: partit
  INTEGER                                         :: my_global_rank, ierror
  INTEGER                                         :: n  
  INTEGER 					  :: status(MPI_STATUS_SIZE,partit%npes) 
  INTEGER                                         :: request(2)
  real(kind=WP)                 		  :: aux(nrecv)
#if defined (__oifs)
  return  !OIFS-FESOM2 coupling uses OASIS3MCT conservative remapping and recieves no net fluxes here.
#endif

  if (action) then
     CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_global_rank, ierror)
     atm_net_fluxes_north=0.
     atm_net_fluxes_south=0.
     if (my_global_rank==target_root) then
	CALL MPI_IRecv(atm_net_fluxes_north(1), nrecv, MPI_DOUBLE_PRECISION, source_root, 111, MPI_COMM_WORLD, request(1), partit%MPIerr)
        CALL MPI_IRecv(atm_net_fluxes_south(1), nrecv, MPI_DOUBLE_PRECISION, source_root, 112, MPI_COMM_WORLD, request(2), partit%MPIerr)
        CALL MPI_Waitall(2, request, status, partit%MPIerr)
     end if
  call MPI_Barrier(partit%MPI_COMM_FESOM, partit%MPIerr)     
  call MPI_AllREDUCE(atm_net_fluxes_north(1), aux, nrecv, MPI_DOUBLE_PRECISION, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  atm_net_fluxes_north=aux
  call MPI_AllREDUCE(atm_net_fluxes_south(1), aux, nrecv, MPI_DOUBLE_PRECISION, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  atm_net_fluxes_south=aux
  end if
END SUBROUTINE net_rec_from_atm
#endif
