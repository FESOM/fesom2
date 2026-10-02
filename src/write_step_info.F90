module write_step_info_module
    USE g_config, only: dt, use_ice, use_icebergs, ib_num, logfile_outfreq, which_ALE, &
            toy_ocean
    USE MOD_MESH
    USE MOD_PARTIT
    USE par_support_module, only: par_ex
    USE MOD_TRACER
    USE MOD_DYN
    USE MOD_ICE
    USE o_PARAM
    USE o_ARRAYS, only: water_flux, heat_flux, pgf_x, pgf_y, Av, Kv, density_dmoc, stress_surf
    USE diagnostics
    USE g_comm_auto
    USE g_support
    USE iceberg_params
    USE io_BLOWUP
    USE g_forcing_arrays
    USE fesom_monitor_module, only: t_monitor, monitor_fill_ocean, monitor_value, monitor_integral, monitor_nonfinite
    USE iceberg_element

    implicit none

    private
    public :: write_step_info, check_blowup, write_enegry_info, &
              plot_fesomlogo, plot_fesomlogo_lildevil, &
              plot_fesomlogo_expl

contains

!
!
!===============================================================================
subroutine write_step_info(istep, outfreq, mon, dynamics, partit)
  implicit none
  integer        , intent(in)            :: istep, outfreq
  type(t_monitor), intent(in)            :: mon
  type(t_dyn)    , intent(in)   , target :: dynamics
  type(t_partit) , intent(inout), target :: partit
  real(kind=WP_full)                     :: int_eta, int_hbar, int_deta, int_dhbar, int_wflux, &
                                            int_hflux, int_temp, int_salt

  if (mod(istep,outfreq)/=0 .or. partit%mype/=0) return

  int_eta   = monitor_integral(mon, 'eta')
  int_hbar  = monitor_integral(mon, 'hbar')
  int_deta  = monitor_integral(mon, 'deta')
  int_dhbar = monitor_integral(mon, 'dhbar')
  int_wflux = monitor_integral(mon, 'wflux')
  int_hflux = monitor_integral(mon, 'hflux')
  int_temp  = monitor_integral(mon, 'temp')
  int_salt  = monitor_integral(mon, 'salt')

  write(*,*) '___CHECK GLOBAL OCEAN VARIABLES --> mstep=',mstep
  write(*,*) '  ___global estimat of eta & hbar____________________'
  write(*,*) '   int(eta), int(hbar)      =', int_eta, int_hbar
  write(*,*) '    --> error(eta-hbar)     =', int_eta-int_hbar
  write(*,*) '   min(eta) , max(eta)      =', mv('eta','min'), mv('eta','max')
  write(*,*) '   max(hbar), max(hbar)     =', mv('hbar','min'), mv('hbar','max')
  write(*,*)
  write(*,*) '   int(deta), int(dhbar)    =', int_deta, int_dhbar
  write(*,*) '    --> error(deta-dhbar)   =', int_deta-int_dhbar
  write(*,*) '    --> error(deta-wflux)   =', int_deta-int_wflux
  write(*,*) '    --> error(dhbar-wflux)  =', int_dhbar-int_wflux
  write(*,*)
  write(*,*) '   -int(wflux)*dt          =', int_wflux*dt*(-1.0)
  write(*,*) '   int(deta )-int(wflux)*dt =', int_deta-int_wflux*dt*(-1.0)
  write(*,*) '   int(dhbar)-int(wflux)*dt =', int_dhbar-int_wflux*dt*(-1.0)
  write(*,*)
  write(*,*) '  ___global min/max/mean  --> mstep=',mstep,'____________'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") '       eta= ', mv('eta','min')   , ' | ', mv('eta','max')   , ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") '      deta= ', mv('deta','min')  , ' | ', mv('deta','max')  , ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") '      hbar= ', mv('hbar','min')  , ' | ', mv('hbar','max')  , ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, ES10.3)") '     wflux= ', mv('wflux','min') , ' | ', mv('wflux','max') , ' | ', int_wflux
  write(*,"(A15, ES10.3, A3, ES10.3, A3, ES10.3)") '     hflux= ', mv('hflux','min') , ' | ', mv('hflux','max') , ' | ', int_hflux
  write(*,"(A15, ES10.3, A3, ES10.3, A3, ES10.3)") '      temp= ', mv('temp','min')  , ' | ', mv('temp','max')  , ' | ', int_temp
  write(*,"(A15, ES10.3, A3, ES10.3, A3, ES10.3)") '      salt= ', mv('salt','min')  , ' | ', mv('salt','max')  , ' | ', int_salt
  if (ldiag_dMOC) then
      write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") '      dens= ', mv('dens','min')-1000.0_WP, ' | ', mv('dens','max')-1000.0_WP, ' | ', 'N.A.'
  end if
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") ' wvel(1,:)= ', mv('wvel1','min') , ' | ', mv('wvel1','max') , ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") ' wvel(2,:)= ', mv('wvel2','min') , ' | ', mv('wvel2','max') , ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") ' uvel(1,:)= ', mv('uvel1','min') , ' | ', mv('uvel1','max') , ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") ' uvel(2,:)= ', mv('uvel2','min') , ' | ', mv('uvel2','max') , ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") ' vvel(1,:)= ', mv('vvel1','min') , ' | ', mv('vvel1','max') , ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") ' vvel(2,:)= ', mv('vvel2','min') , ' | ', mv('vvel2','max') , ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") 'hnode(1,:)= ', mv('hnode1','min'), ' | ', mv('hnode1','max'), ' | ', 'N.A.'
  write(*,"(A15, ES10.3, A3, ES10.3, A3, A10   )") 'hnode(2,:)= ', mv('hnode2','min'), ' | ', mv('hnode2','max'), ' | ', 'N.A.'
  write(*,"(A15, A10   , A3, ES10.3, A3, A10   )") '     cfl_z= ', ' N.A.'   , ' | ', mv('cfl_z','max') , ' | ', 'N.A.'
  write(*,"(A15, A10   , A3, ES10.3, A3, A10   )") '     pgf_x= ', ' N.A.'   , ' | ', mv('pgf_x','max') , ' | ', 'N.A.'
  write(*,"(A15, A10   , A3, ES10.3, A3, A10   )") '     pgf_y= ', ' N.A.'   , ' | ', mv('pgf_y','max') , ' | ', 'N.A.'
  write(*,"(A15, A10   , A3, ES10.3, A3, A10   )") '        Av= ', ' N.A.'   , ' | ', mv('Av','max')    , ' | ', 'N.A.'
  write(*,"(A15, A10   , A3, ES10.3, A3, A10   )") '        Kv= ', ' N.A.'   , ' | ', mv('Kv','max')    , ' | ', 'N.A.'
  if (use_ice)  then
  write(*,"(A15, A10   , A3, ES10.3, A3, A10)")    '     m_ice= ', ' N.A.'   , ' | ', mv('m_ice','max') , ' | ', 'N.A.'
  end if
  !________________________________________________________________________
  ! SSH CG solver. A named quantity in every standard block, like cfl_z
  ! above -- not something that only shows up when it misbehaves. The
  ! cumulative non-convergence count is here so a run that is quietly
  ! stalling is visible in a log people already read.
  ! Skipped on the split-explicit barotropic path, which has no solver.
  if (.not. dynamics%use_ssh_se_subcycl) then
  write(*,*)
  write(*,"(A, I6, A, I6, A, F8.2)")           '   ssh_cg iters= ', dynamics%solverinfo%iters_last, &
       '  | max= ', dynamics%solverinfo%iters_max,                                                  &
       '  | mean= ', real(dynamics%solverinfo%iters_sum)/real(max(dynamics%solverinfo%nsolves,1))
  write(*,"(A, ES10.3, A, ES10.3, A, I6)")     '   ssh_cg rms(r)= ', dynamics%solverinfo%resid_last, &
       '  | rtol= ', dynamics%solverinfo%rtol_last,                                                 &
       '  | nonconv= ', dynamics%solverinfo%nonconv
  end if

contains

  real(kind=WP) function mv(name, what)
    character(len=*), intent(in) :: name, what
    mv = monitor_value(mon, name, what)
  end function mv

end subroutine write_step_info
!
!
!===============================================================================
subroutine check_blowup(istep, mon, ice, dynamics, tracers, partit, mesh)
    implicit none
  
    type(t_monitor), intent(inout)       :: mon
    type(t_ice)   , intent(in)   , target :: ice
    type(t_dyn)   , intent(in)   , target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(in)   , target :: tracers
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                       :: n, nz, istep, ib
    logical                       :: tripped
    integer                       :: el, elidx
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV
    real(kind=WP), dimension(:,:)  , pointer :: Wvel, CFL_z
    real(kind=WP), dimension(:)    , pointer :: ssh_rhs, ssh_rhs_old
    real(kind=WP), dimension(:)    , pointer :: eta_n, d_eta
    real(kind=WP), dimension(:)    , pointer :: u_ice, v_ice
    real(kind=WP), dimension(:)    , pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:)    , pointer :: a_ice_old, m_ice_old, m_snow_old
    real(kind=WP), dimension(:), allocatable, target :: dhbar
    ! Saved, and the global-to-local map is fixed for the run, so it is built
    ! on first use. Each of the five dump sites below used to allocate it
    ! afresh and nothing ever deallocated it, so the second node to blow up
    ! aborted here and took the blowup diagnostic with it.
    integer, dimension(:), save, allocatable :: local_idx_of
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    UV          => dynamics%uv(:,:,:)
    Wvel        => dynamics%w(:,:)
    CFL_z       => dynamics%cfl_z(:,:)
    
    eta_n       => dynamics%eta_n(:)
    u_ice       => ice%uice(:)
    v_ice       => ice%vice(:)
    a_ice       => ice%data(1)%values(:)
    m_ice       => ice%data(2)%values(:)
    m_snow      => ice%data(3)%values(:)
    a_ice_old   => ice%data(1)%values_old(:)
    m_ice_old   => ice%data(2)%values_old(:)
    m_snow_old  => ice%data(3)%values_old(:)
    if ( .not. dynamics%use_ssh_se_subcycl) then 
        d_eta       => dynamics%d_eta(:)
        ssh_rhs     => dynamics%ssh_rhs(:)
        ssh_rhs_old => dynamics%ssh_rhs_old(:)
    else
        allocate(dhbar(myDim_nod2D+eDim_nod2D))
        dhbar = hbar-hbar_old
        d_eta => dhbar
    end if 
    
    !___________________________________________________________________________
    ! Decide from the ocean record of this step (global, so every rank agrees):
    ! eta and the surface layer, temperature and salinity on the wet levels.
    tripped = monitor_nonfinite(mon, 'eta')  > 0     .or. &
              monitor_value(mon, 'eta', 'min') < -10.0 .or. monitor_value(mon, 'eta', 'max') > 10.0 .or. &
              monitor_nonfinite(mon, 'deta') > 0
    if ( .not. trim(which_ALE)=='linfs') then
        tripped = tripped .or. monitor_nonfinite(mon, 'wvel1') > 0 .or. &
                  monitor_nonfinite(mon, 'hnode1') > 0 .or. monitor_value(mon, 'hnode1', 'min') < 0
    end if
    tripped = tripped .or. monitor_nonfinite(mon, 'temp') > 0 .or. &
              monitor_value(mon, 'temp', 'min') < -5.0 .or. monitor_value(mon, 'temp', 'max') > 60 .or. &
              monitor_nonfinite(mon, 'salt') > 0 .or. &
              monitor_value(mon, 'salt', 'min') < 3.0_WP - S_ref_anomaly .or. &
              monitor_value(mon, 'salt', 'max') > 45.0_WP - S_ref_anomaly
    if (.not. tripped) return

    !___________________________________________________________________________
    ! Report every offending node.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz)
    do n=1, myDim_nod2d       
       !___________________________________________________________________
       ! check ssh
       if ( ((eta_n(n) /= eta_n(n)) .or. eta_n(n)<-10.0 .or. eta_n(n)>10.0 .or. (d_eta(n) /= d_eta(n)) ) ) then
!$OMP CRITICAL
          write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
          write(*,*) ' --STOP--> found eta_n become NaN or <-10.0, >10.0'
          write(*,*) 'mype     = ',mype
          write(*,*) 'mstep    = ',istep
          write(*,*) 'node     = ',n
          write(*,*) 'uln, nln    = ',ulevels_nod2D(n), nlevels_nod2D(n)
          write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
          write(*,*)
          write(*,*) 'eta_n(n)    = ',eta_n(n)
          write(*,*) 'd_eta(n)    = ',d_eta(n)
          write(*,*)
          write(*,*) 'zbar_3d_n   = ',zbar_3d_n(:,n)
          write(*,*) 'Z_3d_n      = ',Z_3d_n(:,n)
          write(*,*)
          if ( .not. dynamics%use_ssh_se_subcycl) then 
            write(*,*) 'ssh_rhs = ',ssh_rhs(n),', ssh_rhs_old = ',ssh_rhs_old(n)
          end if
          write(*,*)
          write(*,*) 'hbar = ',hbar(n),', hbar_old = ',hbar_old(n)
          write(*,*)
          write(*,*) 'wflux = ',water_flux(n)
          write(*,*)
          if (.not. toy_ocean) then
            write(*,*) 'u_wind = ',u_wind(n),', v_wind = ',v_wind(n)
            write(*,*)
            do nz=1,nod_in_elem2D_num(n)
                    write(*,*) 'stress_surf(1:2,',nz,') = ',stress_surf(:,nod_in_elem2D(nz,n))
            end do
          end if
          if (use_ice) then
          write(*,*)
          write(*,*) 'm_ice = ',m_ice(n),', m_ice_old = ',m_ice_old(n)
          write(*,*) 'a_ice = ',a_ice(n),', a_ice_old = ',a_ice_old(n)
          end if 
          write(*,*)
          write(*,*) 'Wvel(:, n)  = ',Wvel(ulevels_nod2D(n):nlevels_nod2D(n),n)
          write(*,*)
          write(*,*) 'CFL_z(:,n)  = ',CFL_z(ulevels_nod2D(n):nlevels_nod2D(n),n)
          write(*,*)
          write(*,*) 'hnode(:, n)  = ',hnode(ulevels_nod2D(n):nlevels_nod2D(n),n)
          write(*,*)
          if (use_icebergs) then
            if (.not. allocated(local_idx_of)) then
              allocate(local_idx_of(elem2D))
              call global2local(mesh, partit, local_idx_of, elem2D)
            end if
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                ! global2local zeroes every element this rank does not own, so an
                ! iceberg living on another rank maps to 0 and indexing elem2d_nodes
                ! with it runs off the array.
                if (iceberg_elem(ib) < 1 .or. iceberg_elem(ib) > elem2D) cycle
                if (local_idx_of(iceberg_elem(ib)) == 0) cycle
                if (mesh%elem2d_nodes(1, local_idx_of(iceberg_elem(ib))) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
!$OMP END CRITICAL
       endif
       
       !___________________________________________________________________
       ! check surface vertical velocity --> in case of zlevel and zstar 
       ! vertical coordinate its indicator if Volume is conserved  for 
       ! Wvel(1,n)~maschine preccision
       if ( .not. trim(which_ALE)=='linfs' .and. ( Wvel(1, n) /= Wvel(1, n)  )) then
!$OMP CRITICAL
          write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
          write(*,*) ' --STOP--> found surface layer vertical velocity becomes NaN or >1e-12'
          write(*,*) 'mype     = ',mype
          write(*,*) 'mstep    = ',istep
          write(*,*) 'node     = ',n
          write(*,*) 'uln, nln    = ',ulevels_nod2D(n), nlevels_nod2D(n)
          write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
          write(*,*)
          write(*,*) 'Wvel(1, n)  = ',Wvel(1,n)
          write(*,*) 'Wvel(:, n)  = ',Wvel(:,n)
          write(*,*)
          write(*,*) 'hnode(1, n) = ',hnode(1,n)
          write(*,*) 'hnode(:, n) = ',hnode(:,n)
          write(*,*) 'hflux    = ',heat_flux(n)
          write(*,*) 'wflux    = ',water_flux(n)
          write(*,*)
          write(*,*) 'eta_n    = ',eta_n(n)
          write(*,*) 'd_eta(n)    = ',d_eta(n)
          write(*,*) 'hbar     = ',hbar(n)
          write(*,*) 'hbar_old    = ',hbar_old(n)
          if ( .not. dynamics%use_ssh_se_subcycl) then 
            write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
            write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
          end if
          write(*,*)
          write(*,*) 'CFL_z(:,n)  = ',CFL_z(:,n)
          write(*,*)
          if (use_icebergs) then
            if (.not. allocated(local_idx_of)) then
              allocate(local_idx_of(elem2D))
              call global2local(mesh, partit, local_idx_of, elem2D)
            end if
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                ! global2local zeroes every element this rank does not own, so an
                ! iceberg living on another rank maps to 0 and indexing elem2d_nodes
                ! with it runs off the array.
                if (iceberg_elem(ib) < 1 .or. iceberg_elem(ib) > elem2D) cycle
                if (local_idx_of(iceberg_elem(ib)) == 0) cycle
                if (mesh%elem2d_nodes(1, local_idx_of(iceberg_elem(ib))) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
!$OMP END CRITICAL
       end if ! --> if ( .not. trim(which_ALE)=='linfs' .and. ...
          
       !___________________________________________________________________
       ! check surface layer thinknesss
       if ( .not. trim(which_ALE)=='linfs' .and. ( hnode(1, n) /= hnode(1, n)  .or. hnode(1,n)< 0 )) then
!$OMP CRITICAL
          write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
          write(*,*) ' --STOP--> found surface layer thickness becomes NaN or <0'
          write(*,*) 'mype     = ',mype
          write(*,*) 'mstep    = ',istep
          write(*,*) 'node     = ',n
          write(*,*)
          write(*,*) 'hnode(1, n)  = ',hnode(1,n)
          write(*,*) 'hnode(:, n)  = ',hnode(:,n)
          write(*,*)
          write(*,*) 'eta_n    = ',eta_n(n)
          write(*,*) 'd_eta(n)    = ',d_eta(n)
          write(*,*) 'hbar     = ',hbar(n)
          write(*,*) 'hbar_old    = ',hbar_old(n)
          if ( .not. dynamics%use_ssh_se_subcycl) then 
            write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
            write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
          end if
          write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
          write(*,*)
          if (use_icebergs) then
            if (.not. allocated(local_idx_of)) then
              allocate(local_idx_of(elem2D))
              call global2local(mesh, partit, local_idx_of, elem2D)
            end if
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                ! global2local zeroes every element this rank does not own, so an
                ! iceberg living on another rank maps to 0 and indexing elem2d_nodes
                ! with it runs off the array.
                if (iceberg_elem(ib) < 1 .or. iceberg_elem(ib) > elem2D) cycle
                if (local_idx_of(iceberg_elem(ib)) == 0) cycle
                if (mesh%elem2d_nodes(1, local_idx_of(iceberg_elem(ib))) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
!$OMP END CRITICAL
       end if ! --> if ( .not. trim(which_ALE)=='linfs' .and. ...
          
       
       do nz=ulevels_nod2D(n),nlevels_nod2D(n)-1
          !_______________________________________________________________
          ! check temp
          if ( (tracers%data(1)%values(nz, n) /= tracers%data(1)%values(nz, n)) .or. &
             tracers%data(1)%values(nz, n) < -5.0 .or. tracers%data(1)%values(nz, n)>60) then
!$OMP CRITICAL
             write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
             write(*,*) ' --STOP--> found temperture becomes NaN or <-5.0, >60'
             write(*,*) 'mype     = ',mype
             write(*,*) 'mstep    = ',istep
             write(*,*) 'node     = ',n
             write(*,*) 'lon,lat     = ',geo_coord_nod2D(:,n)/rad
             write(*,*) 'nz       = ',nz
             write(*,*) 'nzmin, nzmax= ',ulevels_nod2D(n),nlevels_nod2D(n)
             write(*,*) 'x=', geo_coord_nod2D(1,n)/rad, ' ; ', 'y=', geo_coord_nod2D(2,n)/rad
             write(*,*) 'temp(nz, n) = ',tracers%data(1)%values(nz, n)
             write(*,*) 'temp(: , n) = ',tracers%data(1)%values(:, n)
             write(*,*) 'temp_old(nz,n)= ',tracers%data(1)%valuesAB(nz, n)
             write(*,*) 'temp_old(: ,n)= ',tracers%data(1)%valuesAB(:, n)
             write(*,*)
             write(*,*) 'hflux    = ',heat_flux(n)
             write(*,*) 'wflux    = ',water_flux(n)
             write(*,*)
             write(*,*) 'eta_n    = ',eta_n(n)
             write(*,*) 'd_eta(n)    = ',d_eta(n)
             write(*,*) 'hbar     = ',hbar(n)
             write(*,*) 'hbar_old    = ',hbar_old(n)
             if ( .not. dynamics%use_ssh_se_subcycl) then 
                write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
                write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
             end if
             write(*,*)
             if (use_ice) then
                write(*,*) 'm_ice    = ',m_ice(n)
                write(*,*) 'm_ice_old   = ',m_ice_old(n)
                write(*,*) 'm_snow      = ',m_snow(n)
                write(*,*) 'm_snow_old  = ',m_snow_old(n)
                write(*,*)
             end if 
             write(*,*) 'hnode    = ',hnode(:,n)
             write(*,*) 'hnode_new   = ',hnode_new(:,n)
             write(*,*)
             write(*,*) 'Kv       = ',Kv(:,n)
             write(*,*)
             write(*,*) 'W          = ',Wvel(:,n)
             write(*,*)
             write(*,*) 'CFL_z(:,n)  = ',CFL_z(:,n)
             write(*,*)
          if (use_icebergs) then
            if (.not. allocated(local_idx_of)) then
              allocate(local_idx_of(elem2D))
              call global2local(mesh, partit, local_idx_of, elem2D)
            end if
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                ! global2local zeroes every element this rank does not own, so an
                ! iceberg living on another rank maps to 0 and indexing elem2d_nodes
                ! with it runs off the array.
                if (iceberg_elem(ib) < 1 .or. iceberg_elem(ib) > elem2D) cycle
                if (local_idx_of(iceberg_elem(ib)) == 0) cycle
                if (mesh%elem2d_nodes(1, local_idx_of(iceberg_elem(ib))) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
             write(*,*)
!$OMP END CRITICAL
          endif ! --> if ( (tracers%data(1)%values(nz, n) /= tracers%data(1)%values(nz, n)) .or. & ...
          
          !_______________________________________________________________
          ! check salt
          ! blowup bounds shifted to anomaly space like the clip in oce_ale_tracer
          ! (S_ref=0 unless use_salt_anomaly)
          if ( (tracers%data(2)%values(nz, n) /= tracers%data(2)%values(nz, n)) .or.  &
             tracers%data(2)%values(nz, n) < 3.0_WP - S_ref_anomaly .or. tracers%data(2)%values(nz, n) > 45.0_WP - S_ref_anomaly ) then
!$OMP CRITICAL
             write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
             write(*,*) ' --STOP--> found salinity becomes NaN or <=3.0, >=45.0'
             write(*,*) 'mype     = ',mype
             write(*,*) 'mstep    = ',istep
             write(*,*) 'node     = ',n
             write(*,*) 'nz       = ',nz
             write(*,*) 'nzmin, nzmax= ',ulevels_nod2D(n),nlevels_nod2D(n)
             write(*,*) 'x=', geo_coord_nod2D(1,n)/rad, ' ; ', 'y=', geo_coord_nod2D(2,n)/rad
             write(*,*) 'salt(nz, n) = ',tracers%data(2)%values(nz, n)
             write(*,*) 'salt(: , n) = ',tracers%data(2)%values(:, n)
             write(*,*)
             write(*,*) 'temp(nz, n) = ',tracers%data(1)%values(nz, n)
             write(*,*) 'temp(: , n) = ',tracers%data(1)%values(:, n)
             write(*,*)
             write(*,*) 'hflux    = ',heat_flux(n)
             write(*,*)
             write(*,*) 'wflux    = ',water_flux(n)
             write(*,*) 'eta_n    = ',eta_n(n)
             write(*,*) 'd_eta(n)    = ',d_eta(n)
             write(*,*) 'hbar     = ',hbar(n)
             write(*,*) 'hbar_old    = ',hbar_old(n)
             if ( .not. dynamics%use_ssh_se_subcycl) then 
                write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
                write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
             end if 
             write(*,*)
             write(*,*) 'hnode    = ',hnode(:,n)
             write(*,*) 'hnode_new   = ',hnode_new(:,n)
             write(*,*)
             write(*,*) 'zbar_3d_n   = ',zbar_3d_n(:,n)
             write(*,*) 'Z_3d_n      = ',Z_3d_n(:,n)
             write(*,*)
             write(*,*) 'Kv       = ',Kv(:,n)
             write(*,*)
             do el=1,nod_in_elem2d_num(n)
                elidx = nod_in_elem2D(el,n)
                 write(*,*) ' elem#=',el,', elemidx=',elidx
                 write(*,*) '    Av =',Av(:,elidx)
             enddo
             write(*,*) 'Wvel     = ',Wvel(:,n)
             write(*,*)
             write(*,*) 'CFL_z(:,n)  = ',CFL_z(:,n)
             write(*,*)
             write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
             write(*,*)
          if (use_icebergs) then
            if (.not. allocated(local_idx_of)) then
              allocate(local_idx_of(elem2D))
              call global2local(mesh, partit, local_idx_of, elem2D)
            end if
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                ! global2local zeroes every element this rank does not own, so an
                ! iceberg living on another rank maps to 0 and indexing elem2d_nodes
                ! with it runs off the array.
                if (iceberg_elem(ib) < 1 .or. iceberg_elem(ib) > elem2D) cycle
                if (local_idx_of(iceberg_elem(ib)) == 0) cycle
                if (mesh%elem2d_nodes(1, local_idx_of(iceberg_elem(ib))) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
!$OMP END CRITICAL
          endif ! --> if ( (tracers%data(2)%values(nz, n) /= tracers%data(2)%values(nz, n)) .or.  & ...
       end do ! --> do nz=1,nlevels_nod2D(n)-1
    end do ! --> do n=1, myDim_nod2d
!$OMP END PARALLEL DO
    !_______________________________________________________________________
    call monitor_fill_ocean(mon, istep, .true., ice, dynamics, tracers, partit, mesh)
    call write_step_info(istep, 1, mon, dynamics, partit)
    if (mype==0) then
        call sleep(1)
        call plot_fesomlogo_expl()
        call plot_fesomlogo_lildevil()

    end if
    call blowup(istep, ice, dynamics, tracers, partit, mesh)
    if (mype==0) write(*,*) ' --> finished writing blow up file'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, abort=1)
end subroutine check_blowup
!===============================================================================
subroutine write_enegry_info(dynamics, partit, mesh)
   IMPLICIT NONE
   type(t_mesh),   intent(in)   , target :: mesh
   type(t_partit), intent(inout), target :: partit
   type(t_dyn)   , intent(in)   , target :: dynamics
   real(kind=WP)                         :: budget(2)
   integer, pointer                      :: mype

   mype            => partit%mype

   if (mype==0) write(*,*) '*******KE budget analysis...*******'
   if (mype==0) write(*,*) '     U     |     V     |     TOTAL'
   call integrate_elem(dynamics%ke_du2(1,:,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_du2(2,:,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. du2=', budget(1),' | ',  budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_pre_xVEL(1,:,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_pre_xVEL(2,:,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. pre=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_nod(dynamics%ke_wrho,  budget(1), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7)") 'w * rho=', budget(1)

   call integrate_elem(dynamics%ke_adv_xVEL(1,:,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_adv_xVEL(2,:,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. adv=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_hvis_xVEL(1,:,:), budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_hvis_xVEL(2,:,:), budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke.  ah=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_vvis_xVEL(1,:,:), budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_vvis_xVEL(2,:,:), budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke.  av=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_cor_xVEL(1,:,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_cor_xVEL(2,:,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. cor=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_wind_xVEL(1,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_wind_xVEL(2,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. wind=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_drag_xVEL(1,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_drag_xVEL(2,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. drag=', budget(1), ' | ', budget(2), ' | ', sum(budget)
   if (mype==0) write(*,*) '***********************************'   
end subroutine write_enegry_info

!
!
!_______________________________________________________________________________
subroutine plot_fesomlogo()
    implicit none 
    character(len=*), parameter :: r = char(27)//'[31m' ! red
    character(len=*), parameter :: c = char(27)//'[36m' ! cyan
    character(len=*), parameter :: z = char(27)//'[0m'    ! reset
    write(*,*)
    write(*,*) c//' .------. ------. ------.  .-----. ,--.   ,-.  .----.  '//z
    write(*,*) c//' |  .--- ́|  .--- ́/    _ / ´  .-.  `|   `. ́   |\_,-.  | '//z
    write(*,*) c//' |  `--. |  `--. \_..`--. |  | |  ||  |`. ́|  |   . ́ . ́ '//z
    write(*,*) c//' |  .-- ́ |  .-- ́ .-._)   \|  | |  ||  |   |  | . ́  /_  '//z
    write(*,*) c//' |  |    |  `---.\       /`  `- ́   ́|  |   |  ||      | '//z
    write(*,*) c//' `-- ́    `------ ́ `----- ́  `----- ́ `-- ́   `-- ́`------ ́ '//z
    write(*,*) '                                          _____           '
    write(*,*) '      ___________                     ,-:` \;´,``-,       '
    write(*,*) '     | .-------. |                  .´-;_,;  `:-;_,`.     '
    write(*,*) '     | |       | |                 /;   `/    ,  _`.-\    ' 
    write(*,*) '     | |       | |                | ´`. (`     /` ` \`|   '
    write(*,*) '     | |__   __| |            ,-C=|:.  `\`-.   \_   / |   '
    write(*,*) '      `---|-|---´          ,-´    |     (   `,  .`\ ;`|   ' 
    write(*,*) '       [====  O]--,    ,--´        \     | .´     `-´/    '
    write(*,*) '     /:::::::::::\ \_,-             `.   ;/        .´     '
    write(*,*) '    /:::::===:::::\                   ``-._____.-´`       '
    write(*,*) '   ´---------------`                                      '
    write(*,*)
end subroutine plot_fesomlogo

!
!
!_______________________________________________________________________________
subroutine plot_fesomlogo_lildevil()
    implicit none 
    character(len=*), parameter :: r = char(27)//'[31m' ! red
    character(len=*), parameter :: c = char(27)//'[36m' ! cyan
    character(len=*), parameter :: z = char(27)//'[0m'    ! reset
    ! character(len=*), parameter :: g = char(27)//'[32m' ! green
    ! character(len=*), parameter :: o = char(27)//'[33m' ! orange
    ! character(len=*), parameter :: b = char(27)//'[34m' ! blue
    ! character(len=*), parameter :: p = char(27)//'[35m' ! purple
    ! write(*,*) '            (`- ́)  _ (`- ́).->          <-. (`- ́)          ' 
    ! write(*,*) '   <-.      ( OO).-/ ( OO)_      .->      \(OO )_         '
    ! write(*,*) '(`- ́)-----.(,------.(_)--\_)(`- ́)----. ,--./  ,-.) .----. '
    ! write(*,*) '(OO|(_\--- ́ |  .--- ́/    _ /( OO).-.  `|   `. ́   |\_,-.  |'
    ! write(*,*) ' / |  `--. (|  `--. \_..`--.( _) | |  ||  |`. ́|  |   . ́ . ́'
    ! write(*,*) ' \_)  .-- ́  |  .-- ́ .-._)   \\|  |)|  ||  |   |  | . ́  /_ '
    ! write(*,*) '  `|  |_)   |  `---.\       / `  `- ́   ́|  |   |  ||      |'
    ! write(*,*) '   `-- ́     `------ ́ `----- ́   `----- ́ `-- ́   `-- ́`------ ́'
    ! write(*,*)
    write(*,*) '                                                          '
    write(*,*) '            '//r//'(`- ́)  _ (`- ́).->          <-. (`- ́)'//z//'          ' 
    write(*,*) '   '//r//'<-.      ( oo).-/ ( oo)_      .->      \(oo )_'//z//'         '
    write(*,*) r//'(`- ́)'//c//'-----.'//r//'('//c//',------.'//r//'(_)'//c//'--'//r//'\_)(`- ́)'//c//'----. ,--.'//r//'/'//c//'  ,-.'//r//')'//c//' .----. '//z
    write(*,*) r//'(oo'//c//'|'//r//'(_\'//c//'--- ́ |  .--- ́/    _ /'//r//'( oo)'//c//'.-.  `|   `. ́   |\_,-.  |'//z
    write(*,*) r//' / '//c//'|  `--. '//r//'('//c//'|  `--. \_..`--.'//r//'( _)'//c//' | |  ||  |`. ́|  |   . ́ . ́'//z
    write(*,*) r//' \_)'//c//'  .-- ́  |  .-- ́ .-._)   \'//r//'\'//c//'|  |'//r//')'//c//'|  ||  |   |  | . ́  /_ '//z
    write(*,*) r//'  `'//c//'|  |'//r//'_)'//c//'   |  `---.\       / `  `- ́   ́|  |   |  ||      |'//z
    write(*,*) c//'   `-- ́     `------ ́ `----- ́   `----- ́ `-- ́   `-- ́`------ ́'//z
    write(*,*) '                                                          '
    write(*,*)
end subroutine plot_fesomlogo_lildevil 

!
!
!_______________________________________________________________________________
subroutine plot_fesomlogo_expl()
    implicit none 
    character(len=*), parameter :: r = char(27)//'[31m' ! red
    character(len=*), parameter :: g = char(27)//'[32m' ! green
    character(len=*), parameter :: o = char(27)//'[33m' ! orange
    character(len=*), parameter :: b = char(27)//'[34m' ! blue
    character(len=*), parameter :: p = char(27)//'[35m' ! purple
    character(len=*), parameter :: c = char(27)//'[36m' ! cyan
    character(len=*), parameter :: z = char(27)//'[0m'  ! reset
    ! write(*,*)
    ! write(*,*) '                                                          '
    ! write(*,*) '        YOUR MODEL BLOW UP, CODE GREMLIN IN ACTION !!!    '
    ! write(*,*) '                              ____                        '
    ! write(*,*) '                       __,-~~/~   `---.                   '
    ! write(*,*) '                     _/_,---(      ,   )                  '
    ! write(*,*) '                 __ /        <   /   )   \___             '
    ! write(*,*) ' - -- ----===;;;`====------------------===;;;===---- -- - '
    ! write(*,*) '.____________.      \/  ~"~"~"~"~"~\~"~)~"/               '
    ! write(*,*) '|Code Gremlin|      (_ (   \  (     >    \)               '
    ! write(*,*) '`----+-------´        \_( _ <         >_>`                 '
    ! write(*,*) '     v                  ~ `-i` ::>|--"                    '
    ! write(*,*) '   (`-´)                    I;|.|.|                       '
    ! write(*,*) '  _( o<)_  __T__           <|i::|i|`                      '
    ! write(*,*) ' (,_ w _,) |TNT|          (` ^`,-* ")                     '
    ! write(*,*) '<-´()^()   |___|\_______.,-#%&(_)%#&#~,.                  '
    ! write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)    '                                                          '
    write(*,*)    '        YOUR MODEL BLOW UP, CODE GREMLIN IN ACTION !!!    '
    write(*,*) b//'                              ____                        '//z
    write(*,*) b//'                       __,-~~/~   `---.                   '//z
    write(*,*) b//'                     _/_,---(      ,   )                  '//z
    write(*,*) b//'                 __ /        <   /   )   \___             '//z
    write(*,*) z//' - -- ----===;;;`====------------------===;;;===---- -- - '//z
    write(*,*) z//'.____________.      \/  ~"~"~"~"~"~\~"~)~"/               '//z
    write(*,*) z//'|Code Gremlin|      (_ (   \  (     >    \)               '//z
    write(*,*) z//'`----+-------´'//o//'        \_( _ <         >_>`                 '//z
    write(*,*) z//'     v        '//o//'           ~ `-i` ::>|--"                    '//z
    write(*,*) r//'   (`-´)  '//z//'              '//o//'    I;|.|.|                       '//z
    write(*,*) r//'  _( o<)_ '//z//' __T__        '//r//'   <|i::|i|`                      '//z
    write(*,*) r//' (,_ w _,)'//z//' |TNT|        '//r//'  (` ^`'//z//',-*'//r//' ")                     '//z
    write(*,*) r//'<-´()^()  '//z//' |___|\_______'//p//'.,-#%&'//z//'(_)'//p//'%#&#~,.                  '//z
    write(*,*)
    write(*,*)
    write(*,*)    ' Things to try (if you havent changed anything in the code)'
    write(*,*)    '     - 1st. increase the step_per_day (namelist.config,    '
    write(*,*)    '       steps_per_day=32(45min),36(40min)...,48(30min)...'
    write(*,*)    '       ...,60(24min)...,72(20min)...,96(15min)...,144(10min)'
    write(*,*)    '       ...,160(9min),180(8min)...,240(6min)...,288(5min)'
    write(*,*)    '       ...,360(4min)...,480(3min)...,720(2min)...,1440(1min))'
    write(*,*)    '     - 2nd. increase slowly background viscosity visc_gamma0 within'
    write(*,*)    '       its bounds (namelist.dyn)'
    write(*,*)    '     - 3nd. increase slowly flow aware viscosity visc_gamma1 within'
    write(*,*)    '       its bounds (namelist.dyn)'
    write(*,*)    '     - 4th. contact developers :-D'  
    write(*,*)
end subroutine plot_fesomlogo_expl
 

end module write_step_info_module
