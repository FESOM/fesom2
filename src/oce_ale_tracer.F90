module oce_ale_tracer_module
    USE g_config
    USE o_PARAM
    USE mod_mesh
    USE MOD_PARTIT
    use par_support_module, only: par_abort
    USE MOD_DYN
    USE MOD_ICE
    USE mod_tracer
    USE g_comm_auto
    USE o_tracers
    USE Toy_Channel_Soufflet
    USE Toy_Channel_Dbgyre
    USE Toy_Neverworld2
    USE o_ARRAYS
    USE g_forcing_arrays
    USE oce_adv_tra_driver_module, only: do_oce_adv_tra
    USE fesom_conservation_module, only: conservation_probe
    USE oce_dyn_module, only: mom_on_device
#if defined (FESOM_PROFILING)
    use fesom_profiler
#endif
    USE diagnostics, only: ldiag_DVD, ldiag_diapmix, density_dmoc_avg, diap_avg_count, &
            dmoc_avg_count, dmoc_is_due, dT_diap, dS_diap, dd_diap
    USE g_forcing_param, only: use_age_tracer
    USE mod_transit
    USE cmor_variables_diag, only: ldiag_cmor, save_cmor_advection
    USE g_clock
    USE o_mixing_KPP_mod
    USE iceberg_params
#if defined(__recom)
    use recom_glovar
    use recom_config
    use recom_ciso
    use o_arrays
#endif
#if defined(__recom)
    use recom_sinking
    use recom_glovar
    use recom_config
    use g_comm_auto
    use g_support
#endif
#if defined (__cvmix)
    use g_cvmix_kpp, only: kpp_nonlcltranspT, kpp_nonlcltranspS, kpp_oblmixc
#endif
#if defined (__recom)
   use recoM_declarations, only: is_erosioninput, is_riverinput
   use recom_glovar
   use recom_config
#endif
#if defined (__ciso)
   use recom_ciso
#endif
    implicit none

    private
    public :: solve_tracers_ale, diff_tracers_ale, &
              diff_ver_part_expl_ale, diff_ver_part_impl_ale, &
              diff_ver_part_redi_expl, diff_part_hor_redi, &
              diff_part_bh, bc_surface, transit_bc_surface, calc_slice

    ! OpenACC: the diffusion of this tracer runs on the device; .false. when a diagnostic or recom
    ! reads del_ttf or values on the host in between (then the whole diffusion runs on the host)
    logical :: diff_on_device = .false.
    ! ... and its implicit vertical diffusion (diff_ver_part_impl_ale): T and S without the
    ! features that need host data inside the column loop (icebergs, depth hosing, toy forcing)
    logical :: impl_on_device = .false.
    ! device copy of T (1) and S (2) equal to the host copy at the end of diff_tracers_ale
    logical :: tr_dev_current(2) = .false.
    ! device-only advecting velocity UV + bolus velocity (GM on the device)
    real(kind=WP), allocatable, target :: uv_tr(:,:,:), w_tr(:,:), we_tr(:,:)
    ! column work arrays of diff_ver_part_impl_ale have this fixed size: as private arrays of
    ! run-time size they would come from the per-thread device heap
    integer, parameter :: impl_nl_max = 128
    ! surface boundary term per node (and per level for the depth hosing tracer), computed
    ! before the column loop of diff_ver_part_impl_ale
    real(kind=WP), allocatable :: impl_bc(:), impl_bcz(:,:)
    ! device copies of the inputs of diff_ver_part_impl_ale, created on the first call
    logical :: impl_inputs_created = .false.
    ! horizontal tracer gradient at nodes in diff_ver_part_redi_expl; allocated (and created
    ! on the device) once: an automatic array in a device data region is not reliably present
    real(kind=WP), allocatable :: tr_xynodes(:,:,:)

contains

!
!
!===============================================================================
! Driving routine    Here with ALE changes!!!
subroutine solve_tracers_ale(ice, dynamics, tracers, partit, mesh)
    implicit none
    type(t_ice)   , intent(in)   , target    :: ice
    type(t_dyn)   , intent(inout), target    :: dynamics
    type(t_tracer), intent(inout), target    :: tracers
    type(t_partit), intent(inout), target    :: partit
    type(t_mesh)  , intent(in)   , target    :: mesh

#if defined(__recom) && defined(__usetp)
! multi FESOM group loop parallelization
    integer             :: num_tracers
    integer             :: tr_num_start_memo

    integer             :: group_i
    integer             :: tr_num_start

    logical             :: has_one_added_tracer
    logical             :: has_one_added_tracer_local_dummy
    integer             :: tr_num_end_local_dummy
    integer             :: tr_num_in_group_local_dummy
    integer             :: tr_num_end
    integer             :: tr_num_in_group_dummy
    integer             :: tr_arr_slice_count_fix_1

    integer             :: Sinkflx_tr_slice_count_fix_1
    integer             :: Benthos_tr_slice_count_fix_1

    integer             :: tr_num_start_local
    integer             :: tr_num_to_send

    logical             :: completed

    logical             :: bBreak
#endif

    !___________________________________________________________________________
    integer                                  :: i, tr_num, node, elem, nzmax, nzmin
    real(kind=WP)                            :: ttf_rhs_bak (mesh%nl-1, partit%myDim_nod2D+partit%eDim_elem2D) ! local variable
    integer                                  :: nz, n, nu1, nl1
    real(kind=WP), dimension(:,:), pointer   :: salt_dev
    !___________________________________________________________________________
    ! pointer on necessary derived types
    logical :: gm_on_device = .false.
    integer :: nlu, nlw
    real(kind=WP), dimension(:,:,:), pointer :: UV, fer_UV, UVadv
    real(kind=WP), dimension(:,:)  , pointer :: Wvel, Wvel_e, Wvel_i, fer_Wvel, Wadv, Weadv
    real(kind=WP_full), dimension(:,:)  , pointer :: del_ttf
    UV     => dynamics%uv(:,:,:)
    UVadv  => UV
    Wadv   => dynamics%w(:,:)
    Weadv  => dynamics%w_e(:,:)
    Wvel   => dynamics%w(:,:)
    Wvel_e => dynamics%w_e(:,:)
    Wvel_i => dynamics%w_i(:,:)
    if (Fer_GM) then
        fer_UV     => dynamics%fer_uv(:,:,:)
        fer_Wvel   => dynamics%fer_w(:,:)
    end if
    del_ttf => tracers%work%del_ttf

#if defined(__recom) && defined(__usetp)
    num_tracers=tracers%num_tracers
#endif

    !___________________________________________________________________________
    if (SPP) then
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call cal_rejected_salt'//achar(27)//'[0m'
        call cal_rejected_salt(ice, partit, mesh)
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call app_rejected_salt'//achar(27)//'[0m'
        call app_rejected_salt(tracers%data(2)%values, partit, mesh)
    end if 

    !___________________________________________________________________________
    ! update 3D velocities with the bolus velocities:
    ! 1. bolus velocities are computed according to GM implementation after R. Ferrari et al., 2010
    ! 2. bolus velocities are used only for advecting tracers and shall be subtracted back afterwards
#ifdef ENABLE_OPENACC
    ! with zstar and the momentum on the device the bolus velocities (computed or uploaded in
    ! vert_vel_ale_acc) are added on the device below and the host copies stay untouched
    gm_on_device = Fer_GM .and. trim(which_ALE)=='zstar' .and. mom_on_device
#endif
    if (Fer_GM .and. .not. gm_on_device) then
        uv_dev_current = .false.
!$OMP PARALLEL DO
        do elem=1, partit%myDim_elem2D+partit%eDim_elem2D
           UV(:, :, elem)    =UV(:, :, elem) + fer_UV(:, :, elem)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do node=1, partit%myDim_nod2D+partit%eDim_nod2D
           Wvel_e(:, node)=Wvel_e(:, node)+fer_Wvel(:, node)
           Wvel  (:, node)=Wvel  (:, node)+fer_Wvel(:, node)
        end do
!$OMP END PARALLEL DO
    end if

    ! Set advective and diffusive components of total tracer fluxes to zero
    ! Before tr_num loop
!#if defined (__recom)           ! not necessarily should belong to recom case
!    tracers%work%tra_advhoriz = 0.0 ! O:G - tra_diag
!    tracers%work%tra_advvert  = 0.0
    ttf_rhs_bak = 0.0
!#endif
    
    !___________________________________________________________________________
    ! loop over all tracers
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("tr_uploads")
#endif
        ! after vert_vel_ale_acc w, w_e, UV, helem, zbar_3d_n and hnode_new are current on the device
        if (.not. vv_on_device) then
            !$ACC UPDATE DEVICE(dynamics%w, dynamics%w_e, dynamics%uv)
            !$ACC UPDATE DEVICE (mesh%helem, mesh%hnode_new, mesh%zbar_3d_n)
        end if
!!!     !$ACC UPDATE DEVICE(tracers%work%fct_ttf_min, tracers%work%fct_ttf_max, tracers%work%fct_plus, tracers%work%fct_minus)
        if (.not. geom_on_device) then
            !$ACC UPDATE DEVICE (mesh%hnode, mesh%z_3d_n)
        end if
#ifdef ENABLE_OPENACC
        if (gm_on_device) then
            ! the advecting velocity UV + bolus velocity goes to a device-only array, so
            ! UV itself stays unchanged on the device (no upload in the next step)
            nlu = size(UV,2)
            nlw = size(Wvel,1)
            if (.not. allocated(uv_tr)) then
                allocate(uv_tr(2, nlu, partit%myDim_elem2D+partit%eDim_elem2D))
                !$ACC ENTER DATA CREATE(uv_tr)
            end if
            !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(nz)
            do elem=1, partit%myDim_elem2D+partit%eDim_elem2D
                do nz=1, nlu
                    uv_tr(1, nz, elem)=UV(1, nz, elem) + fer_UV(1, nz, elem)
                    uv_tr(2, nz, elem)=UV(2, nz, elem) + fer_UV(2, nz, elem)
                end do
            end do
            !$ACC END PARALLEL LOOP
            UVadv => uv_tr
            ! the same for w and w_e: w stays the resolved velocity on the device (output)
            if (.not. allocated(w_tr)) then
                allocate(w_tr(nlw, partit%myDim_nod2D+partit%eDim_nod2D), we_tr(nlw, partit%myDim_nod2D+partit%eDim_nod2D))
                !$ACC ENTER DATA CREATE(w_tr, we_tr)
            end if
            !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(nz)
            do node=1, partit%myDim_nod2D+partit%eDim_nod2D
                do nz=1, nlw
                    we_tr(nz, node)=Wvel_e(nz, node)+fer_Wvel(nz, node)
                    w_tr (nz, node)=Wvel  (nz, node)+fer_Wvel(nz, node)
                end do
            end do
            !$ACC END PARALLEL LOOP
            Wadv  => w_tr
            Weadv => we_tr
        end if
#endif
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("tr_uploads")
#endif
        ! used by the diffusion; computed on the host unless GM/Redi ran on the device
        if (.not. redi_on_device) then
            !$ACC UPDATE DEVICE (Ki, slope_tapered)
        end if
#ifdef ENABLE_OPENACC
        ! inputs of the implicit vertical diffusion, computed on the host (mixing, forcing,
        ! vertical velocity split); unallocated ones (e.g. the KPP arrays of another
        ! scheme) are skipped by the data directives
        if (.not. impl_inputs_created) then
            !$ACC ENTER DATA CREATE(Kv, sw_3d, heat_flux, water_flux, virtual_salt, relax_salt, real_salt_flux)
            !$ACC ENTER DATA CREATE(dynamics%w_i, ghats, blmc)
#if defined (__cvmix)
            !$ACC ENTER DATA CREATE(kpp_nonlcltranspT, kpp_nonlcltranspS, kpp_oblmixc)
#endif
            impl_inputs_created = .true.
        end if
        !$ACC UPDATE DEVICE(sw_3d, heat_flux, water_flux, virtual_salt, relax_salt, real_salt_flux)
        ! with KPP on the device Kv is computed there (mo_convect_acc)
        if (.not. kpp_on_device) then
            !$ACC UPDATE DEVICE(Kv)
        end if
        !$ACC UPDATE DEVICE(dynamics%w_i)
        ! with KPP on the device ghats and blmc are current there
        if (.not. kpp_on_device) then
            !$ACC UPDATE DEVICE(ghats, blmc)
        end if
#if defined (__cvmix)
        !$ACC UPDATE DEVICE(kpp_nonlcltranspT, kpp_nonlcltranspS, kpp_oblmixc)
#endif
#endif

#if defined(__recom) && defined(__usetp)
    call calc_slice(num_tracers, num_fesom_groups, partit%my_fesom_group, tr_num_start, tr_num_end, tr_num_in_group_dummy, has_one_added_tracer)

    tr_arr_slice_count_fix_1 = 1 * (mesh%nl - 1) * (partit%myDim_nod2D + partit%eDim_nod2D)

    Sinkflx_tr_slice_count_fix_1 = 1 * (partit%myDim_nod2D + partit%eDim_nod2D) * bottflx_num
    Benthos_tr_slice_count_fix_1 = 1 * (partit%myDim_nod2D + partit%eDim_nod2D) * benthos_num

    tr_num_start_memo = tr_num_start

    request_count = 0
#endif

#if defined(__recom)
    ! CMIP7 sinking diagnostics (REcoM ExpFlux, BotFlux) are rates of this time step,
    ! summed over the tracers in ver_sinking_recom and ver_sinking_recom_benthos below
    if (Diags) then
        ExpFlux = 0.0_WP
        BotFlux = 0.0_WP
    end if
#endif

#if defined(__recom) && defined(__usetp)
    do tr_num = tr_num_start, tr_num_end
#else
    do tr_num=1, tracers%num_tracers
#endif

#if defined(__recom)
    if(use_MEDUSA) then
            SinkFlx = 0.0d0
#if defined(__usetp)
            SinkFlx_tr(:, :, tr_num) = 0.0d0
#endif !__usetp
    endif
#if defined(__usetp)
    Benthos_tr(:, :, tr_num) = 0.0d0
#endif !__usetp
#endif !__recom

#if defined(__recom)
        SinkingVel1 = 0.0d0 ! OG 16.03.23
        SinkingVel2 = 0.0d0 ! OG 16.03.23
#endif    

        ! do tracer AB (Adams-Bashfort) interpolation only for advectiv part
        ! needed
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call init_tracers_AB'//achar(27)//'[0m'
        ! values were last changed on the host (diffusion, surface fluxes); valuesold
        ! and valuesAB stay on the device, where init_tracers_AB updates them.
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("tr_init_ab")
#endif
        !$ACC UPDATE DEVICE(tracers%data(tr_num)%values)
        call init_tracers_AB(tr_num, tracers, partit, mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("tr_init_ab")
#endif
        if (conservation_freq > 0) call conservation_probe('start', 'a', tr_num, tracers, partit, mesh)
 
        ! advect tracers
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call adv_tracers_ale'//achar(27)//'[0m'
	!here update only those initialized in the init_tracers. (values, valuesAB, edge_up_dn_grad, ...)
        ! it will update del_ttf with contributions from horizontal and vertical advection parts (del_ttf_advhoriz and del_ttf_advvert)
	!$ACC wait(1)
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("tr_adv")
#endif
        call do_oce_adv_tra(dt, UVadv, Wadv, Wvel_i, Weadv, tr_num, dynamics, tracers, partit, mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("tr_adv")
#endif

        !___________________________________________________________________________
        ! update array for total tracer flux del_ttf with the fluxes from horizontal
        ! and vertical advection
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
        do node=1, partit%myDim_nod2d
           tracers%work%del_ttf(:, node)=tracers%work%del_ttf(:, node)+tracers%work%del_ttf_advhoriz(:, node)+tracers%work%del_ttf_advvert(:, node)
        end do
!$OMP END PARALLEL DO
#else
        ! del_ttf stays on the device until the tracer update in diff_tracers_ale
!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(PRESENT)
        do node=1, partit%myDim_nod2d
           do nz=1, mesh%nl-1
              tracers%work%del_ttf(nz, node)=tracers%work%del_ttf(nz, node)+tracers%work%del_ttf_advhoriz(nz, node)+tracers%work%del_ttf_advvert(nz, node)
           end do
        end do
!$ACC END PARALLEL LOOP
#endif

        !___________________________________________________________________________
        ! Save advection snapshot for CMOR tendency decomposition
        if (ldiag_cmor .and. (tr_num <= 2)) then
            !$ACC UPDATE HOST(tracers%work%del_ttf)
            call save_cmor_advection(tr_num, tracers%work%del_ttf(:, 1:partit%myDim_nod2D), partit%myDim_nod2D, mesh%nl-1)
        end if

!if (.FALSE.) then
! O:G - tra_diag
!#if defined (__recom)
!        if (tracers%data(tr_num)%ltra_diag) then
!           do n=1, myDim_nod2D+eDim_nod2D
!              nu1 = ulevels_nod2D(n)
!              nl1 = nlevels_nod2D(n)
!              do nz = nu1, nl1-1
                 ! Horizontal advection part
!                 tracers%work%tra_advhoriz(nz,n,tr_num) = tracers%work%del_ttf_advhoriz(nz,n)
                 ! Vertical advection part
!                 tracers%work%tra_advvert (nz,n,tr_num) = tracers%work%del_ttf_advvert(nz,n)
!              end do
!           end do
!        end if
!#endif
!endif 

        !___________________________________________________________________________
        ! diffuse tracers
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call diff_tracers_ale'//achar(27)//'[0m'
        if (conservation_freq > 0) call conservation_probe('advection', 'a', tr_num, tracers, partit, mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("tr_diff")
#endif
        call diff_tracers_ale(tr_num, dynamics, tracers, ice, partit, mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("tr_diff")
#endif

        !___________________________________________________________________________
        ! Radioactive decay of 14C and 39Ar
        if (tracers%data(tr_num)%ID == 14) tracers%data(tr_num)%values(:,:) = tracers%data(tr_num)%values(:,:) * exp(-decay14 * dt)
        if (tracers%data(tr_num)%ID == 39) tracers%data(tr_num)%values(:,:) = tracers%data(tr_num)%values(:,:) * exp(-decay39 * dt)

!YY: C14 seems to be calculated both in fesom and recom
!YY: decay differently calculated???
#if defined(__ciso)
        ! radioactive decay of 14C
        if (ciso_14 .and. any(c14_tracer_id == tracers%data(tr_num)%ID)) then
          tracers%data(tr_num)%values(:,:) = tracers%data(tr_num)%values(:,:) * (1 - lambda_14 * dt)
        end if    ! ciso & ciso_14
#endif

        !___________________________________________________________________________
        ! relax to salt and temp climatology
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call relax_to_clim'//achar(27)//'[0m'
        ! if ((toy_ocean) .AND. ((tr_num==1) .AND. (TRIM(which_toy)=="soufflet"))) then
        if     ((toy_ocean) .AND. ((TRIM(which_toy)=="soufflet"))) then
            call relax_zonal_temp(tracers%data(1), partit, mesh)
            
        elseif ((toy_ocean) .AND. ((TRIM(which_toy)=="neverworld2"))) then
            call relax_2_tsurf(tracers%data(1), partit, mesh)
            call relax_2_ssurf(tracers%data(2), partit, mesh)

        else
            call relax_to_clim(tr_num, tracers, partit, mesh)
            
        end if
        if (conservation_freq > 0) call conservation_probe('relaxation', 'b', tr_num, tracers, partit, mesh)

#if defined (FESOM_PROFILING)
        call fesom_profiler_start("tr_exchange")
#endif
        call exchange_nod(tracers%data(tr_num)%values(:,:), partit)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("tr_exchange")
#endif
!$OMP BARRIER

!!!        !$ACC UPDATE HOST (tracers%work%fct_ttf_min, tracers%work%fct_ttf_max, tracers%work%fct_plus, tracers%work%fct_minus) &
!!!        !$ACC HOST  (tracers%work%edge_up_dn_grad)

#if defined(__recom) && defined(__usetp)
! broadcast tracer results to fesom groups
        if(num_fesom_groups > 1) then

            do group_i = 0, num_fesom_groups - 1
                call calc_slice(num_tracers, num_fesom_groups, group_i, tr_num_start_local, tr_num_end_local_dummy, tr_num_in_group_local_dummy, has_one_added_tracer_local_dummy)

                tr_num_to_send = tr_num_start_local + (tr_num - tr_num_start_memo)

                if((tr_num == tr_num_end) .and. has_one_added_tracer) then
                    ! skip: if last tracer in group was added to compensate for fragementation it is skipped here and handled after the loop
                else
                    request_count = request_count + 1

! non-blocking communication overlapped with computation in loop
                    call MPI_IBcast(tracers%data(tr_num_to_send)%values(:, :), tr_arr_slice_count_fix_1, MPI_WP, &
                                                group_i, partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, tr_arr_requests(request_count),     partit%MPIerr)

                    if(use_MEDUSA) then
                        call MPI_IBcast(Sinkflx_tr (:, :, tr_num_to_send), Sinkflx_tr_slice_count_fix_1, MPI_WP, &
                                                    group_i, partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, SinkFlx_tr_requests(request_count), partit%MPIerr)
                    endif
                        call MPI_IBcast(Benthos_tr (:, :, tr_num_to_send), Benthos_tr_slice_count_fix_1, MPI_WP, &
                                                    group_i, partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, Benthos_tr_requests(request_count), partit%MPIerr)
                end if
            end do
        end if ! (num_fesom_groups > 1) then
#endif
    end do ! EITHER: tr_num = tr_num_start, tr_num_end OR 1, tracers%num_tracers, depending on __usetp
    
#if defined(__recom) && defined(__usetp)
! if tracer in group was added to compensate for fragmentation its broadcast of the last index is handled here
    if(num_fesom_groups > 1) then
        do group_i = 0, num_fesom_groups - 1
            call calc_slice(num_tracers, num_fesom_groups, group_i, tr_num_start, tr_num_end, tr_num_in_group_dummy, has_one_added_tracer)

            if(has_one_added_tracer) then

                request_count = request_count + 1

                call MPI_IBcast(tracers%data(tr_num_end)%values(:, :), tr_arr_slice_count_fix_1, MPI_WP, &
                                group_i, partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, tr_arr_requests(request_count),     partit%MPIerr)
                if(use_MEDUSA) then
                    call MPI_IBcast(Sinkflx_tr (:, :, tr_num_end), Sinkflx_tr_slice_count_fix_1, MPI_WP, &
                                    group_i, partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, SinkFlx_tr_requests(request_count), partit%MPIerr)
                endif
                    call MPI_IBcast(Benthos_tr (:, :, tr_num_end), Benthos_tr_slice_count_fix_1, MPI_WP, &
                                    group_i, partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, Benthos_tr_requests(request_count), partit%MPIerr)
            end if
        end do
    end if !(num_fesom_groups > 1) then

    if(num_fesom_groups > 1) then
        completed = .false.
        do while (.not. completed)
            call MPI_TESTALL(request_count, tr_arr_requests(:),     completed, MPI_STATUSES_IGNORE, partit%MPIerr)
        end do

        if(use_MEDUSA) then
            completed = .false.
            do while (.not. completed)
                call MPI_TESTALL(request_count, SinkFlx_tr_requests(:), completed, MPI_STATUSES_IGNORE, partit%MPIerr)
            end do
        endif ! (use_MEDUSA) then

            completed = .false.
            do while (.not. completed)
                call MPI_TESTALL(request_count, Benthos_tr_requests(:), completed, MPI_STATUSES_IGNORE, partit%MPIerr)
            end do
    end if ! (num_fesom_groups > 1) then
#endif

#if defined(__recom) && defined(__usetp)
! SinkFlx and Benthos values are buffered per tracer index in the loop above and now summed up to
! avoid non bit identical results regarding global sums when running the tracer loop in parallel
        do tr_num = 1, num_tracers
            if(use_MEDUSA) then
                SinkFlx = SinkFlx + SinkFlx_tr(:, :, tr_num)
            endif
            Benthos = Benthos + Benthos_tr(:, :, tr_num)
!            Sinkingvel1(:,:) = Sinkingvel1(:,:) + Sinkvel1_tr(:, :, tr_num)
!            Sinkingvel2(:,:) = Sinkingvel2(:,:) + Sinkvel2_tr(:, :, tr_num)
        end do

! Each tracer is stepped in one FESOM group only, so the sinking diagnostics hold the
! contributions of this group's tracers: sum them over the groups
        if (Diags .and. num_fesom_groups > 1) then
            call MPI_Allreduce(MPI_IN_PLACE, ExpFlux, size(ExpFlux), MPI_WP, MPI_SUM, &
                               partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, partit%MPIerr)
            call MPI_Allreduce(MPI_IN_PLACE, BotFlux, size(BotFlux), MPI_WP, MPI_SUM, &
                               partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, partit%MPIerr)
        end if
#endif

    !___________________________________________________________________________
    ! 3D restoring for "passive" tracers
    !!!$OMPTODO: add OpenMP later, not needed right now!
    do tr_num=1, ptracers_restore_total
       tracers%data(ptracers_restore(tr_num)%locid)%values(:, ptracers_restore(tr_num)%ind2)=1.0_WP
    end do

    !___________________________________________________________________________
    ! subtract the the bolus velocities back from 3D velocities:
    ! (the device copies of UV and w are unchanged: the sums are in uv_tr, w_tr, we_tr)
    if (Fer_GM .and. .not. gm_on_device) then
        uv_dev_current = .false.
!$OMP PARALLEL DO
        do elem=1, partit%myDim_elem2D+partit%eDim_elem2D
           UV(:, :, elem)    =UV(:, :, elem) - fer_UV(:, :, elem)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do node=1, partit%myDim_nod2D+partit%eDim_nod2D
           Wvel_e(:, node)=Wvel_e(:, node)-fer_Wvel(:, node)
           Wvel  (:, node)=Wvel  (:, node)-fer_Wvel(:, node)
        end do
!$OMP END PARALLEL DO
    end if
    
    ! TODO: do it only when it is coupled to atmosphere 
    !___________________________________________________________________________
    ! to avoid crash with high salinities when coupled to atmosphere
    ! --> if we do only where (tr_arr(:,:,2) < 3._WP ) we also fill up the bottom
    !     topogrpahy with values which are then writte into the output --> thats why
    !     do node=1,.... and tr_arr(node,1:nzmax,2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nzmin, nzmax)
    do node=1,partit%myDim_nod2D+partit%eDim_nod2D
        nzmax=mesh%nlevels_nod2D(node)-1
        nzmin=mesh%ulevels_nod2D(node)
        ! clip bounds shifted to anomaly space (S_ref=0 unless use_salt_anomaly)
        where (tracers%data(2)%values(nzmin:nzmax,node) > 45._WP - S_ref_anomaly)
               tracers%data(2)%values(nzmin:nzmax,node)= 45._WP - S_ref_anomaly
        end where
        where (tracers%data(2)%values(nzmin:nzmax,node) < 3._WP - S_ref_anomaly )
               tracers%data(2)%values(nzmin:nzmax,node) = 3._WP - S_ref_anomaly
        end where
    end do
#ifdef ENABLE_OPENACC
    ! keep the device copies equal to the host (the monitor reads them there)
    ts_on_device = tr_dev_current(1) .and. tr_dev_current(2) .and. ptracers_restore_total==0
    if (ts_on_device) then
        salt_dev => tracers%data(2)%values
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(nz, nzmin, nzmax)
        do node=1, partit%myDim_nod2D+partit%eDim_nod2D
            nzmax=mesh%nlevels_nod2D(node)-1
            nzmin=mesh%ulevels_nod2D(node)
            do nz=nzmin, nzmax
                if (salt_dev(nz,node) > 45._WP - S_ref_anomaly) salt_dev(nz,node) = 45._WP - S_ref_anomaly
                if (salt_dev(nz,node) < 3._WP - S_ref_anomaly)  salt_dev(nz,node) = 3._WP - S_ref_anomaly
            end do
        end do
        !$ACC END PARALLEL LOOP
    end if
#endif
!$OMP END PARALLEL DO

    !---age-code-begin
    if (use_age_tracer) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nzmin, nzmax)
      do node=1,partit%myDim_nod2D+partit%eDim_nod2D
        nzmax=mesh%nlevels_nod2D(node)-1
        nzmin=mesh%ulevels_nod2D(node)
        where (tracers%data(index_age_tracer)%values(nzmin:nzmax,node) < 0._WP )
               tracers%data(index_age_tracer)%values(nzmin:nzmax,node) = 0._WP
        end where
      end do
!$OMP END PARALLEL DO
    end if
    !---age-code-end

    !___________________________________________________________________________
    ! diapycnal mixing diagnostic: count the averaging window, and accumulate the
    ! sigma2 density only on the steps where pressure_bv actually refreshed it
    ! (same dmoc_is_due predicate, which is a pure function of mstep/the clock).
    ! ts_diff2w_diap consumes and resets all three.
    if (ldiag_diapmix) then
        if (dmoc_is_due(mstep)) then
            density_dmoc_avg = density_dmoc_avg + density_dmoc
            dmoc_avg_count   = dmoc_avg_count   + 1.0_WP
        end if
        diap_avg_count   = diap_avg_count   + 1.0_WP
    end if
end subroutine solve_tracers_ale
!
!
!===============================================================================
subroutine diff_tracers_ale(tr_num, dynamics, tracers, ice, partit, mesh)

    implicit none
    integer       , intent(in)   , target :: tr_num
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_ice)   , intent(in)   , target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: n, nzmax, nzmin
    real(kind=WP)                         :: ttf_rhs_bak (mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    integer                               :: nz, nu1, nl1
#if defined(__recom)
    type(tracers_info_type)               :: tracers_info
#endif
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP_full), pointer                :: del_ttf(:,:)
    real(kind=WP_full)                         :: acc
    real(kind=WP)                             :: tnew

    del_ttf => tracers%work%del_ttf
#if defined(__recom)
    dtr_bf         = 0.0_WP
    str_bf         = 0.0_WP
    vert_sink      = 0.0_WP
#endif

    ttf_rhs_bak = 0.0
#ifdef ENABLE_OPENACC
    ! del_ttf arrives on the device from the advection; when a diagnostic, the explicit vertical
    ! scheme or recom reads it on the host, the whole diffusion of this tracer runs there
    ! single precision (use_tracer_carry): the rounding carry of the tracer update is
    ! not on the device, so the diffusion runs on the host
    diff_on_device = .not. (tracers%data(tr_num)%ltra_diag .or. ldiag_diapmix &
                            .or. (.not. tracers%data(tr_num)%i_vert_diff) .or. use_tracer_carry)
#if defined(__recom)
    diff_on_device = .false.
#endif
    if (.not. diff_on_device) then
        ! (the tracer gradients are computed on the device in init_tracers_AB)
        !$ACC UPDATE HOST(tracers%work%del_ttf, tr_xy, tr_z)
    end if
    ! the implicit vertical diffusion runs on the device for T and S unless a feature reads
    ! host-only data inside its column loop
    impl_on_device = diff_on_device .and. (tracers%data(tr_num)%ID==1 .or. tracers%data(tr_num)%ID==2) &
                     .and. .not. (use_icebergs .and. .not. turn_off_hf) .and. .not. toy_ocean &
                     .and. (.not. use_kpp_nonlclflx .or. any(mix_scheme_nmb == (/1, 18, 3, 38/)))
#endif

    if (tracers%data(tr_num)%ltra_diag) then
       call backup_ttf_rhs(del_ttf, ttf_rhs_bak, mesh%ulevels_nod2D, mesh%nlevels_nod2D, partit%myDim_nod2D, partit%eDim_nod2D)
    end if

    !___________________________________________________________________________
    ! diapycnal mixing diagnostic: remember del_ttf before any diffusive operator
    ! has touched it, so that the pure diffusion increment can be extracted below
    if (ldiag_diapmix .and. tr_num<=2) dd_diap = del_ttf

    !___________________________________________________________________________
    ! do horizontal diffusion
    ! write there also horizontal diffusion rhs to del_ttf which is equal the R_T^n
    ! in danilovs srcipt
    ! includes Redi diffusivity if Redi=.true.
    call diff_part_hor_redi(tracers, partit, mesh)  ! seems to be ~9% faster than diff_part_hor
    if (conservation_freq > 0) call conservation_probe('hor_diff', 'a', tr_num, tracers, partit, mesh)

    if (tracers%data(tr_num)%ltra_diag) then
       call store_diag_component(del_ttf, ttf_rhs_bak, mesh%hnode_new, mesh%ulevels_nod2D, mesh%nlevels_nod2D, partit%myDim_nod2D, partit%eDim_nod2D, &
                                 tracers%work%tra_diff_part_hor_redi(:,:,tr_num))
    end if

    if ((.not. tracers%data(tr_num)%i_vert_diff) .and. tracers%data(tr_num)%ltra_diag) then
       call backup_ttf_rhs(del_ttf, ttf_rhs_bak, mesh%ulevels_nod2D, mesh%nlevels_nod2D, partit%myDim_nod2D, partit%eDim_nod2D)
    end if
    !___________________________________________________________________________
    ! do vertical diffusion: explicit
    if (.not. tracers%data(tr_num)%i_vert_diff) call diff_ver_part_expl_ale(tr_num, tracers, partit, mesh)

    ! OG i_vert_diff = TRUE so, we dont call explicit scheme
    ! If we use this, check surface forcing for recom variables (They are not updated)
    if ((.not. tracers%data(tr_num)%i_vert_diff) .and. tracers%data(tr_num)%ltra_diag) then 
       call store_diag_component(del_ttf, ttf_rhs_bak, mesh%hnode_new, mesh%ulevels_nod2D, mesh%nlevels_nod2D, partit%myDim_nod2D, partit%eDim_nod2D, &
                                 tracers%work%tra_diff_part_ver_expl(:,:,tr_num))
    end if

    ! A projection of horizontal Redi diffussivity onto vertical. This par contains horizontal
    ! derivatives and has to be computed explicitly!

    if (tracers%data(tr_num)%ltra_diag .and. Redi) then
       call backup_ttf_rhs(del_ttf, ttf_rhs_bak, mesh%ulevels_nod2D, mesh%nlevels_nod2D, partit%myDim_nod2D, partit%eDim_nod2D)
    end if

    if (Redi) call diff_ver_part_redi_expl(tracers, partit, mesh)
    if (conservation_freq > 0) call conservation_probe('vert_expl', 'a', tr_num, tracers, partit, mesh)

    if (tracers%data(tr_num)%ltra_diag .and. Redi) then
       call store_diag_component(del_ttf, ttf_rhs_bak, mesh%hnode_new, mesh%ulevels_nod2D, mesh%nlevels_nod2D, partit%myDim_nod2D, partit%eDim_nod2D, &
                                 tracers%work%tra_diff_part_ver_redi_expl(:,:,tr_num))
    end if

    !___________________________________________________________________________
    ! diapycnal mixing diagnostic: del_ttf now carries the increment of all
    ! explicit diffusive operators (horizontal + Redi + explicit vertical)
    if (ldiag_diapmix .and. tr_num<=2) dd_diap = del_ttf - dd_diap

!        if (recom_debug .and. mype==0)  print *, tracers%data(tr_num)%ID

#if defined(__recom)
! 1) Remineralization from the benthos
!    Nutrient fluxes come from the bottom boundary
!    Unit [mmol/m2/s]

if (any(recom_remin_tracer_id == tracers%data(tr_num)%ID)) then

! call bottom boundary
      call diff_ver_recom_expl(mesh%nl, &  !--- vert_sink ---
                               mesh%ulevels_nod2D, mesh%nlevels_nod2D, &
                               mesh%nod_in_elem2D_num, &
                               mesh%nod_in_elem2D, mesh%nlevels, mesh%elem_area, &
                               mesh%area, mesh%areasvol, mesh%hnode_new, &
                               tracers%data(tr_num)%ID,  &
                               partit%myDim_nod2D, partit%eDim_nod2D, &
                               partit%mype, partit%MPI_COMM_FESOM, dtr_bf, dt)  !--- str_bf ---
! update tracer fields
        do n=1, partit%myDim_nod2D
            nzmax=mesh%nlevels_nod2D(n)-1
            nzmin=mesh%ulevels_nod2D(n)
!            tr_arr(nzmin:nzmax,n,tr_num)=tr_arr(nzmin:nzmax,n,tr_num)+ &
!                                                dtr_bf(nzmin:nzmax,n)
        tracers%data(tr_num)%values(nzmin:nzmax,n)=tracers%data(tr_num)%values(nzmin:nzmax,n)+ &
                                                dtr_bf(nzmin:nzmax,n)
        end do
end if

! 2) Sinking in water column
!YY: recom_sinking_tracer_id in recom_modules extended for 2. zooplankton
! but not for the combination ciso + 2. zoo!
if (any(recom_sinking_tracer_id == tracers%data(tr_num)%ID)) then

!< activate Ballasting
!< .OG. 04.11.2022

         if (use_ballasting) then
!< get seawater viscosity, seawater_visc_3D
              call get_seawater_viscosity(tr_num, partit%myDim_nod2D, &
                                          mesh%ulevels_nod2D, mesh%nlevels_nod2D, &
                                          tracers%data(1)%values, tracers%data(2)%values) ! seawater_visc_3D

              allocate(tracers_info%ids(tracers%num_tracers))
              allocate(tracers_info%data_pointers(tracers%num_tracers))

              do n = 1, tracers%num_tracers
                tracers_info%ids(n) = tracers%data(n)%id
                tracers_info%data_pointers(n)%tracer_data => tracers%data(n)%values
              end do

!< get particle density of class 1 and 2 !rho_particle1 and rho_particle2
              call get_particle_density(tracers%num_tracers, partit%myDim_nod2D, partit%eDim_nod2D, &
                                        mesh%nl, mesh%ulevels_nod2D, mesh%nlevels_nod2D, &
                                        tracers_info) ! rho_particle = density of particle class 1 and 2

              deallocate(tracers_info%ids)
              deallocate(tracers_info%data_pointers)

!< calculate scaling factors
!< scaling_visc_3D
!< scaling_density1_3D, scaling_density2_3D
              call ballast(partit%myDim_nod2D, mesh%ulevels_nod2D, &
                           mesh%nlevels_nod2D, mesh%geo_coord_nod2D, mesh%z_3d_n,                 &
                           tracers%data(1)%values, tracers%data(2)%values, rad)
        end if

! sinking
        call ver_sinking_recom(tr_num, mesh%nl, &  !--- vert_sink ---
                               mesh%ulevels_nod2D, mesh%nlevels_nod2D, &
                               mesh%zbar_3d_n, mesh%z_3d_n, mesh%nod_in_elem2D_num, &
                               mesh%nod_in_elem2D, mesh%nlevels,       &
                               mesh%area, mesh%areasvol, mesh%hnode, mesh%hnode_new, &
                               tracers%data(tr_num)%ID,  &
                               tracers%data(tr_num)%values(:,:),       &
                               partit%myDim_nod2D, vert_sink, dt)  !--- str_bf ---
! update tracer fields
! sinking into the benthos
        call ver_sinking_recom_benthos(tr_num, mesh%nl, &
                                       mesh%ulevels_nod2D, mesh%nlevels_nod2D, &
                                       mesh%zbar_3d_n, mesh%nod_in_elem2D_num, &
                                       mesh%nod_in_elem2D, mesh%nlevels, mesh%elem_area, &
                                       mesh%area, mesh%areasvol, tracers%data(tr_num)%ID, &
                                       tracers%data(tr_num)%values(:,:),       &
                                       partit%myDim_nod2D, str_bf,             &
                                       partit%mype, partit%MPI_COMM_FESOM,   &
                                       partit%npes, partit%com_nod2D%sPEnum,   &
                                       partit%com_nod2D%rPEnum,  &
                                       partit%s_mpitype_nod2D, partit%r_mpitype_nod2D,              &
                                       partit%s_mpitype_nod3D, partit%r_mpitype_nod3D,              &
                                       partit%com_nod2D%sPE, partit%com_nod2D%rPE,                  &
                                       partit%com_nod2D%req, partit%com_nod2D%nreq, dt)  !--- str_bf ---

! update tracer fields

        do n=1, partit%myDim_nod2D
            nzmax=mesh%nlevels_nod2D(n)-1
            nzmin=mesh%ulevels_nod2D(n)
!            tr_arr(nzmin:nzmax,n,tr_num)=tr_arr(nzmin:nzmax,n,tr_num)+ &
!                                                vert_sink(nzmin:nzmax,n)
!            tr_arr(nzmin:nzmax,n,tr_num)=tr_arr(nzmin:nzmax,n,tr_num)+ &
!                                                str_bf(nzmin:nzmax,n)
        tracers%data(tr_num)%values(nzmin:nzmax,n)=tracers%data(tr_num)%values(nzmin:nzmax,n)+ &
                                                vert_sink(nzmin:nzmax,n)
        tracers%data(tr_num)%values(nzmin:nzmax,n)=tracers%data(tr_num)%values(nzmin:nzmax,n)+ &
                                                str_bf(nzmin:nzmax,n)
        end do
endif
#endif
    !___________________________________________________________________________
    ! Update tracers --> calculate T* see Danilov et al. (2017)
    ! T* =  (dt*R_T^n + h^(n-0.5)*T^(n-0.5))/h^(n+0.5)
#ifdef ENABLE_OPENACC
    if (diff_on_device) then
        ! the same update on the device (no diagnostics on this path); values then goes
        ! to the host for the implicit vertical diffusion
!$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) PRIVATE(nzmin, nzmax)
        do n=1, partit%myDim_nod2D
            nzmax=mesh%nlevels_nod2D(n)-1
            nzmin=mesh%ulevels_nod2D(n)
!$ACC LOOP VECTOR
            do nz=nzmin, nzmax
                del_ttf(nz,n)=del_ttf(nz,n)+tracers%data(tr_num)%values(nz,n)*(mesh%hnode(nz,n)-mesh%hnode_new(nz,n))
                tracers%data(tr_num)%values(nz,n)=tracers%data(tr_num)%values(nz,n)+del_ttf(nz,n)/mesh%hnode_new(nz,n)
            end do
        end do
!$ACC END PARALLEL LOOP
        ! the implicit vertical diffusion reads values: on the host unless it runs on the device
        if (.not. impl_on_device) then
            !$ACC UPDATE HOST(tracers%data(tr_num)%values)
        end if
    else
#endif
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax, acc, tnew)
    do n=1, partit%myDim_nod2D
        nzmax=mesh%nlevels_nod2D(n)-1
        nzmin=mesh%ulevels_nod2D(n)
        if (use_tracer_carry) then
            do nz=nzmin, nzmax
                acc  = tracers%data(tr_num)%values(nz,n)*real(mesh%hnode(nz,n), WP_full) + del_ttf(nz,n) + tracers%work%carry(nz,n,tr_num)
                tnew = real(acc/mesh%hnode_new(nz,n), WP)
                tracers%work%carry(nz,n,tr_num) = real(acc - tnew*real(mesh%hnode_new(nz,n), WP_full), WP)
                del_ttf(nz,n) = del_ttf(nz,n) + tracers%data(tr_num)%values(nz,n)*real(mesh%hnode(nz,n)-mesh%hnode_new(nz,n), WP_full)
                tracers%data(tr_num)%values(nz,n) = tnew
            end do
        else
        del_ttf(nzmin:nzmax,n)=del_ttf(nzmin:nzmax,n)+tracers%data(tr_num)%values(nzmin:nzmax,n)* &
                                    (mesh%hnode(nzmin:nzmax,n)-mesh%hnode_new(nzmin:nzmax,n))
        tracers%data(tr_num)%values(nzmin:nzmax,n)=tracers%data(tr_num)%values(nzmin:nzmax,n)+ &
                                    del_ttf(nzmin:nzmax,n)/mesh%hnode_new(nzmin:nzmax,n)
        end if
        ! WHY NOT ??? --> whats advantage of above --> tested it --> the upper
        ! equation has a 30% smaller nummerical drift
        ! tr_arr(1:nzmax,n,tr_num)=(hnode(1:nzmax,n)*tr_arr(1:nzmax,n,tr_num)+ &
        !                           del_ttf(1:nzmax,n))/hnode_new(1:nzmax,n)
        ! diapycnal mixing diagnostic: same normalisation, turns the explicit
        ! diffusion increment into a tracer change
        if (ldiag_diapmix .and. tr_num<=2) &
            dd_diap(nzmin:nzmax,n)=dd_diap(nzmin:nzmax,n)/mesh%hnode_new(nzmin:nzmax,n)
    end do
!$OMP END PARALLEL DO
#ifdef ENABLE_OPENACC
    end if
#endif

    !___________________________________________________________________________
    if (tracers%data(tr_num)%i_vert_diff) then
        ! do vertical diffusion: implicite

        if (tracers%data(tr_num)%ltra_diag) then
           do n=1, partit%myDim_nod2D+partit%eDim_nod2D
              nu1 = mesh%ulevels_nod2D(n)
              nl1 = mesh%nlevels_nod2D(n)
              do nz = nu1, nl1-1
                 ttf_rhs_bak(nz,n) = tracers%data(tr_num)%values(nz,n)
              end do
           end do
        end if

        ! diapycnal mixing diagnostic: subtract the tracer field before the
        ! implicit solve; adding it back afterwards leaves the implicit vertical
        ! diffusion increment on top of the explicit one already in dd_diap
        if (ldiag_diapmix .and. tr_num<=2) &
            dd_diap(:,1:partit%myDim_nod2D) = dd_diap(:,1:partit%myDim_nod2D) - tracers%data(tr_num)%values(:,1:partit%myDim_nod2D)

        ! (w/out Redi)
        if (conservation_freq > 0) call conservation_probe('divided', 'b', tr_num, tracers, partit, mesh)
        call diff_ver_part_impl_ale(tr_num, dynamics, tracers, ice, partit, mesh)
        if (impl_on_device) then
            !$ACC UPDATE HOST(tracers%data(tr_num)%values)
        end if
        if (conservation_freq > 0) call conservation_probe('vert_impl+bc', 'b', tr_num, tracers, partit, mesh)

        if (ldiag_diapmix .and. tr_num<=2) &
            dd_diap(:,1:partit%myDim_nod2D) = dd_diap(:,1:partit%myDim_nod2D) + tracers%data(tr_num)%values(:,1:partit%myDim_nod2D)

        ! vertical diffusion: implicit
        if (tracers%data(tr_num)%ltra_diag) then
           do n=1, partit%myDim_nod2D+partit%eDim_nod2D
              nu1 = mesh%ulevels_nod2D(n)
              nl1 = mesh%nlevels_nod2D(n)
              do nz = nu1, nl1-1
                 tracers%work%tra_diff_part_ver_impl(nz,n,tr_num) = tracers%data(tr_num)%values(nz,n) - ttf_rhs_bak(nz,n)
                 !if (mype==0)  print *,  tra_diff_part_ver_impl(:,:,tr_num)
              end do
           end do
        end if

    end if

    !___________________________________________________________________________
    ! diapycnal mixing diagnostic: dd_diap now holds the total change of this
    ! tracer over this time step that is due to diffusion. Accumulate it over the
    ! averaging window of ts_diff2w_diap, which consumes and resets it.
    if (ldiag_diapmix) then
        if (tr_num==1) dT_diap(:,1:partit%myDim_nod2D) = dT_diap(:,1:partit%myDim_nod2D) + dd_diap(:,1:partit%myDim_nod2D)
        if (tr_num==2) dS_diap(:,1:partit%myDim_nod2D) = dS_diap(:,1:partit%myDim_nod2D) + dd_diap(:,1:partit%myDim_nod2D)
    end if

    !We DO not set del_ttf to zero because it will not be used in this timestep anymore
    !init_tracers_AB will set it to zero for the next timestep
    if (tracers%data(tr_num)%smooth_bh_tra) then
       call diff_part_bh(tr_num, dynamics, tracers, partit, mesh)  ! alpply biharmonic diffusion (implemented as filter)
    end if
#ifdef ENABLE_OPENACC
    if (tracers%data(tr_num)%ID==1 .or. tracers%data(tr_num)%ID==2) &
        tr_dev_current(tracers%data(tr_num)%ID) = impl_on_device .and. tracers%data(tr_num)%i_vert_diff &
                                                  .and. .not. tracers%data(tr_num)%smooth_bh_tra
#endif
    if (conservation_freq > 0) call conservation_probe('bh_filter', 'b', tr_num, tracers, partit, mesh)

contains

subroutine backup_ttf_rhs(del_ttf, ttf_rhs_bak, ulevels_nod2D, nlevels_nod2D, myDim_nod2D, eDim_nod2D)
    real(kind=WP_full), intent(in)    :: del_ttf(:,:)
    real(kind=WP), intent(inout) :: ttf_rhs_bak(:,:)
    integer      , intent(in)    :: ulevels_nod2D(:)
    integer      , intent(in)    :: nlevels_nod2D(:)
    integer      , intent(in)    :: myDim_nod2D
    integer      , intent(in)    :: eDim_nod2D

    integer :: n, nz, nu1, nl1

    do n=1, myDim_nod2D+eDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        do nz = nu1, nl1-1
            ttf_rhs_bak(nz,n) = del_ttf(nz,n)
        end do
    end do
end subroutine backup_ttf_rhs

subroutine store_diag_component(del_ttf, ttf_rhs_bak, hnode_new, ulevels_nod2D, nlevels_nod2D, myDim_nod2D, eDim_nod2D, target_field)
    real(kind=WP_full), intent(in)    :: del_ttf(:,:)
    real(kind=WP), intent(in)    :: ttf_rhs_bak(:,:)
    real(kind=WP_full), intent(in)  :: hnode_new(:,:)
    integer      , intent(in)    :: ulevels_nod2D(:)
    integer      , intent(in)    :: nlevels_nod2D(:)
    integer      , intent(in)    :: myDim_nod2D
    integer      , intent(in)    :: eDim_nod2D
    real(kind=WP), intent(inout) :: target_field(:,:)

    integer :: n, nz, nu1, nl1

    do n=1, myDim_nod2D+eDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        do nz = nu1, nl1-1
            target_field(nz,n) = (del_ttf(nz,n) - ttf_rhs_bak(nz,n)) / hnode_new(nz,n) ! Unit [Conc]
        end do
    end do
end subroutine store_diag_component

end subroutine diff_tracers_ale
!
!
!===============================================================================
!Vertical diffusive flux(explicit scheme):
subroutine diff_ver_part_expl_ale(tr_num, tracers, partit, mesh)
    implicit none
    integer       , intent(in)   , target :: tr_num
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: n, nz, nl1, ul1
    real(kind=WP)            :: vd_flux(mesh%nl-1)
    real(kind=WP)            :: rdata, flux, rlx
    real(kind=WP)            :: zinv1
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP_full), pointer   :: del_ttf(:,:)
    del_ttf => tracers%work%del_ttf

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nl1, ul1, vd_flux, rdata, flux, rlx, zinv1)
    !___________________________________________________________________________
    do n=1, partit%myDim_nod2D
        nl1=mesh%nlevels_nod2D(n)-1
        ul1=mesh%ulevels_nod2D(n)

        vd_flux=0._WP
        if (tracers%data(tr_num)%ID==1) then
            flux  = -heat_flux(n)/vcpw
            rdata =  Tsurf(n)
            rlx   =  surf_relax_T
        elseif (tracers%data(tr_num)%ID==2) then
            ! use_salt_anomaly background-dilution term (S_ref=0 unless enabled).
            ! See bc_surface (implicit path) below for the derivation and the note
            ! on the ~8e-6/step constancy-error residual it leaves.
            flux  =  virtual_salt(n)+relax_salt(n) + (real_salt_flux(n) + S_ref_anomaly*water_flux(n))*is_nonlinfs
        else
            flux  = 0._WP
            rdata = 0._WP
            rlx=0._WP
        endif
        !_______________________________________________________________________
        !Surface forcing
        vd_flux(ul1)= flux
        do nz=ul1+1,nl1
            !___________________________________________________________________
            zinv1=1.0_WP/(mesh%Z_3d_n(nz-1,n)-mesh%Z_3d_n(nz,n))
            vd_flux(nz) = Kv(nz,n)*(tracers%data(tr_num)%values(nz-1,n)-tracers%data(tr_num)%values(nz,n))*zinv1*mesh%area(nz,n)
        end do

        do nz=ul1,nl1-1
            del_ttf(nz,n) = del_ttf(nz,n) + (vd_flux(nz) - vd_flux(nz+1))/(mesh%zbar_3d_n(nz,n)-mesh%zbar_3d_n(nz+1,n))*dt/mesh%areasvol(nz,n)
        end do
        del_ttf(nl1,n) = del_ttf(nl1,n) + (vd_flux(nl1)/(mesh%zbar_3d_n(nl1,n)-mesh%zbar_3d_n(nl1+1,n)))*dt/mesh%areasvol(nl1,n)

    end do ! --> do n=1, myDim_nod2D
!$OMP END PARALLEL DO
end subroutine diff_ver_part_expl_ale
!
!
!===============================================================================
! vertical diffusivity augmented with Redi contribution [vertical flux of K(3,3)*d_zT]
subroutine diff_ver_part_impl_ale(tr_num, dynamics, tracers, ice, partit, mesh)
    implicit none
    integer       , intent(in)   , target :: tr_num
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh

    type(t_ice)   , intent(in)   , target :: ice

    !___________________________________________________________________________
    real(kind=WP_full)        :: a(impl_nl_max), b(impl_nl_max), c(impl_nl_max), tr(impl_nl_max)
    real(kind=WP_full)        :: cp(impl_nl_max), tp(impl_nl_max)
    integer                  :: nz, n, nzmax, nzmin
    real(kind=WP_full)        :: m, zinv, dz
    real(kind=WP_full)        :: rsss, Ty, Ty1, c1, zinv1, zinv2, v_adv
    real(kind=WP_full)        :: isredi=0._WP_full
    logical                  :: do_wimpl=.true.
    real(kind=WP_full)        :: zbar_n(impl_nl_max), z_n(impl_nl_max)
    real(kind=WP_full)        :: acc
    real(kind=WP)            :: tnew
    integer                  :: tr_id
    logical                  :: sw_dbgyre, hosing_depth
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:), pointer :: trarr
    real(kind=WP), dimension(:,:), pointer :: Wvel_i
    real(kind=WP), dimension(:,:), pointer :: sst, sss ! auxiliary variables needed for transient tracers
    real(kind=WP), dimension(:),   pointer :: a_ice
    trarr  => tracers%data(tr_num)%values(:,:)
    Wvel_i => dynamics%w_i(:,:)

    sst => tracers%data(1)%values(:,:)
    sss => tracers%data(2)%values(:,:)
    a_ice => ice%data(1)%values(:)
    !___________________________________________________________________________
    if ((trim(tracers%data(tr_num)%tra_adv_lim)=='FCT') .OR. (.not. dynamics%use_wsplit)) do_wimpl=.false.
    if (Redi) isredi=1._WP
    Ty    =0.0_WP
    Ty1   =0.0_WP
    if (mesh%nl > impl_nl_max) then
        if (partit%mype==0) write(*,*) 'diff_ver_part_impl_ale: mesh%nl =', mesh%nl, ' > impl_nl_max =', impl_nl_max
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
    end if
    tr_id        = tracers%data(tr_num)%ID
    sw_dbgyre    = use_sw_pene .and. tr_id==1 .and. toy_ocean .and. TRIM(which_toy)=="dbgyre"
    hosing_depth = tr_id==304 .and. use_hosing .and. trim(hosing_mode)=='depth'
    ! fills sw_3d for all nodes once, before the column loop
    if (sw_dbgyre) call cal_shortwave_rad_dbgyre(ice, tracers, partit, mesh)

    !___________________________________________________________________________
    ! surface boundary term (heat, freshwater, salt, relaxation fluxes; transient tracers),
    ! computed per node before the column loop, which only adds it
    if (.not. allocated(impl_bc)) then
        allocate(impl_bc(partit%myDim_nod2D))
        !$ACC ENTER DATA CREATE(impl_bc)
    end if
    if (impl_on_device) then
        ! T and S only (see diff_tracers_ale)
!$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(nzmin)
        do n=1, partit%myDim_nod2D
            nzmin=mesh%ulevels_nod2D(n)
            ! the CASE (1) and CASE (2) terms of bc_surface
            if (tr_id==1) then
                impl_bc(n)=-dt*(heat_flux(n)/vcpw + trarr(nzmin,n)*water_flux(n)*is_nonlinfs)
            else
                impl_bc(n)= dt*(virtual_salt(n) &
                            + relax_salt(n) + (real_salt_flux(n) + S_ref_anomaly*water_flux(n))*is_nonlinfs)
            end if
        end do
!$ACC END PARALLEL LOOP
    elseif (hosing_depth) then
        ! the hosing passive tracer (ID 304) receives its source over the whole column
        if (.not. allocated(impl_bcz)) allocate(impl_bcz(mesh%nl, partit%myDim_nod2D))
        do n=1, partit%myDim_nod2D
            nzmax=mesh%nlevels_nod2D(n)
            nzmin=mesh%ulevels_nod2D(n)
            do nz=nzmin,nzmax
                impl_bcz(nz,n)=bc_surface(n, tr_id, trarr(nz,n), nz, partit, mesh, sst(nz,n), sss(nz,n), a_ice(n))
            end do
        end do
    else
        do n=1, partit%myDim_nod2D
            nzmin=mesh%ulevels_nod2D(n)
            impl_bc(n)=bc_surface(n, tr_id, trarr(nzmin,n), nzmin, partit, mesh, sst(nzmin,n), sss(nzmin,n), a_ice(n))
        end do
    end if
    
    ! solve equation diffusion equation implicite part:
    ! -->   h^(n+0.5)* (T^(n+0.5)-Tstar) = dt*( K_33*d/dz*(T^(n+0.5)-Tstar) + K_33*d/dz*Tstar )
    ! -->   Tnew = T^(n+0.5)-Tstar
    ! -->   h^(n+0.5)* (Tnew) = dt*(K_33*d/dz*Tnew) + K_33*dt*d/dz*Tstar
    ! -->   h^(n+0.5)* (Tnew) = dt*(K_33*d/dz*Tnew) + RHS
    ! -->   solve for T_new
    ! -->   As_1 (Skalar Area), A_1 (Area of edge),               
    !       no Cavity A_1==As_1, yes Cavity A1 !=As_1             
    !                                                             
    !    ----------- zbar_1, A_1                           nvec_up (+1)  
    ! Z_1 o T_1, As_1                                       ^ 
    !    ----------- zbar_2, A_2                       _____|_____    Gaus Theorem:
    ! Z_2 o T_2, As_2                                 /|    |    |\    --> Flux form
    !    ----------- zbar_3, A_3                     * |    |    | *   --> normal vec outwards facing
    ! Z_3 o T_3, As_3                                |\___________/|        
    !    ----------- zbar_4                          | |         | |         
    !        :                                       | |_________| |          
    !                                                |/           \|           
    !                                                *      |      *            
    !                                                 \_____|_____/
    !                                                       |
    !                                                       V nvec_dwn (-1)            
    !                                                           
    ! --> 1st. solve homogenouse part:
    !     f(Tnew) = h^(n+0.5)* (Tnew) - dt*(K_33*dTnew/dz) = 0
    !
    ! --> 2nd. Difference Quotient at Tnew_i in flux form (Gaus Theorem, dont forget normal vectors!!!):
    !        |
    !        +-> Gauss THeorem: int(V', div(F_vec))dV = intcircle(A', F_vec*n_vec)dA
    !        |
    !        +-> du                     = dt*( K_33*d/dz*dTnew )  | *div()_z=d/dz, *int(V',)dV
    !        |   int(V', d/dz *du)dV    = int(V', d/dz *dt*( K_33*d/dz*dTnew ) )dV
    !        |   ...                    = intcircle(A', dt*( K_33*d/dz*dTnew )*n_vec)dA   
    !        |   int(V', d/dz *du)dz*dA = ...
    !        |   int(A', du)dA          = intcircle(A', dt*( K_33*d/dz*dTnew )*n_vec)dA
    !        |   
    !        +-> As_i area of scalar cell = A_i, A_i+1 area of top 
    !        |   and bottom face of scalar cell 
    !        V
    !
    !    As_i*Tnew_i *h_i = dt * [ K_33 * (Tnew_i-1 - Tnew_i  )/(Z_i-1 - Z_i  ) * A_i   * nvec_up(+1)
    !                              +K_33 * (Tnew_i   - Tnew_i+1)/(Z_i   - Z_i+1) * A_i+1 * nvec_dwn(-1) ]
    !    Tnew_i *h_i      = dt * [ K_33 * (Tnew_i-1 - Tnew_i  )/(Z_i-1 - Z_i  ) * A_i  /As_i * nvec_up(+1)
    !                              +K_33 * (Tnew_i   - Tnew_i+1)/(Z_i   - Z_i+1) * A_i+1/As_i * nvec_dwn(-1) ]
    !        |
    !        +-> since we are on scalar cell As_i != A_i != A_i+1 
    !        +-> take into account normal vector direction
    !        V
    !    f(Tnew) = Tnew_i*h_i - dt*K_33 * (Tnew_i-1 - Tnew_i  )/(Z_i-1 - Z_i  ) * A_i  /As_i 
    !                         + dt*K_33 * (Tnew_i   - Tnew_i+1)/(Z_i   - Z_i+1) * A_i+1/As_i
    !            = 0
    !
    ! --> 3rd. solve for coefficents a, b, c:
    !     f(Tnew) = [ a*dTnew_i-1 + b*dTnew_i + c*dTnew_i+1 ]
    !        |
    !        +-> estimate a, b, c by derivation of f(Tnew_i)
    !        |
    !        +-> a = df(Tnew)/dTnew_i-1 = -dt*K_33/(Z_i-1 - Z_i) * A_i/As_i
    !        |
    !        +-> c = df(Tnew)/dTnew_i+1 = -dt * K_33 * 1/(Z_i - Z_i+1) * A_i+1/As_i
    !        |
    !        +-> b = df(Tnew)/dTnew_i   = h_i + dt*K_33/(Z_i-1 - Z_i) * A_i  /As_i
    !                                          + dt*K_33/(Z_i - Z_i+1) * A_i+1/As_i
    !                                   = h_i -(a+c)
    !
    ! --> 4th. solve inhomogenous part:
    !     [ a*dTnew_i-1 + b*dTnew_i + c*dTnew_i+1 ] = RHS/As_i
    !
    !     RHS     = K_33*dt*d/dz*Tstar
    !
    ! --> write as Difference Quotient in flux form
    !     RHS/As_i =  K_33 * dt * (Tstar_i-1 - Tstar_i)/(Z_i-1 - Z_i) * A_i/As_i   * (nvec_up=1)
    !               + K_33 * dt * (Tstar_i - Tstar_i+1)/(Z_i - Z_i+1) * A_i+1/As_i * (nvec_dwn=-1)
    !
    !              =  K_33*dt/(Z_i-1 - Z_i) * A_i/As_i   * Tstar_i-1
    !               - K_33*dt/(Z_i-1 - Z_i) * A_i/As_i   * Tstar_i
    !               - K_33*dt/(Z_i - Z_i+1) * A_i+1/As_i * Tstar_i
    !               + K_33*dt/(Z_i - Z_i+1) * A_i+1/As_i * Tstar_i+1
    !
    !              = -a*Tstar_i-1 + (a+c)*Tstar_i - c * Tstar_i+1
    !              |-> b = h_i - (a+c), a+c = h_i-b
    !
    !              = -a*Tstar_i-1 - (b-h_i)*Tstar_i - c * Tstar_i+1
    !
    ! --> 5th. solve for Tnew_i --> forward sweep algorithm --> see lower
    !  | b_1 c_1 ...            |   |dTnew_1|
    !  | a_2 b_2 c_2 ...        |   |dTnew_2|
    !  |     a_3 b_3 c_3 ...    | * |dTnew_3| = RHS/V_i
    !  |         a_4 b_4 c_4 ...|   |dTnew_4|
    !  |              :         |   |   :   |
    !
    ! --> a = -dt*K_33 / (Z_i-1 - Z_i) * A_i/As_i
    !
    ! --> c = -dt*K_33 / (Z_i - Z_i+1) * A_i+1/As_i
    !
    ! --> b = h^(n+0.5) -[ dt*K_33/(Z_i-1 - Z_i  ) * A_i  /As_i 
    !                     +dt*K_33/(Z_i   - Z_i+1) * A_i+1/As_i ] 
    !       = -(a+c) + h^(n+0.5)
    !
    !___________________________________________________________________________
    ! loop over local nodes

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, nzmax, nzmin, a, b, c, tr, cp, tp, m, zinv, dz, &
!$OMP                                         rsss, Ty, Ty1, c1, zinv1, zinv2, v_adv, zbar_n, z_n, acc, tnew)
!$OMP DO
#else
!$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(impl_on_device) &
!$ACC PRIVATE(nz, nzmax, nzmin, a, b, c, tr, cp, tp, m, zinv, dz, rsss, Ty, Ty1, c1, zinv1, zinv2, v_adv, zbar_n, z_n, acc, tnew)
#endif
    do n=1,partit%myDim_nod2D
        ! initialise
        a  = 0.0_WP
        b  = 0.0_WP
        c  = 0.0_WP
        tr = 0.0_WP
        tp = 0.0_WP
        cp = 0.0_WP
        ! max. number of levels at node n
        nzmax=mesh%nlevels_nod2D(n)
        nzmin=mesh%ulevels_nod2D(n)
        !___________________________________________________________________________
        ! Here can not exchange zbar_n & Z_n with zbar_3d_n & Z_3d_n because
        ! they be calculate from the actualized mesh with hnode_new
        ! calculate new zbar (depth of layers) and Z (mid depths of layers)
        ! depending on layer thinkness over depth at node n
        ! Be carefull here vertical operation have to be done on NEW vertical mesh !!!
        zbar_n=0.0_WP
        Z_n=0.0_WP
        zbar_n(nzmax)=mesh%zbar_n_bot(n)
        Z_n(nzmax-1)=zbar_n(nzmax) + mesh%hnode_new(nzmax-1,n)/2.0_WP
        do nz=nzmax-1,nzmin+1,-1
            zbar_n(nz) = zbar_n(nz+1) + mesh%hnode_new(nz,n)
            Z_n(nz-1)  = zbar_n(nz) + mesh%hnode_new(nz-1,n)/2.0_WP
        end do
        zbar_n(nzmin) = zbar_n(nzmin+1) + mesh%hnode_new(nzmin,n)
        !_______________________________________________________________________
        ! Regular part of coefficients: --> surface layer
        nz=nzmin

        ! 1/dz(nz)
        zinv2=1.0_WP/(Z_n(nz)-Z_n(nz+1))
        zinv=1.0_WP*dt    ! no .../(zbar(1)-zbar(2)) because of  ALE

        ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
        Ty1= (zbar_n(nz+1)-Z_n(   nz+1))*zinv2 *slope_tapered(3,nz  ,n)**2*Ki(nz  ,n) + &
             (Z_n(nz)     -zbar_n(nz+1))*zinv2 *slope_tapered(3,nz+1,n)**2*Ki(nz+1,n)
        Ty1=Ty1*isredi

        ! layer dependent coefficients for for solving dT(1)/dt+d/dz*K_33*d/dz*T(1) = ...
        a(nz)=0.0_WP
        !!PS c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv * (area(nz+1,n)/areasvol(nz,n))
        c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv * mesh%area(nz+1,n)/mesh%areasvol(nz,n)
        b(nz)=-c(nz)+mesh%hnode_new(nz,n)      ! ale
        ! update from the vertical advection --> comes from splitting of vert
        ! velocity into explicite and implicite contribution
        if (do_wimpl) then
            !___________________________________________________________________
            ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for
            ! numerical reasons, to gurante that area/areasvol in case of no
            ! cavity is ==1.0_WP
            v_adv =zinv * ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
            b(nz) =b(nz)+Wvel_i(nz, n)*v_adv

            !!PS v_adv =zinv * ( area(nz+1,n)/areasvol(nz,n) )
            v_adv =zinv * mesh%area(nz+1,n)/mesh%areasvol(nz,n)
            b(nz) =b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
            c(nz) =c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
        end if
        ! backup zinv2 for next depth level
        zinv1=zinv2
        !_______________________________________________________________________
        ! Regular part of coefficients: --> 2nd...nl-2 layer
        do nz=nzmin+1, nzmax-2

            ! 1/dz(nz)
            zinv2=1.0_WP/(Z_n(nz)-Z_n(nz+1))
            ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
            Ty = (zbar_n(nz  )-Z_n(nz     ))*zinv1 *slope_tapered(3,nz-1,n)**2*Ki(nz-1,n)+ &
                 (Z_n(nz-1   )-zbar_n(nz  ))*zinv1 *slope_tapered(3,nz  ,n)**2*Ki(nz  ,n)
            Ty1= (zbar_n(nz+1)-Z_n(nz+1   ))*zinv2 *slope_tapered(3,nz  ,n)**2*Ki(nz  ,n)+ &
                 (Z_n(nz     )-zbar_n(nz+1))*zinv2 *slope_tapered(3,nz+1,n)**2*Ki(nz+1,n)
            Ty =Ty *isredi
            Ty1=Ty1*isredi

            ! layer dependent coefficients for for solving dT(nz)/dt+d/dz*K_33*d/dz*T(nz) = ...
            !___________________________________________________________________
            ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for
            ! numerical reasons, to gurante that area/areasvol in case of no
            ! cavity is ==1.0_WP
            a(nz)=-(Kv(nz,n)  +Ty )*zinv1*zinv * ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
            !!PS c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv * ( area(nz+1,n)/areasvol(nz,n) )
            c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv * mesh%area(nz+1,n)/mesh%areasvol(nz,n)
            b(nz)=-a(nz)-c(nz)+mesh%hnode_new(nz,n)

            ! backup zinv2 for next depth level
            zinv1=zinv2
            ! update from the vertical advection
            if (do_wimpl) then
                !_______________________________________________________________
                ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for
                ! numerical reasons, to gurante that area/areasvol in case of no
                ! cavity is ==1.0_WP
                v_adv=zinv * ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
                a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv
                b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv
                !!PS v_adv=v_adv*areasvol(nz+1,n)/areasvol(nz,n)
                !!PS v_adv=zinv * ( area(nz+1,n)/areasvol(nz,n) )
                v_adv=zinv * mesh%area(nz+1,n)/mesh%areasvol(nz,n)
                b(nz)=b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
                c(nz)=c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
            end if
        end do ! --> do nz=nzmin+1, nzmax-2

        !_______________________________________________________________________
        ! Regular part of coefficients: --> nl-1 layer
        nz=nzmax-1

        zinv=1.0_WP*dt   ! no ... /(zbar(nzmax-1)-zbar(nzmax)) because of ale

        ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
        Ty= (zbar_n(nz)-Z_n(nz)   ) * zinv1 * slope_tapered(3,nz-1,n)**2 * Ki(nz-1,n) + &
            (Z_n(nz-1) -zbar_n(nz)) * zinv1 * slope_tapered(3,nz  ,n)**2 * Ki(nz,n)
        Ty =Ty *isredi
        ! layer dependent coefficients for for solving dT(nz)/dt+d/dz*K_33*d/dz*T(nz) = ...

        !___________________________________________________________________
        ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for
        ! numerical reasons, to gurante that area/areasvol in case of no
        ! cavity is ==1.0_WP
        a(nz)=-(Kv(nz,n)+Ty)*zinv1*zinv* ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
        c(nz)=0.0_WP
        b(nz)=-a(nz)+mesh%hnode_new(nz,n)

        ! update from the vertical advection
        if (do_wimpl) then
            !___________________________________________________________________
            ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for
            ! numerical reasons, to gurante that area/areasvol in case of no
            ! cavity is ==1.0_WP
            v_adv=zinv* ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
            a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv
            b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv
        end if
        
        !_______________________________________________________________________
        ! the rhs (inhomogene part): --> rhs = K_33*dt*d/dz*Tstar --> Tstar...trarr
        ! solve difference quotient for rhs --> tr
        !  RHS at Volume_2:
        !
        !  RHS*V_2 = K_33*dt*(T_1-T_2)/(Z_1-Z_2)*V_2 - K_33*dt*(T_2-T_3)/(Z_2-Z_3)*V_3
        !          = -a*T_1 + (a+c)*T_2 - c*T_3
        !
        ! -+--> tr(1) =(a(1)+c(1))*trarr(1,n)-c(1)*trarr(2,n)
        !  |--> a(1)=0
        nz=nzmin
        dz=mesh%hnode_new(nz,n) ! It would be (zbar(nz)-zbar(nz+1)) if not ALE
        tr(nz)=-(b(nz)-dz)*trarr(nz,n)-c(nz)*trarr(nz+1,n)

        do nz=nzmin+1,nzmax-2
            dz=mesh%hnode_new(nz,n)
            tr(nz)= -a(nz)     * trarr(nz-1,n) &
                    -(b(nz)-dz)* trarr(nz,n) &
                    -c(nz)     * trarr(nz+1,n)
        end do

        nz=nzmax-1
        dz=mesh%hnode_new(nz,n)
        tr(nz)=-a(nz)*trarr(nz-1,n)-(b(nz)-dz)*trarr(nz,n)
        
        !_______________________________________________________________________
        ! Add KPP nonlocal fluxes to the rhs (only T and S currently)
        ! use here blmc or kpp_oblmixc instead of Kv, since Kv already contains
        ! at this point the mixing enhancments from momix, instable
        ! mixing or windmixing which are to much for nonlocal
        ! transports and lead to instability of the model
        if (use_kpp_nonlclflx) then
            if (tr_id==2) then
                rsss=ref_sss
                if (ref_sss_local) rsss=trarr(1,n)
            end if

            !___________________________________________________________________
            ! use fesom1.4 KPP
            if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==18) then
                if     (tr_id==1) then ! temperature
                    ! --> no fluxes to the top out of the surface, no fluxes
                    !     downwards out of the bottom
                    !___surface_________________________________________________
                    nz = nzmin
                    tr(nz)=tr(nz) &
                               +(-MIN(ghats(nz+1,n)*blmc(nz+1,n,2), 1.0_WP)*(mesh%area(nz+1,n)/mesh%areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                    !___bulk____________________________________________________
                    do nz=nzmin+1, nzmax-2
                        tr(nz)=tr(nz) &
                               +( MIN(ghats(nz  ,n)*blmc(nz  ,n,2), 1.0_WP)*(mesh%area(nz  ,n)/mesh%areasvol(nz,n)) &
                                 -MIN(ghats(nz+1,n)*blmc(nz+1,n,2), 1.0_WP)*(mesh%area(nz+1,n)/mesh%areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                    end do
                    !___bottom__________________________________________________
                    nz = nzmax-1
                    tr(nz)=tr(nz) &
                               +( MIN(ghats(nz  ,n)*blmc(nz  ,n,2), 1.0_WP)*(mesh%area(nz  ,n)/mesh%areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt

                elseif (tr_id==2) then ! salinity
                    ! --> no fluxes to the top out of the surface, no fluxes
                    !     downwards out of the bottom
                    !___surface_________________________________________________
                    nz = nzmin
                    tr(nz)=tr(nz) &
                               -(-MIN(ghats(nz+1,n)*blmc(nz+1,n,3), 1.0_WP)*(mesh%area(nz+1,n)/mesh%areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                    !___bulk____________________________________________________
                    do nz=nzmin+1, nzmax-2
                        tr(nz)=tr(nz) &
                               -( MIN(ghats(nz  ,n)*blmc(nz  ,n,3), 1.0_WP)*(mesh%area(nz  ,n)/mesh%areasvol(nz,n)) &
                                 -MIN(ghats(nz+1,n)*blmc(nz+1,n,3), 1.0_WP)*(mesh%area(nz+1,n)/mesh%areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                    end do
                    !___bottom__________________________________________________
                    nz = nzmax-1
                    tr(nz)=tr(nz) &
                               -( MIN(ghats(nz  ,n)*blmc(nz  ,n,3), 1.0_WP)*(mesh%area(nz  ,n)/mesh%areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                end if
#if defined (__cvmix)                   
            !___________________________________________________________________
            ! use cvmix KPP
            elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==38) then
                if     (tr_id==1) then ! temperature
                    !___surface_________________________________________________
                    nz = nzmin
                    tr(nz)=tr(nz) &
                               +(-MIN(kpp_nonlcltranspT(nz+1,n)*kpp_oblmixc(nz+1,n,2), 1.0_WP)*(mesh%area(nz+1,n)/mesh%areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                    !___bulk____________________________________________________
                    do nz=nzmin+1, nzmax-2
                        tr(nz)=tr(nz) &
                               +( MIN(kpp_nonlcltranspT(nz  ,n)*kpp_oblmixc(nz  ,n,2), 1.0_WP)*(mesh%area(nz  ,n)/mesh%areasvol(nz,n)) &
                                 -MIN(kpp_nonlcltranspT(nz+1,n)*kpp_oblmixc(nz+1,n,2), 1.0_WP)*(mesh%area(nz+1,n)/mesh%areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                    end do
                    !___bottom__________________________________________________
                    nz = nzmax-1
                    tr(nz)=tr(nz) &
                               +( MIN(kpp_nonlcltranspT(nz  ,n)*kpp_oblmixc(nz  ,n,2), 1.0_WP)*(mesh%area(nz  ,n)/mesh%areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt

                elseif (tr_id==2) then ! salinity
                    !___surface_________________________________________________
                    nz = nzmin
                    tr(nz)=tr(nz) &
                               -(-MIN(kpp_nonlcltranspS(nz+1,n)*kpp_oblmixc(nz+1,n,3), 1.0_WP)*(mesh%area(nz+1,n)/mesh%areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                    !___bulk____________________________________________________
                    do nz=nzmin+1, nzmax-2
                        tr(nz)=tr(nz) &
                               -( MIN(kpp_nonlcltranspS(nz  ,n)*kpp_oblmixc(nz  ,n,3), 1.0_WP)*(mesh%area(nz  ,n)/mesh%areasvol(nz,n)) &
                                 -MIN(kpp_nonlcltranspS(nz+1,n)*kpp_oblmixc(nz+1,n,3), 1.0_WP)*(mesh%area(nz+1,n)/mesh%areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                    end do
                    !___bottom__________________________________________________
                    nz = nzmax-1
                    tr(nz)=tr(nz) &
                               -( MIN(kpp_nonlcltranspS(nz  ,n)*kpp_oblmixc(nz  ,n,3), 1.0_WP)*(mesh%area(nz  ,n)/mesh%areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                end if
#endif                
            end if
        end if ! --> if (use_kpp_nonlclflx) then
        
        !_______________________________________________________________________
        ! case of activated shortwave penetration into the ocean, ad 3d contribution
        if (use_sw_pene .and. tr_id==1 .and. .not. toy_ocean) then
            do nz=nzmin, nzmax-1
                zinv=1.0_WP*dt  !/(zbar(nz)-zbar(nz+1)) ale!
                !!PS tr(nz)=tr(nz)+(sw_3d(nz, n)-sw_3d(nz+1, n) * ( area(nz+1,n)/areasvol(nz,n)) ) * zinv
                tr(nz)=tr(nz)+(sw_3d(nz, n)-sw_3d(nz+1, n) * mesh%area(nz+1,n)/mesh%areasvol(nz,n)) * zinv
            end do
        elseif (sw_dbgyre) then
            do nz=nzmin, nzmax-1
                zinv=1.0_WP*dt  !/(zbar(nz)-zbar(nz+1)) ale!
                tr(nz)=tr(nz)+(sw_3d(nz, n)-sw_3d(nz+1, n)*mesh%area(nz+1,n)/mesh%area(nz,n))*zinv
            end do

        end if
        
        !_______________________________________________________________________
        ! case of activated shortwave penetration into the ocean, ad 3d contribution
        if (use_icebergs .and. (.not. turn_off_hf) .and. tr_id==1) then
            do nz=nzmin, nzmax-1
                zinv=1.0_WP*dt  !/(zbar(nz)-zbar(nz+1)) ale!
                tr(nz)=tr(nz) + ibhf_n(nz, n) * zinv / vcpw * (mesh%area(nz,n)/mesh%areasvol(nz,n))
            end do
        end if
        
        !_______________________________________________________________________
        !  The first row contains also the boundary condition from heatflux,
        !  freshwaterflux and relaxation terms
        !  --> trarr(1,n)*water_flux(n) : latent heatflux contribution due to
        !      cell volume. If Volume decreases --> temp has to raise, if volume
        !      expended --> temp has to decrease
        !                           (-)   ^                        (-)   ^
        !                            |    |                         |    |
        !   IN MOMENT: heat_flux ~~~~|~~~~|~~~~   ,  water_flux ~~~~|~~~~|~~~~
        !  (BUT CHECK!)              |    |                         |    |
        !                            v   (+)                        v   (+)
        !

        ! The hosing passive tracer (ID 304) receives its source over the whole
        ! column when the anomaly is applied at depth; every other tracer, and
        ! 304 itself in surface mode, gets the usual single surface term.
        if (hosing_depth) then
            do nz=nzmin,nzmax
                tr(nz)= tr(nz)+impl_bcz(nz,n)
            end do
        else
            tr(nzmin)= tr(nzmin)+impl_bc(n)
        end if

        !_______________________________________________________________________
        ! The forward sweep algorithm to solve the three-diagonal matrix
        ! problem
        !
        !  | b_1 c_1 ...            |   |dTnew_1|
        !  | a_2 b_2 c_2 ...        |   |dTnew_2|
        !  |     a_3 b_3 c_3 ...    | * |dTnew_3| = RHS
        !  |         a_4 b_4 c_4 ...|   |dTnew_3|
        !  |              :         |   |   :   |
        !
        ! 1st: define new coefficents:
        !      --> c'_i = c_i/b_i                               ; i=1
        !          c'_i = c_i/(b_i-a_i*c'_i-1)                  ; i = 2,3,...,n-1
        !      --> rhs'_i = rhs_i/b_i                           ; i=1
        !          rhs'_i = (rhs_i-a_i*d'_i-1)/(b_i-a_i*c'_i-1) ; i = 2,3,...,n-1
        !
        ! 2nd: solution is optained by back substitution
        !      --> dTnew_n = rhs'_n
        !      --> dTnew_i = rhs'_i-c'_i*dTnew_i+1 ; i = n-1,n-2,...,1
        !
        ! initialize c-prime and s,t-prime
        cp(nzmin) = c(nzmin)/b(nzmin)
        tp(nzmin) = tr(nzmin)/b(nzmin)

        ! solve for vectors c-prime and t, s-prime
        do nz = nzmin+1,nzmax-1
            m = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
        end do

        ! start with back substitution
        tr(nzmax-1) = tp(nzmax-1)

        ! solve for x from the vectors c-prime and d-prime
        do nz = nzmax-2, nzmin, -1
            tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
        end do
        
        !_______________________________________________________________________
        ! update tracer
        ! tr ... dTnew = T^(n+0.5) - T*
        if (use_tracer_carry) then
            do nz=nzmin,nzmax-1
                acc  = trarr(nz,n) + tr(nz)
                tnew = real(acc, WP)
                tracers%work%carry(nz,n,tr_num) = tracers%work%carry(nz,n,tr_num) &
                                                + real((acc - tnew)*mesh%hnode_new(nz,n), WP)
                trarr(nz,n) = tnew
            end do
        else
        do nz=nzmin,nzmax-1
            ! trarr - before ... T*
            trarr(nz,n)=trarr(nz,n)+tr(nz)
        end do
        end if

    end do ! --> do n=1,myDim_nod2D
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
!$ACC END PARALLEL LOOP
#endif
end subroutine diff_ver_part_impl_ale
!
!
!===============================================================================
subroutine diff_ver_part_redi_expl(tracers, partit, mesh)
    IMPLICIT NONE
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    !___________________________________________________________________________
    integer                  :: n, k, elem, nz
    integer                  :: n2, nl1, ul1, nl2
    real(kind=WP)            :: Tx, Ty, flux, flux_below, zb, zb_below, zn, zn_up
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP_full), pointer   :: del_ttf(:,:)
    del_ttf => tracers%work%del_ttf
    if (.not. allocated(tr_xynodes)) then
        allocate(tr_xynodes(2,mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D))
        !$ACC ENTER DATA CREATE(tr_xynodes)
    end if

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, k, elem, nz, n2, nl1, ul1, nl2, Tx, Ty, flux, flux_below, zb, zb_below, zn, zn_up)
!$OMP DO
#else
!$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) IF(diff_on_device) PRIVATE(nl1, ul1) VECTOR_LENGTH(32)
#endif
    !___________________________________________________________________________
    do n=1, partit%myDim_nod2D
        nl1=mesh%nlevels_nod2D(n)-1
        ul1=mesh%ulevels_nod2D(n)
        !!PS do nz=1, nl1
        ! (OpenACC: the levels over the vector lanes, coalesced column access)
!$ACC LOOP VECTOR PRIVATE(k, elem, Tx, Ty)
        do nz=ul1, nl1
            Tx=0.0_WP
            Ty=0.0_WP
!$ACC LOOP SEQ
            do k=1, mesh%nod_in_elem2D_num(n)
                elem=mesh%nod_in_elem2D(k,n)
                !!PS if(nz.LE.(nlevels(elem)-1)) then
                if( nz.LE.(mesh%nlevels(elem)-1) .and. nz.GE.(mesh%ulevels(elem))) then
                    Tx=Tx+tr_xy(1,nz,elem)*mesh%elem_area(elem)
                    Ty=Ty+tr_xy(2,nz,elem)*mesh%elem_area(elem)
                endif
            end do
            tr_xynodes(1,nz,n)=tx/3.0_WP/mesh%areasvol(nz,n)
            tr_xynodes(2,nz,n)=ty/3.0_WP/mesh%areasvol(nz,n)
        end do
    end do

#ifndef ENABLE_OPENACC
!$OMP END DO
#else
!$ACC END PARALLEL LOOP
#endif
    ! no halo exchange of tr_xynodes is needed !
#ifndef ENABLE_OPENACC
!$OMP DO
#else
!$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(diff_on_device) PRIVATE(nz, nl1, ul1, flux, flux_below, zb, zb_below, zn, zn_up)
#endif
    do n=1, partit%myDim_nod2D
        nl1=mesh%nlevels_nod2D(n)-1
        ul1=mesh%ulevels_nod2D(n)
        ! bottom-up recurrence with scalars instead of per-column work arrays: zb = zbar_n(nz),
        ! zb_below = zbar_n(nz+1), zn = z_n(nz), zn_up = z_n(nz-1), flux_below = vd_flux(nz+1)
        zb_below  = mesh%zbar_n_bot(n)
        zn        = zb_below + mesh%hnode(nl1,n)/2.0_WP
        flux_below= 0.0_WP
        do nz=nl1, ul1+1, -1
            zb    = zb_below + mesh%hnode(nz,n)
            zn_up = zb       + mesh%hnode(nz-1,n)/2.0_WP
            flux  = (zn_up-zb)*(slope_tapered(1,nz-1,n)*tr_xynodes(1,nz-1,n)+slope_tapered(2,nz-1,n)*tr_xynodes(2,nz-1,n))*Ki(nz-1,n)
            flux  = flux+&
                    (zb-zn)  *(slope_tapered(1,nz,n)  *tr_xynodes(1,nz,n)  +slope_tapered(2,nz,n)  *tr_xynodes(2,nz,n))  *Ki(nz,n)
            flux  = flux/(zn_up-zn)*mesh%area(nz,n)
            del_ttf(nz,n) = del_ttf(nz,n) + (flux-flux_below) * dt/mesh%areasvol(nz,n)
            flux_below = flux
            zb_below   = zb
            zn         = zn_up
        enddo
        del_ttf(ul1,n) = del_ttf(ul1,n) + (0.0_WP-flux_below) * dt/mesh%areasvol(ul1,n)
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
!$ACC END PARALLEL LOOP
#endif
end subroutine diff_ver_part_redi_expl
!
!
!===============================================================================
subroutine diff_part_hor_redi(tracers, partit, mesh)
    IMPLICIT NONE
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: edge, n, k
    real(kind=WP)            :: deltaX1, deltaY1, deltaX2, deltaY2
    integer                  :: nl1, ul1, nl2, ul2, nl12, ul12, nz, el(2), elnodes(3), enodes(2)
    real(kind=WP), pointer   :: edge_flux(:,:)
    real(kind=WP)            :: c, Fx, Fy, Tx, Ty, Tx_z, Ty_z, SxTz, SyTz, Tz(2)
    real(kind=WP)            :: Kh, dz
    real(kind=WP)            :: isredi=0._WP
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP_full), pointer   :: del_ttf(:,:)
    del_ttf => tracers%work%del_ttf
    if (.not. allocated(tracers%work%edge_flux)) then
        allocate(tracers%work%edge_flux(mesh%nl-1, partit%myDim_edge2D))
        !$ACC ENTER DATA CREATE(tracers%work%edge_flux)
    end if
    edge_flux => tracers%work%edge_flux

    !___________________________________________________________________________
    if (Redi) isredi=1._WP
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, n, k, deltaX1, deltaY1, deltaX2, deltaY2, &
!$OMP                   nl1, ul1, nl2, ul2, nl12, ul12, nz, el, elnodes, enodes, &
!$OMP                             c, Fx, Fy, Tx, Ty, Tx_z, Ty_z, SxTz, SyTz, Tz, &
!$OMP                                                          Kh, dz)
#endif
    ! Pass 1: the flux through every edge, as seen from its first node (the second node gets minus it).
#ifndef ENABLE_OPENACC
!$OMP DO
#else
!$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(diff_on_device) &
!$ACC PRIVATE(deltaX1, deltaY1, deltaX2, deltaY2, nl1, ul1, nl2, ul2, nl12, ul12, nz, el, elnodes, enodes, &
!$ACC         c, Fx, Fy, Tx, Ty, SxTz, SyTz, Tz, Kh, dz)
#endif
    do edge=1, partit%myDim_edge2D
        !_______________________________________________________________________
        deltaX1=mesh%edge_cross_dxdy(1,edge)
        deltaY1=mesh%edge_cross_dxdy(2,edge)
        el=mesh%edge_tri(:,edge)
        enodes=mesh%edges(:,edge)
        nl1=mesh%nlevels(el(1))-1
        ul1=mesh%ulevels(el(1))
        elnodes=mesh%elem2d_nodes(:,el(1))
        !_______________________________________________________________________
        nl2=0
        ul2=0
        if (el(2)>0) then
            nl2=mesh%nlevels(el(2))-1
            ul2=mesh%ulevels(el(2))
            deltaX2=mesh%edge_cross_dxdy(3,edge)
            deltaY2=mesh%edge_cross_dxdy(4,edge)
        endif
        !_______________________________________________________________________
        nl12=min(nl1,nl2)
        ul12=max(ul1,ul2)
        ! flux of this edge over its column, accumulated per level from (A)-(E) below (no per-edge work array)
        do nz=min(ul1, merge(ul2, ul1, ul2>0)), max(nl1, nl2)
            edge_flux(nz,edge)=0.0_WP
        end do
        !_______________________________________________________________________
        ! (A)
        do nz=ul1, ul12-1
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=mesh%helem(nz, el(1))
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=tr_xy(1,nz,el(1))
            Ty=tr_xy(2,nz,el(1))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=(-deltaX1*Fy+deltaY1*Fx)*dz
            edge_flux(nz,edge) = edge_flux(nz,edge) + c
        end do
        !_______________________________________________________________________
        ! (B)
        if (ul2>0) then
            do nz=ul2,ul12-1
                Kh=sum(Ki(nz, enodes))/2.0_WP
                dz=mesh%helem(nz, el(2))
                Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
                SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
                SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
                Tx=tr_xy(1,nz,el(2))
                Ty=tr_xy(2,nz,el(2))
                Fx=Kh*(Tx+SxTz*isredi)
                Fy=Kh*(Ty+SyTz*isredi)
                c=(deltaX2*Fy-deltaY2*Fx)*dz
                edge_flux(nz,edge) = edge_flux(nz,edge) + c
            end do
        end if
        !_______________________________________________________________________
        ! (C)
        do nz=ul12, nl12
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=sum(mesh%helem(nz, el))/2.0_WP
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=0.5_WP*(tr_xy(1,nz,el(1))+tr_xy(1,nz,el(2)))
            Ty=0.5_WP*(tr_xy(2,nz,el(1))+tr_xy(2,nz,el(2)))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=((deltaX2-deltaX1)*Fy-(deltaY2-deltaY1)*Fx)*dz
            edge_flux(nz,edge) = edge_flux(nz,edge) + c
        enddo
        !_______________________________________________________________________
        ! (D)
        do nz=nl12+1, nl1
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=mesh%helem(nz, el(1))
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=tr_xy(1,nz,el(1))
            Ty=tr_xy(2,nz,el(1))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=(-deltaX1*Fy+deltaY1*Fx)*dz
            edge_flux(nz,edge) = edge_flux(nz,edge) + c
        end do
        !_______________________________________________________________________
        ! (E)
        do nz=nl12+1, nl2
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=mesh%helem(nz, el(2))
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=tr_xy(1,nz,el(2))
            Ty=tr_xy(2,nz,el(2))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=(deltaX2*Fy-deltaY2*Fx)*dz
            edge_flux(nz,edge) = edge_flux(nz,edge) + c
        end do
        !_______________________________________________________________________
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
!$ACC END PARALLEL LOOP
#endif
    ! Pass 2: each node gathers its edges' fluxes (mesh%nod_in_edge2D with sign; the first
    ! node gets +edge_flux, the second -edge_flux) in ascending edge order.
#ifndef ENABLE_OPENACC
!$OMP DO
#else
!$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(diff_on_device) &
!$ACC PRIVATE(k, edge, el, nl1, ul1, nl2, ul2, nl12, ul12, nz)
#endif
    do n=1, partit%myDim_nod2D+partit%eDim_nod2D
        do k=1, mesh%nod_in_edge2D_num(n)
            edge=mesh%nod_in_edge2D(k,n)
            el=mesh%edge_tri(:,edge)
            nl1=mesh%nlevels(el(1))-1
            ul1=mesh%ulevels(el(1))
            nl2=0
            ul2=0
            if (el(2)>0) then
                nl2=mesh%nlevels(el(2))-1
                ul2=mesh%ulevels(el(2))
            endif
            nl12=max(nl1,nl2)
            ul12 = ul1
            if (ul2>0) ul12=min(ul1,ul2)
            do nz=ul12, nl12
                del_ttf(nz,n)=del_ttf(nz,n)+mesh%nod_in_edge2D_sgn(k,n)*edge_flux(nz,edge)*dt/mesh%areasvol(nz,n)
            end do
        end do
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
!$ACC END PARALLEL LOOP
#endif
end subroutine diff_part_hor_redi
!
!
!===============================================================================
SUBROUTINE diff_part_bh(tr_num, dynamics, tracers, partit, mesh)
    IMPLICIT NONE
    integer,        intent(in),    target    :: tr_num
    type(t_dyn)   , intent(inout), target    :: dynamics
    type(t_tracer), intent(inout), target    :: tracers
    type(t_mesh)  , intent(in)   , target    :: mesh
    type(t_partit), intent(inout), target    :: partit
    integer                                  :: n, nz, ed, el(2), en(2), k, elem, nzmin, nzmax
    integer                                  :: elnodes1(3), elnodes2(3)
    real(kind=WP)                            :: u1, v1, len, vi, ww, tt(mesh%nl-1)
    real(kind=WP), pointer                   :: temporary_ttf(:,:)
    real(kind=WP), allocatable, target, save :: bh_ws(:,:)
    real(kind=WP), pointer                   :: edge_flux(:,:)
    real(kind=WP), pointer                   :: UV(:,:,:)
    real(kind=WP), pointer                   :: ttf(:,:)


    UV            => dynamics%uv(:,:,:)
    ttf           => tracers%data(tr_num)%values
    if (.not. allocated(bh_ws)) allocate(bh_ws(size(tracers%work%del_ttf,1), size(tracers%work%del_ttf,2)))
    temporary_ttf => bh_ws
    if (.not. allocated(tracers%work%edge_flux)) allocate(tracers%work%edge_flux(mesh%nl-1, partit%myDim_edge2D))
    edge_flux => tracers%work%edge_flux

!$OMP PARALLEL DO
        do n=1, partit%myDim_nod2D+partit%eDim_nod2D
           temporary_ttf(:, n)=0.0_8
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, ed, el, en, k, elem, nzmin, nzmax, u1, v1, len, vi, tt, ww, &
!$OMP elnodes1, elnodes2)
!$OMP DO
    DO ed=1, partit%myDim_edge2D!+eDim_edge2D
       if (partit%myList_edge2D(ed) > mesh%edge2D_in) cycle
       el=mesh%edge_tri(:,ed)
       en=mesh%edges(:,ed)
       len=sqrt(sum(mesh%elem_area(el)))
       nzmax = minval(mesh%nlevels(el))
       nzmin = maxval(mesh%ulevels(el))
       elnodes1=mesh%elem2d_nodes(:,el(1))
       elnodes2=mesh%elem2d_nodes(:,el(2))
       DO  nz=nzmin, nzmax-1
           u1=maxval(ttf(nz, elnodes1))-minval(ttf(nz, elnodes2))
           v1=minval(ttf(nz, elnodes1))-maxval(ttf(nz, elnodes2))
           vi=u1*u1+v1*v1
           tt(nz)=ttf(nz,en(1))-ttf(nz,en(2))
           vi=sqrt(max(tracers%data(tr_num)%gamma0_tra,            &
                   max(tracers%data(tr_num)%gamma1_tra*sqrt(vi),   &
                       tracers%data(tr_num)%gamma2_tra*     vi)    &
                                                         )*len)
           tt(nz)=tt(nz)*vi
       END DO
       edge_flux(nzmin:nzmax-1,ed)=tt(nzmin:nzmax-1)
    END DO
!$OMP END DO
    ! node gather (first node -tt, second node +tt), ascending edge order
!$OMP DO
    do n=1, partit%myDim_nod2D+partit%eDim_nod2D
       do k=1, mesh%nod_in_edge2D_num(n)
          ed=mesh%nod_in_edge2D(k,n)
          if (partit%myList_edge2D(ed) > mesh%edge2D_in) cycle
          el=mesh%edge_tri(:,ed)
          nzmax = minval(mesh%nlevels(el))
          nzmin = maxval(mesh%ulevels(el))
          do nz=nzmin, nzmax-1
             temporary_ttf(nz,n)=temporary_ttf(nz,n)-mesh%nod_in_edge2D_sgn(k,n)*edge_flux(nz,ed)
          end do
       end do
    end do
!$OMP END DO
!$OMP MASTER
    call exchange_nod(temporary_ttf, partit)
!$OMP END MASTER
!$OMP BARRIER
    ! ===========
    ! Second round:
    ! ===========
!$OMP DO
    DO ed=1, partit%myDim_edge2D!+eDim_edge2D
       if (partit%myList_edge2D(ed)>mesh%edge2D_in) cycle
          el=mesh%edge_tri(:,ed)
          en=mesh%edges(:,ed)
          len=sqrt(sum(mesh%elem_area(el)))
          nzmax = minval(mesh%nlevels(el))
          nzmin = maxval(mesh%ulevels(el))
          elnodes1=mesh%elem2d_nodes(:,el(1))
          elnodes2=mesh%elem2d_nodes(:,el(2))
          DO  nz=nzmin, nzmax-1
              u1=maxval(ttf(nz, elnodes1))-minval(ttf(nz, elnodes2))
              v1=minval(ttf(nz, elnodes1))-maxval(ttf(nz, elnodes2))
              vi=u1*u1+v1*v1
              tt(nz)=temporary_ttf(nz,en(1))-temporary_ttf(nz,en(2))
              vi=sqrt(max(tracers%data(tr_num)%gamma0_tra,            &
                      max(tracers%data(tr_num)%gamma1_tra*sqrt(vi),   &
                          tracers%data(tr_num)%gamma2_tra*     vi)    &
                                                            )*len)
              tt(nz)=-tt(nz)*vi*dt
          END DO
          edge_flux(nzmin:nzmax-1,ed)=tt(nzmin:nzmax-1)
    END DO
!$OMP END DO
!$OMP DO
    do n=1, partit%myDim_nod2D+partit%eDim_nod2D
       do k=1, mesh%nod_in_edge2D_num(n)
          ed=mesh%nod_in_edge2D(k,n)
          if (partit%myList_edge2D(ed) > mesh%edge2D_in) cycle
          el=mesh%edge_tri(:,ed)
          nzmax = minval(mesh%nlevels(el))
          nzmin = maxval(mesh%ulevels(el))
          do nz=nzmin, nzmax-1
             ttf(nz,n)=ttf(nz,n)-mesh%nod_in_edge2D_sgn(k,n)*edge_flux(nz,ed)/mesh%area(nz,n)
          end do
       end do
    end do
!$OMP END DO
!$OMP END PARALLEL
call exchange_nod(ttf, partit)
!$OMP BARRIER
end subroutine diff_part_bh
!
!
!===============================================================================
! this function returns a boundary conditions for a specified tracer ID and surface node
! ID = 0 and 1 are reserved for temperature and salinity
! MB: mesh, sst, sss, and aice are only needed for transient tracers
FUNCTION bc_surface(n, id, sval, nzmin, partit, mesh, sst, sss, aice)
  implicit none

  integer,       intent(in)            :: n, id, nzmin
  real(kind=WP), intent(in)            :: sval, sst, sss, aice
  type(t_partit),intent(inout), target :: partit
  type(t_mesh),  intent(in), target    :: mesh

  REAL(kind=WP)                        :: bc_surface
  character(len=10)                    :: id_string
  real(kind=WP), dimension(:), pointer :: a_ice  !! MB: where is this needed?

  if (use_transit) then
#if defined (__cpl_oasis)
    ! SLP and wind speed in coupled setups. This is a makeshift solution
    ! as long as the true values are not provided by the AGCM / OASIS.
    press_a = mean_slp
    wind_2  = speed_2(stress_atmoce_x(n), stress_atmoce_y(n))
#else
    press_a = press_air(n)
    wind_2  = u_wind(n)**2 + v_wind(n)**2
#endif
    ! Atmospheric input of bomb-14C, CFC-11, CFC-12, and SF6 varies with latitude. To that effect specify
    y_abc = mesh%geo_coord_nod2D(2,n) / rad  ! latitude of atmospheric tracer input
    yy_nh = (10. - y_abc) * 0.05             ! interpolation weight for tropical CFC-12 and SF6 values
  end if ! use_transit

  SELECT CASE (id)
    !___temperature_____________________________________________________________
    CASE (1)
        bc_surface=-dt*(heat_flux(n)/vcpw + sval*water_flux(n)*is_nonlinfs)
    
    !___salinity________________________________________________________________
    CASE (2)
        ! --> real_salt_flux(:): salt flux due to containment/releasing of salt
        !     by forming/melting of sea ice
        ! use_salt_anomaly background-dilution term (S_ref=0 unless enabled, so
        ! this reduces to the original salinity BC). WHY it is needed and where a
        ! residual error comes from:
        !   Storing S = S' + S_ref should be invisible in DP: advection is linear
        !   and, by continuity, adv(S_ref) = -S_ref*dh/dt exactly cancels the
        !   -S_ref*dh/dt from the anomaly transform, leaving the same equation in
        !   S'. That holds ONLY if the free-surface advection preserves a constant
        !   exactly. FESOM's does not: with no correction the DP anomaly-vs-absolute
        !   salt error is ~1.3e-3, a constancy error for the offset S_ref localised
        !   where the surface freshwater flux enters. The S_ref*water_flux term
        !   below cancels ~99.4% of it, leaving a ~8e-6/step residual (the part not
        !   proportional to water_flux: the grid-divergence / FCT-limiter part).
        !   An exact fix would carry S'+S_ref in the surface ADVECTIVE flux across
        !   all schemes, not add a surface source here (a per-layer thickness term
        !   was tried and does not help; it makes it worse combined with this term).
        ! real_salt_flux / relax_salt / virtual_salt are already anomaly-consistent
        ! (real_salt_flux is built from the corrected absolute ice gather; relax_salt
        ! subtracts S_ref from the absolute Ssurf; virtual_salt is 0 in zstar/zlevel).
        bc_surface= dt*(virtual_salt(n) & !--> is zeros for zlevel/zstar
                    + relax_salt(n) + (real_salt_flux(n) + S_ref_anomaly*water_flux(n))*is_nonlinfs)
            
    !___Transient tracers (cases ##6,11,12,14,39)__________________________________
    CASE (6) ! SF6
      if (anthro_transit) then
        ! Select atmospheric input values corresponding to the latitude
        ! Annual values are interpolated to monthly values, this is omitted in the last simulation year
        if (y_abc > 10.)  then       ! Northern hemisphere
          xsf6_a = xsf6_nh(ti_transit)
          if (ti_transit < length_transit) xsf6_a = xsf6_a + month * (xsf6_nh(ti_transit + 1) - xsf6_a) / 12.
        else if (y_abc <- 10.) then  ! Southern hemisphere
          xsf6_a = xsf6_sh(ti_transit)
          if (ti_transit < length_transit) xsf6_a = xsf6_a + month * (xsf6_sh(ti_transit + 1) - xsf6_a) / 12.
        else                         ! Tropical zone, interpolate between NH and SH
          xsf6_a = (1 - yy_nh) * xsf6_nh(ti_transit) + yy_nh * xsf6_sh(ti_transit)
          if (ti_transit < length_transit) &
            xsf6_a = xsf6_a + month * ((1 - yy_nh) * xsf6_nh(ti_transit + 1) + yy_nh * xsf6_sh(ti_transit + 1) - xsf6_a) / 12.
        end if
      else
        !  Constant (global-mean) namelist values are taken
      end if
      ! Local air-sea exchange gas flux of SF6 (m / s):
      bc_surface = dt * (gas_flux("sf6", sst, sss, wind_2, aice, press_a, xsf6_a, sval) - sval * water_flux(n) * is_nonlinfs)
    CASE (11) ! CFC-11
      if (anthro_transit) then
        ! Select atmospheric input values corresponding to the latitude
        ! Annual values are interpolated to monthly values, this is omitted in the last simulation year
        if (y_abc > 10.)  then       ! Northern hemisphere
           xf11_a = xf11_nh(ti_transit)
           if (ti_transit < length_transit) xf11_a = xf11_a + month * (xf11_nh(ti_transit + 1) - xf11_a) / 12.
        else if (y_abc <- 10.) then  ! Southern hemisphere
           xf11_a = xf11_sh(ti_transit)
           if (ti_transit < length_transit) xf11_a = xf11_a + month * (xf11_sh(ti_transit + 1) - xf11_a) / 12.
        else                         ! Tropical zone, interpolate between NH and SH
           xf11_a = (1 - yy_nh) * xf11_nh(ti_transit) + yy_nh * xf11_sh(ti_transit)
           if (ti_transit < length_transit) &
             xf11_a = xf11_a + month * ((1 - yy_nh) * xf11_nh(ti_transit + 1) + yy_nh * xf11_sh(ti_transit + 1) - xf11_a) / 12.
        end if
      else
        ! Constant (global-mean) namelist values are taken
      end if
      ! Local air-sea exchange gas flux of CFC-12 (m / s):
      bc_surface = dt * (gas_flux("f11", sst, sss, wind_2, aice, press_a, xf11_a, sval) - sval * water_flux(n) * is_nonlinfs)
    CASE (12) ! CFC-12
      if (anthro_transit) then
        ! Select atmospheric input values corresponding to the latitude
        ! Annual values are interpolated to monthly values, this is omitted in the last simulation year
        if (y_abc > 10.)  then       ! Northern hemisphere
           xf12_a = xf12_nh(ti_transit)
           if (ti_transit < length_transit) xf12_a = xf12_a + month * (xf12_nh(ti_transit + 1) - xf12_a) / 12.
        else if (y_abc <- 10.) then  ! Southern hemisphere
           xf12_a = xf12_sh(ti_transit)
           if (ti_transit < length_transit) xf12_a = xf12_a + month * (xf12_sh(ti_transit + 1) - xf12_a) / 12.
        else                         ! Tropical zone, interpolate between NH and SH
           xf12_a = (1 - yy_nh) * xf12_nh(ti_transit) + yy_nh * xf12_sh(ti_transit)
           if (ti_transit < length_transit) &
             xf12_a = xf12_a + month * ((1 - yy_nh) * xf12_nh(ti_transit + 1) + yy_nh * xf12_sh(ti_transit + 1) - xf12_a) / 12.
        end if
      else
        ! Constant (global-mean) namelist values are taken
      end if
      ! Local air-sea exchange gas flux of CFC-12 (m / s):
      bc_surface = dt * (gas_flux("f12", sst, sss, wind_2, aice, press_a, xf12_a, sval) - sval * water_flux(n) * is_nonlinfs)
    CASE (14) ! Radiocarbon (more precisely, fractionation-corrected 14C/C):
      if (anthro_transit) then
        ! Select atmospheric input values corresponding to the latitude
        if (y_abc > 30.)  then       ! Northern hemisphere
          r14c_a = r14c_nh(ti_transit)
        else if (y_abc <- 30.) then  ! Southern hemisphere
          r14c_a = r14c_sh(ti_transit)
        else                         ! Tropical zone (values prescribed in input file)
          r14c_a = r14c_tz(ti_transit)
        end if
        xCO2_a = xCO2_ti(ti_transit)
      else if (paleo_transit) then
        r14c_a = r14c_ti(ti_transit)
        xCO2_a = xCO2_ti(ti_transit)
      else
        ! Constant (global-mean) namelist values are taken
      end if
      ! Local isotopic 14CO2/CO2 air-sea exchange flux (m / s)
      bc_surface = dt * (iso_flux("co2", sst, sss, wind_2, aice, press_a, xco2_a, r14c_a, sval, dic_0) - sval * water_flux(n) * is_nonlinfs)
    CASE (39) ! Argon-39 (fractionationation-corrected 39Ar/Ar)
      ! Local isotopic 39Ar/Ar air-sea exchange flux (m / s)
      bc_surface = dt * (iso_flux("arg", sst, sss, wind_2, aice, press_a, xarg_a, r39ar_a, sval, arg_0) - sval * water_flux(n) * is_nonlinfs)
!---  Done with boundary conditions for transient tracers.
#if defined(__recom)
    CASE (1001) ! DIN
        if (use_MEDUSA .and. add_loopback) then  ! OG: add is_MEDUSA_loopback flag is_MEDUSA_loopback flag * lb_flux(n,1)
            bc_surface= dt*(AtmNInput(n) + virtual_din(n) + RiverDIN2D(n)   * is_riverinput &
                                         + ErosionTON2D(n) * is_erosioninput + lb_flux(n,1))
        else
            bc_surface= dt*(AtmNInput(n) + virtual_din(n) + RiverDIN2D(n)   * is_riverinput &
                                         + ErosionTON2D(n) * is_erosioninput)
        end if

    CASE (1002) ! DIC
        if (use_MEDUSA .and. add_loopback) then
            bc_surface= dt*(GloCO2flux_seaicemask(n) + virtual_dic(n) &
                                + RiverDIC2D(n)   * is_riverinput   &
                                + ErosionTOC2D(n) * is_erosioninput &
                                + lb_flux(n,2) + lb_flux(n,5))
        else
            bc_surface= dt*(GloCO2flux_seaicemask(n) + virtual_dic(n) &
                                + RiverDIC2D(n)   * is_riverinput   &
                                + ErosionTOC2D(n) * is_erosioninput)
       end if

    CASE (1003) ! Alk
        if (use_MEDUSA .and. add_loopback) then
            bc_surface= dt*(virtual_alk(n) + relax_alk(n)       &
                            + RiverAlk2D(n) * is_riverinput     &
                            + lb_flux(n,3) + lb_flux(n,5)*2) !CaCO3:Alk burial=1:2
        else
            bc_surface= dt*(virtual_alk(n) + relax_alk(n)       &
                            + RiverAlk2D(n) * is_riverinput)
        end if
    CASE (1004:1010)
        bc_surface=0.0_WP
    CASE (1011) ! DON
        bc_surface= dt*RiverDON2D(n) * is_riverinput
    CASE (1012) ! DOC
        bc_surface= dt*RiverDOC2D(n) * is_riverinput
    CASE (1013:1017)
        bc_surface=0.0_WP
    CASE (1018) ! DSi
        if (use_MEDUSA .and. add_loopback) then
           bc_surface=dt*(virtual_dsi(n) + RiverDSi2D(n)   * is_riverinput &
                        + ErosionTSi2D(n) * is_erosioninput      &
                        + lb_flux(n,4))
        else
            bc_surface=dt*(virtual_dsi(n) + RiverDSi2D(n) * is_riverinput + ErosionTSi2D(n) * is_erosioninput)
        end if

    CASE (1019) ! Fe
        if (useRivFe) then
            bc_surface= dt*(AtmFeInput(n) + virtual_dfe(n) + RiverFe(n))
        else
           bc_surface= dt*(AtmFeInput(n) + virtual_dfe(n))
        end if
    CASE (1020:1021) ! Cal
        bc_surface=0.0_WP
    CASE (1022) ! OXY
        bc_surface= dt*(GloO2flux_seaicemask(n) + virtual_oxy(n))
!        bc_surface=0.0_WP
    CASE (1023:1037)
        bc_surface=0.0_WP  ! OG added bc for recom fields
    CASE (1302) ! Before (1037) ! DIC_13
#if defined (__ciso)
         if (ciso) then
            if (use_MEDUSA .and. add_loopback) then
               bc_surface= dt*(GloCO2flux_seaicemask_13(n) &
                           + lb_flux(n,6) + lb_flux(n,7))
            else
               bc_surface= dt*(GloCO2flux_seaicemask_13(n))
            end if
         else
            bc_surface=0.0_WP
         end if
#endif
    CASE (1305:1321)
         bc_surface=0.0_WP ! organic 13C
    CASE (1402) ! Before (1034) ! DIC_14
#if defined (__ciso)
         if (ciso .and. ciso_14) then
             if (use_MEDUSA .and. add_loopback .and. ciso_organic_14) then
                 bc_surface= dt*(GloCO2flux_seaicemask_14(n) &
                             + lb_flux(n,8) + lb_flux(n,9))
             else
                 bc_surface= dt*GloCO2flux_seaicemask_14(n)
             end if
         else
             bc_surface=0.0_WP
         end if
#endif
    CASE (1405:1421)
         bc_surface=0.0_WP ! organic 14C
#endif
    CASE (101) ! apply boundary conditions to tracer ID=101
        bc_surface= dt*(prec_rain(n))! - real_salt_flux(n)*is_nonlinfs)
!---wiso-code
    CASE (102) ! apply boundary conditions to tracer ID=101 (H218O)
        bc_surface = dt*wiso_flux_oce(n,1)
    CASE (103)  ! apply boundary conditions to tracer ID=102 (HDO)
        bc_surface = dt*wiso_flux_oce(n,2)
    CASE (104)  ! apply boundary conditions to tracer ID=103 (H216O)
        bc_surface = dt*wiso_flux_oce(n,3)
!---wiso-code-end
!---age-code
    CASE (100)
        !bc_surface=-dt*(sval*water_flux(n)*is_nonlinfs)
        bc_surface=0.0_WP
!---age-code-end
    CASE (301)
        bc_surface=0.0_WP
    CASE (302)
        bc_surface=0.0_WP
    CASE (303)
        bc_surface=0.0_WP
    CASE (304) ! hosing passive tracer; nzmin is the level being filled
        if (use_hosing .and. trim(hosing_mode)=='depth') then
            bc_surface= dt*(hosing_flux3D(nzmin,n))
        else
            bc_surface= dt*(hosing_flux(n))
        end if
    CASE (501) ! ice-shelf water due to basal melting
        if (nzmin==1) then
           bc_surface = 0.0_WP
        else
           bc_surface = dt*water_flux(n)*(sval-1.0)
        end if
    CASE DEFAULT
      if (partit%mype==0) then
         write (id_string, "(I3)") id
         if (partit%mype==0) write(*,*) 'invalid ID '//trim(id_string)//' specified in boundary conditions'
         if (partit%mype==0) write(*,*) 'the model will stop!'
      end if
      call par_abort(partit%MPI_COMM_FESOM, partit%mype)
      stop
  END SELECT
  RETURN
END FUNCTION

!===============================================================================
! This function returns a boundary conditions for a specified transient tracer ID and surface node.
! Different to function bc_surface, SST, SSS, and sea ice concentrations are always needed as
! auxiliary variable
FUNCTION transit_bc_surface(n, id, sst, sss, aice, sval, nzmin, partit, mesh)
  implicit none

  integer,       intent(in)            :: n, id, nzmin
  real(kind=WP), intent(in)            :: sst, sss, aice, sval
  type(t_partit),intent(inout), target :: partit
  type(t_mesh),  intent(in), target    :: mesh
  REAL(kind=WP)                        :: transit_bc_surface
  character(len=10)                    :: id_string


  !  --> is_nonlinfs=1.0 for zelvel,zstar ....
  !  --> is_nonlinfs=0.0 for linfs

#if defined (__cpl_oasis)
! SLP and wind speed in coupled setups. This is a makeshift solution
! as long as the true values are not provided by the AGCM / OASIS.
  press_a = mean_slp
  wind_2  = speed_2(stress_atmoce_x(n), stress_atmoce_y(n))
#else
  press_a = press_air(n)
  wind_2  = u_wind(n)**2 + v_wind(n)**2
#endif

! The atmospheric input of bomb 14C, CFC-12, and SF6 depends on latitude. To that effect specify
  y_abc = mesh%geo_coord_nod2D(2,n) / rad  ! latitude of atmospheric tracer input
  yy_nh = (10. - y_abc) * 0.05             ! interpolation weight for tropical tracer values


  SELECT CASE (id)

!   Boundary conditions for additional (transient) tracers (14C, 39Ar, CFC-12, and SF6)
    CASE (14) !   Radiocarbon (more precisely, fractionation-corrected 14C/C):
      if (anthro_transit) then
!       Select atmospheric input values corresponding to the latitude
        if (y_abc > 30.)  then
!         Northern Hemisphere
          r14c_a = r14c_nh(ti_transit)
        else if (y_abc <- 30.) then
!         Southern Hemisphere
          r14c_a = r14c_sh(ti_transit)
        else
!         Tropical zone
          r14c_a = r14c_tz(ti_transit)
        end if
        xCO2_a = xCO2_ti(ti_transit)
      else if (paleo_transit) then
        r14c_a = r14c_ti(ti_transit)
        xCO2_a = xCO2_ti(ti_transit)
      else
!       Constant (global-mean) namelist values are taken
      end if
!     Local isotopic 14CO2/CO2 air-sea exchange flux (in m / s),
!     since F14C is normalized to atmospheric (water) values the isotopic flux has to be
!     corrected for precipitation or evaporation fluxes with different isotopic signatures.
      transit_bc_surface = dt * (iso_flux("co2", sst, sss, wind_2, aice, press_a, xco2_a, r14c_a, sval, dic_0)   &
                                 - sval * water_flux(n) * is_nonlinfs)

    CASE (39) ! Argon-39 (fractionationation-corrected 39Ar/Ar)
!     Local isotopic 39Ar/Ar air-sea exchange flux (in m / s),
!     since F39Ar is normalized to atmospheric (water) values the isotopic flux has to be
!     corrected for precipitation or evaporation fluxes with different isotopic signatures.
      transit_bc_surface = dt * (iso_flux("arg", sst, sss, wind_2, aice, press_a, xarg_a, r39ar_a, sval, arg_0)  &
                                 - sval * water_flux(n) * is_nonlinfs)

    CASE (12) ! CFC-12
      if (anthro_transit) then
!       Select atmospheric input values corresponding to the latitude
!       Annual values are interpolated to monthly values, this is omitted in the last simulation year
        if (y_abc > 10.)  then       ! Northern Hemisphere
!          Northern Hemisphere
           xf12_a = xf12_nh(ti_transit)
           if (ti_transit < length_transit) xf12_a = xf12_a + month * (xf12_nh(ti_transit + 1) - xf12_a) / 12.
        else if (y_abc <- 10.) then
!          Southern Hemisphere
           xf12_a = xf12_sh(ti_transit)
           if (ti_transit < length_transit) xf12_a = xf12_a + month * (xf12_sh(ti_transit + 1) - xf12_a) / 12.
        else
!          Tropical zone, interpolate between NH and SH
           xf12_a = (1 - yy_nh) * xf12_nh(ti_transit) + yy_nh * xf12_sh(ti_transit)
           if (ti_transit < length_transit) &
             xf12_a = xf12_a + month * ((1 - yy_nh) * xf12_nh(ti_transit + 1) + yy_nh * xf12_sh(ti_transit + 1) - xf12_a) / 12.
        end if
      else
!       Constant (global-mean) namelist values are taken
      end if

!     Local air-sea exchange gas flux of CFC-12 (in m / s):
      transit_bc_surface = dt * (gas_flux("f12", sst, sss, wind_2, aice, press_a, xf12_a, sval)  &
                                 - sval * water_flux(n) * is_nonlinfs)

    CASE (6) ! SF6
      if (anthro_transit) then
!       Select atmospheric input values corresponding to the latitude
!       Annual values are interpolated to monthly values, this is omitted in the last simulation year
        if (y_abc > 10.)  then       ! Northern Hemisphere
!         Northern Hemisphere
          xsf6_a = xsf6_nh(ti_transit)
          if (ti_transit < length_transit) xsf6_a = xsf6_a + month * (xsf6_nh(ti_transit + 1) - xsf6_a) / 12.
        else if (y_abc <- 10.) then
!         Southern Hemisphere
          xsf6_a = xsf6_sh(ti_transit)
          if (ti_transit < length_transit) xsf6_a = xsf6_a + month * (xsf6_sh(ti_transit + 1) - xsf6_a) / 12.
        else
!         Tropical zone, interpolate between NH and SH
          xsf6_a = (1 - yy_nh) * xsf6_nh(ti_transit) + yy_nh * xsf6_sh(ti_transit)
          if (ti_transit < length_transit) &
            xsf6_a = xsf6_a + month * ((1 - yy_nh) * xsf6_nh(ti_transit + 1) + yy_nh * xsf6_sh(ti_transit + 1) - xsf6_a) / 12.
        end if
      else
!       Constant (global-mean) namelist values are taken
      end if

!     Local air-sea exchange gas flux of SF6 (in m / s):
      transit_bc_surface = dt * (gas_flux("sf6", sst, sss, wind_2, aice, press_a, xsf6_a, sval)  &
                                 - sval * water_flux(n) * is_nonlinfs)

!   Done with boundary conditions for (transient) tracers.
  END SELECT
  RETURN

END FUNCTION

!===============================================================================
! divide the range specified by indexcount into fesom_group_count equal slices and calculate
! the start_index and end_index for the given fesom_group_id.
! if necessary to compensate for fragmentation, the end index of the first n slices
! might be one higher than for the remaining slices. this is indicated by end_index_is_one_higher
subroutine calc_slice(index_count, fesom_group_count, fesom_group_id, start_index, end_index, index_count_in_group, end_index_is_one_higher)
!   use g_config

    implicit none
    integer, intent(in)      :: index_count
    integer, intent(in)      :: fesom_group_count
    integer, intent(in)      :: fesom_group_id
    integer, intent(out)     :: start_index
    integer, intent(out)     :: end_index
    integer, intent(out)     :: index_count_in_group
    logical, intent(out)     :: end_index_is_one_higher

    integer                  :: group_id_limit_to_adjust_end_index

    index_count_in_group               = index_count / fesom_group_count
    group_id_limit_to_adjust_end_index = mod(index_count, fesom_group_count)
    start_index                        = (fesom_group_id * index_count_in_group) + 1

! adjust loop start and number of loop iterations by 1 if necessary
    if(fesom_group_id < group_id_limit_to_adjust_end_index) then
      start_index = start_index + fesom_group_id
      index_count_in_group = index_count_in_group + 1
      end_index_is_one_higher = .true.
    else
      start_index = start_index + group_id_limit_to_adjust_end_index
      end_index_is_one_higher = .false.
    end if

    end_index  = start_index + index_count_in_group - 1
end subroutine calc_slice

end module oce_ale_tracer_module
