module iceberg_ocean_coupling
 USE MOD_MESH
 use MOD_PARTIT
 use MOD_ICE
 USE MOD_DYN
 use iceberg_params

 public :: reset_ib_fluxes
 public :: prepare_icb2fesom
 public :: icb2fesom
 public :: cap_ibhf_n

 contains


subroutine reset_ib_fluxes()
    use o_param

! kh 18.03.21 not really used here
    use g_config
    use iceberg_params

    wave_erosion_potential  = 0.0
    linit_wave_erosion_pot  = .true.
    ibfwbv  = 0.0
    ibfwb   = 0.0
    ibfwl   = 0.0
    ibfwe   = 0.0
    ibiron  = 0.0
    ibhf    = 0.0
    ibhf_n  = 0.0_WP
end subroutine reset_ib_fluxes


subroutine prepare_icb2fesom(mesh, partit, ib,i_have_element,localelement,depth_ib,height_ib_single)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param
    use MOD_MESH
    use g_config
    use MOD_PARTIT
    use iceberg_params

    implicit none

    logical                 :: i_have_element
    real, intent(in)        :: depth_ib, height_ib_single
    real                    :: lev_low, lev_up
    integer                 :: localelement
    integer                 :: iceberg_node  
    integer, dimension(3)   :: ib_nods_in_ib_elem
    integer                 :: num_ib_nods_in_ib_elem
    real                    :: dz
    real, dimension(:), allocatable      :: tot_area_nods_in_ib_elem
    integer                 :: i, j, k, ib
    integer, dimension(3)   :: idx_d
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    ! Also guard localelement: when an iceberg crosses a PE boundary during
    ! iceberg_step1, i_have_element can remain .true. while local_idx_of()
    ! returned 0 for the new element.  elem2d_nodes(i, 0) is an out-of-bounds
    ! read that produces a garbage node index and may divide by zero area -> NaN.
    if(i_have_element .and. localelement > 0) then
        dz = 0.0
        allocate(tot_area_nods_in_ib_elem(mesh%nl))

        num_ib_nods_in_ib_elem=0
        tot_area_nods_in_ib_elem=0.0
        idx_d = 0

        do i=1,3
            iceberg_node=elem2d_nodes(i,localelement)

            ! LOOP: consider all neighboring pairs (n_up,n_low) of 3D nodes
            ! below n2..
            !innerloop: do k=1, nl+1
            do k=1, nlevels_nod2D(iceberg_node)
                idx_d(i) = k
                lev_up  = mesh%zbar_3d_n(k, iceberg_node)

                if( k==nlevels_nod2D(iceberg_node) ) then
                    lev_low = mesh%zbar_n_bot(iceberg_node)
                else
                    lev_low = mesh%zbar_3d_n(k+1, iceberg_node)
                end if

                if( abs(lev_low)==abs(lev_up) ) then
                    idx_d(i) = idx_d(i) - 1
                    exit
                else if( abs(lev_low)>=abs(depth_ib) ) then
                    exit
                else
                    cycle
                end if
            end do

            if (iceberg_node<=mydim_nod2d) then
                ib_nods_in_ib_elem(i)           = iceberg_node
                num_ib_nods_in_ib_elem          = num_ib_nods_in_ib_elem + 1
                tot_area_nods_in_ib_elem        = tot_area_nods_in_ib_elem + mesh%area(:,iceberg_node)
            else
                ib_nods_in_ib_elem(i)           = 0
            end if
        end do

        do i=1, 3
            iceberg_node=ib_nods_in_ib_elem(i)
            ! Halo nodes are marked 0 above. Test that before indexing, since
            ! Fortran does not promise to short-circuit .or. and Intel does not.
            if (iceberg_node <= 0) cycle
            ! Skip nodes with no ocean column; do NOT skip valid cavity nodes
            if (ulevels_nod2d(iceberg_node) == 0) cycle

            if (iceberg_node>0) then
                ! Guard: tot_area at level 1 is zero when all element nodes are cavity nodes
                ! (mesh%area(1,cavity_node)==0), which would give division by zero / NaN.
                if (tot_area_nods_in_ib_elem(1) > 0.0) then
                    ibfwbv(iceberg_node) = ibfwbv(iceberg_node) - fwbv_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                    ibfwb(iceberg_node) = ibfwb(iceberg_node) - fwb_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                    ibfwl(iceberg_node) = ibfwl(iceberg_node) - fwl_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                    ibfwe(iceberg_node) = ibfwe(iceberg_node) - fwe_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                    ! LA 2026 -- iron released with the meltwater, distributed over
                    ! the element nodes exactly like the freshwater [mol m-2 s-1]
                    if (use_icb_iron) then
                        ibiron(iceberg_node) = ibiron(iceberg_node) - iron_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                    end if
                end if
                !ibhf(iceberg_node) = ibhf(iceberg_node) - hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(1)

                ! Guard against idx_d=0 (can happen if nlevels_nod2D is 0 or loop didn't execute)
                if (idx_d(i) <= 0) cycle

                do j=1,idx_d(i)
                    lev_up  = mesh%zbar_3d_n(j, iceberg_node)
                    if( j==nlevels_nod2D(iceberg_node) ) then
                        lev_low = mesh%zbar_n_bot(iceberg_node)
                    else
                        lev_low = mesh%zbar_3d_n(j+1, iceberg_node)
                    end if
                    dz = abs( lev_low - lev_up )
                    if( (abs(lev_low)>=abs(depth_ib)) .and. (abs(lev_up)<abs(depth_ib)) ) then 
                        dz = abs(abs(lev_up) - abs(depth_ib))
                    end if              
                   
                    ! Also guard tot_area(j): at levels where only some nodes extend,
                    ! the remaining nodes contribute zero area, making tot_area very small.
                    if( abs(depth_ib) > 0.0 .and. tot_area_nods_in_ib_elem(j) > 0.0 ) then
                        ibhf_n(j,iceberg_node) = ibhf_n(j,iceberg_node) & 
                                                    - (hfbv_flux_ib(ib,j)+hfl_flux_ib(ib,j)) & 
                                                    / tot_area_nods_in_ib_elem(j)
                    end if
                end do
                
                if( idx_d(i) > 1 ) then
                    if (tot_area_nods_in_ib_elem(idx_d(i)) > 0.0) &
                        ibhf_n(idx_d(i),iceberg_node) = ibhf_n(idx_d(i),iceberg_node) - 0.5 * hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(idx_d(i))
                    if (tot_area_nods_in_ib_elem(idx_d(i)-1) > 0.0) &
                        ibhf_n(idx_d(i)-1,iceberg_node) = ibhf_n(idx_d(i)-1,iceberg_node) - 0.5 * hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(idx_d(i)-1)
                else if( idx_d(i) == 1 ) then
                    if (tot_area_nods_in_ib_elem(idx_d(i)) > 0.0) &
                        ibhf_n(idx_d(i),iceberg_node) = ibhf_n(idx_d(i),iceberg_node) - hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(idx_d(i))
                end if
                
                ! Wave erosion heat applied at level 1 only on open-ocean nodes (ulevels==1).
                ! Cavity nodes (ulevels>1) skip nz=1 in oce_ale_tracer, so setting ibhf_n(1)
                ! for them wastes the flux and overcounts the heat on the open-ocean node.
                if( height_ib_single .ne. 0.0 .and. tot_area_nods_in_ib_elem(1) > 0.0 &
                        .and. ulevels_nod2d(iceberg_node) == 1 ) then
                    ibhf_n(1,iceberg_node) = ibhf_n(1,iceberg_node) - hfe_flux_ib(ib) &
                            !* abs(height_ib_single) &
                            !* ((abs(height_ib_single)-abs(depth_ib))/abs(height_ib_single)) &
                            / tot_area_nods_in_ib_elem(1)
                end if
            end if
        end do
    end if
end subroutine prepare_icb2fesom

subroutine icb2fesom(mesh, partit, ice)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param

! kh 18.03.21 specification of structure used
    use o_arrays, only: water_flux, heat_flux, wiso_flux_oce, Tclim_ib, is_nonlinfs
    use MOD_MESH
    use g_config
    use MOD_PARTIT
    use MOD_ICE
    use iceberg_params
    integer                                 :: n
    real(kind=WP), dimension(:)  , pointer :: fresh_wa_flux, net_heat_flux
type(t_ice)   , intent(inout), target :: ice
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  if (use_cavity) then
!$OMP PARALLEL DO
    do n=1, myDim_nod2d+eDim_nod2D
        if (ulevels_nod2d(n) > 1) cycle
        if (.not.turn_off_fw) then
                water_flux(n)  = water_flux(n) - (ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n)) !* steps_per_ib_step
        end if
    end do
!$OMP END PARALLEL DO
  else
!$OMP PARALLEL DO
    do n=1, myDim_nod2d+eDim_nod2D
        if (.not.turn_off_fw) then
                water_flux(n)  = water_flux(n) - (ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n)) !* steps_per_ib_step
        end if
    end do
!$OMP END PARALLEL DO
  end if

!---enthalpy-correction-begin
! The ALE kinematic boundary condition in oce_ale_tracer.F90:
!   bc_surface = -dt*(heat_flux(n)/vcpw + sval*water_flux(n)*is_nonlinfs)
! uses sval=T_ocean for the entire water_flux, including the iceberg FW
! contribution. Iceberg meltwater enters at T_melt~=0 degC, not T_ocean.
! When T_ocean < 0 degC the wrong sign of the kinematic FW term creates
! spurious cooling: d(T*h)/dt includes T_ocean*dh instead of T_melt*dh=0.
!
! Fix: add vcpw*T_ocean*fw_iceberg to heat_flux. This cancels the erroneous
! T_ocean*fw_ib term in bc_surface, leaving a net contribution of zero
! (correct for T_melt=0 degC). The sign ensures warming for T_ocean<0 and
! cooling for T_ocean>0, matching the physics of mixing with 0-degC water.
  if (.not. (turn_off_hf .or. turn_off_fw)) then
!$OMP PARALLEL DO
    do n=1, myDim_nod2d+eDim_nod2D
        if (use_cavity .and. ulevels_nod2d(n) > 1) cycle
        heat_flux(n) = heat_flux(n) + vcpw * is_nonlinfs * Tclim_ib(1, n) &
                       * (ibfwb(n) + ibfwl(n) + ibfwe(n) + ibfwbv(n))
    end do
!$OMP END PARALLEL DO
  end if
!---enthalpy-correction-end

!---wiso-code-begin
    if(lwiso) then
      do n=1, myDim_nod2D+eDim_nod2D
      wiso_flux_oce(n,1)=wiso_flux_oce(n,1)+(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))*1000.0*wiso_smow(1)*(1-30.0/1000.0)
      wiso_flux_oce(n,2)=wiso_flux_oce(n,2)+(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))*1000.0*wiso_smow(2)*(1-240.0/1000.0)
      wiso_flux_oce(n,3)=wiso_flux_oce(n,3)+(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))*1000.0
      end do
    end if
!---wiso-code-end
end subroutine icb2fesom

!_______________________________________________________________________________
! Cap the per-layer iceberg heat flux (ibhf_n) so it cannot cool a layer below
! a safe seawater threshold in a single timestep.  Without this cap, shallow
! shelves with several (or a slow-creeping / grounded) icebergs can drive a
! cell's temperature below freezing over many timesteps and trip the
! temperature blow-up check (cf. Chukchi Sea grounding case).
!
! Deliberate trade-off: this caps heat only, not the paired freshwater flux
! (ibfwb/ibfwl/ibfwe/ibfwbv) -- icb2fesom has already consumed those into
! water_flux/heat_flux by the time this runs (it's called right after
! icb2fesom, since ibhf_n itself is only consumed later, in
! oce_ale_tracer.F90). So when the cap binds, the ocean receives the
! iceberg's full freshwater but pays less than the full latent-heat cost for
! it: a bounded, rare-case energy leak rather than a mass leak. Capping the
! freshwater to match would instead leak mass, since the iceberg's own volume
! update in iceberg_newdimensions has already consumed the uncapped melt
! rates by this point -- closing both gaps at once would require moving the
! cap before that commit, per-iceberg rather than per-cell, which cannot see
! multiple icebergs hitting one node in the same step and would still need
! this aggregate cap as a backstop. The diagnostic below quantifies the
! energy leak so it's monitored rather than silently ignored.
!
! Toggle via namelist /icebergs/ l_cap_ibhf_n (default .true.).
!_______________________________________________________________________________
subroutine cap_ibhf_n(tracers, mesh, partit)
    use o_param
    use g_config
    use MOD_MESH
    use MOD_PARTIT
    use MOD_TRACER
    use iceberg_params
    implicit none
    type(t_tracer), intent(in), target :: tracers
    type(t_mesh)  , intent(in), target :: mesh
    type(t_partit), intent(in), target :: partit
    integer        :: n, nz, nzmin, nzmax, ierr_mpi
    real(kind=WP)  :: t_now, hloc, ibhf_n_min, dheat_cell
    ! Conservation diagnostic: heat [W] discarded by the cap this step,
    ! local-PE-only accumulator (n <= myDim_nod2D guard avoids double-counting
    ! halo nodes), then MPI-summed and logged on PE 0.
    real(kind=WP)  :: heat_disc_loc, heat_disc_glob
    integer        :: cells_capped_loc, cells_capped_glob
    ! Lowest temperature we allow iceberg cooling to drive a cell to.
    ! Slightly above the freezing point of seawater at S~34, so the cap kicks
    ! in before we reach the freezing point and protects the blow-up check.
    real(kind=WP), parameter :: ibhf_n_t_threshold = -1.85_WP
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    if (turn_off_hf .or. .not. l_cap_ibhf_n) return

    heat_disc_loc    = 0.0_WP
    cells_capped_loc = 0

!$OMP PARALLEL DO PRIVATE(n, nz, nzmin, nzmax, t_now, hloc, ibhf_n_min, dheat_cell) &
!$OMP             REDUCTION(+: heat_disc_loc, cells_capped_loc)
    do n = 1, myDim_nod2D + eDim_nod2D
        nzmin = ulevels_nod2D(n)
        nzmax = nlevels_nod2D(n)
        do nz = nzmin, nzmax - 1
            if (ibhf_n(nz, n) >= 0.0_WP) cycle   ! only cooling (negative) needs capping
            hloc = hnode(nz, n)
            if (hloc <= 0.0_WP) cycle
            t_now = tracers%data(1)%values(nz, n)
            if (t_now <= ibhf_n_t_threshold) then
                ! Water already at/below the safe threshold; disallow further
                ! iceberg-driven cooling entirely.
                if (n <= myDim_nod2D) then
                    heat_disc_loc    = heat_disc_loc    + (-ibhf_n(nz, n)) * areasvol(nz, n)
                    cells_capped_loc = cells_capped_loc + 1
                end if
                ibhf_n(nz, n) = 0.0_WP
            else
                ! Cooling allowed up to bringing T down to the threshold in
                ! a single timestep:
                !   dT = ibhf_n * dt / (vcpw * h)   (vcpw = rho_w * cp_w)
                ! Most negative ibhf_n that keeps T_new >= T_thresh:
                ibhf_n_min = (ibhf_n_t_threshold - t_now) * vcpw * hloc / dt
                if (ibhf_n(nz, n) < ibhf_n_min) then
                    if (n <= myDim_nod2D) then
                        dheat_cell       = (ibhf_n_min - ibhf_n(nz, n)) * areasvol(nz, n)
                        heat_disc_loc    = heat_disc_loc    + dheat_cell
                        cells_capped_loc = cells_capped_loc + 1
                    end if
                    ibhf_n(nz, n) = ibhf_n_min
                end if
            end if
        end do
    end do
!$OMP END PARALLEL DO

    call MPI_AllReduce(heat_disc_loc,    heat_disc_glob,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr_mpi)
    call MPI_AllReduce(cells_capped_loc, cells_capped_glob, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_FESOM, ierr_mpi)

    if (mype == 0 .and. cells_capped_glob > 0) then
        write(*, '(A, I7, A, ES12.4, A)') &
              " * cap_ibhf_n: capped ", cells_capped_glob, &
              " cells; discarded heat = ", heat_disc_glob * dt, " J/step"
    end if
end subroutine cap_ibhf_n
end module iceberg_ocean_coupling
