module iceberg_ocean_coupling
 USE MOD_MESH
 use MOD_PARTIT
 use MOD_ICE
 USE MOD_DYN
 use iceberg_params

 public :: reset_ib_fluxes
 public :: prepare_icb2fesom
 public :: icb2fesom

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
    ibhf    = 0.0
    ibhf_n  = 0.0_WP
    ! Zero per-iceberg per-level heat flux arrays.  These are only initialised
    ! to 0 at allocation time; without an explicit reset here, stale values
    ! persist whenever ib_n_lvls shrinks (e.g. when an iceberg drifts into a
    ! cavity element and ib_n_lvls drops to 0).  prepare_icb2fesom may then
    ! divide those stale values by a zero node area -> Inf/NaN in ibhf_n.
    hfbv_flux_ib = 0.0
    hfl_flux_ib  = 0.0
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

    if(i_have_element) then                                     ! PE has iceberg
        dz = 0.0                                                ! vertical distance over which heat flux is applied 
        allocate(tot_area_nods_in_ib_elem(mesh%nl))

        num_ib_nods_in_ib_elem=0                                ! number of nodes in this element ???
        tot_area_nods_in_ib_elem=0.0                            ! total area associated to nodes of iceberg element
        idx_d = 0                                               ! index of level directly below or at iceberg base

        ! loop over all three nodes of element
        do i=1,3
            
            ! assign node to iceberg_node
            iceberg_node=elem2d_nodes(i,localelement)

            ! loop over all levels of iceberg_node
            do k=1, nlevels_nod2D(iceberg_node)
                
                idx_d(i) = k                                    ! idx_d holds current level index until ...
                                                                ! ... end of iceberg or bottom topography is reached
                lev_up  = mesh%zbar_3d_n(k, iceberg_node)       ! upper level

                if( k==nlevels_nod2D(iceberg_node) ) then       ! if lower most level is reached ...
                    lev_low = mesh%zbar_n_bot(iceberg_node)     ! ... lower level is equal to bottom topography
                else
                    lev_low = mesh%zbar_3d_n(k+1, iceberg_node) ! ... otherwise, lower level is one level below upper one
                end if

                if( abs(lev_low)==abs(lev_up) ) then            ! if upper level is equal to lower level - when should this happen?
                    idx_d(i) = idx_d(i) - 1                     ! level index is set back by one and exit loop
                    exit
                else if( abs(lev_low)>=abs(depth_ib) ) then     ! if lower level is below iceberg ...
                                                                ! ... i.e. end of iceberg is reached, exit loop
                    exit
                else
                    cycle
                end if
            end do
            ! at the end of the loop, idx_n holds the index of the level ...
            ! ... directly below the iceberg or at iceberg base
            ! -------------------------------------------------------------------

            if (iceberg_node<=mydim_nod2d) then
                ib_nods_in_ib_elem(i)           = iceberg_node
                num_ib_nods_in_ib_elem          = num_ib_nods_in_ib_elem + 1
                ! Only include area for nodes that will actually receive flux
                ! (idx_d>0).  Cavity nodes (idx_d=0) are skipped in the flux
                ! loop below; including their area in the denominator would
                ! under-apply the heat flux by ~N_cavity/N_total per element.
                if (idx_d(i) > 0) &
                    tot_area_nods_in_ib_elem = tot_area_nods_in_ib_elem + mesh%area(:,iceberg_node)
            else
                ib_nods_in_ib_elem(i)           = 0                             ! ... otherwise, add zero to list ib_nods_in_ib_elem
            end if
        end do
        ! at the then of the loop, ib_nods_in_ib_elem holds the nodes of the iceberg elemen - why so complicated?
        ! ... num_ib_nods_in_ib_elem holds the number of nodes assigned to PE
        ! ... tot_area_nods_in_ib_elem holds the total area of the (up to three) nodes assigned to iceberg_elem
        ! -------------------------------------------------------------------

        ! loop over (up to three) nodes of iceberg_elem
        do i=1, 3
            iceberg_node=ib_nods_in_ib_elem(i)
            ! Skip land nodes (ulevels==0) and cavity nodes (ulevels>1).
            ! Surface FW and heat fluxes must only go to open-ocean nodes
            ! (ulevels==1).  In oce_mesh.F90, mesh%area(k,n) is only filled
            ! for k >= ulevels(elem), so mesh%area(1,n)==0 for cavity nodes.
            ! Dividing by tot_area_nods_in_ib_elem(1)==0 gives NaN that
            ! propagates into integrate_nod -> "total iceberg fw flux: NaN".
            if (iceberg_node <= 0 .or. ulevels_nod2d(iceberg_node) /= 1) cycle
            ! Belt-and-suspenders: skip if no surface area was accumulated
            if (tot_area_nods_in_ib_elem(1) <= 0.0) cycle

            ! if iceberg_node in PE, convert freshwater flux to flux density by dividing with tot_area_nods_in_ib_elem ...
            ! ... of upper most level. The total iceberg flux is distributed among up to three nodes and only applied to the ...
            ! ... surface layer.
            if (iceberg_node>0) then
                ibfwbv(iceberg_node) = ibfwbv(iceberg_node) - fwbv_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                ibfwb(iceberg_node) = ibfwb(iceberg_node) - fwb_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                ibfwl(iceberg_node) = ibfwl(iceberg_node) - fwl_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                ibfwe(iceberg_node) = ibfwe(iceberg_node) - fwe_flux_ib(ib) / tot_area_nods_in_ib_elem(1)

                ! Guard against idx_d=0 (can happen if nlevels_nod2D is 0 or loop didn't execute)
                if (idx_d(i) <= 0) cycle

                do j=1,idx_d(i)
                    lev_up  = mesh%zbar_3d_n(j, iceberg_node)           ! upper level
                    if( j==nlevels_nod2D(iceberg_node) ) then           ! if bottom level is reached ...
                        lev_low = mesh%zbar_n_bot(iceberg_node)         ! ... lower level is set to bottom topography
                    else
                        lev_low = mesh%zbar_3d_n(j+1, iceberg_node)     ! ... otherwise, lower level is one level below upper one
                    end if
                    dz = abs( lev_low - lev_up )                        ! distance between lower and upper level
                    
                    ! check if lower level is below iceberg base and upper level is above iceberg base ...
                    if( (abs(lev_low)>=abs(depth_ib)) .and. (abs(lev_up)<abs(depth_ib)) ) then 
                        dz = abs(abs(lev_up) - abs(depth_ib))           ! ... if so, distance dz is calculated as difference ...
                                                                        ! ... between upper level and iceberg base
                    end if              
                   
                    ! Guard against zero node area at level j (can occur in
                    ! mixed cavity/open-ocean elements where not all nodes
                    ! reach depth j, leaving tot_area_nods_in_ib_elem(j)==0).
                    if( abs(depth_ib) > 0.0 .and. tot_area_nods_in_ib_elem(j) > 0.0 ) then
                        ibhf_n(j,iceberg_node) = ibhf_n(j,iceberg_node) & 
                                                    - (hfbv_flux_ib(ib,j)+hfl_flux_ib(ib,j)) & 
                                                    / tot_area_nods_in_ib_elem(j)
                    end if
                end do
                
                ! if idx_d is not only upper most level ...
                if( idx_d(i) > 1 ) then
                    if (tot_area_nods_in_ib_elem(idx_d(i)) > 0.0) &
                        ibhf_n(idx_d(i),iceberg_node) = ibhf_n(idx_d(i),iceberg_node) - 0.5 * hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(idx_d(i))
                    if (tot_area_nods_in_ib_elem(idx_d(i)-1) > 0.0) &
                        ibhf_n(idx_d(i)-1,iceberg_node) = ibhf_n(idx_d(i)-1,iceberg_node) - 0.5 * hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(idx_d(i)-1)
                else
                    
                    ! ... otherwise (should not happen), assign 100% of basal heat flux to level idx_n, i.e. the surface level
                    ibhf_n(idx_d(i),iceberg_node) = ibhf_n(idx_d(i),iceberg_node) - hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(idx_d(i))
                end if
               
                ! if iceberg not yet melted away ...
                if( height_ib_single .ne. 0.0 ) then
                    
                    ! ... assign heat flux due to wave erosion to surface level
                    ibhf_n(1,iceberg_node) = ibhf_n(1,iceberg_node) - hfe_flux_ib(ib) & 
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
                ! LA add heat fluxes here...
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
end module iceberg_ocean_coupling
