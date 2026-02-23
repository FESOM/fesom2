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
end subroutine


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

    if(i_have_element) then
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
                ! Only include area for nodes that will actually receive flux
                ! (idx_d>0).  Cavity nodes (idx_d=0) are skipped in the flux
                ! loop below; including their area in the denominator would
                ! under-apply the heat flux by ~N_cavity/N_total per element.
                if (idx_d(i) > 0) &
                    tot_area_nods_in_ib_elem = tot_area_nods_in_ib_elem + mesh%area(:,iceberg_node)
            else
                ib_nods_in_ib_elem(i)           = 0
            end if
        end do

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

            if (iceberg_node>0) then
                ibfwbv(iceberg_node) = ibfwbv(iceberg_node) - fwbv_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                ibfwb(iceberg_node) = ibfwb(iceberg_node) - fwb_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                ibfwl(iceberg_node) = ibfwl(iceberg_node) - fwl_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
                ibfwe(iceberg_node) = ibfwe(iceberg_node) - fwe_flux_ib(ib) / tot_area_nods_in_ib_elem(1)
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
                   
                    if( abs(depth_ib) > 0.0 ) then
                        ibhf_n(j,iceberg_node) = ibhf_n(j,iceberg_node) & 
                                                    - (hfbv_flux_ib(ib,j)+hfl_flux_ib(ib,j)) & 
                                                    / tot_area_nods_in_ib_elem(j)
                    end if
                end do
                
                if( idx_d(i) > 1 ) then
                    ibhf_n(idx_d(i),iceberg_node) = ibhf_n(idx_d(i),iceberg_node) - 0.5 * hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(idx_d(i))
                    ibhf_n(idx_d(i)-1,iceberg_node) = ibhf_n(idx_d(i)-1,iceberg_node) - 0.5 * hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(idx_d(i)-1)
                else
                    ibhf_n(idx_d(i),iceberg_node) = ibhf_n(idx_d(i),iceberg_node) - hfb_flux_ib(ib) / tot_area_nods_in_ib_elem(idx_d(i))
                end if
                
                if( height_ib_single .ne. 0.0 ) then
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
end module
