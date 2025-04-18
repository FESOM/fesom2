module allocate_icb_interface
  interface
     subroutine allocate_icb(partit, mesh)
        use iceberg_params
        use g_config
        use g_comm
        use g_comm_auto
        use o_param
        use MOD_PARTIT
        use MOD_MESH
        implicit none
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in)   , target :: mesh
     end subroutine allocate_icb
  end interface
end module allocate_icb_interface

subroutine allocate_icb(partit, mesh)
  use iceberg_params
  use g_config
  use g_comm
  use g_comm_auto
  use o_param
  use MOD_PARTIT
  use MOD_MESH

  integer       :: n2
type(t_partit), intent(inout), target :: partit
type(t_mesh), intent(in), target :: mesh
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  n2=myDim_nod2D+eDim_nod2D
  icb_outfreq = step_per_day / steps_per_ib_step

  allocate(ibhf(n2), ibfwb(n2), ibfwl(n2), ibfwe(n2), ibfwbv(n2))
  ibhf      = 0.0
  ibfwb     = 0.0
  ibfwl     = 0.0
  ibfwe     = 0.0
  ibfwbv    = 0.0
  allocate(ibhf_n(mesh%nl, n2))
  ibhf_n    = 0.0_WP

  allocate(wave_erosion_potential(elem2D), linit_wave_erosion_pot(elem2D))
  wave_erosion_potential = 0.0
  linit_wave_erosion_pot = .true.

  allocate(calving_day(ib_num))
  calving_day = 1   !28.0: September 29 for restart in 1 SEP 97 ! 271.0: September 29 for year 1997
  allocate(height_ib(ib_num))
  height_ib = 1.0 ! 250.0 ! 360.0
  allocate(length_ib(ib_num))
  length_ib = 1.0
  allocate(width_ib(ib_num))
  width_ib = 1.0
  allocate(lon_deg(ib_num))
  lon_deg = 0.0
  allocate(lat_deg(ib_num))
  lat_deg = 0.0
  allocate(ini_u(ib_num))
  ini_u = -0.2
  allocate(ini_v(ib_num))
  ini_v = 0.0
  allocate(Co(ib_num))
  Co = 0.85 
  allocate(Ca(ib_num))
  Ca = 0.4
  allocate(Ci(ib_num))
  Ci = 1.0
  allocate(Cdo_skin(ib_num))
  Cdo_skin = 5.0e-3    ! !Cd_oce_ice = 5.0e-3
  allocate(Cda_skin(ib_num))
  Cda_skin = 2.5e-3    ! !similar to Keghouche (2009)
  allocate(l_wave(ib_num))
  l_wave = .false.        ! (use wave radiation force?)
  allocate(conc_sill(ib_num))
  conc_sill =0.90
  allocate(P_sill(ib_num))
  P_sill = 10000.0
  allocate(l_freeze(ib_num))
  l_freeze = .true.       ! (use freezing parametrization?)
  allocate(draft_scale(ib_num))
  draft_scale = 1.0    ! (account for irregularities of draft
  allocate(coriolis_scale(ib_num))
  coriolis_scale = 1.0   ! (scale the body forces, Coriolis and
  allocate(surfslop_scale(ib_num))
  surfslop_scale = 1.0   !  surface slope, by those factors:
  allocate(rho_icb(ib_num))
  rho_icb = 850.0        !Silva et al., Martin
  allocate(rho_h2o(ib_num))
  rho_h2o = 1027.5
  allocate(rho_air(ib_num))
  rho_air = 1.293
  allocate(rho_ice(ib_num))
  rho_ice = 910.0        !910 RT, 945.0 bei Lichey, aus Lemke (1993)
  allocate(u_ib(ib_num))
  allocate(v_ib(ib_num))
  allocate(iceberg_elem(ib_num))
  allocate(find_iceberg_elem(ib_num))
  find_iceberg_elem = .true.
  allocate(f_u_ib_old(ib_num))
  allocate(f_v_ib_old(ib_num))
  allocate(bvl_mean(ib_num))
  bvl_mean = 0.0
  allocate(lvlv_mean(ib_num))
  lvlv_mean = 0.0
  allocate(lvle_mean(ib_num))
  lvle_mean = 0.0
  allocate(lvlb_mean(ib_num))
  lvlb_mean = 0.0 !averaged volume losses
  allocate(fwe_flux_ib(ib_num))
  allocate(fwl_flux_ib(ib_num))
  allocate(fwb_flux_ib(ib_num))
  allocate(fwbv_flux_ib(ib_num))
  allocate(hfe_flux_ib(ib_num))
  allocate(hfl_flux_ib(ib_num,mesh%nl))
  allocate(hfb_flux_ib(ib_num))
  allocate(hfbv_flux_ib(ib_num,mesh%nl))
  allocate(lhfb_flux_ib(ib_num))
  fwe_flux_ib = 0.0
  fwl_flux_ib = 0.0
  fwb_flux_ib = 0.0
  fwbv_flux_ib = 0.0
  hfe_flux_ib = 0.0
  hfl_flux_ib = 0.0
  hfb_flux_ib = 0.0
  hfbv_flux_ib = 0.0
  lhfb_flux_ib = 0.0
  allocate(arr_block(15*ib_num))
  allocate(elem_block(ib_num))
  allocate(pe_block(ib_num))

  allocate(elem_area_glob(elem2D))
  elem_area_glob=0.0
  call gather_elem(elem_area(1:myDim_elem2D), elem_area_glob, partit)
  call MPI_Bcast(elem_area_glob, elem2D, MPI_DOUBLE, 0, MPI_COMM_FESOM, MPIERR)

  allocate(vl_block(4*ib_num))
  allocate(buoy_props(ib_num,13))
  buoy_props = 0.0
  allocate(melted(ib_num))
  melted = .false.
  allocate(grounded(ib_num))
  grounded = .false.
  allocate(scaling(ib_num))
  scaling = 1
end subroutine allocate_icb
