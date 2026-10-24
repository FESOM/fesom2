module icb_allocate_module
    USE iceberg_params
    USE g_config
    USE g_comm
    USE g_comm_auto
    USE o_param
    USE MOD_PARTIT
    USE MOD_MESH

    implicit none

    private
    public :: allocate_icb

contains

subroutine allocate_icb(partit, mesh)

  integer       :: n2
type(t_partit), intent(inout), target :: partit
type(t_mesh), intent(in), target :: mesh
  n2=partit%myDim_nod2D+partit%eDim_nod2D
  icb_outfreq = step_per_day / steps_per_ib_step

  allocate(ibhf(n2), ibfwb(n2), ibfwl(n2), ibfwe(n2), ibfwbv(n2))
  ibhf      = 0.0
  ibfwb     = 0.0
  ibfwl     = 0.0
  ibfwe     = 0.0
  ibfwbv    = 0.0
  allocate(ibiron(n2))
  ibiron    = 0.0
  allocate(ibhf_n(mesh%nl, n2))
  ibhf_n    = 0.0_WP

  allocate(wave_erosion_potential(mesh%elem2D), linit_wave_erosion_pot(mesh%elem2D))
  wave_erosion_potential = 0.0
  linit_wave_erosion_pot = .true.

  allocate(calving_day(ib_num))
  calving_day = 0   !28.0: September 29 for restart in 1 SEP 97 ! 271.0: September 29 for year 1997
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
  allocate(iron_conc_ib(ib_num))
  allocate(iron_flux_ib(ib_num))
  allocate(hfe_flux_ib(ib_num))
  allocate(hfl_flux_ib(ib_num,mesh%nl))
  allocate(hfb_flux_ib(ib_num))
  allocate(hfbv_flux_ib(ib_num,mesh%nl))
  allocate(lhfb_flux_ib(ib_num))
  fwe_flux_ib = 0.0
  fwl_flux_ib = 0.0
  fwb_flux_ib = 0.0
  fwbv_flux_ib = 0.0
  iron_conc_ib = icb_iron_const   ! overwritten by icb_iron.dat / the iron restart
  iron_flux_ib = 0.0
  hfe_flux_ib = 0.0
  hfl_flux_ib = 0.0
  hfb_flux_ib = 0.0
  hfbv_flux_ib = 0.0
  lhfb_flux_ib = 0.0
  allocate(arr_block(16*ib_num))
  allocate(elem_block(ib_num))
  allocate(pe_block(ib_num))

  allocate(elem_area_glob(mesh%elem2D))
  elem_area_glob=0.0
  call gather_elem(mesh%elem_area(1:partit%myDim_elem2D), elem_area_glob, partit)
  call MPI_Bcast(elem_area_glob, mesh%elem2D, MPI_WP, 0, partit%MPI_COMM_FESOM, partit%MPIERR)
  call init_halo_elems(partit, mesh)

  allocate(vl_block(4*ib_num))
  allocate(buoy_props(ib_num,14))
  buoy_props = 0.0
  allocate(melted(ib_num))
  melted = .false.
  allocate(grounded(ib_num))
  grounded = .false.
  allocate(scaling(ib_num))
  scaling = 1
end subroutine allocate_icb

! Corners and cavity rejection of every halo element (eDim and eXDim), filled
! once through the full element halo exchange: the owners compute them for
! their own elements, whose nodes they hold.
subroutine init_halo_elems(partit, mesh)

type(t_partit), intent(inout), target :: partit
type(t_mesh), intent(in), target :: mesh
  real(kind=WP), allocatable :: aux(:,:)
  integer       :: n_my, n_all, n, i, k
  integer       :: elnodes(3)

  n_my  = partit%myDim_elem2D
  n_all = n_my + partit%eDim_elem2D + partit%eXDim_elem2D
  allocate(aux(n_all, 7))
  aux = 0.0_WP
  do n = 1, n_my
     elnodes = mesh%elem2D_nodes(:,n)
     do i = 1, 3
        aux(n, 2*i-1) = mesh%coord_nod2D(1, elnodes(i))
        aux(n, 2*i)   = mesh%coord_nod2D(2, elnodes(i))
     end do
     if (use_cavity) then
        if (all((mesh%cavity_depth(elnodes)/=0.0) .or. (mesh%bc_index_nod2D(elnodes)==0.0))) aux(n, 7) = 1.0_WP
     end if
  end do
  do k = 1, 7
     call exchange_elem(aux(:,k), partit)
  end do

  allocate(halo_elem_xy(2, 3, n_my+1:n_all), halo_elem_reject(n_my+1:n_all))
  do n = n_my+1, n_all
     halo_elem_xy(1,:,n) = aux(n, 1:5:2)
     halo_elem_xy(2,:,n) = aux(n, 2:6:2)
     halo_elem_reject(n) = aux(n, 7) > 0.5_WP
  end do
  deallocate(aux)
end subroutine init_halo_elems

end module icb_allocate_module
