module iceberg_params
use MOD_PARTIT
USE MOD_MESH
use g_config, only: use_cavity, use_cavityonelem
use, intrinsic :: iso_fortran_env, only: real64
USE o_PARAM, only: WP

implicit none
save
  !integer,parameter :: ib_num               ! realistic dataset comprising 6912 icebergs
  real,dimension(:), allocatable:: calving_day  !271.0 	!28.0: September 29 for restart in 1 SEP 97 ! 271.0: September 29 for year 1997
  
  !days since beginning of year (since FESOM first started);
  !2.5 starts an iceberg at 3rd day, 12.00, with respect to model start;
  !For restarts, calving_day is subtracted by the number of modelled days (see iceberg_out subroutine)

  ! ============= REALISTIC INITIAL ICEBERG DISTRIBUTION (SEP/OCT 1997 SNAPSHOT) IN NEAR-COASTAL STRIP AROUND ANTARCTICA  =================== !
  !
  !        from the article 'Near-coastal circum-Antarctic iceberg size distributions determined from Synthetic Aperture Radar images'
  !        by Wesche and Dierking (2014).
  real,dimension(:), allocatable:: height_ib

  !read from file in init_icebergs
  real,dimension(:), allocatable:: length_ib
  real,dimension(:), allocatable:: width_ib
  real,dimension(:), allocatable:: lon_deg
  real,dimension(:), allocatable:: lat_deg
  !in case (l_initial .AND. l_iniuser) = .true. ;
  !initial zonal velocity (positive is to the east):
  real,dimension(:), allocatable:: ini_u 

  !initial meridional velocity (positive is to the north):
  real,dimension(:), allocatable:: ini_v

  ! ================================================================================================================= !
  ! ========= Lichey & Hellmer values ========= !
  real,dimension(:), allocatable:: Co
  real,dimension(:), allocatable:: Ca
  real,dimension(:), allocatable:: Ci
  real,dimension(:), allocatable:: Cdo_skin     ! !Cd_oce_ice = 5.0e-3
  real,dimension(:), allocatable:: Cda_skin     ! !similar to Keghouche (2009)
  ! =========================================== !
  logical,dimension(:), allocatable:: l_wave    ! (use wave radiation force?)
  ! =========================================== !
  real,dimension(:), allocatable:: conc_sill
  real,dimension(:), allocatable:: P_sill
  logical,dimension(:), allocatable:: l_freeze  ! (use freezing parametrization?)
  ! =========================================== !
  logical       :: l_melt = .true.              ! (use melting parametrization?)
  logical       :: l_weeksmellor = .true.       ! (use weeks & mellor stability criterion?)
  logical       :: l_allowgrounding = .true.    ! (are icebergs allowed to ground?)
  real,dimension(:), allocatable:: draft_scale  ! (account for irregularities of draft
  ! =========================================== !
  logical       :: l_tides = .false.            ! (simulate sensitivity to tides? !!check for HLRN-III!!)
  ! =========================================== !
  logical       :: l_initial = .true.           ! (use initial iceberg velocities?)
  logical       :: l_iniuser = .false.          ! (prescribe init. velo or let model decide?)
  ! =========================================== !
  real          :: smallestvol_icb = 10.0       ! (smallest iceberg volume in m^3?)	
  real          :: maxspeed_icb = 3.0           ! (cap iceberg speed at ?? m/s?) !security value !not used
  ! ========= For sensitivity studies ========= !
  real,dimension(:), allocatable:: coriolis_scale   ! (scale the body forces, Coriolis and
  real,dimension(:), allocatable:: surfslop_scale   !  surface slope, by those factors:
  ! =========================================== !
  real,dimension(:), allocatable:: rho_icb          !Silva et al., Martin
  real,dimension(:), allocatable:: rho_h2o
  real,dimension(:), allocatable:: rho_air
  real,dimension(:), allocatable:: rho_ice          !910 RT, 945.0 bei Lichey, aus Lemke (1993)

  !python version
  character(100):: IcebergRestartPath='iceberg.restart'
  character(100):: IcebergRestartPath_ISM='iceberg.restart.ISM'
  character(100):: num_non_melted_icb_file='num_non_melted_icb_file'
  character(100):: file_icb_netcdf='buoys_track.nc' !output file of buoys/icebergs
  character(100):: buoys_xlon_file='icb_longitude.dat'     !buoy position in deg
  character(100):: buoys_ylat_file='icb_latitude.dat'     !buoy position in deg
  character(100):: length_icb_file='icb_length.dat' !iceberg length [m]
  character(100):: width_icb_file='icb_length.dat' !iceberg width [m]
  character(100):: height_icb_file='icb_height.dat' !iceberg height [m]
  character(100):: scaling_file='icb_scaling.dat' !scaling factor
  
  !===== OUTPUT RELATED SETTINGS  =====
  integer :: icb_outfreq           ! 180; for FESOM_dt=2min this is 6 hourly output !120; for FESOM_dt=3min this is 6 hourly output
  logical :: l_geo_out = .true.         ! output in unrotated (.true.) or rotated coordinates
  logical :: ascii_out = .false.        ! old ascii output (slow, more detailed); false: faster nc output
  !===== NUMERICS (DONT HAVE TO BE CHANGED) =====
  logical :: l_semiimplicit = .true.    !false: adams-bashforth for coriolis
  real :: semiimplicit_coeff = 1.0      !1. fully implicit, 0.5 no damping
  real :: AB_coeff = 1.53               !1.5 original AB (amplifying), 1.6 stabilized
  !===== NOTHING MUST BE CHANGED BELOW THIS LINE =====
  real,dimension(:), allocatable:: T_ave
  real,dimension(:), allocatable:: u_ib, v_ib
  integer,dimension(:), allocatable:: iceberg_elem
  logical,dimension(:), allocatable:: find_iceberg_elem
  real,dimension(:), allocatable:: f_u_ib_old, f_v_ib_old
  real,dimension(:), allocatable:: bvl_mean, lvlv_mean, lvle_mean, lvlb_mean !averaged volume losses
  !real,dimension(:), allocatable:: fw_flux_ib, hfb_flux_ib
  real,dimension(:), allocatable:: fwe_flux_ib, fwl_flux_ib, fwb_flux_ib, fwbv_flux_ib
  real,dimension(:), allocatable:: hfe_flux_ib, hfb_flux_ib, lhfb_flux_ib
  real,dimension(:,:), allocatable:: hfl_flux_ib, hfbv_flux_ib
  
  !===== FRESHWATER AND HEAT ARRAYS ON FESOM GRID =====
  real,dimension(:), allocatable:: ibhf    !icb heat flux into ocean 
  real(kind=WP),dimension(:,:), allocatable:: ibhf_n    !icb heat flux into ocean 
  real,dimension(:), allocatable:: ibfwb   !freshwater flux into ocean from basal melting
  real,dimension(:), allocatable:: ibfwbv   !freshwater flux into ocean from basal melting
  real,dimension(:), allocatable:: ibfwl   !freshwater flux into ocean from lateral melting
  real,dimension(:), allocatable:: ibfwe   !freshwater flux into ocean from erosion
  integer,dimension(:), allocatable:: scaling   !scaling factor

  logical,dimension(:), allocatable::   melted  !1 if iceberg melted, 0 otherwise
  logical,dimension(:), allocatable::   grounded    !1 if iceberg grounded, 0 otherwise
  integer   :: num_non_melted_icb = 0 !1 if iceberg melted, 0 otherwise
  !for communication
  real,dimension(:), allocatable:: arr_block
  integer,dimension(:), allocatable:: elem_block
  integer,dimension(:), allocatable:: pe_block
  real(real64), dimension(:), allocatable:: elem_area_glob
  real,dimension(:), allocatable:: vl_block

  !array for output in netcdf
  real,dimension(:,:), allocatable:: buoy_props
  integer:: save_count_buoys
  real:: prev_sec_in_year
!****************************************************************************************************************************
 contains
 ! true if all nodes of the element are either "real" model boundary nodes or shelf nodes
 logical function reject_elem(mesh, partit, elem)
 implicit none
 integer, intent(in) :: elem
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

if (use_cavity) then
! kh 09.08.21 change index_nod2d -> bc_index_nod2d?
 if (.not. use_cavityonelem) then
   reject_elem = all( (mesh%cavity_depth(mesh%elem2D_nodes(:,elem))/=0.0) .OR. (mesh%bc_index_nod2D(mesh%elem2D_nodes(:,elem))==0.0) )
 !else
 end if
else
 reject_elem = all( (mesh%bc_index_nod2D(mesh%elem2D_nodes(:,elem))==0.0) )
endif
 end function reject_elem
 
 ! gives number of "coastal" nodes in cavity setup, i.e. number of nodes that are
 ! either "real" model boundary nodes or shelf nodes
 integer function coastal_nodes(mesh, partit, elem)
 implicit none
 integer, intent(in) :: elem
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

if (use_cavity) then
! kh 09.08.21 change index_nod2d -> bc_index_nod2d?
 if (.not. use_cavityonelem) then
   coastal_nodes = count( (mesh%cavity_depth(mesh%elem2D_nodes(:,elem))/=0.0) .OR. (mesh%bc_index_nod2D(mesh%elem2D_nodes(:,elem))==0.0) )
 !else
 end if
else 
 coastal_nodes = count( (mesh%bc_index_nod2D(mesh%elem2D_nodes(:,elem))==0.0) )
endif
 end function coastal_nodes

end module iceberg_params
