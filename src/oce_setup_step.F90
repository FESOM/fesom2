module oce_initial_state_interface
  interface
    subroutine oce_initial_state(tracers, partit, mesh)
      USE MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      use mod_tracer
      type(t_mesh),   intent(in)  ,  target :: mesh
      type(t_partit), intent(inout), target :: partit
      type(t_tracer), intent(inout), target :: tracers
    end subroutine
  end interface
end module

module tracer_init_interface
  interface
    subroutine tracer_init(tracers, partit, mesh)
      USE MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      use mod_tracer
      type(t_mesh),   intent(in),    target :: mesh
      type(t_partit), intent(inout), target :: partit
      type(t_tracer), intent(inout), target :: tracers
    end subroutine
  end interface
end module

module dynamics_init_interface
  interface
    subroutine dynamics_init(dynamics, partit, mesh)
      USE MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      use MOD_DYN
      type(t_mesh)  , intent(in)   , target :: mesh
      type(t_partit), intent(inout), target :: partit
      type(t_dyn)   , intent(inout), target :: dynamics
    end subroutine
  end interface
end module

module ocean_setup_interface
  interface
    subroutine ocean_setup(dynamics, tracers, partit, mesh)
      USE MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      use mod_tracer
      use MOD_DYN
      type(t_mesh),   intent(in),    target :: mesh
      type(t_partit), intent(inout), target :: partit
      type(t_tracer), intent(inout), target :: tracers
      type(t_dyn), intent(inout), target :: dynamics
    end subroutine
  end interface
end module
module before_oce_step_interface
  interface
    subroutine before_oce_step(dynamics, tracers, partit, mesh)
      USE MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      use mod_tracer
      use MOD_DYN
      type(t_mesh),   intent(in),    target :: mesh
      type(t_partit), intent(inout), target :: partit
      type(t_tracer), intent(inout), target :: tracers
      type(t_dyn), intent(inout), target :: dynamics
    end subroutine
  end interface
end module
!_______________________________________________________________________________
subroutine ocean_setup(dynamics, tracers, partit, mesh)
USE MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
USE MOD_TRACER
USE MOD_DYN
USE o_PARAM
USE o_ARRAYS
USE g_config
USE g_forcing_param, only: use_virt_salt
use g_cvmix_tke
use g_cvmix_idemix
use g_cvmix_pp
use g_cvmix_kpp
use g_cvmix_tidal
use Toy_Channel_Soufflet
use oce_initial_state_interface
use oce_adv_tra_fct_interfaces
IMPLICIT NONE
type(t_mesh),   intent(inout), target :: mesh
type(t_partit), intent(inout), target :: partit
type(t_tracer), intent(inout), target :: tracers
type(t_dyn), intent(inout), target :: dynamics
integer                               :: n
    !___setup virt_salt_flux____________________________________________________
    ! if the ale thinkness remain unchanged (like in 'linfs' case) the vitrual 
    ! salinity flux need to be used
    ! otherwise we set the reference salinity to zero
    if ( .not. trim(which_ALE)=='linfs') then
        use_virt_salt=.false.
        ! this will force the virtual saltinity flux to be zero
        ref_sss_local=.false.
        ref_sss=0._WP
        is_nonlinfs = 1.0_WP
    else
        use_virt_salt=.true.
        is_nonlinfs = 0.0_WP
    end if
   
    !___________________________________________________________________________
    ! initialize arrays for ALE
    if (partit%mype==0) then
       write(*,*) '____________________________________________________________'
       write(*,*) ' --> initialise ALE arrays + sparse SSH stiff matrix'
       write(*,*)
    end if
    call init_ale(partit, mesh)
    call init_stiff_mat_ale(partit, mesh) !!PS test  
    
    !___________________________________________________________________________
    ! initialize arrays from cvmix library for CVMIX_KPP, CVMIX_PP, CVMIX_TKE,
    ! CVMIX_IDEMIX and CVMIX_TIDAL
    ! here translate mix_scheme string into integer --> for later usage only 
    ! integer comparison is required
    select case (trim(mix_scheme))
        case ('KPP'                   ) ; mix_scheme_nmb = 1
        case ('PP'                    ) ; mix_scheme_nmb = 2
        case ('cvmix_KPP'             ) ; mix_scheme_nmb = 3
        case ('cvmix_PP'              ) ; mix_scheme_nmb = 4
        case ('cvmix_TKE'             ) ; mix_scheme_nmb = 5
        case ('cvmix_IDEMIX'          ) ; mix_scheme_nmb = 6
        case ('cvmix_TIDAL'           ) ; mix_scheme_nmb = 7 
        case ('KPP+cvmix_TIDAL'       ) ; mix_scheme_nmb = 17
        case ('PP+cvmix_TIDAL'        ) ; mix_scheme_nmb = 27
        case ('cvmix_KPP+cvmix_TIDAL' ) ; mix_scheme_nmb = 37
        case ('cvmix_PP+cvmix_TIDAL'  ) ; mix_scheme_nmb = 47
        case ('cvmix_TKE+cvmix_IDEMIX') ; mix_scheme_nmb = 56
        case default 
            stop "!not existing mixing scheme!"
            call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    end select

    ! initialise fesom1.4 like KPP
    if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then
        call oce_mixing_kpp_init(partit, mesh)
    ! initialise fesom1.4 like PP
    elseif (mix_scheme_nmb==2 .or. mix_scheme_nmb==27) then
    
    ! initialise cvmix_KPP
    elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then
        call init_cvmix_kpp(partit, mesh)
        
    ! initialise cvmix_PP    
    elseif (mix_scheme_nmb==4 .or. mix_scheme_nmb==47) then
        call init_cvmix_pp(partit, mesh)
        
    ! initialise cvmix_TKE    
    elseif (mix_scheme_nmb==5 .or. mix_scheme_nmb==56) then
        call init_cvmix_tke(partit, mesh)
        
    endif
  
    ! initialise additional mixing cvmix_IDEMIX --> only in combination with 
    ! cvmix_TKE+cvmix_IDEMIX or stand alone for debbuging as cvmix_TKE
    if     (mod(mix_scheme_nmb,10)==6) then
        call init_cvmix_idemix(partit, mesh)
        
    ! initialise additional mixing cvmix_TIDAL --> only in combination with 
    ! KPP+cvmix_TIDAL, PP+cvmix_TIDAL, cvmix_KPP+cvmix_TIDAL, cvmix_PP+cvmix_TIDAL 
    ! or stand alone for debbuging as cvmix_TIDAL   
    elseif (mod(mix_scheme_nmb,10)==7) then
        call init_cvmix_tidal(partit, mesh)
    end if         
    
    !___________________________________________________________________________
    ! set use_density_ref .true. when cavity is used and initialse cavity boundary 
    ! line for the extrapolation of the initialisation
    if (use_cavity .and. .not. use_density_ref) use_density_ref=.true.
    
    ! compute for all cavity points (ulevels_nod2D>1), which is the closest
    ! cavity line point to that point --> use their coordinates and depth -->
    ! use for extrapolation of init state under cavity
    if (use_cavity) call compute_nrst_pnt2cavline(partit, mesh)
      
    if (use_density_ref) call init_ref_density(partit, mesh)
    
    
    !___________________________________________________________________________
    if(partit%mype==0) write(*,*) 'Arrays are set'
        
    !if(open_boundary) call set_open_boundary   !TODO
    
    call oce_adv_tra_fct_init(tracers%work, partit, mesh)
    call muscl_adv_init(tracers%work, partit, mesh) !!PS test
    !=====================
    ! Initialize fields
    ! A user-defined routine has to be called here!
    !=====================
    if (toy_ocean) then  
       SELECT CASE (TRIM(which_toy))
         CASE ("soufflet") !forcing update for soufflet testcase
           if (mod(mstep, soufflet_forc_update)==0) then
              call initial_state_soufflet(dynamics, tracers, partit, mesh)
              call compute_zonal_mean_ini(partit, mesh)  
              call compute_zonal_mean(dynamics, tracers, partit, mesh)
           end if
       END SELECT
    else
       call oce_initial_state(tracers, partit, mesh)   ! Use it if not running tests
    end if

    if (.not.r_restart) then
       do n=1, tracers%num_tracers
          tracers%data(n)%valuesAB=tracers%data(n)%values
       end do
    end if
    
    !___________________________________________________________________________
    ! first time fill up array for hnode & helem
    if (partit%mype==0) then
        write(*,*) '____________________________________________________________'
        write(*,*) ' --> call init_thickness_ale'
        write(*,*)
    end if
    call init_thickness_ale(partit, mesh)
    
    !___________________________________________________________________________
    if(partit%mype==0) write(*,*) 'Initial state'
    if (w_split .and. partit%mype==0) then
        write(*,*) '******************************************************************************'
        write(*,*) 'vertical velocity will be split onto explicit and implicit constitutes;'
        write(*,*) 'maximum allowed CDF on explicit W is set to: ', w_max_cfl
        write(*,*) '******************************************************************************'
    end if
end subroutine ocean_setup
!_______________________________________________________________________________
SUBROUTINE tracer_init(tracers, partit, mesh)
USE MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
USE MOD_TRACER
USE DIAGNOSTICS, only: ldiag_DVD
USE g_ic3d
IMPLICIT NONE
integer        :: elem_size, node_size
integer, save  :: nm_unit  = 104       ! unit to open namelist file, skip 100-102 for cray
integer        :: iost
integer        :: n

integer        :: num_tracers
logical        :: i_vert_diff, smooth_bh_tra
real(kind=WP)  :: gamma0_tra, gamma1_tra, gamma2_tra

type(t_mesh),   intent(in) ,   target               :: mesh
type(t_partit), intent(inout), target               :: partit
type(t_tracer), intent(inout), target               :: tracers
type(nml_tracer_list_type),    target, allocatable  :: nml_tracer_list(:)

namelist /tracer_listsize/ num_tracers
namelist /tracer_list    / nml_tracer_list
namelist /tracer_general / smooth_bh_tra, gamma0_tra, gamma1_tra, gamma2_tra, i_vert_diff

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

! OPEN and read namelist for I/O
open( unit=nm_unit, file='namelist.tra', form='formatted', access='sequential', status='old', iostat=iost )
if (iost == 0) then
   if (mype==0) WRITE(*,*) '     file   : ', 'namelist.tra',' open ok'
else
   if (mype==0) WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist.tra',' ; iostat=',iost
   call par_ex(partit%MPI_COMM_FESOM, partit%mype)
   stop
end if

READ(nm_unit,  nml=tracer_listsize, iostat=iost)
allocate(nml_tracer_list(num_tracers))
READ(nm_unit,  nml=tracer_list,     iostat=iost)
read (nm_unit, nml=tracer_init3d,   iostat=iost)
READ(nm_unit,  nml=tracer_general,  iostat=iost)
close(nm_unit)

do n=1, num_tracers
   if (nml_tracer_list(n)%id==-1) then
      if (mype==0) write(*,*) 'number of tracers will be changed from ', num_tracers, ' to ', n-1, '!'
      num_tracers=n-1
      EXIT
   end if
end do

if (mype==0) write(*,*) 'total number of tracers is: ', num_tracers

elem_size=myDim_elem2D+eDim_elem2D
node_size=myDim_nod2D+eDim_nod2D

tracers%num_tracers=num_tracers

! ================
! Temperature (index=1), Salinity (index=2), etc.
! ================
allocate(tracers%data(num_tracers))
do n=1, tracers%num_tracers
   allocate(tracers%data(n)%values  (nl-1,node_size))
   allocate(tracers%data(n)%valuesAB(nl-1,node_size))
   tracers%data(n)%ID            = nml_tracer_list(n)%id
   tracers%data(n)%tra_adv_hor   = TRIM(nml_tracer_list(n)%adv_hor)
   tracers%data(n)%tra_adv_ver   = TRIM(nml_tracer_list(n)%adv_ver)
   tracers%data(n)%tra_adv_lim   = TRIM(nml_tracer_list(n)%adv_lim)
   tracers%data(n)%tra_adv_ph    = nml_tracer_list(n)%adv_ph
   tracers%data(n)%tra_adv_pv    = nml_tracer_list(n)%adv_pv
   tracers%data(n)%smooth_bh_tra = smooth_bh_tra
   tracers%data(n)%gamma0_tra    = gamma0_tra
   tracers%data(n)%gamma1_tra    = gamma1_tra
   tracers%data(n)%gamma2_tra    = gamma2_tra
   tracers%data(n)%values        = 0.
   tracers%data(n)%valuesAB      = 0.
   tracers%data(n)%i_vert_diff   = i_vert_diff
end do
allocate(tracers%work%del_ttf(nl-1,node_size))
allocate(tracers%work%del_ttf_advhoriz(nl-1,node_size),tracers%work%del_ttf_advvert(nl-1,node_size))
tracers%work%del_ttf          = 0.0_WP
tracers%work%del_ttf_advhoriz = 0.0_WP
tracers%work%del_ttf_advvert  = 0.0_WP
if (ldiag_DVD) then
    allocate(tracers%work%tr_dvd_horiz(nl-1,node_size,2),tracers%work%tr_dvd_vert(nl-1,node_size,2))
    tracers%work%tr_dvd_horiz = 0.0_WP
    tracers%work%tr_dvd_vert  = 0.0_WP
end if
END SUBROUTINE tracer_init
!
!
!_______________________________________________________________________________
SUBROUTINE dynamics_init(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    IMPLICIT NONE
    integer        :: elem_size, node_size
    integer, save  :: nm_unit  = 104       ! unit to open namelist file, skip 100-102 for cray
    integer        :: iost

    integer        :: visc_opt
    real(kind=WP)  :: gamma0_visc, gamma1_visc, gamma2_visc
    real(kind=WP)  :: div_c_visc, leith_c_visc, easybackscat_return
    logical        :: use_ivertvisc
    integer        :: momadv_opt
    logical        :: use_freeslip
    logical        :: use_wsplit
    real(kind=WP)  :: wsplit_maxcfl

    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_dyn)   , intent(inout), target :: dynamics
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

!!PS     ! define dynamics namelist parameter
!!PS     namelist /dynamics_visc    / visc_opt, gamma0_visc, gamma1_visc, gamma2_visc,  &
!!PS                                  div_c_visc, leith_c_visc, use_ivertvisc, easy_bs_return
!!PS     namelist /dynamics_general / momadv_opt, use_freeslip, use_wsplit, wsplit_maxcfl 
!!PS 
!!PS     ! open and read namelist for I/O
!!PS     open(unit=nm_unit, file='namelist.dyn', form='formatted', access='sequential', status='old', iostat=iost )
!!PS     if (iost == 0) then
!!PS         if (mype==0) write(*,*) '     file   : ', 'namelist.dyn',' open ok'
!!PS     else
!!PS         if (mype==0) write(*,*) 'ERROR: --> bad opening file   : ', 'namelist.dyn',' ; iostat=',iost
!!PS         call par_ex(partit%MPI_COMM_FESOM, partit%mype)
!!PS         stop
!!PS     end if
!!PS     read(nm_unit, nml=dynamics_visc   , iostat=iost)
!!PS     read(nm_unit, nml=dynamics_general, iostat=iost)
!!PS     close(nm_unit)

    ! define local vertice & elem array size
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D+eDim_nod2D

    ! allocate data arrays in derived type
    allocate(dynamics%uv(        2, nl-1, elem_size))
    allocate(dynamics%uv_rhs(    2, nl-1, elem_size))
    allocate(dynamics%uv_rhsAB(  2, nl-1, elem_size))
    allocate(dynamics%uvnode(    2, nl-1, node_size))
    allocate(dynamics%uvnode_rhs(2, nl-1, node_size))
    dynamics%uv         = 0.0_WP
    dynamics%uv_rhs     = 0.0_WP
    dynamics%uv_rhsAB   = 0.0_WP
    dynamics%uvnode     = 0.0_WP
    dynamics%uvnode_rhs = 0.0_WP
    
    allocate(dynamics%w(              nl, node_size))
    allocate(dynamics%w_e(            nl, node_size))
    allocate(dynamics%w_i(            nl, node_size))
    allocate(dynamics%cfl_z(          nl, node_size))
    dynamics%w          = 0.0_WP
    dynamics%w_e        = 0.0_WP
    dynamics%w_i        = 0.0_WP
    dynamics%cfl_z      = 0.0_WP
    
    allocate(dynamics%eta_n(      node_size))
    allocate(dynamics%d_eta(      node_size))
    allocate(dynamics%ssh_rhs(    node_size))
    allocate(dynamics%ssh_rhs_old(node_size))
    dynamics%eta_n      = 0.0_WP
    dynamics%d_eta      = 0.0_WP
    dynamics%ssh_rhs    = 0.0_WP
    dynamics%ssh_rhs_old= 0.0_WP    
    
    ! set parameters in derived type
!!PS     dynamics%visc_opt      = visc_opt
!!PS     dynamics%gamma0_visc   = gamma0_visc
!!PS     dynamics%gamma1_visc   = gamma1_visc
!!PS     dynamics%gamma2_visc   = gamma2_visc
!!PS     dynamics%div_c_visc    = div_c_visc
!!PS     dynamics%leith_c_visc  = leith_c_visc
!!PS     dynamics%use_ivertvisc = use_ivertvisc
!!PS     dynamics%momadv_opt    = momadv_opt
!!PS     dynamics%use_freeslip  = use_freeslip
!!PS     dynamics%use_wsplit    = use_wsplit
!!PS     dynamics%wsplit_maxcfl = wsplit_maxcfl

    dynamics%visc_opt      = visc_option
    dynamics%gamma0_visc   = gamma0
    dynamics%gamma1_visc   = gamma1
    dynamics%gamma2_visc   = gamma2
    dynamics%div_c_visc    = Div_c
    dynamics%leith_c_visc  = Leith_c
    dynamics%use_ivertvisc = i_vert_visc
    dynamics%momadv_opt    = mom_adv
    dynamics%use_freeslip  = free_slip
    dynamics%use_wsplit    = w_split
    dynamics%wsplit_maxcfl = w_max_cfl
END SUBROUTINE dynamics_init
!
!
!_______________________________________________________________________________
SUBROUTINE arrays_init(num_tracers, partit, mesh)
USE MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
USE o_ARRAYS
USE o_PARAM
use g_comm_auto
use g_config
use g_forcing_arrays
use o_mixing_kpp_mod ! KPP
USE g_forcing_param, only: use_virt_salt
use diagnostics,     only: ldiag_dMOC, ldiag_DVD
IMPLICIT NONE
integer                               :: elem_size, node_size
integer                               :: n
integer,        intent(in)            :: num_tracers
type(t_mesh),   intent(in),    target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"


elem_size=myDim_elem2D+eDim_elem2D
node_size=myDim_nod2D+eDim_nod2D


! ================
! Velocities
! ================     
!allocate(stress_diag(2, elem_size))!delete me
!!PS allocate(UV(2, nl-1, elem_size))
allocate(UV_rhs(2,nl-1, elem_size))
allocate(UV_rhsAB(2,nl-1, elem_size))
allocate(Visc(nl-1, elem_size))
! ================
! elevation and its rhs
! ================
allocate(eta_n(node_size), d_eta(node_size))
allocate(ssh_rhs(node_size))
! ================
! Monin-Obukhov
! ================
if (use_ice .and. use_momix) allocate(mo(nl,node_size),mixlength(node_size))
if (use_ice .and. use_momix) mixlength=0.
! ================
! Vertical velocity and pressure
! ================
allocate(Wvel(nl, node_size), hpressure(nl,node_size))
allocate(Wvel_e(nl, node_size), Wvel_i(nl, node_size))
allocate(CFL_z(nl, node_size)) ! vertical CFL criteria
allocate(bvfreq(nl,node_size),mixlay_dep(node_size),bv_ref(node_size))
! ================
! Ocean forcing arrays
! ================
allocate(Tclim(nl-1,node_size), Sclim(nl-1, node_size))
allocate(stress_surf(2,myDim_elem2D))    !!! Attention, it is shorter !!! 
allocate(stress_node_surf(2,node_size))
allocate(stress_atmoce_x(node_size), stress_atmoce_y(node_size)) 
allocate(relax2clim(node_size)) 
allocate(heat_flux(node_size), Tsurf(node_size))
allocate(water_flux(node_size), Ssurf(node_size))
allocate(relax_salt(node_size))
allocate(virtual_salt(node_size))

allocate(heat_flux_in(node_size))
allocate(real_salt_flux(node_size)) !PS
! =================
! Arrays used to organize surface forcing
! =================
allocate(Tsurf_t(node_size,2), Ssurf_t(node_size,2))
allocate(tau_x_t(node_size,2), tau_y_t(node_size,2))  

! =================
! All auxiliary arrays
! =================
 
!if(mom_adv==3) then
allocate(vorticity(nl-1,node_size))
vorticity=0.0_WP
!end if

! =================
! Visc and Diff coefs
! =================

allocate(Av(nl,elem_size), Kv(nl,node_size))

Av=0.0_WP
Kv=0.0_WP
if (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then
   allocate(Kv_double(nl,node_size, num_tracers))
   Kv_double=0.0_WP
   !!PS call oce_mixing_kpp_init ! Setup constants, allocate arrays and construct look up table
end if

! =================
! Backscatter arrays
! =================

if(visc_option==8) then

allocate(uke(nl-1,elem_size)) ! Unresolved kinetic energy for backscatter coefficient
allocate(v_back(nl-1,elem_size))  ! Backscatter viscosity
allocate(uke_dis(nl-1,elem_size), uke_back(nl-1,elem_size)) 
allocate(uke_dif(nl-1,elem_size))
allocate(uke_rhs(nl-1,elem_size), uke_rhs_old(nl-1,elem_size))
allocate(UV_dis_tend(2,nl-1,elem_size), UV_back_tend(2,nl-1,elem_size))
allocate(UV_total_tend(2,nl-1,elem_size))

uke=0.0_8
v_back=0.0_8
uke_dis=0.0_8
uke_dif=0.0_8
uke_back=0.0_8
uke_rhs=0.0_8
uke_rhs_old=0.0_8
UV_dis_tend=0.0_8
UV_back_tend=0.0_8
UV_total_tend=0.0_8
end if

!Velocities at nodes
allocate(Unode(2,nl-1,node_size))

! tracer gradients & RHS  
allocate(ttrhs(nl-1,node_size))
allocate(tr_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
allocate(tr_z(nl,myDim_nod2D+eDim_nod2D))

! neutral slope etc. to be used in Redi formulation
allocate(neutral_slope(3, nl-1, node_size))
allocate(slope_tapered(3, nl-1, node_size))
allocate(Ki(nl-1, node_size))

do n=1, node_size
!  Ki(n)=K_hor*area(1,n)/scale_area
   Ki(:,n)=K_hor*(mesh_resolution(n)/100000.0_WP)**2
end do
call exchange_nod(Ki, partit)

neutral_slope=0.0_WP
slope_tapered=0.0_WP

allocate(MLD1(node_size), MLD2(node_size), MLD1_ind(node_size), MLD2_ind(node_size))
if (use_global_tides) then
   allocate(ssh_gp(node_size))
   ssh_gp=0.
end if
! xy gradient of a neutral surface
allocate(sigma_xy(2, nl-1, node_size))
sigma_xy=0.0_WP

! alpha and beta in the EoS
allocate(sw_beta(nl-1, node_size), sw_alpha(nl-1, node_size))
allocate(dens_flux(node_size))
sw_beta  =0.0_WP
sw_alpha =0.0_WP
dens_flux=0.0_WP

if (Fer_GM) then
   allocate(fer_c(node_size),fer_scal(node_size), fer_gamma(2, nl, node_size), fer_K(nl, node_size))
   allocate(fer_wvel(nl, node_size), fer_UV(2, nl-1, elem_size))
   fer_gamma=0.0_WP
   fer_uv=0.0_WP
   fer_wvel=0.0_WP
   fer_K=500._WP
   fer_c=1._WP
   fer_scal = 0.0_WP
end if

if (SPP) then
   allocate(ice_rejected_salt(node_size))
   ice_rejected_salt=0._WP
end if

! =================
! Initialize with zeros 
! =================

!!PS     UV=0.0_WP
    UV_rhs=0.0_WP
    UV_rhsAB=0.0_WP
!
    eta_n=0.0_WP
    d_eta=0.0_WP
    ssh_rhs=0.0_WP
    Wvel=0.0_WP
    Wvel_e	=0.0_WP
    Wvel_i	=0.0_WP
    CFL_z   =0.0_WP
    hpressure=0.0_WP
!
    heat_flux=0.0_WP
    heat_flux_in=0.0_WP
    Tsurf=0.0_WP

    water_flux=0.0_WP
    relax_salt=0.0_WP
    virtual_salt=0.0_WP

    Ssurf=0.0_WP
    
    real_salt_flux=0.0_WP
    
    stress_surf      =0.0_WP
    stress_node_surf =0.0_WP
    stress_atmoce_x  =0.0_WP
    stress_atmoce_y  =0.0_WP
    
    bvfreq=0.0_WP
    mixlay_dep=0.0_WP
    bv_ref=0.0_WP

    MLD1   =0.0_WP
    MLD2   =0.0_WP
    MLD1_ind=0.0_WP
    MLD2_ind=0.0_WP

    relax2clim=0.0_WP

    Tsurf_t=0.0_WP
    Ssurf_t=0.0_WP
    tau_x_t=0.0_WP
    tau_y_t=0.0_WP
    
    ! init field for pressure force 
    allocate(density_ref(nl-1,node_size))
    density_ref = density_0
    allocate(density_m_rho0(nl-1, node_size))
    allocate(density_m_rho0_slev(nl-1, node_size)) !!PS
    if (ldiag_dMOC) then
       allocate(density_dMOC       (nl-1, node_size))
    end if
    allocate(pgf_x(nl-1, elem_size),pgf_y(nl-1, elem_size)) 
    density_m_rho0=0.0_WP
    density_m_rho0_slev=0.0_WP !!PS
    if (ldiag_dMOC) then
       density_dMOC       =0.0_WP
    end if
    pgf_x = 0.0_WP
    pgf_y = 0.0_WP
    
!!PS     ! init dummy arrays
!!PS     allocate(dum_2d_n(node_size), dum_3d_n(nl-1,node_size))
!!PS     allocate(dum_2d_e(elem_size), dum_3d_e(nl-1,elem_size)) 
!!PS     dum_2d_n = 0.0_WP
!!PS     dum_3d_n = 0.0_WP
!!PS     dum_2d_e = 0.0_WP
!!PS     dum_3d_e = 0.0_WP
END SUBROUTINE arrays_init
!
!
!_______________________________________________________________________________
! Here the 3D tracers will be initialized. Initialization strategy depends on a tracer ID.
! ID = 0 and 1 are reserved for temperature and salinity
SUBROUTINE oce_initial_state(tracers, partit, mesh)
USE MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
USE MOD_TRACER
USE o_ARRAYS
USE g_config
USE g_ic3d
  !
  ! reads the initial state or the restart file for the ocean
  !
  implicit none
  integer                  :: i, k, counter, rcounter3, id
  character(len=10)        :: i_string, id_string
  type(t_mesh),   intent(in) ,   target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_tracer), intent(inout), target :: tracers
  real(kind=WP)            :: loc, max_temp, min_temp, max_salt, min_salt

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  if (mype==0) write(*,*) tracers%num_tracers, ' tracers will be used in FESOM'
  if (mype==0) write(*,*) 'tracer IDs are: ', tracers%data(1:tracers%num_tracers)%ID
  !
  ! read ocean state
  ! this must be always done! First two tracers with IDs 0 and 1 are the temperature and salinity.
  if(mype==0) write(*,*) 'read Temperatur climatology from:', trim(filelist(1))
  if(mype==0) write(*,*) 'read Salt       climatology from:', trim(filelist(2))
  call do_ic3d(tracers, partit, mesh)
  
  Tclim=tracers%data(1)%values
  Sclim=tracers%data(2)%values
  Tsurf=Tclim(1,:)
  Ssurf=Sclim(1,:)
  relax2clim=0.0_WP

  ! count the passive tracers which require 3D source (ptracers_restore_total)
  ptracers_restore_total=0
  DO i=3, tracers%num_tracers
     id=tracers%data(i)%ID
     SELECT CASE (id)
     CASE (301)
          ptracers_restore_total=ptracers_restore_total+1
     CASE (302)
          ptracers_restore_total=ptracers_restore_total+1
     CASE (303)
          ptracers_restore_total=ptracers_restore_total+1

     END SELECT
  END DO
  allocate(ptracers_restore(ptracers_restore_total))
  
  rcounter3=0         ! counter for tracers with 3D source
  DO i=3, tracers%num_tracers
     id=tracers%data(i)%ID
     SELECT CASE (id)
       CASE (101)       ! initialize tracer ID=101
         tracers%data(i)%values(:,:)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I3)") i
            write (id_string, "(I3)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (301) !Fram Strait 3d restored passive tracer
         tracers%data(i)%values(:,:)=0.0_WP
         rcounter3    =rcounter3+1
         counter=0
         do k=1, myDim_nod2D+eDim_nod2D
            if     (((geo_coord_nod2D(2,k)>77.5*rad) .and. (geo_coord_nod2D(2,k)<78.*rad))&
               .and.((geo_coord_nod2D(1,k)>0.  *rad) .and. (geo_coord_nod2D(1,k)<10.*rad))) then
               counter=counter+1
            end if
         end do
         allocate(ptracers_restore(rcounter3)%ind2(counter))
         ptracers_restore(rcounter3)%id   =301
         ptracers_restore(rcounter3)%locid=i
         counter=0
         do k=1, myDim_nod2D+eDim_nod2D
            if     (((geo_coord_nod2D(2,k)>77.5*rad) .and. (geo_coord_nod2D(2,k)<78.*rad))&
               .and.((geo_coord_nod2D(1,k)>0.  *rad) .and. (geo_coord_nod2D(1,k)<10.*rad))) then
               counter=counter+1
               ptracers_restore(rcounter3)%ind2(counter)=k
            end if
         end do
         tracers%data(i)%values(:,ptracers_restore(rcounter3)%ind2)=1.
         if (mype==0) then
            write (i_string,  "(I3)") i
            write (id_string, "(I3)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if

       CASE (302) !Bering Strait 3d restored passive tracer
         tracers%data(i)%values(:,:)=0.0_WP
         rcounter3    =rcounter3+1
         counter=0
         do k=1, myDim_nod2D+eDim_nod2D
            if     (((geo_coord_nod2D(2,k)>65.6*rad) .and. (geo_coord_nod2D(2,k)<66.*rad))&
               .and.((geo_coord_nod2D(1,k)>-172.  *rad) .and. (geo_coord_nod2D(1,k)<-166.*rad))) then
               counter=counter+1
            end if
         end do
         allocate(ptracers_restore(rcounter3)%ind2(counter))
         ptracers_restore(rcounter3)%id   =302
         ptracers_restore(rcounter3)%locid=i
         counter=0
         do k=1, myDim_nod2D+eDim_nod2D
            if     (((geo_coord_nod2D(2,k)>65.6*rad) .and. (geo_coord_nod2D(2,k)<66.*rad))&
               .and.((geo_coord_nod2D(1,k)>-172.  *rad) .and. (geo_coord_nod2D(1,k)<-166.*rad))) then
               counter=counter+1
               ptracers_restore(rcounter3)%ind2(counter)=k
            end if
         end do
         tracers%data(i)%values(:,ptracers_restore(rcounter3)%ind2)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I3)") i
            write (id_string, "(I3)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       
      CASE (303) !BSO 3d restored passive tracer
         tracers%data(i)%values(:,:)=0.0_WP
         rcounter3    =rcounter3+1
         counter=0
         do k=1, myDim_nod2D+eDim_nod2D
            if     (((geo_coord_nod2D(2,k)>69.5*rad) .and. (geo_coord_nod2D(2,k)<74.5*rad))&
               .and.((geo_coord_nod2D(1,k)>19.  *rad) .and. (geo_coord_nod2D(1,k)<20.*rad))) then
               counter=counter+1
            end if
         end do
         allocate(ptracers_restore(rcounter3)%ind2(counter))
         ptracers_restore(rcounter3)%id   =303
         ptracers_restore(rcounter3)%locid=i
         counter=0
         do k=1, myDim_nod2D+eDim_nod2D
            if     (((geo_coord_nod2D(2,k)>69.5*rad) .and. (geo_coord_nod2D(2,k)<74.5*rad))&
               .and.((geo_coord_nod2D(1,k)>19.  *rad) .and. (geo_coord_nod2D(1,k)<20.*rad))) then
               counter=counter+1
               ptracers_restore(rcounter3)%ind2(counter)=k
            end if
         end do
         tracers%data(i)%values(:,ptracers_restore(rcounter3)%ind2)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I3)") i
            write (id_string, "(I3)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE DEFAULT
         if (mype==0) then
            write (i_string,  "(I3)") i
            write (id_string, "(I3)") id
            if (mype==0) write(*,*) 'invalid ID '//trim(id_string)//' specified for '//trim(i_string)//' th tracer!!!'
            if (mype==0) write(*,*) 'the model will stop!'
         end if
         call par_ex(partit%MPI_COMM_FESOM, partit%mype)
         stop
     END SELECT
  END DO
end subroutine oce_initial_state
!
!
!==========================================================================
! Here we do things (if applicable) before the ocean timestep will be made
SUBROUTINE before_oce_step(dynamics, tracers, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    USE MOD_DYN
    USE o_ARRAYS
    USE g_config
    USE Toy_Channel_Soufflet
    implicit none
    integer                  :: i, k, counter, rcounter3, id
    character(len=10)        :: i_string, id_string
    type(t_mesh),   intent(in),    target  :: mesh
    type(t_partit), intent(inout), target  :: partit
    type(t_tracer), intent(inout), target  :: tracers
    type(t_dyn), intent(inout), target  :: dynamics

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    if (toy_ocean) then
        SELECT CASE (TRIM(which_toy))
            CASE ("soufflet") !forcing update for soufflet testcase
            if (mod(mstep, soufflet_forc_update)==0) then
                call compute_zonal_mean(dynamics, tracers, partit, mesh)
            end if
        END SELECT
    end if
END SUBROUTINE before_oce_step
