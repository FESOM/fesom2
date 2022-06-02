module array_setup_interface
  interface
    subroutine array_setup(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module oce_initial_state_interface
  interface
    subroutine oce_initial_state(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
!
!
!_______________________________________________________________________________
subroutine ocean_setup(mesh)
USE MOD_MESH
USE o_PARAM
USE g_PARSUP
USE o_ARRAYS
USE g_config
USE g_forcing_param, only: use_virt_salt
use g_cvmix_tke
use g_cvmix_idemix
use g_cvmix_pp
use g_cvmix_kpp
use g_cvmix_tidal
use Toy_Channel_Soufflet
use array_setup_interface
use oce_initial_state_interface
use oce_adv_tra_fct_interfaces
IMPLICIT NONE
type(t_mesh), intent(inout) , target :: mesh
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
    call array_setup(mesh)
    
    !___________________________________________________________________________
    ! initialize arrays for ALE
    if (mype==0) then
       write(*,*) '____________________________________________________________'
       write(*,*) ' --> initialise ALE arrays + sparse SSH stiff matrix'
       write(*,*)
    end if
    call init_ale(mesh)
    call init_stiff_mat_ale(mesh) !!PS test  
    
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
            call par_ex
    end select

    ! initialise fesom1.4 like KPP
    if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then
        call oce_mixing_kpp_init(mesh)
    ! initialise fesom1.4 like PP
    elseif (mix_scheme_nmb==2 .or. mix_scheme_nmb==27) then
    
    ! initialise cvmix_KPP
    elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then
        call init_cvmix_kpp(mesh)
        
    ! initialise cvmix_PP    
    elseif (mix_scheme_nmb==4 .or. mix_scheme_nmb==47) then
        call init_cvmix_pp(mesh)
        
    ! initialise cvmix_TKE    
    elseif (mix_scheme_nmb==5 .or. mix_scheme_nmb==56) then
        call init_cvmix_tke(mesh)
        
    endif
  
    ! initialise additional mixing cvmix_IDEMIX --> only in combination with 
    ! cvmix_TKE+cvmix_IDEMIX or stand alone for debbuging as cvmix_TKE
    if     (mod(mix_scheme_nmb,10)==6) then
        call init_cvmix_idemix(mesh)
        
    ! initialise additional mixing cvmix_TIDAL --> only in combination with 
    ! KPP+cvmix_TIDAL, PP+cvmix_TIDAL, cvmix_KPP+cvmix_TIDAL, cvmix_PP+cvmix_TIDAL 
    ! or stand alone for debbuging as cvmix_TIDAL   
    elseif (mod(mix_scheme_nmb,10)==7) then
        call init_cvmix_tidal(mesh)
    end if         
    
    !___________________________________________________________________________
    ! set use_density_ref .true. when cavity is used and initialse cavity boundary 
    ! line for the extrapolation of the initialisation
    if (use_cavity .and. .not. use_density_ref) use_density_ref=.true.
    
    ! compute for all cavity points (ulevels_nod2D>1), which is the closest
    ! cavity line point to that point --> use their coordinates and depth -->
    ! use for extrapolation of init state under cavity
    if (use_cavity) call compute_nrst_pnt2cavline(mesh)
      
    if (use_density_ref) call init_ref_density(mesh)
    
    
    !___________________________________________________________________________
    if(mype==0) write(*,*) 'Arrays are set'
        
    !if(open_boundary) call set_open_boundary   !TODO
    
    call oce_adv_tra_fct_init(mesh)
    call muscl_adv_init(mesh) !!PS test
    !=====================
    ! Initialize fields
    ! A user-defined routine has to be called here!
    !=====================
    if (toy_ocean) then  
       SELECT CASE (TRIM(which_toy))
         CASE ("soufflet") !forcing update for soufflet testcase
           if (mod(mstep, soufflet_forc_update)==0) then
              call initial_state_soufflet(mesh)
              call compute_zonal_mean_ini(mesh)  
              call compute_zonal_mean(mesh)
           end if
       END SELECT
    else
       call oce_initial_state(mesh)   ! Use it if not running tests
    end if

    if (.not.r_restart) tr_arr_old=tr_arr
    
    !___________________________________________________________________________
    ! first time fill up array for hnode & helem
    if (mype==0) then
        write(*,*) '____________________________________________________________'
        write(*,*) ' --> call init_thickness_ale'
        write(*,*)
    end if
    call init_thickness_ale(mesh)
    
    !___________________________________________________________________________
    if(mype==0) write(*,*) 'Initial state'
    if (w_split .and. mype==0) then
        write(*,*) '******************************************************************************'
        write(*,*) 'vertical velocity will be split onto explicit and implicit constitutes;'
        write(*,*) 'maximum allowed CDF on explicit W is set to: ', w_max_cfl
        write(*,*) '******************************************************************************'
    end if
end subroutine ocean_setup
!
!
!_______________________________________________________________________________
SUBROUTINE array_setup(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
use g_comm_auto
use g_config
use g_forcing_arrays
use o_mixing_kpp_mod ! KPP
USE g_forcing_param, only: use_virt_salt
use diagnostics,     only: ldiag_dMOC, ldiag_DVD
#if defined(__recom)
use REcoM_GloVar
use recom_config
#endif
IMPLICIT NONE
integer     :: elem_size, node_size
integer     :: n
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"


elem_size=myDim_elem2D+eDim_elem2D
node_size=myDim_nod2D+eDim_nod2D


! ================
! Velocities
! ================     
!allocate(stress_diag(2, elem_size))!delete me
allocate(UV(2, nl-1, elem_size))
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
! ================
! Temperature and salinity
! ================
allocate(T_rhs(nl-1, node_size))
allocate(S_rhs(nl-1, node_size))
allocate(tr_arr(nl-1,node_size,num_tracers),tr_arr_old(nl-1,node_size,num_tracers))
allocate(del_ttf(nl-1,node_size))
allocate(del_ttf_advhoriz(nl-1,node_size),del_ttf_advvert(nl-1,node_size))
allocate(dtr_bf(nl-1,node_size)) !jh
allocate(str_bf(nl-1,node_size)) !OG
allocate(vert_sink(nl-1,node_size)) ! OG
allocate(nss(nl-1,node_size)) ! OG

del_ttf          = 0.0_WP
del_ttf_advhoriz = 0.0_WP
del_ttf_advvert  = 0.0_WP
dtr_bf           = 0.0_WP ! jh
str_bf           = 0.0_WP ! OG
vert_sink        = 0.0_WP ! OG
nss              = 0.0_WP ! OG

!!PS allocate(del_ttf_diff(nl-1,node_size))
if (ldiag_DVD) then
    allocate(tr_dvd_horiz(nl-1,node_size,2),tr_dvd_vert(nl-1,node_size,2))
    tr_dvd_horiz = 0.0_WP
    tr_dvd_vert  = 0.0_WP
end if 

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
! ================
! RECOM forcing arrays
! ================
#if defined(__recom)
  if(use_REcoM) then
    !if (restore_alkalinity) then
      allocate(Alk_surf(node_size))
      allocate(relax_alk(node_size))
      allocate(virtual_alk(node_size))
    !endif
  end if
#endif
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
   allocate(Kv_double(nl,node_size,num_tracers))
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
call exchange_nod(Ki)

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

    UV=0.0_WP
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
    T_rhs=0.0_WP
    heat_flux=0.0_WP
    heat_flux_in=0.0_WP
    Tsurf=0.0_WP

    S_rhs=0.0_WP
    water_flux=0.0_WP
    relax_salt=0.0_WP
    virtual_salt=0.0_WP

    Ssurf=0.0_WP
    
    real_salt_flux=0.0_WP
    
    stress_surf      =0.0_WP
    stress_node_surf =0.0_WP
    stress_atmoce_x  =0.0_WP
    stress_atmoce_y  =0.0_WP
    
    tr_arr=0.0_WP
    tr_arr_old=0.0_WP

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
! ================
! RECOM forcing arrays
! ================
#if defined(__recom)
  if(use_REcoM) then
    !if (restore_alkalinity) then
      Alk_surf=0.0_WP
      relax_alk=0.0_WP
      virtual_alk=0.0_WP
    !endif
  end if
#endif        
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
END SUBROUTINE array_setup
!
!
!_______________________________________________________________________________
! Here the 3D tracers will be initialized. Initialization strategy depends on a tracer ID.
! ID = 0 and 1 are reserved for temperature and salinity
SUBROUTINE oce_initial_state(mesh)
USE MOD_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_config
USE g_ic3d
#if defined(__recom)
USE recom_config
USE REcoM_GloVar
use REcoM_ciso
#endif
  !
  ! reads the initial state or the restart file for the ocean
  !
  implicit none
  integer                  :: i, k, counter, rcounter3, id
  character(len=10)        :: i_string, id_string
  type(t_mesh), intent(in) , target :: mesh
  real(kind=WP)            :: loc, max_temp, min_temp, max_salt, min_salt

#include "associate_mesh.h"

  if (mype==0) write(*,*) num_tracers, ' tracers will be used in FESOM'
  if (mype==0) write(*,*) 'tracer IDs are: ', tracer_ID(1:num_tracers)
  !
  ! read ocean state
  ! this must be always done! First two tracers with IDs 0 and 1 are the temperature and salinity.
#if defined(__recom)
  if(use_REcoM) then
    if (mype==0) write(*,*)      
    if (mype==0) print *, achar(27)//'[36m'//'     --> RECOM ON'//achar(27)//'[0m'
    if (ciso) then
      if (mype==0) print *, achar(27)//'[36m'//'     --> CISO ON'//achar(27)//'[0m'
    else
      if (mype==0) print *, achar(27)//'[36m'//'     --> CISO OFF'//achar(27)//'[0m'
    endif
    if (mype==0) then
      write(*,*) '____________________________________________________________'
      write(*,*) ' --> Check namelist.oce for the number of tracers'
      write(*,*)      
    endif

    if(DIC_PI) then 
      filelist(5)='GLODAPv2.2016b.PI_TCO2_fesom2_mmol_fix_z_Fillvalue.nc'
      varlist(5)='PI_TCO2_mmol'
      if (mype==0) then
        write(*,*) '____________________________________________________________'
        write(*,*) ' --> Reading REcoM2 preindustrial DIC for restart'
        write(*,*)
      end if
    end if

    if(mype==0) write(*,*) 'read Iron        climatology from:', trim(filelist(1))
    if(mype==0) write(*,*) 'read Oxygen      climatology from:', trim(filelist(2))
    if(mype==0) write(*,*) 'read Silicate    climatology from:', trim(filelist(3))
    if(mype==0) write(*,*) 'read Alkalinity  climatology from:', trim(filelist(4))
    if(mype==0) write(*,*) 'read DIC         climatology from:', trim(filelist(5))
    if(mype==0) write(*,*) 'read Nitrate     climatology from:', trim(filelist(6))
    if(mype==0) write(*,*) 'read Salt        climatology from:', trim(filelist(7))
    if(mype==0) write(*,*) 'read Temperature climatology from:', trim(filelist(8))
#else
  if(mype==0) write(*,*) 'read Salt       climatology from:', trim(filelist(1))
  if(mype==0) write(*,*) 'read Temperatur climatology from:', trim(filelist(2))
#endif

  call do_ic3d(mesh)
  
  Tclim=tr_arr(:,:,1)
  Sclim=tr_arr(:,:,2)
  Tsurf=tr_arr(1,:,1)
  Ssurf=tr_arr(1,:,2)
  relax2clim=0.0_WP

#if defined(__recom)
    if (restore_alkalinity) then
      if (mype==0)  print *, achar(27)//'[46;1m'//'--> Set surface field for alkalinity restoring'//achar(27)//'[0m'
      Alk_surf=tr_arr(1,:,5)
      if(mype==0) write(*,*),'Alkalinity restoring = true. Field read in.'
    endif
  end if
#endif

  ! count the passive tracers which require 3D source (ptracers_restore_total)
  ptracers_restore_total=0
  DO i=3, num_tracers
     id=tracer_ID(i)
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
  DO i=3, num_tracers
     id=tracer_ID(i)
     if (any(id == idlist)) cycle ! OG recom tracers id's start from 1001
     SELECT CASE (id)
! Read recom variables (hardcoded IDs) OG
       CASE (1004:1017)
         tr_arr(:,:,i)=0.0_WP
        if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (1020:1021)
         tr_arr(:,:,i)=0.0_WP
        if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
   CASE (1023:1028)
         tr_arr(:,:,i)=0.0_WP
        if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
!MB preliminary extension for carbon isotopes which should be refactored in later versions
       CASE (1033:1034)       ! initialize tracer ID=33-34 = DIC_13|14
         tr_arr(:,:,25:26)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (1093:1094)       ! initialize tracer ID=93:94 = phyC_13|14
         tr_arr(:,:,27:28)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (1123:1124)       ! initialize tracer ID=123:124 = detC_13|14
         tr_arr(:,:,29:30)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (1143:1144)       ! initialize tracer ID=143:144 = hetC_13|14
         tr_arr(:,:,31:32)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (1163:1164)       ! initialize tracer ID=163:164 = DOC_13|14
         tr_arr(:,:,33:34)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (1183:1184)       ! initialize tracer ID=183:184 = diaC_13|14
         tr_arr(:,:,35:36)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (1223:1224)       ! initialize tracer ID=223:224 = phyCalc_13|14
         tr_arr(:,:,37:38)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (1233:1234)       ! initialize tracer ID=233:234 = detCalc_13|14
         tr_arr(:,:,39:40)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if

       CASE (101)       ! initialize tracer ID=101
         tr_arr(:,:,i)=0.0_WP
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       CASE (301) !Fram Strait 3d restored passive tracer
         tr_arr(:,:,i)=0.0_WP
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
         tr_arr(:,ptracers_restore(rcounter3)%ind2,i)=1.
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if

       CASE (302) !Bering Strait 3d restored passive tracer
         tr_arr(:,:,i)=0.
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
         tr_arr(:,ptracers_restore(rcounter3)%ind2,i)=1.
         if (mype==0) then
            write (i_string,  "(I4)") i
            write (id_string, "(I4)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
         end if
       
      CASE (303) !BSO 3d restored passive tracer
         tr_arr(:,:,i)=0.
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
         tr_arr(:,ptracers_restore(rcounter3)%ind2,i)=1.
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
         call par_ex
         stop
     END SELECT
  END DO
end subroutine oce_initial_state
!
!
!==========================================================================
! Here we do things (if applicable) before the ocean timestep will be made
SUBROUTINE before_oce_step(mesh)
    USE MOD_MESH
    USE o_ARRAYS
    USE g_PARSUP
    USE g_config
    USE Toy_Channel_Soufflet
    implicit none
    integer                  :: i, k, counter, rcounter3, id
    character(len=10)        :: i_string, id_string
    type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

    if (toy_ocean) then
        SELECT CASE (TRIM(which_toy))
            CASE ("soufflet") !forcing update for soufflet testcase
            if (mod(mstep, soufflet_forc_update)==0) then
                call compute_zonal_mean(mesh)
            end if
        END SELECT
    end if
END SUBROUTINE before_oce_step
