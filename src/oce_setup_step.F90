module oce_initial_state_interface
    interface
        subroutine oce_initial_state(tracers, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in)  ,  target :: mesh
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
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in)  ,  target :: mesh
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
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
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
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout)   , target :: mesh
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
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module
!
!
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
    use g_backscatter
    use Toy_Channel_Soufflet
    use Toy_Channel_Nemo
    use oce_initial_state_interface
    use oce_adv_tra_fct_interfaces
    use init_ale_interface
    use init_thickness_ale_interface
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer                               :: i, n
    
    !___setup virt_salt_flux____________________________________________________
    ! if the ale thinkness remain unchanged (like in 'linfs' case) the vitrual 
    ! salinity flux need to be used
    ! otherwise we set the reference salinity to zero
    if ( .not. trim(which_ALE)=='linfs') then
        use_virt_salt=.false.
        ! this will force the virtual saltinity flux to be zero
        !!PS --> anyway ref_sss or rsss is not used when using zstar 
        !!PS ref_sss_local=.false.
        !!PS ref_sss      = 0.0_WP
        is_nonlinfs  = 1.0_WP
    else
        use_virt_salt=.true.
        is_nonlinfs  = 0.0_WP
    end if
   
    !___________________________________________________________________________
    ! initialize arrays for ALE
    if (partit%mype==0) then
       write(*,*) '____________________________________________________________'
       write(*,*) ' --> initialise ALE arrays + sparse SSH stiff matrix'
       write(*,*)
    end if
    
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_ale'//achar(27)//'[0m'
    call init_ale(dynamics, partit, mesh)
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_stiff_mat_ale'//achar(27)//'[0m'
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
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call oce_mixing_kpp_init'//achar(27)//'[0m'
        call oce_mixing_kpp_init(partit, mesh)
    ! initialise fesom1.4 like PP
    elseif (mix_scheme_nmb==2 .or. mix_scheme_nmb==27) then
    
    ! initialise cvmix_KPP
    elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_cvmix_kpp'//achar(27)//'[0m'
        call init_cvmix_kpp(partit, mesh)
        
    ! initialise cvmix_PP    
    elseif (mix_scheme_nmb==4 .or. mix_scheme_nmb==47) then
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_cvmix_pp'//achar(27)//'[0m'
        call init_cvmix_pp(partit, mesh)
        
    ! initialise cvmix_TKE    
    elseif (mix_scheme_nmb==5 .or. mix_scheme_nmb==56) then
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_cvmix_tke'//achar(27)//'[0m'
        call init_cvmix_tke(partit, mesh)
        
    endif
  
    ! initialise additional mixing cvmix_IDEMIX --> only in combination with 
    ! cvmix_TKE+cvmix_IDEMIX or stand alone for debbuging as cvmix_TKE
    if     (mod(mix_scheme_nmb,10)==6) then
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_cvmix_idemix'//achar(27)//'[0m'
        call init_cvmix_idemix(partit, mesh)
        
    ! initialise additional mixing cvmix_TIDAL --> only in combination with 
    ! KPP+cvmix_TIDAL, PP+cvmix_TIDAL, cvmix_KPP+cvmix_TIDAL, cvmix_PP+cvmix_TIDAL 
    ! or stand alone for debbuging as cvmix_TIDAL   
    elseif (mod(mix_scheme_nmb,10)==7) then
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_cvmix_tidal'//achar(27)//'[0m'
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
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call oce_adv_tra_fct_init'//achar(27)//'[0m'
    call oce_adv_tra_fct_init(tracers%work, partit, mesh)
    call muscl_adv_init(tracers%work, partit, mesh) !!PS test
    !=====================
    ! Initialize fields
    ! A user-defined routine has to be called here!
    !=====================
    if (toy_ocean) then  
       SELECT CASE (TRIM(which_toy))
         CASE ("soufflet") !forcing update for soufflet testcase
         if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call toy_channel'//achar(27)//'[0m'
           if (mod(mstep, soufflet_forc_update)==0) then
              call initial_state_soufflet(dynamics, tracers, partit, mesh)
              call compute_zonal_mean_ini(partit, mesh)  
              call compute_zonal_mean(dynamics, tracers, partit, mesh)
           end if
         CASE("nemo")
             call initial_state_nemo(dynamics, tracers, partit, mesh)
       END SELECT
    else
       if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call oce_initial_state'//achar(27)//'[0m' 
       call oce_initial_state(tracers, partit, mesh)   ! Use it if not running tests
    end if

    if (.not.r_restart) then
       do n=1, tracers%num_tracers
          do i=1, tracers%data(n)%AB_order-1
             tracers%data(n)%valuesold(i,:,:)=tracers%data(n)%values
          end do
       end do
    end if
    
    !___________________________________________________________________________
    ! first time fill up array for hnode & helem
    if (partit%mype==0) then
        write(*,*) '____________________________________________________________'
        write(*,*) ' --> call init_thickness_ale'
        write(*,*)
    end if
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_thickness_ale'//achar(27)//'[0m'
    call init_thickness_ale(dynamics, partit, mesh)
    
    !___________________________________________________________________________
    ! initialise arrays that are needed for backscatter_coef
    if(dynamics%opt_visc==8) call init_backscatter(partit, mesh)
        
    
    !___________________________________________________________________________
    if(partit%mype==0) write(*,*) 'Initial state'
    if (dynamics%use_wsplit .and. partit%mype==0) then
        write(*,*) '******************************************************************************'
        write(*,*) 'vertical velocity will be split onto explicit and implicit constitutes;'
        write(*,*) 'maximum allowed CDF on explicit W is set to: ', dynamics%wsplit_maxcfl
        write(*,*) '******************************************************************************'
    end if

end subroutine ocean_setup
!
!
!_______________________________________________________________________________
SUBROUTINE tracer_init(tracers, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    USE DIAGNOSTICS, only: ldiag_DVD
    USE g_ic3d
    use g_forcing_param, only: use_age_tracer !---age-code
    use g_config, only : lwiso, use_transit   ! add lwiso switch and switch for transient tracers
    use mod_transit, only : index_transit
    IMPLICIT NONE
    type(t_tracer), intent(inout), target               :: tracers
    type(t_partit), intent(inout), target               :: partit
    type(t_mesh),   intent(in) ,   target               :: mesh
    type(nml_tracer_list_type),    target, allocatable  :: nml_tracer_list(:)
    !___________________________________________________________________________
    integer        :: elem_size, node_size
    integer, save  :: nm_unit  = 104       ! unit to open namelist file, skip 100-102 for cray
    integer        :: iost
    integer        :: n
    !___________________________________________________________________________
    ! define tracer namelist parameter
    integer        :: num_tracers
    logical        :: i_vert_diff, smooth_bh_tra
    real(kind=WP)  :: gamma0_tra, gamma1_tra, gamma2_tra
    integer        :: AB_order
    namelist /tracer_listsize/ num_tracers
    namelist /tracer_list    / nml_tracer_list
    namelist /tracer_general / smooth_bh_tra, gamma0_tra, gamma1_tra, gamma2_tra, i_vert_diff, AB_order
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    
    !___________________________________________________________________________
    ! OPEN and read namelist for I/O
    open( unit=nm_unit, file='namelist.tra', form='formatted', access='sequential', status='old', iostat=iost )
    if (iost == 0) then
        if (mype==0) WRITE(*,*) '     file   : ', 'namelist.tra',' open ok'
        else
        if (mype==0) WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist.tra',' ; iostat=',iost
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    end if

    READ(nm_unit, nml=tracer_listsize, iostat=iost)
    allocate(nml_tracer_list(num_tracers))
    READ(nm_unit, nml=tracer_list,     iostat=iost)
    read(nm_unit, nml=tracer_init3d,   iostat=iost)
    READ(nm_unit, nml=tracer_general,  iostat=iost)
    close(nm_unit)

    do n=1, num_tracers
    if (nml_tracer_list(n)%id==-1) then
        if (mype==0) write(*,*) 'number of tracers will be changed from ', num_tracers, ' to ', n-1, '!'
        num_tracers=n-1
        EXIT
    end if
    end do

    !---wiso-code
    !=====================
    ! set necessary water isotope variables
    !=====================
    IF (lwiso) THEN
      ! always assume 3 water isotope tracers in the order H218O, HD16O, H216O
      ! tracers simulated in the model
      nml_tracer_list(num_tracers+1) = nml_tracer_list(1) ! use the same scheme as temperature
      nml_tracer_list(num_tracers+2) = nml_tracer_list(1)
      nml_tracer_list(num_tracers+3) = nml_tracer_list(1)

      nml_tracer_list(num_tracers+1)%id = 101
      nml_tracer_list(num_tracers+2)%id = 102
      nml_tracer_list(num_tracers+3)%id = 103

      index_wiso_tracers(1) = num_tracers+1
      index_wiso_tracers(2) = num_tracers+2
      index_wiso_tracers(3) = num_tracers+3

      num_tracers = num_tracers + 3

      ! tracers initialised from file
      idlist((n_ic3d+1):(n_ic3d+3)) = (/101, 102, 103/)
      filelist((n_ic3d+1):(n_ic3d+3)) = (/'wiso.nc', 'wiso.nc', 'wiso.nc'/)
      varlist((n_ic3d+1):(n_ic3d+3))  = (/'h2o18', 'hDo16', 'h2o16'/)

      n_ic3d = n_ic3d + 3

      if (mype==0) write(*,*) '3 water isotope tracers will be used in FESOM'
    END IF
    !---wiso-code-end

    !---age-code-begin
    if (use_age_tracer) then
      ! add age tracer in the model
      nml_tracer_list(num_tracers+1) = nml_tracer_list(1)
      nml_tracer_list(num_tracers+1)%id = 100
      index_age_tracer = num_tracers+1
      num_tracers = num_tracers + 1

      if (mype==0) write(*,*) '1 water age tracer will be used in FESOM'
    endif
    !---age-code-end

    ! Transient tracers
!! UNDER CONSTRUCTION - Actually we do not want to hardwire the number of transient tracers
    if (use_transit) then
      ! add transient tracers to the model
      nml_tracer_list(num_tracers+1) = nml_tracer_list(1)
      nml_tracer_list(num_tracers+2) = nml_tracer_list(1)
      nml_tracer_list(num_tracers+3) = nml_tracer_list(1)
      nml_tracer_list(num_tracers+4) = nml_tracer_list(1)
      nml_tracer_list(num_tracers+1)%id = 6
      nml_tracer_list(num_tracers+1)%id = 12
      nml_tracer_list(num_tracers+1)%id = 14
      nml_tracer_list(num_tracers+1)%id = 39

      index_transit(1) = num_tracers+1
      index_transit(2) = num_tracers+2
      index_transit(3) = num_tracers+3
      index_transit(4) = num_tracers+4

      num_tracers = num_tracers + 4

      ! tracers initialised from file
      idlist((n_ic3d+1):(n_ic3d+1)) = (/14/)
      filelist((n_ic3d+1):(n_ic3d+1)) = (/'R14C.nc'/)
      varlist((n_ic3d+1):(n_ic3d+1))  = (/'R14C'/)

      if (mype==0) write(*,*) 'XXX Transient tracers will be used in FESOM'
    endif
    ! 'use_transit' end


    if (mype==0) write(*,*) 'total number of tracers is: ', num_tracers

    !___________________________________________________________________________
    ! define local vertice & elem array size + number of tracers
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D+eDim_nod2D
    tracers%num_tracers=num_tracers

    !___________________________________________________________________________
    ! allocate/initialise horizontal velocity arrays in derived type
    ! Temperature (index=1), Salinity (index=2), etc.
    allocate(tracers%data(num_tracers))
    do n=1, tracers%num_tracers
        allocate(tracers%data(n)%values   (                             nl-1, node_size))
        allocate(tracers%data(n)%valuesAB (                             nl-1, node_size))
        tracers%data(n)%AB_order      = AB_order        
        allocate(tracers%data(n)%valuesold(tracers%data(n)%AB_order-1,  nl-1, node_size))
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
        tracers%data(n)%valuesold     = 0.
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
    USE o_param
    IMPLICIT NONE
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_dyn)   , intent(inout), target :: dynamics
    !___________________________________________________________________________
    integer        :: elem_size, node_size
    integer, save  :: nm_unit  = 105       ! unit to open namelist file, skip 100-102 for cray
    integer        :: iost
    !___________________________________________________________________________
    ! define dynamics namelist parameter
    integer        :: opt_visc
    real(kind=WP)  :: visc_gamma0, visc_gamma1, visc_gamma2
    real(kind=WP)  :: visc_easybsreturn
    logical        :: use_ivertvisc=.true.
    integer        :: momadv_opt
    logical        :: use_freeslip =.false.
    logical        :: use_wsplit   =.false.
    logical        :: ldiag_KE     =.false.
    integer        :: AB_order     = 2
    logical        :: check_opt_visc=.true.
    real(kind=WP)  :: wsplit_maxcfl
    namelist /dynamics_visc   / opt_visc, check_opt_visc, visc_gamma0, visc_gamma1, visc_gamma2,  &
                                use_ivertvisc, visc_easybsreturn
    namelist /dynamics_general/ momadv_opt, use_freeslip, use_wsplit, wsplit_maxcfl, ldiag_KE, AB_order
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"   
    
    !___________________________________________________________________________
    ! open and read namelist for I/O
    open(unit=nm_unit, file='namelist.dyn', form='formatted', access='sequential', status='old', iostat=iost )
    if (iost == 0) then
        if (mype==0) write(*,*) '     file   : ', 'namelist.dyn',' open ok'
    else
        if (mype==0) write(*,*) 'ERROR: --> bad opening file   : ', 'namelist.dyn',' ; iostat=',iost
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    end if
    read(nm_unit, nml=dynamics_visc,    iostat=iost)
    read(nm_unit, nml=dynamics_general, iostat=iost)
    close(nm_unit)

    !___________________________________________________________________________
    ! set parameters in derived type
    dynamics%opt_visc          = opt_visc
    dynamics%check_opt_visc    = check_opt_visc
    dynamics%visc_gamma0       = visc_gamma0
    dynamics%visc_gamma1       = visc_gamma1
    dynamics%visc_gamma2       = visc_gamma2
    dynamics%visc_easybsreturn = visc_easybsreturn
    dynamics%use_ivertvisc     = use_ivertvisc
    dynamics%momadv_opt        = momadv_opt
    dynamics%use_freeslip      = use_freeslip
    dynamics%use_wsplit        = use_wsplit
    dynamics%wsplit_maxcfl     = wsplit_maxcfl
    dynamics%ldiag_KE          = ldiag_KE
    dynamics%AB_order          = AB_order
    !___________________________________________________________________________
    ! define local vertice & elem array size
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D+eDim_nod2D
    
    !___________________________________________________________________________
    ! allocate/initialise horizontal velocity arrays in derived type
    allocate(dynamics%uv(        2, nl-1, elem_size))
    allocate(dynamics%uv_rhs(    2, nl-1, elem_size))
    allocate(dynamics%uv_rhsAB(  dynamics%AB_order-1, 2, nl-1, elem_size))
    allocate(dynamics%uvnode(    2, nl-1, node_size))
    dynamics%uv              = 0.0_WP
    dynamics%uv_rhs          = 0.0_WP
    dynamics%uv_rhsAB        = 0.0_WP
    dynamics%uvnode          = 0.0_WP
    if (Fer_GM) then
        allocate(dynamics%fer_uv(2, nl-1, elem_size))
        dynamics%fer_uv      = 0.0_WP
    end if 
    
    !___________________________________________________________________________
    ! allocate/initialise vertical velocity arrays in derived type
    allocate(dynamics%w(              nl, node_size))
    if (dynamics%ldiag_ke) then
       allocate(dynamics%w_old(       nl, node_size))
    end if
    allocate(dynamics%w_e(            nl, node_size))
    allocate(dynamics%w_i(            nl, node_size))
    allocate(dynamics%cfl_z(          nl, node_size))
    dynamics%w               = 0.0_WP
    dynamics%w_e             = 0.0_WP
    dynamics%w_i             = 0.0_WP
    dynamics%cfl_z           = 0.0_WP
    if (Fer_GM) then
        allocate(dynamics%fer_w(      nl, node_size))
        dynamics%fer_w       = 0.0_WP
    end if 
    
    !___________________________________________________________________________
    ! allocate/initialise ssh arrays in derived type
    allocate(dynamics%eta_n(      node_size))
    allocate(dynamics%d_eta(      node_size))
    allocate(dynamics%ssh_rhs(    node_size))
    dynamics%eta_n           = 0.0_WP
    dynamics%d_eta           = 0.0_WP
    dynamics%ssh_rhs         = 0.0_WP
    !!PS     allocate(dynamics%ssh_rhs_old(node_size))
    !!PS     dynamics%ssh_rhs_old= 0.0_WP   

    !___________________________________________________________________________
    ! inititalise working arrays
    allocate(dynamics%work%uvnode_rhs(2, nl-1, node_size))
    allocate(dynamics%work%u_c(nl-1, elem_size))
    allocate(dynamics%work%v_c(nl-1, elem_size))
    dynamics%work%uvnode_rhs = 0.0_WP
    dynamics%work%u_c = 0.0_WP
    dynamics%work%v_c = 0.0_WP
    if (dynamics%opt_visc==5) then
        allocate(dynamics%work%u_b(nl-1, elem_size))
        allocate(dynamics%work%v_b(nl-1, elem_size))
        dynamics%work%u_b = 0.0_WP
        dynamics%work%v_b = 0.0_WP
    end if 
   
    if (dynamics%ldiag_ke) then
       allocate(dynamics%ke_adv    (2, nl-1, elem_size))
       allocate(dynamics%ke_cor    (2, nl-1, elem_size))
       allocate(dynamics%ke_pre    (2, nl-1, elem_size))
       allocate(dynamics%ke_hvis   (2, nl-1, elem_size))
       allocate(dynamics%ke_vvis   (2, nl-1, elem_size))
       allocate(dynamics%ke_umean  (2, nl-1, elem_size))
       allocate(dynamics%ke_u2mean (2, nl-1, elem_size))
       allocate(dynamics%ke_du2    (2, nl-1, elem_size))
       allocate(dynamics%ke_adv_AB (dynamics%AB_order-1, 2, nl-1, elem_size))
       allocate(dynamics%ke_cor_AB (dynamics%AB_order-1, 2, nl-1, elem_size))
       allocate(dynamics%ke_rhs_bak(2, nl-1, elem_size))
       allocate(dynamics%ke_wrho   (nl-1, node_size))
       allocate(dynamics%ke_dW     (nl-1, node_size))
       allocate(dynamics%ke_Pfull  (nl-1, node_size))
       allocate(dynamics%ke_wind   (2, elem_size))
       allocate(dynamics%ke_drag   (2, elem_size))

       allocate(dynamics%ke_pre_xVEL (2, nl-1, elem_size))
       allocate(dynamics%ke_adv_xVEL (2, nl-1, elem_size))
       allocate(dynamics%ke_cor_xVEL (2, nl-1, elem_size))
       allocate(dynamics%ke_hvis_xVEL(2, nl-1, elem_size))
       allocate(dynamics%ke_vvis_xVEL(2, nl-1, elem_size))
       allocate(dynamics%ke_wind_xVEL(2, elem_size))
       allocate(dynamics%ke_drag_xVEL(2, elem_size))
       allocate(dynamics%ke_J(node_size),  dynamics%ke_D(node_size),   dynamics%ke_G(node_size),  &
                dynamics%ke_D2(node_size), dynamics%ke_n0(node_size),  dynamics%ke_JD(node_size), &
                dynamics%ke_GD(node_size), dynamics%ke_swA(node_size), dynamics%ke_swB(node_size))

       dynamics%ke_adv      =0.0_WP
       dynamics%ke_cor      =0.0_WP
       dynamics%ke_pre      =0.0_WP
       dynamics%ke_hvis     =0.0_WP
       dynamics%ke_vvis     =0.0_WP
       dynamics%ke_du2      =0.0_WP
       dynamics%ke_umean    =0.0_WP
       dynamics%ke_u2mean   =0.0_WP
       dynamics%ke_adv_AB   =0.0_WP
       dynamics%ke_cor_AB   =0.0_WP
       dynamics%ke_rhs_bak  =0.0_WP
       dynamics%ke_wrho     =0.0_WP
       dynamics%ke_wind     =0.0_WP
       dynamics%ke_drag     =0.0_WP
       dynamics%ke_pre_xVEL =0.0_WP
       dynamics%ke_adv_xVEL =0.0_WP
       dynamics%ke_cor_xVEL =0.0_WP
       dynamics%ke_hvis_xVEL=0.0_WP
       dynamics%ke_vvis_xVEL=0.0_WP
       dynamics%ke_wind_xVEL=0.0_WP
       dynamics%ke_drag_xVEL=0.0_WP
       dynamics%ke_dW       =0.0_WP
       dynamics%ke_Pfull    =0.0_WP
       dynamics%ke_J        =0.0_WP
       dynamics%ke_D        =0.0_WP
       dynamics%ke_G        =0.0_WP
       dynamics%ke_D2       =0.0_WP
       dynamics%ke_n0       =0.0_WP
       dynamics%ke_JD       =0.0_WP
       dynamics%ke_GD       =0.0_WP
       dynamics%ke_swA      =0.0_WP
       dynamics%ke_swB      =0.0_WP
    end if
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
    integer,        intent(in)            :: num_tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    !___________________________________________________________________________
    integer                               :: elem_size, node_size
    integer                               :: n, nt
    !___________________________________________________________________________
    ! define dynamics namelist parameter
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D+eDim_nod2D

    ! ================
    ! Velocities
    ! ================     
    !allocate(stress_diag(2, elem_size))!delete me
    !!PS allocate(Visc(nl-1, elem_size))
    allocate(t_star(node_size))
    allocate(qsr_c(node_size))
    ! ================
    ! elevation and its rhs
    ! ================

    ! ================
    ! Monin-Obukhov
    ! ================
    if (use_ice .and. use_momix) allocate(mo(nl,node_size),mixlength(node_size))
    if (use_ice .and. use_momix) mixlength=0.
    ! ================
    ! Vertical velocity and pressure
    ! ================
    allocate( hpressure(nl,node_size))
    allocate(bvfreq(nl,node_size),mixlay_dep(node_size),bv_ref(node_size))
    ! ================
    ! Ocean forcing arrays
    ! ================
    allocate(Tclim(nl-1,node_size), Sclim(nl-1, node_size))
    !---
    ! LA: add iceberg tracers 2023-02-08
    allocate(Tclim_ib(nl-1,node_size), Sclim_ib(nl-1, node_size))
    !---
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

    ! tracer gradients & RHS  
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

    allocate(MLD1(node_size), MLD2(node_size), MLD3(node_size))
    allocate(MLD1_ind(node_size), MLD2_ind(node_size), MLD3_ind(node_size))
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
    fer_gamma=0.0_WP
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

    !---wiso-code
    if (lwiso) then
      allocate(tr_arr_ice(node_size,3))  ! add sea ice tracers
      allocate(wiso_flux_oce(node_size,3))
      allocate(wiso_flux_ice(node_size,3))

      ! initialize sea ice isotopes with 0. permill
      ! absolute tracer values are increased by factor 1000. for numerical reasons
      ! (see also routine oce_fluxes in ice_oce_coupling.F90)
      do nt = 1,3
         tr_arr_ice(:,nt)=wiso_smow(nt) * 1000.0_WP
      end do

      ! initialize atmospheric fluxes over open ocean and sea ice
      wiso_flux_oce=0.0_WP
      wiso_flux_ice=0.0_WP
    end if
    !---wiso-code-end

END SUBROUTINE arrays_init
!
!
!_______________________________________________________________________________
! Here the 3D tracers will be initialized. Initialization strategy depends on a tracer ID.
! ID = 0 and 1 are reserved for temperature and salinity
! --> reads the initial state or the restart file for the ocean
SUBROUTINE oce_initial_state(tracers, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    USE o_ARRAYS
    USE g_config
    USE g_ic3d
    ! for additional (transient) tracers:
    use mod_transit, only: id_r14c, id_r39ar, id_f12, id_sf6
    implicit none
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in) ,   target :: mesh
    !___________________________________________________________________________
    integer                  :: i, k, counter, rcounter3, id
    character(len=10)        :: i_string, id_string
    real(kind=WP)            :: loc, max_temp, min_temp, max_salt, min_salt
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    if (mype==0) write(*,*) tracers%num_tracers, ' tracers will be used in FESOM'
    if (mype==0) write(*,*) 'tracer IDs are: ', tracers%data(1:tracers%num_tracers)%ID
    !
    ! read ocean state
    ! this must be always done! First two tracers with IDs 0 and 1 are the temperature and salinity.
    if(mype==0) write(*,*) 'read Temperature climatology from:', trim(filelist(1))
    if(mype==0) write(*,*) 'read Salinity    climatology from:', trim(filelist(2))
    if(any(idlist == 14) .and. mype==0) write(*,*) 'read radiocarbon climatology from:', trim(filelist(3))
    call do_ic3d(tracers, partit, mesh)
    
    Tclim=tracers%data(1)%values
    Sclim=tracers%data(2)%values
    Tsurf=Tclim(1,:)
    Ssurf=Sclim(1,:)
    
    if (use_icebergs) then
      Tclim_ib=tracers%data(1)%values
      Sclim_ib=tracers%data(2)%values
      Tsurf_ib=Tclim(1,:)
      Ssurf_ib=Sclim(1,:)
    end if
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
        !---age-code-begin
        ! FESOM tracers with code id 100 are used as water age
        CASE (100)
          if (mype==0) then
             write (i_string,  "(I3)") i
             write (id_string, "(I3)") id
             write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
             write (*,*) tracers%data(i)%values(1,1)
          end if
        !---age-code-end
        !---wiso-code
        ! FESOM tracers with code id 101, 102, 103 are used as water isotopes
        CASE (101)       ! initialize tracer ID=101 H218O
          if (mype==0) then
             write (i_string,  "(I3)") i
             write (id_string, "(I3)") id
             write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
             write (*,*) tracers%data(i)%values(1,1)
          end if
        CASE (102)       ! initialize tracer ID=102 HD16O
          if (mype==0) then
             write (i_string,  "(I3)") i
             write (id_string, "(I3)") id
             write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
             write (*,*) tracers%data(i)%values(1,1)
          end if
        CASE (103)       ! initialize tracer ID=103 H216O
          if (mype==0) then
             write (i_string,  "(I3)") i
             write (id_string, "(I3)") id
             write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
             write (*,*) tracers%data(i)%values(1,1)
          end if
        !---wiso-code-end

! Transient tracers
       CASE (14)        ! initialize tracer ID=14, fractionation-corrected 14C/C
!        this initialization can be overwritten by calling do_ic3d
!!         if (.not. any(idlist == 14)) then ! CHECK IF THIS LINE IS STILL NECESSARY
         tracers%data(i)%values(:,:) = 0.85
           if (mype==0) then
              write (i_string,  "(I3)") i
              write (id_string, "(I3)") id
              write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
              write (*,*) tracers%data(i)%values(1,1)
           end if
!!         end if
       CASE (39)        ! initialize tracer ID=39, fractionation-corrected 39Ar/Ar
         tracers%data(i)%values(:,:) = 0.85
         if (mype==0) then
            write (i_string,  "(I3)") i
            write (id_string, "(I3)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
            write (*,*) tracers%data(i)%values(1,1)
         end if
       CASE (12)        ! initialize tracer ID=12, CFC-12
         tracers%data(i)%values(:,:) = 0.
         if (mype==0) then
            write (i_string,  "(I3)") i
            write (id_string, "(I3)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
            write (*,*) tracers%data(i)%values(1,1)
         end if
       CASE (6)         ! initialize tracer ID=6, SF6
         tracers%data(i)%values(:,:) = 0.
         if (mype==0) then
            write (i_string,  "(I3)") i
            write (id_string, "(I3)") id
            write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
            write (*,*) tracers%data(i)%values(1,1)
         end if
! Transient tracers end

        !_______________________________________________________________________            
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
            
        !_______________________________________________________________________
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
            
        !_______________________________________________________________________            
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
            
        !_______________________________________________________________________
        CASE (501) ! ice-shelf water due to basal melting
            tracers%data(i)%values(:,:)=0.0_WP
            if (mype==0) then
                write (i_string,  "(I3)") i
                write (id_string, "(I3)") id
                write(*,*) 'initializing '//trim(i_string)//'th tracer with ID='//trim(id_string)
            end if
            
        !_______________________________________________________________________
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
    type(t_dyn)   , intent(inout), target  :: dynamics
    type(t_tracer), intent(inout), target  :: tracers
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    !___________________________________________________________________________
    integer                  :: i, k, counter, rcounter3, id
    character(len=10)        :: i_string, id_string
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    
    !___________________________________________________________________________
    if (toy_ocean) then
        SELECT CASE (TRIM(which_toy))
            CASE ("soufflet") !forcing update for soufflet testcase
            if (mod(mstep, soufflet_forc_update)==0) then
                call compute_zonal_mean(dynamics, tracers, partit, mesh)
            end if
        END SELECT
    end if
END SUBROUTINE before_oce_step
