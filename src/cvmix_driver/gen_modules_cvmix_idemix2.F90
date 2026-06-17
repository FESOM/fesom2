!
!
!===============================================================================
! module interface to FESOM2.0 for the CVMIX IDEMIX2 extension for the calculation 
! of the internal wave energy and its dissipation in Turbulent Kinetic Energy
! vertical mixing scheme --> it is based one its pyOM2 implementation by Carsten Eden
!
! @see Olbers D., Eden C.:
!       A Global Model for the Diapycnal Diffusivity Induced Internal Gravity Waves.
!       J. Phys. Oceanogr., 43, 1759-1779. doi: 10.1175/JPO-D-12-0207.1, 2013.
! @see Eden C., Czeschel L., Olbers D.:
!       Towards Energetically Consistent Ocean Models. 
!       J. Phys. Oceanogr., 44, 3160-3184, doi: 10.1175/JPO-D-13-0260.1, 2014.
! written by Patrick Scholz, 04.06.2026
module g_cvmix_idemix2
    
    !___________________________________________________________________________
    ! module calls from cvmix library
    use cvmix_idemix2, only : cvmix_idemix2_init                            , &
                              cvmix_idemix2_compute_param                   , & 
                              cvmix_idemix2_compute_compart_groupvel        , & 
                              cvmix_idemix2_compute_compart_interact_tscale , &
                              cvmix_idemix2_compute_M2_dissipation          , &
                              cvmix_idemix2_compute_vert_struct_fct         , & 
                              cvmix_idemix2_compute_vdiff_vdiss_Eiw         , &
                              cvmix_idemix2_compute_Eiw_waveinteract        , &
                              cvmix_idemix2_compute_Eiw_diss2KvAv
                                             
    use cvmix_put_get, only : cvmix_put
    use cvmix_kinds_and_types 
    
    !___________________________________________________________________________
    ! module calls from FESOM
    use g_config , only: dt, flag_debug, logfile_outfreq
    use o_param           
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_arrays
    use o_tracers
    use g_comm_auto 
    use g_interp
    use g_read_other_NetCDF
    use g_dist2coast, only: compute_dist2coast
    use find_up_downwind_triangles_interface
    use par_support_interfaces, only: init_mpi_types_fbin
    use g_support,              only: smooth_nod
    implicit none
    public
    
    !___________________________________________________________________________
    ! OCECTL/CVMIX_IDEMIX_PARAM namelist parameters
    ! time scale for vertical symmetrisation (sec)
    ! real(kind=WP) :: idemix_tau_v = 86400.0 ! old
    real(kind=WP)      :: idemix2_tau_v = 172800.0  ! from Pollman et al. (2017), use 2days, 2*86400sec
    
    ! time scale for horizontal symmetrisation, only necessary for lateral diffusion (sec), 
    ! use 15days, 15*86400sec
    real(kind=WP)      :: idemix2_tau_h = 1296000.0
    
    ! constant of order one derived from the shape of the spectrum in m space (dimensionless)
    real(kind=WP)      :: idemix2_gamma = 1.570
    
    ! spectral bandwidth in modes (dimensionless), or mode number scale, it describes over 
    ! over how many equivalent modes the energy is spread
    ! real(kind=WP) :: idemix_jstar = 10.0 ! old 
    real(kind=WP)      :: idemix2_jstar = 5.0    ! from Pollman et al. (2017)
    
    ! dissipation parameter (dimensionless)
    ! real(kind=WP) :: idemix_mu0   = 1.33333333 ! old 
    real(kind=WP)      :: idemix2_mu0   = 0.33333333   ! from Pollman et al. (2017), use 2days

    ! use superbee-spectral advection scheme and Adams-Bashfort timestepping 
!     logical            :: idemix2_enable_superbee_adv = .true.
    logical            :: idemix2_enable_AB  = .true.
    real(kind=WP)      :: idemix2_AB_epsilon = 0.1_WP
    ! number of spectral bins used for the M2 tidal and near-inertial wave (niw) components
    ! e.g 50+2 in loop used as do fbin=2,51 ... 1 and 52 serve as spectral boudary condition
    ! The first and last bins are used for boundary conditions in the spectral space.
    ! They ensure smooth transitions and prevent numerical instabilities at the edges of the spectral domain.
    ! Numerical Stability: The advection scheme (e.g., Superbee) requires ghost cells or boundary values.
    ! Skipping the first and last bins avoids out-of-bounds errors when calculating gradients.
    integer            :: idemix2_nfbin=52
    
    ! enable adding of idemix full horizontal tendency from 
    ! div(grad(Eiw*v0)*v0*tauh) diffusion term --> explicit
    logical            :: idemix2_enable_hor_diff_expl = .false.
    
    ! enable idemix1 functionality of homogenous diffusion into all directions
    logical            :: idemix2_enable_hor_diff_impl_iter = .false.
    integer            :: idemix2_hor_diff_niter       = 5   ! from Pollman et al. (2017)
    
    ! define shelf is defined as distance from coast (default 300km)
    real(kind=WP)      :: idemix2_shelf_dist = 300.0e3 
    
    ! scal down baroclinic wave speed when not all modes are used
    real(kind=WP)      :: idemix2_scal_cn = 1.0
    
    ! minimum cn floor: prevents cg_compart→0 in convective columns near the M2 critical latitude (~74.5°N).
    ! Recommended value: 0.1 m/s (symmetric with the cstar=max(1e-2,...) floor already applied for Eiw).
    real(kind=WP)      :: idemix2_cn_min  = 0.0
    
    !___________________________________________________________________________
    ! enable M2 tidal component as significant source of internal wave energy in IDEMIX2, 
    ! lower frequency modes
    ! Frequency: Semi-diurnal (twice daily, ~12.42 hour period)
    ! Source   : Barotropic tidal flow over rough topography
    ! Energy   : Higher energy, lower frequency
    ! Impact   : Dominant in regions with strong tidal currents
    logical            :: idemix2_enable_M2    = .true.
    
    ! filelocation for idemix2 M2 forcing (Summed anisotropic M2-tide generation modes 1-2 (W/m2))
    character(MAX_PATH):: idemix2_M2forc_file = './idemix2_IT_forc_M2modes1+2_aniso.nc'
    character(MAX_PATH):: idemix2_M2forc_vname= 'Flux_to2'
    character(MAX_PATH):: idemix2_M2forc_zname= 'phi'
    
    !___________________________________________________________________________
    ! enable Near-Inertial Wave (NIW) component of the IDEMIX module. NIWs are generated by 
    ! wind stress and have frequencies close to the local inertial frequency.
    ! Frequency: Near the local inertial frequency (f, typically 1-4 cycles per day at mid-latitudes)
    ! Source   : Wind forcing, especially during storms
    ! Energy   : Lower energy, higher frequency than M2
    ! Impact   : Important in the upper ocean, especially after storms
    logical            :: idemix2_enable_niw   = .false.
    
    ! amount of surface for2cing that is used
    real(kind=WP)      :: idemix2_fniw_usage   = 0.2
    
    ! filelocation for idemix2 NIW forcing (near inertial waves)
    character(MAX_PATH):: idemix2_niwforc_file = './dummy.nc'
    character(MAX_PATH):: idemix2_niwforc_vname= 'dummy'
    
    !___________________________________________________________________________
    ! filelocation for idemix2 bottom forcing
    logical            :: idemix2_enable_bot   = .true.
    
    character(MAX_PATH):: idemix2_botforc_file = './tidal_energy_gx1v6_20090205_rgrid.nc'
    character(MAX_PATH):: idemix2_botforc_vname= 'wave_dissipation'
    ! total global Energy input that should be conserved if 0.0 no conservation is applied
    real(kind=WP)      :: idemix2_botforc_Etot = 0.0_WP ! units W
    
    !___________________________________________________________________________
    ! maximum cutoff value for interal wave energy --> is there for stability
    real(kind=WP)      :: idemix2_Eiw_maxthresh= 0.1
    ! Option A: smooth IDEMIX2 internal fields to improve long-term stability.
    ! All flags default to .false. so Option B (TKE-interface smoothing in
    ! gen_modules_cvmix_tke.F90) is tested first.
    logical            :: idemix2_smooth_Eiw_diss = .false.   ! smooth iwe2_E_iw_diss after vertical solve
    logical            :: idemix2_smooth_Eiw      = .false.   ! smooth iwe2_E_iw(:,:,tip1) after wave-wave
    logical            :: idemix2_smooth_alpha_c  = .false.   ! smooth iwe2_alpha_c after parameter loop
    integer            :: idemix2_smooth_niter    = 1
    
    !___________________________________________________________________________
    ! switch on extended diagnostic
    logical            :: idemix2_diag_Ecompart= .false.
    logical            :: idemix2_diag_Eiw     = .false.
    logical            :: idemix2_diag_WWI     = .false.
    
!     !___________________________________________________________________________
!     ! enable lee wave source of internal wave energy in IDEMIX2, 
!     ! Lee Wave Generation: Generated when mean flow encounters topography Energy 
!     ! transfers from mean flow to internal waves
!     logical            :: idemix2_enable_leew   = .false.
!     
!     ! filelocation for idemix2 M2 forcing (Summed anisotropic M2-tide generation modes 1-2 (W/m2))
!     character(MAX_PATH):: idemix2_leewforc_file = './idemix2_lee_forc_Eden.nc'
!     character(MAX_PATH):: idemix2_leewforc_vname= 'C_lee'
    
    !___________________________________________________________________________
    ! filelocation for idemix2 Root Mean Square Topographic Height forcing. Measures the roughness or 
    ! variability of seafloor topography. Units: Meters (m), Role: Represents the standard deviation 
    ! of seafloor height variations. Used to calculate energy transfer from tides/NIW to internal waves.
    character(MAX_PATH):: idemix2_hrmsforc_file = './idemix2_forcing_t-scattering_Goff2023_1deg.nc'
    character(MAX_PATH):: idemix2_hrmsforc_vname= 'HRMS'
    
    ! filelocation for idemix2 Characteristic horizontal length scale of topographic features.
    ! Units: Meters (m), Role: Represents the dominant wavelength of seafloor roughness. Used to 
    ! normalize the topographic forcing.
    character(MAX_PATH):: idemix2_hlamforc_file = './idemix2_forcing_t-scattering_1deg.nc'
    character(MAX_PATH):: idemix2_hlamforc_vname= 'LAMBDA_G10'
    
    namelist /param_idemix2/ idemix2_tau_v, idemix2_tau_h, idemix2_gamma, idemix2_jstar, idemix2_mu0, idemix2_scal_cn, idemix2_cn_min, &
                             idemix2_enable_AB, idemix2_AB_epsilon, idemix2_nfbin, & ! idemix2_enable_superbee_adv
                             idemix2_enable_hor_diff_expl, idemix2_enable_hor_diff_impl_iter, idemix2_hor_diff_niter, &
                             idemix2_shelf_dist, &
                             idemix2_botforc_Etot, &
                             idemix2_enable_M2  , idemix2_M2forc_file  , idemix2_M2forc_vname  , idemix2_M2forc_zname, &
                             idemix2_enable_niw , idemix2_niwforc_file , idemix2_niwforc_vname , idemix2_fniw_usage, &
!                              idemix2_enable_leew, idemix2_leewforc_file, idemix2_leewforc_vname, &
                             idemix2_enable_bot , idemix2_botforc_file , idemix2_botforc_vname , &
                             idemix2_hrmsforc_file, idemix2_hrmsforc_vname, &
                             idemix2_hlamforc_file, idemix2_hlamforc_vname, &
                             idemix2_Eiw_maxthresh, &
                             idemix2_smooth_Eiw_diss, idemix2_smooth_Eiw, idemix2_smooth_alpha_c, idemix2_smooth_niter, &
                             idemix2_diag_Ecompart, idemix2_diag_Eiw, idemix2_diag_WWI
    
    !___________________________________________________________________________
    ! CVMIX-IDEMIX variables
    real(kind=WP), allocatable, dimension(:)    :: iwe2_phit, iwe2_phiu, iwe2_dphit, iwe2_dphiu
    
    !
    ! --- M2 related global variables ---
    real(kind=WP)                               :: iwe2_omega_M2 
    real(kind=WP), allocatable, dimension(:,:,:):: iwe2_M2_uv, iwe2_E_M2, iwe2_E_M2_divh, iwe2_E_M2_divs
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_M2_w, iwe2_fM2, iwe2_E_M2_struct 
    real(kind=WP), allocatable, dimension(:)    :: w_M2_e, iwe2_alpha_M2_c, iwe2_M2_tau    
    
    ! optional diagnostic
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_E_M2_dt  ,                     &
                                                   iwe2_E_M2_advh, iwe2_E_M2_advs    , &
                                                   iwe2_E_M2_diss, iwe2_E_M2_diss_wwi, &
                                                   iwe2_E_M2_forc,                     &
                                                   iwe2_E_M2_refl
    
    
    ! --- niw related global variables --- 
    real(kind=WP), allocatable, dimension(:)    :: iwe2_omega_niw
    real(kind=WP), allocatable, dimension(:,:,:):: iwe2_niw_uv, iwe2_E_niw, iwe2_E_niw_divh, iwe2_E_niw_divs
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_niw_w, iwe2_fniw, iwe2_E_niw_struct
    real(kind=WP), allocatable, dimension(:)    :: w_niw_e, iwe2_niw_tau
    
    ! optional diagnostic
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_E_niw_dt  ,                      & 
                                                   iwe2_E_niw_advh, iwe2_E_niw_advs    , &
                                                   iwe2_E_niw_diss,                      &
                                                   iwe2_E_niw_forc,                      &
                                                   iwe2_E_niw_refl
    
    ! --- Eiw - internal wave energy related variables ---
    real(kind=WP), allocatable, dimension(:,:,:):: iwe2_E_iw
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_E_iw_diss
    real(kind=WP), allocatable, dimension(:)    :: iwe_E_iw_vint
    
    ! optional diagnostic
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_E_iw_dt, iwe2_E_iw_fbot, & 
                                                   iwe2_E_iw_vdif, iwe2_E_iw_hdif
    real(kind=WP), allocatable, dimension(:)    :: iwe2_E_iw_fsrf
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_E_iw_diss_M2, iwe2_E_iw_diss_niw
    
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_Av
    
    ! --- general idemix variable ---
    real(kind=WP), allocatable, dimension(:)    :: iwe2_cn
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_c0, iwe2_v0, iwe2_alpha_c
    
    ! --- forcing realted variables ---
    real(kind=WP), allocatable, dimension(:)    :: iwe2_topo_hrms, iwe2_topo_hlam, iwe2_topo_dist
    real(kind=WP), allocatable, dimension(:)    :: iwe2_fleew, iwe2_fsrf, iwe2_fbot_n
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_fbot
    
    ! --- support variables ---
    real(kind=WP), allocatable, dimension(:)    :: vol_nodB2T
    real(kind=WP), allocatable, dimension(:,:)  :: vol_wcelli
    real(kind=WP), allocatable, dimension(:)    :: iwe2_grady_coriol, aux
    real(kind=WP), allocatable, dimension(:,:,:):: iwe2_gradxy_e, iwe2_gradxy_n
    real(kind=WP), allocatable, dimension(:,:  ):: iwe2_flx_uv, iwe2_flx_w
    integer                                     :: iwe2_ti=1, iwe2_tip1=2

    ! reflective coastal BC: for each edge touching exactly one coast node,
    ! iwe2_bc_reflect_bin stores the mirror spectral bin (0 = edge not flagged),
    ! iwe2_bc_reflect_sgn = +1 if ednodes(1) is the coast node, -1 if ednodes(2) is.
    integer      , allocatable, dimension(:)    :: iwe2_refl_bin
    integer      , allocatable, dimension(:)    :: iwe2_refl_sgn
    ! per-timestep scratch: blocked coast flux accumulated per (mirror bin, interior node),
    ! filled by apply_reflect_bc_spctrl, applied to Edivh in hsintegrate_Ecompart
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_refl_src
    ! mask: .true. for own nodes that are interior to at least one coast edge
    logical,       allocatable, dimension(:)    :: iwe2_refl_node
    ! pre-computed periodic wrap indices for cross-spectral superbee stencil
    ! (branch-free inner loop in adv_Ecompart_crss_spctrl_superbee)
    integer,       allocatable, dimension(:)    :: iwe2_idxp1, iwe2_idxp2, iwe2_idxm1
        
    ! load variables from CVMix list
    type(cvmix_data_type)                       :: CVMix_vars

    contains
    
    
    
    !
    !
    !
    !===========================================================================
    ! allocate and initialize IDEMIX variables --> call initialisation 
    ! routine from cvmix library
    subroutine init_cvmix_idemix2(partit, mesh)
        implicit none
        type(t_mesh),   intent(inout), target :: mesh
        type(t_partit), intent(inout), target :: partit
        
        character(len=cvmix_strlen) :: nmlfile
        logical                     :: file_exist=.False.
        integer                     :: node_size, elem_size, elem, node, nfbin, &
                                       fbin_i, elnodes(3), nzmax, nlu, nln, nz
        real(kind=WP)               :: loc_Etot=0.0_WP, glb_Etot=0.0_WP
        real(kind=WP)               :: t0, t1
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        !_______________________________________________________________________
        if(mype==0) then
            write(*,*) '____________________________________________________________'
            write(*,*) ' --> initialise IDEMIX2'
            write(*,*)
        end if
                !_______________________________________________________________________
        ! read cvmix namelist file 
        nmlfile ='namelist.cvmix'    ! name of ocean namelist file
        ! check if cvmix namelist file exists if not use default values 
        file_exist=.False.
        inquire(file=trim(nmlfile),exist=file_exist) 
        if (file_exist) then
            open(20,file=trim(nmlfile))
                read(20,nml=param_idemix2)
            close(20)
        else
            write(*,*) '     could not find namelist.cvmix, will use default values !'    
        end if    
        
        !_______________________________________________________________________
        if (mype==0) then
            write(*,*) "     IDEMIX2:"
            write(*,*) "     ├> parameter                     "
            write(*,*) "     │  ┌──────────────────────────── "
            write(*,*) "     │  ├> idemix2_tau_v            = ", idemix2_tau_v
            write(*,*) "     │  ├> idemix2_tau_h            = ", idemix2_tau_h
            write(*,*) "     │  ├> idemix2_gamma            = ", idemix2_gamma
            write(*,*) "     │  ├> idemix2_jstar            = ", idemix2_jstar
            write(*,*) "     │  ├> idemix2_mu0              = ", idemix2_mu0
            write(*,*) "     │  ├> idemix2_scal_cn          = ", idemix2_scal_cn
            write(*,*) "     │  ├> idemix2_cn_min           = ", idemix2_cn_min
            write(*,*) "     │  ├> idemix2_Eiw_maxthresh    = ", idemix2_Eiw_maxthresh
            write(*,*) "     │  ├> idemix2...hor_diff_expl  = ", idemix2_enable_hor_diff_expl
            write(*,*) "     │  └> idemix2...hor_diff_impl_iter= ", idemix2_enable_hor_diff_impl_iter
            write(*,*) "     │     └> idemix2_hor_diff_niter= ", idemix2_hor_diff_niter
            write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_AB_timestep         = ", idemix2_enable_AB
            write(*,*) "     │  └> idemix2_AB_epsilon       = ", idemix2_AB_epsilon
            write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_nfbin               = ", idemix2_nfbin
            write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_enable_M2           = ", idemix2_enable_M2
            write(*,*) "     │  ┌──────────────────────────── "
            write(*,*) "     │  └> idemix2_M2forc_file      = ", trim(idemix2_M2forc_file)
            write(*,*) "     │     ├> idemix2_M2forc_vname  = ", trim(idemix2_M2forc_vname)
            write(*,*) "     │     └> idemix2_M2forc_zname  = ", trim(idemix2_M2forc_zname)
            write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_enable_niw          = ", idemix2_enable_niw 
            write(*,*) "     │  ┌──────────────────────────── "
            write(*,*) "     │  ├> idemix2_fniw_usage       = ", idemix2_fniw_usage 
            write(*,*) "     │  └> idemix2_niwforc_file     = ", trim(idemix2_niwforc_file)
            write(*,*) "     │     └> idemix2_niwforc_vname = ", trim(idemix2_niwforc_vname)
            write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_enable_bot          = ", idemix2_enable_bot
            write(*,*) "     │  ┌──────────────────────────── "
            write(*,*) "     │  └> idemix2_botforc_file     = ", trim(idemix2_botforc_file)
            write(*,*) "     │     ├> idemix2_botforc_vname = ", trim(idemix2_botforc_vname)
            write(*,*) "     │     └> idemix2_botforc_Etot  = ", idemix2_botforc_Etot
            write(*,*) "     │                                "
!             write(*,*) "     ├> idemix2_enable_leew         = ", idemix2_enable_leew
!             write(*,*) "     │  └> idemix2_leewforc_file    = ", trim(idemix2_leewforc_file)
!             write(*,*) "     │     └> idemix2_leewforc_vname= ", trim(idemix2_leewforc_vname)
!             write(*,*) "     │                                "
            write(*,*) "     ├> topographic height forcing    "
            write(*,*) "     │  ┌──────────────────────────── "
            write(*,*) "     │  ├> idemix2_hrmsforc_file    = ", trim(idemix2_hrmsforc_file)
            write(*,*) "     │  │  └> idemix2_hrmsforc_vname= ", trim(idemix2_hrmsforc_vname)
            write(*,*) "     │  └> idemix2_hlamforc_file    = ", trim(idemix2_hlamforc_file)
            write(*,*) "     │     └> idemix2_hlamforc_vname= ", trim(idemix2_hlamforc_vname)
            write(*,*) "     │                                "
            WRITE(*,*) "     ├> idemix2_shelf_dist          = ", idemix2_shelf_dist
            write(*,*) "     │                                "
            write(*,*) "     ├> smoothing switches            "
            write(*,*) "     │  ┌──────────────────────────── "
            write(*,*) "     │  ├> idemix2_smooth_Eiw_diss  = ", idemix2_smooth_Eiw_diss 
            write(*,*) "     │  ├> idemix2_smooth_Eiw       = ", idemix2_smooth_Eiw
            write(*,*) "     │  ├> idemix2_smooth_alpha_c   = ", idemix2_smooth_alpha_c
            write(*,*) "     │  └> idemix2_smooth_niter     = ", idemix2_smooth_niter
            write(*,*) "     │                                "
            write(*,*) "     ├> diagnostic switches           "
            write(*,*) "     │  ┌──────────────────────────── "
            write(*,*) "     │  ├> idemix2_diag_Ecompart    = ", idemix2_diag_Ecompart 
            write(*,*) "     │  ├> idemix2_diag_Eiw         = ", idemix2_diag_Eiw
            write(*,*) "     │  └> idemix2_diag_WWI         = ", idemix2_diag_WWI
            write(*,*)
            write(*,*) "     IDEMIX2 inputs:"
        end if
        
        !_______________________________________________________________________
        ! allocate + initialse idemix arrays --> with size myDim_nod2D+eDim_nod2D
        node_size=myDim_nod2D+eDim_nod2D
        elem_size=myDim_elem2D+eDim_elem2D
        nfbin    =idemix2_nfbin
        
        ! Build MPI datatypes for spectral bin halo exchange
        call init_mpi_types_fbin(idemix2_nfbin, partit)
        
        ! 3d vertical group velocities (m/s) for the continuous internal wave spectrum
        ! 3d horizontal group velocities (m/s) for the continuous internal wave spectrum
        allocate(iwe2_c0(nl,node_size), iwe2_v0(nl,node_size))
        iwe2_c0(:,:)         = 0.0_WP
        iwe2_v0(:,:)         = 0.0_WP
        
        ! baroclionic gravity wave speed
        allocate(iwe2_cn(node_size))
        iwe2_cn(:)           = 0.0_WP
        
        ! Eiw 3d enery dissipation coefficient 
        allocate(iwe2_alpha_c(nl,node_size))
        iwe2_alpha_c(:,:)    = 0.0_WP
        
        ! initialise Eiw - internal wave energy variables
        ! index (..., 1:3) timestep index E^(n-1), E^(n), E^(n+1)
        allocate(iwe2_E_iw(   nl, node_size, 2))
        iwe2_E_iw(     :,:,:)= 0.0_WP
        
        allocate(iwe2_E_iw_diss(nl, node_size))
        iwe2_E_iw_diss(  :,:)= 0.0_WP ! Eiw production from disspation 
        
        ! Diagnostics
        if (idemix2_diag_Eiw) then
            allocate( iwe2_E_iw_fbot(nl, node_size) &
                    , iwe2_E_iw_vdif(nl, node_size) &
                    , iwe2_E_iw_hdif(nl, node_size) &
                    , iwe2_E_iw_fsrf(    node_size) &
                    )
            iwe2_E_iw_fbot(  :,:)= 0.0_WP
            iwe2_E_iw_hdif(  :,:)= 0.0_WP
            iwe2_E_iw_vdif(  :,:)= 0.0_WP
            iwe2_E_iw_fsrf(    :)= 0.0_WP
        end if
        if (idemix2_diag_Eiw .or. idemix2_diag_WWI) then
            allocate( iwe2_E_iw_dt(  nl, node_size))
            iwe2_E_iw_dt(    :,:)= 0.0_WP
        end if
        
        ! initialise M2 variables
        if (idemix2_enable_M2) then
            ! M2 energy dissipation
            allocate(iwe2_alpha_M2_c(node_size))
            iwe2_alpha_M2_c(:)    = 0.0_WP
            
            ! M2 dissipation timescale
            allocate(iwe2_M2_tau(node_size))
            iwe2_M2_tau(:)        = 0.0_WP
            
            ! M2 cross spectral propagation at elements
            allocate(w_M2_e(nfbin))
            w_M2_e(:)             = 0.0_WP
            
            ! M2 zonal/merid group velocity, and cross spectral propagation
            ! at vertices
            allocate(iwe2_M2_uv(2, nfbin,elem_size))
            allocate(iwe2_M2_w(    nfbin,node_size))
            iwe2_M2_uv(:,:,:)     = 0.0_WP
            iwe2_M2_w(:,:)        = 0.0_WP
            
            ! M2 forcing
            allocate(iwe2_fM2(nfbin, node_size))
            iwe2_fM2(:,:)         = 0.0_WP
            
            ! M2 wave energy, and divergence of M2 wave energy
            ! index (..., 1:3) timestep index E^(n-1), E^(n), E^(n+1)
            allocate(  iwe2_E_M2(     nfbin, node_size, 2) &
                     , iwe2_E_M2_divh(nfbin, node_size, 2) &
                     , iwe2_E_M2_divs(nfbin, node_size, 2) &
                    )
            iwe2_E_M2(      :,:,:)= 0.0_WP
            iwe2_E_M2_divh( :,:,:)= 0.0_WP
            iwe2_E_M2_divs( :,:,:)= 0.0_WP
            
            ! structure function for M2 energy
            allocate(iwe2_E_M2_struct(nl, node_size))
            iwe2_E_M2_struct(:,:) = 0.0_WP
            
            ! Diagnostics
            if (idemix2_diag_Ecompart) then
                ! diagnostic for M2 spectral energy advection
                allocate( iwe2_E_M2_advh(nfbin, node_size) &
                        , iwe2_E_M2_advs(nfbin, node_size) &
                        , iwe2_E_M2_diss(nfbin, node_size) &
                        , iwe2_E_M2_forc(nfbin, node_size) &
                        , iwe2_E_M2_refl(nfbin, node_size) )
                iwe2_E_M2_advh(:,:) = 0.0_WP
                iwe2_E_M2_advs(:,:) = 0.0_WP
                iwe2_E_M2_diss(:,:) = 0.0_WP
                iwe2_E_M2_forc(:,:) = 0.0_WP
                iwe2_E_M2_refl(:,:) = 0.0_WP
            end if
            
            if (idemix2_diag_Ecompart .or. idemix2_diag_WWI) then
                ! diagnostic for M2 spectral energy advection
                allocate( iwe2_E_M2_dt(  nfbin, node_size) )
                iwe2_E_M2_dt(  :,:) = 0.0_WP
            end if 
            
            if (idemix2_diag_WWI) then
                ! diagnostic for M2 wave-wave interaction
                allocate( iwe2_E_iw_diss_M2( nl   , node_size   ) &
                        , iwe2_E_M2_diss_wwi(nfbin, node_size) )
                iwe2_E_iw_diss_M2( :,:) = 0.0_WP
                iwe2_E_M2_diss_wwi(:,:) = 0.0_WP
            end if
        end if

        ! initialise niw variables
        if (idemix2_enable_niw) then
            ! niw frequency
            allocate(iwe2_omega_niw(node_size))
            iwe2_omega_niw(:)      = 0.0_WP
            
            ! niw dissipation timescale
            allocate(iwe2_niw_tau(node_size))
            iwe2_niw_tau(:)        = 0.0_WP
            
            ! niw cross spectral propagation at elements
            allocate(w_niw_e(nfbin))
            w_niw_e(:)             = 0.0_WP
            
            ! niw zonal/merid group velocity, and cross spectral propagation
            ! at vertices
            allocate(iwe2_niw_uv(2, nfbin,elem_size))
            allocate(iwe2_niw_w(    nfbin,node_size))
            iwe2_niw_uv(:,:,:)     = 0.0_WP
            iwe2_niw_w(:,:)        = 0.0_WP
            
            ! niw forcing
            allocate(iwe2_fniw(nfbin, node_size))
            iwe2_fniw(:,:)         = 0.0_WP
            
            ! niw wave energy, and divergence of niw wave energy
            ! index (..., 1:3) timestep index E^(n-1), E^(n), E^(n+1)
            allocate(  iwe2_E_niw(     nfbin, node_size, 2) &
                     , iwe2_E_niw_divh(nfbin, node_size, 2) &
                     , iwe2_E_niw_divs(nfbin, node_size, 2) &
                     )
            iwe2_E_niw(      :,:,:)= 0.0_WP
            iwe2_E_niw_divh( :,:,:)= 0.0_WP
            iwe2_E_niw_divs( :,:,:)= 0.0_WP
            
            ! structure function for niw energy
            allocate(iwe2_E_niw_struct(nl, node_size))
            iwe2_E_niw_struct(:,:) = 0.0_WP
            
            ! Diagnostics
            if (idemix2_diag_Ecompart) then
                ! diagnostic for niw spectral energy advection
                allocate( iwe2_E_niw_advh(nfbin, node_size) &
                        , iwe2_E_niw_advs(nfbin, node_size) &
                        , iwe2_E_niw_diss(nfbin, node_size) &
                        , iwe2_E_niw_forc(nfbin, node_size) &
                        , iwe2_E_niw_refl(nfbin, node_size) )
                iwe2_E_niw_advh(:,:) = 0.0_WP
                iwe2_E_niw_advs(:,:) = 0.0_WP
                iwe2_E_niw_diss(:,:) = 0.0_WP
                iwe2_E_niw_forc(:,:) = 0.0_WP
                iwe2_E_niw_refl(:,:) = 0.0_WP
            end if
            
            if (idemix2_diag_Ecompart .or. idemix2_diag_WWI) then
                ! diagnostic for niw spectral energy advection
                allocate( iwe2_E_niw_dt(  nfbin, node_size))
                iwe2_E_niw_dt(  :,:) = 0.0_WP
            end if
            
            if (idemix2_diag_WWI) then
                ! diagnostic for niw wave-wave interaction
                allocate( iwe2_E_iw_diss_niw(nl, node_size) )
                iwe2_E_iw_diss_niw(:,:) = 0.0_WP
            end if
        end if 
        
        ! reflected boundary condition
        allocate(iwe2_refl_bin(myDim_edge2D))
        allocate(iwe2_refl_sgn(myDim_edge2D))
        allocate(iwe2_refl_src(idemix2_nfbin, myDim_nod2D))
        allocate(iwe2_refl_node(myDim_nod2D+eDim_nod2D))
        iwe2_refl_bin(:)  = 0
        iwe2_refl_sgn(:)  = 0
        iwe2_refl_src(:,:)= 0.0_WP
        iwe2_refl_node(:) = .false.
        
        ! allocate 1d Spectral space coordinates
        allocate(  iwe2_phit( nfbin)              &  
                 , iwe2_phiu( nfbin)              &  
                 , iwe2_dphit(nfbin)              &
                 , iwe2_dphiu(nfbin)              &
                )
        iwe2_phit(:)         = 0.0_WP
        iwe2_phiu(:)         = 0.0_WP
        iwe2_dphit(:)        = 0.0_WP
        iwe2_dphiu(:)        = 0.0_WP
                
        ! internal wave related vertical viscosity and diffusivity
        if(mix_scheme_nmb==7) then 
            allocate(iwe2_Av(nl,elem_size))
            iwe2_Av(:,:)     = 0.0_WP
        endif 
        
        ! forcing fields, M2 tidal forcing (spectral) and NIW forcing (spectral)
        allocate(  iwe2_fbot_n(    node_size)     &
                 , iwe2_fbot(nl, node_size)       &
                 , iwe2_fsrf(    node_size)       &
!                  , iwe2_fleew(   myDim_elem2D)    &
                 )
        iwe2_fbot_n(:)       = 0.0_WP
        iwe2_fbot(:,:)       = 0.0_WP
        iwe2_fsrf(  :)       = 0.0_WP
!         iwe2_fleew( :)       = 0.0_WP
        
        ! Topographic height and Topographic wavelength
        allocate(  iwe2_topo_hrms( node_size)   &
                 , iwe2_topo_hlam( node_size)   & 
                 , iwe2_topo_dist( node_size)   &
                )
        iwe2_topo_hrms(:)    = 0.0_WP
        iwe2_topo_hlam(:)    = 0.0_WP
        iwe2_topo_dist(:)    = 0.0_WP
        
        ! support gradient variables 
        allocate(  iwe2_grady_coriol(      elem_size)   &   
                 , iwe2_gradxy_e(2, nfbin, elem_size)   &
                 , iwe2_gradxy_n(2, nfbin, node_size)   &
                 )
        iwe2_grady_coriol(:) = 0.0_WP
        iwe2_gradxy_e(:,:,:) = 0.0_WP
        iwe2_gradxy_n(:,:,:) = 0.0_WP
        
        ! support horizontal edge flux and vertical flux and accumulated divergence variables 
        allocate(  iwe2_flx_uv(nfbin, partit%myDim_edge2D) & 
                 , iwe2_flx_w( nfbin, partit%myDim_nod2D ) & 
                 )
        iwe2_flx_uv(:,:)     = 0.0
        iwe2_flx_w( :,:)     = 0.0
        
        ! support inverse volume of wcell
        allocate(vol_wcelli(nl, node_size))
        vol_wcelli(:,:)     = 0.0_WP
        allocate(vol_nodB2T(node_size))
        vol_nodB2T(:)     = 0.0_WP
        
        !_______________________________________________________________________
        ! width of spectral frequency bins
        iwe2_dphit(:)   = 2.0_WP*pi/(nfbin-2)
        iwe2_dphiu(:)   = iwe2_dphit(:)
        ! iwe2_phiu(k) represents the center of the k-th spectral bin
        ! iwe2_phit(k) represents the edge of the k-th spectral bin
        iwe2_phit(1) = -iwe2_dphit(1)
        iwe2_phiu(1 )= iwe2_phit(1)+iwe2_dphit(1)/2.0_WP
        do fbin_i=2,nfbin
            iwe2_phit(fbin_i) = iwe2_phit(fbin_i-1) + iwe2_dphit(fbin_i)
            iwe2_phiu(fbin_i) = iwe2_phiu(fbin_i-1) + iwe2_dphiu(fbin_i)
        end do

        !_______________________________________________________________________
        ! pre-compute periodic wrap indices for cross-spectral superbee stencil
        allocate(iwe2_idxp1(nfbin), iwe2_idxp2(nfbin), iwe2_idxm1(nfbin))
        do fbin_i = 1, nfbin
            iwe2_idxp1(fbin_i) = fbin_i+1; if (iwe2_idxp1(fbin_i)>nfbin) iwe2_idxp1(fbin_i)=iwe2_idxp1(fbin_i)-(nfbin-2)
            iwe2_idxp2(fbin_i) = fbin_i+2; if (iwe2_idxp2(fbin_i)>nfbin) iwe2_idxp2(fbin_i)=iwe2_idxp2(fbin_i)-(nfbin-2)
            iwe2_idxm1(fbin_i) = fbin_i-1; if (iwe2_idxm1(fbin_i)<1    ) iwe2_idxm1(fbin_i)=nfbin-2
        end do

        !_______________________________________________________________________
        ! pre-compute reflective coastal BC lookup tables
        call init_reflect_bc( iwe2_refl_bin  & ! OUT: mirror spectral bin index per coast edge (0=not coast)
                            , iwe2_refl_sgn  & ! OUT: sign per coast edge (+1 if ednodes(1) is coast, -1 otherwise)
                            , iwe2_refl_src  & ! OUT: per-timestep reflect source accumulator (m²/s³)
                            , iwe2_refl_node & ! OUT: logical mask, .true. for coast-adjacent nodes
                            , partit         & ! IN
                            , mesh           & ! IN
                            )

        !_______________________________________________________________________
        ! read idemix M2 forcing from cfsr data --> file
        if (idemix2_enable_M2) then
            t0=MPI_Wtime()
            ! omega_M2 is fixed, as the M2 tidal forcing is independent of the local inertial frequency.
            ! Physical Implication: At high latitudes where |f| > ω_M2, the M2 tide becomes evanescent 
            ! (cannot propagate as internal waves). This is handled in the code by checking ω_M2 > |f| 
            ! for M2 tidal forcing when computing iwe_M2_tau, if  ω_M2 < |f| --> iwe_M2_tau=0.0
            iwe2_omega_M2 =  2*pi/( 12*60*60 + 25.2 *60 )   ! M2 frequency in 1/s
            
            file_exist=.False.
            inquire(file=trim(idemix2_M2forc_file),exist=file_exist) 
            if (file_exist) then
                if (mype==0) write(*,*) '     ├> read IDEMIX2 M2 wave forcing'
                call read_other_NetCDF_3d(                             &
                                          trim(idemix2_M2forc_file)  , & 
                                          trim(idemix2_M2forc_vname) , & 
                                          trim(idemix2_M2forc_zname) , &
                                          iwe2_fM2                   , &  
                                          .true.                     , & !.true. -> on node
                                          partit                     , &
                                          mesh                         &
                                          )
                do fbin_i = 1, nfbin                          
                    iwe2_fM2(fbin_i,:) = iwe2_fM2(fbin_i,:)/density_0/iwe2_dphit(fbin_i)
                end do    
            else
                if (mype==0) then
                    print *, achar(27)//'[33m'
                    write(*,*) '____________________________________________________________________'
                    write(*,*) ' ERROR: IDEMIX2 M2 forcing file not found! Cant apply IDEMIX2'
                    write(*,*) '        vertical mixing parameterisation! '
                    write(*,*) '        ├> file: ', trim(idemix2_M2forc_file)
                    write(*,*) '        └> check: namelist.cvmix, idemix2_M2forc_file &  '
                    write(*,*) '____________________________________________________________________'
                    print *, achar(27)//'[0m'
                    write(*,*)
                    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
                end if
            end if 
            
            t1=MPI_Wtime()
            if (mype==0) write(*,*) '     │  └> elapsed time:', t1-t0
        end if ! --> if (idemix2_enable_M2) then
        
        !_______________________________________________________________________
        ! read idemix niw or srf forcing from cfsr data --> file 
        t0=MPI_Wtime()
        file_exist=.False.
        inquire(file=trim(idemix2_niwforc_file),exist=file_exist) 
        if (file_exist) then
            call read_other_NetCDF(                                    &
                                    trim(idemix2_niwforc_file)       , & 
                                    trim(idemix2_niwforc_vname)      , &
                                    1                                , & 
                                    iwe2_fsrf                        , & 
                                    .true.                           , & ! NN for missing values
                                    .true.                           , & ! interpolate to vertices  
                                    partit                           , & 
                                    mesh                               &
                                    )
            
            
            ! only 20% (idemix2_fniw_usage) from the surface forcing goes into internal waves
            ! iwe2_fsrf becomes surface forcing variable when idemix2_enable_niw=.false.
            ! in this case standard idemix1 behaviour 
            iwe2_fsrf = max(0.0_WP,iwe2_fsrf)
            iwe2_fsrf = iwe2_fsrf/density_0*idemix2_fniw_usage
            
            ! check for total tidal energy that is infused through the surface, see how 
            ! much is lossed during interpolation and compare with value of the 
            ! original files
            loc_Etot = 0.0_WP
            do node=1, myDim_nod2D
                loc_Etot = loc_Etot + area(1,node)*iwe2_fsrf(node)*density_0                    
            end do
            call MPI_AllREDUCE(loc_Etot, glb_Etot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
            
            if (idemix2_enable_niw) then
                if (mype==0) write(*,*) '     ├> read IDEMIX2 niw wave forcing'
                do fbin_i =2, nfbin-1
                    iwe2_fniw(fbin_i, :) = iwe2_fsrf(:)/(2.0_WP*pi)
                end do 
                iwe2_fsrf(:) = 0.0_WP
            else
                if (mype==0) write(*,*) '     ├> read IDEMIX2 surface wave forcing'
            end if
            if (mype==0) write(*,*) '     ├> IDEMIX2 total srf. energy Etot_srf =', glb_Etot*1.0e-12, ' TW'
            
        else
            if (mype==0) then
                print *, achar(27)//'[33m'
                write(*,*) '____________________________________________________________________'
                write(*,*) ' ERROR: IDEMIX2 NIW/Srf. forcing file not found! '
                write(*,*) '        Required for srf. energy input even when idemix2_enable_niw=.false.!'
                write(*,*) '        Cant apply IDEMIX2 vertical mixing parameterisation! '
                write(*,*) '        ├> file: ', trim(idemix2_niwforc_file)
                write(*,*) '        └> check: namelist.cvmix, idemix2_niwforc_file &  '
                write(*,*) '____________________________________________________________________'
                print *, achar(27)//'[0m'
                write(*,*)
                call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
            end if
        end if 
        
        
        ! NIWs are super-inertial (ω_niw > |f|), which is essential for their propagation as internal waves.
        ! The factor 1.05 ensures the frequency is 5% above the local inertial frequency.
        ! Why 1.05? Represents the minimum frequency for NIWs to propagate. Accounts for the Doppler shift 
        ! from background flows. Physical Context: NIWs are generated by wind events and have frequencies 
        ! close to f. They can only propagate as internal waves if ω > |f|.
        if (idemix2_enable_niw) then
            do node = 1, myDim_nod2D+eDim_nod2D
                iwe2_omega_niw(node) = max(1.0e-8_WP, abs( 1.05 * mesh%coriolis_node(node) ) )
            end do
        end if 
        t1=MPI_Wtime()
        if (mype==0) write(*,*) '     │  └> elapsed time:', t1-t0
        
        !_______________________________________________________________________
        ! read idemix bottom near tidal forcing from cesm data set --> file 
        ! from N. Brüggemann interpoalted to regular grid
        if (idemix2_enable_bot) then
            t0=MPI_Wtime()
            file_exist=.False.
            inquire(file=trim(idemix2_botforc_file),exist=file_exist) 
            if (file_exist) then
                if (mype==0) write(*,*) '     ├> read IDEMIX2 near tidal bottom forcing'
                call read_other_NetCDF(trim(idemix2_botforc_file), trim(idemix2_botforc_vname), 1, iwe2_fbot_n, .true., .true., partit, mesh)                !                                                                                                |       |
                !                                                   .true.=NN for missing values ----------------+       |
                !                                                   .true.=interpolate to vertices ----------------------+

                ! make sure forcing is all positive no numerical negative values
                iwe2_fbot_n = max(0.0_WP, iwe2_fbot_n)

                ! check for total tidal energy that is infused through the bottom
                loc_Etot = 0.0_WP
                do node=1, myDim_nod2D
                    loc_Etot = loc_Etot + area(1,node)*iwe2_fbot_n(node)
                end do
                call MPI_AllREDUCE(loc_Etot, glb_Etot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
                if (mype==0) write(*,*) "     │  └> IDEMIX2 total tidal energy Etot_bot =", glb_Etot*1.0e-12, ' TW'

                ! normalize total tidal energy at bottom with respect to the total
                ! tidal energy that is e.g in the original forcing files to accomodate
                ! non concerving losses during interpolation. This is only done when
                ! in namelist.cvmix: idemix2_botforc_Etot \= 0.0_WP
                if (idemix2_botforc_Etot /= 0.0_WP) then
                    iwe2_fbot_n = iwe2_fbot_n * idemix2_botforc_Etot/glb_Etot

                    loc_Etot = 0.0_WP
                    do node=1, myDim_nod2D
                        loc_Etot = loc_Etot + area(1,node)*iwe2_fbot_n(node)
                    end do
                    call MPI_AllREDUCE(loc_Etot, glb_Etot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
                    if (mype==0) write(*,*) "     │  └> IDEMIX2 Etot_bot after normalizing =", glb_Etot*1.0e-12, ' TW'
                end if

                ! divide by density_0 --> convert from W/m^2 to m^3/s^3
                iwe2_fbot_n = iwe2_fbot_n/density_0
                
            else
                if (mype==0) then
                    print *, achar(27)//'[33m'
                    write(*,*) '____________________________________________________________________'
                    write(*,*) ' ERROR: IDEMIX2 bottom forcing file not found! Cant apply IDEMIX'
                    write(*,*) '        vertical mixing parameterisation! '
                    write(*,*) '        ├> file: ', trim(idemix2_botforc_file)
                    write(*,*) '        └> check: namelist.cvmix, idemix2_botforc_file'
                    write(*,*) '____________________________________________________________________'
                    print *, achar(27)//'[0m'
                    write(*,*)
                    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
                end if 
            end if 
            t1=MPI_Wtime()
            if (mype==0) write(*,*) '     │  └> elapsed time:', t1-t0
        end if !--> if (idemix2_enable_bot) then
        
!         !_______________________________________________________________________
!         ! read Lee-Wave forcing 
!         if (idemix2_enable_leew) then 
!             t0=MPI_Wtime()
!             file_exist=.False.
!             inquire(file=trim(idemix2_leewforc_file),exist=file_exist) 
!             if (file_exist) then
!                 if (mype==0) write(*,*) '     ├> read IDEMIX2 lee wave forcing --> add to bottom forcing'
!                 call read_other_NetCDF(trim(idemix2_leewforc_file), trim(idemix2_leewforc_vname), 1, iwe2_fleew, .true., .false., partit, mesh)
!                 
!                 iwe2_fbot_n = iwe2_fbot_n + iwe2_fleew
!                 !                      |
!                 !                      +-> no division by density_0, Lee wave Forcing
!                 !                          already in units of m^3/s^3
!             else
!                 if (mype==0) then
!                     print *, achar(27)//'[33m'
!                     write(*,*) '____________________________________________________________________'
!                     write(*,*) ' ERROR: IDEMIX2 Lee-Wave forcing file not found! Cant apply IDEMIX2'
!                     write(*,*) '        vertical mixing parameterisation! '
!                     write(*,*) '        ├> file: ', trim(idemix2_hrmsforc_file)
!                     write(*,*) '        └> check: namelist.cvmix, idemix2_leewforc_file'
!                     write(*,*) '____________________________________________________________________'
!                     print *, achar(27)//'[0m'
!                     write(*,*)
!                     call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
!                 end if
!             end if 
!         end if ! --> if (idemix2_enable_leew) then 
        
        !_______________________________________________________________________
        ! read idemix HRMS and HLAM forcing --> file 
        if (idemix2_enable_M2 .or. idemix2_enable_niw ) then
            t0=MPI_Wtime()
            ! topo_hrms (Root Mean Square Topographic Height), Definition: Measures the roughness or variability of seafloor topography.
            ! Units: Meters (m), Role: Represents the standard deviation of seafloor height variations. Used to calculate 
            ! energy transfer from tides/NIW to internal waves.
            file_exist=.False.
            inquire(file=trim(idemix2_hrmsforc_file),exist=file_exist) 
            if (file_exist) then
                if (mype==0) write(*,*) '     ├> read IDEMIX2 HRMS forcing'
                call read_other_NetCDF(trim(idemix2_hrmsforc_file), trim(idemix2_hrmsforc_vname), 1, iwe2_topo_hrms, .true., .true., partit, mesh)
            else
                if (mype==0) then
                    print *, achar(27)//'[33m'
                    write(*,*) '____________________________________________________________________'
                    write(*,*) ' ERROR: IDEMIX2 HRMS forcing file not found! Cant apply IDEMIX2'
                    write(*,*) '        vertical mixing parameterisation! '
                    write(*,*) '        ├> file: ', trim(idemix2_hrmsforc_file)
                    write(*,*) '        └> check: namelist.cvmix, idemix2_hrmsforc_file'
                    write(*,*) '____________________________________________________________________'
                    print *, achar(27)//'[0m'
                    write(*,*)
                    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
                end if
            end if 
            
            ! topo_hlam (Topographic Wavelength), Definition: Characteristic horizontal length scale of topographic features.
            ! Units: Meters (m), Role: Represents the dominant wavelength of seafloor roughness. Used to normalize the topographic forcing.
            file_exist=.False.
            inquire(file=trim(idemix2_hlamforc_file),exist=file_exist) 
            if (file_exist) then
                if (mype==0) write(*,*) '     ├> read IDEMIX2 HLAM forcing'
                call read_other_NetCDF(trim(idemix2_hlamforc_file), trim(idemix2_hlamforc_vname), 1, iwe2_topo_hlam, .true., .true., partit, mesh)
            else
                if (mype==0) then
                    print *, achar(27)//'[33m'
                    write(*,*) '____________________________________________________________________'
                    write(*,*) ' ERROR: IDEMIX2 HLAM forcing file not found! Cant apply IDEMIX2'
                    write(*,*) '        vertical mixing parameterisation! '
                    write(*,*) '        ├> file: ', trim(idemix2_hlamforc_file)
                    write(*,*) '        └> check: namelist.cvmix, idemix2_hlamforc_file'
                    write(*,*) '____________________________________________________________________'
                    print *, achar(27)//'[0m'
                    write(*,*)
                    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
                end if
            end if
            
            ! In the M2 tidal and NIW energy calculations:
            !  --> fxc = topo_hrms(i,j)**2 * 2*pi / (1d-12 + topo_lam(i,j))
            ! This term estimates the energy flux from barotropic tides to internal waves. The ratio topo_hrms²/topo_lam represents 
            ! the topographic "steepness" affecting wave generation.
            ! Physical Interpretation: 
            !  --> High topo_hrms and small topo_lam: Indicates rough, steep topography that efficiently generates internal waves. 
            !  --> Low topo_hrms and large topo_lam: Indicates smooth, flat topography with minimal wave generation.
            ! These parameters are crucial for accurately representing the conversion of barotropic tidal energy to internal 
            ! waves over rough topography.
            
            t1=MPI_Wtime()
            if (mype==0) write(*,*) '     │  └> elapsed time:', t1-t0
        
        end if ! --> if (idemix2_enable_M2 .or. idemix2_enable_niw ) then
        
        !_______________________________________________________________________
        ! distribute nodal bottom flux across staircase depth levels using cap
        ! area fractions: iwe2_fbot(nz,n) = fbot_n * (area(nz-1,n)-area(nz,n)) / area(1,n)
        do node = 1, myDim_nod2D+eDim_nod2D
            nlu = nlevels_nod2D_min(node)
            nln = nlevels_nod2D(node)
            do nz = nlu, nln
                iwe2_fbot(nz, node) = iwe2_fbot_n(node) * (area(nz-1, node) - area(nz, node)) / area(1, node)
            end do
        end do
        
        !_______________________________________________________________________
        ! compute centroid distance from nearset coastal point together with
        ! idemix2_shelf_dist defines what is shelf and what not  
        if (idemix2_enable_M2 .or. idemix2_enable_niw) then 
            t0=MPI_Wtime() 
            call compute_dist2coast(iwe2_topo_dist, mesh, partit, 'node')
            t1=MPI_Wtime()
            if (partit%mype==0) write(*,*) '        └> dist2coast elapsed time:', t1-t0
        end if 
        !_______________________________________________________________________
        ! compute d/dy of coriolis 
        do elem = 1, myDim_elem2D
            elnodes = elem2d_nodes(:,elem)
            iwe2_grady_coriol(elem) = sum(gradient_sca(4:6, elem)*mesh%coriolis_node(elnodes))
        end do       
        
        !_______________________________________________________________________
        ! compute scalar cell volume from top to bottom 
        call compute_vol_nodB2T_fix(vol_nodB2T, mesh, partit)
        
        !_______________________________________________________________________
        ! initialise IDEMIX parameters
        call cvmix_idemix2_init(                                 &
              tau_v               = idemix2_tau_v                & !IN: vertical dissipation timescale (s)
            , tau_h               = idemix2_tau_h                & !IN: horizontal diffusion timescale (s)
            , gamma               = idemix2_gamma                & !IN: spectral bandwidth parameter
            , jstar               = idemix2_jstar                & !IN: GM76 spectral peak mode number
            , mu0                 = idemix2_mu0                  & !IN: wave-wave interaction coefficient
            , nfbin               = idemix2_nfbin                & !IN: number of spectral angular bins
            , shelf_dist          = idemix2_shelf_dist           & !IN: shelf distance threshold (m)
            , enable_M2           = idemix2_enable_M2            & !IN: activate M2 tidal compartment
            , enable_niw          = idemix2_enable_niw           & !IN: activate NIW compartment
            , enable_AB_timestep  = idemix2_enable_AB            & !IN: use Adams-Bashforth 2nd-order time stepping
            , enable_hor_diffusion= idemix2_enable_hor_diff_expl & !IN: enable horizontal Laplacian diffusion
            , enable_hor_diff_iter= idemix2_enable_hor_diff_impl_iter & !IN: use iterative horizontal diffusion
            , hor_diff_niter      = idemix2_hor_diff_niter       & !IN: number of horizontal diffusion iterations
            )
                                
        if (partit%mype==0) write(*,*)                  
    end subroutine init_cvmix_idemix2



    !
    !
    !
    !===========================================================================
    ! pre-compute reflective coastal BC lookup tables iwe2_bc_reflect_bin and
    ! iwe2_bc_reflect_sgn.  Must be called after iwe2_phit is initialised.
    ! 
    ! SETUP (pre-computed once in init_reflect_bc)
    ! ─────────────────────────────────────────────
    ! LAND │  n_coast ──────edge──────  n_int
    !         │    sgn=+1                (interior)
    !         │
    !         │  iwe2_bc_reflect_bin(edge) = kk   (mirror spectral bin)
    !         │  iwe2_bc_reflect_sgn(edge) = +1   (coast = ednodes(1))
    subroutine init_reflect_bc(   refl_bin   &
                                , refl_sgn   &
                                , refl_src   &
                                , refl_node  &
                                , partit     &
                                , mesh       &
                               )
        implicit none
        type(t_mesh)  , intent(in), target :: mesh
        type(t_partit), intent(in), target :: partit
        integer       , intent(inout)      :: refl_bin(partit%myDim_edge2D)
        integer       , intent(inout)      :: refl_sgn(partit%myDim_edge2D)
        real(kind=WP) , intent(inout)      :: refl_src(idemix2_nfbin, partit%myDim_nod2D)
        logical       , intent(inout)      :: refl_node(partit%myDim_nod2D+partit%eDim_nod2D)
        
        integer       :: edge, n1, n2, n_int, node, kk_bc, fbinj, nfbin
        real(kind=WP) :: lon1, lat1, lon2, lat2, lat_mean
        real(kind=WP) :: edge_dx, edge_dy, phi_mirror, dist, dist_min
        logical, allocatable :: is_coast_node(:)
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        nfbin = idemix2_nfbin
    
        ! > refl_bin(edge) stores the MIRROR SPECTRAL BIN for a flagged
        ! coast edge.  0 = edge not flagged, BC skipped.
        ! 
        ! > HOW IT IS COMPUTED (init_reflect_bc)
        ! For each edge with exactly one coast node:
        !   1. Compute direction FROM coast node BACK INTO the ocean:
        !        edge_dx = (lon_interior - lon_coast) * cos(lat_mean)
        !        edge_dy =  lat_interior - lat_coast
        !   2. phi_mirror = atan2(edge_dy, edge_dx)  mapped to [0, 2*pi)
        !   3. kk = interior bin (2..nfbin-1) with angle closest to phi_mirror
        !   4. refl_bin(edge) = kk
        ! 
        ! > HOW IT IS USED (apply_reflect_bc_spctrl):
        ! ALL toward-coast bins at this edge (any fbini where sgn*flux > 0)
        ! are blocked and their energy reflected into the SAME single bin kk:
        !   reflect_src(kk, n_int) += sgn * outflux   (for all toward-coast fbini)
        ! One kk per edge — an approximation based on the edge-normal geometry.
        ! refl_bin(:) = 0
        
        ! > refl_sgn(edge) encodes WHICH of the two edge nodes is the coast
        ! node, and from that, WHICH flux sign means "toward coast". Each edge has 
        ! two endpoint nodes: edges(1, edge) and edges(2, edge). The sign convention 
        ! in the FV scheme is:
        !   flux(fbini, edge) > 0  =  transport from edges(2) -> edges(1)
        !   flux(fbini, edge) < 0  =  transport from edges(1) -> edges(2)
        !   sgn = +1  -->  coast node = edges(1, edge)
        !                  interior   = edges(2, edge)
        !                  flux toward coast when: flux > 0
        !   sgn = -1  -->  coast node = edges(2, edge)
        !                  interior   = edges(1, edge)
        !                  flux toward coast when: flux < 0
        ! 
        ! > The single condition:
        !     sgn * flux(fbini, edge) > 0
        ! 
        ! > captures both cases -- it is true whenever the flux is directed toward
        ! whichever node is the coast node, regardless of edge orientation.
        ! It also determines the interior node n_int:
        !     n_int = merge(edges(2, edge), edges(1, edge), sgn == +1)
        ! And the reflected energy magnitude is always:
        !     sgn * outflux = |outflux| > 0
        ! 
        ! in both orientations (positive for sgn=+1 with outflux>0,
        !                       positive for sgn=-1 with outflux<0).
        ! refl_sgn(:) = 0
        
        ! > refl_src is a 2D scratch array (nfbin, nod2D) used to implement
        ! reflective boundary conditions at coastlines for horizontal NIW energy advection.
        ! 1. During apply_reflect_bc_spctrl:
        !    For every coast-adjacent edge flagged in iwe2_bc_reflect_bin, any flux
        !    directed toward the coast node is zeroed out and the blocked energy is
        !    accumulated into refl_src(kk, n_int) at the interior node (n_int),
        !    with the sign corrected so energy is conserved.
        ! 2. After advection in hsintegrate_Ecompart:
        !    The accumulated blocked flux is re-injected as a divergence source
        !    Edivh += refl_src / vol_s, effectively bouncing the energy back
        !    into the interior in bin kk (the mirror/reflected propagation direction).
        !    
        ! > In short: it is a flux-blocking + re-injection buffer that prevents NIW energy
        ! from leaking into land, redirecting it instead into the spectrally mirrored
        ! propagation direction at the last ocean node.        
        ! refl_src(:,:) = 0.0_WP
        
        ! a node is a coast node if it belongs to at least one boundary edge
        allocate(is_coast_node(myDim_nod2D+eDim_nod2D))
        is_coast_node(:)       = .false.
        do edge = 1, myDim_edge2D
            if (edge_tri(2, edge) == 0) then
                is_coast_node(edges(1, edge)) = .true.
                is_coast_node(edges(2, edge)) = .true.
            end if
        end do

        do edge = 1, myDim_edge2D
            if (edge_tri(2, edge) == 0) cycle   ! skip coast-parallel boundary edges
            n1 = edges(1, edge)
            n2 = edges(2, edge)
            ! only process edges with exactly one coast node
            if (is_coast_node(n1) .eqv. is_coast_node(n2)) cycle
            
            ! geographic coordinates of the two nodes (lon/lat in radians)
            lon1 = geo_coord_nod2D(1, n1);  lat1 = geo_coord_nod2D(2, n1)
            lon2 = geo_coord_nod2D(1, n2);  lat2 = geo_coord_nod2D(2, n2)
            lat_mean = 0.5_WP*(lat1 + lat2)
            
            ! edge direction FROM coast node BACK TO interior node (the reflection direction)
            if (is_coast_node(n2)) then
                ! coast = n2, interior = n1: reflect direction is n2 -> n1
                edge_dx = (lon1 - lon2) * cos(lat_mean)
                edge_dy =  lat1 - lat2
                refl_sgn(edge) = -1  ! negative flux = transport n1->n2 = INTO coast
            else
                ! coast = n1, interior = n2: reflect direction is n1 -> n2
                edge_dx = (lon2 - lon1) * cos(lat_mean)
                edge_dy =  lat2 - lat1
                refl_sgn(edge) = +1  ! positive flux = transport n2->n1 = INTO coast
            end if
            
            ! mirror angle = direction away from coast, mapped to [0, 2*pi)
            phi_mirror = atan2(edge_dy, edge_dx)
            if (phi_mirror < 0.0_WP) phi_mirror = phi_mirror + 2.0_WP*pi
            
            ! find nearest interior spectral bin (skip ghost cells 1 and nfbin)
            kk_bc    = 2
            dist_min = (iwe2_phit(2) - phi_mirror)**2
            do fbinj = 3, nfbin-1
                dist = (iwe2_phit(fbinj) - phi_mirror)**2
                if (dist < dist_min) then
                    dist_min = dist
                    kk_bc    = fbinj
                end if
            end do
            refl_bin(edge) = kk_bc
        end do

        ! create logical array if node is a reflected boudnary node
        refl_node(:) = .false.
        do edge = 1, myDim_edge2D
            if (refl_bin(edge) == 0) cycle
            n_int = merge(edges(2,edge), edges(1,edge), refl_sgn(edge) == +1)
            if (n_int <= myDim_nod2D) refl_node(n_int) = .true.
        end do

        deallocate(is_coast_node)
          
        if (mype==0) write(*,*) '     ├> IDEMIX2 reflective coastal BC lookup table built'
    end subroutine init_reflect_bc



    !
    !
    !
    !===========================================================================
    ! calculate IDEMIX2 internal wave energy and its dissipation
    subroutine calc_cvmix_idemix2(istep, partit, mesh)
        implicit none
        integer,        intent(in)            :: istep
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(inout), target :: partit
        integer       :: node, elem, edge, node_size, elem_size, k, fbini, nfbin
        integer       :: nz, nln, nl1, nl2, nl12, nu1, nu2, nu12, uln, iter
        integer       :: elnodes(3), el(2), ednodes(2)
        real(kind=WP) :: lat_n_deg, lat_e_deg
        real(kind=WP) :: cn, cn_e, cn_gradx_e, cn_grady_e, omega_niw_e
        real(kind=WP) :: inv_area, area_third
        logical       :: topo_shelf=.False.
        ! Timing variables
        real(kind=WP) :: t_start, t_end, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
        real(kind=WP) :: time_params, time_groupvel, time_M2_integ, time_niw_integ
        real(kind=WP) :: time_Eiw_vdiff=0.0, time_Eiw_hdiff=0.0, time_waveint=0.0, time_Kv_Av=0.0, time_total=0.0

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        node_size  = myDim_nod2D
        elem_size  = myDim_elem2D
        nfbin      = idemix2_nfbin
        
        !_______________________________________________________________________
        ! Start total timing
        t_start = MPI_Wtime()
        
        !_______________________________________________________________________
        ! do tau timestepping indices shift (2-slot: ti and tip1 alternate 1<->2)
        iwe2_ti   = iwe2_tip1
        iwe2_tip1 = 3 - iwe2_ti
        
        !_______________________________________________________________________
        ! compute scalar cell volume from top to bottom 
        ! --> if we take not the ALE variying ocean layer into account --> slight 
        !     numerical drift
        ! --> if take the fixed standard levels (zbar) to compute the ocean volume 
        !     than only normal numerical drift, ocean volume can be computed at 
        !     initialisation point
        ! call compute_vol_nodB2T_hnode(vol_nodB2T, mesh, partit)
        
        ! compute inverse cell volume for horizontal diffusion; uses Z_3d_n and
        ! hnode which change every timestep under ALE, so must be updated here
        call compute_vol_wcell(vol_wcelli, mesh, partit)
        
        !_______________________________________________________________________
        ! 1. Compute IDEMIX2 parameters (cn, alpha_c, v0, c0, tau, struct functions)
        t0 = MPI_Wtime()
!$OMP PARALLEL DO PRIVATE(node, lat_n_deg, nln, uln, nz, cn, topo_shelf)
        do node=1, myDim_nod2D+eDim_nod2D
            !___________________________________________________________________
            uln = ulevels_nod2D(node)
            nln = nlevels_nod2D(node)
            
            !___________________________________________________________________
            ! re-initialse cross spectral velocites, for later accumulation through 
            ! elemental values
            if (idemix2_enable_M2 ) iwe2_M2_w( :,node) = 0.0_WP
            if (idemix2_enable_niw) iwe2_niw_w(:,node) = 0.0_WP
            
            !___________________________________________________________________
            lat_n_deg = geo_coord_nod2D(2,node) / rad
            
            !___________________________________________________________________
            ! compute baroclinic gravity wave speed            
            cn=0.0_WP
            do nz=uln, nln-1
                cn=cn+hnode(nz,node)*(sqrt(max(bvfreq(nz,node), 0._WP)) + sqrt(max(bvfreq(nz+1,node), 0._WP)))/2._WP
            end do
            cn = cn/pi * idemix2_scal_cn
            cn = max(cn, idemix2_cn_min)
            iwe2_cn(node)=cn
            
            !
            !
            !___________________________________________________________________
            ! 1st. compute idemix2 parameter over the vertical water column and 
            ! local horizontal parameters on elements !!!:
            call cvmix_idemix2_compute_param(            &
                  nlev     = nln-uln+1                   & !IN:  number of active levels
                , coriolis = mesh%coriolis_node(   node) & !IN:  Coriolis parameter (1/s)
                , Nsqr     = bvfreq(      uln:nln, node) & !IN:  buoyancy frequency squared (1/s²)
                , cn       = cn                          & !IN:  first baroclinic mode speed (m/s)
                , alpha_c  = iwe2_alpha_c(uln:nln, node) & !OUT: wave-wave dissipation rate (1/s)
                , c0       = iwe2_c0(     uln:nln, node) & !OUT: vertical representative group speed (m/s)
                , v0       = iwe2_v0(     uln:nln, node) & !OUT: horizontal representative group speed (m/s)
                )
            
            !
            !
            !___________________________________________________________________
            ! 3rd. compute dissipation time scales and rates
            if (idemix2_enable_M2 .or. idemix2_enable_niw) then 
                if (iwe2_topo_dist(node) <= idemix2_shelf_dist) then
                    topo_shelf = .True.
                else    
                    topo_shelf = .False.
                end if
                if (idemix2_enable_M2) then 
                    call cvmix_idemix2_compute_compart_interact_tscale(&
                          dtime        = dt                          & !IN:  time step (s)
                        , lat          = lat_n_deg                   & !IN:  latitude (degrees)
                        , coriolis     = mesh%coriolis_node(   node) & !IN:  Coriolis parameter (1/s)
                        , cn           = cn                          & !IN:  first baroclinic mode speed (m/s)
                        , zbottom      = -zbar_n_bot(          node) & !IN:  bottom depth, positive (m)
                        , topo_hrms    = iwe2_topo_hrms(       node) & !IN:  RMS topographic height (m)
                        , topo_hlam    = iwe2_topo_hlam(       node) & !IN:  characteristic topographic length scale (m)
                        , topo_shelf   = topo_shelf                  & !IN:  .true. if node is on shelf
                        , omega_compart= iwe2_omega_M2               & !IN:  M2 tidal frequency (1/s)
                        , tau_compart  = iwe2_M2_tau(          node) & !OUT: M2 bottom scattering timescale (1/s)
                        , compart_name = 'M2'                        & !IN:  compartment name
                        )
                    call cvmix_idemix2_compute_M2_dissipation(       &
                          lat          = lat_n_deg                   & !IN:  latitude (degrees)
                        , cn           = cn                          & !IN:  first baroclinic mode speed (m/s)
                        , zbottom      = -zbar_n_bot(          node) & !IN:  bottom depth, positive (m)
                        , alpha_M2_c   = iwe2_alpha_M2_c(      node) & !OUT: M2 continuous (PSI) dissipation rate (1/s)
                        )
                end if 
                if (idemix2_enable_niw) then 
                    call cvmix_idemix2_compute_compart_interact_tscale(&
                          dtime        = dt                          & !IN:  time step (s)
                        , lat          = lat_n_deg                   & !IN:  latitude (degrees)
                        , coriolis     = mesh%coriolis_node(   node) & !IN:  Coriolis parameter (1/s)
                        , cn           = cn                          & !IN:  first baroclinic mode speed (m/s)
                        , zbottom      = -zbar_n_bot(          node) & !IN:  bottom depth, positive (m)
                        , topo_hrms    = iwe2_topo_hrms(       node) & !IN:  RMS topographic height (m)
                        , topo_hlam    = iwe2_topo_hlam(       node) & !IN:  characteristic topographic length scale (m)
                        , topo_shelf   = topo_shelf                  & !IN:  .true. if node is on shelf
                        , omega_compart= iwe2_omega_niw(       node) & !IN:  NIW frequency (lat-dependent, 1/s)
                        , tau_compart  = iwe2_niw_tau(         node) & !OUT: NIW bottom scattering timescale (1/s)
                        , compart_name = 'niw'                       & !IN:  compartment name
                        )
                end if  
            end if 
            
            !
            !
            !___________________________________________________________________
            ! 4th. compute structure function for M2 and NIW energy
            if ((idemix2_enable_M2) .and. (idemix2_enable_niw)) then 
                call cvmix_idemix2_compute_vert_struct_fct(          &
                      nlev        = nln-uln+1                        & !IN:    number of active levels
                    , dzw         = hnode(          uln:nln-1, node) & !IN:    layer thickness (m)
                    , coriolis    = mesh%coriolis_node(        node) & !IN:    Coriolis parameter (1/s)
                    , cn          = cn                               & !IN:    first baroclinic mode speed (m/s)
                    , Nsqr        = bvfreq(           uln:nln, node) & !IN:    buoyancy frequency squared (1/s²)
                    , omega_M2    = iwe2_omega_M2                    & !IN:    M2 tidal frequency (1/s)
                    , omega_niw   = iwe2_omega_niw(            node) & !IN:    NIW frequency, lat-dependent (1/s)
                    , E_struct_M2 = iwe2_E_M2_struct( uln:nln, node) & !INOUT: M2 energy vertical structure function
                    , E_struct_niw= iwe2_E_niw_struct(uln:nln, node) & !INOUT: NIW energy vertical structure function
                    )
                    
            else if (idemix2_enable_M2) then 
                call cvmix_idemix2_compute_vert_struct_fct(          &
                      nlev        = nln-uln+1                        & !IN:    number of active levels
                    , dzw         = hnode(          uln:nln-1, node) & !IN:    layer thickness (m)
                    , coriolis    = mesh%coriolis_node(        node) & !IN:    Coriolis parameter (1/s)
                    , cn          = cn                               & !IN:    first baroclinic mode speed (m/s)
                    , Nsqr        = bvfreq(           uln:nln, node) & !IN:    buoyancy frequency squared (1/s²)
                    , omega_M2    = iwe2_omega_M2                    & !IN:    M2 tidal frequency (1/s)
                    , E_struct_M2 = iwe2_E_M2_struct( uln:nln, node) & !INOUT: M2 energy vertical structure function
                    )
                    
            else if (idemix2_enable_niw) then 
                call cvmix_idemix2_compute_vert_struct_fct(          &
                      nlev        = nln-uln+1                        & !IN:    number of active levels
                    , dzw         = hnode(          uln:nln-1, node) & !IN:    layer thickness (m)
                    , coriolis    = mesh%coriolis_node(        node) & !IN:    Coriolis parameter (1/s)
                    , cn          = cn                               & !IN:    first baroclinic mode speed (m/s)
                    , Nsqr        = bvfreq(           uln:nln, node) & !IN:    buoyancy frequency squared (1/s²)
                    , omega_niw   = iwe2_omega_niw(            node) & !IN:    NIW frequency, lat-dependent (1/s)
                    , E_struct_niw= iwe2_E_niw_struct(uln:nln, node) & !INOUT: NIW energy vertical structure function
                    )
                    
            end if         
        end do ! --> do node=1, myDim_nod2D+eDim_nod2D
!$OMP END PARALLEL DO
        t1 = MPI_Wtime()
        time_params = t1 - t0

        !___________________________________________________________________________
        if (idemix2_enable_M2 .or. idemix2_enable_niw) then 
            !_______________________________________________________________________
            ! 2. Compute idemix2 group velocites for M2 an NIW
            t1 = MPI_Wtime()
            do elem = 1, myDim_elem2D
                nln        = nlevels(elem)-1
                uln        = ulevels(elem)
                elnodes    = elem2d_nodes(:,elem)
                lat_e_deg  = sum(geo_coord_nod2D(2,elnodes)) / (3.0_WP*rad)
                area_third = elem_area(elem) / 3.0_WP
                if (idemix2_enable_M2)  w_M2_e(:)  = 0.0_WP
                if (idemix2_enable_niw) w_niw_e(:) = 0.0_WP
                
                !___________________________________________________________________
                ! compute baroclionic gravity wave speed on elements
                cn_e        = sum(iwe2_cn(elnodes))/3.0
                
                ! compute gradient of  baroclionic gkdot_x_M2 = sqrt(fxa)/omega_M2*cn_gradxravity wave speed on elements
                cn_gradx_e  = sum(gradient_sca(1:3,elem)*iwe2_cn(elnodes))
                cn_grady_e  = sum(gradient_sca(4:6,elem)*iwe2_cn(elnodes))
                                
                ! average iwe2_omega_niw to elem
                if (idemix2_enable_niw) omega_niw_e = sum(iwe2_omega_niw(elnodes))/3.0
                
                !___________________________________________________________________
                if (idemix2_enable_M2) then 
                    call cvmix_idemix2_compute_compart_groupvel(     &
                          nfbin        = idemix2_nfbin               & !IN:  number of spectral angular bins
                        , coriolis     = mesh%coriolis(        elem) & !IN:  Coriolis parameter at element (1/s)
                        , coriol_grady = iwe2_grady_coriol(    elem) & !IN:  meridional gradient of Coriolis (1/(m·s))
                        , coslat       = mesh%elem_cos(        elem) & !IN:  cosine of latitude at element
                        , cn           = cn_e                        & !IN:  first baroclinic mode speed at element (m/s)
                        , cn_gradx     = cn_gradx_e                  & !IN:  zonal gradient of cn at element (1/s)
                        , cn_grady     = cn_grady_e                  & !IN:  meridional gradient of cn at element (1/s)
                        , phit         = iwe2_phit                   & !IN:  angular bin edges (rad)
                        , phiu         = iwe2_phiu                   & !IN:  angular bin centres (rad)
                        , omega_compart= iwe2_omega_M2               & !IN:  M2 tidal frequency (1/s)
                        , u_compart    = iwe2_M2_uv(1,      :, elem) & !OUT: zonal group velocity per bin (m/s)
                        , v_compart    = iwe2_M2_uv(2,      :, elem) & !OUT: meridional group velocity per bin (m/s)
                        , w_compart    = w_M2_e(            :      ) & !OUT: cross-spectral propagation rate per bin (1/s)
                        )
                    ! --> here w_M2_e and w_niw_e are still on elements but for the proper 
                    !     finite volume advection implementation we need them on nodes !
                    do k=1,3
                        iwe2_M2_w(:, elnodes(k)) = iwe2_M2_w(:, elnodes(k)) + w_M2_e(:)*area_third
                    end do
                    
                end if
                if (idemix2_enable_niw) then 
                    call cvmix_idemix2_compute_compart_groupvel(      &
                          nfbin        = idemix2_nfbin                & !IN:  number of spectral angular bins
                        , coriolis     = mesh%coriolis(        elem)  & !IN:  Coriolis parameter at element (1/s)
                        , coriol_grady = iwe2_grady_coriol(    elem)  & !IN:  meridional gradient of Coriolis (1/(m·s))
                        , coslat       = mesh%elem_cos(        elem)  & !IN:  cosine of latitude at element
                        , cn           = cn_e                         & !IN:  first baroclinic mode speed at element (m/s)
                        , cn_gradx     = cn_gradx_e                   & !IN:  zonal gradient of cn at element (1/s)
                        , cn_grady     = cn_grady_e                   & !IN:  meridional gradient of cn at element (1/s)
                        , phit         = iwe2_phit                    & !IN:  angular bin edges (rad)
                        , phiu         = iwe2_phiu                    & !IN:  angular bin centres (rad)
                        , omega_compart= omega_niw_e                  & !IN:  NIW frequency, elem-averaged (1/s)
                        , u_compart    = iwe2_niw_uv(1,      :, elem) & !OUT: zonal group velocity per bin (m/s)
                        , v_compart    = iwe2_niw_uv(2,      :, elem) & !OUT: meridional group velocity per bin (m/s)
                        , w_compart    = w_niw_e(            :      ) & !OUT: cross-spectral propagation rate per bin (1/s)
                        )
                    
                    ! Mask out NIW velocities near equator where physics is singular
                    if (abs(lat_e_deg) < 5.0_WP) then  ! ~5° latitude
                        iwe2_niw_uv( :, :, elem) = 0.0_WP
                        w_niw_e(           :   ) = 0.0_WP
                    end if
                    
                    do k=1,3
                        iwe2_niw_w(:, elnodes(k)) = iwe2_niw_w(:, elnodes(k)) + w_niw_e(:)*area_third
                    end do
                end if
                ! --> at the end we still need to normalize iwe2_M2_w and iwe2_niw_w 
                !     with the scalararea!
            end do ! --> do elem = 1, myDim_elem2D
            t2 = MPI_Wtime()
            time_groupvel = t2 - t1
            
            ! finalize elem2node averaging of iwe2_M2_w and iwe2_niw_w
            ! cross spectral exachange has to be related to nodes, since general advection 
            ! is related to nodes independent of the volume
            if (idemix2_enable_M2) then 
                call exchange_elem_fbin(iwe2_M2_uv, partit)
                call exchange_nod_fbin(iwe2_M2_w, partit)
    !$OMP PARALLEL DO PRIVATE(node, inv_area)
                do node=1, myDim_nod2D+eDim_nod2D
                    inv_area = 1.0_WP / area(1, node)
                    iwe2_M2_w(:, node) = iwe2_M2_w(:, node) * inv_area
                end do
    !$OMP END PARALLEL DO
            end if

            if (idemix2_enable_niw) then
                call exchange_elem_fbin(iwe2_niw_uv, partit)
                call exchange_nod_fbin(iwe2_niw_w, partit)
    !$OMP PARALLEL DO PRIVATE(node, inv_area)
                do node=1, myDim_nod2D+eDim_nod2D
                    inv_area = 1.0_WP / area(1, node)
                    iwe2_niw_w(:, node) = iwe2_niw_w(:, node) * inv_area
                end do
    !$OMP END PARALLEL DO
            end if
            
            
            !_______________________________________________________________________
            ! 3. Horizontal spectral integrate energy compartment E_M2^(n+1) and 
            ! E_niw^(n+1) equation. 
            ! dE/dt = -div(vec_u * E) - tau*E + Forc    
            ! E^(n+1) = E^n + dt*( -div(vec_u^n*E^n) - tau*E^n + Forc^n)
            t3 = MPI_Wtime()
            if (idemix2_enable_M2) then
                if (idemix2_diag_Ecompart) then
                    call hsintegrate_Ecompart(              &
                        iwe2_ti, iwe2_tip1                  & ! IN:    current and next time indices
                        , 'M2'                              & ! IN:    compartment name
                        , iwe2_E_M2                         & ! INOUT: energy array [2,nfbin,node]
                        , iwe2_E_M2_divh                    & ! INOUT: horizontal divergence flux accumulator
                        , iwe2_E_M2_divs                    & ! INOUT: cross-spectral divergence accumulator
                        , iwe2_M2_uv                        & ! IN:    horizontal group velocities on elements
                        , iwe2_M2_w                         & ! INOUT: cross-spectral propagation speed on nodes
                        , iwe2_fM2                          & ! IN:    M2 tidal surface forcing
                        , iwe2_M2_tau                       & ! IN:    M2 dissipation timescale per node
                        , iwe2_dphit                        & ! IN:    angular bin widths
                        , iwe2_gradxy_e                     & ! INOUT: energy gradient on elements (work array)
                        , iwe2_gradxy_n                     & ! INOUT: energy gradient on nodes (work array)
                        , iwe2_flx_uv                       & ! INOUT: horizontal edge flux (work array)
                        , iwe2_flx_w                        & ! INOUT: cross-spectral flux (work array)
                        , iwe2_refl_bin                     & ! IN:    mirror bin lookup table
                        , iwe2_refl_sgn                     & ! IN:    coast edge sign
                        , iwe2_refl_src                     & ! INOUT: reflect source accumulator
                        , iwe2_refl_node                    & ! IN:    coast-adjacent node mask
                        , vol_nodB2T                        & ! IN:    scalar cell volumes bottom-to-top
                        , partit                            & ! INOUT
                        , mesh                              & ! IN
                        , idemix2_enable_AB                 & ! IN:    enable AB2 time stepping
                        , iwe2_E_M2_dt                      & ! INOUT: optional diagnostic: total tendency
                        , iwe2_E_M2_advh                    & ! INOUT: optional diagnostic: horiz. advection
                        , iwe2_E_M2_advs                    & ! INOUT: optional diagnostic: spectral advection
                        , iwe2_E_M2_diss                    & ! INOUT: optional diagnostic: dissipation
                        , iwe2_E_M2_forc                    & ! INOUT: optional diagnostic: forcing
                        , iwe2_E_M2_refl                    & ! INOUT: optional diagnostic: coast reflection
                        )
                else
                    call hsintegrate_Ecompart(              &
                        iwe2_ti, iwe2_tip1                 & ! IN:    current and next time indices
                        , 'M2'                              & ! IN:    compartment name
                        , iwe2_E_M2                         & ! INOUT: energy array [2,nfbin,node]
                        , iwe2_E_M2_divh                    & ! INOUT: horizontal divergence flux accumulator
                        , iwe2_E_M2_divs                    & ! INOUT: cross-spectral divergence accumulator
                        , iwe2_M2_uv                        & ! IN:    horizontal group velocities on elements
                        , iwe2_M2_w                         & ! INOUT: cross-spectral propagation speed on nodes
                        , iwe2_fM2                          & ! IN:    M2 tidal surface forcing
                        , iwe2_M2_tau                       & ! IN:    M2 dissipation timescale per node
                        , iwe2_dphit                        & ! IN:    angular bin widths
                        , iwe2_gradxy_e                     & ! INOUT: energy gradient on elements (work array)
                        , iwe2_gradxy_n                     & ! INOUT: energy gradient on nodes (work array)
                        , iwe2_flx_uv                       & ! INOUT: horizontal edge flux (work array)
                        , iwe2_flx_w                        & ! INOUT: cross-spectral flux (work array)
                        , iwe2_refl_bin                     & ! IN:    mirror bin lookup table
                        , iwe2_refl_sgn                     & ! IN:    coast edge sign
                        , iwe2_refl_src                     & ! INOUT: reflect source accumulator
                        , iwe2_refl_node                    & ! IN:    coast-adjacent node mask
                        , vol_nodB2T                        & ! IN:    scalar cell volumes bottom-to-top
                        , partit                            & ! INOUT
                        , mesh                              & ! IN
                        , idemix2_enable_AB                 & ! IN:    enable AB2 time stepping
                        )
                end if
                ! Compute global total NIW energy for conservation check
                ! call check_global_energy(iwe2_E_M2, iwe2_tip1, iwe2_dphit, vol_nodB2T, iwe2_fM2, iwe2_M2_tau, partit, mesh, 'M2')
                ! call check_flux_conservation(  iwe2_E_niw_divh, iwe2_E_niw_divs, iwe2_ti, iwe2_dphit, vol_nodB2T, partit, mesh, 'M2')
            end if
            t4 = MPI_Wtime()
            time_M2_integ = t4 - t3
            
            t4 = MPI_Wtime()
            if (idemix2_enable_niw) then
                ! create complete Ecomaprt diagnostic
                if (idemix2_diag_Ecompart) then
                    call hsintegrate_Ecompart(              &
                        iwe2_ti, iwe2_tip1                 & ! IN:    current and next time indices
                        , 'NIW'                             & ! IN:    compartment name
                        , iwe2_E_niw                        & ! INOUT: energy array [2,nfbin,node]
                        , iwe2_E_niw_divh                   & ! INOUT: horizontal divergence flux accumulator
                        , iwe2_E_niw_divs                   & ! INOUT: cross-spectral divergence accumulator
                        , iwe2_niw_uv                       & ! IN:    horizontal group velocities on elements
                        , iwe2_niw_w                        & ! INOUT: cross-spectral propagation speed on nodes
                        , iwe2_fniw                         & ! IN:    NIW surface forcing
                        , iwe2_niw_tau                      & ! IN:    NIW dissipation timescale per node
                        , iwe2_dphit                        & ! IN:    angular bin widths
                        , iwe2_gradxy_e                     & ! INOUT: energy gradient on elements (work array)
                        , iwe2_gradxy_n                     & ! INOUT: energy gradient on nodes (work array)
                        , iwe2_flx_uv                       & ! INOUT: horizontal edge flux (work array)
                        , iwe2_flx_w                        & ! INOUT: cross-spectral flux (work array)
                        , iwe2_refl_bin                     & ! IN:    mirror bin lookup table
                        , iwe2_refl_sgn                     & ! IN:    coast edge sign
                        , iwe2_refl_src                     & ! INOUT: reflect source accumulator
                        , iwe2_refl_node                    & ! IN:    coast-adjacent node mask
                        , vol_nodB2T                        & ! IN:    scalar cell volumes bottom-to-top
                        , partit                            & ! INOUT
                        , mesh                              & ! IN
                        , idemix2_enable_AB                 & ! IN:    enable AB2 time stepping
                        , iwe2_E_niw_dt                     & ! INOUT: optional diagnostic: total tendency
                        , iwe2_E_niw_advh                   & ! INOUT: optional diagnostic: horiz. advection
                        , iwe2_E_niw_advs                   & ! INOUT: optional diagnostic: spectral advection
                        , iwe2_E_niw_diss                   & ! INOUT: optional diagnostic: dissipation
                        , iwe2_E_niw_forc                   & ! INOUT: optional diagnostic: forcing
                        , iwe2_E_niw_refl                   & ! INOUT: optional diagnostic: coast reflection
                        )
                else
                call hsintegrate_Ecompart(                  &
                        iwe2_ti, iwe2_tip1                 & ! IN:    current and next time indices
                        , 'NIW'                             & ! IN:    compartment name
                        , iwe2_E_niw                        & ! INOUT: energy array [2,nfbin,node]
                        , iwe2_E_niw_divh                   & ! INOUT: horizontal divergence flux accumulator
                        , iwe2_E_niw_divs                   & ! INOUT: cross-spectral divergence accumulator
                        , iwe2_niw_uv                       & ! IN:    horizontal group velocities on elements
                        , iwe2_niw_w                        & ! INOUT: cross-spectral propagation speed on nodes
                        , iwe2_fniw                         & ! IN:    NIW surface forcing
                        , iwe2_niw_tau                      & ! IN:    NIW dissipation timescale per node
                        , iwe2_dphit                        & ! IN:    angular bin widths
                        , iwe2_gradxy_e                     & ! INOUT: energy gradient on elements (work array)
                        , iwe2_gradxy_n                     & ! INOUT: energy gradient on nodes (work array)
                        , iwe2_flx_uv                       & ! INOUT: horizontal edge flux (work array)
                        , iwe2_flx_w                        & ! INOUT: cross-spectral flux (work array)
                        , iwe2_refl_bin                     & ! IN:    mirror bin lookup table
                        , iwe2_refl_sgn                     & ! IN:    coast edge sign
                        , iwe2_refl_src                     & ! INOUT: reflect source accumulator
                        , iwe2_refl_node                    & ! IN:    coast-adjacent node mask
                        , vol_nodB2T                        & ! IN:    scalar cell volumes bottom-to-top
                        , partit                            & ! INOUT
                        , mesh                              & ! IN
                        , idemix2_enable_AB                 & ! IN:    enable AB2 time stepping
                        )
                end if 
                ! Compute global total NIW energy for conservation check
                ! call check_global_energy(iwe2_E_niw, iwe2_tip1, iwe2_dphit, vol_nodB2T, iwe2_fniw, iwe2_niw_tau, partit, mesh, 'niw')
                ! call check_flux_conservation(  iwe2_E_niw_divh, iwe2_E_niw_divs, iwe2_ti, iwe2_dphit, vol_nodB2T, partit, mesh, 'niw')
            end if
            t5 = MPI_Wtime()
            time_niw_integ = t5 - t4
        
        end if ! --> if (idemix2_enable_M2 .or. idemix2_enable_niw) then 
        
        
        !_______________________________________________________________________
        ! 4. Integrate IDEMIX equation vertical, solve vertical diffusion and
        t6 = MPI_Wtime()
        ! dissipation part implicitly
        ! Eiw^(t+1) = Eiw^(t) + dt*[  d/dz( c_0 * tau_v * d/dz(c_0*E_iw))^(t+1)
        !                           - alpha_c*Eiw^(t+1)
        !                           + Forc^(t) ]
!$OMP PARALLEL DO PRIVATE(node, uln, nln)
        do node = 1, myDim_nod2D+eDim_nod2D
            uln = ulevels_nod2D(node)
            nln = nlevels_nod2D(node)
            if (idemix2_diag_Eiw) then 
                call cvmix_idemix2_compute_vdiff_vdiss_Eiw(             &
                      nlev     = nln-uln+1                              & !IN:    number of active levels
                    , dzw      = hnode(           uln:nln-1, node)      & !IN:    layer thickness (m)
                    , dt       = dt                                     & !IN:    time step (s)
                    , c0       = iwe2_c0(         uln:nln  , node)      & !IN:    first baroclinic mode speed (m/s)
                    , alpha_c  = iwe2_alpha_c(    uln:nln  , node)      & !IN:    wave-wave dissipation rate (1/s)
                    , fsrf     = iwe2_fsrf(                  node)      & !IN:    surface IW energy flux (m²/s³)
                    , fbot     = iwe2_fbot(       uln:nln  , node)      & !IN:    bottom tidal energy flux (m²/s³)
                    , Eiw_maxthresh = idemix2_Eiw_maxthresh             & !IN:
                    , Eiw_old  = iwe2_E_iw(uln:nln  , node, iwe2_ti)    & !IN:    IW energy at current time (m²/s²)
                    , Eiw_new  = iwe2_E_iw(uln:nln  , node, iwe2_tip1)  & !OUT:   IW energy at next time (m²/s²)
                    , Eiw_diss = iwe2_E_iw_diss(  uln:nln  , node)      & !OUT:   IW dissipation rate for TKE coupling (m²/s³)
                    , Eiw_dt   = iwe2_E_iw_dt(    uln:nln  , node)      & !OUT:   optional diagnostic: total tendency
                    , Eiw_vdif = iwe2_E_iw_vdif(  uln:nln  , node)      & !OUT:   optional diagnostic: vertical diffusion tendency
                    , Eiw_srf  = iwe2_E_iw_fsrf(             node)      & !OUT:   optional diagnostic: surface forcing
                    , Eiw_bot  = iwe2_E_iw_fbot(  uln:nln  , node)      & !OUT:   optional diagnostic: bottom forcing
                )
            else
                call cvmix_idemix2_compute_vdiff_vdiss_Eiw(             &
                      nlev     = nln-uln+1                              & !IN:    number of active levels
                    , dzw      = hnode(           uln:nln-1, node)      & !IN:    layer thickness (m)
                    , dt       = dt                                     & !IN:    time step (s)
                    , c0       = iwe2_c0(         uln:nln  , node)      & !IN:    first baroclinic mode speed (m/s)
                    , alpha_c  = iwe2_alpha_c(    uln:nln  , node)      & !IN:    wave-wave dissipation rate (1/s)
                    , fsrf     = iwe2_fsrf(                  node)      & !IN:    surface IW energy flux (m²/s³)
                    , fbot     = iwe2_fbot(       uln:nln  , node)      & !IN:    bottom tidal energy flux (m²/s³)
                    , Eiw_maxthresh = idemix2_Eiw_maxthresh             & !IN:
                    , Eiw_old  = iwe2_E_iw(uln:nln, node, iwe2_ti)      & !IN:    IW energy at current time (m²/s²)
                    , Eiw_new  = iwe2_E_iw(uln:nln, node, iwe2_tip1)    & !OUT:   IW energy at next time (m²/s²)
                    , Eiw_diss = iwe2_E_iw_diss(  uln:nln  , node)      & !OUT:   IW dissipation rate for TKE coupling (m²/s³)
                    )
            end if     
        end do ! --> do node = 1, myDim_nod2D
!$OMP END PARALLEL DO
        t7 = MPI_Wtime()
        time_Eiw_vdiff = t7 - t6

        !_______________________________________________________________________
        ! Exchange E_iw(..., tp1) after vertical diffusion - needed for horizontal diffusion gradients
        if (idemix2_enable_hor_diff_expl .or. idemix2_enable_hor_diff_impl_iter) then
            call exchange_nod(iwe2_E_iw(:, :, iwe2_tip1), partit)
        end if 
        
        !_______________________________________________________________________
        ! 5+6. Horizontal Eiw diffusion (explicit single-pass or iterative implicit)
        ! Eiw^(t+1) += div_h( v_0 * tau_h * grad_h(v_0 * Eiw) )
        t7 = MPI_Wtime()
        ! Diagnostic: capture negated pre-diffusion iwe2_tip1 state.
        ! The before/after approach is used for both modes because the iterative
        ! implicit call aliases Eiw_old and Eiw to the same array, making any
        ! in-subroutine Eiw-Eiw_old difference zero (Fortran aliasing bug).
        if ((idemix2_enable_hor_diff_expl .or. idemix2_enable_hor_diff_impl_iter) .and. idemix2_diag_Eiw) then
            iwe2_E_iw_hdif = -iwe2_E_iw(:, :, iwe2_tip1)
        end if

        if (idemix2_enable_hor_diff_expl) then
            ! Exchange E_iw(...,ti) halo - needed for explicit horizontal diffusion gradients
            call exchange_nod(iwe2_E_iw(:, :, iwe2_ti), partit)
            call compute_hdiff_Eiw(         &
                  iwe2_E_iw(:, :, iwe2_ti)  & ! IN:    IW energy at current time
                , iwe2_E_iw(:, :, iwe2_tip1)& ! INOUT: IW energy updated by horizontal diffusion
                , iwe2_v0                   & ! IN:    first baroclinic mode group speed (m/s)
                , 1                         & ! IN:    number of horizontal diffusion iterations
                , partit                    & ! INOUT
                , mesh                      & ! IN
                )
        end if

        if (idemix2_enable_hor_diff_impl_iter) then
            do iter=1, idemix2_hor_diff_niter
                call compute_hdiff_Eiw(             &
                      iwe2_E_iw(:, :, iwe2_tip1)    & ! IN:    IW energy from previous iteration
                    , iwe2_E_iw(:, :, iwe2_tip1)    & ! INOUT: IW energy updated in-place (iterative)
                    , iwe2_v0                       & ! IN:    first baroclinic mode group speed (m/s)
                    , idemix2_hor_diff_niter        & ! IN:    total number of diffusion iterations
                    , partit                        & ! INOUT
                    , mesh                          & ! IN
                    )
                call exchange_nod(iwe2_E_iw(:, :, iwe2_tip1), partit)
            end do
        end if

        ! Diagnostic: iwe2_E_iw_hdif = (Eiw_after - Eiw_before) / dt  [m²/s³]
        if ((idemix2_enable_hor_diff_expl .or. idemix2_enable_hor_diff_impl_iter) .and. idemix2_diag_Eiw) then
            do node = 1, myDim_nod2D
                uln = ulevels_nod2D(node)
                nln = nlevels_nod2D(node)
                iwe2_E_iw_hdif(uln:nln, node) = (iwe2_E_iw_hdif(uln:nln, node) &
                                                + iwe2_E_iw(uln:nln, node, iwe2_tip1)) / dt
                iwe2_E_iw_dt(uln:nln, node)   =  iwe2_E_iw_dt(uln:nln, node) &
                                                + iwe2_E_iw_hdif(uln:nln, node)
            end do
        end if
        t8 = MPI_Wtime()
        time_Eiw_hdiff = t8 - t7
        
        !_______________________________________________________________________
        ! 7. compute wave-wave interaction
        t8 = MPI_Wtime() 
        if (idemix2_enable_M2 .or. idemix2_enable_niw) then
!$OMP PARALLEL DO PRIVATE(node, uln, nln)
            do node = 1, myDim_nod2D+eDim_nod2D
                uln = ulevels_nod2D(node)
                nln = nlevels_nod2D(node)
                !_______________________________________________________________
                if (idemix2_enable_M2 .and. idemix2_enable_niw) then
                    if (idemix2_diag_WWI) then
                        call cvmix_idemix2_compute_Eiw_waveinteract(                  &
                              nlev          = nln-uln+1                               & !IN:    number of active levels
                            , nfbin         = idemix2_nfbin                           & !IN:    number of spectral angular bins
                            , dzw           = hnode(                uln:nln-1, node)  & !IN:    layer thickness (m)
                            , dphi          = iwe2_dphit                              & !IN:    angular bin widths
                            , dt            = dt                                      & !IN:    time step (s)
                            , flag_posdef   = .True.                                  & !IN:    enforce positive-definiteness
                            , E_iw_old      = iwe2_E_iw( uln:nln, node, iwe2_ti)      & !IN:    IDEMIX1 IW energy at current time (m²/s²)
                            , E_iw_new      = iwe2_E_iw( uln:nln, node, iwe2_tip1)    & !INOUT: IDEMIX1 IW energy updated by wave-wave interaction
                            , E_M2_old      = iwe2_E_M2( :      , node, iwe2_ti)      & !IN:    M2 compartment energy at current time (m²/s²)
                            , E_M2_new      = iwe2_E_M2( :      , node, iwe2_tip1)    & !INOUT: M2 compartment energy updated
                            , E_M2_struct   = iwe2_E_M2_struct(       :      , node)  & !IN:    M2 vertical structure function
                            , alpha_M2_c    = iwe2_alpha_M2_c(                 node)  & !IN:    M2 continuous dissipation rate
                            , tau_M2        = iwe2_M2_tau(                     node)  & !IN:    M2 interaction timescale (s)
                            , E_niw_old     = iwe2_E_niw(:      , node, iwe2_ti)      & !IN:    NIW compartment energy at current time (m²/s²)
                            , E_niw_new     = iwe2_E_niw(:      , node, iwe2_tip1)    & !INOUT: NIW compartment energy updated
                            , E_niw_struct  = iwe2_E_niw_struct(      :      , node)  & !IN:    NIW vertical structure function
                            , tau_niw       = iwe2_niw_tau(                    node)  & !IN:    NIW interaction timescale (s)
                            , E_iw_dt       = iwe2_E_iw_dt(           :      , node)  & !INOUT: optional diagnostic: IDEMIX1 tendency from wave-wave
                            , E_iw_diss_M2  = iwe2_E_iw_diss_M2(      :      , node)  & !INOUT: optional diagnostic: IDEMIX1 dissipation from M2
                            , E_iw_diss_niw = iwe2_E_iw_diss_niw(     :      , node)  & !INOUT: optional diagnostic: IDEMIX1 dissipation from NIW
                            , E_M2_dt       = iwe2_E_M2_dt(           :      , node)  & !INOUT: optional diagnostic: M2 tendency from wave-wave
                            , E_M2_diss_wwi = iwe2_E_M2_diss_wwi(     :      , node)  & !INOUT: optional diagnostic: M2 dissipation from wave-wave
                            )
                    else
                        call cvmix_idemix2_compute_Eiw_waveinteract(                  &
                              nlev          = nln-uln+1                               & !IN:    number of active levels
                            , nfbin         = idemix2_nfbin                           & !IN:    number of spectral angular bins
                            , dzw           = hnode(                uln:nln-1, node)  & !IN:    layer thickness (m)
                            , dphi          = iwe2_dphit                              & !IN:    angular bin widths
                            , dt            = dt                                      & !IN:    time step (s)
                            , flag_posdef   = .True.                                  & !IN:    enforce positive-definiteness
                            , E_iw_old      = iwe2_E_iw( uln:nln, node, iwe2_ti)      & !IN:    IDEMIX1 IW energy at current time (m²/s²)
                            , E_iw_new      = iwe2_E_iw( uln:nln, node, iwe2_tip1)    & !INOUT: IDEMIX1 IW energy updated by wave-wave interaction
                            , E_M2_old      = iwe2_E_M2( :      , node, iwe2_ti)      & !IN:    M2 compartment energy at current time (m²/s²)
                            , E_M2_new      = iwe2_E_M2( :      , node, iwe2_tip1)    & !INOUT: M2 compartment energy updated
                            , E_M2_struct   = iwe2_E_M2_struct(       :      , node)  & !IN:    M2 vertical structure function
                            , alpha_M2_c    = iwe2_alpha_M2_c(                 node)  & !IN:    M2 continuous dissipation rate
                            , tau_M2        = iwe2_M2_tau(                     node)  & !IN:    M2 interaction timescale (s)
                            , E_niw_old     = iwe2_E_niw(:      , node, iwe2_ti)      & !IN:    NIW compartment energy at current time (m²/s²)
                            , E_niw_new     = iwe2_E_niw(:      , node, iwe2_tip1)    & !INOUT: NIW compartment energy updated
                            , E_niw_struct  = iwe2_E_niw_struct(      :      , node)  & !IN:    NIW vertical structure function
                            , tau_niw       = iwe2_niw_tau(                    node)  & !IN:    NIW interaction timescale (s)
                            )
                    end if
                    
                elseif (idemix2_enable_M2) then
                    if (idemix2_diag_WWI) then
                        call cvmix_idemix2_compute_Eiw_waveinteract(                  &
                              nlev          = nln-uln+1                               & !IN:    number of active levels
                            , nfbin         = idemix2_nfbin                           & !IN:    number of spectral angular bins
                            , dzw           = hnode(                uln:nln-1, node)  & !IN:    layer thickness (m)
                            , dphi          = iwe2_dphit                              & !IN:    angular bin widths
                            , dt            = dt                                      & !IN:    time step (s)
                            , flag_posdef   = .True.                                  & !IN:    enforce positive-definiteness
                            , E_iw_old      = iwe2_E_iw( uln:nln, node, iwe2_ti)      & !IN:    IDEMIX1 IW energy at current time (m²/s²)
                            , E_iw_new      = iwe2_E_iw( uln:nln, node, iwe2_tip1)    & !INOUT: IDEMIX1 IW energy updated by wave-wave interaction
                            , E_M2_old      = iwe2_E_M2( :      , node, iwe2_ti)      & !IN:    M2 compartment energy at current time (m²/s²)
                            , E_M2_new      = iwe2_E_M2( :      , node, iwe2_tip1)    & !INOUT: M2 compartment energy updated
                            , E_M2_struct   = iwe2_E_M2_struct(       :      , node)  & !IN:    M2 vertical structure function
                            , alpha_M2_c    = iwe2_alpha_M2_c(                 node)  & !IN:    M2 continuous dissipation rate
                            , tau_M2        = iwe2_M2_tau(                     node)  & !IN:    M2 interaction timescale (s)
                            , E_iw_dt       = iwe2_E_iw_dt(           :      , node)  & !INOUT: optional diagnostic: IDEMIX1 tendency from wave-wave
                            , E_iw_diss_M2  = iwe2_E_iw_diss_M2(      :      , node)  & !INOUT: optional diagnostic: IDEMIX1 dissipation from M2
                            , E_M2_dt       = iwe2_E_M2_dt(           :      , node)  & !INOUT: optional diagnostic: M2 tendency from wave-wave
                            , E_M2_diss_wwi = iwe2_E_M2_diss_wwi(     :      , node)  & !INOUT: optional diagnostic: M2 dissipation from wave-wave
                            )
                    else
                        call cvmix_idemix2_compute_Eiw_waveinteract(                  &
                              nlev          = nln-uln+1                               & !IN:    number of active levels
                            , nfbin         = idemix2_nfbin                           & !IN:    number of spectral angular bins
                            , dzw           = hnode(                uln:nln-1, node)  & !IN:    layer thickness (m)
                            , dphi          = iwe2_dphit                              & !IN:    angular bin widths
                            , dt            = dt                                      & !IN:    time step (s)
                            , flag_posdef   = .True.                                  & !IN:    enforce positive-definiteness
                            , E_iw_old      = iwe2_E_iw( uln:nln, node, iwe2_ti)      & !IN:    IDEMIX1 IW energy at current time (m²/s²)
                            , E_iw_new      = iwe2_E_iw( uln:nln, node, iwe2_tip1)    & !INOUT: IDEMIX1 IW energy updated by wave-wave interaction
                            , E_M2_old      = iwe2_E_M2( :      , node, iwe2_ti)      & !IN:    M2 compartment energy at current time (m²/s²)
                            , E_M2_new      = iwe2_E_M2( :      , node, iwe2_tip1)    & !INOUT: M2 compartment energy updated
                            , E_M2_struct   = iwe2_E_M2_struct(       :      , node)  & !IN:    M2 vertical structure function
                            , alpha_M2_c    = iwe2_alpha_M2_c(                 node)  & !IN:    M2 continuous dissipation rate
                            , tau_M2        = iwe2_M2_tau(                     node)  & !IN:    M2 interaction timescale (s)
                            )
                    end if
                elseif (idemix2_enable_niw) then
                    if (idemix2_diag_WWI) then
                        call cvmix_idemix2_compute_Eiw_waveinteract(                  &
                              nlev          = nln-uln+1                               & !IN:    number of active levels
                            , nfbin         = idemix2_nfbin                           & !IN:    number of spectral angular bins
                            , dzw           = hnode(                uln:nln-1, node)  & !IN:    layer thickness (m)
                            , dphi          = iwe2_dphit                              & !IN:    angular bin widths
                            , dt            = dt                                      & !IN:    time step (s)
                            , flag_posdef   = .True.                                  & !IN:    enforce positive-definiteness
                            , E_iw_old      = iwe2_E_iw( uln:nln, node, iwe2_ti)      & !IN:    IDEMIX1 IW energy at current time (m²/s²)
                            , E_iw_new      = iwe2_E_iw( uln:nln, node, iwe2_tip1)    & !INOUT: IDEMIX1 IW energy updated by wave-wave interaction
                            , E_niw_old     = iwe2_E_niw(:      , node, iwe2_ti)      & !IN:    NIW compartment energy at current time (m²/s²)
                            , E_niw_new     = iwe2_E_niw(:      , node, iwe2_tip1)    & !INOUT: NIW compartment energy updated
                            , E_niw_struct  = iwe2_E_niw_struct(      :      , node)  & !IN:    NIW vertical structure function
                            , tau_niw       = iwe2_niw_tau(                    node)  & !IN:    NIW interaction timescale (s)
                            , E_iw_dt       = iwe2_E_iw_dt(           :      , node)  & !INOUT: optional diagnostic: IDEMIX1 tendency from wave-wave
                            , E_iw_diss_niw = iwe2_E_iw_diss_niw(     :      , node)  & !INOUT: optional diagnostic: IDEMIX1 dissipation from NIW
                            )
                    else
                        call cvmix_idemix2_compute_Eiw_waveinteract(                  &
                              nlev          = nln-uln+1                               & !IN:    number of active levels
                            , nfbin         = idemix2_nfbin                           & !IN:    number of spectral angular bins
                            , dzw           = hnode(                uln:nln-1, node)  & !IN:    layer thickness (m)
                            , dphi          = iwe2_dphit                              & !IN:    angular bin widths
                            , dt            = dt                                      & !IN:    time step (s)
                            , flag_posdef   = .True.                                  & !IN:    enforce positive-definiteness
                            , E_iw_old      = iwe2_E_iw( uln:nln, node, iwe2_ti)      & !IN:    IDEMIX1 IW energy at current time (m²/s²)
                            , E_iw_new      = iwe2_E_iw( uln:nln, node, iwe2_tip1)    & !INOUT: IDEMIX1 IW energy updated by wave-wave interaction
                            , E_niw_old     = iwe2_E_niw(:      , node, iwe2_ti)      & !IN:    NIW compartment energy at current time (m²/s²)
                            , E_niw_new     = iwe2_E_niw(:      , node, iwe2_tip1)    & !INOUT: NIW compartment energy updated
                            , E_niw_struct  = iwe2_E_niw_struct(      :      , node)  & !IN:    NIW vertical structure function
                            , tau_niw       = iwe2_niw_tau(                    node)  & !IN:    NIW interaction timescale (s)
                            )
                    end if
                end if ! -->  (idemix2_enable_M2 .and. idemix2_enable_niw) then
            end do ! --> for node = 1, myDim_nod2D
!$OMP END PARALLEL DO
        end if ! --> if (idemix2_enable_M2 .or . idemix2_enable_niw) then
        t9 = MPI_Wtime()
        time_waveint = t9 - t8

        !_______________________________________________________________________
        ! 8. write IDEMIX2 diffusivities and viscositie to FESOM only when IDEMIX2 is
        t9 = MPI_Wtime() 
        ! used alone --> mostly for debuging --> otherwise TKE Av and Kv are use
        if(mix_scheme_nmb==7) then 
            
            !___________________________________________________________________
            ! write out diffusivity --> convert from elem to vertices
!$OMP PARALLEL DO PRIVATE(node, uln, nln)
            do node=1, myDim_nod2D
                uln = ulevels_nod2D(node)
                nln = nlevels_nod2D(node)
                !_______________________________________________________________
                ! convert Eiw disspation into Kv and Av on vertices 
                call cvmix_idemix2_compute_Eiw_diss2KvAv(         &
                          nlev    = nln-uln+1                     & !IN:    number of active levels
                        , Eiw_diss= iwe2_E_iw_diss(uln:nln, node) & !IN:    IW dissipation rate (m²/s³)
                        , Nsqr    = bvfreq(        uln:nln, node) & !IN:    buoyancy frequency squared (1/s²)
                        , KappaH  = Kv(            uln:nln, node) & !INOUT: diapycnal diffusivity (m²/s)
                        , KappaM  = iwe2_Av(       uln:nln, node) & !INOUT: vertical viscosity (m²/s)
                        )
            end do
!$OMP END PARALLEL DO
            call exchange_nod(iwe2_Av, partit)    
            call exchange_nod(Kv , partit)
            
            !___________________________________________________________________
            ! more idemix2 viscosity from vertices to elements
            do elem=1,myDim_elem2D+eDim_elem2D
                uln = ulevels(elem)
                nln = nlevels(elem)
                elnodes = elem2d_nodes(:,elem)
                do nz=uln, nln
                    Av(nz, elem) = sum(iwe2_Av(nz, elnodes))/3.0_WP
                end do
            end do
            
        end if
        t10 = MPI_Wtime()
        time_Kv_Av = t10 - t9
        
        !_______________________________________________________________________
        ! Timing report
        t_end = MPI_Wtime()
        time_total = t_end - t_start
        
        if (mype==0 .and. mod(istep, logfile_outfreq)==0) then
            write(*,*)
            write(*,'(A)') ' ___IDEMIX2 EXECUTION TIMES_____________________________'
            write(*,'(A, ES10.3, A, F6.2, A)') '     Params (cn,tau,struct)    : ', time_params,    ' s  (', 100.0*time_params/time_total,    '%)'
            write(*,'(A, ES10.3, A, F6.2, A)') '     Group velocities (u,v,w)  : ', time_groupvel,  ' s  (', 100.0*time_groupvel/time_total,  '%)'
            write(*,'(A, ES10.3, A, F6.2, A)') '     M2 spectral integration   : ', time_M2_integ,  ' s  (', 100.0*time_M2_integ/time_total,  '%)'
            write(*,'(A, ES10.3, A, F6.2, A)') '     NIW spectral integration  : ', time_niw_integ, ' s  (', 100.0*time_niw_integ/time_total, '%)'
            write(*,'(A, ES10.3, A, F6.2, A)') '     Eiw vertical diffusion    : ', time_Eiw_vdiff, ' s  (', 100.0*time_Eiw_vdiff/time_total, '%)'
            write(*,'(A, ES10.3, A, F6.2, A)') '     Eiw horizontal diffusion  : ', time_Eiw_hdiff, ' s  (', 100.0*time_Eiw_hdiff/time_total, '%)'
            write(*,'(A, ES10.3, A, F6.2, A)') '     Wave-wave interaction     : ', time_waveint,   ' s  (', 100.0*time_waveint/time_total,   '%)'
            write(*,'(A, ES10.3, A, F6.2, A)') '     Kv/Av output              : ', time_Kv_Av,     ' s  (', 100.0*time_Kv_Av/time_total,     '%)'
            write(*,'(A)') ' _______________________________________'
            write(*,'(A, ES10.3,A)')        '     IDEMIX2 TOTAL time        : ', time_total, ' s'
            write(*,*)
            write(*,*)
        end if
        
    end subroutine calc_cvmix_idemix2    

    
    
    !
    !
    !___________________________________________________________________________
    ! horizontal superbee advection of spectral bins
    subroutine adv_Ecompart_hor_spctrl_superbee(&
                  vel                           & 
                , ttf                           &
                , ttf_grad_n                    &
                , flux                          &
                , partit                        &
                , mesh                          &
                , flag_2ndord_time              &
                , flag_posdef                   &
                )

        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit)  , intent(inout), target :: partit
        type(t_mesh)    , intent(in)   , target :: mesh
        real(kind=WP)   , intent(in)            :: ttf(          idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(in)            :: ttf_grad_n(2, idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(in)            :: vel(2,        idemix2_nfbin, partit%myDim_elem2D+partit%eDim_elem2D)
        real(kind=WP)   , intent(inout)         :: flux(         idemix2_nfbin, partit%myDim_edge2D)
        logical         , intent(in)            :: flag_2ndord_time
        logical         , intent(in)            :: flag_posdef
        
        !___LOCAL VARIABLES_____________________________________________________
        real(kind=WP)                           :: dx1, dy1, dx2, dy2, dxdy12(2), dt_over_edlen, n_x, n_y, n_len  
        real(kind=WP)                           :: u1, u2, v1, v2, Ue, CFL, dh1, dh2
        real(kind=WP)                           :: R, ttfp2, ttfm1, ttf0, ttfp1, dttf0p1, Tmean1, Tmean2, Cr, vflux, vfabs
        integer                                 :: el(2), el2, ednodes(2), edge, fbini, nfbin, nln, uln, nz
        logical, parameter                      :: flag_volfix=.True.
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        
        !_______________________________________________________________________
        nfbin=idemix2_nfbin
        
        ! this advection does !!! NOT !!! go over the vertical dimension it goes 
        ! over the domain of the spectral bins 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, ednodes, el, el2, fbini, nln, &
!$OMP                                  dx1, dy1, dx2, dy2, dxdy12, dh1, dh2, n_x, n_y, n_len, dt_over_edlen, &
!$OMP                                  u1, v1, u2, v2, Ue, CFL, vflux, vfabs, &
!$OMP                                  ttf0, ttfp1, dttf0p1, ttfp2, ttfm1, R, Cr, Tmean1, Tmean2)
!$OMP DO        
        do edge=1, myDim_edge2D
            !___________________________________________________________________
            ! local indice of nodes that span up edge ed
            ednodes= edges(:,edge)
            
            ! local index of element that contribute to edge
            el     = edge_tri(:,edge)
            
            ! edge_cross_dxdy(1:2,ed)... dx,dy distance from center of edge to
            ! element centroid el(1) --> needed to calc flux perpedicular to edge from elem el(1)
            dx1    = edge_cross_dxdy(1,edge)
            dy1    = edge_cross_dxdy(2,edge)
            
            ! length of edge dx dy
            dxdy12 = edge_dxdy(:,edge)*r_earth  
            
            ! compute elemental local bottom depth for volume flux computation 
            if (flag_volfix) then 
                nln = nlevels(el(1))
                dh1 = abs(zbar(nln))
            else
                uln = ulevels(el(1))
                nln = nlevels(el(1))-1
                dh1  = 0.0_WP
                do nz=uln,nln
                    dh1 = dh1 + helem(nz, el(1))
                end do
            end if 
            
            !___________________________________________________________________
            ! same parameter but for other element el(2) that contributes to edge ed
            ! if el(2)==0 than edge is boundary edge
            dx2    = 0.0_WP
            dy2    = 0.0_WP
            dh2    = 0.0_WP
            el2    = el(1) ! --> just a dummy index here is vel(el2) are multiplied with 0.0_WP
            if(el(2)>0) then
                el2      = el(2) ! --> dummy replaced with valid el(2), dx2 and dy2
                dx2      = edge_cross_dxdy(3,edge)
                dy2      = edge_cross_dxdy(4,edge)
                
                ! compute elemental local bottom depth for volume flux computation 
                if (flag_volfix) then 
                    nln = nlevels(el2)
                    dh2 = abs(zbar(nln))
                else
                    uln = ulevels(el2)
                    nln = nlevels(el2)-1
                    do nz=uln,nln
                        dh2 = dh2 + helem(nz, el2)
                    end do
                end if 
                
                dxdy12(1)= dxdy12(1) * (elem_cos(el(1)) + elem_cos(el2)) * 0.5_WP
                
                ! compute mean edge-segment normal vector
                !                         (dx1,dy1)
                !                             ^ 
                !                             │
                !                             ├──> (dy1,-dx1)  
                !                             │
                !                ●━━━━━━━━━━━━┿━━━━━━━━━━►● 
                !                             │
                !               (dy2,-dx2) <──┼──> (-dy2,dx2)  
                !                             │ 
                !                             v  
                !                         (dx2,dy2)
                ! 
                n_x      = (  dy1  + (-dy2))/2
                n_y      = ((-dx1) +   dx2 )/2
                
            else
                dxdy12(1)= dxdy12(1) * elem_cos(el(1))   
                ! compute edge-segment normal vector
                n_x      =  dy1
                n_y      = -dx1
            end if
            dt_over_edlen = dt/sqrt(dxdy12(1)**2 + dxdy12(2)**2)
            
            ! Normalize normal vector to unit length for CFL calculation
            n_len = sqrt(n_x**2 + n_y**2)
            n_x   = n_x / n_len
            n_y   = n_y / n_len
            
            !___________________________________________________________________
            !loop over spectral bins
            do fbini = 2, nfbin-1
                !_______________________________________________________________
                ! compute volume flux across the segments from el(1)
                ! NOTE: NO dphi here - matches pyOM2 (dphi only in cross-spectral divergence)
                u1     = vel(1, fbini, el(1))
                v1     = vel(2, fbini, el(1))
                u2     = vel(1, fbini, el2  )
                v2     = vel(2, fbini, el2  )
                
                vflux  =  (-v1*dx1 + u1*dy1) * dh1 &
                        + ( v2*dx2 - u2*dy2) * dh2 ! NO * dphi(fbini)
                
                ! compute approximated edge centered and along edge projected 
                ! mean velocity --> need to add component to make second order
                ! in time
                Ue    = sqrt((0.5_WP*( u1 + u2 )*n_x)**2 +  (0.5_WP*( v1 + v2 )*n_y)**2)
                CFL   = min(1.0_WP, Ue*dt_over_edlen)
                
                !_______________________________________________________________
                ! compute tracer difference allong edge
                ttf0    = ttf(fbini, ednodes(1))
                ttfp1   = ttf(fbini, ednodes(2))
                dttf0p1 = (ttfp1-ttf0)
                
                !_______________________________________________________________
                ! tracer Slope Ratio Calculation for upwind augmented point
                !
                !                             o
                !            ○._            .´ `.            .-○
                !            |  `._       .´  ^ dx1,dy1   .-´  |
                !            |     `._  .´    ├─>  `.  .-´     |
                !            |   □    `●━━━━━━┿━━━━━►●´----□---|-->●ttf_p2
                !            |      .-´│`.  <─┤    .´│`._      |
                !            |   .-´  dx2,dy2 v  .´  │   `._   |
                !            ○.-´      │    `. .´    │      `._○ 
                !                      │      o      │
                !                      ├------------>┤
                !                      │    dxdy12   │
                !                      │             │
                !                      v             v
                !                    ttf_0      ttf_grad_n_p1
                !                        
                ttfp2   =   ttf0                                                & 
                          + 2.0_WP * dxdy12(1)*ttf_grad_n(1, fbini, ednodes(2)) &
                          + 2.0_WP * dxdy12(2)*ttf_grad_n(2, fbini, ednodes(2))  
                
                ! considering here we want to advect energy we can ensure the 
                ! variable to be positiv definit
                ttfp2   = merge(max(0.0_WP,ttfp2), ttfp2, flag_posdef)
                
                ! compute tracer slope 
                R       = (ttfp1-ttfp2)/(-dttf0p1+small)
                
                ! apply superbee limiter
                Cr      = slimiter_superbee(R)  ! DEBUG: force pure upwind (no reconstruction)
               
                ! construct edge centered tracer value
                ! T_(i+0.5) = T_i+1 - [ (T_i+1 - T_i)/dx ]_Limited * dx/2
                !  --> this is seconds order in space, but first order in time 
                !  
                ! T_(i+0.5) = T_i+1 - [ (T_i+1 - T_i)/dx ]_Limited * (dx/2 - u_ed*n_ed*dt/2)
                ! T_(i+0.5) = T_i+1 - [ (T_i+1 - T_i)/dx ]_Limited * dx/2 *(1 - CFL_h)
                !  --> becomes second order in space and time  (Direct space-time scheme)
                !  --> u_ed  = (vel(el(1)) + vel(el(2)))/2
                !  --> n_ed  = (n_1 + n_2)/2
                !  --> CFL_h = u_ed*n_ed*dt/dx
                Tmean2  = ttfp1 + 0.5_WP * Cr * (-dttf0p1) * (1.0_WP-merge(CFL, 0.0_WP, flag_2ndord_time)) 
                
                !_______________________________________________________________
                ! tracer Slope Ratio Calculation for downwind augmented point
                !
                !                             o
                !            ○._            .´ `.            .-○
                !            |  `._       .´  ^ dx1,dy1   .-´  |
                !            |     `._  .´    ├─>  `.  .-´     |
                !  ttf_m1●<--|---□----`●━━━━━━┿━━━━━►●´    □   |
                !            |      .-´│`.  <─┤    .´│`._      |
                !            |   .-´  dx2,dy2 v  .´  │   `._   |
                !            ○.-´      │    `. .´    │      `._○ 
                !                      │      o      │
                !                      ├------------>┤
                !                      │    dxdy12   │
                !                      │             │
                !                      v             v
                !                 ttf_grad_n_0     ttf_p1 
                !                        
                ttfm1   =   ttfp1                                               & 
                          - 2.0_WP * dxdy12(1)*ttf_grad_n(1, fbini, ednodes(1)) &
                          - 2.0_WP * dxdy12(2)*ttf_grad_n(2, fbini, ednodes(1))  
                
                ! considering here we want to advect energy we can ensure the 
                ! variable to be positiv definit
                ttfm1   = merge(max(0.0_WP,ttfm1), ttfm1, flag_posdef)
                
                ! compute tracer slope 
                R       = (ttf0-ttfm1)/(dttf0p1+small)
                
                ! apply superbee limiter
                Cr      = slimiter_superbee(R)  ! DEBUG: force pure upwind (no reconstruction)
                
                ! construct edge centered tracer value
                ! T_(i+0.5) = T_i + [ (T_i+1 - T_i)/dx ]_Limited * dx/2
                !  --> this is seconds order in space, but first order in time 
                !  
                ! T_(i+0.5) = T_i + [ (T_i+1 - T_i)/dx ]_Limited * (dx/2 - u_ed*n_ed*dt/2)
                ! T_(i+0.5) = T_i + [ (T_i+1 - T_i)/dx ]_Limited * dx/2 *(1 - CFL_h)
                !  --> becomes second order in space and time  (Direct space-time scheme)
                !  --> u_ed  = (vel(el(1)) + vel(el(2)))/2
                !  --> n_ed  = (n_1 + n_2)/2
                !  --> CFL_h = u_ed*n_ed*dt/dx
                Tmean1  = ttf0 + 0.5_WP * Cr * (dttf0p1) * (1.0_WP-merge(CFL, 0.0_WP, flag_2ndord_time))
                
                !_______________________________________________________________
                ! Horizontal flux
                vfabs = abs(vflux)
                flux(fbini, edge) = -0.5_WP*(                    &
                                           (vflux+vfabs)*Tmean1+ &
                                           (vflux-vfabs)*Tmean2  &
                                          )

            end do ! --> for fbini = 1, nfbin 
        end do !--> do edge=1, myDim_edge2D
!$OMP END DO
!$OMP END PARALLEL
    end subroutine adv_Ecompart_hor_spctrl_superbee



    !
    !
    !___________________________________________________________________________
    ! reflective coastal boundary condition for horizontal spectral energy flux
    !
    ! For every interior edge touching exactly one coast node, any flux directed
    ! TOWARD that coast node is blocked and re-injected at the same edge into the
    ! mirror spectral bin (bin whose propagation angle points back into the ocean).
    ! This conserves total energy: the divergence seen by adv_Ecompart_flx2tra_spctrl
    ! has zero net flux at coastal edges.
    !
    ! Sign convention (matches adv_Ecompart_flx2tra_spctrl):
    !   flux > 0  ->  transport from ednodes(2) to ednodes(1)
    !   flux < 0  ->  transport from ednodes(1) to ednodes(2)
    !
    ! iwe2_bc_reflect_sgn(edge) = +1  when ednodes(1) is the coast node
    !                           = -1  when ednodes(2) is the coast node
    ! Transport is INTO the coast when  sgn * flux > 0.
    !___________________________________________________________________________
    !
    ! EACH TIMESTEP — 4 steps
    ! ────────────────────────
    !
    ! STEP 1 │  adv_Ecompart_hor_spctrl_superbee
    !        │
    !        │  Computes flux(fbini, edge) for all bins.
    !        │  In bin fbini, wave energy at n_int propagates toward coast:
    !        │
    !        │    n_coast <━━━━━━━━━━━━━━  n_int
    !        │              flux(fbini) > 0
    !        │              (positive = toward n_coast)
    !        │
    !        │  Superbee enforces positivity on its own fluxes here.
    !
    ! STEP 2 │  apply_reflect_bc_spctrl
    !        │
    !        │  For each flagged edge where sgn * flux(fbini) > 0:
    !        │
    !        │    flux(fbini, edge) = 0          <- coast-directed flux blocked
    !        │    reflect_src(kk, n_int) += sgn * outflux   <- energy saved for n_int
    !        │
    !        │    n_coast   ................  n_int
    !        │              flux = 0  (blocked)    reflect_src(kk) = |outflux|
    !        │
    !        │  Coast node is NEVER written to.
    !        │  sgn * outflux = |outflux| > 0 regardless of edge orientation.
    !
    ! STEP 3 │  adv_Ecompart_flx2tra_spctrl
    !        │
    !        │  Accumulates flux divergence into Edivh, divides by vol_s.
    !        │  Bin fbini flux was zeroed -> n_int keeps its energy in fbini,
    !        │                               n_coast gets nothing.
    !        │
    !        │    Edivh(fbini, n_int)  = 0  (no export toward coast)
    !        │    Edivh(fbini, n_coast)= 0  (nothing received)
    !
    ! STEP 4 │  inject reflect_src into Edivh  (in hsintegrate_Ecompart)
    !        │
    !        │    Edivh(kk, n_int) += reflect_src(kk, n_int) / vol_s(n_int)
    !        │
    !        │  Energy appears as a SOURCE at n_int in the mirror bin kk --
    !        │  the reflected wave now propagates AWAY from the coast:
    !        │
    !        │    n_coast               n_int ━━━━━━━━━━━━>
    !        │                                bin kk  (away from coast)
    !
    ! ENERGY BALANCE
    ! ──────────────
    ! n_coast : gained 0 (fbini blocked) + 0 (not in kk path) = 0  (unchanged)
    ! n_int   : gained 0 in fbini (blocked)
    !           gained |outflux|/vol_s in kk (reflect_src)
    !           net total = 0  ->  spectral redistribution only     (conserved)
    ! fbini --> energy blocked at coast, stays at n_int
    ! kk    <-- same energy re-injected at n_int in mirror direction
    !
    subroutine apply_reflect_bc_spctrl(   flux      &
                                        , refl_src  &
                                        , refl_bin  &
                                        , refl_sgn  &
                                        , vol_s     &
                                        , partit    & 
                                        , mesh)
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in   ), target :: mesh
        real(kind=WP) , intent(inout)         :: flux(idemix2_nfbin, partit%myDim_edge2D)
        integer       , intent(in   )         :: refl_bin(   partit%myDim_edge2D )
        integer       , intent(in   )         :: refl_sgn(   partit%myDim_edge2D ) 
        real(kind=WP) , intent(out  )         :: refl_src(idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP) , intent(in   )         :: vol_s(partit%myDim_nod2D+partit%eDim_nod2D)
        !___LOCAL VARIABLES_____________________________________________________
        integer       :: edge, fbini, kk, nfbin, sgn, n_int
        real(kind=WP) :: outflux, ivols_int
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

        nfbin = idemix2_nfbin

        refl_src(:,:) = 0.0_WP

        ! Zero every coast-directed flux so no energy reaches the coast node.
        ! For each blocked flux in bin fbini, accumulate two terms into reflect_src:
        !   +sgn*outflux in mirror bin kk    (energy gained in reflected direction)
        !   -sgn*outflux in source bin fbini (energy lost from incident direction)
        ! Net: energy is redistributed fbini→kk at n_int with zero total change.
        ! Applied to Edivh after adv_Ecompart_flx2tra_spctrl in hsintegrate_Ecompart.
!$OMP PARALLEL DO PRIVATE(edge, fbini, kk, sgn, n_int, outflux)
        do edge = 1, myDim_edge2D
            kk  = refl_bin(edge)
            ! edge is not connected to coast
            if (kk == 0) cycle

            !   sgn = +1  -->  coast node = edges(1, edge)
            !                  interior   = edges(2, edge)
            !   sgn = -1  -->  coast node = edges(2, edge)
            !                  interior   = edges(1, edge)
            sgn   = refl_sgn(edge)
            
            ! which one is interior coastal edge node
            n_int = merge(edges(2, edge), edges(1, edge), sgn == +1)

            ! local scaler volume
            ivols_int = 1/vol_s(n_int)
            
            do fbini = 2, nfbin-1
                outflux = flux(fbini, edge)*ivols_int
                if (sgn * outflux > 0.0_WP) then
                    ! block flux towards coast by setting to zero
                    flux(fbini, edge) = 0.0_WP
                    if (n_int <= myDim_nod2D) then
!$OMP ATOMIC
                        refl_src(kk,    n_int) = refl_src(kk,    n_int) + sgn * outflux
!$OMP ATOMIC
                        refl_src(fbini, n_int) = refl_src(fbini, n_int) - sgn * outflux
                    end if
                end if
            end do ! --> do fbini = 2, nfbin-1
        end do ! --> do edge = 1, myDim_edge2D
!$OMP END PARALLEL DO

    end subroutine apply_reflect_bc_spctrl



    !
    !
    !___________________________________________________________________________
    ! cross spectral superbee advection across spectral bins
    subroutine adv_Ecompart_crss_spctrl_superbee( &
                  cs                              & ! cross-spectral exchange rate
                , ttf                             & ! Energy compartment @ node
                , dphi                            & ! width spectral bin
                , flux                            & ! cross-spectral flux
                , partit                          &
                , mesh                            &
                , flag_2ndord_time                &
                , flag_posdef                     &
                ! flag_wlimitcfl                  &
                )
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit),intent(inout), target :: partit
        type(t_mesh),  intent(in)   , target :: mesh
        real(kind=WP), intent(inout)         :: cs(  idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP), intent(in)            :: ttf( idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP), intent(in)            :: dphi(idemix2_nfbin)
        real(kind=WP), intent(inout)         :: flux(idemix2_nfbin, partit%myDim_nod2D)
        logical      , intent(in)            :: flag_2ndord_time
        logical      , intent(in)            :: flag_posdef
        ! logical      , intent(in)            :: flag_wlimitcfl
        
        !___LOCAL VARIABLES_____________________________________________________
        integer                              :: node, fbini, nfbin
        real(kind=WP)                        :: R, ttf0, ttfp1, dttf0p1, Tmean1, Tmean2, Cr, vflux, vfabs
        real(kind=WP)                        :: CFL !, cflmax=0.50_WP
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

        !_______________________________________________________________________
        nfbin=idemix2_nfbin
        
        !_______________________________________________________________________
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, fbini, &
!$OMP                                  ttf0, ttfp1, dttf0p1, R, Cr, Tmean1, Tmean2, &
!$OMP                                  CFL, vflux, vfabs)
!$OMP DO
#else
!$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
        !_______________________________________________________________________
        ! this advection does !!! NOT !!! go over the vertical dimension it goes 
        ! over the domain of the spectral bins 
        do node = 1, myDim_nod2D
            
            !___________________________________________________________________
            !loop over spectral bins 
!$ACC LOOP VECTOR     
            ! Compute flux at all interfaces (1 to nfbin-1)
            ! Bins 1 and nfbin are boundary bins (should have zero energy)
            do fbini = 1, nfbin
                ! wrap indices pre-computed in init (branch-free, module arrays)

                CFL     = abs(cs(fbini, node)) * dt / dphi(fbini)

                !_______________________________________________________________
                ! compute tracer difference
                ttf0    = ttf(fbini              , node)
                ttfp1   = ttf(iwe2_idxp1(fbini)  , node)
                dttf0p1 = ttfp1 - ttf0

                !_______________________________________________________________
                ! tracer Slope Ratio Calculation for upwind point
                ! compute tracer slope
                R      = (ttfp1-ttf(iwe2_idxp2(fbini),node))/(-dttf0p1+small)

                ! apply superbee limiter
                Cr     = slimiter_superbee(R)

                ! construct edge centered tracer value
                ! T_(i+0.5) = T_i+1 - [ (T_i+1 - T_i)/dz ]_Limited * dz/2
                !  --> this is seconds order in space, but first order in time
                !
                ! T_(i+0.5) = T_i+1 - [ (T_i+1 - T_i)/dz ]_Limited * (dz/2 -|W|*dt/2)
                ! T_(i+0.5) = T_i+1 - [ (T_i+1 - T_i)/dz ]_Limited * dz/2 *(1 - CFL)
                !  --> CFL = W*dt/dx
                Tmean2 = ttfp1 + 0.5_WP*Cr*(-dttf0p1) * (1.0_WP-merge(CFL, 0.0_WP, flag_2ndord_time))
                Tmean2 = merge(max(0.0_WP, Tmean2), Tmean2, flag_posdef)

                !_______________________________________________________________
                ! tracer Slope Ratio Calculation for downwind point
                ! compute tracer slope
                R      = (ttf0-ttf(iwe2_idxm1(fbini),node))/(dttf0p1+small)

                ! apply superbee limiter
                Cr     = slimiter_superbee(R)

                ! construct edge centered tracer value
                ! T_(i+0.5) = T_i + [ (T_i+1 - T_i)/dz ]_Limited * dz/2
                !  --> this is seconds order in space, but first order in time
                !
                ! T_(i+0.5) = T_i + [ (T_i+1 - T_i)/dz ]_Limited * (dz/2 -|W|*dt/2)
                ! T_(i+0.5) = T_i + [ (T_i+1 - T_i)/dz ]_Limited * dz/2 *(1 - CFL)
                !  --> CFL = W*dt/dx
                Tmean1 = ttf0 + 0.5_WP*Cr*(dttf0p1) * (1.0_WP-merge(CFL, 0.0_WP, flag_2ndord_time))
                Tmean1 = merge(max(0.0_WP, Tmean1), Tmean1, flag_posdef)

                !_______________________________________________________________
                ! cross spectral flux matching pyOM2 (NO negative sign, NO area)
                vflux = cs(fbini, node)
                vfabs = abs(vflux)
                flux(fbini, node) = -0.5_WP*(                     &
                                            (vflux+vfabs)*Tmean1+ &
                                            (vflux-vfabs)*Tmean2  &
                                            )
                
            end do !--> do fbini = 1, nfbin
!$ACC END LOOP
        end do ! --> do node = 1, myDim_nod2D
        
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
    !$ACC END PARALLEL LOOP
#endif
    end subroutine adv_Ecompart_crss_spctrl_superbee
    
    
    
    !
    !
    !___________________________________________________________________________
    ! superbee slope limiter
    pure elemental real(kind=WP) function slimiter_superbee(R) result(Cr)
        real(kind=WP), intent(in) :: R
        Cr = max( 0._WP, min(1.0_WP, 2.0_WP*R), min(2.0_WP, R) )
    end function slimiter_superbee
    
    
    
    !
    !
    !___________________________________________________________________________
    ! build tracer divergence from fluxes
    ! 
    ! dE/dt = -div(vec_u * E) + Forc    | *int(dV)
    ! E^(n+1)*V^(n+1) - E^n*V^n = dt*( -int_circ( vec_u^n*vec_n*E^n*dA) + Forc*V^n)
    ! |
    ! |--> V^(n+1) = V^(n)
    ! v
    ! E^(n+1) = E^(n) + dt/V^n * ( - Flx_hv + Forc*V^n)
    ! |--> Flx_hv = Flx_h + Flx_v
    subroutine adv_Ecompart_flx2tra_spctrl(&
                  flx_h                   &
                , flx_v                   &
                , refl_src                &
                , refl_node               &
                , div_h                   &
                , div_v                   &
                , dphit                   &
                , vol_s                   &
                , partit                  &
                , mesh                    &
                    )
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in   ), target :: mesh
        real(kind=WP) , intent(in)            :: flx_h(   idemix2_nfbin, partit%myDim_edge2D) 
        real(kind=WP) , intent(in)            :: flx_v(   idemix2_nfbin, partit%myDim_nod2D) 
        real(kind=WP) , intent(inout)         :: refl_src(idemix2_nfbin, partit%myDim_nod2D ) 
        logical       , intent(in)            :: refl_node(              partit%myDim_nod2D+partit%eDim_nod2D ) 
        
        real(kind=WP) , intent(inout)         :: div_h(   idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP) , intent(inout)         :: div_v(   idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP) , intent(in   )         :: dphit(   idemix2_nfbin)
        real(kind=WP) , intent(in   )         :: vol_s(partit%myDim_nod2D+partit%eDim_nod2D)
        !___LOCAL VARIABLES_____________________________________________________
        integer                               :: node, edge, fbini, nfbin, ednodes(2)
        real(kind=WP)                         :: inv_dphi(idemix2_nfbin), ivols_lcl
        
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        
        nfbin = idemix2_nfbin
        !
        !
        !_______________________________________________________________________
#ifndef ENABLE_OPENACC
        !$OMP PARALLEL DO COLLAPSE(2)
#else
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
        do node=1, myDim_nod2D+eDim_nod2D
            do fbini=1, nfbin
                div_h(fbini, node)=0.0_WP
                div_v(fbini, node)=0.0_WP
            end do
        end do
#ifndef ENABLE_OPENACC
        !$OMP END PARALLEL DO
#else
        !$ACC END PARALLEL LOOP
#endif 
        do fbini=1, nfbin
            inv_dphi(fbini) = 1 / dphit(fbini)
        end do
        !
        !
        !_______________________________________________________________________
        ! Vertical part
#ifndef ENABLE_OPENACC
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, edge, ednodes, fbini)
       !$OMP DO
#else
       !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
        do node=1, myDim_nod2d
            !$ACC LOOP VECTOR
            do fbini=2,nfbin-1
                ! Cross-spectral divergence: (flux(fbini) - flux(fbini-1))/dphi
                ! flux(fbini) is at interface between fbini and fbini+1
                ! flux(fbini-1) is at interface between fbini-1 and fbini
                div_v(fbini,node) = div_v(fbini,node) + (flx_v(fbini,node)-flx_v(fbini-1,node)) * inv_dphi(fbini)
            end do
            !$ACC END LOOP
        end do
#ifndef ENABLE_OPENACC
       !$OMP END DO
#else
       !$ACC END PARALLEL LOOP
#endif

        !
        !
        !_______________________________________________________________________
        ! Horizontal part
#ifndef ENABLE_OPENACC
       !$OMP DO
#else
#   if !defined(DISABLE_OPENACC_ATOMICS)
       !$ACC PARALLEL LOOP GANG PRIVATE(enodes, el) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#   else
       !$ACC UPDATE SELF(dttf_h, flux_h)
#   endif
#endif
        do edge=1, myDim_edge2D 
            ednodes  = edges(:,edge)
#ifndef ENABLE_OPENACC
#   if defined(_OPENMP) && !defined(__openmp_reproducible)
           call omp_set_lock(partit%plock(ednodes(1)))
#   else
           !$OMP ORDERED
#   endif 
#else
#   if !defined(DISABLE_OPENACC_ATOMICS)
           !$ACC LOOP VECTOR
#   endif
#endif      
            !
            !___________________________________________________________________
            ! horizontal flux contribution to ednode_1
            do fbini=2, nfbin-1
#if !defined(DISABLE_OPENACC_ATOMICS)
               !$ACC ATOMIC UPDATE
#endif
                ! Horizontal divergence: flux/area (NO inv_dphi - that's only for cross-spectral!)
                div_h(fbini, ednodes(1))=div_h(fbini, ednodes(1))+flx_h(fbini, edge)
                
                
#ifndef ENABLE_OPENACC
#   if defined(_OPENMP)  && !defined(__openmp_reproducible)
            ! if there is now __openmp or __openmp_reproducible assemble horizontal
            ! divergence in one loop
            end do
            
            !___________________________________________________________________
            call omp_unset_lock(partit%plock(ednodes(1)))
            call omp_set_lock  (partit%plock(ednodes(2)))
            
            !
            !___________________________________________________________________
            ! horizontal flux contribution to ednode_2
            do fbini=2, nfbin-1
#   endif    
#else
#   if !defined(DISABLE_OPENACC_ATOMICS)
                !$ACC ATOMIC UPDATE
#   endif                
#endif
                ! Horizontal divergence: flux/area (NO inv_dphi - that's only for cross-spectral!)
                div_h(fbini,ednodes(2))=div_h(fbini,ednodes(2))-flx_h(fbini,edge)
            end do
            
#ifndef ENABLE_OPENACC
#   if defined(_OPENMP)  && !defined(__openmp_reproducible)
           call omp_unset_lock(partit%plock(ednodes(2)))
#   else
            !$OMP END ORDERED
#   endif
#else
#   if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC END LOOP
#   endif
#endif
        end do

#ifndef ENABLE_OPENACC
   !$OMP END DO
   !$OMP END PARALLEL
#else
#   if !defined(DISABLE_OPENACC_ATOMICS)
   !$ACC END PARALLEL LOOP
#   else
   !$ACC UPDATE DEVICE(dttf_h)
#   endif
#endif
        !_______________________________________________________________________
        ! normalize with total volume of scalar cell from top to bottom
#ifndef ENABLE_OPENACC
        ! precompute reciprocal outside inner loop: FP division is ~6x more expensive
        ! than multiply on CPU, and the outer loop already saturates thread count
        !$OMP PARALLEL DO PRIVATE(ivols_lcl)
        do node=1, myDim_nod2d! +eDim_nod2D
            ivols_lcl = 1/vol_s(node)
            do fbini=2,nfbin-1
                div_h(fbini, node) = div_h(fbini,node)*ivols_lcl
            end do
        end do
        !$OMP END PARALLEL DO
#else
        ! on GPU, division throughput is far less asymmetric than on CPU and memory
        ! bandwidth dominates; inline division to keep loops perfectly nested for
        ! COLLAPSE(2), which gives a flat warp-filling iteration space
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
        do node=1, myDim_nod2d! +eDim_nod2D
            do fbini=2,nfbin-1
                div_h(fbini, node) = div_h(fbini,node)/vol_s(node)
            end do
        end do
        !$ACC END PARALLEL LOOP
#endif

        !_______________________________________________________________________
        ! inject blocked coast fluxes as source at interior nodes in mirror bin
#ifndef ENABLE_OPENACC
        !$OMP PARALLEL DO
        do node = 1, myDim_nod2D
            if (.not. refl_node(node)) cycle
            do fbini = 2, nfbin-1
                div_h(fbini, node) = div_h(fbini, node) + refl_src(fbini, node)
            end do
        end do
        !$OMP END PARALLEL DO
#else
        ! cycle on outer loop is semantically wrong in a COLLAPSE(2) combined iteration
        ! space (would cycle fbini, not node); condition moved inside inner loop —
        ! a predicated assignment is GPU-friendly and restores perfect nesting
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
        do node = 1, myDim_nod2D
            do fbini = 2, nfbin-1
                if (refl_node(node)) div_h(fbini, node) = div_h(fbini, node) + refl_src(fbini, node)
            end do
        end do
        !$ACC END PARALLEL LOOP
#endif

    end subroutine adv_Ecompart_flx2tra_spctrl
    
    
    
    !
    !
    !___________________________________________________________________________
    ! integrate equation for wave energy compartments
    subroutine hsintegrate_Ecompart(      &
                ti, tip1                  &
                , Ename                   &
                , E                       & 
                , Edivh                   &
                , Edivs                   &
                , vel                     &
                , cs                      & ! cross-spectral exchange rate phi_dot
                , forc                    &
                , tauE                    &
                , dphit                   &
                , gradxy_e                &
                , gradxy_n                &
                , flx_uv                  &
                , flx_cs                  &
                , refl_bin                &
                , refl_sgn                &
                , refl_src                &
                , refl_node               &
                , vol_s                   &
                , partit                  &
                , mesh                    &
                , flag_AB2                & ! to 2nd order Adams-Bashfort on top 
                , Edt                     &
                , Eadvh                   &
                , Eadvs                   &
                , Ediss                   &
                , Eforc                   &
                , Erefl                   &
                )
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit)  , intent(inout), target  :: partit
        type(t_mesh)    , intent(in   ), target  :: mesh
        integer         , intent(in   )          :: ti, tip1
        character(len=*), intent(in   )          :: Ename
        real(kind=WP)   , intent(inout)          :: E(    idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D, 2)
        real(kind=WP)   , intent(inout)          :: Edivh(idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D, 2)
        real(kind=WP)   , intent(inout)          :: Edivs(idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D, 2)
        real(kind=WP)   , intent(in   )          :: vel(     2, idemix2_nfbin, partit%myDim_elem2D+partit%eDim_elem2D)
        real(kind=WP)   , intent(inout)          :: cs(         idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(in   )          :: forc(       idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(in   )          :: tauE(                      partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(inout)          :: gradxy_e(2, idemix2_nfbin, partit%myDim_elem2D+partit%eDim_elem2D)
        real(kind=WP)   , intent(inout)          :: gradxy_n(2, idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(in   )          :: dphit(      idemix2_nfbin)
        real(kind=WP)   , intent(inout)          :: flx_uv(     idemix2_nfbin, partit%myDim_edge2D)
        real(kind=WP)   , intent(inout)          :: flx_cs(     idemix2_nfbin, partit%myDim_nod2D )
        integer         , intent(in   )          :: refl_bin(   partit%myDim_edge2D )
        integer         , intent(in   )          :: refl_sgn(   partit%myDim_edge2D ) 
        real(kind=WP)   , intent(inout)          :: refl_src(   idemix2_nfbin, partit%myDim_nod2D ) 
        logical         , intent(in   )          :: refl_node(                 partit%myDim_nod2D +partit%eDim_nod2D) 
        real(kind=WP)   , intent(in   )          :: vol_s(                     partit%myDim_nod2D +partit%eDim_nod2D)
        logical         , intent(in)             :: flag_AB2
        
        real(kind=WP)   , intent(inout), optional:: Edt(        idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP)   , intent(inout), optional:: Eadvh(      idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP)   , intent(inout), optional:: Eadvs(      idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP)   , intent(inout), optional:: Ediss(      idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP)   , intent(inout), optional:: Eforc(      idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP)   , intent(inout), optional:: Erefl(      idemix2_nfbin, partit%myDim_nod2D)
        
        !___LOCAL VARIABLES_____________________________________________________
        integer                                  :: node, fbini, nfbin
        real(kind=WP)                            :: inv_denom
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"        
        !_______________________________________________________________________
        nfbin=idemix2_nfbin
        iwe2_flx_uv(:,:)     = 0.0
        iwe2_flx_w( :,:)     = 0.0
        
        !_______________________________________________________________________
        ! Exchange current timelevel (contiguous with timelevel-last layout)
        call exchange_nod_fbin(E(:, :, ti), partit)

        !_______________________________________________________________________
        ! compute gradient of iwe2_E on nodes
        call tracer_gradient_elements(E(:,:,ti), gradxy_e, partit, mesh, .False.)
        call exchange_elem_fbin(gradxy_e, partit)
        call interp_e2n(gradxy_e, gradxy_n, mesh, partit, .False.)
        call exchange_nod_fbin(gradxy_n, partit)

        !_______________________________________________________________________
        ! compute horizontal superbee advected tracer flux
        call adv_Ecompart_hor_spctrl_superbee(            &
                                          vel             &
                                        , E(:,:,ti)       &
                                        , gradxy_n        &
                                        , flx_uv          &
                                        , partit          &
                                        , mesh            &
                                        , .True.          & ! do 2nd order in space & time
                                        , .True.          & ! enforce positive definit
                                        )

        ! apply reflective coastal BC: redirect coast-directed fluxes into mirror bins
        call apply_reflect_bc_spctrl(     flx_uv          &
                                        , refl_src        &
                                        , refl_bin        &
                                        , refl_sgn        &
                                        , vol_s           &
                                        , partit          &
                                        , mesh            &
                                        )

        !_______________________________________________________________________
        ! compute cross spectral superbee advected tracer flux
        call adv_Ecompart_crss_spctrl_superbee(           &
                                          cs              &
                                        , E(:,:,ti)       &
                                        , dphit           &
                                        , flx_cs          &
                                        , partit          &
                                        , mesh            &
                                        , .True.          & ! do 2nd order in space & time
                                        , .True.          & ! enforce positive definit
                                        )

        !_______________________________________________________________________
        ! compute flux divergence
        call adv_Ecompart_flx2tra_spctrl(                 &
                                         flx_uv           &
                                        ,flx_cs           &
                                        ,refl_src         &
                                        ,refl_node        &
                                        ,Edivh(:,:,ti)    &
                                        ,Edivs(:,:,ti)    &
                                        ,dphit            &
                                        ,vol_s            &
                                        ,partit           &
                                        ,mesh             &
                                        )

        !_______________________________________________________________________
        ! Eden and Olbers 2014, eq. 2
        ! integrate E_M2^(n+1) = E_M2^(n) + dt*[ 
        !                        - div(c_uv*E_M2^n)         --> (the minus is hidden in the flux computation)
        !                        - d/dphi (dphit/dt*E_M2^n) --> (the minus is hidden in the flux computation)
        !                        + forc
        !                        + W + T]
        !                      
        ! eq. 7 & 8: Wave-Wave (T) and Topographic interaction
        ! T_M2  + W_M2  =  alpha_c_M2*E_M2*int(Eiw*dz) - 1/tau_M2 * E_M2
        !                  └───────────┬─────────────┘
        !                              └> this part added in wave-wave-interaction
        ! T_niw + W_niw =  -1/tau_niw * E_niw
        ! --> here variable tauE is already 1/tau
        
        !_______________________________________________________________________
        ! Adams-Bashforth 2nd order + implicit dissipation
        if (flag_AB2) then
            do node= 1, myDim_nod2d
                inv_denom = 1.0_WP / (1.0_WP + dt*tauE(node))
                do fbini=2,nfbin-1
                    E(fbini, node, tip1) = ( E(fbini, node, ti) &
                                            + dt * ( &
                                                    + (1.5_WP+idemix2_AB_epsilon)*Edivh(fbini, node, ti  ) &
                                                    - (0.5_WP+idemix2_AB_epsilon)*Edivh(fbini, node, tip1) &
                                                    + (1.5_WP+idemix2_AB_epsilon)*Edivs(fbini, node, ti  ) &
                                                    - (0.5_WP+idemix2_AB_epsilon)*Edivs(fbini, node, tip1) &
                                                    +                              forc(fbini, node)       &
                                                ) &
                                            ) * inv_denom
                end do
                E(    1, node, tip1) = E(nfbin-1, node, tip1)
                E(nfbin, node, tip1) = E(      2, node, tip1)
            end do
        !_______________________________________________________________________
        ! Forward Euler + implicit dissipation
        else
            do node= 1, myDim_nod2d
                inv_denom = 1.0_WP / (1.0_WP + dt*tauE(node))
                do fbini=2,nfbin-1
                    E(fbini, node, tip1) = ( E(fbini, node, ti) &
                                            + dt * ( &
                                                    + Edivh(fbini, node, ti) & ! div(c_uv*E^n)
                                                    + Edivs(fbini, node, ti) & ! d/dphi (w*E^n)
                                                    + forc( fbini, node)     & ! forc
                                                   ) &
                                            ) * inv_denom
                end do
                E(    1, node, tip1) = E(nfbin-1, node, tip1)
                E(nfbin, node, tip1) = E(      2, node, tip1)
            end do
        end if

        !_______________________________________________________________________
        ! additional diagnostics
        if (present(Edt))   Edt(:,:)   =  (E(:, 1:myDim_nod2D, tip1)-E(:, 1:myDim_nod2D, ti)) / dt

        if (present(Eadvh)) then
            if (flag_AB2) then
                Eadvh(:,:) =  ((1.5+idemix2_AB_epsilon)*Edivh(:, 1:myDim_nod2D, ti  ) &
                            - (0.5+idemix2_AB_epsilon)*Edivh(:, 1:myDim_nod2D, tip1))
            else
                Eadvh(:,:) = Edivh(:, 1:myDim_nod2D, ti)
            end if
        end if

        if (present(Eadvs)) then
            if (flag_AB2) then
                Eadvs(:,:) =  ((1.5+idemix2_AB_epsilon)*Edivs(:, 1:myDim_nod2D, ti  ) &
                            - (0.5+idemix2_AB_epsilon)*Edivs(:, 1:myDim_nod2D, tip1))
            else
                Eadvs(:,:) = Edivs(:, 1:myDim_nod2D, ti)
            end if
        end if

        if (present(Ediss)) then
            do node= 1, myDim_nod2d
                Ediss(:,node) =  - tauE(node)*E(:, node, tip1)
            end do
        end if

        if (present(Eforc)) Eforc(:,:) = forc(:, 1:myDim_nod2D)

        if (present(Erefl)) Erefl(:,:) = refl_src(:, 1:myDim_nod2D)

    end subroutine hsintegrate_Ecompart
    
    
    !
    !
    !___________________________________________________________________________
    ! compute inverse volume of w scalar cell 
    ! ~~~~~~~~~o zbar_1~~~~~~~~~~
    !          |
    !          x Z_1 -------->^
    !          |              |
    !          o zbar_2, w_2  | h_wcell * areasvol = Vol_wcell --> Vol_wcelli = 1/Vol_wcell
    !          |              | 
    !          x Z_2 -------->v
    !          |
    !          o zbar_3
    !          :
    !          :
    !          :
    ! ---------o zbar_n----------
    ! ///////////////////////////
    subroutine compute_vol_wcell(vol_wcelli, mesh, partit)
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit)  , intent(inout), target :: partit
        type(t_mesh)    , intent(in   ), target :: mesh
        real(kind=WP)   , intent(inout)         :: vol_wcelli(:,:)
        !___LOCAL VARIABLES_____________________________________________________
        integer                                 :: node, nz, nln, uln
        real(kind=WP)                           :: dz_trr(mesh%nl)
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"   
        do node = 1, myDim_nod2D + eDim_nod2D
            ! number of above bottom levels at node
            nln = nlevels_nod2D(node)-1
            uln = ulevels_nod2D(node)
                
            ! thickness of mid-level to mid-level interface at node
            dz_trr             = 0.0_WP
            dz_trr(uln+1:nln)  = Z_3d_n(uln:nln-1,node)-Z_3d_n(uln+1:nln,node)
            dz_trr(uln)        = hnode(uln,node)/2.0_WP
            dz_trr(nln+1)      = hnode(nln,node)/2.0_WP
                
            ! surface cell 
            vol_wcelli(uln,node)   = 1/areasvol(uln,node)/dz_trr(uln)
            do nz=uln+1,nln
                ! inverse volumne centered around full depth levels
                vol_wcelli(nz,node)= 1/areasvol(nz-1,node)/dz_trr(nz)
            end do 
            ! bottom cell 
            vol_wcelli(nln+1,node) = 1/areasvol(nln,node)/dz_trr(nln+1)
            
        end do !-->do node = 1,node_size
    end subroutine compute_vol_wcell
    
    
    
    !
    !
    !___________________________________________________________________________
    ! compute volume of scalar cell from bottom to top
    ! ~~~~~~~~~o zbar_1~~~~~~~~~~
    !          |
    !          x Z_1 
    !          |             
    !          o zbar_2 
    !          |     
    !          x Z_2 
    !          |
    !          o zbar_3
    !          :
    !          :
    !          :
    ! ---------o zbar_n----------
    ! ///////////////////////////
    subroutine compute_vol_nodB2T_fix(vol_nodB2T, mesh, partit)
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit)  , intent(inout), target :: partit
        type(t_mesh)    , intent(in   ), target :: mesh
        real(kind=WP)   , intent(inout)         :: vol_nodB2T(partit%myDim_nod2D+partit%eDim_nod2D)
        !___LOCAL VARIABLES_____________________________________________________
        integer                                 :: node, elem, nz, nln, uln, elnodes(3)
!         real(kind=WP)                           :: lcl_sumvol, glb_sumvol, vol
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"   
        vol_nodB2T = 0.0_WP
!         lcl_sumvol = 0.0_WP
        do node = 1, myDim_nod2D
            nln = nlevels_nod2D(node)-1
            uln = ulevels_nod2D(node)
            do nz=uln,nln
!                 lcl_sumvol       = lcl_sumvol       + areasvol(nz,node)*abs(zbar(nz)-zbar(nz+1))
                vol_nodB2T(node) = vol_nodB2T(node) + areasvol(nz,node)*abs(zbar(nz)-zbar(nz+1))
            end do
        end do !-->do node = 1,node_size
        ! call MPI_Allreduce(lcl_sumvol, glb_sumvol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        ! if (partit%mype == 0) then
        !     write(*,*) ' debug: vol_node_fix: ', glb_sumvol
        ! end if
        call exchange_nod(vol_nodB2T, partit)
        
        ! lcl_sumvol = 0.0_WP
        ! do elem = 1, myDim_elem2D
        !     ! Nodes are uniquely partitioned (e.g., via METIS) — each node belongs 
        !     ! to exactly one rank Elements are then distributed so that each rank 
        !     ! has all elements touching its owned nodes — but boundary elements 
        !     ! necessarily appear on multiple ranks (they share vertices across partitions)
        !     elnodes=mesh%elem2D_nodes(:,elem)
        !     if (elnodes(1) > myDim_nod2D) cycle
        !     nln = nlevels(elem)-1
        !     uln = ulevels(elem)
        !     do nz=uln,nln
        !         lcl_sumvol       = lcl_sumvol       + elem_area(elem)*abs(zbar(nz)-zbar(nz+1))
        !     end do
        ! end do !-->do node = 1,node_size
        ! call MPI_Allreduce(lcl_sumvol, glb_sumvol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        ! if (partit%mype == 0) then
        !     write(*,*) ' debug: vol_elem_fix: ', glb_sumvol
        ! end if
        
    end subroutine compute_vol_nodB2T_fix
    
    
    
    !
    !
    !___________________________________________________________________________
    ! compute volume of scalar cell from bottom to top
    ! ~~~~~~~~~o zbar_1~~~~~~~~~~
    !          |
    !          x Z_1 
    !          |             
    !          o zbar_2 -------+
    !          |               |
    !          x Z_2           |--------> hnode(2, node)
    !          |               |
    !          o zbar_3--------+
    !          :
    !          :
    !          :
    ! ---------o zbar_n----------
    ! ///////////////////////////
    subroutine compute_vol_nodB2T_hnode(vol_nodB2T, mesh, partit)
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit)  , intent(inout), target :: partit
        type(t_mesh)    , intent(in   ), target :: mesh
        real(kind=WP)   , intent(inout)         :: vol_nodB2T(partit%myDim_nod2D+partit%eDim_nod2D)
        !___LOCAL VARIABLES_____________________________________________________
        integer                                 :: node, elem, nz, nln, uln
!         real(kind=WP)                           :: lcl_sumvol, glb_sumvol, vol
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"   
        vol_nodB2T = 0.0_WP
        ! lcl_sumvol = 0.0_WP
        do node = 1, myDim_nod2D
            nln = nlevels_nod2D(node)-1
            uln = ulevels_nod2D(node)
            do nz=uln,nln
                ! lcl_sumvol       = lcl_sumvol       + areasvol(nz,node)*hnode(nz,node)
                vol_nodB2T(node) = vol_nodB2T(node) + areasvol(nz,node)*hnode(nz,node)
            end do
        end do !-->do node = 1,node_size
        ! call MPI_Allreduce(lcl_sumvol, glb_sumvol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        ! if (partit%mype == 0) then
        !     write(*,*) ' debug: vol_node_hnode: ', glb_sumvol
        ! end if
        call exchange_nod(vol_nodB2T, partit)
        
    end subroutine compute_vol_nodB2T_hnode
    
    
    
    !
    !
    !___________________________________________________________________________
    ! --> add contribution from horizontal wave propagation
    ! Since IDEMIX is used in a global model configuration the vertical 
    ! internal wave mixing (call cvmix_coeffs_idemix) have to be extended by 
    ! a lateral diffusion term (see. Olbers D., Eden C., 2013, A Global Model 
    ! for the Diapycnal Diffusivity Induced Internal Gravity Waves...)
    !
    ! diffusion term = div_h( v_0 * idemix_tau_h * grad_h(v_0*E_iw) )
    ! 
    ! use Gaussian integral satz ... int(div vec_A)dV = ringint(A*vec_n)dA
    !                                    div vec_A    = 1/V * sum_i=1...nface( A_i*vec_n_i)*A_i
    subroutine compute_hdiff_Eiw(         &
                Eiw_old                 , &
                Eiw                     , &
                v0                      , &
                n_hor_iter              , &
                partit                  , &
                mesh                      &
                )
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit)  , intent(inout), target :: partit
        type(t_mesh)    , intent(in   ), target :: mesh
        real(kind=WP)   , intent(in   )         :: Eiw_old(:,:)
        real(kind=WP)   , intent(inout)         :: Eiw(:,:)
        real(kind=WP)   , intent(in   )         :: v0(:,:)
        integer         , intent(in   )         :: n_hor_iter
        !___LOCAL VARIABLES_____________________________________________________
        integer                                 :: edge, ednodes(2), edel(2),          &
                                                   elnodes1(3), elnodes2(3),           & 
                                                   node, nz,                           &
                                                   nl1, nl2, nl12,                     &
                                                   nu1, nu2, nu12
        real(kind=WP)                           :: dz_trr1(mesh%nl), dz_trr2(mesh%nl), &
                                                   grad_v0Eiw(2)
        real(kind=WP)                           :: dx1, dx2, dy1, dy2, vflux
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h" 
        !_______________________________________________________________________
        ! calculate horizontal diffusion term div( grad(Eiw*v0) * tau_h*v0) 
        do edge=1,myDim_edge2D
            !___________________________________________________________________
            ! ednode ... vertices that form the edge
            ednodes  = edges(:,edge)
            
            ! deltaX1,deltaY1 ... cross edge dx and dy for the two elements 
            ! that contribute to the edge
            !                 o
            !                / \
            !               /   \       + ... triangle centroids
            !              /  +  \      x ... edge mid points
            !             /   ^   \
            !            /    |----\--------------->edge_cross_dxdy(1:2,edge)
            ! ednodes(1)o-----x-----o ednodes(2)
            !            \    |----/--------------->edge_cross_dxdy(3:4,edge)
            !             \   v   /
            !              \  +  /     
            !               \   /
            !                \ /
            !                 o
            !                 
            dx1      = edge_cross_dxdy(1,edge)
            dy1      = edge_cross_dxdy(2,edge)
            ! el ... elements that contribute to edge
            edel     = edge_tri(:,edge)
            ! elnodes1 ... nodes that contribute to element el(1) 
            elnodes1 = elem2d_nodes(:,edel(1))
            ! nl1 ... number of layers at element el(1)
            nl1      = nlevels(edel(1))
            ! nu1 ... upper index of ocean default = 1 but can consider cavity !=1
            nu1      = ulevels(edel(1))
                    
            ! thickness of mid-level to mid-level interface of element el(1)
            dz_trr1         = 0.0_WP
            dz_trr1(nu1)    = helem(nu1  ,edel(1))/2.0_WP
            do nz = nu1+1, nl1-1
                dz_trr1(nz) = helem(nz-1 ,edel(1))/2.0_WP + helem(nz,edel(1))/2.0_WP
            end do
            dz_trr1(nl1)    = helem(nl1-1,edel(1))/2.0_WP
                    
            !___________________________________________________________________
            ! the same as above but for el(2)--> if el(2)==0 than this edge 
            ! is a boundary edge and el(2) does not exist
            nl2=0
            nu2=0
            if (edel(2)>0) then 
                dx2      = edge_cross_dxdy(3,edge)
                dy2      = edge_cross_dxdy(4,edge)
                elnodes2 = elem2d_nodes(:,edel(2))
                nl2      = nlevels(edel(2))
                nu2      = ulevels(edel(2))
                        
                ! thickness of mid-level to mid-level interface of element el(2)
                dz_trr2         = 0.0_WP
                dz_trr2(nu2)    = helem(nu2  ,edel(2))/2.0_WP
                do nz = nu2+1, nl2-1
                    dz_trr2(nz) = helem(nz-1 ,edel(2))/2.0_WP + helem(nz,edel(2))/2.0_WP
                end do
                dz_trr2(nl2)    = helem(nl2-1,edel(2))/2.0_WP
            endif
            
            !___________________________________________________________________
            ! common depth levels between elements edel(1) and edel(2)
            nl12=min(nl1,nl2)
            nu12=max(nu1,nu2)
            
            !___________________________________________________________________
            ! (A) goes only into this loop when the edge has only facing element
            ! el(1) --> so the edge is a boundary edge --> this is for ocean 
            ! surface in case of cavity
            do nz=nu1,nu12-1
                !_______________________________________________________________
                ! --> calc: grad_h(v_0*E_iw)
                ! calculate flux from el(1) with respect to edge mid 
                ! point
                grad_v0Eiw(1) = sum(gradient_sca(1:3,edel(1))*v0(nz,elnodes1)*Eiw_old(nz,elnodes1))
                grad_v0Eiw(2) = sum(gradient_sca(4:6,edel(1))*v0(nz,elnodes1)*Eiw_old(nz,elnodes1))
                
                ! calculate flux --> grad_v0Eiw is now located on elements
                vflux = (grad_v0Eiw(1)*dy1-grad_v0Eiw(2)*dx1)*dz_trr1(nz)
                
                !_______________________________________________________________
                ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                ! multiply vflux with iwe_v0 interpolate to the edge-
                ! mid point 
                vflux = vflux * (v0(nz,ednodes(1))+v0(nz,ednodes(2)))*0.5_WP
                
                !_______________________________________________________________
                ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                ! sum fluxes over the surface --> gaussian integral satz
                !  --> keep in mind vol_wcelli is the inverse volume !!!
                Eiw(nz,ednodes(1)) = Eiw(nz,ednodes(1)) + dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(1))*vflux
                Eiw(nz,ednodes(2)) = Eiw(nz,ednodes(2)) - dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(2))*vflux
            end do !-->do nz=nu1,nu12-1
            
            !___________________________________________________________________
            ! (B) goes only into this loop when the edge has only facing elemenmt
            ! el(2) --> so the edge is a boundary edge --> this is for ocean 
            ! surface in case of cavity
            if (nu2 > 0) then 
                do nz=nu2,nu12-1
                    !___________________________________________________________
                    ! --> calc: grad_h(v_0*E_iw)
                    ! first calculate flux from el(1) with respect to edge mid 
                    ! point
                    grad_v0Eiw(1) = sum(gradient_sca(1:3,edel(2))*v0(nz,elnodes2)*Eiw_old(nz,elnodes2))
                    grad_v0Eiw(2) = sum(gradient_sca(4:6,edel(2))*v0(nz,elnodes2)*Eiw_old(nz,elnodes2))
                    
                    ! calculate flux 
                    vflux = -(grad_v0Eiw(1)*dy2-grad_v0Eiw(2)*dx2)*dz_trr2(nz)
                    !       |--> minus sign comes from the fact that the the 
                    !            normal vectors (dx1,dy1) and (dx2,dy2) face 
                    !            in opposite direction (Right-Hand-Rule)
                    
                    !___________________________________________________________
                    ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                    ! multiply vflux with iwe_v0 interpolate to the edge-
                    ! mid point 
                    vflux = vflux * (v0(nz,ednodes(1))+v0(nz,ednodes(2)))*0.5_WP
                    
                    !___________________________________________________________
                    ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                    ! sum fluxes over the surface --> gaussian integral satz
                    Eiw(nz,ednodes(1)) = Eiw(nz,ednodes(1)) + dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(1))*vflux
                    Eiw(nz,ednodes(2)) = Eiw(nz,ednodes(2)) - dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(2))*vflux
                    
                end do !-->do nz=nu2,nu12-1
            end if 
            
            !___________________________________________________________________
            ! (C) goes only into this loop when the edge has two facing elements
            ! --> so the edge is not a boundary edge
            do nz=nu12,nl12
                !_______________________________________________________________
                ! --> calc: grad_h(v_0*E_iw)
                ! calculate grad(iwe*iwe_v0) for el(1)
                grad_v0Eiw(1) = sum(gradient_sca(1:3,edel(1))*v0(nz,elnodes1)*Eiw_old(nz,elnodes1))
                grad_v0Eiw(2) = sum(gradient_sca(4:6,edel(1))*v0(nz,elnodes1)*Eiw_old(nz,elnodes1))
                
                ! calculate flux --> grad_v0Eiw is now located on elements
                vflux = (grad_v0Eiw(1)*dy1-grad_v0Eiw(2)*dx1)*dz_trr1(nz)
                
                ! calculate grad(iwe*iwe_v0) for el(2) and average for el(1)
                ! and el(2)
                grad_v0Eiw(1) = sum(gradient_sca(1:3,edel(2))*v0(nz,elnodes2)*Eiw_old(nz,elnodes2))
                grad_v0Eiw(2) = sum(gradient_sca(4:6,edel(2))*v0(nz,elnodes2)*Eiw_old(nz,elnodes2))
                
                ! calculate flux 
                vflux = vflux + -(grad_v0Eiw(1)*dy2-grad_v0Eiw(2)*dx2)*dz_trr2(nz)
                    
                !_______________________________________________________________
                ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                ! multiply vflux with iwe_v0 interpolate to the edge-
                ! mid point 
                vflux = vflux * (v0(nz,ednodes(1))+v0(nz,ednodes(2)))*0.5_WP
                
                !_______________________________________________________________
                ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                ! sum fluxes over the surface --> gaussian integral satz
                Eiw(nz,ednodes(1)) = Eiw(nz,ednodes(1)) + dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(1))*vflux
                Eiw(nz,ednodes(2)) = Eiw(nz,ednodes(2)) - dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(2))*vflux
            end do !-->do nz=1,n2
            
            !___________________________________________________________________
            ! (D) goes only into this loop when the edge has only facing element
            ! el(1) --> so the edge is a boundary edge
            do nz=nl12+1,nl1
                !_______________________________________________________________
                ! --> calc: grad_h(v_0*E_iw)
                ! calculate flux from el(1) with respect to edge mid 
                ! point
                grad_v0Eiw(1) = sum(gradient_sca(1:3,edel(1))*v0(nz,elnodes1)*Eiw_old(nz,elnodes1))
                grad_v0Eiw(2) = sum(gradient_sca(4:6,edel(1))*v0(nz,elnodes1)*Eiw_old(nz,elnodes1))
                
                ! calculate flux 
                vflux = (grad_v0Eiw(1)*dy1-grad_v0Eiw(2)*dx1)*dz_trr1(nz)
                
                !_______________________________________________________________
                ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                ! multiply vflux with iwe_v0 interpolate to the edge-
                ! mid point 
                vflux = vflux * (v0(nz,ednodes(1))+v0(nz,ednodes(2)))*0.5_WP
                
                !_______________________________________________________________
                ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                ! sum fluxes over the surface --> gaussian integral satz
                Eiw(nz,ednodes(1)) = Eiw(nz,ednodes(1)) + dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(1))*vflux
                Eiw(nz,ednodes(2)) = Eiw(nz,ednodes(2)) - dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(2))*vflux
            end do !-->do nz=nl12+1,nl1
            
            !___________________________________________________________________
            ! (E) goes only into this loop when the edge has only facing elemenmt
            ! el(2) --> so the edge is a boundary edge
            do nz=nl12+1,nl2
                !_______________________________________________________________
                ! --> calc: grad_h(v_0*E_iw)
                ! first calculate flux from el(1) with respect to edge mid 
                ! point
                grad_v0Eiw(1) = sum(gradient_sca(1:3,edel(2))*v0(nz,elnodes2)*Eiw_old(nz,elnodes2))
                grad_v0Eiw(2) = sum(gradient_sca(4:6,edel(2))*v0(nz,elnodes2)*Eiw_old(nz,elnodes2))
                
                ! calculate flux 
                vflux = -(grad_v0Eiw(1)*dy2-grad_v0Eiw(2)*dx2)*dz_trr2(nz)
                
                !_______________________________________________________________
                ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                ! multiply vflux with iwe_v0 interpolate to the edge-
                ! mid point 
                vflux = vflux * (v0(nz,ednodes(1))+v0(nz,ednodes(2)))*0.5_WP
                
                !_______________________________________________________________
                ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                ! sum fluxes over the surface --> gaussian integral satz
                Eiw(nz,ednodes(1)) = Eiw(nz,ednodes(1)) + dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(1))*vflux
                Eiw(nz,ednodes(2)) = Eiw(nz,ednodes(2)) - dt*idemix2_tau_h/n_hor_iter*vol_wcelli(nz,ednodes(2))*vflux
                
            end do !-->do nz=nl12+1,nl1
        end do !-->do edge=1,myDim_edge2D
        
    end subroutine compute_hdiff_Eiw



    !
    !
    !___________________________________________________________________________
    ! Check if velocity field is divergence-free using gradient_vec operator
    subroutine check_Ecompart_div_uv(vel_uv, div_uv, partit, mesh, compartment_name)
        implicit none
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        real(kind=WP),  intent(in)            :: vel_uv(2, idemix2_nfbin, partit%myDim_elem2D+partit%eDim_elem2D)
        real(kind=WP),  intent(inout)         :: div_uv(idemix2_nfbin, partit%myDim_elem2D+partit%eDim_elem2D)
        character(len=*), intent(in)          :: compartment_name
        
        integer                               :: elem, fbini, j, neighbor_elem
        real(kind=WP)                         :: local_div_sum, global_div_sum, global_div_max
        real(kind=WP)                         :: u_neighbor(idemix2_nfbin), v_neighbor(idemix2_nfbin)
        
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        
        ! Compute divergence of velocity field using gradient_vec operator
        ! gradient_vec(1:3, elem) are coefficients for ∂/∂x
        ! gradient_vec(4:6, elem) are coefficients for ∂/∂y
        do elem = 1, myDim_elem2D
        
            ! Initialize divergence arrays
            div_uv(:, elem) = 0.0_WP
            
            do j = 1, 3
                neighbor_elem = mesh%elem_neighbors(j, elem)
                if (neighbor_elem > 0) then
                    u_neighbor(:) = vel_uv(1, :, neighbor_elem)
                    v_neighbor(:) = vel_uv(2, :, neighbor_elem)
                else
                    u_neighbor(:) = vel_uv(1, :, elem)
                    v_neighbor(:) = vel_uv(2, :, elem)
                end if 
                
                do fbini = 2, idemix2_nfbin-1
                    ! Compute Total divergence
                    ! Compute gradients using gradient_vec coefficients
                    ! ∂u/∂x = Σ(gradient_vec(1:3, elem) * u_neighbor)
                    ! ∂v/∂y = Σ(gradient_vec(4:6, elem) * v_neighbor)
                    div_uv(fbini, elem) = div_uv(fbini, elem) &
                                          + (mesh%gradient_vec(j  , elem)*u_neighbor(fbini)) &
                                          + (mesh%gradient_vec(j+3, elem)*v_neighbor(fbini))
                end do
            end do
        end do
        
    end subroutine check_Ecompart_div_uv

    
    
    !
    !
    !___________________________________________________________________________
    ! Check flux conservation: sum of all divergences weighted by area should be zero
    subroutine check_flux_conservation(div_h, div_s, ti, dphit, vol_s, partit, mesh, compartment_name)
        implicit none
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        real(kind=WP),  intent(in)            :: div_h(idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D, 3)
        real(kind=WP),  intent(in)            :: div_s(idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D, 3)
        integer,        intent(in)            :: ti
        real(kind=WP),  intent(in)            :: dphit(idemix2_nfbin)
        real(kind=WP),  intent(in)            :: vol_s(partit%myDim_nod2D+partit%eDim_nod2D)
        character(len=*), intent(in)          :: compartment_name
        
        integer                               :: node, fbini, lcl_nedge, glb_nedge
        real(kind=WP)                         :: lcl_divh_sum, lcl_divs_sum
        real(kind=WP)                         :: glb_divh_sum, glb_divs_sum
        
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        
        ! Compute local sum: div * area * dphi
        ! For a closed domain with no sources/sinks, this should sum to zero
        lcl_divh_sum = 0.0_WP
        do node = 1, myDim_nod2D + eDim_nod2D  ! Full extended domain: telescoping holds locally per rank
            do fbini = 2, idemix2_nfbin-1  ! Exclude ghost bins
                lcl_divh_sum = lcl_divh_sum + div_h(fbini, node, ti)*vol_s(node)
            end do
        end do
!         write(*,*) lcl_divh_sum, mype

        lcl_divs_sum = 0.0_WP
        do node = 1, myDim_nod2D  ! Full extended domain: telescoping holds locally per rank
            do fbini = 2, idemix2_nfbin-1  ! Exclude ghost bins
                lcl_divs_sum = lcl_divs_sum + div_s(fbini, node, ti)*dphit(fbini)
            end do
        end do
        
        ! Sum across all MPI ranks
        glb_divh_sum = 0.0_WP
        glb_divs_sum = 0.0_WP
        call MPI_Allreduce(lcl_divh_sum, glb_divh_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        call MPI_Allreduce(lcl_divs_sum, glb_divs_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        
        if (mype == 0) then
            write(*,*) ' debug: flux total conservation check (',trim(compartment_name),'):',&
                       '  Σ(div_h*vols) = ', glb_divh_sum, &
                       ', Σ(div_s*dphi) = ', glb_divs_sum
        end if
        
    end subroutine check_flux_conservation



    !
    !
    !___________________________________________________________________________
    ! Compute and print global total NIW energy for conservation check
    subroutine check_global_energy(Ecompart, ti, dphit, vol_s, forc, tauE, partit, mesh, compartment_name)
        implicit none
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        real(kind=WP),  intent(in)            :: Ecompart(idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D, 3)
        real(kind=WP),  intent(in)            :: vol_s(partit%myDim_nod2D+partit%eDim_nod2D)
        integer,        intent(in)            :: ti
        real(kind=WP),  intent(in)            :: dphit(idemix2_nfbin)
        real(kind=WP),  intent(in)            :: forc(idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP),  intent(in)            :: tauE(partit%myDim_nod2D+partit%eDim_nod2D)
        character(len=*), intent(in)          :: compartment_name
        
        integer                               :: node, fbini
        real(kind=WP)                         :: lcl_energy, lcl_vol, glb_energy, glb_vol
        real(kind=WP)                         :: lcl_forc, glb_forc, lcl_diss, glb_diss
        real(kind=WP)                         :: lcl_maxE, glb_maxE, lcl_minE, glb_minE
        
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        
        lcl_energy = 0.0_WP
        lcl_vol    = 0.0_WP
        lcl_forc   = 0.0_WP
        lcl_diss   = 0.0_WP
        lcl_maxE   = -1.0e33
        lcl_minE   =  1.0e33
        do node = 1, myDim_nod2D  ! Only owned nodes
            lcl_vol = lcl_vol + vol_s(node)
            do fbini = 2, idemix2_nfbin-1  ! Exclude ghost bins
                lcl_energy = lcl_energy + Ecompart(fbini, node, ti) * dphit(fbini) * vol_s(node)
                lcl_forc   = lcl_forc   + forc(fbini, node)         * dphit(fbini) * vol_s(node)
                lcl_diss   = lcl_diss   + tauE(node) * Ecompart(fbini, node, ti) * dphit(fbini) * vol_s(node)
                lcl_maxE   = max(lcl_maxE, Ecompart(fbini, node, ti))
                lcl_minE   = min(lcl_minE, Ecompart(fbini, node, ti))
            end do
        end do
        
        ! Sum across all MPI ranks
        glb_energy = 0.0_WP
        glb_vol    = 0.0_WP
        glb_forc   = 0.0_WP
        glb_diss   = 0.0_WP
        call MPI_Allreduce(lcl_energy, glb_energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        call MPI_Allreduce(lcl_vol   , glb_vol   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        call MPI_Allreduce(lcl_forc  , glb_forc  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        call MPI_Allreduce(lcl_diss  , glb_diss  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        call MPI_Allreduce(lcl_minE  , glb_minE  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        call MPI_Allreduce(lcl_maxE  , glb_maxE  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        
        
        ! Print on rank 0
        if (mype == 0) then
            write(*,*) ' debug: global Σ(E_',trim(compartment_name),'*dphi*vols) = ', glb_energy, ' J', &
                             ', global min(E_',trim(compartment_name),') = ', glb_minE, &
                             ', global max(E_',trim(compartment_name),') = ', glb_maxE
            write(*,*) ' debug: global Σ(forc_',trim(compartment_name),'*dphi*vols) = ', glb_forc, ' W', &
                             ', global Σ(tau*E_',trim(compartment_name),'*dphi*vols) = ', glb_diss, ' W', &
                             ', global forc_',trim(compartment_name),'/diss_',trim(compartment_name),' = ', glb_forc/glb_diss, ' W'
        end if
        
    end subroutine check_global_energy
    
    
    !
    !___________________________________________________________________________
    ! Compute and print maximum CFL numbers for horizontal and cross-spectral advection
    subroutine check_max_cfl(vel, cs, dphit, partit, mesh, compartment_name)
        implicit none
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        real(kind=WP),  intent(in)            :: vel(2, idemix2_nfbin, partit%myDim_elem2D+partit%eDim_elem2D)
        real(kind=WP),  intent(in)            :: cs(idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP),  intent(in)            :: dphit(idemix2_nfbin)
        character(len=*), intent(in)          :: compartment_name
        
        integer                               :: edge, fbini, nfbin, node, el(2), el2, ednodes(2)
        real(kind=WP)                         :: dx1, dy1, dx2, dy2, dxdy12(2), edlen, n_x, n_y, n_len
        real(kind=WP)                         :: u1, u2, v1, v2, Ue, CFL_h, CFL_cs
        real(kind=WP)                         :: lcl_cfl_h, glb_cfl_h, lcl_cfl_cs, glb_cfl_cs
        
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        
        nfbin = idemix2_nfbin
        lcl_cfl_h  = 0.0_WP
        lcl_cfl_cs = 0.0_WP
        
        ! Max horizontal CFL over all edges
        do edge = 1, myDim_edge2D
            ednodes = edges(:,edge)
            el      = edge_tri(:,edge)
            
            dx1     = edge_cross_dxdy(1,edge)
            dy1     = edge_cross_dxdy(2,edge)
            dxdy12  = edge_dxdy(:,edge)*r_earth
            
            el2 = el(1)
            dx2 = 0.0_WP
            dy2 = 0.0_WP
            if (el(2) > 0) then
                el2     = el(2)
                dx2     = edge_cross_dxdy(3,edge)
                dy2     = edge_cross_dxdy(4,edge)
                dxdy12(1) = dxdy12(1) * (elem_cos(el(1)) + elem_cos(el2)) * 0.5_WP
                n_x     = (  dy1  + (-dy2))/2
                n_y     = ((-dx1) +   dx2 )/2
            else
                dxdy12(1) = dxdy12(1) * elem_cos(el(1))
                n_x     =  dy1
                n_y     = -dx1
            end if
            edlen = sqrt(dxdy12(1)**2 + dxdy12(2)**2)
            n_len = sqrt(n_x**2 + n_y**2)
            if (n_len > 0.0_WP) then
                n_x = n_x / n_len
                n_y = n_y / n_len
            end if
            
            do fbini = 2, nfbin-1
                u1  = vel(1, fbini, el(1))
                v1  = vel(2, fbini, el(1))
                u2  = vel(1, fbini, el2)
                v2  = vel(2, fbini, el2)
                Ue  = sqrt((0.5_WP*(u1+u2)*n_x)**2 + (0.5_WP*(v1+v2)*n_y)**2)
                CFL_h = Ue * dt / edlen
                lcl_cfl_h = max(lcl_cfl_h, CFL_h)
            end do
        end do
        
        ! Max cross-spectral CFL over all nodes and bins
        do node = 1, myDim_nod2D
            do fbini = 2, nfbin-1
                CFL_cs = abs(cs(fbini, node)) * dt / dphit(fbini)
                lcl_cfl_cs = max(lcl_cfl_cs, CFL_cs)
            end do
        end do
        
        ! Global max across MPI ranks
        glb_cfl_h  = 0.0_WP
        glb_cfl_cs = 0.0_WP
        call MPI_Allreduce(lcl_cfl_h , glb_cfl_h , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        call MPI_Allreduce(lcl_cfl_cs, glb_cfl_cs, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        
        if (mype == 0) then
            write(*,*) ' debug: CFL_max_horizontal(',trim(compartment_name),') = ', glb_cfl_h, &
                       ', CFL_max_cross_spectral(',trim(compartment_name),') = ', glb_cfl_cs
!             write(*,*) ' debug: dt_max_FE(',trim(compartment_name),') = ', dt/max(glb_cfl_h, glb_cfl_cs), ' s', &
!                        ', dt_max_AB2(',trim(compartment_name),') = ', 0.5_WP*dt/max(glb_cfl_h, glb_cfl_cs), ' s'
        end if
        
    end subroutine check_max_cfl


!     !===========================================================================
!     ! Compute velocity flux divergence consistent with the edge-based Superbee
!     ! flux scheme. For tracer h=1, the Superbee flux reduces to -vflux (since
!     ! gradients are zero for a constant field). This gives the numerical velocity
!     ! divergence as seen by the flux operator.
!     !
!     ! Used for compatible transport: E(n+1) = [E(n)+dt*Edivh] / [1+dt*div_vel]
!     ! which preserves constant fields exactly while maintaining conservation.
!     subroutine compute_vel_div_consistent(vel, div_vel, partit, mesh)
!         implicit none
!         type(t_partit), intent(inout), target :: partit
!         type(t_mesh),   intent(in),    target :: mesh
!         real(kind=WP),  intent(in)            :: vel(2, idemix2_nfbin, partit%myDim_elem2D+partit%eDim_elem2D)
!         real(kind=WP),  intent(inout)         :: div_vel(idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D)
!         
!         integer        :: edge, fbini, nfbin, node, el(2), el2, ednodes(2)
!         real(kind=WP)  :: dx1, dy1, dx2, dy2, u1, v1, u2, v2, vflux
!         
! #include "../associate_part_def.h"
! #include "../associate_mesh_def.h"
! #include "../associate_part_ass.h"
! #include "../associate_mesh_ass.h"
!         
!         nfbin = idemix2_nfbin
!         div_vel = 0.0_WP
!         
!         ! Edge loop: compute vflux using SAME geometry as Superbee routine
!         do edge = 1, myDim_edge2D
!             ednodes = edges(:, edge)
!             el      = edge_tri(:, edge)
!             
!             dx1 = edge_cross_dxdy(1, edge)
!             dy1 = edge_cross_dxdy(2, edge)
!             
!             dx2 = 0.0_WP
!             dy2 = 0.0_WP
!             el2 = el(1)
!             if (el(2) > 0) then
!                 el2 = el(2)
!                 dx2 = edge_cross_dxdy(3, edge)
!                 dy2 = edge_cross_dxdy(4, edge)
!             end if
!             
!             do fbini = 2, nfbin-1
!                 u1 = vel(1, fbini, el(1))
!                 v1 = vel(2, fbini, el(1))
!                 u2 = vel(1, fbini, el2)
!                 v2 = vel(2, fbini, el2)
!                 
!                 vflux = (-v1*dx1 + u1*dy1) + (v2*dx2 - u2*dy2)
!                 
!                 ! For tracer h=1: flux_h1 = -vflux
!                 ! Assembly matches Superbee: div(node1) += flux, div(node2) -= flux
!                 div_vel(fbini, ednodes(1)) = div_vel(fbini, ednodes(1)) - vflux
!                 div_vel(fbini, ednodes(2)) = div_vel(fbini, ednodes(2)) + vflux
!             end do
!         end do
!         
!         ! Assemble halo contributions to owner (same as for div_h)
!         call assemble_nod2D_fbin(div_vel, partit)
!         ! Distribute complete values to halos
!         call exchange_nod_fbin(div_vel, partit)
!         ! Normalize by area (same as for div_h)
!         do node = 1, myDim_nod2D + eDim_nod2D
!             do fbini = 2, nfbin-1
!                 div_vel(fbini, node) = div_vel(fbini, node) / areasvol(1, node)
!             end do
!         end do
!         
!     end subroutine compute_vel_div_consistent


!_______________________________________________________________________________
! Stage current state into slot 1 before writing restart (call before rotation).
! The two cases are mutually exclusive and complementary:
!   tip1==2: fresh E is in slot 2 → copy to slot 1; div already in slot 1 (ti=1)
!   tip1==1: fresh E already in slot 1; current div is in slot 2 (ti=2) → copy to slot 1
!   
!   the output stream is fixed on location iwe2_E_iw(:, :, 1) thats why we need it to
!   update  when iwe2_tip1==2 and its a restart writing moment
subroutine prepare_idemix2_restart()
    if (iwe2_tip1 == 2) then
        iwe2_E_iw(:, :, 1) = iwe2_E_iw(:, :, 2)
        if (allocated(iwe2_E_M2))  iwe2_E_M2(:, :, 1)  = iwe2_E_M2(:, :, 2)
        if (allocated(iwe2_E_niw)) iwe2_E_niw(:, :, 1) = iwe2_E_niw(:, :, 2)
    else
        if (allocated(iwe2_E_M2)) then
            iwe2_E_M2_divh(:, :, 1) = iwe2_E_M2_divh(:, :, 2)
            iwe2_E_M2_divs(:, :, 1) = iwe2_E_M2_divs(:, :, 2)
        end if
        if (allocated(iwe2_E_niw)) then
            iwe2_E_niw_divh(:, :, 1) = iwe2_E_niw_divh(:, :, 2)
            iwe2_E_niw_divs(:, :, 1) = iwe2_E_niw_divs(:, :, 2)
        end if
    end if
end subroutine prepare_idemix2_restart


!_______________________________________________________________________________
! After loading restart: slot 1 holds saved E and div.
! E: fill slot 2 for AB bootstrap (iwe2_ti=1,iwe2_tip1=2 at restart init).
! Div: saved div goes to slot 2 (=tip1 at init, read as "previous" on first step);
!      slot 1 (=ti) zeroed — it will be overwritten fresh in the first step.
subroutine apply_idemix2_restart()
    iwe2_E_iw(:, :, 2) = iwe2_E_iw(:, :, 1)
    if (allocated(iwe2_E_M2)) then
        iwe2_E_M2(:, :, 2)       = iwe2_E_M2(:, :, 1)
        iwe2_E_M2_divh(:, :, 2)  = iwe2_E_M2_divh(:, :, 1)
        iwe2_E_M2_divh(:, :, 1)  = 0.0_WP
        iwe2_E_M2_divs(:, :, 2)  = iwe2_E_M2_divs(:, :, 1)
        iwe2_E_M2_divs(:, :, 1)  = 0.0_WP
    end if
    if (allocated(iwe2_E_niw)) then
        iwe2_E_niw(:, :, 2)      = iwe2_E_niw(:, :, 1)
        iwe2_E_niw_divh(:, :, 2) = iwe2_E_niw_divh(:, :, 1)
        iwe2_E_niw_divh(:, :, 1) = 0.0_WP
        iwe2_E_niw_divs(:, :, 2) = iwe2_E_niw_divs(:, :, 1)
        iwe2_E_niw_divs(:, :, 1) = 0.0_WP
    end if
end subroutine apply_idemix2_restart


end module g_cvmix_idemix2

