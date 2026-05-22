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
! written by Patrick Scholz, 10.05.2019
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
    use g_config , only: dt, flag_debug
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
    logical            :: idemix2_enable_AB_timestep  = .true.
    integer            :: idemix2_AB_epsilon = 0.1_WP
    ! number of spectral bins used for the M2 tidal and near-inertial wave (niw) components
    ! e.g 50+2 in loop used as do fbin=2,51 ... 1 and 52 serve as spectral boudary condition
    ! The first and last bins are used for boundary conditions in the spectral space.
    ! They ensure smooth transitions and prevent numerical instabilities at the edges of the spectral domain.
    ! Numerical Stability: The advection scheme (e.g., Superbee) requires ghost cells or boundary values.
    ! Skipping the first and last bins avoids out-of-bounds errors when calculating gradients.
    integer            :: idemix2_nfbin=52
    
    ! enable adding of idemix full horizontal tendency from 
    ! div(grad(Eiw*v0)*v0*tauh) diffusion term
    logical            :: idemix2_enable_hor_diffusion = .false.
    
    ! enable idemix1 functionality of homogenous diffusion into all directions
    logical            :: idemix2_enable_hor_diff_iter = .false.
    integer            :: idemix2_hor_diff_niter       = 5   ! from Pollman et al. (2017)
    
    ! define shelf is defined as distance from coast (default 300km)
    real(kind=WP)      :: idemix2_shelf_dist = 300.0e3 
    
    ! scal down baroclinic wave speed when not all modes are used
    real(kind=WP)      :: idemix2_scal_cn = 1.0
    
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
    
    namelist /param_idemix2/ idemix2_tau_v, idemix2_tau_h, idemix2_gamma, idemix2_jstar, idemix2_mu0, idemix2_scal_cn, &
                             idemix2_enable_AB_timestep, idemix2_nfbin, & ! idemix2_enable_superbee_adv
                             idemix2_enable_hor_diffusion, idemix2_enable_hor_diff_iter, idemix2_hor_diff_niter, &
                             idemix2_shelf_dist, &
                             idemix2_botforc_Etot, &
                             idemix2_enable_M2  , idemix2_M2forc_file  , idemix2_M2forc_vname  , idemix2_M2forc_zname, &
                             idemix2_enable_niw , idemix2_niwforc_file , idemix2_niwforc_vname , idemix2_fniw_usage, &
!                              idemix2_enable_leew, idemix2_leewforc_file, idemix2_leewforc_vname, &
                             idemix2_enable_bot , idemix2_botforc_file , idemix2_botforc_vname , &
                             idemix2_hrmsforc_file, idemix2_hrmsforc_vname, &
                             idemix2_hlamforc_file, idemix2_hlamforc_vname
    
    !___________________________________________________________________________
    ! CVMIX-IDEMIX variables
    real(kind=WP), allocatable, dimension(:)    :: iwe2_phit, iwe2_phiu, iwe2_dphit, iwe2_dphiu
    
    ! M2 related global variables
    real(kind=WP)                               :: iwe2_omega_M2 
    real(kind=WP), allocatable, dimension(:,:,:):: iwe2_M2_uv, iwe2_E_M2, iwe2_E_M2_divh, iwe2_E_M2_divs
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_M2_w, iwe2_fM2, iwe2_E_M2_struct 
    real(kind=WP), allocatable, dimension(:)    :: w_M2_e, iwe2_alpha_M2_c, iwe2_M2_tau    
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_E_M2_dt, iwe2_E_M2_advh, iwe2_E_M2_advs, iwe2_E_M2_diss, iwe2_E_M2_forc
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_M2_div_uv

    ! niw related global variables
    real(kind=WP), allocatable, dimension(:)    :: iwe2_omega_niw
    real(kind=WP), allocatable, dimension(:,:,:):: iwe2_niw_uv, iwe2_E_niw, iwe2_E_niw_divh, iwe2_E_niw_divs
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_niw_w, iwe2_fniw, iwe2_E_niw_struct
    real(kind=WP), allocatable, dimension(:)    :: w_niw_e, iwe2_niw_tau
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_E_niw_dt, iwe2_E_niw_advh, iwe2_E_niw_advs, iwe2_E_niw_diss, iwe2_E_niw_forc
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_niw_div_uv
    
    ! general idemix variable
    real(kind=WP), allocatable, dimension(:)    :: iwe2_cn
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_c0, iwe2_v0, iwe2_alpha_c
    
    ! forcing realted variables
    real(kind=WP), allocatable, dimension(:)    :: iwe2_topo_hrms, iwe2_topo_hlam, iwe2_topo_dist
    real(kind=WP), allocatable, dimension(:)    :: iwe2_fbot_e, iwe2_fleew, iwe2_fsrf
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_fbot
    
    ! Eiw - internal wave energy related variables
    real(kind=WP), allocatable, dimension(:,:,:):: iwe2_E_iw
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_E_iw_diss, iwe2_E_iw_dt, iwe2_E_iw_fbot, & 
                                                   iwe2_E_iw_vdif, iwe2_E_iw_hdif
    real(kind=WP), allocatable, dimension(:)    :: iwe2_E_iw_fsrf, iwe_E_iw_vint
    
    real(kind=WP), allocatable, dimension(:,:)  :: iwe2_Av
    
    ! support variables 
    real(kind=WP), allocatable, dimension(:)    :: vol_nodB2T
    real(kind=WP), allocatable, dimension(:,:)  :: vol_wcelli
    integer      , allocatable, dimension(:,:)  :: edge_up_dn_tri      
    real(kind=WP), allocatable, dimension(:)    :: iwe2_grady_coriol, aux
    real(kind=WP), allocatable, dimension(:,:,:):: iwe2_gradxy_e, iwe2_gradxy_n
    real(kind=WP), allocatable, dimension(:,:  ):: iwe2_flx_uv, iwe2_flx_w
    integer                                     :: iwe2_taum1=1, iwe2_tau=2, iwe2_taup1=3, otaum1
    
    ! real(kind=WP), allocatable, dimension(:)    :: iwe2_cn_e, iwe2_cn_gradx, iwe2_cn_grady
    
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
                                       fbin_i, elnodes(3), nzmax
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
            write(*,*) "     IDEMIX2 parameters:"
            write(*,*) "     ├> idemix2_tau_v               = ", idemix2_tau_v
            write(*,*) "     ├> idemix2_tau_h               = ", idemix2_tau_h
            write(*,*) "     ├> idemix2_gamma               = ", idemix2_gamma
            write(*,*) "     ├> idemix2_jstar               = ", idemix2_jstar
            write(*,*) "     ├> idemix2_mu0                 = ", idemix2_mu0
            write(*,*) "     │                                "
!             write(*,*) "     ├> idemix2_superbee_adv        = ", idemix2_enable_superbee_adv
            write(*,*) "     ├> idemix2_AB_timestep         = ", idemix2_enable_AB_timestep
            write(*,*) "     ├> idemix2_nfbin               = ", idemix2_nfbin
            write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_enable_M2           = ", idemix2_enable_M2
            write(*,*) "     │  └> idemix2_M2forc_file      = ", trim(idemix2_M2forc_file)
            write(*,*) "     │     ├> idemix2_M2forc_vname  = ", trim(idemix2_M2forc_vname)
            write(*,*) "     │     └> idemix2_M2forc_zname  = ", trim(idemix2_M2forc_zname)
            write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_enable_niw          = ", idemix2_enable_niw 
            write(*,*) "     │  ├> idemix2_fniw_usage       = ", idemix2_fniw_usage 
            write(*,*) "     │  └> idemix2_niwforc_file     = ", trim(idemix2_niwforc_file)
            write(*,*) "     │     └> idemix2_niwforc_vname = ", trim(idemix2_niwforc_vname)
            write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_enable_bot          = ", idemix2_enable_bot
            write(*,*) "     │  └> idemix2_botforc_file     = ", trim(idemix2_botforc_file)
            write(*,*) "     │     ├> idemix2_botforc_vname = ", trim(idemix2_botforc_vname)
            write(*,*) "     │     └> idemix2_botforc_Etot  = ", idemix2_botforc_Etot
            write(*,*) "     │                                "
!             write(*,*) "     ├> idemix2_enable_leew         = ", idemix2_enable_leew
!             write(*,*) "     │  └> idemix2_leewforc_file    = ", trim(idemix2_leewforc_file)
!             write(*,*) "     │     └> idemix2_leewforc_vname= ", trim(idemix2_leewforc_vname)
!             write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_hrmsforc_file       = ", trim(idemix2_hrmsforc_file)
            write(*,*) "     │  └> idemix2_hrmsforc_vname   = ", trim(idemix2_hrmsforc_vname)
            write(*,*) "     │                                "
            write(*,*) "     ├> idemix2_hlamforc_file       = ", trim(idemix2_hlamforc_file)
            write(*,*) "     │  └> idemix2_hlamforc_vname   = ", trim(idemix2_hlamforc_vname)
            write(*,*) "     │                                "
            WRITE(*,*) "     └> idemix2_shelf_dist         = ", idemix2_shelf_dist
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
        
!         allocate(iwe2_cn_e(elem_size))
!         allocate(iwe2_cn_gradx(elem_size))
!         allocate(iwe2_cn_grady(elem_size))
!         iwe2_cn_e(:)         = 0.0_WP
!         iwe2_cn_gradx(:)     = 0.0_WP
!         iwe2_cn_grady(:)     = 0.0_WP
        
        ! Eiw 3d enery dissipation coefficient 
        allocate(iwe2_alpha_c(nl,node_size))
        iwe2_alpha_c(:,:)    = 0.0_WP
        
        ! initialise Eiw - internal wave energy variables
        ! index (1:3,...) timestep index E^(n-1), E^(n), E^(n+1)
        allocate(iwe2_E_iw(   3, nl, node_size))
        allocate(iwe2_E_iw_diss(nl, node_size), iwe2_E_iw_dt(nl, node_size), iwe2_E_iw_fbot(nl, node_size))
        allocate(iwe2_E_iw_vdif(nl, node_size), iwe2_E_iw_hdif(nl, node_size))
        allocate(iwe2_E_iw_fsrf(node_size))
        iwe2_E_iw(     :,:,:)= 0.0_WP
        
        ! Eiw production terms: 
        iwe2_E_iw_dt(    :,:)= 0.0_WP ! total Eiw production
        iwe2_E_iw_diss(  :,:)= 0.0_WP ! Eiw production from disspation 
        iwe2_E_iw_fbot(  :,:)= 0.0_WP ! Eiw production from bottom forcing 
        iwe2_E_iw_hdif(  :,:)= 0.0_WP ! Eiw production from horizontal diffusion
        iwe2_E_iw_vdif(  :,:)= 0.0_WP ! Eiw production from vertical diffusion
        iwe2_E_iw_fsrf(    :)= 0.0_WP ! Eiw production from surface forcing 
        
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
            
            allocate(iwe2_M2_div_uv(    nfbin,elem_size))
            iwe2_M2_div_uv        = 0.0_WP
            
            ! M2 forcing
            allocate(iwe2_fM2(nfbin, node_size))
            iwe2_fM2(:,:)         = 0.0_WP
            
            ! M2 wave energy, and divergence of M2 wave energy 
            ! index (1:3,...) timestep index E^(n-1), E^(n), E^(n+1)
            allocate(  iwe2_E_M2(     3, nfbin, node_size) &
                     , iwe2_E_M2_divh(3, nfbin, node_size) &
                     , iwe2_E_M2_divs(3, nfbin, node_size) &
                    )
            iwe2_E_M2(      :,:,:)= 0.0_WP
            iwe2_E_M2_divh( :,:,:)= 0.0_WP
            iwe2_E_M2_divs( :,:,:)= 0.0_WP
            
            ! structure function for M2 energy 
            allocate(iwe2_E_M2_struct(nl, node_size))
            iwe2_E_M2_struct(:,:) = 0.0_WP
            
            ! diagnostic for M2 spectral energy advection 
            allocate(  iwe2_E_M2_dt(  nfbin, myDim_nod2D)  &
                     , iwe2_E_M2_advh(nfbin, myDim_nod2D)  &
                     , iwe2_E_M2_advs(nfbin, myDim_nod2D)  &
                     , iwe2_E_M2_diss(nfbin, myDim_nod2D)  &
                     , iwe2_E_M2_forc(nfbin, myDim_nod2D)  &
                     )
            iwe2_E_M2_dt(  :,:)   = 0.0_WP
            iwe2_E_M2_advh(:,:)   = 0.0_WP
            iwe2_E_M2_advs(:,:)   = 0.0_WP
            iwe2_E_M2_diss(:,:)   = 0.0_WP
            iwe2_E_M2_forc(:,:)   = 0.0_WP
            
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
            
            allocate(iwe2_niw_div_uv(    nfbin,elem_size))
            iwe2_niw_div_uv        = 0.0_WP
            
            ! niw forcing
            allocate(iwe2_fniw(nfbin, node_size))
            iwe2_fniw(:,:)         = 0.0_WP
            
            ! niw wave energy, and divergence of niw wave energy 
            ! index (1:3,...) timestep index E^(n-1), E^(n), E^(n+1)
            allocate(  iwe2_E_niw(     3, nfbin, node_size)   &
                     , iwe2_E_niw_divh(3, nfbin, node_size)   &
                     , iwe2_E_niw_divs(3, nfbin, node_size)   &
                     )
            iwe2_E_niw(      :,:,:)= 0.0_WP
            iwe2_E_niw_divh( :,:,:)= 0.0_WP
            iwe2_E_niw_divs( :,:,:)= 0.0_WP
            
            ! structure function for niw energy 
            allocate(iwe2_E_niw_struct(nl, node_size))
            iwe2_E_niw_struct(:,:) = 0.0_WP
        
            ! diagnostic for niw spectral energy advection 
            allocate(  iwe2_E_niw_dt(  nfbin, myDim_nod2D)   &
                     , iwe2_E_niw_advh(nfbin, myDim_nod2D)   &
                     , iwe2_E_niw_advs(nfbin, myDim_nod2D)   &
                     , iwe2_E_niw_diss(nfbin, myDim_nod2D)   &
                     , iwe2_E_niw_forc(nfbin, myDim_nod2D)   &
                     )
            iwe2_E_niw_dt(  :,:)   = 0.0_WP
            iwe2_E_niw_advh(:,:)   = 0.0_WP
            iwe2_E_niw_advs(:,:)   = 0.0_WP
            iwe2_E_niw_diss(:,:)   = 0.0_WP
            iwe2_E_niw_forc(:,:)   = 0.0_WP
            
        end if 
        
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
        allocate(iwe2_Av(nl,elem_size))
        iwe2_Av(:,:)         = 0.0_WP

        ! forcing fields, M2 tidal forcing (spectral) and NIW forcing (spectral)
        allocate(  iwe2_fbot_e(  myDim_elem2D)    &
                 , iwe2_fbot(nl, node_size)       &
                 , iwe2_fsrf(    node_size)     &
!                  , iwe2_fleew(   myDim_elem2D)    &                 
                 )
        iwe2_fbot_e(:)       = 0.0_WP
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
                end if
            end if 
            if (.not. file_exist) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
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
                                    .false.                          , & ! NN for missing values
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
                write(*,*) ' ERROR: IDEMIX2 NIW forcing file not found! Cant apply IDEMIX2'
                write(*,*) '        vertical mixing parameterisation! '
                write(*,*) '        ├> file: ', trim(idemix2_niwforc_file)
                write(*,*) '        └> check: namelist.cvmix, idemix2_niwforc_file &  '
                write(*,*) '____________________________________________________________________'
                print *, achar(27)//'[0m'
                write(*,*)
            end if
        end if 
        if (.not. file_exist) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        
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
                call read_other_NetCDF(trim(idemix2_botforc_file), trim(idemix2_botforc_vname), 1, iwe2_fbot_e, .true., .false., partit, mesh)
                !                                                                                                  |           
                !                                 .false.=interpolate on element centroids instead of vertices <---+   
                
                ! make sure forcing is all positive no numerical negative values
                iwe2_fbot_e = max(0.0_WP,iwe2_fbot_e)
                
                
                ! check for total tidal energy that is infused through the bottom, see how 
                ! much is lossed during interpolation and compare with value of the 
                ! original files
                loc_Etot = 0.0_WP
                do elem=1, myDim_elem2D
                    ! REMEMBER!!!: the partition on elements is not unique there are 
                    ! elements that belong to two CPUs. For unique elements the index
                    ! of the First trinagle node must  be <= myDim_nod2D
                    if (elem2D_nodes(1,elem)<=myDim_nod2D) then
                        loc_Etot = loc_Etot + elem_area(elem)*iwe2_fbot_e(elem)
                    end if     
                end do
                call MPI_AllREDUCE(loc_Etot, glb_Etot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
                if (mype==0) write(*,*) "     │  └> IDEMIX2 total tidal energy Etot_bot =", glb_Etot*1.0e-12, ' TW'
                
                ! normalize total tidal energy at bottom with respect to the total 
                ! tidal energy that is e.g in the original forcing files to accomodate 
                ! non concerving losses during interpolation. This is only done when 
                ! in namelist.cvmix: idemix2_botforc_Etot \= 0.0_WP
                if (idemix2_botforc_Etot /= 0.0_WP) then
                    iwe2_fbot_e = iwe2_fbot_e * idemix2_botforc_Etot/glb_Etot
                    
                    loc_Etot = 0.0_WP
                    do elem=1, myDim_elem2D
                        ! REMEMBER!!!: the partition on elements is not unique there are 
                        ! elements that belong to two CPUs. For unique elements the index
                        ! of the First trinagle node must  be <= myDim_nod2D
                        if (elem2D_nodes(1,elem)<=myDim_nod2D) then
                            loc_Etot = loc_Etot + elem_area(elem)*iwe2_fbot_e(elem)
                        end if     
                    end do
                    call MPI_AllREDUCE(loc_Etot, glb_Etot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
                    if (mype==0) write(*,*) "     │  └> IDEMIX2 Etot_bot after normalizing =", glb_Etot*1.0e-12, ' TW'
                end if 
                
                ! divide by density_0 --> convert from W/m^2 to m^3/s^3
                iwe2_fbot_e  = iwe2_fbot_e/density_0
                
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
                end if 
            end if 
            if (.not. file_exist) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
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
!                 iwe2_fbot_e = iwe2_fbot_e + iwe2_fleew
!                 !                           |
!                 !                           +-> no dividion by density_0, Lee wave Forcing
!                 !                               already in units of m^3/s^3
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
!                 end if
!             end if 
!             if (.not. file_exist) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
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
                end if
            end if 
            if (.not. file_exist) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
            
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
            if (.not. file_exist) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
            
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
        ! convert elemental bot forcing into nodal bottom forcing 
        do elem = 1, myDim_elem2D
            nzmax = nlevels(elem)
            elnodes = elem2d_nodes(:,elem)
            iwe2_fbot(nzmax, elnodes(:)) = iwe2_fbot(nzmax, elnodes(:)) + iwe2_fbot_e(elem)/3.0
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
        ! compute d/dy of coriolis 
        call find_up_downwind_triangles(partit, mesh, mesh%edge_up_dn_tri)
        
        !_______________________________________________________________________
        ! compute scalar cell volume from top to bottom 
        ! --> if take the fixed standard levels (zbar) to compute the ocean volume 
        !     than only normal numerical drift, ocean volume can be computed at 
        !     initialisation point
        call compute_vol_nodB2T_fix(vol_nodB2T, mesh, partit)
        
        !_______________________________________________________________________
        ! initialise IDEMIX parameters
        call cvmix_idemix2_init(                                 &
              tau_v               = idemix2_tau_v                &
            , tau_h               = idemix2_tau_h                &
            , gamma               = idemix2_gamma                &
            , jstar               = idemix2_jstar                &
            , mu0                 = idemix2_mu0                  &
            , nfbin               = idemix2_nfbin                &
            , shelf_dist          = idemix2_shelf_dist           &
            , enable_M2           = idemix2_enable_M2            &
            , enable_niw          = idemix2_enable_niw           &
            , enable_AB_timestep  = idemix2_enable_AB_timestep   &
            , enable_hor_diffusion= idemix2_enable_hor_diffusion &
            , enable_hor_diff_iter= idemix2_enable_hor_diff_iter &
            , hor_diff_niter      = idemix2_hor_diff_niter       &
            ! enable_superbee_adv = idemix2_enable_superbee_adv  &
            )
                                
        if (partit%mype==0) write(*,*)                  
    end subroutine init_cvmix_idemix2
    
    
    
    !
    !
    !
    !===========================================================================
    ! calculate IDEMIX2 internal wave energy and its dissipation
    subroutine calc_cvmix_idemix2(partit, mesh)   
        implicit none
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(inout), target :: partit
        integer       :: node, elem, edge, node_size, elem_size, k, fbini, nfbin, nodeH_size, elemH_size
        integer       :: nz, nln, nl1, nl2, nl12, nu1, nu2, nu12, uln, iter  
        integer       :: elnodes(3), el(2), ednodes(2) 
        real(kind=WP) :: lat_n_deg, lat_e_deg
        real(kind=WP) :: cn, cn_e, cn_gradx_e, cn_grady_e, omega_niw_e
        logical       :: topo_shelf=.False.
        logical       :: debug=.false.

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        node_size  = myDim_nod2D
        nodeH_size = myDim_nod2D+eDim_nod2D
        elem_size  = myDim_elem2D
        elemH_size = myDim_elem2D+eDim_elem2D
        nfbin      = idemix2_nfbin
        !_______________________________________________________________________
        ! do tau timestepping indices shift
        otaum1= iwe2_taum1 
        iwe2_taum1 = iwe2_tau
        iwe2_tau   = iwe2_taup1
        iwe2_taup1 = otaum1;
        ! if (mype==0) write(*,*) " --> taum1, tau, taup1: ", iwe2_taum1, iwe2_tau, iwe2_taup1
        
        !_______________________________________________________________________
        ! compute scalar cell volume from top to bottom 
        ! --> if we take not the ALE variying ocean layer into account --> slight 
        !     numerical drift
        ! --> if take the fixed standard levels (zbar) to compute the ocean volume 
        !     than only normal numerical drift, ocean volume can be computed at 
        !     initialisation point
        ! call compute_vol_nodB2T_hnode(vol_nodB2T, mesh, partit)
        
        !_______________________________________________________________________
        do node=1, myDim_nod2D+eDim_nod2D
            !___________________________________________________________________
            ! re-initialse cross spectral velocites, for later accumulation through 
            ! elemental values
            if (idemix2_enable_M2 ) iwe2_M2_w( :,node) = 0.0_WP
            if (idemix2_enable_niw) iwe2_niw_w(:,node) = 0.0_WP
            
            !___________________________________________________________________
            lat_n_deg = geo_coord_nod2D(2,node) * 180.0/pi   
            
            !___________________________________________________________________
            ! compute baroclinic gravity wave speed
            nln = mesh%nl
            uln = 1
            do nz=1, nod_in_elem2D_num(node)
                nln=min(nln, nlevels(nod_in_elem2D(nz, node)))
                uln=max(uln, ulevels(nod_in_elem2D(nz, node)))
            end do
            cn=0.0_WP
            do nz=uln, nln-1
                cn=cn+hnode(nz,node)*sqrt(max(bvfreq(nz,node), 0._WP) + max(bvfreq(nz+1,node), 0._WP))/2._WP
            end do
            cn = cn/pi * idemix2_scal_cn
            iwe2_cn(node)=cn
            
            !
            !
            !___________________________________________________________________
            ! 1st. compute idemix2 parameter over the vertical water column and 
            ! local horizontal parameters on elements !!!:
            !  IN VARIABLES:
            ! --------------
            ! nlev      ... number of mid-depth levels
            ! coriolis  ... coriolis parameter 
            ! Nsqr      ... squared buoyancy frequency
            ! cn        ... baroclionic gravity wave speed
            ! 
            ! OUT VARIABLES:
            ! --------------
            ! alpha_c   ... 3d enery dissipation coefficient 
            ! v0        ... 3d horizontal group velocities (m/s) for the continuous 
            !               internal wave spectrum
            ! c0        ... 3d vertical group velocities (m/s) for the continuous 
            !               internal wave spectrum
            call cvmix_idemix2_compute_param(            &
                  nlev     = nln-uln+1                   & !IN 
                , coriolis = mesh%coriolis_node(   node) & !IN 
                , Nsqr     = bvfreq(      uln:nln, node) & !IN
                , cn       = cn                          & !IN
                , alpha_c  = iwe2_alpha_c(uln:nln, node) & !OUT
                , c0       = iwe2_c0(     uln:nln, node) & !OUT
                , v0       = iwe2_v0(     uln:nln, node) & !OUT
                )
            
            !
            !
            !___________________________________________________________________
            ! 3rd. compute dissipation time scales and rates
            !  IN VARIABLES:
            ! --------------
            ! dtime     ... time step
            ! lat       ... elem lat coordinates
            ! coriolis  ... coriolis parameter 
            ! omega_M2  ... M2 frequency (fixed)
            ! omega_niw ... NIW frequency (lat dependent)
            ! cn        ... baroclionic gravity wave speed
            ! zbottom   ... bottom depth at nodes (depth must be positive, *(-1))
            ! topo_hrms ... root mean square topographic height
            ! topo_hlam ... characteristic horiz. length scale of topographic features
            ! topo_shelf... .true. if  point is considered on shelf
            ! 
            ! OUT VARIABLES:
            ! --------------
            ! tau_M2    ... M2 Tidal Dissipation Timescale     
            ! tau_niw   ... NIW Dissipation Timescale
            ! alpha_M2_cont... M2 Continuous Dissipation Rate (Background M2 tidal 
            !                  dissipation rate
            if (idemix2_enable_M2 .or. idemix2_enable_niw) then 
                if (iwe2_topo_dist(node) <= idemix2_shelf_dist) then
                    topo_shelf = .True.
                else    
                    topo_shelf = .False.
                end if
                if (idemix2_enable_M2) then 
                    call cvmix_idemix2_compute_compart_interact_tscale(&
                          dtime        = dt                          & !IN 
                        , lat          = lat_n_deg                   & !IN
                        , coriolis     = mesh%coriolis_node(   node) & !IN 
                        , cn           = cn                          & !IN
                        , zbottom      = -zbar_n_bot(          node) & !IN
                        , topo_hrms    = iwe2_topo_hrms(       node) & !IN
                        , topo_hlam    = iwe2_topo_hlam(       node) & !IN
                        , topo_shelf   = topo_shelf                  & !IN
                        , omega_compart= iwe2_omega_M2               & !IN
                        , tau_compart  = iwe2_M2_tau(          node) & !OUT
                        , compart_name = 'M2'                        &
                        )
                    call cvmix_idemix2_compute_M2_dissipation(       &
                          lat          = lat_n_deg                   & !IN
                        , cn           = cn                          & !IN
                        , zbottom      = -zbar_n_bot(          node) & !IN
                        , alpha_M2_c   = iwe2_alpha_M2_c(      node) & !OUT
                        )    
                end if 
                if (idemix2_enable_niw) then 
                    call cvmix_idemix2_compute_compart_interact_tscale(&
                          dtime        = dt                          & !IN 
                        , lat          = lat_n_deg                   & !IN
                        , coriolis     = mesh%coriolis_node(   node) & !IN 
                        , cn           = cn                          & !IN
                        , zbottom      = -zbar_n_bot(          node) & !IN
                        , topo_hrms    = iwe2_topo_hrms(       node) & !IN
                        , topo_hlam    = iwe2_topo_hlam(       node) & !IN
                        , topo_shelf   = topo_shelf                  & !IN
                        , omega_compart= iwe2_omega_niw(       node) & !IN
                        , tau_compart  = iwe2_niw_tau(         node) & !OUT
                        , compart_name = 'niw'                       &
                        )
                end if  
            end if 
            
            !
            !
            !___________________________________________________________________
            ! 4th. compute structure function for M2 and NIW energy
            !  IN VARIABLES:
            ! --------------
            ! nlev      ... number of full depth levels
            ! dzw       ... layyer thickness, distance between full depth levels
            ! coriolis  ... coriolis parameter 
            ! omega_M2  ... M2 frequency (fixed)
            ! omega_niw ... NIW frequency (lat dependent)
            ! cn        ... baroclionic gravity wave speed
            ! Nsqr      ... squared buoyancy frequency
            ! 
            ! OUT VARIABLES:
            ! --------------
            ! E_struct_M2 ... structure function for M2 energy     
            ! E_struct_niw... structure function for NIW energy
            if ((idemix2_enable_M2) .and. (idemix2_enable_niw)) then 
                call cvmix_idemix2_compute_vert_struct_fct(          &
                      nlev        = nln-uln+1                        & !IN 
                    , dzw         = hnode(          uln:nln-1, node) & !IN
                    , coriolis    = mesh%coriolis_node(        node) & !IN 
                    , cn          = cn                               & !IN
                    , Nsqr        = bvfreq(           uln:nln, node) & !IN
                    , omega_M2    = iwe2_omega_M2                    & !IN
                    , omega_niw   = iwe2_omega_niw(            node) & !IN
                    , E_struct_M2 = iwe2_E_M2_struct( uln:nln, node) & !OUT
                    , E_struct_niw= iwe2_E_niw_struct(uln:nln, node) & !OUT
                    )
                    
            else if (idemix2_enable_M2) then 
                call cvmix_idemix2_compute_vert_struct_fct(          &
                      nlev        = nln-uln+1                        & !IN 
                    , dzw         = hnode(          uln:nln-1, node) & !IN
                    , coriolis    = mesh%coriolis_node(        node) & !IN 
                    , cn          = cn                               & !IN
                    , Nsqr        = bvfreq(           uln:nln, node) & !IN
                    , omega_M2    = iwe2_omega_M2                    & !IN
                    , E_struct_M2 = iwe2_E_M2_struct( uln:nln, node) & !OUT
                    )
                    
            else if (idemix2_enable_niw) then 
                call cvmix_idemix2_compute_vert_struct_fct(          &
                      nlev        = nln-uln+1                        & !IN 
                    , dzw         = hnode(          uln:nln-1, node) & !IN
                    , coriolis    = mesh%coriolis_node(        node) & !IN 
                    , cn          = cn                               & !IN
                    , Nsqr        = bvfreq(           uln:nln, node) & !IN
                    , omega_niw   = iwe2_omega_niw(            node) & !IN
                    , E_struct_niw= iwe2_E_niw_struct(uln:nln, node) & !OUT
                    )
                    
            end if         
        end do ! --> do node=1, myDim_nod2D+eDim_nod2D
        
        
        
        !_______________________________________________________________________
        ! 5th. compute idemix2 group velocites for M2 an NIW
        do elem = 1, myDim_elem2D
            nln        = nlevels(elem)-1
            uln        = ulevels(elem)
            elnodes    = elem2d_nodes(:,elem)
            lat_e_deg  = sum(geo_coord_nod2D(2,elnodes))/3.0 * 180.0/pi    
            if (idemix2_enable_M2)  w_M2_e(:)  = 0.0_WP
            if (idemix2_enable_niw) w_niw_e(:) = 0.0_WP
            
            !___________________________________________________________________
            ! compute baroclionic gravity wave speed on elements
            cn_e        = sum(iwe2_cn(elnodes))/3.0
            
            ! compute gradient of  baroclionic gkdot_x_M2 = sqrt(fxa)/omega_M2*cn_gradxravity wave speed on elements
            cn_gradx_e  = sum(gradient_sca(1:3,elem)*iwe2_cn(elnodes))
            cn_grady_e  = sum(gradient_sca(4:6,elem)*iwe2_cn(elnodes))
            
!             iwe2_cn_e(elem)     = cn_e
!             iwe2_cn_gradx(elem) = cn_gradx_e
!             iwe2_cn_grady(elem) = cn_grady_e
            
            ! average to elem
            if (idemix2_enable_niw) omega_niw_e = sum(iwe2_omega_niw(elnodes))/3.0
            
            !
            !
            !___________________________________________________________________
            !  IN VARIABLES:
            ! --------------
            ! nlev      ... number of mid-depth levels
            ! nfbin     ... number of spectral bins used for the M2 tidal and 
            !               near-inertial wave (niw)
            ! dtime     ... time step
            ! coriolis  ... coriolis parameter 
            ! grady_coriol ... lat gradient of coriolis paramter
            ! coslat    ... cosine of latitude @elem
            ! Nsqr      ... squared buoyancy frequency
            ! omega_M2  ... M2 frequency (fixed)
            ! omega_niw ... NIW frequency (lat dependent)
            ! cn        ... baroclionic gravity wave speed
            ! cn_gradx  ... zonal gradient of  baroclionic gravity wave speed
            ! cn_grady  ... merid gradient of  baroclionic gravity wave speed
            ! phit      ... edge of the k-th spectral bin
            ! phiu      ... center of the k-th spectral bin
            ! 
            ! OUT VARIABLES:
            ! --------------
            ! alpha_c   ... 3d enery dissipation coefficient 
            ! c0        ... 3d horizontal group velocities (m/s) for the continuous 
            !              z internal wave spectrum
            ! v0        ... 3d vertical group velocities (m/s) for the continuous 
            !               internal wave spectrum
            ! cg_M2     ... 2d group velocity of M2 internal tidal waves 
            ! cg_niw    ... 2d group velocity of near-inertial waves (NIW)
            ! u_M2      ... zonal component of M2 internal tide group velocity
            ! v_M2      ... meridional component of M2 internal tide group velocity
            ! w_M2      ... cross spectral propagation of M2 internal tide group velocity
            ! u_niw     ... zonal component of NIW internal group velocity
            ! v_niw     ... meridional component of NIW internal group velocity
            ! w_niw     ... cross spectral propagation of NIW internal group velocity
            if (idemix2_enable_M2) then 
                call cvmix_idemix2_compute_compart_groupvel(     &
                      nfbin        = idemix2_nfbin               & !IN
                    , coriolis     = mesh%coriolis(        elem) & !IN 
                    , coriol_grady = iwe2_grady_coriol(    elem) & !IN
                    , coslat       = mesh%elem_cos(        elem) & !IN
                    , cn           = cn_e                        & !IN
                    , cn_gradx     = cn_gradx_e                  & !IN
                    , cn_grady     = cn_grady_e                  & !IN
                    , phit         = iwe2_phit                   & !IN 
                    , phiu         = iwe2_phiu                   & !IN 
                    , omega_compart= iwe2_omega_M2               & !IN
                    , u_compart    = iwe2_M2_uv(1,      :, elem) & !OUT
                    , v_compart    = iwe2_M2_uv(2,      :, elem) & !OUT
                    , w_compart    = w_M2_e(            :      ) & !OUT
                    )
                ! --> here w_M2_e and w_niw_e are still on elements but for the proper 
                !     finite volume advection implementation we need them on nodes !
                do k=1,3
                    iwe2_M2_w(:, elnodes(k)) = iwe2_M2_w(:, elnodes(k)) + w_M2_e(:)*elem_area(elem)/3.0_WP
                end do
                
            end if
            if (idemix2_enable_niw) then 
                call cvmix_idemix2_compute_compart_groupvel(      &
                      nfbin        = idemix2_nfbin                & !IN
                    , coriolis     = mesh%coriolis(        elem)  & !IN 
                    , coriol_grady = iwe2_grady_coriol(    elem)  & !IN
                    , coslat       = mesh%elem_cos(        elem)  & !IN
                    , cn           = cn_e                         & !IN
                    , cn_gradx     = cn_gradx_e                   & !IN
                    , cn_grady     = cn_grady_e                   & !IN
                    , phit         = iwe2_phit                    & !IN 
                    , phiu         = iwe2_phiu                    & !IN 
                    , omega_compart= omega_niw_e                  & !IN
                    , u_compart    = iwe2_niw_uv(1,      :, elem) & !OUT
                    , v_compart    = iwe2_niw_uv(2,      :, elem) & !OUT
                    , w_compart    = w_niw_e(            :      ) & !OUT
                    )
                
                ! Mask out NIW velocities near equator where physics is singular
                if (abs(lat_e_deg) < 5.0_WP) then  ! ~5° latitude
                    iwe2_niw_uv( :, :, elem) = 0.0_WP
                    w_niw_e(           :   ) = 0.0_WP
                end if
                
                do k=1,3
                    iwe2_niw_w(:, elnodes(k)) = iwe2_niw_w(:, elnodes(k)) + w_niw_e(:)*elem_area(elem)/3.0_WP
                end do
            end if
            ! --> at the end we still need to normalize iwe2_M2_w and iwe2_niw_w 
            !     with the scalararea!
        end do ! --> do elem = 1, myDim_elem2D         
        
        ! finalize elem2node averaging of iwe2_M2_w and iwe2_niw_w
        ! cross spectral exachange has to be related to nodes, since general advection 
        ! is related to nodes independent of the volume
        if (idemix2_enable_M2) then 
            call exchange_elem_fbin(iwe2_M2_uv, partit)
            call exchange_nod_fbin(iwe2_M2_w, partit)
            do node=1, myDim_nod2D+eDim_nod2D
                iwe2_M2_w( :, node) = iwe2_M2_w(:, node)/area(1, node)
            end do
        end if
        
        if (idemix2_enable_niw) then 
            call exchange_elem_fbin(iwe2_niw_uv, partit)
            call exchange_nod_fbin(iwe2_niw_w, partit)
            do node=1, myDim_nod2D+eDim_nod2D
                iwe2_niw_w( :, node) = iwe2_niw_w(:, node)/area(1, node)
            end do            
        end if
        
        
        !_______________________________________________________________________
        ! 6th. horizontal spectral integrate energy compartment E_M2^(n+1) and 
        ! E_niw^(n+1) equation. 
        ! dE/dt = -div(vec_u * E) - tau*E + Forc    
        ! E^(n+1) = E^n + dt*( -div(vec_u^n*E^n) - tau*E^n + Forc^n)
        if (idemix2_enable_M2) then 
            call hsintegrate_Ecompart(              &
                  iwe2_taum1, iwe2_tau, iwe2_taup1  &
                , 'M2'                              &  
                , iwe2_E_M2                         & 
                , iwe2_E_M2_divh                    &
                , iwe2_E_M2_divs                    &
                , iwe2_M2_uv                        &
                , iwe2_M2_w                         &
                , iwe2_fM2                          &
                , iwe2_M2_tau                       &
                , iwe2_dphit                        &
                , iwe2_gradxy_e                     &
                , iwe2_gradxy_n                     &
                , iwe2_flx_uv                       &
                , iwe2_flx_w                        &
                , vol_nodB2T                        &
                , partit                            &
                , mesh                              &
                , .True.                            & ! flag_AB2, do 2nd order adams-bashfort in time
                , iwe2_E_M2_dt                      & ! optional: diagnostic
                , iwe2_E_M2_advh                    & ! optional: diagnostic
                , iwe2_E_M2_advs                    & ! optional: diagnostic
                , iwe2_E_M2_diss                    & ! optional: diagnostic
                , iwe2_E_M2_forc                    & ! optional: diagnostic
                )
                
            ! Compute global total NIW energy for conservation check
            call check_global_energy(iwe2_E_M2, iwe2_taup1, iwe2_dphit, vol_nodB2T, iwe2_fM2, iwe2_M2_tau, partit, mesh, 'M2')
!             call check_flux_conservation(  iwe2_E_niw_divh, iwe2_E_niw_divs, iwe2_tau, iwe2_dphit, vol_nodB2T, partit, mesh, 'M2')
        end if 
        
        if (idemix2_enable_niw) then
            call hsintegrate_Ecompart(              &
                  iwe2_taum1, iwe2_tau, iwe2_taup1  & 
                , 'NIW'                             &
                , iwe2_E_niw                        & 
                , iwe2_E_niw_divh                   &
                , iwe2_E_niw_divs                   &
                , iwe2_niw_uv                       &
                , iwe2_niw_w                        &
                , iwe2_fniw                         &
                , iwe2_niw_tau                      &
                , iwe2_dphit                        &
                , iwe2_gradxy_e                     &
                , iwe2_gradxy_n                     &
                , iwe2_flx_uv                       &
                , iwe2_flx_w                        &
                , vol_nodB2T                        &
                , partit                            &
                , mesh                              &
                , .True.                            & ! flag_AB2, do 2nd order adams-bashfort in time
                , iwe2_E_niw_dt                     & ! optional: diagnostic
                , iwe2_E_niw_advh                   & ! optional: diagnostic
                , iwe2_E_niw_advs                   & ! optional: diagnostic
                , iwe2_E_niw_diss                   & ! optional: diagnostic
                , iwe2_E_niw_forc                   & ! optional: diagnostic
                )
            
            ! Compute global total NIW energy for conservation check
            call check_global_energy(iwe2_E_niw, iwe2_taup1, iwe2_dphit, vol_nodB2T, iwe2_fniw, iwe2_niw_tau, partit, mesh, 'niw')
!             call check_max_cfl(iwe2_niw_uv, iwe2_niw_w, iwe2_dphit, partit, mesh, 'niw')
!             call check_flux_conservation(  iwe2_E_niw_divh, iwe2_E_niw_divs, iwe2_tau, iwe2_dphit, vol_nodB2T, partit, mesh, 'niw')
        end if 
        
        
        
        !_______________________________________________________________________
        ! 7th. Integrate IDEMIX equation vertical, solve vertical diffusion and 
        ! dissipation part implicitly
        ! Eiw^(t+1) = Eiw^(t) + dt*[  d/dz( c_0 * tau_v * d/dz(c_0*E_iw))^(t+1) 
        !                           - alpha_c*Eiw^(t+1)
        !                           + Forc^(t) ]
        do node = 1, myDim_nod2D
            uln = ulevels_nod2D(node)
            nln = nlevels_nod2D(node)
            call cvmix_idemix2_compute_vdiff_vdiss_Eiw(              &
                  nlev     = nln-uln+1                               & 
                , dzw      = hnode(           uln:nln-1, node)       &
                , dt       = dt                                      &
                , c0       = iwe2_c0(         uln:nln  , node)       &
                , alpha_c  = iwe2_alpha_c(    uln:nln  , node)       &
                , fsrf     = iwe2_fsrf(                  node)       &
                , fbot     = iwe2_fbot(       uln:nln  , node)       &
                , Eiw_old  = iwe2_E_iw(iwe2_tau  , uln:nln  , node)  &
                , Eiw_new  = iwe2_E_iw(iwe2_taup1, uln:nln  , node)  &
                , Eiw_diss = iwe2_E_iw_diss(  uln:nln  , node)       &
                , Eiw_dt   = iwe2_E_iw_dt(    uln:nln  , node)       & ! optional: diagnostic
                , Eiw_vdif = iwe2_E_iw_vdif(  uln:nln  , node)       & ! optional: diagnostic
                , Eiw_srf  = iwe2_E_iw_fsrf(             node)       & ! optional: diagnostic
                , Eiw_bot  = iwe2_E_iw_fbot(  uln:nln  , node)       & ! optional: diagnostic
                )
        end do ! --> do node = 1, myDim_nod2D
        
        ! 8th. add lateral diffusion term (see. Olbers D., Eden C., 2013, A Global Model 
        ! for the Diapycnal Diffusivity Induced Internal Gravity Waves...)
        ! Eiw^(t+1) = Eiw^(t+1) + div_h( v_0 * tau_h * grad_h(v_0*E_iw^(t)) )
        if (idemix2_enable_hor_diffusion .or. idemix2_enable_hor_diff_iter) then 
            iwe2_E_iw_hdif = 0.0_WP
        end if 
        if (idemix2_enable_hor_diffusion) then 
            call compute_hdiff_Eiw(                 &  
                  iwe2_E_iw(iwe2_tau  , :, :)       &
                , iwe2_E_iw(iwe2_taup1, :, :)       &
                , iwe2_v0                           &
                , 1                                 &
                , partit                            &
                , mesh                              &
                , Eiw_dt   = iwe2_E_iw_dt(    :, :) & ! optional: diagnostic
                , Eiw_hdif = iwe2_E_iw_hdif(  :, :) & ! optional: diagnostic                   
                )
        end if
        
        
        
        !_______________________________________________________________________
        ! 9th. add tendency due to lateral diffusion with iterative method in case of 
        ! high resolution
        if (idemix2_enable_hor_diff_iter) then
            do iter=1, idemix2_hor_diff_niter
                call compute_hdiff_Eiw(                 &
                      iwe2_E_iw(iwe2_taup1, :, :)       &
                    , iwe2_E_iw(iwe2_taup1, :, :)       &
                    , iwe2_v0                           &
                    , idemix2_hor_diff_niter            &
                    , partit                            &
                    , mesh                              &
                    , Eiw_dt   = iwe2_E_iw_dt(    :, :) & ! optional: diagnostic
                    , Eiw_hdif = iwe2_E_iw_hdif(  :, :) & ! optional: diagnostic
                    )
            end do
        end if 
        
        
        
        !_______________________________________________________________________
        ! 10th. compute wave-wave interaction 
        if (idemix2_enable_M2 .or. idemix2_enable_niw) then 
            do node = 1, myDim_nod2D+eDim_nod2D
                uln = ulevels_nod2D(node)
                nln = nlevels_nod2D(node)
                !_______________________________________________________________
                if     (idemix2_enable_M2 .and. idemix2_enable_niw) then 
                    call cvmix_idemix2_compute_Eiw_waveinteract(                      &
                              nlev         = nln-uln+1                                & 
                            , nfbin        = idemix2_nfbin                            &
                            , dzw          = hnode(                uln:nln-1, node)   &
                            , dphi         = iwe2_dphit                               &
                            , dt           = dt                                       &
                            , E_iw_old     = iwe2_E_iw( iwe2_tau  ,  uln:nln, node)   &
                            , E_iw_new     = iwe2_E_iw( iwe2_taup1,  uln:nln, node)   &
                            , E_M2_old     = iwe2_E_M2( iwe2_tau  ,  :      , node)   &
                            , E_M2_new     = iwe2_E_M2( iwe2_taup1,  :      , node)   &
                            , E_M2_struct  = iwe2_E_M2_struct(       :      , node)   &
                            , alpha_M2_c   = iwe2_alpha_M2_c(                 node)   &
                            , tau_M2       = iwe2_M2_tau(                     node)   &
                            , E_niw_old    = iwe2_E_niw(iwe2_tau  ,  :      , node)   &
                            , E_niw_new    = iwe2_E_niw(iwe2_taup1,  :      , node)   &
                            , E_niw_struct = iwe2_E_niw_struct(      :      , node)   &
                            , tau_niw      = iwe2_niw_tau(                    node)   &
                            )
                elseif (idemix2_enable_M2) then 
                    call cvmix_idemix2_compute_Eiw_waveinteract(                      &
                              nlev         = nln-uln+1                                & 
                            , nfbin        = idemix2_nfbin                            &
                            , dzw          = hnode(                uln:nln-1, node)   &
                            , dphi         = iwe2_dphit                               &
                            , dt           = dt                                       &
                            , E_iw_old     = iwe2_E_iw( iwe2_tau  ,  uln:nln, node)   &
                            , E_iw_new     = iwe2_E_iw( iwe2_taup1,  uln:nln, node)   &
                            , E_M2_old     = iwe2_E_M2( iwe2_tau  ,  :      , node)   &
                            , E_M2_new     = iwe2_E_M2( iwe2_taup1,  :      , node)   &
                            , E_M2_struct  = iwe2_E_M2_struct(       :      , node)   &
                            , alpha_M2_c   = iwe2_alpha_M2_c(                 node)   &
                            , tau_M2       = iwe2_M2_tau(                     node)   &
                            )
                elseif (idemix2_enable_niw) then 
                    call cvmix_idemix2_compute_Eiw_waveinteract(                      &
                              nlev         = nln-uln+1                                & 
                            , nfbin        = idemix2_nfbin                            &
                            , dzw          = hnode(                uln:nln-1, node)   &
                            , dphi         = iwe2_dphit                               &
                            , dt           = dt                                       &
                            , E_iw_old     = iwe2_E_iw( iwe2_tau  ,  uln:nln, node)   &
                            , E_iw_new     = iwe2_E_iw( iwe2_taup1,  uln:nln, node)   &
                            , E_niw_old    = iwe2_E_niw(iwe2_tau  ,  :      , node)   &
                            , E_niw_new    = iwe2_E_niw(iwe2_taup1,  :      , node)   &
                            , E_niw_struct = iwe2_E_niw_struct(      :      , node)   &
                            , tau_niw      = iwe2_niw_tau(                    node)   &
                            )
                end if ! -->  (idemix2_enable_M2 .and. idemix2_enable_niw) then 
            end do ! --> for node = 1, myDim_nod2D
        end if ! --> if (idemix2_enable_M2 .or . idemix2_enable_niw) then 



        !_______________________________________________________________________
        ! 11th. write IDEMIX2 diffusivities and viscositie to FESOM only when IDEMIX2 is 
        ! used alone --> mostly for debuging --> otherwise TKE Av and Kv are use
        if(mix_scheme_nmb==7) then 
            
            !___________________________________________________________________
            ! write out diffusivity --> convert from elem to vertices
            do node=1, myDim_nod2D+eDim_nod2D
                uln = ulevels_nod2D(node)
                nln = nlevels_nod2D(node)
                !_______________________________________________________________
                ! convert Eiw disspation into Kv and Av on vertices 
                call cvmix_idemix2_compute_Eiw_diss2KvAv(         &
                          nlev    = nln-uln+1                     &
                        , Eiw_diss= iwe2_E_iw_diss(uln:nln, node) &
                        , Nsqr    = bvfreq(        uln:nln, node) &
                        , KappaH  = Kv(            uln:nln, node) &
                        , KappaM  = iwe2_Av(       uln:nln, node) &
                        )
            end do        
                
            !___________________________________________________________________
            ! more idemix2 viscosity from vertices to elements
            do elem=1,myDim_elem2D
                uln = ulevels(elem)
                nln = nlevels(elem)
                elnodes = elem2d_nodes(:,elem)
                do nz=uln, nln
                    Av(nz, elem) = sum(iwe2_Av(nz, elnodes))/3.0_WP
                end do
            end do
            call exchange_elem(Av, partit)
            call exchange_nod(Kv , partit)
        end if 
        
    end subroutine calc_cvmix_idemix2    

    
    
    !
    !
    !___________________________________________________________________________
    ! horizontal superbee advection of spectral bins
    subroutine adv_Ecompart_hor_spctrl_superbee(   &
                vel                         , & 
                ttf                         , &
                ttf_grad_n                  , &
                flux                        , &
                partit                      , &
                mesh                        , &
                flag_2ndord_time            , &
                flag_posdef                   &
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
        real(kind=WP)                           :: dx1, dy1, dx2, dy2, dxdy12(2), edlen, n_x, n_y, n_len  
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
!$OMP                                  dx1, dy1, dx2, dy2, dxdy12, dh1, dh2, n_x, n_y, n_len, edlen, &
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
            edlen = sqrt(dxdy12(1)**2 + dxdy12(2)**2)
            
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
                CFL   = min(1.0_WP, Ue*dt/edlen)
                
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
    ! cross spectral superbee advection across spectral bins
    subroutine adv_Ecompart_crss_spctrl_superbee( &
                cs                              , & ! cross-spectral exchange rate
                ttf                             , & ! Energy compartment @ node
                dphi                            , & ! width spectral bin
                flux                            , & ! cross-spectral flux
                partit                          , &
                mesh                            , &
                flag_2ndord_time                , &
                flag_posdef                       &
                ! flag_wlimitcfl                    &
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
        integer                              :: node, fbini, nfbin, idxp1, idxp2, idxm1
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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, fbini, idxp1, idxp2, idxm1, &
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
                !_______________________________________________________________
                ! Periodic wrapping for gradient stencil
                idxp1   = fbini+1 ; if (idxp1>nfbin) idxp1=idxp1-(nfbin-2)
                idxp2   = fbini+2 ; if (idxp2>nfbin) idxp2=idxp2-(nfbin-2)
                idxm1   = fbini-1 ; if (idxm1<1    ) idxm1=nfbin-2
                
                ! ! limit W cross spectral exchange rate to the CFL limit 
                ! cs(fbini, node) = merge(  min(abs(cs(fbini, node)), cflmax*dphi(fbini)/dt)*sign(1.0_WP, cs(fbini, node)) &
                !                         , cs(fbini, node) &
                !                         , flag_wlimitcfl &
                !                        )
                
                CFL     = abs(cs(fbini, node)) * dt / dphi(fbini)
                
                !_______________________________________________________________
                ! compute tracer difference 
                ttf0    = ttf(fbini  , node)
                ttfp1   = ttf(idxp1  , node)
                dttf0p1 = ttfp1 - ttf0
                
                !_______________________________________________________________
                ! tracer Slope Ratio Calculation for upwind point
                ! compute tracer slope 
                R      = (ttfp1-ttf(idxp2,node))/(-dttf0p1+small)
                
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
                R      = (ttf0-ttf(idxm1,node))/(dttf0p1+small)
                
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
    subroutine adv_Ecompart_flx2tra_spctrl(    &
                flx_h                   , &
                flx_v                   , &
                div_h                   , &
                div_v                   , &
                dphit                   , &
                vol_s                   , &
                partit                  , &
                mesh                      &
                    )
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in   ), target :: mesh
        real(kind=WP) , intent(in)            :: flx_h( idemix2_nfbin, partit%myDim_edge2D) 
        real(kind=WP) , intent(in)            :: flx_v( idemix2_nfbin, partit%myDim_nod2D) 
        real(kind=WP) , intent(inout)         :: div_h( idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP) , intent(inout)         :: div_v( idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP) , intent(in   )         :: dphit( idemix2_nfbin)
        real(kind=WP) , intent(in   )         :: vol_s(partit%myDim_nod2D+partit%eDim_nod2D)
        !___LOCAL VARIABLES_____________________________________________________
        integer                               :: node, edge, fbini, nfbin, ednodes(2)
        real(kind=WP)                         :: inv_dphi(idemix2_nfbin)
        
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
        
        ! normalize with total volume of scalar cell from top to bottom
#ifndef ENABLE_OPENACC
        !$OMP PARALLEL DO COLLAPSE(2)
#else
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
        do node=1, myDim_nod2d+eDim_nod2D
            do fbini=2,nfbin-1
                div_h(fbini, node) = div_h(fbini,node)/vol_s(node)
            end do
        end do
#ifndef ENABLE_OPENACC
        !$OMP END PARALLEL DO
#else
        !$ACC END PARALLEL LOOP
#endif

    end subroutine adv_Ecompart_flx2tra_spctrl
    
    
    
    !
    !
    !___________________________________________________________________________
    ! integrate equation for wave energy compartments
    subroutine hsintegrate_Ecompart(      & 
                tim1, ti, tip1            &
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
                , vol_s                   &
                , partit                  &
                , mesh                    &
                , flag_AB2                & ! to 2nd order Adams-Bashfort on top 
                , Edt                     &
                , Eadvh                   &
                , Eadvs                   &
                , Ediss                   &
                , Eforc                   &
                )
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit)  , intent(inout), target  :: partit
        type(t_mesh)    , intent(in   ), target  :: mesh
        integer         , intent(in   )          :: ti, tim1, tip1
        character(len=*), intent(in   )          :: Ename
        real(kind=WP)   , intent(inout)          :: E(       3, idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(inout)          :: Edivh(   3, idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(inout)          :: Edivs(   3, idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(in   )          :: vel(     2, idemix2_nfbin, partit%myDim_elem2D+partit%eDim_elem2D)
        real(kind=WP)   , intent(inout)          :: cs(         idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(in   )          :: forc(       idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(in   )          :: tauE(                      partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(inout)          :: gradxy_e(2, idemix2_nfbin, partit%myDim_elem2D+partit%eDim_elem2D)
        real(kind=WP)   , intent(inout)          :: gradxy_n(2, idemix2_nfbin, partit%myDim_nod2D +partit%eDim_nod2D)
        real(kind=WP)   , intent(in   )          :: dphit(      idemix2_nfbin)
        real(kind=WP)   , intent(inout)          :: flx_uv(     idemix2_nfbin, partit%myDim_edge2D)
        real(kind=WP)   , intent(inout)          :: flx_cs(     idemix2_nfbin, partit%myDim_nod2D ) 
        real(kind=WP)   , intent(in   )          :: vol_s(                     partit%myDim_nod2D +partit%eDim_nod2D)
        logical         , intent(in)             :: flag_AB2
        
        
        real(kind=WP)   , intent(inout), optional:: Edt(        idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP)   , intent(inout), optional:: Eadvh(      idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP)   , intent(inout), optional:: Eadvs(      idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP)   , intent(inout), optional:: Ediss(      idemix2_nfbin, partit%myDim_nod2D)
        real(kind=WP)   , intent(inout), optional:: Eforc(      idemix2_nfbin, partit%myDim_nod2D)
        
        !___LOCAL VARIABLES_____________________________________________________
        integer                                  :: node, fbini, nfbin
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"        
        !_______________________________________________________________________
        nfbin=idemix2_nfbin
        iwe2_flx_uv(:,:)     = 0.0
        iwe2_flx_w( :,:)     = 0.0
        
        !_______________________________________________________________________
        ! Exchange all 3 time levels (minimal overhead, avoids array slicing issues)
        call exchange_nod_fbin(E, partit)
        
        !_______________________________________________________________________
        ! compute gradient of iwe2_E on nodes
        call tracer_gradient_elements(E(ti,:,:), gradxy_e, partit, mesh, .False.)
        call exchange_elem_fbin(gradxy_e, partit)
        call interp_e2n(gradxy_e, gradxy_n, mesh, partit, .False.)
        call exchange_nod_fbin(gradxy_n, partit)
        
        !_______________________________________________________________________
        ! compute horizontal superbee advected tracer flux 
        call adv_Ecompart_hor_spctrl_superbee(            &
                                          vel             &
                                        , E(ti,:,:)       &
                                        , gradxy_n        &
                                        , flx_uv          &
                                        , partit          &
                                        , mesh            &
                                        , .True.          & ! do 2nd order in space & time
                                        , .True.          & ! enforce positive definit 
                                        )
                                        
        ! compute cross spectral superbee advected tracer flux                                 
        call adv_Ecompart_crss_spctrl_superbee(           &  
                                          cs              &
                                        , E(ti,:,:)       & 
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
                                        ,Edivh(ti,:,:)    &
                                        ,Edivs(ti,:,:)    &
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
                do fbini=2,nfbin-1
                    E(tip1, fbini, node) = ( E(ti, fbini, node) &
                                            + dt * ( &
                                                    + (1.5_WP+idemix2_AB_epsilon)*Edivh(ti  , fbini, node) &
                                                    - (0.5_WP+idemix2_AB_epsilon)*Edivh(tim1, fbini, node) &
                                                    + (1.5_WP+idemix2_AB_epsilon)*Edivs(ti  , fbini, node) &
                                                    - (0.5_WP+idemix2_AB_epsilon)*Edivs(tim1, fbini, node) &
                                                    +                              forc(      fbini, node) &
                                                ) &
                                            ) / (1.0_WP + dt*tauE(node))  
                end do
                E(tip1,     1, node) = E(tip1, nfbin-1, node)
                E(tip1, nfbin, node) = E(tip1,       2, node)
            end do    
        !_______________________________________________________________________
        ! Forward Euler + implicit dissipation
        else
            do node= 1, myDim_nod2d
                do fbini=2,nfbin-1
                    E(tip1, fbini, node) = ( E(ti, fbini, node) &
                                            + dt * ( &
                                                    + Edivh(ti, fbini, node) & ! div(c_uv*E^n)
                                                    + Edivs(ti, fbini, node) & ! d/dphi (w*E^n)
                                                    + forc(     fbini, node) & ! forc
                                                   ) &
                                            ) / (1.0_WP + dt*tauE(node)) 
                end do
                E(tip1,     1, node) = E(tip1, nfbin-1, node)
                E(tip1, nfbin, node) = E(tip1,       2, node)
            end do
        end if     
        
        !_______________________________________________________________________
        ! additional diagnostics
        if (present(Edt))   Edt(:,:)   =  (E(tip1, :, 1:myDim_nod2D)-E(ti, :, 1:myDim_nod2D))
        
        if (present(Eadvh)) then 
            if (flag_AB2) then    
                Eadvh(:,:) =  ((1.5+idemix2_AB_epsilon)*Edivh(ti  , :, 1:myDim_nod2D) &
                             - (0.5+idemix2_AB_epsilon)*Edivh(tim1, :, 1:myDim_nod2D))
            else                             
                Eadvh(:,:) = Edivh(ti, :, 1:myDim_nod2D)                             
            end if 
        end if 
        
        if (present(Eadvs)) then 
            if (flag_AB2) then    
                Eadvs(:,:) =  ((1.5+idemix2_AB_epsilon)*Edivs(ti  , :, 1:myDim_nod2D) &
                             - (0.5+idemix2_AB_epsilon)*Edivs(tim1, :, 1:myDim_nod2D))
            else                 
                Eadvs(:,:) = Edivs(ti, :, 1:myDim_nod2D)
            end if
        end if     
        
        if (present(Ediss)) then
            do node= 1, myDim_nod2d
                Ediss(:,node) =  - tauE(node)*E(tip1, :, node)                                      
            end do
        end if     
        
        if (present(Eforc)) Eforc(:,:) = forc(:, 1:myDim_nod2D)
        
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
        real(kind=WP)                           :: lcl_sumvol, glb_sumvol, vol
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"   
        vol_nodB2T = 0.0_WP
        lcl_sumvol = 0.0_WP
        do node = 1, myDim_nod2D
            nln = nlevels_nod2D(node)-1
            uln = ulevels_nod2D(node)
            do nz=uln,nln
                lcl_sumvol       = lcl_sumvol       + areasvol(nz,node)*abs(zbar(nz)-zbar(nz+1))
                vol_nodB2T(node) = vol_nodB2T(node) + areasvol(nz,node)*abs(zbar(nz)-zbar(nz+1))
            end do
        end do !-->do node = 1,node_size
        call MPI_Allreduce(lcl_sumvol, glb_sumvol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        if (partit%mype == 0) then
            write(*,*) ' debug: vol_node_fix: ', glb_sumvol
        end if
        call exchange_nod(vol_nodB2T, partit)
        
        lcl_sumvol = 0.0_WP
        do elem = 1, myDim_elem2D
            ! Nodes are uniquely partitioned (e.g., via METIS) — each node belongs 
            ! to exactly one rank Elements are then distributed so that each rank 
            ! has all elements touching its owned nodes — but boundary elements 
            ! necessarily appear on multiple ranks (they share vertices across partitions)
            elnodes=mesh%elem2D_nodes(:,elem)
            if (elnodes(1) > myDim_nod2D) cycle
            nln = nlevels(elem)-1
            uln = ulevels(elem)
            do nz=uln,nln
                lcl_sumvol       = lcl_sumvol       + elem_area(elem)*abs(zbar(nz)-zbar(nz+1))
            end do
        end do !-->do node = 1,node_size
        call MPI_Allreduce(lcl_sumvol, glb_sumvol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
        if (partit%mype == 0) then
            write(*,*) ' debug: vol_elem_fix: ', glb_sumvol
        end if
        
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
                mesh                    , &
                Eiw_dt                  , &
                Eiw_hdif                  &
                )
        implicit none
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_partit)  , intent(inout), target :: partit
        type(t_mesh)    , intent(in   ), target :: mesh
        real(kind=WP)   , intent(in   )         :: Eiw_old(:,:)
        real(kind=WP)   , intent(inout)         :: Eiw(:,:)
        real(kind=WP)   , intent(in   )         :: v0(:,:)
        integer         , intent(in   )         :: n_hor_iter
        real(kind=WP)   , intent(inout),optional:: Eiw_dt(:,:)
        real(kind=WP)   , intent(inout),optional:: Eiw_hdif(:,:)
        
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
        
        !_______________________________________________________________________
        ! add to total change in Eiw 
        do node = 1, myDim_nod2D
            nu1 = ulevels_nod2D(node)
            nl1 = nlevels_nod2D(node)
            do nz = nu1, nl1
                Eiw_hdif(nz, node) = Eiw_hdif(nz, node) + (Eiw(nz, node)-Eiw_old(nz, node))/dt
                Eiw_dt(nz, node)   = Eiw_dt(nz, node) + (Eiw(nz, node)-Eiw_old(nz, node))/dt
            end do        
        end do
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
        real(kind=WP),  intent(in)            :: div_h(3, idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP),  intent(in)            :: div_s(3, idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D)
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
                lcl_divh_sum = lcl_divh_sum + div_h(ti, fbini, node)*vol_s(node)
            end do
        end do
!         write(*,*) lcl_divh_sum, mype
        
        lcl_divs_sum = 0.0_WP
        do node = 1, myDim_nod2D  ! Full extended domain: telescoping holds locally per rank
            do fbini = 2, idemix2_nfbin-1  ! Exclude ghost bins
                lcl_divs_sum = lcl_divs_sum + div_s(ti, fbini, node)*dphit(fbini) 
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
        real(kind=WP),  intent(in)            :: Ecompart(3, idemix2_nfbin, partit%myDim_nod2D+partit%eDim_nod2D)
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
                lcl_energy = lcl_energy + Ecompart(ti, fbini, node) * dphit(fbini) * vol_s(node)
                lcl_forc   = lcl_forc   + forc(fbini, node)         * dphit(fbini) * vol_s(node)
                lcl_diss   = lcl_diss   + tauE(node) * Ecompart(ti, fbini, node) * dphit(fbini) * vol_s(node)
                lcl_maxE   = max(lcl_maxE, Ecompart(ti, fbini, node))
                lcl_minE   = min(lcl_minE, Ecompart(ti, fbini, node))
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
                             ', global forc/diss = ', glb_forc/glb_diss, ' W'
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


end module g_cvmix_idemix2

