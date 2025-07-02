!=======================================================================
!
! This module defines and and initializes the namelists
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
! Adapted for Icepack 1.5.0 by F. Kauker (frank.kauker@awi.de)
!
!=======================================================================

submodule (icedrv_main) icedrv_set

  use icepack_intfc, only: &
       icepack_init_parameters, icepack_init_tracer_flags, icepack_init_tracer_sizes,     &
       icepack_init_tracer_indices, icepack_query_parameters, icepack_query_tracer_flags, &
       icepack_query_tracer_sizes, icepack_query_tracer_indices, icepack_warnings_flush,  &
       icepack_warnings_aborted
  use icedrv_system, only: icedrv_system_abort, icedrv_system_init
      
contains 

  module subroutine set_icepack(flag_debug, ice, partit)
    !
    ! originally 'subroutine input_data' in icedrv_init.F90
    !
    use MOD_ICE
    
    implicit none

    logical (kind=log_kind), intent(in) :: flag_debug
    type(t_partit), intent(inout), target :: partit
    type(t_ice), intent(inout), target :: ice

    ! local variables
    character(len=char_len) :: nml_filename, diag_filename
    character(len=*), parameter :: subname = '(set_icepack)'
    logical(kind=log_kind) :: tr_pond, wave_spec
    integer(kind=int_kind) :: &
         nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd,  &
         nt_iso, nt_isosno, nt_isoice, nt_aero, nt_fsd, nt_FY, ntrcr, nt_iage, nt_smice,   &
         nt_smliq, nt_rhos, nt_rsnw, nml_error, diag_error, mpi_error, n                      
    real(kind=dbl_kind) :: rplvl, rptopo, puny

    ! env namelist

    logical(kind=log_kind) :: conserv_check   ! run conservation checks and abort if checks fail
                                              ! originally in setup_nml (frank.kauker@awi.de)
    integer(kind=int_kind) :: trbgcz          ! set to 1 for zbgc tracers (needs TRBGCS = 0 and TRBRI = 1)
    integer(kind=int_kind) :: trbri           ! set to 1 for brine height tracer
    integer(kind=int_kind) :: trage           ! set to 1 for ice age tracer
    integer(kind=int_kind) :: trfy            ! set to 1 for first-year ice area tracer
    integer(kind=int_kind) :: trlvl           ! set to 1 for level and deformed ice tracers
    integer(kind=int_kind) :: trpnd           ! set to 1 for melt pond tracers
    integer(kind=int_kind) :: trsnow          ! set to 1 for advanced snow physics tracers 
    integer(kind=int_kind) :: trbgcs          ! set to 1 for skeletal layer tracers (needs TRBGCZ = 0)

    ! tracer namelist
    
    logical (kind=log_kind)  :: tr_iage
    logical (kind=log_kind)  :: tr_FY
    logical (kind=log_kind)  :: tr_lvl
    logical (kind=log_kind)  :: tr_pond_topo
    logical (kind=log_kind)  :: tr_pond_lvl
    logical (kind=log_kind)  :: tr_snow
    logical (kind=log_kind)  :: tr_aero
    logical (kind=log_kind)  :: tr_fsd
    logical (kind=log_kind)  :: tr_iso

    ! grid namelist
 
    integer (kind=int_kind)  :: kcatbound

    ! thermo namelist
    
    character (len=char_len) :: conduct
    integer (kind=int_kind)  :: kitd
    integer (kind=int_kind)  :: ktherm
    real (kind=dbl_kind)     :: a_rapid_mode
    real (kind=dbl_kind)     :: Rac_rapid_mode  
    real (kind=dbl_kind)     :: aspect_rapid_mode 
    real (kind=dbl_kind)     :: dSdt_slow_mode    
    real (kind=dbl_kind)     :: phi_c_slow_mode   
    real (kind=dbl_kind)     :: phi_i_mushy       
    real (kind=dbl_kind)     :: hi_min
    real (kind=dbl_kind)     :: ksno
    real (kind=dbl_kind)     :: Tliquidus_max
    real (kind=dbl_kind)     :: floediam
    real (kind=dbl_kind)     :: hfrazilmin

    ! dynamics namelist
    
    integer (kind=int_kind)  :: kstrength       
    integer (kind=int_kind)  :: krdg_partic     
    integer (kind=int_kind)  :: krdg_redist     
    real (kind=dbl_kind)     :: mu_rdg          
    real (kind=dbl_kind)     :: Cf

    ! shortwave namelist

    logical (kind=log_kind)  :: sw_redist  
    character (len=char_len) :: shortwave      
    character (len=char_len) :: albedo_type
    character (len=char_len) :: snw_ssp_table
    real (kind=dbl_kind)     :: albicev        
    real (kind=dbl_kind)     :: albicei         
    real (kind=dbl_kind)     :: albsnowv        
    real (kind=dbl_kind)     :: albsnowi      
    real (kind=dbl_kind)     :: albocn      
    real (kind=dbl_kind)     :: ahmax         
    real (kind=dbl_kind)     :: R_ice         
    real (kind=dbl_kind)     :: R_pnd         
    real (kind=dbl_kind)     :: R_snw         
    real (kind=dbl_kind)     :: dT_mlt        
    real (kind=dbl_kind)     :: rsnw_mlt       
    real (kind=dbl_kind)     :: kalg           
    real (kind=dbl_kind)     :: sw_frac
    real (kind=dbl_kind)     :: sw_dtemp

    ! ponds namelist
    
    character (len=char_len) :: frzpnd
    real (kind=dbl_kind)     :: hp1 
    real (kind=dbl_kind)     :: hs0 
    real (kind=dbl_kind)     :: hs1
    real (kind=dbl_kind)     :: dpscale         
    real (kind=dbl_kind)     :: rfracmin      
    real (kind=dbl_kind)     :: rfracmax     
    real (kind=dbl_kind)     :: pndaspect

     ! forcing namelist
    
    logical (kind=log_kind)  :: formdrag
    logical (kind=log_kind)  :: calc_strair     
    logical (kind=log_kind)  :: calc_Tsfc
    logical (kind=log_kind)  :: calc_dragio
    logical (kind=log_kind)  :: highfreq        
    logical (kind=log_kind)  :: update_ocn_f
    logical (kind=log_kind)  :: l_mpond_fresh
    logical (kind=log_kind)  :: oceanmixed_ice
    character (len=char_len) :: atmbndy
    character (len=char_len) :: fbot_xfer_type
    character (len=char_len) :: cpl_frazil
    character (len=char_len) :: congel_freeze
    character (len=char_len) :: tfrz_option
    character (len=char_len) :: saltflux_option
    character (len=char_len) :: wave_spec_type
    integer (kind=int_kind)  :: natmiter        
    real (kind=dbl_kind)     :: ustar_min      
    real (kind=dbl_kind)     :: emissivity      
    real (kind=dbl_kind)     :: dragio      
    real (kind=dbl_kind)     :: atmiter_conv ! Flux convergence tolerance
    real (kind=dbl_kind)     :: ice_ref_salinity ! Ice reference salinity for fluxes

    ! to snychronize density definition between icepack and fesom2, becomes
    ! crucial for waterflux computation and volume conservation
    real (kind=dbl_kind)     :: rhoice, rhosno, rhowat, rhofwt ! why here? frank.kauker@awi.de
    
    ! snow namelist
    
    logical (kind=log_kind)  :: use_smliq_pnd
    logical (kind=log_kind)  :: snwgrain
    character (len=char_len) :: snwredist
    character (len=char_len) :: snw_aging_table 
    real (kind=dbl_kind)     :: rsnw_fall
    real (kind=dbl_kind)     :: rsnw_tmax
    real (kind=dbl_kind)     :: rhosnew
    real (kind=dbl_kind)     :: rhosmin
    real (kind=dbl_kind)     :: rhosmax
    real (kind=dbl_kind)     :: windmin
    real (kind=dbl_kind)     :: drhosdwind
    real (kind=dbl_kind)     :: snwlvlfac

    ! zbgc namelist

    logical (kind=log_kind)  :: tr_brine
    logical (kind=log_kind)  :: tr_zaero
    logical (kind=log_kind)  :: modal_aero
    logical (kind=log_kind)  :: skl_bgc
    logical (kind=log_kind)  :: z_tracers
    logical (kind=log_kind)  :: dEdd_algae
    logical (kind=log_kind)  :: solve_zbgc
    logical (kind=log_kind)  :: restore_bgc
    logical (kind=log_kind)  :: scale_bgc
    logical (kind=log_kind)  :: solve_zsal
    logical (kind=log_kind)  :: tr_bgc_Nit
    logical (kind=log_kind)  :: tr_bgc_C
    logical (kind=log_kind)  :: tr_bgc_chl
    logical (kind=log_kind)  :: tr_bgc_Am
    logical (kind=log_kind)  :: tr_bgc_Sil
    logical (kind=log_kind)  :: tr_bgc_DMS
    logical (kind=log_kind)  :: tr_bgc_PON
    logical (kind=log_kind)  :: tr_bgc_hum
    logical (kind=log_kind)  :: tr_bgc_DON
    logical (kind=log_kind)  :: tr_bgc_Fe
    character (char_len)     :: bgc_flux_type

    !-----------------------------------------------------------------
    ! Namelist variables
    !-----------------------------------------------------------------

    ! env_nml is not used in 1.5.0 because it causes a lot of trouble (see definition of n_bgc, nltrcr, max_nsw, max_ntrcr)
    ! All variables have to be pre-defined in this file and cannot be changed after compilation (fkauker@awi.de)

    namelist / grid_nml /                                               &
         kcatbound

    namelist /thermo_nml/ &
         kitd,           ktherm,          ksno,     conduct,             &
         a_rapid_mode,   Rac_rapid_mode,  aspect_rapid_mode,             &
         dSdt_slow_mode, phi_c_slow_mode, phi_i_mushy,                   &
         floediam,       hfrazilmin,      Tliquidus_max,    hi_min

    namelist /dynamics_nml/ &
         kstrength,      krdg_partic,    krdg_redist,    mu_rdg,         &
         Cf

    namelist /shortwave_nml/ &
         shortwave,      albedo_type,     albocn,                        &
         albicev,        albicei,         albsnowv,      albsnowi,       &
         ahmax,          R_ice,           R_pnd,         R_snw,          &
         sw_redist,      sw_frac,         sw_dtemp,                      &
         dT_mlt,         rsnw_mlt,        kalg,          snw_ssp_table

    namelist /ponds_nml/ &
         hs0,            dpscale,         frzpnd,                        &
         rfracmin,       rfracmax,        pndaspect,     hs1,            &
         hp1
    
    namelist /snow_nml/ &
         snwredist,      snwgrain,       rsnw_fall,     rsnw_tmax,      &
         rhosnew,        rhosmin,        rhosmax,       snwlvlfac,      &
         windmin,        drhosdwind,     use_smliq_pnd, snw_aging_table
    
    namelist /forcing_nml/ &
         atmbndy,         calc_strair,     calc_Tsfc,       &
         update_ocn_f,    l_mpond_fresh,   ustar_min,       &
         fbot_xfer_type,  oceanmixed_ice,  emissivity,      &
         formdrag,        highfreq,        natmiter,        &
         atmiter_conv,    calc_dragio,     congel_freeze,   &
         tfrz_option,     saltflux_option, ice_ref_salinity, &
         wave_spec_type,  cpl_frazil

    nml_filename = 'namelist.icepack'        ! name of icepack namelist file

    !-----------------------------------------------------------------
    ! Derived types used by icepack
    !-----------------------------------------------------------------
    call icedrv_system_init(partit)
    p_partit => partit
    nx         = p_partit%myDim_nod2D  + p_partit%eDim_nod2D
    nx_elem    = p_partit%myDim_elem2D + p_partit%eDim_elem2D
    nx_nh      = p_partit%myDim_nod2D
    nx_elem_nh = p_partit%myDim_elem2D
    mype       = p_partit%mype

    !-----------------------------------------------------------------
    ! !!! THE VARIABLES below are originally defined in
    ! icedrv_domain_size.F90 with the help of environment variables.
    ! Because a consistency check follows afterwards these variables
    ! cannot be changed in namelist.icepack when Icepack is coupled to
    ! FESOM2 (frank.kauker@awi.de)
    !----------------------------------------------------------------- 
    ncat      = 5          ! number of categories
    nfsd      = 1          ! number of floe size categories
    nilyr     = 4          ! number of ice layers per category
    nslyr     = 4          ! number of snow layers per category
    n_iso     = 0          ! number of  water isotopes (snow, sea ice),
                           ! name changed in this routine because because n_iso is used in fesom as well 
    n_aero    = 0          ! number of aerosols in use
    n_zaero   = 0          ! number of z aerosols in use
    n_algae   = 0          ! number of algae in use
    n_doc     = 0          ! number of DOC pools in use
    n_dic     = 0          ! number of DIC pools in use
    n_don     = 0          ! number of DON pools in use
    n_fed     = 0          ! number of Fe  pools in use dissolved Fe
    n_fep     = 0          ! number of Fe  pools in use particulate Fe
    nfreq     = 25         ! number of wave frequencies ! HARDWIRED FOR NOW
    nblyr     = 4          ! number of bio/brine layers per category
                           ! maximum number of biology tracers + aerosols
    ndtd      = 1          ! dynamic time steps per thermodynamic time step
    conserv_check = .false.! if .true., run conservation checks and abort if checks fail
                           ! originally in setup_nml (frank.kauker@awi.de)
                           ! *** add to kscavz in icepack_zbgc_shared.F90
    trbgcz    = 0          ! set to 1 for zbgc tracers (needs trbgcs = 0 and trbri = 1)
    trbri     = 0          ! set to 1 for brine height tracer
    trage     = 0          ! set to 1 for ice age tracer
    trfy      = 0          ! set to 1 for first-year ice area tracer
    trlvl     = 1          ! set to 1 for level and deformed ice tracers
    trpnd     = 0          ! set to 1 for melt pond tracers
    trbgcs    = 0          ! set to 1 for skeletal layer tracers (needs

    n_bgc     = (n_algae*2 + n_doc + n_dic + n_don + n_fed + n_fep +n_zaero &
         + 8)         ! nit, am, sil, dmspp, dmspd, dms, pon, humic
    nltrcr    = (n_bgc*trbgcz)*trbri ! number of zbgc (includes zaero)
                                     ! and zsalinity tracers
    max_nsw   = (nilyr+nslyr+2) & ! total chlorophyll plus aerosols
         * (n_zaero)     ! number of tracers active in shortwave calculation
    max_ntrcr =   1    & ! 1 = surface temperature
         + nilyr       & ! ice salinity
         + nilyr       & ! ice enthalpy
         + nslyr       & ! snow enthalpy
                         ! optional tracers:
         + nfsd        & ! number of floe size categories
         + trage       & ! age
         + trfy        & ! first-year area
         + trlvl*2     & ! level/deformed ice
         + trpnd*3     & ! ponds
         + trsnow*4*nslyr & ! snow redistribution/metamorphism
         + n_iso*2     & ! number of isotopes (in ice and snow)
         + n_aero*4    & ! number of aerosols * 4 aero layers
         + trbri       & ! brine height
         + trbgcs*n_bgc                 & ! skeletal layer BGC
         + n_bgc*trbgcz*trbri*(nblyr+3) & ! zbgc (off if TRBRI=0)
         + n_bgc*trbgcz                 & ! mobile/stationary phase tracer
         + 1             ! for unused tracer flags
    
    !-----------------------------------------------------------------
    ! query Icepack default values
    !-----------------------------------------------------------------

    call icepack_query_parameters(ustar_min_out=ustar_min, Cf_out=Cf, &
         albicev_out=albicev, albicei_out=albicei, ksno_out = ksno,   &
         albsnowv_out=albsnowv, albsnowi_out=albsnowi, hi_min_out=hi_min, &
         natmiter_out=natmiter, ahmax_out=ahmax, shortwave_out=shortwave, &
         atmiter_conv_out = atmiter_conv, calc_dragio_out=calc_dragio, &
         albedo_type_out=albedo_type, R_ice_out=R_ice, R_pnd_out=R_pnd, &
         R_snw_out=R_snw, dT_mlt_out=dT_mlt, rsnw_mlt_out=rsnw_mlt, &
         kstrength_out=kstrength, krdg_partic_out=krdg_partic, &
         krdg_redist_out=krdg_redist, mu_rdg_out=mu_rdg, &
         atmbndy_out=atmbndy, calc_strair_out=calc_strair, &
         formdrag_out=formdrag, highfreq_out=highfreq, &
         emissivity_out=emissivity, &
         kitd_out=kitd, kcatbound_out=kcatbound, hs0_out=hs0, &
         dpscale_out=dpscale, frzpnd_out=frzpnd, &
         rfracmin_out=rfracmin, rfracmax_out=rfracmax, &
         pndaspect_out=pndaspect, hs1_out=hs1, hp1_out=hp1, &
         ktherm_out=ktherm, calc_Tsfc_out=calc_Tsfc, &
         floediam_out=floediam, hfrazilmin_out=hfrazilmin, &
         update_ocn_f_out = update_ocn_f, cpl_frazil_out = cpl_frazil, &
         conduct_out=conduct, a_rapid_mode_out=a_rapid_mode, &
         Rac_rapid_mode_out=Rac_rapid_mode, &
         aspect_rapid_mode_out=aspect_rapid_mode, &
         dSdt_slow_mode_out=dSdt_slow_mode, &
         phi_c_slow_mode_out=phi_c_slow_mode, Tliquidus_max_out=Tliquidus_max, &
         phi_i_mushy_out=phi_i_mushy, conserv_check_out=conserv_check, &
         tfrz_option_out=tfrz_option, saltflux_option_out=saltflux_option, &
         ice_ref_salinity_out=ice_ref_salinity, kalg_out=kalg, &
         fbot_xfer_type_out=fbot_xfer_type, puny_out=puny, &
         congel_freeze_out=congel_freeze, &
         wave_spec_type_out=wave_spec_type, dragio_out=dragio, &
         sw_redist_out=sw_redist, sw_frac_out=sw_frac, sw_dtemp_out=sw_dtemp, &
         snwredist_out=snwredist, use_smliq_pnd_out=use_smliq_pnd, &
         snwgrain_out=snwgrain, rsnw_fall_out=rsnw_fall, rsnw_tmax_out=rsnw_tmax, &
         rhosnew_out=rhosnew, rhosmin_out = rhosmin, rhosmax_out=rhosmax, &
         windmin_out=windmin, drhosdwind_out=drhosdwind, snwlvlfac_out=snwlvlfac, &
         snw_aging_table_out=snw_aging_table)

    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)
          
    !-----------------------------------------------------------------
    ! Synchronize values between Icepack and FESOM2 values
    !-----------------------------------------------------------------
    ! Make the namelists.ice and namelist.icepack consistent (icepack wins
    ! over fesom)
    ice%cd_oce_ice = dragio
    ice%thermo%albw = albocn
    ! in terms of density definition fesom wins over icepack, otherwise
    ! we can get in trouble with the waterflux computation and thus the 
    ! volume conservation under zstar
    !!! make sure that this values are not overwritten in 'namelist.icepack' !!! frank.kauker@awi.de
    rhoice = ice%thermo%rhoice 
    rhosno = ice%thermo%rhosno
    rhowat = ice%thermo%rhowat
    rhofwt = ice%thermo%rhofwt
    
    !-----------------------------------------------------------------
    ! other default values 
    !-----------------------------------------------------------------

    l_mpond_fresh = .false.     ! logical switch for including meltpond freshwater
                                ! flux feedback to ocean model
    oceanmixed_ice  = .false.   ! if true, use internal ocean mixed layer
    wave_spec_type  = 'none'    ! type of wave spectrum forcing

    ! extra tracers - have to be consistent to integer values in 'env namelist - STANDARD VALUES' (fkauker@awi.de]
    tr_iage         = .false.   ! chronological ice age
    tr_FY           = .false.   ! area of first year ice
    tr_lvl          = .true.    ! level ice 
    tr_pond_topo    = .false.   ! explicit melt ponds (topographic)
    tr_pond_lvl     = .false.   ! level-ice melt ponds
    tr_snow         = .false.   ! snow redistribution or metamorphosis tracers
    tr_aero         = .false.   ! aerosols
    tr_fsd          = .false.   ! floe size distribution
    tr_iso          = .false.   ! water isotope tracers
    
    !-----------------------------------------------------------------
    ! open error and diagnostic file
    !-----------------------------------------------------------------

    if (partit%mype == 0) write(*,*) 'icedrv_set: Error output will be in file '
    if (partit%mype == 0) write(*,*) '       icepack.errors'

    diag_filename = 'icepack.errors'
    open (ice_stderr, file=diag_filename, status='unknown', iostat=diag_error)
    if (diag_error /= 0) then
       if (partit%mype == 0) write(nu_diag,*) subname,' Error while opening error file'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif

    if (partit%mype == 0) write(*,*) 'icedrv_set: Diagnostic output will be in file '
    if (partit%mype == 0) write(*,*) '       icepack.diagnostics'

    diag_filename = 'icepack.diagnostics'
    open (nu_diag, file=diag_filename, status='unknown', iostat=diag_error)
    if (diag_error /= 0) then
       if (partit%mype == 0) write(nu_diag,*) subname,' Error while opening diagnostic file'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif

    !-----------------------------------------------------------------
    ! read from input file
    !-----------------------------------------------------------------
    
    open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: open file '// &
            trim(nml_filename), &
            file=__FILE__, line=__LINE__)
    endif

    rewind(unit=nu_nml, iostat=nml_error)
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: grid_nml rewind ', &
            file=__FILE__, line=__LINE__)
    endif
    nml_error =  1
    if (partit%mype == 0) write(nu_diag,*) subname,' Reading grid_nml'
    do while (nml_error > 0)
       read(nu_nml, nml=grid_nml,iostat=nml_error)
    end do
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: grid_nml reading ', &
            file=__FILE__, line=__LINE__)
    endif

    rewind(unit=nu_nml, iostat=nml_error)
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: thermo_nml rewind ', &
            file=__FILE__, line=__LINE__)
    endif
    nml_error =  1
    if (partit%mype == 0) write(nu_diag,*) subname,' Reading thermo_nml'
    do while (nml_error > 0)
       read(nu_nml, nml=thermo_nml,iostat=nml_error)
    end do
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: thermo_nml reading ', &
            file=__FILE__, line=__LINE__)
    endif

    rewind(unit=nu_nml, iostat=nml_error)
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: shortwave_nml rewind ', &
            file=__FILE__, line=__LINE__)
    endif
    nml_error =  1
    if (partit%mype == 0) write(nu_diag,*) subname,' Reading shortwave_nml'
    do while (nml_error > 0)
       read(nu_nml, nml=shortwave_nml,iostat=nml_error)
    end do
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: shortwave_nml reading ', &
            file=__FILE__, line=__LINE__)
    endif

    rewind(unit=nu_nml, iostat=nml_error)
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: ponds_nml rewind ', &
            file=__FILE__, line=__LINE__)
    endif
    nml_error =  1
    if (partit%mype == 0) write(nu_diag,*) subname,' Reading ponds_nml'
    do while (nml_error > 0)
       read(nu_nml, nml=ponds_nml,iostat=nml_error)
    end do
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: ponds_nml reading ', &
            file=__FILE__, line=__LINE__)
    endif

    rewind(unit=nu_nml, iostat=nml_error)
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: snow_nml rewind ', &
            file=__FILE__, line=__LINE__)
    endif
    nml_error =  1
    if (partit%mype == 0) write(nu_diag,*) subname,' Reading snow_nml'
    do while (nml_error > 0)
       read(nu_nml, nml=snow_nml,iostat=nml_error)
    end do
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: snow_nml reading ', &
            file=__FILE__, line=__LINE__)
    endif

    rewind(unit=nu_nml, iostat=nml_error)
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: forcing_nml rewind ', &
            file=__FILE__, line=__LINE__)
    endif
    nml_error =  1
    if (partit%mype == 0) write(nu_diag,*) subname,' Reading forcing_nml'
    do while (nml_error > 0)
       read(nu_nml, nml=forcing_nml,iostat=nml_error)
    end do
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: forcing_nml reading ', &
            file=__FILE__, line=__LINE__)
    endif

    rewind(unit=nu_nml, iostat=nml_error)
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: dynamics_nml rewind ', &
            file=__FILE__, line=__LINE__)
    endif
    nml_error =  1
    if (partit%mype == 0) write(nu_diag,*) subname,' Reading dynamics_nml'
    do while (nml_error > 0)
       read(nu_nml, nml=dynamics_nml,iostat=nml_error)
    end do
    if (nml_error /= 0) then
       call icedrv_system_abort(string=subname//'ERROR: dynamics_nml reading ', &
            file=__FILE__, line=__LINE__)
    endif
    
    close(nu_nml)
    
    if (partit%mype == 0) write(nu_diag,*) subname,' Reading namelist completed'

    !-----------------------------------------------------------------
    ! set up diagnostics output and resolve conflicts
    !-----------------------------------------------------------------
    
    if (partit%mype == 0) write(nu_diag,*) ' '
    if (partit%mype == 0) write(nu_diag,*) '-----------------------------------'
    if (partit%mype == 0) write(nu_diag,*) '  ICEPACK model diagnostic output  '
    if (partit%mype == 0) write(nu_diag,*) '-----------------------------------'
    if (partit%mype == 0) write(nu_diag,*) ' '

    if (ice%whichEVP == 1 .or. ice%whichEVP == 2) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: whichEVP = 1 or 2'
       if (partit%mype == 0) write (nu_diag,*) 'Adaptive or Modified EVP formulations'
       if (partit%mype == 0) write (nu_diag,*) 'are not allowed when using Icepack (yet).'
       if (partit%mype == 0) write (nu_diag,*) 'Standard EVP will be used instead'
       if (partit%mype == 0) write (nu_diag,*) '         whichEVP = 0'
       ice%whichEVP = 0
    endif
    
    if (ncat == 1 .and. kitd == 1) then
       if (partit%mype == 0) write (nu_diag,*) 'Remapping the ITD is not allowed for ncat=1.'
       if (partit%mype == 0) write (nu_diag,*) 'Use kitd = 0 (delta function ITD) with kcatbound = 0'
       if (partit%mype == 0) write (nu_diag,*) 'or for column configurations use kcatbound = -1'
       if (partit%mype == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (ncat /= 1 .and. kcatbound == -1) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: ITD required for ncat > 1'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting kitd and kcatbound to default values'
       kitd = 1
       kcatbound = 0
    endif
    
    rplvl  = c0
    rptopo = c0
    if (tr_pond_lvl ) rplvl  = c1
    if (tr_pond_topo) rptopo = c1
    
    tr_pond = .false. ! explicit melt ponds
    if (rplvl + rptopo > puny) tr_pond = .true.
    if (rplvl + rptopo > c1 + puny) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Must use only one melt pond scheme'
       if (partit%mype == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (tr_pond_lvl .and. .not. tr_lvl) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: tr_pond_lvl=T but tr_lvl=F'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting tr_lvl=T'
       tr_lvl = .true.
    endif
    
    if (tr_pond_lvl .and. abs(hs0) > puny) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: tr_pond_lvl=T and hs0/=0'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting hs0=0'
       hs0 = c0
    endif
    
    if (trim(shortwave(1:4)) /= 'dEdd' .and. tr_pond .and. calc_tsfc) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Must use dEdd shortwave'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: with tr_pond and calc_tsfc=T.'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting shortwave = dEdd'
       shortwave = 'dEdd'
    endif

    if (snwredist(1:3) == 'ITD' .and. .not. tr_snow) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: snwredist on but tr_snow=F'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (snwredist(1:4) == 'bulk' .and. .not. tr_lvl) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: snwredist=bulk but tr_lvl=F'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (snwredist(1:6) == 'ITDrdg' .and. .not. tr_lvl) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: snwredist=ITDrdg but tr_lvl=F'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (use_smliq_pnd .and. .not. snwgrain) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: use_smliq_pnd = T but'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: snow metamorphosis not used'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (use_smliq_pnd .and. .not. tr_snow) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: use_smliq_pnd = T but'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: snow tracers are not active'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (snwgrain .and. .not. tr_snow) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: snwgrain=T but tr_snow=F'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (trim(snw_aging_table) /= 'test') then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: snw_aging_table /= test'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: netcdf not available'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (tr_iso .and. n_iso==0) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: isotopes activated but'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: not allocated in tracer array.'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Activate in compilation script.'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (tr_aero .and. n_aero==0) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: aerosols activated but'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: not allocated in tracer array.'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Activate in compilation script.'
       if (partit%mype == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (tr_aero .and. trim(shortwave(1:4)) /= 'dEdd') then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: aerosols activated but dEdd'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: shortwave is not.'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting shortwave = dEdd'
       shortwave = 'dEdd'
    endif

    if (snwgrain .and. trim(shortwave(1:4)) /= 'dEdd') then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: snow grain radius activated but'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: dEdd shortwave is not.'
    endif
    
    if (snwredist(1:4) /= 'none' .and. trim(shortwave(1:4)) /= 'dEdd') then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: snow redistribution activated but'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: dEdd shortwave is not.'
    endif
    
    rfracmin = min(max(rfracmin,c0),c1)
    rfracmax = min(max(rfracmax,c0),c1)
    
    if (ktherm == 2 .and. .not. calc_Tsfc) then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: ktherm = 2 and calc_Tsfc = F'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting calc_Tsfc = T'
       calc_Tsfc = .true.
    endif
    
    if (ktherm == 1 .and. trim(tfrz_option) /= 'linear_salt') then
       if (partit%mype == 0) write (nu_diag,*) &
            'WARNING: ktherm = 1 and tfrz_option = ',trim(tfrz_option)
       if (partit%mype == 0) write (nu_diag,*) &
            'WARNING: For consistency, set tfrz_option = linear_salt'
    endif
    
    if (ktherm == 2 .and. trim(tfrz_option) /= 'mushy') then
       if (partit%mype == 0) write (nu_diag,*) &
            'WARNING: ktherm = 2 and tfrz_option = ',trim(tfrz_option)
       if (partit%mype == 0) write (nu_diag,*) &
            'WARNING: For consistency, set tfrz_option = mushy'
    endif

    if (ktherm == 1 .and. trim(saltflux_option) /= 'constant') then
       if (mype == 0) write (nu_diag,*) &
            'WARNING: ktherm = 1 and saltflux_option = ',trim(saltflux_option)
       if (mype == 0) write (nu_diag,*) &
            'WARNING: For consistency, set saltflux_option = constant'
    endif
    
    if (ktherm == 0) then
       if (mype == 0) write (nu_diag,*) 'WARNING: ktherm = 0 zero-layer thermodynamics'
       if (mype == 0) write (nu_diag,*) 'WARNING: has been deprecated'
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
    
    if (formdrag) then
       if (trim(atmbndy) == 'constant') then
          if (partit%mype == 0) write (nu_diag,*) 'WARNING: atmbndy = constant not allowed with formdrag'
          if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting atmbndy = similarity'
          atmbndy = 'similarity'
       endif
    
       if (.not. calc_strair) then
          if (partit%mype == 0) write (nu_diag,*) 'WARNING: formdrag=T but calc_strair=F'
          if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting calc_strair=T'
          calc_strair = .true.
       endif
    
       if (.not. tr_lvl) then
          if (partit%mype == 0) write (nu_diag,*) 'WARNING: formdrag=T but tr_lvl=F'
          if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting tr_lvl=T'
          tr_lvl = .true.
          !max_ntrcr = max_ntrcr + 2 ! tr_lvl brings two more tracers
          ! (double counting? - not present in original code frank.kauker@awi.de)
       endif
    endif

    if (trim(fbot_xfer_type) == 'Cdn_ocn' .and. .not. formdrag)  then
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: formdrag=F but fbot_xfer_type=Cdn_ocn'
       if (partit%mype == 0) write (nu_diag,*) 'WARNING: Setting fbot_xfer_type = constant'
       fbot_xfer_type = 'constant'
    endif
    
    wave_spec = .false.
    if (tr_fsd .and. (trim(wave_spec_type) /= 'none')) wave_spec = .true.
    if (tr_fsd .and. (trim(wave_spec_type) == 'none')) then
       if (mype == 0) write (nu_diag,*) 'WARNING: tr_fsd=T but wave_spec=F - not recommended'
    endif

    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call fesom_to_icepack_para'//achar(27)//'[0m'
    call fesom_to_icepack_para() ! no atmoshperic data yet frank.kauker@awi.de                                         
    !-----------------------------------------------------------------
    ! spew
    !-----------------------------------------------------------------
           
    if (partit%mype == 0) then
       write(nu_diag,*) ' Document ice_in namelist parameters:'
       write(nu_diag,*) ' ==================================== '
       write(nu_diag,*) ' '
       write(nu_diag,1010) ' conserv_check             = ', conserv_check
       write(nu_diag,1020) ' kitd                      = ', kitd
       write(nu_diag,1020) ' kcatbound                 = ', kcatbound
       write(nu_diag,1020) ' ndtd                      = ', ndtd
       write(nu_diag,1020) ' kstrength                 = ', kstrength
       write(nu_diag,1020) ' krdg_partic               = ', krdg_partic
       write(nu_diag,1020) ' krdg_redist               = ', krdg_redist
       if (krdg_redist == 1) &
            write(nu_diag,1000) ' mu_rdg                    = ', mu_rdg
       if (kstrength == 1) &
            write(nu_diag,1000) ' Cf                        = ', Cf
       write(nu_diag,1000) ' ksno                      = ', ksno
       write(nu_diag,1030) ' shortwave                 = ', trim(shortwave)
       write(nu_diag,1000) ' -------------------------------'
       write(nu_diag,1000) ' BGC coupling is switched OFF '
       write(nu_diag,1000) ' not implemented with FESOM2  '
       write(nu_diag,1000) ' -------------------------------'
       if (trim(shortwave(1:4)) == 'dEdd') then
          write(nu_diag,1000) ' R_ice                     = ', R_ice
          write(nu_diag,1000) ' R_pnd                     = ', R_pnd
          write(nu_diag,1000) ' R_snw                     = ', R_snw
          write(nu_diag,1000) ' dT_mlt                    = ', dT_mlt
          write(nu_diag,1000) ' rsnw_mlt                  = ', rsnw_mlt
          write(nu_diag,1000) ' kalg                      = ', kalg
          write(nu_diag,1000) ' hp1                       = ', hp1
          write(nu_diag,1000) ' hs0                       = ', hs0
       else
          write(nu_diag,1030) ' albedo_type               = ', trim(albedo_type)
          write(nu_diag,1000) ' albicev                   = ', albicev
          write(nu_diag,1000) ' albicei                   = ', albicei
          write(nu_diag,1000) ' albsnowv                  = ', albsnowv
          write(nu_diag,1000) ' albsnowi                  = ', albsnowi
          write(nu_diag,1000) ' albocn                    = ', albocn
          write(nu_diag,1000) ' ahmax                     = ', ahmax
       endif
       if (trim(shortwave) == 'dEdd_snicar_ad') then
          write(nu_diag,1030) ' snw_ssp_table             = ', trim(snw_ssp_table)
       endif
       write(nu_diag,1010) ' sw_redist                 = ', sw_redist
       write(nu_diag,1005) ' sw_frac                   = ', sw_frac
       write(nu_diag,1005) ' sw_dtemp                  = ', sw_dtemp
       write(nu_diag,1000) ' rhos                      = ', rhosno
       write(nu_diag,1000) ' rhoi                      = ', rhoice
       write(nu_diag,1000) ' rhow                      = ', rhowat
       write(nu_diag,1000) ' rhofwt                    = ', rhofwt
       write(nu_diag,1000) ' rfracmin                  = ', rfracmin
       write(nu_diag,1000) ' rfracmax                  = ', rfracmax
       if (tr_pond_lvl) then
          write(nu_diag,1000) ' hs1                       = ', hs1
          write(nu_diag,1000) ' dpscale                   = ', dpscale
          write(nu_diag,1030) ' frzpnd                    = ', trim(frzpnd)
       endif
       if (tr_pond .and. .not. tr_pond_lvl) &
          write(nu_diag,1000) ' pndaspect                 = ', pndaspect
       if (tr_snow) then
          write(nu_diag,1030) ' snwredist                 = ', trim(snwredist)
          write(nu_diag,1010) ' snwgrain                  = ', snwgrain
          write(nu_diag,1010) ' use_smliq_pnd             = ', use_smliq_pnd
          write(nu_diag,1030) ' snw_aging_table           = ', trim(snw_aging_table)
          write(nu_diag,1000) ' rsnw_fall                 = ', rsnw_fall
          write(nu_diag,1000) ' rsnw_tmax                 = ', rsnw_tmax
          write(nu_diag,1000) ' rhosnew                   = ', rhosnew
          write(nu_diag,1000) ' rhosmin                   = ', rhosmin
          write(nu_diag,1000) ' rhosmax                   = ', rhosmax
          write(nu_diag,1000) ' windmin                   = ', windmin
          write(nu_diag,1000) ' drhosdwind                = ', drhosdwind
          write(nu_diag,1000) ' snwlvlfac                 = ', snwlvlfac
       endif
       write(nu_diag,1020)    ' ktherm                    = ', ktherm
       if (ktherm == 1) &
          write(nu_diag,1030)     ' conduct                   = ', trim(conduct)
       write(nu_diag,1005)    ' emissivity                = ', emissivity
       write(nu_diag,1005)    ' zlvl_v                    = ', zlvl_v
       write(nu_diag,1005)    ' zlvl_s                    = ', zlvl_s
       if (ktherm == 2) then
          write(nu_diag,1005) ' a_rapid_mode              = ', a_rapid_mode
          write(nu_diag,1005) ' Rac_rapid_mode            = ', Rac_rapid_mode
          write(nu_diag,1005) ' aspect_rapid_mode         = ', aspect_rapid_mode
          write(nu_diag,1005) ' dSdt_slow_mode            = ', dSdt_slow_mode
          write(nu_diag,1005) ' phi_c_slow_mode           = ', phi_c_slow_mode
          write(nu_diag,1005) ' phi_i_mushy               = ', phi_i_mushy
          write(nu_diag,1005) ' Tliquidus_max             = ', Tliquidus_max
       endif
       write(nu_diag,1030) ' atmbndy                   = ', trim(atmbndy)
       write(nu_diag,1010) ' formdrag                  = ', formdrag
       write(nu_diag,1010) ' highfreq                  = ', highfreq
       write(nu_diag,1020) ' natmiter                  = ', natmiter
       write(nu_diag,1005) ' atmiter_conv              = ', atmiter_conv
       write(nu_diag,1010) ' calc_strair               = ', calc_strair
       write(nu_diag,1010) ' calc_Tsfc                 = ', calc_Tsfc    
       write(nu_diag,1010) ' calc_dragio               = ', calc_dragio
       write(nu_diag,1005) ' floediam                  = ', floediam
       write(nu_diag,1005) ' hfrazilmin                = ', hfrazilmin
       write(nu_diag,1030) ' cpl_frazil                = ', trim(cpl_frazil)
       write(nu_diag,1010) ' update_ocn_f              = ', update_ocn_f
       write(nu_diag,1010) ' wave_spec                 = ', wave_spec
       if (wave_spec) &
          write(nu_diag,1030) ' wave_spec_type            = ', trim(wave_spec_type)
       write(nu_diag,1005) ' dragio                    = ', dragio
       write(nu_diag,1010) ' l_mpond_fresh             = ', l_mpond_fresh
       write(nu_diag,1005) ' ustar_min                 = ', ustar_min
       write(nu_diag,1005) ' hi_min                    = ', hi_min
       write(nu_diag,1030) ' fbot_xfer_type            = ', trim(fbot_xfer_type)
       write(nu_diag,1010) ' oceanmixed_ice            = ', oceanmixed_ice
       write(nu_diag,1030) ' tfrz_option               = ', trim(tfrz_option)
       write(nu_diag,*)    ' saltflux_option           = ', trim(saltflux_option)
       if (trim(saltflux_option)=='constant') &
          write(nu_diag,1005) ' ice_ref_salinity          = ', ice_ref_salinity
       write(nu_diag,1030) ' congel_freeze             = ', trim(congel_freeze)
       ! tracers
       write(nu_diag,1010) ' tr_iage                   = ', tr_iage
       write(nu_diag,1010) ' tr_FY                     = ', tr_FY
       write(nu_diag,1010) ' tr_lvl                    = ', tr_lvl
       write(nu_diag,1010) ' tr_pond_lvl               = ', tr_pond_lvl
       write(nu_diag,1010) ' tr_pond_topo              = ', tr_pond_topo
       write(nu_diag,1010) ' tr_snow                   = ', tr_snow
       write(nu_diag,1010) ' tr_aero                   = ', tr_aero
       write(nu_diag,1010) ' tr_fsd                    = ', tr_fsd
       write(nu_diag,1010) ' tr_iso                    = ', tr_iso
    endif

    !-----------------------------------------------------------------
    ! Compute number of tracers
    !-----------------------------------------------------------------

    nt_Tsfc = 1           ! index tracers, starting with Tsfc = 1
    ntrcr = 1             ! count tracers, starting with Tsfc = 1
    
    nt_qice = ntrcr + 1
    ntrcr = ntrcr + nilyr ! qice in nilyr layers
    
    nt_qsno = ntrcr + 1
    ntrcr = ntrcr + nslyr ! qsno in nslyr layers
    
    nt_sice = ntrcr + 1
    ntrcr = ntrcr + nilyr ! sice in nilyr layers
    
    nt_iage = max_ntrcr
    if (tr_iage) then
       ntrcr = ntrcr + 1
       nt_iage = ntrcr   ! chronological ice age
    endif
 
    nt_FY = max_ntrcr
    if (tr_FY) then
       ntrcr = ntrcr + 1
       nt_FY = ntrcr     ! area of first year ice
    endif
 
    nt_alvl = max_ntrcr
    nt_vlvl = max_ntrcr
    if (tr_lvl) then
       ntrcr = ntrcr + 1
       nt_alvl = ntrcr   ! area of level ice
       ntrcr = ntrcr + 1
       nt_vlvl = ntrcr   ! volume of level ice
    endif
 
    nt_apnd = max_ntrcr
    nt_hpnd = max_ntrcr
    nt_ipnd = max_ntrcr
    if (tr_pond) then            ! all explicit melt pond schemes
       ntrcr = ntrcr + 1
       nt_apnd = ntrcr
       ntrcr = ntrcr + 1
       nt_hpnd = ntrcr
       if (tr_pond_lvl) then
          ntrcr = ntrcr + 1    ! refrozen pond ice lid thickness
          nt_ipnd = ntrcr      ! on level-ice ponds (if frzpnd='hlid')
       endif
       if (tr_pond_topo) then
          ntrcr = ntrcr + 1    !
          nt_ipnd = ntrcr      ! refrozen pond ice lid thickness
       endif
    endif
 
    nt_smice = max_ntrcr
    nt_smliq = max_ntrcr
    nt_rhos  = max_ntrcr
    nt_rsnw  = max_ntrcr
    if (tr_snow) then
       nt_smice = ntrcr + 1
       ntrcr = ntrcr + nslyr     ! mass of ice in nslyr snow layers                     
       nt_smliq = ntrcr + 1
       ntrcr = ntrcr + nslyr     ! mass of liquid in nslyr snow layers
       nt_rhos = ntrcr + 1
       ntrcr = ntrcr + nslyr     ! snow density in nslyr layers
       nt_rsnw = ntrcr + 1
       ntrcr = ntrcr + nslyr     ! snow grain radius in nslyr layers
    endif

    nt_fsd = max_ntrcr
    if (tr_fsd) then
       nt_fsd = ntrcr + 1       ! floe size distribution
       ntrcr = ntrcr + nfsd
    end if
    
    nt_isosno = max_ntrcr
    nt_isoice = max_ntrcr
    if (tr_iso) then            ! isotopes
       nt_isosno = ntrcr + 1
       ntrcr = ntrcr + n_iso    ! n_iso species in snow
       nt_isoice = ntrcr + 1
       ntrcr = ntrcr + n_iso    ! n_iso species in ice
    end if

    nt_aero = max_ntrcr - 4*n_aero
    if (tr_aero) then
       nt_aero = ntrcr + 1
       ntrcr = ntrcr + 4*n_aero ! 4 dEdd layers, n_aero species
    endif
 
    if (ntrcr > max_ntrcr-1) then
       if (partit%mype == 0) write(ice_stderr,*) 'max_ntrcr-1 < number of namelist tracers'
       if (partit%mype == 0) write(ice_stderr,*) 'max_ntrcr-1 = ',max_ntrcr-1,' ntrcr = ',ntrcr
       call icepack_warnings_flush(ice_stderr)
       call icedrv_system_abort(file=__FILE__,line=__LINE__)
    endif
 
    if (partit%mype == 0) then 
       write(nu_diag,1010) ' tr_iage                   = ', tr_iage

       write(nu_diag,*)    ' '
       write(nu_diag,1020) ' max_ntrcr                 = ', max_ntrcr
       write(nu_diag,1020) ' ntrcr                     = ', ntrcr
       write(nu_diag,*)    ' '
       write(nu_diag,1020) ' nt_qice                   = ', nt_qice
       write(nu_diag,1020) ' nt_qsno                   = ', nt_qsno
       write(nu_diag,1020) ' nt_sice                   = ', nt_sice
       write(nu_diag,*)' '
       write(nu_diag,1020) ' ncat                      = ', ncat
       write(nu_diag,1020) ' nilyr                     = ', nilyr
       write(nu_diag,1020) ' nslyr                     = ', nslyr
       write(nu_diag,1020) ' nblyr                     = ', nblyr
       write(nu_diag,1020) ' nfsd                      = ', nfsd
       write(nu_diag,1020) ' n_iso                     = ', n_iso
       write(nu_diag,1020) ' n_aero                    = ', n_aero
       write(nu_diag,*)' '
       if (formdrag) then
          if (nt_apnd==0) then
             write(ice_stderr,*)'ERROR: nt_apnd:',nt_apnd
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          elseif (nt_hpnd==0) then
             write(ice_stderr,*)'ERROR: nt_hpnd:',nt_hpnd
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          elseif (nt_ipnd==0) then
             write(ice_stderr,*)'ERROR: nt_ipnd:',nt_ipnd
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          elseif (nt_alvl==0) then
             write(ice_stderr,*)'ERROR: nt_alvl:',nt_alvl
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          elseif (nt_vlvl==0) then
             write(ice_stderr,*)'ERROR: nt_vlvl:',nt_vlvl
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
       endif
    endif ! mype == 0

1000 format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
1005 format (a30,2x,f10.6) ! float
1010 format (a30,2x,l6)    ! logical
1020 format (a30,2x,i6)    ! integer
1030 format (a30,   a8)    ! character
1040 format (a30,2x,6i6)   ! integer
1050 format (a30,2x,6a6)   ! character

    !-----------------------------------------------------------------
    ! set Icepack values
    !-----------------------------------------------------------------  

    call icepack_init_parameters(ustar_min_in=ustar_min, Cf_in=Cf, &
         albicev_in=albicev, albicei_in=albicei, ksno_in=ksno, &
         albsnowv_in=albsnowv, albsnowi_in=albsnowi, hi_min_in=hi_min, &
         natmiter_in=natmiter, ahmax_in=ahmax, shortwave_in=shortwave, &
         atmiter_conv_in = atmiter_conv, calc_dragio_in=calc_dragio, &
         albedo_type_in=albedo_type, R_ice_in=R_ice, R_pnd_in=R_pnd, &
         R_snw_in=R_snw, dT_mlt_in=dT_mlt, rsnw_mlt_in=rsnw_mlt, &
         kstrength_in=kstrength, krdg_partic_in=krdg_partic, &
         krdg_redist_in=krdg_redist, mu_rdg_in=mu_rdg, &
         atmbndy_in=atmbndy, calc_strair_in=calc_strair, &
         formdrag_in=formdrag, highfreq_in=highfreq, &
         emissivity_in=emissivity, &
         kitd_in=kitd, kcatbound_in=kcatbound, hs0_in=hs0, &
         dpscale_in=dpscale, frzpnd_in=frzpnd, &
         rfracmin_in=rfracmin, rfracmax_in=rfracmax, &
         pndaspect_in=pndaspect, hs1_in=hs1, hp1_in=hp1, &
         floediam_in=floediam, hfrazilmin_in=hfrazilmin, &
         ktherm_in=ktherm, calc_Tsfc_in=calc_Tsfc, &
         conduct_in=conduct, a_rapid_mode_in=a_rapid_mode, &
         update_ocn_f_in=update_ocn_f, cpl_frazil_in=cpl_frazil, &
         Rac_rapid_mode_in=Rac_rapid_mode, &
         aspect_rapid_mode_in=aspect_rapid_mode, &
         dSdt_slow_mode_in=dSdt_slow_mode, &
         phi_c_slow_mode_in=phi_c_slow_mode, Tliquidus_max_in=Tliquidus_max, &
         phi_i_mushy_in=phi_i_mushy, conserv_check_in=conserv_check, &
         tfrz_option_in=tfrz_option, saltflux_option_in=saltflux_option, &
         ice_ref_salinity_in=ice_ref_salinity, kalg_in=kalg, &
         fbot_xfer_type_in=fbot_xfer_type, &
         congel_freeze_in=congel_freeze, &
         wave_spec_type_in=wave_spec_type, wave_spec_in=wave_spec, &
         sw_redist_in=sw_redist, sw_frac_in=sw_frac, sw_dtemp_in=sw_dtemp, &
         snwredist_in=snwredist, use_smliq_pnd_in=use_smliq_pnd, &
         snw_aging_table_in=snw_aging_table, &
         snwgrain_in=snwgrain, rsnw_fall_in=rsnw_fall, rsnw_tmax_in=rsnw_tmax, &
         rhosnew_in=rhosnew, rhosmin_in=rhosmin, rhosmax_in=rhosmax, &
         windmin_in=windmin, drhosdwind_in=drhosdwind, snwlvlfac_in=snwlvlfac)
    call icepack_init_tracer_sizes(ntrcr_in=ntrcr, &
         ncat_in=ncat, nilyr_in=nilyr, nslyr_in=nslyr, nblyr_in=nblyr, &
         nfsd_in=nfsd, n_iso_in=n_iso, n_aero_in=n_aero)
    call icepack_init_tracer_flags(tr_iage_in=tr_iage, &
         tr_FY_in=tr_FY, tr_lvl_in=tr_lvl, tr_aero_in=tr_aero, &
         tr_iso_in=tr_iso, tr_snow_in=tr_snow, &
         tr_pond_in=tr_pond, &
         tr_pond_lvl_in=tr_pond_lvl, &
         tr_pond_topo_in=tr_pond_topo, tr_fsd_in=tr_fsd)
    call icepack_init_tracer_indices(nt_Tsfc_in=nt_Tsfc, &
         nt_sice_in=nt_sice, nt_qice_in=nt_qice, &
         nt_qsno_in=nt_qsno, nt_iage_in=nt_iage, &
         nt_fy_in=nt_fy, nt_alvl_in=nt_alvl, nt_vlvl_in=nt_vlvl, &
         nt_apnd_in=nt_apnd, nt_hpnd_in=nt_hpnd, nt_ipnd_in=nt_ipnd, &
         nt_smice_in=nt_smice, nt_smliq_in=nt_smliq, &
         nt_rhos_in=nt_rhos, nt_rsnw_in=nt_rsnw, &
         nt_aero_in=nt_aero, nt_fsd_in=nt_fsd, &
         nt_isosno_in=nt_isosno, nt_isoice_in=nt_isoice)
    
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    call mpi_barrier(p_partit%mpi_comm_fesom, mpi_error)
    
  end subroutine set_icepack

  !=======================================================================

  module subroutine set_grid_icepack(mesh)

    use mod_mesh

    implicit none
         
    integer (kind=int_kind)                      :: i
    real(kind=dbl_kind)                          :: puny
    real(kind=dbl_kind), dimension(:,:), pointer :: coord_nod2D  
    character(len=*), parameter                  :: subname = '(init_grid_icepack)'
    type(t_mesh), intent(in), target             :: mesh
                    
                     
    !-----------------------------------------------------------------
    ! query Icepack values
    !-----------------------------------------------------------------
    
    call icepack_query_parameters(puny_out=puny)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)
    
    coord_nod2D(1:2,1:nx) => mesh%coord_nod2D    
    
    !-----------------------------------------------------------------
    ! create hemisphereic masks
    !-----------------------------------------------------------------
    
    lmask_n(:) = .false.
    lmask_s(:) = .false.
    
    do i = 1, nx
       if (coord_nod2D(2,i) >= -puny) lmask_n(i) = .true. ! N. Hem.
       if (coord_nod2D(2,i) <  -puny) lmask_s(i) = .true. ! S. Hem.
    enddo
    
    !-----------------------------------------------------------------
    ! longitudes and latitudes
    !-----------------------------------------------------------------
    
    lon_val(:) = coord_nod2D(1,:)
    lat_val(:) = coord_nod2D(2,:)
    
  end subroutine set_grid_icepack
  
  !=======================================================================
  
end submodule icedrv_set







