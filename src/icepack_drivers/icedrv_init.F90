!=======================================================================
!
! This module defines and and initializes the namelists
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!
!=======================================================================

      module icedrv_init

      use icedrv_kinds
      use icedrv_constants, only: nu_diag, ice_stdout, nu_diag_out, nu_nml
      use icedrv_constants, only: c0, c1, c2, c3, p2, p5, puny
      use icepack_intfc,    only: icepack_init_parameters
      use icepack_intfc,    only: icepack_init_fsd
      use icepack_intfc,    only: icepack_init_tracer_flags
      use icepack_intfc,    only: icepack_init_tracer_sizes
      use icepack_intfc,    only: icepack_init_tracer_indices
      use icepack_intfc,    only: icepack_init_trcr
      use icepack_intfc,    only: icepack_query_parameters
      use icepack_intfc,    only: icepack_query_tracer_flags
      use icepack_intfc,    only: icepack_query_tracer_sizes
      use icepack_intfc,    only: icepack_query_tracer_indices
      use icepack_intfc,    only: icepack_warnings_flush
      use icepack_intfc,    only: icepack_warnings_aborted
      use icedrv_system,    only: icedrv_system_abort 

      implicit none
      private
      public :: read_namelist_icepack
      
      contains 

      subroutine read_namelist_icepack(icepack_settings)

          use icedrv_settings

          implicit none

          character(len=char_len)    :: nml_filename
          integer  (kind=int_kind)   :: nml_error,    & ! namelist i/o error flag
                                        n               ! loop index
          type(t_icepack_settings), intent(inout), target :: icepack_settings

#include "associate_icepack_settings.h"

          !-----------------------------------------------------------------
          ! Namelist definition
          !-----------------------------------------------------------------
 
          nml_filename = 'namelist.icepack'        ! name of icepack namelist file
   
          namelist / env_nml /                                                &
             nicecat, nfsdcat, nicelyr, nsnwlyr, ntraero, trzaero, tralg,     &
             trdoc,   trdic,   trdon,   trfed,   trfep,   nbgclyr, trbgcz,    &
             trzs,    trbri,   trage,   trfy,    trlvl,   trpnd,   trbgcs

          namelist /thermo_nml/                                               &
             kitd,           ktherm,          conduct,                        &
             a_rapid_mode,   Rac_rapid_mode,  aspect_rapid_mode,              &
             dSdt_slow_mode, phi_c_slow_mode, phi_i_mushy

          namelist /dynamics_nml/                                              &
             kstrength,      krdg_partic,    krdg_redist,    mu_rdg,           &
             Cf

          namelist /shortwave_nml/                                             &
             shortwave,      albedo_type,                                      &
             albicev,        albicei,         albsnowv,      albsnowi,         &
             ahmax,          R_ice,           R_pnd,         R_snw,&
             dT_mlt,         rsnw_mlt,        kalg

          namelist /ponds_nml/                                                 &
             hs0,            dpscale,         frzpnd,                          &
             rfracmin,       rfracmax,        pndaspect,     hs1,              &
             hp1

          namelist /tracer_nml/                                                &
             tr_iage,      tr_FY,        tr_lvl,       tr_pond_cesm,           &
             tr_pond_lvl,  tr_pond_topo, tr_aero,      tr_fsd

          namelist /forcing_nml/                                               &
             atmbndy,         calc_strair,     calc_Tsfc,                      &
             update_ocn_f,    l_mpond_fresh,   ustar_min,                      &
             fbot_xfer_type,  oceanmixed_ice,  emissivity,                     &
             formdrag,        highfreq,        natmiter,                       &
             tfrz_option,     wave_spec_type

          !-----------------------------------------------------------------
          ! env namelist - STANDARD VALUES
          !-----------------------------------------------------------------

          nicecat   = 5           ! number of ice thickness categories
          nicelyr   = 4           ! number of vertical layers in the ice
          nsnwlyr   = 4           ! number of vertical layers in the snow
          nfsdcat   = 1           ! number of floe size categories
          ntraero   = 0           ! number of aerosol tracers (up to max_aero in ice_domain_size.f90)
          trzaero   = 0           ! number of z aerosol tracers (up to max_aero = 6)
          tralg     = 0           ! number of algal tracers (up to max_algae = 3)
          trdoc     = 0           ! number of dissolve organic carbon (up to max_doc = 3)
          trdic     = 0           ! number of dissolve inorganic carbon (up to max_dic = 1)
          trdon     = 0           ! number of dissolve organic nitrogen (up to max_don = 1)
          trfed     = 0           ! number of dissolved iron tracers (up to max_fe  = 2)
          trfep     = 0           ! number of particulate iron tracers (up to max_fe  = 2)
          nbgclyr   = 4           ! number of zbgc layers 
          trbgcz    = 0           ! set to 1 for zbgc tracers (needs trbgcs = 0 and trbri = 1)
          trzs      = 0           ! set to 1 for zsalinity tracer (needs trbri = 1)
          trbri     = 0           ! set to 1 for brine height tracer 
          trage     = 0           ! set to 1 for ice age tracer
          trfy      = 0           ! set to 1 for first-year ice area tracer
          trlvl     = 0           ! set to 1 for level and deformed ice tracers
          trpnd     = 0           ! set to 1 for melt pond tracers
          trbgcs    = 0           ! set to 1 for skeletal layer tracers (needs

          !-----------------------------------------------------------------
          ! Read namelist env_nml
          !-----------------------------------------------------------------

          open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
          if (nml_error /= 0) then
             nml_error = -1
          else
             nml_error =  1
          end if

          if (nml_error > 0)
              print*,'Reading namelist file   ',nml_filename

              print*,'Reading setup_nml'
              read(nu_nml, nml=env_nml,iostat=nml_error)
          end if

          if (nml_error == 0) close(nu_nml)
          if (nml_error /= 0) then
             write(ice_stdout,*) 'error reading namelist'
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
             close(nu_nml)
          end if

          !-----------------------------------------------------------------
          ! Derived quantities used by the icepack model
          !-----------------------------------------------------------------

          ncat      = nicecat    ! number of categories
          nfsd      = nfsdcat    ! number of floe size categories
          nilyr     = nicelyr    ! number of ice layers per category
          nslyr     = nsnwlyr    ! number of snow layers per category
          n_aero    = ntraero    ! number of aerosols in use
          n_zaero   = trzaero    ! number of z aerosols in use
          n_algae   = tralg      ! number of algae in use
          n_doc     = trdoc      ! number of DOC pools in use
          n_dic     = trdic      ! number of DIC pools in use
          n_don     = trdon      ! number of DON pools in use
          n_fed     = trfed      ! number of Fe  pools in use dissolved Fe
          n_fep     = trfep      ! number of Fe  pools in use particulate Fe
          nfreq     = 25         ! number of wave frequencies ! HARDWIRED FOR NOW
          nblyr     = nbgclyr    ! number of bio/brine layers per category
                                 ! maximum number of biology tracers +
                                 ! aerosols
                                 ! *** add to kscavz in
                                 ! icepack_zbgc_shared.F90
          n_bgc     = (n_algae*2 + n_doc + n_dic + n_don + n_fed + n_fep +n_zaero &
                    + 8)         ! nit, am, sil, dmspp, dmspd, dms, pon, humic
          nltrcr    = (n_bgc*trbgcz+trzs)*trbri ! number of zbgc (includes zaero)
                                                ! and zsalinity tracers
          max_nsw   = (nilyr+nslyr+2) & ! total chlorophyll plus aerosols
                    * (1+trzaero)       ! number of tracers active in shortwave calculation
          max_ntrcr =   1         & ! 1 = surface temperature
                    + nilyr       & ! ice salinity
                    + nilyr       & ! ice enthalpy
                    + nslyr       & ! snow enthalpy
                                    !!!!! optional tracers:
                    + nfsd        & ! number of floe size categories
                    + trage       & ! age
                    + trfy        & ! first-year area
                    + trlvl*2     & ! level/deformed ice
                    + trpnd*3     & ! ponds
                    + n_aero*4    & ! number of aerosols * 4 aero layers
                    + trbri       & ! brine height
                    + trbgcs*n_bgc                 & ! skeletal layer BGC
                    + trzs  *trbri* nblyr          & ! zsalinity  (off if TRBRI=0)
                    + n_bgc*trbgcz*trbri*(nblyr+3) & ! zbgc (off if TRBRI=0)
                    + n_bgc*trbgcz                 & ! mobile/stationary phase tracer
                    + 1             ! for unused tracer flags

          !-----------------------------------------------------------------
          ! query Icepack default values
          !-----------------------------------------------------------------

           call icepack_query_parameters(ustar_min_out=ustar_min, Cf_out=Cf,     &
                albicev_out=albicev, albicei_out=albicei,                        &
                albsnowv_out=albsnowv, albsnowi_out=albsnowi,                    &
                natmiter_out=natmiter, ahmax_out=ahmax, shortwave_out=shortwave, &
                albedo_type_out=albedo_type, R_ice_out=R_ice, R_pnd_out=R_pnd,   &
                R_snw_out=R_snw, dT_mlt_out=dT_mlt, rsnw_mlt_out=rsnw_mlt,       &
                kstrength_out=kstrength, krdg_partic_out=krdg_partic,            &
                krdg_redist_out=krdg_redist, mu_rdg_out=mu_rdg,                  &
                atmbndy_out=atmbndy, calc_strair_out=calc_strair,                &
                formdrag_out=formdrag, highfreq_out=highfreq,                    &
                emissivity_out=emissivity,                                       &
                kitd_out=kitd, kcatbound_out=kcatbound, hs0_out=hs0,             &
                dpscale_out=dpscale, frzpnd_out=frzpnd,                          &
                rfracmin_out=rfracmin, rfracmax_out=rfracmax,                    &
                pndaspect_out=pndaspect, hs1_out=hs1, hp1_out=hp1,               &
                ktherm_out=ktherm, calc_Tsfc_out=calc_Tsfc,                      &
                update_ocn_f_out = update_ocn_f,                                 &
                conduct_out=conduct, a_rapid_mode_out=a_rapid_mode,              &
                Rac_rapid_mode_out=Rac_rapid_mode,                               &
                aspect_rapid_mode_out=aspect_rapid_mode,                         &
                dSdt_slow_mode_out=dSdt_slow_mode,                               &
                phi_c_slow_mode_out=phi_c_slow_mode,                             &
                phi_i_mushy_out=phi_i_mushy,                                     &
                tfrz_option_out=tfrz_option, kalg_out=kalg,                      &
                fbot_xfer_type_out=fbot_xfer_type, puny_out=puny,                &
                wave_spec_type_out=wave_spec_type)
          call icepack_warnings_flush(nu_diag)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__, line=__LINE__)

          !-----------------------------------------------------------------
          ! other default values 
          !-----------------------------------------------------------------

          ndtd = 1                    ! dynamic time steps per thermodynamic time step
          l_mpond_fresh = .false.     ! logical switch for including meltpond freshwater
                                      ! flux feedback to ocean model
          oceanmixed_ice  = .false.   ! if true, use internal ocean mixed layer
          wave_spec_type  = 'none'    ! type of wave spectrum forcing

          tr_iage      = .false.      ! ice age
          tr_FY        = .false.      ! ice age
          tr_lvl       = .false.      ! level ice 
          tr_pond_cesm = .false.      ! CESM melt ponds
          tr_pond_lvl  = .false.      ! level-ice melt ponds
          tr_pond_topo = .false.      ! explicit melt ponds (topographic)
          tr_aero      = .false.      ! aerosols
          tr_fsd       = .false.      ! floe size distribution


          !-----------------------------------------------------------------
          ! read from input file
          !-----------------------------------------------------------------

          open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
          if (nml_error /= 0) then
             nml_error = -1
          else
             nml_error =  1
          endif
      
          do while (nml_error > 0)
             print*,'Reading namelist file   ',nml_filename

             print*,'Reading setup_nml'
             read(nu_nml, nml=setup_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             print*,'Reading grid_nml'
             read(nu_nml, nml=grid_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             print*,'Reading tracer_nml'
             read(nu_nml, nml=tracer_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             print*,'Reading thermo_nml'
             read(nu_nml, nml=thermo_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             print*,'Reading shortwave_nml'
             read(nu_nml, nml=shortwave_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             print*,'Reading ponds_nml'
             read(nu_nml, nml=ponds_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             print*,'Reading forcing_nml'
             read(nu_nml, nml=forcing_nml,iostat=nml_error)
             if (nml_error /= 0) exit
          end do
          if (nml_error == 0) close(nu_nml)
          if (nml_error /= 0) then
             write(ice_stdout,*) 'error reading namelist'
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
          close(nu_nml)

      end subroutine read_namelist_icepack
!=======================================================================

      end module icedrv_init







