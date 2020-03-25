!=======================================================================
!
! Defines and initializes namelists
!
! author L. Zampieri
!
!=======================================================================

      module icedrv_namelist

      use icedrv_kinds

      implicit none
      private
      public :: namelist_settings
      
      contains 

      subroutine namelist_settings(icepack_settings)

          use icedrv_settings

          implicit none
          character(len=100)                              :: nmlfile
          type(t_icepack_settings), intent(inout), target :: icepack_settings

#include "associate_icepack_settings.h"

          ! Standard values
          
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

          ! Read namelist

          namelist / env_nml /  nicecat, nfsdcat, nicelyr, nsnwlyr, ntraero, trzaero, tralg,     &
                                trdoc,   trdic,   trdon,   trfed,   trfep,   nbgclyr, trbgcz,    &
                                trzs,    trbri,   trage,   trfy,    trlvl,   trpnd,   trbgcs

          nmlfile ='namelist.icepack'        ! name of icepack namelist file
          open (10,file=nmlfile)
          read (10,NML=env_nml)
          close (10)

          ! Derived quantities used by icepack

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

      end subroutine namelist_settings
!=======================================================================

      end module icedrv_namelist







