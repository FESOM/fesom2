!=======================================================================
!
! Defines the domain size, number of categories and layers.
!
! author L. Zampieri
!
!=======================================================================

      module icedrv_domain_size

      use icedrv_kinds

!=======================================================================

          implicit none
      
          public                                                             &
                  nx, ncat, nfsd, nilyr, nslyr, n_aero, n_zaero, n_algae,    &
                  n_doc, n_dic, n_don, n_fed, n_fep, nblyr, n_bgc, nltrcr,   &
                  max_nsw, max_ntrcr, nfreq, ndtd

          private

          ! setting variables used by the model

          integer (kind=int_kind), save  :: nx                   ! number of grid cells and gost cells for each mesh partition
          integer (kind=int_kind), save  :: ncat                 ! number of categories in use
          integer (kind=int_kind), save  :: nfsd                 ! number of floe size categories in use
          integer (kind=int_kind), save  :: nilyr                ! number of ice layers per category in use
          integer (kind=int_kind), save  :: nslyr                ! number of snow layers per category in use
          integer (kind=int_kind), save  :: n_aero               ! number of aerosols in use
          integer (kind=int_kind), save  :: n_zaero              ! number of z aerosols in use
          integer (kind=int_kind), save  :: n_algae              ! number of algae in use
          integer (kind=int_kind), save  :: n_doc                ! number of DOC pools in use
          integer (kind=int_kind), save  :: n_dic                ! number of DIC pools in use
          integer (kind=int_kind), save  :: n_don                ! number of DON pools in use
          integer (kind=int_kind), save  :: n_fed                ! number of Fe pools in use dissolved Fe
          integer (kind=int_kind), save  :: n_fep                ! number of Fe pools in use particulate Fe
          integer (kind=int_kind), save  :: nblyr                ! number of bio/brine layers per category
          integer (kind=int_kind), save  :: n_bgc                ! nit, am, sil, dmspp, dmspd, dms, pon, humic
          integer (kind=int_kind), save  :: nltrcr               ! number of zbgc (includes zaero) and zsalinity tracers
          integer (kind=int_kind), save  :: max_nsw              ! number of tracers active in shortwave calculation
          integer (kind=int_kind), save  :: max_ntrcr            ! number of tracers in total
          integer (kind=int_kind), save  :: nfreq                ! number of wave frequencies ! HARDWIRED FOR NOW
          integer (kind=int_kind), save  :: ndtd                 ! dynamic time steps per thermodynamic time step

!=======================================================================    

      end module icedrv_domain_size

!=======================================================================

