!=======================================================================
!
! Defines the domain size, number of categories and layers.
!
! author L. Zampieri
!a
!=======================================================================

      module icedrv_settings

      use icedrv_kinds

!=======================================================================

      implicit none

      type t_icepack_settings

          ! env namelist

          integer (kind=int_kind)  :: nicecat              ! number of ice thickness categories
          integer (kind=int_kind)  :: nfsdcat              ! number of floe size categories
          integer (kind=int_kind)  :: nicelyr              ! number of vertical layers in the ice
          integer (kind=int_kind)  :: nsnwlyr              ! number of vertical layers in the snow
          integer (kind=int_kind)  :: ntraero              ! number of aerosol tracers (up to max_aero in ice_domain_size.F90)
          integer (kind=int_kind)  :: trzaero              ! number of z aerosol tracers (up to max_aero = 6)
          integer (kind=int_kind)  :: tralg                ! number of algal tracers (up to max_algae = 3)
          integer (kind=int_kind)  :: trdoc                ! number of dissolve organic carbon (up to max_doc = 3)
          integer (kind=int_kind)  :: trdic                ! number of dissolve inorganic carbon (up to max_dic = 1)
          integer (kind=int_kind)  :: trdon                ! number of dissolve organic nitrogen (up to max_don = 1)
          integer (kind=int_kind)  :: trfed                ! number of dissolved iron tracers (up to max_fe  = 2)
          integer (kind=int_kind)  :: trfep                ! number of particulate iron tracers (up to max_fe  = 2)
          integer (kind=int_kind)  :: nbgclyr              ! number of zbgc layers
          integer (kind=int_kind)  :: trbgcz               ! set to 1 for zbgc tracers (needs TRBGCS = 0 and TRBRI = 1)
          integer (kind=int_kind)  :: trzs                 ! set to 1 for zsalinity tracer (needs TRBRI = 1)
          integer (kind=int_kind)  :: trbri                ! set to 1 for brine height tracer
          integer (kind=int_kind)  :: trage                ! set to 1 for ice age tracer
          integer (kind=int_kind)  :: trfy                 ! set to 1 for first-year ice area tracer
          integer (kind=int_kind)  :: trlvl                ! set to 1 for level and deformed ice tracers
          integer (kind=int_kind)  :: trpnd                ! set to 1 for melt pond tracers
          integer (kind=int_kind)  :: trbgcs               ! set to 1 for skeletal layer tracers (needs TRBGCZ = 0)

          ! setting variables used by the model

          integer (kind=int_kind)  :: ncat                 ! number of categories in use
          integer (kind=int_kind)  :: nfsd                 ! number of floe size categories in use
          integer (kind=int_kind)  :: nilyr                ! number of ice layers per category in use
          integer (kind=int_kind)  :: nslyr                ! number of snow layers per category in use
          integer (kind=int_kind)  :: n_aero               ! number of aerosols in use
          integer (kind=int_kind)  :: n_zaero              ! number of z aerosols in use
          integer (kind=int_kind)  :: n_algae              ! number of algae in use
          integer (kind=int_kind)  :: n_doc                ! number of DOC pools in use
          integer (kind=int_kind)  :: n_dic                ! number of DIC pools in use
          integer (kind=int_kind)  :: n_don                ! number of DON pools in use
          integer (kind=int_kind)  :: n_fed                ! number of Fe pools in use dissolved Fe
          integer (kind=int_kind)  :: n_fep                ! number of Fe pools in use particulate Fe
          integer (kind=int_kind)  :: nblyr                ! number of bio/brine layers per category
          integer (kind=int_kind)  :: n_bgc                ! nit, am, sil, dmspp, dmspd, dms, pon, humic
          integer (kind=int_kind)  :: nltrcr               ! number of zbgc (includes zaero) and zsalinity tracers
          integer (kind=int_kind)  :: max_nsw              ! number of tracers active in shortwave calculation
          integer (kind=int_kind)  :: max_ntrcr            ! number of tracers in total
          integer (kind=int_kind)  :: nfreq                ! number of wave frequencies ! HARDWIRED FOR NOW

          ! tracer namelist
        
          logical (kind=log_kind)  :: tr_iage
          logical (kind=log_kind)  :: tr_FY
          logical (kind=log_kind)  :: tr_lvl
          logical (kind=log_kind)  :: tr_pond_cesm
          logical (kind=log_kind)  :: tr_pond_topo
          logical (kind=log_kind)  :: tr_pond_lvl
          logical (kind=log_kind)  :: tr_aero
          logical (kind=log_kind)  :: tr_fsd

          ! thermo namelist

          integer (kind=int_kind)  :: kitd
          integer (kind=int_kind)  :: ktherm
          character (len=char_len) :: conduct
          real (kind=dbl_kind)     :: a_rapid_mode
          real (kind=dbl_kind)     :: Rac_rapid_mode  
          real (kind=dbl_kind)     :: aspect_rapid_mode 
          real (kind=dbl_kind)     :: dSdt_slow_mode    
          real (kind=dbl_kind)     :: phi_c_slow_mode   
          real (kind=dbl_kind)     :: phi_i_mushy       

          ! dynamics namelist
          
          integer (kind=int_kind)  :: kstrength       
          integer (kind=int_kind)  :: krdg_partic     
          integer (kind=int_kind)  :: krdg_redist     
          integer (kind=int_kind)  :: mu_rdg          
          real (kind=dbl_kind)     :: Cf              
          
          ! shortwave namelist
      
          character (len=char_len) :: shortwave      
          character (len=char_len) :: albedo_type    
          real (kind=dbl_kind)     :: albicev        
          real (kind=dbl_kind)     :: albicei         
          real (kind=dbl_kind)     :: albsnowv        
          real (kind=dbl_kind)     :: albsnowi      
          real (kind=dbl_kind)     :: ahmax         
          real (kind=dbl_kind)     :: R_ice         
          real (kind=dbl_kind)     :: R_pnd         
          real (kind=dbl_kind)     :: R_snw         
          real (kind=dbl_kind)     :: dT_mlt        
          real (kind=dbl_kind)     :: rsnw_mlt       
          real (kind=dbl_kind)     :: kalg           

          ! ponds namelist

          real (kind=dbl_kind)     :: hp1 
          real (kind=dbl_kind)     :: hs0 
          real (kind=dbl_kind)     :: hs1
          real (kind=dbl_kind)     :: dpscale         
          character (len=char_len) :: frzpnd         
          real (kind=dbl_kind)     :: rfracmin      
          real (kind=dbl_kind)     :: rfracmax     
          real (kind=dbl_kind)     :: pndaspect   

          ! forcing namelist

          logical (kind=log_kind)  :: formdrag        
          character (len=char_len) :: atmbndy         
          logical (kind=log_kind)  :: calc_strair     
          logical (kind=log_kind)  :: calc_Tsfc       
          logical (kind=log_kind)  :: highfreq        
          integer (kind=int_kind)  :: natmiter        
          real (kind=dbl_kind)     :: ustar_min      
          real (kind=dbl_kind)     :: emissivity      
          character (len=char_len) :: fbot_xfer_type  
          logical (kind=log_kind)  :: update_ocn_f    
          logical (kind=log_kind)  :: l_mpond_fresh   
          character (len=char_len) :: tfrz_option     
          logical (kind=log_kind)  :: oceanmixed_ice   
          character (len=char_len) :: wave_spec_type   

      end type t_icepack_settings

!=======================================================================    

      end module icedrv_settings

!=======================================================================

