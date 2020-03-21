integer (kind=int_kind), pointer :: nicecat              ! number of ice thickness categories
integer (kind=int_kind), pointer :: nfsdcat              ! number of floe size categories
integer (kind=int_kind), pointer :: nicelyr              ! number of vertical layers in the ice
integer (kind=int_kind), pointer :: nsnwlyr              ! number of vertical layers in the snow
integer (kind=int_kind), pointer :: ntraero              ! number of aerosol tracers (up to max_aero in ice_domain_size.F90)
integer (kind=int_kind), pointer :: trzaero              ! number of z aerosol tracers (up to max_aero = 6)
integer (kind=int_kind), pointer :: tralg                ! number of algal tracers (up to max_algae = 3)
integer (kind=int_kind), pointer :: trdoc                ! number of dissolve organic carbon (up to max_doc = 3)
integer (kind=int_kind), pointer :: trdic                ! number of dissolve inorganic carbon (up to max_dic = 1)
integer (kind=int_kind), pointer :: trdon                ! number of dissolve organic nitrogen (up to max_don = 1)
integer (kind=int_kind), pointer :: trfed                ! number of dissolved iron tracers (up to max_fe  = 2)
integer (kind=int_kind), pointer :: trfep                ! number of particulate iron tracers (up to max_fe  = 2)
integer (kind=int_kind), pointer :: nbgclyr              ! number of zbgc layers
integer (kind=int_kind), pointer :: trbgcz               ! set to 1 for zbgc tracers (needs TRBGCS = 0 and TRBRI = 1)
integer (kind=int_kind), pointer :: trzs                 ! set to 1 for zsalinity tracer (needs TRBRI = 1)
integer (kind=int_kind), pointer :: trbri                ! set to 1 for brine height tracer
integer (kind=int_kind), pointer :: trage                ! set to 1 for ice age tracer
integer (kind=int_kind), pointer :: trfy                 ! set to 1 for first-year ice area tracer
integer (kind=int_kind), pointer :: trlvl                ! set to 1 for level and deformed ice tracers
integer (kind=int_kind), pointer :: trpnd                ! set to 1 for melt pond tracers
integer (kind=int_kind), pointer :: trbgcs               ! set to 1 for skeletal layer tracers (needs TRBGCZ = 0)


integer (kind=int_kind), pointer :: ncat                 ! number of categories in use
integer (kind=int_kind), pointer :: nfsd                 ! number of floe size categories in use
integer (kind=int_kind), pointer :: nilyr                ! number of ice layers per category in use
integer (kind=int_kind), pointer :: nslyr                ! number of snow layers per category in use
integer (kind=int_kind), pointer :: n_aero               ! number of aerosols in use
integer (kind=int_kind), pointer :: n_zaero              ! number of z aerosols in use 
integer (kind=int_kind), pointer :: n_algae              ! number of algae in use 
integer (kind=int_kind), pointer :: n_doc                ! number of DOC pools in use
integer (kind=int_kind), pointer :: n_dic                ! number of DIC pools in use
integer (kind=int_kind), pointer :: n_don                ! number of DON pools in use
integer (kind=int_kind), pointer :: n_fed                ! number of Fe  pools in use dissolved Fe
integer (kind=int_kind), pointer :: n_fep                ! number of Fe  pools in use particulate Fe
integer (kind=int_kind), pointer :: nblyr                ! number of bio/brine layers per category
integer (kind=int_kind), pointer :: n_bgc                ! nit, am, sil, dmspp, dmspd, dms, pon, humic
integer (kind=int_kind), pointer :: nltrcr               ! number of zbgc (includes zaero) and zsalinity tracers
integer (kind=int_kind), pointer :: max_nsw              ! number of tracers active in shortwave calculation
integer (kind=int_kind), pointer :: max_ntrcr            ! number of tracers in total
integer (kind=int_kind), pointer :: nfreq                ! number of wave frequencies ! HARDWIRED FOR NOW

nicecat => icepack_settings%nicecat
nfsdcat => icepack_settings%nfsdcat
nicelyr => icepack_settings%nicelyr
nsnwlyr => icepack_settings%nsnwlyr
ntraero => icepack_settings%ntraero
trzaero => icepack_settings%trzaero
tralg   => icepack_settings%tralg
trdoc   => icepack_settings%trdoc
trdic   => icepack_settings%trdic
trdon   => icepack_settings%trdon
trfed   => icepack_settings%trfed
trfep   => icepack_settings%trfep
nbgclyr => icepack_settings%nbgclyr
trbgcz  => icepack_settings%trbgcz 
trzs    => icepack_settings%trzs
trbri   => icepack_settings%trbri
trage   => icepack_settings%trage
trfy    => icepack_settings%trfy
trlvl   => icepack_settings%trlvl
trpnd   => icepack_settings%trpnd
trbgcs  => icepack_settings%trbgcs

ncat    => icepack_settings%ncat
nfsd    => icepack_settings%nfsd
nilyr   => icepack_settings%nilyr
nslyr   => icepack_settings%nslyr
n_aero  => icepack_settings%n_aero
n_zaero => icepack_settings%n_zaero
n_algae => icepack_settings%n_algae
n_doc   => icepack_settings%n_doc
n_dic   => icepack_settings%n_dic
n_don   => icepack_settings%n_don
n_fed   => icepack_settings%n_fed
n_fep   => icepack_settings%n_fep
nblyr   => icepack_settings%nblyr
n_bgc   => icepack_settings%n_bgc
nltrcr  => icepack_settings%nltrcr
max_nsw => icepack_settings%max_nsw
nfreq   => icepack_settings%nfreq
max_ntrcr => icepack_settings%max_ntrcr






























