! env namelist

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

! setting variables used by the model

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

! tracer namelist

logical (kind=log_kind), pointer :: tr_iage
logical (kind=log_kind), pointer :: tr_FY
logical (kind=log_kind), pointer :: tr_lvl
logical (kind=log_kind), pointer :: tr_pond_cesm
logical (kind=log_kind), pointer :: tr_pond_topo
logical (kind=log_kind), pointer :: tr_pond_lvl
logical (kind=log_kind), pointer :: tr_aero
logical (kind=log_kind), pointer :: tr_fsd

! thermo namelist

integer (kind=int_kind),  pointer :: kitd
integer (kind=int_kind),  pointer :: ktherm
character (len=char_len), pointer :: conduct
real (kind=dbl_kind),     pointer :: a_rapid_mode
real (kind=dbl_kind),     pointer :: Rac_rapid_mode
real (kind=dbl_kind),     pointer :: aspect_rapid_mode
real (kind=dbl_kind),     pointer :: dSdt_slow_mode
real (kind=dbl_kind),     pointer :: phi_c_slow_mode
real (kind=dbl_kind),     pointer :: phi_i_mushy

! dynamics namelist

integer (kind=int_kind), pointer  :: kstrength
integer (kind=int_kind), pointer  :: krdg_partic
integer (kind=int_kind), pointer  :: krdg_redist
integer (kind=int_kind), pointer  :: mu_rdg
real (kind=dbl_kind),    pointer  :: Cf

! shortwave namelist

character (len=char_len), pointer :: shortwave
character (len=char_len), pointer :: albedo_type
real (kind=dbl_kind),     pointer :: albicev
real (kind=dbl_kind),     pointer :: albicei
real (kind=dbl_kind),     pointer :: albsnowv
real (kind=dbl_kind),     pointer :: albsnowi
real (kind=dbl_kind),     pointer :: ahmax
real (kind=dbl_kind),     pointer :: R_ice
real (kind=dbl_kind),     pointer :: R_pnd
real (kind=dbl_kind),     pointer :: R_snw
real (kind=dbl_kind),     pointer :: dT_mlt
real (kind=dbl_kind),     pointer :: rsnw_mlt
real (kind=dbl_kind),     pointer :: kalg

! forcing namelist

logical (kind=log_kind),  pointer :: formdrag
character (len=char_len), pointer :: atmbndy
logical (kind=log_kind),  pointer :: calc_strair
logical (kind=log_kind),  pointer :: calc_Tsfc
logical (kind=log_kind),  pointer :: highfreq
integer (kind=int_kind),  pointer :: natmiter
real (kind=dbl_kind),     pointer :: ustar_min
real (kind=dbl_kind),     pointer :: emissivity
character (len=char_len), pointer :: fbot_xfer_type
logical (kind=log_kind),  pointer :: update_ocn_f
logical (kind=log_kind),  pointer :: l_mpond_fresh
character (len=char_len), pointer :: tfrz_option
logical (kind=log_kind),  pointer :: oceanmixed_ice
character (len=char_len), pointer :: wave_spec_type

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

tr_iage      => icepack_settings%tr_iage
tr_FY        => icepack_settings%tr_FY
tr_lvl       => icepack_settings%tr_lvl
tr_pond_cesm => icepack_settings%tr_pond_cesm
tr_pond_topo => icepack_settings%tr_pond_topo
tr_pond_lvl  => icepack_settings%tr_pond_lvl
tr_aero      => icepack_settings%tr_aero
tr_fsd       => icepack_settings%tr_fsd

kitd               => icepack_settings%kitd
ktherm             => icepack_settings%ktherm
conduct            => icepack_settings%conduct
a_rapid_mode       => icepack_settings%a_rapid_mode
Rac_rapid_mode     => icepack_settings%Rac_rapid_mode
aspect_rapid_mode  => icepack_settings%aspect_rapid_mode
dSdt_slow_mode     => icepack_settings%dSdt_slow_mode
phi_c_slow_mode    => icepack_settings%phi_c_slow_mode
phi_i_mushy        => icepack_settings%phi_i_mushy

kstrength   => icepack_settings%kstrength
krdg_partic => icepack_settings%krdg_partic
krdg_redist => icepack_settings%krdg_redist
mu_rdg      => icepack_settings%mu_rdg
Cf          => icepack_settings%Cf

shortwave   => icepack_settings%shortwave
albedo_type => icepack_settings%albedo_type
albicev     => icepack_settings%albicev
albicei     => icepack_settings%albicei
albsnowv    => icepack_settings%albsnowv
albsnowi    => icepack_settings%albsnowi
ahmax       => icepack_settings%ahmax
R_ice       => icepack_settings%R_ice
R_pnd       => icepack_settings%R_pnd
R_snw       => icepack_settings%R_snw
dT_mlt      => icepack_settings%dT_mlt
rsnw_mlt    => icepack_settings%rsnw_mlt
kalg        => icepack_settings%kalg

formdrag       => icepack_settings%formdrag
atmbndy        => icepack_settings%atmbndy
calc_strair    => icepack_settings%calc_strair
calc_Tsfc      => icepack_settings%calc_Tsfc
highfreq       => icepack_settings%highfreq
natmiter       => icepack_settings%natmiter
ustar_min      => icepack_settings%ustar_min
emissivity     => icepack_settings%emissivity
fbot_xfer_type => icepack_settings%fbot_xfer_type
update_ocn_f   => icepack_settings%update_ocn_f
l_mpond_fresh  => icepack_settings%l_mpond_fresh
tfrz_option    => icepack_settings%tfrz_option
oceanmixed_ice => icepack_settings%oceanmixed_ice
wave_spec_type => icepack_settings%wave_spec_type



























