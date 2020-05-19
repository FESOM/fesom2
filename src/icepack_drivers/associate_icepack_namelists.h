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
integer (kind=int_kind), pointer :: ndtd                 ! dynamic time steps per thermodynamic time step

! tracer namelist

logical (kind=log_kind), pointer :: tr_iage
logical (kind=log_kind), pointer :: tr_FY
logical (kind=log_kind), pointer :: tr_lvl
logical (kind=log_kind), pointer :: tr_pond_cesm
logical (kind=log_kind), pointer :: tr_pond_topo
logical (kind=log_kind), pointer :: tr_pond_lvl
logical (kind=log_kind), pointer :: tr_aero
logical (kind=log_kind), pointer :: tr_fsd

! grid namelist

integer (kind=int_kind), pointer  :: kcatbound

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
real (kind=dbl_kind),    pointer  :: mu_rdg
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

! pond namelist

real (kind=dbl_kind),     pointer :: hp1
real (kind=dbl_kind),     pointer :: hs0
real (kind=dbl_kind),     pointer :: hs1
real (kind=dbl_kind),     pointer :: dpscale
real (kind=dbl_kind),     pointer :: rfracmin
real (kind=dbl_kind),     pointer :: rfracmax
real (kind=dbl_kind),     pointer :: pndaspect
character (len=char_len), pointer :: frzpnd

!=====================================================
!=====================================================

nicecat => icepack_namelists%nicecat
nfsdcat => icepack_namelists%nfsdcat
nicelyr => icepack_namelists%nicelyr
nsnwlyr => icepack_namelists%nsnwlyr
ntraero => icepack_namelists%ntraero
trzaero => icepack_namelists%trzaero
tralg   => icepack_namelists%tralg
trdoc   => icepack_namelists%trdoc
trdic   => icepack_namelists%trdic
trdon   => icepack_namelists%trdon
trfed   => icepack_namelists%trfed
trfep   => icepack_namelists%trfep
nbgclyr => icepack_namelists%nbgclyr
trbgcz  => icepack_namelists%trbgcz 
trzs    => icepack_namelists%trzs
trbri   => icepack_namelists%trbri
trage   => icepack_namelists%trage
trfy    => icepack_namelists%trfy
trlvl   => icepack_namelists%trlvl
trpnd   => icepack_namelists%trpnd
trbgcs  => icepack_namelists%trbgcs

tr_iage      => icepack_namelists%tr_iage
tr_FY        => icepack_namelists%tr_FY
tr_lvl       => icepack_namelists%tr_lvl
tr_pond_cesm => icepack_namelists%tr_pond_cesm
tr_pond_topo => icepack_namelists%tr_pond_topo
tr_pond_lvl  => icepack_namelists%tr_pond_lvl
tr_aero      => icepack_namelists%tr_aero
tr_fsd       => icepack_namelists%tr_fsd

kcatbound    => icepack_namelists%kcatbound

kitd               => icepack_namelists%kitd
ktherm             => icepack_namelists%ktherm
conduct            => icepack_namelists%conduct
a_rapid_mode       => icepack_namelists%a_rapid_mode
Rac_rapid_mode     => icepack_namelists%Rac_rapid_mode
aspect_rapid_mode  => icepack_namelists%aspect_rapid_mode
dSdt_slow_mode     => icepack_namelists%dSdt_slow_mode
phi_c_slow_mode    => icepack_namelists%phi_c_slow_mode
phi_i_mushy        => icepack_namelists%phi_i_mushy

kstrength   => icepack_namelists%kstrength
krdg_partic => icepack_namelists%krdg_partic
krdg_redist => icepack_namelists%krdg_redist
mu_rdg      => icepack_namelists%mu_rdg
Cf          => icepack_namelists%Cf

shortwave   => icepack_namelists%shortwave
albedo_type => icepack_namelists%albedo_type
albicev     => icepack_namelists%albicev
albicei     => icepack_namelists%albicei
albsnowv    => icepack_namelists%albsnowv
albsnowi    => icepack_namelists%albsnowi
ahmax       => icepack_namelists%ahmax
R_ice       => icepack_namelists%R_ice
R_pnd       => icepack_namelists%R_pnd
R_snw       => icepack_namelists%R_snw
dT_mlt      => icepack_namelists%dT_mlt
rsnw_mlt    => icepack_namelists%rsnw_mlt
kalg        => icepack_namelists%kalg

formdrag       => icepack_namelists%formdrag
atmbndy        => icepack_namelists%atmbndy
calc_strair    => icepack_namelists%calc_strair
calc_Tsfc      => icepack_namelists%calc_Tsfc
highfreq       => icepack_namelists%highfreq
natmiter       => icepack_namelists%natmiter
ustar_min      => icepack_namelists%ustar_min
emissivity     => icepack_namelists%emissivity
fbot_xfer_type => icepack_namelists%fbot_xfer_type
update_ocn_f   => icepack_namelists%update_ocn_f
l_mpond_fresh  => icepack_namelists%l_mpond_fresh
tfrz_option    => icepack_namelists%tfrz_option
oceanmixed_ice => icepack_namelists%oceanmixed_ice
wave_spec_type => icepack_namelists%wave_spec_type

hp1            => icepack_namelists%hp1
hs0            => icepack_namelists%hs0
hs1            => icepack_namelists%hs1
dpscale        => icepack_namelists%dpscale
rfracmin       => icepack_namelists%rfracmin
rfracmax       => icepack_namelists%rfracmax
pndaspect      => icepack_namelists%pndaspect
frzpnd         => icepack_namelists%frzpnd
