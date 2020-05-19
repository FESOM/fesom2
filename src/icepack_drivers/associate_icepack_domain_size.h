! setting variables used by the model

integer (kind=int_kind), pointer :: nx                   ! number of grid cells and gost cells for each mesh partition            
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

!=====================================================
!=====================================================

nx      => icepack_domain_size%nx
ncat    => icepack_domain_size%ncat
nfsd    => icepack_domain_size%nfsd
nilyr   => icepack_domain_size%nilyr
nslyr   => icepack_domain_size%nslyr
n_aero  => icepack_domain_size%n_aero
n_zaero => icepack_domain_size%n_zaero
n_algae => icepack_domain_size%n_algae
n_doc   => icepack_domain_size%n_doc
n_dic   => icepack_domain_size%n_dic
n_don   => icepack_domain_size%n_don
n_fed   => icepack_domain_size%n_fed
n_fep   => icepack_domain_size%n_fep
nblyr   => icepack_domain_size%nblyr
n_bgc   => icepack_domain_size%n_bgc
nltrcr  => icepack_domain_size%nltrcr
max_nsw => icepack_domain_size%max_nsw
nfreq   => icepack_domain_size%nfreq
max_ntrcr => icepack_domain_size%max_ntrcr

