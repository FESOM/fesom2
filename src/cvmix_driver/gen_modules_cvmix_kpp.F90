!
!
!===============================================================================
! @see  Large W.G., McWilliams J.C., Doney S.C. -> KPD94
!       Oceanic Vertical Mixing: A Review and a Model with a Nonlocal
!       Boundary Layer Parameterizations.
!       Rev. of Geophys., XX,XX–XX. doi: , 1994.
! written by Patrick Scholz, 14.05.2019  --> based on the MPIOM-CVMIX interface 
! provided by Nils Brüggeman, Oliver Gutjahr
module g_cvmix_kpp
    !___________________________________________________________________________
    ! module calls from cvmix library
    use cvmix_kinds_and_types
    
    USE cvmix_kpp, only : cvmix_init_kpp, cvmix_put_kpp, CVmix_get_kpp_real, &
                          cvmix_coeffs_kpp, cvmix_kpp_compute_OBL_depth,     &
                          cvmix_kpp_compute_turbulent_scales,                &   
                          cvmix_kpp_compute_bulk_Richardson,                 & 
                          cvmix_kpp_compute_unresolved_shear,                & 
                          cvmix_kpp_compute_kOBL_depth,                      &
                          cvmix_kpp_compute_StokesXi,                        &
                          cvmix_kpp_ustokes_SL_model
    
    !___________________________________________________________________________
    ! module calls from FESOM
    use g_config
    use o_param           
    USE MOD_ICE
    USE MOD_DYN
    USE mod_tracer
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_arrays
    use g_comm_auto 
    use g_forcing_arrays
    use g_support
    use o_mixing_KPP_mod
    implicit none
    
    !___Parameter for the init of KPP___________________________________________
    ! kpp_use_fesomkpp ... general flag to use fesom1.4 flavour of KPP mixing 
    ! which corresbonds to MOM5 if false use MOM6 flavour of KPP mixing 
    
    logical           :: kpp_use_fesomkpp = .false. 
    
    ! Critical bulk Richardson number,(OBL_depth = where bulk Ri = Ri_crit)
    real(kind=WP)     :: kpp_Rib_crit=0.3  
    
    ! von Karman constant (dimensionless)
    real(kind=WP)     :: kpp_vonKarman=0.40   
    
    ! If non-zero, set the minimum depth for the OBL (m)
    real(kind=WP)     :: kpp_minOBLdepth=0.0
    
    ! Min for the squared unresolved velocity used in Rib CVMix calculation (m2/s2)
    real(kind=WP)     :: kpp_minVtsqr=1.0e-10
    
    ! Fraction of OBL depth considered in the surface layer (nondim) (epsilon=0.1)
    real(kind=WP)     :: kpp_surf_layer_ext=0.10
    
    ! If non-zero, is a distance from the bottom that the OBL can not penetrate through (m)
    real(kind=WP)     :: kpp_deepOBLoffset = 0.0
    
    ! Parameter for multiplying by non-local term
    real(kind=WP)     :: kpp_cs2=6.32739901508
    
    ! By default, cvmix use sigma construction from Danabasoglu et al. when computing
    ! turbulent scales. Set l_LMD_ws = .true. to use Large et al. construction.
    logical           :: kpp_use_LMDws = .true.
    
    ! If true, add enhanced diffusivity at base of boundary layer (in the original
    logical           :: kpp_use_enhanceKv = .true. 
    
    ! If true, compute Ekman depth limit for OBLdepth 
    logical           :: kpp_use_compEkman = .true.  
    
    ! If true, compute Monin-Obukhov limit for OBLdepth
    logical           :: kpp_use_monob     = .true.
        
    logical           :: kpp_cs_is_one     = .false.                                !
    
    ! If true use horizontal smoothing on boundary layer mixing coefficient as 
    ! its done in the original kpp of fesom2.0
    logical           :: kpp_use_smoothblmc= .true.
    
    ! how offen should smoothing be applied
    integer           :: kpp_smoothblmc_nmb= 3
    
    ! Sets MOM6 method for using shortwave radiation in surface buoyancy flux
    ! --> all   ... all shortwave pentration is used
    ! --> mxl   ... only shortwave pentration with obldepth is used -->  comes 
    !               closest what FESOM1.4 was doing
    ! --> lvl1  ... only shortwave pentration of first layer is used
    ! --> fesom ... is not in MOM6 but uses the method from fesom1.4
    character(len=20) :: kpp_sw_method ="mxl"    
 
    ! interpolation type used to interpolate bulk Richardson number
    ! --> determining OBL depth: linear,quadratic,cubic
    ! --> tests have shown that "linear" and "cubic" lead to very similar Kv profiles,
    !     with "linear" producing slightly higher max. Kv values (+3%) compared to "cubic".
    !     "quad" produces an up to 20% weaker maximum diffusivity with an up to 500m
    !     shallower Kv profile. Using "cubic" together with core2 grid leads to
    !     slightly colder Arctic Ocean (~0.2°C) when compared to "linear". Using 
    !     "linear" reduces the cold blop anomaly by around 0.8°C  compared to "cubic" 
    character(len=20) :: kpp_interptype_ri = "linear"             
    
    ! interpolation type used to interpolate diff and visc at OBL_depth: 
    ! linear,quadratic,cubic,LMD94 (default for cvmix)
    character(len=20) :: kpp_interptype_atobl = "LMD94"             
    
    ! Method used in CVMix for setting diffusivity and NLT profile functions:
    ! SimpleShapes      = sigma*(1-sigma)^2 for both diffusivity and NLT, Shape 
    !                     functions for both the gradient and nonlocal term vanish
    !                     at interface  
    ! MatchGradient     = Shape functions for nonlocal term vanishes at interface, 
    !                     but gradient term matches interior values, sigma*(1-sigma)^2 
    !                     for NLT; diffusivity profile from matching
    ! MatchBoth         = Shape function for both the gradient and nonlocal terms
    !                     match interior values
    ! ParabolicNonLocal = Shape function for the nonlocal term is (1-sigma)^2, 
    !                     gradient term is sigma*(1-sigma)^2
    ! --> MOM6 recommends parabolic shape function for the nonlocal transport 
    ! --> tests have shown that "ParabolicNonLocal" leads to deepest and highest Kv 
    !     values, "SimpleShapes" lead to intermediate Kv values but with shallowest 
    !     extension, "Matchboth" lead to weakest maximum Kv values with intermediate
    !     vertical extension and  "MatchGradient" leads to intermediate Kv values and 
    !     intermediate vertical extension
    character(len=20) :: kpp_matchtechc = "ParabolicNonLocal"  
    
    ! MOM6 over-ride of CVMix Non Local Transport (NLT) shape function if not 
    ! NLT_shape="cvmix" Over-rides the result from CVMix.  Allowed values are:
    ! cvmix     - Uses the profiles from CVmix specified by MATCH_TECHNIQUE (default)
    ! linear    - A linear profile, 1-sigma
    ! parabolic - A parablic profile, (1-sigma)^2
    ! cubic     - A cubic profile, (1-sigma)^2(1+2*sigma)
    ! cubic_LMD - The original KPP profile
    character(len=20) :: kpp_nlt_shape = "cvmix"  
        
    ! Ri-number dependet mixing scheme below the OBL: 'PP' or 'KPP'
    character(len=10) :: kpp_internalmix= "KPP"      

    !___Windstress below ice____________________________________________________
    ! If True, reduce the wind stress (ustar) under sea ice. --> took this from 
    ! mpiom not tested yet
    logical           :: kpp_reduce_tauuice = .false. 
    
    !___Stokes Similarty package________________________________________________
    ! If true, use Stokes Similarty package (i.e. include wave‐related / Stokes drift 
    ! effects in the surface layer). Triggers usage of additional routines 
    ! that alter the shape functions, or mixing formulations, incorporating wave / 
    ! Stokes drift effects consistent with Monin–Obukhov similarity theory 
    ! (MOST). The code logic probably augments or replaces parts of the standard boundary 
    ! (layer similarity (or nonlocal mixing) using a Stokes‐drift‐aware correction.
    logical           :: kpp_use_StokesMOST= .true.
    
    ! approximate proportionality between surface wind velocity and stokes velocity
    ! U_stokes ~ kpp_A_stokes * U_wind
    real(kind=WP)     :: kpp_A_stokes      = 0.005 ! a
    
    !___Langmuir option____________________________________________________
    ! Option of Langmuir enhanced mixing apply an enhancement factor to the
    ! turbulent velocity scale
    ! LWF16     -  MixingCoefEnhancement = Langmuir_EFactor
    ! RWHGK16   -  MixingCoefEnhancement = cvmix_one + ShapeNoMatchAtS/NMshapeMax * &
    !                                      (Langmuir_EFactor - cvmix_one)
    ! NONE      -  Langmuir switched off, MixingCoefEnhancement=1 
    character(len=20) :: kpp_langmuir_mixing= "LWF16"
        
    ! Option of Langmuir turbulence enhanced entrainment - modify the unresolved shear
    ! LWF16     -  Li Q., Webb A., Fox-Kemper B., Craig A., Danabasoglu G., 
    !              Large W., Vertenstein M., 2016, Langmuir mixing effects on 
    !              global climate: WAVEWATCH III in CESM, Ocean Modelling 103 (2016) 145–160
    !              
    ! LF17      -  Li Q., Fox-Kemper B., Breivik O., Webb A., 2017, Statistical 
    !              models of global Langmuir mixing, Ocean Modelling 113 (2017) 95–114
    !              
    ! RWHGK16   -  Reichl B., Wang D., Hara T., Ginis I. and Kukulka T, 2016, Impact 
    !              of Sea-State-Dependent Langmuir Turbulence on the Ocean
    !              Response to a Tropical Cyclone, Mon. Wea. Rev., 144
    !              
    ! NONE      -  
    character(len=20) :: kpp_langmuir_entrainment= "LF17"
    
    !___Mixing below OBL________________________________________________________
    ! Parameters to run shear-dependent LM94 scheme below the mixed layer
    ! leading coefficient of shear mixing formula, units: m^2/s: default= 5e-3  
    real(kind=WP)     :: kpp_Av0 = 5.0e-3 
    real(kind=WP)     :: kpp_Kv0 = 5.0e-3 
    
    ! critical shear  Richardson number value, units: unitless (0.7 in LMD94)
    real(kind=WP)     :: kpp_Ri0 = 0.7
    
    ! Exponent of unitless factor of diffusities,units:unitless (3 in LMD94)
    real(kind=WP)     :: kpp_loc_exp = 3.0
    
    ! Parameter in case of PP mixing below the OBL
    real(kind=WP)     :: kpp_pp_Av0     = 0.01
    real(kind=WP)     :: kpp_pp_alpha   = 5.0
    real(kind=WP)     :: kpp_pp_loc_exp = 2.0

    !___Background Viscosity/Diffusivity________________________________________
    ! If True use non constant background diffusivity of Qiang from FESOM1.4
    logical           :: kpp_use_nonconstKvb = .true.
    
    ! Values for const. background viscosity and diffusivity
    real(kind=WP)     :: kpp_Avbckg     = 1.0e-4
    real(kind=WP)     :: kpp_Kvbckg     = 1.0e-5
    
    real(kind=WP), parameter :: kpp_epsln   = 1.0e-40_WP ! a small value
    
    !___________________________________________________________________________
    ! set which are namelist parameter
    namelist /param_kpp/ kpp_Rib_crit, kpp_vonKarman, kpp_surf_layer_ext,        &
                         kpp_minVtsqr, kpp_use_enhanceKv, kpp_use_compEkman,     &
                         kpp_use_monob, kpp_cs_is_one, kpp_interptype_ri,        &
                         kpp_interptype_atobl, kpp_minOBLdepth, kpp_matchtechc,  &
                         kpp_Avbckg, kpp_Kvbckg, kpp_internalmix, kpp_Av0,       &
                         kpp_Kv0, kpp_Ri0, kpp_loc_exp, kpp_pp_Av0,              & 
                         kpp_use_nonconstKvb, kpp_pp_alpha, kpp_pp_loc_exp,      &
                         kpp_reduce_tauuice, kpp_use_smoothblmc,                 &
                         kpp_smoothblmc_nmb, kpp_use_fesomkpp, kpp_deepOBLoffset,&
                         kpp_langmuir_mixing, kpp_langmuir_entrainment,          &
                         kpp_use_StokesMOST, kpp_A_stokes,                       &
                         kpp_use_LMDws, kpp_sw_method, kpp_nlt_shape
    
    !___________________________________________________________________________
    ! 1d arrays
    !-----------
    ! bulk richardson number
    real(kind=WP), allocatable, dimension(:)    :: kpp_bulkRi
    ! shear richardson number
    real(kind=WP), allocatable, dimension(:)    :: kpp_shearRi
    ! squared velocity shear referenced to the surface
    real(kind=WP), allocatable, dimension(:)    :: kpp_dvsurf2
    ! buoyancy difference between surface and layer
    real(kind=WP), allocatable, dimension(:)    :: kpp_dbsurf
    ! turbulent velocity scale for scalars (tracer) and momentum
    real(kind=WP), allocatable, dimension(:)    :: kpp_ws_cntr

    ! 2d arrays
    !-----------
    ! depth of oceananic boundary layer (OBL)
    real(kind=WP), allocatable, dimension(:)    :: kpp_obldepth, kpp_nzobldepth
    ! surface buoyancy flux
    real(kind=WP), allocatable, dimension(:)    :: kpp_sbuoyflx, kpp_buoyflx_nl
    
    ! stokes, langmuir related varaibles 
    real(kind=WP), allocatable, dimension(:)    :: kpp_stokesXi_z, kpp_stokesVt_z, kpp_stokesXi, &
                                                   kpp_EFactor, kpp_LaSL
                                                   
    real(kind=WP), allocatable, dimension(:)    :: kpp_uS_t, kpp_vS_t, &
                                                   kpp_uS_c, kpp_vS_c, &
                                                   kpp_uS_m, kpp_vS_m
    
    
    ! 3d arrays
    !-----------
    ! nodal viscosities/diffusivities
    real(kind=WP), allocatable, dimension(:,:)  :: kpp_Av, kpp_Kv
    ! non local transport for temp. and salt
    real(kind=WP), allocatable, dimension(:,:)  :: kpp_nonlcltranspT, kpp_nonlcltranspS
    ! obl mixing coefficient for momentum, temp/salt diffusivity
    real(kind=WP), allocatable, dimension(:,:,:):: kpp_oblmixc
    
    
    type(cvmix_global_params_type)              :: CVmix_params_in ! --> is neede for CVmix_params_in%Gravity
    
    contains
    !
    !
    !
    !===========================================================================
    ! allocate and initialize CVMIX KPP variables --> call initialisation 
    ! routine from cvmix library
    subroutine init_cvmix_kpp(partit, mesh)
        implicit none
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(inout), target :: partit
        character(len=MAX_PATH) :: nmlfile
        logical            :: nmlfile_exist=.False.
        integer            :: node_size
        integer fileunit
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        !_______________________________________________________________________
        if(mype==0) then
            write(*,*) '____________________________________________________________'
            write(*,*) ' --> initialise CVMIX_KPP'
            write(*,*)
        end if
        
        !_______________________________________________________________________
        ! allocate + initialse kpp arrays --> with size myDim_nod2D+eDim_nod2D
        node_size=myDim_nod2D+eDim_nod2D
        
        ! allocate 1D depth variable 
        allocate(kpp_bulkRi(nl))
        allocate(kpp_shearRi(nl))
        allocate(kpp_dvsurf2(nl-1))
        allocate(kpp_dbsurf(nl-1))
        allocate(kpp_ws_cntr(nl-1))
        !!PS allocate(kpp_w_cntr(nl-1))
        kpp_bulkRi        = 0.0_WP
        kpp_shearRi       = 0.0_WP
        kpp_dvsurf2       = 0.0_WP
        kpp_dbsurf        = 0.0_WP
        kpp_ws_cntr       = 0.0_WP
        !!PS kpp_wm_cntr       = 0.0_WP
        
        ! allocate horizontal 2D variable 
        allocate(kpp_obldepth(node_size),kpp_nzobldepth(node_size))
        allocate(kpp_sbuoyflx(node_size))
        kpp_obldepth      = 0.0_WP
        kpp_nzobldepth    = 0.0_WP
        kpp_sbuoyflx      = 0.0_WP          ! total surface buoyancy flux (includes solar + non-solar)
        
        allocate(kpp_buoyflx_nl(nl))
        kpp_buoyflx_nl    =0.0_WP
        
        
        ! allocate 3D variable 
        allocate(kpp_Av(nl,node_size),kpp_Kv(nl,node_size))
        allocate(kpp_nonlcltranspT(nl,node_size),kpp_nonlcltranspS(nl,node_size))
        kpp_Av            = 0.0_WP
        kpp_Kv            = 0.0_WP
        kpp_nonlcltranspT = 0.0_WP
        kpp_nonlcltranspS = 0.0_WP
        
        allocate(kpp_oblmixc(nl,node_size,3))
        kpp_oblmixc       = 0.0_WP
        
        allocate(kpp_uS_t(nl-1), kpp_vS_t(nl-1)) ! top
        allocate(kpp_uS_c(nl-1), kpp_vS_c(nl-1)) ! center   
        allocate(kpp_uS_m(nl-1), kpp_vS_m(nl-1)) ! mean
        kpp_uS_t         = 0.0_WP
        kpp_vS_t         = 0.0_WP
        kpp_uS_c         = 0.0_WP
        kpp_vS_c         = 0.0_WP
        kpp_uS_m         = 0.0_WP
        kpp_vS_m         = 0.0_WP

        allocate(kpp_stokesXi_z(nl-1))
        allocate(kpp_stokesVt_z(nl-1))
        kpp_stokesXi_z = 0.0_WP
        kpp_stokesVt_z = 0.0_WP
        
        allocate(kpp_EFactor(node_size))
        allocate(kpp_LaSL(node_size))
        allocate(kpp_stokesXi(node_size))
        kpp_EFactor       = 0.0_WP ! Langmuir enhancement factor for entrainment 
        kpp_LaSL          = 0.0_WP ! surface layer averaged Langmuir number (units: none)
        kpp_stokesXi   = 0.0_WP
        
        
        !_______________________________________________________________________
        ! read cvmix namelist file 
        nmlfile ='namelist.cvmix'    ! name of ocean namelist file
        ! check if cvmix namelist file exists if not use default values 
        inquire(file=trim(nmlfile),exist=nmlfile_exist) 
        if (nmlfile_exist) then
            open(newunit=fileunit,file=trim(nmlfile))
                read(fileunit,nml=param_kpp)
            close(fileunit)
        else
            write(*,*) '     could not find namelist.cvmix, will use default values !'
        end if
        
        !_______________________________________________________________________
        ! write info to log file 
        if (mype==0) then
            write(*,*) "     kpp_Rib_crit        = ", kpp_Rib_crit
            write(*,*) "     kpp_minOBLdepth     = ", kpp_minOBLdepth
            write(*,*) "     kpp_deepOBLoffset   = ", kpp_deepOBLoffset
            write(*,*) "     kpp_minVtsqr        = ", kpp_minVtsqr
            write(*,*) "     kpp_vonKarman       = ", kpp_vonKarman
            write(*,*) "     kpp_surf_layer_ext  = ", kpp_surf_layer_ext
            write(*,*) "     kpp_interptype_ri   = ", kpp_interptype_ri
            write(*,*) "     kpp_interptype_atobl= ", kpp_interptype_atobl
            write(*,*) "     kpp_use_LMDws       = ", kpp_use_LMDws
            write(*,*) "     kpp_sw_method       = ", kpp_sw_method
            write(*,*) "     kpp_nlt_shape       = ", kpp_nlt_shape
            write(*,*) "     kpp_use_compEkman   = ", kpp_use_compEkman
            write(*,*) "     kpp_use_monob       = ", kpp_use_monob
            write(*,*) "     kpp_matchtechc      = ", kpp_matchtechc
            write(*,*) "     kpp_use_enhanceKv   = ", kpp_use_enhanceKv
            write(*,*) "     kpp_cs_is_one       = ", kpp_cs_is_one
            write(*,*) "     kpp_interptype_ri   = ", kpp_interptype_ri
            write(*,*) "     kpp_matchtechc      = ", kpp_matchtechc
            write(*,*) "     kpp_internalmix     = ", kpp_internalmix
            write(*,*) "     kpp_reduce_tauuice  = ", kpp_reduce_tauuice
            write(*,*) "     kpp_langmuir_mixing = ", kpp_langmuir_mixing
            write(*,*) "     kpp_langmuir_entrainment = ", kpp_langmuir_entrainment
            if (kpp_internalmix .eq. 'KPP') then 
                write(*,*) "     kpp_Av0             = ", kpp_Av0
                write(*,*) "     kpp_Kv0             = ", kpp_Kv0
                write(*,*) "     kpp_Ri0             = ", kpp_Ri0
                write(*,*) "     kpp_loc_exp         = ", kpp_loc_exp
            else
                write(*,*) "     kpp_pp_Av0          = ", kpp_pp_Av0
                write(*,*) "     kpp_pp_alpha        = ", kpp_pp_alpha
                write(*,*) "     kpp_pp_loc_exp      = ", kpp_pp_loc_exp
            end if
            write(*,*) "     kpp_use_nonconstKvb = ", kpp_use_nonconstKvb
            if (kpp_use_nonconstKvb .eqv. .false.) then 
                write(*,*) "     kpp_Avbckg         = ", kpp_Avbckg
                write(*,*) "     kpp_Kvbckg         = ", kpp_Kvbckg
            end if
        end if
        
        !_______________________________________________________________________
        ! Initialise CVMIX
         ! call the cvmix subroutine to initialise all required namelists
        call cvmix_init_kpp(Ri_crit                  = kpp_Rib_crit,            & 
                            minOBLdepth              = kpp_minOBLdepth,         & 
                            minVtsqr                 = kpp_minVtsqr,            &
                            vonKarman                = kpp_vonKarman,           & 
                            surf_layer_ext           = kpp_surf_layer_ext,      &
                            interp_type              = kpp_interptype_ri,       &
                            interp_type2             = kpp_interptype_atobl,    &
                            lEkman                   = kpp_use_compEkman,       &
                            lStokesMOST              = kpp_use_StokesMOST,      &
                            lMonOb                   = kpp_use_monob,           &
                            MatchTechnique           = kpp_matchtechc,          &
                            lenhanced_diff           = kpp_use_enhanceKv,       &
                            lnonzero_surf_nonlocal   = kpp_cs_is_one,           &
                            Langmuir_mixing_str      = kpp_langmuir_mixing,     &
                            Langmuir_entrainment_str = kpp_langmuir_entrainment,&
                            l_LMD_ws                 = kpp_use_LMDws            & 
                            )
        
    end subroutine init_cvmix_kpp
    !
    !
    !
    !===========================================================================
    ! calculate PP vertrical mixing coefficients from CVMIX library
    subroutine calc_cvmix_kpp(ice, dynamics, tracers, partit, mesh)
        type(t_ice)   , intent(in),    target :: ice
        type(t_dyn)   , intent(in),    target :: dynamics
        type(t_tracer), intent(in),    target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        !_______________________________________________________________________
        integer       :: node, elem, nz, nln, nun,  elnodes(3), aux_nz, sld_nz
        real(kind=WP) :: vshear2, dz2, aux, aux_wm(mesh%nl), aux_ws(mesh%nl)
        real(kind=WP) :: aux_coeff, sigma, stable
        real(kind=WP) :: ustar, aux_surfbuoyflx_nl(mesh%nl), wind_norm
        
        integer       :: nzsfc, nztmp
        real(kind=WP) :: sldepth, sfc_temp, sfc_salt, sfc_u, sfc_v, htot, delh, rho_sfc, rho_nz, oblguess
        real(kind=WP) :: rhopot, bulk_0, bulk_pz, bulk_pz2
        real(kind=WP) :: sfc_rhopot, sfc_bulk_0, sfc_bulk_pz, sfc_bulk_pz2
        real(kind=WP) :: uS_sld_t, vS_sld_t, uS_sld_c, vS_sld_c, uS_sld_m, vS_sld_m, uv_SLmean
         !_______________________________________________________________________
        ! pointer on necessary derived types
        real(kind=WP), dimension(:)    , pointer :: a_ice
        real(kind=WP), dimension(:,:)  , pointer :: temp, salt
        real(kind=WP), dimension(:,:,:), pointer :: UVnode
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
        temp=>tracers%data(1)%values(:,:)
        salt=>tracers%data(2)%values(:,:)
        UVnode=>dynamics%uvnode(:,:,:)
        a_ice        => ice%data(1)%values(:)
        
        !_______________________________________________________________________
        kpp_Av = 0.0_WP
        kpp_Kv = 0.0_WP
        do node= 1, myDim_nod2D
            !___________________________________________________________________
            nln = nlevels_nod2D(node)-1
            nun = ulevels_nod2D(node)
            
            !___________________________________________________________________
            ! initialide 1d arrays
            kpp_dvsurf2 = 0.0_WP
            kpp_dbsurf  = 0.0_WP
            kpp_ws_cntr = 0.0_WP
            kpp_bulkRi  = 0.0_WP
            kpp_shearRi = 0.0_WP
            
            !___________________________________________________________________
            !   ||    ||    ||    ||    ||    ||    ||    ||    ||    ||    ||  
            !  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_ 
            !  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  / 
            !   \/    \/    \/    \/    \/    \/    \/    \/    \/    \/    \/  
            !___________________________________________________________________
            ! 1) compute: - vertical shear with respect to surface and 
            !             - shear Richardson number
            !             - bouyancy flux with respect to surface
            !             - friction velocity (ustar) at surface (m/s)
            
            !___2D Quantities___________________________________________________
            ! calculate surface bouyancy flux after eq. A2c & A2d & A3b & A3d 
            ! in Large et al. 1994
            ! kpp_sbuoyflx includes:    Q_sensibel
            !                           Q_latent
            !                           Q_shortwave (surface part)
            !                           Q_longwave (surface part)
            !                           W_freshwater
            kpp_sbuoyflx(node) = -g * &
                                    (sw_alpha(nun,node)*heat_flux( node) / vcpw + &   !heat_flux & water_flux: positive up
                                     sw_beta( nun,node)*water_flux(node)*salt(nun,node))
            
            kpp_buoyflx_nl = 0.0_WP
            if (use_sw_pene) then
                ! coeffcient to transfer SW temp flux into buoyancy flux
                aux_coeff       = g*sw_alpha(nun,node)
                
                do nz = nun, nln
                    ! sw_3d is the temperature through the full depth levels into/
                    ! out off the tacervolume 
                    ! --> (sw_3d(1,node)-sw_3d(nz,node)) contains all the penetrated 
                    ! heatflux until level nz
                    kpp_buoyflx_nl(nz) = kpp_sbuoyflx(node)+aux_coeff*(sw_3d(nun,node)-sw_3d(nz+1,node))
                                                   ! look in oce_shortwave_pene.F90<--|
                                                   ! substract swsurf from surface 
                                                   ! heat flux 
                end do
            else
                kpp_buoyflx_nl(nun:nln) = kpp_sbuoyflx(node)
            end if  
            
            ! calculate friction velocity (ustar) at surface (m/s)
            ustar = sqrt( sqrt( stress_node_surf(1,node)**2 + stress_node_surf(2,node)**2 )*density_0_r ) ! @ the surface (eqn. 2)
            
            ! reduce friction velocity under ice --> approach take from the mpiom
            ! interface --> haven't been tested before with FESOM
            if (kpp_reduce_tauuice) then
                ustar = ustar*(1.0_WP-a_ice(node))**2
            end if 
            
            !___3D Quantities___________________________________________________
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> call shear variables'//achar(27)//'[0m'
            if (kpp_use_fesomkpp) then
                !!PS do nz=2, nln
                do nz=nun+1, nln
                    !___________________________________________________________
                    ! calculate squared velocity shear referenced to the surface
                    ! --> cvmix wants to have it with  respect to the midlevel rather than full levels
                    !!PS kpp_dvsurf2(nz) = ((UVnode(1,nz-1,node)+UVnode(1,nz,node))*0.5 - UVnode( 1,1,node) )**2 + &
                    !!PS                   ((UVnode(2,nz-1,node)+UVnode(2,nz,node))*0.5 - UVnode( 2,1,node) )**2
                    kpp_dvsurf2(nz) = ((UVnode(1,nz-1,node)+UVnode(1,nz,node))*0.5 - UVnode( 1,nun,node) )**2 + &
                                      ((UVnode(2,nz-1,node)+UVnode(2,nz,node))*0.5 - UVnode( 2,nun,node) )**2
                    !___________________________________________________________
                    ! calculate shear Richardson number Ri = N^2/(du/dz)^2
                    dz2     = (Z_3d_n( nz-1,node)-Z_3d_n( nz,node))**2
                    vshear2 = (UVnode(1,nz-1,node)-UVnode(1,nz,node))**2 + &
                              (UVnode(2,nz-1,node)-UVnode(2,nz,node))**2 
                    vshear2 = vshear2/dz2
                    kpp_shearRi(nz) = max(bvfreq(nz,node),0.0_WP)/(vshear2+kpp_epsln)
                    
                    !_______________________________________________________________
                    ! buoyancy difference with respect to the surface --> computed in
                    ! oce_ale_pressure_bf.F90 --> subroutine pressure_bv 
                    ! --> dbsfc(nz,node)
                    !!PS call densityJM_components(temp(1,node), salt(1,node), sfc_bulk_0, sfc_bulk_pz, sfc_bulk_pz2, sfc_rhopot, mesh)
                    call densityJM_components(temp(nun,node), salt(nun,node), sfc_bulk_0, sfc_bulk_pz, sfc_bulk_pz2, sfc_rhopot, mesh)
                    call densityJM_components(temp(nz, node), salt(nz, node), bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)                    
                    rho_nz  = bulk_0   + Z_3d_n(nz,node)*(bulk_pz   + Z_3d_n(nz,node)*bulk_pz2)
                    rho_nz  = rho_nz*rhopot/(rho_nz+0.1_WP*Z_3d_n(nz,node))-density_0
                    rho_sfc = sfc_bulk_0   + Z_3d_n(nz,node)*(sfc_bulk_pz   + Z_3d_n(nz,node)*sfc_bulk_pz2)
                    rho_sfc = rho_sfc*sfc_rhopot/(rho_sfc+0.1_WP*Z_3d_n(nz,node))-density_0
                    kpp_dbsurf(nz) = -g * density_0_r *(rho_sfc-rho_nz)
                end do 
                
            ! --> MOM6 style --> see MOM_CVMix_KPP.F90 
            else
                do nz=nun, nln
                    !___________________________________________________________
                    ! Calculate the surface layer depth, averaged surface layer 
                    ! quantities
                    sldepth = kpp_surf_layer_ext*max( max(-Z_3d_n(nz,node),-zbar_3d_n(nun+1,node)),kpp_minOBLdepth )
                    nzsfc = nz
                    do nztmp = nun, nz
                        if (-zbar_3d_n(nztmp+1,node)>=sldepth) then
                            nzsfc = nztmp
                            exit
                        end if                        
                    end do
                    
                    if (kpp_use_StokesMOST) then
                        call compute_stokes_velocities_MOM6style( &
                            zbar_3d_n(nz,node),          &
                            zbar_3d_n(nz+1,node),        &
                            kpp_A_stokes,                &
                            u_wind(node), v_wind(node) , &
                            kpp_uS_t(nz), kpp_vS_t(nz) , &
                            kpp_uS_c(nz), kpp_vS_c(nz) , &
                            kpp_uS_m(nz), kpp_vS_m(nz))    
                                                    
                        call compute_stokes_velocities_MOM6style(& 
                            zbar_3d_n(nzsfc,node),       &
                            sldepth,                     &
                            kpp_A_stokes,                &
                            u_wind(node), v_wind(node),  &
                            uS_sld_t, vS_sld_t , &
                            uS_sld_c, vS_sld_c , &
                            uS_sld_m, vS_sld_m)                                   
                        
                        
                        call cvmix_kpp_compute_StokesXi (&
                            zbar_3d_n(nun:nln+1,node),   & ! full depth levels
                            Z_3d_n(nun:nln,node),        & ! mid depth levels
                            nzsfc,                       & ! cell index of Surface Layer Depth
                            sldepth,                     & ! surface layer depth > 0
                            kpp_buoyflx_nl(nz),          & ! surfce buoyancy flux (m2/s3) total
                            ustar,                       & ! turbulent friction velocity at surface (m/s), 
                            mesh%coriolis_node(node),    & ! Coriolis parameter (1/s) dim=1 
                            UVnode(1,nun:nln,node),      & ! zonal Eulerian mean horizontal velocity components
                            UVnode(2,nun:nln,node),      & ! merid Eulerian mean horizontal velocity components
                            !___stokes velocities___________________________________________
                            kpp_uS_t(nun:nln),           & ! zonal Surface Stokes drift velocity (Wave-induced drift)
                            kpp_VS_t(nun:nln),           & ! merid Surface Stokes drift velocity (Wave-induced drift)
                            kpp_uS_m(nun:nln),           & ! 
                            kpp_VS_m(nun:nln),           & ! 
                            uS_sld_t,                    & ! 
                            vS_sld_t,                    & ! 
                            uS_sld_m,                    & ! 
                            vS_sld_m,                    & ! 
                            !___outputs_____________________________________________________
                            kpp_stokesXi_z(nz)           & ! (out) Stokes Similartiy parameter 
                            )
                            
                        kpp_stokesVt_z(nz) = 0.0! kpp_stokesXi_z(nz)
                        
                        ! kpp_stokesXi_z: It represents the strength of Langmuir turbulence — 
                        ! how much the interaction between wind-driven shear and Stokes 
                        ! drift enhances vertical mixing.
                        
                        ! kpp_stokesVt_z: It represents the vortex force or momentum 
                        ! tendency term arising from the interaction between the Stokes 
                        ! drift and planetary rotation (the Coriolis effect).
                        ! carries the vertical structure of the Stokes drift contribution 
                        ! to momentum via the vortex force
                        ! --> switched off because in the KPP-only setup they’re not yet 
                        !     applying the explicit Stokes–Coriolis vortex term (which requires 
                        !     a coupled momentum solver and consistent Stokes profile).
                        
                        ! average quantities over surface layers
                        sfc_temp = 0.0_WP
                        sfc_salt = 0.0_WP
                        sfc_u    = 0.0_WP
                        sfc_v    = 0.0_WP
                        htot     = 0.0_WP
                        do nztmp = nun, nzsfc
                            delh     = min( max(0.0_WP,sldepth-htot), hnode(nztmp,node) )
                            htot     = htot+delh
                            sfc_temp = sfc_temp + temp(nztmp,node)*delh
                            sfc_salt = sfc_salt + salt(nztmp,node)*delh
                            sfc_u    = sfc_u    + (UVnode(1,nztmp,node)+kpp_uS_m(nztmp)) *delh
                            sfc_v    = sfc_v    + (UVnode(2,nztmp,node)+kpp_vS_m(nztmp)) *delh
                        end do
                        sfc_temp = sfc_temp/htot
                        sfc_salt = sfc_salt/htot
                        sfc_u    = sfc_u/htot
                        sfc_v    = sfc_v/htot
                        
                        !___________________________________________________________
                        ! calculate vertical shear between present layer and surface
                        ! averaged sfc_u and sfc_v
                        kpp_dvsurf2(nz) = (UVnode(1,nz,node)+kpp_uS_c(nz)-sfc_u)**2 + &
                                          (UVnode(2,nz,node)+kpp_vS_c(nz)-sfc_v)**2
                    
                    else ! kpp_use_StokesMOST ==.false.
                    
                        kpp_stokesXi_z = 0.0_WP
                        kpp_stokesVt_z = 0.0_WP
                        
                        ! average quantities over surface layers
                        sfc_temp = 0.0_WP
                        sfc_salt = 0.0_WP
                        sfc_u    = 0.0_WP
                        sfc_v    = 0.0_WP
                        htot     = 0.0_WP
                        do nztmp = nun, nzsfc
                            delh     = min( max(0.0_WP,sldepth-htot), hnode(nztmp,node) )
                            htot     = htot+delh
                            sfc_temp = sfc_temp + temp(nztmp,node)*delh
                            sfc_salt = sfc_salt + salt(nztmp,node)*delh
                            sfc_u    = sfc_u    + UVnode(1,nztmp,node) *delh
                            sfc_v    = sfc_v    + UVnode(2,nztmp,node) *delh
                        end do
                        sfc_temp = sfc_temp/htot
                        sfc_salt = sfc_salt/htot
                        sfc_u    = sfc_u/htot
                        sfc_v    = sfc_v/htot
                        
                        !___________________________________________________________
                        ! calculate vertical shear between present layer and surface
                        ! averaged sfc_u and sfc_v
                        kpp_dvsurf2(nz) = (UVnode(1,nz,node)-sfc_u)**2 + &
                                          (UVnode(2,nz,node)-sfc_v)**2
                        
                    end if ! --> if (kpp_use_StokesMOST) then
                    
                    !___________________________________________________________
                    ! calculate buoyancy difference between the surface averaged 
                    ! and the grid points blow 
                    ! --> bring density of surface point adiabatically to the same 
                    !     depth level as the deep point --> than calculate bouyancy 
                    !     difference
                    call densityJM_components(sfc_temp, sfc_salt, sfc_bulk_0, sfc_bulk_pz, sfc_bulk_pz2, sfc_rhopot, mesh)
                    call densityJM_components(temp(nz,node), salt(nz,node), bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)                    
                    rho_nz  = bulk_0   + Z_3d_n(nz,node)*(bulk_pz   + Z_3d_n(nz,node)*bulk_pz2)
                    rho_nz  = rho_nz*rhopot/(rho_nz+0.1_WP*Z_3d_n(nz,node))-density_0
                    rho_sfc = sfc_bulk_0   + Z_3d_n(nz,node)*(sfc_bulk_pz   + Z_3d_n(nz,node)*sfc_bulk_pz2)
                    rho_sfc = rho_sfc*sfc_rhopot/(rho_sfc+0.1_WP*Z_3d_n(nz,node))-density_0
                    kpp_dbsurf(nz) = -g * density_0_r *(rho_sfc-rho_nz)
                end do ! --> do nz=1, nln   
                
                !!PS do nz=2, nln 
                do nz=nun+1, nln 
                    !___________________________________________________________
                    ! calculate shear Richardson number Ri = N^2/(du/dz)^2 for 
                    ! mixing parameterisation below ocean boundary layer 
                    dz2     = (Z_3d_n( nz-1,node)-Z_3d_n( nz,node))**2
                    vshear2 = (UVnode(1,nz-1,node)-UVnode(1,nz,node))**2 + &
                              (UVnode(2,nz-1,node)-UVnode(2,nz,node))**2 
                    vshear2 = vshear2/dz2
                    kpp_shearRi(nz) = max(bvfreq(nz,node),0.0_WP)/(vshear2+kpp_epsln)
                end do ! --> do nz=1, nln
            end if
            
            
            
            !___________________________________________________________________
            !   ||    ||    ||    ||    ||    ||    ||    ||    ||    ||    ||  
            !  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_ 
            !  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  / 
            !   \/    \/    \/    \/    \/    \/    \/    \/    \/    \/    \/  
            !___________________________________________________________________
            ! 2) Calculate Richardson-dependent internal diffusivities, viscosities 
            ! (either 'PP' or 'KPP')
            ! --> After communication with CVMIX people, interior diffusivites 
            !     should be computed prior to the computation of KPP diffusivities
            !___________________________________________________________________
            ! --> PP parameterisation after Pacanowski and Philander 1981
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc internal mixing'//achar(27)//'[0m'
            if (kpp_internalmix .eq. 'PP') then
                do nz = nun+1, nln
                    kpp_Av(nz,node) = kpp_pp_Av0     /(1.0_WP+kpp_pp_alpha*kpp_shearRi(nz))**kpp_pp_loc_exp + kpp_Avbckg
                    kpp_Kv(nz,node) = kpp_Av(nz,node)/(1.0_WP+kpp_pp_alpha*kpp_shearRi(nz)) ! + background diffusivity
                end do
            
            !___________________________________________________________________
            ! --> KPP shear parameterization below the mixed layer after (Large et
            ! al.,1994)     
            elseif (kpp_internalmix .eq. 'KPP') then
                do nz = nun+1, nln
                    ! to avoid if conditions: if Ri>0, if 0<Ri<Ri0, if Ri>Ri0
                    ! see Large et al., 1994 eq. 28a, 28b, 28c 
                    aux = max(kpp_shearRi(nz),0.0_WP)
                    aux = min(aux/kpp_Ri0,1.0_WP)
                    aux = (1.0_WP-aux**2.0)
                    aux = aux**kpp_loc_exp
                    kpp_Av(nz,node) = kpp_Av0 * aux + kpp_Avbckg
                    kpp_Kv(nz,node) = kpp_Kv0 * aux ! + background diffusivity comes in next step
                end do 
                
            else
                write(*,*) " --> Error: this kpp_internalmix scheme is not supported"
                write(*,*) "     for the mixing below the OBL, either KPP or PP !"
                call par_ex(partit%MPI_COMM_FESOM, partit%mype)
            end if 
            
            
            !___________________________________________________________________
            ! 3) Set Background diffusivities either to const. background 
            ! diffusivity (kpp_Kvbckg) or non-const. background diffusivity 
            ! of Qiang from FESOM1.4
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc background diff'//achar(27)//'[0m'
            if (kpp_use_nonconstKvb) then
                do nz = nun+1, nln
                    call Kv0_background_qiang( &
                        aux,geo_coord_nod2D(2,node)/rad,abs(zbar_3d_n(nz,node))&
                        )
                    kpp_Kv(nz,node) = kpp_Kv(nz,node) + aux
                end do 
            else
                do nz = nun+1, nln
                    kpp_Kv(nz,node) = kpp_Kv(nz,node) + kpp_Kvbckg
                end do 
            end if
            
            
            !___________________________________________________________________
            ! 4) Compute turbulent scales 
            ! --> make and educated guess for w_s since we need to compute bulkRi
            ! --> w_s is re-computed in call cvmix_coeffs_kpp_low under consideration 
            !     StokesXi when CVmix_kpp_params_in%lStokesMOST = .True.
            ! --> sigma_coord must be here 1.0!!! --> see MOM6 --> parameterizations/vertical/MOM_CVMix_KPP.F90:1343
            call cvmix_kpp_compute_turbulent_scales(           &
                sigma_coord     = 1.0_WP            ,          & ! (in)  sigma: Normalized surface layer depth
                OBL_depth       = abs(Z_3d_n(nun:nln,node)),   & ! (in)  Assume OBL depth (m) =  mid-depth level
                surf_buoy_force = kpp_buoyflx_nl(nun:nln),     & ! (in)  surfce buoyancy flux (m2/s3) consider sw_pene
                surf_fric_vel   = ustar,                       & ! (in)  turbulent friction velocity at surface (m/s)
                xi              = kpp_stokesVt_z(nun:nln),    & ! (in)  Stokes parameter xi= Ps/(PU+PS+PB)
                w_s             = kpp_ws_cntr(nun:nln)         & ! (out) Turbulent velocity scale profile (m/s) for skalars
                ) 
            ! --> need w_s to compute cvmix_kpp_compute_bulk_Richardson(...)
            
            
            
            !___________________________________________________________________
            ! 5) Make a first guess about the OBLdepth and compute Langmuir 
            ! number and Langmuir enhancement factor 
            !ATTENTION! here make a hirst educated guess, and use the OBL from the previous 
            ! time step. the OBL of the actual timestep can only be computed after bulkRI 
            ! is known
            oblguess = max(kpp_obldepth(node), kpp_minOBLdepth)
            
            ! compute norm of wind velocity 
            wind_norm = sqrt(u_wind(node)**2 + v_wind(node)**2)
            
            ! computes the surface layer averaged Stokes drift, given
            ! the 10-meter wind (m/s) and the boundary layer depth (m).
            uv_SLmean = cvmix_kpp_ustokes_SL_model(wind_norm, oblguess, CVmix_params_in)
            
            ! Copute Langmuir enhance,ent factor & SL langmuir number
            call compute_Efactor(               &
                kpp_langmuir_entrainment,       &
                ustar,                          & ! friction velocity
                wind_norm,                      & ! surface wind velocity 
                uv_SLmean,                      & ! 
                kpp_LaSL(node),                 & ! Surface Layer Langmuir Number
                kpp_EFactor(node)               & ! Langmuir Entrainment Enhancement Factor
                )             
            
            
            
            !___________________________________________________________________
            ! 6) Calculate Bulk Richardson number at centers of cell 
            ! --> assuming OBLdepth is at mid-depth level interface Z_3d_n
            ! Hereafter kpp_bulkRi(k) is known for each column, then CVMix interpolates 
            ! to find the actual OBLdepth. This approach avoids need to iterate
            ! on the OBLdepth calculation. It follows the approach that used in MOM6
            ! |
            ! +-> extent of oceanic OBL depth depends on surrface forcing, oceanic 
            ! |   bouyancy and velocity profile 
            ! +-> define Bulk Richardson number relative to surface
            !     Rib = ( buoy_surf-buoy(z) )*z/( |v_surf-v(z)|^2 + v_t^2(z) )
            !     --> v_t = turbulent velocity shear, needs to be parameterized
            !     --> v_t = z * ws * N * v_tc
            !                             |-> v_tc = Cv * sqrt(0.2/Cs/epsilon) / kpp_vonKarman^2 / kpp_Rib_crit
            kpp_bulkRi(nun:nln+1) = cvmix_kpp_compute_bulk_Richardson(   &
                zt_cntr         = Z_3d_n(     nun:nln  ,node),  & ! (in) Depth of cell center (m)
                delta_buoy_cntr = kpp_dbsurf( nun:nln),         & ! (in) Bulk buoyancy difference, Br-B(z) (1/s)
                delta_Vsqr_cntr = kpp_dvsurf2(nun:nln),         & ! (in) Square of resolved velocity difference (m2/s2)
                ws_cntr         = kpp_ws_cntr(nun:nln),         & ! (in) Turbulent velocity scale profile (m/s)
                Nsqr_iface      = bvfreq(     nun:nln+1,node),  & ! (in) Buoyancy frequency (1/s)
                bfsfc           = kpp_buoyflx_nl(nun:nln),      & ! (in) surface buoyancy flux (units: m^2/s^3)
                ustar           = ustar,                        & ! (in) friction velocity (units: m/s)
                EFactor         = kpp_EFactor(node),            & ! (in) Langmuir enhancement factor for entrainment (units: none)
                LaSL            = kpp_LaSL(node)                & ! (in) surface layer averaged Langmuir number (units: none)
                )    
            
            
            !___________________________________________________________________
            ! 7) Compute depth of oceanic boundary layer (kpp_obldepth)
            ! |
            ! +-> kpp_obldepth h is smallest value of z at which this bulk Richardson 
            ! |    number Rib(z) is equal a critical value Ric (kpp_Rib_crit)
            ! +-> if kpp_use_compEkman=.True. ... kpp_obldepth can not be deeper than 
            ! !   ekmann layer depth (h_ekmann = 0.7 * ustar/coriolis)
            ! +-> if kpp_use_monob=.True. ... kpp_obldepth can not be deeper than monin
            !     obukov length (h_monob = ustar^3/( kpp_sbuoyflx*kpp_vonKarman )
            !
            ! Attention !!!: in case of kpp_use_compEkman=.True. cvmix DOES NOT 
            ! check if surface buoyancy forcing is stable, it limits to ekman depth
            ! everywhere and eradicates possible obldepths from bulkRi number 
            ! --> thats why kpp_use_compEkman = .False.
            ! --> check if its still the case
            call cvmix_kpp_compute_OBL_depth(            &
                Ri_bulk    = kpp_bulkRi(nun:nln+1),      & ! (in) Bulk Richardson number dim=(ke+1)
                zw_iface   = zbar_3d_n( nun:nln+1,node), & ! (in) Height of interfaces (m) dim=(ke+1)
                OBL_depth  = kpp_obldepth(node),         & ! (out) OBL depth (m) dim=1
                kOBL_depth = kpp_nzobldepth(node),       & ! (out) level (+fraction) of OBL extent dim=1
                zt_cntr    = Z_3d_n(    nun:nln  ,node), & ! (in) Height of cell centers (m) dim=(ke)
                surf_fric  = ustar,                      & ! (in) Turbulent friction velocity at surface (m/s) dim=1
                surf_buoy  = aux_surfbuoyflx_nl(nun:nln),& ! (in) Buoyancy flux at surface (m2/s3) dim=1
                Coriolis   = mesh%coriolis_node(node),   & ! (in) Coriolis parameter (1/s) dim=1
                Xi         = kpp_stokesXi_z(nun:nln)     &
                )    
            kpp_nzobldepth(node) = kpp_nzobldepth(node) + nun - 1
            
            !___safty switches for kpp_nzobldepth_______________________________
            ! A hack from MOM6 to avoid KPP reaching the bottom. It was needed during
            ! development because KPP was unable to handle vanishingly small layers
            ! near the bottom
            if (kpp_deepOBLoffset > 0.0) then
                ! aux ... zBottomMinusOffset
                aux = zbar_3d_n(nln+1,node) + min(kpp_deepOBLoffset,-0.1_WP*zbar_3d_n(nln+1,node))
                kpp_obldepth(node) = min( kpp_obldepth(node), -aux )
            end if 
            
            ! no shallower than top layer
            !!PS kpp_obldepth(node)  = max( kpp_obldepth(node), abs(zbar_3d_n(2, node)) )
            kpp_obldepth(node)  = max( kpp_obldepth(node), abs(zbar_3d_n(nun+1, node)) )
            
            ! no deeper than bottom layer 
            !!PS kpp_obldepth(node)  = min( kpp_obldepth(node), abs(zbar_3d_n(nlevels_nod2D(node),node)) )
            kpp_obldepth(node)  = min( kpp_obldepth(node), abs(zbar_3d_n(nln+1,node)) )
            
            ! safety for kOBL
            if (kpp_nzobldepth(node) > nln+1) then
                kpp_nzobldepth(node) = nln+1
            end if 
            
            
            !___________________________________________________________________
            ! 8) Now were OBLdepth is known, recompute SLdepth (surface layer depth)
            ! and recompute StokeXi similarity parameter and Langmuir enhancment factor 
            if (kpp_use_StokesMOST) then 
                !___________________________________________________________________
                ! re compute surface layer depth (SLDEPTH) now that know the actual 
                ! OBLDepth
                ! In the Large et al. (1994) K-Profile Parameterization (KPP), the 
                ! surface layer depth! hs (often called the surface layer thickness) 
                ! is a single scalar quantity per column, defined as:
                !
                !                           h_s = ϵ*h_b
                !
                ! where:
                ! hb ... boundary layer depth (OBL depth)
                ! ϵ  ... a constant fraction, typically 0.1
                !
                ! It is used to compute the averaged surface layer quantities 
                ! (mean T, S, U, V, etc.) for the bulk Richardson number test.
                ! The surface layer in KPP is the thin top region where fluxes are 
                ! applied and where surface averages are taken — not something that 
                ! varies with depth.
                sldepth = kpp_surf_layer_ext*max(kpp_obldepth(node), kpp_minOBLdepth)
                nzsfc = kpp_obldepth(node)
                do nztmp = nun, int(kpp_obldepth(node))
                    if (-zbar_3d_n(nztmp+1,node)>=sldepth) then
                        nzsfc = nztmp
                        exit
                    end if                        
                end do
                
                call compute_stokes_velocities_MOM6style(& 
                    zbar_3d_n(nzsfc,node),       &
                    sldepth,                     &
                    kpp_A_stokes,                &
                    u_wind(node), v_wind(node),  &
                    uS_sld_t, vS_sld_t ,         &
                    uS_sld_c, vS_sld_c ,         &
                    uS_sld_m, vS_sld_m) 
                        
                call cvmix_kpp_compute_StokesXi (&
                    zbar_3d_n(nun:nln+1,node),   & ! full depth levels
                    Z_3d_n(nun:nln,node),        & ! mid depth levels
                    nzsfc,                       & ! cell index of Surface Layer Depth
                    sldepth,                     & ! surface layer depth > 0
                    kpp_buoyflx_nl(kpp_nzobldepth(node)),          & ! surfce buoyancy flux (m2/s3) total
                    ustar,                       & ! turbulent friction velocity at surface (m/s), 
                    mesh%coriolis_node(node),    & ! Coriolis parameter (1/s) dim=1 
                    UVnode(1,nun:nln,node),      & ! zonal Eulerian mean horizontal velocity components
                    UVnode(2,nun:nln,node),      & ! merid Eulerian mean horizontal velocity components
                    !___stokes velocities___________________________________________
                    kpp_uS_t(nun:nln),           & ! zonal Surface Stokes drift velocity (Wave-induced drift)
                    kpp_VS_t(nun:nln),           & ! merid Surface Stokes drift velocity (Wave-induced drift)
                    kpp_uS_m(nun:nln),           & ! 
                    kpp_VS_m(nun:nln),           & ! 
                    uS_sld_t,                    & ! 
                    vS_sld_t,                    & ! 
                    uS_sld_m,                    & ! 
                    vS_sld_m,                    & ! 
                    !___outputs_____________________________________________________
                    kpp_stokesXi(node) )           ! (out) Stokes Similartiy parameter 
            end if ! --> if (kpp_use_StokesMOST) then 
            
            ! computes the surface layer averaged Stokes drift, given
            ! the 10-meter wind (m/s) and the boundary layer depth (m).
            uv_SLmean = cvmix_kpp_ustokes_SL_model(wind_norm, kpp_obldepth(node), CVmix_params_in)
            
            ! Copute Langmuir enhance,ent factor & SL langmuir number
            call compute_Efactor(               &
                kpp_langmuir_entrainment,       &
                ustar,                          & ! friction velocity
                wind_norm,                      & ! surface wind velocity 
                uv_SLmean,                      & ! 
                kpp_LaSL(node),                 & ! Surface Layer Langmuir Number
                kpp_EFactor(node)               & ! Langmuir Entrainment Enhancement Factor
                )    
            
            !___________________________________________________________________
            ! 9) Call CVMix/KPP to obtain OBL diffusivities, viscosities and non-
            ! local transports
            !___________________________________________________________________
            ! the dealing of shortwave penetration comes from the original kpp
            ! parameterisation of FESOM1.4 & FESOM2.0
            ! --> interpolate contribution that comes from shortwave penetration
            ! to the depth of the obldepth
            aux_surfbuoyflx_nl(1) = kpp_sbuoyflx(node)
            if (use_sw_pene .and. kpp_use_fesomkpp) then
                aux_nz = int(kpp_nzobldepth(node))
                ! take only penetrated shortwave radiation heatflux into account 
                ! that reached until the obldepth --> do linear interpolation 
                aux_surfbuoyflx_nl(1) = aux_surfbuoyflx_nl(1) + &
                                        aux_coeff * &
                                        ( sw_3d(nun,node) - &
                                            ( sw_3d(aux_nz-1, node) + &
                                                ( sw_3d(aux_nz, node) - sw_3d(aux_nz-1, node) ) &
                                                * ( kpp_obldepth(node) + zbar_3d_n( aux_nz-1,node) ) &
                                                / ( zbar_3d_n(aux_nz-1,node) - zbar_3d_n(aux_nz,node) ) ) )
                                                
            ! MOM6 provides different option how buoyancy flux is influenced by 
            ! short wave penetration flux
            ! --> mxl comes closest to what FESOM1.4 was doing 
            elseif (use_sw_pene .and. (.not. kpp_use_fesomkpp)) then
                if     (trim(kpp_sw_method) == 'all')  then
                    aux_surfbuoyflx_nl(1) = aux_surfbuoyflx_nl(1)+aux_coeff*sw_3d(1,node) 
                elseif (trim(kpp_sw_method) == 'mxl')  then
                    aux_surfbuoyflx_nl(1) = aux_surfbuoyflx_nl(1)+aux_coeff*(sw_3d(1,node)-sw_3d(int(kpp_nzobldepth(node))+1,node))
                elseif (trim(kpp_sw_method) == 'lvl1') then
                    aux_surfbuoyflx_nl(1) = aux_surfbuoyflx_nl(1)+aux_coeff*(sw_3d(1,node)-sw_3d(2,node))
                ! --> 'fesom' not part of mom6    
                elseif (trim(kpp_sw_method) == 'fesom') then
                    aux_surfbuoyflx_nl(1) = aux_surfbuoyflx_nl(1) + &
                                        aux_coeff * &
                                        ( sw_3d(nun,node) - &
                                            ( sw_3d(aux_nz-1, node) + &
                                                ( sw_3d(aux_nz, node) - sw_3d(aux_nz-1, node) ) &
                                                * ( kpp_obldepth(node) + zbar_3d_n( aux_nz-1,node) ) &
                                                / ( zbar_3d_n(aux_nz-1,node) - zbar_3d_n(aux_nz,node) ) ) )
                end if
            end if 
            
            !___________________________________________________________________
            ! 10) write kpp_Av and kpp_Kv that contain only the ocean internal 
            ! mixing + background mixing into KdiffT, KdiffS, Kvisc --> are only 
            ! required for  kpp_matchtechc= 'MatchGradient' or 'MatchBoth'  as
            ! well as for the diffusive enhancement at the last OBL layer
            !!PS if ( .not. trim(kpp_matchtechc)=='MatchBoth') then 
                kpp_oblmixc(:,node,1) = kpp_Av(:,node)
                kpp_oblmixc(:,node,2) = kpp_Kv(:,node)
                kpp_oblmixc(:,node,3) = kpp_Kv(:,node)
            !!PS else
            !!PS     ! --> this part is done in MOM6 and MPIOM but i dont realy see 
            !!PS     ! the reason for doing this --> need to test this !!! -->
            !!PS     ! this leads to slightly lower and slightly less deeper
            !!PS     ! diffusivity values --> i would not use it !      
            !!PS     kpp_oblmixc(:,node,1) = 0.0_WP
            !!PS     kpp_oblmixc(:,node,2) = 0.0_WP
            !!PS     kpp_oblmixc(:,node,3) = 0.0_WP
            !!PS end if 
            
            !___________________________________________________________________
            ! 11) compute the turbulent diffusion coefficients
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc kpp coeff'//achar(27)//'[0m'
            call cvmix_coeffs_kpp(                     &
                Mdiff_out = kpp_oblmixc(:,node,1),     & ! (inout) new_Mdiff: Total viscosity (m2/s)
                Tdiff_out = kpp_oblmixc(:,node,2),     & ! (inout) new_Tdiff: Total heat diffusivity (m2/s)
                Sdiff_out = kpp_oblmixc(:,node,3),     & ! (inout) new_Sdiff: Total salt diffusivity (m2/s)
                zw        = zbar_3d_n(:,node),         & ! (in)    zw_iface: Height of interfaces (m)
                zt        = Z_3d_n(:,node),            & ! (in)    zt_cntr: Height of level centers (m)
                old_Mdiff = kpp_oblmixc(:,node,1),     & ! (in)    MDiff_iface: Original viscosity (m2/s)
                old_Tdiff = kpp_oblmixc(:,node,2),     & ! (in)    Tdiff_iface: Original heat diffusivity (m2/s)
                old_Sdiff = kpp_oblmixc(:,node,3),     & ! (in)    Sdiff_iface: Original salt diffusivity (m2/s)
                OBL_depth = kpp_obldepth(node),        & ! (in)    BoundaryLayerDepth: OBL depth (m)
                kOBL_depth= kpp_nzobldepth(node),      & ! (in)    kOBL_depth: level (+fraction) of OBL extent
                Tnonlocal = kpp_nonlcltranspT(:,node), & ! (out)   kpp_Tnonlocal_iface: Non-local heat transport (non-dimensional)
                Snonlocal = kpp_nonlcltranspS(:,node), & ! (out)   kkp_Snonlocal_iface: Non-local salt transport (non-dimensional)
                surf_fric = ustar,                 & ! (in)    SurfaceFriction:Turbulent friction velocity at surface (m/s)
                surf_buoy = aux_surfbuoyflx_nl(1),     & ! (in)    SurfaceBuoynacyForcing: Buoyancy flux at surface (m2/s3)
                nlev      = nlevels_nod2D(node)-1,     & ! (in)    nlev: Number of levels to compute coeffs for
                max_nlev  = nl-1,                      & ! (in)    max_lev: maximum vertical levels
                Langmuir_EFactor = kpp_EFactor(node),  & ! (in)    Langmuir enhancement factor
                StokesXI         = kpp_stokesXI(node)  & !
                )
                
            
            !___________________________________________________________________
            ! 12) computation of nonlocal transport terms done in cvmix_coeffs_kpp
            ! --> the diffusive enhancement of the last OBL layer is applied at 
            !     the end of call cvmix_coeffs_kpp (kpp_use_enhanceKv=.true.) 
            !     --> After Large et al. it is used to remove bias and dampen 
            !     oscillation, see Appendix D and eq. D6 
            
            !___________________________________________________________________
            ! MOM6: Over-write CVMix NLT shape function with one of the 
            !       following choices. The CVMix code has yet to update for 
            !       these options, so we compute it here. Note that 
            !       nonLocalTrans = Cs * G(sigma) (LMD94 notation), with Cs =
            !       6.32739901508.
            !       Start do-loop at k=2, since k=1 is ocean surface (sigma=0)
            !       and we do not wish to double-count the surface forcing.
            !       Only compute nonlocal transport for 0 <= sigma <= 1.
            !       MOM6 recommended shape is the parabolic; it gives deeper 
            !       boundary layer and no spurious extrema.
            ! --> I dont think this part is neccessary either since the treatment 
            !     of parabolic non-local transport seems to be already included 
            !     in the cvmix used heren, kept only for comparability reasons 
            !     between FESOM2 and MOM6/MPIOM --> keep kpp_nlt_shape='cvmix'
            if (aux_surfbuoyflx_nl(1) < 0.0_WP) then
                if (kpp_nlt_shape == "cubic") then
                    !!PS do nz = 2,nln+1
                    do nz = nun+1,nln+1
                        sigma = min(1.0_wp,-zbar_3d_n(nz,node)/kpp_obldepth(node))
                        kpp_nonlcltranspT(nz,node) = (1.0_wp - sigma)**2 * (1.0_wp + 2.0_wp*sigma)!*cs2
                        kpp_nonlcltranspS(nz,node) = kpp_nonlcltranspT(nz,node)
                    end do
                elseif (kpp_nlt_shape == "parabolic") then
                    !!PS do nz = 2,nln+1
                    do nz = nun+1,nln+1
                        sigma = min(1.0_wp,-zbar_3d_n(nz,node)/kpp_obldepth(node))
                        kpp_nonlcltranspT(nz,node) = (1.0_wp - sigma)**2 !*cs2
                        kpp_nonlcltranspS(nz,node) = kpp_nonlcltranspT(nz,node)
                    end do
                elseif (kpp_nlt_shape == "linear") then
                    !!PS do nz = 2,nln+1
                    do nz = nun+1,nln+1
                        sigma = min(1.0_wp,-zbar_3d_n(nz,node)/kpp_obldepth(node))
                        kpp_nonlcltranspT(nz,node) = (1.0_wp - sigma)!*cs2
                        kpp_nonlcltranspS(nz,node) = kpp_nonlcltranspT(nz,node)
                    end do
                elseif (kpp_nlt_shape == "cubic_LMD") then
                    !!PS do nz = 2,nln+1
                    do nz = nun+1,nln+1
                        sigma = min(1.0_wp,-zbar_3d_n(nz,node)/kpp_obldepth(node))
                        kpp_nonlcltranspT(nz,node) = kpp_cs2 * sigma*(1.0_wp -sigma)**2
                        kpp_nonlcltranspS(nz,node) = kpp_nonlcltranspT(nz,node)
                    end do 
                end if     
            end if 
        end do !--> do node= 1,myDim_nod2D   
        
        
        !_______________________________________________________________________
        ! 13) horizontal smoothing of the OBL mixing coefficient --> approach from
        ! original kpp parameterisation of FESOM1.4 & FESOM2.0
        !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc smooth kpp_oblmixc'//achar(27)//'[0m'    
        if (kpp_use_smoothblmc .and. kpp_use_fesomkpp) then
            call exchange_nod(kpp_oblmixc(:,:,1), partit)
            call exchange_nod(kpp_oblmixc(:,:,2), partit)
            call exchange_nod(kpp_oblmixc(:,:,3), partit)
            do nz=1, 3
                !_______________________________________________________________
                ! all loops go over myDim_nod2D so no halo information --> for 
                ! smoothing haloinfo is required --> therefor exchange_nod
                call smooth_nod(kpp_oblmixc(:,:,nz), kpp_smoothblmc_nmb, partit, mesh)
            end do
        end if 
        
        
        !_______________________________________________________________________
        !   ||    ||    ||    ||    ||    ||    ||    ||    ||    ||    ||  
        !  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_ 
        !  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  / 
        !   \/    \/    \/    \/    \/    \/    \/    \/    \/    \/    \/  
        !_______________________________________________________________________
        ! 14) combine KPP OBL mixing coefficents with internal mixing coefficients
        ! --> ...min(kpp_oblmixc(nz,node,2),kpp_oblmixc(nz,node,3))... comes from MOM6
        !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc combine OBL & internal mixing'//achar(27)//'[0m'    
        do node= 1, myDim_nod2D
            ! FESOM14 approach
            if (kpp_use_fesomkpp) then
                !!PS do nz = 2,int(kpp_nzobldepth(node))
                do nz = ulevels_nod2D(node)+1,int(kpp_nzobldepth(node))
                    kpp_Av(nz,node) = max(kpp_Av(nz,node),kpp_oblmixc(nz,node,1))
                    kpp_Kv(nz,node) = max(kpp_Kv(nz,node),kpp_oblmixc(nz,node,2))
                end do
            ! MOM6/MPIOM approach
            else
                !!PS do nz = 2,nlevels_nod2D(node)-1
                do nz = ulevels_nod2D(node)+1,nlevels_nod2D(node)-1
                    kpp_Av(nz,node) = kpp_oblmixc(nz,node,1)
                    kpp_Kv(nz,node) = min(kpp_oblmixc(nz,node,2),kpp_oblmixc(nz,node,3))
                end do
            end if 
        end do
        
        
        !_______________________________________________________________________
        ! write out diffusivities to FESOM2.0 --> diffusivities remain on nodes
        call exchange_nod(kpp_Kv, partit)
        Kv = kpp_Kv
        
        
        !_______________________________________________________________________
        ! write out viscosities to FESOM2.0 --> viscosities for FESOM2.0 are 
        ! defined on elements --> interpolate therefor from nodes to elements
        call exchange_nod(kpp_Av, partit)
        Av = 0.0_WP
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            !!PS do nz=2,nlevels(elem)-1
            do nz=ulevels(elem),nlevels(elem)-1
                Av(nz,elem) = sum(kpp_Av(nz,elnodes))/3.0_WP    ! (elementwise)                
            end do
        end do
        
        
    end subroutine calc_cvmix_kpp    
    
    
    
    !
    !___________________________________________________________________________
    ! Since we dont have a fully coupled wave model yet, approximate computation 
    ! of Stokes drift velocities based on surface wind.
    ! This routine provides a simple parameterization of the Stokes drift velocity
    ! profile following the approach of Breivik et al. (2016, *J. Phys. Oceanogr.*)
    ! using a Phillips spectrum approximation. It is intended as a placeholder
    ! until a fully coupled wave model is available.
    ! The Stokes drift velocity describes the net Lagrangian drift induced by
    ! surface gravity waves. In deep water, the vertical structure of the Stokes
    ! drift decays exponentially with depth:
    !
    !    uS(z) = uS0 * exp(-2*k_p *|z - z_sfc}| )
    !
    ! where:
    !   - uS0 is the surface Stokes drift velocity,
    !   - k_p is the wavenumber corresponding to the spectral peak,
    !   - z is depth of layer,  .
    !
    ! The peak wavenumber and frequency are estimated from the wind speed using
    ! a fetch-limited empirical relationship:
    !
    !     T_p = alpha * UV10, f_p = 1 / T_p, 
    !     k_p = (2 * pi * f_p)^2 / g
    !
    ! with a typical choice alpha = 0.8.
    !
    ! The magnitude of the surface Stokes drift is approximated as:
    !
    !     |uS0| = C * UV10^2
    !
    ! where C ~ 0.01 (tunable between 0.005–0.02), following Ardhuin et al. (2009).
    subroutine compute_stokes_velocities_MOM6style(ztop   , zbot   , &
                                                   Astokes,          &
                                                   uwind  , vwind  , &
                                                   uS_top , vS_top , &
                                                   uS_mid , vS_mid , &     
                                                   uS_mean, vS_mean)
        implicit none
        real(kind=WP), intent(in)  :: ztop, zbot
        real(kind=WP), intent(in)  :: Astokes, uwind, vwind
        real(kind=WP), intent(out) :: uS_top, vS_top, uS_mid, vS_mid, uS_mean, vS_mean
        
        real(WP) :: uv10, Tp, fp, kp, uS0_mag, uS0, vS0, C, fexp
        
        uv10 = sqrt(uwind**2 + vwind**2) 
        
        !_______________________________________________________________________
        ! Empirical relationships
        ! Peak period Tp ≈ α * U10  (fetch-limited wind-sea)
        ! A common choice: Tp ≈ 0.8 * U10  (Tp in seconds, U10 in m/s)
        Tp = 0.8 * uv10
        fp = 1.0 / Tp
        kp = (2.0 * pi * fp)**2 / g     ! deep-water dispersion relation

        ! Empirical coefficient linking U10 to surface Stokes drift (Ardhuin 2009)
        ! uS0 ≈ C * U10^2  (typical C = 0.01, can tune between 0.005–0.02)
        !!PS Astokes = 0.01
        uS0_mag = Astokes * uv10**2
        uS0 = uS0_mag * uwind/uv10
        vS0 = uS0_mag * vwind/uv10
  
        ! Compute depth-dependent values
        ! Interface (bottom of cell)
        fexp = exp(-2.0 * kp * abs(zbot))
        uS_top = uS0 * fexp
        vS_top = vS0 * fexp

        ! Midpoint
        fexp = -exp(-2.0*kp * abs((ztop + zbot)*0.5_WP))
        uS_mid = uS0 * fexp
        vS_mid = vS0 * fexp

        ! Layer mean (analytic integration)
        fexp = exp(-2.0 * kp * abs(ztop)) - exp(-2.0 * kp * abs(zbot))
        uS_mean = 0.5 * uS0 * fexp / (kp * abs(ztop - zbot))
        vS_mean = 0.5 * vS0 * fexp / (kp * abs(ztop - zbot))
        
    end subroutine compute_stokes_velocities_MOM6style
    
    
    !
    !
    !___________________________________________________________________________
    ! Langmuir entrainment enhancement factor (EFactor)
    !
    ! Derivation / justification (compact):
    !   - Use the surface-layer Langmuir number LaSL = sqrt(u_star / uS_SL).
    !   - Define s = LaSL^{-2} = uS_SL / u_star (this is CVMix's `lasl_sqr_i`).
    !   - Empirically fit LES results (Li et al. 2016/2017 and followups) with a
    !     simple quadratic correction to E^2:
    !        E^2 = 1 + a*s + b*s^2
    !     so that E -> 1 when s -> 0 and the curve can match LES at moderate/large s.
    !   - CVMix / Li-style coefficients (from LES fits) use
    !        a = 1 / (1.5)^2
    !        b = 1 / (5.4)^4
    !     giving:
    !        E = sqrt( 1 + (1/1.5^2) * s + (1/5.4^4) * s^2 )
    !
    !   - Physically: E multiplies the shear-related turbulent velocity scale (not
    !     the convective component) to represent Langmuir-enhanced shear mixing.
    !
    ! References:
    !   Li & Fox-Kemper (2017), J. Phys. Oceanogr. — LES-based study of Langmuir
    !     effects and KPP modifications (see section on enhancement factor and
    !     fitting).  (primary LES reference used in CVMix).    
    !   - CVMix implementation (the `cvmix_kpp_EFactor_model` / `cvmix_kpp_ustokes_SL_model` routines) — the code you inspected uses the same `lasl_sqr_i` → `EFactor` mapping. :contentReference[oaicite:5]{index=5}
    subroutine compute_Efactor( method,                 &
                                ustar,                  & ! friction velocity
                                wind_norm,              & ! norm surface wind velocity 
                                SL_norm,                & ! norm stokes drift averaged over SL
                                LaSL,                   & ! Surface Layer Langmuir Number
                                EFactor)                  ! Langmuir Entrainment Enhancement Factor
        implicit none 
        character    , intent(in)  :: method
        real(kind=WP), intent(in)  :: ustar
        real(kind=WP), intent(in)  :: wind_norm
        real(kind=WP), intent(in)  :: SL_norm
        real(kind=WP), intent(out) :: LaSL
        real(kind=WP), intent(out) :: EFactor
        
        real(kind=WP)              :: SL_eps=1.0e-8
        
        if (wind_norm .gt. 0.0_WP .and. ustar .gt. 0.0_WP) then 
        
            ! Compute the Surface Layer Langmuir Number
            if (SL_norm < SL_eps) then 
                LaSL    = sqrt(ustar / SL_eps)
            else    
                LaSL    = sqrt(ustar / SL_norm)
            end if
                
            ! Compute the Langmuir Entrainment Enhancement Factor (EFactor)
            ! Once LaSL is known, EFactor represents the enhancement of entrainment 
            ! and mixing due to Langmuir turbulence.
            if (method=='LF16') then 
                ! According to Li et al. (2016) --> Mc Williams and Sullivan (2000) proposed 
                EFactor =  sqrt(1.0_wp + 0.08 / LaSL**4)
            
            else ! (kpp_langmuir_entrainment=='LF17') then 
                ! According to Li et al. (2017) --> Harcourt and D’Asaro (2008) and
                ! Van Roekel et al. (2012) proposed
                EFactor =  sqrt(1.0_wp + (1.0_wp / 1.5_wp**2) / LaSL**2 + &
                                         (1.0_wp / 5.4_wp**4) / LaSL**4)
            end if 
            
        else
        
            EFactor = 1.0_WP
            LaSL    = 0.0_WP
            
        end if    
        
        !_______________________________________________________________________
        return
            
    end subroutine compute_Efactor          
        
end module g_cvmix_kpp
