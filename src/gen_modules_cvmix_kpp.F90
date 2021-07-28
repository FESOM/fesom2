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
    USE cvmix_kpp, only : cvmix_init_kpp, cvmix_put_kpp, CVmix_get_kpp_real, &
                          cvmix_coeffs_kpp, cvmix_kpp_compute_OBL_depth,     &
                          cvmix_kpp_compute_turbulent_scales,                &   
                          cvmix_kpp_compute_bulk_Richardson,                 & 
                          cvmix_kpp_compute_unresolved_shear,                & 
                          cvmix_kpp_params_type,                             & 
                          cvmix_kpp_compute_kOBL_depth                    
    
    !___________________________________________________________________________
    ! module calls from FESOM
    use g_config
    use o_param           
    use mod_mesh
    use g_parsup
    use o_arrays
    use g_comm_auto 
    use i_arrays
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
    ! cvmix_kpp.F90 there is a bug in the enhancement --> i fixed that in the file 
    ! cvmix_kpp_fix.F90)
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
                         kpp_use_LMDws, kpp_sw_method,kpp_nlt_shape
    
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
    real(kind=WP), allocatable, dimension(:)    :: kpp_obldepth,kpp_nzobldepth
    ! surface buoyancy flux
    real(kind=WP), allocatable, dimension(:)    :: kpp_sbuoyflx
    
    ! 3d arrays
    !-----------
    ! nodal viscosities/diffusivities
    real(kind=WP), allocatable, dimension(:,:)  :: kpp_Av, kpp_Kv
    ! non local transport for temp. and salt
    real(kind=WP), allocatable, dimension(:,:)  :: kpp_nonlcltranspT, kpp_nonlcltranspS
    ! obl mixing coefficient for momentum, temp/salt diffusivity
    real(kind=WP), allocatable, dimension(:,:,:):: kpp_oblmixc
    
    contains
    !
    !
    !
    !===========================================================================
    ! allocate and initialize CVMIX KPP variables --> call initialisation 
    ! routine from cvmix library
    subroutine init_cvmix_kpp(mesh)
        implicit none
        character(len=MAX_PATH) :: nmlfile
        logical            :: nmlfile_exist=.False.
        integer            :: node_size
        type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
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
        kpp_bulkRi        = 0.0_WP
        kpp_shearRi       = 0.0_WP
        kpp_dvsurf2       = 0.0_WP
        kpp_dbsurf        = 0.0_WP
        kpp_ws_cntr       = 0.0_WP
        
        ! allocate horizontal 2D variable 
        allocate(kpp_obldepth(node_size),kpp_nzobldepth(node_size))
        allocate(kpp_sbuoyflx(node_size))
        kpp_obldepth      = 0.0_WP
        kpp_nzobldepth    = 0.0_WP
        kpp_sbuoyflx      = 0.0_WP
        
        ! allocate 3D variable 
        allocate(kpp_Av(nl,node_size),kpp_Kv(nl,node_size))
        allocate(kpp_nonlcltranspT(nl,node_size),kpp_nonlcltranspS(nl,node_size))
        kpp_Av            = 0.0_WP
        kpp_Kv            = 0.0_WP
        kpp_nonlcltranspT = 0.0_WP
        kpp_nonlcltranspS = 0.0_WP
        
        allocate(kpp_oblmixc(nl,node_size,3))
        kpp_oblmixc       = 0.0_WP
        
        !_______________________________________________________________________
        ! read cvmix namelist file 
        nmlfile ='namelist.cvmix'    ! name of ocean namelist file
        ! check if cvmix namelist file exists if not use default values 
        inquire(file=trim(nmlfile),exist=nmlfile_exist) 
        if (nmlfile_exist) then
            open(20,file=trim(nmlfile))
                read(20,nml=param_kpp)
            close(20)
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
        call cvmix_init_kpp(Ri_crit        = kpp_Rib_crit,             & 
                            minOBLdepth    = kpp_minOBLdepth,         & 
                            minVtsqr       = kpp_minVtsqr,            &
                            vonKarman      = kpp_vonKarman,           & 
                            surf_layer_ext = kpp_surf_layer_ext,      &
                            interp_type    = kpp_interptype_ri,       &
                            interp_type2   = kpp_interptype_atobl,    &
                            lEkman         = kpp_use_compEkman,       &
                            lMonOb         = kpp_use_monob,           &
                            MatchTechnique = kpp_matchtechc,          &
                            lenhanced_diff = kpp_use_enhanceKv,       &
                            l_LMD_ws       = kpp_use_LMDws,           & 
                            lnonzero_surf_nonlocal = kpp_cs_is_one)
    end subroutine init_cvmix_kpp
    !
    !
    !
    !===========================================================================
    ! calculate PP vertrical mixing coefficients from CVMIX library
    subroutine calc_cvmix_kpp(mesh)
        type(t_mesh), intent(in)  , target :: mesh
        integer       :: node, elem, nz, nln, nun,  elnodes(3), aux_nz
        real(kind=WP) :: vshear2, dz2, aux, aux_wm(mesh%nl), aux_ws(mesh%nl)
        real(kind=WP) :: aux_coeff, sigma, stable
        real(kind=WP) :: aux_ustar, aux_surfbuoyflx_nl(mesh%nl)
        
        integer       :: nzsfc, nztmp
        real(kind=WP) :: sldepth, sfc_temp, sfc_salt, sfc_u, sfc_v, htot, delh, rho_sfc, rho_nz
        real(kind=WP) :: rhopot, bulk_0, bulk_pz, bulk_pz2
        real(kind=WP) :: sfc_rhopot, sfc_bulk_0, sfc_bulk_pz, sfc_bulk_pz2
#include "associate_mesh.h"
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
            
            !___3D Quantities___________________________________________________
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> call shear variables'//achar(27)//'[0m'
            if (kpp_use_fesomkpp) then
                !!PS do nz=2, nln
                do nz=nun+1, nln
                    !___________________________________________________________
                    ! calculate squared velocity shear referenced to the surface
                    ! --> cvmix wants to have it with  respect to the midlevel rather than full levels
                    !!PS kpp_dvsurf2(nz) = ((Unode(1,nz-1,node)+Unode(1,nz,node))*0.5 - Unode( 1,1,node) )**2 + &
                    !!PS                   ((Unode(2,nz-1,node)+Unode(2,nz,node))*0.5 - Unode( 2,1,node) )**2
                    kpp_dvsurf2(nz) = ((Unode(1,nz-1,node)+Unode(1,nz,node))*0.5 - Unode( 1,nun,node) )**2 + &
                                      ((Unode(2,nz-1,node)+Unode(2,nz,node))*0.5 - Unode( 2,nun,node) )**2
                    !___________________________________________________________
                    ! calculate shear Richardson number Ri = N^2/(du/dz)^2
                    dz2     = (Z_3d_n( nz-1,node)-Z_3d_n( nz,node))**2
                    vshear2 = (Unode(1,nz-1,node)-Unode(1,nz,node))**2 + &
                              (Unode(2,nz-1,node)-Unode(2,nz,node))**2 
                    vshear2 = vshear2/dz2
                    kpp_shearRi(nz) = max(bvfreq(nz,node),0.0_WP)/(vshear2+kpp_epsln)
                    
                    !_______________________________________________________________
                    ! buoyancy difference with respect to the surface --> computed in
                    ! oce_ale_pressure_bf.F90 --> subroutine pressure_bv 
                    ! --> dbsfc(nz,node)
                    !!PS call densityJM_components(tr_arr(1,node,1), tr_arr(1,node,2), sfc_bulk_0, sfc_bulk_pz, sfc_bulk_pz2, sfc_rhopot, mesh)
                    call densityJM_components(tr_arr(nun,node,1), tr_arr(nun,node,2), sfc_bulk_0, sfc_bulk_pz, sfc_bulk_pz2, sfc_rhopot, mesh)
                    call densityJM_components(tr_arr(nz,node,1), tr_arr(nz,node,2), bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)                    
                    rho_nz  = bulk_0   + Z_3d_n(nz,node)*(bulk_pz   + Z_3d_n(nz,node)*bulk_pz2)
                    rho_nz  = rho_nz*rhopot/(rho_nz+0.1_WP*Z_3d_n(nz,node))-density_0
                    rho_sfc = sfc_bulk_0   + Z_3d_n(nz,node)*(sfc_bulk_pz   + Z_3d_n(nz,node)*sfc_bulk_pz2)
                    rho_sfc = rho_sfc*sfc_rhopot/(rho_sfc+0.1_WP*Z_3d_n(nz,node))-density_0
                    kpp_dbsurf(nz) = -g * density_0_r *(rho_sfc-rho_nz)
                end do 
            else
                !!PS do nz=1, nln
                do nz=nun, nln
                    !___________________________________________________________
                    ! Calculate the surface layer depth, averaged surface layer 
                    ! quantities
                    
                    !!PS sldepth = kpp_surf_layer_ext*max( max(-Z_3d_n(nz,node),-zbar_3d_n(2,node)),kpp_minOBLdepth )
                    sldepth = kpp_surf_layer_ext*max( max(-Z_3d_n(nz,node),-zbar_3d_n(nun+1,node)),kpp_minOBLdepth )
                    nzsfc = nz
                    !!PS do nztmp = 1, nz
                    do nztmp = nun, nz
                        if (-zbar_3d_n(nztmp+1,node)>=sldepth) then
                            nzsfc = nztmp
                            exit
                        end if                        
                    end do
                    
                    ! average quantities over surface layers
                    sfc_temp = 0.0_WP
                    sfc_salt = 0.0_WP
                    sfc_u    = 0.0_WP
                    sfc_v    = 0.0_WP
                    htot     = 0.0_WP
                    !!PS do nztmp = 1, nzsfc
                    do nztmp = nun, nzsfc
                        delh     = min( max(0.0_WP,sldepth-htot), hnode(nztmp,node) )
                        htot     = htot+delh
                        sfc_temp = sfc_temp + tr_arr(nztmp,node,1)*delh
                        sfc_salt = sfc_salt + tr_arr(nztmp,node,2)*delh
                        sfc_u    = sfc_u    + Unode(1,nztmp,node) *delh
                        sfc_v    = sfc_v    + Unode(2,nztmp,node) *delh
                    end do
                    sfc_temp = sfc_temp/htot
                    sfc_salt = sfc_salt/htot
                    sfc_u    = sfc_u/htot
                    sfc_v    = sfc_v/htot
                    
                    !___________________________________________________________
                    ! calculate vertical shear between present layer and surface
                    ! averaged sfc_u and sfc_v
                    kpp_dvsurf2(nz) = (Unode(1,nz,node)-sfc_u)**2 + &
                                      (Unode(2,nz,node)-sfc_v)**2
                    
                    !___________________________________________________________
                    ! calculate buoyancy difference between the surface averaged 
                    ! and the grid points blow 
                    ! --> bring density of surface point adiabatically to the same 
                    !     depth level as the deep point --> than calculate bouyancy 
                    !     difference
                    call densityJM_components(sfc_temp, sfc_salt, sfc_bulk_0, sfc_bulk_pz, sfc_bulk_pz2, sfc_rhopot, mesh)
                    call densityJM_components(tr_arr(nz,node,1), tr_arr(nz,node,2), bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)                    
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
                    vshear2 = (Unode(1,nz-1,node)-Unode(1,nz,node))**2 + &
                              (Unode(2,nz-1,node)-Unode(2,nz,node))**2 
                    vshear2 = vshear2/dz2
                    kpp_shearRi(nz) = max(bvfreq(nz,node),0.0_WP)/(vshear2+kpp_epsln)
                end do ! --> do nz=1, nln
            end if
            
            !___2D Quantities___________________________________________________
            ! calculate surface bouyancy flux after eq. A2c & A2d & A3b & A3d 
            ! in Large et al. 1994
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> call surface buyflux[0m'
            !!PS kpp_sbuoyflx(node) = -g * &
            !!PS                         (sw_alpha(1,node)*heat_flux( node) / vcpw + &   !heat_flux & water_flux: positive up
            !!PS                          sw_beta( 1,node)*water_flux(node)*tr_arr(1,node,2))
            kpp_sbuoyflx(node) = -g * &
                                    (sw_alpha(nun,node)*heat_flux( node) / vcpw + &   !heat_flux & water_flux: positive up
                                     sw_beta( nun,node)*water_flux(node)*tr_arr(nun,node,2))
            
            
            ! calculate friction velocity (ustar) at surface (m/s)
            aux_ustar = sqrt( sqrt( stress_atmoce_x(node)**2 + stress_atmoce_y(node)**2 )*density_0_r ) ! @ the surface (eqn. 2)
            
            ! reduce friction velocity under ice --> approach take from the mpiom
            ! interface --> haven't been tested before with FESOM
            if (kpp_reduce_tauuice) then
                aux_ustar = aux_ustar*(1.0_WP-a_ice(node))**2
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
                !!PS do nz = 2, nln
                do nz = nun+1, nln
                    kpp_Av(nz,node) = kpp_pp_Av0     /(1.0_WP+kpp_pp_alpha*kpp_shearRi(nz))**kpp_pp_loc_exp + kpp_Avbckg
                    kpp_Kv(nz,node) = kpp_Av(nz,node)/(1.0_WP+kpp_pp_alpha*kpp_shearRi(nz)) ! + background diffusivity
                end do
            
            !___________________________________________________________________
            ! --> KPP shear parameterization below the mixed layer after (Large et
            ! al.,1994)     
            elseif (kpp_internalmix .eq. 'KPP') then
                !!PS do nz = 2, nln
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
                call par_ex
            end if 
            
            !___________________________________________________________________
            ! 3) Set Background diffusivities either to const. background 
            ! diffusivity (kpp_Kvbckg) or non-const. background diffusivity 
            ! of Qiang from FESOM1.4
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc background diff'//achar(27)//'[0m'
            if (kpp_use_nonconstKvb) then
                !!PS do nz = 2, nln
                do nz = nun+1, nln
                    call Kv0_background_qiang( &
                        aux,geo_coord_nod2D(2,node)/rad,abs(zbar_3d_n(nz,node))&
                        )
                    kpp_Kv(nz,node) = kpp_Kv(nz,node) + aux
                end do 
            else
                !!PS do nz = 2, nln
                do nz = nun+1, nln
                    kpp_Kv(nz,node) = kpp_Kv(nz,node) + kpp_Kvbckg
                end do 
            end if
            
            !___________________________________________________________________
            !   ||    ||    ||    ||    ||    ||    ||    ||    ||    ||    ||  
            !  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_ 
            !  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  / 
            !   \/    \/    \/    \/    \/    \/    \/    \/    \/    \/    \/  
            !___________________________________________________________________
            ! 4) Calculate the turbulent velocity scales w_s for scalars at 
            ! the cell centered mid-depth levels
            ! --> more equivalent to what fesom_kpp is doing --> depending only on
            ! surface bouyancy flux calcualted directly from heat, freshwater flux
            ! and contribution from shortwave penetration 
            ! eq. A2c & A2d & A3b & A3d in Large et al. 1994
            ! --> avoid if condition in depth loop --> should be a tiny bit faster 
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc turbulent vel. scale'//achar(27)//'[0m'
            aux_surfbuoyflx_nl = 0.0_WP
            !!PS aux_surfbuoyflx_nl(1:nln) = kpp_sbuoyflx(node)
            aux_surfbuoyflx_nl(nun:nln) = kpp_sbuoyflx(node)
            if (use_sw_pene) then
                ! coeffcient to transfer SW temp flux into buoyancy flux
                !!PS aux_coeff       = g*sw_alpha(1,node)
                aux_coeff       = g*sw_alpha(nun,node)
                
                !!PS do nz = 1, nln
                do nz = nun, nln
                    ! sw_3d is the temperature through the full depth levels into/
                    ! out off the tacervolume 
                    ! --> (sw_3d(1,node)-sw_3d(nz,node)) contains all the penetrated 
                    ! heatflux until level nz
                    aux_surfbuoyflx_nl(nz) = kpp_sbuoyflx(node)+aux_coeff*(sw_3d(nun,node)-sw_3d(nz+1,node))
                                                   ! look in oce_shortwave_pene.F90<--|
                                                   ! substract swsurf from surface 
                                                   ! heat flux 
                end do    
            end if     
                
            ! REMEMBER !!!: kpp_use_LMDws==.true. use Large et al way of computing 
            ! sigma ... 
            ! if ((surf_buoy_force.ge.cvmix_zero) .and. l_LMD_ws) then
            !    sigma_loc(:) = sigma_coord(:)
            ! else
            !    sigma_loc(:) = min(surf_layer_ext, sigma_coord(:))
            ! end if
            ! --> FESOM1.4/MOM5 approach
            if (kpp_use_LMDws .and. kpp_use_fesomkpp) then
                stable = 0.5_WP + SIGN( 0.5_WP, kpp_sbuoyflx(node) )
                sigma  = stable + ( 1.0_WP - stable ) * kpp_surf_layer_ext
            ! --> MOM6 approach
            else
                sigma = kpp_surf_layer_ext
            end if 
            
            !!PS call cvmix_kpp_compute_turbulent_scales(         &
            !!PS     sigma_coord     = sigma             ,        & ! (in)  sigma: Normalized surface layer depth
            !!PS     OBL_depth       = abs(Z_3d_n(1:nln,node)),   & ! (in)  Assume OBL depth (m) =  mid-depth level
            !!PS     surf_buoy_force = aux_surfbuoyflx_nl(1:nln), & ! (in)  surfce buoyancy flux (m2/s3) consider sw_pene
            !!PS     surf_fric_vel   = aux_ustar,                 & ! (in)  turbulent friction velocity at surface (m/s)
            !!PS     w_s             = kpp_ws_cntr(1:nln)    & ! (out) Turbulent velocity scale profile (m/s) for skalars
            !!PS    ) 
            call cvmix_kpp_compute_turbulent_scales(         &
                sigma_coord     = sigma             ,        & ! (in)  sigma: Normalized surface layer depth
                OBL_depth       = abs(Z_3d_n(nun:nln,node)),   & ! (in)  Assume OBL depth (m) =  mid-depth level
                surf_buoy_force = aux_surfbuoyflx_nl(nun:nln), & ! (in)  surfce buoyancy flux (m2/s3) consider sw_pene
                surf_fric_vel   = aux_ustar,                 & ! (in)  turbulent friction velocity at surface (m/s)
                w_s             = kpp_ws_cntr(nun:nln)    & ! (out) Turbulent velocity scale profile (m/s) for skalars
                )     
            
            !___________________________________________________________________
            ! 5) Calculate Bulk Richardson number at centers of cell 
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
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc bulk richards'//achar(27)//'[0m'
            !!PS kpp_bulkRi(1:nln+1) = cvmix_kpp_compute_bulk_Richardson(   &
            !!PS     zt_cntr         = Z_3d_n(     1:nln  ,node),  & ! (in) Depth of cell center (m)
            !!PS     delta_buoy_cntr = kpp_dbsurf( 1:nln),         & ! (in) Bulk buoyancy difference, Br-B(z) (1/s)
            !!PS     delta_Vsqr_cntr = kpp_dvsurf2(1:nln),         & ! (in) Square of resolved velocity difference (m2/s2)
            !!PS     ws_cntr         = kpp_ws_cntr(1:nln),         & ! (in) Turbulent velocity scale profile (m/s)
            !!PS     Nsqr_iface      = bvfreq(     1:nln+1,node)   & ! (in) Buoyancy frequency (1/s)
            !!PS     )
            kpp_bulkRi(nun:nln+1) = cvmix_kpp_compute_bulk_Richardson(   &
                zt_cntr         = Z_3d_n(     nun:nln  ,node),  & ! (in) Depth of cell center (m)
                delta_buoy_cntr = kpp_dbsurf( nun:nln),         & ! (in) Bulk buoyancy difference, Br-B(z) (1/s)
                delta_Vsqr_cntr = kpp_dvsurf2(nun:nln),         & ! (in) Square of resolved velocity difference (m2/s2)
                ws_cntr         = kpp_ws_cntr(nun:nln),         & ! (in) Turbulent velocity scale profile (m/s)
                Nsqr_iface      = bvfreq(     nun:nln+1,node)   & ! (in) Buoyancy frequency (1/s)
                )    
            
            !___________________________________________________________________
            ! 6) Compute depth of oceanic boundary layer (kpp_obldepth)
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
            aux_surfbuoyflx_nl(1) = kpp_sbuoyflx(node)
            if (use_sw_pene) then
                ! take total heatflux from shortwave radiation into account --> 
                ! here only needed to calculate monin-obukov mixing length
                !!PS aux_surfbuoyflx_nl(1) = aux_surfbuoyflx_nl(1)+aux_coeff*sw_3d(1,node)
                aux_surfbuoyflx_nl(1) = aux_surfbuoyflx_nl(1)+aux_coeff*sw_3d(nun,node) 
            end if 
            
            !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc obl depth'//achar(27)//'[0m'
            !!PS call cvmix_kpp_compute_OBL_depth(          &
            !!PS     Ri_bulk    = kpp_bulkRi(1:nln+1),      & ! (in) Bulk Richardson number dim=(ke+1)
            !!PS     zw_iface   = zbar_3d_n( 1:nln+1,node), & ! (in) Height of interfaces (m) dim=(ke+1)
            !!PS     zt_cntr    = Z_3d_n(    1:nln  ,node), & ! (in) Height of cell centers (m) dim=(ke)
            !!PS     surf_fric  = aux_ustar,                & ! (in) Turbulent friction velocity at surface (m/s) dim=1
            !!PS     surf_buoy  = aux_surfbuoyflx_nl(1),    & ! (in) Buoyancy flux at surface (m2/s3) dim=1
            !!PS     Coriolis   = coriolis_node(node),      & ! (in) Coriolis parameter (1/s) dim=1
            !!PS     OBL_depth  = kpp_obldepth(node),       & ! (out) OBL depth (m) dim=1
            !!PS     kOBL_depth = kpp_nzobldepth(node)      & ! (out) level (+fraction) of OBL extent dim=1
            !!PS     )
            call cvmix_kpp_compute_OBL_depth(          &
                Ri_bulk    = kpp_bulkRi(nun:nln+1),      & ! (in) Bulk Richardson number dim=(ke+1)
                zw_iface   = zbar_3d_n( nun:nln+1,node), & ! (in) Height of interfaces (m) dim=(ke+1)
                zt_cntr    = Z_3d_n(    nun:nln  ,node), & ! (in) Height of cell centers (m) dim=(ke)
                surf_fric  = aux_ustar,                & ! (in) Turbulent friction velocity at surface (m/s) dim=1
                surf_buoy  = aux_surfbuoyflx_nl(1),    & ! (in) Buoyancy flux at surface (m2/s3) dim=1
                Coriolis   = coriolis_node(node),      & ! (in) Coriolis parameter (1/s) dim=1
                OBL_depth  = kpp_obldepth(node),       & ! (out) OBL depth (m) dim=1
                kOBL_depth = kpp_nzobldepth(node)      & ! (out) level (+fraction) of OBL extent dim=1
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
            
            ! model level of OBL depth (note: float number)
            kpp_nzobldepth(node)= cvmix_kpp_compute_kOBL_depth(zbar_3d_n(:,node), &
                                                               Z_3d_n(:,node),    &
                                                               kpp_obldepth(node) )
            
            ! safety for kOBL
            !!PS if (kpp_nzobldepth(node) > nlevels_nod2D(node)) then
            !!PS     kpp_nzobldepth(node) = nlevels_nod2D(node)
            !!PS end if 
            if (kpp_nzobldepth(node) > nln+1) then
                kpp_nzobldepth(node) = nln+1
            end if 
            
            !___________________________________________________________________
            ! 7) Call CVMix/KPP to obtain OBL diffusivities, viscosities and non-
            ! local transports
            !___________________________________________________________________
            ! the dealing of shortwave penetration comes from the original kpp
            ! parameterisation of FESOM1.4 & FESOM2.0
            ! --> interpolate contribution that comes from shortwave penetration
            ! to the depth of the obldepth
            aux_surfbuoyflx_nl(1) = kpp_sbuoyflx(node)
            if (use_sw_pene .and. kpp_use_fesomkpp .eqv. .true.) then
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
            elseif (use_sw_pene .and. kpp_use_fesomkpp .eqv. .false.) then
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
            ! write kpp_Av and kpp_Kv that contain only the ocean internal 
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
            !  compute the turbulent diffusion coefficients
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
                surf_fric = aux_ustar,                 & ! (in)    SurfaceFriction:Turbulent friction velocity at surface (m/s)
                surf_buoy = aux_surfbuoyflx_nl(1),     & ! (in)    SurfaceBuoynacyForcing: Buoyancy flux at surface (m2/s3)
                nlev      = nlevels_nod2D(node)-1,     & ! (in)    nlev: Number of levels to compute coeffs for
                max_nlev  = nl-1,                      & ! (in)    max_lev: maximum vertical levels
                w_s       = aux_ws(:),                 & ! (out)   turbulent velocity scale for tracer
                w_m       = aux_wm(:)                  & ! (out)   turbulent velocity scale for momentum
                )
                !!PS w_s       = kpp_ws(:,node),            & ! (out)   turbulent velocity scale for tracer
                !!PS w_m       = kpp_wm(:,node)             & ! (out)   turbulent velocity scale for momentum
                !!PS )
            !___________________________________________________________________
            ! --> computation of nonlocal transport terms done in cvmix_coeffs_kpp
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
        ! 8) horizontal smoothing of the OBL mixing coefficient --> approach from
        ! original kpp parameterisation of FESOM1.4 & FESOM2.0
        !!PS if (flag_debug .and. mype==0)  print *, achar(27)//'[35m'//'         --> calc smooth kpp_oblmixc'//achar(27)//'[0m'    
        if (kpp_use_smoothblmc .and. kpp_use_fesomkpp) then
            call exchange_nod(kpp_oblmixc(:,:,1))
            call exchange_nod(kpp_oblmixc(:,:,2))
            call exchange_nod(kpp_oblmixc(:,:,3))
            do nz=1, 3
                !_______________________________________________________________
                ! all loops go over myDim_nod2D so no halo information --> for 
                ! smoothing haloinfo is required --> therefor exchange_nod
                call smooth_nod(kpp_oblmixc(:,:,nz), kpp_smoothblmc_nmb, mesh)
            end do
        end if 
        
        !_______________________________________________________________________
        !   ||    ||    ||    ||    ||    ||    ||    ||    ||    ||    ||  
        !  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_  _||_ 
        !  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  / 
        !   \/    \/    \/    \/    \/    \/    \/    \/    \/    \/    \/  
        !_______________________________________________________________________
        ! 9) combine KPP OBL mixing coefficents with internal mixing coefficients
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
        call exchange_nod(kpp_Kv)
        Kv = kpp_Kv
        
        !_______________________________________________________________________
        ! write out viscosities to FESOM2.0 --> viscosities for FESOM2.0 are 
        ! defined on elements --> interpolate therefor from nodes to elements
        call exchange_nod(kpp_Av)
        Av = 0.0_WP
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            !!PS do nz=2,nlevels(elem)-1
            do nz=ulevels(elem),nlevels(elem)-1
                Av(nz,elem) = sum(kpp_Av(nz,elnodes))/3.0_WP    ! (elementwise)                
            end do
        end do
    end subroutine calc_cvmix_kpp    
end module g_cvmix_kpp
