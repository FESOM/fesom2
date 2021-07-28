!
!
!===============================================================================
! module interface to FESOM2.0 for the CVMIX TKE vertical mixing scheme --> Its based
! on the module interface for MPIOM programed by Nils Brüggeman & Oliver Gutjahr
! This module calls subroutines from the CVMix library for the
! calculation of vertical mixing: TKE scheme 
!
! @see  Gaspar, P., Y. Grégoris, and J.-M. Lefevre, 1990
!       J. Geophys. Res., 95(C9), 16179–16193, doi:10.1029/JC095iC09p16179.
!
! @see  Blanke, B., P. Delecluse
!       J. Phys. Oceanogr., 23, 1363–1388. doi:10.1175/1520-0485(1993)023<1363:VOTTAO>2.0.CO;2
!
! written by Patrick Scholz, 10.05.2019
module g_cvmix_tke
    !___________________________________________________________________________
    ! module calls from cvmix library
    use cvmix_tke,         only: init_tke, cvmix_coeffs_tke
    use cvmix_put_get,     only: cvmix_put
    use cvmix_kinds_and_types
    use g_cvmix_idemix,    only: iwe, iwe_Tdis, iwe_alpha_c
    
    !___________________________________________________________________________
    ! module calls from FESOM
    use g_config , only: dt
    use o_param           
    use mod_mesh
    use g_parsup
    use o_arrays
    use g_comm_auto 
    implicit none
    
    !___________________________________________________________________________
    ! CVMIX-TKE namelist parameters
    real(kind=WP) :: tke_c_k        = 0.1              
    real(kind=WP) :: tke_c_eps      = 0.7               
    real(kind=WP) :: tke_alpha      = 30.0              
    real(kind=WP) :: tke_mxl_min    = 1.0e-8            
    real(kind=WP) :: tke_kappaM_min = 0.0                
    real(kind=WP) :: tke_kappaM_max = 100.0              
    ! real(kind=WP) :: cd             = 3.75 ! for Dirichlet boundary conditions
    real(kind=WP) :: tke_cd         = 1.0  ! for Neumann boundary conditions 
    real(kind=WP) :: tke_surf_min   = 1.0e-4             
    real(kind=WP) :: tke_min        = 1.0e-6
    
    ! tke_mxl_choice ... Can only be 1 or 2, choice of calculation of mixing 
    ! length; currently only Blanke, B., P. Delecluse option is implemented
    integer       :: tke_mxl_choice = 2 
    
    logical       :: tke_only             = .true.             
    logical       :: use_ubound_dirichlet = .false.           
    logical       :: use_lbound_dirichlet = .false.
    
    ! apply time relaxation to avo/dvo
    ! FIXME: nils: Do we need that
    logical       :: timerelax_tke = .false.
    
    real(kind=WP) :: relne=0.4       ! percentage of new value
    real(kind=WP) :: relax=0.6       ! percentage of old value 1-relne
    
    namelist /param_tke/ tke_c_k, tke_c_eps, tke_alpha, tke_mxl_min, tke_kappaM_min, tke_kappaM_max, &
                         tke_cd, tke_surf_min, tke_min, tke_mxl_choice, &
                         use_ubound_dirichlet, use_lbound_dirichlet, & 
                         timerelax_tke, relne, relax
     
    !___________________________________________________________________________
    ! CVMIX-TKE 3D variables
    
    ! dissipation of turbulent kinetic energy
    real(kind=WP), allocatable, dimension(:,:) :: tke_diss 
    
    ! forcing by internal wave energy (from idemix)
    real(kind=WP), allocatable, dimension(:,:) :: tke_in3d_iwe, tke_in3d_iwdis, tke_in3d_iwealphac
    
    ! dummy varibles
    real(kind=WP), allocatable, dimension(:,:) :: cvmix_dummy_1, cvmix_dummy_2, cvmix_dummy_3
    
    ! TKE diagnostics from Nils 
    real(kind=WP), allocatable, dimension(:,:) :: tke
    real(kind=WP), allocatable, dimension(:,:) :: tke_Lmix
    real(kind=WP), allocatable, dimension(:,:) :: tke_Pr  
    real(kind=WP), allocatable, dimension(:,:) :: tke_Tbpr
    real(kind=WP), allocatable, dimension(:,:) :: tke_Tspr
    real(kind=WP), allocatable, dimension(:,:) :: tke_Tdif
    real(kind=WP), allocatable, dimension(:,:) :: tke_Tdis
    real(kind=WP), allocatable, dimension(:,:) :: tke_Twin
    real(kind=WP), allocatable, dimension(:,:) :: tke_Tiwf
    real(kind=WP), allocatable, dimension(:,:) :: tke_Tbck
    real(kind=WP), allocatable, dimension(:,:) :: tke_Ttot                              
    real(kind=WP), allocatable, dimension(:,:) :: tke_Av, tke_Kv
    real(kind=WP), allocatable, dimension(:)   :: tke_Av_old, tke_Kv_old, tke_old
    
    ! forcing square of vertical velocity shear
!!PS     real(kind=WP), allocatable, dimension(:,:) :: tke_in3d_vshear2,tke_in3d_bvfreq2
    
    !! THIS should come from IDEMIX module when its writen
    !real(kind=WP), allocatable, dimension(:,:) :: iwe, iwe_Tdis, iwe_alpha_c
    
    !___________________________________________________________________________
    ! CVMIX-TKE 2D variables
    ! surface forcing of TKE (unit?)
    real(kind=WP), allocatable, dimension(:) :: tke_forc2d_normstress     
    ! buoyancy forcing at surface from heat/freshwater (unit?)
    real(kind=WP), allocatable, dimension(:) :: tke_forc2d_rhosurf         
    ! bottom friction (unit?)
    real(kind=WP), allocatable, dimension(:) :: tke_forc2d_botfrict

    ! nils
    integer :: tstep_count
  
    contains
    !
    !
    !
    !===========================================================================
    ! allocate and initialize TKE 2D and 3D variables --> call initialisation 
    ! routine from cvmix library
    subroutine init_cvmix_tke(mesh)
        implicit none
        character(len=cvmix_strlen)       :: nmlfile
        logical                  :: nmlfile_exist=.False.
        integer                  :: node_size
        type(t_mesh), intent(in), target :: mesh
#include "associate_mesh.h"
        !_______________________________________________________________________
        if(mype==0) then
            write(*,*) '____________________________________________________________'
            write(*,*) ' --> initialise CVMIX_TKE'
            write(*,*)
        end if
        
        !_______________________________________________________________________
        ! allocate + initialse tke arrays --> with size myDim_nod2D+eDim_nod2D
        node_size=myDim_nod2D+eDim_nod2D
        ! initialise TKE viscosities and diffusivities
        allocate(tke_Av(nl,node_size),tke_Kv(nl,node_size))
        tke_Av        = 0.0_WP
        tke_Kv        = 0.0_WP
        
        ! initialize 2D TKE input fields 
        allocate(tke_forc2d_normstress(node_size))
        allocate(tke_forc2d_botfrict(node_size))
        allocate(tke_forc2d_rhosurf(node_size))
        tke_forc2d_normstress= 0.0_WP
        tke_forc2d_botfrict  = 0.0_WP ! --> Nils B.: has yet no function is for later extension
        tke_forc2d_rhosurf   = 0.0_WP ! --> Nils B.: has yet no function is for later extension
        
        ! initialize 2D TKE input fields that come from IDEMIX
        allocate(tke_in3d_iwe(nl, node_size))
        allocate(tke_in3d_iwdis(nl, node_size))
        allocate(tke_in3d_iwealphac(nl, node_size))
        tke_in3d_iwe       = 0.0_WP
        tke_in3d_iwdis     = 0.0_WP
        tke_in3d_iwealphac = 0.0_WP
        
        ! initialise dummy values for debugging
        allocate(cvmix_dummy_1(nl,node_size))
        allocate(cvmix_dummy_2(nl,node_size))
        allocate(cvmix_dummy_3(nl,node_size))
        cvmix_dummy_1  = 0.0_WP
        cvmix_dummy_2  = 0.0_WP
        cvmix_dummy_3  = 0.0_WP
        
        ! initialise TKE diagnostic
        allocate(tke(nl,node_size))
        tke            = 0.0_WP
        
        allocate(tke_Pr(nl,node_size))
        allocate(tke_Lmix(nl,node_size))
        tke_Pr         = 0.0_WP
        tke_Lmix       = 0.0_WP
        
        allocate(tke_Tbpr(nl,node_size))
        allocate(tke_Tspr(nl,node_size))
        allocate(tke_Tdif(nl,node_size))
        allocate(tke_Tdis(nl,node_size))
        allocate(tke_Tiwf(nl,node_size))
        allocate(tke_Twin(nl,node_size))
        allocate(tke_Tbck(nl,node_size))
        allocate(tke_Ttot(nl,node_size))
        tke_Tbpr       = 0.0_WP
        tke_Tspr       = 0.0_WP
        tke_Tdif       = 0.0_WP
        tke_Tdis       = 0.0_WP
        tke_Twin       = 0.0_WP
        tke_Tiwf       = 0.0_WP
        tke_Tbck       = 0.0_WP
        tke_Ttot       = 0.0_WP
        
        ! nils (for debugging)
        tstep_count = 0
        
        !_______________________________________________________________________
        ! read cvmix namelist file 
        nmlfile ='namelist.cvmix'    ! name of ocean namelist file
        ! check if cvmix namelist file exists if not use default values 
        inquire(file=trim(nmlfile),exist=nmlfile_exist) 
        if (nmlfile_exist) then
            open(20,file=trim(nmlfile))
                read(20,nml=param_tke)
            close(20)
        else
            write(*,*) '     could not find namelist.cvmix, will use default values !'    
        end if
        
        !_______________________________________________________________________
        if(mix_scheme_nmb==56) tke_only=.False.
        
        if (mype==0) then
            write(*,*) "     tke_only       = ", tke_only
            write(*,*) "     tke_c_k        = ", tke_c_k
            write(*,*) "     tke_c_eps      = ", tke_c_eps
            write(*,*) "     tke_alpha      = ", tke_alpha
            write(*,*) "     tke_cd         = ", tke_cd
            write(*,*) "     tke_surf_min   = ", tke_surf_min
            write(*,*) "     tke_min        = ", tke_min
            write(*,*) "     tke_kappaM_min = ", tke_kappaM_min
            write(*,*) "     tke_kappaM_max = ", tke_kappaM_max
            write(*,*) "     tke_mxl_choice = ", tke_mxl_choice
            write(*,*)
        end if
        
        !_______________________________________________________________________
        ! call tke initialisation routine from cvmix library
        call init_tke(c_k            = tke_c_k,            &
                      c_eps          = tke_c_eps,          &
                      cd             = tke_cd,             &
                      alpha_tke      = tke_alpha,          &
                      mxl_min        = tke_mxl_min,        &
                      kappaM_min     = tke_kappaM_min,     &
                      kappaM_max     = tke_kappaM_max,     &
                      tke_mxl_choice = tke_mxl_choice,     &
                      use_ubound_dirichlet = use_ubound_dirichlet, &
                      use_lbound_dirichlet = use_lbound_dirichlet, &
                      only_tke       = tke_only,           &
                      tke_min        = tke_min,            &
                      tke_surf_min   = tke_surf_min    )
    end subroutine init_cvmix_tke
    !
    !
    !
    !===========================================================================
    ! calculate TKE vertical mixing coefficients from CVMIX library
    subroutine calc_cvmix_tke(mesh)
        implicit none
        type(t_mesh), intent(in), target :: mesh
        integer       :: node, elem, nelem, nz, nln, nun, elnodes(3), node_size
        real(kind=WP) :: tvol
        real(kind=WP) :: dz_trr(mesh%nl), bvfreq2(mesh%nl), vshear2(mesh%nl)
        real(kind=WP) :: tke_Av_old(mesh%nl), tke_Kv_old(mesh%nl), tke_old(mesh%nl)
        
#include "associate_mesh.h"

        node_size = myDim_nod2D
        !_______________________________________________________________________
        ! calculate all neccessary forcing for TKE 
        tke_forc2d_normstress = 0.0_WP
        tke_forc2d_botfrict   = 0.0_WP
        tke_forc2d_rhosurf    = 0.0_WP
        
        ! load things from idemix when selected
        if (.not. tke_only) then
            tke_in3d_iwe       = iwe
            tke_in3d_iwdis     = -iwe_Tdis
            tke_in3d_iwealphac = iwe_alpha_c
        endif
        
        !_______________________________________________________________________
        do node = 1,node_size
            !___________________________________________________________________
            ! number of above bottom levels at node
            nln = nlevels_nod2D(node)-1
            nun = ulevels_nod2D(node)
            
            !___________________________________________________________________
            ! calcualte for TKE surface momentum forcing --> norm of nodal 
            ! surface wind stress --> tke_forc2d_normstress --> interpolate from elements
            ! to nodes
            tke_forc2d_normstress(node) =   sqrt( &
                                                stress_node_surf(1,node)**2 + &
                                                stress_node_surf(2,node)**2 &
                                                )/density_0
            
            !___________________________________________________________________
            ! calculate for TKE 3D vertical velocity shear
            vshear2=0.0_WP
            do nz=nun+1,nln
                vshear2(nz)=(( Unode(1, nz-1, node) - Unode(1, nz, node))**2 + &
                             ( Unode(2, nz-1, node) - Unode(2, nz, node))**2)/ &
                             ((Z_3d_n(nz-1,node)-Z_3d_n(nz,node))**2)
            end do 
            
            !___________________________________________________________________
            ! calculate for TKE bottom friction and surface density forcing --> 
            ! NIELS: bottom friction variable should be included at a later stage 
            ! same as tke_forc2d_rhosurf
            ! tke_forc2d_botfrict(node) = 0.0_WP
            ! tke_forc2d_rhosurf(node)  = 0.0_WP
            
            !___________________________________________________________________
            ! calculate for TKE square of Brünt-Väisälä frequency, be aware that
            ! bvfreq contains already the squared brünt Väisälä frequency ...
            bvfreq2        = 0.0_WP
            !!PS bvfreq2(2:nln) = bvfreq(2:nln,node)
            bvfreq2(nun+1:nln) = bvfreq(nun+1:nln,node)
            
            !___________________________________________________________________
            ! dz_trr distance between tracer points, surface and bottom dz_trr is half 
            ! the layerthickness ...
            dz_trr            = 0.0_WP
            dz_trr(nun+1:nln) = abs(Z_3d_n(nun:nln-1,node)-Z_3d_n(nun+1:nln,node))
            dz_trr(nun)       = hnode(nun,node)/2.0_WP
            dz_trr(nln+1)     = hnode(nln,node)/2.0_WP
            
            !___________________________________________________________________
            ! main cvmix call to calculate tke
            tke_Av_old = tke_Av(:,node)
            tke_Kv_old = tke_Kv(:,node)
            tke_old    = tke(:,node)
            
            call cvmix_coeffs_tke(&
                ! parameter
                dzw          = hnode(nun:nln,node),               & ! distance between layer interface --> hnode
                dzt          = dz_trr(nun:nln+1),                   & ! distnace between tracer points
!                 nlev         = nln,                         &
                nlev         = nln-nun+1,                         &
                max_nlev     = nl-1,                        &
                dtime        = dt,                          &
                rho_ref      = density_0,                   &
                grav         = g,                           &
                ! essentials
                tke_new      = tke(       nun:nln+1,node),                 & ! out--> turbulent kinetic energy
                KappaM_out   = tke_Av(    nun:nln+1,node),              & ! out
                KappaH_out   = tke_Kv(    nun:nln+1,node),              & ! out
                tke_old      = tke_old(   nun:nln+1),                  & ! in --> turbulent kinetic energy previous time step
                old_KappaM   = tke_Av_old(nun:nln+1),               & ! in
                old_KappaH   = tke_Kv_old(nun:nln+1),               & ! in
                Ssqr         = vshear2(   nun:nln+1),                  & ! in --> square vert. vel. shear
                Nsqr         = bvfreq2(   nun:nln+1),                  & ! in --> square brunt Väisälä freq
                alpha_c      = tke_in3d_iwealphac(nun:nln+1,node),  & ! in for IDEMIX Ri
                E_iw         = tke_in3d_iwe(nun:nln+1,node),        & ! in for IDEMIX Ri
                ! forcing
                forc_tke_surf= tke_forc2d_normstress(   node), & ! in --> wind stress  
                forc_rho_surf= tke_forc2d_rhosurf(      node), & ! in
                bottom_fric  = tke_forc2d_botfrict(     node), & ! in
                iw_diss      = tke_in3d_iwdis(nun:nln+1,node), & ! in
                ! diagnostics
                tke_Tbpr     = tke_Tbpr(nun:nln+1,node),            & ! buoyancy production
                tke_Tspr     = tke_Tspr(nun:nln+1,node),            & ! shear production 
                tke_Tdif     = tke_Tdif(nun:nln+1,node),            & ! vertical diffusion d/dz(k d/dz)TKE
                tke_Tdis     = tke_Tdis(nun:nln+1,node),            & ! dissipation
                tke_Twin     = tke_Twin(nun:nln+1,node),            & ! wind forcing
                tke_Tiwf     = tke_Tiwf(nun:nln+1,node),            & ! internal wave forcing when idemix is used
                tke_Tbck     = tke_Tbck(nun:nln+1,node),            & ! background forcing only active if IDEMIX is not active, forcing that results from resetting TKE to minimum background TKE value
                tke_Ttot     = tke_Ttot(nun:nln+1,node),            & ! sum of all terms
                tke_Lmix     = tke_Lmix(nun:nln+1,node),            & ! mixing length scale of the TKE scheme
                tke_Pr       = tke_Pr(  nun:nln+1,node),              & ! Prantl number
                ! debugging
                cvmix_int_1  = cvmix_dummy_1(nun:nln+1,node),        & !
                cvmix_int_2  = cvmix_dummy_2(nun:nln+1,node),        & !
                cvmix_int_3  = cvmix_dummy_3(nun:nln+1,node),        & !
                i = 1,                                       &
                j = 1,                                       &
                tstep_count = tstep_count                    &
                )
            
            tke_Av(nln+1,node)=0.0_WP
            tke_Kv(nln+1,node)=0.0_WP
            tke_Av(nun  ,node)=0.0_WP
            tke_Kv(nun  ,node)=0.0_WP
            
        end do !--> do node = 1,node_size
        
        !_______________________________________________________________________
        ! write out diffusivity
        call exchange_nod(tke_Kv)
        Kv = tke_Kv
            
        !_______________________________________________________________________
        ! write out viscosity -->interpolate therefor from nodes to elements
        call exchange_nod(tke_Av) !Warning: don't forget to communicate before averaging on elements!!!
        Av = 0.0_WP
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            do nz=ulevels(elem)+1,nlevels(elem)-1
                Av(nz,elem) = sum(tke_Av(nz,elnodes))/3.0_WP    ! (elementwise)                
            end do
        end do
    end subroutine calc_cvmix_tke
end module g_cvmix_tke
