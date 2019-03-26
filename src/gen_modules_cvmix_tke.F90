!
!
!===============================================================================
! module interface to FESOM2.0 for the CVMIX TKE vertical mixing scheme --> Its based
! on the module interface for MPIOM programed by Nils Brüggeman & Oliver Gutjahr
! This module calls subroutines from the CVMix library for the
! calculation of vertical mixing: TKE scheme 
!
! @see  Gaspar, P., Y. Grégoris, and J.-M. Lefevre
!       J. Geophys. Res., 95(C9), 16179–16193, doi:10.1029/JC095iC09p16179.
!
! @see  Blanke, B., P. Delecluse
!       J. Phys. Oceanogr., 23, 1363–1388. doi:10.1175/1520-0485(1993)023<1363:VOTTAO>2.0.CO;2
!
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
    use o_mesh
    use g_parsup
    use o_arrays
    use g_comm_auto 
    implicit none
    
    !___________________________________________________________________________
    ! CVMIX-TKE namelist parameters
    real(kind=WP) :: c_k          = 0.1              
    real(kind=WP) :: c_eps        = 0.7               
    real(kind=WP) :: alpha_tke    = 30.0              
    real(kind=WP) :: mxl_min      = 1.0e-8            
    real(kind=WP) :: kappaM_min   = 0.0                
    real(kind=WP) :: kappaM_max   = 100.0              
    ! real(kind=WP) :: cd         = 3.75 ! for Dirichlet boundary conditions
    real(kind=WP) :: cd           = 1.0  ! for Neumann boundary conditions 
    real(kind=WP) :: tke_surf_min = 1.0e-4             
    real(kind=WP) :: tke_min      = 1.0e-6
    
    ! tke_mxl_choice ... Can only be 1 or 2, choice of calculation of mixing 
    ! length; currently only Blanke, B., P. Delecluse option is implemented
    integer       :: tke_mxl_choice = 2 
    
    logical       :: l_tke_active         = .true.             
    logical       :: only_tke             = .true.             
    logical       :: use_ubound_dirichlet = .false.           
    logical       :: use_lbound_dirichlet = .false.
    
    ! apply time relaxation to avo/dvo
    ! FIXME: nils: Do we need that
    logical       :: timerelax_tke = .false.           
    
    real(kind=WP) :: relne=0.4       ! percentage of new value
    real(kind=WP) :: relax=0.6       ! percentage of old value 1-relne
    
    namelist /param_tke/ c_k, c_eps, alpha_tke, mxl_min, kappaM_min, kappaM_max, &
                         cd, tke_surf_min, tke_min, tke_mxl_choice, &
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
    real(kind=WP), allocatable, dimension(:,:) :: tke_Av, tke_Av_old
    real(kind=WP), allocatable, dimension(:,:) :: tke_Kv, tke_Kv_old
    
    ! forcing square of vertical velocity shear
    real(kind=WP), allocatable, dimension(:,:) :: in3d_vshear2
    
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
    subroutine init_cvmix_tke
        character(len=100) :: nmlfile
        logical            :: nmlfile_exist=.False.
        integer            :: node_size
        
        !_______________________________________________________________________
        ! allocate + initialse all tke arrays
        node_size=myDim_nod2D+eDim_nod2D
        
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
        
        ! initialise TKE viscosities and diffusivities
        allocate(tke_Av(nl,node_size),tke_Av_old(nl,node_size))
        allocate(tke_Kv(nl,node_size),tke_Kv_old(nl,node_size))
        tke_Av        = 0.0_WP
        tke_Kv        = 0.0_WP
        tke_Av_old    = 0.0_WP
        tke_Kv_old    = 0.0_WP
        
        ! nils (for debugging)
        tstep_count = 0
        
!!PS         !_______________________________________________________________________
!!PS         ! THIS SHOULD COME FROM IDEMIX --> MUST BE FIXED WHEN IDEMIX MOUDLE IS 
!!PS         ! WRITTEN 
!!PS         allocate(iwe(nl,node_size))
!!PS         allocate(iwe_Tdis(nl,node_size))
!!PS         allocate(iwe_alpha_c(nl,node_size))
!!PS         iwe              = 0.0_WP
!!PS         iwe_Tdis         = 0.0_WP
!!PS         iwe_alpha_c      = 0.0_WP
        
        !_______________________________________________________________________
        ! read cvmix namelist file 
        nmlfile ='namelist.cvmix'    ! name of ocean namelist file
        ! check if cvmix namelist file exists if not use default values 
        inquire(file=trim(nmlfile),exist=nmlfile_exist) 
        if (nmlfile_exist) then
            open(20,file=trim(nmlfile))
                read(20,nml=param_tke)
            close(20)
        end if  
        
        !_______________________________________________________________________
        if(trim(mix_scheme)=='cvmix_TKE+IDEMIX') only_tke=.False.
        
        !_______________________________________________________________________
        ! call tke initialisation routine from cvmix library
        call init_tke(c_k            = c_k,            &
                      c_eps          = c_eps,          &
                      cd             = cd,             &
                      alpha_tke      = alpha_tke,      &
                      mxl_min        = mxl_min,        &
                      kappaM_min     = kappaM_min,     &
                      kappaM_max     = kappaM_max,     &
                      tke_mxl_choice = tke_mxl_choice, &
                      use_ubound_dirichlet = use_ubound_dirichlet, &
                      use_lbound_dirichlet = use_lbound_dirichlet, &
                      only_tke       = only_tke,       &
                      tke_min        = tke_min,        &
                      tke_surf_min   = tke_surf_min    )
            
    end subroutine init_cvmix_tke
    !
    !
    !
    !===========================================================================
    ! calculate TKE vertrical mixing coefficients from CVMIX library
    subroutine calc_cvmix_tke
        
        integer       :: node, elem, nelem, nz, nln, elnodes(3)
        real(kind=WP) :: tvol
        real(kind=WP) :: dz_trr(nl), bvfreq2(nl), vshear2(nl)
        
        !_______________________________________________________________________
        ! calculate all neccessary forcing for TKE 
        tke_forc2d_normstress = 0.0_WP
        tke_forc2d_botfrict   = 0.0_WP
        tke_forc2d_rhosurf    = 0.0_WP
        in3d_vshear2          = 0.0_WP
        
        ! load things from idemix when selected
        if (.not. only_tke) then
            tke_in3d_iwe       = iwe
            tke_in3d_iwdis     = -iwe_Tdis
            tke_in3d_iwealphac = iwe_alpha_c
        endif
        
        !_______________________________________________________________________
        do node = 1,myDim_nod2D
            !___________________________________________________________________
            ! calcualte for TKE surface momentum forcing --> norm of nodal 
            ! surface wind stress --> tke_forc2d_normstress --> interpolate from elements
            ! to nodes
            tvol = 0.0_WP
            do nelem=1,nod_in_elem2D_num(node)
                elem = nod_in_elem2D(nelem,node)
                tvol = tvol + elem_area(elem)
                tke_forc2d_normstress(node) = tke_forc2d_normstress(node) &
                                         + sqrt(stress_surf(1,elem)**2 + stress_surf(1,elem)**2)*elem_area(elem)
            end do !--> do nelem=1,nod_in_elem2D_num(node)
            tke_forc2d_normstress(node) = tke_forc2d_normstress(node)/tvol
                
            !___________________________________________________________________
            ! calculate for TKE 3D vertical velocity shear
            nln = nlevels_nod2D(node)-1
            vshear2=0.0_WP
            do nz=2,nln
                vshear2(nz)=(( Unode(1, nz-1, node) - Unode(1, nz, node))**2 + &
                             ( Unode(2, nz-1, node) - Unode(2, nz, node))**2)/ &
                             ((Z_3d_n(nz-1,node)-Z_3d_n(nz,node))**2)
            end do 
            ! in3d_vshear2(1)    = in3d_vshear2(2)
            ! in3d_vshear2(nln+1)= in3d_vshear2(nln)
            
            !___________________________________________________________________
            ! calculate for TKE bottom friction --> after niels bottom friction 
            ! variable should be included at a later stage same as tke_forc2d_rhosurf
            ! tke_forc2d_botfrict(node) = 0.0_WP
            
            !___________________________________________________________________
            ! calculate for TKE square of Brünt-Väisälä frequency
            bvfreq2          = 0.0_WP
            bvfreq2(1:nln+1) = bvfreq(1:nln+1,node)**2
            
            !___________________________________________________________________
            ! dz_trr distance between tracer points, surface and bottom dz_trr is half 
            ! the layerthickness ...
            dz_trr         = 0.0_WP
            dz_trr(2:nln)  = Z_3d_n(1:nln-1,node)-Z_3d_n(2:nln,node)
            dz_trr(1)      = hnode(1,node)/2.0_WP
            dz_trr(nln+1)  = hnode(nln,node)/2.0_WP
            
            !___________________________________________________________________
            ! main cvmix call to calculate tke
            tke_Av_old(:,node) = tke_Av(:,node)
            tke_Kv_old(:,node) = tke_Kv(:,node)
            call cvmix_coeffs_tke(&
                ! parameter
                dzw          = hnode(1:nln,node),             & ! distance between layer interface --> hnode
                dzt          = dz_trr(1:nln+1),               & ! distnace between tracer points
                nlev         = nln,                           &
                max_nlev     = nl,                            &
                dtime        = dt,                            &
                rho_ref      = density_0,                     &
                grav         = g,                             &
                ! essentials
                tke_new      = tke(1:nln+1,node),             & ! out--> turbulent kinetic energy
                KappaM_out   = tke_Av(1:nln+1,node),          & ! out
                KappaH_out   = tke_Kv(1:nln+1,node),          & ! out
                tke_old      = tke(1:nln+1,node),             & ! in --> turbulent kinetic energy previous time step
                Ssqr         = vshear2(1:nln+1),              & ! in --> square vert. vel. shear
                Nsqr         = bvfreq2(1:nln+1),              & ! in --> square brunt Väisälä freq
                old_kappaM   = tke_Av_old(1:nln+1,node),      & ! in
                old_KappaH   = tke_Kv_old(1:nln+1,node),      & ! in
                alpha_c      = tke_in3d_iwealphac(:,node),    & ! in for IDEMIX Ri
                E_iw         = tke_in3d_iwe(:,node),          & ! in for IDEMIX Ri
                ! forcing
                forc_tke_surf= tke_forc2d_normstress(node),   & ! in --> wind stress  
                forc_rho_surf= tke_forc2d_rhosurf(node),      & ! in
                bottom_fric  = tke_forc2d_botfrict(node),     & ! in
                iw_diss      = tke_in3d_iwdis(:,node),        & ! in
                ! diagnostics
                tke_Tbpr     = tke_Tbpr(1:nln+1,node),        & ! buoyancy production
                tke_Tspr     = tke_Tspr(1:nln+1,node),        & ! shear production 
                tke_Tdif     = tke_Tdif(1:nln+1,node),        & ! vertical diffusion d/dz(k d/dz)TKE
                tke_Tdis     = tke_Tdis(1:nln+1,node),        & ! dissipation
                tke_Twin     = tke_Twin(1:nln+1,node),        & ! wind forcing
                tke_Tiwf     = tke_Tiwf(1:nln+1,node),        & ! internal wave forcing when idemix is used
                tke_Tbck     = tke_Tbck(1:nln+1,node),        & ! background forcing only active if IDEMIX is not active, forcing that results from resetting TKE to minimum background TKE value
                tke_Ttot     = tke_Ttot(1:nln+1,node),        & ! sum of all terms
                tke_Lmix     = tke_Lmix(1:nln+1,node),        & ! mixing length scale of the TKE scheme
                tke_Pr       = tke_Pr(1:nln+1,node),          & ! Prantl number
               ! debugging
               cvmix_int_1  = cvmix_dummy_1(1:nln+1,node),       & !
               cvmix_int_2  = cvmix_dummy_2(1:nln+1,node),       & !
               cvmix_int_3  = cvmix_dummy_3(1:nln+1,node),       & !
               i = 1,                                            &
               j = 1,                                            &
               tstep_count = tstep_count                         &
               )
      
        end do !--> do node = 1,myDim_nod2D
        
        !_______________________________________________________________________
        ! Time relaxation of avo/dvo
        ! FIXME: nils: Why should we want to do this???
        !if ( timerelax_tke) THEN
        !  
        !    !___________________________________________________________________
        !    ! write out diffusivity
        !    Kv = relax*tke_Kv_old + relne*tke_Kv
        !    
        !    !___________________________________________________________________
        !    ! write out viscosity -->interpolate therefor from nodes to elements
        !    call exchange_nod(tke_Av) !Warning: don't forget to communicate before averaging on elements!!!
        !    do elem=1, myDim_elem2D
        !        elnodes=elem2D_nodes(:,elem)
        !        do nz=1,nlevels(elem)-1
        !            Av(nz,elem) = sum(relax*tke_Av_old(nz,elnodes)+relne*tke_Av(nz,elnodes))/3.0_WP    ! (elementwise)                
        !        end do
        !        Av(nlevels(elem),elem ) = Av(nlevels(elem)-1,elem )
        !    end do
        !    
        !else
            !___________________________________________________________________
            ! write out diffusivity
            Kv = tke_Kv
            
            !___________________________________________________________________
            ! write out viscosity -->interpolate therefor from nodes to elements
            call exchange_nod(tke_Av) !Warning: don't forget to communicate before averaging on elements!!!
            do elem=1, myDim_elem2D
                elnodes=elem2D_nodes(:,elem)
                do nz=1,nlevels(elem)-1
                    Av(nz,elem) = sum(tke_Av(nz,elnodes))/3.0_WP    ! (elementwise)                
                end do
                Av(nlevels(elem),elem ) = Av(nlevels(elem)-1,elem )
            end do
        !end
        
    end subroutine calc_cvmix_tke
end module g_cvmix_tke
