!
!
!===============================================================================
! module interface to FESOM2.0 for the CVMIX IDEMIX extension for the calculation 
! of hte internal wave energy and its dissipationof in Turbulent Kinetic Energy
! vertical mixing scheme --> Its based on the module interface for MPIOM programed 
! by Nils Brüggeman & Oliver Gutjahr, This module calls subroutines from the 
! CVMix library 
!
! @see Olbers D., Eden C.:
!       A Global Model for the Diapycnal Diffusivity Induced Internal Gravity Waves.
!       J. Phys. Oceanogr., 43, 1759-1779. doi: 10.1175/JPO-D-12-0207.1, 2013.
! @see Eden C., Czeschel L., Olbers D.:
!       Towards Energetically Consistent Ocean Models. 
!       J. Phys. Oceanogr., 44, 3160-3184, doi: 10.1175/JPO-D-13-0260.1, 2014.
!!PS !
!!PS !
!!PS !===============================================================================
!!PS module g_cvmix_idemix_param  
!!PS     use o_param, only: WP
!!PS     !___________________________________________________________________________
!!PS     ! OCECTL/CVMIX_IDEMIX_PARAM namelist parameters
!!PS     ! time scale for vertical symmetrisation (sec)
!!PS     real(kind=WP) :: tau_v = 86400.0
!!PS     
!!PS     ! time scale for horizontal symmetrisation, only necessary for lateral diffusion (sec)
!!PS     real(kind=WP) :: tau_h = 1296000.0
!!PS     
!!PS     ! constant of order one derived from the shape of the spectrum in m space (dimensionless)
!!PS     real(kind=WP) :: gamma = 1.570
!!PS     
!!PS     ! spectral bandwidth in modes (dimensionless)
!!PS     real(kind=WP) :: jstar = 10.0
!!PS     
!!PS     ! dissipation parameter (dimensionless)
!!PS     real(kind=WP) :: mu0   = 1.33333333
!!PS     
!!PS     integer       :: n_hor_iwe_prop_iter = 1
!!PS     
!!PS     namelist /cvmix_idemix/ tau_v, tau_h, gamma, jstar, mu0, n_hor_iwe_prop_iter
!!PS end module g_cvmix_idemix_param    
!
!
!===============================================================================
module g_cvmix_idemix
    
    !___________________________________________________________________________
    ! module calls from cvmix library
    use cvmix_idemix,  only :  init_idemix, calc_idemix_v0, cvmix_coeffs_idemix                 
    use cvmix_put_get, only : cvmix_put
    use cvmix_kinds_and_types 
    
    !___________________________________________________________________________
    ! module calls from FESOM
    use g_config , only: dt
    use o_param           
    use o_mesh
    use g_parsup
    use o_arrays
    use g_comm_auto 
    implicit none
    public
    
    !___________________________________________________________________________
    ! OCECTL/CVMIX_IDEMIX_PARAM namelist parameters
    ! time scale for vertical symmetrisation (sec)
    real(kind=WP) :: tau_v = 86400.0
    
    ! time scale for horizontal symmetrisation, only necessary for lateral diffusion (sec)
    real(kind=WP) :: tau_h = 1296000.0
    
    ! constant of order one derived from the shape of the spectrum in m space (dimensionless)
    real(kind=WP) :: gamma = 1.570
    
    ! spectral bandwidth in modes (dimensionless)
    real(kind=WP) :: jstar = 10.0
    
    ! dissipation parameter (dimensionless)
    real(kind=WP) :: mu0   = 1.33333333
    
    integer       :: n_hor_iwe_prop_iter = 1
    
    namelist /param_idemix/ tau_v, tau_h, gamma, jstar, mu0, n_hor_iwe_prop_iter
    
    !___________________________________________________________________________
    ! CVMIX-IDEMIX variables
    
    ! diagnostic fields
    real(kind=WP), allocatable, dimension(:,:) :: iwe
    real(kind=WP), allocatable, dimension(:,:) :: iwe_diss
    real(kind=WP), allocatable, dimension(:,:) :: iwe_alpha_c
    real(kind=WP), allocatable, dimension(:,:) :: iwe_c0
    real(kind=WP), allocatable, dimension(:,:) :: iwe_v0
    real(kind=WP), allocatable, dimension(:,:) :: cvmix_dummy_1
    real(kind=WP), allocatable, dimension(:,:) :: cvmix_dummy_2
    real(kind=WP), allocatable, dimension(:,:) :: cvmix_dummy_3
    real(kind=WP), allocatable, dimension(:,:) :: iwe_Ttot
    real(kind=WP), allocatable, dimension(:,:) :: iwe_Tdif
    real(kind=WP), allocatable, dimension(:,:) :: iwe_Thdi
    real(kind=WP), allocatable, dimension(:,:) :: iwe_Tdis
    real(kind=WP), allocatable, dimension(:,:) :: iwe_Tsur
    real(kind=WP), allocatable, dimension(:,:) :: iwe_Tbot
    real(kind=WP), allocatable, dimension(:,:) :: iwe_Av
    real(kind=WP), allocatable, dimension(:,:) :: iwe_Kv
    
    ! forcing square of Brunt Väisälä frequency
    real(kind=WP), allocatable, dimension(:,:) :: in3d_bvfreq2
    
    real(kind=WP), allocatable, dimension(:,:) :: weto_w
    real(kind=WP), allocatable, dimension(:,:) :: vol_wcell
    real(kind=WP), allocatable, dimension(:,:) :: vol_wcelli
    real(kind=WP), allocatable, dimension(:,:) :: Fx
    real(kind=WP), allocatable, dimension(:,:) :: Fy

    real(kind=WP), allocatable, dimension(:)   :: forc_iw_bottom_2D
    real(kind=WP), allocatable, dimension(:)   :: forc_iw_surface_2D

    real(kind=WP), dimension(:,:), pointer     :: dummy_3D
    real(kind=WP), dimension(:,:), pointer     :: dummy_2D

    ! load variables from CVMix list
    type(cvmix_data_type) :: CVMix_vars

    ! nils
    integer :: tstep_count

    contains
    !
    !
    !
    !===========================================================================
    ! allocate and initialize IDEMIX variables --> call initialisation 
    ! routine from cvmix library
    subroutine init_cvmix_idemix
        
        character(len=100) :: nmlfile
        logical            :: nmlfile_exist=.False.
        integer            :: node_size
        
        !_______________________________________________________________________
        ! allocate + initialse all idemix arrays
        node_size=myDim_nod2D+eDim_nod2D
        
        ! initialise 3D  IDEMIX fields
        allocate(iwe(nl,node_size))
        allocate(iwe_diss(nl,node_size))
        iwe(:,:)         = 0.0_WP
        iwe_diss(:,:)    = 0.0_WP
        
        allocate(cvmix_dummy_1(nl,node_size))
        allocate(cvmix_dummy_2(nl,node_size))
        allocate(cvmix_dummy_3(nl,node_size))
        cvmix_dummy_1(:,:) = 0.0_WP
        cvmix_dummy_2(:,:) = 0.0_WP
        cvmix_dummy_3(:,:) = 0.0_WP
        
        ! diagnostic 
        allocate(iwe_Ttot(nl,node_size))
        allocate(iwe_Tdif(nl,node_size))
        allocate(iwe_Thdi(nl,node_size))
        allocate(iwe_Tdis(nl,node_size))
        allocate(iwe_Tsur(nl,node_size))
        allocate(iwe_Tbot(nl,node_size))
        iwe_Ttot(:,:)    = 0.0_WP
        iwe_Tdif(:,:)    = 0.0_WP
        iwe_Thdi(:,:)    = 0.0_WP
        iwe_Tdis(:,:)    = 0.0_WP
        iwe_Tsur(:,:)    = 0.0_WP
        iwe_Tbot(:,:)    = 0.0_WP
        
        ! internal wave related vertical viscosity and diffusivity
        allocate(iwe_Av(nl,node_size))
        allocate(iwe_Kv(nl,node_size))
        iwe_Av(:,:)     = 0.0_WP
        iwe_Kv(:,:)     = 0.0_WP
        
        allocate(iwe_c0(nl,node_size))
        allocate(iwe_v0(nl,node_size))
        allocate(iwe_alpha_c(nl,node_size))
        iwe_c0(:,:)      = 0.0_WP
        iwe_v0(:,:)      = 0.0_WP
        iwe_alpha_c(:,:) = 0.0_WP
        
        allocate(weto_w(nl,node_size))
        allocate(vol_wcell(nl,node_size))
        allocate(vol_wcelli(nl,node_size))
        weto_w(:,:)      = 0.0_WP
        vol_wcell(:,:)   = 0.0_WP
        vol_wcelli(:,:)  = 0.0_WP
        
        allocate(Fx(nl,node_size))
        allocate(Fy(nl,node_size))
        Fx(:,:)          = 0.0_WP
        Fy(:,:)          = 0.0_WP
        
        ! 2D
        allocate(forc_iw_bottom_2D(node_size))
        allocate(forc_iw_surface_2D(node_size))
        forc_iw_bottom_2D(:) = 0.0_WP
        forc_iw_surface_2D(:)= 0.0_WP
        
        ! nils (for debugging)
        tstep_count = 0
        
        !_______________________________________________________________________
        ! read cvmix namelist file 
        nmlfile ='namelist.cvmix'    ! name of ocean namelist file
        ! check if cvmix namelist file exists if not use default values 
        inquire(file=trim(nmlfile),exist=nmlfile_exist) 
        if (nmlfile_exist) then
            open(20,file=trim(nmlfile))
                read(20,nml=param_idemix)
            close(20)
        end if    
        
        !_______________________________________________________________________
        ! initialise IDEMIX parameters
        call init_idemix(tau_v,tau_h,gamma,jstar,mu0)! ,handle_old_vals)! ,idemix_userdef_constants)
        
    end subroutine init_cvmix_idemix
    !
    !
    !
    !===========================================================================
    ! calculate IDEMIX internal wave energy and its dissipation
    subroutine calc_cvmix_idemix
    
        integer       :: node, elem, edge
        integer       :: nz, nln, nl1, nl2, nl12
        integer       :: elnodes1(3), elnodes2(3), el(2), ednodes(2) 
        real(kind=WP) :: dz_trr(nl), bvfreq2(nl), vflux, dz_el, aux, cflfac
        real(kind=WP) :: grad_v0Eiw(2), deltaX1, deltaY1, deltaX2, deltaY2
        logical       :: debug=.false.
        
        ! nils
        tstep_count = tstep_count + 1
        
        !_______________________________________________________________________
        ! Bottom and surface forcing
        if (tstep_count==1) then
            ! convert from W/m^2 to m^3/s^3
            forc_iw_bottom_2D  = forc_iw_bottom_2D/density_0
            ! only 20% of the niw-input are available to penetrate into the deeper ocean
            forc_iw_surface_2D = forc_iw_surface_2D/density_0 * 0.2 
        end if
        
        !_______________________________________________________________________
        do node = 1,myDim_nod2D
            nln = nlevels_nod2D(node)-1
            
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
            ! main call to calculate idemix
            call cvmix_coeffs_idemix(&
                ! parameter
                dzw             = hnode(1:nln,node),          &
                dzt             = dz_trr(1:nln+1),            &
                nlev            = nln,                        &
                max_nlev        = nl,                         & 
                dtime           = dt,                         &
                coriolis        = coriolis_node(node),        &
                ! essentials 
                iwe_old         = iwe(1:nln+1,node),          & ! in
                iwe_new         = iwe(1:nln+1,node),          & ! out
                forc_iw_surface = forc_iw_surface_2D(node),   & ! in
                forc_iw_bottom  = forc_iw_bottom_2D(node),    & ! in
                ! FIXME: nils: better output IDEMIX Ri directly
                alpha_c         = iwe_alpha_c(1:nln+1,node),  & ! out (for Ri IMIX)
                ! only for Osborn shortcut 
                ! FIXME: nils: put this to cvmix_tke
                KappaM_out      = iwe_Av(1:nln+1,node),      & ! out
                KappaH_out      = iwe_Kv(1:nln+1,node),      & ! out
                Nsqr            = bvfreq2(1:nln+1),           & ! in
                ! diagnostics
                iwe_Ttot        = iwe_Ttot(1:nln+1,node),     &
                iwe_Tdif        = iwe_Tdif(1:nln+1,node),     &
                iwe_Thdi        = iwe_Thdi(1:nln+1,node),     &
                iwe_Tdis        = iwe_Tdis(1:nln+1,node),     &
                iwe_Tsur        = iwe_Tsur(1:nln+1,node),     &
                iwe_Tbot        = iwe_Tbot(1:nln+1,node),     &
                c0              = iwe_c0(1:nln+1,node),       &
                v0              = iwe_v0(1:nln+1,node),       &
                ! debugging
                debug           = debug,                      &
                !i = i,                                        &
                !j = j,                                        &
                !tstep_count = tstep_count,                    &
                cvmix_int_1     = cvmix_dummy_1(1:nln+1,node),&
                cvmix_int_2     = cvmix_dummy_2(1:nln+1,node),&
                cvmix_int_3     = cvmix_dummy_3(1:nln+1,node) &
                )
            
        end do !-->do node = 1,myDim_nod2D
        
        !_______________________________________________________________________
        ! --> add contribution from horizontal wave propagation
        ! Since IDEMIX is used in a global model configuration the vertical 
        ! internal wave mixing (call cvmix_coeffs_idemix) have to be extended by 
        ! a lateral diffusion term (see. Olbers D., Eden C., 2013, A Global Model 
        ! for the Diapycnal Diffusivity Induced Internal Gravity Waves...)
        !
        ! diffusion term = div_h( v_0 * tau_h * grad_h(v_0*E_iw) )
        ! 
        ! use Gaussian integral satz ... int(div vec_A)dV = ringint(A*vec_n)dA
        !                                    div vec_A    = 1/V * sum_i=1...nface( A_i*vec_n_i)*A_i
        !
        if (n_hor_iwe_prop_iter>0) then
        
            ! make boundary exchange for iwe, and iwe_v0 --> for propagation need
            ! to calculate edge contribution that crosses the halo
            call exchange_nod(iwe)
            
            ! temporarily store old iwe values for diag
            iwe_Thdi = iwe
            
            !___________________________________________________________________
            ! calculate inverse volume and restrict iwe_v0 to fullfill stability 
            ! criterium --> CFL
            ! CFL Diffusion : CFL = v0^2 * dt/dx^2, CFL < 0.5
            ! --> limit v0 to CFL=0.2
            ! --> v0 = sqrt(CFL * dx^2 / dt)
            cflfac = 0.2_WP
            ! |--> FROM NILS: fac=0.2 ist geschätzt. Würde ich erstmal so 
            !      probieren. Der kommt aus dem stabilitätskriterium für Diffusion 
            !      (ähnlich berechnet wie das CFL Kriterium nur halt für den 
            !      Diffusions anstatt für den Advektionsterm). Normalerweise 
            !      sollte der Grenzwert aber nicht zu oft auftreten. Ich hatte 
            !      mal damit rum-experimentiert, aber letztendlich war die Lösung 
            !      das Iterativ zu machen und ggf. n_hor_iwe_prop_iter zu erhöhen. 
            !      Du kannst IDEMIX erstmal ohne den Term ausprobieren und sehen, 
            !      ob es läuft, dann kannst du den dazuschalten und hoffen, dass 
            !      es nicht explodiert. Eigentlich sollte der Term alles glatter 
            !      machen, aber nahe der ML kann der schon Probleme machen.
            do node = 1,myDim_nod2D
                nln = nlevels_nod2D(node)-1
                dz_trr         = 0.0_WP
                dz_trr(2:nln)  = Z_3d_n(1:nln-1,node)-Z_3d_n(2:nln,node)
                dz_trr(1)      = hnode(1,node)/2.0_WP
                dz_trr(nln+1)  = hnode(nln,node)/2.0_WP
                do nz=1,nln+1
                    ! inverse volumne 
                    vol_wcelli(nz,node) = area_inv(nz,node)/dz_trr(nz)
                    
                    ! restrict iwe_v0
                    aux = sqrt(cflfac*(area(nz,node)/pi*4.0_WP)/(tau_h*dt/n_hor_iwe_prop_iter))
                    !                  `--------+-------------´
                    !                           |-> comes from mesh_resolution=sqrt(area(1, :)/pi)*2._WP
                    iwe(nz,node) = min(iwe(nz,node),aux)
                end do 
            end do !-->do node = 1,myDim_nod2D
            call exchange_nod(vol_wcelli)
            call exchange_nod(iwe_v0)
            
            !___________________________________________________________________
            ! calculate horizontal diffusion term for internal wave energy
            do edge=1,myDim_edge2D
                !_______________________________________________________________
                ! deltaX1,deltaY1 ... cross edge dx and dy for the two elements 
                ! that contribute to the edge
                !                 o
                !                / \
                !               /   \       + ... triangle centroids
                !              /  +  \      x ... edge mid points
                !             /   ^   \
                !            /    |----\--------------->edge_cross_dxdy(1:2,edge)
                ! ednodes(1)o-----x-----o ednodes(2)
                !            \    |----/--------------->edge_cross_dxdy(3:4,edge)
                !             \   v   /
                !              \  +  /     
                !               \   /
                !                \ /
                !                 o
                deltaX1  = edge_cross_dxdy(1,edge)
                deltaY1  = edge_cross_dxdy(2,edge)
                ! ednode ... vertices that form the edgeiwe_Thdi
                ednodes  = edges(:,edge)
                ! el ... elements that contribute to edge
                el       = edge_tri(:,edge)
                ! elnodes1 ... nodes that contribute to element el(1) 
                elnodes1 = elem2d_nodes(:,el(1))
                ! nl1 ... number of layers at element el(1)
                nl1      = nlevels(el(1))
                
                !_______________________________________________________________
                ! the same as above but for el(2)--> if el(2)==0 than this edge 
                ! is a boundary edge
                nl2=0
                if (el(2)>0) then 
                    deltaX2=edge_cross_dxdy(3,edge)
                    deltaY2=edge_cross_dxdy(4,edge)
                    elnodes1= elem2d_nodes(:,el(2))
                    nl2=nlevels(el(2))
                endif
                
                !_______________________________________________________________
                ! goes only into this loop when the edge as two facing elements
                ! --> so the edge is not a boundary edge
                nl12=min(nl1,nl2)
                do nz=1,nl12
                    !___________________________________________________________
                    ! --> calc: grad_h(v_0*E_iw)
                    ! calculate grad(iwe*iwe_v0) for el(1)
                    grad_v0Eiw(1) = sum(gradient_sca(1:3,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    grad_v0Eiw(2) = sum(gradient_sca(4:6,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    
                    ! calculate grad(iwe*iwe_v0) for el(2) and average for el(1)
                    ! and el(2)
                    grad_v0Eiw(1) = (grad_v0Eiw(1) + sum(gradient_sca(1:3,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2)))*0.5_WP
                    grad_v0Eiw(2) = (grad_v0Eiw(2) + sum(gradient_sca(4:6,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2)))*0.5_WP
                    dz_el         = sum(helem(nz, el))*0.5_WP
                    
                    ! calculate flux 
                    vflux = ((deltaX2-deltaX1)*grad_v0Eiw(2)-(deltaY2-deltaY1)*grad_v0Eiw(1))*dz_el
                    
                    !___________________________________________________________
                    ! --> calc: v_0*tau_h* grad_h(v_0*E_iw)
                    ! multiply vflux with iwe_v0 interpolate to the edge-
                    ! mid point 
                    vflux = vflux * (iwe_v0(nz,ednodes(1))+iwe_v0(nz,ednodes(2)))*0.5_WP
                    
                    !___________________________________________________________
                    ! --> calc: div(v_0*tau_h* grad_h(v_0*E_iw))
                    ! sum fluxes over the surface --> gaussian integral satz
                    iwe(nz,ednodes(1)) = iwe(nz,ednodes(1)) + dt*tau_h/n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(1))*vflux
                    iwe(nz,ednodes(2)) = iwe(nz,ednodes(2)) - dt*tau_h/n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(2))*vflux
                end do !-->do nz=1,n2
                
                !_______________________________________________________________
                ! goes only into this loop when the edge has only facing element
                ! el(1) --> so the edge is a boundary edge
                do nz=nl12+1,nl1
                    !___________________________________________________________
                    ! --> calc: grad_h(v_0*E_iw)
                    ! calculate flux from el(1) with respect to edge mid 
                    ! point
                    grad_v0Eiw(1) = sum(gradient_sca(1:3,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    grad_v0Eiw(2) = sum(gradient_sca(4:6,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    dz_el         = helem(nz, el(1))
                    
                    ! calculate flux 
                    vflux = (grad_v0Eiw(1)*deltaY1-grad_v0Eiw(2)*deltaX1)*dz_el
                    
                    !___________________________________________________________
                    ! --> calc: v_0*tau_h* grad_h(v_0*E_iw)
                    ! multiply vflux with iwe_v0 interpolate to the edge-
                    ! mid point 
                    vflux = vflux * (iwe_v0(nz,ednodes(1))+iwe_v0(nz,ednodes(2)))*0.5_WP
                    
                    !___________________________________________________________
                    ! --> calc: div(v_0*tau_h* grad_h(v_0*E_iw))
                    ! sum fluxes over the surface --> gaussian integral satz
                    iwe(nz,ednodes(1)) = iwe(nz,ednodes(1)) + dt*tau_h/n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(1))*vflux
                    iwe(nz,ednodes(2)) = iwe(nz,ednodes(2)) - dt*tau_h/n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(2))*vflux
                end do !-->do nz=nl12+1,nl1
                
                !_______________________________________________________________
                ! goes only into this loop when the edge has only facing elemenmt
                ! el(2) --> so the edge is a boundary edge
                do nz=nl12+1,nl2
                    !___________________________________________________________
                    ! --> calc: grad_h(v_0*E_iw)
                    ! first calculate flux from el(1) with respect to edge mid 
                    ! point
                    grad_v0Eiw(1) = sum(gradient_sca(1:3,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2))
                    grad_v0Eiw(2) = sum(gradient_sca(4:6,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2))
                    dz_el         = helem(nz, el(2))
                    
                    ! calculate flux 
                    vflux = -(grad_v0Eiw(1)*deltaY2-grad_v0Eiw(2)*deltaX2)*dz_el
                    !       |--> minus sign comes from the fact that the the 
                    !            normal vectors (dx1,dy1) and (dx2,dy2) face 
                    !            in opposite direction (Right-Hand-Rule)
                    
                    !___________________________________________________________
                    ! --> calc: v_0*tau_h* grad_h(v_0*E_iw)
                    ! multiply vflux with iwe_v0 interpolate to the edge-
                    ! mid point 
                    vflux = vflux * (iwe_v0(nz,ednodes(1))+iwe_v0(nz,ednodes(2)))*0.5_WP
                    
                    !___________________________________________________________
                    ! --> calc: div(v_0*tau_h* grad_h(v_0*E_iw))
                    ! sum fluxes over the surface --> gaussian integral satz
                    iwe(nz,ednodes(1)) = iwe(nz,ednodes(1)) + dt*tau_h/n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(1))*vflux
                    iwe(nz,ednodes(2)) = iwe(nz,ednodes(2)) - dt*tau_h/n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(2))*vflux
                end do !-->do nz=nl12+1,nl1
            end do !-->do edge=1,myDim_edge2D
            
            ! diagnostic: add horizontal propgation to the total production rate
            ! of internal wave energy iwe_Tot
            iwe_Thdi = (iwe - iwe_Thdi)/dt
            iwe_Ttot = iwe_Ttot + iwe_Thdi
        end if !-->if (n_hor_iwe_prop_iter>0) then
        
        !_______________________________________________________________________
        ! write IDEMIX diffusivities and viscositie to FESOM only when IDEMIX is 
        ! used alone --> mostly for debuging --> otherwise TKE Av and Kv are use
        if(trim(mix_scheme)=='cvmix_IDEMIX') then 
            !___________________________________________________________________
            ! write out diffusivity
            Kv = iwe_Kv
            
            !___________________________________________________________________
            ! write out viscosity -->interpolate therefor from nodes to elements
            call exchange_nod(iwe_Av) !Warning: don't forget to communicate before averaging on elements!!!
            do elem=1, myDim_elem2D
                elnodes1=elem2D_nodes(:,elem)
                do nz=1,nlevels(elem)-1
                    Av(nz,elem) = sum(iwe_Av(nz,elnodes1))/3.0_WP    ! (elementwise)                
                end do
                Av(nlevels(elem),elem ) = Av(nlevels(elem)-1,elem )
            end do
        end if 
    end subroutine calc_cvmix_idemix
end module g_cvmix_idemix
