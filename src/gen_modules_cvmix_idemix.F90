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
! written by Patrick Scholz, 10.05.2019
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
    use mod_mesh
    use g_parsup
    use o_arrays
    use g_comm_auto 
    use g_read_other_NetCDF
    implicit none
    public
    
    !___________________________________________________________________________
    ! OCECTL/CVMIX_IDEMIX_PARAM namelist parameters
    ! time scale for vertical symmetrisation (sec)
    real(kind=WP) :: idemix_tau_v = 86400.0
    
    ! time scale for horizontal symmetrisation, only necessary for lateral diffusion (sec)
    real(kind=WP) :: idemix_tau_h = 1296000.0
    
    ! constant of order one derived from the shape of the spectrum in m space (dimensionless)
    real(kind=WP) :: idemix_gamma = 1.570
    
    ! spectral bandwidth in modes (dimensionless)
    real(kind=WP) :: idemix_jstar = 10.0
    
    ! dissipation parameter (dimensionless)
    real(kind=WP) :: idemix_mu0   = 1.33333333

    ! amount of surface forcing that is used
    real(kind=WP) :: idemix_sforcusage = 0.2
    
    
    integer       :: idemix_n_hor_iwe_prop_iter = 1
    
    ! filelocation for idemix surface forcing
    character(MAX_PATH):: idemix_surforc_file = './fourier_smooth_2005_cfsr_inert_rgrid.nc'

    ! filelocation for idemix bottom forcing
    character(MAX_PATH):: idemix_botforc_file = './tidal_energy_gx1v6_20090205_rgrid.nc'

    namelist /param_idemix/ idemix_tau_v, idemix_tau_h, idemix_gamma, idemix_jstar, idemix_mu0, idemix_n_hor_iwe_prop_iter, &
                            idemix_sforcusage, idemix_surforc_file, idemix_botforc_file
                            
                            
    
    !___________________________________________________________________________
    ! CVMIX-IDEMIX variables
    
    ! diagnostic fields
    real(kind=WP), allocatable, dimension(:,:) :: iwe
    real(kind=WP), allocatable, dimension(:)   :: iwe_old
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
    
    real(kind=WP), allocatable, dimension(:,:) :: vol_wcelli
    
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
    subroutine init_cvmix_idemix(mesh)
        implicit none
        character(len=cvmix_strlen) :: nmlfile
        logical                  :: file_exist=.False.
        integer                  :: node_size

        type(t_mesh), intent(in), target :: mesh

#include "associate_mesh.h"     
        !_______________________________________________________________________
        if(mype==0) then
            write(*,*) '____________________________________________________________'
            write(*,*) ' --> initialise IDEMIX'
            write(*,*)
        end if
            
        !_______________________________________________________________________
        ! allocate + initialse idemix arrays --> with size myDim_nod2D+eDim_nod2D
        node_size=myDim_nod2D+eDim_nod2D
        allocate(iwe(nl,node_size))
        iwe(:,:)         = 0.0_WP
        allocate(iwe_v0(nl,node_size))
        iwe_v0(:,:)      = 0.0_WP
        allocate(vol_wcelli(nl,node_size))
        vol_wcelli(:,:)  = 0.0_WP
        
        ! internal wave related vertical viscosity and diffusivity
        allocate(iwe_Av(nl,node_size))
        allocate(iwe_Kv(nl,node_size))
        iwe_Av(:,:)     = 0.0_WP
        iwe_Kv(:,:)     = 0.0_WP
        
        allocate(iwe_old(nl))
        allocate(iwe_diss(nl,node_size))
        iwe_old(:)       = 0.0_WP
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
        
        allocate(iwe_c0(nl,node_size))
        allocate(iwe_alpha_c(nl,node_size))
        iwe_c0(:,:)      = 0.0_WP
        iwe_alpha_c(:,:) = 0.0_WP
        
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
        file_exist=.False.
        inquire(file=trim(nmlfile),exist=file_exist) 
        if (file_exist) then
            open(20,file=trim(nmlfile))
                read(20,nml=param_idemix)
            close(20)
        else
            write(*,*) '     could not find namelist.cvmix, will use default values !'    
        end if    
        
        !_______________________________________________________________________
        if (mype==0) then
            write(*,*) "     idemix_tau_v        = ", idemix_tau_v
            write(*,*) "     idemix_tau_h        = ", idemix_tau_h
            write(*,*) "     idemix_gamma        = ", idemix_gamma
            write(*,*) "     idemix_jstar        = ", idemix_jstar
            write(*,*) "     idemix_mu0          = ", idemix_mu0
            write(*,*) "     idemix_n_hor_iwe_...= ", idemix_n_hor_iwe_prop_iter
            write(*,*) "     idemix_surforc_file = ", idemix_surforc_file
            write(*,*) "     idemix_botforc_file = ", idemix_botforc_file
            write(*,*)
        end if
        
        !_______________________________________________________________________
        ! read idemix surface near inertial wave forcing from cfsr data --> file 
        ! from N. Brüggemann interpoalted to regular grid
        file_exist=.False.
        inquire(file=trim(idemix_surforc_file),exist=file_exist) 
        if (file_exist) then
            if (mype==0) write(*,*) ' --> read IDEMIX near inertial wave surface forcing'
            call read_other_NetCDF(idemix_surforc_file, 'var706', 1, forc_iw_surface_2D, .true., mesh) 
            ! only 20% of the niw-input are available to penetrate into the deeper ocean
            forc_iw_surface_2D = forc_iw_surface_2D/density_0 * idemix_sforcusage 
            
        else
            if (mype==0) then
                write(*,*) '____________________________________________________________________'
                write(*,*) ' ERROR: IDEMIX surface forcing file not found! Cant apply IDEMIX'
                write(*,*) '        vertical mixing parameterisation! '
                write(*,*) '        --> check your namelist.cvmix, idemix_surforc_file &  '
                write(*,*) '            idemix_botforc_file'
                write(*,*) '____________________________________________________________________'
            end if
            call par_ex(0)
        end if 
        
        !_______________________________________________________________________
        ! read idemix bottom near tidal forcing from cesm data set --> file 
        ! from N. Brüggemann interpoalted to regular grid
        file_exist=.False.
        inquire(file=trim(idemix_surforc_file),exist=file_exist) 
        if (file_exist) then
            if (mype==0) write(*,*) ' --> read IDEMIX near tidal bottom forcing'
            call read_other_NetCDF(idemix_botforc_file, 'wave_dissipation', 1, forc_iw_bottom_2D, .true., mesh) 
            ! convert from W/m^2 to m^3/s^3
            forc_iw_bottom_2D  = forc_iw_bottom_2D/density_0
            
        else
            if (mype==0) then
                write(*,*) '____________________________________________________________________'
                write(*,*) ' ERROR: IDEMIX bottom forcing file not found! Cant apply IDEMIX'
                write(*,*) '        vertical mixing parameterisation! '
                write(*,*) '        --> check your namelist.cvmix, idemix_surforc_file &  '
                write(*,*) '            idemix_botforc_file'
                write(*,*) '____________________________________________________________________'
            end if 
            call par_ex(0)
        end if 
        
        !_______________________________________________________________________
        ! initialise IDEMIX parameters
        call init_idemix(idemix_tau_v,idemix_tau_h,idemix_gamma,idemix_jstar,idemix_mu0)! ,handle_old_vals)! ,idemix_userdef_constants)
    end subroutine init_cvmix_idemix
    !
    !
    !
    !===========================================================================
    ! calculate IDEMIX internal wave energy and its dissipation
    subroutine calc_cvmix_idemix(mesh)
        implicit none
        type(t_mesh), intent(in), target :: mesh
        integer       :: node, elem, edge, node_size
        integer       :: nz, nln, nl1, nl2, nl12, nu1, nu2, nu12, uln 
        integer       :: elnodes1(3), elnodes2(3), el(2), ednodes(2) 
        real(kind=WP) :: dz_trr(mesh%nl), dz_trr2(mesh%nl), bvfreq2(mesh%nl), vflux, dz_el, aux, cflfac
        real(kind=WP) :: grad_v0Eiw(2), deltaX1, deltaY1, deltaX2, deltaY2
        logical       :: debug=.false.

#include "associate_mesh.h"   
        ! nils
        tstep_count = tstep_count + 1
          
        !_______________________________________________________________________
        node_size = myDim_nod2D
        do node = 1,node_size
            nln = nlevels_nod2D(node)-1
            uln = ulevels_nod2D(node)
            
            !___________________________________________________________________
            ! calculate for TKE square of Brünt-Väisälä frequency, be aware that
            ! bvfreq contains already the squared brünt Väisälä frequency ...
            bvfreq2        = 0.0_WP
            bvfreq2(uln:nln) = bvfreq(uln:nln,node)
            
            !___________________________________________________________________
            ! dz_trr distance between tracer points, surface and bottom dz_trr is half 
            ! the layerthickness ...
            dz_trr         = 0.0_WP
            dz_trr(uln+1:nln)  = abs(Z_3d_n(uln:nln-1,node)-Z_3d_n(uln+1:nln,node))
            dz_trr(uln)      = hnode(uln,node)/2.0_WP
            dz_trr(nln+1)  = hnode(nln,node)/2.0_WP
            
            !___________________________________________________________________
            ! main call to calculate idemix
            iwe_old = iwe(:,node)
            
            call cvmix_coeffs_idemix(&
                ! parameter
                dzw             = hnode(:,node),                &
                dzt             = dz_trr(:),                    &
                nlev            = nln,                          &
                max_nlev        = nl-1,                         &  
                dtime           = dt,                           &
                coriolis        = coriolis_node(node),          &
                ! essentials 
                iwe_new         = iwe(:,node),                  & ! out
                iwe_old         = iwe_old(:),                   & ! in
                forc_iw_surface = forc_iw_surface_2D(node),     & ! in
                forc_iw_bottom  = forc_iw_bottom_2D(node),      & ! in
                ! FIXME: nils: better output IDEMIX Ri directly
                alpha_c         = iwe_alpha_c(:,node),          & ! out (for Ri IMIX)
                ! only for Osborn shortcut 
                ! FIXME: nils: put this to cvmix_tke
                KappaM_out      = iwe_Av(:,node),               & ! out
                KappaH_out      = iwe_Kv(:,node),               & ! out
                Nsqr            = bvfreq2(:),                   & ! in
                ! diagnostics
                iwe_Ttot        = iwe_Ttot(:,node),             &
                iwe_Tdif        = iwe_Tdif(:,node),             &
                iwe_Thdi        = iwe_Thdi(:,node),             &
                iwe_Tdis        = iwe_Tdis(:,node),             &
                iwe_Tsur        = iwe_Tsur(:,node),             &
                iwe_Tbot        = iwe_Tbot(:,node),             &
                c0              = iwe_c0(:,node),               &
                v0              = iwe_v0(:,node),               &
                ! debugging
                debug           = debug,                        &
                !i = i,                                        &
                !j = j,                                        &
                !tstep_count = tstep_count,                    &
                cvmix_int_1     = cvmix_dummy_1(:,node),        &
                cvmix_int_2     = cvmix_dummy_2(:,node),        &
                cvmix_int_3     = cvmix_dummy_3(:,node)         &
                )
            
        end do !-->do node = 1,node_size
            
        !_______________________________________________________________________
        ! --> add contribution from horizontal wave propagation
        ! Since IDEMIX is used in a global model configuration the vertical 
        ! internal wave mixing (call cvmix_coeffs_idemix) have to be extended by 
        ! a lateral diffusion term (see. Olbers D., Eden C., 2013, A Global Model 
        ! for the Diapycnal Diffusivity Induced Internal Gravity Waves...)
        !
        ! diffusion term = div_h( v_0 * idemix_tau_h * grad_h(v_0*E_iw) )
        ! 
        ! use Gaussian integral satz ... int(div vec_A)dV = ringint(A*vec_n)dA
        !                                    div vec_A    = 1/V * sum_i=1...nface( A_i*vec_n_i)*A_i
        !
        if (idemix_n_hor_iwe_prop_iter>0) then
        
            ! make boundary exchange for iwe, and iwe_v0 --> for propagation need
            ! to calculate edge contribution that crosses the halo
            call exchange_nod(iwe)
            
            !___________________________________________________________________
            ! calculate inverse volume and restrict iwe_v0 to fullfill stability 
            ! criterium --> CFL
            ! CFL Diffusion : CFL = v0^2 * dt/dx^2, CFL < 0.5
            ! --> limit v0 to CFL=0.2
            ! --> v0 = sqrt(CFL * dx^2 / dt)
            cflfac = 0.2_WP
            ! |--> FROM NILS: "fac=0.2 ist geschätzt. Würde ich erstmal so 
            !      probieren. Der kommt aus dem stabilitätskriterium für Diffusion 
            !      (ähnlich berechnet wie das CFL Kriterium nur halt für den 
            !      Diffusions anstatt für den Advektionsterm). Normalerweise 
            !      sollte der Grenzwert aber nicht zu oft auftreten. Ich hatte 
            !      mal damit rum-experimentiert, aber letztendlich war die Lösung 
            !      das Iterativ zu machen und ggf. idemix_n_hor_iwe_prop_iter zu erhöhen. 
            !      Du kannst IDEMIX erstmal ohne den Term ausprobieren und sehen, 
            !      ob es läuft, dann kannst du den dazuschalten und hoffen, dass 
            !      es nicht explodiert. Eigentlich sollte der Term alles glatter 
            !      machen, aber nahe der ML kann der schon Probleme machen".
            do node = 1,node_size
            
                ! temporarily store old iwe values for diag
                iwe_Thdi(:,node) = iwe(:,node)
                
                ! number of above bottom levels at node
                nln = nlevels_nod2D(node)-1
                uln = ulevels_nod2D(node)
                
                ! thickness of mid-level to mid-level interface at node
                dz_trr         = 0.0_WP
                dz_trr(uln+1:nln)  = Z_3d_n(uln:nln-1,node)-Z_3d_n(uln+1:nln,node)
                dz_trr(uln)      = hnode(uln,node)/2.0_WP
                dz_trr(nln+1)  = hnode(nln,node)/2.0_WP
                
                ! surface cell 
                vol_wcelli(uln,node) = 1/(areasvol(uln,node)*dz_trr(uln))
                aux = sqrt(cflfac*(area(uln,node)/pi*4.0_WP)/(idemix_tau_h*dt/idemix_n_hor_iwe_prop_iter))
                iwe_v0(uln,node) = min(iwe_v0(uln,node),aux)
                
                ! bulk cells
                !!PS do nz=2,nln
                do nz=uln+1,nln
                    ! inverse volumne 
                    vol_wcelli(nz,node) = 1/(areasvol(nz-1,node)*dz_trr(nz))
                    
                    ! restrict iwe_v0
                    aux = sqrt(cflfac*(area(nz-1,node)/pi*4.0_WP)/(idemix_tau_h*dt/idemix_n_hor_iwe_prop_iter))
                    !                  `--------+-------------´
                    !                           |-> comes from mesh_resolution=sqrt(area(1, :)/pi)*2._WP
                    iwe_v0(nz,node) = min(iwe_v0(nz,node),aux)
                end do 
                
                ! bottom cell 
                vol_wcelli(nln+1,node) = 1/(areasvol(nln,node)*dz_trr(nln+1))
                aux = sqrt(cflfac*(area(nln,node)/pi*4.0_WP)/(idemix_tau_h*dt/idemix_n_hor_iwe_prop_iter))
                iwe_v0(nln+1,node) = min(iwe_v0(nln+1,node),aux)
                
            end do !-->do node = 1,node_size
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
                ! nu1 ... upper index of ocean default = 1 but can consider cavity !=1
                nu1      = ulevels(el(1))
                
                ! thickness of mid-level to mid-level interface of element el(1)
                dz_trr         = 0.0_WP
                dz_trr(1)      = helem(1,el(1))/2.0_WP
                !!PS do nz=2,nl1-1
                do nz=nu1+1,nl1-1
                    dz_trr(nz) = helem(nz-1,el(1))/2.0_WP + helem(nz,el(1))/2.0_WP
                end do
                dz_trr(nl1)    = helem(nl1-1,el(1))/2.0_WP
                
                !_______________________________________________________________
                ! the same as above but for el(2)--> if el(2)==0 than this edge 
                ! is a boundary edge and el(2) does not exist
                nl2=0
                nu2=0
                if (el(2)>0) then 
                    deltaX2  = edge_cross_dxdy(3,edge)
                    deltaY2  = edge_cross_dxdy(4,edge)
                    elnodes2 = elem2d_nodes(:,el(2))
                    nl2      = nlevels(el(2))
                    nu2      = ulevels(el(2))
                    
                    ! thickness of mid-level to mid-level interface of element el(2)
                    dz_trr2         = 0.0_WP
                    dz_trr2(1)      = helem(1,el(2))/2.0_WP
                    !!PS do nz=2,nl2-1
                    do nz=nu2+1,nl2-1
                        dz_trr2(nz) = helem(nz-1,el(2))/2.0_WP + helem(nz,el(2))/2.0_WP
                    end do
                    dz_trr2(nl2)    = helem(nl2-1,el(2))/2.0_WP
                endif
                
                !_______________________________________________________________
                nl12=min(nl1,nl2)
                nu12=max(nu1,nu2)
                
                !_______________________________________________________________
                ! (A) goes only into this loop when the edge has only facing element
                ! el(1) --> so the edge is a boundary edge --> this is for ocean 
                ! surface in case of cavity
                do nz=nu1,nu12-1
                    !___________________________________________________________
                    ! --> calc: grad_h(v_0*E_iw)
                    ! calculate flux from el(1) with respect to edge mid 
                    ! point
                    grad_v0Eiw(1) = sum(gradient_sca(1:3,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    grad_v0Eiw(2) = sum(gradient_sca(4:6,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    dz_el         = dz_trr(nz)
                    
                    ! calculate flux 
                    vflux = (grad_v0Eiw(1)*deltaY1-grad_v0Eiw(2)*deltaX1)*dz_el
                    
                    !___________________________________________________________
                    ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                    ! multiply vflux with iwe_v0 interpolate to the edge-
                    ! mid point 
                    vflux = vflux * (iwe_v0(nz,ednodes(1))+iwe_v0(nz,ednodes(2)))*0.5_WP
                    
                    !___________________________________________________________
                    ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                    ! sum fluxes over the surface --> gaussian integral satz
                    iwe(nz,ednodes(1)) = iwe(nz,ednodes(1)) + dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(1))*vflux
                    iwe(nz,ednodes(2)) = iwe(nz,ednodes(2)) - dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(2))*vflux
                end do !-->do nz=nu1,nu12-1
                
                !_______________________________________________________________
                ! (B) goes only into this loop when the edge has only facing elemenmt
                ! el(2) --> so the edge is a boundary edge --> this is for ocean 
                ! surface in case of cavity
                if (nu2 > 0) then 
                    do nz=nu2,nu12-1
                        !___________________________________________________________
                        ! --> calc: grad_h(v_0*E_iw)
                        ! first calculate flux from el(1) with respect to edge mid 
                        ! point
                        grad_v0Eiw(1) = sum(gradient_sca(1:3,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2))
                        grad_v0Eiw(2) = sum(gradient_sca(4:6,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2))
                        dz_el         = dz_trr2(nz)
                        
                        ! calculate flux 
                        vflux = -(grad_v0Eiw(1)*deltaY2-grad_v0Eiw(2)*deltaX2)*dz_el
                        !       |--> minus sign comes from the fact that the the 
                        !            normal vectors (dx1,dy1) and (dx2,dy2) face 
                        !            in opposite direction (Right-Hand-Rule)
                        
                        !___________________________________________________________
                        ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                        ! multiply vflux with iwe_v0 interpolate to the edge-
                        ! mid point 
                        vflux = vflux * (iwe_v0(nz,ednodes(1))+iwe_v0(nz,ednodes(2)))*0.5_WP
                        
                        !___________________________________________________________
                        ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                        ! sum fluxes over the surface --> gaussian integral satz
                        iwe(nz,ednodes(1)) = iwe(nz,ednodes(1)) + dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(1))*vflux
                        iwe(nz,ednodes(2)) = iwe(nz,ednodes(2)) - dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(2))*vflux
                        
                    end do !-->do nz=nu2,nu12-1
                end if 
                !_______________________________________________________________
                ! (C) goes only into this loop when the edge has two facing elements
                ! --> so the edge is not a boundary edge
                !!PS do nz=1,nl12
                do nz=nu12,nl12
                    !___________________________________________________________
                    ! --> calc: grad_h(v_0*E_iw)
                    ! calculate grad(iwe*iwe_v0) for el(1)
                    grad_v0Eiw(1) = sum(gradient_sca(1:3,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    grad_v0Eiw(2) = sum(gradient_sca(4:6,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    
                    ! calculate grad(iwe*iwe_v0) for el(2) and average for el(1)
                    ! and el(2)
                    grad_v0Eiw(1) = (grad_v0Eiw(1) + sum(gradient_sca(1:3,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2)))*0.5_WP
                    grad_v0Eiw(2) = (grad_v0Eiw(2) + sum(gradient_sca(4:6,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2)))*0.5_WP
                    
                    ! mean mid level to midlevel tickness between el(1) and el(2)
                    dz_el         = (dz_trr(nz)+dz_trr2(nz))*0.5_WP
                    
                    ! calculate flux 
                    vflux = ((deltaX2-deltaX1)*grad_v0Eiw(2)-(deltaY2-deltaY1)*grad_v0Eiw(1))*dz_el
                    
                    !___________________________________________________________
                    ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                    ! multiply vflux with iwe_v0 interpolate to the edge-
                    ! mid point 
                    vflux = vflux * (iwe_v0(nz,ednodes(1))+iwe_v0(nz,ednodes(2)))*0.5_WP
                    
                    !___________________________________________________________
                    ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                    ! sum fluxes over the surface --> gaussian integral satz
                    iwe(nz,ednodes(1)) = iwe(nz,ednodes(1)) + dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(1))*vflux
                    iwe(nz,ednodes(2)) = iwe(nz,ednodes(2)) - dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(2))*vflux
                end do !-->do nz=1,n2
                
                !_______________________________________________________________
                ! (D) goes only into this loop when the edge has only facing element
                ! el(1) --> so the edge is a boundary edge
                do nz=nl12+1,nl1
                    !___________________________________________________________
                    ! --> calc: grad_h(v_0*E_iw)
                    ! calculate flux from el(1) with respect to edge mid 
                    ! point
                    grad_v0Eiw(1) = sum(gradient_sca(1:3,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    grad_v0Eiw(2) = sum(gradient_sca(4:6,el(1))*iwe_v0(nz,elnodes1)*iwe(nz,elnodes1))
                    dz_el         = dz_trr(nz)
                    
                    ! calculate flux 
                    vflux = (grad_v0Eiw(1)*deltaY1-grad_v0Eiw(2)*deltaX1)*dz_el
                    
                    !___________________________________________________________
                    ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                    ! multiply vflux with iwe_v0 interpolate to the edge-
                    ! mid point 
                    vflux = vflux * (iwe_v0(nz,ednodes(1))+iwe_v0(nz,ednodes(2)))*0.5_WP
                    
                    !___________________________________________________________
                    ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                    ! sum fluxes over the surface --> gaussian integral satz
                    iwe(nz,ednodes(1)) = iwe(nz,ednodes(1)) + dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(1))*vflux
                    iwe(nz,ednodes(2)) = iwe(nz,ednodes(2)) - dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(2))*vflux
                end do !-->do nz=nl12+1,nl1
                
                !_______________________________________________________________
                ! (E) goes only into this loop when the edge has only facing elemenmt
                ! el(2) --> so the edge is a boundary edge
                do nz=nl12+1,nl2
                    !___________________________________________________________
                    ! --> calc: grad_h(v_0*E_iw)
                    ! first calculate flux from el(1) with respect to edge mid 
                    ! point
                    grad_v0Eiw(1) = sum(gradient_sca(1:3,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2))
                    grad_v0Eiw(2) = sum(gradient_sca(4:6,el(2))*iwe_v0(nz,elnodes2)*iwe(nz,elnodes2))
                    dz_el         = dz_trr2(nz)
                    
                    ! calculate flux 
                    vflux = -(grad_v0Eiw(1)*deltaY2-grad_v0Eiw(2)*deltaX2)*dz_el
                    !       |--> minus sign comes from the fact that the the 
                    !            normal vectors (dx1,dy1) and (dx2,dy2) face 
                    !            in opposite direction (Right-Hand-Rule)
                    
                    !___________________________________________________________
                    ! --> calc: v_0*idemix_tau_h* grad_h(v_0*E_iw)
                    ! multiply vflux with iwe_v0 interpolate to the edge-
                    ! mid point 
                    vflux = vflux * (iwe_v0(nz,ednodes(1))+iwe_v0(nz,ednodes(2)))*0.5_WP
                    
                    !___________________________________________________________
                    ! --> calc: div(v_0*idemix_tau_h* grad_h(v_0*E_iw))
                    ! sum fluxes over the surface --> gaussian integral satz
                    iwe(nz,ednodes(1)) = iwe(nz,ednodes(1)) + dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(1))*vflux
                    iwe(nz,ednodes(2)) = iwe(nz,ednodes(2)) - dt*idemix_tau_h/idemix_n_hor_iwe_prop_iter*vol_wcelli(nz,ednodes(2))*vflux
                    
                end do !-->do nz=nl12+1,nl1
            end do !-->do edge=1,myDim_edge2D
            
            !___________________________________________________________________
            ! diagnostic: add horizontal propgation to the total production rate
            ! of internal wave energy iwe_Tot
            do node = 1,node_size
                iwe_Thdi(:,node) = ( iwe(:,node) - iwe_Thdi(:,node) )/dt
                iwe_Ttot(:,node) = iwe_Ttot(:,node) + iwe_Thdi(:,node)
            end do
        end if !-->if (idemix_n_hor_iwe_prop_iter>0) then
        
        !_______________________________________________________________________
        ! write IDEMIX diffusivities and viscositie to FESOM only when IDEMIX is 
        ! used alone --> mostly for debuging --> otherwise TKE Av and Kv are use
        if(mix_scheme_nmb==6) then 
            !___________________________________________________________________
            ! write out diffusivity
            call exchange_nod(iwe_Kv)
            Kv = iwe_Kv
            
            !___________________________________________________________________
            ! write out viscosity -->interpolate therefor from nodes to elements
            call exchange_nod(iwe_Av) !Warning: don't forget to communicate before averaging on elements!!!
            do elem=1, myDim_elem2D
                elnodes1=elem2D_nodes(:,elem)
                !!PS do nz=1,nlevels(elem)-1
                do nz=ulevels(elem),nlevels(elem)-1
                    Av(nz,elem) = sum(iwe_Av(nz,elnodes1))/3.0_WP    ! (elementwise)                
                end do
            end do
        end if 
    end subroutine calc_cvmix_idemix
end module g_cvmix_idemix
