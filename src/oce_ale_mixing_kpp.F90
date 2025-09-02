MODULE o_mixing_KPP_mod
  ! Original numerical algorithm by Bill Large at NCAR, 1994
  ! Equation numbers in the code refer to the Large etal paper. 
  !
  ! Modified from FESOM1.4 (Qiang's code) to FESOM2.0 by ozgur 
  !                                                  Oct. 2016
  ! checked by ??
  !---------------------------------------------------------------  
  USE o_PARAM
  USE MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_TRACER
  USE MOD_DYN
  USE o_ARRAYS
  USE g_config
  USE g_forcing_arrays
  USE g_comm
  USE g_support
  IMPLICIT NONE
  
  private 

  public oce_mixing_kpp_init
  public oce_mixing_kpp

  public hbl, ghats, blmc, Bo,dbsfc !!!!! ghats is nonzero only for scalars in unstable forcing conditions

  private bldepth    ! only within o_mixing_kpp_mod
  private wscale
  private ri_iwmix
  private ddmix
  private blmix_kpp
  private enhance

  real(KIND=WP), dimension(:),     allocatable  :: bfsfc    ! surface buoyancy forcing    (m^2/s^3)
  real(KIND=WP), dimension(:),     allocatable  :: caseA    ! = 1 in case A; =0 in case B
  real(KIND=WP), dimension(:),     allocatable  :: stable   ! = 1 in stable forcing; =0 in unstable
  real(KIND=WP), dimension(:,:),   allocatable  :: dkm1     ! boundary layer diff_cbt at kbl-1 level
  real(KIND=WP), dimension(:,:,:), allocatable  :: blmc     ! boundary layer mixing coefficients
  real(KIND=WP), dimension(:),     allocatable  :: ustar    ! surface friction velocity       (m/s)
  real(KIND=WP), dimension(:),     allocatable  :: Bo       ! surface turb buoy. forcing  (m^2/s^3)
  real(KIND=WP), dimension(:,:),   allocatable  :: dVsq     ! (velocity shear re sfc)^2   (m/s)^2
  real(KIND=WP), dimension(:,:),   allocatable  :: dbsfc 

  integer, dimension(:), allocatable            :: kbl      ! index of first grid level below hbl
  real(KIND=WP), dimension(:,:), allocatable    :: ghats    ! nonlocal transport (s/m^2)
  real(KIND=WP), dimension(:), allocatable      :: hbl      ! boundary layer depth

  real(KIND=WP), parameter :: epsln             = 1.0e-40_WP ! a small value

  real(KIND=WP), parameter :: epsilon_kpp       = 0.1_WP
  real(KIND=WP), parameter :: vonk              = 0.4_WP
  real(KIND=WP), parameter :: conc1             = 5.0_WP

  real(KIND=WP), parameter :: zmin            = -4.e-7_WP  ! m3/s3 limit for lookup table of wm and ws
  real(KIND=WP), parameter :: zmax            = 0.0_WP     ! m3/s3 limit for lookup table of wm and ws
  real(KIND=WP), parameter :: umin            = 0.0_WP     ! m/s limit for lookup table of wm and ws
  real(KIND=WP), parameter :: umax            = 0.04_WP    ! m/s limit for lookup table of wm and ws

  real(KIND=WP) :: cg                ! non-dimensional coefficient for counter-gradient term
  real(KIND=WP) :: Vtc              ! non-dimensional coefficient for velocity 
                                    ! scale of turbulant velocity shear        
                                    ! (=function of concv,concs,epsilon_kpp,vonk,Ricr)

  real(KIND=WP) :: deltaz           ! delta zehat in table
  real(KIND=WP) :: deltau           ! delta ustar in table

  integer, parameter                        :: nni = 890         ! number of values for zehat in the look up table
  integer, parameter                        :: nnj = 480         ! number of values for ustar in the look up table
  real(KIND=WP), dimension(0:nni+1,0:nnj+1) :: wmt ! lookup table for wm, the turbulent velocity scale for momentum
  real(KIND=WP), dimension(0:nni+1,0:nnj+1) :: wst ! lookup table for ws, the turbulent velocity scale scalars
  logical                    :: smooth_blmc      =.true.
  logical                    :: smooth_hbl       =.false.
  logical                    :: smooth_Ri_hor    =.false.
  logical                    :: smooth_Ri_ver    =.false.
  logical                    :: limit_hbl_ekmmob =.false. !.true.

  integer                                   :: n          ! to perform loop iterations 

contains

  !#######################################################################
  ! <SUBROUTINE NAME="oce_mixing_kpp_init">
  !
  ! Initialization for the KPP vertical mixing scheme
  !
  !     output:
  !       Vtc = non-dimensional constant used in calc. bulk Ri              
  !       cg  = constant used in calc.nonlocal transport term                
  !       wmt = turbulent velocity scale for momentum                         
  !       wst = turbulent velocity scale for scaler                          

!      *******************************************************************
!       Call this routine under oce_setup_step.F90 (array_setup subroutine)
!       Change the size of Kv under oce_modules.F90 (o_ARRAYS module)
!     KPP: Kv(nl,node_size,2) for T and S and Av(nl,node_size) 
!     PP:  Kv(nl,node_size) and Av(nl,elem_size)
!      *******************************************************************

  subroutine oce_mixing_kpp_init(partit, mesh)

     IMPLICIT NONE

     real(KIND=WP), parameter :: cstar  = 10.0_WP   ! proportionality coefficient for nonlocal transport
     real(KIND=WP), parameter :: conam  = 1.257_WP
     real(KIND=WP), parameter :: concm  = 8.380_WP 
     real(KIND=WP), parameter :: conc2  = 16.0_WP
     real(KIND=WP), parameter :: zetam  = -0.2_WP
     real(KIND=WP), parameter :: conas  = -28.86_WP
     real(KIND=WP), parameter :: concs  = 98.96_WP
     real(KIND=WP), parameter :: conc3  = 16.0_WP
     real(KIND=WP), parameter :: zetas  = -1.0_WP

     real(KIND=WP) :: zehat ! = zeta * ustar**3
     real(KIND=WP) :: zeta  ! = stability parameter d/L
     real(KIND=WP) :: usta

     integer :: i, j

     type(t_mesh),   intent(in),    target :: mesh
     type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

     allocate (    ghats    ( nl-1, myDim_nod2D+eDim_nod2D        ))   ! nonlocal transport (s/m^2)
     allocate (    hbl      (       myDim_nod2D+eDim_nod2D        ))   ! boundary layer depth
     allocate (    bfsfc    (       myDim_nod2D+eDim_nod2D        ))   ! surface buoyancy forcing    (m^2/s^3)
     allocate (    caseA    (       myDim_nod2D+eDim_nod2D        ))   ! = 1 in case A; =0 in case B
     allocate (    stable   (       myDim_nod2D+eDim_nod2D        ))   ! = 1 in stable forcing; =0 in unstable
     allocate (    dkm1     (       myDim_nod2D+eDim_nod2D,  3    ))   ! boundary layer diff at kbl-1 level
     allocate (    blmc     ( nl,   myDim_nod2D+eDim_nod2D,  3    ))   ! boundary layer mixing coefficients
     allocate (    ustar    (       myDim_nod2D+eDim_nod2D        ))   ! surface friction velocity       (m/s)
     allocate (    Bo       (       myDim_nod2D+eDim_nod2D        ))   ! surface turb buoy. forcing  (m^2/s^3)
     allocate (    dVsq     ( nl,   myDim_nod2D+eDim_nod2D        ))   ! (velocity shear re sfc)^2   (m/s)^2
     allocate (    dbsfc    ( nl,   myDim_nod2D+eDim_nod2D        ))   ! buoyancy re sfc
     allocate (    kbl      (       myDim_nod2D+eDim_nod2D        ))   ! index of first grid level below hbl

!$OMP PARALLEL DO
     DO n=1, myDim_nod2D+eDim_nod2D
        ghats  (:, n    ) = 0.0_WP
        hbl    (   n    ) = 0.0_WP
        bfsfc  (   n    ) = 0.0_WP
        caseA  (   n    ) = 0.0_WP
        stable (   n    ) = 0.0_WP
        dkm1   (   n,  :) = 0.0_WP
        blmc   (:, n,  :) = 0.0_WP
        ustar  (   n    ) = 0.0_WP
        Bo     (   n    ) = 0.0_WP
        dVsq   (:, n    ) = 0.0_WP
        dbsfc  (:, n    ) = 0.0_WP
        kbl    (   n    ) = 0.0_WP
    END DO
!$OMP END PARALLEL DO
!      *******************************************************************
!       Initialize some constants for kmix subroutines, and initialize
!       for kmix subroutine "wscale" the 2D-lookup table for wm and ws
!       as functions of ustar and zetahat (=vonk*sigma*hbl*bfsfc).
!      *******************************************************************
!      *******************************************************************
!       Define some non-dimensional constants  (recall epsilon_kpp=0.1)
!      *******************************************************************

!      Vtc used in eqn. 23
       Vtc = concv * sqrt(0.2_WP/concs/epsilon_kpp) / vonk**2 / Ricr
 
!      *******************************************************************
!       The nonlocal transport term is nonzero ONLY FOR SCALARS in
!       unstable (convective) forcing conditions where it has been 
!       successfully parametrized by Mailhot and Benoit [1982] 
!       as ghats = cs * surface kinematic salinity flux / (wst * hbl)
!       cs = cstar * vonk * (concs * vonk * epsilon_kpp) ** ( 1/3 )
!       cg = cs in eqn. 20
!      *******************************************************************

       cg = cstar * vonk * (concs * vonk * epsilon_kpp)**(1._WP/3._WP)

!      *******************************************************************
!       Construct the wm and ws lookup tables (eqn. 13 & B1)
!      *******************************************************************

     deltaz = (zmax-zmin)/real(nni+1,WP) 
     deltau = (umax-umin)/real(nnj+1,WP)

     do i=0,nni+1
       zehat = deltaz*(i) + zmin
       do j=0,nnj+1
          usta = deltau*(j) + umin
          zeta = zehat/(usta**3+epsln)

          if(zehat >= 0._WP) then
             wmt(i,j) = vonk*usta/(1.+conc1*zeta)
             wst(i,j) = wmt(i,j)
          else
             if(zeta > zetam) then
                wmt(i,j) = vonk* usta * (1._WP-conc2*zeta)**(1._WP/4._WP)
             else
                wmt(i,j) = vonk* (conam*usta**3-concm*zehat)**(1._WP/3._WP)
             endif
             if(zeta > zetas) then
                wst(i,j) = vonk* usta * (1._WP-conc3*zeta)**(1._WP/2._WP)
             else
                wst(i,j) = vonk* (conas*usta**3-concs*zehat)**(1._WP/3._WP)
             endif
          endif
       enddo
    enddo
  end subroutine oce_mixing_kpp_init

  !#######################################################################
  ! <SUBROUTINE NAME="oce_mixing_kpp">
  !
  ! This subroutine computes the vertical diffusivity and viscosity according
  ! to the KPP scheme of Large etal. In brief, the scheme does the 
  ! following:
  !
  ! --Compute interior mixing everywhere:                               
  !   interior mixing gets computed due to constant internal wave 
  !   background activity ("visc_cbu_iw" and "diff_cbt_iw"). ??????
  !   Mixing is enhanced in places of static instability (local Ri < 0).
  !   Additionally, mixing can be enhanced by contribution from shear 
  !   instability which is a function of the local Ri.
  !
  ! --Double diffusion:
  !   Interior mixing can be enhanced by double diffusion due to salt
  !   fingering and diffusive convection.
  !
  ! --Boundary layer:
  !
  !   (A) Boundary layer depth:
  !       at every gridpoint the depth of the oceanic boundary layer 
  !       ("hbl") gets computed by evaluating bulk richardson numbers.
  !
  !   (B) Boundary layer mixing:
  !       within the boundary layer, above hbl, vertical mixing is 
  !       determined by turbulent surface fluxes, and interior mixing at
  !       the lower boundary, i.e. at hbl.
  !
  ! outputs
  !
  !  hbl   = boundary layer depth (meters)
  !  ghats = nonlocal transport coefficient (s/m^2)
  !  viscA = viscosity coefficient  (m^2/s) 
  !  diffK = diffusion coefficient (m^2/s) 
  !
  !---------------------------------------------------------------  
  SUBROUTINE oce_mixing_KPP(viscAE, diffK, dynamics, tracers, partit, mesh)

     IMPLICIT NONE

!      *******************************************************************
!     Define allocatble arrays under oce_modules.F90
!     Allocate arrays under oce_setup_step.F90
!      *******************************************************************
     type(t_mesh),   intent(in),    target :: mesh
     type(t_partit), intent(inout), target :: partit
     type(t_tracer), intent(in),    target :: tracers
     type(t_dyn)   , intent(in),    target :: dynamics
     integer                    :: node, kn, elem, elnodes(3)
     integer                    :: nz, ns, j, q, lay, lay_mi, nzmin, nzmax
     real(KIND=WP)              :: smftu, smftv, aux, vol
     real(KIND=WP)              :: dens_up, minmix
     real(KIND=WP)              :: u_loc, v_loc
     real(kind=WP)              :: usurf, vsurf
     real(kind=WP)              :: rhopot, bulk, pz
     real(kind=WP)              :: bulk_0, bulk_pz, bulk_pz2
     real(kind=WP)              :: rho_surf, rho_insitu
     real(KIND=WP), dimension(mesh%nl, partit%myDim_elem2D+partit%eDim_elem2D), intent(inout)                     :: viscAE!for momentum (elements)
     real(KIND=WP), dimension(mesh%nl, partit%myDim_nod2D +partit%eDim_nod2D)                                     :: viscA !for momentum (nodes)
     real(KIND=WP), dimension(mesh%nl, partit%myDim_nod2D +partit%eDim_nod2D, tracers%num_tracers), intent(inout) :: diffK !for T and S
     real(kind=WP), dimension(:,:,:), pointer :: UVnode
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UVnode=>dynamics%uvnode(:,:,:)

!$OMP PARALLEL DO
     DO node=1, myDim_nod2D+eDim_nod2D
        ViscA(:, node) = 0.0_WP
     END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax, usurf, vsurf, u_loc, v_loc)
     DO node=1, myDim_nod2D !+eDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)

!      *******************************************************************
!       Eqn. 21
!       dVsq: Velocity shear referenced to surface (defined @ Z layer)
!       Compute squared vertical difference of velocity
!     (Vr-Vd)**2 
!        
!       dbsfc: Buoyancy difference between surface and grid points below 
!       i.e., Compute buoyancy difference with respect to "zref",i.e. the 
!       surface(m/s2) (defined @ zbar layer)
!       dbsfc = g * [ drho{1,k}/rho{1,k} - drho{k,k}/rho{k,k} ]
!     (Br-Bd) 
!
!       Squared buoyancy frequency which is required in this module,
!       is calculated in a separate routine
!      *******************************************************************

         ! Surface layer is our reference dVsq(m2/s2) & dbsfc(m/s2)
         dVsq (nzmin,node) = 0.0_WP  
         dbsfc(nzmin,node) = 0.0_WP 

         ! Surface velocity
         usurf = UVnode(1,nzmin,node)   
         vsurf = UVnode(2,nzmin,node)
         DO nz=nzmin+1, nzmax-1

         ! Squared velocity shear referenced to surface (@ Z)
            u_loc = 0.5_WP * ( UVnode(1,nz-1,node) + UVnode(1,nz,node) )
            v_loc = 0.5_WP * ( UVnode(2,nz-1,node) + UVnode(2,nz,node) )           
            dVsq(nz,node) = ( usurf - u_loc )**2 + ( vsurf - v_loc )**2
         ! dbsfc (buoyancy difference with respect to the surface (m/s2)) is now computed in oce_ale_pressure_bv.F90
         END DO
         dVsq ( nzmax, node ) = dVsq  ( nzmax-1, node )
     END DO
!$OMP END PARALLEL DO

!      *******************************************************************
!       compute thermal and haline expansion coefficients (without factor of rho).
!       thermal expansion coefficient without 1/rho factor       (kg/m3/C)
!             talpha= -d(rho{k,k})/d(T(k))
!       salt expansion coefficient without 1/rho factor          (kg/m3/PSU)
!             sbeta = d(rho{k,k})/d(S(k))

!       sw_alpha(nz,n) and sw_beta(nz,n) are computed @ Z where T and S are defined
!       You have to call sw_alpha_beta here if you DO NOT use Ferrari Gent McWilliams scheme 
!       Reason: oce_timestep(n) is called after subroutine oce_mixing_(K)PP 
!       where compute_sigma_xy -> sw_alpha_beta is called (Fer_GM should be set to true)
!      *******************************************************************
!       friction velocity, turbulent sfc buoyancy forcing
!       ustar = sqrt( sqrt( stress_atmoce_x^2 + stress_atmoce_y^2 ) / rho ) (m/s)
!       bo =  -g * ( Talpha*heat_flux/vcpw + Sbeta * salinity*water_flux ) (m^2/s^3)
!      *******************************************************************

  DO node=1, myDim_nod2D !+eDim_nod2D
     nzmin = ulevels_nod2D(node)
!!PS      ustar(node) = sqrt( sqrt( stress_atmoce_x(node)**2 + stress_atmoce_y(node)**2 )*density_0_r ) ! @ the surface (eqn. 2)
     ustar(node) = sqrt( sqrt( stress_node_surf(1,node)**2 + stress_node_surf(2,node)**2 )*density_0_r ) ! @ the surface (eqn. 2)
    
! Surface buoyancy forcing (eqns. A2c & A2d & A3b & A3d)
     !!PS Bo(node)  = -g * ( sw_alpha(1,node) * heat_flux(node)  / vcpw             &   !heat_flux & water_flux: positive up
     !!PS                  + sw_beta (1,node) * water_flux(node) * tr_arr(1,node,2))
     Bo(node)  = -g * ( sw_alpha(nzmin,node) * heat_flux(node)  / vcpw             &   !heat_flux & water_flux: positive up
                      + sw_beta (nzmin,node) * water_flux(node) * tracers%data(2)%values(nzmin,node)) 
  END DO
      
! compute interior mixing coefficients everywhere, due to constant 
! internal wave activity, static instability, and local shear 
! instability.
     CALL ri_iwmix(viscA, diffK, dynamics, tracers, partit, mesh)
! add double diffusion
     IF (double_diffusion) then
        CALL ddmix(diffK, tracers, partit, mesh)
     END IF

! boundary layer mixing coefficients: diagnose new b.l. depth
     CALL bldepth(partit, mesh)
! boundary layer diffusivities
     CALL blmix_kpp(viscA, diffK, partit, mesh)
! enhance diffusivity at interface kbl - 1
     CALL enhance(viscA, diffK, partit, mesh)
    
     if (smooth_blmc) then
        call exchange_nod(blmc(:,:,1), partit)
        call exchange_nod(blmc(:,:,2), partit)
        call exchange_nod(blmc(:,:,3), partit)
        do j=1, 3
          !_____________________________________________________________________  
          ! all loops go over myDim_nod2D so no halo information --> for smoothing 
          ! haloinfo is required --> therefor exchange_nod
         call smooth_nod(blmc(:,:,j), 3, partit, mesh)
        end do
     end if
!$OMP BARRIER
! then combine blmc and viscA/diffK
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax)
     DO node=1, myDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)
        DO nz=nzmin+1,nzmax-1
          IF (nz < kbl(node)) then ! within the bounday layer
             viscA(nz,node  ) = MAX(viscA(nz,node  ), blmc(nz,node,1))
             diffK(nz,node,1) = MAX(diffK(nz,node,1), blmc(nz,node,2))
             diffK(nz,node,2) = MAX(diffK(nz,node,2), blmc(nz,node,3))                    
          ELSE
             ghats(nz,node)=0.0_WP    ! outside the boundary layer set nonlocal terms to zero
          ENDIF
        END DO
     END DO    
!$OMP END PARALLEL DO
     !_____________________________________________________________________________
     ! do all node loops only over myDim_nod2D --> therefore do an halo exchange 
     ! only at the end should save some time
     call exchange_nod(diffK(:,:,1), partit)
     call exchange_nod(diffK(:,:,2), partit)
     call exchange_nod(ghats, partit)
  
     ! OVER ELEMENTS 
     call exchange_nod(viscA, partit) !Warning: don't forget to communicate before averaging on elements!!!
!$OMP BARRIER
     minmix=3.0e-3_WP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, elnodes, nz, nzmin, nzmax)
     DO elem    = 1, myDim_elem2D
        elnodes = elem2D_nodes(:,elem)
        nzmin   = ulevels(elem)
        nzmax   = nlevels(elem)
        DO nz=nzmin,nzmax-1
           viscAE(nz, elem) = SUM(viscA(nz, elnodes))/3.0_WP    ! (elementwise)                
        END DO
        viscAE( nlevels(elem), elem ) = viscAE( nlevels(elem)-1, elem)
        ! Set the mixing coeff. in the first layer above some limiting value
        ! this is very helpful to avoid huge surface velocity when vertical
        ! viscosity is very small derived from the KPP scheme.
        ! I strongly recommend this trick, at least in the current FESOM version.    
        IF (viscAE(nzmin,elem) < minmix) viscAE(nzmin,elem) = minmix    
     END DO
!$OMP END PARALLEL DO
  END SUBROUTINE oce_mixing_kpp


  !#######################################################################
  ! <SUBROUTINE NAME="bldepth">
  !
  ! The oceanic planetray boundary layer depth, hbl, is determined as
  ! the shallowest depth where the bulk richardson number is
  ! equal to the critical value, Ricr.
  !
  ! Bulk Richardson numbers are evaluated by computing velocity and
  ! buoyancy differences between values at level k and surface
  ! reference values.
  !
  ! In this configuration, the reference values are equal to the
  ! mean values between the first two levels (FESOMv2 first level
  ! @Z). When using a very fine vertical grid, these values should be 
  ! computed as the vertical average of velocity and buoyancy from 
  ! the surface down to epsilon_kpp*zcoord(k).
  !
  ! When the bulk richardson number at k exceeds Ricr, hbl is
  ! linearly interpolated between grid levels zcoord(k) and zcoord(k-1).
  !
  ! The water column and the surface forcing are diagnosed for 
  ! stable/ustable forcing conditions, and where hbl is relative 
  ! to grid points (caseA), so that conditional branches can be 
  ! avoided in later subroutines.
  !
  !
  !  input         
  !      real ustar(t2d)      = surface friction velocity          (m/s)      
  !      real dVsq(t3d)       = (velocity shear re sfc)^2          (m/s)^2
  !      real Bo(t2d)         = surface turbulent buoyancy forcing (m^2/s^3)  
  !      real sw_3d(t3d)      = radiative buoyancy forcing         (C m/s)          
  !      real f_c(t2d)        = Coriolis parameter                 (1/s)            
  !
  !  output
  !      real hbl(t2d)        ! boundary layer depth               (m)      
  !      real bfsfc(t2d)      ! Bo+radiation absorbed (m^2/s^3)      
  !      real stable(t2d)     ! =1 in stable forcing; =0 unstable          
  !      real caseA(t2d)      ! =1 in case A, =0 in case B                
  !      integer kbl(t2d)     ! index of first grid level below hbl        
  !
  SUBROUTINE bldepth(partit, mesh)

     IMPLICIT NONE

     real(KIND=WP)            :: Ritop, bvsq, Vtsq
     real(KIND=WP)            :: hekman, hmonob, hlimit
     real(KIND=WP)            :: Rib_km1, Rib_k, coeff_sw, zk, zkm1
     real(KIND=WP)            :: sigma, zehat, wm, ws, dzup, dzupE(3)
     integer                  :: node, nk, nz, elem, elnodes(3), nzmin, nzmax

     real(KIND=WP), parameter :: cekman = 0.7_WP  ! constant for Ekman depth
     real(KIND=WP), parameter :: cmonob = 1.0_WP  ! constant for Monin-Obukhov depth

     type(t_mesh),   intent(in),    target :: mesh
     type(t_partit), intent(inout), target :: partit

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

!$OMP PARALLEL DO
     ! Initialize hbl and kbl to bottomed out values
     DO node=1, myDim_nod2D !+eDim_nod2D
        ! Index of first grid level below hbl
        kbl(node) = nlevels_nod2D(node)      
        ! Boundary layer depth
        hbl(node) = ABS( zbar_3d_n( nlevels_nod2d(node),node ) )
     END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax, coeff_sw, Rib_km1, zk, zkm1, sigma, zehat, &
!$OMP                                  wm, ws, bvsq, Vtsq, Ritop, Rib_k, dzup, hekman, hmonob, hlimit)
!$OMP DO
     DO node=1, myDim_nod2D !+eDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)

        IF (use_sw_pene)  THEN                
           coeff_sw = g * sw_alpha(nzmin,node)  ! @ the surface @ Z (m/s2/K)
        END IF
 
        Rib_km1 = 0.0_WP
        bfsfc(node) = Bo(node)
        
        DO nz=nzmin+1,nzmax

           zk   = ABS( zbar_3d_n(nz,  node) )  
           zkm1 = ABS( zbar_3d_n(nz-1,node) ) 

          ! bfsfc = Bo + sw contribution
           IF (use_sw_pene)  THEN     
              bfsfc(node) = Bo(node) + coeff_sw * ( sw_3d(nzmin,node) - sw_3d(nz,node) )  ! coeff_sw [m/s2/K] sw_3d: [K m/s], positive downward
           END IF 

           stable(node) = 0.5_WP + SIGN( 0.5_WP, bfsfc(node) )
           sigma  = stable(node) + ( 1.0_WP - stable(node) ) * epsilon_kpp

          !-----------------------------------------------------------------------
          ! compute velocity scales at sigma, for z=-zbar(nz)     :
          !-----------------------------------------------------------------------

           zehat= vonk * sigma * zk * bfsfc(node)
           call wscale(zehat, ustar(node), wm, ws)    

          !-----------------------------------------------------------------------
          ! compute the turbulent shear contribution to Rib
          ! eqn. (23)
          !-----------------------------------------------------------------------        

!           IF (nz<nk) THEN
!              bvsq=( bvfreq(nz-1,node) + bvfreq(nz,node) )*0.5         ! bvfreq squared value is calculated @ zbar and bvsq @ Z
!           ELSE
              bvsq = bvfreq(nz ,node)                                   ! bvfreq squared and bvsq are calculated @ zbar
!           END IF
 
           Vtsq = zk * ws * SQRT(ABS(bvsq)) * Vtc

          !-----------------------------------------------------------------------
          !  compute bulk Richardson number at new level
          !  eqn. (21)
          !-----------------------------------------------------------------------

           Ritop = zk    *   dbsfc( nz ,node )                
           Rib_k = Ritop /  (dVsq ( nz ,node ) + Vtsq + epsln )  ! dbsfc (m/s2) and dVsq @ Z are calculated under oce_mixing_kpp
           dzup  = zk    -   zkm1

           IF (Rib_k > Ricr) THEN
             ! linearly interpolate to find hbl where Rib = Ricr
              hbl(node) = zkm1 + dzup*(Ricr-Rib_km1)/(Rib_k-Rib_km1+epsln)
              kbl(node) = nz
              EXIT
           ELSE
              Rib_km1=Rib_k
           END IF

       !-----------------------------------------------------------------------
       ! find stability and buoyancy forcing for boundary layer
       !-----------------------------------------------------------------------

           IF (use_sw_pene)  THEN     
       ! Linear interpolation of sw_3d to depth of hbl
              bfsfc(node) = Bo(node) + & 
                            coeff_sw * &
                            ( sw_3d(nzmin,node) - &
                                            ( sw_3d(nz-1,node) + &
                                                               ( sw_3d(nz,node) - sw_3d(nz-1,node) ) * ( hbl(node) - zkm1 ) / dzup &
                                            ) &
                            )
              stable(node) = 0.5_WP + SIGN( 0.5_WP, bfsfc(node) )
              bfsfc (node) = bfsfc(node) + stable(node) * epsln  ! ensures bfsfc never=0
           END IF 
        END DO
       !-----------------------------------------------------------------------
       !        check hbl limits for hekman or hmonob
       !        eqn. (24)
       !-----------------------------------------------------------------------

        !!PS IF (bfsfc(node) > 0.0_WP) THEN
        IF (bfsfc(node) > 0.0_WP .and. nzmin==1) THEN
                                          !-> no ekman or monin-obukov when there is cavity  
           hekman = cekman * ustar(node) / MAX( ABS (mesh%coriolis_node(node) ), epsln)
           hmonob = cmonob * ustar(node) * ustar(node) * ustar(node)     &
                /vonk / (bfsfc(node) + epsln) 
           hlimit = stable(node) * AMIN1( hekman, hmonob )  
           hbl(node) = AMIN1( hbl(node), hlimit )
           hbl(node) = MAX( hbl(node), ABS(zbar_3d_n(2,node)) ) 
        END IF
  END DO
!$OMP END DO
!$OMP END PARALLEL

  if (smooth_hbl) then
    call exchange_nod(hbl, partit)
    call smooth_nod(hbl, 3, partit, mesh)
  end if
!$OMP BARRIER

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax, dzup, coeff_sw)
  DO node=1, myDim_nod2D
       nzmax = nlevels_nod2D(node)
       nzmin = ulevels_nod2D(node)
       !-----------------------------------------------------------------------
       !     find new kbl 
       !-----------------------------------------------------------------------
        !!PS kbl(node) = nk
        kbl(node) = nzmax
        !!PS DO nz=2,nk
        DO nz=nzmin+1,nzmax
           IF (ABS(zbar_3d_n(nz,node)) > hbl(node)) THEN
              kbl(node) = nz
              EXIT
           END IF
        END DO

       !-----------------------------------------------------------------------
       !     find stability and buoyancy forcing for final hbl values
       !-----------------------------------------------------------------------      
        IF (use_sw_pene)  THEN     
           coeff_sw = g * sw_alpha(nzmin,node)  ! @ the surface @ Z (m/s2/K)
           ! Linear interpolation of sw_3d to depth of hbl
           bfsfc(node) = Bo(node) + & 
                         coeff_sw * &
                         ( sw_3d(nzmin,node) - &
                                         ( sw_3d(kbl(node)-1, node) + &
                                                                    ( sw_3d(kbl(node), node) - sw_3d(kbl(node)-1, node) ) &
                                                                    * ( hbl(node) + zbar_3d_n( kbl(node)-1,node) ) &
                                                                    / ( zbar_3d_n( kbl(node)-1,node) - zbar_3d_n(kbl(node),node) ) ) )
           stable(node)=0.5_WP + SIGN(0.5_WP, bfsfc(node))
           bfsfc(node) =bfsfc(node) + stable(node) * epsln 
        END IF
       !-----------------------------------------------------------------------
       !     determine caseA and caseB
       !     (if hbl is below (deeper than) the mid point of level kbl
       !     then caseA=0  else caseA=1)
       !-----------------------------------------------------------------------
        dzup         = zbar_3d_n(kbl(node)-1,node) - zbar_3d_n(kbl(node),node)
        caseA(node)  = 0.5_WP + SIGN( 0.5_WP, ABS( zbar_3d_n(kbl(node),node) ) - 0.5_WP * dzup - hbl(node) )
  END DO
!$OMP END PARALLEL DO
  END SUBROUTINE bldepth


  !#######################################################################
  ! <SUBROUTINE NAME="wscale">
  !
  ! Compute turbulent velocity scales.
  ! Use a 2D-lookup table for wm and ws as functions of ustar and
  ! zetahat (=vonk*sigma*hbl*bfsfc=zeta*ustar**3).
  !
  ! Note: the lookup table is only used for unstable conditions
  ! (zehat <= 0), in the stable domain wm (=ws) gets computed
  ! directly.
  !
  !  input
  !      real zehat    zeta*us**3
  !      real us       utar(nod), surface friction velocity    (m/s)    
  !  output                                                             
  !      real wm, ws   turbulent velocity scales 
  !
  SUBROUTINE wscale(zehat, us, wm, ws)

     IMPLICIT NONE

     real(KIND=WP), intent(in)     :: zehat, us
     real(KIND=WP), intent(out)    :: wm, ws
     real(KIND=WP)                 :: zdiff, udiff, zfrac, ufrac, fzfrac
     real(KIND=WP)                 :: wam, wbm, was, wbs, u3
     integer                       :: iz, izp1, ju, jup1

    !-----------------------------------------------------------------------
    !     use lookup table for zehat < zmax only; otherwise use
    !     stable formulae
    !-----------------------------------------------------------------------

    ! zehat = vonk * sigma * abs(depth) * bfsfc

     IF (zehat <= zmax) THEN

        zdiff = zehat-zmin
        !zdiff  = MAX(zehat-zmin, 0.)
        iz    = INT( zdiff/deltaz)
        iz    = MIN( iz , nni )
        iz    = MAX( iz , 0  )
        izp1  = iz + 1

        udiff = us-umin
        ju    = INT( MIN(udiff/deltau,real(nnj,WP)))
        ju    = MAX( ju , 0  )
        jup1  = ju+1

        zfrac = zdiff/deltaz - real(iz,WP)
        ufrac = udiff/deltau - real(ju,WP)

        fzfrac= 1._WP-zfrac
        wam   = (fzfrac)  * wmt(iz,jup1) + zfrac * wmt(izp1,jup1)
        wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac * wmt(izp1,ju  )
        wm    = (1._WP-ufrac)* wbm          + ufrac * wam

        was   = (fzfrac)  * wst(iz,jup1) + zfrac * wst(izp1,jup1)
        wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac * wst(izp1,ju  )
        ws    = (1._WP-ufrac)* wbs          + ufrac * was

     ELSE
        u3    = us*us*us
        wm    = vonk * us * u3 / ( u3 + conc1*zehat + epsln )
        ws    = wm
     ENDIF

  END SUBROUTINE wscale

  !#######################################################################
  ! <SUBROUTINE NAME="ri_iwmix">
  !
  ! Compute interior viscosity and diffusivity due 
  ! to shear instability (dependent on a local richardson number),
  ! to background internal wave activity, and 
  ! to static instability (local richardson number < 0).                    
  ! outputs:
  !    visc = viscosity coefficient (m**2/s)       
  !    diff = diffusion coefficient (m**2/s)     
  !
  subroutine ri_iwmix(viscA, diffK, dynamics, tracers, partit, mesh)
     IMPLICIT NONE
     type(t_mesh),   intent(in),    target :: mesh
     type(t_partit), intent(inout), target :: partit
     type(t_tracer), intent(in),    target :: tracers
     type(t_dyn), intent(in),    target :: dynamics
     integer                     :: node, nz, mr, nzmin, nzmax
     real(KIND=WP) , parameter   :: Riinfty = 0.8_WP                ! local Richardson Number limit for shear instability (LMD 1994 uses 0.7)
     real(KIND=WP)               :: ri_prev, tmp
     real(KIND=WP)               :: Rigg, ratio, frit
     real(KIND=WP)               :: dz_inv, shear, aux, dep, lat, Kv0_b

     real(KIND=WP), dimension(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D             ),         intent(inout) :: viscA !for momentum (nodes)
     real(KIND=WP), dimension(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D ,tracers%num_tracers), intent(inout) :: diffK !for T and S

! Put them under the namelist.oce
     logical                     :: smooth_richardson_number = .false.
     integer                     :: num_smoothings = 1              ! for vertical smoothing of Richardson number
     real(kind=WP), dimension(:,:,:), pointer :: UVnode
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
     UVnode=>dynamics%uvnode(:,:,:)

     ! Compute Richardson number and store it as diffK to save memory
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax, dz_inv, shear, mr, ri_prev, tmp)
     DO node=1, myDim_nod2D! +eDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)
        !!PS DO nz=2,nlevels_nod2d(node)-1
        DO nz=nzmin+1,nzmax-1
           dz_inv = 1.0_WP / (Z_3d_n(nz-1,node)-Z_3d_n(nz,node))  ! > 0
           shear  = ( UVnode(1, nz-1, node) - UVnode(1, nz, node) )**2 + &
                    ( UVnode(2, nz-1, node) - UVnode(2, nz, node) )**2 
           shear  = shear * dz_inv * dz_inv
           diffK(nz,node,1) = MAX( bvfreq(nz,node), 0.0_WP ) / (shear + epsln)  ! To avoid NaNs at start
        END DO                                                                  ! minimum Richardson number is 0
!      *******************************************************************
!       No need to set Richardson number for the surface and bottom layers
!       diffK @ zbar. Model do not use these levels !!!!!!!       
!      *******************************************************************
        diffK(nzmin,node,1)=diffK(nzmin+1,node,1)
        diffK(nzmax,node,1)=diffK(nzmax-1,node,1)

        ! smooth Richardson number in the vertical using a 1-2-1 filter
        IF(smooth_Ri_ver .and. nzmax > 2) then
           DO mr=1, num_smoothings
              ri_prev = 0.25_WP * diffK(1, node, 1)
              DO nz=nzmin+1, nzmax-1
                 tmp = diffK(nz,node,1)
                 diffK(nz,node,1) = ri_prev + 0.5_WP * diffK(nz,node,1) + 0.25_WP * diffK(nz+1,node,1)
                 ri_prev = 0.25_WP * tmp
              END DO
           END DO
        END IF
     END DO
!$OMP END PARALLEL DO
     IF (smooth_Ri_hor) then
       call smooth_nod(diffK(:,:,1), 3, partit, mesh)
     END IF
!$OMP BARRIER
    !___________________________________________________________________________
    ! compute viscA and diffK
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax, Rigg, ratio, frit, Kv0_b)
    do node=1, myDim_nod2D
       nzmin = ulevels_nod2D(node)
       nzmax = nlevels_nod2D(node)
       do nz=nzmin+1, nzmax-1
          !___________________________________________________________________
          ! evaluate function of Ri# for shear instability eqn. (28b&c)
          Rigg  = AMAX1( diffK(nz,node,1) , 0.0_WP)
          ratio = AMIN1( Rigg/Riinfty , 1.0_WP )
          frit  = (1.0_WP - ratio*ratio)
          frit  = frit*frit*frit
          !___________________________________________________________________
          ! viscosity
          viscA(nz,node) =  visc_sh_limit * frit + A_ver  ! A_ver= 1.e-4 Vertical harm. visc.
            
          !___________________________________________________________________
          ! diffusivity
          ! set constant background diffusivity with namelist value K_ver
          if (Kv0_const) then
             diffK(nz,node,1) = diff_sh_limit * frit + K_ver                
             ! set latitudinal and depth dependent background diffusivity after Qiangs
             ! FESOM1.4 approach
             else
             ! --> see in oce_ale_mixing_pp.F90 --> there are different 
             ! schemes of the vertical background diffusivity possible strongly
             ! depending on purpos and tuning especially with arctic focus
             call Kv0_background_qiang(Kv0_b,geo_coord_nod2D(2,node)/rad,abs(zbar_3d_n(nz,node)))
                  diffK(nz,node,1) =  diff_sh_limit * frit + Kv0_b
           end if 
           diffK(nz,node,2) = diffK(nz,node,1)  
       end do ! --> do nz=2,nlevels_nod2d(node)-1
        
        !_______________________________________________________________________
        !!! No need to set surface and bottom diffusivity. diffK @ zbar      !!!
        !!! Model do not use these levels !!!!!!!                            !!!
        viscA( nzmin, node    ) = viscA( nzmin+1, node    )
        diffK( nzmin, node, 1 ) = diffK( nzmin+1, node, 1 )
        diffK( nzmin, node, 2 ) = diffK( nzmin+1, node, 2 )

        viscA( nzmax, node    ) = viscA( nzmax-1, node    )  
        diffK( nzmax, node, 1 ) = diffK( nzmax-1, node, 1 )
        diffK( nzmax, node, 2 ) = diffK( nzmax-1, node, 2 )
        
    end do !-->do node=1, myDim_nod2D+eDim_nod2D
!$OMP END PARALLEL DO
  end subroutine ri_iwmix

 !#######################################################################
  ! <SUBROUTINE NAME="ddmix">
  !
  ! Rrho dependent interior flux parameterization.
  ! Add double-diffusion diffusivities to Ri-mix values at blending
  ! interface and below. 
  ! Salt fingering code modified july 2003
  ! by stephen.griffies@noaa.gov based on NCAR CCSM2.x
  !
  ! output: update diffu
  !
  subroutine ddmix(diffK, tracers, partit, mesh)

     IMPLICIT NONE
     type(t_mesh),   intent(in),    target :: mesh
     type(t_partit), intent(inout), target :: partit
     type(t_tracer), intent(in),    target :: tracers
     real(KIND=WP), parameter       :: Rrho0               = 1.9_WP          ! limit for double diffusive density ratio
     real(KIND=WP), parameter       :: dsfmax              = 1.e-4_WP        ! (m^2/s) max diffusivity in case of salt fingering
     real(KIND=WP), parameter       :: viscosity_molecular = 1.5e-6_WP       ! (m^2/s)

     integer                        :: node, nz, nzmin, nzmax
     real(KIND=WP)                  :: alphaDT, betaDS
     real(KIND=WP)                  :: diffdd, Rrho, prandtl

     real(KIND=WP), dimension(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D, 2), intent(inout)   :: diffK ! for T and S

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax, alphaDT, betaDS, Rrho, diffdd, prandtl)
     DO node=1, myDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)
        DO nz=nzmin+1,nzmax-1
           alphaDT = sw_alpha(nz-1,node) * tracers%data(1)%values(nz-1,node)
           betaDS  = sw_beta (nz-1,node) * tracers%data(2)%values(nz-1,node)

           IF (alphaDT > betaDS .and. betaDS > 0.0_WP) THEN

!      *******************************************************************
!       salt fingering case,  eqn. (31)
!       subtropical and tropical thermoclines where warm salty 
!       water overlies cold freshwater
!      *******************************************************************
              Rrho   = MIN(alphaDT / betaDS, Rrho0)

              diffdd = 1.0_WP -( (Rrho-1.0_WP) / (Rrho0-1.0_WP) )
              diffdd = dsfmax * diffdd * diffdd * diffdd
 
              diffK(nz,node,1) = diffK(nz,node,1) + 0.7_WP * diffdd         ! for temperature
              diffK(nz,node,2) = diffK(nz,node,2) + diffdd                  ! for salinity

           ELSE IF ( alphaDT < 0.0_WP .and. alphaDT > betaDS ) then

!      *******************************************************************
!       diffusive convection eqn. (32)
!       staircases occur primarily in the arctic and adjacent regions
!       where cold freshwater is over warm salty water
!      *******************************************************************
   
              Rrho    = alphaDT / betaDS 
              diffdd  = viscosity_molecular * 0.909_WP * exp( 4.6_WP * exp( -0.54_WP * (1.0_WP / Rrho - 1.0_WP) ) )  

            ! eqn. (34)
              prandtl = 0.15_WP * Rrho
              IF (Rrho > 0.5_WP) prandtl = (1.85_WP - 0.85_WP / Rrho) * Rrho
              diffK(nz,node,1) = diffK(nz,node,1) + diffdd               ! for temperature
              diffK(nz,node,2) = diffK(nz,node,2) + prandtl * diffdd     ! for salinity
           ENDIF
        END DO

!      *******************************************************************
!       No need to set surface and bottom diffusivity. diffK @ zbar
!       Model do not use these levels !!!!!!!       
!      *******************************************************************
        diffK( nzmin, node, 1 ) = diffK( nzmin+1, node, 1 )
        diffK( nzmin, node, 2 ) = diffK( nzmin+1, node, 2 )
        diffK( nzmax, node, 1 ) = diffK( nzmax-1, node, 1 )
        diffK( nzmax, node, 2 ) = diffK( nzmax-1, node, 2 )
     END DO
!$OMP END PARALLEL DO
  end subroutine ddmix

 !#######################################################################
  ! <SUBROUTINE NAME="blmix_kpp">
  !
  ! Mixing coefficients within boundary layer depend on surface
  ! forcing and the magnitude and gradient of interior mixing below
  ! the boundary layer ("matching").
  !
  !     inputs:
  !
  !      real ustar(2d)    ! surface friction velocity         (m/s) 
  !      real bfsfc(2d)    ! surface buoyancy forcing     (m^2/s^3)  
  !      real hbl(2d)      ! boundary layer depth              (m)   
  !      real stable(2d)   ! = 1 in stable forcing                   
  !      real caseA(2d)    ! = 1 in case A                           
  !      integer kbl(2d)   ! index of first grid level below hbl
  !
  !     outputs:
  !
  !      real dkm1(2d,3) = boundary layer diff_cbt at kbl-1 level 
  !      real blmc(3d,3) = boundary layer mixing coeff.(m**2/s)   
  !      real ghats(3d)  = nonlocal scalar transport              
  !
  subroutine blmix_kpp(viscA,diffK, partit, mesh)

     IMPLICIT NONE
     type(t_mesh),   intent(in),    target :: mesh
     type(t_partit), intent(inout), target :: partit
     integer           :: node, nz, kn, knm1, knp1, nl1, nu1
     real(KIND=WP)     :: delhat, R, dvdzup, dvdzdn
     real(KIND=WP)     :: viscp, difsp, diftp, visch, difsh, difth, f1
     real(KIND=WP)     :: sig, a1, a2, a3, Gm, Gs, Gt
     real(KIND=WP)     :: sigma, zehat, wm, ws 
     real(KIND=WP)     :: gat1m, gat1t, gat1s, dat1m, dat1s, dat1t

     real(KIND=WP)     :: dthick(mesh%nl), diff_col(mesh%nl,3)

     real(KIND=WP), dimension(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D    ), intent(inout) :: viscA ! for momentum (nodes)
     real(KIND=WP), dimension(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D, 2 ), intent(inout) :: diffK ! for T and S

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
!$OMP PARALLEL DO
     DO node=1, myDim_nod2D+eDim_nod2D
        blmc   (:, node,  :) = 0.0_WP
     END DO
!$OMP END PARALLEL DO
!    *******************************************************************
!     Kv over the NODE 
!    *******************************************************************
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, nz, kn, knm1, knp1, nl1, nu1, delhat, R, dvdzup, dvdzdn, viscp, difsp, diftp, visch, difsh, difth, f1, sig, &
!$OMP                                            a1, a2, a3, Gm, Gs, Gt, sigma, zehat, wm, ws, gat1m, gat1t, gat1s, dat1m, dat1s, dat1t, dthick, diff_col)
!$OMP DO
     DO node=1, myDim_nod2D
        nl1=nlevels_nod2d(node)
        nu1=ulevels_nod2d(node)

        if(nl1     < 3) cycle  ! a temporary solution
        if(nl1-nu1 < 2) cycle

        dthick(nu1+1:nl1-1)=0.5_WP*(hnode(nu1:nl1-2,node)+hnode(nu1+1:nl1-1,node) )
        dthick(nu1)=hnode(nu1,node)*0.5_WP
        dthick(nl1)=hnode(nl1-1,node)*0.5_WP
        
        diff_col(nu1:nl1-1,1)=viscA(nu1:nl1-1,node)
        diff_col(nu1:nl1-1,2:3)=diffK(nu1:nl1-1,node,:)
        diff_col(nl1,:)=diff_col(nl1-1,:)

!      *******************************************************************
!       Compute velocity scales at hbl. (recall epsilon_kpp=0.1)
!      *******************************************************************

        sigma = stable(node) * 1.0_WP + (1.0_WP-stable(node)) * epsilon_kpp
        zehat= vonk * sigma * hbl(node) * bfsfc(node)
        call wscale(zehat, ustar(node), wm, ws)
        
        kn = INT(caseA(node)+epsln) *(kbl(node) -1) +   &
             (1-INT(caseA(node)+epsln)) * kbl(node)
        kn   = MIN(kn,nl1-1)
        knm1 = MAX(kn-1,nu1)
        knp1 = MIN(kn+1,nl1)

!      *******************************************************************
!       Find the interior viscosities and derivatives at hbl(i) 
!       eqn. (18)
!      *******************************************************************
        delhat = ABS(Z_3d_n(kn,node))-hbl(node)
        R      = 1.0_WP - delhat / dthick(kn)

        dvdzup = (diff_col(knm1,1) - diff_col(kn,1))/dthick(kn)
        dvdzdn = (diff_col(kn,1) - diff_col(knp1,1))/dthick(knp1)
        viscp  = 0.5_WP * ( (1.0_WP - R) * (dvdzup + ABS(dvdzup))+         &
            R  * (dvdzdn + abs(dvdzdn)) )
        
        
        dvdzup = (diff_col(knm1,3) - diff_col(kn,3))/dthick(kn)
        dvdzdn = (diff_col(kn,3) - diff_col(knp1,3))/dthick(knp1)     
        difsp  = 0.5_WP * ( (1.0_WP - R) * (dvdzup + ABS(dvdzup))+         &
             R  * (dvdzdn + ABS(dvdzdn)) )

        dvdzup = (diff_col(knm1,2) - diff_col(kn,2))/dthick(kn)
        dvdzdn = (diff_col(kn,2) - diff_col(knp1,2))/dthick(knp1)
        
        diftp  = 0.5_WP * ( (1.0_WP - R) * (dvdzup + ABS(dvdzup))+         &
             R  * (dvdzdn + ABS(dvdzdn)) )

        visch  = diff_col(kn,1) + viscp * delhat
        difsh  = diff_col(kn,3) + difsp * delhat
        difth  = diff_col(kn,2) + diftp * delhat

        f1 = stable(node) * conc1 * bfsfc(node) / (ustar(node)**4+epsln)

        gat1m = visch / (hbl(node) + epsln) / (wm + epsln)
        dat1m = -viscp / (wm+epsln) + f1 * visch
        dat1m = min(dat1m, 0.0_WP) 

        gat1s = difsh  / (hbl(node) + epsln) / (ws + epsln)
        dat1s = -difsp / (ws+epsln) + f1 * difsh 
        dat1s = min(dat1s, 0.0_WP) 

        gat1t = difth /  (hbl(node) + epsln) / (ws + epsln)
        dat1t = -diftp / (ws+epsln) + f1 * difth 
        dat1t = min(dat1t, 0.0_WP) 
        
        !!PS DO nz=2,nlevels_nod2d(node)-1
        DO nz=nu1+1,nl1-1

           if (nz >= kbl(node)) exit

!      *******************************************************************
!       Compute turbulent velocity scales on the interfaces
!      *******************************************************************
           sig   = ABS(Z_3d_n(nz,node)) / (hbl(node)+epsln)
           sigma = stable(node) * sig                          &
                + (1.0_WP - stable(node)) * AMIN1(sig, epsilon_kpp)
           zehat= vonk * sigma * hbl(node) * bfsfc(node)
           call wscale(zehat, ustar(node), wm, ws)

!      *******************************************************************
!       Compute the dimensionless shape functions at the interfaces
!       eqn. (11)
!      *******************************************************************

           a1 = sig    - 2.0_WP
           a2 = 3.0_WP - 2.0_WP * sig
           a3 = sig    - 1.0_WP

           Gm = a1 + a2 * gat1m + a3 * dat1m
           Gs = a1 + a2 * gat1s + a3 * dat1s
           Gt = a1 + a2 * gat1t + a3 * dat1t

!      *******************************************************************
!       Compute boundary layer diffusivities at the interfaces
!       eqn. (10)
!      *******************************************************************

           blmc(nz,node,1) = hbl(node) * wm * sig * (1.0_WP + sig * Gm) 
           blmc(nz,node,2) = hbl(node) * ws * sig * (1.0_WP + sig * Gt) 
           blmc(nz,node,3) = hbl(node) * ws * sig * (1.0_WP + sig * Gs) 
        
!      *******************************************************************
!       Nonlocal transport term = ghats * <ws>o (eqn. 20)
!      *******************************************************************

           ghats(nz,node) = (1.0_WP - stable(node)) * cg    &
                      / (ws * hbl(node) + epsln)

        END DO

!      *******************************************************************
!       Find diffusivities at kbl-1 grid level 
!      *******************************************************************

        sig   =  ABS(zbar_3d_n(kbl(node)-1,node))  / (hbl(node) + epsln)
        sigma =  stable(node) * sig                  &
              + (1.0_WP - stable(node)) * MIN(sig, epsilon_kpp)
        zehat= vonk * sigma * hbl(node) * bfsfc(node)
        call wscale(zehat, ustar(node), wm, ws)

        a1 = sig    - 2.0_WP
        a2 = 3.0_WP - 2.0_WP * sig
        a3 = sig    - 1.0_WP

        Gm = a1 + a2 * gat1m + a3 * dat1m
        Gs = a1 + a2 * gat1s + a3 * dat1s
        Gt = a1 + a2 * gat1t + a3 * dat1t
        
        dkm1(node,1) = hbl(node) * wm * sig * (1.0_WP + sig * Gm)
        dkm1(node,2) = hbl(node) * ws * sig * (1.0_WP + sig * Gt)
        dkm1(node,3) = hbl(node) * ws * sig * (1.0_WP + sig * Gs)
     END DO
!$OMP END DO
!$OMP END PARALLEL
  end subroutine blmix_kpp

  !#######################################################################
  ! <SUBROUTINE NAME="enhance">
  !
  ! Enhance the diffusivity at the kbl-.5 interface
  !
  ! input
  !      integer kbl(n2)   =  grid above hbl                     
  !      real hbl(n2)      =  boundary layer depth (m)           
  !      real dkm1(n2,3)   =  bl diffusivity at kbl-1 grid level 
  !      real caseA(n2)    =  1 in caseA, = 0 in case B
  !
  ! input/output
  !      real ghats(ij_bounds,nk)  =  nonlocal transport     (s/m**2)
  !      modified ghats at kbl(i)-1 interface        
  ! output
  !      real blmc(n3,3) = enhanced boundary layer mixing coefficient
  !
  subroutine enhance(viscA, diffK, partit, mesh)
     IMPLICIT NONE
     type(t_mesh),   intent(in),    target :: mesh
     type(t_partit), intent(inout), target :: partit
     real(KIND=WP),  dimension(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D),   intent(inout) :: viscA !for momentum (nodes)
     real(kind=WP),  dimension(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D,2), intent(inout) :: diffK !for T and S
     integer           :: nz, node, k
     real(kind=WP)     :: delta, dkmp5, dstar

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nz, node, k, delta, dkmp5, dstar)
     DO node=1, myDim_nod2D

        k     = kbl(node) - 1
        delta = (hbl(node) + zbar_3d_n(k,node)) / (zbar_3d_n(k,node) - zbar_3d_n(k+1,node))

       ! momentum
        dkmp5 = caseA(node) * viscA(k,node)      &
              + ( 1.0_WP - caseA(node) )    * blmc( k, node, 1 )
        dstar = ( 1.0_WP - delta )**2 * dkm1( node, 1 ) + delta**2 * dkmp5
        blmc( k, node, 1 ) = (1.0_WP - delta) * viscA(k,node)  &
              + delta * dstar
        
       ! temperature:
        dkmp5 = caseA(node) * diffK(k,node,1)  &
              + ( 1.0_WP - caseA(node) ) * blmc( k, node, 2 )
        dstar = ( 1.0_WP - delta )**2 * dkm1( node, 2 ) + delta**2 * dkmp5    
        blmc( k, node, 2 ) = ( 1.0_WP - delta ) * diffK( k, node, 1)  &
             + delta * dstar
             
       ! salinity:   
        dkmp5 = caseA(node) * diffK(k,node,2)  &
              + ( 1.0_WP - caseA(node) ) * blmc( k, node, 3 )
        dstar = ( 1.0_WP - delta )**2 * dkm1( node, 3 ) + delta**2 * dkmp5
        blmc( k, node, 3 ) = ( 1.0_WP - delta ) * diffK( k, node, 2 )  &
             + delta * dstar
             
        ghats(k,node) = (1.0_WP-caseA(node)) * ghats(k,node) ! plot ghats
     END DO
!$OMP END PARALLEL DO
  end subroutine enhance
END MODULE o_mixing_KPP_mod
