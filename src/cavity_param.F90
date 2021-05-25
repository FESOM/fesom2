!
!
!_______________________________________________________________________________
! compute boundary line off cavity --> boundary line is defined my nodes 
! that have at least one cavity nodes as nearest neighbour.
! Than compute for all cavity points (ulevels_nod2D>1), which is the closest
! cavity line point to that point --> use their coordinates and depth
subroutine compute_nrst_pnt2cavline(mesh)
    use MOD_MESH
    use o_PARAM , only: WP
    use o_ARRAYS, only: Z_3d_n
    use g_PARSUP
    implicit none

    type(t_mesh), intent(inout) , target :: mesh
    integer                                             :: node, kk, elnodes(3), gnode, aux_idx
    integer,       allocatable, dimension(:)            :: cavl_idx, lcl_cavl_idx
    real(kind=WP), allocatable, dimension(:)            :: cavl_lon, cavl_lat, cavl_dep,lcl_cavl_lon, lcl_cavl_lat, lcl_cavl_dep
    real(kind=WP)                                       :: aux_x, aux_y, aux_d, aux_dmin
    
#include "associate_mesh.h"  

    !___________________________________________________________________________
    if (mype==0) write(*,*) ' --> compute cavity line '
    allocate(lcl_cavl_idx(mesh%nod2d), cavl_idx(mesh%nod2d))
    allocate(lcl_cavl_lon(mesh%nod2d), cavl_lon(mesh%nod2d))
    allocate(lcl_cavl_lat(mesh%nod2d), cavl_lat(mesh%nod2d))
    allocate(lcl_cavl_dep(mesh%nod2d), cavl_dep(mesh%nod2d))
    lcl_cavl_idx = 0
    lcl_cavl_lon = 0.0
    lcl_cavl_lat = 0.0
    lcl_cavl_dep = 0.0
    cavl_idx     = 0
    cavl_lon     = 0.0
    cavl_lat     = 0.0
    cavl_dep     = 0.0
    
    do node=1,myDim_nod2d  ! should not include eDim_nod2d 
        !_______________________________________________________________________
        ! node is a cavity point
        if (ulevels_nod2D(node)>1) then 
            cycle
        !_______________________________________________________________________    
        ! node is not a cavity point --> but should have a neighbour that is cavity
        ! point
        else
            ! loop over neighbouring elements to verttice node
            do kk=1,nod_in_elem2D_num(node)
                elnodes = elem2D_nodes(:,nod_in_elem2D(kk,node))
                ! found non cavity point with a cavity point as neighbour
                if (any(ulevels_nod2D(elnodes)>1)) then
                    gnode = myList_nod2D(node)
                    lcl_cavl_idx(gnode) = 1
                    lcl_cavl_lon(gnode) = geo_coord_nod2D(1,node)
                    lcl_cavl_lat(gnode) = geo_coord_nod2D(2,node) 
                    lcl_cavl_dep(gnode) = Z_3d_n(nlevels_nod2D(node)-1,node)
                    exit
                end if 
            end do
        end if 
    end do
    
    !___________________________________________________________________________
    ! wait until all cpus are finished
    call MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
    
    ! mpi reduce from local to global 
    call MPI_AllREDUCE(lcl_cavl_idx, cavl_idx, nod2d, MPI_INTEGER         , MPI_SUM, MPI_COMM_WORLD, MPIerr)
    call MPI_AllREDUCE(lcl_cavl_lon, cavl_lon, nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
    call MPI_AllREDUCE(lcl_cavl_lat, cavl_lat, nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
    call MPI_AllREDUCE(lcl_cavl_dep, cavl_dep, nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
  
    !___________________________________________________________________________
    ! deallocate local arrays
    deallocate(lcl_cavl_idx, lcl_cavl_lon, lcl_cavl_lat, lcl_cavl_dep)

    !___________________________________________________________________________
    ! now search for all cavity points (ulevels_nod2d>1), which is the closest 
    ! cavity line points with respect to that point
    if (mype==0) write(*,*) ' --> compute nearest cavity line point with respect to cavity point'
    allocate(mesh%cavity_nrst_cavlpnt_xyz(3,myDim_nod2d+eDim_nod2d))
    mesh%cavity_nrst_cavlpnt_xyz = 0.0_WP
    do node=1,myDim_nod2d !+eDim_nod2d
        !_______________________________________________________________________
        ! --> node is not a cavity point --> cycle over it and go to next one
        if(ulevels_nod2D(node)==1) cycle 
        
        !_______________________________________________________________________
        ! do global search over all cavity line points to search for closest one
        aux_x=geo_coord_nod2D(1,node)
        aux_y=geo_coord_nod2D(2,node)
        aux_dmin = 3000.0e3_WP
        do gnode=1,nod2d
            if (cavl_idx(gnode)/=1) cycle
            call dist_on_earth(aux_x, aux_y, cavl_lon(gnode), cavl_lat(gnode), aux_d)
            if (aux_d<aux_dmin) then
                aux_idx=gnode
                aux_dmin=aux_d
            end if
        end do
        
        !_______________________________________________________________________
        ! write out coordinates and bottom depth for that closest cavity
        ! line point
        mesh%cavity_nrst_cavlpnt_xyz(1,node)=cavl_lon(aux_idx)
        mesh%cavity_nrst_cavlpnt_xyz(2,node)=cavl_lat(aux_idx)
        mesh%cavity_nrst_cavlpnt_xyz(3,node)=cavl_dep(aux_idx)
    end do
    
    !___________________________________________________________________________
    ! deallocate global arrays
    deallocate(cavl_idx, cavl_lon, cavl_lat, cavl_dep)
    
end subroutine compute_nrst_pnt2cavline
!
!
!_______________________________________________________________________________
! The three-equation model of ice-shelf ocean interaction (Hellmer et al., 1997)
! Code derived from BRIOS subroutine iceshelf (which goes back to H.Hellmer's 2D ice shelf model code)
! adjusted for use in FESOM by Ralph Timmermann, 16.02.2011
! Reviewed by ?
! adapted by P. SCholz for FESOM2.0
subroutine cavity_heat_water_fluxes_3eq(mesh)
    use MOD_MESH
    use o_PARAM , only: density_0, WP
    use o_ARRAYS, only: heat_flux, water_flux, tr_arr, Z_3d_n, Unode, density_m_rho0,density_ref
    use i_ARRAYS, only: net_heat_flux, fresh_wa_flux
    use g_PARSUP
    implicit none
    
    !___________________________________________________________________________
    type(t_mesh), intent(inout) , target :: mesh
    real (kind=WP)  :: temp,sal,tin,zice
    real (kind=WP)  :: rhow, rhor, rho
    real (kind=WP)  :: gats1, gats2, gas, gat
    real (kind=WP)  :: ep1,ep2,ep3,ep4,ep5,ep31
    real (kind=WP)  :: ex1,ex2,ex3,ex4,ex5,ex6
    real (kind=WP)  :: vt1,sr1,sr2,sf1,sf2,tf1,tf2,tf,sf,seta,re
    integer        :: node, nzmax, nzmin
    
    !___________________________________________________________________________
    real(kind=WP),parameter ::  rp =   0.                        !reference pressure
    real(kind=WP),parameter ::  a   = -0.0575                    !Foldvik&Kvinge (1974)
    real(kind=WP),parameter ::  b   =  0.0901
    real(kind=WP),parameter ::  c   =  7.61e-4

    real(kind=WP),parameter ::  pr  =  13.8                      !Prandtl number      [dimensionless]
    real(kind=WP),parameter ::  sc  =  2432.                     !Schmidt number      [dimensionless]
    real(kind=WP),parameter ::  ak  =  2.50e-3                   !dimensionless drag coeff.
    real(kind=WP),parameter ::  sak1=  sqrt(ak)
    real(kind=WP),parameter ::  un  =  1.95e-6                   !kinematic viscosity [m2/s]
    real(kind=WP),parameter ::  pr1 =  pr**(2./3.)               !Jenkins (1991)
    real(kind=WP),parameter ::  sc1 =  sc**(2./3.)

    real(kind=WP),parameter ::  tob=  -20.                       !temperatur at the ice surface
    real(kind=WP),parameter ::  rhoi=  920.                      !mean ice density
    real(kind=WP),parameter ::  cpw =  4180.0                    !Barnier et al. (1995)
    real(kind=WP),parameter ::  lhf =  3.33e+5                   !latent heat of fusion
    real(kind=WP),parameter ::  tdif=  1.54e-6                   !thermal conductivity of ice shelf !RG4190 / RG44027
    real(kind=WP),parameter ::  atk =  273.15                    !0 deg C in Kelvin
    real(kind=WP),parameter ::  cpi =  152.5+7.122*(atk+tob)     !Paterson:"The Physics of Glaciers"

    real(kind=WP),parameter ::  L    = 334000.                   ! [J/Kg]

    ! hemw = helium content of the glacial meltwater
    ! oomw = isotopic fractionation due to melting
    ! oofw = isotopic fractionation due to freezing
    !      hemw=  4.02*14.
    !      oomw= -30.
    !      oofw= -2.5
    
#include "associate_mesh.h"  

    !___________________________________________________________________________
    do node=1,myDim_nod2D !+eDim_nod2D  
        nzmin = ulevels_nod2D(node)
        if(nzmin==1) cycle ! if no cavity skip that node
        
        !_______________________________________________________________________
        temp = tr_arr(nzmin, node,1)
        sal  = tr_arr(nzmin, node,2)
        zice = Z_3d_n(nzmin, node)  !(<0)
        
        !_______________________________________________________________________
        ! Calculate the in-situ temperature tin
        !call potit(s(i,j,N,lrhs)+35.0,t(i,j,N,lrhs),-zice(i,j),rp,tin)
        call potit(sal,temp,abs(zice),rp,tin)
        
        !_______________________________________________________________________
        ! Calculate or prescribe the turbulent heat and salt transfer coeff. GAT and GAS
        ! velocity-dependent approach of Jenkins (1991)
        !rt      vt1  = 0.25*sqrt((u(i,j,N,lrhs)+u(i+1,j,N,lrhs))**2
        !rt     &                +(v(i,j,N,lrhs)+v(i,j+1,N,lrhs))**2)
        ! if(vt1.eq.0.) vt1=0.001
        !rt      re   = Hz_r(i,j,N)*ds/un        !Reynolds number
        
        vt1  = sqrt(Unode(1,nzmin,node)*Unode(1,nzmin,node)+Unode(2,nzmin,node)*Unode(2,nzmin,node))
        vt1  = max(vt1,0.001_WP)
        !vt1  = max(vt1,0.005) ! CW
        re   = 10._WP/un                   !vt1*re (=velocity times length scale over kinematic viscosity) is the Reynolds number
        
        gats1= sak1*vt1
        gats2= 2.12_WP*log(gats1*re)-9._WP
        gat  = gats1/(gats2+12.5_WP*pr1)
        gas  = gats1/(gats2+12.5_WP*sc1)
            
        !RG3417 gat  = 1.00e-4   ![m/s]  RT: to be replaced by velocity-dependent equations later
        !RG3417 gas  = 5.05e-7   ![m/s]  RT: to be replaced by velocity-dependent equations later
            
        !_______________________________________________________________________
        ! Calculate
        ! density in the boundary layer: rhow
        ! and interface pressure pg [dbar],
        ! Solve a quadratic equation for the interface salinity sb 
        ! to determine the melting/freezing rate seta.
        ! call fcn_density(temp,sal,zice,rho)
        ! rhow = rho !fcn_density returns full in-situ density now!
        ! in previous FESOM version, which has density anomaly from fcn_density, 
        ! so used density_0+rho was rhow= rho0+rho(i,j,N) in BRIOS
        rhow = density_m_rho0(nzmin,node) + density_ref(nzmin,node)
        rhor = rhoi/rhow
        
        ep1  = cpw*gat
        ep2  = cpi*gas
        ep3  = lhf*gas
        ep31 = -rhor*cpi*tdif/zice   !RG4190 / RG44027  ! CW
        ep4  = b+c*zice
        ep5  = gas/rhor
        
        ! negative heat flux term in the ice (due to -kappa/D)
        !    ex1 = a*(ep1-ep2)
        !    ex2 = ep1*(ep4-tin)+ep2*(tob+a*sal-ep4)-ep3
        !    ex3 = sal*(ep2*(ep4-tob)+ep3)
        !    ex4 = ex2/ex1
        !    ex5 = ex3/ex1
        
        !_______________________________________________________________________
        ! RT RG4190/RG44027:
        ! In case of melting ice account for changing temperature gradient, 
        ! i.e. switch from heat conduction to heat capacity approach
        tf = a*sal+ep4
        if(tin.lt.tf) then
            !freezing
            ex1 = a*(ep1+ep31)
            ex2 = ep1*(tin-ep4)+ep3+ep31*(tob-ep4)      ! heat conduction
            ex3 = ep3*sal
            ex6 = 0.5_WP
        else
            !melting
            ex1 = a*(ep1-ep2)
            ex2 = ep1*(ep4-tin)+ep2*(tob+a*sal-ep4)-ep3   ! heat capacity
            ex3 = sal*(ep2*(ep4-tob)+ep3)
            ex6 = -0.5_WP
        endif
        
        !_______________________________________________________________________
        !RT RG4190-
        ex4 = ex2/ex1
        ex5 = ex3/ex1
        
        sr1 = 0.25_WP*ex4*ex4-ex5
        !sr2 = -0.5*ex4
        sr2 = ex6*ex4   ! modified for RG4190 / RG44027 ! CW
        sf1 = sr2+sqrt(sr1)
        tf1 = a*sf1+ep4
        sf2 = sr2-sqrt(sr1)
        tf2 = a*sf2+ep4
        
        !_______________________________________________________________________
        ! Salinities < 0 psu are not defined, therefore pick the positive of the 
        ! two solutions:
        if(sf1.gt.0.) then
            tf = tf1
            sf = sf1
        else
            tf = tf2
            sf = sf2
        endif
        
        !_______________________________________________________________________
        ! Calculate the melting/freezing rate [m/s]
        ! seta = ep5*(1.0-sal/sf)     !rt thinks this is not needed
        !rt  t_surf_flux(i,j)=gat*(tf-tin)
        !rt  s_surf_flux(i,j)=gas*(sf-(s(i,j,N,lrhs)+35.0))
        
        heat_flux(node)  = rhow*cpw*gat*(tin-tf)      ! [W/m2]  ! positive for upward
        water_flux(node) =          gas*(sf-sal)/sf   ! [m/s]   !
        
        !      qo=-rhor*seta*oofw
        !      if(seta.le.0.) then
        !         qc=rhor*seta*hemw
        !         qo=rhor*seta*oomw
        !      endif
        
        ! write(*,'(a10,i10,9f10.3)') 'ice shelf',n,zice,rhow,temp,sal,tin,tf,sf,heat_flux(n),water_flux(n)*86400.*365.
        
        !for saving to output:
        net_heat_flux(node)=-heat_flux(node)   ! positive down
        fresh_wa_flux(node)=-water_flux(node)
    end do
end subroutine cavity_heat_water_fluxes_3eq
!
!
!_______________________________________________________________________________
! Compute the heat and freshwater fluxes under ice cavity using simple 2equ.
! Coded by Adriana Huerta-Casas
! Reviewed by Qiang Wang
subroutine cavity_heat_water_fluxes_2eq(mesh)
    use MOD_MESH
    use o_PARAM , only: WP
    use o_ARRAYS, only: heat_flux, water_flux, tr_arr, Z_3d_n
    use i_ARRAYS, only: net_heat_flux, fresh_wa_flux
    use g_PARSUP
    implicit none

    type(t_mesh), intent(inout) , target :: mesh
    integer        :: node, nzmin
    real(kind=WP)   :: gama, L, aux
    real(kind=WP)   :: c2, c3, c4, c5, c6
    real(kind=WP)   :: t_i, s_i, p, t_fz
    
#include "associate_mesh.h"  

    !___________________________________________________________________________
    ! parameter for computing heat and water fluxes
    gama = 1.0e-4_WP     ! heat exchange velocity [m/s]
    L    = 334000._WP    ! water to ice latent heat [J/Kg], same as used by the ice model

    ! parameter for computing freezing temperature (UNESCO 1983 equ.)
    c3 = 1.710523e-3_WP
    c4 = -2.154996e-4_WP
    c5 = -0.0575_WP
    c6 = -7.53e-4_WP
    
    !___________________________________________________________________________
    do node=1,myDim_nod2D      
        nzmin = ulevels_nod2D(node)
        if(nzmin==1) cycle
        t_i  = tr_arr(nzmin,node,1)
        s_i  = tr_arr(nzmin,node,2)
        t_fz = c3*(s_i**(3./2.)) + c4*(s_i**2) + c5*s_i + c6*abs(Z_3d_n(nzmin,node))
        
        heat_flux(node)=vcpw*gama*(t_i - t_fz)  ! Hunter2006 used cpw=3974J/Kg (*rhowat)
        water_flux(node) = -1.0*heat_flux(node)/(L*1000.0)  
        
        !for saving to output:
        net_heat_flux(node)=-heat_flux(node)
        fresh_wa_flux(node)=-water_flux(node)
    end do
end subroutine cavity_heat_water_fluxes_2eq
!
!
!_______________________________________________________________________________
! Compute the momentum fluxes under ice cavity
! Moved to this separated routine by Qiang, 20.1.2012
subroutine cavity_momentum_fluxes(mesh)
    use MOD_MESH
    use o_PARAM , only: density_0, C_d, WP
    use o_ARRAYS, only: UV, Unode, stress_surf, stress_node_surf
    use i_ARRAYS, only: u_w, v_w
    use g_PARSUP   
    implicit none
    
    !___________________________________________________________________________
    type(t_mesh), intent(inout) , target :: mesh
    integer        :: elem, elnodes(3), nzmin, node
    real(kind=WP)  :: aux

#include "associate_mesh.h" 

    !___________________________________________________________________________
    do elem=1,myDim_elem2D
        elnodes= elem2D_nodes(:,elem)
        nzmin  = ulevels(elem)
        if(nzmin==1) cycle   
        
        ! momentum stress:
        ! need to check the sensitivity to the drag coefficient
        ! here I use the bottom stress coefficient, which is 3e-3, for this FO2 work.
        aux=sqrt(UV(1,nzmin,elem)**2+UV(2,nzmin,elem)**2)*density_0*C_d 
        stress_surf(1,elem)=-aux*UV(1,nzmin,elem)
        stress_surf(2,elem)=-aux*UV(2,nzmin,elem)
    end do
    
    !___________________________________________________________________________
    do node=1,myDim_nod2D+eDim_nod2D   
        nzmin  = ulevels_nod2d(node)
        if(nzmin==1) cycle   
        
        ! momentum stress:
        ! need to check the sensitivity to the drag coefficient
        ! here I use the bottom stress coefficient, which is 3e-3, for this FO2 work.
        aux=sqrt(Unode(1,nzmin,node)**2+Unode(2,nzmin,node)**2)*density_0*C_d 
        stress_node_surf(1,node)=-aux*Unode(1,nzmin,node)
        stress_node_surf(2,node)=-aux*Unode(2,nzmin,node)
    end do
end subroutine cavity_momentum_fluxes
!
!
!_______________________________________________________________________________
subroutine cavity_ice_clean_vel(mesh)
    use MOD_MESH
    use i_ARRAYS, only: U_ice, V_ice
    use g_PARSUP   
    implicit none
    type(t_mesh), intent(inout) , target :: mesh
    integer        :: node
    
#include "associate_mesh.h" 

    do node=1,myDim_nod2d+eDim_nod2d           
        if(ulevels_nod2D(node)>1) then
            U_ice(node)=0._WP
            V_ice(node)=0._WP
        endif
    enddo
end subroutine cavity_ice_clean_vel
!
!
!_______________________________________________________________________________
subroutine cavity_ice_clean_ma(mesh)
    use MOD_MESH
    use i_ARRAYS, only: m_ice, m_snow, a_ice
    use g_PARSUP   
    implicit none
    type(t_mesh), intent(inout) , target :: mesh
    integer        :: node
    
#include "associate_mesh.h" 

    do node=1,myDim_nod2d+eDim_nod2d           
        if(ulevels_nod2D(node)>1) then
            m_ice(node) =0.0_WP
            m_snow(node)=0.0_WP
            a_ice(node) =0.0_WP
        endif
    enddo
end subroutine cavity_ice_clean_ma
!
!
!_______________________________________________________________________________
subroutine dist_on_earth(lon1, lat1, lon2, lat2, dist)
  ! distance on the earth between two points
  ! input: lon1 lat2 and lon2 lat2 in radian
  ! output: dist in m
  use o_param
  implicit none
  real(kind=WP)  :: lon1, lat1, lon2, lat2, alpha1, dist

  alpha1=acos(cos(lat1)*cos(lat2)*cos(lon1-lon2)+sin(lat1)*sin(lat2))
  dist=r_earth*abs(alpha1)
end subroutine dist_on_earth
!
!
!_______________________________________________________________________________
! Berechnet aus dem Salzgehalt[psu] (SALZ), der pot. Temperatur[oC]
! (PT) und dem Referenzdruck[dbar] (REFPRES) die in-situ Temperatur
! [oC] (TIN) bezogen auf den in-situ Druck[dbar] (PRES) mit Hilfe
! eines Iterationsverfahrens aus.
subroutine potit(salz,pt,pres,rfpres,tin)
  use o_PARAM , only: WP
    integer iter
    real(kind=WP) :: salz,pt,pres,rfpres,tin
    real(kind=WP) :: epsi, pt1,ptd,pttmpr

     real(kind=WP), parameter :: tpmd=0.001_WP

    epsi = 0.
    do iter=1,100
        tin  = pt+epsi
        pt1  = pttmpr(salz,tin,pres,rfpres)
        ptd  = pt1-pt
        if(abs(ptd).lt.tpmd) return
        epsi = epsi-ptd
    enddo
    write(6,*) ' WARNING!'
    write(6,*) ' in-situ temperature calculation has not converged.'
    stop
    return
end subroutine potit
!
!
!_______________________________________________________________________________
! Berechnet aus dem Salzgehalt/psu (SALZ), der in-situ Temperatur/degC
! (TEMP) und dem in-situ Druck/dbar (PRES) die potentielle Temperatur/
! degC (PTTMPR) bezogen auf den Referenzdruck/dbar (RFPRES). Es wird
! ein Runge-Kutta Verfahren vierter Ordnung verwendet.
! Checkwert: PTTMPR = 36.89073 DegC
!       fuer SALZ   =    40.0 psu
!            TEMP   =    40.0 DegC
!            PRES   = 10000.000 dbar
!            RFPRES =     0.000 dbar
real(kind=WP) function pttmpr(salz,temp,pres,rfpres)
  use o_PARAM , only: WP

    real(kind=WP) :: salz,temp,pres,rfpres
    real(kind=WP) :: p,t,dp,dt,q
    real(kind=WP) :: adlprt
    real(kind=WP), parameter :: ct2  =  0.29289322_WP
    real(kind=WP), parameter :: ct3  =  1.707106781_WP
    real(kind=WP), parameter :: cq2a =  0.58578644_WP
    real(kind=WP), parameter :: cq2b =  0.121320344_WP
    real(kind=WP), parameter :: cq3a =  3.414213562_WP
    real(kind=WP), parameter :: cq3b = -4.121320344_WP    
    

    p  = pres
    t  = temp
    dp = rfpres-pres
    dt = dp*adlprt(salz,t,p)
    t  = t +0.5*dt
    q = dt
    p  = p +0.5*dp
    dt = dp*adlprt(salz,t,p)
    t  = t + ct2*(dt-q)
    q  = cq2a*dt + cq2b*q
    dt = dp*adlprt(salz,t,p)
    t  = t + ct3*(dt-q)
    q  = cq3a*dt + cq3b*q
    p  = rfpres
    dt = dp*adlprt(salz,t,p)

    pttmpr = t + (dt-q-q)/6.0

end function pttmpr
!
!
!_______________________________________________________________________________
! Berechnet aus dem Salzgehalt/psu (SALZ), der in-situ Temperatur/degC
! (TEMP) und dem in-situ Druck/dbar (PRES) den adiabatischen Temperatur-
! gradienten/(K Dbar^-1) ADLPRT.
! Checkwert: ADLPRT =     3.255976E-4 K dbar^-1
!       fuer SALZ   =    40.0 psu
!            TEMP   =    40.0 DegC
!            PRES   = 10000.000 dbar
real(kind=WP) function adlprt(salz,temp,pres)
  
  use o_PARAM , only: WP
  real(kind=WP) :: salz,temp,pres, ds
  real(kind=WP), parameter :: s0 = 35.0
  real(kind=WP), parameter :: a0 =  3.5803E-5
  real(kind=WP), parameter :: a1 =  8.5258E-6
  real(kind=WP), parameter :: a2 = -6.8360E-8
  real(kind=WP), parameter :: a3 =  6.6228E-10
  real(kind=WP), parameter :: b0 =  1.8932E-6
  real(kind=WP), parameter :: b1 = -4.2393E-8
  real(kind=WP), parameter :: c0 =  1.8741E-8
  real(kind=WP), parameter :: c1 = -6.7795E-10
  real(kind=WP), parameter :: c2 =  8.7330E-12
  real(kind=WP), parameter :: c3 = -5.4481E-14
  real(kind=WP), parameter :: d0 = -1.1351E-10
  real(kind=WP), parameter :: d1 =  2.7759E-12
  real(kind=WP), parameter :: e0 = -4.6206E-13
  real(kind=WP), parameter :: e1 =  1.8676E-14
  real(kind=WP), parameter :: e2 = -2.1687E-16

  ds = salz-s0
  adlprt = ( ( (e2*temp + e1)*temp + e0 )*pres                     &
       + ( (d1*temp + d0)*ds                                  &
       + ( (c3*temp + c2)*temp + c1 )*temp + c0 ) )*pres   &
       + (b1*temp + b0)*ds +  ( (a3*temp + a2)*temp + a1 )*temp + a0
end function adlprt
