subroutine pressure_bv
! fill in the hydrostatic pressure and the Brunt-Vaisala frequency 
! in a single pass the using split form of the equation of state
! as proposed by NR
USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_config
use i_arrays
IMPLICIT NONE

real(kind=WP)         :: dz_inv, bv,  a, rho_up, rho_dn, t, s
integer               :: node, nz, nl1
real(kind=WP)         :: rhopot(nl), bulk_0(nl), bulk_pz(nl), bulk_pz2(nl), rho(nl)
real(kind=WP)         :: bulk_up, bulk_dn, smallvalue, buoyancy_crit
 

  smallvalue=1.0e-20
  buoyancy_crit=0.0003

 ! Screen salinity
 a=0.0_8
 DO node=1, myDim_nod2D+eDim_nod2D
     DO nz=1,nlevels_nod2d(node)-1
     a=min(a,tr_arr(nz,node,2))
     enddo
 enddo
 if(a<0.) then
   write (*,*)'s<0 happens!', a
   pe_status=1
   DO node=1, myDim_nod2D+eDim_nod2D
      DO nz=1, nlevels_nod2d(node)-1
         if (tr_arr(nz, node, 2) < 0) write (*,*) 'the model blows up at n=', mylist_nod2D(node), ' ; ', 'nz=', nz
      END DO
   END DO	
 endif
 DO node=1, myDim_nod2D+eDim_nod2D
     nl1=nlevels_nod2d(node)-1
     DO nz=1, nl1
        t=tr_arr(nz, node,1)
        s=tr_arr(nz, node,2)
        call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz))
        rho(nz)= bulk_0(nz)   + Z(nz)*(bulk_pz(nz)   + Z(nz)*bulk_pz2(nz))
        rho(nz)=rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*Z(nz))-density_0
     END DO
    ! -------
    ! Pressure
    ! -------
    hpressure(1, node)=-Z(1)*rho(1)*g
    DO nz=2, nl1
	a=0.5_WP*g*(rho(nz-1)*(zbar(nz-1)-zbar(nz))+rho(nz)*(zbar(nz)-zbar(nz+1)))
	hpressure(nz, node)=hpressure(nz-1, node)+a
    END DO
    ! -------   
    ! BV frequency:  bvfreq(nl,:), squared value is stored   
    ! -------
     DO nz=2,nl1
        bulk_up = bulk_0(nz-1) + zbar(nz)*(bulk_pz(nz-1) + zbar(nz)*bulk_pz2(nz-1)) 
        bulk_dn = bulk_0(nz)   + zbar(nz)*(bulk_pz(nz)   + zbar(nz)*bulk_pz2(nz))
        rho_up = bulk_up*rhopot(nz-1) / (bulk_up + 0.1*zbar(nz))  
        rho_dn = bulk_dn*rhopot(nz)   / (bulk_dn + 0.1*zbar(nz))  
        dz_inv=1.0_WP/(Z(nz-1)-Z(nz))  
        bvfreq(nz,node)  = -g*dz_inv*(rho_up-rho_dn)/density_0
     END DO
     bvfreq(1,node)=bvfreq(2,node)
     bvfreq(nl1+1,node)=bvfreq(nl1,node) 
     ! ------
     ! The mixed layer depth 
     ! ------
     ! mixlay_depth    
     ! bv_ref
  END DO
     ! BV is defined on full levels except for the first and the last ones.
end subroutine pressure_bv
! ===========================================================================
SUBROUTINE densityJM_local(t, s, pz, rho_out)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP, only: par_ex,pe_status
IMPLICIT NONE

  !
  ! - calculates in-situ density as a function of potential temperature
  !   (relative to the surface)
  !   using the Jackett and McDougall equation of state
  !   (Copyright (c) 1992, CSIRO, Australia)
  ! - has been derived from the SPEM subroutine rhocal
  !
  !---------------------------------------------------------------------------

  real(kind=WP), intent(IN)  :: t,s,pz
  real(kind=WP), intent(OUT) :: rho_out                 
  real(kind=WP)              :: rhopot, bulk
  real(kind=WP)              :: bulk_0, bulk_pz, bulk_pz2
  !compute secant bulk modulus

  call densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot)

  bulk = bulk_0 + pz*(bulk_pz + pz*bulk_pz2) 

  rho_out = bulk*rhopot / (bulk + 0.1*pz) - density_0

end subroutine densityJM_local
		
! ===========================================================================
SUBROUTINE densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP, only: par_ex,pe_status
IMPLICIT NONE

  !
  ! - calculates in-situ density as a function of potential temperature
  !   (relative to the surface)
  !   using the Jackett and McDougall equation of state
  !   (Copyright (c) 1992, CSIRO, Australia)
  ! - has been derived from the SPEM subroutine rhocal
  !
  ! Ralph Timmermann, August 2005
  !---------------------------------------------------------------------------
  ! N. Rakowski 2014 the split form
  !---------------------------------------------------------------------------
  real(kind=WP), intent(IN)  :: t,s
  real(kind=WP), intent(OUT) :: bulk_0, bulk_pz, bulk_pz2, rhopot
  real(kind=WP)              :: s_sqrt

  real(kind=WP), parameter   :: a0    = 19092.56,     at   = 209.8925
  real(kind=WP), parameter   :: at2   = -3.041638,    at3  = -1.852732e-3
  real(kind=WP), parameter   :: at4   = -1.361629e-5
  real(kind=WP), parameter   :: as    = 104.4077,     ast  = -6.500517
  real(kind=WP), parameter   :: ast2  = .1553190,     ast3 = 2.326469e-4 
  real(kind=WP), parameter   :: ass   = -5.587545,    asst = 0.7390729 
  real(kind=WP), parameter   :: asst2 = -1.909078e-2
  real(kind=WP), parameter   :: ap    = -4.721788e-1, apt  = -1.028859e-2
  real(kind=WP), parameter   :: apt2  = 2.512549e-4,  apt3 = 5.939910e-7 
  real(kind=WP), parameter   :: aps   = 1.571896e-2,  apst = 2.598241e-4 
  real(kind=WP), parameter   :: apst2 = -7.267926e-6, apss = -2.042967e-3
  real(kind=WP), parameter   :: ap2   = 1.045941e-5,  ap2t = -5.782165e-10 
  real(kind=WP), parameter   :: ap2t2 = 1.296821e-7
  real(kind=WP), parameter   :: ap2s  = -2.595994e-7,ap2st = -1.248266e-9 
  real(kind=WP), parameter   :: ap2st2= -3.508914e-9
  
  real(kind=WP), parameter   :: b0 = 999.842594,    bt  = 6.793952e-2
  real(kind=WP), parameter   :: bt2 = -9.095290e-3, bt3 = 1.001685e-4
  real(kind=WP), parameter   :: bt4 = -1.120083e-6, bt5 = 6.536332e-9
  real(kind=WP), parameter   :: bs = 0.824493,      bst = -4.08990e-3
  real(kind=WP), parameter   :: bst2 = 7.64380e-5,  bst3 = -8.24670e-7		
  real(kind=WP), parameter   :: bst4 = 5.38750e-9	
  real(kind=WP), parameter   :: bss = -5.72466e-3,  bsst = 1.02270e-4
  real(kind=WP), parameter   :: bsst2 = -1.65460e-6,bss2 = 4.8314e-4

  !compute secant bulk modulus

  s_sqrt = sqrt(s)

  bulk_0 =  a0      + t*(at   + t*(at2  + t*(at3 + t*at4)))      &
          + s* (as  + t*(ast  + t*(ast2 + t*ast3))               &
               + s_sqrt*(ass  + t*(asst + t*asst2)))                

  bulk_pz =  ap  + t*(apt  + t*(apt2 + t*apt3))                  &
                  + s*(aps + t*(apst + t*apst2) + s_sqrt*apss)

  bulk_pz2 = ap2 + t*(ap2t + t*ap2t2)		                 &
                + s *(ap2s + t*(ap2st + t*ap2st2))

  rhopot =  b0 + t*(bt + t*(bt2 + t*(bt3  + t*(bt4  + t*bt5))))	 &
               + s*(bs + t*(bst + t*(bst2 + t*(bst3 + t*bst4)))  &
                  + s_sqrt*(bss + t*(bsst + t*bsst2))            &
                       + s* bss2)

end subroutine densityJM_components
! ===================================================================
function ptheta(s,t,p,pr)
  ! Compute local potential temperature at pr
  ! using bryden 1973 polynomial for adiabatic lapse rate
  ! and runge-kutta 4-th order integration algorithm.
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! fofonoff,n.,1977,deep-sea res.,24,489-491
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       reference prs   pr       decibars
  !       potential tmp.  theta    deg celsius
  ! checkvalue: theta= 36.89073 c,s=40 (ipss-78),t=40 deg c,
  ! p=10000 decibars,pr=0 decibars
  !
  ! Coded by ??
  ! Reviewed by ??
  !--------------------------------------------------------
  
  implicit none
  real*8 			:: ptheta, s, t, p, pr
  real*8 			:: h, xk, q
  real*8, external	        :: atg

  h = pr - p
  xk = h*atg(s,t,p)
  t = t + 0.5*xk
  q = xk
  p = p + 0.5*h
  xk = h*atg(s,t,p)
  t = t + 0.29289322*(xk-q)
  q = 0.58578644*xk + 0.121320344*q
  xk = h*atg(s,t,p)
  t = t + 1.707106781*(xk-q)
  q = 3.414213562*xk - 4.121320344*q
  p = p + 0.5*h
  xk = h*atg(s,t,p)
  ptheta = t + (xk-2.0*q)/6.0
  return
end function ptheta
!
!-----------------------------------------------------------------
!
function atg(s,t,p)
  ! adiabatic temperature gradient deg c per decibar
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       adiabatic       atg      deg. c/decibar
  ! checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
  ! t=40 deg c,p0=10000 decibars
  !
  ! Coded by ??
  ! Reviewed by ??
  !--------------------------------------------------------
  
  implicit none
  real*8  atg, s, t, p, ds

  ds = s - 35.0
  atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p   &
       +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t        &
       +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p             &
       +(-4.2393e-8*t+1.8932e-6)*ds                          &
       +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5

  return
end function atg

!======================================================================
