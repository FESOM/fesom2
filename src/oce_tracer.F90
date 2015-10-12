! Contains subroutines involving tracers
!	pressure
!	densityJM
!	pp_diff (Pakanowsky-Philander vert. diffusion)
!	tracer_gradient_nodes  (Horizontal gradients at nodes)
!	tracer_gradient_elements (Horizontal gradients at elements)
!	Solve_tracers - main routine
!	init_tracers
!	init_tracers_AB
!	adv_tracers - advection part driver
!	diff_tracers - diffusion part driver
!	diff_part_hor
!	diff_part_ver
!	ver_redi_gm
!	diff_ver_part_expl
!	diff_ver_part_impl
!	relax_to_clim !Updates tracer field and relaxes it to climatalogy if needed
!===========================================================================

SUBROUTINE pressure
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
IMPLICIT NONE

INTEGER    nz, n, nzmax
REAL(KIND=WP)   a, rho(nl)

 DO n=1,myDim_nod2D+eDim_nod2D  !! m=1,myDim_nod2D+eDim_nod2D
                                !! n=myList_nod2D(m)
 nzmax=nlevels_nod2d(n)
 ! Compute density in the vertical column
 DO nz=1,nzmax-1
 call densityJM(tr_arr(nz,n,1), tr_arr(nz,n,2), Z(nz), rho(nz))
 END DO
 hpressure(1, n)=-Z(1)*rho(1)*g

 DO nz=2, nzmax-1
	a=0.5_WP*g*(rho(nz-1)*(zbar(nz-1)-zbar(nz))+rho(nz)*(zbar(nz)-zbar(nz+1)))
	hpressure(nz,n)=hpressure(nz-1,n)+a
 END DO
 END DO

END SUBROUTINE pressure
!===========================================================================
SUBROUTINE densityJM(t, s, pz, rho_out)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP, only: par_ex,pe_status, myDim_nod2D, eDim_nod2D, mylist_nod2D
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

  real*8, intent(IN)	:: t,s,pz
  real*8, intent(OUT)	:: rho_out                 
  real*8        	:: rhopot, bulk
  integer 		:: n, nz
  !compute secant bulk modulus

  if (s.lt.0.) then
     write (*,*)'s<0 happens!',t,s
     pe_status=1
     DO n=1, myDim_nod2D+eDim_nod2D
         DO nz=1, nlevels_nod2d(n)-1
            if (tr_arr(nz, n, 2) < 0) write (*,*) 'the model blows up at n=', mylist_nod2D(n), ' ; ', 'nz=', nz
	 END DO
     END DO	
  endif

  bulk = 19092.56 + t*(209.8925 			&
       - t*(3.041638 - t*(-1.852732e-3			&
       - t*(1.361629e-5))))				&
       + s*(104.4077 - t*(6.500517			&
       -  t*(.1553190 - t*(-2.326469e-4))))		&
       + sqrt(s**3)*(-5.587545				&
       + t*(0.7390729 - t*(1.909078e-2)))		&
       - pz *(4.721788e-1 + t*(1.028859e-2		&
       + t*(-2.512549e-4 - t*(5.939910e-7))))		&
       - pz*s*(-1.571896e-2				&
       - t*(2.598241e-4 + t*(-7.267926e-6)))		&
       - pz*sqrt(s**3)					&
       *2.042967e-3 + pz*pz*(1.045941e-5		&
       - t*(5.782165e-10 - t*(1.296821e-7)))		&
       + pz*pz*s					&
       *(-2.595994e-7					&
       + t*(-1.248266e-9 + t*(-3.508914e-9)))

  rhopot = ( 999.842594				&
       + t*( 6.793952e-2		        &
       + t*(-9.095290e-3			&
       + t*( 1.001685e-4			&
       + t*(-1.120083e-6			&
       + t*( 6.536332e-9)))))			&
       + s*( 0.824493				&
       + t *(-4.08990e-3			&
       + t *( 7.64380e-5			&
       + t *(-8.24670e-7			&
       + t *( 5.38750e-9)))))			&
       + sqrt(s**3)*(-5.72466e-3		&
       + t*( 1.02270e-4			&
       + t*(-1.65460e-6)))			&
       + 4.8314e-4*s**2)
  rho_out = rhopot / (1.0 + 0.1*pz/bulk) - density_0
end subroutine densityJM
!
!----------------------------------------------------------------------------
!
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
SUBROUTINE solve_tracers
use g_parsup
use o_PARAM, only: tracer_adv,num_tracers
use o_arrays
use o_mesh
use g_comm_auto
use o_tracers
IMPLICIT NONE
integer :: tr_num
!real(kind=8)      :: t6, t2, t3, t4, t5
!1-Temperature
!2-Salinity
!>2-Rest
do tr_num=1,num_tracers
!   t2=MPI_Wtime()
select case (tracer_adv)
case(1,2) !Miura, Quadratic reconstr.
   call init_tracers(tr_num)
case(3,4,5) !MUSCL(+FCT)
   call init_tracers_AB(tr_num)
CASE DEFAULT !unknown
	if (mype==0) write(*,*) 'Unknown advection type. Check your namelists.'
        call par_ex(1)
END SELECT
!  t3=MPI_Wtime()
   call adv_tracers(tr_num)
!  t4=MPI_Wtime()
   call diff_tracers(tr_num)
!  t5=MPI_Wtime()
   call relax_to_clim(tr_num)
!  t6=MPI_Wtime()
!   if (mype==0) then
!    write(*,*) 'INIT_TRACERS',tr_num, t3-t2
!    write(*,*) 'DIFF_TRACERS',tr_num, t4-t3
!    write(*,*) 'ADV_TRACERS',tr_num, t5-t3
!    write(*,*) 'RELAX_TO_CLIM',tr_num, t6-t5
!   endif
   call exchange_nod(tr_arr(:,:,tr_num))
enddo

!where (tr_arr(:, :, 2) < 20.)
!tr_arr(:, :, 2)=20.
!end where

end subroutine solve_tracers
