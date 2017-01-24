subroutine fcn_density(t,s,z,rho)
  !
  ! - calculates insitu density as a function of potential temperature
  !   (t is relative to the surface)
  !   using the Jackett and McDougall equation of state (1992?)
  ! Qiang 02,07,2010:  Should this be updated (1995 or 2003)? The current 
  !   version is also different to the international equation of state 
  !   (Unesco 1983). What is the exact reference for this version then? 

  implicit none

  real(kind=8), intent(IN)       :: t, s, z
  real(kind=8), intent(OUT)      :: rho                 
  real(kind=8)                   :: rhopot, bulk

     bulk = 19092.56 + t*(209.8925 				&
          - t*(3.041638 - t*(-1.852732e-3			&
          - t*(1.361629e-5))))				&
          + s*(104.4077 - t*(6.500517			&
          -  t*(.1553190 - t*(-2.326469e-4))))		&
          + sqrt(s**3)*(-5.587545				&
          + t*(0.7390729 - t*(1.909078e-2)))		&
          - z *(4.721788e-1 + t*(1.028859e-2		&
          + t*(-2.512549e-4 - t*(5.939910e-7))))		&
          - z*s*(-1.571896e-2				&
          - t*(2.598241e-4 + t*(-7.267926e-6)))		&
          - z*sqrt(s**3)					&
          *2.042967e-3 + z*z*(1.045941e-5			&
          - t*(5.782165e-10 - t*(1.296821e-7)))		&
          + z*z*s						&
          *(-2.595994e-7					&
          + t*(-1.248266e-9 + t*(-3.508914e-9)))

     rhopot = ( 999.842594					&
          + t*( 6.793952e-2			&
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
     rho = rhopot / (1.0 + 0.1*z/bulk)
end subroutine fcn_density
!
!----------------------------------------------------------------------------
!
subroutine compute_density
  use o_MESH
  use o_PARAM
  use o_array
  use g_PARFE
  implicit none
  !
  integer         :: n2, n, k
  real(kind=8)    :: z
  !
  do n2=1, nod2d
     do k=1,num_layers_below_nod2d(n2)+1
        n=nod3d_below_nod2d(k,n2)
        z=min(coord_nod3D(3,n), 0.0) 
        call fcn_density(tracer(n,1), tracer(n,2), z, density_insitu(n))
     end do
  end do
end subroutine compute_density
!
!----------------------------------------------------------------------------
!
subroutine compute_pressure
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use o_MATRICES
  use g_PARFE

  implicit none
  !
  integer         :: i, n, node_lo, node_hi
  real(kind=8)    :: dens, denscor, z_up, z_lo, wd_ice_eff
  
  do n=1, nod2d
     node_hi = nod3D_below_nod2D(1,n)
     z_up=0.0 ! already re-checked, qiang, 15.03.2011
     hpressure(node_hi)=0.0
     do i=2, num_layers_below_nod2D(n)+1
        node_lo = nod3D_below_nod2D(i,n) 
        dens=0.5_8*(density_insitu(node_hi)+density_insitu(node_lo))
        z_lo=coord_nod3d(3,node_lo)
        hpressure(node_lo)=hpressure(node_hi)+(z_up-z_lo)*dens*rho0r
        node_hi=node_lo
        z_up=z_lo
     end do
  end do
end subroutine compute_pressure
