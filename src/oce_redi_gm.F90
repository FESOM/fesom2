! A set of routines to implement Redi and GM diffusivity tensors.
! Tapering and expansion/contraction is adapted from FESOM as 
! implemented by Q. Wang.
! sergey.danilov@awi.de 2012
!Contains:
!	sw_alpha_beta
!	compute_neutral_slope
!	redi_gm_coef 
!----------------------------------------------------------------------------
!
subroutine sw_alpha_beta(TF1,SF1)
  ! DESCRIPTION:
  !   A function to calculate the thermal expansion coefficient
  !   and saline contraction coefficient. (elementwise)
  !
  ! INPUT:
  !   tracer(:,2) = salinity              [psu      (PSS-78)]
  !   tracer(:,1) = potential temperature [degree C (ITS-90)]
  !   z           = pressure (or -depth)  [db]
  !
  ! OUTPUT:
  !   sw_alpha = Thermal expansion coeff (alpha) [degree_C.^-1]
  !   sw_beta  = Saline contraction coeff (beta) [psu.^-1]
  !
  ! Qiang Wang, 25,11,2004
  !
  ! REFERENCE:
  !    McDougall, T.J. 1987.  Neutral Surfaces
  !    Journal of Physical Oceanography, vol 17, 1950-1964,
  !-----------------------------------------------------------------
  ! CHECK VALUE:
  !    sw_beta=0.72088e-3 psu^-1 @ S=40.0psu, ptmp=10.0C (ITS-90), p=4000db
  !    a_over_b=0.34765 psu*C^-1 @ S=40.0psu, ptmp=10.0C, p=4000db
  !-----------------------------------------------------------------
  use o_mesh
  use o_arrays
  use g_parsup
  use o_param
  implicit none
  !
  integer        :: n, nz
  real(kind=WP)  :: t1,t1_2,t1_3,t1_4,p1,p1_2,p1_3,s1,s35,s35_2 
  real(kind=WP)  :: a_over_b    
  real(kind=WP)  :: TF1(nl-1, myDim_nod2D+eDim_nod2D),SF1(nl-1, myDim_nod2D+eDim_nod2D)

  do n = 1,myDim_nod2d
     do nz=1, nlevels_nod2D(n)-1
     
     t1 = TF1(nz,n)*1.00024_WP
     s1 = SF1(nz,n)
     p1 = abs(Z(nz)) 
     
     t1_2 = t1*t1
     t1_3 = t1_2*t1
     t1_4 = t1_3*t1
     p1_2 = p1*p1
     p1_3 = p1_2*p1
     s35 = s1-35.0_8
     s35_2 = s35*s35

     ! calculate beta
     sw_beta(nz,n) = 0.785567e-3 - 0.301985e-5*t1 &
          + 0.555579e-7*t1_2 - 0.415613e-9*t1_3 &
          + s35*(-0.356603e-6 + 0.788212e-8*t1 &
          + 0.408195e-10*p1 - 0.602281e-15*p1_2) &
          + s35_2*(0.515032e-8) & 
          + p1*(-0.121555e-7 + 0.192867e-9*t1 - 0.213127e-11*t1_2) &
          + p1_2*(0.176621e-12 - 0.175379e-14*t1) &
          + p1_3*(0.121551e-17)

     ! calculate the thermal expansion / saline contraction ratio
     a_over_b = 0.665157e-1 + 0.170907e-1*t1 &
          - 0.203814e-3*t1_2 + 0.298357e-5*t1_3 &
          - 0.255019e-7*t1_4 &
          + s35*(0.378110e-2 - 0.846960e-4*t1 &
          - 0.164759e-6*p1 - 0.251520e-11*p1_2) &
          + s35_2*(-0.678662e-5) &
          + p1*(0.380374e-4 - 0.933746e-6*t1 + 0.791325e-8*t1_2) &
          + p1_2*t1_2*(0.512857e-12) &
          - p1_3*(0.302285e-13)

     ! calculate alpha
     sw_alpha(nz,n) = a_over_b*sw_beta(nz,n)
   end do
 end do     
end subroutine sw_alpha_beta
!
!====================================================================
!
subroutine compute_neutral_slope(TF1,SF1)
  !--------------------------------------------------------------------
  ! DESCRIPTION:
  !   A function to calculate GM_rhs on nodes
  !
  ! INPUT:
  !   SF          = salinity              [psu      (PSS-78)]
  !   TF          = potential temperature [degree C (ITS-90)]
  !
  ! OUTPUT:
  ! based on thermal expansion and saline contraction coefficients computes:
  !   density horizontal gradient on nodes scaled with g/rho_0 (GM_rhs)
  !   squared boyancy (GM_N2)
  !-------------------------------------------------------------------
  use o_mesh
  use o_param
  use o_arrays
  use g_parsup
  use g_comm
  implicit none
  !
  real(kind=WP), intent(IN)   :: TF1(nl-1, myDim_nod2D+eDim_nod2D), SF1(nl-1, myDim_nod2D+eDim_nod2D)
  real(kind=WP)               :: tx, ty, sx, sy, eps, dz, vol
  real(kind=WP)               :: sigma_x, sigma_y, sigma_z
  integer                     :: n, nz, elnodes(3),elem, k
  !
  ! calculate alpha and beta(elementwise)
  eps=1.0e-8
  
  call sw_alpha_beta(TF1, SF1)

  ! =========
  ! Vertical derivatives
  ! on full levels
  ! =========
  DO n=1,myDim_nod2D+eDim_nod2D
     DO nz=2,nlevels_nod2D(n)-1
        dz=Z(nz-1)-Z(nz)
        tz(nz,n)=(TF1(nz-1,n)-TF1(nz,n))/dz
        sz(nz,n)=(SF1(nz-1,n)-SF1(nz,n))/dz
     END DO
        tz(1,n)=0.0_WP
	tz(1,n)=0.0_WP
	nz=nlevels_nod2D(n)
        tz(nz,n)=0.0_WP
        sz(nz,n)=0.0_WP
   END DO

  DO n=1, myDim_nod2D
     DO nz=1, nlevels_nod2D(n)-1
        vol=0.0_WP
	tx=0.0_WP
	ty=0.0_WP
	sx=0.0_WP
	sy=0.0_WP
	DO k=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(k, n)
           if(nlevels(elem)-1 < nz) cycle
           tvol=tvol+elem_area(elem)
	   elnodes=elem2D_nodes(:,elem)

           tx=tx+sum(gradient_sca(1:3,elem)*TF1(nz,elnodes))*elem_area(elem)
	   ty=ty+sum(gradient_sca(4:6,elem)*TF1(nz,elnodes))*elem_area(elem)
           sx=sx+sum(gradient_sca(1:3,elem)*SF1(nz,elnodes))*elem_area(elem)
	   sy=sy+sum(gradient_sca(4:6,elem)*SF1(nz,elnodes))*elem_area(elem)
        END DO 
    
	  sigma_xyz(1,nz,n)=(-sw_alpha(nz,n)*tx+sw_beta(nz,n)*sx)/tvol
	  sigma_xyz(2,nz,n)=(-sw_alpha(nz,n)*ty+sw_beta(nz,n)*sy)/tvol
	  sigma_xyz(3,nz,n)=0.5_WP*(-sw_alpha(nz,n)*(tz(nz,n)+tz(nz+1,n))+ &
	             sw_beta(nz,n)*(sz(nz,n)+sz(nz+1,n)))
	  sigma_xyz(3,nz,n)=sigma_xyz(3,nz,n)+sign(eps, sigma_z)
     END DO
  END DO 
  call exchange_nod(neutral_slope)
end subroutine compute_neutral_slope
