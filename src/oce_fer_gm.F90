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
subroutine compute_sigma_xy(TF1,SF1)
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
  use g_comm_auto
  implicit none
  !
  real(kind=WP), intent(IN)   :: TF1(nl-1, myDim_nod2D+eDim_nod2D), SF1(nl-1, myDim_nod2D+eDim_nod2D)
  real(kind=WP)               :: tx, ty, sx, sy, dz, vol
  integer                     :: n, nz, elnodes(3),elem, k
  !
  ! calculate alpha and beta(elementwise)
  call sw_alpha_beta(TF1, SF1) ! this is already called under oce_mixing_kpp if you are not calling KPP it should be here

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
           vol=vol+elem_area(elem)
	   elnodes=elem2D_nodes(:,elem)

           tx=tx+sum(gradient_sca(1:3,elem)*TF1(nz,elnodes))*elem_area(elem)
	   ty=ty+sum(gradient_sca(4:6,elem)*TF1(nz,elnodes))*elem_area(elem)
           sx=sx+sum(gradient_sca(1:3,elem)*SF1(nz,elnodes))*elem_area(elem)
	   sy=sy+sum(gradient_sca(4:6,elem)*SF1(nz,elnodes))*elem_area(elem)

        END DO 
        
	  sigma_xy(1,nz,n)=(-sw_alpha(nz,n)*tx+sw_beta(nz,n)*sx)/vol*density_0
	  sigma_xy(2,nz,n)=(-sw_alpha(nz,n)*ty+sw_beta(nz,n)*sy)/vol*density_0
     END DO
  END DO 

  call exchange_nod(sigma_xy)
end subroutine compute_sigma_xy
!===========================================================================
subroutine fer_solve_Gamma
   USE o_MESH
   USE o_PARAM
   USE o_ARRAYS, ONLY: sigma_xy, bvfreq, fer_gamma, fer_c, fer_K
   USE g_PARSUP
   USE g_CONFIG
   use g_comm_auto
   IMPLICIT NONE

   integer                         :: nz, n, nzmax
   real*8                          :: zinv1,zinv2, zinv, m, r
   real*8                          :: a(nl), b(nl), c(nl)
   real*8                          :: cp(nl), tp(2,nl)
   real*8, dimension(:,:), pointer :: tr
   DO n=1,myDim_nod2D
          tr=>fer_gamma(:,:,n)
!         nzmax=nlevels_nod2D(n)
          nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
          ! The first row
          c(1)=0.0_WP
          a(1)=0.0_WP
          b(1)=1.0_WP

          zinv2=1.0_WP/(zbar(1)-zbar(2))
          DO nz=2, nzmax-1
             zinv1=zinv2
	     zinv2=1.0_WP/(zbar(nz)-zbar(nz+1))
	     zinv =1.0_WP/(Z(nz-1)-Z(nz))
             a(nz)= fer_c(n)*zinv1*zinv
             c(nz)= fer_c(n)*zinv2*zinv
             b(nz)=-a(nz)-c(nz)-max(bvfreq(nz,n), 1.e-12)
          END DO
          ! The last row
          nz=nzmax
          ! The first row
          c(nz)=0.0_WP
          a(nz)=0.0_WP
          b(nz)=1.0_WP
          ! ===========================================
          ! The rhs:
          tr(:, 1)=0.
          tr(:, nzmax)=0.
          DO nz=2, nzmax-1
             r=g/density_0
             tr(1, nz)=r*0.5_WP*sum(sigma_xy(1,nz-1:nz,n))*fer_K(n)
             tr(2, nz)=r*0.5_WP*sum(sigma_xy(2,nz-1:nz,n))*fer_K(n)
          END DO 
         ! =============================================
          ! The sweep algorithm
          ! initialize c-prime and s,t-prime
          cp(1) = c(1)/b(1)
          tp(:,1) = tr(:,1)/b(1)
! solve for vectors c-prime and t, s-prime
          DO nz = 2, nzmax
           m = b(nz)-cp(nz-1)*a(nz)
           cp(nz) = c(nz)/m
           tp(:,nz) = (tr(:,nz)-tp(:,nz-1)*a(nz))/m
          END DO
! initialize x
          tr(:,nzmax) = tp(:,nzmax)
         ! solve for x from the vectors c-prime and d-prime
          do nz = nzmax-1, 1, -1
             tr(:,nz) = tp(:,nz)-cp(nz)*tr(:,nz+1)
          end do
   END DO   !!! cycle over nodes

   call exchange_nod(fer_gamma)
   END subroutine fer_solve_Gamma
!====================================================================
subroutine fer_gamma2vel
  USE o_MESH
  USE o_PARAM
  USE o_ARRAYS, ONLY: fer_gamma, fer_uv
  USE g_PARSUP
  USE g_CONFIG
  use g_comm_auto
  IMPLICIT NONE

   integer                         :: nz, el, elnod(3)
   real*8                          :: zinv
   real*8                          :: onethird=1._WP/3._WP

   DO el=1, myDim_elem2D
      elnod=elem2D_nodes(:,el)
      DO nz=1, nlevels(el)-1 !nl-1
         zinv=onethird/(zbar(nz)-zbar(nz+1))
         fer_uv(1,nz,el)=sum(fer_gamma(1,nz,elnod)-fer_gamma(1,nz+1,elnod))*zinv
         fer_uv(2,nz,el)=sum(fer_gamma(2,nz,elnod)-fer_gamma(2,nz+1,elnod))*zinv
      END DO
   END DO
   call exchange_elem(fer_uv)
!  call exchange_elem(fer_uv(1,:,:))
!  call exchange_elem(fer_uv(2,:,:))
end subroutine fer_gamma2vel
!====================================================================
subroutine fer_compute_C_K
  USE o_MESH
  USE o_PARAM
  USE o_ARRAYS, ONLY: fer_c, fer_k, bvfreq, coriolis_node
  USE g_PARSUP
  USE g_CONFIG
  use g_comm_auto
  IMPLICIT NONE

   integer                         :: n, nz, nzmax
   real*8                          :: reso, c1, rosb, scaling, rr_ratio
   real*8                          :: d, dz
   real*8                          :: x0=1.5, sigma=.15 ! Fermi function parameters to cut off GM where Rossby radius is resolved
   real*8                          :: c_min=0.5, f_min=1.e-6, r_max=200000.
   DO n=1, myDim_nod2D
      c1=0._wp
      nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
!     reso=sqrt(area(1, n)/pi)*2._WP
      reso=mesh_resolution(n)
      d=0._WP
      DO nz=1, nzmax-1
	 dz=(zbar(nz)-zbar(nz+1))
         d=d+dz
         c1=c1+dz*(sqrt(max(bvfreq(nz,n), 0._WP))+sqrt(max(bvfreq(nz+1,n), 0._WP)))/2.
      END DO
      c1=max(c_min, c1/pi) !ca. first baroclinic gravity wave speed limited from below by c_min
      rosb=min(c1/max(abs(coriolis_node(n)), f_min), r_max)
      rr_ratio=min(reso/rosb, 5._WP)
      scaling=1._WP/(1._WP+exp(-(rr_ratio-x0)/sigma))
      fer_c(n)=c1
      fer_k(n)=1500._WP*scaling
   END DO
   call exchange_nod(fer_c)
   call exchange_nod(fer_k)
end subroutine fer_compute_C_K
!====================================================================
