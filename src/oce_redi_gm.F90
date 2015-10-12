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
  !   A function to calculate neutral slopes (elementwise)
  !
  ! INPUT:
  !   sw_alpha    = Thermal expansion coeff (alpha) [degree_C.^-1]
  !   sw_beta     = Saline contraction coeff (beta) [psu.^-1]
  !   SF          = salinity              [psu      (PSS-78)]
  !   TF          = potential temperature [degree C (ITS-90)]
  !
  ! OUTPUT:
  !   neutral_slope = Neutral_slope on control volumes
  !
  !-------------------------------------------------------------------
  use o_mesh
  use o_param
  use o_arrays
  use g_parsup
  use g_comm
  implicit none
  !
  real(kind=WP)               :: tx, ty, sx, sy, eps, tvol
  real(kind=WP)               :: tsv(2,nl,myDim_nod2D+eDIm_nod2D)
  integer                     :: n, nz, elnodes(3),elem, k,i
  real(kind=WP)  :: TF1(nl-1, myDim_nod2D+eDim_nod2D),SF1(nl-1, myDim_nod2D+eDim_nod2D)
  real(kind=WP)  :: etemp(3,nl-1, myDim_elem2D)
  !
  ! calculate alpha and beta(elementwise)
  eps=1.0e-8
  
call compute_bvfreq(TF1,SF1)
call compute_mixlay(TF1,SF1)
call calc_bv_ref


  call sw_alpha_beta(TF1,SF1)


  ! =========
  ! Vertical derivatives
  ! on full levels
  ! =========
  DO n=1,myDim_nod2D+eDim_nod2D
     DO nz=2,nlevels_nod2D(n)-1
        tvol=Z(nz-1)-Z(nz)
        tsv(1,nz,n)=(TF1(nz-1,n)-TF1(nz,n))/tvol
        tsv(2,nz,n)=(SF1(nz-1,n)-SF1(nz,n))/tvol
     END DO
        tsv(:,1,n)=0.0_WP
	nz=nlevels_nod2D(n)
	tsv(:,nz,n)=0.0_WP
   END DO

  DO n=1, myDim_nod2D
     DO nz=1, nlevels_nod2D(n)-1
        tvol=0.0_WP
	tx=0.0_WP
	ty=0.0_WP
	sx=0.0_WP
	sy=0.0_WP
	DO k=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(k,n)
           if(nlevels(elem)-1<nz) cycle
           tvol=tvol+elem_area(elem)
           tx=tx+tt_xy_stored(1,nz,elem,1)*elem_area(elem)
	   ty=ty+tt_xy_stored(2,nz,elem,1)*elem_area(elem)
	   sx=sx+tt_xy_stored(1,nz,elem,2)*elem_area(elem)
	   sy=sy+tt_xy_stored(2,nz,elem,2)*elem_area(elem)
        END DO 
    
	  tx=(-sw_alpha(nz,n)*tx+sw_beta(nz,n)*sx)/tvol
	  ty=(-sw_alpha(nz,n)*ty+sw_beta(nz,n)*sy)/tvol
	  sx=0.5_WP*(-sw_alpha(nz,n)*(tsv(1,nz,n)+tsv(1,nz+1,n))+ &
	             sw_beta(nz,n)*(tsv(2,nz,n)+tsv(2,nz+1,n)))
	  sy=sx+sign(eps,sx)
	!if(sx>=0) then 
	!sx=sx+eps
	!else
	!sx=sx-eps
	!end if
				
	neutral_slope(1,nz,n)=-tx/sy
	neutral_slope(2,nz,n)=-ty/sy
	neutral_slope(3,nz,n)=sqrt(tx*tx+ty*ty)/abs(sy)
	
     END DO
  END DO 
  call exchange_n3D(neutral_slope, 3, nl-1)
end subroutine compute_neutral_slope
!
!=================================================================
!
subroutine redi_gm_coef
USE o_ARRAYS
USE o_MESH
USE g_PARSUP
Implicit none
!variables used for Redi/GM scheme
  real(kind=WP)                 :: S(3), fcn1,fcn2, lambd, Z_mean, depth_scale
  real(kind=WP)                 :: S_cr, S_d, c_speed, factor, fc,factor1
  integer                       :: nz, n   
  real*8 :: bv_min=0.2
  real*8 :: t1
  real*8 :: redi_factor=1.
  Kd=0.

!  S_cr=1.0e-2 ;S_d=1.0e-3 !default
  S_cr=4.0e-3 ;S_d=1.0e-3 !test

  c_speed=2.0_WP
  
  ! S_cr and S_d are used to compute taper function fcn1. 
  t1=0.1
  DO n=1, myDim_nod2D+eDim_nod2D
     DO nz=1,nl-1

        !factor=sqrt(area(nz,n)/scale_area)
	! K_GM scaling
	if (sqrt(area(nz,n))>=5e4) then
		factor=1.
	elseif (sqrt(area(nz,n))>=25e3) then
		factor=t1+(1.-t1)*(sqrt(area(nz,n))*1e-3-25)/25.
        else
		factor=t1*(area(nz,n)/625e6)
        endif

	if (bvfreq(nz,n)<=0d0) then
	   factor1=1.
        else
	   factor1=min(max(bv_min,bvfreq(nz,n)/bv_ref(n)),1.)
        endif


        factor=factor*factor1*area(nz,n)/scale_area
        !factor=factor*factor1
        ! neutral slope
        S = neutral_slope(:,nz,n)   ! S(1:3): Sx,Sy and |S|
        ! prepare for tapering
        ! define fcn1 and fcn2 which are required for taper
        ! fcn1, hyperbolic tangent
        fcn1 = 0.5_WP*(1.0_WP + tanh((S_cr - S(3))/S_d))

        ! if (nz,n) is near the surface 
        ! we need tapering function fcn2, a sine function of depth.
	! It is not needed for dm95 scheme
        if (taper_scheme/='odm95') then  
           ! the limited first baroclinic Rossby radius
           ! [following Large et al(1997)] to handle singularity
           ! at the equator

	   fc=abs(coriolis_node(n)) 
	   lambd=min(max(15000._WP, c_speed/fc),100000.0_WP)
           !critical depth, above which sine tapering is necessary.
           depth_scale = lambd*S(3)
           !depth_scale = max(lambd*S(3),mixlay_dep(n))

           !if in the surface layer we add tapering function f2
           fcn2 = 1.0_WP
           if (-Z(nz) < depth_scale)  then
              fcn2 = 0.5_WP*(1.0_WP+sin(-pi*Z(nz)/depth_scale-pi/2.0_WP))
           end if
        end if
	!taper functions fcn1 and fcn2 have been calculated.
        !
        !Second, select one tapering scheme 
        !('odm95','ldd97','g2004',or just for our test '4test'):
        !
        ! (I) first optional scheme: oDM95
        ! ************************************
        ! Danabasoglu,G. and J. C.McWilliams, 1995:
        ! Sensitivity of the global ocean circulation to 
        ! parameterizations of mesoscale tracer transports
        ! Journal of Climate, 8, 2967-2987
        !
        ! For steep slope regions:
        !    Exponential taper applied to both the neutral and GM operator
        !    parts.
        !
        ! This method was similar to Gerdes 1991. In this scheme
        ! all components of the neutral mixing 
        ! tensor are rapidly scaled to zero for S>Smax, including the
	! (3,3) term. 
        !
        if (taper_scheme == 'odm95') then  
           Kd(1,nz,n)= K_hor*fcn1*factor
           Kd(2,nz,n)= K_hor*fcn1*factor
           Kd(3,nz,n)= K_GM*fcn1*factor
        end if  ! tapering_scheme=='odm95'
        !
        ! (II) second optional scheme: LDD97
        ! ************************************
        ! Large, W.G. et al. 1997
        ! Sensitivity to surface forcing and boundary layer mixing in a global 
        !  ocean model: annual-mean climatology
        ! JPO, 27, 2418-2447
        !
        ! For steep slope regions:
        !    Exponential taper applied to both neutral operator and GM operator
        ! For near surface part:
        !    Sine taper also applied to both neutral operator and GM operator. 
        ! 
        ! This method was based on GM95, adding a sin tapering. 
        ! (used in Large et al s OGCMs) Also as an optional scheme in MITgcm.
        !
        if (taper_scheme == 'ldd97') then  
           Kd(1,nz,n) = K_hor*fcn1*fcn2*factor
           Kd(2,nz,n) = Kd(1,nz,n)
           Kd(3,nz,n) = K_GM*fcn1*fcn2*factor
        end if  ! tapering_scheme=='ldd97'

        ! (III) fesom 
        ! ***********************************
        ! For steep slope region:
        ! a) no taper applied to diagonal piece of horizontal neutral operator
        ! b) hyperbolic tangent(exponential) taper applied to off-diagonal piece of
        !    horizontal operator and to diagonal and off-diagonal piece of vertical
        !    neutral diffusion operator. a)+b) means we transfer the tracer diffusion
        !    to a horizontal-vertical manner in regions of steep neutral slopes.
        ! c) Exponential taper applied to GM operator.
        ! For surface layer with small slope:
        ! a) sine taper applied to both neutral operator and GM operator, except the
        !    diagonal piece of the horizontal diffusion.
        !
        ! In one word, here we use ldd97, but always keep the diagonal
	! part of the horizontal diffusion following the suggestion of
	! Griffies (2004). 
        !
	
        if (taper_scheme == 'fesom') then
           ! diffusion part:
           Kd(1,nz,n) = K_redi*factor!*fcn1                    ! diag Hor
           Kd(4,nz,n) = Kd(1,nz,n)*fcn2                        ! diag Ver

	   Kd(2,nz,n) = K_GM*fcn1*fcn2*factor                   ! off diag
           ! skewsion part:
           Kd(3,nz,n) = K_GM*fcn1*fcn2*factor                   ! GM
        end if  ! taper_scheme=='fesom' 
   END DO
  END DO 	
 
end subroutine redi_gm_coef

subroutine compute_bvfreq(TF,SF )
use o_mesh
use o_arrays
use g_parsup
real*8, dimension(nl-1, myDim_nod2D+eDim_nod2D) :: TF,SF
    real*8    :: dz, aux, bf_prev, tmp, fx, dens_up,dens_low
    integer :: i, k, kn, mr, nodup, nodlo
    real*8    :: drhodz_min               = 1.e-10
    fx=g/density_0
    do i=1,myDim_nod2D+eDim_nod2D
       kn=nlevels_nod2D(i)-2
       do k=1, kn
          nodup=k
          nodlo=k+1
          call densityJM(TF(nodup,i),SF(nodup,i),Z(nodlo),dens_up)
          call densityJM(TF(nodlo,i),SF(nodlo,i),Z(nodlo),dens_low)
          dz=Z(nodup)-Z(nodlo)
          aux=-(dens_up-dens_low) / dz
          bvfreq(k,i)=fx*aux
       end do
       kn=nlevels_nod2D(i)-1
       bvfreq(k,i)=bvfreq(k-1,i)

!       if(smooth_bf .and. kn>2) then
!          do mr=1,num_121_smoothings
!             bf_prev=0.25*bvfreq(nod3d_below_nod2d(1,i))
!             do k=2,kn-1
!                nodup=nod3d_below_nod2d(k,i)
!                nodlo=nod3d_below_nod2d(k+1,i)
!                tmp=bvfreq(nodup)
!                bvfreq(nodup)=bf_prev+0.5*bvfreq(nodup)+0.25*bvfreq(nodlo)
!                bf_prev=0.25*tmp
!             end do
!          end do
!       end if

    end do

  end subroutine compute_bvfreq


subroutine compute_mixlay(TF,SF)
  use o_MESH
  use g_PARsup
  use o_arrays
  use g_config
  implicit none

  real*8, dimension(nl-1, myDim_nod2D+eDim_nod2D) :: TF,SF
  integer         :: m, n2, n3, k, row
  real(kind=8)    :: dens_surf, dens, md, bf, bf_up 
  real(kind=8)    :: buoyancy_crit, smallvalue

  smallvalue=1.0e-20
  buoyancy_crit=0.0003

  !mixed layer depth
  do n2=1,myDim_nod2d+eDim_nod2D
     row=n2
     !call fcn_dens0(tracer(row,1),tracer(row,2),dens_surf)
     call densityJM(TF(1,row),SF(1,row),0.,dens_surf)
     md=0.0
     bf_up=0.0
     do k=2,nlevels_nod2D(n2)-1
        call densityJM(TF(1,row),SF(1,row),Z(k),dens)!_surf)
        bf=g*(dens-dens_surf)/dens
        if(bf>=buoyancy_crit) then
           md=md+(Z(k)-md)/(bf-bf_up+smallvalue)*(buoyancy_crit-bf_up)
           exit
        else
           md=Z(k)
           bf_up=bf
        end if
     end do
     mixlay_dep(n2)=abs(md)     
  end do

end subroutine compute_mixlay

subroutine calc_bv_ref
use o_arrays
use g_parsup
use o_mesh
implicit none
real*8 :: small_number=1d-16
integer n,nz
bv_ref=0d0
do n=1,myDim_nod2D+eDim_nod2D
   do nz=1,nlevels_nod2D(n)-1
	if ((abs(Z(nz))>=mixlay_dep(n)) .and. ( bvfreq(nz,n)>=small_number )) then
           bv_ref(n)=bvfreq(nz,n)
           exit
        endif
   enddo
enddo
end subroutine calc_bv_ref
