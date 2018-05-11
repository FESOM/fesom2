! ====================================================================
! New evp implementation following Bouillion et al. 2013
! and Kimmritz et al. 2015 (mEVP) and Kimmritz et al. 2016 (aEVP)
! ====================================================================  
subroutine stress_tensor_m
  ! Internal stress tensor
  ! New implementation following Boullion et al, Ocean Modelling 2013.
  ! SD, 30.07.2014
  !===================================================================
  use o_param
  use i_param
  use o_mesh
  use g_config
  use i_arrays
  use g_parsup
  implicit none

  integer        :: elem, elnodes(3)
  real(kind=8)   :: dx(3), dy(3), msum, asum
  real(kind=8)   :: eps11, eps12, eps22, eps1, eps2, pressure, delta
  real(kind=8)   :: val3, meancos, usum, vsum, vale
  real(kind=8)   :: det1, det2, r1, r2, r3, si1, si2
  
  val3=1.0_8/3.0_8
  vale=1.0_8/(ellipse**2)
  det2=1.0_8/(1.0_8+alpha_evp)
  det1=alpha_evp*det2
   do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)

     msum=sum(m_ice(elnodes))*val3
     if(msum<=0.01) cycle !DS
     asum=sum(a_ice(elnodes))*val3
     
     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)     
      ! METRICS:
          vsum=sum(v_ice_aux(elnodes))
          usum=sum(u_ice_aux(elnodes))
          meancos=metric_factor(elem)
      !  
      ! ====== Deformation rate tensor on element elem:
     eps11=sum(dx*u_ice_aux(elnodes))
     eps11=eps11-val3*vsum*meancos                !metrics
     eps22=sum(dy*v_ice_aux(elnodes))
     eps12=0.5_8*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes))
     eps12=eps12+0.5_8*val3*usum*meancos          !metrics 
     
      ! ======= Switch to eps1,eps2
     eps1=eps11+eps22
     eps2=eps11-eps22   
     
      ! ====== moduli:
     delta=eps1**2+vale*(eps2**2+4.0_8*eps12**2)
     delta=sqrt(delta)
    
     pressure=pstar*msum*exp(-c_pressure*(1.0_8-asum))/(delta+delta_min)
    
        r1=pressure*(eps1-delta) 
        r2=pressure*eps2*vale
        r3=pressure*eps12*vale
        si1=sigma11(elem)+sigma22(elem)
        si2=sigma11(elem)-sigma22(elem)

        si1=det1*si1+det2*r1
        si2=det1*si2+det2*r2
        sigma12(elem)=det1*sigma12(elem)+det2*r3
        sigma11(elem)=0.5_8*(si1+si2)
        sigma22(elem)=0.5_8*(si1-si2)
  end do
 ! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of 
 ! Boullion et al Ocean Modelling 2013, but in an implicit mode:
 ! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
 ! and similarly for si2 and sigma12
  
end subroutine stress_tensor_m

!
! ==================================================================
! 
subroutine ssh2rhs
  ! Compute the contribution from the elevation to the rhs
  ! S.D. 30.07.2014
  use o_param
  use i_param
  use o_mesh
  use g_config
  use i_arrays
  use g_parsup
  implicit none
  
  integer       :: row, elem, elnodes(3)
  real(kind=8)  :: dx(3), dy(3), vol
  real(kind=8)  :: val3, meancos, aa, bb
  real(kind=8)  :: elevation_elem(3)

  val3=1.0_8/3.0_8
  
  ! use rhs_m and rhs_a for storing the contribution from elevation:
  do row=1, myDim_nod2d 
     rhs_a(row)=0.0
     rhs_m(row)=0.0
  end do

  do elem=1,myDim_elem2d         
     elnodes=elem2D_nodes(:,elem)
     vol=elem_area(elem)
     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)     
     elevation_elem=elevation(elnodes)
     bb=g*val3*vol
     aa=bb*sum(dx*elevation_elem)
     bb=bb*sum(dy*elevation_elem)
        rhs_a(elnodes)=rhs_a(elnodes)-aa	    
        rhs_m(elnodes)=rhs_m(elnodes)-bb
  end do
  
  
  
end subroutine ssh2rhs
!
!===================================================================
!
subroutine stress2rhs_m

  ! add internal stress to the rhs
  ! SD, 30.07.2014
  !-----------------------------------------------------------------  
  use o_param
  use i_param
  use i_therm_param
  use o_mesh
  use g_config
  use i_arrays
  use g_parsup
  implicit none
  
  integer       :: k, row, elem, elnodes(3)
  real(kind=8)  :: dx(3), dy(3), vol
  real(kind=8)  :: val3, mf, aa, bb
  real(kind=8)  :: mass, cluster_area, elevation_elem(3)

  val3=1.0_8/3.0_8
  
  do row=1, myDim_nod2d 
     u_rhs_ice(row)=0.0
     v_rhs_ice(row)=0.0
  end do

  do elem=1,myDim_elem2d         
     elnodes=elem2D_nodes(:,elem)
     if(sum(a_ice(elnodes)) < 0.01) cycle !DS
     
     vol=elem_area(elem)
     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)     
     mf=metric_factor(elem)                               !metrics

     do k=1,3
        row=elnodes(k)
        u_rhs_ice(row)=u_rhs_ice(row) - vol* &
             (sigma11(elem)*dx(k)+sigma12(elem)*dy(k))    &
       -vol*sigma12(elem)*val3*mf                         !metrics
        v_rhs_ice(row)=v_rhs_ice(row) - vol* &
             (sigma12(elem)*dx(k)+sigma22(elem)*dy(k))    &
       +vol*sigma11(elem)*val3*mf                         ! metrics
     end do
  end do
  
  do row=1, myDim_nod2d               
     mass=(m_ice(row)*rhoice+m_snow(row)*rhosno)
     mass=mass/(1.0_8+mass*mass)
        u_rhs_ice(row)=(u_rhs_ice(row)*mass + rhs_a(row))/area(1,row) 
        v_rhs_ice(row)=(v_rhs_ice(row)*mass + rhs_m(row))/area(1,row) 
  end do

end subroutine stress2rhs_m
!
!===================================================================
!
subroutine EVPdynamics_m
  ! assemble rhs and solve for ice velocity
  ! New implementation based on Bouillion et al. Ocean Modelling 2013
  ! SD 30.07.14
  !---------------------------------------------------------

  use o_param
  use i_param
  use i_therm_param
  use o_mesh
  use g_config
  use i_arrays
  use o_arrays
  use g_parsup
  use g_comm_auto

  implicit none
  integer         :: steps, shortstep, i, ed
  real(kind=8)    :: rdt, drag, det, fc
  real(kind=8)    :: inv_thickness(myDim_nod2D), umod, rhsu, rhsv
  logical         :: ice_el(myDim_elem2D)

!NR for stress_tensor_m
  integer        :: el, elnodes(3)
  real(kind=8)   :: dx(3), dy(3), msum, asum
  real(kind=8)   :: eps11, eps12, eps22, eps1, eps2, pressure, pressure_fac(myDim_elem2D), delta
  real(kind=8)   :: val3, meancos, vale
  real(kind=8)   :: det1, det2, r1, r2, r3, si1, si2

!NR for stress2rhs_m  
  integer       :: k, row
  real(kind=8)  :: vol
  real(kind=8)  :: mf, aa, bb
  real(kind=8)  :: mass(myDim_nod2D)


  
  val3=1.0_8/3.0_8
  vale=1.0_8/(ellipse**2)
  det2=1.0_8/(1.0_8+alpha_evp)
  det1=alpha_evp*det2
  rdt=ice_dt/(1.0*evp_rheol_steps)
  steps=evp_rheol_steps
  
  u_ice_aux=u_ice    ! Initialize solver variables
  v_ice_aux=v_ice

!NR inlined, to have all initialization in one place.
!  call ssh2rhs
  
  ! use rhs_m and rhs_a for storing the contribution from elevation:
  do row=1, myDim_nod2d 
     rhs_a(row)=0.0
     rhs_m(row)=0.0
  end do

  do el=1,myDim_elem2d         
     elnodes=elem2D_nodes(:,el)
     vol=elem_area(el)
     dx=gradient_sca(1:3,el)
     dy=gradient_sca(4:6,el)
     bb=g*val3*vol
     aa=bb*sum(dx*elevation(elnodes))
     bb=bb*sum(dy*elevation(elnodes))
     rhs_a(elnodes)=rhs_a(elnodes)-aa	    
     rhs_m(elnodes)=rhs_m(elnodes)-bb
  end do

! precompute thickness (the inverse is needed) and mass (scaled by area)
  do i=1,myDim_nod2D
     inv_thickness(i) = 0._8
     if (a_ice(i) >= 0.01_8) then
        inv_thickness(i) = (rhoice*m_ice(i)+rhosno*m_snow(i))/a_ice(i)
        inv_thickness(i) = 1.0_8/max(inv_thickness(i), 9.0_8)  ! Limit the mass

        mass(i) = (m_ice(i)*rhoice+m_snow(i)*rhosno)
        mass(i) = mass(i)/((1.0_8+mass(i)*mass(i))*area(1,i))

        ! scale rhs_a, rhs_m, too.
        rhs_a(i) = rhs_a(i)/area(1,i) 
        rhs_m(i) = rhs_m(i)/area(1,i) 

     endif
  enddo

! precompute pressure factor
  do el=1,myDim_elem2D
     elnodes=elem2D_nodes(:,el)

     pressure_fac(el) = 0._8
     ice_el(el) = .false.
     msum=sum(m_ice(elnodes))*val3
     if(msum > 0.01) then
        ice_el(el) = .true.
        asum=sum(a_ice(elnodes))*val3     
     
        pressure_fac(el) = pstar*msum*exp(-c_pressure*(1.0_8-asum))
     endif
  end do

  do row=1, myDim_nod2d 
     u_rhs_ice(row)=0.0
     v_rhs_ice(row)=0.0
  end do

  do shortstep=1, steps

!NR inlining, to make it easier to have local arrays and fuse loops
!NR    call stress_tensor_m
  ! Internal stress tensor
  ! New implementation following Boullion et al, Ocean Modelling 2013.
  ! SD, 30.07.2014
  !===================================================================
 
   do el=1,myDim_elem2D
     elnodes=elem2D_nodes(:,el)

     if(.not. ice_el(el)) cycle !DS

!     if(sum(a_ice(elnodes)) < 0.01) cycle !DS
     
     dx=gradient_sca(1:3,el)
     dy=gradient_sca(4:6,el)     
      ! METRICS:
     meancos=metric_factor(el)
      !  
      ! ====== Deformation rate tensor on element elem:
     eps11 = sum(dx*u_ice_aux(elnodes)) - val3*sum(v_ice_aux(elnodes))*meancos                !metrics
     eps22 = sum(dy*v_ice_aux(elnodes))
     eps12 = 0.5_8*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes)) &
            +0.5_8*val3*sum(u_ice_aux(elnodes))*meancos          !metrics 
     
      ! ======= Switch to eps1,eps2
     eps1=eps11+eps22
     eps2=eps11-eps22   
     
      ! ====== moduli:
     delta=eps1**2+vale*(eps2**2+4.0_8*eps12**2)
     delta=sqrt(delta)
    
     pressure = pressure_fac(el)/(delta+delta_min)
    
     r1 = pressure*(eps1-delta) 
     r2 = pressure*eps2*vale
     r3 = pressure*eps12*vale
     si1=sigma11(el)+sigma22(el)
     si2=sigma11(el)-sigma22(el)

     si1 = det1*si1+det2*r1
     si2 = det1*si2+det2*r2
     sigma12(el) = det1*sigma12(el)+det2*r3
     sigma11(el) = 0.5_8*(si1+si2)
     sigma22(el) = 0.5_8*(si1-si2)
!  end do   ! fuse loops
 ! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of 
 ! Boullion et al Ocean Modelling 2013, but in an implicit mode:
 ! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
 ! and similarly for si2 and sigma12

!NR inlining  call stress2rhs_m
  ! add internal stress to the rhs
  ! SD, 30.07.2014
  !-----------------------------------------------------------------  

!  do el=1,myDim_elem2d         
!     elnodes=elem2D_nodes(:,el)
     
     vol=elem_area(el)

     do k=1,3
        row=elnodes(k)
        u_rhs_ice(row)=u_rhs_ice(row) - vol* &
             (sigma11(el)*dx(k)+sigma12(el)*dy(k))    &
             -vol*sigma12(el)*val3*meancos                         !metrics

        v_rhs_ice(row)=v_rhs_ice(row) - vol* &
             (sigma12(el)*dx(k)+sigma22(el)*dy(k))    &
             +vol*sigma11(el)*val3*meancos                         ! metrics
     end do
  end do
  
  do i=1, myDim_nod2d 
     if (a_ice(i) >= 0.01) then                   ! Skip if ice is absent              

        u_rhs_ice(i) = u_rhs_ice(i)*mass(i) + rhs_a(i)
        v_rhs_ice(i) = v_rhs_ice(i)*mass(i) + rhs_m(i)
!  end do   !NR fuse loops
 !============= stress2rhs_m ends ======================

!     do i=1,myDim_nod2D

        umod = sqrt((u_ice_aux(i)-u_w(i))**2+(v_ice_aux(i)-v_w(i))**2)
        drag = rdt*Cd_oce_ice*umod*density_0*inv_thickness(i)

        !rhs for water stress, air stress, and u_rhs_ice/v (internal stress + ssh)
        rhsu = u_ice(i)+drag*u_w(i)+rdt*(inv_thickness(i)*stress_atmice_x(i)+u_rhs_ice(i)) + beta_evp*u_ice_aux(i)
        rhsv = v_ice(i)+drag*v_w(i)+rdt*(inv_thickness(i)*stress_atmice_y(i)+v_rhs_ice(i)) + beta_evp*v_ice_aux(i)

        !solve (Coriolis and water stress are treated implicitly)        
        det = bc_index_nod2D(i) / ((1.0_8+beta_evp+drag)**2 + (rdt*coriolis_node(i))**2)

        u_ice_aux(i) = det*((1.0+beta_evp+drag)*rhsu+fc*rhsv)
        v_ice_aux(i) = det*((1.0+beta_evp+drag)*rhsv-fc*rhsu)
        end if
     end do

     do  ed=1, myDim_edge2D
         ! boundary conditions
         if (myList_edge2D(ed) > edge2D_in) then
            u_ice_aux(edges(1:2,ed))=0.0_WP
            v_ice_aux(edges(1:2,ed))=0.0_WP
         endif
     end do    
     call exchange_nod_begin(u_ice_aux, v_ice_aux)

     do row=1, myDim_nod2d 
        u_rhs_ice(row)=0.0
        v_rhs_ice(row)=0.0
     end do

     call exchange_nod_end
  end do

  where (a_ice < 0.01) ! Added 28.10.14 for full compatibility with the VP solver 
        u_ice_aux=0._WP
        v_ice_aux=0._WP
  end where

  u_ice=u_ice_aux
  v_ice=v_ice_aux

end subroutine EVPdynamics_m
!
!
!
! ====================================================================
! aEVP implementation: Similar to mEVP, but alpha is variable.
! The subroutines involved are with _a.
! ====================================================================
!
subroutine find_alpha_field_a
  ! EVP stability parameter alpha is computed at each element
  ! aEVP implementation
  ! SD, 13.02.2017
  ! ==================================================================
  use o_param
  use i_param
  use i_therm_param
  use o_mesh
  use g_config
  use i_arrays
  use g_parsup
  implicit none

  integer        :: elem, elnodes(3)
  real(kind=8)   :: dx(3), dy(3), msum, asum
  real(kind=8)   :: eps11, eps12, eps22, eps1, eps2, pressure, delta
  real(kind=8)   :: val3, meancos, usum, vsum, vale
  val3=1.0_8/3.0_8
  vale=1.0_8/(ellipse**2)
   do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)

     msum=sum(m_ice(elnodes))*val3
     if(msum<=0.01) cycle !DS
     asum=sum(a_ice(elnodes))*val3
     
     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)     
     ! METRICS:
     vsum=sum(v_ice_aux(elnodes))
     usum=sum(u_ice_aux(elnodes))
     meancos=metric_factor(elem)
     !  
     ! ====== Deformation rate tensor on element elem:
     eps11=sum(dx*u_ice_aux(elnodes))
     eps11=eps11-val3*vsum*meancos                !metrics
     eps22=sum(dy*v_ice_aux(elnodes))
     eps12=0.5_8*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes))
     eps12=eps12+0.5_8*val3*usum*meancos          !metrics 
     
      ! ======= Switch to eps1,eps2
     eps1=eps11+eps22
     eps2=eps11-eps22   
     
      ! ====== moduli:
     delta=eps1**2+vale*(eps2**2+4.0_8*eps12**2)
     delta=sqrt(delta)
         
     pressure=pstar*exp(-c_pressure*(1.0_8-asum))/(delta+delta_min) ! no multiplication
                                                                    ! with thickness (msum) 
     alpha_evp_array(elem)=max(50.0,sqrt(c_aevp*pressure/rhoice/elem_area(elem)))  
      ! /voltriangle(elem) for FESOM1.4
      ! We do not allow alpha to be too small!
   end do
  end subroutine find_alpha_field_a  
! ====================================================================

subroutine stress_tensor_a
  ! Internal stress tensor
  ! New implementation following Boullion et al, Ocean Modelling 2013.
  ! and Kimmritz et al., Ocean Modelling 2016
  ! SD, 14.02.2017
  !===================================================================
  use o_param
  use i_param
  use o_mesh
  use g_config
  use i_arrays
  use g_parsup
  implicit none

  integer        :: elem, elnodes(3)
  real(kind=8)   :: dx(3), dy(3), msum, asum
  real(kind=8)   :: eps11, eps12, eps22, eps1, eps2, pressure, delta
  real(kind=8)   :: val3, meancos, usum, vsum, vale
  real(kind=8)   :: det1, det2, r1, r2, r3, si1, si2
  
  val3=1.0_8/3.0_8
  vale=1.0_8/(ellipse**2)
   do elem=1,myDim_elem2D
     det2=1.0_8/(1.0_8+alpha_evp_array(elem))     ! Take alpha from array
     det1=alpha_evp_array(elem)*det2
  
     elnodes=elem2D_nodes(:,elem)

     msum=sum(m_ice(elnodes))*val3
     if(msum<=0.01) cycle !DS
     asum=sum(a_ice(elnodes))*val3
     
     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)     
     ! METRICS:
     vsum=sum(v_ice_aux(elnodes))
     usum=sum(u_ice_aux(elnodes))
     meancos=metric_factor(elem)
     !  
     ! ====== Deformation rate tensor on element elem:
     eps11=sum(dx*u_ice_aux(elnodes))
     eps11=eps11-val3*vsum*meancos                !metrics
     eps22=sum(dy*v_ice_aux(elnodes))
     eps12=0.5_8*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes))
     eps12=eps12+0.5_8*val3*usum*meancos          !metrics 
     
      ! ======= Switch to eps1,eps2
     eps1=eps11+eps22
     eps2=eps11-eps22   
     
      ! ====== moduli:
     delta=eps1**2+vale*(eps2**2+4.0_8*eps12**2)
     delta=sqrt(delta)
    
     pressure=pstar*msum*exp(-c_pressure*(1.0_8-asum))/(delta+delta_min)
    
        r1=pressure*(eps1-delta) 
        r2=pressure*eps2*vale
        r3=pressure*eps12*vale
        si1=sigma11(elem)+sigma22(elem)
        si2=sigma11(elem)-sigma22(elem)

        si1=det1*si1+det2*r1
        si2=det1*si2+det2*r2
        sigma12(elem)=det1*sigma12(elem)+det2*r3
        sigma11(elem)=0.5_8*(si1+si2)
        sigma22(elem)=0.5_8*(si1-si2)
  end do
 ! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of 
 ! Boullion et al Ocean Modelling 2013, but in an implicit mode:
 ! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
 ! and similarly for si2 and sigma12
  
end subroutine stress_tensor_a
!
!===================================================================
!
subroutine EVPdynamics_a
  ! assemble rhs and solve for ice velocity
  ! New implementation based on Bouillion et al. Ocean Modelling 2013
  ! and Kimmritz et al., Ocean Modelling  2016 
  ! SD 14.02.17
  !---------------------------------------------------------

use o_mesh
use o_param
use o_mesh
use i_arrays
USE o_arrays
use i_param
use o_PARAM
use i_therm_param
use g_parsup
use g_comm_auto

  implicit none
  integer         :: steps, shortstep, i, ed
  real(kind=8)    :: rdt, drag, det, fc
  real(kind=8)    :: thickness, inv_thickness, umod, rhsu, rhsv
  REAL(kind=8)    :: t0,t1, t2, t3, t4, t5, t00, txx
 
  steps=evp_rheol_steps
  rdt=ice_dt/(1.0*evp_rheol_steps)  
  u_ice_aux=u_ice    ! Initialize solver variables
  v_ice_aux=v_ice
  call ssh2rhs
 
  do shortstep=1, steps 
     call stress_tensor_a
     call stress2rhs_m    ! _m=_a, so no _m version is the only one!
     do i=1,myDim_nod2D 
         thickness=(rhoice*m_ice(i)+rhosno*m_snow(i))/max(a_ice(i),0.01)
         thickness=max(thickness, 9.0)   ! Limit if it is too small (0.01 m)
         inv_thickness=1.0_8/thickness

         umod=sqrt((u_ice_aux(i)-u_w(i))**2+(v_ice_aux(i)-v_w(i))**2)
         drag=rdt*Cd_oce_ice*umod*density_0*inv_thickness

         !rhs for water stress, air stress, and u_rhs_ice/v (internal stress + ssh)
         rhsu=u_ice(i)+drag*u_w(i)+rdt*(inv_thickness*stress_atmice_x(i)+u_rhs_ice(i))
         rhsv=v_ice(i)+drag*v_w(i)+rdt*(inv_thickness*stress_atmice_y(i)+v_rhs_ice(i))

         rhsu=beta_evp*u_ice_aux(i)+rhsu
	 rhsv=beta_evp*v_ice_aux(i)+rhsv
         !solve (Coriolis and water stress are treated implicitly)
         fc=rdt*coriolis_node(i)
         det=(1.0_8+beta_evp+drag)**2+fc**2
         det=bc_index_nod2D(i)/det
         u_ice_aux(i)=det*((1.0+beta_evp_array(i)+drag)*rhsu+fc*rhsv)
         v_ice_aux(i)=det*((1.0+beta_evp_array(i)+drag)*rhsv-fc*rhsu)
     end do
     do  ed=1, myDim_edge2D
         ! boundary conditions
         if (myList_edge2D(ed) > edge2D_in) then
            u_ice_aux(edges(1:2,ed))=0.0_WP
            v_ice_aux(edges(1:2,ed))=0.0_WP
         endif
     end do    
     call exchange_nod(u_ice_aux, v_ice_aux)
  end do
  
    do i=1, myDim_nod2D+eDim_nod2D   ! Added 28.10.14 for full compatibility with 
                                     ! the VP solver 
    if (a_ice(i)<=0.01) then
       u_ice_aux(i)=0.0
       v_ice_aux(i)=0.0
    end if
    end do 
    
    u_ice=u_ice_aux
    v_ice=v_ice_aux
 
  call find_alpha_field_a             ! alpha_evp_array is initialized with alpha_evp;
                                      ! At this stage we already have non-trivial velocities. 
  call find_beta_field_a
end subroutine EVPdynamics_a
!
! =================================================================
!
 
subroutine find_beta_field_a 
! beta_evp_array is defined at nodes, and this is the only 
! reason we need it in addition to alpha_evp_array (we work with 
! alpha=beta, and keep different names for generality; mEVP can work with 
! alpha \ne beta, but not aEVP).

use o_mesh
use o_param
USE i_param
use i_arrays
use g_parsup 
Implicit none
integer :: n
    DO n=1, myDim_nod2D
       ! ==============
       ! FESOM1.4 and stand-alone FESIM
       ! beta_evp_array(n) =  maxval(alpha_evp_array(nod_in_elem2D(n)%addresses(1:nod_in_elem2D(n)%nmb)))
       ! ==============
       ! FESOM2.0
       beta_evp_array(n) =  maxval(alpha_evp_array(nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
    END DO

end subroutine find_beta_field_a
! 
! ================================================================
!
