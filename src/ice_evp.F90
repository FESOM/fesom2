!
! Contains: Three routines of EVP dynamics+auxiliary routine that computes 
! velocity gradients. The driving routine is EVPdynamics.
!
!===========================================================================
SUBROUTINE ice_vel_gradients
! Compute derivatives of velocity by least square interpolation.
! For the no-slip case, it is assumed that velocity at 
! the boundary edge == 0. For the free-slip case, there are only 2
! neighbours
USE o_MESH
USE i_ARRAYS
USE i_PARAM
USE g_PARSUP
USE g_comm
IMPLICIT NONE
real(kind=WP)    :: x(3), y(3), u, v, r1, r2
real(kind=WP)    :: zc(2), un, grad_aux(4)
real(kind=WP)    :: u_aux, v_aux
integer         :: elem, el, j, nz, m
   DO elem=1, myDim_elem2D     !! P elem=1,elem2D
      grad_aux=0.0_WP
      DO j=1,3
	 el=elem_neighbors(j,elem)
         if(el>0) then
	   u=U_ice(el)-U_ice(elem)
           v=V_ice(el)-V_ice(elem)
         	grad_aux(1)=grad_aux(1)+gradient_vec(j,elem)*u
		grad_aux(2)=grad_aux(2)+gradient_vec(j+3,elem)*u
		grad_aux(3)=grad_aux(3)+gradient_vec(j,elem)*v
		grad_aux(4)=grad_aux(4)+gradient_vec(j+3,elem)*v
	  else
	  ! ===============
	  ! Boundary element
	  ! ===============
	  if(ice_free_slip) then
	       zc=edge_dxdy(:,elem_edges(j,elem))
	       zc(1)=zc(1)*elem_cos(elem)
	       un=-2*(U_ice(elem)*zc(2)-V_ice(elem)*zc(1))/sum(zc*zc)
               u=U_ice(elem)+un*zc(2)
               v=V_ice(elem)-un*zc(1)
	  else     ! noslip
	       u=-U_ice(elem)
	       v=-V_ice(elem)     
	  end if 
		u=u-U_ice(elem)
                v=v-V_ice(elem)
         	grad_aux(1)=grad_aux(1)+gradient_vec(j,elem)*u
		grad_aux(2)=grad_aux(2)+gradient_vec(j+3,elem)*u
		grad_aux(3)=grad_aux(3)+gradient_vec(j,elem)*v
		grad_aux(4)=grad_aux(4)+gradient_vec(j+3,elem)*v
	  end if
       END DO   ! cycle over neighbor elements
 	 ice_grad_vel(:,elem)=grad_aux(:)
	 
  END DO	  
 call exchange_e2D(ice_grad_vel, 4)

END SUBROUTINE ice_vel_gradients
!===========================================================================
subroutine stress_tensor
! EVP rheology. The routine computes stress tensor components based on ice 
! velocity field. They are stored as elemental arrays (sigma11, sigma22 and
! sigma12). 
use o_param
use i_param
use o_mesh
use i_arrays
use g_parsup
USE g_CONFIG
use g_comm_auto
implicit none

real(kind=WP)   :: eps11, eps12, eps22, eta, xi, ice_strength, delta, aa
integer         :: elem, elnodes(3)
real(kind=WP)   :: asum, msum, vale, dx(3), dy(3)
real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2, dte 
real(kind=WP)   :: zeta, delta_inv

  vale=1.0_WP/(ellipse**2)
   
  dte=ice_dt/(1.0_WP*evp_rheol_steps)
  det1=1.0_WP+0.5_WP*Tevp_inv*dte
  det2=1.0_WP+0.5_WP*Tevp_inv*dte*ellipse**2 
     
  det1=1.0_WP/det1
  det2=1.0_WP/det2

  call ice_vel_gradients
  
  do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
      ! ===== Check if there is ice on elem
        aa=product(m_ice(elnodes))*product(a_ice(elnodes))
        if (aa==0) CYCLE     ! There is no ice in elem
      ! =====	
      ! ===== Deformation rate tensor on element elem:
     eps11=ice_grad_vel(1,elem)-V_ice(elem)*metric_factor(elem)   
     eps22=ice_grad_vel(4,elem)   
    
     eps12=0.5_WP*(ice_grad_vel(2,elem)+ice_grad_vel(3,elem)+ &
                  U_ice(elem)*metric_factor(elem))      
      ! ===== moduli:
     delta=(eps11**2+eps22**2)*(1.0_WP+vale)+4.0_WP*vale*eps12**2 + &
            2.0_WP*eps11*eps22*(1.0_WP-vale)
     delta=sqrt(delta)
     msum=sum(m_ice(elnodes))/3.0_WP
     asum=sum(a_ice(elnodes))/3.0_WP
     
      ! ===== Hunke and Dukowicz c*h*p*
     ice_strength=pstar*(msum)*exp(-c_pressure*(1.0_WP-asum))
      ! =======================================
      ! ===== Here the EVP rheology piece starts
      ! =======================================
     ice_strength=0.5_WP*ice_strength
      ! ===== viscosity zeta should exceed zeta_min
      ! (done via limiting delta from above)
      
      !if(delta>pressure/zeta_min) delta=pressure/zeta_min
           !It does not work properly by 
	   !creating response where ice_strength is small
           ! Uncomment and test if necessary
      
      ! ===== if delta is too small or zero, viscosity will too large (unlimited)
      ! (limit delta_inv)
     delta_inv=1.0_WP/max(delta,delta_min) 
     zeta=ice_strength*delta_inv			     
      ! ===== Limiting pressure/Delta  (zeta): it may still happen that pressure/Delta 
      ! is too large in some regions and CFL criterion is violated.
      ! The regularization below was introduced by Hunke, 
      ! but seemingly is not used in the current CICE. 
      ! Without it divergence and zeta can be noisy (but code 
      ! remains stable), using it reduces viscosities too strongly.
      ! It is therefore commented
      
      !if (zeta>Clim_evp*voltriangle(elem)) then
      !zeta=Clim_evp*voltriangle(elem)
      !end if 
      
      ice_strength=ice_strength*Tevp_inv
      zeta=zeta*Tevp_inv
      				     
     r1=zeta*(eps11+eps22) - ice_strength
     r2=zeta*(eps11-eps22)
     r3=zeta*eps12
     si1=sigma11(elem)+sigma22(elem)
     si2=sigma11(elem)-sigma22(elem)
     
     si1=det1*(si1+dte*r1)
     si2=det2*(si2+dte*r2)
     sigma12(elem)=det2*(sigma12(elem)+dte*r3)
     sigma11(elem)=0.5_WP*(si1+si2)
     sigma22(elem)=0.5_WP*(si1-si2)
    end do
     call exchange_elem(sigma11)
     call exchange_elem(sigma12)
     call exchange_elem(sigma22)
end subroutine stress_tensor

!===================================================================
subroutine stress2rhs
! EVP implementation:
! Computes the divergence of stress tensor and puts the result into the
! rhs vectors 

USE o_MESH
USE o_PARAM
USE i_PARAM
USE i_therm_param
USE i_arrays
USE g_PARSUP
 

IMPLICIT NONE
INTEGER      :: elem, ed, elnodes(3), el(2)  
REAL(kind=8) :: mass, uc, vc, xe, ye



 DO elem=1, myDim_elem2D
     U_rhs_ice(elem)=0.0_WP
     V_rhs_ice(elem)=0.0_WP
 END DO
 
 ! Stress divergence
 DO  ed=1,myDim_edge2D+eDim_edge2D   
   el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)
   ye=edge_dxdy(2,ed)         ! xe, ye are in radians!
   if(myList_edge2D(ed)<=edge2D_in) then     ! elements on both sides
   uc=0.5_WP*r_earth*(ye*(sigma11(el(1))+sigma11(el(2)))- &
                  xe*(elem_cos(el(1))*sigma12(el(1))+ &
		  elem_cos(el(2))*sigma12(el(2))))
   vc=0.5_WP*r_earth*(ye*(sigma12(el(1))+sigma12(el(2)))- &
                  xe*(elem_cos(el(1))*sigma22(el(1))+ &
		  elem_cos(el(2))*sigma22(el(2))))
   U_rhs_ice(el(1))=U_rhs_ice(el(1))+uc
   U_rhs_ice(el(2))=U_rhs_ice(el(2))-uc
   V_rhs_ice(el(1))=V_rhs_ice(el(1))+vc
   V_rhs_ice(el(2))=V_rhs_ice(el(2))-vc
   else
   uc=r_earth*(ye*sigma11(el(1))- xe*elem_cos(el(1))*sigma12(el(1)))
   vc=r_earth*(ye*sigma12(el(1))- xe*elem_cos(el(1))*sigma22(el(1)))
   U_rhs_ice(el(1))=U_rhs_ice(el(1))+uc
   V_rhs_ice(el(1))=V_rhs_ice(el(1))+vc
   end if
 END DO
 do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     mass=elem_area(elem)*(rhoice*sum(m_ice(elnodes))+rhosno*sum(m_snow(elnodes)))/3.0_WP
     if (mass>0) then
     U_rhs_ice(elem)=U_rhs_ice(elem)/mass - &
                     g*sum(gradient_sca(1:3,elem)*elevation(elnodes)) 
     V_rhs_ice(elem)=V_rhs_ice(elem)/mass - &
                     g*sum(gradient_sca(4:6,elem)*elevation(elnodes)) 
     else
     U_rhs_ice(elem)=0.0_WP 
     V_rhs_ice(elem)=0.0_WP
     end if
 END DO
end subroutine stress2rhs
!===================================================================
subroutine EVPdynamics
! EVP implementation. Does subcycling and boundary conditions.  
USE o_MESH
USE o_PARAM
USE i_ARRAYS
USE i_PARAM
USE i_therm_param
USE g_PARSUP
USE o_ARRAYS
USE g_CONFIG
use g_comm_auto
 
IMPLICIT NONE
integer          :: steps, shortstep
real(kind=WP)    :: rdt, am, sm, mm
real(kind=WP)    :: drag, inv_mass, det, umod, rhsu, rhsv
integer          :: elem, elnodes(3)
real(kind=WP)    :: ax, ay, yx, tx, ty
real(kind=WP)    :: t0,t1,t2,t3,t4
    rdt=ice_dt/(1.0*evp_rheol_steps)
    ax=cos(theta_io)
    ay=sin(theta_io)

do shortstep=1, 2 !evp_rheol_steps 
 t0=MPI_Wtime()     
 call stress_tensor
 t1=MPI_Wtime()     
! if(shortstep<3) write(*,*) 'stress ', mype, minval(sigma11(1:myDim_elem2D))
 call stress2rhs
 t2=MPI_Wtime()     
! if(shortstep<3)  write(*,*) 'rhs  ', mype, maxval(U_rhs_ice(1:myDim_elem2D)),minval(U_rhs_ice(1:myDim_elem2D))
 do elem=1,myDim_elem2D 
    elnodes=elem2D_nodes(:,elem)
    am=sum(a_ice(elnodes))/3.0_WP
    mm=sum(m_ice(elnodes))/3.0_WP
    sm=sum(m_snow(elnodes))/3.0_WP
    if (am < 0.01) cycle               ! Skip if ice is absent
    inv_mass=(rhoice*mm+rhosno*sm)/am
    inv_mass=max(inv_mass, 9.0)        ! Limit the mass 
                                       ! if it is too small
    inv_mass=1.0/inv_mass
    umod=sqrt((U_ice(elem)-U_w(elem))**2+(V_ice(elem)-V_w(elem))**2)
    drag=Cd_oce_ice*umod*density_0*inv_mass
    tx=sum(stress_atmice_x(elnodes))/3.0_WP
    ty=sum(stress_atmice_y(elnodes))/3.0_WP
    rhsu=U_ice(elem)+rdt*(drag*(ax*U_w(elem)-ay*V_w(elem))+ &
                    inv_mass*tx+U_rhs_ice(elem))
    rhsv=V_ice(elem)+rdt*(drag*(ax*V_w(elem)+ay*U_w(elem))+ &
                    inv_mass*ty+V_rhs_ice(elem))

    det=(1.+ax*drag*rdt)**2+(rdt*coriolis(elem)+rdt*ay*drag)**2
    det=1.0_WP/det
    U_ice(elem)=det*((1.0+ax*drag*rdt)*rhsu+rdt*(coriolis(elem)+ay*drag)*rhsv)
    V_ice(elem)=det*((1.0+ax*drag*rdt)*rhsv-rdt*(coriolis(elem)+ay*drag)*rhsu)
 end do
 t3=MPI_Wtime()     
 call exchange_elem(U_ice)
 call exchange_elem(V_ice)    
 t4=MPI_Wtime()     
 END DO
 if (mype==1) write(*,*) 'EVP Timing:',t1-t0,t2-t0,t3-t0,t4-t0
end subroutine EVPdynamics
!===================================================================
subroutine stress_tensor_n
! EVP rheology. The routine computes stress tensor components based on ice 
! velocity field. They are stored as elemental arrays (sigma11, sigma22 and
! sigma12). The ocean velocity is at nodal locations.
use o_param
use i_param
use o_mesh
use i_arrays
use g_parsup
USE g_CONFIG
implicit none

real(kind=WP)   :: eps11, eps12, eps22, eta, xi, ice_strength, delta, aa
integer         :: elem, elnodes(3)
real(kind=WP)   :: asum, msum, vale, dx(3), dy(3)
real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2, dte 
real(kind=WP)   :: zeta, delta_inv, d1, d2

  vale=1.0_WP/(ellipse**2)
   
  dte=ice_dt/(1.0_WP*evp_rheol_steps)
  det1=1.0_WP+0.5_WP*Tevp_inv*dte
  det2=1.0_WP+0.5_WP*Tevp_inv*dte*ellipse**2 
     
  det1=1.0_WP/det1
  det2=1.0_WP/det2

  do elem=1,myDim_elem2D
     
     elnodes=elem2D_nodes(:,elem)
      ! ===== Check if there is ice on elem
        aa=product(m_ice(elnodes))*product(a_ice(elnodes))
        if (aa==0) CYCLE     ! There is no ice in elem
      ! =====	
      ! ===== Deformation rate tensor on element elem:
           !du/dx
       d1=sum(gradient_sca(1:3, elem)*U_ice(elnodes))
       d2=sum(V_ice(elnodes))/3.0_WP
     eps11=d1-d2*metric_factor(elem)   
     eps22=sum(gradient_sca(4:6, elem)*V_ice(elnodes))
       d1=sum(gradient_sca(4:6, elem)*U_ice(elnodes)) + &
          sum(gradient_sca(1:3, elem)*V_ice(elnodes))
       d2=sum(U_ice(elnodes))/3.0_WP	   
     eps12=0.5_WP*(d1+ d2*metric_factor(elem))      
      ! ===== moduli:
     delta=(eps11**2+eps22**2)*(1.0_WP+vale)+4.0_WP*vale*eps12**2 + &
            2.0_WP*eps11*eps22*(1.0_WP-vale)
     delta=sqrt(delta)
     msum=sum(m_ice(elnodes))/3.0_WP
     asum=sum(a_ice(elnodes))/3.0_WP
     
      ! ===== Hunke and Dukowicz c*h*p*
     ice_strength=pstar*(msum)*exp(-c_pressure*(1.0_WP-asum))
      ! =======================================
      ! ===== Here the EVP rheology piece starts
      ! =======================================
     ice_strength=0.5_WP*ice_strength
      ! ===== viscosity zeta should exceed zeta_min
      ! (done via limiting delta from above)
      
      !if(delta>pressure/zeta_min) delta=pressure/zeta_min
           !It does not work properly by 
	   !creating response where ice_strength is small
           ! Uncomment and test if necessary
      
      ! ===== if delta is too small or zero, viscosity will too large (unlimited)
      ! (limit delta_inv)
     delta_inv=1.0_WP/max(delta,delta_min) 
     zeta=ice_strength*delta_inv			     
      ! ===== Limiting pressure/Delta  (zeta): it may still happen that pressure/Delta 
      ! is too large in some regions and CFL criterion is violated.
      ! The regularization below was introduced by Hunke, 
      ! but seemingly is not used in the current CICE. 
      ! Without it divergence and zeta can be noisy (but code 
      ! remains stable), using it reduces viscosities too strongly.
      ! It is therefore commented
      
      !if (zeta>Clim_evp*voltriangle(elem)) then
      !zeta=Clim_evp*voltriangle(elem)
      !end if 
      
      ice_strength=ice_strength*Tevp_inv
      zeta=zeta*Tevp_inv
      				     
     r1=zeta*(eps11+eps22) - ice_strength
     r2=zeta*(eps11-eps22)
     r3=zeta*eps12
     si1=sigma11(elem)+sigma22(elem)
     si2=sigma11(elem)-sigma22(elem)
     
     si1=det1*(si1+dte*r1)
     si2=det2*(si2+dte*r2)
     sigma12(elem)=det2*(sigma12(elem)+dte*r3)
     sigma11(elem)=0.5_WP*(si1+si2)
     sigma22(elem)=0.5_WP*(si1-si2)
    end do
    
end subroutine stress_tensor_n
!===================================================================
subroutine stress2rhs_e
! EVP implementation:
! Computes the divergence of stress tensor and puts the result into the
! rhs vectors. Velocity is at nodes. 
! The divergence is computed in a cysly over edges. It is slower that the
! approach in stress2rhs_e inherited from FESOM


USE o_MESH
USE o_PARAM
USE i_PARAM
USE i_therm_param
USE i_arrays
USE g_PARSUP


IMPLICIT NONE
INTEGER      :: n, elem, ed, elnodes(3), el(2), ednodes(2)  
REAL(kind=8) :: mass, uc, vc,  deltaX1, deltaX2, deltaY1, deltaY2

 DO n=1, myDim_nod2D
     U_rhs_ice(n)=0.0
     V_rhs_ice(n)=0.0
 END DO
 
 ! Stress divergence
 DO  ed=1,myDim_edge2D
   ednodes=edges(:,ed) 
   el=edge_tri(:,ed)
   deltaX1=edge_cross_dxdy(1,ed)
   deltaY1=edge_cross_dxdy(2,ed)
   deltaX2=edge_cross_dxdy(3,ed)
   deltaY2=edge_cross_dxdy(4,ed)
   
   if(myList_edge2D(ed)>edge2D_in) cycle     
   ! elements on both sides
   uc=-sigma12(el(1))*deltaX1+sigma11(el(1))*deltaY1
   vc=-sigma22(el(1))*deltaX1+sigma12(el(1))*deltaY1
   uc=uc+sigma12(el(2))*deltaX2-sigma11(el(2))*deltaY2
   vc=vc+sigma22(el(2))*deltaX2-sigma12(el(2))*deltaY2
   U_rhs_ice(ednodes(1))=U_rhs_ice(ednodes(1))+uc
   U_rhs_ice(ednodes(2))=U_rhs_ice(ednodes(2))-uc
   V_rhs_ice(ednodes(1))=V_rhs_ice(ednodes(1))+vc
   V_rhs_ice(ednodes(2))=V_rhs_ice(ednodes(2))-vc
   END DO
 
 DO n=1, myDim_nod2D
      mass=area(1,n)*(rhoice*m_ice(n)+rhosno*m_snow(n)) 
      if(mass>0.) then 
         U_rhs_ice(n)=U_rhs_ice(n)/mass
         V_rhs_ice(n)=V_rhs_ice(n)/mass
      else
         U_rhs_ice(n)=0.0_WP
         V_rhs_ice(n)=0.0_WP
      end if
 END DO
 !
 ! elevation gradient contribution      
 !
 do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     uc=elem_area(elem)*g*sum(gradient_sca(1:3,elem)*elevation(elnodes))/3.0_WP
     vc=elem_area(elem)*g*sum(gradient_sca(4:6,elem)*elevation(elnodes))/3.0_WP
     U_rhs_ice(elnodes)=U_rhs_ice(elnodes) - uc/area(1,elnodes)
     V_rhs_ice(elnodes)=V_rhs_ice(elnodes) - vc/area(1,elnodes)
 END DO
end subroutine stress2rhs_e
!===================================================================
subroutine stress2rhs_n
! EVP implementation:
! Computes the divergence of stress tensor and puts the result into the
! rhs vectors 

USE o_MESH
USE o_PARAM
USE i_PARAM
USE i_THERM_PARAM
USE g_PARSUP
USE i_arrays

IMPLICIT NONE
INTEGER      :: n, elem, elnodes(3), k, ed, ednodes(2) 
REAL(kind=8) :: mass, aa
REAL(kind=8) :: elevation_elem(3)
REAL(kind=8) :: dx(3), dy(3), val3

val3=1/3.0_WP
 DO  n=1, myDim_nod2D
     U_rhs_ice(n)=0.0
     V_rhs_ice(n)=0.0
     rhs_a(n)=0.0       ! these are used as temporal storage here
     rhs_m(n)=0.0       ! for the contribution due to ssh
 END DO
 
 do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
      ! ===== Skip if ice is absent
     aa=product(m_ice(elnodes))*product(a_ice(elnodes))
     if (aa==0) CYCLE
      ! =====

     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)
     elevation_elem=elevation(elnodes)
     
     DO k=1,3
        n=elnodes(k)
        U_rhs_ice(n)=U_rhs_ice(n) - elem_area(elem) * &
             (sigma11(elem)*dx(k)+sigma12(elem)*(dy(k)) &
             +sigma12(elem)*val3*metric_factor(elem))            !metrics
        V_rhs_ice(n)=V_rhs_ice(n) - elem_area(elem) * &
             (sigma12(elem)*dx(k)+sigma22(elem)*dy(k) &   
             -sigma11(elem)*val3*metric_factor(elem))
     END DO
      ! use rhs_m and rhs_a for storing the contribution from elevation:
     aa=9.81*elem_area(elem)/3.0_WP
     DO k=1,3
        n=elnodes(k)
        rhs_a(n)=rhs_a(n)-aa*sum(dx*elevation_elem)	    
        rhs_m(n)=rhs_m(n)-aa*sum(dy*elevation_elem)
     END DO     
 end do 
 
  DO n=1, myDim_nod2D
     mass=area(1,n)*(m_ice(n)*rhoice+m_snow(n)*rhosno)
     
     if (mass.ne.0) then
     U_rhs_ice(n)=U_rhs_ice(n)/mass + rhs_a(n)/area(1,n) 
     V_rhs_ice(n)=V_rhs_ice(n)/mass + rhs_m(n)/area(1,n)
     else
     U_rhs_ice(n)=0.  
     V_rhs_ice(n)=0.
     end if
  END DO 
  
end subroutine stress2rhs_n
!===================================================================
!===================================================================
subroutine EVPdynamics_n
! EVP implementation. Does subcycling and boundary conditions.  
! Velocities at nodes
USE o_MESH
USE o_PARAM
USE i_ARRAYS
USE i_PARAM
USE i_therm_param
USE g_PARSUP
USE o_ARRAYS
USE g_CONFIG
USE g_comm_auto

IMPLICIT NONE
integer          :: steps, shortstep
real(kind=WP)    :: rdt, am, sm, mm
real(kind=WP)    :: drag, inv_mass, det, umod, rhsu, rhsv
integer          :: n, ed, ednodes(2)
real(kind=WP)    :: ax, ay, yx, tx, ty
    rdt=ice_dt/(1.0*evp_rheol_steps)
    ax=cos(theta_io)
    ay=sin(theta_io)
    
    
do shortstep=1, evp_rheol_steps 
 
 call stress_tensor_n
 
 !if(shortstep<3) write(*,*) mype,'stress ', minval(sigma11(1:myDim_elem2D))
 call stress2rhs_n
 
 !if(shortstep<3)  write(*,*) mype, 'rhs  ', maxval(U_rhs_ice(1:myDim_nod2D))
 do n=1,myDim_nod2D 
    am=a_ice(n)
    mm=m_ice(n)
    sm=m_snow(n)
    if (am < 0.01) cycle               ! Skip if ice is absent
    inv_mass=(rhoice*mm+rhosno*sm)/am
    inv_mass=max(inv_mass, 9.0)        ! Limit the mass 
                                       ! if it is too small
    inv_mass=1.0/inv_mass
    umod=sqrt((U_ice(n)-U_w(n))**2+(V_ice(n)-V_w(n))**2)
    drag=Cd_oce_ice*umod*density_0*inv_mass
    tx=stress_atmice_x(n)
    ty=stress_atmice_y(n)
    rhsu=U_ice(n)+rdt*(drag*(ax*U_w(n)-ay*V_w(n))+ &
                    inv_mass*tx+U_rhs_ice(n))
    rhsv=V_ice(n)+rdt*(drag*(ax*V_w(n)+ay*U_w(n))+ &
                    inv_mass*ty+V_rhs_ice(n))

    det=(1.+ax*drag*rdt)**2+(rdt*coriolis_node(n)+rdt*ay*drag)**2
    det=1.0_WP/det
    U_ice(n)=det*((1.0+ax*drag*rdt)*rhsu+rdt*(coriolis_node(n)+ay*drag)*rhsv)
    V_ice(n)=det*((1.0+ax*drag*rdt)*rhsv-rdt*(coriolis_node(n)+ay*drag)*rhsu)
 end do
 DO  ed=1,myDim_edge2D
   ! boundary conditions
   ednodes=edges(:,ed) 
   if(myList_edge2D(ed)<=edge2D_in) cycle 
   U_ice(ednodes(1))=0.0_WP
   U_ice(ednodes(2))=0.0_WP
   V_ice(ednodes(1))=0.0_WP
   V_ice(ednodes(2))=0.0_WP
 end do  
 
 call exchange_nod(U_ice)
 call exchange_nod(V_ice)    
 !if(shortstep<5)  write(*,*) mype, 'vel  ', maxval(U_ice(1:myDim_nod2D))
 END DO
 
end subroutine EVPdynamics_n
