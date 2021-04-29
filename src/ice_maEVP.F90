module ice_maEVP_interfaces
  interface
    subroutine ssh2rhs(mesh)
      use mod_mesh
      type(t_mesh), intent(in), target  :: mesh
    end subroutine

    subroutine stress_tensor_a(mesh)
      use mod_mesh
      type(t_mesh), intent(in), target  :: mesh
    end subroutine

    subroutine stress2rhs_m(mesh)
      use mod_mesh
      type(t_mesh), intent(in), target  :: mesh
    end subroutine

    subroutine find_alpha_field_a(mesh)
      use mod_mesh
      type(t_mesh), intent(in), target  :: mesh
    end subroutine

    subroutine find_beta_field_a(mesh)
      use mod_mesh
      type(t_mesh), intent(in), target  :: mesh
    end subroutine
  end interface  
end module

! ====================================================================
! New evp implementation following Bouillion et al. 2013
! and Kimmritz et al. 2015 (mEVP) and Kimmritz et al. 2016 (aEVP)
! ====================================================================  
subroutine stress_tensor_m(mesh)
  ! Internal stress tensor
  ! New implementation following Boullion et al, Ocean Modelling 2013.
  ! SD, 30.07.2014
  !===================================================================
  use o_param
  use i_param
  use mod_mesh
  use g_config
  use i_arrays
  use g_parsup

#if defined (__icepack)
use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
#endif

  implicit none

  integer         :: elem, elnodes(3)
  real(kind=WP)   :: dx(3), dy(3), msum, asum
  real(kind=WP)   :: eps1, eps2, pressure, delta
  real(kind=WP)   :: val3, meancos, usum, vsum, vale
  real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2
  type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"

  val3=1.0_WP/3.0_WP
  vale=1.0_WP/(ellipse**2)
  det2=1.0_WP/(1.0_WP+alpha_evp)
  det1=alpha_evp*det2
   do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     !_______________________________________________________________________
	 ! if element has any cavity node skip it 
     if (ulevels(elem) > 1) cycle
     
     msum=sum(m_ice(elnodes))*val3
     if(msum<=0.01_WP) cycle !DS
     asum=sum(a_ice(elnodes))*val3
     
     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)     
      ! METRICS:
          vsum=sum(v_ice_aux(elnodes))
          usum=sum(u_ice_aux(elnodes))
          meancos=metric_factor(elem)
      !  
      ! ====== Deformation rate tensor on element elem:
     eps11(elem)=sum(dx*u_ice_aux(elnodes))
     eps11(elem)=eps11(elem)-val3*vsum*meancos                !metrics
     eps22(elem)=sum(dy*v_ice_aux(elnodes))
     eps12(elem)=0.5_WP*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes))
     eps12(elem)=eps12(elem)+0.5_WP*val3*usum*meancos          !metrics 
     
      ! ======= Switch to eps1,eps2
     eps1=eps11(elem)+eps22(elem)
     eps2=eps11(elem)-eps22(elem)   
     
      ! ====== moduli:
     delta=eps1**2+vale*(eps2**2+4.0_WP*eps12(elem)**2)
     delta=sqrt(delta)
    
#if defined (__icepack)
     pressure = sum(strength(elnodes))*val3/max(delta,delta_min)
#else
     pressure=pstar*msum*exp(-c_pressure*(1.0_WP-asum))/max(delta,delta_min)
#endif
    
        r1=pressure*(eps1-max(delta,delta_min))
        r2=pressure*eps2*vale
        r3=pressure*eps12(elem)*vale
        si1=sigma11(elem)+sigma22(elem)
        si2=sigma11(elem)-sigma22(elem)

        si1=det1*si1+det2*r1
        si2=det1*si2+det2*r2
        sigma12(elem)=det1*sigma12(elem)+det2*r3
        sigma11(elem)=0.5_WP*(si1+si2)
        sigma22(elem)=0.5_WP*(si1-si2)

#if defined (__icepack)
        rdg_conv_elem(elem)  = -min((eps11(elem)+eps22(elem)),0.0_WP)
        rdg_shear_elem(elem) = 0.5_WP*(delta - abs(eps11(elem)+eps22(elem)))
#endif

  end do
 ! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of 
 ! Boullion et al Ocean Modelling 2013, but in an implicit mode:
 ! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
 ! and similarly for si2 and sigma12
end subroutine stress_tensor_m

!
! ==================================================================
! 
subroutine ssh2rhs(mesh)
  ! Compute the contribution from the elevation to the rhs
  ! S.D. 30.07.2014
  use o_param
  use i_param
  use mod_mesh
  use g_config
  use i_arrays
  use g_parsup
  use i_therm_param
  implicit none
  
  integer                  :: row, elem, elnodes(3), n
  real(kind=WP)            :: dx(3), dy(3), vol
  real(kind=WP)            :: val3, meancos, aa, bb, p_ice(3)
  type(t_mesh), intent(in) , target :: mesh
  
#include "associate_mesh.h"

  val3=1.0_WP/3.0_WP
  
  ! use rhs_m and rhs_a for storing the contribution from elevation:
  do row=1, myDim_nod2d 
     rhs_a(row)=0.0_WP
     rhs_m(row)=0.0_WP
  end do
  
  !_____________________________________________________________________________
  ! use floating sea ice for zlevel and zstar
  if (use_floatice .and.  .not. trim(which_ale)=='linfs') then
    do elem=1,myDim_elem2d         
        elnodes=elem2D_nodes(:,elem)
        !_______________________________________________________________________
        ! if element has any cavity node skip it 
        if (ulevels(elem) > 1) cycle
        
        !_______________________________________________________________________
        vol=elem_area(elem)
        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)     
        
        !_______________________________________________________________________
        ! add pressure gradient from sea ice --> in case of floating sea ice
        p_ice=(rhoice*m_ice(elnodes)+rhosno*m_snow(elnodes))*inv_rhowat
        do n=1,3
            p_ice(n)=min(p_ice(n),max_ice_loading)
        end do
        
        !_______________________________________________________________________
        bb=g*val3*vol
        aa=bb*sum(dx*(elevation(elnodes)+p_ice))
        bb=bb*sum(dy*(elevation(elnodes)+p_ice))
        rhs_a(elnodes)=rhs_a(elnodes)-aa    
        rhs_m(elnodes)=rhs_m(elnodes)-bb
    end do
  else
    do elem=1,myDim_elem2d         
        elnodes=elem2D_nodes(:,elem)
        !_______________________________________________________________________
        ! if element has any cavity node skip it 
        if (ulevels(elem) > 1) cycle
        
        vol=elem_area(elem)
        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)     
        bb=g*val3*vol
        aa=bb*sum(dx*elevation(elnodes))
        bb=bb*sum(dy*elevation(elnodes))
        rhs_a(elnodes)=rhs_a(elnodes)-aa   
        rhs_m(elnodes)=rhs_m(elnodes)-bb
    end do
  end if 
end subroutine ssh2rhs
!
!===================================================================
!
subroutine stress2rhs_m(mesh)

  ! add internal stress to the rhs
  ! SD, 30.07.2014
  !-----------------------------------------------------------------  
  use o_param
  use i_param
  use i_therm_param
  use mod_mesh
  use g_config
  use i_arrays
  use g_parsup
  implicit none
  
  integer                  :: k, row, elem, elnodes(3)
  real(kind=WP)            :: dx(3), dy(3), vol
  real(kind=WP)            :: val3, mf, aa, bb
  real(kind=WP)            :: mass, cluster_area, elevation_elem(3)
  type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

  val3=1.0_WP/3.0_WP
  
  do row=1, myDim_nod2d 
     u_rhs_ice(row)=0.0_WP
     v_rhs_ice(row)=0.0_WP
  end do

  do elem=1,myDim_elem2d         
     elnodes=elem2D_nodes(:,elem)
     !_______________________________________________________________________
     ! if element has any cavity node skip it 
     if (ulevels(elem) > 1) cycle

     if(sum(a_ice(elnodes)) < 0.01_WP) cycle !DS
     
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
     !_________________________________________________________________________
     ! if cavity node skip it 
     if ( ulevels_nod2d(row)>1 ) cycle
     
     mass=(m_ice(row)*rhoice+m_snow(row)*rhosno)
     mass=mass/(1.0_WP+mass*mass)
     u_rhs_ice(row)=(u_rhs_ice(row)*mass + rhs_a(row))/area(1,row) 
     v_rhs_ice(row)=(v_rhs_ice(row)*mass + rhs_m(row))/area(1,row) 
  end do
end subroutine stress2rhs_m
!
!===================================================================
!
subroutine EVPdynamics_m(mesh)
  ! assemble rhs and solve for ice velocity
  ! New implementation based on Bouillion et al. Ocean Modelling 2013
  ! SD 30.07.14
  !---------------------------------------------------------

  use o_param
  use i_param
  use i_therm_param
  use mod_mesh
  use g_config
  use i_arrays
  use o_arrays
  use g_parsup
  use g_comm_auto

#if defined (__icepack)
  use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
  use icedrv_main,   only: icepack_to_fesom
#endif

  implicit none
  integer          :: steps, shortstep, i, ed,n
  real(kind=WP)    :: rdt, drag, det
  real(kind=WP)    :: inv_thickness(myDim_nod2D), umod, rhsu, rhsv
  logical          :: ice_el(myDim_elem2D), ice_nod(myDim_nod2D)

!NR for stress_tensor_m
  integer         :: el, elnodes(3)
  real(kind=WP)   :: dx(3), dy(3), msum, asum
  real(kind=WP)   :: eps1, eps2, pressure, pressure_fac(myDim_elem2D), delta
  real(kind=WP)   :: val3, meancos, vale
  real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2

!NR for stress2rhs_m  
  integer        :: k, row
  real(kind=WP)  :: vol
  real(kind=WP)  :: mf,aa, bb,p_ice(3)
  real(kind=WP)  :: mass(myDim_nod2D)
  type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"

  val3=1.0_WP/3.0_WP
  vale=1.0_WP/(ellipse**2)
  det2=1.0_WP/(1.0_WP+alpha_evp)
  det1=alpha_evp*det2
  rdt=ice_dt
  steps=evp_rheol_steps
  
  u_ice_aux=u_ice    ! Initialize solver variables
  v_ice_aux=v_ice

#if defined (__icepack)
  a_ice_old(:)  = a_ice(:)
  m_ice_old(:)  = a_ice(:)
  m_snow_old(:) = m_snow(:)

  call icepack_to_fesom (nx_in=(myDim_nod2D+eDim_nod2D), &
                         aice_out=a_ice,                 &
                         vice_out=m_ice,                 &
                         vsno_out=m_snow)
#endif

!NR inlined, to have all initialization in one place.
!  call ssh2rhs
  
  ! use rhs_m and rhs_a for storing the contribution from elevation:
  do row=1, myDim_nod2d 
     rhs_a(row)=0.0_WP
     rhs_m(row)=0.0_WP
  end do
  
  !_____________________________________________________________________________
  ! use floating sea ice for zlevel and zstar
  if (use_floatice .and.  .not. trim(which_ale)=='linfs') then
    do el=1,myDim_elem2d         
        elnodes=elem2D_nodes(:,el)
        
        !_______________________________________________________________________
        ! if element has any cavity node skip it 
        if (ulevels(el) > 1) cycle
        
        !_______________________________________________________________________
        vol=elem_area(el)
        dx=gradient_sca(1:3,el)
        dy=gradient_sca(4:6,el)     
        
        !_______________________________________________________________________
        ! add pressure gradient from sea ice --> in case of floating sea ice
        p_ice=(rhoice*m_ice(elnodes)+rhosno*m_snow(elnodes))*inv_rhowat
        do n=1,3
            p_ice(n)=min(p_ice(n),max_ice_loading)
        end do
        
        !_______________________________________________________________________
        bb=g*val3*vol
        aa=bb*sum(dx*(elevation(elnodes)+p_ice))
        bb=bb*sum(dy*(elevation(elnodes)+p_ice))
        rhs_a(elnodes)=rhs_a(elnodes)-aa    
        rhs_m(elnodes)=rhs_m(elnodes)-bb
    end do
  !_____________________________________________________________________________
  ! use levitating sea ice for linfs, zlevel and zstar  
  else
    do el=1,myDim_elem2d         
        elnodes=elem2D_nodes(:,el)
        !_______________________________________________________________________
        ! if element has any cavity node skip it 
        if (ulevels(el) > 1)  cycle
        
        vol=elem_area(el)
        dx=gradient_sca(1:3,el)
        dy=gradient_sca(4:6,el)
        bb=g*val3*vol
        aa=bb*sum(dx*elevation(elnodes))
        bb=bb*sum(dy*elevation(elnodes))
        rhs_a(elnodes)=rhs_a(elnodes)-aa    
        rhs_m(elnodes)=rhs_m(elnodes)-bb
    end do
  end if

! precompute thickness (the inverse is needed) and mass (scaled by area)
  do i=1,myDim_nod2D
     inv_thickness(i) = 0._WP
     mass(i) = 0._WP
     ice_nod(i) = .false.
     !_________________________________________________________________________
     ! if cavity ndoe skip it 
     if ( ulevels_nod2d(i)>1 ) cycle
     
     if (a_ice(i) >= 0.01_WP) then
        inv_thickness(i) = (rhoice*m_ice(i)+rhosno*m_snow(i))/a_ice(i)
        inv_thickness(i) = 1.0_WP/max(inv_thickness(i), 9.0_WP)  ! Limit the mass

        mass(i) = (m_ice(i)*rhoice+m_snow(i)*rhosno)
        mass(i) = mass(i)/((1.0_WP+mass(i)*mass(i))*area(1,i))

        ! scale rhs_a, rhs_m, too.
        rhs_a(i) = rhs_a(i)/area(1,i) 
        rhs_m(i) = rhs_m(i)/area(1,i) 

        ice_nod(i) = .true.
     endif
  enddo

! precompute pressure factor
  do el=1,myDim_elem2D
     elnodes=elem2D_nodes(:,el)

     pressure_fac(el) = 0._WP
     ice_el(el) = .false.
     
     !_______________________________________________________________________
     ! if element has any cavity node skip it 
     if (ulevels(el) > 1) cycle
     
     msum=sum(m_ice(elnodes))*val3
     if(msum > 0.01) then
        ice_el(el) = .true.
        asum=sum(a_ice(elnodes))*val3          
        pressure_fac(el) = det2*pstar*msum*exp(-c_pressure*(1.0_WP-asum))
     endif
  end do

  do row=1, myDim_nod2d 
     u_rhs_ice(row)=0.0_WP
     v_rhs_ice(row)=0.0_WP
  end do

!=======================================
! Ice EVPdynamics Iteration main loop:
!=======================================

#if defined (__icepack)
  rdg_conv_elem(:)  = 0.0_WP
  rdg_shear_elem(:) = 0.0_WP
#endif

  do shortstep=1, steps

!NR inlining, to make it easier to have local arrays and fuse loops
!NR    call stress_tensor_m
  ! Internal stress tensor
  ! New implementation following Boullion et al, Ocean Modelling 2013.
  ! SD, 30.07.2014
  !===================================================================
 
   do el=1,myDim_elem2D
     !__________________________________________________________________________
     if (ulevels(el)>1) cycle

     !__________________________________________________________________________
     if(ice_el(el)) then
     
        elnodes=elem2D_nodes(:,el)
        dx=gradient_sca(1:3,el)
        dy=gradient_sca(4:6,el)     
        ! METRICS:
        meancos = val3*metric_factor(el)
        !  
        ! ====== Deformation rate tensor on element elem:
        eps11(el) = sum(dx(:)*u_ice_aux(elnodes)) - sum(v_ice_aux(elnodes))*meancos                !metrics
        eps22(el) = sum(dy(:)*v_ice_aux(elnodes))
        eps12(el) = 0.5_WP*(sum(dy(:)*u_ice_aux(elnodes) + dx(:)*v_ice_aux(elnodes)) &
                         +sum(u_ice_aux(elnodes))*meancos )          !metrics 
        
        ! ======= Switch to eps1,eps2
        eps1 = eps11(el) + eps22(el)
        eps2 = eps11(el) - eps22(el)   
        
        ! ====== moduli:
        delta = sqrt(eps1**2+vale*(eps2**2+4.0_WP*eps12(el)**2))
        
        pressure = pressure_fac(el)/(delta+delta_min)
        
!        si1 = det1*(sigma11(el)+sigma22(el)) + pressure*(eps1-delta) 
!        si2 = det1*(sigma11(el)-sigma22(el)) + pressure*eps2*vale
!        sigma11(el) = 0.5_WP*(si1+si2)
!        sigma22(el) = 0.5_WP*(si1-si2)
!NR directly insert si1, si2 cancels some operations and should increase accuracy
        sigma12(el) = det1*sigma12(el) +       pressure*eps12(el)*vale
        sigma11(el) = det1*sigma11(el) + 0.5_WP*pressure*(eps1 - delta + eps2*vale)
        sigma22(el) = det1*sigma22(el) + 0.5_WP*pressure*(eps1 - delta - eps2*vale)

#if defined (__icepack)
        rdg_conv_elem(el)  = -min((eps11(el)+eps22(el)),0.0_WP)
        rdg_shear_elem(el) = 0.5_WP*(delta - abs(eps11(el)+eps22(el)))
#endif

        !  end do   ! fuse loops
        ! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of 
        ! Boullion et al Ocean Modelling 2013, but in an implicit mode:
        ! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
        ! and similarly for si2 and sigma12

        !NR inlining  call stress2rhs_m
        ! add internal stress to the rhs
        ! SD, 30.07.2014
  !-----------------------------------------------------------------  
        if (elnodes(1) <= myDim_nod2D) then
           u_rhs_ice(elnodes(1)) = u_rhs_ice(elnodes(1)) - elem_area(el)* &
                (sigma11(el)*dx(1)+sigma12(el)*(dy(1) + meancos))                         !metrics 
           v_rhs_ice(elnodes(1)) = v_rhs_ice(elnodes(1)) - elem_area(el)* &
                (sigma12(el)*dx(1)+sigma22(el)*dy(1) - sigma11(el)*meancos)               ! metrics                                              
        end if

        if (elnodes(2) <= myDim_nod2D) then
           u_rhs_ice(elnodes(2)) = u_rhs_ice(elnodes(2)) - elem_area(el)* &
                (sigma11(el)*dx(2)+sigma12(el)*(dy(2) + meancos))                         !metrics 
           v_rhs_ice(elnodes(2)) = v_rhs_ice(elnodes(2)) - elem_area(el)* &
                (sigma12(el)*dx(2)+sigma22(el)*dy(2) - sigma11(el)*meancos)               ! metrics                                              
        end if

        if (elnodes(3) <= myDim_nod2D) then
           u_rhs_ice(elnodes(3)) = u_rhs_ice(elnodes(3)) - elem_area(el)* &
                (sigma11(el)*dx(3)+sigma12(el)*(dy(3) + meancos))                         !metrics 
           v_rhs_ice(elnodes(3)) = v_rhs_ice(elnodes(3)) - elem_area(el)* &
                (sigma12(el)*dx(3)+sigma22(el)*dy(3) - sigma11(el)*meancos)               ! metrics                                              
        end if
     end if
  end do ! --> do el=1,myDim_elem2D
  
  do i=1, myDim_nod2d 
     !__________________________________________________________________________
     if (ulevels_nod2D(i)>1) cycle
     
     !__________________________________________________________________________
     if (ice_nod(i)) then                   ! Skip if ice is absent              

        u_rhs_ice(i) = u_rhs_ice(i)*mass(i) + rhs_a(i)
        v_rhs_ice(i) = v_rhs_ice(i)*mass(i) + rhs_m(i)

 ! end do   !NR fuse loops
 !============= stress2rhs_m ends ======================

 !    do i=1,myDim_nod2D
 
        umod = sqrt((u_ice_aux(i)-u_w(i))**2+(v_ice_aux(i)-v_w(i))**2)
        drag = rdt*Cd_oce_ice*umod*density_0*inv_thickness(i)

        !rhs for water stress, air stress, and u_rhs_ice/v (internal stress + ssh)
        rhsu = u_ice(i)+drag*u_w(i)+rdt*(inv_thickness(i)*stress_atmice_x(i)+u_rhs_ice(i)) + beta_evp*u_ice_aux(i)
        rhsv = v_ice(i)+drag*v_w(i)+rdt*(inv_thickness(i)*stress_atmice_y(i)+v_rhs_ice(i)) + beta_evp*v_ice_aux(i)

        !solve (Coriolis and water stress are treated implicitly)        
        det = bc_index_nod2D(i) / ((1.0_WP+beta_evp+drag)**2 + (rdt*coriolis_node(i))**2)

        u_ice_aux(i) = det*((1.0_WP+beta_evp+drag)*rhsu +rdt*coriolis_node(i)*rhsv)
        v_ice_aux(i) = det*((1.0_WP+beta_evp+drag)*rhsv -rdt*coriolis_node(i)*rhsu)
        end if
     end do ! --> do i=1, myDim_nod2d 

    !___________________________________________________________________________ 
    ! apply sea ice velocity boundary condition 
    do ed=1,myDim_edge2D
        !_______________________________________________________________________
        ! apply coastal sea ice velocity boundary conditions
        if(myList_edge2D(ed) > edge2D_in) then
            u_ice_aux(edges(:,ed))=0.0_WP
            v_ice_aux(edges(:,ed))=0.0_WP
        end if
        
        !_______________________________________________________________________
        ! apply sea ice velocity boundary conditions at cavity-ocean edge
        if (use_cavity) then 
            if ( (ulevels(edge_tri(1,ed))>1) .or. &
                 ( edge_tri(2,ed)>0 .and. ulevels(edge_tri(2,ed))>1) ) then 
                u_ice_aux(edges(1:2,ed))=0.0_WP
                v_ice_aux(edges(1:2,ed))=0.0_WP
            end if 
        end if 
    end do ! --> do ed=1,myDim_edge2D
    
    !___________________________________________________________________________
    call exchange_nod_begin(u_ice_aux, v_ice_aux)

    do row=1, myDim_nod2d 
        u_rhs_ice(row)=0.0_WP
        v_rhs_ice(row)=0.0_WP
    end do

    call exchange_nod_end
    
  end do ! --> do shortstep=1, steps

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
subroutine find_alpha_field_a(mesh)
  ! EVP stability parameter alpha is computed at each element
  ! aEVP implementation
  ! SD, 13.02.2017
  ! ==================================================================
  use o_param
  use i_param
  use i_therm_param
  use mod_mesh
  use g_config
  use i_arrays
  use g_parsup

#if defined (__icepack)
  use icedrv_main,   only: strength
#endif

  implicit none

  integer                  :: elem, elnodes(3)
  real(kind=WP)            :: dx(3), dy(3), msum, asum
  real(kind=WP)            :: eps1, eps2, pressure, delta
  real(kind=WP)            :: val3, meancos, usum, vsum, vale
  type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

  val3=1.0_WP/3.0_WP
  vale=1.0_WP/(ellipse**2)
   do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     !_______________________________________________________________________
     ! if element has any cavity node skip it 
     if (ulevels(elem) > 1) cycle
        
     msum=sum(m_ice(elnodes))*val3
     if(msum<=0.01_WP) cycle !DS
     asum=sum(a_ice(elnodes))*val3
     
     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)     
     ! METRICS:
     vsum=sum(v_ice_aux(elnodes))
     usum=sum(u_ice_aux(elnodes))
     meancos=metric_factor(elem)
     !  
     ! ====== Deformation rate tensor on element elem:
     eps11(elem)=sum(dx*u_ice_aux(elnodes))
     eps11(elem)=eps11(elem)-val3*vsum*meancos                !metrics
     eps22(elem)=sum(dy*v_ice_aux(elnodes))
     eps12(elem)=0.5_WP*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes))
     eps12(elem)=eps12(elem)+0.5_WP*val3*usum*meancos          !metrics 
     
      ! ======= Switch to eps1,eps2
     eps1=eps11(elem)+eps22(elem)
     eps2=eps11(elem)-eps22(elem)   
     
      ! ====== moduli:
     delta=eps1**2+vale*(eps2**2+4.0_WP*eps12(elem)**2)
     delta=sqrt(delta)
         
#if defined (__icepack)
     pressure = sum(strength(elnodes))*val3/(delta+delta_min)/msum
#else
     pressure = pstar*exp(-c_pressure*(1.0_WP-asum))/(delta+delta_min) ! no multiplication
                                                                       ! with thickness (msum)
#endif
     !adjust c_aevp such, that alpha_evp_array and beta_evp_array become in acceptable range
     alpha_evp_array(elem)=max(50.0_WP,sqrt(ice_dt*c_aevp*pressure/rhoice/elem_area(elem)))
     ! /voltriangle(elem) for FESOM1.4
     ! We do not allow alpha to be too small!
   end do
  end subroutine find_alpha_field_a  
! ====================================================================

subroutine stress_tensor_a(mesh)
  ! Internal stress tensor
  ! New implementation following Boullion et al, Ocean Modelling 2013.
  ! and Kimmritz et al., Ocean Modelling 2016
  ! SD, 14.02.2017
  !===================================================================
  use o_param
  use i_param
  use mod_mesh
  use g_config
  use i_arrays
  use g_parsup

#if defined (__icepack)
  use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
#endif

  implicit none

  integer                   :: elem, elnodes(3)
  real(kind=WP)             :: dx(3), dy(3), msum, asum
  real(kind=WP)             :: eps1, eps2, pressure, delta
  real(kind=WP)             :: val3, meancos, usum, vsum, vale
  real(kind=WP)             :: det1, det2, r1, r2, r3, si1, si2
  type(t_mesh), intent(in)  , target :: mesh

#include "associate_mesh.h"
  
  val3=1.0_WP/3.0_WP
  vale=1.0_WP/(ellipse**2)
   do elem=1,myDim_elem2D
     !__________________________________________________________________________
     ! if element has any cavity node skip it 
     if (ulevels(elem) > 1) cycle
     
     !__________________________________________________________________________
     det2=1.0_WP/(1.0_WP+alpha_evp_array(elem))     ! Take alpha from array
     det1=alpha_evp_array(elem)*det2
  
     elnodes=elem2D_nodes(:,elem)
     
     msum=sum(m_ice(elnodes))*val3
     if(msum<=0.01_WP) cycle !DS
     asum=sum(a_ice(elnodes))*val3
     
     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)     
     ! METRICS:
     vsum=sum(v_ice_aux(elnodes))
     usum=sum(u_ice_aux(elnodes))
     meancos=metric_factor(elem)
     !  
     ! ====== Deformation rate tensor on element elem:
     eps11(elem)=sum(dx*u_ice_aux(elnodes))
     eps11(elem)=eps11(elem)-val3*vsum*meancos                !metrics
     eps22(elem)=sum(dy*v_ice_aux(elnodes))
     eps12(elem)=0.5_WP*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes))
     eps12(elem)=eps12(elem)+0.5_WP*val3*usum*meancos          !metrics 
     
      ! ======= Switch to eps1,eps2
     eps1=eps11(elem)+eps22(elem)
     eps2=eps11(elem)-eps22(elem)   
     
      ! ====== moduli:
     delta=eps1**2+vale*(eps2**2+4.0_WP*eps12(elem)**2)
     delta=sqrt(delta)
   
#if defined (__icepack)
     pressure = sum(strength(elnodes))*val3/(delta+delta_min)
#else
     pressure=pstar*msum*exp(-c_pressure*(1.0_WP-asum))/(delta+delta_min)
#endif
    
        r1=pressure*(eps1-delta) 
        r2=pressure*eps2*vale
        r3=pressure*eps12(elem)*vale
        si1=sigma11(elem)+sigma22(elem)
        si2=sigma11(elem)-sigma22(elem)

        si1=det1*si1+det2*r1
        si2=det1*si2+det2*r2
        sigma12(elem)=det1*sigma12(elem)+det2*r3
        sigma11(elem)=0.5_WP*(si1+si2)
        sigma22(elem)=0.5_WP*(si1-si2)

#if defined (__icepack)
        rdg_conv_elem(elem)  = -min((eps11(elem)+eps22(elem)),0.0_WP)
        rdg_shear_elem(elem) = 0.5_WP*(delta - abs(eps11(elem)+eps22(elem)))
#endif

  end do
 ! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of 
 ! Boullion et al Ocean Modelling 2013, but in an implicit mode:
 ! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
 ! and similarly for si2 and sigma12
end subroutine stress_tensor_a
!
!===================================================================
!
subroutine EVPdynamics_a(mesh)
  ! assemble rhs and solve for ice velocity
  ! New implementation based on Bouillion et al. Ocean Modelling 2013
  ! and Kimmritz et al., Ocean Modelling  2016 
  ! SD 14.02.17
  !---------------------------------------------------------

use o_param
use mod_mesh
use i_arrays
USE o_arrays
use i_param
use o_PARAM
use i_therm_param
use g_parsup
use g_config, only: use_cavity
use g_comm_auto
use ice_maEVP_interfaces

#if defined (__icepack)
  use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem
#endif

  implicit none
  integer          :: steps, shortstep, i, ed
  real(kind=WP)    :: rdt, drag, det, fc
  real(kind=WP)    :: thickness, inv_thickness, umod, rhsu, rhsv
  REAL(kind=WP)    :: t0,t1, t2, t3, t4, t5, t00, txx
  type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"

  steps=evp_rheol_steps
  rdt=ice_dt
  u_ice_aux=u_ice    ! Initialize solver variables
  v_ice_aux=v_ice
  call ssh2rhs(mesh)

#if defined (__icepack)
  rdg_conv_elem(:)  = 0.0_WP
  rdg_shear_elem(:) = 0.0_WP
#endif
 
  do shortstep=1, steps 
     call stress_tensor_a(mesh)
     call stress2rhs_m(mesh)    ! _m=_a, so no _m version is the only one!
     do i=1,myDim_nod2D 
     
         !_______________________________________________________________________
         ! if element has any cavity node skip it 
         if (ulevels_nod2d(i)>1) cycle
        
         thickness=(rhoice*m_ice(i)+rhosno*m_snow(i))/max(a_ice(i),0.01_WP)
         thickness=max(thickness, 9.0_WP)   ! Limit if it is too small (0.01 m)
         inv_thickness=1.0_WP/thickness

         umod=sqrt((u_ice_aux(i)-u_w(i))**2+(v_ice_aux(i)-v_w(i))**2)
         drag=rdt*Cd_oce_ice*umod*density_0*inv_thickness

         !rhs for water stress, air stress, and u_rhs_ice/v (internal stress + ssh)
         rhsu=u_ice(i)+drag*u_w(i)+rdt*(inv_thickness*stress_atmice_x(i)+u_rhs_ice(i))
         rhsv=v_ice(i)+drag*v_w(i)+rdt*(inv_thickness*stress_atmice_y(i)+v_rhs_ice(i))

         rhsu=beta_evp_array(i)*u_ice_aux(i)+rhsu
         rhsv=beta_evp_array(i)*v_ice_aux(i)+rhsv
         !solve (Coriolis and water stress are treated implicitly)
         fc=rdt*coriolis_node(i)
         det=(1.0_WP+beta_evp_array(i)+drag)**2+fc**2
         det=bc_index_nod2D(i)/det
         u_ice_aux(i)=det*((1.0_WP+beta_evp_array(i)+drag)*rhsu+fc*rhsv)
         v_ice_aux(i)=det*((1.0_WP+beta_evp_array(i)+drag)*rhsv-fc*rhsu)
     end do
     
    !___________________________________________________________________________ 
    ! apply sea ice velocity boundary condition 
    do ed=1,myDim_edge2D
        !_______________________________________________________________________
        ! apply coastal sea ice velocity boundary conditions
        if(myList_edge2D(ed) > edge2D_in) then
            u_ice_aux(edges(:,ed))=0.0_WP
            v_ice_aux(edges(:,ed))=0.0_WP
        end if
        
        !_______________________________________________________________________
        ! apply sea ice velocity boundary conditions at cavity-ocean edge
        if (use_cavity) then 
            if ( (ulevels(edge_tri(1,ed))>1) .or. &
                 ( edge_tri(2,ed)>0 .and. ulevels(edge_tri(2,ed))>1) ) then 
                u_ice_aux(edges(1:2,ed))=0.0_WP
                v_ice_aux(edges(1:2,ed))=0.0_WP
            end if 
        end if 
    end do ! --> do ed=1,myDim_edge2D
    
     call exchange_nod(u_ice_aux, v_ice_aux)
  end do
     
    u_ice=u_ice_aux
    v_ice=v_ice_aux
 
  call find_alpha_field_a(mesh)             ! alpha_evp_array is initialized with alpha_evp;
                                      ! At this stage we already have non-trivial velocities. 
  call find_beta_field_a(mesh)
end subroutine EVPdynamics_a
!
! =================================================================
!
subroutine find_beta_field_a(mesh)
! beta_evp_array is defined at nodes, and this is the only 
! reason we need it in addition to alpha_evp_array (we work with 
! alpha=beta, and keep different names for generality; mEVP can work with 
! alpha \ne beta, but not aEVP).

use mod_mesh
use o_param
USE i_param
use i_arrays
use g_parsup 
Implicit none
integer :: n

type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"

    DO n=1, myDim_nod2D
       !_______________________________________________________________________
       ! if element has any cavity node skip it 
       if (ulevels_nod2d(n)>1) cycle
       
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
