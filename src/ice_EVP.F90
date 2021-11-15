module ice_EVP_interfaces
    interface
        subroutine stress_tensor(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(in)   , target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
        
        subroutine stress2rhs(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(in)   , target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface  
end module
module ice_EVPdynamics_interfaces
    interface
        subroutine EVPdynamics(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(in)   , target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface  
end module
!
! Contains routines of EVP dynamics
!
!_______________________________________________________________________________
! EVP rheology. The routine computes stress tensor components based on ice 
! velocity field. They are stored as elemental arrays (sigma11, sigma22 and
! sigma12). The ocean velocity is at nodal locations.
subroutine stress_tensor(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use i_param
    USE g_CONFIG
#if defined (__icepack)
    use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
#endif
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=WP)   :: eta, xi, delta, aa
    integer         :: el, elnodes(3)
    real(kind=WP)   :: asum, msum, vale, dx(3), dy(3)
    real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2, dte 
    real(kind=WP)   :: zeta, delta_inv, d1, d2
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: ice_strength
    real(kind=WP), dimension(:), pointer  :: U_ice, V_ice 
    real(kind=WP), dimension(:), pointer  :: eps11, eps12, eps22
    real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    ice_strength => ice%work%ice_strength
    U_ice        => ice%uvice(1, :)
    V_ice        => ice%uvice(2, :)
    eps11        => ice%work%eps11(:)
    eps12        => ice%work%eps12(:)
    eps22        => ice%work%eps22(:)
    sigma11      => ice%work%sigma11(:)
    sigma12      => ice%work%sigma12(:)
    sigma22      => ice%work%sigma22(:)
    
    !___________________________________________________________________________
    vale = 1.0_WP/(ice%ellipse**2)
    dte  = ice%ice_dt/(1.0_WP*ice%evp_rheol_steps)
    det1 = 1.0_WP/(1.0_WP + 0.5_WP*ice%Tevp_inv*dte)
    det2 = 1.0_WP/(1.0_WP + 0.5_WP*ice%Tevp_inv*dte) !*ellipse**2 
     
    !___________________________________________________________________________
    do el=1,myDim_elem2D
        !_______________________________________________________________________
        ! if element contains cavity node skip it 
        if (ulevels(el) > 1) cycle
        
        !_______________________________________________________________________
        ! Check if there is ice on elem
        ! There is no ice in elem 
        ! if (any(m_ice(elnodes)<= 0.) .or. any(a_ice(elnodes) <=0.)) CYCLE     
        if (ice_strength(el) > 0.) then
        
            !___________________________________________________________________
            ! Deformation rate tensor on element elem:
            !du/dx
            eps11(el) = sum(gradient_sca(1:3,el)*U_ice(elem2D_nodes(1:3,el))) &
                        - metric_factor(el) * sum(V_ice(elem2D_nodes(1:3,el)))/3.0_WP
            eps22(el) = sum(gradient_sca(4:6, el)*V_ice(elem2D_nodes(1:3,el)))
            eps12(el) = 0.5_WP*(sum(gradient_sca(4:6,el)*U_ice(elem2D_nodes(1:3,el))) &
                        + sum(gradient_sca(1:3,el)*V_ice(elem2D_nodes(1:3,el))) &
                        + metric_factor(el) * sum(U_ice(elem2D_nodes(1:3,el)))/3.0_WP)
            !___________________________________________________________________
            ! moduli:
            delta = sqrt((eps11(el)*eps11(el) + eps22(el)*eps22(el))*(1.0_WP+vale) + &
                          4.0_WP*vale*eps12(el)*eps12(el) + &
                          2.0_WP*eps11(el)*eps22(el)*(1.0_WP-vale))
            
            !___________________________________________________________________
            ! ===== Here the EVP rheology piece starts
            ! ===== viscosity zeta should exceed zeta_min
            ! (done via limiting delta from above)
            
            !if(delta>pressure/zeta_min) delta=pressure/zeta_min
                !It does not work properly by 
            !creating response where ice_strength is small
                ! Uncomment and test if necessary
            
            ! ===== if delta is too small or zero, viscosity will too large (unlimited)
            ! (limit delta_inv)
            delta_inv = 1.0_WP/max(delta,ice%delta_min) 
            zeta = ice_strength(el)*delta_inv
            
            ! ===== Limiting pressure/Delta  (zeta): it may still happen that pressure/Delta 
            ! is too large in some regions and CFL criterion is violated.
            ! The regularization below was introduced by Hunke, 
            ! but seemingly is not used in the current CICE. 
            ! Without it divergence and zeta can be noisy (but code 
            ! remains stable), using it reduces viscosities too strongly.
            ! It is therefore commented
        
            !if (zeta>Clim_evp*voltriangle(el)) then
            !zeta=Clim_evp*voltriangle(el)
            !end if 
        
            zeta = zeta*ice%Tevp_inv
                            
            r1  = zeta*(eps11(el)+eps22(el)) - ice_strength(el)*ice%Tevp_inv
            r2  = zeta*(eps11(el)-eps22(el))*vale
            r3  = zeta*eps12(el)*vale
            
            si1 = det1*(sigma11(el) + sigma22(el) + dte*r1)
            si2 = det2*(sigma11(el) - sigma22(el) + dte*r2)
            
            sigma12(el) = det2*(sigma12(el)+dte*r3)
            sigma11(el) = 0.5_WP*(si1+si2)
            sigma22(el) = 0.5_WP*(si1-si2)

#if defined (__icepack)
            rdg_conv_elem(el)  = -min((eps11(el)+eps22(el)),0.0_WP)
            rdg_shear_elem(el) = 0.5_WP*(delta - abs(eps11(el)+eps22(el)))
#endif
        endif
    end do
end subroutine stress_tensor
!
!
!_______________________________________________________________________________
! ! subroutine stress_tensor_no1(ice_strength, partit, mesh)
! ! ! EVP rheology. The routine computes stress tensor components based on ice 
! ! ! velocity field. They are stored as elemental arrays (sigma11, sigma22 and
! ! ! sigma12). The ocean velocity is at nodal locations.
! ! use o_param
! ! use i_param
! ! use i_arrays
! ! USE g_CONFIG
! ! USE MOD_MESH
! ! USE MOD_PARTIT
! ! USE MOD_PARSUP
! ! implicit none
! ! type(t_mesh),   intent(in),    target :: mesh
! ! type(t_partit), intent(inout), target :: partit
! ! real(kind=WP), intent(in) :: ice_strength(partit%mydim_elem2D)
! ! real(kind=WP)   :: eta, xi, delta, aa
! ! integer         :: el, elnodes(3)
! ! real(kind=WP)   :: asum, msum, vale, dx(3), dy(3)
! ! real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2, dte 
! ! real(kind=WP)   :: zeta, delta_inv, d1, d2
! ! 
! ! #include "associate_part_def.h"
! ! #include "associate_mesh_def.h"
! ! #include "associate_part_ass.h"
! ! #include "associate_mesh_ass.h"
! ! 
! !   vale = 1.0_WP/(ellipse**2)
! !    
! !   dte  = ice%ice_dt/(1.0_WP*ice%evp_rheol_steps)
! !   det1 = 1.0_WP/(1.0_WP + 0.5_WP*ice%Tevp_inv*dte)
! !   det2 = 1.0_WP/(1.0_WP + 0.5_WP*ice%Tevp_inv*dte) !*ellipse**2 
! !      
! ! 
! !   do el=1,myDim_elem2D
! !      !__________________________________________________________________________
! !      ! if element contains cavity node skip it 
! !      if (ulevels(el) > 1) cycle
! !       ! ===== Check if there is ice on elem
! ! 
! !      ! There is no ice in elem 
! !      ! if (any(m_ice(elnodes)<= 0.) .or. any(a_ice(elnodes) <=0.)) CYCLE     
! !      if (ice_strength(el) > 0.) then
! !       ! =====	
! !       ! ===== Deformation rate tensor on element elem:
! !            !du/dx
! ! 
! !         eps11(el) = sum(mesh%gradient_sca(1:3,el)*U_ice(mesh%elem2D_nodes(1:3,el))) &
! !                -mesh% metric_factor(el) * sum(V_ice(mesh%elem2D_nodes(1:3,el)))/3.0_WP
! ! 
! !         eps22(el) = sum(mesh%gradient_sca(4:6, el)*V_ice(mesh%elem2D_nodes(1:3,el)))
! ! 
! !         eps12(el) = 0.5_WP*(sum(mesh%gradient_sca(4:6,el)*U_ice(mesh%elem2D_nodes(1:3,el))) &
! !                       + sum(mesh%gradient_sca(1:3,el)*V_ice(mesh%elem2D_nodes(1:3,el))) &
! !                        + mesh%metric_factor(el) * sum(U_ice(mesh%elem2D_nodes(1:3,el)))/3.0_WP)
! !         ! ===== moduli:
! !         delta = sqrt((eps11(el)*eps11(el) + eps22(el)*eps22(el))*(1.0_WP+vale) + 4.0_WP*vale*eps12(el)*eps12(el) + &
! !                               2.0_WP*eps11(el)*eps22(el)*(1.0_WP-vale))
! ! 
! !        ! =======================================
! !        ! ===== Here the EVP rheology piece starts
! !        ! =======================================
! ! 
! !       ! ===== viscosity zeta should exceed zeta_min
! !       ! (done via limiting delta from above)
! !       
! !       !if(delta>pressure/zeta_min) delta=pressure/zeta_min
! !            !It does not work properly by 
! ! 	   !creating response where ice_strength is small
! !            ! Uncomment and test if necessary
! !       
! !       ! ===== if delta is too small or zero, viscosity will too large (unlimited)
! !       ! (limit delta_inv)
! !         delta_inv = 1.0_WP/max(delta,delta_min)
! !         
! ! !!PS         delta_inv = delta/(delta+delta_min)
! !         
! !         zeta = ice_strength(el)*delta_inv			     
! !       ! ===== Limiting pressure/Delta  (zeta): it may still happen that pressure/Delta 
! !       ! is too large in some regions and CFL criterion is violated.
! !       ! The regularization below was introduced by Hunke, 
! !       ! but seemingly is not used in the current CICE. 
! !       ! Without it divergence and zeta can be noisy (but code 
! !       ! remains stable), using it reduces viscosities too strongly.
! !       ! It is therefore commented
! !       
! !       !if (zeta>Clim_evp*voltriangle(el)) then
! !       !zeta=Clim_evp*voltriangle(el)
! !       !end if 
! !       
! !         zeta = zeta*Tevp_inv
! !         
! !         r1  = zeta*(eps11(el)+eps22(el)) - ice_strength(el)*ice%Tevp_inv
! !         r2  = zeta*(eps11(el)-eps22(el))*vale
! !         r3  = zeta*eps12(el)*vale
! !         
! !         si1 = det1*(sigma11(el) + sigma22(el) + dte*r1)
! !         si2 = det2*(sigma11(el) - sigma22(el) + dte*r2)
! !         
! !         sigma12(el) = det2*(sigma12(el)+dte*r3)
! !         sigma11(el) = 0.5_WP*(si1+si2)
! !         sigma22(el) = 0.5_WP*(si1-si2)
! !      endif
! !   end do
! ! end subroutine stress_tensor_no1
!
!
!_______________________________________________________________________________
! ! subroutine stress2rhs_e(ice, partit, mesh)
! !     ! EVP implementation:
! !     ! Computes the divergence of stress tensor and puts the result into the
! !     ! rhs vectors. Velocity is at nodes. 
! !     ! The divergence is computed in a cysly over edges. It is slower that the
! !     ! approach in stress2rhs_e inherited from FESOM
! !     USE MOD_ICE
! !     USE MOD_PARTIT
! !     USE MOD_PARSUP
! !     USE MOD_MESH
! !     USE o_PARAM
! !     USE i_PARAM
! !     USE i_therm_param
! !     USE i_arrays
! !     use g_config, only: use_cavity
! !     IMPLICIT NONE
! !     type(t_mesh),   intent(in),    target :: ice
! !     type(t_partit), intent(inout), target :: partit
! !     type(t_mesh),   intent(in),    target :: mesh
! !     !___________________________________________________________________________
! !     integer       :: n, elem, ed, elnodes(3), el(2), ednodes(2)  
! !     real(kind=WP) :: mass, uc, vc,  deltaX1, deltaX2, deltaY1, deltaY2
! !     !___________________________________________________________________________
! !     ! pointer on necessary derived types
! !     real(kind=WP), dimension(:), pointer  :: inv_areamass, elevation
! !     real(kind=WP), dimension(:), pointer  :: U_rhs_ice, V_rhs_ice
! !     real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
! ! #include "associate_part_def.h"
! ! #include "associate_mesh_def.h"
! ! #include "associate_part_ass.h"
! ! #include "associate_mesh_ass.h"
! !     inv_areamass    => ice%work%inv_areamass
! !     U_rhs_ice       => ice%uvice_rhs(1, :)
! !     V_rhs_ice       => ice%uvice_rhs(2, :)
! !     sigma11         => ice%work%sigma11(:)
! !     sigma12         => ice%work%sigma12(:)
! !     sigma22         => ice%work%sigma22(:)
! !     elevation       => ice%srfoce_ssh(:)
! !     !___________________________________________________________________________
! !     DO n=1, myDim_nod2D
! !         U_rhs_ice(n)=0.0_WP
! !         V_rhs_ice(n)=0.0_WP
! !     END DO
! !     
! !     !___________________________________________________________________________
! !     ! Stress divergence
! !     DO  ed=1,myDim_edge2D
! !         ednodes=edges(:,ed) 
! !         el=edge_tri(:,ed)
! !         if(myList_edge2D(ed)>edge2D_in) cycle    
! !         !___________________________________________________________________________
! !         ! stress boundary condition at ocean cavity boundary edge ==0
! !         if (use_cavity) then 
! !             if ( (ulevels(el(1))>1) .or.  ( el(2)>0 .and. ulevels(el(2))>1) ) cycle
! !         end if 
! !         !___________________________________________________________________________
! !         ! elements on both sides
! !         uc = - sigma12(el(1))*edge_cross_dxdy(1,ed) + sigma11(el(1))*edge_cross_dxdy(2,ed) &
! !             + sigma12(el(2))*edge_cross_dxdy(3,ed) - sigma11(el(2))*edge_cross_dxdy(4,ed)
! !         
! !         vc = - sigma22(el(1))*edge_cross_dxdy(1,ed) + sigma12(el(1))*edge_cross_dxdy(2,ed) &
! !             + sigma22(el(2))*edge_cross_dxdy(3,ed) - sigma12(el(2))*edge_cross_dxdy(4,ed)
! !         !___________________________________________________________________________
! !         U_rhs_ice(ednodes(1)) = U_rhs_ice(ednodes(1)) + uc
! !         U_rhs_ice(ednodes(2)) = U_rhs_ice(ednodes(2)) - uc
! !         V_rhs_ice(ednodes(1)) = V_rhs_ice(ednodes(1)) + vc
! !         V_rhs_ice(ednodes(2)) = V_rhs_ice(ednodes(2)) - vc
! !     END DO
! !     
! !     !___________________________________________________________________________
! !     do n=1, myDim_nod2D
! !         !_______________________________________________________________________
! !         ! if cavity node skip it 
! !         if ( ulevels_nod2d(n) > 1 ) cycle
! !         
! !         !_______________________________________________________________________
! !         mass = area(1,n)*(rhoice*m_ice(n)+rhosno*m_snow(n)) 
! !         if(mass > 1.e-3_WP) then 
! !             U_rhs_ice(n) = U_rhs_ice(n) * inv_areamass
! !             V_rhs_ice(n) = V_rhs_ice(n) * inv_areamass
! !         else
! !             U_rhs_ice(n)=0.0_WP
! !             V_rhs_ice(n)=0.0_WP
! !         end if
! !     end do
! !     
! !     !___________________________________________________________________________
! !     ! elevation gradient contribution      
! !     do elem=1,myDim_elem2D
! !         !_______________________________________________________________________
! !         ! if element contains cavity node skip it 
! !         if (ulevels(elem) > 1) cycle
! !         
! !         !_______________________________________________________________________
! !         elnodes=elem2D_nodes(:,elem)
! !         uc=elem_area(elem)*g*sum(gradient_sca(1:3,elem)*elevation(elnodes))/3.0_WP
! !         vc=elem_area(elem)*g*sum(gradient_sca(4:6,elem)*elevation(elnodes))/3.0_WP
! !         U_rhs_ice(elnodes)=U_rhs_ice(elnodes) - uc/area(1,elnodes)
! !         V_rhs_ice(elnodes)=V_rhs_ice(elnodes) - vc/area(1,elnodes)
! !     end do
! ! end subroutine stress2rhs_e
!
!
!_______________________________________________________________________________
! EVP implementation: Computes the divergence of stress tensor and puts the 
! result into the rhs vectors 
subroutine stress2rhs(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE o_PARAM
    USE i_PARAM
    USE i_THERM_PARAM
    USE i_arrays, only:
    IMPLICIT NONE
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: n, el,  k
    real(kind=WP)                         :: val3
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: ice_strength, inv_areamass
    real(kind=WP), dimension(:), pointer  :: U_rhs_ice, V_rhs_ice
    real(kind=WP), dimension(:), pointer  :: rhs_a, rhs_m
    real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    ice_strength    => ice%work%ice_strength
    inv_areamass    => ice%work%inv_areamass
    U_rhs_ice       => ice%uvice_rhs(1, :)
    V_rhs_ice       => ice%uvice_rhs(2, :)
    rhs_a           => ice%data(1)%values_rhs(:)
    rhs_m           => ice%data(2)%values_rhs(:)
    sigma11         => ice%work%sigma11(:)
    sigma12         => ice%work%sigma12(:)
    sigma22         => ice%work%sigma22(:)
    
    !___________________________________________________________________________
    val3=1.0_WP/3.0_WP
    do  n=1, myDim_nod2D
        U_rhs_ice(n)=0.0_WP
        V_rhs_ice(n)=0.0_WP
    end do
    
    !___________________________________________________________________________
    do el=1,myDim_elem2D
        !_______________________________________________________________________
        ! Skip if ice is absent
        ! if (any(m_ice(elnodes)<= 0.) .or. any(a_ice(elnodes) <=0.)) CYCLE 
        !_______________________________________________________________________
        ! if element contains cavity node skip it 
        if (ulevels(el) > 1) cycle
        
        !_______________________________________________________________________
        if (ice_strength(el) > 0._WP) then
            do k=1,3
                U_rhs_ice(elem2D_nodes(k,el)) = U_rhs_ice(elem2D_nodes(k,el)) &
                - elem_area(el) * &
                    (sigma11(el)*gradient_sca(k,el) + sigma12(el)*gradient_sca(k+3,el) &
                    +sigma12(el)*val3*metric_factor(el))            !metrics
                V_rhs_ice(elem2D_nodes(k,el)) = V_rhs_ice(elem2D_nodes(k,el)) & 
                    - elem_area(el) * &
                    (sigma12(el)*gradient_sca(k,el) + sigma22(el)*gradient_sca(k+3,el) &   
                    -sigma11(el)*val3*metric_factor(el))
            end do
        endif
    end do 
    
    !___________________________________________________________________________
    do n=1, myDim_nod2D
        !_______________________________________________________________________
        ! if cavity node skip it 
        if (ulevels_nod2d(n)>1) cycle
        
        !_______________________________________________________________________
        if (inv_areamass(n) > 0._WP) then
            U_rhs_ice(n) = U_rhs_ice(n)*inv_areamass(n) + rhs_a(n)
            V_rhs_ice(n) = V_rhs_ice(n)*inv_areamass(n) + rhs_m(n)
        else
            U_rhs_ice(n) = 0._WP
            V_rhs_ice(n) = 0._WP
        endif
    end do
end subroutine stress2rhs
!
!
!_______________________________________________________________________________
! EVP implementation. Does subcycling and boundary conditions.  
! Velocities at nodes
subroutine EVPdynamics(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE i_therm_param
    USE o_ARRAYS, only: coriolis_node
    USE g_CONFIG, only: max_ice_loading, use_cavity, use_floatice, which_ale, flag_debug
    USE g_comm_auto
    use ice_EVP_interfaces
#if defined (__icepack)
    USE icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
    USE icedrv_main,   only: icepack_to_fesom   
#endif
    IMPLICIT NONE
    type(t_ice)   , intent(in)   , target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer         :: steps, shortstep
    real(kind=WP)   :: rdt, asum, msum, r_a, r_b
    real(kind=WP)   :: drag, det, umod, rhsu, rhsv
    integer         :: n, ed, ednodes(2), el,  elnodes(3)
    real(kind=WP)   :: ax, ay, aa, elevation_dx, elevation_dy
    real(kind=WP)   :: elevation_elem(3), p_ice(3)
    integer         :: use_pice
    real(kind=WP)   :: eta, xi, delta
    integer         :: k
    real(kind=WP)   :: vale, dx(3), dy(3), val3
    real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2, dte 
    real(kind=WP)   :: zeta, delta_inv, d1, d2
    integer         :: elem
    real(kind=WP)   :: mass, uc, vc,  deltaX1, deltaX2, deltaY1, deltaY2
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer :: ice_strength, inv_areamass, inv_mass
    real(kind=WP), dimension(:), pointer :: stress_atmice_x, stress_atmice_y
    real(kind=WP), dimension(:), pointer :: U_ice, V_ice, U_ice_old, V_ice_old, &
                                            U_w, V_w, U_rhs_ice, V_rhs_ice
    real(kind=WP), dimension(:), pointer :: a_ice_old, m_ice_old, m_snow_old
    real(kind=WP), dimension(:), pointer :: rhs_a, rhs_m, elevation
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice           => ice%data(1)%values(:)
    m_ice           => ice%data(2)%values(:)
    m_snow          => ice%data(3)%values(:)
    a_ice_old       => ice%data(1)%values_old(:)
    m_ice_old       => ice%data(2)%values_old(:)
    m_snow_old      => ice%data(3)%values_old(:)
    rhs_a           => ice%data(1)%values_rhs(:)
    rhs_m           => ice%data(2)%values_rhs(:)
    elevation       => ice%srfoce_ssh(:)
    U_ice           => ice%uvice(1, :)
    V_ice           => ice%uvice(2, :)
    U_rhs_ice       => ice%uvice_rhs(1, :)
    V_rhs_ice       => ice%uvice_rhs(2, :)
    U_ice_old       => ice%uvice_old(1, :)
    V_ice_old       => ice%uvice_old(2, :)
    U_w             => ice%srfoce_uv(1, :)
    V_w             => ice%srfoce_uv(2, :)
    stress_atmice_x => ice%stress_atmice_xy(1, :)
    stress_atmice_y => ice%stress_atmice_xy(2, :)
    ice_strength    => ice%work%ice_strength
    inv_areamass    => ice%work%inv_areamass
    inv_mass        => ice%work%inv_mass
    
    !___________________________________________________________________________
    ! If Icepack is used, always update the tracers
#if defined (__icepack)
    a_ice_old(:)  = a_ice(:)
    m_ice_old(:)  = m_ice(:)
    m_snow_old(:) = m_snow(:)
    call icepack_to_fesom (nx_in=(myDim_nod2D+eDim_nod2D), &
                            aice_out=a_ice,                &
                            vice_out=m_ice,                &
                            vsno_out=m_snow)
#endif

    !___________________________________________________________________________
    rdt=ice%ice_dt/(1.0*ice%evp_rheol_steps)
    ax=cos(ice%theta_io)
    ay=sin(ice%theta_io)
    
    ! Precompute values that are never changed during the iteration
    inv_areamass =0.0_WP
    inv_mass     =0.0_WP
    rhs_a        =0.0_WP
    rhs_m        =0.0_WP
    do n=1,myDim_nod2D 
        !_______________________________________________________________________
        ! if cavity node skip it 
        if (ulevels_nod2d(n)>1) cycle
        
        !_______________________________________________________________________
        if ((rhoice*m_ice(n)+rhosno*m_snow(n)) > 1.e-3_WP) then
            inv_areamass(n) = 1._WP/(area(1,n)*(rhoice*m_ice(n)+rhosno*m_snow(n))) 
        else
            inv_areamass(n) = 0._WP
        endif
        
        if (a_ice(n) < 0.01_WP) then
            ! Skip if ice is absent
            inv_mass(n) = 0._WP
        else
            inv_mass(n) = (rhoice*m_ice(n)+rhosno*m_snow(n))/a_ice(n)
            ! Limit the mass if it is too small
            inv_mass(n) = 1.0_WP/max(inv_mass(n), 9.0_WP)        
        endif
        !!PS rhs_a(n)=0.0_WP       ! these are used as temporal storage here
        !!PS rhs_m(n)=0.0_WP       ! for the contribution due to ssh
    enddo
    
    !___________________________________________________________________________
    use_pice=0
    if (use_floatice .and.  .not. trim(which_ale)=='linfs') use_pice=1
    if ( .not. trim(which_ALE)=='linfs') then
        ! for full free surface include pressure from ice mass
        ice_strength=0.0_WP
        do el = 1,myDim_elem2D
            elnodes = elem2D_nodes(:,el)
            !___________________________________________________________________
            ! if element has any cavity node skip it 
            if (ulevels(el) > 1) cycle
            
            !___________________________________________________________________
            if (any(m_ice(elnodes)<=0._WP) .or. any(a_ice(elnodes)<=0._WP)) then
                ! There is no ice in elem
                ice_strength(el) = 0._WP
                
            !___________________________________________________________________
            else
                msum = sum(m_ice(elnodes))/3.0_WP
                asum = sum(a_ice(elnodes))/3.0_WP
                
                !_______________________________________________________________
                ! Hunke and Dukowicz c*h*p*
#if defined (__icepack)
                ice_strength(el) = ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
#else
                ice_strength(el) = ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
#endif
                ice_strength(el) = 0.5_WP*ice_strength(el)
                
                !_______________________________________________________________
                ! use rhs_m and rhs_a for storing the contribution from elevation:
                aa = 9.81_WP*elem_area(el)/3.0_WP
                
                !_______________________________________________________________
                ! add and limit pressure from ice weight in case of floating ice
                ! like in FESOM 1.4
                p_ice=(rhoice*m_ice(elnodes)+rhosno*m_snow(elnodes))*inv_rhowat
                do n=1,3
                    p_ice(n)=min(p_ice(n),max_ice_loading)
                end do
                !!PS p_ice= 0.0_WP
                
                !_______________________________________________________________
                elevation_elem = elevation(elnodes)
                elevation_dx   = sum(gradient_sca(1:3,el)*(elevation_elem+p_ice*use_pice))   
                elevation_dy   = sum(gradient_sca(4:6,el)*(elevation_elem+p_ice*use_pice))
                
                !_______________________________________________________________
                rhs_a(elnodes) = rhs_a(elnodes)-aa*elevation_dx
                rhs_m(elnodes) = rhs_m(elnodes)-aa*elevation_dy
            end if
        enddo
    !___________________________________________________________________________    
    else
        ! for linear free surface
        ice_strength=0.0_WP
        do el = 1,myDim_elem2D
            elnodes = elem2D_nodes(:,el)
            !___________________________________________________________________
            ! if element has any cavity node skip it 
            if (ulevels(el) > 1) cycle
            
            !___________________________________________________________________
            if (any(m_ice(elnodes)<=0._WP) .or. any(a_ice(elnodes)<=0._WP)) then
                ! There is no ice in elem
                ice_strength(el) = 0._WP
            else
                msum = sum(m_ice(elnodes))/3.0_WP
                asum = sum(a_ice(elnodes))/3.0_WP
                
                ! ===== Hunke and Dukowicz c*h*p*
#if defined (__icepack)
                ice_strength(el) = ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
#else
                ice_strength(el) = ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
#endif
                ice_strength(el) = 0.5_WP*ice_strength(el)
                
                ! use rhs_m and rhs_a for storing the contribution from elevation:
                elevation_dx = sum(gradient_sca(1:3,el)*elevation(elnodes))
                elevation_dy = sum(gradient_sca(4:6,el)*elevation(elnodes))
                aa = 9.81_WP*elem_area(el)/3.0_WP
                rhs_a(elnodes) = rhs_a(elnodes)-aa*elevation_dx
                rhs_m(elnodes) = rhs_m(elnodes)-aa*elevation_dy
            end if
        enddo
    endif ! --> if ( .not. trim(which_ALE)=='linfs') then
 
    !___________________________________________________________________________
    do n=1,myDim_nod2D 
        if (ulevels_nod2d(n)>1) cycle
        rhs_a(n) = rhs_a(n)/area(1,n)
        rhs_m(n) = rhs_m(n)/area(1,n)
    enddo
    ! End of Precomputing

    !___________________________________________________________________________
    ! And the ice stepping starts
#if defined (__icepack)
    rdg_conv_elem(:)  = 0.0_WP
    rdg_shear_elem(:) = 0.0_WP
#endif
    do shortstep=1, ice%evp_rheol_steps 
        !_______________________________________________________________________
        !!PS if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call stress_tensor'//achar(27)//'[0m'
        call stress_tensor(ice, partit, mesh)
        !!PS if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call stress2rhs'//achar(27)//'[0m'
        call stress2rhs(ice, partit, mesh) 
        
        !_______________________________________________________________________
        U_ice_old = U_ice !PS
        V_ice_old = V_ice !PS
        do n=1,myDim_nod2D 
            !___________________________________________________________________
            ! if cavity node skip it 
            if ( ulevels_nod2d(n)>1 ) cycle
            
            !___________________________________________________________________
            if (a_ice(n) >= 0.01_WP) then               ! Skip if ice is absent
                umod = sqrt((U_ice(n)-U_w(n))**2+(V_ice(n)-V_w(n))**2)
                drag = ice%cd_oce_ice*umod*density_0*inv_mass(n)
                rhsu = U_ice(n) +rdt*(drag*(ax*U_w(n) - ay*V_w(n))+ &
                       inv_mass(n)*stress_atmice_x(n) + U_rhs_ice(n))
                rhsv = V_ice(n) +rdt*(drag*(ax*V_w(n) + ay*U_w(n))+ &
                       inv_mass(n)*stress_atmice_y(n) + V_rhs_ice(n))
                
                r_a = 1._WP + ax*drag*rdt
                r_b = rdt*(coriolis_node(n) + ay*drag)
                det = 1.0_WP/(r_a*r_a + r_b*r_b)
                
                U_ice(n) = det*(r_a*rhsu +r_b*rhsv)
                V_ice(n) = det*(r_a*rhsv -r_b*rhsu)
            else  ! Set velocities to 0 if ice is absent 
                U_ice(n) = 0.0_WP
                V_ice(n) = 0.0_WP
            end if
        end do
        
        !_______________________________________________________________________
        ! apply sea ice velocity boundary condition 
        DO  ed=1,myDim_edge2D
            !___________________________________________________________________
            ! apply coastal sea ice velocity boundary conditions
            if(myList_edge2D(ed) > edge2D_in) then
                U_ice(edges(1:2,ed))=0.0_WP
                V_ice(edges(1:2,ed))=0.0_WP
            endif
            
            !___________________________________________________________________
            ! apply sea ice velocity boundary conditions at cavity-ocean edge
            if (use_cavity) then 
                if ( (ulevels(edge_tri(1,ed))>1) .or. &
                    ( edge_tri(2,ed)>0 .and. ulevels(edge_tri(2,ed))>1) ) then 
                    U_ice(edges(1:2,ed))=0.0_WP
                    V_ice(edges(1:2,ed))=0.0_WP
                end if 
            end if 
        end do
        call exchange_nod(U_ice, V_ice, partit)   
    end do
end subroutine EVPdynamics
