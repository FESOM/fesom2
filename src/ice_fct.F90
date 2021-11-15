module ice_fct_interfaces
  interface
    subroutine ice_mass_matrix_fill(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine

    subroutine ice_solve_high_order(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine

    subroutine ice_solve_low_order(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine
    
    subroutine ice_fem_fct(tr_array_id, ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        integer   :: tr_array_id
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine
  end interface
end module

module ice_TG_rhs_div_interfaces
    interface
        subroutine ice_TG_rhs_div(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine    
    end interface
end module

module ice_update_for_div_interface
    interface
        subroutine ice_update_for_div(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine    
    end interface
end module

module ice_fct_solve_interface
    interface
        subroutine ice_fct_solve(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine    
    end interface
end module
! ! ! 
! ! ! This file collect subroutines implementing FE-FCT
! ! ! advection scheme by Loehner et al.
! ! ! There is a tunable parameter ice_gamma_fct.
! ! ! Increasing it leads to positivity preserving solution.
! ! 
! ! ! Driving routine is fct_ice_solve. It calles other routines 
! ! ! that do low-order and figh order solutions and then combine them in a flux
! ! ! corrected way. Taylor-Galerkin scheme is used as a high-order one.
! ! 
! ! ! The code is adapted from  FESOM
! ! !
! ! ! =====================================================================
! ! subroutine ice_TG_rhs(partit, mesh)
! !   use MOD_MESH
! !   USE MOD_PARTIT
! !   USE MOD_PARSUP
! !   use i_Arrays
! !   use i_PARAM
! !   use o_PARAM
! !   USE g_CONFIG
! !   implicit none 
! !   real(kind=WP)   :: diff, entries(3),  um, vm, vol, dx(3), dy(3) 
! !   integer         :: n, q, row, elem, elnodes(3)
! !   type(t_partit), intent(inout), target :: partit
! !   type(t_mesh),   intent(in),    target :: mesh
! ! 
! ! #include "associate_part_def.h"
! ! #include "associate_mesh_def.h"
! ! #include "associate_part_ass.h"
! ! #include "associate_mesh_ass.h"
! ! 
! !     ! Taylor-Galerkin (Lax-Wendroff) rhs
! !     DO row=1, myDim_nod2D
! !         rhs_m(row)=0._WP
! !         rhs_a(row)=0._WP
! !         rhs_ms(row)=0._WP        
! ! #if defined (__oifs)
! !         ths_temp(row)=0._WP
! ! #endif /* (__oifs) */
! !     END DO
! !     
! !     ! Velocities at nodes
! !     do elem=1,myDim_elem2D          !assembling rhs over elements
! !         elnodes=elem2D_nodes(:,elem)
! !         !_______________________________________________________________________
! !         ! if cavity element skip it 
! !         if (ulevels(elem)>1) cycle
! !         
! !         !derivatives
! !         dx=gradient_sca(1:3,elem)
! !         dy=gradient_sca(4:6,elem)
! !         vol=elem_area(elem)
! !         !um=sum(U_ice(elnodes))/3.0_WP
! !         !vm=sum(V_ice(elnodes))/3.0_WP
! !         um=sum(U_ice(elnodes))
! !         vm=sum(V_ice(elnodes))
! !      
! !         !diffusivity
! !         diff=ice_diff*sqrt(elem_area(elem)/scale_area)
! !         DO n=1,3
! !             row=elnodes(n)
! !             DO q = 1,3 
! !                 !entries(q)= vol*dt*((dx(n)*um+dy(n)*vm)/3.0_WP - &
! !                 !            diff*(dx(n)*dx(q)+ dy(n)*dy(q))- &
! !                 !	       0.5*dt*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))) 
! !                 entries(q)= vol*ice%ice_dt*((dx(n)*(um+u_ice(elnodes(q)))+ &
! !                             dy(n)*(vm+v_ice(elnodes(q))))/12.0_WP - &
! !                             diff*(dx(n)*dx(q)+ dy(n)*dy(q))- &
! !                             0.5_WP*ice%ice_dt*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0_WP)    
! !             END DO
! !             rhs_m(row)=rhs_m(row)+sum(entries*m_ice(elnodes))
! !             rhs_a(row)=rhs_a(row)+sum(entries*a_ice(elnodes))
! !             rhs_ms(row)=rhs_ms(row)+sum(entries*m_snow(elnodes))
! ! #if defined (__oifs)
! !             rhs_temp(row)=rhs_temp(row)+sum(entries*ice_temp(elnodes))
! ! #endif /* (__oifs) */
! !         END DO
! !     end do
! ! end subroutine ice_TG_rhs   
!
!----------------------------------------------------------------------------
!
! subroutine ice_fct_init(partit, mesh)
!   use o_PARAM
!   use MOD_MESH
!   USE MOD_PARTIT
!   USE MOD_PARSUP
!   use i_ARRAYS
!   use ice_fct_interfaces
!   implicit none
!   integer   :: n_size
!   type(t_partit), intent(inout), target :: partit
!   type(t_mesh),   intent(in),    target :: mesh
! 
! #include "associate_part_def.h"
! #include "associate_mesh_def.h"
! #include "associate_part_ass.h"
! #include "associate_mesh_ass.h"
! 
!   
!   n_size=myDim_nod2D+eDim_nod2D
!   
!   ! Initialization of arrays necessary to implement FCT algorithm
!   allocate(m_icel(n_size), a_icel(n_size), m_snowl(n_size))  ! low-order solutions
!   m_icel=0.0_WP
!   a_icel=0.0_WP 
!   m_snowl=0.0_WP
! #if defined (__oifs)
!   allocate(m_templ(n_size))  
!   allocate(dm_temp(n_size))  
! #endif /* (__oifs) */
!   allocate(icefluxes(myDim_elem2D,3))
!   allocate(icepplus(n_size), icepminus(n_size))
!   icefluxes = 0.0_WP
!   icepplus = 0.0_WP
!   icepminus= 0.0_WP
!   
! #if defined (__oifs)
!   m_templ=0.0_WP
!   dm_temp=0.0_WP
! #endif /* (__oifs) */
!   
!   allocate(dm_ice(n_size), da_ice(n_size), dm_snow(n_size))  ! increments of high
!   dm_ice = 0.0_WP                                            ! order solutions
!   da_ice = 0.0_WP
!   dm_snow = 0.0_WP
!   
!   
! end subroutine ice_fct_init
!
!
!_______________________________________________________________________________
subroutine ice_fct_solve(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use g_config, only: flag_debug
    use ice_fct_interfaces
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in),    target :: mesh
    !___________________________________________________________________________
    ! pointer on necessary derived types
    
    !___________________________________________________________________________
    ! Driving routine
    ! uses arrays of low-order solutions as temp storage. It should preceed the 
    ! call of low order solution.  
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call ice_solve_high_order...'//achar(27)//'[0m' 
    call ice_solve_high_order(ice, partit, mesh)   
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call ice_solve_low_order...'//achar(27)//'[0m' 
    call ice_solve_low_order(ice, partit, mesh)
    
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call ice_fem_fct...'//achar(27)//'[0m' 
    call ice_fem_fct(1, ice, partit, mesh)    ! m_ice
    call ice_fem_fct(2, ice, partit, mesh)    ! a_ice
    call ice_fem_fct(3, ice, partit, mesh)    ! m_snow
#if defined (__oifs)
    call ice_fem_fct(4, ice, partit, mesh)    ! ice_temp
#endif /* (__oifs) */

end subroutine ice_fct_solve
!
!
!_______________________________________________________________________________
!============================
! Low-order solution
!============================
!
! It is assumed that m_ice, a_ice and m_snow from the previous time step 
! are known at 1:myDim_nod2D+eDim_nod2D.
! We add diffusive contribution to the rhs. The diffusion operator
! is implemented as the difference between the consistent and lumped mass
! matrices acting on the field from the previous time step. The consistent 
! mass matrix on the lhs is replaced with the lumped one.       
subroutine ice_solve_low_order(ice, partit, mesh) 
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use i_PARAM
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer       :: row, clo, clo2, cn, location(100)
    real(kind=WP) :: gamma
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer :: a_icel, m_icel, m_snowl
    real(kind=WP), dimension(:), pointer :: rhs_a, rhs_m, rhs_ms
    real(kind=WP), dimension(:), pointer :: mass_matrix
#if defined (__oifs)
    real(kind=WP), dimension(:), pointer :: ice_temp, rhs_temp, m_templ
#endif /* (__oifs) */    
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice       => ice%data(1)%values(:)
    a_icel      => ice%data(1)%valuesl(:)
    rhs_a       => ice%data(1)%values_rhs(:)
    m_ice       => ice%data(2)%values(:)
    m_icel      => ice%data(2)%valuesl(:)
    rhs_m       => ice%data(2)%values_rhs(:)
    m_snow      => ice%data(3)%values(:)
    m_snowl     => ice%data(3)%valuesl(:)
    rhs_ms      => ice%data(3)%values_rhs(:)
#if defined (__oifs)
    ice_temp    => ice%data(4)%values(:)
    m_templ     => ice%data(4)%valuesl(:)
    rhs_temp    => ice%data(4)%values_rhs(:)
#endif /* (__oifs) */ 
    mass_matrix => ice%work%fct_massmatrix(:)
    
    !___________________________________________________________________________
    ! Added diffusivity parameter Adjust it to ensure posivity of solution    
    gamma=ice%ice_gamma_fct         
    
    !___________________________________________________________________________
    do row=1,myDim_nod2D
        !_______________________________________________________________________
        ! if there is cavity no ice fxt low order
        if (ulevels_nod2D(row)>1) cycle
        
        !_______________________________________________________________________
        clo=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)+1
        clo2=ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(1)
        cn=clo2-clo+1
        location(1:cn)=nn_pos(1:cn, row)
        m_icel(row)=(rhs_m(row)+gamma*sum(mass_matrix(clo:clo2)* &
                    m_ice(location(1:cn))))/area(1,row) + &
                    (1.0_WP-gamma)*m_ice(row)
        a_icel(row)=(rhs_a(row)+gamma*sum(mass_matrix(clo:clo2)* &
                    a_ice(location(1:cn))))/area(1,row) + &
                    (1.0_WP-gamma)*a_ice(row)
        m_snowl(row)=(rhs_ms(row)+gamma*sum(mass_matrix(clo:clo2)* &
                    m_snow(location(1:cn))))/area(1,row) + &
                    (1.0_WP-gamma)*m_snow(row)
#if defined (__oifs)
        m_templ(row)=(rhs_temp(row)+gamma*sum(mass_matrix(clo:clo2)* &
                  ice_temp(location(1:cn))))/area(1,row) + &
                  (1.0_WP-gamma)*ice_temp(row)
#endif /* (__oifs) */
    end do
    ! Low-order solution must be known to neighbours
    call exchange_nod(m_icel,a_icel,m_snowl, partit)
#if defined (__oifs)
    call exchange_nod(m_templ, partit)
#endif /* (__oifs) */
end subroutine ice_solve_low_order     
!
!
!_______________________________________________________________________________
subroutine ice_solve_high_order(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: n,i,clo,clo2,cn,location(100),row
    real(kind=WP)                         :: rhs_new
    integer                               :: num_iter_solve=3
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer :: a_icel, m_icel, m_snowl
    real(kind=WP), dimension(:), pointer :: rhs_a, rhs_m, rhs_ms
    real(kind=WP), dimension(:), pointer :: da_ice, dm_ice, dm_snow
    real(kind=WP), dimension(:), pointer :: mass_matrix
#if defined (__oifs)
    real(kind=WP), dimension(:), pointer :: dm_temp, rhs_temp, m_templ
#endif /* (__oifs) */
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_icel      => ice%data(1)%valuesl(:)
    rhs_a       => ice%data(1)%values_rhs(:)
    da_ice      => ice%data(1)%dvalues(:)
    m_icel      => ice%data(2)%valuesl(:)
    rhs_m       => ice%data(2)%values_rhs(:)
    dm_ice      => ice%data(2)%dvalues(:)
    m_snowl     => ice%data(3)%valuesl(:)
    rhs_ms      => ice%data(3)%values_rhs(:)
    dm_snow     => ice%data(3)%dvalues(:)
#if defined (__oifs)
    m_templ     => ice%data(4)%valuesl(:)
    rhs_temp    => ice%data(4)%values_rhs(:)
    dm_temp     => ice%data(4)%dvalues(:)
#endif /* (__oifs) */ 
    mass_matrix => ice%work%fct_massmatrix(:)

    !___________________________________________________________________________
    ! Does Taylor-Galerkin solution
    ! --> the first approximation
    do row=1,myDim_nod2D
        !_______________________________________________________________________
        ! if cavity node skip it 
        if (ulevels_nod2d(row)>1) cycle
        
        dm_ice(row)=rhs_m(row)/area(1,row)
        da_ice(row)=rhs_a(row)/area(1,row)
        dm_snow(row)=rhs_ms(row)/area(1,row)
#if defined (__oifs)
        dm_temp(row)=rhs_temp(row)/area(1,row)
#endif /* (__oifs) */
    end do
    call exchange_nod(dm_ice, da_ice, dm_snow, partit)
#if defined (__oifs)
    call exchange_nod(dm_temp, partit)
#endif /* (__oifs) */
    
    !___________________________________________________________________________
    ! --> iterate 
    do n=1,num_iter_solve-1
        do row=1,myDim_nod2D
            ! if cavity node skip it 
            if (ulevels_nod2d(row)>1) cycle
            
            clo=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)+1
            clo2=ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(1)
            cn=clo2-clo+1
            location(1:cn)=nn_pos(1:cn,row)
            rhs_new     = rhs_m(row) - sum(mass_matrix(clo:clo2)*dm_ice(location(1:cn)))
            m_icel(row) = dm_ice(row)+rhs_new/area(1,row)
            rhs_new     = rhs_a(row) - sum(mass_matrix(clo:clo2)*da_ice(location(1:cn)))
            a_icel(row) = da_ice(row)+rhs_new/area(1,row)
            rhs_new     = rhs_ms(row) - sum(mass_matrix(clo:clo2)*dm_snow(location(1:cn)))
            m_snowl(row)= dm_snow(row)+rhs_new/area(1,row) 
#if defined (__oifs)
            rhs_new     = rhs_temp(row) - sum(mass_matrix(clo:clo2)*dm_temp(location(1:cn)))
            m_templ(row)= dm_temp(row)+rhs_new/area(1,row)
#endif /* (__oifs) */
        end do
        do row=1,myDim_nod2D
            ! if cavity node skip it 
            if (ulevels_nod2d(row)>1) cycle
            dm_ice(row)=m_icel(row)
            da_ice(row)=a_icel(row)
            dm_snow(row)=m_snowl(row)
#if defined (__oifs)
            dm_temp(row)=m_templ(row)
#endif /* (__oifs) */
        end do
        call exchange_nod(dm_ice, da_ice, dm_snow, partit)
#if defined (__oifs)
        call exchange_nod(dm_temp, partit)
#endif /* (__oifs) */
    end do
end subroutine ice_solve_high_order
!
!
!_______________________________________________________________________________
! Flux corrected transport algorithm for tracer advection
!
! It is based on Loehner et al. (Finite-element flux-corrected 
! transport (FEM-FCT) for the Euler and Navier-Stokes equation, 
! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
! Turek. (kuzmin@math.uni-dortmund.de) 
subroutine ice_fem_fct(tr_array_id, ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use i_param
    use o_PARAM
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in),    target :: mesh
    !___________________________________________________________________________
    integer   :: tr_array_id
    integer   :: icoef(3,3),n,q, elem,elnodes(3),row
    real(kind=WP)   :: vol, flux, ae, gamma
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:)  , pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:)  , pointer :: a_icel, m_icel, m_snowl
    real(kind=WP), dimension(:)  , pointer :: da_ice, dm_ice, dm_snow
    real(kind=WP), dimension(:)  , pointer :: tmax, tmin, icepplus, icepminus
    real(kind=WP), dimension(:,:), pointer :: icefluxes
#if defined (__oifs)
    real(kind=WP), dimension(:)  , pointer :: dm_temp, ice_temp, m_templ
#endif /* (__oifs) */
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice       => ice%data(1)%values(:)
    a_icel      => ice%data(1)%valuesl(:)
    da_ice      => ice%data(1)%dvalues(:)
    m_ice       => ice%data(2)%values(:)
    m_icel      => ice%data(2)%valuesl(:)
    dm_ice      => ice%data(2)%dvalues(:)
    m_snow      => ice%data(3)%values(:)
    m_snowl     => ice%data(3)%valuesl(:)
    dm_snow     => ice%data(3)%dvalues(:)
#if defined (__oifs)
    ice_temp    => ice%data(4)%values(:)
    m_templ     => ice%data(4)%valuesl(:)
    dm_temp     => ice%data(4)%dvalues(:)
#endif /* (__oifs) */ 
    tmax        => ice%work%fct_tmax
    tmin        => ice%work%fct_tmin
    icepplus    => ice%work%fct_plus
    icepminus   => ice%work%fct_minus
    icefluxes   => ice%work%fct_fluxes
    
    !___________________________________________________________________________
    ! It should coinside with gamma in ts_solve_low_order
    gamma=ice%ice_gamma_fct          
  
    !___________________________________________________________________________
    ! Compute elemental antidiffusive fluxes to nodes
    ! This is the most unpleasant part --- 
    ! it takes memory and time. For every element 
    ! we need its antidiffusive contribution to 
    ! each of its 3 nodes
    tmax = 0.0_WP
    tmin = 0.0_WP
    
    ! Auxiliary elemental operator (mass matrix- lumped mass matrix)
    icoef=1
    do n=1,3   ! three upper nodes
        ! Cycle over rows  row=elnodes(n)
        icoef(n,n)=-2
    end do	    
    
    !___________________________________________________________________________
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        elnodes=elem2D_nodes(:,elem)
        
        !_______________________________________________________________________
        ! if cavity cycle over
        if(ulevels(elem)>1) cycle !LK89140
        
        !_______________________________________________________________________
        vol=elem_area(elem)
        if (tr_array_id==1) then
            do q=1,3       
            icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*m_ice(elnodes) + &
                            dm_ice(elnodes)))*(vol/area(1,elnodes(q)))/12.0_WP
            end do
        end if
        
        if (tr_array_id==2) then
            do q=1,3       
            icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*a_ice(elnodes) + &
                            da_ice(elnodes)))*(vol/area(1,elnodes(q)))/12.0_WP
            end do
        end if
        
        if (tr_array_id==3) then
            do q=1,3       
            icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*m_snow(elnodes) + &
                            dm_snow(elnodes)))*(vol/area(1,elnodes(q)))/12.0_WP
            end do
        end if
        
#if defined (__oifs)
        if (tr_array_id==4) then
            do q=1,3
            icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*ice_temp(elnodes) + &
                            dm_temp(elnodes)))*(vol/area(1,elnodes(q)))/12.0_WP
            end do
        end if
#endif /* (__oifs) */
    end do
     
    !___________________________________________________________________________
    ! Screening the low-order solution
    ! TO BE ADDED IF FOUND NECESSARY
    ! Screening means comparing low-order solutions with the
    ! solution on the previous time step and using whichever 
    ! is greater/smaller in computations of max/min below
    
    !___________________________________________________________________________
    ! Cluster min/max
    if (tr_array_id==1) then
        do row=1, myDim_nod2D
            if (ulevels_nod2d(row)>1) cycle
            n=nn_num(row)
            tmax(row)=maxval(m_icel(nn_pos(1:n,row)))
            tmin(row)=minval(m_icel(nn_pos(1:n,row)))
                ! Admissible increments
            tmax(row)=tmax(row)-m_icel(row)
            tmin(row)=tmin(row)-m_icel(row)
        end do
    end if
    
    if (tr_array_id==2) then
        do row=1, myDim_nod2D
            if (ulevels_nod2d(row)>1) cycle
            n=nn_num(row)
            tmax(row)=maxval(a_icel(nn_pos(1:n,row)))
            tmin(row)=minval(a_icel(nn_pos(1:n,row)))
                ! Admissible increments
            tmax(row)=tmax(row)-a_icel(row)
            tmin(row)=tmin(row)-a_icel(row)
        end do
    end if
 
    if (tr_array_id==3) then
        do row=1, myDim_nod2D
            if (ulevels_nod2d(row)>1) cycle
            n=nn_num(row)
            tmax(row)=maxval(m_snowl(nn_pos(1:n,row)))
            tmin(row)=minval(m_snowl(nn_pos(1:n,row)))
                ! Admissible increments
            tmax(row)=tmax(row)-m_snowl(row)
            tmin(row)=tmin(row)-m_snowl(row)
        end do
    end if

#if defined (__oifs)
    if (tr_array_id==4) then
        do row=1, myDim_nod2D
            if (ulevels_nod2d(row)>1) cycle
            n=nn_num(row)
            tmax(row)=maxval(m_templ(nn_pos(1:n,row)))
            tmin(row)=minval(m_templ(nn_pos(1:n,row)))
            ! Admissible increments
            tmax(row)=tmax(row)-m_templ(row)
            tmin(row)=tmin(row)-m_templ(row)
        end do
    end if
#endif /* (__oifs) */
 
    !___________________________________________________________________________
    ! Sums of positive/negative fluxes to node row
    icepplus=0._WP
    icepminus=0._WP
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        elnodes=elem2D_nodes(:,elem)
        
        !_______________________________________________________________________
        ! if cavity cycle over
        if(ulevels(elem)>1) cycle !LK89140
        
        !_______________________________________________________________________
        do q=1,3
            n=elnodes(q) 
            flux=icefluxes(elem,q)
            if (flux>0) then
                icepplus(n)=icepplus(n)+flux
            else
                icepminus(n)=icepminus(n)+flux	  
            end if
        end do  
    end do   
        
    !___________________________________________________________________________
    ! The least upper bound for the correction factors
    do n=1,myDim_nod2D
        !_______________________________________________________________________
        ! if cavity cycle over
        if(ulevels_nod2D(n)>1) cycle !LK89140
        
        flux=icepplus(n)
        if (abs(flux)>0) then
            icepplus(n)=min(1.0_WP,tmax(n)/flux)
        else
            icepplus(n)=0._WP
        end if
        
        flux=icepminus(n)
        if (abs(flux)>0) then
            icepminus(n)=min(1.0_WP,tmin(n)/flux)
        else
            icepminus(n)=0._WP
        end if
    end do
    ! pminus and pplus are to be known to neighbouting PE
    call exchange_nod(icepminus, icepplus, partit)
    
    !___________________________________________________________________________
    ! Limiting
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        elnodes=elem2D_nodes(:,elem)
        
        !_______________________________________________________________________
        ! if cavity cycle over
        if(ulevels(elem)>1) cycle !LK89140
        
        !_______________________________________________________________________
        ae=1.0_WP
        do q=1,3
            n=elnodes(q)  
            flux=icefluxes(elem,q)
            if(flux>=0._WP) ae=min(ae,icepplus(n))
            if(flux<0._WP) ae=min(ae,icepminus(n))
        end do
        icefluxes(elem,:)=ae*icefluxes(elem,:)
    end do   
  
    !_______________________________________________________________________
    ! Update the solution 
    if(tr_array_id==1) then
        do n=1,myDim_nod2D
            if(ulevels_nod2D(n)>1) cycle !LK89140
            m_ice(n)=m_icel(n)
        end do      
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            
            !___________________________________________________________________
            ! if cavity cycle over
            if(ulevels(elem)>1) cycle !LK89140
            
            do q=1,3
                n=elnodes(q)  
                m_ice(n)=m_ice(n)+icefluxes(elem,q)
            end do
        end do   
    end if
    
    if(tr_array_id==2) then
        do n=1,myDim_nod2D
            if(ulevels_nod2D(n)>1) cycle !LK89140
            a_ice(n)=a_icel(n)
        end do      
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            
            !___________________________________________________________________
            ! if cavity cycle over
            if(ulevels(elem)>1) cycle !LK89140
            
            do q=1,3
                n=elnodes(q)  
                a_ice(n)=a_ice(n)+icefluxes(elem,q)
            end do
        end do   
    end if
    
    if(tr_array_id==3) then
        do n=1,myDim_nod2D
            if(ulevels_nod2D(n)>1) cycle !LK89140
            m_snow(n)=m_snowl(n)
        end do      
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            
            !___________________________________________________________________
            ! if cavity cycle over
            if(ulevels(elem)>1) cycle !LK89140
            
            do q=1,3
                n=elnodes(q)  
                m_snow(n)=m_snow(n)+icefluxes(elem,q)
            end do
        end do   
    end if
    
#if defined (__oifs)
    if(tr_array_id==4) then
        do n=1,myDim_nod2D
            if(ulevels_nod2D(n)>1) cycle !LK89140
            ice_temp(n)=m_templ(n)
        end do
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            !___________________________________________________________________
            ! if cavity cycle over
            if(ulevels(elem)>1) cycle !LK89140
            
            do q=1,3
            n=elnodes(q)
            ice_temp(n)=ice_temp(n)+icefluxes(elem,q)
            end do
        end do
    end if
#endif /* (__oifs) */
    call exchange_nod(m_ice, a_ice, m_snow, partit)
#if defined (__oifs)
    call exchange_nod(ice_temp, partit)
#endif /* (__oifs) */    
end subroutine ice_fem_fct
!
!
!_______________________________________________________________________________
subroutine ice_mass_matrix_fill(ice, partit, mesh)
! Used in ice_fct inherited from FESOM
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use i_PARAM
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                             :: n, n1, n2, row
    integer                             :: elem, elnodes(3), q, offset, col, ipos 
    integer, allocatable                :: col_pos(:)
    real(kind=WP)                       :: aa
    integer                             :: flag=0,iflag=0
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer :: mass_matrix
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    mass_matrix => ice%work%fct_massmatrix(:)
    
    !___________________________________________________________________________
    ! a)
    allocate(col_pos(myDim_nod2D+eDim_nod2D))
    
    DO elem=1,myDim_elem2D
        elnodes=elem2D_nodes(:,elem) 
        
        !_______________________________________________________________________
        do n=1,3
            row=elnodes(n)
            if(row>myDim_nod2D) cycle
            !___________________________________________________________________
            ! Global-to-local neighbourhood correspondence  
            DO q=1,nn_num(row)
                col_pos(nn_pos(q,row))=q
            END DO 
            offset=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)
            DO q=1,3 
                col=elnodes(q)
                !_______________________________________________________________
                ! if element is cavity cycle over
                if(ulevels(elem)>1) cycle
                
                ipos=offset+col_pos(col)
                mass_matrix(ipos)=mass_matrix(ipos)+elem_area(elem)/12.0_WP
                if(q==n) then                     
                    mass_matrix(ipos)=mass_matrix(ipos)+elem_area(elem)/12.0_WP
                end if
            END DO
        end do
    END DO
  
    ! TEST: area==sum of row entries in mass_matrix:
    DO q=1,myDim_nod2D
        !___________________________________________________________________
        ! if cavity cycle over
        if(ulevels_nod2d(q)>1) cycle
        
        !_______________________________________________________________________
        offset=ssh_stiff%rowptr(q)-ssh_stiff%rowptr(1)+1
        n=ssh_stiff%rowptr(q+1)-ssh_stiff%rowptr(1)
        aa=sum(mass_matrix(offset:n))  
        !!PS if(abs(area(1,q)-aa)>.1_WP) then
        if(abs(area(ulevels_nod2d(q),q)-aa)>.1_WP) then
            iflag=q
            flag=1
        endif
    END DO
    if(flag>0) then
        offset=ssh_stiff%rowptr(iflag)-ssh_stiff%rowptr(1)+1
        n=ssh_stiff%rowptr(iflag+1)-ssh_stiff%rowptr(1)
        aa=sum(mass_matrix(offset:n))
        write(*,*) '#### MASS MATRIX PROBLEM', mype, iflag, aa, area(1,iflag), ulevels_nod2D(iflag)
    endif
    deallocate(col_pos)
end subroutine ice_mass_matrix_fill
!
!
!_______________________________________________________________________________
subroutine ice_TG_rhs_div(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use i_PARAM
    use o_PARAM
    USE g_CONFIG
    implicit none 
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=WP)            :: diff, entries(3),  um, vm, vol, dx(3), dy(3) 
    integer                  :: n, q, row, elem, elnodes(3)
    real(kind=WP)            :: c1, c2, c3, c4, cx1, cx2, cx3, cx4, entries2(3) 
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer :: rhs_a, rhs_m, rhs_ms 
    real(kind=WP), dimension(:), pointer :: rhs_adiv, rhs_mdiv, rhs_msdiv
    real(kind=WP), dimension(:), pointer :: u_ice, v_ice
#if defined (__oifs)
    real(kind=WP), dimension(:), pointer :: rhs_temp, rhs_tempdiv
#endif /* (__oifs) */
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice       => ice%data(1)%values(:)
    m_ice       => ice%data(2)%values(:)
    m_snow      => ice%data(3)%values(:)
    u_ice       => ice%uvice(1, :)
    v_ice       => ice%uvice(2, :)
    rhs_a       => ice%data(1)%values_rhs(:)
    rhs_m       => ice%data(2)%values_rhs(:)
    rhs_ms      => ice%data(3)%values_rhs(:)
    rhs_adiv    => ice%data(1)%values_div_rhs(:)
    rhs_mdiv    => ice%data(2)%values_div_rhs(:)
    rhs_msdiv   => ice%data(3)%values_div_rhs(:)
#if defined (__oifs)
    rhs_temp    => ice%data(4)%values_rhs(:)
    rhs_tempdiv => ice%data(4)%values_div_rhs(:)
#endif /* (__oifs) */    
    
    !___________________________________________________________________________
    ! Computes the rhs in a Taylor-Galerkin way (with upwind type of 
    ! correction for the advection operator)
    ! In this version I tr to split divergent term off, so that FCT works without it.
    do row=1, myDim_nod2D
        rhs_m(row)      = 0.0_WP
        rhs_a(row)      = 0.0_WP
        rhs_ms(row)     = 0.0_WP
#if defined (__oifs)
        rhs_temp(row)   = 0.0_WP
#endif /* (__oifs) */
        rhs_mdiv(row)   = 0.0_WP
        rhs_adiv(row)   = 0.0_WP
        rhs_msdiv(row)  = 0.0_WP
#if defined (__oifs)
        rhs_tempdiv(row)= 0.0_WP        
#endif /* (__oifs) */
    end do
    
    !___________________________________________________________________________
    do elem=1,myDim_elem2D          !assembling rhs over elements
        elnodes=elem2D_nodes(:,elem)
        !_______________________________________________________________________
        ! if cavity element skip it 
        if (ulevels(elem)>1) cycle
        
        !derivatives
        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)
        vol=elem_area(elem)
        um=sum(u_ice(elnodes))
        vm=sum(v_ice(elnodes))
        !_______________________________________________________________________
        ! this is exact computation (no assumption of u=const on elements used 
        ! in the standard version)
        c1=(um*um+sum(u_ice(elnodes)*u_ice(elnodes)))/12.0_WP 
        c2=(vm*vm+sum(v_ice(elnodes)*v_ice(elnodes)))/12.0_WP
        c3=(um*vm+sum(v_ice(elnodes)*u_ice(elnodes)))/12.0_WP
        c4=sum(dx*u_ice(elnodes)+dy*v_ice(elnodes))
        !_______________________________________________________________________
        do n=1,3
            row=elnodes(n)
            do q = 1,3 
                entries(q)= vol*ice%ice_dt*((1.0_WP-0.5_WP*ice%ice_dt*c4)*(dx(n)*(um+u_ice(elnodes(q)))+ &
                            dy(n)*(vm+v_ice(elnodes(q))))/12.0_WP - &
                            0.5_WP*ice%ice_dt*(c1*dx(n)*dx(q)+c2*dy(n)*dy(q)+c3*(dx(n)*dy(q)+dx(q)*dy(n))))
                        !um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0)
                entries2(q)=0.5_WP*ice%ice_dt*(dx(n)*(um+u_ice(elnodes(q)))+ &
                            dy(n)*(vm+v_ice(elnodes(q)))-dx(q)*(um+u_ice(row))- &
                            dy(q)*(vm+v_ice(row)))  
            end do
            !___________________________________________________________________
            cx1=vol*ice%ice_dt*c4*(sum(m_ice(elnodes))+m_ice(elnodes(n))+sum(entries2*m_ice(elnodes)))/12.0_WP
            cx2=vol*ice%ice_dt*c4*(sum(a_ice(elnodes))+a_ice(elnodes(n))+sum(entries2*a_ice(elnodes)))/12.0_WP
            cx3=vol*ice%ice_dt*c4*(sum(m_snow(elnodes))+m_snow(elnodes(n))+sum(entries2*m_snow(elnodes)))/12.0_WP
#if defined (__oifs)
            cx4=vol*ice%ice_dt*c4*(sum(ice_temp(elnodes))+ice_temp(elnodes(n))+sum(entries2*ice_temp(elnodes)))/12.0_WP
#endif /* (__oifs) */
            !___________________________________________________________________
            rhs_m(row)=rhs_m(row)+sum(entries*m_ice(elnodes))+cx1
            rhs_a(row)=rhs_a(row)+sum(entries*a_ice(elnodes))+cx2
            rhs_ms(row)=rhs_ms(row)+sum(entries*m_snow(elnodes))+cx3
#if defined (__oifs)
            rhs_temp(row)=rhs_temp(row)+sum(entries*ice_temp(elnodes))+cx4
#endif /* (__oifs) */
            !___________________________________________________________________
            rhs_mdiv(row)=rhs_mdiv(row)-cx1
            rhs_adiv(row)=rhs_adiv(row)-cx2
            rhs_msdiv(row)=rhs_msdiv(row)-cx3
#if defined (__oifs)
            rhs_tempdiv(row)=rhs_tempdiv(row)-cx4
#endif /* (__oifs) */
        end do
    end do
end subroutine ice_TG_rhs_div 
!
!
!_______________________________________________________________________________
subroutine ice_update_for_div(ice, partit, mesh)
    USE MOD_ICE    
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use i_PARAM
    use o_PARAM
    USE g_CONFIG
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target   :: ice
    type(t_partit), intent(inout), target   :: partit
    type(t_mesh)  , intent(in)   , target   :: mesh
    !___________________________________________________________________________
    integer                                 :: n,i,clo,clo2,cn,location(100),row
    real(kind=WP)                           :: rhs_new
    integer                                 :: num_iter_solve=3
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer :: a_icel, m_icel, m_snowl
    real(kind=WP), dimension(:), pointer :: rhs_a, rhs_m, rhs_ms 
    real(kind=WP), dimension(:), pointer :: rhs_adiv, rhs_mdiv, rhs_msdiv
    real(kind=WP), dimension(:), pointer :: da_ice, dm_ice, dm_snow
    real(kind=WP), dimension(:), pointer :: mass_matrix
#if defined (__oifs)
    real(kind=WP), dimension(:), pointer :: ice_temp, m_templ, rhs_temp, rhs_tempdiv, dm_temp
#endif /* (__oifs) */
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice       => ice%data(1)%values(:)
    a_icel      => ice%data(1)%valuesl(:)
    rhs_a       => ice%data(1)%values_rhs(:)
    rhs_adiv    => ice%data(1)%values_div_rhs(:)
    da_ice      => ice%data(1)%dvalues(:)
    
    m_ice       => ice%data(2)%values(:)
    m_icel      => ice%data(2)%valuesl(:)
    rhs_m       => ice%data(2)%values_rhs(:)
    rhs_mdiv    => ice%data(2)%values_div_rhs(:)
    dm_ice      => ice%data(2)%dvalues(:)
    
    m_snow      => ice%data(3)%values(:)
    m_snowl     => ice%data(3)%valuesl(:)
    rhs_ms      => ice%data(3)%values_rhs(:)
    rhs_msdiv   => ice%data(3)%values_div_rhs(:)
    dm_snow     => ice%data(3)%dvalues(:)
#if defined (__oifs)
    ice_temp    => ice%data(4)%values(:)
    m_templ     => ice%data(4)%valuesl(:)
    rhs_temp    => ice%data(4)%values_rhs(:)
    rhs_tempdiv => ice%data(4)%values_div_rhs(:)
    dm_temp     => ice%data(4)%dvalues(:)
#endif /* (__oifs) */    
    mass_matrix => ice%work%fct_massmatrix(:)
    
    !___________________________________________________________________________
    ! Does Taylor-Galerkin solution
    ! --> the first approximation
    do row=1,myDim_nod2D
        ! if cavity node skip it 
        if (ulevels_nod2d(row)>1) cycle
        dm_ice(row) =rhs_mdiv(row) /area(1,row)
        da_ice(row) =rhs_adiv(row) /area(1,row)
        dm_snow(row)=rhs_msdiv(row)/area(1,row)
#if defined (__oifs)
        dm_temp(row)=rhs_tempdiv(row)/area(1,row)
#endif /* (__oifs) */
    end do
    call exchange_nod(dm_ice, partit)
    call exchange_nod(da_ice, partit)
    call exchange_nod(dm_snow, partit)
#if defined (__oifs)
    call exchange_nod(dm_temp, partit)
#endif /* (__oifs) */

    !___________________________________________________________________________
    ! --> iterate 
    do n=1,num_iter_solve-1
        do row=1,myDim_nod2D
            ! if cavity node skip it 
            if (ulevels_nod2d(row)>1) cycle
            
            clo=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)+1
            clo2=ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(1)
            cn=clo2-clo+1
            location(1:cn)=nn_pos(1:cn, row)
            
            rhs_new     = rhs_mdiv(row) - sum(mass_matrix(clo:clo2)*dm_ice(location(1:cn)))
            m_icel(row) = dm_ice(row)+rhs_new/area(1,row)
            rhs_new     = rhs_adiv(row) - sum(mass_matrix(clo:clo2)*da_ice(location(1:cn)))
            a_icel(row) = da_ice(row)+rhs_new/area(1,row)
            rhs_new     = rhs_msdiv(row) - sum(mass_matrix(clo:clo2)*dm_snow(location(1:cn)))
            m_snowl(row)= dm_snow(row)+rhs_new/area(1,row)
#if defined (__oifs)
            rhs_new     = rhs_tempdiv(row) - sum(mass_matrix(clo:clo2)*dm_temp(location(1:cn)))
            m_templ(row)= dm_temp(row)+rhs_new/area(1,row)
#endif /* (__oifs) */
        end do
        do row=1,myDim_nod2D
            ! if cavity node skip it 
            if (ulevels_nod2d(row)>1) cycle
            dm_ice(row) = m_icel(row)
            da_ice(row) = a_icel(row)
            dm_snow(row)= m_snowl(row)
#if defined (__oifs)
            dm_temp(row)= m_templ(row)
#endif /* (__oifs) */
        end do
        call exchange_nod(dm_ice, partit)
        call exchange_nod(da_ice, partit)
        call exchange_nod(dm_snow, partit)
#if defined (__oifs)
        call exchange_nod(dm_temp, partit)
#endif /* (__oifs) */
    end do
    m_ice   = m_ice+dm_ice
    a_ice   = a_ice+da_ice
    m_snow  = m_snow+dm_snow
#if defined (__oifs)
    ice_temp= ice_temp+dm_temp
#endif /* (__oifs) */
end subroutine ice_update_for_div
! =============================================================
