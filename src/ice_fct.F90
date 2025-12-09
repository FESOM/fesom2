module ice_fct_interfaces
  interface
    subroutine ice_mass_matrix_fill(mesh)
      use MOD_MESH
      type(t_mesh), intent(in)              , target :: mesh
    end subroutine

    subroutine ice_solve_high_order(mesh)
      use MOD_MESH
      type(t_mesh), intent(in)              , target :: mesh
    end subroutine

    subroutine ice_solve_low_order(mesh)
      use MOD_MESH
      type(t_mesh), intent(in)              , target :: mesh
    end subroutine
    
    subroutine ice_fem_fct(tr_array_id, mesh)
      use MOD_MESH
      integer   :: tr_array_id
      type(t_mesh), intent(in)              , target :: mesh
    end subroutine
  end interface
end module

! 
! This file collect subroutines implementing FE-FCT
! advection scheme by Loehner et al.
! There is a tunable parameter ice_gamma_fct.
! Increasing it leads to positivity preserving solution.

! Driving routine is fct_ice_solve. It calles other routines 
! that do low-order and figh order solutions and then combine them in a flux
! corrected way. Taylor-Galerkin scheme is used as a high-order one.

! The code is adapted from  FESOM
!
! =====================================================================
subroutine ice_TG_rhs(mesh)
  use MOD_MESH
  use i_Arrays
  use i_PARAM
  use g_PARSUP
  use o_PARAM
  USE g_CONFIG
  implicit none 
  real(kind=WP)   :: diff, entries(3),  um, vm, vol, dx(3), dy(3) 
  integer         :: n, q, row, elem, elnodes(3)
  type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"

    ! Taylor-Galerkin (Lax-Wendroff) rhs
    DO row=1, myDim_nod2D
        rhs_m(row)=0._WP
        rhs_a(row)=0._WP
        rhs_ms(row)=0._WP        
#if defined (__oifs)
        ths_temp(row)=0._WP
#endif /* (__oifs) */
    END DO
    
    ! Velocities at nodes
    do elem=1,myDim_elem2D          !assembling rhs over elements
        elnodes=elem2D_nodes(:,elem)
        !_______________________________________________________________________
        ! if cavity element skip it 
        if (ulevels(elem)>1) cycle
        
        !derivatives
        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)
        vol=elem_area(elem)
        !um=sum(U_ice(elnodes))/3.0_WP
        !vm=sum(V_ice(elnodes))/3.0_WP
        um=sum(U_ice(elnodes))
        vm=sum(V_ice(elnodes))
     
        !diffusivity
        diff=ice_diff*sqrt(elem_area(elem)/scale_area)
        DO n=1,3
            row=elnodes(n)
            DO q = 1,3 
                !entries(q)= vol*dt*((dx(n)*um+dy(n)*vm)/3.0_WP - &
                !            diff*(dx(n)*dx(q)+ dy(n)*dy(q))- &
                !	       0.5*dt*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))) 
                entries(q)= vol*ice_dt*((dx(n)*(um+u_ice(elnodes(q)))+ &
                            dy(n)*(vm+v_ice(elnodes(q))))/12.0_WP - &
                            diff*(dx(n)*dx(q)+ dy(n)*dy(q))- &
                            0.5_WP*ice_dt*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0_WP)    
            END DO
            rhs_m(row)=rhs_m(row)+sum(entries*m_ice(elnodes))
            rhs_a(row)=rhs_a(row)+sum(entries*a_ice(elnodes))
            rhs_ms(row)=rhs_ms(row)+sum(entries*m_snow(elnodes))
#if defined (__oifs)
            rhs_temp(row)=rhs_temp(row)+sum(entries*ice_temp(elnodes))
#endif /* (__oifs) */
        END DO
    end do
end subroutine ice_TG_rhs   
!
!----------------------------------------------------------------------------
!
subroutine ice_fct_init(mesh)
  use o_PARAM
  use MOD_MESH
  use i_ARRAYS
  use g_PARSUP
  use ice_fct_interfaces
  implicit none
  integer   :: n_size
  type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"

  
  n_size=myDim_nod2D+eDim_nod2D
  
  ! Initialization of arrays necessary to implement FCT algorithm
  allocate(m_icel(n_size), a_icel(n_size), m_snowl(n_size))  ! low-order solutions
  m_icel=0.0_WP
  a_icel=0.0_WP 
  m_snowl=0.0_WP
#if defined (__oifs)
  allocate(m_templ(n_size))  
  allocate(dm_temp(n_size))  
#endif /* (__oifs) */
  allocate(icefluxes(myDim_elem2D,3))
  allocate(icepplus(n_size), icepminus(n_size))
  icefluxes = 0.0_WP
  icepplus = 0.0_WP
  icepminus= 0.0_WP
  
#if defined (__oifs)
  m_templ=0.0_WP
  dm_temp=0.0_WP
#endif /* (__oifs) */
  
  allocate(dm_ice(n_size), da_ice(n_size), dm_snow(n_size))  ! increments of high
  dm_ice = 0.0_WP                                            ! order solutions
  da_ice = 0.0_WP
  dm_snow = 0.0_WP
  
  ! Fill in  the mass matrix    
  call ice_mass_matrix_fill(mesh)
  if (mype==0) write(*,*) 'Ice FCT is initialized' 
end subroutine ice_fct_init
!
!----------------------------------------------------------------------------
!
subroutine ice_fct_solve(mesh)
  use MOD_MESH
  use ice_fct_interfaces
  implicit none
  type(t_mesh), intent(in)              , target :: mesh
  ! Driving routine
  call ice_solve_high_order(mesh)   ! uses arrays of low-order solutions as temp
                                    ! storage. It should preceed the call of low
                                    ! order solution.  
  call ice_solve_low_order(mesh)

  call ice_fem_fct(1, mesh)    ! m_ice
  call ice_fem_fct(2, mesh)    ! a_ice
  call ice_fem_fct(3, mesh)    ! m_snow
#if defined (__oifs)
  call ice_fem_fct(4, mesh)    ! ice_temp
#endif /* (__oifs) */

end subroutine ice_fct_solve
!
!
!_______________________________________________________________________________
subroutine ice_solve_low_order(mesh)
 
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
    
    use MOD_MESH
    use o_MESH
    use i_ARRAYS
    use i_PARAM
    use g_PARSUP
    use g_comm_auto
    implicit none
    integer       :: row, clo, clo2, cn, location(100)
    real(kind=WP) :: gamma
    type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"

    gamma=ice_gamma_fct         ! Added diffusivity parameter
                                ! Adjust it to ensure posivity of solution    
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
    call exchange_nod(m_icel,a_icel,m_snowl)

#if defined (__oifs)
    call exchange_nod(m_templ)
#endif /* (__oifs) */


end subroutine ice_solve_low_order     
!
!
!_______________________________________________________________________________
subroutine ice_solve_high_order(mesh)

  use MOD_MESH
  use O_MESH
  use i_ARRAYS
  use g_PARSUP
  use o_PARAM
  use g_comm_auto
  implicit none
  !
  integer                              :: n,i,clo,clo2,cn,location(100),row
  real(kind=WP)                        :: rhs_new
  integer                              :: num_iter_solve=3
  type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"
  ! Does Taylor-Galerkin solution
  !
  !the first approximation
  do row=1,myDim_nod2D
     !___________________________________________________________________________
     ! if cavity node skip it 
     if (ulevels_nod2d(row)>1) cycle
     
     dm_ice(row)=rhs_m(row)/area(1,row)
     da_ice(row)=rhs_a(row)/area(1,row)
     dm_snow(row)=rhs_ms(row)/area(1,row)
#if defined (__oifs)
     dm_temp(row)=rhs_temp(row)/area(1,row)
#endif /* (__oifs) */
  end do

  call exchange_nod(dm_ice, da_ice, dm_snow)

#if defined (__oifs)
     call exchange_nod(dm_temp)
#endif /* (__oifs) */
  !iterate 
  do n=1,num_iter_solve-1
     do row=1,myDim_nod2D
        !___________________________________________________________________________
        ! if cavity node skip it 
        if (ulevels_nod2d(row)>1) cycle
        
        clo=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)+1
        clo2=ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(1)
        cn=clo2-clo+1
        location(1:cn)=nn_pos(1:cn,row)
        rhs_new=rhs_m(row) - sum(mass_matrix(clo:clo2)*dm_ice(location(1:cn)))
        m_icel(row)=dm_ice(row)+rhs_new/area(1,row)
        rhs_new=rhs_a(row) - sum(mass_matrix(clo:clo2)*da_ice(location(1:cn)))
        a_icel(row)=da_ice(row)+rhs_new/area(1,row)
        rhs_new=rhs_ms(row) - sum(mass_matrix(clo:clo2)*dm_snow(location(1:cn)))
        m_snowl(row)=dm_snow(row)+rhs_new/area(1,row) 
#if defined (__oifs)
        rhs_new=rhs_temp(row) - sum(mass_matrix(clo:clo2)*dm_temp(location(1:cn)))
        m_templ(row)=dm_temp(row)+rhs_new/area(1,row)
#endif /* (__oifs) */
     end do
     do row=1,myDim_nod2D
        !_______________________________________________________________________
        ! if cavity node skip it 
        if (ulevels_nod2d(row)>1) cycle
        
        dm_ice(row)=m_icel(row)
        da_ice(row)=a_icel(row)
        dm_snow(row)=m_snowl(row)
#if defined (__oifs)
        dm_temp(row)=m_templ(row)
#endif /* (__oifs) */
     end do
     call exchange_nod(dm_ice, da_ice, dm_snow)

#if defined (__oifs)
     call exchange_nod(dm_temp)
#endif /* (__oifs) */

  end do
end subroutine ice_solve_high_order
!
!
!_______________________________________________________________________________
subroutine ice_fem_fct(tr_array_id, mesh)
    ! Flux corrected transport algorithm for tracer advection
    !
    ! It is based on Loehner et al. (Finite-element flux-corrected 
    ! transport (FEM-FCT) for the Euler and Navier-Stokes equation, 
    ! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
    ! Turek. (kuzmin@math.uni-dortmund.de) 
    !

    use MOD_MESH
    use O_MESH
    use i_arrays
    use i_param
    use o_PARAM
    use g_PARSUP
    use g_comm_auto
    implicit none

    integer   :: tr_array_id
    integer   :: icoef(3,3),n,q, elem,elnodes(3),row
    real(kind=WP), allocatable, dimension(:) :: tmax, tmin 
    real(kind=WP)   :: vol, flux, ae, gamma
    type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"
  
    gamma=ice_gamma_fct        ! It should coinside with gamma in 
                             ! ts_solve_low_order  
  
    !==========================
    ! Compute elemental antidiffusive fluxes to nodes
    !==========================
    ! This is the most unpleasant part --- 
    ! it takes memory and time. For every element 
    ! we need its antidiffusive contribution to 
    ! each of its 3 nodes
    
    allocate(tmax(myDim_nod2D), tmin(myDim_nod2D))
    tmax = 0.0_WP
    tmin = 0.0_WP
    
    ! Auxiliary elemental operator (mass matrix- lumped mass matrix)
    icoef=1
    do n=1,3   ! three upper nodes
        ! Cycle over rows  row=elnodes(n)
        icoef(n,n)=-2
    end do	    
    
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
     
    !==========================   
    ! Screening the low-order solution
    !==========================
    ! TO BE ADDED IF FOUND NECESSARY
    ! Screening means comparing low-order solutions with the
    ! solution on the previous time step and using whichever 
    ! is greater/smaller in computations of max/min below
    
    !==========================
    ! Cluster min/max
    !==========================
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
 
    !=========================
    ! Sums of positive/negative fluxes to node row
    !=========================
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
        
    !========================
    ! The least upper bound for the correction factors
    !========================
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
    call exchange_nod(icepminus, icepplus)
    
    !======================== 
    ! Limiting
    !======================== 
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
  
    !==========================
    ! Update the solution 
    !==========================
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
    
    call exchange_nod(m_ice, a_ice, m_snow)

#if defined (__oifs)
    call exchange_nod(ice_temp)
#endif /* (__oifs) */    

    deallocate(tmin, tmax)
end subroutine ice_fem_fct
!
!
!_______________________________________________________________________________
SUBROUTINE ice_mass_matrix_fill(mesh)
! Used in ice_fct inherited from FESOM
  use MOD_MESH
  use O_MESH
  use i_PARAM
  use i_ARRAYS
  use g_PARSUP
  !
  implicit none
  integer                             :: n, n1, n2, row

  integer                             :: elem, elnodes(3), q, offset, col, ipos 
  integer, allocatable                :: col_pos(:)
  real(kind=WP)                       :: aa
  integer                             :: flag=0,iflag=0
  type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"
    !
    ! a)
    allocate(mass_matrix(sum(nn_num(1:myDim_nod2D))))
    mass_matrix =0.0_WP
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
END SUBROUTINE ice_mass_matrix_fill
!
!=========================================================
!
subroutine ice_TG_rhs_div(mesh)
  use MOD_MESH
  use i_Arrays
  use i_PARAM
  use g_PARSUP
  use o_PARAM
  USE g_CONFIG
  implicit none 
  real(kind=WP)            :: diff, entries(3),  um, vm, vol, dx(3), dy(3) 
  integer                  :: n, q, row, elem, elnodes(3)
  real(kind=WP)            :: c1, c2, c3, c4, cx1, cx2, cx3, cx4, entries2(3) 
  type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

 ! Computes the rhs in a Taylor-Galerkin way (with upwind type of 
 ! correction for the advection operator)
 ! In this version I tr to split divergent term off, so that FCT works without it.

  DO row=1, myDim_nod2D
                  !! row=myList_nod2D(m)
     rhs_m(row)=0.0_WP
     rhs_a(row)=0.0_WP
     rhs_ms(row)=0.0_WP
#if defined (__oifs)
     rhs_temp(row)=0.0_WP
#endif /* (__oifs) */
     rhs_mdiv(row)=0.0_WP
     rhs_adiv(row)=0.0_WP
     rhs_msdiv(row)=0.0_WP
#if defined (__oifs)
     rhs_tempdiv(row)=0.0_WP        
#endif /* (__oifs) */
  END DO
  do elem=1,myDim_elem2D          !assembling rhs over elements
     elnodes=elem2D_nodes(:,elem)
     !___________________________________________________________________________
     ! if cavity element skip it 
     if (ulevels(elem)>1) cycle
     
      !derivatives
     dx=gradient_sca(1:3,elem)
     dy=gradient_sca(4:6,elem)
     vol=elem_area(elem)
     um=sum(u_ice(elnodes))
     vm=sum(v_ice(elnodes))
      ! this is exact computation (no assumption of u=const on elements used 
      ! in the standard version)
     c1=(um*um+sum(u_ice(elnodes)*u_ice(elnodes)))/12.0_WP 
     c2=(vm*vm+sum(v_ice(elnodes)*v_ice(elnodes)))/12.0_WP
     c3=(um*vm+sum(v_ice(elnodes)*u_ice(elnodes)))/12.0_WP
     c4=sum(dx*u_ice(elnodes)+dy*v_ice(elnodes))
     DO n=1,3
        row=elnodes(n)
!!PS         if(ulevels_nod2D(row)>1) cycle !LK89140
        DO q = 1,3 
            entries(q)= vol*ice_dt*((1.0_WP-0.5_WP*ice_dt*c4)*(dx(n)*(um+u_ice(elnodes(q)))+ &
                        dy(n)*(vm+v_ice(elnodes(q))))/12.0_WP - &
                        0.5_WP*ice_dt*(c1*dx(n)*dx(q)+c2*dy(n)*dy(q)+c3*(dx(n)*dy(q)+dx(q)*dy(n))))
                       !um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0)
            entries2(q)=0.5_WP*ice_dt*(dx(n)*(um+u_ice(elnodes(q)))+ &
                        dy(n)*(vm+v_ice(elnodes(q)))-dx(q)*(um+u_ice(row))- &
                        dy(q)*(vm+v_ice(row)))  
        END DO
        cx1=vol*ice_dt*c4*(sum(m_ice(elnodes))+m_ice(elnodes(n))+sum(entries2*m_ice(elnodes)))/12.0_WP
        cx2=vol*ice_dt*c4*(sum(a_ice(elnodes))+a_ice(elnodes(n))+sum(entries2*a_ice(elnodes)))/12.0_WP
        cx3=vol*ice_dt*c4*(sum(m_snow(elnodes))+m_snow(elnodes(n))+sum(entries2*m_snow(elnodes)))/12.0_WP
#if defined (__oifs)
        cx4=vol*ice_dt*c4*(sum(ice_temp(elnodes))+ice_temp(elnodes(n))+sum(entries2*ice_temp(elnodes)))/12.0_WP
#endif /* (__oifs) */

        rhs_m(row)=rhs_m(row)+sum(entries*m_ice(elnodes))+cx1
        rhs_a(row)=rhs_a(row)+sum(entries*a_ice(elnodes))+cx2
        rhs_ms(row)=rhs_ms(row)+sum(entries*m_snow(elnodes))+cx3
#if defined (__oifs)
        rhs_temp(row)=rhs_temp(row)+sum(entries*ice_temp(elnodes))+cx4
#endif /* (__oifs) */
        
        rhs_mdiv(row)=rhs_mdiv(row)-cx1
        rhs_adiv(row)=rhs_adiv(row)-cx2
        rhs_msdiv(row)=rhs_msdiv(row)-cx3
#if defined (__oifs)
        rhs_tempdiv(row)=rhs_tempdiv(row)-cx4
#endif /* (__oifs) */

     END DO
  end do
end subroutine ice_TG_rhs_div 
!
!
!_______________________________________________________________________________
subroutine ice_update_for_div(mesh)
    use MOD_MESH
    use O_MESH
    use i_Arrays
    use i_PARAM
    use g_PARSUP
    use o_PARAM
    USE g_CONFIG
    use g_comm_auto
    implicit none
    !
    integer                                 :: n,i,clo,clo2,cn,location(100),row
    real(kind=WP)                           :: rhs_new
    integer                                 :: num_iter_solve=3
    type(t_mesh), intent(in)                , target :: mesh

#include "associate_mesh.h"

    ! Does Taylor-Galerkin solution
    !
    !the first approximation
    do row=1,myDim_nod2D
                 !! row=myList_nod2D(m)
        !___________________________________________________________________________
        ! if cavity node skip it 
        if (ulevels_nod2d(row)>1) cycle
        dm_ice(row) =rhs_mdiv(row) /area(1,row)
        da_ice(row) =rhs_adiv(row) /area(1,row)
        dm_snow(row)=rhs_msdiv(row)/area(1,row)
#if defined (__oifs)
        dm_temp(row)=rhs_tempdiv(row)/area(1,row)
#endif /* (__oifs) */
    end do
    call exchange_nod(dm_ice)
    call exchange_nod(da_ice)
    call exchange_nod(dm_snow)
#if defined (__oifs)
    call exchange_nod(dm_temp)
#endif /* (__oifs) */

    !iterate 
    do n=1,num_iter_solve-1
        do row=1,myDim_nod2D
            !___________________________________________________________________________
            ! if cavity node skip it 
            if (ulevels_nod2d(row)>1) cycle
                  !! row=myList_nod2D(m)
            clo=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)+1
            clo2=ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(1)
            cn=clo2-clo+1
            location(1:cn)=nn_pos(1:cn, row)
            rhs_new=rhs_mdiv(row) - sum(mass_matrix(clo:clo2)*dm_ice(location(1:cn)))
            m_icel(row)=dm_ice(row)+rhs_new/area(1,row)
            rhs_new=rhs_adiv(row) - sum(mass_matrix(clo:clo2)*da_ice(location(1:cn)))
            a_icel(row)=da_ice(row)+rhs_new/area(1,row)
            rhs_new=rhs_msdiv(row) - sum(mass_matrix(clo:clo2)*dm_snow(location(1:cn)))
            m_snowl(row)=dm_snow(row)+rhs_new/area(1,row)
#if defined (__oifs)
            rhs_new=rhs_tempdiv(row) - sum(mass_matrix(clo:clo2)*dm_temp(location(1:cn)))
            m_templ(row)=dm_temp(row)+rhs_new/area(1,row)
#endif /* (__oifs) */
        end do
        do row=1,myDim_nod2D
            !___________________________________________________________________________
            ! if cavity node skip it 
            if (ulevels_nod2d(row)>1) cycle
                  !! row=myList_nod2D(m)
            dm_ice(row)=m_icel(row)
            da_ice(row)=a_icel(row)
            dm_snow(row)=m_snowl(row)
#if defined (__oifs)
            dm_temp(row)=m_templ(row)
#endif /* (__oifs) */
        end do
        call exchange_nod(dm_ice)
        call exchange_nod(da_ice)
        call exchange_nod(dm_snow)
#if defined (__oifs)
        call exchange_nod(dm_temp)
#endif /* (__oifs) */
    end do
    m_ice=m_ice+dm_ice
    a_ice=a_ice+da_ice
    m_snow=m_snow+dm_snow
#if defined (__oifs)
    ice_temp=ice_temp+dm_temp
#endif /* (__oifs) */

end subroutine ice_update_for_div
! =============================================================
