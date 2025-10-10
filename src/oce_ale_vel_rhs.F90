module oce_ale_vel_rhs_module
    USE MOD_ICE
    USE MOD_DYN
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_ARRAYS
    use o_PARAM
    use g_CONFIG
    use g_forcing_param
    use g_forcing_arrays
    use g_comm_auto
    use g_sbf
    
    implicit none
    
    private
    public :: compute_vel_rhs, momentum_adv_scalar

contains

!
!
!_______________________________________________________________________________
subroutine compute_vel_rhs(ice, dynamics, partit, mesh)
    USE MOD_ICE
    USE MOD_DYN
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_ARRAYS, only: ssh_gp, pgf_x, pgf_y
    use o_PARAM
    use g_CONFIG
    use g_forcing_param, only: use_virt_salt
    use g_forcing_arrays, only: press_air
    use g_comm_auto
    use g_sbf, only: l_mslp
    use momentum_adv_scalar_transpv_interface
    implicit none 
    type(t_ice)   , intent(inout), target :: ice
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: elem, elnodes(3), nz, nzmax, nzmin 
    real(kind=WP)            :: ff, mm 
    real(kind=WP)            :: Fx, Fy, pre(3)
    logical, save            :: lfirst=.true.
    real(kind=WP)            :: p_ice(3), p_air(3), p_eta(3)
    integer                  :: use_pice
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:),   pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:,:,:), pointer :: UV_rhsAB
    real(kind=WP), dimension(:)    ,   pointer :: eta_n
    real(kind=WP), dimension(:)    ,   pointer :: m_ice, m_snow, a_ice
    real(kind=WP)                  ,   pointer :: rhoice, rhosno, inv_rhowat
    real(kind=WP)                              :: ab1, ab2, ab3 !Adams-Bashforth coefficients
    real(kind=WP), dimension(:,:,:),   pointer :: UVh
    real(kind=WP), dimension(:,:)  ,   pointer :: UVBT_4AB
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV        => dynamics%uv(:,:,:)
    UV_rhs    => dynamics%uv_rhs(:,:,:)
    UV_rhsAB  => dynamics%uv_rhsAB(:,:,:,:)
    eta_n     => dynamics%eta_n(:)
    m_ice     => ice%data(2)%values(:)
    m_snow    => ice%data(3)%values(:)
    rhoice    => ice%thermo%rhoice
    rhosno    => ice%thermo%rhosno
    inv_rhowat=> ice%thermo%inv_rhowat
    
    ! if split-explicite ssh subcycling is used
    if (dynamics%use_ssh_se_subcycl) then
        UVh      => dynamics%se_uvh
        UVBT_4AB => dynamics%se_uvBT_4AB
    end if 
    
    !___________________________________________________________________________
    use_pice=0
    if (use_floatice .and.  .not. trim(which_ale)=='linfs') use_pice=1
    if ((toy_ocean)  .and. (trim(which_toy)=="soufflet"))   use_pice=0

    IF     (dynamics%AB_order==2)  THEN
            ab1=-(0.5_WP+epsilon)
            ab2= (1.5_WP+epsilon)
            ab3=  0.0_WP
    ELSEIF (dynamics%AB_order==3) THEN
            ab1=  5.0_WP/12.0_WP
            ab2=-16.0_WP/12.0_WP
            ab3= 23.0_WP/12.0_WP
    ELSE 
       write(*,*) 'unsuppported AB scheme for momentum, use 2 or 3'
       call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
       stop
    END IF 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax, elnodes, ff, mm, Fx, Fy, pre, p_ice, p_air, p_eta)
    do elem=1, myDim_elem2D
        nzmax = nlevels(elem)
        nzmin = ulevels(elem)
        !_______________________________________________________________________
        ! Take care of the AB part
        ! 2nd order Adams-Bashfort
        if (dynamics%AB_order==2) then
            do nz=nzmin,nzmax-1
                UV_rhs(1,nz,elem)=ab1*UV_rhsAB(1,1,nz,elem)
                UV_rhs(2,nz,elem)=ab1*UV_rhsAB(1,2,nz,elem)
                !                        |
                !                        V
                ! Here Adams-Bashfort from previouse time n-0.5
                ! Thats why (0.5_WP+epsilon)*... to achieve second
                ! order AB2 --> fAB2 = 3/2*fab_n - 1/2*fab_n-1   
            end do
            if (dynamics%ldiag_ke) then
                do nz=nzmin,nzmax-1
                    dynamics%ke_adv(:,nz,elem)=ab1*dynamics%ke_adv_AB(1, :,nz,elem)
                    dynamics%ke_cor(:,nz,elem)=ab1*dynamics%ke_cor_AB(1, :,nz,elem)
                end do
            end if
        
        ! 3rd order Adams-Bashfort
        elseif (dynamics%AB_order==3) then
            do nz=nzmin,nzmax-1
                UV_rhs(1,nz,elem)=ab1*UV_rhsAB(2,1,nz,elem)+ab2*UV_rhsAB(1,1,nz,elem)
                UV_rhs(2,nz,elem)=ab1*UV_rhsAB(2,2,nz,elem)+ab2*UV_rhsAB(1,2,nz,elem)
            end do
            if (dynamics%ldiag_ke) then
                do nz=nzmin,nzmax-1
                    dynamics%ke_adv(:,nz,elem)=ab1*dynamics%ke_adv_AB(2,:,nz,elem)+ab2*dynamics%ke_adv_AB(1,:,nz,elem)
                    dynamics%ke_cor(:,nz,elem)=ab1*dynamics%ke_cor_AB(2,:,nz,elem)+ab2*dynamics%ke_cor_AB(1,:,nz,elem)
                end do
            end if      
        end if 
        
        !_______________________________________________________________________
        ! Sea level and pressure contribution   -\nabla(\eta +hpressure/rho_0)
        ! and the Coriolis force + metric terms
        elnodes=elem2D_nodes(:,elem)
        
        !  p_eta=g*eta_n(elnodes)*(1-theta)        !! this place needs update (1-theta)!!!
        p_eta = g*eta_n(elnodes)   
        
        ff    = mesh%coriolis(elem)*elem_area(elem)
        !mm=metric_factor(elem)*elem_area(elem)
        
        !_______________________________________________________________________
        ! contribution from sea level pressure
        if (l_mslp) then
            p_air = press_air(elnodes)/1000
        else                   !|-> convert press_air from: Pa--> bar)
            p_air = 0.0_WP
        end if
        
        !_______________________________________________________________________
        ! in case of ALE zlevel and zstar add pressure from ice to atmospheric pressure
        ! to account for floating ice
        if (use_pice > 0) then
            p_ice = 0.0_WP
            p_ice = (m_ice(elnodes)*rhoice+m_snow(elnodes)*rhosno)*inv_rhowat
            ! limit maximum ice loading like in FESOM1.4
            p_ice = g*min(p_ice,max_ice_loading)
        else
            p_ice = 0.0_WP
        endif
        
        !_______________________________________________________________________
        ! apply pressure gradient force, as well as contributions from gradient of 
        ! the sea surface height as well as ice pressure in case of floating sea ice
        ! to velocity rhs
        pre = -(p_eta+p_ice+p_air)
        if (use_global_tides) then
           pre=pre-ssh_gp(elnodes)+0.1*p_eta
        end if
        Fx  = sum(gradient_sca(1:3, elem)*pre)
        Fy  = sum(gradient_sca(4:6, elem)*pre)
        
        !_______________________________________________________________________
        ! when ssh split-explicite subcycling method is setted use transport velocities
        ! u*h, v*h instead of u,v 
        if (.not. dynamics%use_ssh_se_subcycl) then
            do nz=nzmin, nzmax-1
                ! add pressure gradient terms
                UV_rhs(  1, nz, elem) = UV_rhs(1, nz, elem) + (Fx-pgf_x(nz, elem))*elem_area(elem) 
                UV_rhs(  2, nz, elem) = UV_rhs(2, nz, elem) + (Fy-pgf_y(nz, elem))*elem_area(elem)
                
                ! add coriolis force, initialise AB2 array of actual timestep 
                ! with coriolis term
                if     (dynamics%AB_order==2) then
                    UV_rhsAB(1,1,nz,elem) = UV(2,nz,elem)*ff! + mm*UV(1,nz,elem)*UV(2,nz,elem)
                    UV_rhsAB(1,2,nz,elem) =-UV(1,nz,elem)*ff! - mm*UV(1,nz,elem)*UV(2,nz,elem)
                elseif (dynamics%AB_order==3) then 
                    UV_rhsAB(2,1,nz,elem) = UV_rhsAB(1,1,nz,elem)
                    UV_rhsAB(2,2,nz,elem) = UV_rhsAB(1,2,nz,elem)
                    UV_rhsAB(1,1,nz,elem) = UV(2,nz,elem)*ff
                    UV_rhsAB(1,2,nz,elem) =-UV(1,nz,elem)*ff
                end if 
            end do
        else
            UVBT_4AB(1:2, elem) = 0.0_WP
            do nz=nzmin, nzmax-1
                ! add pressure gradient terms
                UV_rhs(  1, nz, elem) = UV_rhs(1, nz, elem) + (Fx-pgf_x(nz, elem))*elem_area(elem)*helem(nz,elem)
                UV_rhs(  2, nz, elem) = UV_rhs(2, nz, elem) + (Fy-pgf_y(nz, elem))*elem_area(elem)*helem(nz,elem)
                
                ! add coriolis force, initialise AB2 array of actual timestep 
                ! with coriolis term
                if     (dynamics%AB_order==2) then
                    UV_rhsAB(1,1,nz,elem) = UVh(2,nz,elem)*ff! + mm*UV(1,nz,elem)*UV(2,nz,elem)
                    UV_rhsAB(1,2,nz,elem) =-UVh(1,nz,elem)*ff! - mm*UV(1,nz,elem)*UV(2,nz,elem)
                elseif (dynamics%AB_order==3) then 
                    UV_rhsAB(2,1,nz,elem) = UV_rhsAB(1,1,nz,elem)
                    UV_rhsAB(2,2,nz,elem) = UV_rhsAB(1,2,nz,elem)
                    UV_rhsAB(1,1,nz,elem) = UVh(2,nz,elem)*ff
                    UV_rhsAB(1,2,nz,elem) =-UVh(1,nz,elem)*ff
                end if
                
                ! compute barotropic velocity for adams-bashfort time stepping
                ! UVBT_4AB(1:2, elem)--> actual timestep, 
                ! UVBT_4AB(3:4, elem)--> previous timestep (is setted in 
                ! call compute_BC_BT_SE_vtransp) 
                UVBT_4AB(1, elem)     = UVBT_4AB(1, elem) + UVh(1, nz, elem)  ! 
                UVBT_4AB(2, elem)     = UVBT_4AB(2, elem) + UVh(2, nz, elem)  !
            end do    
        end if 
        
        !_______________________________________________________________________
        if (dynamics%ldiag_ke) then
            do nz=nzmin,nzmax-1
                dynamics%ke_pre(1,nz,elem)= (Fx-pgf_x(nz,elem))*dt!*elem_area(elem) !not to divide it aterwards (at the end of this subroutine)
                dynamics%ke_pre(2,nz,elem)= (Fy-pgf_y(nz,elem))*dt!*elem_area(elem) !but account for DT here              
                
                if (dynamics%AB_order==3) then
                    dynamics%ke_cor_AB(2,1,nz,elem) = dynamics%ke_cor_AB(1,1,nz,elem)
                    dynamics%ke_cor_AB(2,2,nz,elem) = dynamics%ke_cor_AB(1,2,nz,elem)
                    
                    dynamics%ke_adv_AB(2,1,nz,elem)= dynamics%ke_adv_AB(1,1,nz,elem)
                    dynamics%ke_adv_AB(2,2,nz,elem)= dynamics%ke_adv_AB(1,2,nz,elem)
                end if
                
                dynamics%ke_cor_AB(1,1,nz,elem)= UV(2,nz,elem)*ff
                dynamics%ke_cor_AB(1,2,nz,elem)=-UV(1,nz,elem)*ff
                
                dynamics%ke_adv_AB(1,1,nz,elem)= 0.0_WP
                dynamics%ke_adv_AB(1,2,nz,elem)= 0.0_WP
            end do
        end if

    end do
!$OMP END PARALLEL DO
    
    !___________________________________________________________________________
    ! advection --> add momentum advection to actual timerstep adams-bashfort 
    ! array UV_rhsAB
    if (dynamics%momadv_opt==1) then
       if (mype==0) write(*,*) 'in moment not adapted mom_adv advection typ for ALE, check your namelist'
       call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
    elseif (dynamics%momadv_opt==2) then
        if (.not. dynamics%use_ssh_se_subcycl) then
            call momentum_adv_scalar(dynamics, partit, mesh)
        else
            call momentum_adv_scalar_transpv(dynamics, partit, mesh)
        end if     
    end if
    
    !___________________________________________________________________________
    ! Update the rhs
    IF (dynamics%AB_order==2) THEN
        ff=ab2
    ELSEIF (dynamics%AB_order==3) THEN
        ff=ab3
    END IF
    
    if (lfirst.and.(.not.r_restart)) then
        ff=1.0_WP
        lfirst=.false.
    end if
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax)
    do elem=1, myDim_elem2D
        nzmin = ulevels(elem)
        nzmax = nlevels(elem)
        do nz=nzmin,nzmax-1
            UV_rhs(1,nz,elem)=dt*(UV_rhs(1,nz,elem)+UV_rhsAB(1,1,nz,elem)*ff)/elem_area(elem)
            UV_rhs(2,nz,elem)=dt*(UV_rhs(2,nz,elem)+UV_rhsAB(1,2,nz,elem)*ff)/elem_area(elem)
            !                       |                   |
            !                       V                   V
            !        fAB = (f_pgf - 1/2*fab_n-1)    +3/2*fab_n
            !        
            ! until here: UV_rhs = dt*[ (R_advec + R_coriolis)^n + R_pressure] 
            ! --> horizontal viscosity contribution still missing is added in 
            !     call viscosity_filter      
            ! --> vertical viscosity contribution still missing is added in 
            !     call impl_vert_visc_ale      
        end do
    end do
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    if (dynamics%ldiag_ke) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax)
        do elem=1, myDim_elem2D
            nzmin = ulevels(elem)
            nzmax = nlevels(elem)
            do nz=nzmin,nzmax-1
                dynamics%ke_adv(:,nz,elem)=dt*(dynamics%ke_adv(:,nz,elem)+dynamics%ke_adv_AB(1,:,nz,elem)*ff)/elem_area(elem)
                dynamics%ke_cor(:,nz,elem)=dt*(dynamics%ke_cor(:,nz,elem)+dynamics%ke_cor_AB(1,:,nz,elem)*ff)/elem_area(elem)
            end do
        end do
!$OMP END PARALLEL DO
    end if
    
    ! =======================  
    ! U_rhs contains all contributions to velocity from old time steps   
    ! =======================
END SUBROUTINE compute_vel_rhs

!
!
!_______________________________________________________________________________
! Momentum advection on scalar control volumes with ALE adaption--> exchange zinv(nz)
! against hnode(nz,node)
subroutine momentum_adv_scalar(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    USE o_PARAM
    use g_comm_auto
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: n, nz, el1, el2
    integer                  :: nl1, nl2, ul1, ul2, nod(2), el, ed, k, nle, ule
    real(kind=WP)            :: un1(1:mesh%nl-1), un2(1:mesh%nl-1)
    real(kind=WP)            :: wu(1:mesh%nl), wv(1:mesh%nl)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:),   pointer :: UV, UVnode_rhs
    real(kind=WP), dimension(:,:,:,:), pointer :: UV_rhsAB
    real(kind=WP), dimension(:,:)    , pointer :: Wvel_e
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV        =>dynamics%uv(:,:,:)
    UV_rhsAB  =>dynamics%uv_rhsAB(:,:,:,:)
    UVnode_rhs=>dynamics%work%uvnode_rhs(:,:,:)
    Wvel_e    =>dynamics%w_e(:,:)

    !___________________________________________________________________________
    ! 1st. compute vertical momentum advection component: w * du/dz, w*dv/dz
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, el1, el2, nl1, nl2, ul1, ul2, nod, el, ed, k, nle, ule, un1, un2, wu, wv)
!$OMP DO
    do n=1,myDim_nod2d
        nl1 = nlevels_nod2D(n)-1
        ul1 = ulevels_nod2D(n)
        wu(1:nl1+1) = 0._WP
        wv(1:nl1+1) = 0._WP
        
        !_______________________________________________________________________
        ! loop over adjacent elements of vertice n 
        do k=1,nod_in_elem2D_num(n)
            el = nod_in_elem2D(k,n)
            !___________________________________________________________________
            nle = nlevels(el)-1
            ule = ulevels(el)
            
            !___________________________________________________________________
            ! accumulate horizontal velocities at full depth levels (top and 
            ! bottom faces of prism) 
            ! account here also for boundary condition below cavity --> 
            ! horizontal velocity at cavity-ocean interce ule (if ule>1) must be  
            ! zero ???
            if (ule==1) then
                wu(ule) = wu(ule) + UV(1,ule,el)*elem_area(el)
                wv(ule) = wv(ule) + UV(2,ule,el)*elem_area(el)
            end if 
            
            ! interpolate horizontal velocity from mid-depth levels to full
            ! depth levels of upper and lower prism faces and average over adjacent
            ! elements of vertice n
            wu(ule+1:nle) = wu(ule+1:nle) + 0.5_WP*(UV(1,ule+1:nle,el)+UV(1,ule:nle-1,el))*elem_area(el)
            wv(ule+1:nle) = wv(ule+1:nle) + 0.5_WP*(UV(2,ule+1:nle,el)+UV(2,ule:nle-1,el))*elem_area(el)
        enddo
        
        !_______________________________________________________________________
        ! multiply w*du and w*dv
        wu(ul1:nl1) = wu(ul1:nl1)*Wvel_e(ul1:nl1,n)
        wv(ul1:nl1) = wv(ul1:nl1)*Wvel_e(ul1:nl1,n)
        
        !_______________________________________________________________________
        ! compute w*du/dz, w*dv/dz
        do nz=ul1,nl1
            ! Here 1/3 because 1/3 of the area is related to the node --> comes from
            ! averaging the elemental velocities
            UVnode_rhs(1,nz,n) = - (wu(nz) - wu(nz+1) ) / (3._WP*hnode(nz,n)) 
            UVnode_rhs(2,nz,n) = - (wv(nz) - wv(nz+1) ) / (3._WP*hnode(nz,n)) 
            
        enddo
        
        !_______________________________________________________________________
        ! To get a clean checksum, set the remaining values to zero
        UVnode_rhs(1:2,nl1+1:nl-1,n) = 0._WP
        UVnode_rhs(1:2,1:ul1-1   ,n) = 0._WP
    end do
!$OMP END DO

    !___________________________________________________________________________
    ! 2nd. compute horizontal advection component: u*du/dx, u*dv/dx & v*du/dy, v*dv/dy
    ! loop over triangle edges
!$OMP DO
    do ed=1, myDim_edge2D
        nod = edges(:,ed)   
        el1 = edge_tri(1,ed)   
        el2 = edge_tri(2,ed)
        nl1 = nlevels(el1)-1
        ul1 = ulevels(el1)
        
        !_______________________________________________________________________
        ! compute horizontal normal velocity with respect to the edge from triangle 
        ! centroid towards triangel edge mid-pointe for element el1
        !                     .o.    
        !                   ./   \.                
        !                 ./  el1  \.   
        !               ./     x     \. 
        !             ./       |-------\.-----------------edge_cross_dxdy(1:2,ed) --> (dx,dy)
        !            /         |->n_vec  \
        !    nod(1) o----------O----------o nod(2)   
        !            \.        |->n_vec ./
        !              \.      |------./------------------edge_cross_dxdy(3:4,ed) --> (dx,dy)
        !                \.    x    ./
        !                  \. el2 ./
        !                    \. ./  
        !                      °
        un1(ul1:nl1) =   UV(2,ul1:nl1,el1)*edge_cross_dxdy(1,ed)   &
                       - UV(1,ul1:nl1,el1)*edge_cross_dxdy(2,ed)  
                       
        !_______________________________________________________________________
        ! compute horizontal normal velocity with respect to the edge from triangle 
        ! centroid towards triangel edge mid-pointe for element el2 when it is valid
        ! --> if its a boundary triangle el2 will be not valid
        if (el2>0) then ! --> el2 is valid element
            nl2 = nlevels(el2)-1
            ul2 = ulevels(el2)
            
            un2(ul2:nl2) = - UV(2,ul2:nl2,el2)*edge_cross_dxdy(3,ed) &
                           + UV(1,ul2:nl2,el2)*edge_cross_dxdy(4,ed)
            
            ! fill with zeros to combine the loops
            ! Usually, no or only a very few levels have to be filled. In this case, 
            ! computing "zeros" is cheaper than the loop overhead.
            un1(nl1+1:max(nl1,nl2)) = 0._WP
            un2(nl2+1:max(nl1,nl2)) = 0._WP
            un1(1:ul1-1)            = 0._WP
            un2(1:ul2-1)            = 0._WP

#if defined(__openmp_reproducible)
!$OMP ORDERED
#endif
            
            ! first edge node
            ! Do not calculate on Halo nodes, as the result will not be used. 
            ! The "if" is cheaper than the avoided computiations.
            if (nod(1) <= myDim_nod2d) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_set_lock(partit%plock(nod(1)))
#endif
                do nz=min(ul1,ul2), max(nl1,nl2)
                    ! add w*du/dz+(u*du/dx+v*du/dy) & w*dv/dz+(u*dv/dx+v*dv/dy)
                    UVnode_rhs(1,nz,nod(1)) = UVnode_rhs(1,nz,nod(1)) + un1(nz)*UV(1,nz,el1) + un2(nz)*UV(1,nz,el2) 
                    UVnode_rhs(2,nz,nod(1)) = UVnode_rhs(2,nz,nod(1)) + un1(nz)*UV(2,nz,el1) + un2(nz)*UV(2,nz,el2)
                end do
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_unset_lock(partit%plock(nod(1)))
#endif
            endif
            
            ! second edge node
            if (nod(2) <= myDim_nod2d) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_set_lock(partit%plock(nod(2)))
#endif
                do nz=min(ul1,ul2), max(nl1,nl2)
                    ! add w*du/dz+(u*du/dx+v*du/dy) & w*dv/dz+(u*dv/dx+v*dv/dy)
                    UVnode_rhs(1,nz,nod(2)) = UVnode_rhs(1,nz,nod(2)) - un1(nz)*UV(1,nz,el1) - un2(nz)*UV(1,nz,el2)
                    UVnode_rhs(2,nz,nod(2)) = UVnode_rhs(2,nz,nod(2)) - un1(nz)*UV(2,nz,el1) - un2(nz)*UV(2,nz,el2)
                end do
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_unset_lock(partit%plock(nod(2)))
#endif
            endif
            
        else  ! el2 is not a valid element --> ed is a boundary edge, there is only the contribution from el1
            ! first edge node
            if (nod(1) <= myDim_nod2d) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_set_lock(partit%plock(nod(1)))
#endif
                do nz=ul1, nl1
                    ! add w*du/dz+(u*du/dx+v*du/dy) & w*dv/dz+(u*dv/dx+v*dv/dy)
                    UVnode_rhs(1,nz,nod(1)) = UVnode_rhs(1,nz,nod(1)) + un1(nz)*UV(1,nz,el1)
                    UVnode_rhs(2,nz,nod(1)) = UVnode_rhs(2,nz,nod(1)) + un1(nz)*UV(2,nz,el1)
                end do ! --> do nz=ul1, nl1
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_unset_lock(partit%plock(nod(1)))
#endif
            endif 
            
            ! second edge node
            if  (nod(2) <= myDim_nod2d) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_set_lock(partit%plock(nod(2)))
#endif
                do nz=ul1, nl1
                    ! add w*du/dz+(u*du/dx+v*du/dy) & w*dv/dz+(u*dv/dx+v*dv/dy)
                    UVnode_rhs(1,nz,nod(2)) = UVnode_rhs(1,nz,nod(2)) - un1(nz)*UV(1,nz,el1)
                    UVnode_rhs(2,nz,nod(2)) = UVnode_rhs(2,nz,nod(2)) - un1(nz)*UV(2,nz,el1)
                end do ! --> do nz=ul1, nl1
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_unset_lock(partit%plock(nod(2)))
#endif
            endif
        endif ! --> if (el2>0) then

#if defined(__openmp_reproducible)
!$OMP END ORDERED
#endif

    end do ! --> do ed=1, myDim_edge2D
!$OMP END DO

    !___________________________________________________________________________
    ! divide total nodal advection by scalar area
!$OMP DO
    do n=1,myDim_nod2d
        nl1 = nlevels_nod2D(n)-1
        ul1 = ulevels_nod2D(n)
        UVnode_rhs(1,ul1:nl1,n) = UVnode_rhs(1,ul1:nl1,n) *areasvol_inv(ul1:nl1,n)
        UVnode_rhs(2,ul1:nl1,n) = UVnode_rhs(2,ul1:nl1,n) *areasvol_inv(ul1:nl1,n)
    end do !-->do n=1,myDim_nod2d
!$OMP END DO
    !___________________________________________________________________________
!$OMP MASTER
    call exchange_nod(UVnode_rhs, partit)
!$OMP END MASTER
!$OMP BARRIER
    !___________________________________________________________________________
    ! convert total nodal advection from vertice --> elements
!$OMP DO
    do el=1, myDim_elem2D
        nl1 = nlevels(el)-1
        ul1 = ulevels(el)
        UV_rhsAB(1,1:2,ul1:nl1,el) = UV_rhsAB(1,1:2,ul1:nl1,el) &
                + elem_area(el)*(UVnode_rhs(1:2,ul1:nl1,elem2D_nodes(1,el)) &
                + UVnode_rhs(1:2,ul1:nl1,elem2D_nodes(2,el)) & 
                + UVnode_rhs(1:2,ul1:nl1,elem2D_nodes(3,el))) / 3.0_WP     
    
    end do ! --> do el=1, myDim_elem2D
!$OMP END DO

    if (dynamics%ldiag_ke) then !we repeat the computation here and there are multiple ways to speed it up
!$OMP DO
       do el=1, myDim_elem2D
          nl1 = nlevels(el)-1
          ul1 = ulevels(el)
          dynamics%ke_adv_AB(1,1:2,ul1:nl1,el) = dynamics%ke_adv_AB(1,1:2,ul1:nl1,el) &
                + elem_area(el)*(UVnode_rhs(1:2,ul1:nl1,elem2D_nodes(1,el)) &
                + UVnode_rhs(1:2,ul1:nl1,elem2D_nodes(2,el)) & 
                + UVnode_rhs(1:2,ul1:nl1,elem2D_nodes(3,el))) / 3.0_WP     
       end do
!$OMP END DO
    end if
!$OMP END PARALLEL
end subroutine momentum_adv_scalar

end module oce_ale_vel_rhs_module

