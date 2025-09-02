!
!
!_______________________________________________________________________________
module momentum_adv_scalar_transpv_interface
    interface
        subroutine momentum_adv_scalar_transpv(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine momentum_adv_scalar_transpv
    end interface
end module momentum_adv_scalar_transpv_interface
!
!
!_______________________________________________________________________________
module impl_vert_visc_ale_vtransp_interface
    interface
        subroutine impl_vert_visc_ale_vtransp(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine impl_vert_visc_ale_vtransp
    end interface
end module impl_vert_visc_ale_vtransp_interface
!
!
!_______________________________________________________________________________
module compute_ssh_split_explicit_interface
    interface
        subroutine compute_BT_rhs_SE_vtransp(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine compute_BT_rhs_SE_vtransp

        subroutine compute_BT_step_SE_ale(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine compute_BT_step_SE_ale
        
        subroutine update_trim_vel_ale_vtransp(mode, dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        integer       , intent(in)            :: mode
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine update_trim_vel_ale_vtransp
        
    end interface
end module compute_ssh_split_explicit_interface
!
!
!_______________________________________________________________________________
! Transports are used instead of velocities, Urhs, Vrhs are also for transports. 
subroutine momentum_adv_scalar_transpv(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE o_PARAM
    USE g_comm_auto
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: node, elem, ed, nz
    integer                  :: nl1, ul1, nl2, ul2, nl12, ul12
    real(kind=WP)            :: uv12, uv1, uv2, qc, qu, qd, num_ord=0.95_WP
    integer                  :: ednodes(2), edelem(2)
    real(kind=WP)            :: wu(mesh%nl), wv(mesh%nl), un1(mesh%nl), un2(mesh%nl)

    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:,:), pointer :: UV_rhsAB
    real(kind=WP), dimension(:,:,:)  , pointer :: UV, UVnode_rhs, UVnode, UVh
    real(kind=WP), dimension(:,:)    , pointer :: Wvel_e
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV        =>dynamics%uv(:,:,:)
    UV_rhsAB  =>dynamics%uv_rhsAB(:,:,:,:)
    UVnode_rhs=>dynamics%work%uvnode_rhs(:,:,:)
    UVnode    =>dynamics%uvnode(:,:,:)
    Wvel_e    =>dynamics%w_e(:,:)
    UVh       =>dynamics%se_uvh(:,:,:)
    !___________________________________________________________________________
    ! 1st. compute vertical momentum advection component: w * du/dz, w*dv/dz
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, elem, ed, nz, nl1, ul1, nl2, ul2, nl12, ul12, &
!$OMP                                  uv1, uv2, uv12, qc, qu, qd, wu, wv, &
!$OMP                                  ednodes, edelem, un1, un2)
!$OMP DO    
    do node=1, myDim_nod2D
        ul1 = ulevels_nod2D(node)
        nl1 = nlevels_nod2D(node)
        wu(1:nl1) = 0.0_WP
        wv(1:nl1) = 0.0_WP
        UVnode_rhs(:,:,node) = 0.0_WP
        !_______________________________________________________________________
        ! surface
        nz=ul1
        wu(nz)=UVnode(1, nz, node)*Wvel_e(nz, node)*area(nz, node)
        wv(nz)=UVnode(2, nz, node)*Wvel_e(nz, node)*area(nz, node)
        
        !_______________________________________________________________________
        ! subsurface --> centered 2nd order
        nz=ul1+1          ! Central differences at the second and last but one levels
        wu(nz)=0.5_WP*(UVnode(1, nz, node)+UVnode(1, nz-1, node))*Wvel_e(nz, node)*area(nz, node)
        wv(nz)=0.5_WP*(UVnode(2, nz, node)+UVnode(2, nz-1, node))*Wvel_e(nz, node)*area(nz, node)
        
        !_______________________________________________________________________
        ! bulk --> centered 4th order (num_ord=1), upwind 3rd order (num_ord=0)
        do nz=ul1+2, nl1-2
            qc   = (UVnode(1, nz-1, node)-UVnode(1, nz  , node))/(hnode(nz-1, node)+hnode(nz  , node))      ! factor 2 is accounted for later   
            qu   = (UVnode(1, nz  , node)-UVnode(1, nz+1, node))/(hnode(nz  , node)+hnode(nz+1, node))
            qd   = (UVnode(1, nz-2, node)-UVnode(1, nz-1, node))/(hnode(nz-2, node)+hnode(nz-1, node))
            
            uv1  =  UVnode(1, nz  , node)+(2*qc+qu)*hnode(nz  , node)/6.0_WP                    ! Gradient reconstruction 2(2qc+qu)(h/2)(1/6)            
            uv2  =  UVnode(1, nz-1, node)-(2*qc+qd)*hnode(nz-1, node)/6.0_WP
            uv12 = (Wvel_e(nz, node)+abs(Wvel_e(nz, node)))*uv1+ &
                   (Wvel_e(nz, node)-abs(Wvel_e(nz, node)))*uv2
            wu(nz)=0.5_WP*(num_ord*(uv1+uv2)*Wvel_e(nz, node)+(1.0_WP-num_ord)*uv12)*area(nz, node)
            
            qc   = (UVnode(2, nz-1, node)-UVnode(2, nz  , node))/(hnode(nz-1, node)+hnode(nz  , node))
            qu   = (UVnode(2, nz  , node)-UVnode(2, nz+1, node))/(hnode(nz  , node)+hnode(nz+1, node))
            qd   = (UVnode(2, nz-2, node)-UVnode(2, nz-1, node))/(hnode(nz-2, node)+hnode(nz-1, node))
            
            uv1  =  UVnode(2, nz  , node)+(2*qc+qu)*hnode(nz  , node)/6.0_WP
            uv2  =  UVnode(2, nz-1, node)-(2*qc+qd)*hnode(nz-1, node)/6.0_WP
            uv12 = (Wvel_e(nz, node)+abs(Wvel_e(nz, node)))*uv1+ &
                   (Wvel_e(nz, node)-abs(Wvel_e(nz, node)))*uv2
            wv(nz)=0.5_WP*(num_ord*(uv1+uv2)*Wvel_e(nz, node)+(1.0_WP-num_ord)*uv12)*area(nz, node)
        end do ! --> do nz=ul1+2, nl1-2
     
        !_______________________________________________________________________
        ! one layer above bottom --> centered 2nd order
        nz=nl1-1
        wu(nz)=0.5_WP*(UVnode(1, nz, node)+UVnode(1, nz-1, node))*Wvel_e(nz, node)*area(nz, node)
        wv(nz)=0.5_WP*(UVnode(2, nz, node)+UVnode(2, nz-1, node))*Wvel_e(nz, node)*area(nz, node)
        
        !_______________________________________________________________________
        ! bottom layer --> boundary condition
        nz=nl1
        wu(nz)=0.0_WP
        wv(nz)=0.0_WP
        
        !_______________________________________________________________________
        ! set to the rhs for transports, not velocities!!! --> No division by h  
        do nz=ul1, nl1-1
            UVnode_rhs(1, nz, node)= UVnode_rhs(1, nz, node) - (wu(nz)-wu(nz+1))
            UVnode_rhs(2, nz, node)= UVnode_rhs(2, nz, node) - (wv(nz)-wv(nz+1))
        end do
        
    end do ! --> do node=1, myDim_nod2D
!$OMP END DO


    !___________________________________________________________________________
    ! 2nd. compute horizontal advection component: u*du/dx, u*dv/dx & v*du/dy, v*dv/dy
    ! loop over triangle edges
!$OMP DO    
    do ed=1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        ednodes = edges(:,ed)   
        
        ! local index of element that contribute to edge
        edelem  = edge_tri(1:2,ed)   
        
        !_______________________________________________________________________
        ! index off surface layer in case of cavity !=1 and index of mid depth 
        ! bottom layer
        ul1     = ulevels(edelem(1))
        nl1     = nlevels(edelem(1))-1
        
        !_______________________________________________________________________
        !NR --> Natalja Style
        un1          = 0.0_WP
        un1(ul1:nl1) = (  UVh(2, ul1:nl1, edelem(1))*edge_cross_dxdy(1,ed) & 
                        - UVh(1, ul1:nl1, edelem(1))*edge_cross_dxdy(2,ed)) 
        
        !_______________________________________________________________________
        ! if edelem(2)==0 than edge is boundary edge
        if(edelem(2)>0) then
            ul2 = ulevels(edelem(2))
            nl2 = nlevels(edelem(2))-1
            
            !___________________________________________________________________
            !NR --> Natalja Style
            un2          = 0.0_WP
            un2(ul2:nl2) = -(  UVh(2, ul2:nl2, edelem(2))*edge_cross_dxdy(3,ed) & 
                             - UVh(1, ul2:nl2, edelem(2))*edge_cross_dxdy(4,ed)) 
            
            !___________________________________________________________________
            ! nl12 ... minimum number of layers -1 between element edelem(1) & edelem(2) that 
            ! contribute to edge ed
            ! nu12 ... upper index of layers between element edelem(1) & edelem(2) that 
            ! contribute to edge ed
            ! be carefull !!! --> if ed is a boundary edge than edelem(1)~=0 and edelem(2)==0
            !                     that means nl1>0, nl2==0, nl12=min(nl1,nl2)=0 !!!
            ul12 = max(ul1, ul2)
            nl12 = min(nl1, nl2)
            
            !___________________________________________________________________
            ! ensure openmp numerical reproducability
#if defined(__openmp_reproducible)
!$OMP ORDERED
#endif
            !___________________________________________________________________
            !NR add contribution to first edge node --> ednodes(1)
            !NR Do not calculate on Halo nodes, as the result will not be used. 
            !NR The "if" is cheaper than the avoided computiations.
            if (ednodes(1) <= myDim_nod2d) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_set_lock(partit%plock(ednodes(1)))
#endif          
                ! cavity domain where only edelem(1) exist
                do nz=ul1 , ul12-1
                    UVnode_rhs(1, nz, ednodes(1)) = UVnode_rhs(1, nz, ednodes(1)) + un1(nz)*UV(1, nz, edelem(1))
                    UVnode_rhs(2, nz, ednodes(1)) = UVnode_rhs(2, nz, ednodes(1)) + un1(nz)*UV(2, nz, edelem(1))
                end do
                ! cavity domain where only edelem(2) exist
                do nz=ul2 , ul12-1
                    UVnode_rhs(1, nz, ednodes(1)) = UVnode_rhs(1, nz, ednodes(1)) + un2(nz)*UV(1, nz, edelem(2))
                    UVnode_rhs(2, nz, ednodes(1)) = UVnode_rhs(2, nz, ednodes(1)) + un2(nz)*UV(2, nz, edelem(2))
                end do
                ! bulk domain where edelem(1) and edelem(2) exist
                do nz=ul12, nl12
                    UVnode_rhs(1, nz, ednodes(1)) = UVnode_rhs(1, nz, ednodes(1)) + (un1(nz) + un2(nz))*(UV(1, nz, edelem(1))+UV(1, nz, edelem(2)))*0.5_WP
                    UVnode_rhs(2, nz, ednodes(1)) = UVnode_rhs(2, nz, ednodes(1)) + (un1(nz) + un2(nz))*(UV(2, nz, edelem(1))+UV(2, nz, edelem(2)))*0.5_WP
                end do
                ! bottom domain where only edelem(1) exist
                do nz=nl12+1, nl1
                    UVnode_rhs(1, nz, ednodes(1)) = UVnode_rhs(1, nz, ednodes(1)) + un1(nz)*UV(1, nz, edelem(1))
                    UVnode_rhs(2, nz, ednodes(1)) = UVnode_rhs(2, nz, ednodes(1)) + un1(nz)*UV(2, nz, edelem(1))
                end do
                ! bottom domain where only edelem(2) exist
                do nz=nl12+1, nl2
                    UVnode_rhs(1, nz, ednodes(1)) = UVnode_rhs(1, nz, ednodes(1)) + un2(nz)*UV(1, nz, edelem(2))
                    UVnode_rhs(2, nz, ednodes(1)) = UVnode_rhs(2, nz, ednodes(1)) + un2(nz)*UV(2, nz, edelem(2))
                end do
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(ednodes(1)))
#endif
            end if
            !___________________________________________________________________
            !NR add contribution to second edge node --> ednodes(2)
            !NR Do not calculate on Halo nodes, as the result will not be used. 
            !NR The "if" is cheaper than the avoided computiations.
            if (ednodes(2) <= myDim_nod2d) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_set_lock(partit%plock(ednodes(2)))
#endif 
                ! cavity domain where only edelem(1) exist
                do nz=ul1 , ul12-1
                    UVnode_rhs(1, nz, ednodes(2)) = UVnode_rhs(1, nz, ednodes(2)) - un1(nz)*UV(1, nz, edelem(1))
                    UVnode_rhs(2, nz, ednodes(2)) = UVnode_rhs(2, nz, ednodes(2)) - un1(nz)*UV(2, nz, edelem(1))
                end do
                ! cavity domain where only edelem(2) exist
                do nz=ul2 , ul12-1
                    UVnode_rhs(1, nz, ednodes(2)) = UVnode_rhs(1, nz, ednodes(2)) - un2(nz)*UV(1, nz, edelem(2))
                    UVnode_rhs(2, nz, ednodes(2)) = UVnode_rhs(2, nz, ednodes(2)) - un2(nz)*UV(2, nz, edelem(2))
                end do
                ! bulk domain where edelem(1) and edelem(2) exist
                do nz=ul12, nl12
                    UVnode_rhs(1, nz, ednodes(2)) = UVnode_rhs(1, nz, ednodes(2)) - (un1(nz) + un2(nz))*(UV(1, nz, edelem(1))+UV(1, nz, edelem(2)))*0.5_WP
                    UVnode_rhs(2, nz, ednodes(2)) = UVnode_rhs(2, nz, ednodes(2)) - (un1(nz) + un2(nz))*(UV(2, nz, edelem(1))+UV(2, nz, edelem(2)))*0.5_WP
                end do
                ! bottom domain where only edelem(1) exist
                do nz=nl12+1, nl1
                    UVnode_rhs(1, nz, ednodes(2)) = UVnode_rhs(1, nz, ednodes(2)) - un1(nz)*UV(1, nz, edelem(1))
                    UVnode_rhs(2, nz, ednodes(2)) = UVnode_rhs(2, nz, ednodes(2)) - un1(nz)*UV(2, nz, edelem(1))
                end do
                ! bottom domain where only edelem(2) exist
                do nz=nl12+1, nl2
                    UVnode_rhs(1, nz, ednodes(2)) = UVnode_rhs(1, nz, ednodes(2)) - un2(nz)*UV(1, nz, edelem(2))
                    UVnode_rhs(2, nz, ednodes(2)) = UVnode_rhs(2, nz, ednodes(2)) - un2(nz)*UV(2, nz, edelem(2))
                end do
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(ednodes(2)))
#endif
            end if 
        
        !_______________________________________________________________________
        ! if edelem(2)==0 than edge is boundary edge
        else ! --> if(edelem(2)>0) then
            !___________________________________________________________________
            !NR add contribution to first edge node --> ednodes(1)
            !NR Do not calculate on Halo nodes, as the result will not be used. 
            !NR The "if" is cheaper than the avoided computiations.
            if (ednodes(1) <= myDim_nod2d) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_set_lock(partit%plock(ednodes(1)))
#endif          
                ! bulk domain where only edelem(1) exist
                do nz=ul1 , nl1
                    UVnode_rhs(1, nz, ednodes(1)) = UVnode_rhs(1, nz, ednodes(1)) + un1(nz)*UV(1, nz, edelem(1))
                    UVnode_rhs(2, nz, ednodes(1)) = UVnode_rhs(2, nz, ednodes(1)) + un1(nz)*UV(2, nz, edelem(1))
                end do
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(ednodes(1)))
#endif
            end if
            !___________________________________________________________________
            !NR add contribution to second edge node --> ednodes(2)
            !NR Do not calculate on Halo nodes, as the result will not be used. 
            !NR The "if" is cheaper than the avoided computiations.
            if (ednodes(2) <= myDim_nod2d) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_set_lock(partit%plock(ednodes(2)))
#endif 
                ! bulk domain where only edelem(1) exist
                do nz=ul1 , nl1
                    UVnode_rhs(1, nz, ednodes(2)) = UVnode_rhs(1, nz, ednodes(2)) - un1(nz)*UV(1, nz, edelem(1))
                    UVnode_rhs(2, nz, ednodes(2)) = UVnode_rhs(2, nz, ednodes(2)) - un1(nz)*UV(2, nz, edelem(1))
                end do
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(ednodes(2)))
#endif
            end if 
        end if ! --> if(edelem(2)>0) then
        
#if defined(__openmp_reproducible)
!$OMP END ORDERED
#endif

    end do ! --> do ed=1, myDim_edge2D  
!$OMP END DO 

    !___________________________________________________________________________
    ! divide total nodal momentum advection by scalar area
!$OMP DO    
    do node=1,myDim_nod2d
        ul1 = ulevels_nod2D(node)
        nl1 = nlevels_nod2D(node)-1
        UVnode_rhs(1, ul1:nl1, node) = UVnode_rhs(1, ul1:nl1, node) * areasvol_inv(ul1:nl1, node)
        UVnode_rhs(2, ul1:nl1, node) = UVnode_rhs(2, ul1:nl1, node) * areasvol_inv(ul1:nl1, node)
    end do ! --> do node=1,myDim_nod2d
!$OMP END DO

    !___________________________________________________________________________
!$OMP MASTER    
    call exchange_nod(UVnode_rhs, partit)
!$OMP END MASTER
!$OMP BARRIER

    !___________________________________________________________________________
    ! convert total nodal advection from vertice --> elements
!$OMP DO    
    do elem=1, myDim_elem2D
        ul1 = ulevels(elem)
        nl1 = nlevels(elem)-1
        UV_rhsAB(1, 1:2, ul1:nl1, elem) = UV_rhsAB(1, 1:2, ul1:nl1, elem) + elem_area(elem)* &
                ( UVnode_rhs(1:2, ul1:nl1, elem2D_nodes(1, elem)) &
                + UVnode_rhs(1:2, ul1:nl1, elem2D_nodes(2, elem)) & 
                + UVnode_rhs(1:2, ul1:nl1, elem2D_nodes(3, elem))) / 3.0_WP     
    end do ! --> do edelem=1, myDim_elem2D
!$OMP END DO    
    
    !___________________________________________________________________________
    ! for energz diagnostic 
    if (dynamics%ldiag_ke) then !we repeat the computation here and there are multiple ways to speed it up
!$OMP DO
        do elem=1, myDim_elem2D
            ul1 = ulevels(elem)
            nl1 = nlevels(elem)-1
            dynamics%ke_adv_AB(1, 1, ul1:nl1, elem) = dynamics%ke_adv_AB(1, 1, ul1:nl1, elem) + elem_area(elem)* &
                    ( UVnode_rhs(1, ul1:nl1, elem2D_nodes(1, elem)) &
                    + UVnode_rhs(1, ul1:nl1, elem2D_nodes(2, elem)) & 
                    + UVnode_rhs(1, ul1:nl1, elem2D_nodes(3, elem))) / 3.0_WP / helem(ul1:nl1, elem) 
            dynamics%ke_adv_AB(1, 2, ul1:nl1, elem) = dynamics%ke_adv_AB(1, 2, ul1:nl1, elem) + elem_area(elem)* &
                    ( UVnode_rhs(2, ul1:nl1, elem2D_nodes(1, elem)) &
                    + UVnode_rhs(2, ul1:nl1, elem2D_nodes(2, elem)) & 
                    + UVnode_rhs(2, ul1:nl1, elem2D_nodes(3, elem))) / 3.0_WP / helem(ul1:nl1, elem)            
        end do
!$OMP END DO
    end if
!$OMP END PARALLEL
end subroutine momentum_adv_scalar_transpv

!
!
!_______________________________________________________________________________
!SD Transport velocity version (T in the name)
!SD Solve  U**-U*=dt*D_vert u**, where U*=Uh+rhs (rhs contains dt/elem_area)
!SD Express through Delta U=U**-U* and then introduce
!SD Delta u= Delta U/helem; after Delta u is found, rhs=rhs +h*Delta u.
!SD The rhs accumulated previously is not affected in this step, this is different
!SD from FESOM, and I would recomment this variant.
!
!
! solve equation vertical viscosity implicite part:
! --> U^(n+0.5,**) - U^(n+0.5,*) = dt*( Av * d/dz * u^(n+0.5,**)  )|^t_b
!        |
!        +-> du = u^(n+0.5,**) - u^(n+0.5,*) = dU/h^(*) 
!        |
!        +-> du * h^(*)   = dU
!        |   u^(n+0.5,**) = du + u^(n+0.5,*)   
!        V
!     du * h^(*) = dt*( Av*d/dz*du )|^t_b + dt*( Av*d/dz*u^(n+0.5,*) )|^t_b
!     du - dt/h^(*)*( Av*d/dz*du ) = dt/h^(*)*( Av*d/dz*u^(n+0.5,*) )   
! --> solve for du
!                                                 ^
!                                                /|\ nvec_up (+1)
!                                                 |
!       ----------- zbar_1, A_1              *----|----*
!   Z_1 o T_1, Ac_1                          |\   |  ./|
!       ----------- zbar_2, A_2              | \   ./  |   Gaus Theorem:
!   Z_2 o T_2, Ac_2                          |  \ /    |    --> Flux form
!       ----------- zbar_3, A_3              |   |     |    --> normal vec outwards facing
!   Z_3 o T_3, Ac_3                          *---|-----*
!       ----------- zbar_4                    \  | | ./
!           :                                  \ | |/
!                                               \|/|
!                                                * |
!                                                  V nvec_dwn (-1)
!
! --> 1st. solve homogenouse part:
!     f(du) = du - dt/h^(*)*( Av*d/dz*du ) = 0
!
! --> 2nd. compute difference quotient at du_i using Gauss-Theorem --> flux form
!        |
!        +-> Gauss THeorem: int(V', div(F_vec))dV = intcircle(A', F_vec*n_vec)dA
!        |
!        +-> du                     = dt/h^(*)*( Av*d/dz*du )  | *div()_z=d/dz, *int(V',)dV
!        |   int(V', d/dz *du)dV    = int(V', d/dz *dt/h^(*)*( Av*d/dz*du ) )dV
!        |   ...                    = intcircle(A', dt/h^(*)*( Av*d/dz*du )*n_vec)dA   
!        |   int(V', d/dz *du)dz*dA = ...
!        |   int(A', du)dA          = intcircle(A', dt/h^(*)*( Av*d/dz*du )*n_vec)dA
!        |   
!        +-> Ac_i area of vector cell = are of triangle, A_i, A_i+1 area of top 
!        |   and bottom face
!        V
!         
!     du_i * Ac_i = dt * [ Av_i  /h_i * (du_i-1 - du_i  )/(Z_i-1 - Z_i  ) * A_i   * nvec_up(+1)  
!                         +Av_i+1/h_i * (du_i   - du_i+1)/(Z_i   - Z_i+1) * A_i+1 * nvec_dwn(-1) ] 
!        |
!        +-> since we are on triangles Ac_i = A_i = A_i+1 --> can kick out A_i/Ac_i
!        |   and A_i+1/Ac_i
!        +-> take into account normal vector direction
!        V
!     f(du_i) = du_i - dt*Av_i  /h_i * (du_i-1 - du_i  )/(Z_i-1 - Z_i  )
!                    + dt*Av_i+1/h_i * (du_i   - du_i+1)/(Z_i   - Z_i+1)
!             = 0
!
! --> 3rd. solve for coefficents a, b, c (homogenous part):
!     f(du_i) = [ a*du_i-1 + b*du_i + c*du_i+1 ]
!        |
!        +-> estimate a, b, c by derivation of f(du_i)
!        |
!        +-> a = d[f(du_i)]/d[du_i-1] = - dt*Av_i/h_i / (Z_i-1 - Z_i) 
!        |
!        +-> c = d[f(du_i)]/d[du_i+1] = - dt*Av_i+1/h_i / (Z_i - Z_i+1)
!        |
!        +-> b = d[f(du_i)]/d[du_i]   = 1 + dt*Av_i/h_i / (Z_i-1 - Z_i) + dt*Av_i+1/h_i / (Z_i - Z_i+1)
!                                       1 - a - c
!
! --> 4th. solve inhomogenous part:
!     [ a*du_i-1 + b*du_i + c*du_i+1 ] = RHS/A_i
!
!     RHS/A_i = dt* [ Av_i  /h_i * (u^(n+0.5,*)_i-1 - u^(n+0.5,*)_i  )/(Z_i-1 - Z_i  ) * A_i  /A_i * nvec_up(+1) 
!             |      +Av_i+1/h_i * (u^(n+0.5,*)_i   - u^(n+0.5,*)_i+1)/(Z_i   - Z_i+1) * A_i+1/A_i * nvec_dwn(-1) ]
!             |
!             +-> since we are on triangles A_i = A_i+1 --> can kick out A_i
!             +-> take into account normal vector direction
!             V
!             = -a*u^(n+0.5,*)_i-1 + (a+c)*u^(n+0.5,*)_i - c*u^(n+0.5,*)_i+1
!
! --> 5th. solve for du_i --> forward sweep algorithm --> see lower
!     | b_1 c_1 ...            |   |du_1|
!     | a_2 b_2 c_2 ...        |   |du_2|
!     |     a_3 b_3 c_3 ...    | * |du_3| = RHS/A_i
!     |         a_4 b_4 c_4 ...|   |du_4|
!     |              :         |   | :  |
!
subroutine impl_vert_visc_ale_vtransp(dynamics, partit, mesh)
    USE MOD_MESH
    USE o_PARAM
    USE o_ARRAYS, only: Av, stress_surf
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE g_CONFIG, only: dt
    IMPLICIT NONE
    !___________________________________________________________________________
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    real(kind=WP)              ::  a(mesh%nl-1), b(mesh%nl-1), c(mesh%nl-1), ur(mesh%nl-1), vr(mesh%nl-1)
    real(kind=WP)              ::  cp(mesh%nl-1), up(mesh%nl-1), vp(mesh%nl-1), uu(mesh%nl-1), vv(mesh%nl-1)  
    integer                    ::  elem, nz, nzmax, nzmin, elnodes(3)
    real(kind=WP)              ::  zinv, m, friction, wu, wd
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs, UVh
    real(kind=WP), dimension(:,:)  , pointer :: Wvel_i
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV     =>dynamics%uv(:,:,:)
    UV_rhs =>dynamics%uv_rhs(:,:,:)
    Wvel_i =>dynamics%w_i(:,:)
    UVh    =>dynamics%se_uvh(:,:,:)
    
    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nz, nzmin, nzmax, elem, elnodes, &
!$OMP                                  a, b, c, m, ur, vr, cp, up, vp,  & 
!$OMP                                  zinv, friction, wu, wd, uu, vv)
    
    !___________________________________________________________________________
!$OMP DO
    do elem=1,myDim_elem2D
        elnodes= elem2D_nodes(:,elem)
        nzmin  = ulevels(elem)
        nzmax  = nlevels(elem)
        uu     = 0.0_WP
        vv     = 0.0_WP
        
        !_______________________________________________________________________
        ! New velocities: compute u_k^(n+1/2, *) = u_k^(n-1/2) + UV_rhs/helem 
        ! UV_rhs is here in terms of transport velocity, we need real velocity 
        ! therefor divide with helem
        uu(nzmin:nzmax-1)=(UVh(1, nzmin:nzmax-1, elem)+UV_rhs(1, nzmin:nzmax-1, elem))/helem(nzmin:nzmax-1, elem)   !! u*=U*/h
        vv(nzmin:nzmax-1)=(UVh(2, nzmin:nzmax-1, elem)+UV_rhs(2, nzmin:nzmax-1, elem))/helem(nzmin:nzmax-1, elem) 
                
        !_______________________________________________________________________
        ! Operator + rhs
        ! Regular part of coefficients:
        
        !_____1st layer (surface)_____
        nz    = nzmin
        zinv  = 2.0_WP*dt/helem(nz, elem)
        c( nz)= -Av(nz+1,elem)/(helem(nz,elem)+helem(nz+1,elem))*zinv     
        a( nz)= 0.0_WP
        b( nz)= -c(nz)+1.0_WP
        
        ! Update from the vertical advection
        wu    = sum(Wvel_i(nz  , elnodes))/3._WP
        wd    = sum(Wvel_i(nz+1, elnodes))/3._WP
        b( nz)= b(nz)+wu*zinv
        b( nz)= b(nz)-min(0._WP, wd)*zinv
        c( nz)= c(nz)-max(0._WP, wd)*zinv
        
        ur(nz)= +c(nz)*uu(nz)-c(nz)*uu(nz+1)
        vr(nz)= +c(nz)*vv(nz)-c(nz)*vv(nz+1)
        
        !_____bulk layers_____
        do nz=nzmin+1, nzmax-2
            zinv  = 2.0_WP*dt/helem(nz, elem)     
            a( nz)= -Av(nz  , elem)/(helem(nz-1, elem)+helem(nz  , elem))*zinv
            c( nz)= -Av(nz+1, elem)/(helem(nz  , elem)+helem(nz+1, elem))*zinv
            b( nz)= -a(nz)-c(nz)+1.0_WP
            
            ! Update from the vertical advection
            wu=sum(Wvel_i(nz  , elnodes))/3._WP
            wd=sum(Wvel_i(nz+1, elnodes))/3._WP
            a( nz)= a(nz)+min(0._WP, wu)*zinv
            b( nz)= b(nz)+max(0._WP, wu)*zinv
            
            ur(nz)= -a(nz)*uu(nz-1)+(a(nz)+c(nz))*uu(nz)-c(nz)*uu(nz+1)               !! add dt*D_vert u*
            vr(nz)= -a(nz)*vv(nz-1)+(a(nz)+c(nz))*vv(nz)-c(nz)*vv(nz+1)               !! to the rhs (because of Delta u)
        end do
        
        !_____bottom layer_____
        nz    = nzmax-1
        zinv  = 2.0_WP*dt/helem(nz,elem)
        a( nz)= -Av(nz, elem)/(helem(nz-1, elem)+helem(nz,elem))*zinv
        b( nz)= -a(nz)+1.0_WP
        c( nz)= 0.0_WP
        
        ! Update from the vertical advection
        wu    = sum(Wvel_i(nz-1, elnodes))/3._WP
        a( nz)= a(nz)+min(0._WP, wu)*zinv
        b( nz)= b(nz)+max(0._WP, wu)*zinv
        
        ur(nz)= -a(nz)*uu(nz-1)+a(nz)*uu(nz)
        vr(nz)= -a(nz)*vv(nz-1)+a(nz)*vv(nz)

        !_______________________________________________________________________
        ! Forcing to  the rhs:
        ! The first row contains surface forcing
        nz    = nzmin
        zinv  = 1.0_WP*dt/helem(nz, elem)
        ur(nz)= ur(nz) + zinv*stress_surf(1, elem)/density_0
        vr(nz)= vr(nz) + zinv*stress_surf(2, elem)/density_0
        if (dynamics%ldiag_ke) then
           dynamics%ke_wind(1,elem)=stress_surf(1, elem)/density_0*dt
           dynamics%ke_wind(2,elem)=stress_surf(2, elem)/density_0*dt
        end if
        
        ! The last row contains bottom friction
        nz    = nzmax-1
        zinv  = 1.0_WP*dt/helem(nz, elem)
        friction=-C_d*sqrt(UV(1, nz, elem)**2 + UV(2, nz, elem)**2)
        ur(nz)= ur(nz) + zinv*friction*UV(1, nz, elem)
        vr(nz)= vr(nz) + zinv*friction*UV(2, nz, elem)
        if (dynamics%ldiag_ke) then
           dynamics%ke_drag(1, elem)=friction*UV(1, nz, elem)*dt
           dynamics%ke_drag(2, elem)=friction*UV(2, nz, elem)*dt
        end if
        
        !_______________________________________________________________________
        ! The sweep algorithm --> solve for du_i
        ! initialize c-prime and s,t-prime
        cp(nzmin) =  c(nzmin)/b(nzmin)
        up(nzmin) = ur(nzmin)/b(nzmin)
        vp(nzmin) = vr(nzmin)/b(nzmin) 
        
        ! sove for vectors c-prime and t, s-prime
        do nz=nzmin+1, nzmax-1
            m      = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            up(nz) = (ur(nz)-up(nz-1)*a(nz))/m
            vp(nz) = (vr(nz)-vp(nz-1)*a(nz))/m
        end do
    
        ! initialize x
        ur(nzmax-1) = up(nzmax-1)
        vr(nzmax-1) = vp(nzmax-1)
        
        ! solve for x from the vectors c-prime and d-prime
        do nz=nzmax-2, nzmin, -1
            ur(nz) = up(nz)-cp(nz)*ur(nz+1)
            vr(nz) = vp(nz)-cp(nz)*vr(nz+1)
        end do
        
        !_______________________________________________________________________
        ! RHS update
        ! ur and vr are in terms of du velocity difference, UV_rhs are in terms 
        ! of transport velocity
        do nz=nzmin, nzmax-1
            UV_rhs(1, nz, elem)=UV_rhs(1, nz, elem)+ur(nz)*helem(nz, elem)             !! The rhs is for transport
            UV_rhs(2, nz, elem)=UV_rhs(2, nz, elem)+vr(nz)*helem(nz, elem)
        end do
        
        !_______________________________________________________________
     
    end do ! --> do elem=1,myDim_elem2D
!$OMP END DO
!$OMP END PARALLEL
end subroutine impl_vert_visc_ale_vtransp
!
!
!_______________________________________________________________________________
!SD Compute vertical integral of transport velocity rhs omitting the contributions from
!SD the elevation and Coriolis. The elevation and Coriolis are accounted for 
!SD explicitly in BT equations, and should therefore be removed from the vertically integrated rhs.
subroutine compute_BT_rhs_SE_vtransp(dynamics, partit, mesh)
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE MOD_DYN
    USE g_config, only: dt, r_restart
    USE g_comm_auto
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout)   , target :: mesh
    !___________________________________________________________________________
    real(kind=WP)       :: vert_sum_u, vert_sum_v, Fx, Fy, ab1, ab2, hh
    integer             :: elem, nz, nzmin, nzmax, elnodes(3), ed, el(2)
    logical, save       :: sfirst
    real(kind=WP)       :: update_ubt, update_vbt, vi, len
    
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV_rhs
    real(kind=WP), dimension(:)    , pointer :: eta_n
    real(kind=WP), dimension(:,:)  , pointer :: UVBT_4AB, UVBT_rhs
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    eta_n     =>dynamics%eta_n(:)
    UV_rhs    =>dynamics%uv_rhs(:,:,:)
    UVBT_rhs  =>dynamics%se_uvBT_rhs(:,:)
    UVBT_4AB  =>dynamics%se_uvBT_4AB(:,:)
    
    !___________________________________________________________________________
    ab1=(1.5_WP+epsilon)
    if(sfirst .and. (.not. r_restart)) then
        ab1=1.0_WP
        sfirst=.false.
    end if
    ab2=1.0_WP-ab1

    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem, elnodes, nz, nzmin, nzmax, &
!$OMP                                  Fx, Fy, hh, vert_sum_u, vert_sum_v)
!$OMP DO
    do elem=1,  myDim_elem2D
        elnodes= elem2D_nodes(:,elem)
        nzmin  = ulevels(elem)
        nzmax  = nlevels(elem)-1
        
        !_______________________________________________________________________
        ! vertically integrate UV_rhs --> for barotropic equatiobn 
        vert_sum_u=0.0_WP
        vert_sum_v=0.0_WP
        do nz=nzmin, nzmax
            vert_sum_u=vert_sum_u + UV_rhs(1, nz, elem)   
            vert_sum_v=vert_sum_v + UV_rhs(2, nz, elem)   
        end do
        
        !_______________________________________________________________________
        ! Remove the contribution from the elevation will be accounted explicitely 
        ! for in the barotropic equation
        Fx = g*dt*sum(gradient_sca(1:3,elem)*eta_n(elnodes))
        Fy = g*dt*sum(gradient_sca(4:6,elem)*eta_n(elnodes))
        
        ! total ocean depth H
        !PS hh = sum(helem(nzmin:nzmax, elem))
        !PS hh = -zbar_e_bot(elem) + sum(eta_n(elnodes))/3.0_WP
        hh = -zbar_e_bot(elem) + sum(eta_n(elnodes))/3.0_WP
        
        vert_sum_u=vert_sum_u + Fx*hh
        vert_sum_v=vert_sum_v + Fy*hh
        
        !_______________________________________________________________________
        ! Remove the contribution from the Coriolis will be accounted explicitely 
        ! for in the barotropic equation
        ! UVBT_rhs ... baroclinic forcing term in barotropic equation R_b
        ! --> d/dt*U_bt + f*e_z x U_bt + g*H* grad(eta) = R_bt
        ! --> from AB2 in oce_ale_vel_rhs.F90 --> ab2=-0.5 (from previouse time
        !     step) --> ab1=1.5 or 0.0 (first tstep)
        UVBT_rhs(1, elem)=vert_sum_u - dt*mesh%coriolis(elem)* &
                                       (ab1*UVBT_4AB(2,elem)+ab2*UVBT_4AB(4,elem))         ! for AB-interpolated
        UVBT_rhs(2, elem)=vert_sum_v + dt*mesh%coriolis(elem)* &
                                       (ab1*UVBT_4AB(1,elem)+ab2*UVBT_4AB(3,elem))    
        
        !_______________________________________________________________________
        ! save actual of vertical integrated transport velocity UVBT_4AB(1:2,elem)  
        ! timestep for next time loop step
        ! UVBT_4AB ... is U_bt in barotropic equation
        ! --> d/dt*U_bt + f*e_z x U_bt + g*H*grad_H(eta) = R_bt
        UVBT_4AB(3:4,elem)=UVBT_4AB(1:2,elem)    
        
    end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine compute_BT_rhs_SE_vtransp
!
!
!_______________________________________________________________________________
! Barotropic time stepping with Forward-Backward dissipative method 
! (Demange et al. 2019) is used where eta and U_bt = sum_k(U_k) = sum_k(u_k*h_k) 
! are estimated from the equations ...
! --> d/dt*U_bt + f*e_z x U_bt + g*H*grad_H(eta) = R_bt
! --> d/dt*eta + div_h(U_bt) + FW = 0
!
!
! --> forward-backward dissipative time stepping by Demagne et al. 2019
! --> equation (6) in T. Banerjee et al.,Split-Explicite external
!     mode solver in FESOM2, 
!     
! Ubt^(n+(m+1)/M) = Ubt^(n+(m)/M) - dt/M*[ 
!                                  + 0.5*f*e_z x (Ubt^(n+(m+1)/M) + Ubt^(n+(m)/M))
!                                  - h*H^m*grad_H*eta^((n+m)/M)
!                                  - Rbt-->UVBT_rhs ] 
! 
! eta^(n+(m+1)/M) = eta^(n+(m)/M) - dt/M * div_H * [(1+theta)*Ubt^(n+(m+1)/M) - theta*Ubt^(n+(m)/M)]
!
subroutine compute_BT_step_SE_ale(dynamics, partit, mesh)
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE MOD_DYN
    USE g_comm_auto
    USE g_config,  only: dt, which_ALE, use_cavity_fw2press
    USE g_support, only: integrate_nod
    use o_ARRAYS,  only: water_flux
    IMPLICIT NONE
    !___________________________________________________________________________
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout)   , target :: mesh
    !___________________________________________________________________________
    real(kind=WP)                   :: dtBT, BT_inv, hh, ff, rx, ry, a, b, d, c1, c2, ax, ay
    real(kind=WP)                   :: deltaX1, deltaY1, deltaX2, deltaY2, thetaBT
    integer                         :: step, elem, edge, node, elnodes(3), ednodes(2), edelem(2), nzmax
    real(kind=WP)                   :: update_ubt, update_vbt, vi, len
    
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:)    , pointer :: eta_n, bottomdrag
    real(kind=WP), dimension(:,:)  , pointer :: UVBT_rhs, UVBT, UVBT_theta, UVBT_mean, UVBT_12, UVBT_harmvisc
    real(kind=WP), dimension(:,:,:), pointer :: UV
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    eta_n     =>dynamics%eta_n(:)
    UVBT_rhs  =>dynamics%se_uvBT_rhs(:,:)
    UVBT      =>dynamics%se_uvBT(:,:)
    UVBT_theta=>dynamics%se_uvBT_theta(:,:)
    UVBT_mean =>dynamics%se_uvBT_mean(:,:)
    UVBT_12   =>dynamics%se_uvBT_12(:,:)
    if (dynamics%se_bottdrag) then
        UV           =>dynamics%uv(:,:,:)
        bottomdrag   =>dynamics%se_uvBT_stab_bdrag(:)
    end if 
    if (dynamics%se_visc) then 
        UVBT_harmvisc=>dynamics%se_uvBT_stab_hvisc(:,:)
    end if 
    
    !___________________________________________________________________________
    ! Dissipation parameter of FB dissipative method 0.14 is the default value 
    ! from Demange et al.
    thetaBT= dynamics%se_BTtheta
    
    ! BTsteps should be 30 or 40.
    dtBT   = dt/dynamics%se_BTsteps                
    BT_inv = 1.0_WP/(1.0_WP*dynamics%se_BTsteps)
    
    !___SPLIT-EXPLICITE STABILIZATION___________________________________________
    ! trim R (UVBT_rhs):
    ! UVBT_rhs --> UVBT_rhs - div_h Ah H^n div_h(Ubt^n/H^n) + Cd*|Ubot|* Ubt^n/H^n 
    ! The intention here is to approximately remove the term thatwill be
    ! added later to the barotropic momentum equation (see subroutine 
    ! compute_BT_step_SE_ale)
    ! --> use only harmonmic viscosity operator applied to the barotropic
    !     velocity
    if (dynamics%se_visc) then 
        !_______________________________________________________________________
        ! remove viscosity
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, edelem, ednodes, hh, len, &
!$OMP                                  vi, update_ubt, update_vbt)
!$OMP DO        
        do edge=1, myDim_edge2D+eDim_edge2D
                
            ! if ed is an outer boundary edge, skip it
            if(myList_edge2D(edge)>edge2D_in) cycle
                
            ! elem indices that participate in edge
            edelem  = edge_tri(:,edge)
            ednodes = edges(:,edge) 
            
            ! total ocean depth H
            !PS nzmax = minval(nlevels(edelem))
            !PS hh    = -zbar(nzmax)
            !PS hh      = minval(-zbar_e_bot(edelem))
            hh    = -sum(zbar_e_bot(edelem))*0.5_WP + sum(hbar(ednodes))*0.5_WP
            
            len     = sqrt(sum(elem_area(edelem)))
            update_ubt=(UVBT(1, edelem(1))-UVBT(1, edelem(2)))/hh
            update_vbt=(UVBT(2, edelem(1))-UVBT(2, edelem(2)))/hh
            vi=update_ubt*update_ubt + update_vbt*update_vbt
            vi=-dt*sqrt(max(dynamics%se_visc_gamma0,           &
                        max(dynamics%se_visc_gamma1*sqrt(vi),  &
                            dynamics%se_visc_gamma2*vi)        &
                    )*len)
            update_ubt=update_ubt*vi
            update_vbt=update_vbt*vi
            
            !___________________________________________________________________
#if defined(_OPENMP) && !defined(__openmp_reproducible)
            call omp_set_lock(partit%plock(edelem(1)))
#else
!$OMP ORDERED
#endif
            UVBT_rhs(1, edelem(1))=UVBT_rhs(1, edelem(1))-update_ubt/elem_area(edelem(1))*hh
            UVBT_rhs(2, edelem(1))=UVBT_rhs(2, edelem(1))-update_vbt/elem_area(edelem(1))*hh
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(edelem(1)))
            call omp_set_lock  (partit%plock(edelem(2)))
#endif
            UVBT_rhs(1, edelem(2))=UVBT_rhs(1, edelem(2))+update_ubt/elem_area(edelem(2))*hh
            UVBT_rhs(2, edelem(2))=UVBT_rhs(2, edelem(2))+update_vbt/elem_area(edelem(2))*hh
#if defined(_OPENMP) && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(edelem(2)))
#else
!$OMP END ORDERED
#endif
        end do ! --> do edge=1, myDim_edge2D+eDim_edge2D
!$OMP END DO       
!$OMP END PARALLEL 
    end if ! --> if (dynamics%se_visc) then     
    
    !___________________________________________________________________________
    ! remove bottom drag
    if (dynamics%se_bottdrag) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem, elnodes, nzmax, hh)
!$OMP DO    
        do elem=1, myDim_elem2D
            elnodes= elem2D_nodes(:,elem)
            nzmax  = nlevels(elem)
            
            ! total ocean depth H
            !PS hh     = -zbar(nzmax)+sum(eta_n(elnodes))/3.0_WP
            !PS hh     = -zbar(nzmax)
            !PS hh     = -zbar_e_bot(elem)
            hh     = -zbar_e_bot(elem) + sum(hbar(elnodes))/3.0_WP
            
            bottomdrag(elem) = dt*C_d*sqrt(UV(1, nzmax-1, elem)**2 + UV(2, nzmax-1, elem)**2)
            UVBT_rhs(1, elem)=UVBT_rhs(1, elem) + bottomdrag(elem)*UVBT(1, elem)/hh
            UVBT_rhs(2, elem)=UVBT_rhs(2, elem) + bottomdrag(elem)*UVBT(2, elem)/hh
        end do
!$OMP END DO       
!$OMP END PARALLEL        
    end if ! --> if (dynamics%se_bottdrag) then
    
    !___________________________________________________________________________
    ! initialise UVBT_mean with zeros --> OMP style
!$OMP PARALLEL DO
    do elem=1, myDim_elem2D+eDim_elem2D
        UVBT_mean(:, elem) = 0.0_WP
    end do
!$OMP END PARALLEL DO    
   
    !___________________________________________________________________________
    ! eta_n   elevation used in BT stepping, it is just a copy of eta_n
    ! UBT and VBT are transport velocities
    do step=1, dynamics%se_BTsteps
        !#######################################################################
        !##########    Dissipative forward--backward time stepping    ##########
        !#######################################################################
        
        !_______________________________________________________________________
        ! compute harmonic viscosity for stability
        if (dynamics%se_visc) then 
            
            !___________________________________________________________________
            ! initialise UVBT_harmvisc with zeros --> OMP style
!$OMP PARALLEL DO
            do elem=1, myDim_elem2D+eDim_elem2D
                UVBT_harmvisc(:, elem) = 0.0_WP
            end do
!$OMP END PARALLEL DO
            
            !___________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, edelem, ednodes, hh, len, &
!$OMP                                  vi, update_ubt, update_vbt)
!$OMP DO                  
            do edge=1, myDim_edge2D+eDim_edge2D
                    
                ! if ed is an outer boundary edge, skip it
                if(myList_edge2D(edge)>edge2D_in) cycle
                    
                ! elem indices that participate in edge
                edelem  = edge_tri(:, edge)
                ednodes = edges(:, edge)
                
                ! total ocean depth H
                !PS nzmax = minval(nlevels(edelem))
                !PS hh    = -zbar(nzmax)
                !PS hh    = minval(-zbar_e_bot(edelem))
                hh      = -sum(zbar_e_bot(edelem))*0.5_WP + sum(hbar(ednodes))*0.5_WP
                
                len     = sqrt(sum(elem_area(edelem)))
                update_ubt=(UVBT(1, edelem(1))-UVBT(1, edelem(2)))/hh
                update_vbt=(UVBT(2, edelem(1))-UVBT(2, edelem(2)))/hh
                vi=update_ubt*update_ubt + update_vbt*update_vbt
                vi=dt*sqrt(max(dynamics%se_visc_gamma0,           &
                           max(dynamics%se_visc_gamma1*sqrt(vi),   &
                               dynamics%se_visc_gamma2*vi)         &
                        )*len)
                update_ubt=update_ubt*vi
                update_vbt=update_vbt*vi
                
                !_______________________________________________________________
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                call omp_set_lock(partit%plock(edelem(1)))
#else
!$OMP ORDERED
#endif
                UVBT_harmvisc(1, edelem(1))=UVBT_harmvisc(1, edelem(1))-update_ubt/elem_area(edelem(1))*hh
                UVBT_harmvisc(2, edelem(1))=UVBT_harmvisc(2, edelem(1))-update_vbt/elem_area(edelem(1))*hh
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(edelem(1)))
                call omp_set_lock  (partit%plock(edelem(2)))
#endif
                UVBT_harmvisc(1, edelem(2))=UVBT_harmvisc(1, edelem(2))+update_ubt/elem_area(edelem(2))*hh
                UVBT_harmvisc(2, edelem(2))=UVBT_harmvisc(2, edelem(2))+update_vbt/elem_area(edelem(2))*hh
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(edelem(2)))
#else
!$OMP END ORDERED
#endif
            end do ! --> do edge=1, myDim_edge2D+eDim_edge2D
!$OMP END DO       
!$OMP END PARALLEL 
        end if ! -> if (dynamics%se_visc) then 
        
        !_______________________________________________________________________
        ! Advance velocities. I use SI stepping for the Coriolis
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem, elnodes, hh, ff, rx, ry, a, b,  &
!$OMP                                  d, ax, ay)
!$OMP DO           
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            !___________________________________________________________________
            ! compute term: 
            ! AAA = - dt/M*[ + 0.5*f*e_z x (Ubt^(n+(m+1)/M) + Ubt^(n+(m)/M))
            !          - h*H^m*grad_H*eta^((n+m)/M)
            !          - Rbt-->UVBT_rhs ] 
            ! total ocean depth H
            !PS hh = -zbar(nlevels(elem))+sum(eta_n(elnodes))/3.0_WP ! Total fluid depth
            !PS hh = -zbar(nlevels(elem)) ! Total fluid depth
            !PS hh = -zbar_e_bot(elem)
            hh = -zbar_e_bot(elem) + sum(hbar(elnodes))/3.0_WP
            
            ff = mesh%coriolis(elem)
            
            rx =   dtBT*(-g*hh*sum(gradient_sca(1:3,elem)*eta_n(elnodes)) + ff*UVBT(2, elem))  & 
                 + BT_inv*UVBT_rhs(1, elem)                                                   !& 
!PS                  + BT_inv*(UVBT_harmvisc(1, elem) - bottomdrag(elem)*UVBT(1, elem)/hh) ! <-- stabilization terms
                 
            ry =   dtBT*(-g*hh*sum(gradient_sca(4:6,elem)*eta_n(elnodes)) - ff*UVBT(1, elem))  &
                 + BT_inv*UVBT_rhs(2, elem)                                                   !&
!PS                  + BT_inv*(UVBT_harmvisc(2, elem) - bottomdrag(elem)*UVBT(2, elem)/hh) ! <-- stabilization terms
            if (dynamics%se_visc) then 
                rx = rx + BT_inv*UVBT_harmvisc(1, elem) 
                ry = ry + BT_inv*UVBT_harmvisc(2, elem) 
            end if 
            if (dynamics%se_bottdrag) then
                rx = rx - BT_inv*bottomdrag(elem)*UVBT(1, elem)/hh
                ry = ry - BT_inv*bottomdrag(elem)*UVBT(2, elem)/hh
            end if 
            
            ! compute new velocity Ubt^(n+(m+1)/M), Vbt^(n+(m+1)/M) considering 
            ! in terms of Increments (deltaU) and semi-Implicit Coriolis
            ! (We do it here based on increments since it saves us some significant
            ! digits for the accuracy)
            ! Increments:
            ! deltaU                       = Ubt^(n+(m+1)/M)-Ubt^(n+m/M)
            ! Ubt^(n+(m+1)/M)+Ubt^(n+m/M)  = Ubt^(n+(m+1)/M)-Ubt^(n+m/M)+2*Ubt^(n+m/M)
            !                              = deltaU + 2*Ubt^(n+m/M)
            !
            ! Ubt^(n+(m+1)/M)-Ubt^(n+m/M)  = - dt/M*[ + 0.5*f*e_z x (Ubt^(n+(m+1)/M) + Ubt^(n+m/M))
            !                                - h*H^m*grad_H*eta^(n+m/M)
            !                                - Rbt-->UVBT_rhs ] 
            !                                + dt/M*[grad_h*A_h*H^n*grad_h(Ubt^(n+m/M)/H^n)]
            !                                - dt/M*[cd*|Ubot|*Ubt^(n+(m+1)/M)/H^n]
            !
            ! --> a = dt/(2*M)*ff
            ! --> b = dt/M*cd*|Ubot|/H^n  
            !
            ! deltaU + b*deltaU - a*deltaV   =  dt/M*f*Vbt^(n+m/M) - h*H^m*gradx_H*eta^(n+m/M) + Rbtx
            !                                   + dt/M*[grad_h*A_h*H^n*grad_h(Ubt^(n+m/M)/H^n)]
            !                                   + dt/M*[cd*|Ubot|*Ubt^(n+m/M)/H^n]
            ! deltaV + b*deltaV + a*deltaU   = -dt/M*f*Ubt^(n+m/M) - h*H^m*grady_H*eta^(n+m/M) + Rbty
            !                                   + dt/M*[grad_h*A_h*H^n*grad_h(Vbt^(n+m/M)/H^n)]
            !                                   + dt/M*[cd*|Ubot|VUbt^(n+m/M)/H^n]
            !                                \________________________v___________________________/
            !                                                      Rx, Ry 
            !
            ! | 1+b -a  | * | deltaU | = MAT* | deltaU | = | Rx |
            ! |  a  1+b |   | deltaV | =      | deltaV | = | Ry |
            ! --> b' = b+1  
            ! --> d  = 1/(b'^2 + a^2) ; inv(MAT) = d* |  b'  a  |
            !                                         | -a   b' |
            !
            ! | deltaU | = inv(MAT) * | Rx | 
            ! | deltaV | =            | Ry |
            !
            ! Ubt^(n+(m+1)/M) = Ubt^(n+m/M) + d*( b'*Rx + a *Ry) 
            ! Vbt^(n+(m+1)/M) = Vbt^(n+m/M) + d*(-a *Rx + b'*Ry) 
            !
            ! Semi-Implicit Coriolis
            a  = dtBT*ff*0.5_WP
            if ( (dynamics%se_bdrag_si) .and. (dynamics%se_bottdrag) ) then 
                b  = 1.0_WP+BT_inv*bottomdrag(elem)/hh
            else
                b  = 1.0_WP
            end if 
            d  = 1.0_WP/(b*b + a*a)
            ax = d*(  b*rx + a*ry )
            ay = d*( -a*rx + b*ry )
            
            !___________________________________________________________________
            ! compute new velocities Ubt^(n+(m+1)/M) at barotropic time step (n+(m+1)/M) ...
            ! Ubt^(n+(m+1)/M)   = Ubt^(n+(m)/M) + AAA
            ! equation (6) in T. Banerjee et al.,Split-Explicite external
            ! mode solver in FESOM2, 
            ! compute barotropic velocity at time step (n+(m+1)/M)
            UVBT(      1, elem) = UVBT(1, elem) + ax   
            UVBT(      2, elem) = UVBT(2, elem) + ay
            
            !___________________________________________________________________
            ! velocities for dissipative time stepping of thickness equation
            ! compute: [(1+theta)*Ubt^(n+(m+1)/M) - theta*Ubt^(n+(m)/M)] by ...
            ! 
            ! --> AAA = Ubt^(n+(m+1)/M) - Ubt^(n+(m)/M)    | *theta
            !     AAA * theta = (Ubt^(n+(m+1)/M) - Ubt^(n+(m)/M) )*theta
            ! --> Ubt^(n+(m+1)/M) + AAA*theta = Ubt^(n+(m+1)/M) + (Ubt^(n+(m+1)/M) - Ubt^(n+(m)/M) )*theta
            !     (1+theta)*Ubt^(n+(m+1)/M) + theta*Ubt^(n+(m)/M)
            UVBT_theta(1, elem) = UVBT(1, elem) + thetaBT*ax   ! New vel*(1+thetaBT)-old vel *thetaBT
            UVBT_theta(2, elem) = UVBT(2, elem) + thetaBT*ay
            
            !___________________________________________________________________
            ! Mean BT velocity to trim 3D velocity in tracers, equation (10) in 
            ! T. Banerjee et al.,Split-Explicite external mode solver in FESOM2, 
            UVBT_mean( 1, elem) = UVBT_mean( 1, elem) + UVBT_theta(1, elem)*BT_inv
            UVBT_mean( 2, elem) = UVBT_mean( 2, elem) + UVBT_theta(2, elem)*BT_inv
            
        end do
!$OMP END DO       
!$OMP END PARALLEL 
        
        !_______________________________________________________________________
!$OMP MASTER
        call exchange_elem_begin(UVBT, partit)
        call exchange_elem_end(partit)
!$OMP END MASTER
!$OMP BARRIER

        !_______________________________________________________________________
        ! Store mid-step velocity (to trim 3D velocities in momentum)
        if(step==dynamics%se_BTsteps/2) then
            UVBT_12=UVBT
        end if 
        
        !_______________________________________________________________________
        ! Advance thickness
        ! compute: dt/M * div_H * [(1+theta)*Ubt^(n+(m+1)/M) - theta*Ubt^(n+(m)/M)]
        ! and advance ssh --> eta^(n+(m+1)/M) = eta^(n+(m)/M) - dt/M * div_H * [...]
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, ednodes, edelem, c1, c2, &
!$OMP                                  deltaX1, deltaX2, deltaY1, deltaY2)
!$OMP DO
        do edge=1, myDim_edge2D
            ednodes = edges(:,edge)
            edelem  = edge_tri(:,edge)
            
            !___________________________________________________________________
            ! compute divergence div_H * [(1+theta)*Ubt^(n+(m+1)/M) - theta*Ubt^(n+(m)/M)]
            deltaX1 = edge_cross_dxdy(1,edge)
            deltaY1 = edge_cross_dxdy(2,edge)
            c1 = UVBT_theta(2, edelem(1))*deltaX1 - UVBT_theta(1, edelem(1))*deltaY1
            c2 = 0.0_WP
            if(edelem(2)>0) then
                deltaX2=edge_cross_dxdy(3,edge)
                deltaY2=edge_cross_dxdy(4,edge)
                c2=-(UVBT_theta(2, edelem(2))*deltaX2 - UVBT_theta(1, edelem(2))*deltaY2)
            end if
            
            !___________________________________________________________________
            ! advance ssh --> eta^(n+(m+1)/M) = eta^(n+(m)/M) - dt/M * div_H * [...]
            ! equation (6) in T. Banerjee et al.,Split-Explicite external
            ! mode solver in FESOM2, 
#if defined(_OPENMP) && !defined(__openmp_reproducible)
            call omp_set_lock(partit%plock(ednodes(1)))
#else
!$OMP ORDERED
#endif            
            eta_n(ednodes(1))=eta_n(ednodes(1)) + (c1+c2)*dtBT/areasvol(1,ednodes(1))
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(ednodes(1)))
            call omp_set_lock  (partit%plock(ednodes(2)))
#endif            
            eta_n(ednodes(2))=eta_n(ednodes(2)) - (c1+c2)*dtBT/areasvol(1,ednodes(2))
#if defined(_OPENMP) && !defined(__openmp_reproducible) 
            call omp_unset_lock(partit%plock(ednodes(2)))
#else
!$OMP END ORDERED
#endif            
        end do
!$OMP END DO       
!$OMP END PARALLEL 

        !_______________________________________________________________________
        ! Apply freshwater boundary condition 
        if ( .not. trim(which_ALE)=='linfs') then
!$OMP PARALLEL DO
            do node=1,myDim_nod2D
                eta_n(node)=eta_n(node) - dtBT*water_flux(node)
            end do
!$OMP END PARALLEL DO            
        end if 
        
        !_______________________________________________________________________
!$OMP MASTER        
        call exchange_nod(eta_n, partit)
!$OMP END MASTER
!$OMP BARRIER

    end do ! --> do step=1, dynamics%se_BTsteps
    
    !___________________________________________________________________________
    hbar_old = hbar
    hbar     = eta_n
end subroutine compute_BT_step_SE_ale
!
!
!_______________________________________________________________________________
! Trim U and Uh to be consistent with BT transport
subroutine update_trim_vel_ale_vtransp(mode, dynamics, partit, mesh)
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE MOD_DYN
    use g_comm
    IMPLICIT NONE
    !___________________________________________________________________________
    integer       , intent(in)            :: mode
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout)   , target :: mesh
    integer                               :: elem, nz, nzmin, nzmax
    real(kind=WP)                         :: ubar, vbar, hh_inv
    real(kind=WP)                         :: usum(2), udiff(2)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UVh, UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: UVBT_mean, UVBT_12
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV        =>dynamics%uv(:,:,:)
    UVh       =>dynamics%se_uvh(:,:,:)
    UV_rhs    =>dynamics%uv_rhs(:,:,:)
    UVBT_mean =>dynamics%se_uvBT_mean(:,:)
    UVBT_12   =>dynamics%se_uvBT_12(:,:)
    
    !___________________________________________________________________________
    !
    ! Transport version
    !
    !___________________________________________________________________________
    ! Trim full velocity to ensure that its vertical sum equals to the mean BT velocity
    ! The trimmed Uh,Vh are consistent with new total height defined by eta_n   
    if (mode==1) then 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax, ubar, vbar, hh_inv)
!$OMP DO        
        do elem=1, myDim_elem2D
            nzmin  = ulevels(elem)
            nzmax  = nlevels(elem)-1
            
            !___________________________________________________________________
            ubar   = 0.0_WP
            vbar   = 0.0_WP
            hh_inv = 0.0_WP
            do nz=nzmin, nzmax-1
                !PS !___________________________________________________________
                !PS ! Update transport velocities: U^(n+1/2,**) = U^(n+1/2,*) + U_rhs
                !PS UVh(1, nz, elem)=UVh(1, nz, elem)+UV_rhs(1, nz, elem)   
                !PS UVh(2, nz, elem)=UVh(2, nz, elem)+UV_rhs(2, nz, elem)
                !PS 
                !PS !___________________________________________________________
                !PS ! vertically integrate updated transport velocity: sum(k, U_k^(n+1/2,**) )
                !PS ubar  = ubar+UVh(1, nz, elem)
                !PS vbar  = vbar+UVh(2, nz, elem)
                !PS hh_inv= hh_inv+helem(nz,elem)
                !_______________________________________________________________
                ! vertically integrate updated transport velocity: sum(k, U_k^(n+1/2,**) )
                ! --> the actual update of the transport velocity is done after the 
                !     if (dynamics%ldiag_ke) block
                ubar  = ubar+UVh(1, nz, elem) + UV_rhs(1, nz, elem)   
                vbar  = vbar+UVh(2, nz, elem) + UV_rhs(2, nz, elem)   
                hh_inv= hh_inv+helem(nz,elem)
            end do ! --> do nz=nzmin, nzmax
            ! inverse total heigh over which the trimming part is distributed -->
            ! here exclude bottom cell
            hh_inv=1.0_WP/hh_inv
            
            ! sum transport over entire water column including bottom cell need 
            ! the trimming with the barotropic transport
            nz    = nzmax
            ubar  = ubar+UVh(1, nz, elem) + UV_rhs(1, nz, elem)   
            vbar  = vbar+UVh(2, nz, elem) + UV_rhs(2, nz, elem)
            
            !___________________________________________________________________
            if (dynamics%ldiag_ke) then
                do nz=nzmin, nzmax
                    ! U_(n+1)   - U_n   =  Urhs_n |* (U_(n+1)+U_n)
                    ! U_(n+1)^2 - U_n^2 =  Urhs_n * (U_(n+1)+U_n) 
                    !                                  | 
                    !                                  +-> U_(n+1) = U_n+Urhs_n
                    ! U_(n+1)^2 - U_n^2 =  Urhs_n * (2*U_n + Urhs) 
                    !                          |            |
                    !                          v            v 
                    !                        udiff         usum 
                    if (nz==nzmax) then 
                        usum(1)  = 2.0_WP*UVh(1,nz,elem) + UV_rhs(1, nz, elem)
                        usum(2)  = 2.0_WP*UVh(2,nz,elem) + UV_rhs(2, nz, elem)
                    else
                        usum(1)  = 2.0_WP*UVh(1,nz,elem) + UV_rhs(1, nz, elem) + (UVBT_mean(1, elem)-ubar)*helem(nz,elem)*hh_inv
                        usum(2)  = 2.0_WP*UVh(2,nz,elem) + UV_rhs(2, nz, elem) + (UVBT_mean(2, elem)-vbar)*helem(nz,elem)*hh_inv
                    end if     
                    ! transform: transport vel --> velocity
                    usum(1)  = usum(1)/helem(nz,elem)
                    usum(2)  = usum(2)/helem(nz,elem)
                    
                    udiff(1) = UV_rhs(1, nz, elem) + (UVBT_mean(1, elem)-ubar)*helem(nz,elem)*hh_inv
                    udiff(2) = UV_rhs(2, nz, elem) + (UVBT_mean(2, elem)-vbar)*helem(nz,elem)*hh_inv
                    ! transform: transport vel --> velocity
                    udiff(1) = udiff(1)/helem(nz,elem)
                    udiff(2) = udiff(2)/helem(nz,elem)
                    
                    ! (U_(n+1)^2 - U_n^2)/2 = usum*udiff/2
                    dynamics%ke_du2(      :,nz,elem)  = usum*udiff/2.0_WP
                    
                    dynamics%ke_pre_xVEL( :,nz,elem)  = usum*dynamics%ke_pre (:,nz,elem)/2.0_WP
                    dynamics%ke_adv_xVEL( :,nz,elem)  = usum*dynamics%ke_adv (:,nz,elem)/2.0_WP
                    dynamics%ke_cor_xVEL( :,nz,elem)  = usum*dynamics%ke_cor (:,nz,elem)/2.0_WP
                    dynamics%ke_hvis_xVEL(:,nz,elem)  = usum*dynamics%ke_hvis(:,nz,elem)/2.0_WP
                    dynamics%ke_vvis_xVEL(:,nz,elem)  = usum*dynamics%ke_vvis(:,nz,elem)/2.0_WP
                    
                    ! U_(n+0.5)    = U_n + 0.5*Urhs
                    dynamics%ke_umean(    :,nz,elem)  = usum/2.0_WP
                    ! U_(n+0.5)^2
                    dynamics%ke_u2mean(   :,nz,elem)  = (usum*usum)/4.0_WP
                    
                    if (nz==nzmin) then
                        dynamics%ke_wind_xVEL(:,elem) = usum*dynamics%ke_wind(:,elem)/2.0_WP
                    end if
                    
                    if (nz==nzmax) then
                        dynamics%ke_drag_xVEL(:,elem) = usum*dynamics%ke_drag(:,elem)/2.0_WP
                    end if
                end do ! --> do nz=nzmin, nzmax
            end if ! --> if (dynamics%ldiag_ke) then
            
            !___________________________________________________________________
            ! finalize horizontal transport by making vertically integrated 
            ! transport equal to the value obtained from the barotropic solution
            ! 
            ! --> equation (11) in T. Banerjee et al.,Split-Explicite external
            !     mode solver in FESOM2, 
            ! U_k^(n+1/2) =  U^(n+1/2,**) 
            !              + [<<Ubar>>^(n+1/2) - sum(k, U_k^(n+1/2,**)] * h_k^(n+1/2)/sum(k,h_k^(n+1/2))
            do nz=nzmin, nzmax-1
                !PS UVh(1, nz, elem)= UVh(1, nz, elem)+(UVBT_mean(1, elem)-ubar)*helem(nz,elem)*hh_inv
                !PS UVh(2, nz, elem)= UVh(2, nz, elem)+(UVBT_mean(2, elem)-vbar)*helem(nz,elem)*hh_inv
                UVh(1, nz, elem)= UVh(1, nz, elem) + UV_rhs(1, nz, elem) + (UVBT_mean(1, elem)-ubar)*helem(nz,elem)*hh_inv
                UVh(2, nz, elem)= UVh(2, nz, elem) + UV_rhs(2, nz, elem) + (UVBT_mean(2, elem)-vbar)*helem(nz,elem)*hh_inv
                UV( 1, nz, elem)= UVh(1, nz, elem)/helem(nz,elem)  ! velocities are still needed    
                UV( 2, nz, elem)= UVh(2, nz, elem)/helem(nz,elem) 
            end do ! --> do nz=nzmin, nzmax
            
            ! apply no trimming in the bottom cell
            nz = nzmax
            UVh(1, nz, elem)= UVh(1, nz, elem) + UV_rhs(1, nz, elem)
            UVh(2, nz, elem)= UVh(2, nz, elem) + UV_rhs(2, nz, elem)
            UV( 1, nz, elem)= UVh(1, nz, elem)/helem(nz,elem)  ! velocities are still needed    
            UV( 2, nz, elem)= UVh(2, nz, elem)/helem(nz,elem) 
            
        end do ! --> do elem=1, myDim_elem2D
!$OMP END DO       
!$OMP END PARALLEL    
        
        !_______________________________________________________________________
!$OMP MASTER        
        call exchange_elem(UVh, partit)  ! This exchange can be avoided, but test first.
        call exchange_elem(UV, partit)
!$OMP END MASTER
!$OMP BARRIER

    !___________________________________________________________________________
    ! Trim to make velocity consistent with BT velocity at n+1/2
    ! Velocity will be used to advance momentum
    elseif (mode==2) then 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax, ubar, vbar, hh_inv)
!$OMP DO      
        do elem=1, myDim_elem2D    
            nzmin  = ulevels(elem)
            nzmax  = nlevels(elem)-1
            
            !___________________________________________________________________
            ubar   = 0.0_WP
            vbar   = 0.0_WP
            hh_inv = 0.0_WP    
            do nz=nzmin, nzmax-1
                ubar=ubar+UVh(1, nz, elem)
                vbar=vbar+UVh(2, nz, elem)
                hh_inv=hh_inv+helem(nz,elem)
            end do
            nz=nzmax
            ubar=ubar+UVh(1, nz, elem)
            vbar=vbar+UVh(2, nz, elem)
            
            !___________________________________________________________________
            hh_inv=1.0_WP/hh_inv    ! Postpone 2nd order, just use available thickness 
            do nz=nzmin, nzmax-1
                UVh(1, nz, elem)= UVh(1, nz, elem)+(UVBT_12(1, elem)-ubar)*helem(nz,elem)*hh_inv
                UVh(2, nz, elem)= UVh(2, nz, elem)+(UVBT_12(2, elem)-vbar)*helem(nz,elem)*hh_inv
                UV( 1, nz, elem)= UVh(1, nz, elem)/helem(nz,elem)  ! velocities are still needed    
                UV( 2, nz, elem)= UVh(2, nz, elem)/helem(nz,elem)  ! to compute momentum advection
            end do
            nz=nzmax
            UVh(1, nz, elem)= UVh(1, nz, elem)
            UVh(2, nz, elem)= UVh(2, nz, elem)
            UV( 1, nz, elem)= UVh(1, nz, elem)/helem(nz,elem)  ! velocities are still needed    
            UV( 2, nz, elem)= UVh(2, nz, elem)/helem(nz,elem)  ! to compute momentum advection
        end do
!$OMP END DO       
!$OMP END PARALLEL   

        !_______________________________________________________________________
!$OMP MASTER        
        call exchange_elem(UVh, partit)    ! 
        call exchange_elem(UV , partit)   ! Check if this is needed 
!$OMP END MASTER
!$OMP BARRIER
    
    !___________________________________________________________________________
    ! skip trimming only adnvance velocity
    ! Velocity will be used to advance momentum
    elseif (mode==3) then 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax, ubar, vbar, hh_inv)
!$OMP DO      
        do elem=1, myDim_elem2D    
            nzmin  = ulevels(elem)
            nzmax  = nlevels(elem)-1
            
            !___________________________________________________________________
            do nz=nzmin, nzmax
                UVh(1, nz, elem)= UVh(1, nz, elem) + UV_rhs(1, nz, elem)
                UVh(2, nz, elem)= UVh(2, nz, elem) + UV_rhs(2, nz, elem)
                UV( 1, nz, elem)= UVh(1, nz, elem)/helem(nz,elem)  ! velocities are still needed    
                UV( 2, nz, elem)= UVh(2, nz, elem)/helem(nz,elem)  ! to compute momentum advection
            end do
        end do
!$OMP END DO       
!$OMP END PARALLEL   

        !_______________________________________________________________________
!$OMP MASTER        
        call exchange_elem(UVh, partit)    ! 
        call exchange_elem(UV , partit)   ! Check if this is needed 
!$OMP END MASTER
!$OMP BARRIER

    end if ! --> if (mode==1) then
end subroutine update_trim_vel_ale_vtransp
!
!
!_______________________________________________________________________________
! Trim U and Uh to be consistent with BT transport
subroutine compute_thickness_zstar(dynamics, partit, mesh)
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE MOD_DYN
    use g_comm
    implicit none
    !___________________________________________________________________________
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout)   , target :: mesh
    integer                               :: node, elem, nz, nzmin, nzmax, elnodes(3) 
    real(kind=WP)                         :: hh_inv
    
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer :: eta_n
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    eta_n => dynamics%eta_n(:)
    
    !___________________________________________________________________________
    ! leave bottom layer again untouched
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax, hh_inv)
!$OMP DO    
    do node=1, myDim_nod2D+eDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D_min(node)-1
        !PS hh_inv=-1.0_WP/zbar(nzmax)
        hh_inv=-1.0_WP/zbar_3d_n(nzmax, node)
        do nz=nzmin, nzmax-1
            hnode_new(nz,node)=(zbar(nz)-zbar(nz+1))*(1.0_WP+hh_inv*eta_n(node))
        end do
    end do ! --> do node=1, myDim_nod2D+eDim_nod2D
!$OMP END DO       
!$OMP END PARALLEL 
    
!PS     !___________________________________________________________________________
!PS     ! --> update mean layer thinkness at element
!PS !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem, elnodes, nz, nzmin, nzmax)
!PS !$OMP DO        
!PS     do elem=1, myDim_elem2D
!PS         nzmin = ulevels(elem)
!PS         nzmax = nlevels(elem)-1
!PS         elnodes=elem2D_nodes(:,elem)
!PS         do nz=nzmin, nzmax-1
!PS             helem(nz,elem)=sum(hnode_new(nz,elnodes))/3.0_WP
!PS         end do
!PS     end do
!PS !$OMP END DO       
!PS !$OMP END PARALLEL 
!PS 
!PS     !___________________________________________________________________________
!PS !$OMP MASTER
!PS     call exchange_elem(helem, partit)
!PS !$OMP END MASTER
!PS !$OMP BARRIER

end subroutine compute_thickness_zstar


