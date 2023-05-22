!
!
!_______________________________________________________________________________
module momentum_adv_scalar_4splitexpl_interface
    interface
        subroutine momentum_adv_scalar_4splitexpl(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module
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
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module
!
!
!_______________________________________________________________________________
module compute_ssh_split_explicit_interface
    interface
        subroutine compute_BC_BT_SE_vtransp(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine

        subroutine compute_BT_step_SE_ale(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
        
        subroutine update_trim_vel_ale_vtransp(mode, dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        integer       , intent(in)            :: mode
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

!
!
!_______________________________________________________________________________
! Transports are used instead of velocities, Urhs, Vrhs are also for transports. 
subroutine momentum_adv_scalar_4splitexpl(dynamics, partit, mesh)
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
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhsAB, UVnode_rhs, UVnode, UVh
    real(kind=WP), dimension(:,:)  , pointer :: Wvel_e
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV        =>dynamics%uv(:,:,:)
    UV_rhsAB  =>dynamics%uv_rhsAB(:,:,:)
    UVnode_rhs=>dynamics%work%uvnode_rhs(:,:,:)
    UVnode    =>dynamics%uvnode(:,:,:)
    Wvel_e    =>dynamics%w_e(:,:)
    UVh       =>dynamics%se_uvh(:,:,:)
    !___________________________________________________________________________
    ! 1st. compute vertical momentum advection component: w * du/dz, w*dv/dz
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, elem, ed, nz, nl1, ul1, nl2, ul2, nl12, ul12, &
!$OMP                                  uv1, uv2, uv12, qc, qu, qd, wu, wv, &
!$OMP                                  ednodes, edelem, dx1, dy1, dx2, dy2, un, uu, vv)

!$OMP DO    
    do node=1, myDim_nod2D
        nl1 = nlevels_nod2D(node)
        ul1 = ulevels_nod2D(node)
        
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
            
            uv1  = UVnode(1, nz  , node)+(2*qc+qu)*hnode(nz  , node)/6.0_WP                    ! Gradient reconstruction 2(2qc+qu)(h/2)(1/6)            
            uv2  = UVnode(1, nz-1, node)-(2*qc+qd)*hnode(nz-1, node)/6.0_WP
            uv12 = (Wvel_e(nz, node)+abs(Wvel_e(nz, node)))*uv1+ &
                   (Wvel_e(nz, node)-abs(Wvel_e(nz, node)))*uv2
            wu(nz)=0.5_WP*(num_ord*(uv1+uv2)*Wvel_e(nz, node)+(1.0_WP-num_ord)*uv12)*area(nz, node)
            
            qc   = (UVnode(2, nz-1, node)-UVnode(2, nz  , node))/(hnode(nz-1, node)+hnode(nz  , node))
            qu   = (UVnode(2, nz  , node)-UVnode(2, nz+1, node))/(hnode(nz  , node)+hnode(nz+1, node))
            qd   = (UVnode(2, nz-2, node)-UVnode(2, nz-1, node))/(hnode(nz-2, node)+hnode(nz-1, node))
            
            uv1  = UVnode(2, nz  , node)+(2*qc+qu)*hnode(nz  , node)/6.0_WP
            uv2  = UVnode(2, nz-1, node)-(2*qc+qd)*hnode(nz-1, node)/6.0_WP
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
        do nz=1, nl1-1
            UVnode_rhs(1, nz, node)= -(wu(nz)-wu(nz+1))
            UVnode_rhs(2, nz, node)= -(wv(nz)-wv(nz+1))
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
        nl1     = nlevels(edelem(1))-1
        ul1     = ulevels(edelem(1))
        
        !_______________________________________________________________________
        !NR --> Natalja Style
        un1          = 0.0_WP
        un1(ul1:nl1) = (  UVh(2, ul1:nl1, edelem(1))*edge_cross_dxdy(1,ed) & 
                        - UVh(1, ul1:nl1, edelem(1))*edge_cross_dxdy(2,ed)) 
        
        !_______________________________________________________________________
        ! if edelem(2)==0 than edge is boundary edge
        if(edelem(2)>0) then
            nl2 = nlevels(edelem(2))-1
            ul2 = ulevels(edelem(2))

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
                do nz=ul1 , nl1-1
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

!PS         !_______________________________________________________________________
!PS         ! nl12 ... minimum number of layers -1 between element edelem(1) & edelem(2) that 
!PS         ! contribute to edge ed
!PS         ! nu12 ... upper index of layers between element edelem(1) & edelem(2) that 
!PS         ! contribute to edge ed
!PS         ! be carefull !!! --> if ed is a boundary edge than edelem(1)~=0 and edelem(2)==0
!PS         !                     that means nl1>0, nl2==0, nl12=min(nl1,nl2)=0 !!!
!PS         ul12 = max(ul1, ul2)
!PS         nl12 = min(nl1, nl2)
!PS         
!PS         !_______________________________________________________________________
!PS         ! (A1) goes only into this loop when the edge has only facing element
!PS         ! edelem(1) --> so the edge is a boundary edge --> this is for ocean 
!PS         ! surface in case of cavity
!PS         do nz=nu1, nu12-1
!PS             un= (UVh(2, nz, edelem(1))*x1 - UVh(1, nz, edelem(1))*y1) 
!PS             uu=un*UV(1, nz, edelem(1)) ! the momentum to be carried depends on velocities
!PS             vv=un*UV(2, nz, edelem(1))
!PS             UVnode_rhs(1, nz, ednodes(1))=UVnode_rhs(1, nz, ednodes(1))+uu
!PS             UVnode_rhs(1, nz, ednodes(2))=UVnode_rhs(1, nz, ednodes(2))-uu
!PS             UVnode_rhs(2, nz, ednodes(1))=UVnode_rhs(2, nz, ednodes(1))+vv
!PS             UVnode_rhs(2, nz, ednodes(2))=UVnode_rhs(2, nz, ednodes(2))-vv
!PS         end do ! --> do nz=nu1, nu12-1
!PS         
!PS         !_______________________________________________________________________
!PS         ! (A2) goes only into this loop when the edge has a facing elemenmt
!PS         ! edelem(2) --> so the edge is a boundary edge --> this is for ocean 
!PS         ! surface in case of cavity
!PS         if (nu2 > 0) then 
!PS             do nz=nu2, nu12-1
!PS                 un=-(UVh(2, nz, edelem(2))*x2- UVh(1, nz, edelem(2))*y2)
!PS                 uu=un*UV(1, nz, edelem(2))
!PS                 vv=un*UV(2, nz, edelem(2))
!PS                 UVnode_rhs(1, nz, ednodes(1))=UVnode_rhs(1, nz, ednodes(1))+uu
!PS                 UVnode_rhs(1, nz, ednodes(2))=UVnode_rhs(1, nz, ednodes(2))-uu
!PS                 UVnode_rhs(2, nz, ednodes(1))=UVnode_rhs(2, nz, ednodes(1))+vv
!PS                 UVnode_rhs(2, nz, ednodes(2))=UVnode_rhs(2, nz, ednodes(2))-vv
!PS             end do ! --> do nz=nu2, nu12-1
!PS         end if 
!PS         
!PS         !_______________________________________________________________________
!PS         ! (B) Both segments
!PS         ! loop over depth layers from shared upper layer index nu12 to shared 
!PS         ! lower layer index nl12
!PS         do nz=nu12, nl12
!PS             un=  (UVh(2, nz, edelem(1))*x1 - UVh(1, nz, edelem(1))*y1) &
!PS                 -(UVh(2, nz, edelem(2))*x2 - UVh(1, nz, edelem(2))*y2)
!PS             uu=un*(UV(1, nz, edelem(1)) + UV(1, nz, edelem(2)))*0.5_WP! the momentum to be carried depends on velocities
!PS             vv=un*(UV(2, nz, edelem(1)) + UV(2, nz, edelem(2)))*0.5_WP
!PS             UVnode_rhs(1, nz, ednodes(1))=UVnode_rhs(1, nz, ednodes(1))+uu
!PS             UVnode_rhs(1, nz, ednodes(2))=UVnode_rhs(1, nz, ednodes(2))-uu
!PS             UVnode_rhs(2, nz, ednodes(1))=UVnode_rhs(2, nz, ednodes(1))+vv
!PS             UVnode_rhs(2, nz, ednodes(2))=UVnode_rhs(2, nz, ednodes(2))-vv
!PS         end do ! --> do nz=nu12, nl12
!PS         
!PS         !_______________________________________________________________________
!PS         ! (C1) remaining segments from the shared lower lyer index nl12 to bottom 
!PS         ! of element edelem(1)
!PS         do nz=nl12+1, nl1
!PS             un= (UVh(2, nz, edelem(1))*x1 - UVh(1, nz, edelem(1))*y1) 
!PS             uu=un*UV(1, nz, edelem(1)) ! the momentum to be carried depends on velocities
!PS             vv=un*UV(2, nz, edelem(1))
!PS             UVnode_rhs(1, nz, ednodes(1))=UVnode_rhs(1, nz, ednodes(1))+uu
!PS             UVnode_rhs(1, nz, ednodes(2))=UVnode_rhs(1, nz, ednodes(2))-uu
!PS             UVnode_rhs(2, nz, ednodes(1))=UVnode_rhs(2, nz, ednodes(1))+vv
!PS             UVnode_rhs(2, nz, ednodes(2))=UVnode_rhs(2, nz, ednodes(2))-vv
!PS         end do ! --> do nz=nl12+1, nl1
!PS         
!PS         !_______________________________________________________________________
!PS         ! (C2) remaining segments from the shared lower lyer index nl12 to bottom 
!PS         ! of element edelem(1)
!PS         do nz=nl12+1, nl2
!PS             un=-(UVh(2, nz, edelem(2))*x2- UVh(1, nz, edelem(2))*y2)
!PS             uu=un*UV(1, nz, edelem(2))
!PS             vv=un*UV(2, nz, edelem(2))
!PS             UVnode_rhs(1, nz, ednodes(1))=UVnode_rhs(1, nz, ednodes(1))+uu
!PS             UVnode_rhs(1, nz, ednodes(2))=UVnode_rhs(1, nz, ednodes(2))-uu
!PS             UVnode_rhs(2, nz, ednodes(1))=UVnode_rhs(2, nz, ednodes(1))+vv
!PS             UVnode_rhs(2, nz, ednodes(2))=UVnode_rhs(2, nz, ednodes(2))-vv
!PS         end do ! --> do nz=nl12+1, nl2
!PS     end do ! --> do ed=1, myDim_edge2D  
!PS !$OMP END DO


    !___________________________________________________________________________
    ! divide total nodal momentum advection by scalar area
!$OMP DO    
    do node=1,myDim_nod2d
        nl1 = nlevels_nod2D(node)-1
        ul1 = ulevels_nod2D(node)
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
        nl1 = nlevels(elem)-1
        ul1 = ulevels(elem)
        UV_rhsAB(1:2, ul1:nl1, elem) = UV_rhsAB(1:2, ul1:nl1, elem) + elem_area(elem)* &
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
          nl1 = nlevels(elem)-1
          ul1 = ulevels(elem)
          dynamics%ke_adv_AB(1:2, ul1:nl1, elem) = dynamics%ke_adv_AB(1:2, ul1:nl1, elem) + elem_area(elem)* &
                ( UVnode_rhs(1:2, ul1:nl1, elem2D_nodes(1, elem)) &
                + UVnode_rhs(1:2, ul1:nl1, elem2D_nodes(2, elem)) & 
                + UVnode_rhs(1:2, ul1:nl1, elem2D_nodes(3, elem))) / 3.0_WP     
       end do
!$OMP END DO
    end if
!$OMP END PARALLEL
end subroutine

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
!
! --> 1st. solve homogenouse part:
!     f(du) = du - dt/h^(*)*( Av*d/dz*du ) = 0
!
! --> 2nd. compute difference quotient at du_i using Gauss-Theorem --> flux form
!     V_i * du_i = dt * [ Av_i/h_i * (du_i-1 - du_i)/(Z_i-1 - Z_i) * A_i * nvec_up(+1)  
!                        +Av_i+1/h_i * (du_i - du_i+1)/(Z_i - Z_i+1) * A_i+1 * nvec_dwn(-1) ] 
!
!     f(du_i) = du_i - dt*Av_i/h_i   * (du_i-1 - du_i)/(Z_i-1 - Z_i) * A_i/V_i+1
!                    + dt*Av_i+1/h_i * (du_i - du_i+1)/(Z_i - Z_i+1) * A_i+1/V_i+1 
!             = 0
!
! --> 3rd. solve for coefficents a, b, c (homogenous part):
!     f(du_i) = [ a*du_i-1 + b*du_i + c*du_i+1 ]
!        |
!        +-> estimate a, b, c by derivation of d(du_i)
!        |
!        +-> a = d[f(du_i)]/d[du_i-1] = - dt*Av_i/h_i / (Z_i-1 - Z_i) 
!        |
!        +-> c = d[f(du_i)]/d[du_i+1] = - dt*Av_i+1/h_i / (Z_i - Z_i+1)
!        |
!        +-> b = d[f(du_i)]/d[du_i]   = 1 + dt*Av_i/h_i / (Z_i-1 - Z_i) + dt*Av_i+1/h_i / (Z_i - Z_i+1)
!                                       1 - a - c
!
! --> 4th. solve inhomogenous part:
!     [ a*du_i-1 + b*du_i + c*du_i+1 ] = RHS/V_i
!
!     RHS/V_i = dt* [ Av_i  /h_i * (u^(n+0.5,*)_i-1 - u^(n+0.5,*)_i)/(Z_i-1 - Z_i) * A_i   * nvec_up(+1) 
!             |      +Av_i+1/h_i * (u^(n+0.5,*)_i - u^(n+0.5,*)_i+1)/(Z_i - Z_i+1) * A_i+1 * nvec_dwn(-1) ]
!             |  
!             V
!             = -a*u^(n+0.5,*)_i-1 + (a+c)*u^(n+0.5,*)_i - c*u^(n+0.5,*)_i+1
!
! --> 5th. solve for du_i --> forward sweep algorithm --> see lower
!     | b_1 c_1 ...            |   |du_1|
!     | a_2 b_2 c_2 ...        |   |du_2|
!     |     a_3 b_3 c_3 ...    | * |du_3| = RHS/V_i
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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a, b, c, ur, vr, cp, up, vp, elem, & 
!$OMP                                  nz, nzmin, nzmax, elnodes, &
!$OMP                                  zinv, m, friction, wu, wv, uu, vv)
    
    !___________________________________________________________________________
!$OMP DO
    do elem=1,myDim_elem2D
        elnodes= elem2D_nodes(:,elem)
        nzmin  = ulevels(elem)
        nzmax  = nlevels(elem)
        
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
subroutine compute_BC_BT_SE_vtransp(dynamics, partit, mesh)
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE MOD_DYN
    USE g_config, only: dt, r_restart
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=WP)       :: vert_sum_u, vert_sum_v, Fx, Fy, ab1, ab2, hh
    integer             :: elem, nz, nzmin, nzmax, elnodes(4)
    logical, save       :: sfirst
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV_rhs
    real(kind=WP), dimension(:)    , pointer :: eta_n
    real(kind=WP), dimension(:,:)  , pointer :: UVhBT_4AB, UVhBC_rhs
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    eta_n     =>dynamics%eta_n(:)
    UV_rhs    =>dynamics%uv_rhs(:,:,:)
    UVhBC_rhs =>dynamics%se_uvhBC_rhs(:,:)
    UVhBT_4AB =>dynamics%se_uvhBT_4AB(:,:)
    
    !___________________________________________________________________________
    ab1=(1.5_WP+epsilon)
    if(sfirst .and. (.not. r_restart)) then
        ab1=1.0_WP
        sfirst=.false.
    end if
    ab2=1.0_WP-ab1

    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax, Fx, Fy, & 
!$OMP                                  vert_sum_u, vert_sum_v, hh)

    !___________________________________________________________________________
!$OMP DO
    do elem=1,  myDim_elem2D
        elnodes= elem2D_nodes(:,elem)
        nzmin  = ulevels(elem)
        nzmax  = nlevels(elem)-1
        
        !_______________________________________________________________________
        Fx=g*dt*sum(gradient_sca(1:4,elem)*eta_n(elnodes))
        Fy=g*dt*sum(gradient_sca(5:8,elem)*eta_n(elnodes))
        
        !_______________________________________________________________________
        ! Sum 3D rhs   
        vert_sum_u=0.0_WP
        vert_sum_v=0.0_WP
        do nz=nzmin, nzmax
            vert_sum_u=vert_sum_u + UV_rhs(1, nz, elem)   
            vert_sum_v=vert_sum_v + UV_rhs(2, nz, elem)   
        end do
        
        !_______________________________________________________________________
        ! Remove the contribution from the elevation
        hh=sum(helem(nzmin:nzmax, elem))
        vert_sum_u=vert_sum_u + Fx*hh
        vert_sum_v=vert_sum_v + Fy*hh
        
        !_______________________________________________________________________
        ! Remove the contribution from the Coriolis
        UVhBC_rhs(1, elem)=vert_sum_u - dt*mesh%coriolis(elem)* &
                                       (ab1*UVhBT_4AB(2,elem)+ab2*UVhBT_4AB(4,elem))         ! for AB-interpolated
        UVhBC_rhs(2, elem)=vert_sum_v + dt*mesh%coriolis(elem)* &
                                       (ab1*UVhBT_4AB(1,elem)+ab2*UVhBT_4AB(3,elem))    
        
        !_______________________________________________________________________
        ! save actual of vertical integrated transport velocity UVhBT_4AB(1:2,elem)  
        ! timestep for next time loop step
        UVhBT_4AB(3:4,elem)=UVhBT_4AB(1:2,elem)    
    end do
!PS !$OMP END DO
!PS !$OMP END PARALLEL
end subroutine compute_BC_BT_SE_vtransp
!
!
!_______________________________________________________________________________
! Barotropic time stepping with Forward-Backward dissipative method
! (Demange et al. 2019)
subroutine compute_BT_step_SE_ale(dynamics, partit, mesh)
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE MOD_DYN
    use g_comm_auto
    USE g_config, only: dt
    IMPLICIT NONE
    !___________________________________________________________________________
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=WP)                   :: dtBT, BT_inv, hh, f, rx, ry, a, d, c1, c2, ax, ay
    real(kind=WP)                   :: deltaX1, deltaY1, deltaX2, deltaY2, thetaBT
    integer                         :: step, elem, elnodes(3), edge, ednodes(2), edelem(2)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:)    , pointer :: eta_n
    real(kind=WP), dimension(:,:)  , pointer :: UVhBC_rhs, UVBT, UVBT_theta, UVBT_mean, UVBT_12
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    eta_n     =>dynamics%eta_n(:)
    UVhBC_rhs =>dynamics%se_uvhBC_rhs(:,:)
    UVBT      =>dynamics%se_uvBT(:,:)
    UVBT_theta=>dynamics%se_uvBT_theta(:,:)
    UVBT_mean =>dynamics%se_uvBT_mean(:,:)
    UVBT_12   =>dynamics%se_uvBT_12(:,:)
    
    !___________________________________________________________________________
    ! Dissipation parameter of FB dissipative method 0.14 is the default value 
    ! from Demange et al.
    thetaBT= 0.14_WP
    
    ! BTsteps should be 30 or 40.
    dtBT   = dt/dynamics%splitexpl_BTsteps                
    BT_inv = 1.0_WP/(1.0_WP*dynamics%splitexpl_BTsteps)
    
    !___________________________________________________________________________
    ! eta_n   elevation used in BT stepping, it is just a copy of eta_n
    ! UBT and VBT are transport velocities
    UVBT_mean= 0.0_WP
    do step=1, dynamics%splitexpl_BTsteps
        !
        ! Dissipative forward--backward time stepping
        !
        !_______________________________________________________________________
        ! Advance velocities. I use SI stepping for the Coriolis
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            !PS hh = -zbar(nlevels(elem))+sum(eta_n(elnodes))/3.0_WP 
            hh = -zbar_e_bot(elem)+sum(eta_n(elnodes))/3.0_WP ! Total fluid depth
            f  = mesh%coriolis(elem)
            rx = dtBT*(-g*hh*sum(gradient_sca(1:3,elem)*eta_n(elnodes)) + f*UVBT(1, elem)) + BT_inv*UVhBC_rhs(1, elem)
            ry = dtBT*(-g*hh*sum(gradient_sca(4:6,elem)*eta_n(elnodes)) - f*UVBT(2, elem)) + BT_inv*UVhBC_rhs(2, elem)
            
            ! Semi-Implicit Coriolis
            a  = dtBT*f*0.5_WP
            d  = 1.0_WP/(1.0_WP+a*a)
            ax = d*(rx+a*ry)
            ay = d*(-a*rx+ry)
            UVBT(      1, elem) = UVBT(      1, elem) + ax   ! New velocity
            UVBT(      2, elem) = UVBT(      2, elem) + ay
            
            ! velocities for dissipative time stepping of thickness equation
            UVBT_theta(1, elem) = UVBT_theta(1, elem) + thetaBT*ax   ! New vel*(1+thetaBT)-old vel *thetaBT
            UVBT_theta(2, elem) = UVBT_theta(2, elem) + thetaBT*ay
            
            ! Mean BT velocity to trim 3D velocity in tracers
            UVBT_mean( 1, elem) = UVBT_mean( 1, elem) + UVBT_theta(1, elem)*BT_inv
            UVBT_mean( 2, elem) = UVBT_mean( 2, elem) + UVBT_theta(2, elem)*BT_inv
        end do
        
        !_______________________________________________________________________
        ! Store mid-step velocity (to trim 3D velocities in momentum)
        if(step==dynamics%splitexpl_BTsteps/2) then
            UVBT_12=UVBT
        end if 
        
        !_______________________________________________________________________
        ! Advance thickness
        do edge=1, myDim_edge2D
            ednodes = edges(:,edge)
            edelem  = edge_tri(:,edge)
            deltaX1 = edge_cross_dxdy(1,edge)
            deltaY1 = edge_cross_dxdy(2,edge)
            
            c1 = UVBT_theta(2, edelem(1))*deltaX1 - UVBT_theta(1, edelem(1))*deltaY1
            c2 = 0.0_WP
            if(edelem(2)>0) then
                deltaX2=edge_cross_dxdy(3,edge)
                deltaY2=edge_cross_dxdy(4,edge)
                c2=-(UVBT_theta(2, edelem(2))*deltaX2 - UVBT_theta(1, edelem(2))*deltaY2)
            end if
            eta_n(ednodes(1))=eta_n(ednodes(1))+(c1+c2)*dtBT/area(1,ednodes(1))
            eta_n(ednodes(2))=eta_n(ednodes(2))-(c1+c2)*dtBT/area(1,ednodes(2))
        end do
        
        !_______________________________________________________________________
        call exchange_nod(eta_n, partit)
        
    end do ! --> do step=1, dynamics%splitexpl_BTsteps
    
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
    use g_comm_auto
    IMPLICIT NONE
    !___________________________________________________________________________
    integer       , intent(in)            :: mode
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    integer                               :: i,elem, elnodes(4), nz, m
    real(kind=WP)                         :: eta(4), ff
    real(kind=WP)                         :: Fx, Fy, ubar, vbar, hh_inv
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UVh, UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: UVBT_mean, UVBT_12
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UVh       =>dynamics%se_uvh(:,:,:)
    UV_rhs    =>dynamics%uv_rhs(:,:,:)
    UVBT_mean =>dynamics%se_uvBT_mean(:,:)
    UVBT_12   =>dynamics%se_uvBT_12(:,:)
    if (mode==2) then 
        UV    =>dynamics%uv(:,:,:)
    end if
    
    !___________________________________________________________________________
    !
    ! Transport version
    !
    !___________________________________________________________________________
    ! Trim full velocity to ensure that its vertical sum equals to the mean BT velocity
    ! The trimmed Uh,Vh are consistent with new total height defined by eta_n   
    if (mode==1) then              
        do elem=1, myDim_elem2D
            ubar=0.0_WP
            vbar=0.0_WP
            hh_inv=0.0_WP
            do nz=1, nlevels(elem)-1
                ! Update transport velocities
                UVh(1, nz, elem)=UVh(1, nz, elem)+UV_rhs(1, nz, elem)   
                UVh(2, nz, elem)=UVh(2, nz, elem)+UV_rhs(2, nz, elem)
                ubar=ubar+UVh(1, nz, elem)
                vbar=vbar+UVh(2, nz, elem)
                hh_inv=hh_inv+helem(nz,elem)
            end do
            
            hh_inv=1.0_WP/hh_inv     ! Postpone AB and 2nd order, just use available thickness
            DO nz=1, nlevels(elem)-1
                UVh(1, nz, elem)= UVh(1, nz, elem)+(UVBT_mean(1, elem)-ubar)*helem(nz,elem)*hh_inv
                UVh(2, nz, elem)= UVh(2, nz, elem)+(UVBT_mean(2, elem)-vbar)*helem(nz,elem)*hh_inv
                UV( 1, nz, elem)= UVh(1, nz, elem)/helem(nz,elem)  ! velocities are still needed    
                UV( 2, nz, elem)= UVh(2, nz, elem)/helem(nz,elem) 
            end do
        end do
        call exchange_elem(UVh, partit)  ! This exchange can be avoided, but test first.
        
    !___________________________________________________________________________
    ! Trim to make velocity consistent with BT velocity at n+1/2
    ! Velocity will be used to advance momentum
    else                         
        do elem=1, myDim_elem2D    
            ubar=0.0_WP
            vbar=0.0_WP
            hh_inv=0.0_WP    
            do nz=1, nlevels(elem)-1
                ubar=ubar+UVh(1, nz, elem)
                vbar=vbar+UVh(2, nz, elem)
                hh_inv=hh_inv+helem(nz,elem)
            end do
            
            hh_inv=1.0_WP/hh_inv    ! Postpone 2nd order, just use available thickness 
            do nz=1, nlevels(elem)-1
                UVh(1, nz, elem)= UVh(1, nz, elem)+(UVBT_12(1, elem)-ubar)*helem(nz,elem)*hh_inv
                UVh(2, nz, elem)= UVh(2, nz, elem)+(UVBT_12(2, elem)-vbar)*helem(nz,elem)*hh_inv
                UV( 1, nz, elem)= UVh(1, nz, elem)/helem(nz,elem)  ! velocities are still needed    
                UV( 2, nz, elem)= UVh(2, nz, elem)/helem(nz,elem)  ! to compute momentum advection
            end do
        end do
        call exchange_elem(UVh, partit)    ! 
        call exchange_elem(UV , partit)   ! Check if this is needed 
    end if ! --> if (mode==1) then
end subroutine update_trim_vel_ale_vtransp



