module momentum_adv_scalar_4splitexpl_interface
    interface
        subroutine momentum_adv_scalar_4splitexpl(dynamics, partit, mesh)
        use mod_mesh
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
!   Transports are used instead of velocities, Urhs, Vrhs are also for transports. 
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









