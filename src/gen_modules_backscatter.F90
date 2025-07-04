module g_backscatter

    !___________________________________________________________________________
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN

    !___________________________________________________________________________
    USE o_ARRAYS, only: bvfreq

    !___________________________________________________________________________
    USE o_param
    USE g_CONFIG
    USE g_comm_auto
    USE g_support
    USE g_rotate_grid
    IMPLICIT NONE
    
    !___________________________________________________________________________
    ! allocate backscatter arrays
    real(kind=WP), allocatable, dimension(:,:)  :: v_back
    real(kind=WP), allocatable, dimension(:,:)  :: uke, uke_back, uke_dis, uke_dif, uke_adv
    real(kind=WP), allocatable, dimension(:,:)  :: uke_rhs, uke_rhs_old
    real(kind=WP), allocatable, dimension(:,:)  :: UV_dis_posdef_b2, UV_dis_posdef, UV_back_posdef                                                
    real(kind=WP), allocatable, dimension(:,:,:):: UV_back, UV_dis
    real(kind=WP), allocatable, dimension(:,:,:):: UV_dis_tend, UV_total_tend, UV_back_tend

    contains 
    !
    !
    !___________________________________________________________________________
    ! allocate/initialise backscatter arrays
    subroutine init_backscatter(dynamics, partit, mesh)
        implicit none 
        integer                               :: elem_size
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
        
        elem_size = myDim_elem2D + eDim_elem2D
        allocate(v_back(          nl-1, elem_size)) ! Backscatter viscosity
        allocate(uke(             nl-1, elem_size)) ! Unresolved kinetic energy for backscatter coefficient
        allocate(uke_dis(         nl-1, elem_size)) 
        allocate(uke_back(        nl-1, elem_size)) 
        allocate(uke_dif(         nl-1, elem_size))
        allocate(uke_rhs(         nl-1, elem_size))
        allocate(uke_rhs_old(     nl-1, elem_size))
        allocate(UV_dis(       2, nl-1, elem_size))
        allocate(UV_back(      2, nl-1, elem_size))
        allocate(UV_dis_tend(  2, nl-1, elem_size))
        allocate(UV_back_tend( 2, nl-1, elem_size))
        allocate(UV_total_tend(2, nl-1, elem_size))
        uke           = 0.0_WP
        v_back        = 0.0_WP
        uke_dis       = 0.0_WP
        uke_dif       = 0.0_WP
        uke_back      = 0.0_WP
        uke_rhs       = 0.0_WP
        uke_rhs_old   = 0.0_WP
        UV_dis        = 0.0_WP     
        UV_dis_tend   = 0.0_WP
        UV_back       = 0.0_WP
        UV_back_tend  = 0.0_WP
        UV_total_tend = 0.0_WP
        
        if (dynamics%uke_advection) then
            allocate(uke_adv(    nl-1, elem_size))
            uke_adv     = 0.0_WP
        endif

    end subroutine init_backscatter

    !
    !
    !_______________________________________________________________________________
    subroutine visc_filt_dbcksc(dynamics, partit, mesh)
        IMPLICIT NONE
        
        real(kind=WP)  :: u1, v1, le(2), len, len2, crosslen, vi, uke1 
        integer       :: nz, ed, el(2)
        real(kind=WP)  , allocatable  :: uke_d(:,:)
        !!PS real(kind=WP)  , allocatable  :: UV_back(:,:,:), UV_dis(:,:,:)
        real(kind=WP)  , allocatable  :: uuu(:)
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        real(kind=WP) , dimension(:,:,:), pointer :: UV, UV_rhs
        real(kind=WP) , dimension(:,:)  , pointer :: U_c, V_c
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
        
        UV     => dynamics%uv(:,:,:)
        UV_rhs => dynamics%uv_rhs(:,:,:)
        U_c    => dynamics%work%u_c(:,:)
        V_c    => dynamics%work%v_c(:,:)
        
        ! An analog of harmonic viscosity operator.  
        ! It adds to the rhs(0) Visc*(u1+u2+u3-3*u0)/area
        ! on triangles, which is Visc*Laplacian/4 on equilateral triangles. 
        ! The contribution from boundary edges is neglected (free slip). 
        ! Filter is applied twice. 
        ed=myDim_elem2D+eDim_elem2D
        !!PS allocate(UV_back(2,nl-1,ed), UV_dis(2,nl-1, ed)) 
        allocate(uke_d(nl-1,ed)) 
        allocate(uuu(ed))
        UV_back= 0.0_WP
        UV_dis = 0.0_WP
        uke_d  = 0.0_WP
        U_c    = 0.0_WP
        V_c    = 0.0_WP
        
        DO ed=1, myDim_edge2D+eDim_edge2D
            if(myList_edge2D(ed)>edge2D_in) cycle
            el=edge_tri(:,ed)
            !New viscosity lines
            len=sqrt(sum(elem_area(el)))
            DO  nz=1,minval(nlevels(el))-1
                u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
                v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
                ! New viscosity lines
                vi=u1*u1+v1*v1
                vi=sqrt(max(dynamics%visc_gamma0, max(dynamics%visc_gamma1*sqrt(vi), dynamics%visc_gamma2*vi))*len)
                u1=u1*vi
                v1=v1*vi
                
                U_c(nz,el(1))=U_c(nz,el(1))-u1
                U_c(nz,el(2))=U_c(nz,el(2))+u1
                V_c(nz,el(1))=V_c(nz,el(1))-v1
                V_c(nz,el(2))=V_c(nz,el(2))+v1
            END DO 
        END DO
        
        call exchange_elem(U_c, partit)
        call exchange_elem(V_c, partit) 
        
        DO ed=1, myDim_edge2D+eDim_edge2D
            if(myList_edge2D(ed)>edge2D_in) cycle
            el=edge_tri(:,ed)
            le=edge_dxdy(:,ed)
            le(1)=le(1)*sum(elem_cos(el))*0.25_WP
            len=sqrt(le(1)**2+le(2)**2)*r_earth
            !Weighting new viscosity 
            len2=sqrt(sum(elem_area(el)))
            le(1)=edge_cross_dxdy(1,ed)-edge_cross_dxdy(3,ed)
            le(2)=edge_cross_dxdy(2,ed)-edge_cross_dxdy(4,ed)
            crosslen=sqrt(le(1)**2+le(2)**2) 
            
            DO  nz=1,minval(nlevels(el))-1
                vi=dt*len*(v_back(nz,el(1))+v_back(nz,el(2)))/crosslen
                !if(mype==0) write(*,*) 'vi ', vi , ' and ed' , ed
                !if(mype==0) write(*,*) 'dt*len/crosslen ', dt*len/crosslen, ' and ed' , ed
                !vi=max(vi,0.005*len*dt) ! This helps to reduce noise in places where 
                                        ! Visc is small and decoupling might happen 
                !Backscatter contribution
                u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))*vi
                v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))*vi
                
                !UKE diffusion
                vi=dt*len*(dynamics%K_back*sqrt(elem_area(el(1))/scale_area)+dynamics%K_back*sqrt(elem_area(el(2))/scale_area))/crosslen
                uke1=(uke(nz,el(1))-uke(nz,el(2)))*vi
                
                UV_back(1,nz,el(1))=UV_back(1,nz,el(1))-u1/elem_area(el(1))
                UV_back(1,nz,el(2))=UV_back(1,nz,el(2))+u1/elem_area(el(2))
                UV_back(2,nz,el(1))=UV_back(2,nz,el(1))-v1/elem_area(el(1))
                UV_back(2,nz,el(2))=UV_back(2,nz,el(2))+v1/elem_area(el(2))  
                
                !Correct scaling for the diffusion?
                uke_d(nz,el(1))=uke_d(nz,el(1))-uke1/elem_area(el(1))
                uke_d(nz,el(2))=uke_d(nz,el(2))+uke1/elem_area(el(2))
                
                !New viscosity param part
                u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
                v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
                vi=u1*u1+v1*v1
                vi=-dt*sqrt(max(dynamics%visc_gamma0, max(dynamics%visc_gamma1*sqrt(vi), dynamics%visc_gamma2*vi))*len2)
                u1=vi*(U_c(nz,el(1))-U_c(nz,el(2)))
                v1=vi*(V_c(nz,el(1))-V_c(nz,el(2)))
                
                UV_dis(1,nz,el(1))=UV_dis(1,nz,el(1))-u1/elem_area(el(1))
                UV_dis(1,nz,el(2))=UV_dis(1,nz,el(2))+u1/elem_area(el(2))
                UV_dis(2,nz,el(1))=UV_dis(2,nz,el(1))-v1/elem_area(el(1))
                UV_dis(2,nz,el(2))=UV_dis(2,nz,el(2))+v1/elem_area(el(2))
                
            END DO 
        END DO
        call exchange_elem(UV_back, partit)
        
        DO  nz=1, nl-1
            uuu=0.0_WP
            uuu=UV_back(1,nz,:)
            call smooth_elem(uuu,dynamics%smooth_back_tend, partit, mesh)
            UV_back(1,nz,:)=uuu
            uuu=0.0_WP
            uuu=UV_back(2,nz,:)
            call smooth_elem(uuu,dynamics%smooth_back_tend, partit, mesh)
            UV_back(2,nz,:)=uuu 
        END DO
        
        DO ed=1, myDim_elem2D
            DO  nz=1,nlevels(ed)-1
            UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+UV_dis(1,nz,ed)+UV_back(1,nz,ed)
            UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+UV_dis(2,nz,ed)+UV_back(2,nz,ed)               
            END DO 
        END DO
        
        UV_dis_tend=UV_dis!+UV_back
        UV_total_tend=UV_dis+UV_back
        UV_back_tend=UV_back
        uke_dif=uke_d
        
        call uke_update(dynamics, partit, mesh)
        
        !!PS deallocate(UV_dis,UV_back) 
        deallocate(uke_d)
        deallocate(uuu)
    end subroutine visc_filt_dbcksc


    subroutine backscatter_adv(dynamics, partit, mesh)
        IMPLICIT NONE
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(inout), target :: partit
        integer                               :: n, nz, el1, el2, ul1
        integer                               :: nl1, nl2, nod(2), el, ed, k, nle
        real(kind=WP), dimension(:,:,:), pointer :: UV
        real(kind=WP), dimension(:,:),   pointer :: Wvel_e
        real(kind=WP), allocatable, dimension(:) :: contr1, contr2, wuke
        real(kind=WP), allocatable, dimension(:,:)  :: advnode
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
        Wvel_e    =>dynamics%w_e(:,:)
        UV        =>dynamics%UV(:,:,:)
        
        allocate(contr1(1:nl-1))
        allocate(contr2(1:nl-1))
        allocate(wuke(1:nl-1))
        allocate(advnode(nl-1, myDim_nod2d+eDim_nod2D))

        advnode=0.0_WP
        contr1=0.0_WP
        contr2=0.0_WP
        wuke=0.0_WP
         
        !vertical advection
        do n=1,myDim_nod2d
           nl1 = nlevels_nod2D(n)-1
           
           !find average of uke between vertical layers
           do k=1,nod_in_elem2D_num(n)
               el = nod_in_elem2D(k,n)
               nle = nlevels(el)-1
           
               wuke(1) = wuke(1) + uke(1,el)*elem_area(el)
               wuke(2:nle) = wuke(2:nle) + 0.5_WP*(uke(2:nle,el)+uke(1:nle-1,el))*elem_area(el)      

           end do

           wuke(1:nl) = wuke(1:nl1)*Wvel_e(1:nl1,n)

           do nz=1,nl1
              ! Here 1/3 because 1/3 of the area is related to the node
              advnode(nz,n) = - (wuke(nz) - wuke(nz+1) ) / (3._WP*hnode(nz,n)) 
           end do
           
           advnode(nl1+1:nl-1,n) = 0._WP

        end do

        !horizontal advection
        DO ed=1, myDim_edge2D
           nod = edges(:,ed) 
           el1  = edge_tri(1,ed)   
           el2  = edge_tri(2,ed)
           nl1 = nlevels(el1)-1

           !calculate fluxes throught the each edge 
           !for the edges that are on the boundary and have only one element 
           contr1(1:nl1) = UV(2,1:nl1,el1)*edge_cross_dxdy(1,ed) - UV(1,1:nl1,el1)*edge_cross_dxdy(2,ed)  
           
           ! if we have the other element as well
           if (el2>0) then
              nl2 = nlevels(el2)-1
              
              contr2(1:nl2) = - UV(2,1:nl2,el2)*edge_cross_dxdy(3,ed) + UV(1,1:nl2,el2)*edge_cross_dxdy(4,ed)
              
              contr1(nl1+1:max(nl1,nl2)) = 0._WP
              contr2(nl2+1:max(nl1,nl2)) = 0._WP
              
              if (nod(1) <= myDim_nod2d) then
                 do nz=1, max(nl1,nl2)
                    advnode(nz,nod(1)) = advnode(nz,nod(1)) + contr1(nz)*uke(nz,el1) + contr2(nz)*uke(nz,el2)             
                 end do
              endif
              
              if (nod(2) <= myDim_nod2d) then
                 do nz=1, max(nl1,nl2)
                    advnode(nz,nod(2)) = advnode(nz,nod(2)) - contr1(nz)*uke(nz,el1) - contr2(nz)*uke(nz,el2)
                 end do
              endif
           
           else  ! ed is a boundary edge, there is only the contribution from el1
              if (nod(1) <= myDim_nod2d) then
                 do nz=1, nl1
                    advnode(nz,nod(1)) = advnode(nz,nod(1)) + contr1(nz)*uke(nz,el1)
                 end do
              endif
              ! second edge node
              if  (nod(2) <= myDim_nod2d) then
                 do nz=1, nl1
                    advnode(nz,nod(2)) = advnode(nz,nod(2)) - contr1(nz)*uke(nz,el1)
                 end do
              endif
           endif
        end do

        ! Multiply by the segment area
        do n=1,myDim_nod2d
           nl1 = nlevels_nod2D(n)-1
           advnode(1:nl1,n) = advnode(1:nl1,n) * area_inv(1:nl1,n)
        end do

        call exchange_nod(advnode, partit)

        !Calculate advection as an average across the nodes
        do n=1,myDim_nod2d
           nl1 = nlevels_nod2D(n)-1
           advnode(1:nl1,n) = advnode(1:nl1,n) *areasvol_inv(1:nl1,n)
        end do

        do el=1, myDim_elem2D
           nl1 = nlevels(el)-1
           uke_adv(1:nl1,el) = uke_adv(1:nl1,el) &
                + elem_area(el)*(advnode(1:nl1,elem2D_nodes(1,el)) &
                + advnode(1:nl1,elem2D_nodes(2,el)) & 
                + advnode(1:nl1,elem2D_nodes(3,el)))/ 3.0_WP

            do nz=1,nlevels(el)-1   
                uke_adv(nz,el)=dt*uke_adv(nz,el)/elem_area(el)
            end do
        end do

        deallocate(contr1)
        deallocate(contr2)
        deallocate(wuke)
        deallocate(advnode)
    end subroutine backscatter_adv

    !
    !
    !_______________________________________________________________________________
    subroutine backscatter_coef(dynamics, partit, mesh)
        IMPLICIT NONE
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(inout), target :: partit
        integer                               :: elem, nz
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

        !Potentially add the Rossby number scaling to the script...
        !check if sign is right! Different in the Jansen paper
        !Also check with the normalization by area; as before we use element length sqrt(2*elem_area(ed))
        v_back=0.0_WP
        DO  elem=1, myDim_elem2D 
            DO  nz=1,nlevels(elem)-1
                !v_back(1,ed)=c_back*sqrt(2.0_WP*elem_area(ed))*sqrt(max(2.0_WP*uke(1,ed),0.0_WP))*(3600.0_WP*24.0_WP/tau_c)*4.0_WP/sqrt(2.0_WP*elem_area(ed))**2 !*sqrt(max(2.0_WP*uke(1,ed),0.0_WP))
                !v_back(nz,elem)=-c_back*sqrt(4._WP/sqrt(3.0_WP)*elem_area(elem))*sqrt(max(2.0_WP*uke(nz,elem),0.0_WP)) !Is the scaling correct
                v_back(nz,elem)=min(-dynamics%c_back*sqrt(elem_area(elem))*sqrt(max(2.0_WP*uke(nz,elem),0.0_WP)),0.2*elem_area(elem)/dt) !Is the scaling correct
                !Scaling by sqrt(2*elem_area) or sqrt(elem_area)?
            END DO
        END DO
        call exchange_elem(v_back, partit)
    end subroutine backscatter_coef
    !
    !
    !_______________________________________________________________________________
    subroutine uke_update(dynamics, partit, mesh)
        IMPLICIT NONE

        !I had to change uke(:) to uke(:,:) to make output and restart work!!
        !Why is it necessary to implement the length of the array? It doesn't work without!
        !integer, intent(in)        :: t_levels
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh

        real(kind=WP)     :: hall, h1_eta, hnz, vol
        integer           :: elnodes(3), nz, ed, edi, node, j, elem, q
        real(kind=WP), allocatable  :: uuu(:), work_array(:), U_work(:,:), V_work(:,:), rosb_array(:,:), work_uv(:)
        integer           :: kk, nzmax, el
        real(kind=WP)     :: c1, rosb, vel_u, vel_v, vel_uv, scaling, reso
        real*8            :: c_min=0.5, f_min=1.e-6, r_max=200000., ex, ey, a1, a2, len_reg, dist_reg(2) ! Are those values still correct?
        real(kind=WP), dimension(:,:,:), pointer :: UV
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
        UV => dynamics%uv(:,:,:)
        
        !rosb_dis=1._WP !Should be variable to control how much of the dissipated energy is backscattered
        !rossby_num=2
        
        ed=myDim_elem2D+eDim_elem2D
        allocate(uuu(ed)) 
        
        uke_back=0.0_WP
        uke_dis=0.0_WP
        DO ed=1, myDim_elem2D
            DO nz=1, nlevels(ed)-1  
                uke_dis(nz,ed) =(UV(1,nz,ed)*UV_dis_tend( 1,nz,ed)+UV(2,nz,ed)*UV_dis_tend( 2,nz,ed))   
                uke_back(nz,ed)=(UV(1,nz,ed)*UV_back_tend(1,nz,ed)+UV(2,nz,ed)*UV_back_tend(2,nz,ed))
            END DO
        END DO
        
        DO  nz=1,nl-1
            uuu=0.0_WP
            uuu=uke_back(nz,:)
            call smooth_elem(uuu,dynamics%smooth_back, partit, mesh) !3) ?
            uke_back(nz,:)=uuu
        END DO
        
        !Timestepping use simple backward timestepping; all components should have dt in it, unless they need it twice
        !Amplitudes should be right given the correction of the viscosities; check for all, also for biharmonic
        !uke(1,ed)=uke(1,ed)-uke_dis(1,ed)-uke_back(1,ed)+uke_dif(1,ed)
        ed=myDim_elem2D+eDim_elem2D
        allocate(U_work(nl-1,myDim_nod2D+eDim_nod2D),V_work(nl-1,myDim_nod2D+eDim_nod2D))
        allocate(work_uv(myDim_nod2D+eDim_nod2D))
        allocate(rosb_array(nl-1,ed))
        call exchange_elem(UV, partit)
        rosb_array=0._WP
        DO nz=1, nl-1
            work_uv=0._WP
            DO node=1, myDim_nod2D
                vol=0._WP
                U_work(nz,node)=0._WP 
                V_work(nz,node)=0._WP 
                DO j=1, nod_in_elem2D_num(node)
                    elem=nod_in_elem2D(j, node)
                    U_work(nz,node)=U_work(nz,node)+UV(1,nz,elem)*elem_area(elem)
                    V_work(nz,node)=V_work(nz,node)+UV(2,nz,elem)*elem_area(elem)
                    vol=vol+elem_area(elem)
                END DO
                U_work(nz,node)=U_work(nz,node)/vol
                V_work(nz,node)=U_work(nz,node)/vol
            END DO
            work_uv=U_work(nz,:)
            call exchange_nod(work_uv, partit)
            U_work(nz,:)=work_uv
            work_uv=V_work(nz,:)
            call exchange_nod(work_uv, partit)
            V_work(nz,:)=work_uv    
        END DO
        
        DO el=1,myDim_elem2D
            DO nz=1, nlevels(el)-1     
                rosb_array(nz,el)=sqrt((sum(gradient_sca(1:3,el)*U_work(nz,elem2D_nodes(1:3,el)))-&
                    sum(gradient_sca(4:6, el)*V_work(nz,elem2D_nodes(1:3,el))))**2+&
                    (sum(gradient_sca(4:6, el)*U_work(nz,elem2D_nodes(1:3,el)))+&
                    sum(gradient_sca(1:3, el)*V_work(nz,elem2D_nodes(1:3,el))))**2)
                ! hall=hall+hnz
            END DO
            ! rosb_array(el)=rosb_array(el)/hall
        END DO
        
        DO ed=1, myDim_elem2D
            scaling=1._WP
            IF(dynamics%uke_scaling) then
                reso=sqrt(elem_area(ed)*4._wp/sqrt(3._wp))
                rosb=0._wp
                elnodes=elem2D_nodes(:, ed)   
                DO kk=1,3
                    c1=0._wp
                    nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(elnodes(kk)), elnodes(kk))), 1)
                    !Vertical average; same scaling in the vertical
                    DO nz=1, nzmax-1
                        c1=c1+hnode_new(nz,elnodes(kk))*(sqrt(max(bvfreq(nz,elnodes(kk)), 0._WP))+sqrt(max(bvfreq(nz+1,elnodes(kk)), 0._WP)))/2.
                    END DO
                    c1=max(c_min, c1/pi) !ca. first baroclinic gravity wave speed limited from below by c_min
                    !Cutoff K_GM depending on (Resolution/Rossby radius) ratio
                    rosb=rosb+min(c1/max(abs(mesh%coriolis_node(elnodes(kk))), f_min), r_max)
                END DO
                rosb=rosb/3._WP
                scaling=1._WP/(1._WP+(dynamics%uke_scaling_factor*reso/rosb))!(4._wp*reso/rosb))
            END IF
            
            DO nz=1, nlevels(ed)-1  
                elnodes=elem2D_nodes(:,ed)
                rosb_array(nz,ed)=rosb_array(nz,ed)/max(abs(sum(mesh%coriolis_node(elnodes(:)))), f_min)
                uke_dis(nz,ed)=scaling*1._WP/(1._WP+rosb_array(nz,ed)/dynamics%rosb_dis)*uke_dis(nz,ed)
            END DO
        END DO
        
        deallocate(U_work, V_work)
        deallocate(rosb_array)
        deallocate(work_uv)
        
        call exchange_elem(uke_dis, partit)
        DO nz=1, nl-1  
            uuu=uke_dis(nz,:)
            call smooth_elem(uuu,dynamics%smooth_dis, partit, mesh)
            uke_dis(nz,:)=uuu
        END DO

        if (dynamics%uke_advection) then
            call backscatter_adv(dynamics, partit, mesh)
            call exchange_elem(uke_adv, partit)

            DO ed=1, myDim_elem2D
                DO  nz=1,nlevels(ed)-1
                uke_rhs_old(nz,ed)=uke_rhs(nz,ed)
                uke_rhs(nz,ed)=-uke_dis(nz,ed)-uke_back(nz,ed)+uke_dif(nz,ed)-uke_adv(nz, ed)
                uke(nz,ed)=uke(nz,ed)+1.5_WP*uke_rhs(nz,ed)-0.5_WP*uke_rhs_old(nz,ed)
                END DO
            END DO
        else
            DO ed=1, myDim_elem2D
                DO  nz=1,nlevels(ed)-1
                uke_rhs_old(nz,ed)=uke_rhs(nz,ed)
                uke_rhs(nz,ed)=-uke_dis(nz,ed)-uke_back(nz,ed)+uke_dif(nz,ed)
                uke(nz,ed)=uke(nz,ed)+1.5_WP*uke_rhs(nz,ed)-0.5_WP*uke_rhs_old(nz,ed)
                END DO
            END DO
        end if
        
        call exchange_elem(uke, partit)
        deallocate(uuu)
        
    end subroutine uke_update
end module g_backscatter

