
module compute_CFLz_interface
    interface
        subroutine compute_CFLz(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine compute_CFLz
    end interface
end module compute_CFLz_interface

module compute_Wvel_split_interface
    interface
        subroutine compute_Wvel_split(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine compute_Wvel_split
    end interface
end module compute_Wvel_split_interface

module compute_vert_vel_transpv_interface
    interface        
        subroutine compute_vert_vel_transpv(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine compute_vert_vel_transpv
    end interface
end module compute_vert_vel_transpv_interface

module oce_ale_interfaces
    interface
        subroutine init_bottom_elem_thickness(partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine init_bottom_elem_thickness

        subroutine init_bottom_node_thickness(partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine init_bottom_node_thickness
        
        subroutine init_surface_elem_depth(partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine init_surface_elem_depth

        subroutine init_surface_node_depth(partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine init_surface_node_depth

        subroutine impl_vert_visc_ale(dynamics, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine impl_vert_visc_ale

        subroutine update_stiff_mat_ale(partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine update_stiff_mat_ale

        subroutine compute_ssh_rhs_ale(dynamics, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine compute_ssh_rhs_ale

        subroutine solve_ssh_ale(dynamics, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine solve_ssh_ale

        subroutine compute_hbar_ale(dynamics, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine compute_hbar_ale

        subroutine vert_vel_ale(dynamics, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine vert_vel_ale

        subroutine update_thickness_ale(partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine update_thickness_ale
    end interface
end module oce_ale_interfaces

module init_ale_interface
    interface
        subroutine init_ale(dynamics, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine init_ale
    end interface
end module init_ale_interface

module init_thickness_ale_interface
    interface
        subroutine init_thickness_ale(dynamics, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine init_thickness_ale
    end interface
end module init_thickness_ale_interface

module oce_timestep_ale_interface
    interface
        subroutine oce_timestep_ale(n, ice, dynamics, tracers, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        use MOD_DYN
        use MOD_ICE
        integer       , intent(in)            :: n
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_ice), intent(inout), target :: ice
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine oce_timestep_ale
    end interface
end module oce_timestep_ale_interface


! CONTENT:
! ------------
!    subroutine ale_init
!    subroutine init_bottom_elem_thickness
!    subroutine init_bottom_elem_thickness
!    subroutine init_thickness_ale
!    subroutine update_thickness_ale
!    subroutine restart_thickness_ale
!    subroutine init_stiff_mat_ale
!    subroutine update_stiff_mat_ale
!    subroutine compute_ssh_rhs_ale
!    subroutine compute_hbar_ale
!    subroutine vert_vel_ale
!    subroutine solve_ssh_ale
!    subroutine impl_vert_visc_ale
!    subroutine oce_timestep_ale
!
! initially written by Sergey Danilov, adapted and expanded by Patrick Scholz and
! Dmitry Sidorenko, 14.02.2019
!    
!===============================================================================
! allocate & initialise arrays for Arbitrary-Langrangian-Eularian (ALE) method
subroutine init_ale(dynamics, partit, mesh)
    USE o_PARAM
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE o_ARRAYS
!    USE g_config, only: which_ale, use_cavity, use_partial_cell

! kh 18.03.21
    USE g_config, only: which_ale, use_cavity, use_partial_cell, ib_async_mode

    USE g_forcing_param, only: use_virt_salt
    use oce_ale_interfaces
    Implicit NONE
     
! kh 18.03.21
    integer             :: i, j

    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh
    !___________________________________________________________________________
    integer                               :: n, nzmax, nzmin, elnodes(3), elem
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
nl => mesh%nl

    !___allocate________________________________________________________________
    ! hnode and hnode_new: layer thicknesses at nodes. 
    allocate(mesh%hnode(1:nl-1, myDim_nod2D+eDim_nod2D))
    allocate(mesh%hnode_new(1:nl-1, myDim_nod2D+eDim_nod2D))
    ! ssh_rhs_old: auxiliary array to store an intermediate part of the rhs computations.
    allocate(dynamics%ssh_rhs_old(myDim_nod2D+eDim_nod2D))
    dynamics%ssh_rhs_old = 0.0_WP
    
    ! hbar, hbar_old: correspond to the elevation, but on semi-integer time steps.
    allocate(mesh%hbar(myDim_nod2D+eDim_nod2D))
    allocate(mesh%hbar_old(myDim_nod2D+eDim_nod2D))
    
    ! helem: layer thickness at elements. It is interpolated from hnode.
    allocate(mesh%helem(1:nl-1, myDim_elem2D+eDim_nod2D))
    
    ! dhe: The increment of total fluid depth on elements. It is used to update the matrix
    ! of the ssh operator.      
    allocate(mesh%dhe(myDim_elem2D))
    
    allocate(mesh%zbar_3d_n(nl,myDim_nod2D+eDim_nod2D))
    
    ! Z_n: mid depth of layers due to ale thinkness variactions at ervery node n 
!    allocate(mesh%Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D)) 
! kh 18.03.21
    if (ib_async_mode == 0) then
        allocate(mesh%Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D)) 
        allocate(mesh%Z_3d_n_ib(nl-1,myDim_nod2D+eDim_nod2D)) 
        Z_3d_n(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)              => mesh%Z_3d_n(:,:)
        Z_3d_n_ib(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)           => mesh%Z_3d_n_ib(:,:)
        !allocate(Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D))
        !allocate(Z_3d_n_ib(nl-1,myDim_nod2D+eDim_nod2D))
    else
! kh 18.03.21 support "first touch" idea
!$omp parallel sections num_threads(2)
!$omp section
        allocate(mesh%Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D)) 
        Z_3d_n(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)              => mesh%Z_3d_n(:,:)
        !allocate(Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D))
        do i = 1, myDim_nod2D+eDim_nod2D
            do j = 1, nl-1
                Z_3d_n(j, i) = 0._WP
            end do
        end do
!$omp section
        allocate(mesh%Z_3d_n_ib(nl-1,myDim_nod2D+eDim_nod2D)) 
        Z_3d_n_ib(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)           => mesh%Z_3d_n_ib(:,:)
        !allocate(Z_3d_n_ib(nl-1,myDim_nod2D+eDim_nod2D))
        do i = 1, myDim_nod2D+eDim_nod2D
            do j = 1, nl-1
                Z_3d_n_ib(j, i) = 0._WP
            end do
        end do
!$omp end parallel sections
    end if
    
    ! bottom_elem_tickness: changed bottom layer thinkness due to partial cells
    allocate(mesh%bottom_elem_thickness(myDim_elem2D+eDim_nod2D))
    allocate(mesh%zbar_e_bot(myDim_elem2D+eDim_elem2D)) 
    allocate(mesh%zbar_e_srf(myDim_elem2D+eDim_elem2D)) 
    
    ! also change bottom thickness at nodes due to partial cell --> bottom 
    ! thickness at nodes is the volume weighted mean of sorounding elemental
    ! thicknesses
    allocate(mesh%bottom_node_thickness(myDim_nod2D+eDim_nod2D))
    allocate(mesh%zbar_n_bot(myDim_nod2D+eDim_nod2D)) 
    allocate(mesh%zbar_n_srf(myDim_nod2D+eDim_nod2D)) 

    ! reassociate after the allocation (no pointer exists before)
#include "associate_mesh_ass.h"
 
    !___initialize______________________________________________________________
    hbar      = 0.0_WP
    hbar_old  = 0.0_WP
    dhe       = 0.0_WP
    hnode     = 0.0_WP
    hnode_new = 0.0_WP
    helem     = 0.0_WP
    !___________________________________________________________________________
    ! calculate thickness of partial bottom layer cells as well as depth depth
    ! of partial cell bootom layer
    zbar_n_bot = 0.0
    zbar_e_bot = 0.0
    call init_bottom_elem_thickness(partit, mesh)
    call init_bottom_node_thickness(partit, mesh)
    
    ! compute depth of partial cell ocean-cavity interface
    zbar_n_srf = zbar(1)
    zbar_e_srf = zbar(1)
    call init_surface_elem_depth(partit, mesh)
    call init_surface_node_depth(partit, mesh)
    
    !___________________________________________________________________________
    ! initialise 3d field of depth levels and mid-depth levels
    zbar_3d_n  = 0.0_WP
    Z_3d_n     = 0.0_WP
    do n=1,myDim_nod2D+eDim_nod2D 
        ! max. number of levels at node n
        nzmin=ulevels_nod2D(n)
        nzmax=nlevels_nod2D(n)
        
        !_______________________________________________________________________
        ! create dummy zbar full depth levels within cavity --> need to compute 
        ! cavity pressure b oundary condition
        zbar_3d_n(1:nzmin-1,n)=zbar(1:nzmin-1);
        
        !_______________________________________________________________________
        ! in case of partial cells and use_cavity surface depth is different from zbar(nzmin)
        zbar_3d_n(nzmin,n)=zbar_n_srf(n);
        
        zbar_3d_n(nzmin+1:nzmax-1,n)=zbar(nzmin+1:nzmax-1);
        
        ! in case of partial cells bottom depth is different from zbar(nzmax)
        zbar_3d_n(nzmax,n)=zbar_n_bot(n);
        
        !_______________________________________________________________________
        ! create dummy Z mid depth levels within cavity  --> need to compute 
        ! cavity pressure b oundary condition
        Z_3d_n(1:nzmin-1,n)=Z(1:nzmin-1);
        
        !_______________________________________________________________________
        ! in case of partial cells bottom mid depth is different from Z(nzmax-1)
        Z_3d_n(nzmin,n) =zbar_3d_n(nzmin,n)+(zbar_3d_n(nzmin+1,n)-zbar_n_srf(n))/2;
        
        Z_3d_n(nzmin+1:nzmax-2,n) =Z(nzmin+1:nzmax-2);
        
        ! in case of partial cells bottom mid depth is different from Z(nzmax-1)
        Z_3d_n(nzmax-1,n) =zbar_3d_n(nzmax-1,n)+(zbar_n_bot(n)-zbar_3d_n(nzmax-1,n))/2;
        
    end do

end subroutine init_ale
!
!
!===============================================================================
subroutine init_bottom_elem_thickness(partit, mesh)
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_ARRAYS
    use g_config,only: use_partial_cell, partial_cell_thresh, use_depthonelem
    use g_comm_auto
    use g_support
    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer       :: elem, elnodes(3), nle
    real(kind=WP) :: dd
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    
    !___________________________________________________________________________
    ! If we use partial cells, the thickness of bottom cell is adjusted.
    ! The adjustment is limited. It cannot be more than + (1/2) of the deeper
    ! layer, nor -(1/2) of the current layer. 
    if(use_partial_cell) then 
        !Adjust the thickness of elemental bottom cells
        do elem=1, myDim_elem2D
            
            !___________________________________________________________________
            if (use_depthonelem) then 
                dd=depth(elem)
            else
                elnodes=elem2D_nodes(:,elem) 
                ! elemental topographic depth
                dd=sum(depth(elnodes))/3.0_WP
            end if 
            
            ! number of full depth levels at elem
            nle=nlevels(elem)
            
            !___________________________________________________________________
            ! Only apply Partial Cells when the initial full cell bottom
            ! layer thickness is above the treshhold partial_cell_thresh
            if (zbar(nle-1)-zbar(nle)<=partial_cell_thresh) then
                zbar_e_bot(elem) = zbar(nle)
                bottom_elem_thickness(elem)=zbar(nle-1)-zbar_e_bot(elem)
                cycle
            end if 
                
            !___________________________________________________________________
            ! if topographic depth dd is deeper than depth of deepest full cell 
            ! depth level zbar(nle)
            !       : 
            !       : 
            ! ______________ zbar(nle-1)--------->+---->+
            !                                     |     |
            !                                     |     |
            ! -------------- Z(nle-1)             |--case1--> dd1=
            !                                     |     |
            !                                     |     |
            ! ______________ zbar(nle)            |     |--case2--> dd1 = 
            ! / / / / / / /                       |     |
            !  / / o dd case1 ------------------->+     |
            ! -------------- Z(nle)(mid-depth)--------->+
            !  / / / / / / /
            ! / /  o dd case2
            !  / / / / / /
            if(dd<zbar(nle)) then 
                if(nle==nl) then
                    zbar_e_bot(elem) = max(dd,zbar(nle)+(zbar(nle)-Z(nle-1)))
                    
                else
                    ! case 1 : max(Z(nle),dd) = dd
                    ! case 2 : max(Z(nle),dd) = Z(nle)
                    zbar_e_bot(elem) = max(Z(nle),dd)
                end if
                bottom_elem_thickness(elem)=zbar(nle-1)-zbar_e_bot(elem)
            !___________________________________________________________________
            ! if topographic depth dd is shallower than depth of deepest full cell 
            ! depth level zbar(nle)
            !        : 
            !        : 
            ! ______________ zbar(nle-1)--------->+---->+
            !                                     |--dd case1--> zbar_n_bot=
            !      o dd case1                     |     |
            ! -------------- Z(nle-1)(mid-depth)->+     |--dd case 2 --> zbar_n_bot=
            !      o dd case2 ------------------------->+
            ! ______________ zbar(nle) 
            ! / / / / / / / 
            !  / / / / / / /
            ! / / / / / / / 
            else
!!PS                 !_______________________________________________________________
!!PS                 ! if a thicker partial bottom layer thickness is more realistic than
!!PS                 ! always apply it, BUT when a thinner bottom layer thickness is more 
!!PS                 ! realistic than only apply it when the initial full cell bottom
!!PS                 ! layer thickness is above the treshhold partial_cell_thresh to 
!!PS                 ! not allow already thin layers to become even thinner
!!PS                 if (zbar(nle-1)-zbar(nle)<=partial_cell_thresh) then
!!PS                     zbar_e_bot(elem) = zbar(nle)
!!PS                     bottom_elem_thickness(elem)=zbar(nle-1)-zbar_e_bot(elem)
!!PS                     cycle
!!PS                 end if     
                
                ! case 1 : min(Z(nle-1),dd) = Z(nle-1)
                ! case 2 : min(Z(nle-1),dd) = dd
                zbar_e_bot(elem) = min(Z(nle-1),dd)
                bottom_elem_thickness(elem)=zbar(nle-1)-zbar_e_bot(elem)
            end if  
        end do ! --> do elem=1, myDim_elem2D
        
    !___________________________________________________________________________
    ! use full bottom cells
    else
        do elem=1, myDim_elem2D
            nle=nlevels(elem)
            bottom_elem_thickness(elem)=zbar(nle-1)-zbar(nle)
            zbar_e_bot(elem) = zbar(nle)
        end do
    end if 
    
    !___________________________________________________________________________
    call exchange_elem(zbar_e_bot, partit)
    call exchange_elem(bottom_elem_thickness, partit)
    
end subroutine init_bottom_elem_thickness
!
!
!===============================================================================
subroutine init_bottom_node_thickness(partit, mesh)
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_ARRAYS
    use g_config,only: use_partial_cell
    use g_comm_auto
    use g_support
    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer       :: node, nln, elem, elemi, nelem
    real(kind=WP) :: dd 
    real(kind=WP) :: hnbot, tvol 
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    
    !___________________________________________________________________________
    ! If we use partial cells, the thickness of bottom cell is adjusted.
    ! The adjustment is limited. It cannot be more than + (1/2) of the deeper
    ! layer, nor -(1/2) of the current layer. 
    if(use_partial_cell) then 
        !Adjust the thickness of nodal bottom cells
        do node=1, myDim_nod2D
            
!!PS             !___________________________________________________________________
!!PS             ! nodal topographic depth must be as least as deep as deepest bottom depth 
!!PS             ! of the sorounding elements          
!!PS             dd = depth(node)
!!PS             
!!PS             ! number of full depth levels at node
!!PS             nln  = nlevels_nod2D(node)
!!PS             
!!PS             !___________________________________________________________________
!!PS             ! if topographic depth dd is deeper than depth of deepest full cell 
!!PS             ! depth level zbar(nle)
!!PS             !       : 
!!PS             !       : 
!!PS             ! ______________ zbar(nle-1)--------->+---->+
!!PS             !                                     |     |
!!PS             !                                     |     |
!!PS             ! -------------- Z(nle-1)             |--case1--> zbar_n_bot=
!!PS             !                                     |     |
!!PS             !                                     |     |
!!PS             ! ______________ zbar(nle)            |     |--case2--> zbar_n_bot = 
!!PS             ! / / / / / / /                       |     |
!!PS             !  / / o dd case1 ------------------->+     |
!!PS             ! -------------- Z(nle)(mid-depth)--------->+
!!PS             !  / / / / / / /
!!PS             ! / /  o dd case2
!!PS             !  / / / / / /
!!PS             if(dd<zbar(nln)) then 
!!PS                 if(nln==nl) then
!!PS                     zbar_n_bot(node) = max(dd,zbar(nln)+(zbar(nln)-Z(nln-1)))
!!PS                     
!!PS                 else
!!PS                     ! case 1 : max(Z(nle),dd) = dd
!!PS                     ! case 2 : max(Z(nle),dd) = Z(nle)
!!PS                     zbar_n_bot(node) = max(Z(nln),dd)
!!PS                 end if
!!PS             !___________________________________________________________________
!!PS             ! if topographic depth dd is shallower than depth of deepest full cell 
!!PS             ! depth level zbar(nle)
!!PS             !        : 
!!PS             !        : 
!!PS             ! ______________ zbar(nle-1)--------->+---->+
!!PS             !                                     |--dd case1--> zbar_n_bot=
!!PS             !      o dd case1                     |     |
!!PS             ! -------------- Z(nle-1)(mid-depth)->+     |--dd case 2 --> zbar_n_bot=
!!PS             !      o dd case2 ------------------------->+
!!PS             ! ______________ zbar(nle) 
!!PS             ! / / / / / / / 
!!PS             !  / / / / / / /
!!PS             ! / / / / / / / 
!!PS             else
!!PS                 ! case 1 : min(Z(nle-1),dd) = Z(nle-1)
!!PS                 ! case 2 : min(Z(nle-1),dd) = dd
!!PS                 zbar_n_bot(node) = min(Z(nln-1),dd)
!!PS                 
!!PS             end if        
            !___________________________________________________________________
            ! compute vertice partial bottom depth from the deepest sorounding
            ! elemental partial bottom depths
            nln   = nlevels_nod2D(node)
            nelem = nod_in_elem2d_num(node)
            zbar_n_bot(node)           = minval(zbar_e_bot(nod_in_elem2d(1:nelem,node))) 
            bottom_node_thickness(node)= zbar(nln-1)-zbar_n_bot(node)
        end do ! --> do node=1, myDim_nod2D+eDim_nod2D
        
    !___________________________________________________________________________
    ! use full bottom cells
    else
        do node=1,myDim_nod2D
            nln = nlevels_nod2D(node)
            zbar_n_bot(node)           = zbar(nln)
            bottom_node_thickness(node)= zbar(nln-1)-zbar_n_bot(node)
        end do
    end if ! --> if(use_partial_cell) then 

    !___________________________________________________________________________
    call exchange_nod(zbar_n_bot, partit)
    call exchange_nod(bottom_node_thickness, partit)
    
end subroutine init_bottom_node_thickness
!
!
!===============================================================================
subroutine init_surface_elem_depth(partit, mesh)
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_ARRAYS
    use g_config,only: use_cavity, use_cavity_partial_cell, cavity_partial_cell_thresh, use_cavityonelem
    use g_comm_auto
    use g_support
    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer       :: elem, elnodes(3), ule
    real(kind=WP) :: dd
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    
    if (use_cavity) then 
        
        !_______________________________________________________________________
        ! If we use partial cells and cavity, the thickness of surface cell is adjusted.
        ! The adjustment is limited. It cannot be more than + (1/2) of the deeper
        ! layer, nor -(1/2) of the current layer. 
        ! Adjust the thickness of elemental surface cells under the cavity
        do elem=1, myDim_elem2D
            !___________________________________________________________________
            ule=ulevels(elem)    
            if (ule==1) cycle
            
            !___________________________________________________________________
            ! elemental cavity depth
            if (use_cavity_partial_cell) then 
            
                !_______________________________________________________________
                if (use_cavityonelem) then 
                    dd=cavity_depth(elem)
                else
                    elnodes=elem2D_nodes(:,elem) 
                    ! elemental cavity depth
                    dd=sum(cavity_depth(elnodes))/3.0_WP
                end if 
                
                !_______________________________________________________________
                ! Only apply Surface Partial Cells when the initial full cell surface
                ! layer thickness is above the treshhold cavity_partial_cell_thresh
                if (zbar(ule)-zbar(ule+1)<=cavity_partial_cell_thresh) then
                    zbar_e_srf(elem) = zbar(ule)
                    cycle
                end if         
                
                if(dd<zbar(ule)) then 
!!PS                     !_______________________________________________________________
!!PS                     ! if a thicker partial surface layer thickness is more realistic than
!!PS                     ! always apply it, BUT when a thinner surface layer thickness is more 
!!PS                     ! realistic than only apply it when the initial full cell surface
!!PS                     ! layer thickness is above the treshhold cavity_partial_cell_thresh to 
!!PS                     ! not allow already thin layers to become even thinner
!!PS                     if (zbar(ule)-zbar(ule+1)<=cavity_partial_cell_thresh) then
!!PS                         zbar_e_srf(elem) = zbar(ule)
!!PS                     else
!!PS                         zbar_e_srf(elem) = max(Z(ule),dd)
!!PS                     end if
                    zbar_e_srf(elem) = max(Z(ule),dd)
                else
                    zbar_e_srf(elem) = min(Z(ule-1),dd)
                end if
            else
                zbar_e_srf(elem) = zbar(ule)
            end if 
                
        end do ! --> do elem=1, myDim_elem2D
        
        !_______________________________________________________________________
        call exchange_elem(zbar_e_srf, partit)
    end if 
end subroutine init_surface_elem_depth
!
!
!===============================================================================
subroutine init_surface_node_depth(partit, mesh)
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_ARRAYS
    use g_config,only:  use_cavity, use_cavity_partial_cell
    use g_comm_auto
    use g_support
    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer       :: node, uln, nelem, elemi
    real(kind=WP) :: dd 
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    
    !___________________________________________________________________________
    if (use_cavity) then 
        !___________________________________________________________________________
        ! If we use partial cells and cavity, the thickness of surface cell is adjusted.
        ! The adjustment is limited. It cannot be more than + (1/2) of the deeper
        ! layer, nor -(1/2) of the current layer. 
        !Adjust the thickness of nodal surface cells
        do node=1, myDim_nod2D
            !___________________________________________________________________
            ! number of full depth levels at node
            uln  = ulevels_nod2D(node)
            if (uln==1) cycle
            
            !___________________________________________________________________
            ! nodal cavity depth  
            if (use_cavity_partial_cell) then 
!!PS                 dd = cavity_depth(node)
!!PS                 if(dd<zbar(uln)) then 
!!PS                     zbar_n_srf(node) = max(Z(uln),dd)
!!PS                 else
!!PS                     zbar_n_srf(node) = min(Z(uln-1),dd)
!!PS                 end if  
                
                nelem =  nod_in_elem2d_num(node)
                zbar_n_srf(node)=maxval(zbar_e_srf(nod_in_elem2d(1:nelem,node))) 
                
            else
                zbar_n_srf(node) = zbar(uln)
                
            end if 
        end do ! --> do node=1, myDim_nod2D+eDim_nod2D
        
        !_______________________________________________________________________
        call exchange_nod(zbar_n_srf, partit)
    end if 
end subroutine init_surface_node_depth
!
!
!===============================================================================
! initialize thickness arrays based on the current hbar 
subroutine init_thickness_ale(dynamics, partit, mesh)
! For z-star case: we stretch scalar thicknesses (nodal) 
! through nlevels_nod2D_min -2 layers. Layer nlevels_nod2D_min-1
! should not be touched if partial cell is implemented (it is).
! In lower layers scalar prisms are modified by the bottom.  
! Important: nlevels_nod2D_min has to be allocated and filled. 
    use g_config,only: dt, which_ale
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use g_comm_auto
    implicit none
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer       :: n, nz, elem, elnodes(3), nzmin, nzmax
    real(kind=WP) :: dd 
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer :: ssh_rhs_old, eta_n
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    ssh_rhs_old=>dynamics%ssh_rhs_old(:)
    eta_n      =>dynamics%eta_n(:)
    
    !___________________________________________________________________________
    
    if(mype==0) then
        write(*,*) '____________________________________________________________'
        write(*,*) ' --> initialise ALE layerthicknesses, depth levels and middepth levels'
        write(*,*)
    end if
    !___________________________________________________________________________
    ! Fill in ssh_rhs_old
    !!PS ssh_rhs_old=(hbar-hbar_old)*area(1,:)/dt
    do n=1,myDim_nod2D+eDim_nod2D
        ssh_rhs_old(n)=(hbar(n)-hbar_old(n))*areasvol(ulevels_nod2D(n),n)/dt  ! --> TEST_cavity
    end do
    
    ! -->see equation (14) FESOM2:from finite elements to finie volume
    eta_n=alpha*hbar_old+(1.0_WP-alpha)*hbar   
    
    if     (trim(which_ale)=='linfs') then
        !_______________________________________________________________________
        ! >->->->->->->->->->->->     linear-free-surface    <-<-<-<-<-<-<-<-<-<
        !_______________________________________________________________________
        ! no layer thickness variation in any layer
        do n=1,myDim_nod2D+eDim_nod2D
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D(n)-1
            !!PS hnode(1,n)=(zbar(1)-zbar(2))
            !!PS do nz=2,nlevels_nod2D(n)-2
            do nz=nzmin,nzmax-1
                !!PS hnode(nz,n)=(zbar(nz)-zbar(nz+1))
                hnode(nz,n)=(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
            end do      
            
            ! set bottom node thickness
            !!PS hnode(nlevels_nod2D(n)-1,n)=bottom_node_thickness(n)
            hnode(nzmax,n)=bottom_node_thickness(n)
            
            !!PS do nz=nlevels_nod2D(n),nl-1
            !!PS --> can skip this, hnode(:,:) is initialised with 0.0_WP
            !!PS do nz=nzmax+1,nl-1
            !!PS     hnode(nz,n)=0.0_WP
            !!PS end do
        end do
        
        do elem=1,myDim_elem2D
            nzmin = ulevels(elem)
            nzmax = nlevels(elem)-1
            
            !!PS do nz=1,nlevels(elem)-2
            helem(nzmin,elem)=(zbar_e_srf(elem)-zbar(nzmin+1))
            do nz = nzmin+1, nzmax-1
                helem(nz,elem)=(zbar(nz)-zbar(nz+1))
            end do
            
            ! set bottom elem thickness
            !!PS helem(nlevels(elem)-1,elem)=bottom_elem_thickness(elem)
            helem(nzmax,elem)=bottom_elem_thickness(elem)
            
            !!PS do nz=nlevels(elem),nl-1
            !!PS --> can skip this, helem(:,:) is initialised with 0.0_WP
            !!PS do nz=nzmax+1,nl-1
            !!PS     helem(nz,elem)=0.0_WP
            !!PS end do
        end do
        
    elseif (trim(which_ale)=='zlevel') then
        !_______________________________________________________________________
        ! >->->->->->->->->->->->->->->     z-level    <-<-<-<-<-<-<-<-<-<-<-<-<
        !_______________________________________________________________________
        ! --> include all ssh variations into the top layer 
        do n=1,myDim_nod2D+eDim_nod2D
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D(n)-1
            
            ! put all ssh variation (hbar) into first layer 
            !!PS hnode(1,n)=hbar(n)+(zbar(1)-zbar(2))
            if (nzmin == 1) then 
                ! only allow open ocean to move with ssh
                !!PS hnode(nzmin,n)=hbar(n)+(zbar(nzmin)-zbar(nzmin+1))
                hnode(nzmin,n)=hbar(n)+(zbar_3d_n(nzmin,n)-zbar_3d_n(nzmin+1,n))
            else
                ! in case of cavity no movement with ssh, cavity-ocean boundary is fixed
                hnode(nzmin,n)=(zbar_3d_n(nzmin,n)-zbar_3d_n(nzmin+1,n))
            endif 
            
            ! leave lower levels untouched
            !!PS do nz=2,nlevels_nod2D(n)-2
            do nz=nzmin+1,nzmax-1
                !!PS hnode(nz,n)=(zbar(nz)-zbar(nz+1))
                hnode(nz,n)=(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
            end do 
            
            ! set bottom node thickness
            !!PS hnode(nlevels_nod2D(n)-1,n)=bottom_node_thickness(n)
            hnode(nzmax,n)=bottom_node_thickness(n)
            
            !!PS --> can skip this, hnode(:,:) is initialised with 0.0_WP
            !!PS ! layer thickness of bottom layer equal 0
            !!PS do nz=nlevels_nod2D(n),nl-1
            !!PS     hnode(nz,n)=0.0_WP
            !!PS end do
        end do
        
        do elem=1,myDim_elem2D
            nzmin = ulevels(elem)
            nzmax = nlevels(elem)-1
            
            elnodes=elem2D_nodes(:,elem) 
            
            ! interpolated ssh variation at element elem
            dhe(elem)=sum(hbar(elnodes))/3.0_WP
            
            ! store elemtal ssh varition only in first layer
            !!PS helem(1,elem)=dhe(elem)+(zbar(1)-zbar(2))
            if (nzmin==1) then 
                helem(nzmin,elem)=dhe(elem)+(zbar_e_srf(elem)-zbar(nzmin+1))
            else
                helem(nzmin,elem)=(zbar_e_srf(elem)-zbar(nzmin+1))
            end if 
            
            ! lower layers leave untouched 
            !!PS do nz=2,nlevels(elem)-2
            do nz=nzmin+1,nzmax-1
                helem(nz,elem)=(zbar(nz)-zbar(nz+1))
            end do
            
            ! elemental bottom layer thickness
            !!PS helem(nlevels(elem)-1,elem)=bottom_elem_thickness(elem)
            helem(nzmax,elem)=bottom_elem_thickness(elem)
            
            !!PS --> can skip this, helem(:,:) is initialised with 0.0_WP
            !!PS ! fill thickness below bottom layer
            !!PS do nz=nlevels(elem),nl-1
            !!PS     helem(nz,elem)=0.0_WP
            !!PS end do
            
        end do
        
    elseif (trim(which_ale)=='zstar' ) then
        !_______________________________________________________________________
        ! >->->->->->->->->->->->->->->     z-star     <-<-<-<-<-<-<-<-<-<-<-<-<
        !_______________________________________________________________________
        ! --> calcualte layer thinkness at depth layer and node
        do n=1, myDim_nod2D+eDim_nod2D
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D(n)-1
            
            if (nzmin==1) then 
                ! depth anomaly until the last minus one level where the scalar 
                ! prism is not intersected with bottom.
                !!PS dd=zbar(1)-zbar(nlevels_nod2D_min(n)-1)
                dd=zbar(nzmin)-zbar(nlevels_nod2D_min(n)-1)  
                
                ! calc layer thinkness for depth layer nz and node n. distribute 
                ! hbar surface elevation linear over verical column
                !!PS do nz=1,nlevels_nod2D_min(n)-2
                do nz=nzmin,nlevels_nod2D_min(n)-2
                    hnode(nz,n)=(zbar(nz)-zbar(nz+1))*(1.0_WP+hbar(n)/dd)
                end do
                
                ! do not distribute hbar into cells that intersect somehow with 
                ! bottom layer 
                !!PS do nz=nlevels_nod2D_min(n)-1, nlevels_nod2D(n)-1
                do nz=nlevels_nod2D_min(n)-1, nzmax-1
                    hnode(nz,n)=(zbar(nz)-zbar(nz+1))
                end do
            else
                ! in case of cavity dont distribute ssh --> cavity-ocean boudnary
                ! is fixed
                do nz=nzmin,nzmax-1
                    hnode(nz,n)=(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
                end do 
            end if 
            ! set bottom node thickness
            hnode(nzmax,n)=bottom_node_thickness(n)
            
            !!PS --> can skip this, hnode(:,:) is initialised with 0.0_WP
            !!PS ! layer thickness of bottom layer equal 0
            !!PS do nz=nlevels_nod2D(n),nl-1
            !!PS     hnode(nz,n)=0.0_WP
            !!PS end do
        end do
        
        !_______________________________________________________________________
        ! --> calculate mean layer thinkness at element
        do elem=1, myDim_elem2D
            nzmin = ulevels(elem)
            nzmax = nlevels(elem)-1
            
            elnodes=elem2D_nodes(:, elem)
            
            ! interpolated ssh variation at element elem
            if (nzmin==1) then 
                dhe(elem)=sum(hbar(elnodes))/3.0_WP
                helem(nzmin,elem)=sum(hnode(nzmin,elnodes))/3.0_WP
            else
                dhe = 0.0_WP
                helem(nzmin,elem)=zbar_e_srf(elem)-zbar(nzmin+1)
            end if 
            
            !!PS do nz=1,nlevels(elem)-2
            !!PS do nz=nzmin,nzmax-1
            do nz=nzmin+1,nzmax-1
                helem(nz,elem)=sum(hnode(nz,elnodes))/3.0_WP
            end do
            
            ! elemental bottom layer thickness
            !!PS helem(nlevels(elem)-1,elem)=bottom_elem_thickness(elem)
            helem(nzmax,elem)=bottom_elem_thickness(elem)
            
            !!PS --> can skip this, helem(:,:) is initialised with 0.0_WP 
            !!PS ! fill thickness below bottom layer
            !!PS do nz=nlevels(elem),nl-1
            !!PS     helem(nz,elem)=0.0_WP
            !!PS end do
        end do
    else
        if (mype==0) then
        write(*,*)
        write(*,*) '____________________________________________________________'
        write(*,*) 'The vertical ALE discretisation ', which_ale,' is currently not supported!!!'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
        end if 
    endif
    
    !___________________________________________________________________________
    hnode_new=hnode  ! Should be initialized, because only variable part is updated.
    call exchange_elem(helem, partit)
    !!PS call check_total_volume(partit, mesh)
    
end subroutine init_thickness_ale
!
!
!===============================================================================
! update thickness arrays based on the current hbar 
subroutine update_thickness_ale(partit, mesh)
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_ARRAYS
    use g_config,only: which_ale,lzstar_lev,min_hnode
    use diagnostics, only: ldiag_DVD 
    use g_comm_auto

    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer :: n, nz, elem, elnodes(3),nzmax, nzmin
    integer      , dimension(:), allocatable :: idx
    real(kind=WP)                            :: rescue_hnode_old(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    ! >->->->->->->->->->->->->->->       z-level      <-<-<-<-<-<-<-<-<-<-<-<-<
    !___________________________________________________________________________
    if     (trim(which_ale)=='zlevel') then
        
        !_______________________________________________________________________
        ! idx is only needed for local star case to estimate over how much 
        ! depth layers hnode, depthlevel and mid-depthlevel need to be updated
        allocate(idx(lzstar_lev))
        
        ! if lzstar_lev=4 --> idx = /1,2,3,4/
        idx = (/(nz, nz=1, lzstar_lev, 1)/)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax, elem, elnodes)
!$OMP DO        
        !_______________________________________________________________________
        do elem=1,myDim_elem2D
            elnodes=elem2D_nodes(:, elem)
            nzmin = ulevels(elem)
            nzmax = nlevels(elem)-1
            
            !___________________________________________________________________
            ! if there is a cavity surface layer thickness is not update, its 
            ! kept fixed 
            if (nzmin > 1) cycle
            
            !___________________________________________________________________
            ! actualize elemental layer thinkness in first lzstar_lev layers
            if (any(hnode_new(nzmin+1:nzmin+lzstar_lev-1,elnodes(1)) - hnode(nzmin+1:nzmin+lzstar_lev-1,elnodes(1))/=0.0_WP) .or. &
                any(hnode_new(nzmin+1:nzmin+lzstar_lev-1,elnodes(2)) - hnode(nzmin+1:nzmin+lzstar_lev-1,elnodes(2))/=0.0_WP) .or. &
                any(hnode_new(nzmin+1:nzmin+lzstar_lev-1,elnodes(3)) - hnode(nzmin+1:nzmin+lzstar_lev-1,elnodes(3))/=0.0_WP)      &
                ) then
                ! --> case local zstar
                ! try to limitate over how much layers i realy need to distribute
                ! the change in ssh, so that the next loops run only over the 
                ! nesseccary levels and not over all lzstar_lev levels
                nz    = max(nzmin, maxval(pack(nzmin+idx-1,hnode_new(nzmin+1:nzmin+lzstar_lev-1, elnodes(1)) - hnode(nzmin+1:nzmin+lzstar_lev-1, elnodes(1))/=0.0_WP)))
                nz    = max(nz   , maxval(pack(nzmin+idx-1,hnode_new(nzmin+1:nzmin+lzstar_lev-1, elnodes(2)) - hnode(nzmin+1:nzmin+lzstar_lev-1, elnodes(2))/=0.0_WP)))
                nz    = max(nz   , maxval(pack(nzmin+idx-1,hnode_new(nzmin+1:nzmin+lzstar_lev-1, elnodes(3)) - hnode(nzmin+1:nzmin+lzstar_lev-1, elnodes(3))/=0.0_WP)))
                nzmax = min(nz   , nzmax-1)
                do nz=nzmin,nzmax
                    helem(nz,elem)=sum(hnode_new(nz,elnodes))/3.0_WP
                end do    
            !___________________________________________________________________
            ! only actualize elemental layer thickness in first layer 
            else
                ! --> case normal zlevel
                helem(nzmin,elem)=sum(hnode_new(nzmin,elnodes))/3.0_WP
            end if
        end do
!$OMP END DO 
        !_______________________________________________________________________
!$OMP DO 
        do n=1,myDim_nod2D+eDim_nod2D
            nzmin = ulevels_nod2D_max(n)
            nzmax = nlevels_nod2D_min(n)-1
            
            !___________________________________________________________________
            ! if there is a cavity surface layer thickness is not update, its 
            ! kept fixed 
            if (nzmin > 1) cycle
            
            !___________________________________________________________________
            ! actualize layer thinkness in first lzstar_lev layers
            if ( (any(hnode_new(nzmin+1:nzmin+lzstar_lev-1,n)-hnode(nzmin+1:nzmin+lzstar_lev-1,n)/=0.0_WP)) ) then
                ! --> case local zstar 
                ! try to limitate over how much layers i realy need to distribute
                ! the change in ssh, so that the next loops run only over the 
                ! nesseccary levels and not over all lzstar_lev levels
                !!PS nz = max(1,maxval(pack(idx,hnode_new(1:lzstar_lev,n)-hnode(1:lzstar_lev,n)/=0.0_WP)))
                nz = max(nzmin,maxval(pack(nzmin+idx-1,hnode_new(nzmin:nzmin+lzstar_lev-1,n)-hnode(nzmin:nzmin+lzstar_lev-1,n)/=0.0_WP)))
                
                ! nlevels_nod2D_min(n)-1 ...would be hnode of partial bottom cell but this
                ! one is not allowed to change so go until nlevels_nod2D_min(n)-2
                nzmax = min(nz,nzmax-1)
                ! do not touch zbars_3d_n that are involved in the bottom cell !!!!
                ! this ones are set up during initialisation and are not touched afterwards
                ! --> nlevels_nod2D_min(n),nlevels_nod2D_min(n)-1
                do nz=nzmax,nzmin,-1
                    ! --> case normal zlevel
                    ! at this point hnode = hode^(n) becomes hnode^(n+1) but i need the 
                    ! hnode^(n) value for the DVD analysis in call diagnostic() therefore 
                    ! i need to rescue this value before it gets dumped 
                    rescue_hnode_old(nz) = hnode(nz,n) 
                
                    hnode(nz,n)     = hnode_new(nz,n)
                    zbar_3d_n(nz,n) = zbar_3d_n(nz+1,n)+hnode_new(nz,n)
                    Z_3d_n(nz,n)    = zbar_3d_n(nz+1,n)+hnode_new(nz,n)/2.0_WP
                end do    
                
                ! save hnode^(n) now in the variable hnode_new --> be carefull this 
                ! might be confusing here for later ! --> it will only happen when 
                ! the DVD diagnostic is active !
                if (ldiag_DVD) then 
                    do nz=nzmax,nzmin,-1
                        hnode_new(nz,n) = rescue_hnode_old(nz)
                    end do
                end if 
                    
            !___________________________________________________________________
            ! only actualize layer thinkness in first layer 
            else
                ! --> case normal zlevel
                ! at this point hnode = hode^(n) becomes hnode^(n+1) but i need the 
                ! hnode^(n) value for the DVD analysis in call diagnostic() therefore 
                ! i need to rescue this value before it gets dumped 
                rescue_hnode_old(nzmin) = hnode(nzmin,n) 
                
                hnode(nzmin,n)    = hnode_new(nzmin,n)
                zbar_3d_n(nzmin,n)= zbar_3d_n(nzmin+1,n)+hnode_new(nzmin,n)
                Z_3d_n(nzmin,n)   = zbar_3d_n(nzmin+1,n)+hnode_new(nzmin,n)/2.0_WP
                
                ! save hnode^(n) now in the variable hnode_new --> be carefull this 
                ! might be confusing here for later ! --> it will only happen when 
                ! the DVD diagnostic is active !
                if (ldiag_DVD) then 
                    hnode_new(nzmin,n) = rescue_hnode_old(nzmin)
                end if 
               
            end if
        end do
!$OMP END DO
!$OMP END PARALLEL        
        !_______________________________________________________________________
        deallocate(idx)
        
    !___________________________________________________________________________
    ! >->->->->->->->->->->->->->->       z-star       <-<-<-<-<-<-<-<-<-<-<-<-<
    !___________________________________________________________________________
    elseif (trim(which_ale)=='zstar' ) then
        
        ! --> update layer thinkness, depth layer  and mid-depth layer at node
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax)
        do n=1, myDim_nod2D+eDim_nod2D
            ! actualize 3d depth levels and mid-depth levels from bottom to top
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D_min(n)-2
            
            !___________________________________________________________________
            ! if there is a cavity layer thickness is not updated, its 
            ! kept fixed 
            if (nzmin > 1) cycle
            
            !___________________________________________________________________
            ! do not touch zbars_3d_n that are involved in the bottom cell !!!!
            ! --> nlevels_nod2D_min(n),nlevels_nod2D_min(n)-1
            do nz=nzmax, nzmin,-1
                ! at this point hnode = hode^(n) becomes hnode^(n+1) but i need the 
                ! hnode^(n) value for the DVD analysis in call diagnostic() therefore 
                ! i need to rescue this value before it gets dumped 
                rescue_hnode_old(nz) = hnode(nz,n) 
                
                hnode(nz,n)     = hnode_new(nz,n)
                zbar_3d_n(nz,n) = zbar_3d_n(nz+1,n) + hnode_new(nz,n)
                Z_3d_n(nz,n)    = zbar_3d_n(nz+1,n) + hnode_new(nz,n)/2.0_WP
            end do
            
            !___________________________________________________________________
            ! save hnode^(n) now in the variable hnode_new --> be carefull this 
            ! might be confusing here for later ! --> it will only happen when 
            ! the DVD diagnostic is active !
            if (ldiag_DVD) then 
                do nz=nzmax, nzmin,-1
                    hnode_new(nz,n) = rescue_hnode_old(nz)
                end do
            end if 
            
        end do
!$OMP END PARALLEL DO
        !_______________________________________________________________________
        ! --> update mean layer thinkness at element
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, elnodes, nz, nzmin, nzmax)
        do elem=1, myDim_elem2D
            nzmin = ulevels(elem)
            nzmax = nlevels(elem)-1
            !___________________________________________________________________
            ! if there is a cavity layer thickness is not updated, its 
            ! kept fixed 
            if (nzmin > 1) cycle
            
            !___________________________________________________________________
            elnodes=elem2D_nodes(:, elem)
            do nz=nzmin, nzmax-1
                helem(nz,elem)=sum(hnode(nz,elnodes))/3.0_WP
            end do
        end do
!$OMP END PARALLEL DO
    endif
    
    !___________________________________________________________________________
!$OMP MASTER
    call exchange_elem(helem, partit)
!$OMP END MASTER
!$OMP BARRIER

end subroutine update_thickness_ale
!
!
!===============================================================================
! update thickness arrays based on the current hbar 
subroutine restart_thickness_ale(partit, mesh)
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_ARRAYS
    use g_config,only: which_ale,lzstar_lev,min_hnode
    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer :: n, nz, elem, elnodes(3), nzmax, nzmin, lcl_lzstar_lev
    integer      , dimension(:), allocatable :: idx
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    if(mype==0) then
        write(*,*) '____________________________________________________________'
        write(*,*) ' --> restart ALE layerthicknesses, depth levels and middepth levels'
        write(*,*)
    end if
        
    !___________________________________________________________________________
    ! >->->->->->->->->->->->->       z-star/zlevel       <-<-<-<-<-<-<-<-<-<-<-
    !___________________________________________________________________________
    if (trim(which_ale)=='zstar' .or. trim(which_ale)=='zlevel' ) then
        ! restart depthlevels (zbar_3d_n) and mitdpethlevels (Z_3d_n)
        ! dont forget also at restart zbar_3d_n and Z_3d_n are first initialised 
        ! and filled up in ale_init there bottom depth zbar_3d_n(nlevels_nod2d) 
        ! ist set according if there are partial cells or not 
        do n=1, myDim_nod2D+eDim_nod2D
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D(n)-1
            
            !___________________________________________________________________
            ! if there is a cavity layer thickness is not updated, its 
            ! kept fixed 
            if (nzmin > 1) cycle
            
            !___________________________________________________________________
            ! be sure that bottom layerthickness uses partial cell layer thickness
            ! in case its activated, especially when you make a restart from a non 
            ! partiall cell runs towards a simulation with partial cells
            hnode(nzmax,n) = bottom_node_thickness(n)
            
            !___________________________________________________________________
            do nz=nzmax-1,nzmin,-1
                zbar_3d_n(nz,n) =zbar_3d_n(nz+1,n) + hnode(nz,n)
                Z_3d_n(nz,n)    =zbar_3d_n(nz+1,n) + hnode(nz,n)/2.0_WP
            end do
        end do
        
        !_______________________________________________________________________
        ! restart element layer thinkness (helem) and the increment of total 
        ! fluid depth on elements (dhe)
        do elem=1, myDim_elem2D
            nzmin = ulevels(elem)
            nzmax = nlevels(elem)-1
            
            !___________________________________________________________________
            ! if there is a cavity layer thickness is not updated, its 
            ! kept fixed 
            if (nzmin > 1) cycle
            
            !___________________________________________________________________
            elnodes=elem2D_nodes(:, elem)
            do nz=nzmin,nzmax-1
                helem(nz,elem)=sum(hnode(nz,elnodes))/3.0_WP
            end do
            
            !___________________________________________________________________
            ! be sure elemental bottom thickness has partial cells in it, when 
            ! its used after restart
            helem(nzmax,elem)=bottom_elem_thickness(elem)
            
            !___________________________________________________________________
            ! for the first time steps of a restart or initialisation dhe must 
            ! be the absolute value of the elemental sea surface height, afterwards
            ! when the stiffness matrix is update dhe becomes the anomaly between 
            ! old and new elemental sea surface height (dhe(elem)=sum(hbar(elnodes)
            ! -hbar_old(elnodes))/3.0_WP in subroutine compute_hbar_ale)
            dhe(elem)=sum(hbar(elnodes))/3.0_WP
            
        end do
        
    !___________________________________________________________________________
    ! >->->->->->->->->->->->->->->       linfs      <-<-<-<-<-<-<-<-<-<-<-->->-
    !___________________________________________________________________________
    ! for the case you make a restart from zstar --> linfs, hnode that comes with
    ! the restart needs to be overwritten with the standard layer thicknesses. 
    ! Since zbar_3d_n and Z_3d_n are initialised with standard levels
    else ! --> which_ale == 'linfs'
        do n=1,myDim_nod2D+eDim_nod2D
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D(n)-1
            do nz=nzmin,nzmax-1
                hnode(nz,n)=(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
            end do      
            hnode(nzmax,n)=bottom_node_thickness(n)
        end do
        
        do elem=1,myDim_elem2D
            nzmin = ulevels(elem)
            nzmax = nlevels(elem)-1
            helem(nzmin,elem)=(zbar_e_srf(elem)-zbar(nzmin+1))
            do nz = nzmin+1, nzmax-1
                helem(nz,elem)=(zbar(nz)-zbar(nz+1))
            end do
            helem(nzmax,elem)=bottom_elem_thickness(elem)
        end do    
    endif

end subroutine restart_thickness_ale
!
!
!===============================================================================
! Stiffness matrix for the elevation
! 
! We first use local numbering and assemble the matrix
! Having completed this step we substitute global contiguous numbers.
!
! Our colind cannot be used to access local node neighbors
! This is a reminder from FESOM
!        do q=1, nghbr_nod2D(row)%nmb
!           col_pos(nghbr_nod2D(row)%addresses(q)) = q
!        enddo
!       ipos=sshstiff%rowptr(row)+col_pos(col)-sshstiff%rowptr(1)
! 
! To achive it we should use global arrays n_num and n_pos.
! Reserved for future. 
subroutine init_stiff_mat_ale(partit, mesh)
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_CONFIG
    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh
    !___________________________________________________________________________
    integer                             :: n, n1, n2, i, j, row, ed, fileID
    integer                             :: elnodes(3), el(2)
    integer                             :: npos(3), offset, nini, nend
    real(kind=WP)                       :: factor, fy(3)
    integer, allocatable                :: n_num(:), n_pos(:,:), pnza(:), rpnza(:)
    integer, allocatable                :: mapping(:)
    character*10                        :: npes_string
    character(MAX_PATH)                 :: dist_mesh_dir, file_name
    real(kind=WP)                       :: t0, t1
    integer                             :: ierror              ! MPI, return error code
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    t0=MPI_Wtime()
    if (mype==0) then
        write(*,*) '____________________________________________________________'
        write(*,*) ' --> initialise ssh operator using unperturbed ocean depth'
    endif
    
    !___________________________________________________________________________
    ! a) pre allocate stiff_mat: dimenssion,  reduced row vector ... and number 
    ! of neighbors and idices of neighbors
    ssh_stiff%dim=nod2D
    allocate(mesh%ssh_stiff%rowptr(myDim_nod2D+1), mesh%ssh_stiff%rowptr_loc(myDim_nod2D+1))
    ssh_stiff%rowptr(1)=1               ! This has to be updated as
                                        ! contiguous numbering is required
    
    allocate(n_num(myDim_nod2D+eDim_nod2D),n_pos(12,myDim_nod2D))
    n_pos=0
    
    !___________________________________________________________________________
    ! b) Neighbourhood information
    do n=1,myDim_nod2d
        ! number of neigbouring nodes to node n
        n_num(n)=1
        ! position indices of neigbouring nodes, first entry in n_pos-->n_pos(1,n)
        ! is the number of the node itself --> n_pos(1,n)=n
        n_pos(1,n)=n
    end do   
    
    ! determine neighbourhood of nodes via loop over 2d edge array (array adges ...
    ! contains only unique edges)
    do n=1, myDim_edge2D
        n1=edges(1,n)
        n2=edges(2,n)
        ! ... if(n1<=myDim_nod2D) --> because dont use extended nodes
        if (n1<=myDim_nod2D) then
            n_pos(n_num(n1)+1,n1)=n2
            n_num(n1)=n_num(n1)+1
        end if
        ! ... if(n2<=myDim_nod2D) --> because dont use extended nodes
        if (n2<=myDim_nod2D) then
            n_pos(n_num(n2)+1,n2)=n1
            n_num(n2)=n_num(n2)+1
        end if
    end do
    ! n_num...contains the number of neighbors, n_pos...contains their indices 

    !___________________________________________________________________________
    ! c) fill up reduced row vector: indice entry where sparse entrys switch to 
    ! next row, remember ssh_stiff%rowptr(1) is initialised with on and 
    ! ssh_stiff%rowptr is allocated with a length of myDim_nod2D+1
    do n=1,myDim_nod2D
        ssh_stiff%rowptr(n+1) = ssh_stiff%rowptr(n)+n_num(n)
    end do
    
    !___________________________________________________________________________
    ! d) how many nonzero entries sparse matrix has? --> Corresponds to last
    ! entry in ssh_stiff%rowptr with index myDim_nod2D+1 minus 1
    ssh_stiff%nza = ssh_stiff%rowptr(myDim_nod2D+1)-1
    
    ! allocate column as colind value array of sparse matrix, have length of nonzero 
    ! entrie of sparse matrix 
    allocate(mesh%ssh_stiff%colind(ssh_stiff%nza), mesh%ssh_stiff%colind_loc(ssh_stiff%nza))
    allocate(mesh%ssh_stiff%values(ssh_stiff%nza))
    ssh_stiff%values=0.0_WP  
    
    !___________________________________________________________________________
    ! e) fill up sparse matrix column index 
    do n=1,myDim_nod2D
        ! for every node points n, estimate start (nini) and end (nend) indices of neighbouring nodes
        ! in sparse matrix
        nini = ssh_stiff%rowptr(n) 
        nend = ssh_stiff%rowptr(n+1)- 1
        ! fill colind with local indices location from n_pos
        ssh_stiff%colind(nini:nend) = n_pos(1:n_num(n),n)
    end do
    ssh_stiff%colind_loc=ssh_stiff%colind
    ssh_stiff%rowptr_loc=ssh_stiff%rowptr

    !!! Thus far everything is in local numbering.    !!!
    !!! We will update it later when the values are   !!!
    !!! filled in                                     !!!
    
    !___________________________________________________________________________
    ! f) fill in  M/dt-alpha*theta*g*dt*\nabla(H\nabla\eta))
    !
    ! is here reinitialised to become auxilary array to switch from local to 
    ! global node indices
    n_num=0 
    
    ! 1st do secod term of lhs od equation (18) of "FESOM2 from finite element to finite volumes"
    ! stiffness part
    factor = g*dt*alpha*theta
    
    ! loop over edges
    do ed=1,myDim_edge2D   !! Attention
        
        ! el ... which two elements contribute to edge 
        el=edge_tri(:,ed)
        ! loop over two triangle elements
        do i=1,2  ! Two elements related to the edge
                ! It should be just grad on elements 
            
            if (el(i)<1) cycle ! if el(i)<1, it means its an outer boundary edge this
                            ! has only one triangle element to which it contribute
                            
            ! which three nodes span up triangle el(i)
            ! elnodes ... node indices 
            elnodes=elem2D_nodes(:,el(i))
            
            ! calc value for stiffness matrix something like H*div --> zbar is maximum depth(m)
            ! at element el(i)
            ! Attention: here corrected with bottom depth of partial cells !!!
            
            !!PS fy(1:3) = (zbar_e_bot(el(i)))* & !-> cavity
            !!PS fy(1:3) = (zbar_e_bot(el(i))-zbar(ulevels(el(i))))* &
            fy(1:3) = (zbar_e_bot(el(i))-zbar_e_srf(el(i)))* &
                      ( gradient_sca(1:3,el(i)) * edge_cross_dxdy(2*i  ,ed)   &
                       -gradient_sca(4:6,el(i)) * edge_cross_dxdy(2*i-1,ed) )
            
            if(i==2) fy=-fy
            
            ! who is node point 1 of edge ed --> connected with row index of sparse matrix
            row=edges(1,ed)
            if (row <= myDim_nod2D) then
                !n... loop over neighbouring nodes to 1st node point of edge ed 
                DO n = SSH_stiff%rowptr(row), SSH_stiff%rowptr(row+1)-1
                    ! SSH_stiff%colind(n) contains still local indices location
                    ! n is in global indexing
                    ! n_num(SSH_stiff%colind(n)) = n --> with n_num connect local index location 
                    ! SSH_stiff%colind(n) with global index n
                    n_num(SSH_stiff%colind(n)) = n
                END DO
                !npos contains global index location of local nodes that spanup elemental
                ! el(i) which contributes to edge ed
                ! npos = [1 x 3]
                npos = n_num(elnodes)
                ! fill sparse martix value array with stiffness info
                SSH_stiff%values(npos) = SSH_stiff%values(npos) + fy*factor
            endif
            
            ! who is node point 2 of edge ed
            row=edges(2,ed)
            if (row <= myDim_nod2D) then
                ! same like for row=edges(1,ed)
                DO n = SSH_stiff%rowptr(row), SSH_stiff%rowptr(row+1)-1
                    n_num(SSH_stiff%colind(n)) = n
                END DO
                npos = n_num(elnodes)
                SSH_stiff%values(npos) = SSH_stiff%values(npos) - fy*factor
            endif
        end do
    end do
    
    ! 2nd do first term of lhs od equation (18) of "FESOM2 from finite element to finite volumes"
    ! Mass matrix part
    do row=1, myDim_nod2D
        ! if cavity no time derivative for eta in case of rigid lid approximation at
        ! thee cavity-ocean interface, which means cavity-ocean interface is not allowed 
        ! to move vertically.
        if (ulevels_nod2D(row)>1) cycle
        offset = ssh_stiff%rowptr(row)
        SSH_stiff%values(offset) = SSH_stiff%values(offset)+ areasvol(ulevels_nod2D(row),row)/dt
    end do
    deallocate(n_pos,n_num)
    
    !___________________________________________________________________________
    ! g) Global contiguous numbers:
    ! Now we need to exchange between PE to know their 
    ! numbers of non-zero entries (nza):
    allocate(pnza(npes), rpnza(npes))    
    pnza(1:npes)=0
    rpnza=0
    
    ! number of nonzero entries at every CPU
    pnza(mype+1)=ssh_stiff%nza
    call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
    !collect this number from all CPUs into rpnza
    call MPI_AllREDUCE( pnza, rpnza, &
        npes, MPI_INTEGER,MPI_SUM, &
        MPI_COMM_FESOM, MPIerr)
        
    if (mype==0) then
        offset=0
    else
        ! calculate offset for all cpus mype~=0
        offset=sum(rpnza(1:mype))  ! This is offset for mype 
    end if
    
    !--> make sparse matrix row pointers from local to global by nonzero entrie 
    !    offset
    ssh_stiff%rowptr=ssh_stiff%rowptr+offset   ! pointers are global
    
    !___________________________________________________________________________
    ! replace local nza with a global one
    ssh_stiff%nza=sum(rpnza(1:npes))
    deallocate(rpnza, pnza)
    
    !___________________________________________________________________________
    ! colindices are now local to PE. We need to make them local contiguous
    ! (i) global natural: 
    do n=1,ssh_stiff%rowptr(myDim_nod2D+1)-ssh_stiff%rowptr(1)  
        ! myList_nod2D ... contains global node index of every meshpoit that belongs 
        ! to a CPU
        ! ssh_stiff%colind(n) ... contains local node index (1:myDim_nod2d)
        ! myList_nod2D(ssh_stiff%colind(n))  ... converts local index to global index
        ssh_stiff%colind(n)=myList_nod2D(ssh_stiff%colind(n))    
    end do
    
    !___________________________________________________________________________
    allocate(mapping(nod2d))
    ! 0 proc reads the data in chunks and distributes it between other procs
    write(npes_string,"(I10)") npes
    dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'
    file_name=trim(dist_mesh_dir)//'rpart.out'
    fileID=10
    if (mype==0) then
        write(*,*) '     > in stiff_mat_ale, reading ', trim(file_name)
        open(fileID, file=trim(file_name))
        ! n ... how many cpus
        read(fileID, *) n      
        ! 1st part of rpart.out: mapping(1:npes) = how many 2d node points owns every CPU
        read(fileID, *) mapping(1:npes)
        ! 2nd part of rpart.out: mapping for contigous numbering of how the 2d mesh points are
        ! locate on the CPUs: e.g node point 1, lies on CPU 2 and is there the 5th node point. 
        ! If CPU1 owns in total 5000 node points that is the mapping for the node point 1 =5005
        read(fileID, *) mapping
        close(fileID) 
    end if
    call MPI_BCast(mapping, nod2D, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
    
    !___________________________________________________________________________
    ! (ii) global PE contiguous: 
    do n=1,ssh_stiff%rowptr(myDim_nod2D+1)-ssh_stiff%rowptr(1)  
        ! convert global mesh node point numbering to global numbering of how the single 
        ! node points are contigous located on the CPUs
        ssh_stiff%colind(n)=mapping(ssh_stiff%colind(n))
    end do    
    
    !___________________________________________________________________________
    deallocate(mapping)
    t1=MPI_Wtime()
    if (mype==0) then
        write(*,*) '     took: ',t1-t0, ' sec'
        write(*,*) 
    endif

end subroutine init_stiff_mat_ale
!
!
!===============================================================================
! Update ssh stiffness matrix for a new elevation
! --> update second term on lhs of equation 18 in Danilov et al 2017
!
! 1/tau*(eta^(n+1)-eta^n) - alpha*theta*g*tau*grad*(int_H^hbar (grad*(eta^(n+1)-eta^n))*dz) = ...
!
!   1/tau*(eta^(n+1)-eta^n)                                          |
!                                                                    |--> this part is done in using the unpertubed bottom depth 
!   - alpha*theta*g*tau*grad*(int_H^0(grad*(eta^(n+1)-eta^n))*dz)    |    in the initialisation of the stiff matrix
!
!   - alpha*theta*g*tau*grad*(int_0^hbar(grad*(eta^(n+1)-eta^n))*dz) |--> the update from pertubations in the bottom depth
!                                                                         due to changes in ssh is done here 
!   = ssh_rhs                                                             in the update of the stiff matrix
!
subroutine update_stiff_mat_ale(partit, mesh)
    use g_config,only: dt
    use o_PARAM
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_ARRAYS
    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer                             :: n, i, j, k, row, ed, n2
    integer                             :: enodes(2), elnodes(3), el(2)
    integer                             :: elem, npos(3), offset, nini, nend
    real(kind=WP)                       :: factor 
    real(kind=WP)                       :: fx(3), fy(3)
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    ! update secod term of lhs od equation (18) of "FESOM2 from finite element 
    ! to finite volumes" --> stiff matrix part
    ! loop over lcal edges
    factor=g*dt*alpha*theta

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, i, j, k, row, ed, n2, enodes, elnodes, el, elem, npos, offset, nini, nend, fx, fy)
    do ed=1,myDim_edge2D   !! Attention
        ! enodes ... local node indices of nodes that edge ed
        enodes=edges(:,ed)        
        ! el ... local element indices of the two elments that contribute to edge
        ! el(1) or el(2) < 0 than edge is boundary edge
        el=edge_tri(:,ed)
        !_______________________________________________________________________
        do j=1,2 
            ! row ... local indice od edge node 1 or 2
            row=enodes(j)
            if(row>myDim_nod2D) cycle    !! Attention
            
            !___________________________________________________________________
            ! sparse indice offset for node with index row
            offset=SSH_stiff%rowptr(row)-ssh_stiff%rowptr(1)           
            !___________________________________________________________________
            do i=1, 2  ! Two elements related to the edge
                        ! It should be just grad on elements 
                ! elem ... local element index to calc grad on that element
                elem=el(i)                
                if(elem<1) cycle                
                ! elnodes ... local node indices of nodes that form element elem
                elnodes=elem2D_nodes(:,elem)
                ! we have to put it here for OMP compatibility. The MPI version might become a bit slower :(
                ! loop over number of neghbouring nodes of node-row
                do k=1, 3
                   do n=1, SSH_stiff%rowptr(row+1)-SSH_stiff%rowptr(row)
                      ! npos ... global sparse matrix indices of local mesh points elnodes
                      if (elnodes(k)==nn_pos(n, row)) then
                          npos(k)=offset+n !start with the next k
                          EXIT
                      end if
                   end do
                end do                       
                ! here update of second term on lhs of eq. 18 in Danilov etal 2017
                ! --> in the initialisation of the stiff matrix the integration went 
                !     over the unperturbed ocean depth using -zbar_e_bot 
                ! --> here this therm is now updated with the actual change in ssh 
                !     interpolated to the element dhe
                ! calculate: - alpha*theta*g*tau*grad*(int_0^hbar(grad*(eta^(n+1)-eta^n))*dz)
                fy(1:3) = -dhe(elem)* &
                         ( gradient_sca(1:3,el(i)) * edge_cross_dxdy(2*i  ,ed)   &
                          -gradient_sca(4:6,el(i)) * edge_cross_dxdy(2*i-1,ed) )
                          
                if(i==2) fy=-fy
                if(j==2) fy=-fy
                
                ! In the computation above, I've used rules from ssh_rhs (where it is 
                ! on the rhs. So the sign is changed in the expression below.
                ! npos... sparse matrix indices position of node points elnodes
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                   call omp_set_lock  (partit%plock(row)) ! it shall be sufficient to block writing into the same row of SSH_stiff
#else
!$OMP ORDERED
#endif
                   SSH_stiff%values(npos)=SSH_stiff%values(npos) + fy*factor
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                   call omp_unset_lock(partit%plock(row))
#else
!$OMP END ORDERED
#endif                
            end do ! --> do i=1,2
        end do ! --> do j=1,2 
    end do ! --> do ed=1,myDim_edge2D 
!$OMP END PARALLEL DO
!DS this check will work only on 0pe because SSH_stiff%rowptr contains global pointers
!if (mype==0) then
!do row=1, myDim_nod2D
!   nini=SSH_stiff%rowptr(row)
!   nend=SSH_stiff%rowptr(row+1)-1
!   factor=sum(SSH_stiff%values(nini:nend))/area(1,row)*dt
!   if (abs(factor-1._WP)>1.e-7) then
!      write(*,*) 'ssh_stiff mype/row/sum(vals)=', mylist_nod2D(row), factor
!   end if
!end do
!end if
!DS

end subroutine update_stiff_mat_ale
!
!
!===============================================================================
! compute new ssh_rhs which is first term on rhs of eq.18 in S. Danilov et al.: 
!"FESOM2: from finite elements to finite volumes"
!
! ssh_rhs = alpha * grad[  int_hbot^hbar(n+0.5)( u^n+deltau)dz + W(n+0.5) ]
! In the semiimplicit method: 
! ssh_rhs=-alpha*\nabla\int(U_n+U_rhs)dz-(1-alpha)*...
! see "FESOM2: from finite elements to finte volumes, S. Danilov..." eq. (11) rhs
subroutine compute_ssh_rhs_ale(dynamics, partit, mesh)
    use g_config,only: which_ALE, dt, use_cavity_fw2press
    use MOD_MESH
    use o_ARRAYS, only: water_flux
    use o_PARAM
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use g_comm_auto
    implicit none
    type(t_mesh)  , intent(inout), target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_dyn)   , intent(inout), target :: dynamics
    !___________________________________________________________________________
    integer       :: ed, el(2), enodes(2), nz, n, nzmin, nzmax
    real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2 
    real(kind=WP) :: dumc1_1, dumc1_2, dumc2_1, dumc2_2 !!PS
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:)    , pointer :: ssh_rhs, ssh_rhs_old
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV         => dynamics%uv(:,:,:)
    UV_rhs     => dynamics%uv_rhs(:,:,:)
    ssh_rhs    => dynamics%ssh_rhs(:)
    ssh_rhs_old=> dynamics%ssh_rhs_old(:)

    !___________________________________________________________________________
    ! loop over local edges
!$OMP PARALLEL DO
    do n=1, myDim_nod2D+eDim_nod2D
       ssh_rhs(n)=0.0_WP
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ed, el, enodes, n, nz, nzmin, nzmax, c1, c2, deltaX1, deltaX2, deltaY1, deltaY2, &
!$OMP                                                                                 dumc1_1, dumc1_2, dumc2_1, dumc2_2)
!$OMP DO
    do ed=1, myDim_edge2D      
        ! local indice of nodes that span up edge ed
        enodes=edges(:,ed)
        ! local index of element that contribute to edge
        el=edge_tri(:,ed)
        
        !_______________________________________________________________________
        ! calc depth integral: alpha*\nabla\int(U_n+U_rhs)dz for el(1)
        c1=0.0_WP
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,ed) 
        deltaY1=edge_cross_dxdy(2,ed)
        
        nzmin = ulevels(el(1))
        nzmax = nlevels(el(1))-1
        do nz=nzmin, nzmax
            c1=c1+alpha*((UV(2,nz,el(1))+UV_rhs(2,nz,el(1)))*deltaX1- &
                         (UV(1,nz,el(1))+UV_rhs(1,nz,el(1)))*deltaY1)*helem(nz,el(1))
        end do
        
        !_______________________________________________________________________
        ! if ed is not a boundary edge --> calc depth integral: 
        ! alpha*\nabla\int(U_n+U_rhs)dz for el(2) 
        c2=0.0_WP
        if(el(2)>0) then
            ! edge_cross_dxdy(3:4,ed)... dx,dy distance from element centroid el(2) to 
            ! center of edge --> needed to calc flux perpedicular to edge from elem el(2)
            deltaX2=edge_cross_dxdy(3,ed)
            deltaY2=edge_cross_dxdy(4,ed)
            nzmin = ulevels(el(2))
            nzmax = nlevels(el(2))-1
            do nz=nzmin, nzmax
                c2=c2-alpha*((UV(2,nz,el(2))+UV_rhs(2,nz,el(2)))*deltaX2- &
                             (UV(1,nz,el(2))+UV_rhs(1,nz,el(2)))*deltaY2)*helem(nz,el(2))
            end do
        end if
        
        !_______________________________________________________________________
        ! calc netto "flux"
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call   omp_set_lock(partit%plock(enodes(1)))
#else
!$OMP ORDERED
#endif
        ssh_rhs(enodes(1))=ssh_rhs(enodes(1))+(c1+c2)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(enodes(1)))
        call   omp_set_lock(partit%plock(enodes(2)))
#endif
        ssh_rhs(enodes(2))=ssh_rhs(enodes(2))-(c1+c2)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(enodes(2)))
#else
!$OMP END ORDERED
#endif

    end do
!$OMP END DO

    !___________________________________________________________________________
    ! take into account water flux
    ! at this point: ssh_rhs     = -alpha * nabla*int(u^n + deltau dz)
    !                ssh_rhs_old = - nabla*int(u^n dz) - water_flux*area
    !                
    ! (eta_(n+1)-eta_n)/dt = ssh_rhs - alpha*water_flux*area + (1-alpha)*ssh_rhs_old
    !                      = ssh_rhs
    !                      
    ! shown in eq (12) rhs of "FESOM2: from finite elements to finte volumes, S. Danilov..." eq. (11) rhs
    if ( .not. trim(which_ALE)=='linfs') then
!$OMP DO
        do n=1,myDim_nod2D
            nzmin = ulevels_nod2D(n)
            if (ulevels_nod2D(n)>1) then
                ! use_cavity_fw2press=.true.: adds freshwater under the cavity thereby 
                ! increasing the local pressure
                if (use_cavity_fw2press) then
                    ssh_rhs(n)=ssh_rhs(n)-alpha*water_flux(n)*areasvol(nzmin,n)
                end if
            else
                ssh_rhs(n)=ssh_rhs(n)-alpha*water_flux(n)*areasvol(nzmin,n)+(1.0_WP-alpha)*ssh_rhs_old(n)
            end if 
        end do
!$OMP END DO
    else
!$OMP DO
        do n=1,myDim_nod2D
            ssh_rhs(n)=ssh_rhs(n)+(1.0_WP-alpha)*ssh_rhs_old(n)
        end do
!$OMP END DO
    end if
!$OMP END PARALLEL
    call exchange_nod(ssh_rhs, partit)
!$OMP BARRIER
end subroutine compute_ssh_rhs_ale
!
!
!===============================================================================
! compute old ssh_rhs which is second term on rhs of eq.18...
!
! ssh_rhs_old = grad[  int_hbot^hbar(n-0.5)(u^n)dz + W(n-0.5) ]
!
! as well as new ale surface elevation hbar from eq. 13
!
! hbar(n+0.5) = hbar(n-0.5) - tau*ssh_rhs_old
!
! in S. Danilov et al.: "FESOM2: from finite elements to finite volumes"
!
! see "FESOM2: from finite elements to finte volumes, S. Danilov..." 
! hbar(n+1)-hbar(n)=tau*ssh_rhs_old 
! ssh_rhs_old=-\nabla\int(U_n)dz-water_flux*area (if free surface)
! Find new elevation hbar
subroutine compute_hbar_ale(dynamics, partit, mesh)
    use g_config,only: dt, which_ALE, use_cavity
    use MOD_MESH
    use o_ARRAYS, only: water_flux
    use o_PARAM
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use g_comm_auto
    implicit none
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh
    !___________________________________________________________________________
    integer       :: ed, el(2), enodes(2), elem, elnodes(3), n, nz, nzmin, nzmax
    real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2 
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV
    real(kind=WP), dimension(:)    , pointer :: ssh_rhs, ssh_rhs_old
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV         => dynamics%uv(:,:,:)
    ssh_rhs    => dynamics%ssh_rhs(:)
    ssh_rhs_old=> dynamics%ssh_rhs_old(:)

    !___________________________________________________________________________
    ! compute the rhs

!$OMP PARALLEL DO
    do n=1, myDim_nod2D+eDim_nod2D
       ssh_rhs_old(n)=0.0_WP
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ed, el, enodes, elem, elnodes, n, nz, nzmin, nzmax, &
!$OMP                                            c1, c2, deltaX1, deltaX2, deltaY1, deltaY2)
!$OMP DO
    do ed=1, myDim_edge2D                     
        ! local indice of nodes that span up edge ed
        enodes=edges(:,ed)
        ! local index of element that contribute to edge
        el=edge_tri(:,ed)
        
        !_______________________________________________________________________
        ! cal depth integal: \nabla\int(U_n)dz for el(1)
        c1=0.0_WP
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,ed)
        deltaY1=edge_cross_dxdy(2,ed)
        
        nzmin = ulevels(el(1))
        nzmax = nlevels(el(1))-1
        !!PS do nz=1, nlevels(el(1))-1
        do nz=nzmin, nzmax 
            c1=c1+(UV(2,nz,el(1))*deltaX1-UV(1,nz,el(1))*deltaY1)*helem(nz,el(1))
        end do
        !_______________________________________________________________________
        ! if ed is not a boundary edge --> calc depth integral: \nabla\int(U_n)dz 
        ! for el(2)
        c2=0.0_WP
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,ed)
            deltaY2=edge_cross_dxdy(4,ed)
            nzmin = ulevels(el(2))
            nzmax = nlevels(el(2))-1
            !!PS do nz=1, nlevels(el(2))-1
            do nz=nzmin, nzmax
                c2=c2-(UV(2,nz,el(2))*deltaX2-UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
            end do
        end if
        !_______________________________________________________________________
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call   omp_set_lock(partit%plock(enodes(1)))
#else
!$OMP ORDERED
#endif
        ssh_rhs_old(enodes(1))=ssh_rhs_old(enodes(1))+(c1+c2)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(enodes(1)))
        call   omp_set_lock(partit%plock(enodes(2)))
#endif
        ssh_rhs_old(enodes(2))=ssh_rhs_old(enodes(2))-(c1+c2)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(enodes(2)))
#else
!$OMP END ORDERED
#endif
    end do
!$OMP END DO
!$OMP END PARALLEL

    !___________________________________________________________________________
    ! take into account water flux
    if (.not. trim(which_ALE)=='linfs') then
!$OMP PARALLEL DO
        do n=1,myDim_nod2D
            if (ulevels_nod2D(n)>1) cycle
            ssh_rhs_old(n)=ssh_rhs_old(n)-water_flux(n)*areasvol(ulevels_nod2D(n),n)
        end do
!$OMP END PARALLEL DO
        call exchange_nod(ssh_rhs_old, partit)
!$OMP BARRIER
    end if 
    !___________________________________________________________________________
    ! update the thickness
!$OMP PARALLEL DO
    do n=1, myDim_nod2D+eDim_nod2D
       hbar_old(n)=hbar(n)
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do n=1,myDim_nod2D
        if (ulevels_nod2D(n) > 1) cycle ! --> if cavity node hbar == hbar_old
        hbar(n)=hbar_old(n)+ssh_rhs_old(n)*dt/areasvol(ulevels_nod2D(n),n)
    end do
!$OMP END PARALLEL DO
    call exchange_nod(hbar, partit)
!$OMP BARRIER
    !___________________________________________________________________________
    ! fill the array for updating the stiffness matrix
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, elnodes)
    do elem=1,myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        if (ulevels(elem) > 1) then 
            dhe(elem) = 0.0_WP
        else
            dhe(elem) = sum(hbar(elnodes)-hbar_old(elnodes))/3.0_WP
        endif 
    end do
!$OMP END PARALLEL DO
end subroutine compute_hbar_ale
!
!
!===============================================================================
! calculate vertical velocity from eq.3 in S. Danilov et al. : FESOM2: from 
! finite elements to finite volumes. 
!
!   dh_k/dt + grad(u*h)_k + (w^t - w^b) + water_flux_k=1 = 0
!
!   w^t = w^b - dh_k/dt - grad(u*h)_k - water_flux=1
!   --> do cumulativ summation from bottom to top
!
! > for linfs : dh_k/dt = 0
! > for zlevel: dh_k/dt_k=1 != 0
! > for zstar : dh_k/dt_k=1...kbot-1 != 0
!
subroutine vert_vel_ale(dynamics, partit, mesh)
    use g_config,only: dt, which_ALE, min_hnode, lzstar_lev, flag_warn_cflz
    use MOD_MESH
    use o_ARRAYS, only: water_flux
    use o_PARAM
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use g_comm_auto
    use io_RESTART !!PS
    use g_forcing_arrays !!PS
    use compute_Wvel_split_interface
    use compute_CFLz_interface
    implicit none
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh
    !___________________________________________________________________________
    integer       :: el(2), enodes(2), n, nz, ed, nzmin, nzmax, uln1, uln2, nln1, nln2
    real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, dd, dd1, dddt, cflmax
    ! still to be understood but if you allocate these arrays statically the results will be different:
    real(kind=WP) :: c1(mesh%nl-1), c2(mesh%nl-1)
    ! --> zlevel with local zstar
    real(kind=WP) :: dhbar_total, dhbar_rest, distrib_dhbar_int
    real(kind=WP), dimension(:), allocatable :: max_dhbar2distr, cumsum_maxdhbar, distrib_dhbar
    integer      , dimension(:), allocatable :: idx
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, fer_UV
    real(kind=WP), dimension(:,:)  , pointer :: Wvel, Wvel_e, Wvel_i, CFL_z, fer_Wvel
    real(kind=WP), dimension(:)    , pointer :: ssh_rhs, ssh_rhs_old
    real(kind=WP), dimension(:)    , pointer :: eta_n, d_eta
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV          => dynamics%uv(:,:,:)
    Wvel        => dynamics%w(:,:)
    Wvel_e      => dynamics%w_e(:,:)
    Wvel_i      => dynamics%w_i(:,:)
    CFL_z       => dynamics%cfl_z(:,:)
    eta_n       => dynamics%eta_n(:)
    if (Fer_GM) then
        fer_UV  => dynamics%fer_uv(:,:,:)
        fer_Wvel=> dynamics%fer_w(:,:)
    end if
    
    if ( .not. dynamics%use_ssh_se_subcycl) then
        d_eta       => dynamics%d_eta(:)
        ssh_rhs     => dynamics%ssh_rhs(:)
        ssh_rhs_old => dynamics%ssh_rhs_old(:)
    else
    
    end if 
    
    !___________________________________________________________________________
    ! Contributions from levels in divergence
!$OMP PARALLEL DO
    DO n=1, myDim_nod2D+eDim_nod2D
       Wvel(:, n)=0.0_WP
       if (Fer_GM) then
          fer_Wvel(:, n)=0.0_WP
       end if
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ed, enodes, el, deltaX1, deltaY1, nz, nzmin, nzmax, deltaX2, deltaY2, c1, c2)
!$OMP DO
    do ed=1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        enodes=edges(:,ed)   
        
        ! local index of element that contribute to edge
        el=edge_tri(:,ed)
        
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,ed)
        deltaY1=edge_cross_dxdy(2,ed)
        
        !_______________________________________________________________________
        ! calc div(u_vec*h) for every layer 
        ! do it with gauss-law: int( div(u_vec)*dV) = int( u_vec * n_vec * dS )
        nzmin = ulevels(el(1))
        nzmax = nlevels(el(1))-1
! we introduced c1 & c2 as arrays here to avoid deadlocks when in OpenMP mode
        do nz = nzmax, nzmin, -1
            ! --> h * u_vec * n_vec
            ! --> e_vec = (dx,dy), n_vec = (-dy,dx);
            ! --> h * u*(-dy) + v*dx
            c1(nz)=( UV(2,nz,el(1))*deltaX1 - UV(1,nz,el(1))*deltaY1 )*helem(nz,el(1))
            ! inflow(outflow) "flux" to control volume of node enodes1
            ! is equal to outflow(inflow) "flux" to control volume of node enodes2
            if (Fer_GM) then
                c2(nz)=(fer_UV(2,nz,el(1))*deltaX1- fer_UV(1,nz,el(1))*deltaY1)*helem(nz,el(1))
            end if
        end do
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_set_lock  (partit%plock(enodes(1)))
#else
!$OMP ORDERED
#endif
        Wvel       (nzmin:nzmax, enodes(1))= Wvel    (nzmin:nzmax, enodes(1))+c1(nzmin:nzmax)
        if (Fer_GM) then
           fer_Wvel(nzmin:nzmax, enodes(1))= fer_Wvel(nzmin:nzmax, enodes(1))+c2(nzmin:nzmax)
        end if
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(enodes(1)))
        call omp_set_lock  (partit%plock(enodes(2)))
#endif
        Wvel       (nzmin:nzmax, enodes(2))= Wvel    (nzmin:nzmax, enodes(2))-c1(nzmin:nzmax)
        if (Fer_GM) then
           fer_Wvel(nzmin:nzmax, enodes(2))= fer_Wvel(nzmin:nzmax, enodes(2))-c2(nzmin:nzmax)
        end if
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(enodes(2)))
#else
!$OMP END ORDERED
#endif
        !_______________________________________________________________________
        ! if ed is not a boundary edge --> calc div(u_vec*h) for every layer
        ! for el(2)
        c1 = 0.0_WP
        c2 = 0.0_WP
        if(el(2)>0)then
            deltaX2=edge_cross_dxdy(3,ed)
            deltaY2=edge_cross_dxdy(4,ed)
            nzmin = ulevels(el(2))
            nzmax = nlevels(el(2))-1   
            do nz = nzmax, nzmin, -1
                c1(nz)=-(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
                if (Fer_GM) then
                    c2(nz)=-(fer_UV(2,nz,el(2))*deltaX2-fer_UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
                end if
            end do
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_set_lock  (partit%plock(enodes(1)))
#else
!$OMP ORDERED
#endif
            Wvel       (nzmin:nzmax, enodes(1))= Wvel    (nzmin:nzmax, enodes(1))+c1(nzmin:nzmax)
            if (Fer_GM) then
               fer_Wvel(nzmin:nzmax, enodes(1))= fer_Wvel(nzmin:nzmax, enodes(1))+c2(nzmin:nzmax)
            end if
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(enodes(1)))
            call omp_set_lock  (partit%plock(enodes(2)))
#endif
            Wvel       (nzmin:nzmax, enodes(2))= Wvel    (nzmin:nzmax, enodes(2))-c1(nzmin:nzmax)
            if (Fer_GM) then
               fer_Wvel(nzmin:nzmax, enodes(2))= fer_Wvel(nzmin:nzmax, enodes(2))-c2(nzmin:nzmax)
            end if
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(enodes(2)))
#else
!$OMP END ORDERED
#endif
        end if !~-> if(el(2)>0)then       
    end do ! --> do ed=1, myDim_edge2D
!$OMP END DO
!$OMP END PARALLEL
    ! |
    ! |
    ! +--> until here Wvel contains the thickness divergence div(u)
    !___________________________________________________________________________
    ! cumulative summation of div(u_vec*h) vertically
    ! W_k = W_k+1 - div(h_k*u_k)
    ! W_k ... vertical flux trough 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax)
    do n=1, myDim_nod2D
        nzmin = ulevels_nod2D(n)
        nzmax = nlevels_nod2d(n)-1
        
        do nz=nzmax,nzmin,-1
            Wvel(nz,n)=Wvel(nz,n)+Wvel(nz+1,n)
            if (Fer_GM) then 
                fer_Wvel(nz,n)=fer_Wvel(nz,n)+fer_Wvel(nz+1,n)
            end if
        end do
        
        !_______________________________________________________________________
        if ( any(Wvel(nzmin:nzmax, n)/=Wvel(nzmin:nzmax, n)) ) then
            write(*,*) ' --> subroutine vert_vel_ale --> found Nan in Wvel after cumulativ summation(...)'
            write(*,*) ' mype =', mype
            write(*,*) ' node =', n
            write(*,*) ' Wvel(nzmin:nzmax, n)=', Wvel(nzmin:nzmax, n) 
        end if
    end do
!$OMP END PARALLEL DO
    !___________________________________________________________________________
    ! divide with depth dependent cell area to convert from Vertical flux to 
    ! physical vertical velocities in units m/s
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax)
    do n=1, myDim_nod2D
        nzmin = ulevels_nod2D(n)
        nzmax = nlevels_nod2d(n)-1
        
        do nz=nzmin,nzmax
            Wvel(nz,n)=Wvel(nz,n)/area(nz,n)
            if (Fer_GM) then 
                fer_Wvel(nz,n)=fer_Wvel(nz,n)/area(nz,n)
            end if
            
        end do
        
        !_______________________________________________________________________
        if ( any(Wvel(nzmin:nzmax, n)/=Wvel(nzmin:nzmax, n)) ) then
            write(*,*) ' --> subroutine vert_vel_ale --> found Nan in Wvel after divide with area'
            write(*,*) ' mype =', mype
            write(*,*) ' node =', n
            write(*,*) ' Wvel(nzmin:nzmax, n)=', Wvel(nzmin:nzmax, n) 
        end if
    end do
!$OMP END PARALLEL DO
    ! |
    ! |--> (A) linear free surface: dh/dt=0 ; W_t-W_b = -div(hu)
    ! |
    ! |        Full free surface cases:
    ! |        ------------------------
    ! |--> (B) ZLEVEL:  k=1 --> W_t-W_b = -div(hu) - dh/dt - Wflux_k1
    ! |                 k>1 --> W_t-W_b = -div(hu)
    ! |
    ! |--> (C) ZSTAR:  W_t-W_b = -div(hu) - dh/dt - Wflux_k1
    !                                         |
    !                                         |--> (dh/dt)_k = 1/H*dh/dt
    
    !___________________________________________________________________________
    ! Correct for free surface (zlevel and zstar)
    if(trim(which_ALE)=='zlevel') then
        !_______________________________________________________________________
        ! Update the upper level
        ! water_flux is positive if out of the ocean
        ! Wvel(1,n) should be 0 up to machine precision only then volume is 
        ! conserved --> at this place should be checked.
        
        !_______________________________________________________________________
        ! idx is only needed for local star case to estimate over how much 
        ! depth layers change in ssh needs to be distributed
        !!PS allocate(max_dhbar2distr(nl-1),distrib_dhbar(nl-1),idx(nl-1),cumsum_maxdhbar(nl-1))
        !!PS idx = (/(nz,nz=1,nl-1,1)/)
       allocate(max_dhbar2distr(lzstar_lev), distrib_dhbar(lzstar_lev), idx(lzstar_lev), cumsum_maxdhbar(lzstar_lev))
       idx = (/(nz, nz=1, lzstar_lev, 1)/)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax, dhbar_total, max_dhbar2distr, cumsum_maxdhbar, &
!$OMP                                                    distrib_dhbar, dhbar_rest, distrib_dhbar_int)
!$OMP DO
        do n=1, myDim_nod2D
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D_min(n)-1
            
            !_______________________________________________________________________
            ! compute new surface vertical velocity and layer thickness only when 
            ! there is no cavity --> when there is cavity treat like linfs
            if (nzmin==1) then 
                !___________________________________________________________________
                ! total ssh change to distribute
                dhbar_total = hbar(n)-hbar_old(n)
                
                !___________________________________________________________________
                ! if new surface layerthickness at node n is smaller than the initial 
                ! layerthickness*min_hnode than go from zlevel to local zstar approach
                ! over the first lzstar_lev layers.
                ! --> otherwise it can happen, especially with floating ice, that 
                !     layerthickness becomes to small or even negativ and model 
                !     blows up
                !!PS if (dhbar_total<0.0_WP .and. hnode(1,n)+dhbar_total<=(zbar(1)-zbar(2))*min_hnode ) then
                if (dhbar_total < 0.0_WP .and. hnode(nzmin,n)+dhbar_total<=(zbar(nzmin)-zbar(nzmin+1))*min_hnode ) then 
                    ! --> do local zstar case 
                    !_______________________________________________________________
                    ! max_dhbar2distr ... how much negative ssh change can be maximal 
                    ! distributed per layer (must be negativ, if positive or ==0 
                    ! layer reached already minimum layerthickness)
                    max_dhbar2distr = 0.0_WP
                    !max_dhbar2distr = (zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1))*min_hnode - hnode(1:lzstar_lev,n);
                    max_dhbar2distr = (zbar(nzmin:nzmin+lzstar_lev-1)-zbar(nzmin+1:nzmin+lzstar_lev-1+1))*min_hnode - hnode(nzmin:nzmin+lzstar_lev-1,n);
                    where (max_dhbar2distr>=0.0_WP) max_dhbar2distr=0.0_WP
                    
                    !_______________________________________________________________
                    ! if vertical CFL criteria at a certain node is at its limit 
                    ! don't take away further layer thickness --> take it than better 
                    ! from a deeper layer
                    !!PS where (CFL_z(1:lzstar_lev,n)>=0.95_WP) max_dhbar2distr=0.0_WP
                    where (CFL_z(nzmin:nzmin+lzstar_lev-1,n)>=0.95_WP) max_dhbar2distr=0.0_WP
                    
                    !_______________________________________________________________
                    ! try to limitate over how much layers i realy need to distribute
                    ! the change in ssh, so that the next loops run only over the 
                    ! nesseccary levels and not over all lzstar_lev levels
                    ! --> do this with cumulativ summation of maximum dhbar that can 
                    !     be distributed per layer. Than search index where this
                    !     cumulativ sum is larger than dhbar_total
                    cumsum_maxdhbar(1)            =  max_dhbar2distr(1)
                    cumsum_maxdhbar(2:lzstar_lev) = (/(max_dhbar2distr(nz)+max_dhbar2distr(nz-1),nz=2,lzstar_lev,1)/)
                    nz = minval(pack(idx,cumsum_maxdhbar<dhbar_total))
                    nz = min(nz,lzstar_lev)
                    
                    !_______________________________________________________________
                    ! calc array for distribution of ssh change over layers
                    distrib_dhbar = 0.0_WP
                    dhbar_rest    = dhbar_total
                    
                    ! nlevels_nod2D_min(n)-1 ...would be hnode of partial bottom 
                    ! cell but this one is not allowed to change so go until 
                    ! nlevels_nod2D_min(n)-2
                    !!PS nzmax = min(nz,nlevels_nod2D_min(n)-2)
                    nzmax = min(nz,nzmax-1)
                    do nz=1,nzmax
                        distrib_dhbar(nz) = max(dhbar_rest,max_dhbar2distr(nz))    
                        dhbar_rest        = dhbar_rest - distrib_dhbar(nz)
                        dhbar_rest        = min(0.0_WP,dhbar_rest)
                    end do
                    
                    !_______________________________________________________________
                    if ( abs(sum(distrib_dhbar)-dhbar_total)>1.0e-10 ) then
                        write(*,*) " --> problem <-- with conservation of dhbar distribution over depth"
                        write(*,*) "                 there are not enough layers to distribute "
                        write(*,*) "                 all change in ssh "
                        write(*,*) "                  > mype        =",mype
                        write(*,*) "                  > node        =",n
                        write(*,*) "                  > mstep       =",mstep
                        write(*,*) "                  > dhbar_total =",dhbar_total
                        write(*,*) "                  > dhbar_rest  =",dhbar_rest
                        write(*,*) "                  > lzstar_lev  =",lzstar_lev
                        write(*,*) "                  > nzmax       =",nzmax
                        write(*,*) "                  > max_dhbar2distr=",max_dhbar2distr
                        write(*,*) "                  > hnode_min=",(zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1))*min_hnode
                        write(*,*) "                  > hnode_now=",hnode(1:lzstar_lev,n)
                        
                    end if 
                    
                    !_______________________________________________________________
                    distrib_dhbar_int = 0.0_WP
                    do nz=nzmax,1,-1
                        !___________________________________________________________
                        ! --> integrate ssh distribution from down to up
                        distrib_dhbar_int = distrib_dhbar_int + distrib_dhbar(nz)
                        
                        !___________________________________________________________
                        ! --> distribute change in ssh over layers in hnode and Wvel
                        !!PS Wvel(nz,n)        = Wvel(nz,n) - distrib_dhbar_int/dt
                        !!PS hnode_new(nz,n)   = hnode(nz,n)+ distrib_dhbar(nz)
                        Wvel(nzmin+nz-1,n)        = Wvel(nzmin+nz-1,n) - distrib_dhbar_int/dt
                        hnode_new(nzmin+nz-1,n)   = hnode(nzmin+nz-1,n)+ distrib_dhbar(nz)
                    end do
                    
                !___________________________________________________________________
                ! in case local zstar was applied must allow the mesh in case of 
                ! positive ssh change to return to the normal zlevel case, that means
                ! to first "refill" the subsurface layerthickness and with the rest 
                ! than the surface layerthickness
                !!PS elseif (dhbar_total>0.0_WP .and. & 
                !!PS         any(hnode(2:lzstar_lev,n)/=(zbar(2:lzstar_lev)-zbar(3:lzstar_lev+1))) &
                !!PS         ) then
                elseif (dhbar_total>0.0_WP .and. & 
                        any(hnode(nzmin+1:nzmin+lzstar_lev-1,n)/=(zbar(nzmin+1:nzmin+lzstar_lev-1)-zbar(nzmin+2:nzmin+lzstar_lev-1+1))) &
                        ) then
                    ! --> do return to zlevel
                    !_______________________________________________________________
                    ! max_dhbar2distr ... how much positive ssh change must be 
                    ! distributed in the subsurface layers to be able to return to 
                    ! the init layerthickness
                    max_dhbar2distr   = 0.0_WP
                    !!PS max_dhbar2distr   = (zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1)) - hnode(1:lzstar_lev,n);
                    max_dhbar2distr   = (zbar(nzmin:nzmin+lzstar_lev-1)-zbar(nzmin+1:nzmin+lzstar_lev-1+1)) - hnode(nzmin:nzmin+lzstar_lev-1,n);
                    ! there is no limitation in the surface layer how much positive 
                    ! ssh change can be put there (1000.0_WP is just an arbitrary 
                    ! high value that should no be reached by dhbar_total)
                    max_dhbar2distr(1)= 1000.0_WP
                    
                    !_______________________________________________________________
                    ! try to limitate over how much layers i realy need to distribute
                    ! the change in ssh, so that the next loops run only over the 
                    ! nesseccary levels and not over all lzstar_lev levels
                    !!PS nz = maxval(pack(idx,hnode(1:lzstar_lev,n)/=(zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1))))
                    nz = maxval(pack(idx,hnode(nzmin:nzmin+lzstar_lev-1,n)/=(zbar(nzmin:nzmin+lzstar_lev-1)-zbar(nzmin+1:nzmin+lzstar_lev-1+1))))
                    
                    ! nlevels_nod2D_min(n)-1 ...would be hnode of partial bottom 
                    ! cell but this one is not allowed to change so go until 
                    ! nlevels_nod2D_min(n)-2
                    !!PS nzmax = min(nz,nlevels_nod2D_min(n)-2)
                    nzmax = min(nz,nzmax-1)
                    
                    !_______________________________________________________________
                    ! calc array for distribution of ssh change over layers
                    dhbar_rest        = dhbar_total
                    distrib_dhbar     = 0.0_WP
                    distrib_dhbar_int = 0.0_WP
                    !!PS do nz=nzmax,1,-1
                    do nz=nzmax,1,-1
                        !___________________________________________________________
                        distrib_dhbar(nz) = min(dhbar_rest,max_dhbar2distr(nz))    
                        dhbar_rest        = dhbar_rest - distrib_dhbar(nz)
                        dhbar_rest        = max(0.0_WP,dhbar_rest)
                        
                        !___________________________________________________________
                        ! --> integrate ssh distribution from down to up
                        distrib_dhbar_int = distrib_dhbar_int + distrib_dhbar(nz)
                        
                        !___________________________________________________________
                        ! --> distribute change in ssh over layers in hnode and Wvel
                        Wvel(     nzmin+nz-1,n) = Wvel( nzmin+nz-1,n) - distrib_dhbar_int/dt
                        hnode_new(nzmin+nz-1,n) = hnode(nzmin+nz-1,n) + distrib_dhbar(nz)
                        
                    end do
                    
                !___________________________________________________________________
                else
                    ! --> do normal zlevel case
                    ! only distribute change in ssh for Wvel and hnode_new into the 
                    ! surface layer
                    !!PS Wvel(1,n)      = Wvel(1,n) -dhbar_total/dt
                    !!PS hnode_new(1,n) = hnode(1,n)+dhbar_total
                    Wvel(nzmin,n)      = Wvel(nzmin,n)  - dhbar_total/dt
                    hnode_new(nzmin,n) = hnode(nzmin,n) + dhbar_total
                    
                end if ! --> if (dhbar_total<0 .and. hnode(1,n)+dhbar_total<=... ) then 
                
                !___________________________________________________________________
                ! Add surface fresh water flux as upper boundary condition for continutity
                !!PS Wvel(1,n) = Wvel(1,n)-water_flux(n)
                Wvel(nzmin,n) = Wvel(nzmin,n)-water_flux(n)
            
            end if ! --> if (nzmin==1) then 
            
        end do ! --> do n=1, myDim_nod2D
!$OMP END DO
!$OMP END PARALLEL
        !_______________________________________________________________________
        deallocate(max_dhbar2distr,distrib_dhbar,idx,cumsum_maxdhbar)
        
    !___________________________________________________________________________
    elseif (trim(which_ALE)=='zstar') then
        ! distribute total change in ssh (hbar(n)-hbar_old(n)) over all layers 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax, dd, dd1, dddt)
        do n=1, myDim_nod2D
            nzmin = ulevels_nod2D(n)
            !!PS nzmin = ulevels_nod2D_max(n)
            nzmax = nlevels_nod2D_min(n)-1
            
            !_______________________________________________________________________
            ! compute new surface vertical velocity and layer thickness only when 
            ! there is no cavity --> when there is cavity treat like linfs
            if (nzmin==1) then 
                
                !___________________________________________________________________
                ! --> be careful Sergey suggest in his paper to use the unperturbed
                !     ocean levels NOT the actual one !!! but spoke with Sergey its not 
                !     so important which to use as long as it is consistent and 
                !     volume is conserved
                !!PS dd1=zbar_3d_n(nlevels_nod2D_min(n)-1,n)
                dd1=zbar_3d_n(nzmax,n)
                
                ! This is the depth the stretching is applied (area(nz,n)=area(1,n))
                !!ps dd=zbar_3d_n(1,n)-dd1    
                dd=zbar_3d_n(nzmin,n)-dd1  
                
                ! how much of (hbar(n)-hbar_old(n)) is distributed into each layer
                ! 1/H*dhbar
                dd=(hbar(n)-hbar_old(n))/dd
                
                !___________________________________________________________________
                ! 1/H*dhbar/dt
                dddt=dd/dt
                
                !___________________________________________________________________
                !!PS do nz=1,nlevels_nod2D_min(n)-2
                do nz=nzmin, nzmax-1
                    ! why  *(zbar(nz)-dd1) ??? 
                    ! because here Wvel_k = SUM_k:kmax(div(h_k*v_k))/V_k
                    ! but Wvel_k = Wvel_k+1 - div(h_k*v_k) - h⁰_k/H*dhbar/dt
                    !                |--> Wvel_k+1 = Wvel_k+2 - div(h_k+1*v_k+1) - h⁰_k+1/H*dhbar/dt
                    !                                  |--> Wvel_k+2 = Wvel_k+3 - div(h_k+2*v_k+2) - h⁰_k+2/H*dhbar/dt
                    !
                    ! Wvel_k             = SUM_i=k:kmax(div(h_i*v_i)) + 1/H*dhbar/dt*SUM_i=k:kmax(h⁰_k)
                    ! SUM_i=k:kmax(h⁰_k) = (zbar(nz)-dd1)
                    ! --> this strange term zbar_3d_n(nz,n)-dd1)*dddt --> comes from 
                    !     the vertical integration bottom to top of Wvel
                    Wvel(nz,n)    =Wvel(nz,n) -(zbar_3d_n(nz,n)-dd1)*dddt
                    
                    hnode_new(nz,n)=hnode(nz,n)+(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))*dd
                end do
                
                !___________________________________________________________________
                ! Add surface fresh water flux as upper boundary condition for 
                ! continutity
                Wvel(nzmin,n)=Wvel(nzmin,n)-water_flux(n) 
                
                !_______________________________________________________________________
                if ( any(Wvel(nzmin:nzmax, n)/=Wvel(nzmin:nzmax, n)) ) then
                    write(*,*) ' --> subroutine vert_vel_ale --> found Nan in Wvel after ssh is distributed'
                    write(*,*) ' mype =', mype
                    write(*,*) ' node =', n
                    write(*,*) ' Wvel(nzmin:nzmax, n)=', Wvel(nzmin:nzmax, n) 
                end if
                
            endif ! --> if (nzmin==1) then
            
        end do ! --> do n=1, myDim_nod2D
!$OMP END PARALLEL DO
        ! The implementation here is a bit strange, but this is to avoid 
        ! unnecessary multiplications and divisions by area. We use the fact 
        ! that we apply stretching only over the part of the column
        ! where area(nz,n)=area(1,n)
    endif ! --> if(trim(which_ALE)=='....') then

!$OMP PARALLEL DO
        do n=1, myDim_nod2D+eDim_nod2D
            if (any( hnode_new(:,n) < 0.0_WP)) then
                write(*,*) ' --> fatal problem <--: layerthickness of a layer became smaller than zero'
                write(*,*) "          mype = ", mype
                write(*,*) "         mstep = ", mstep
                write(*,*) "          node = ", n
                write(*,*) 'glon,glat      = ',geo_coord_nod2D(:,n)/rad
                write(*,*)
                write(*,*) 'water_flux     = ', water_flux(n)
                write(*,*)
                write(*,*) "eta_n          = ", eta_n(n)
                if ( .not. dynamics%use_ssh_se_subcycl) then
                    write(*,*) "d_eta          = ", d_eta(n)
                    write(*,*) "ssh_rhs        = ", ssh_rhs(n)
                    write(*,*) "ssh_rhs_old    = ", ssh_rhs_old(n)
                end if 
                write(*,*)
                write(*,*) "zbar_3d_n(1,n) = ", zbar_3d_n(1,n)
                write(*,*) "dd1            = ", dd1
                write(*,*) "nlevels_nod2D_min(n)-1 = ",nlevels_nod2D_min(n)-1
                write(*,*)
                write(*,*) "dhbar/H        = ", dd
                write(*,*) "1/H*dhbar/dt   = ", dddt
                write(*,*) "hbar(n)        = ", hbar(n)
                write(*,*) "hbar_old(n)    = ", hbar_old(n)
                write(*,*) "hbar(n)-hbar_old(n) = ", hbar(n)-hbar_old(n)
                write(*,*)
                write(*,*) "hnode_new(:,n) = ", hnode_new(:,n)
                write(*,*) "hnode(:,n)     = ", hnode(:,n)
                write(*,*)
                write(*,*) "zbar_3d_n(:,n) = ", zbar_3d_n(:,n)
                write(*,*) "Z_3d_n(:,n)    = ", Z_3d_n(:,n)
                write(*,*)
                write(*,*) "zbar_n_bot(n)  = ", zbar_n_bot(n)
                write(*,*) "bottom_node_thickness(n) = ", bottom_node_thickness(n)
                write(*,*)
            end if 
        end do
!$OMP END PARALLEL DO
    !___________________________________________________________________________
    call exchange_nod(Wvel, partit)
    call exchange_nod(hnode_new, partit)   ! Or extend cycles above  
    if (Fer_GM) call exchange_nod(fer_Wvel, partit)
!$OMP BARRIER    

    !___________________________________________________________________________
    ! compute vertical CFL_z criteria
    call compute_CFLz(dynamics, partit, mesh)
    
    !___________________________________________________________________________
    ! compute implicite explicite splitting of vetical velocity Wvel according 
    ! to CFL_z criteria
    call compute_Wvel_split(dynamics, partit, mesh)
    
end subroutine vert_vel_ale

!
!
!_______________________________________________________________________________
! calculate vertical velocity from eq.3 in S. Danilov et al. : FESOM2: from 
! finite elements to finite volumes. 
!
!   dh_k/dt + grad(u*h)_k + (w^t - w^b) + water_flux_k=1 = 0
!
!   w^t = w^b - dh_k/dt - grad(u*h)_k - water_flux=1
!   --> do cumulativ summation from bottom to top
subroutine compute_vert_vel_transpv(dynamics, partit, mesh)
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE MOD_DYN
    use o_ARRAYS, only: water_flux
    use g_config, only: dt, which_ale
    use g_comm_auto
    use compute_Wvel_split_interface
    use compute_CFLz_interface
    implicit none
    !___________________________________________________________________________
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    integer                               :: node, elem, nz, ed, nzmin, nzmax, ednodes(2), edelem(2) 
    real(kind=WP)                         :: hh_inv, deltaX1, deltaX2, deltaY1, deltaY2 
    real(kind=WP)                         :: c1(mesh%nl-1), c2(mesh%nl-1)
    
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UVh, fer_UV
    real(kind=WP), dimension(:,:)  , pointer :: Wvel, Wvel_e, fer_Wvel
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UVh         => dynamics%se_uvh(:,:,:)
    Wvel        => dynamics%w(:,:)
    Wvel_e      => dynamics%w_e(:,:)
    if (Fer_GM) then
        fer_UV  => dynamics%fer_uv(:,:,:)
        fer_Wvel=> dynamics%fer_w(:,:)
    end if
    
    !___________________________________________________________________________
!$OMP PARALLEL DO
    do node=1, myDim_nod2D+eDim_nod2D
        Wvel(:, node)=0.0_WP
        if (Fer_GM) then
            fer_Wvel(:, node)=0.0_WP
        end if
    end do ! --> do node=1, myDim_nod2D+eDim_nod2D
!$OMP END PARALLEL DO

    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ed, ednodes, edelem, nz, nzmin, nzmax, & 
!$OMP                                  deltaX1, deltaY1, deltaX2, deltaY2, c1, c2)
!$OMP DO
    do ed=1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        ednodes=edges(:,ed)   
        
        ! local index of element that contribute to edge
        edelem=edge_tri(:,ed)
        
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid edelem(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem edelem(1)
        deltaX1=edge_cross_dxdy(1,ed)
        deltaY1=edge_cross_dxdy(2,ed)
        
        !_______________________________________________________________________
        ! calc div(u_vec*h) for every layer 
        ! do it with gauss-law: int( div(u_vec)*dV) = int( u_vec * n_vec * dS )
        nzmin = ulevels(edelem(1))
        nzmax = nlevels(edelem(1))-1
        ! we introduced c1 & c2 as arrays here to avoid deadlocks when in OpenMP mode
        do nz = nzmax, nzmin, -1
            ! --> h * u_vec * n_vec
            ! --> e_vec = (dx,dy), n_vec = (-dy,dx);
            ! --> h * u*(-dy) + v*dx
            c1(nz)=( UVh(2, nz, edelem(1))*deltaX1 - UVh(1, nz, edelem(1))*deltaY1 )
            ! inflow(outflow) "flux" to control volume of node enodes1
            ! is equal to outflow(inflow) "flux" to control volume of node enodes2
        end do ! --> do nz=nzmax,nzmin,-1
        if (Fer_GM) then
            do nz = nzmax, nzmin, -1
                c2(nz)=(fer_UV(2, nz, edelem(1))*deltaX1 - fer_UV(1, nz, edelem(1))*deltaY1)*helem(nz, edelem(1))
            end do
        end if 
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_set_lock  (partit%plock(ednodes(1)))
#else
!$OMP ORDERED
#endif        
        Wvel(nzmin:nzmax, ednodes(1))= Wvel(nzmin:nzmax, ednodes(1))+c1(nzmin:nzmax)
        if (Fer_GM) then
            fer_Wvel(nzmin:nzmax, ednodes(1))= fer_Wvel(nzmin:nzmax, ednodes(1))+c2(nzmin:nzmax)
        end if
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(ednodes(1)))
        call omp_set_lock  (partit%plock(ednodes(2)))
#endif        
        Wvel(nzmin:nzmax, ednodes(2))= Wvel(nzmin:nzmax, ednodes(2))-c1(nzmin:nzmax)
        if (Fer_GM) then
            fer_Wvel(nzmin:nzmax, ednodes(2))= fer_Wvel(nzmin:nzmax, ednodes(2))-c2(nzmin:nzmax)
        end if
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(ednodes(2)))
#else
!$OMP END ORDERED
#endif
        
        !_______________________________________________________________________
        ! if ed is not a boundary edge --> calc div(u_vec*h) for every layer
        ! for edelem(2)
        c1 = 0.0_WP
        c2 = 0.0_WP
        if(edelem(2)>0)then
            deltaX2=edge_cross_dxdy(3,ed)
            deltaY2=edge_cross_dxdy(4,ed)
            nzmin = ulevels(edelem(2))
            nzmax = nlevels(edelem(2))-1   
            do nz = nzmax, nzmin, -1
                c1(nz)=-(UVh(2, nz, edelem(2))*deltaX2 - UVh(1, nz, edelem(2))*deltaY2)
            end do ! --> do nz=nzmax,nzmin,-1
            if (Fer_GM) then
                do nz = nzmax, nzmin, -1
                    c2(nz)=-(fer_UV(2, nz, edelem(2))*deltaX2-fer_UV(1, nz, edelem(2))*deltaY2)*helem(nz, edelem(2))
                end do ! --> do nz=nzmax,nzmin,-1
            end if 
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_set_lock  (partit%plock(ednodes(1)))
#else
!$OMP ORDERED
#endif                
            Wvel(nzmin:nzmax, ednodes(1))= Wvel(nzmin:nzmax, ednodes(1))+c1(nzmin:nzmax)
            if (Fer_GM) then
                fer_Wvel(nzmin:nzmax, ednodes(1))= fer_Wvel(nzmin:nzmax, ednodes(1))+c2(nzmin:nzmax)
            end if 
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(ednodes(1)))
            call omp_set_lock  (partit%plock(ednodes(2)))
#endif            
            Wvel(nzmin:nzmax, ednodes(2))= Wvel(nzmin:nzmax, ednodes(2))-c1(nzmin:nzmax)
            if (Fer_GM) then
                fer_Wvel(nzmin:nzmax, ednodes(2))= fer_Wvel(nzmin:nzmax, ednodes(2))-c2(nzmin:nzmax)
            end if
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(ednodes(2)))
#else
!$OMP END ORDERED
#endif            
        end if !--> if(edelem(2)>0)then
        
    end do ! --> do ed=1, myDim_edge2D
!$OMP END DO    
!$OMP END PARALLEL

    !___________________________________________________________________________
    ! add the contribution from -dh/dt * area
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax)   
    do node=1, myDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2d(node)-1
        do nz=nzmax,nzmin,-1
            Wvel(nz, node)=Wvel(nz, node)-(hnode_new(nz, node)-hnode(nz, node))*area(nz, node)/dt 
        end do ! --> do nz=nzmax,nzmin,-1
    end do ! --> do node=1, myDim_nod2D
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    ! Sum up to get W*area
    ! cumulative summation of div(u_vec*h) vertically
    ! W_k = W_k+1 - div(h_k*u_k)
    ! W_k ... vertical flux troughdo node=1, myDim_nod2D
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax)       
    do node=1, myDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2d(node)-1
        do nz=nzmax,nzmin,-1
            Wvel(nz, node)=Wvel(nz, node)+Wvel(nz+1, node)  
            if (Fer_GM) then 
                fer_Wvel(nz, node)=fer_Wvel(nz, node)+fer_Wvel(nz+1, node)
            end if
        end do ! --> do nz=nzmax,nzmin,-1
    end do ! --> do node=1, myDim_nod2D
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    ! divide with depth dependent cell area to convert from Vertical flux to 
    ! physical vertical velocities in units m/s
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax)           
    do node=1, myDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2d(node)-1
        do nz=nzmin,nzmax
            Wvel(nz, node)=Wvel(nz, node)/area(nz, node)
            if (Fer_GM) then 
                fer_Wvel(nz, node)=fer_Wvel(nz, node)/area(nz, node)
            end if
        end do ! --> do nz=nzmax,nzmin,-1
    end do ! --> do node=1, myDim_nod2D
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    ! Add surface fresh water flux as upper boundary condition for 
    ! continutity
    if (.not. (trim(which_ale)=='linfs' )) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nzmin)       
        do node=1, myDim_nod2D
            nzmin = ulevels_nod2D(node)
            if (nzmin==1) Wvel(nzmin, node)=Wvel(nzmin, node)-water_flux(node) 
        end do ! --> do node=1, myDim_nod2D
!$OMP END PARALLEL DO        
    end if


    !___________________________________________________________________________
!$OMP MASTER    
    call exchange_nod(Wvel, partit)
    if (Fer_GM) call exchange_nod(fer_Wvel, partit)
!$OMP END MASTER    
!$OMP BARRIER

    !___________________________________________________________________________
    ! compute vertical CFL_z criteria
    call compute_CFLz(dynamics, partit, mesh)
    
    !___________________________________________________________________________
    ! compute implicite explicite splitting of vetical velocity Wvel according 
    ! to CFL_z criteria
    call compute_Wvel_split(dynamics, partit, mesh)
    
end subroutine compute_vert_vel_transpv
!
!
!_______________________________________________________________________________
! compute vertical CFL_z criteria and print out warning when critical value over
! stepped
subroutine compute_CFLz(dynamics, partit, mesh)
    use g_config, only: dt, flag_warn_cflz
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use o_PARAM
    use g_comm_auto
    implicit none
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh
    !___________________________________________________________________________
    integer       :: node, nz, nzmin, nzmax
    real(kind=WP) :: cflmax, c1, c2
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:)  , pointer :: Wvel, CFL_z
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    Wvel  => dynamics%w(:,:)
    CFL_z => dynamics%cfl_z(:,:)

    !___________________________________________________________________________
    ! calc vertical CFL criteria for debugging purpose and vertical Wvel splitting
!$OMP PARALLEL DO
    do node=1, myDim_nod2D+eDim_nod2D
       CFL_z(1,node)=0._WP
    end do
!$OMP END PARALLEL DO

    !___________________________________________________________________________
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax, c1, c2)
    do node=1, myDim_nod2D+eDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)-1
        do nz=nzmin,nzmax
            c1=abs(Wvel(nz,node)  *dt/hnode_new(nz,node)) !c1->c1(1) is made for the sake of reproducibility with the master branch (rounding error)
            c2=abs(Wvel(nz+1,node)*dt/hnode_new(nz,node)) !otherwise just add these terms (c(1) & c(2)) to CFL_z, respectively!
            ! strong condition:
            ! total volume change induced by the vertical motion
            ! no matter, upwind or downwind !
            CFL_z(nz,  node)=CFL_z(nz,node)+c1
            CFL_z(nz+1,node)=               c2
        end do ! --> do nz=nzmin,nzmax
    end do ! --> do node=1, myDim_nod2D+eDim_nod2D
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    cflmax=0.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node) REDUCTION(max:cflmax)
!$OMP DO
    do node=1, myDim_nod2D+eDim_nod2D
       cflmax=max(cflmax, maxval(CFL_z(:, node)))
    end do ! --> do node=1, myDim_nod2D+eDim_nod2D
!$OMP END DO
!$OMP END PARALLEL

    !___________________________________________________________________________
    if (cflmax > 1.0_WP .and. flag_warn_cflz) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax)
        do node=1, myDim_nod2D
            nzmin = ulevels_nod2D(node)
            nzmax = nlevels_nod2D(node)-1
            do nz=nzmin, nzmax
                if (abs(CFL_z(nz,node)-cflmax) < 1.e-12 .and. CFL_z(nz,node) > 1.75_WP .and. CFL_z(nz,node)<=2.5_WP ) then
                    print '(A, A, F4.2, A, I6, A, F7.2,A,F6.2, A, I3,I3,I3)', achar(27)//'[33m'//' --> WARNING CFLz>1.75:'//achar(27)//'[0m',&
                          'CFLz_max=',cflmax,',mstep=',mstep,',glon/glat=',geo_coord_nod2D(1,node)/rad,'/',geo_coord_nod2D(2,node)/rad,&
                          ',nz/nzmin/nzmax=',nz,nzmin,nzmax
                elseif (abs(CFL_z(nz,node)-cflmax) < 1.e-12 .and. CFL_z(nz,node) > 2.5_WP) then          
                    print '(A, A, F4.2, A, I6, A, F7.2,A,F6.2, A, I3,I3,I3)', achar(27)//'[31m'//' --> WARNING CFLz>2.5:'//achar(27)//'[0m',&
                          'CFLz_max=',cflmax,',mstep=',mstep,',glon/glat=',geo_coord_nod2D(1,node)/rad,'/',geo_coord_nod2D(2,node)/rad,&
                          ',nz/nzmin/nzmax=',nz,nzmin,nzmax
                    !!PS write(*,*) '***********************************************************'
                    !!PS write(*,*) 'max. CFL_z = ', cflmax, ' mype = ', mype
                    !!PS write(*,*) 'mstep      = ', mstep
                    !!PS write(*,*) 'glon, glat = ', geo_coord_nod2D(:,node)/rad
                    !!PS write(*,*) '2D node    = ', myList_nod2D(node)
                    !!PS write(*,*) 'nz         = ', nz
                    !!PS write(*,*) '***********************************************************'
                end if
            end do ! --> do nz=nzmin,nzmax
        end do ! --> do node=1, myDim_nod2D
!$OMP END PARALLEL DO
    end if ! --> if (cflmax > 1.0_WP .and. flag_warn_cflz) then
end subroutine compute_CFLz
!
!
!_______________________________________________________________________________ 
subroutine compute_Wvel_split(dynamics, partit, mesh)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use o_PARAM
    use g_comm_auto
    implicit none
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh
    !___________________________________________________________________________
    integer                               :: node, nz, nzmin, nzmax
    real(kind=WP)                         :: dd
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:)  , pointer :: Wvel, Wvel_e, Wvel_i, CFL_z
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    Wvel   => dynamics%w(    :,:)
    Wvel_e => dynamics%w_e(  :,:)
    Wvel_i => dynamics%w_i(  :,:)
    CFL_z  => dynamics%cfl_z(:,:)

    !___________________________________________________________________________
    ! Split implicit vertical velocity onto implicit and explicit components using CFL criteria:
    ! wsplit_maxcfl constrains the allowed explicit w according to the CFL at this place
    ! wsplit_maxcfl=1 means   w_exp  is cut at at the maximum of allowed CFL
    ! wsplit_maxcfl=0 means   w_exp  is zero (everything computed implicitly)
    ! wsplit_maxcfl=inf menas w_impl is zero (everything computed explicitly)
    ! a guess for optimal choice of wsplit_maxcfl would be 0.95
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax, dd)
    do node=1, myDim_nod2D+eDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)
        do nz=nzmin, nzmax
            Wvel_e(nz, node)=Wvel(nz, node)
            Wvel_i(nz, node)=0.0_WP
            if (dynamics%use_wsplit  .and. (CFL_z(nz, node) > dynamics%wsplit_maxcfl)) then
                dd=max((CFL_z(nz, node)-dynamics%wsplit_maxcfl), 0.0_WP)/max(dynamics%wsplit_maxcfl, 1.e-12)
                Wvel_e(nz, node)=(1.0_WP/(1.0_WP+dd))*Wvel(nz, node) !explicit part =1. if dd=0.
                Wvel_i(nz, node)=(dd    /(1.0_WP+dd))*Wvel(nz, node) !implicit part =1. if dd=inf
            end if
        end do ! --> do nz=nzmin,nzmax
    end do ! --> do node=1, myDim_nod2D+eDim_nod2D
!$OMP END PARALLEL DO
end subroutine compute_Wvel_split
!
!
!
!===============================================================================
! solve  eq.18 in S. Danilov et al. : FESOM2: from finite elements to finite volumes. 
! for (eta^(n+1)-eta^n) = d_eta
subroutine solve_ssh_ale(dynamics, partit, mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use g_comm_auto
    use g_config, only: which_ale
    use ssh_solve_preconditioner_interface
    use ssh_solve_cg_interface
    implicit none
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    logical, save        :: lfirst=.true.
    integer              :: n
    
    ! pointer on necessary derived types
    real(kind=WP), pointer  :: droptol, soltol
    integer, pointer  :: maxiter, restart, lutype, fillin, ident
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    ident   => dynamics%solverinfo%ident
    maxiter => dynamics%solverinfo%maxiter
    restart => dynamics%solverinfo%restart
    lutype  => dynamics%solverinfo%lutype
    fillin  => dynamics%solverinfo%fillin
    droptol => dynamics%solverinfo%droptol
    soltol  => dynamics%solverinfo%soltol

    if (lfirst) call ssh_solve_preconditioner(dynamics%solverinfo, partit, mesh)
    call ssh_solve_cg(dynamics%d_eta, dynamics%ssh_rhs, dynamics%solverinfo, partit, mesh)
    call exchange_nod(dynamics%d_eta, partit) !is this required after calling psolve ?
    lfirst=.false.

end subroutine solve_ssh_ale
!
!
!===============================================================================
subroutine impl_vert_visc_ale(dynamics, partit, mesh)
    USE MOD_MESH
    USE o_PARAM
    USE o_ARRAYS, only: Av, stress_surf
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE g_CONFIG !,only: dt
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    real(kind=WP)              ::  a(mesh%nl-1), b(mesh%nl-1), c(mesh%nl-1), ur(mesh%nl-1), vr(mesh%nl-1)
    real(kind=WP)              ::  cp(mesh%nl-1), up(mesh%nl-1), vp(mesh%nl-1)
    integer                    ::  nz, elem, nzmin, nzmax, elnodes(3)
    real(kind=WP)              ::  zinv, m, friction, wu, wd
    real(kind=WP)              ::  zbar_n(mesh%nl), Z_n(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: Wvel_i
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV     =>dynamics%uv(:,:,:)
    UV_rhs =>dynamics%uv_rhs(:,:,:)
    Wvel_i =>dynamics%w_i(:,:)

    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a, b, c, ur, vr, cp, up, vp, elem, nz, nzmin, nzmax, elnodes, &
!$OMP                                                          zinv, m, friction, wu, wd, zbar_n, Z_n)

!$OMP DO
    DO elem=1,myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        nzmin = ulevels(elem)
        nzmax = nlevels(elem)
        
        !___________________________________________________________________________
        ! Here can not exchange zbar_n & Z_n with zbar_3d_n & Z_3d_n because 
        ! they run over elements here 
        zbar_n=0.0_WP
        Z_n   =0.0_WP
        ! in case of partial cells zbar_n(nzmax) is not any more at zbar(nzmax), 
        ! zbar_n(nzmax) is now zbar_e_bot(elem), 
        zbar_n(nzmax)=zbar_e_bot(elem)
        Z_n(nzmax-1)=zbar_n(nzmax) + helem(nzmax-1,elem)/2.0_WP
        !!PS do nz=nzmax-1,2,-1
        do nz=nzmax-1,nzmin+1,-1
            zbar_n(nz) = zbar_n(nz+1) + helem(nz,elem)
            Z_n(nz-1) = zbar_n(nz) + helem(nz-1,elem)/2.0_WP
        end do
        !!PS zbar_n(1) = zbar_n(2) + helem(1,elem)
        zbar_n(nzmin) = zbar_n(nzmin+1) + helem(nzmin,elem)
        
        !___________________________________________________________________________
        ! Operator
        ! Regular part of coefficients:
        !!PS do nz=2, nzmax-2
        do nz=nzmin+1, nzmax-2
            zinv=1.0_WP*dt/(zbar_n(nz)-zbar_n(nz+1))
            a(nz)=-Av(nz,elem)/(Z_n(nz-1)-Z_n(nz))*zinv
            c(nz)=-Av(nz+1,elem)/(Z_n(nz)-Z_n(nz+1))*zinv
            b(nz)=-a(nz)-c(nz)+1.0_WP
            ! Update from the vertical advection
            wu=sum(Wvel_i(nz,   elnodes))/3._WP
            wd=sum(Wvel_i(nz+1, elnodes))/3._WP
            a(nz)=a(nz)+min(0._WP, wu)*zinv
            b(nz)=b(nz)+max(0._WP, wu)*zinv
            
            b(nz)=b(nz)-min(0._WP, wd)*zinv
            c(nz)=c(nz)-max(0._WP, wd)*zinv
        end do
        ! The last row
        zinv=1.0_WP*dt/(zbar_n(nzmax-1)-zbar_n(nzmax))
        a(nzmax-1)=-Av(nzmax-1,elem)/(Z_n(nzmax-2)-Z_n(nzmax-1))*zinv
        b(nzmax-1)=-a(nzmax-1)+1.0_WP
        c(nzmax-1)=0.0_WP
        
        ! Update from the vertical advection
        wu=sum(Wvel_i(nzmax-1, elnodes))/3._WP
        a(nzmax-1)=a(nzmax-1)+min(0._WP, wu)*zinv
        b(nzmax-1)=b(nzmax-1)+max(0._WP, wu)*zinv
        
        ! The first row
        !!PS zinv=1.0_WP*dt/(zbar_n(1)-zbar_n(2))
        !!PS c(1)=-Av(2,elem)/(Z_n(1)-Z_n(2))*zinv
        !!PS a(1)=0.0_WP
        !!PS b(1)=-c(1)+1.0_WP
        zinv=1.0_WP*dt/(zbar_n(nzmin)-zbar_n(nzmin+1))
        c(nzmin)=-Av(nzmin+1,elem)/(Z_n(nzmin)-Z_n(nzmin+1))*zinv
        a(nzmin)=0.0_WP
        b(nzmin)=-c(nzmin)+1.0_WP
        
        ! Update from the vertical advection
        !!PS wu=sum(Wvel_i(1, elnodes))/3._WP
        !!PS wd=sum(Wvel_i(2, elnodes))/3._WP
        wu=sum(Wvel_i(nzmin, elnodes))/3._WP
        wd=sum(Wvel_i(nzmin+1, elnodes))/3._WP
        
        !!PS b(1)=b(1)+wu*zinv
        !!PS b(1)=b(1)-min(0._WP, wd)*zinv
        !!PS c(1)=c(1)-max(0._WP, wd)*zinv
        b(nzmin)=b(nzmin)+wu*zinv
        b(nzmin)=b(nzmin)-min(0._WP, wd)*zinv
        c(nzmin)=c(nzmin)-max(0._WP, wd)*zinv
        
        ! ===========================
        ! The rhs:
        ! ===========================
        !!PS ur(1:nzmax-1)=UV_rhs(1,1:nzmax-1,elem)
        !!PS vr(1:nzmax-1)=UV_rhs(2,1:nzmax-1,elem)
        ur(nzmin:nzmax-1)=UV_rhs(1,nzmin:nzmax-1,elem)
        vr(nzmin:nzmax-1)=UV_rhs(2,nzmin:nzmax-1,elem)
        
        ! The first row contains surface forcing
        !!PS ur(1)= ur(1)+zinv*stress_surf(1,elem)/density_0
        !!PS vr(1)= vr(1)+zinv*stress_surf(2,elem)/density_0
        ur(nzmin)= ur(nzmin)+zinv*stress_surf(1,elem)/density_0
        vr(nzmin)= vr(nzmin)+zinv*stress_surf(2,elem)/density_0

        if (dynamics%ldiag_ke) then
           dynamics%ke_wind(1,elem)=stress_surf(1,elem)/density_0*dt
           dynamics%ke_wind(2,elem)=stress_surf(2,elem)/density_0*dt
        end if
        
        ! The last row contains bottom friction
        zinv=1.0_WP*dt/(zbar_n(nzmax-1)-zbar_n(nzmax))
        !!PS friction=-C_d*sqrt(UV(1,nlevels(elem)-1,elem)**2+ &
        !!PS             UV(2,nlevels(elem)-1,elem)**2)
        
        if ((toy_ocean) .AND. (TRIM(which_toy)=="dbgyre")) then
           friction=-C_d

        else if ((toy_ocean) .AND. (TRIM(which_toy)=="soufflet")) then
           friction=-C_d

        else
           friction=-C_d*sqrt(UV(1,nzmax-1,elem)**2+ &
                           UV(2,nzmax-1,elem)**2)
        end if

        ur(nzmax-1)=ur(nzmax-1)+zinv*friction*UV(1,nzmax-1,elem)
        vr(nzmax-1)=vr(nzmax-1)+zinv*friction*UV(2,nzmax-1,elem)

        if (dynamics%ldiag_ke) then
           dynamics%ke_drag(1,elem)=friction*UV(1,nzmax-1,elem)*dt
           dynamics%ke_drag(2,elem)=friction*UV(2,nzmax-1,elem)*dt
        end if

        ! Model solves for the difference to the timestep N and therefore we need to 
        ! update the RHS for advective and diffusive contributions at the previous time step
        !!PS do nz=2, nzmax-2
        do nz=nzmin+1, nzmax-2
            ur(nz)=ur(nz)-a(nz)*UV(1,nz-1,elem)-(b(nz)-1.0_WP)*UV(1,nz,elem)-c(nz)*UV(1,nz+1,elem)
            vr(nz)=vr(nz)-a(nz)*UV(2,nz-1,elem)-(b(nz)-1.0_WP)*UV(2,nz,elem)-c(nz)*UV(2,nz+1,elem)
        end do
        !!PS ur(1)=ur(1)-(b(1)-1.0_WP)*UV(1,1,elem)-c(1)*UV(1,2,elem)
        !!PS vr(1)=vr(1)-(b(1)-1.0_WP)*UV(2,1,elem)-c(1)*UV(2,2,elem)
        ur(nzmin)=ur(nzmin)-(b(nzmin)-1.0_WP)*UV(1,nzmin,elem)-c(nzmin)*UV(1,nzmin+1,elem)
        vr(nzmin)=vr(nzmin)-(b(nzmin)-1.0_WP)*UV(2,nzmin,elem)-c(nzmin)*UV(2,nzmin+1,elem)
        
        ur(nzmax-1)=ur(nzmax-1)-a(nzmax-1)*UV(1,nzmax-2,elem)-(b(nzmax-1)-1.0_WP)*UV(1,nzmax-1,elem)
        vr(nzmax-1)=vr(nzmax-1)-a(nzmax-1)*UV(2,nzmax-2,elem)-(b(nzmax-1)-1.0_WP)*UV(2,nzmax-1,elem)
        
        ! ===========================
        ! The sweep algorithm
        ! ===========================
        ! initialize c-prime and s,t-prime
        !!PS cp(1) = c(1)/b(1)
        !!PS up(1) = ur(1)/b(1)
        !!PS vp(1) = vr(1)/b(1)
        cp(nzmin) = c(nzmin)/b(nzmin)
        up(nzmin) = ur(nzmin)/b(nzmin)
        vp(nzmin) = vr(nzmin)/b(nzmin)
        
        ! solve for vectors c-prime and t, s-prime
        !!PS do nz = 2,nzmax-1
        do nz = nzmin+1,nzmax-1
            m = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            up(nz) = (ur(nz)-up(nz-1)*a(nz))/m
            vp(nz) = (vr(nz)-vp(nz-1)*a(nz))/m
        enddo
        ! initialize x
        ur(nzmax-1) = up(nzmax-1)
        vr(nzmax-1) = vp(nzmax-1)
        
        ! solve for x from the vectors c-prime and d-prime
        !!PS do nz = nzmax-2, 1, -1
        do nz = nzmax-2, nzmin, -1
            ur(nz) = up(nz)-cp(nz)*ur(nz+1)
            vr(nz) = vp(nz)-cp(nz)*vr(nz+1)
        end do
        
        ! ===========================
        ! RHS update
        ! ===========================
        !!PS do nz=1,nzmax-1
        do nz=nzmin,nzmax-1
            UV_rhs(1,nz,elem)=ur(nz)
            UV_rhs(2,nz,elem)=vr(nz)
        end do
    end do   !!! cycle over elements
!$OMP END DO
!$OMP END PARALLEL
end subroutine impl_vert_visc_ale
!
!
!===============================================================================
subroutine oce_timestep_ale(n, ice, dynamics, tracers, partit, mesh)
    use g_config
    use MOD_MESH
    use MOD_TRACER
    use MOD_DYN
    USE MOD_ICE
    use o_ARRAYS
    use o_PARAM
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_comm_auto
    use io_RESTART !PS
    use o_mixing_KPP_mod
#if defined (__cvmix)       
    use g_cvmix_tke
    use g_cvmix_idemix
    use g_cvmix_pp
    use g_cvmix_kpp
    use g_cvmix_tidal
#endif    
    use Toy_Channel_Soufflet
    use oce_ale_interfaces
    use compute_vert_vel_transpv_interface
    use compute_ssh_split_explicit_interface
    use pressure_bv_interface
    use pressure_force_4_linfs_interface
    use pressure_force_4_zxxxx_interface
    use compute_vel_rhs_interface
    use solve_tracers_ale_interface
    use write_step_info_interface
    use check_blowup_interface
    use fer_solve_interface
    use impl_vert_visc_ale_vtransp_interface
#if defined (FESOM_PROFILING)
    use fesom_profiler
#endif
    
    IMPLICIT NONE
    integer       , intent(in)            :: n
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    type(t_ice)   , intent(inout), target :: ice
    !___________________________________________________________________________
    real(kind=8)      :: t0,t1, t2, t30, t3, t4, t5, t6, t7, t8, t9, t10, loc, glo
    integer           :: node
    integer           :: nz, elem, nzmin, nzmax !for KE diagnostic
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer :: eta_n
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    eta_n => dynamics%eta_n(:)
    
    !___________________________________________________________________________
    t0=MPI_Wtime()
!PS     water_flux = 0.0_WP
!PS     heat_flux  = 0.0_WP
!PS     stress_surf= 0.0_WP
!PS     stress_node_surf= 0.0_WP

#if defined (FESOM_PROFILING)
    call fesom_profiler_start("oce_mix_pres")
#endif
    !___________________________________________________________________________
    ! calculate equation of state, density, pressure and mixed layer depths
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call pressure_bv'//achar(27)//'[0m'
    call pressure_bv(tracers, partit, mesh)            !!!!! HeRE change is made. It is linear EoS now.

    !___________________________________________________________________________
    ! calculate calculate pressure gradient force
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call pressure_force_4_...'//achar(27)//'[0m'
    if (trim(which_ale)=='linfs') then
        call pressure_force_4_linfs(tracers, partit, mesh)
    else  
        call pressure_force_4_zxxxx(tracers, partit, mesh)
    end if
    
    !___________________________________________________________________________
    ! check validity of visc_opt=5 selection
    ! --> need to know buoyancy frequency to do so.
    ! --> only check on the first timestep 
    if (n==1) call check_viscopt(dynamics, partit, mesh)
    
    !___________________________________________________________________________
    ! calculate alpha and beta
    ! it will be used for KPP, Redi, GM etc. Shall we keep it on in general case?
    call sw_alpha_beta(tracers%data(1)%values, tracers%data(2)%values, partit, mesh)

    ! computes the xy gradient of a neutral surface; will be used by Redi, GM etc.
    call compute_sigma_xy(tracers%data(1)%values,tracers%data(2)%values, partit, mesh)

    ! compute both: neutral slope and tapered neutral slope. Can be later combined with compute_sigma_xy
    ! will be primarily used for computing Redi diffusivities. etc?
    call compute_neutral_slope(partit, mesh)
 
    !___________________________________________________________________________
    call status_check(partit)
    !___________________________________________________________________________
    ! >>>>>>                                                             <<<<<<
    ! >>>>>>    calculate vertical mixing coefficients for tracer (Kv)   <<<<<<
    ! >>>>>>    and momentum (Av) using using FESOM implementations      <<<<<<
    ! >>>>>>    for PP and KPP as well as mixing schemes provided by     <<<<<<
    ! >>>>>>    by the CVMIX library                                     <<<<<< 
    ! >>>>>>                                                             <<<<<<
    !___________________________________________________________________________
    
    !___EXTENSION OF MIXING SCHEMES_____________________________________________
    ! add CVMIX IDEMIX (internal wave energy) parameterisation for 
    ! vertical mixing (dissipation of energy by internal gravity waves) 
    ! extension from Olbers and Eden, 2013, "A global Model for the diapycnal 
    ! diffusivity induced by internal gravity waves" --> use together with 
    ! cvmix_TKE (mix_scheme='cvmix_TKE+cvmix_IDEMIX') or standalone for debbuging 
    ! (mix_scheme='cvmix_IDEMIX') --> If idemix is used together with tke it needs 
    ! to be called prior to tke
    ! for debugging
#if defined (__cvmix)       
    if  (mod(mix_scheme_nmb,10)==6) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call calc_cvmix_idemix'//achar(27)//'[0m'
        call calc_cvmix_idemix(partit, mesh)
    end if 
#endif    

    !___MAIN MIXING SCHEMES_____________________________________________________
    ! use FESOM2.0 tuned k-profile parameterization for vertical mixing 
    if (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call oce_mixing_KPP'//achar(27)//'[0m' 
        call oce_mixing_KPP(Av, Kv_double, dynamics, tracers, partit, mesh)
!$OMP PARALLEL DO
        do node=1, myDim_nod2D+eDim_nod2D
           Kv(:, node)=Kv_double(:, node, 1)
        end do
!$OMP END PARALLEL DO

        call mo_convect(ice, partit, mesh)
        
    ! use FESOM2.0 tuned pacanowski & philander parameterization for vertical 
    ! mixing     
    else if(mix_scheme_nmb==2 .or. mix_scheme_nmb==27) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call oce_mixing_PP'//achar(27)//'[0m' 
        call oce_mixing_PP(dynamics, partit, mesh)
        call mo_convect(ice, partit, mesh)
#if defined (__cvmix)           
    ! use CVMIX KPP (Large at al. 1994) 
    else if(mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call calc_cvmix_kpp'//achar(27)//'[0m'
        call calc_cvmix_kpp(ice, dynamics, tracers, partit, mesh)
        call mo_convect(ice, partit, mesh)
        
    ! use CVMIX PP (Pacanowski and Philander 1981) parameterisation for mixing
    ! based on Richardson number Ri = N^2/(du/dz)^2, using Brunt Väisälä frequency
    ! N^2 and vertical horizontal velocity shear dui/dz
    else if(mix_scheme_nmb==4 .or. mix_scheme_nmb==47) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call calc_cvmix_pp'//achar(27)//'[0m'
        call calc_cvmix_pp(dynamics, partit, mesh)
        call mo_convect(ice, partit, mesh)
        
    ! use CVMIX TKE (turbulent kinetic energy closure) parameterisation for 
    ! vertical mixing with or without the IDEMIX (dissipation of energy by 
    ! internal gravity waves) extension from Olbers and Eden, 2013, "A global 
    ! Model for the diapycnal diffusivity induced by internal gravity waves" 
    else if(mix_scheme_nmb==5 .or. mix_scheme_nmb==56) then    
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call calc_cvmix_tke'//achar(27)//'[0m'
        call calc_cvmix_tke(dynamics, partit, mesh)
        call mo_convect(ice, partit, mesh)
#endif    
    end if   

#if defined (__cvmix)       
    !___EXTENSION OF MIXING SCHEMES_____________________________________________
    ! add CVMIX TIDAL mixing scheme of Simmons et al. 2004 "Tidally driven mixing 
    ! in a numerical model of the ocean general circulation", ocean modelling to 
    ! the already computed viscosities/diffusivities of KPP, PP, cvmix_KPP or 
    ! cvmix_PP --> use standalone for debugging --> needs to be called after main
    ! mixing schemes
    if ( mod(mix_scheme_nmb,10)==7) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call calc_cvmix_tidal'//achar(27)//'[0m'
        call calc_cvmix_tidal(partit, mesh)
        
    end if
#endif    
    t1=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_mix_pres")
    call fesom_profiler_start("oce_dyn_momentum")
#endif    
    
    !___________________________________________________________________________
    ! add contribution from momentum advection, coriolis and pressure gradient |
    ! force to UV_rhs
    ! UV_rhs = dt*[ (R_advec + R_coriolis)_AB2^n + R_pressure^n ]
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call compute_vel_rhs'//achar(27)//'[0m'
    call compute_vel_rhs(ice, dynamics, partit, mesh)
     
    !___________________________________________________________________________
    ! Energy diagnostic contribution
    if (dynamics%ldiag_ke) then
        ! if use solver
        if (.not. dynamics%use_ssh_se_subcycl) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax)
            do elem=1, myDim_elem2D
                nzmax = nlevels(elem)
                nzmin = ulevels(elem)
                do nz=nzmin,nzmax-1
                    dynamics%ke_rhs_bak(:,nz,elem)=dynamics%UV_rhs(:,nz,elem)
                end do
            end do
!$OMP END PARALLEL DO
        
        ! if use splitexpl subcycl. UV_rhs in units of transport --> therefor *1/helem
        else
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax)
            do elem=1, myDim_elem2D
                nzmax = nlevels(elem)
                nzmin = ulevels(elem)
                do nz=nzmin,nzmax-1
                    dynamics%ke_rhs_bak(:,nz,elem)=dynamics%UV_rhs(:,nz,elem)/helem(nz,elem) 
                end do
            end do
!$OMP END PARALLEL DO        
        end if 
    end if
    
    !___________________________________________________________________________
    ! add contribution from horizontal viscosity to UV_rhs
    ! UV_rhs = dt*[ (R_advec + R_coriolis)^n + R_pressure + R_hviscos] 
    ! UV_rhs = UV_rhs + dt*R_hviscos
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call viscosity_filter'//achar(27)//'[0m'
    call viscosity_filter(dynamics%opt_visc, dynamics, partit, mesh)
    
    !___________________________________________________________________________
    ! Energy diagnostic contribution
    if (dynamics%ldiag_ke) then
        ! if use solver
        if (.not. dynamics%use_ssh_se_subcycl) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax)
            do elem=1, myDim_elem2D
                nzmax = nlevels(elem)
                nzmin = ulevels(elem)
                do nz=nzmin,nzmax-1
                    dynamics%ke_hvis(:,nz,elem)=dynamics%UV_rhs(:,nz,elem) - dynamics%ke_rhs_bak(:,nz,elem)
                end do
            end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax)
            do elem=1, myDim_elem2D
                nzmax = nlevels(elem)
                nzmin = ulevels(elem)
                do nz=nzmin,nzmax-1
                    dynamics%ke_rhs_bak(:,nz,elem)=dynamics%UV_rhs(:,nz,elem)
                end do
            end do
!$OMP END PARALLEL DO

        ! if use splitexpl subcycl. UV_rhs in units of transport --> therefor *1/helem
        else
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax)
            do elem=1, myDim_elem2D
                nzmax = nlevels(elem)
                nzmin = ulevels(elem)
                do nz=nzmin,nzmax-1
                    dynamics%ke_hvis(:,nz,elem)=dynamics%UV_rhs(:,nz,elem)/helem(nz,elem) - dynamics%ke_rhs_bak(:,nz,elem)
                end do
            end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax)
            do elem=1, myDim_elem2D
                nzmax = nlevels(elem)
                nzmin = ulevels(elem)
                do nz=nzmin,nzmax-1
                    dynamics%ke_rhs_bak(:,nz,elem)=dynamics%UV_rhs(:,nz,elem)/helem(nz,elem) 
                end do
            end do
!$OMP END PARALLEL DO        
        
        end if 
    end if
    
    !___________________________________________________________________________
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call impl_vert_visc_ale'//achar(27)//'[0m'
    if(dynamics%use_ivertvisc) then
        if ( .not. dynamics%use_ssh_se_subcycl ) then
            call impl_vert_visc_ale(dynamics,partit, mesh)
        else
            call impl_vert_visc_ale_vtransp(dynamics, partit, mesh)
        end if 
    end if
    
    !___________________________________________________________________________
    if (dynamics%ldiag_ke) then
        ! if use solver
        if (.not. dynamics%use_ssh_se_subcycl) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, nz, nzmin, nzmax)
            do elem=1, myDim_elem2D
                nzmax = nlevels(elem)
                nzmin = ulevels(elem)
                do nz=nzmin,nzmax-1
                    dynamics%ke_vvis(:,nz,elem)=dynamics%UV_rhs(:,nz,elem) - dynamics%ke_rhs_bak(:,nz,elem)
                end do
            end do
!$OMP END PARALLEL DO

        ! if use splitexpl subcycl. UV_rhs in units of transport --> therefor *1/helem
        else
            do elem=1, myDim_elem2D
                nzmax = nlevels(elem)
                nzmin = ulevels(elem)
                do nz=nzmin,nzmax-1
                    dynamics%ke_vvis(:,nz,elem)=dynamics%UV_rhs(:,nz,elem)/helem(nz,elem) - dynamics%ke_rhs_bak(:,nz,elem)
                end do
            end do
        end if 
    end if
    t2=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_dyn_momentum")
    call fesom_profiler_start("oce_ssh_solve")
#endif        
    !___________________________________________________________________________
    ! >->->->->->->->->->->->->     ALE-part starts     <-<-<-<-<-<-<-<-<-<-<-<-
    !___________________________________________________________________________
    
    !___________________________________________________________________________
    ! Compute SSH via solver
    ! Update stiffness matrix by dhe=hbar(n+1/2)-hbar(n-1/2) on elements, only
    ! needed for zlevel and zstar
    if (.not. dynamics%use_ssh_se_subcycl) then
        if (.not. trim(which_ale)=='linfs') call update_stiff_mat_ale(partit, mesh)
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call compute_ssh_rhs_ale'//achar(27)//'[0m'
        ! ssh_rhs=-alpha*\nabla\int(U_n+U_rhs)dz-(1-alpha)*...
        ! see "FESOM2: from finite elements to finte volumes, S. Danilov..." eq. (18) rhs
        call compute_ssh_rhs_ale(dynamics, partit, mesh)

        ! Take updated ssh matrix and solve --> new ssh!
        t30=MPI_Wtime() 
        call solve_ssh_ale(dynamics, partit, mesh)
        
        if ((toy_ocean) .AND. (TRIM(which_toy)=="soufflet")) then
            if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call relax_zonal_vel'//achar(27)//'[0m'
            call relax_zonal_vel(dynamics, partit, mesh)
        end if     
        t3=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_ssh_solve")
    call fesom_profiler_start("oce_vel_update")
#endif 

        ! estimate new horizontal velocity u^(n+1)
        ! u^(n+1) = u* + [-g * tau * theta * grad(eta^(n+1)-eta^(n)) ]
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call update_vel'//achar(27)//'[0m'
        ! ke will be computed inside there if dynamics%ldiag_ke is .TRUE.
        call update_vel(dynamics, partit, mesh)
        
        ! --> eta_(n) --> eta_(n+1) = eta_(n) + deta = eta_(n) + (eta_(n+1) + eta_(n))
        t4=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_vel_update")
    call fesom_profiler_start("oce_hbar_calc")
#endif 
        
        ! Update to hbar(n+3/2) and compute dhe to be used on the next step
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call compute_hbar_ale'//achar(27)//'[0m'
        call compute_hbar_ale(dynamics, partit, mesh)

        !___________________________________________________________________________
        ! - Current dynamic elevation alpha*hbar(n+1/2)+(1-alpha)*hbar(n-1/2)
        !   equation (14) Danlov et.al "the finite volume sea ice ocean model FESOM2
        !   ...if we do it here we don't need to write hbar_old into a restart file...
        ! - where(ulevels_nod2D==1) is used here because of the rigid lid 
        !   approximation under the cavity 
        ! - at points in the cavity the time derivative term in ssh matrix will be 
        !   omitted; and (14) will not be applied at cavity points. Additionally,
        !   since there is no real elevation, but only surface pressure, there is 
        !   no layer motion under the cavity. In this case the ice sheet acts as a 
        !   rigid lid.
!$OMP PARALLEL DO
        do node=1, myDim_nod2D+eDim_nod2D
            if (ulevels_nod2D(node)==1) eta_n(node)=alpha*hbar(node)+(1.0_WP-alpha)*hbar_old(node)
        end do
!$OMP END PARALLEL DO
        ! --> eta_(n)
        ! call zero_dynamics !DS, zeros several dynamical variables; to be used for testing new implementations!
        t5=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_hbar_calc")
    call fesom_profiler_start("oce_gm_redi")
#endif 
    
    !___________________________________________________________________________
    ! Compute SSH via split-explicite subcycling
    else    
        ! Compute vertical integral of transport velocity rhs omitting the contributions from
        ! the elevation and Coriolis. 
        t30=MPI_Wtime()
        call compute_BT_rhs_SE_vtransp(dynamics, partit, mesh)
        
        ! Do barotropic step, get eta_{n+1} and BT transport 
        call compute_BT_step_SE_ale(dynamics, partit, mesh)
        t3=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_ssh_solve")
    call fesom_profiler_start("oce_vel_update")
#endif        
        ! Trim U to be consistent with BT transport
        call update_trim_vel_ale_vtransp(1, dynamics, partit, mesh) 
        t4=MPI_Wtime()
        t5=t4
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_vel_update")
    call fesom_profiler_start("oce_hbar_calc")
    call fesom_profiler_end("oce_hbar_calc")
    call fesom_profiler_start("oce_gm_redi")
#endif
    end if ! --> if (.not. dynamics%use_ssh_se_subcycl) then
    
    !___________________________________________________________________________
    ! Do horizontal and vertical scaling of GM/Redi  diffusivity 
    if (Fer_GM .or. Redi) then
        call init_Redi_GM(partit, mesh)
    end if
    
    ! Implementation of Gent & McWiliams parameterization after R. Ferrari et al., 2010
    ! does not belong directly to ALE formalism
    if (Fer_GM) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call fer_solve_Gamma'//achar(27)//'[0m'
        call fer_solve_Gamma(partit, mesh)
        call fer_gamma2vel(dynamics, partit, mesh)
    end if
    t6=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_gm_redi")
    call fesom_profiler_start("oce_vert_vel")
#endif 
 
    !___________________________________________________________________________
    ! keep the old vertical velocity for computation of the mean between the timesteps (is used in compute_ke_wrho)
    if (dynamics%ldiag_ke) then
!$OMP PARALLEL DO
        do node=1, myDim_nod2D+eDim_nod2D
        dynamics%w_old(:, node)=dynamics%w(:, node)
        end do
!$OMP END PARALLEL DO
    end if
    
    !___________________________________________________________________________
    ! The main step of ALE procedure --> this is were the magic happens --> here 
    ! is decided how change in hbar is distributed over the vertical layers
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call vert_vel_ale'//achar(27)//'[0m'
    if ( .not. dynamics%use_ssh_se_subcycl) then
        call vert_vel_ale(dynamics, partit, mesh)
    else
        if (trim(which_ale)=='zstar' ) then
            call compute_thickness_zstar(dynamics, partit, mesh)
        else
            hnode_new = hnode
        end if 
        call compute_vert_vel_transpv(dynamics, partit, mesh)
    end if    
    t7=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_vert_vel")
    call fesom_profiler_start("oce_tracer_solve")
#endif   
     
    !___________________________________________________________________________
    ! energy diagnostic computation
    if (dynamics%ldiag_ke) then
       call compute_ke_wrho(dynamics, partit, mesh)
       call compute_apegen (dynamics, tracers, partit, mesh)
    end if
 
    !___________________________________________________________________________
    ! solve tracer equation
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call solve_tracers_ale'//achar(27)//'[0m'
    call solve_tracers_ale(ice, dynamics, tracers, partit, mesh)
    t8=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_tracer_solve")
    call fesom_profiler_start("oce_thickness_update")
#endif 
     
    !___________________________________________________________________________
    ! Update hnode=hnode_new, helem
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call update_thickness_ale'//achar(27)//'[0m'
    call update_thickness_ale(partit, mesh)
    t9=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_thickness_update")
    call fesom_profiler_start("oce_blowup_check")
#endif 
    
!PS     !___________________________________________________________________________
!PS     ! Trim to make velocity consistent with BT velocity at n+1/2
!PS     ! Velocity will be used to advance momentum
!PS     if (dynamics%use_ssh_se_subcycl) then
!PS         call update_trim_vel_ale_vtransp(2, dynamics, partit, mesh)   
!PS     end if
    
    !___________________________________________________________________________
    ! write out global fields for debugging
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call write_step_info'//achar(27)//'[0m'
    call write_step_info(n,logfile_outfreq, ice, dynamics, tracers, partit, mesh)
    
    !___________________________________________________________________________
    ! write energy diagnostic info (dynamics%ldiag_ke = .true.)
    if ( (dynamics%ldiag_ke) .and. (mod(n,logfile_outfreq)==0) ) then
        call write_enegry_info(dynamics, partit, mesh)
    end if
    
    ! check model for blowup --> ! write_step_info and check_blowup require 
    ! togeather around 2.5% of model runtime
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call check_blowup'//achar(27)//'[0m'
    call check_blowup(n, ice, dynamics, tracers, partit, mesh)
    t10=MPI_Wtime()
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("oce_blowup_check")
#endif

    !___________________________________________________________________________
    ! write out execution times for ocean step parts
    rtime_oce          = rtime_oce + (t10-t0)-(t10-t9)
    rtime_oce_mixpres  = rtime_oce_mixpres + (t1-t0)
    rtime_oce_dyn      = rtime_oce_dyn + (t2-t1)+(t7-t6)+(t4-t3)
    rtime_oce_dynssh   = rtime_oce_dynssh + (t3-t2)+(t5-t4)
    rtime_oce_solvessh = rtime_oce_solvessh + (t3-t30)
    rtime_oce_GMRedi   = rtime_oce_GMRedi + (t6-t5)
    rtime_oce_solvetra = rtime_oce_solvetra + (t8-t7)
    rtime_tot          = rtime_tot + (t10-t0)-(t10-t9)
    if(mod(n,logfile_outfreq)==0 .and. mype==0) then  
        write(*,*) '___ALE OCEAN STEP EXECUTION TIMES______________________'
        write(*,"(A, ES10.3)") '     Oce. Mix,Press.. :', t1-t0
        write(*,"(A, ES10.3)") '     Oce. Dynamics    :', t2-t1
        write(*,"(A, ES10.3)") '     Oce. Update Vel. :', t4-t3
        write(*,"(A, ES10.3)") '     Oce. Fer-GM.     :', t6-t5
        write(*,*) '    _______________________________'
        write(*,"(A, ES10.3)") '     ALE-Solve SSH    :', t3-t2
        write(*,"(A, ES10.3)") '     ALE-Calc. hbar   :', t5-t4
        write(*,"(A, ES10.3)") '     ALE-Update+W     :', t7-t6
        write(*,"(A, ES10.3)") '     ALE-Solve Tracer :', t8-t7
        write(*,"(A, ES10.3)") '     ALE-Update hnode :', t9-t8
        write(*,*) '    _______________________________'
        write(*,"(A, ES10.3)") '     check for blowup :', t10-t9
        write(*,*) '    _______________________________'
        write(*,"(A, ES10.3)") '     Oce. TOTAL       :', t10-t0
        write(*,*)
        write(*,*)
    end if    
end subroutine oce_timestep_ale

