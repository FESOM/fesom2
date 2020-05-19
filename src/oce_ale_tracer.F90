module diff_part_hor_redi_interface
  interface
    subroutine diff_part_hor_redi(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module adv_tracers_muscle_ale_interface
  interface
    subroutine adv_tracers_muscle_ale(ttfAB, num_ord, do_Xmoment, mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh), intent(in) , target :: mesh
      integer                  :: do_Xmoment
      real(kind=WP)            :: ttfAB(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP)            :: num_ord
    end subroutine
  end interface
end module
module adv_tracers_vert_ppm_ale_interface
  interface
    subroutine adv_tracers_vert_ppm_ale(ttf, do_Xmoment, mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh), intent(in) , target :: mesh
      integer                  :: do_Xmoment
      real(kind=WP)            :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    end subroutine
  end interface
end module
module adv_tracers_ale_interface
  interface
    subroutine adv_tracers_ale(tr_num, mesh)
      use mod_mesh
      integer :: tr_num
      type(t_mesh), intent(in) , target :: mesh    
    end subroutine
  end interface
end module
module diff_ver_part_expl_ale_interface
  interface
    subroutine diff_ver_part_expl_ale(tr_num, mesh)
      use MOD_MESH
      type(t_mesh), intent(in) , target :: mesh    
      integer                  :: tr_num
    end subroutine
  end interface
end module
module diff_ver_part_redi_expl_interface
  interface
    subroutine diff_ver_part_redi_expl(mesh)
      use MOD_MESH
      type(t_mesh), intent(in) , target :: mesh    
    end subroutine
  end interface
end module
module diff_ver_part_impl_ale_interface
  interface
    subroutine diff_ver_part_impl_ale(tr_num, mesh)
      use MOD_MESH
      type(t_mesh), intent(in) , target :: mesh
      integer                  :: tr_num
    end subroutine
  end interface
end module
module diff_tracers_ale_interface
  interface
    subroutine diff_tracers_ale(tr_num, mesh)
      use mod_mesh
      integer, intent(in)      :: tr_num
      type(t_mesh), intent(in) , target :: mesh
    end subroutine
  end interface
end module

!
!
!===============================================================================
! Driving routine    Here with ALE changes!!!
subroutine solve_tracers_ale(mesh)
    use g_config, only: flag_debug
    use g_parsup
    use o_PARAM, only: tracer_adv, num_tracers, SPP, Fer_GM
    use o_arrays
    use mod_mesh
    use g_comm_auto
    use o_tracers
    use adv_tracers_ale_interface
    use diff_tracers_ale_interface
    
    implicit none
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: tr_num, node, nzmax
    real(kind=WP)            :: aux_tr(mesh%nl-1,myDim_nod2D+eDim_nod2D)

#include "associate_mesh.h"
    !___________________________________________________________________________
    if (SPP) call cal_rejected_salt(mesh)
    if (SPP) call app_rejected_salt(mesh)
    
    !___________________________________________________________________________
    ! update 3D velocities with the bolus velocities:
    ! 1. bolus velocities are computed according to GM implementation after R. Ferrari et al., 2010
    ! 2. bolus velocities are used only for advecting tracers and shall be subtracted back afterwards
    if (Fer_GM) then
        UV    =UV    +fer_UV
        Wvel_e=Wvel_e+fer_Wvel
        Wvel  =Wvel  +fer_Wvel
    end if
    !___________________________________________________________________________
    ! loop over all tracers 
    do tr_num=1,num_tracers
        ! do tracer AB (Adams-Bashfort) interpolation only for advectiv part 
        ! needed
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call init_tracers_AB'//achar(27)//'[0m'
        call init_tracers_AB(tr_num, mesh)
        
        ! advect tracers
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call adv_tracers_ale'//achar(27)//'[0m'
        call adv_tracers_ale(tr_num, mesh)
        
        ! diffuse tracers 
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call diff_tracers_ale'//achar(27)//'[0m'
        call diff_tracers_ale(tr_num, mesh)
        
        ! relax to salt and temp climatology
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call relax_to_clim'//achar(27)//'[0m'
        call relax_to_clim(tr_num, mesh)
        
        call exchange_nod(tr_arr(:,:,tr_num))
    end do
    
    !___________________________________________________________________________
    do tr_num=1, ptracers_restore_total           
        tr_arr(:,ptracers_restore(tr_num)%ind2,ptracers_restore(tr_num)%locid)=1.0_WP
    end do
    
    !___________________________________________________________________________
    ! subtract the the bolus velocities back from 3D velocities:
    if (Fer_GM) then
        UV    =UV    -fer_UV
        Wvel_e=Wvel_e-fer_Wvel
        Wvel  =Wvel  -fer_Wvel
    end if
    
    !___________________________________________________________________________
    ! to avoid crash with high salinities when coupled to atmosphere
    ! --> if we do only where (tr_arr(:,:,2) < 3._WP ) we also fill up the bottom 
    !     topogrpahy with values which are then writte into the output --> thats why
    !     do node=1,.... and tr_arr(node,1:nzmax,2)
    do node=1,myDim_nod2D+eDim_nod2D
        nzmax=nlevels_nod2D(node)-1
        where (tr_arr(1:nzmax,node,2) > 45._WP)
            tr_arr(1:nzmax,node,2)=45._WP
        end where

        where (tr_arr(1:nzmax,node,2) < 3._WP )
            tr_arr(1:nzmax,node,2)=3._WP
        end where
    end do
end subroutine solve_tracers_ale
!
!
!===============================================================================
subroutine adv_tracers_ale(tr_num, mesh)
    use g_config, only: flag_debug
    use g_parsup
    use mod_mesh
    use o_PARAM, only: tracer_adv
    use o_arrays
    use diagnostics, only: ldiag_DVD, compute_diag_dvd_2ndmoment_klingbeil_etal_2014, & 
                           compute_diag_dvd_2ndmoment_burchard_etal_2008, compute_diag_dvd
    use adv_tracers_muscle_ale_interface
    use adv_tracers_vert_ppm_ale_interface
    implicit none
    integer :: tr_num, node, nz
    type(t_mesh), intent(in) , target :: mesh    
    ! del_ttf ... initialised and setted to zero in call init_tracers_AB(tr_num)
    ! --> del_ttf ... equivalent to R_T^n in Danilov etal FESOM2: "from finite element
    !     to finite volume". At the end R_T^n should contain all advection therms and 
    !     the terms due to diffusion.
    ! del_ttf=0d0
    
    !___________________________________________________________________________
    ! if ldiag_DVD=.true. --> compute tracer second moments for the calcualtion 
    ! of discret variance decay
    if (ldiag_DVD .and. tr_num <= 2) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call compute_diag_dvd_2ndmoment'//achar(27)//'[0m'
        call compute_diag_dvd_2ndmoment_klingbeil_etal_2014(tr_num, mesh)
        !!PS call compute_diag_dvd_2ndmoment_burchard_etal_2008(tr_num)
    end if    
    
    !___________________________________________________________________________
    ! horizontal ale tracer advection 
    ! here --> add horizontal advection part to del_ttf(nz,n) = del_ttf(nz,n) + ...
    del_ttf_advhoriz = 0.0_WP
    del_ttf_advvert  = 0.0_WP
    select case (tracer_adv)
        case(1) !MUSCL
            ! --> tr_arr_old ... AB interpolated tracer from call init_tracers_AB(tr_num)
            if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call adv_tracers_muscle_ale'//achar(27)//'[0m'
            call adv_tracers_muscle_ale(tr_arr_old(:,:,tr_num), .25_WP, 1, mesh)
            !                                                      |    | 
            !             fraction of fourth-order contribution <--'    |
            !                              1st tracer moment is used <--'
            
            if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call adv_tracers_vert_ppm_ale'//achar(27)//'[0m'
            call adv_tracers_vert_ppm_ale(tr_arr(:,:,tr_num), 1, mesh)
            !                                                 | 
            !                    1st tracer moment is used <--'
            
        case(2) !MUSCL+FCT(3D)
            if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call adv_tracer_fct_ale'//achar(27)//'[0m'
            call adv_tracer_fct_ale(tr_arr_old(:,:,tr_num),tr_arr(:,:,tr_num), 1.0_WP, 1, mesh)
            !                                                                     |    | 
            !                            fraction of fourth-order contribution <--'    | 
            !                                             1st tracer moment is used <--'
            
        case default !unknown
            if (mype==0) write(*,*) 'Unknown ALE advection type. Check your namelists.'
            call par_ex(1)
    end select
    
    !___________________________________________________________________________
    ! update array for total tracer flux del_ttf with the fluxes from horizontal
    ! and vertical advection
    del_ttf=del_ttf+del_ttf_advhoriz+del_ttf_advvert
    
    !___________________________________________________________________________
    ! compute discrete variance decay after Burchard and Rennau 2008
    if (ldiag_DVD .and. tr_num <= 2) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call compute_diag_dvd'//achar(27)//'[0m'
        call compute_diag_dvd(tr_num, mesh)
    end if     
    
end subroutine adv_tracers_ale
!
!
!===============================================================================
! Horizontal ALE advection based on the gradient reconstruction.
! The last argument 'num_ord<1' defines the share of the
! 4th order centered contribution, and (1-num_ord) is done with 3rd order upwind.
! Dissipation comes only from the first part. num_ord=0.75--0.85 is 
! recommended if stable. 
! It is assumed that velocity is at n+1/2, where n is time step, hence only tracer field 
! is AB2 interpolated to n+1/2. 
! ttfAB --> corresponds to array tr_arr_old(:,:,tr_num) which is created by routine 
!             call init_tracers_AB(tr_num)
!             tr_arr_old(:,:,tr_num)=-(0.5+epsilon)*tr_arr_old(:,:,tr_num)+(1.5+epsilon)*tr_arr(:,:,tr_num)
subroutine adv_tracers_muscle_ale(ttfAB, num_ord, do_Xmoment, mesh)
    use MOD_MESH
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_CONFIG
    use g_comm_auto
    implicit none
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: el(2), enodes(2), n, nz, ed
    integer                  :: nl1, nl2, n2
    integer                  :: do_Xmoment !--> = [1,2] compute 1st or 2nd moment of tracer transport
    real(kind=WP)            :: c1, deltaX1, deltaY1, deltaX2, deltaY2, vflux=0.0_WP
    real(kind=WP)            :: c_lo(2)
    real(kind=WP)            :: Tmean1, Tmean2, a
    real(kind=WP)            :: ttfAB(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP)            :: num_ord

#include "associate_mesh.h"

    !___________________________________________________________________________
    ! Horizontal advection
    ! loop over loval edges 
    do ed=1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        enodes=edges(:,ed)   
        
        ! local index of element that contribute to edge
        el=edge_tri(:,ed)
        
        ! number of layers -1 at elem el(1)
        nl1=nlevels(el(1))-1
        
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,ed)
        deltaY1=edge_cross_dxdy(2,ed)
        a=r_earth*elem_cos(el(1))
        
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        deltaX2=0.0_WP
        deltaY2=0.0_WP
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,ed)
            deltaY2=edge_cross_dxdy(4,ed)
            ! number of layers -1 at elem el(2)
            nl2=nlevels(el(2))-1
            a=0.5_WP*(a+r_earth*elem_cos(el(2)))
        end if     
        
        ! n2 ... minimum number of layers -1 between element el(1) & el(2) that 
        ! contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        !                     that means nl1>0, nl2==0, n2=min(nl1,nl2)=0 !!!
        n2=min(nl1,nl2)
        
        !_______________________________________________________________________
        ! Both segments
        ! loop over depth layers from top to n2
        ! be carefull !!! --> if ed is a boundary edge, el(2)==0 than n2=0 so 
        !                     you wont enter in this loop
        do nz=1, n2
            !___________________________________________________________________
            ! MUSCL-type reconstruction
            ! check if upwind or downwind triagle is necessary
            !
            ! cross product between velocity vector and cross vector edge-elem-center
            !        o
            !       / \
            !      /   \
            !     /     \
            !    /   x   \
            !   /    |____\_____(dx,dy)
            !  o-----v-----o
            !  1   edge    2
            ! cross product > 0 --> angle vec_v and (dx,dy) --> [0   180] --> upwind triangle
            ! cross product < 0 --> angle vec_v and (dx,dy) --> [180 360] --> downwind triangle 
            !
            !       o                  o      !     o                  o
            !      / \                / \     !    / \                / \
            !     /   \    \ vec_v   /   \    !   /   \        /     /   \
            !    /  up \    \       / dn  \   !  /  up \      /     / dn  \
            !   o-------o----+---->o-------o  ! o-------o----+---->o-------o
            !           1   /      2          !         1     \vec_v
            !              /vec_v             !                \
            !   --> downwind triangle         ! --> upwind triangle 
            !  
            !  edge_up_dn_grad(1,nz,edge) ... gradTR_x upwind
            !  edge_up_dn_grad(2,nz,edge) ... gradTR_x downwind
            !  edge_up_dn_grad(3,nz,edge) ... gradTR_y upwind
            !  edge_up_dn_grad(4,nz,edge) ... gradTR_y downwind
            
            !___________________________________________________________________
            ! in case there is either no upwind or downwind triangle, than 
            ! calculate only low order solution 
            ! no upwind   triangle --> c_lo(1)=0, otherwise = 1
            ! no downwind triangle --> c_lo(2)=0, otherwise = 1
            ! (real(...) --> convert from integer to float)
            c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
            c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
            
            !___________________________________________________________________
            ! use downwind triangle to interpolate Tracer to edge center with 
            ! fancy scheme --> Linear upwind reconstruction
            ! T_n+0.5 = T_n+1 - 1/2*deltax*GRADIENT
            ! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
            ! T_n+0.5 = T_n+1 - 2/6*(T_n+1-T_n) + 1/6*gradT_down
            ! --> edge_up_dn_grad ... contains already elemental tracer gradient 
            !     of up and dn wind triangle
            ! --> Tmean2 ... edge center interpolated Tracer using tracer
            !     gradient info from upwind triangle
            Tmean2=ttfAB(nz, enodes(2))- &
                    (2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
                     edge_dxdy(1,ed)*a*edge_up_dn_grad(2,nz,ed)+ &
                     edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(4,nz,ed) &
                    )/6.0_WP*c_lo(2)
            
            ! use upwind triangle to interpolate Tracer to edge center with 
            ! fancy scheme --> Linear upwind reconstruction
            ! T_n+0.5 = T_n + 1/2*deltax*GRADIENT
            ! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
            ! T_n+0.5 = T_n + 2/6*(T_n+1-T_n) + 1/6*gradT_down
            ! --> Tmean1 ... edge center interpolated Tracer using tracer
            !     gradient info from downwind triangle
            Tmean1=ttfAB(nz, enodes(1))+ &
                    (2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
                     edge_dxdy(1,ed)*a*edge_up_dn_grad(1,nz,ed)+ &
                     edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(3,nz,ed) &
                    )/6.0_WP*c_lo(1)
                    
            !___________________________________________________________________
            ! volume flux along the edge segment ed
            ! netto volume flux along segment that comes from edge node 1 and 2
            !
            !                         
            !                         C1 (centroid el(1)) --> (u1,v1)
            !                         x
            !                         ^ 
            !               (dx1,dy1) |       
            !                         |---> vec_n1 (dy1,-dx1)--> project vec_u1 onto vec_n1 --> -v1*dx1+u1*dy1 -->
            !                         |                                                                          |
            !    enodes(1) o----------O---------o enodes(2)                                                      |-> calculate volume flux out of/in
            !          vflux_________/|                                                                          |   the volume of enode1(enode2) through
            !                     <---| vec_n2 (dy2,-dx2)--> project vec_u2 onto vec_n2 --> -v2*dx2+u2*dy2     -->   sections of dx1,dy1 and dx2,dy2
            !               (dx2,dy2) |                                                                              --> vflux 
            !                         v
            !                         x
            !                         C2 (centroid el(2)) --> (u2,v2)   
            
            ! here already assumed that ed is NOT! a boundary edge so el(2) should exist
            vflux=(-UV(2,nz,el(1))*deltaX1 + UV(1,nz,el(1))*deltaY1)*helem(nz,el(1)) &
                  +(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
            
            ! WHY vflux = (...vflux(el1)...) !!!-!!! (vflux(...el(2))...) ???
            ! --> Right-Hand-Rule: vectors vec_n1 and vec_n2 face in oppostite
            !     direction, the minus sign switches the direction of vec_n2
            !___________________________________________________________________
            ! tracer flux upwind
            ! if vflux (+) --> c1 = 2*vflux*Tmean1
            ! if vflux (-) --> c1 = -2*vflux*Tmean2
            ! so only use upwind tracer flux !!!!!!!!!!!!!
            c1=(vflux+abs(vflux))*( Tmean1**do_Xmoment )+ &
               (vflux-abs(vflux))*( Tmean2**do_Xmoment )
            !                                   ||
            !                                  _||_
            !                                  \  / 
            !                                   \/
            ! for the calculation of the discrete variance decay (DVD) diagnostic of 
            ! Klingbeil et al, 2014, Quantification of spurious dissipation and mixing 
            ! - Discrete variance decay in a Finite-Volume framework
            ! --> need the second moments of the tracers times flux at the control 
            !     volume interface
            ! --> if : do_Xmoment==1 --> 1st tracer moment  
            ! --> if : do_Xmoment==2 --> 2nd tracer moment  
            
            !___________________________________________________________________
            ! combined with centered
            ! num_ord is the fraction of fourth-order contribution in the HO solution
            ! (1-num_ord) is done with 3rd order upwind
!!PS             c1=-0.5_WP*( (1.0_WP-num_ord)*c1+vflux*num_ord*( (Tmean1+Tmean2)**do_Xmoment ) )
            c1=-0.5_WP*(1.0_WP-num_ord)*c1 - vflux*num_ord*( (0.5_WP*(Tmean1+Tmean2))**do_Xmoment ) 
            !                                                         |____________|
            !                                                            v
            !                                           dont use fourth order solution
            !                                           if its at the boundary
            
            !___________________________________________________________________
            ! write horizontal ale advection into rhs
            del_ttf_advhoriz(nz,enodes(1))=del_ttf_advhoriz(nz,enodes(1))+c1*dt/area(nz,enodes(1))
            del_ttf_advhoriz(nz,enodes(2))=del_ttf_advhoriz(nz,enodes(2))-c1*dt/area(nz,enodes(2))  
            
        end do ! --> do nz=1, n2
        
        !_______________________________________________________________________
        ! remaining segments on the left or on the right
        if(nl1>nl2) then 
            ! be carefull !!! --> if ed is a boundary edge, el(2)==0 than nl1>0 
            !                     and nl2==0, n2=0, so for boundary edges you will 
            !                     skip the previouse do loop and always end up 
            !                     in this part of the if condition
            do nz=1+n2,nl1
                !_______________________________________________________________
                ! check if upwind or downwind triangle exist, decide if high or 
                ! low order solution is calculated c_lo=1 --> high order, c_lo=0-->low order
                c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
                c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
                
                Tmean2=ttfAB(nz, enodes(2))- &
                        (2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
                        edge_dxdy(1,ed)*a*edge_up_dn_grad(2,nz,ed)+ &
                        edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(4,nz,ed))/6.0_WP*c_lo(2)
                
                Tmean1=ttfAB(nz, enodes(1))+ &
                        (2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
                        edge_dxdy(1,ed)*a*edge_up_dn_grad(1,nz,ed)+ &
                        edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(3,nz,ed))/6.0_WP*c_lo(1)
                        
                !_______________________________________________________________
                ! volume flux across the segments
                vflux=(-UV(2,nz,el(1))*deltaX1 + UV(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
                
                !_______________________________________________________________
                ! tracer flux upwind
                c1=(vflux+abs(vflux))*( Tmean1**do_Xmoment ) + &
                   (vflux-abs(vflux))*( Tmean2**do_Xmoment )
                
                !_______________________________________________________________
                ! combined with centered
                ! num_ord is the fraction of fourth-order contribution in the HO solution
                ! (1-num_ord) is done with 3rd order upwind
!!PS                 c1=-0.5_WP*( (1.0_WP-num_ord)*c1+vflux*num_ord*( (Tmean1+Tmean2)**do_Xmoment ) )
                c1=-0.5_WP*(1.0_WP-num_ord)*c1 - vflux*num_ord*( (0.5_WP*(Tmean1+Tmean2))**do_Xmoment )
                !                                                         |____________|
                !                                                            v
                !                                           dont use fourth order solution
                !                                           if its at the boundary
                !_______________________________________________________________
                ! write horizontal ale advection into rhs
                del_ttf_advhoriz(nz,enodes(1))=del_ttf_advhoriz(nz,enodes(1))+c1*dt/area(nz,enodes(1))
                del_ttf_advhoriz(nz,enodes(2))=del_ttf_advhoriz(nz,enodes(2))-c1*dt/area(nz,enodes(2)) 
                
            end do ! --> do nz=1+n2,nl1
        else
            do nz=n2+1,nl2
                !_______________________________________________________________
                ! check if upwind or downwind triangle exist, decide if high or 
                ! low order solution is calculated c_lo=1 --> high order, c_lo=0-->low order
                c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
                c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
                
                Tmean2=ttfAB(nz, enodes(2))- &
                        (2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
                        edge_dxdy(1,ed)*a*edge_up_dn_grad(2,nz,ed)+ &
                        edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(4,nz,ed))/6.0_WP*c_lo(2)
                
                Tmean1=ttfAB(nz, enodes(1))+ &
                        (2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
                        edge_dxdy(1,ed)*a*edge_up_dn_grad(1,nz,ed)+ &
                        edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(3,nz,ed))/6.0_WP*c_lo(1)
                        
                !_______________________________________________________________
                ! volume flux across the segments
                vflux=(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
                
                !_______________________________________________________________
                ! tracer flux upwind
                c1=(vflux+abs(vflux))*( Tmean1**do_Xmoment )+&
                   (vflux-abs(vflux))*( Tmean2**do_Xmoment )
                
                !_______________________________________________________________
                ! combined with centered
                ! num_ord is the fraction of fourth-order contribution in the HO solution
                ! (1-num_ord) is done with 3rd order upwind
!!PS                 c1=-0.5_WP*( (1.0_WP-num_ord)*c1+vflux*num_ord*( (Tmean1+Tmean2)**do_Xmoment ) )
                c1=-0.5_WP*(1.0_WP-num_ord)*c1 - vflux*num_ord*( (0.5_WP*(Tmean1+Tmean2))**do_Xmoment )
                !                                                        |_____________|
                !                                                           v
                !                                           dont use fourth order solution
                !                                           if its at the boundary
                !_______________________________________________________________
                ! write horizontal ale advection into rhs
                del_ttf_advhoriz(nz,enodes(1))=del_ttf_advhoriz(nz,enodes(1))+c1*dt/area(nz,enodes(1))
                del_ttf_advhoriz(nz,enodes(2))=del_ttf_advhoriz(nz,enodes(2))-c1*dt/area(nz,enodes(2))  
                
            end do ! --> do nz=n2+1,nl2
        end if ! --> if(nl1>nl2) then
    end do ! --> do ed=1, myDim_edge2D
end subroutine adv_tracers_muscle_ale
!
!
!===============================================================================
! Vertical ALE advection with PPM reconstruction (5th order)
subroutine adv_tracers_vert_ppm_ale(ttf, do_Xmoment, mesh)
    use g_config
    use MOD_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_forcing_arrays
    implicit none
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: n, nz, nzmax
    integer                  :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP)            :: tvert(mesh%nl), tv(mesh%nl), aL, aR, aj, x
    real(kind=WP)            :: dzjm1, dzj, dzjp1, dzjp2, deltaj, deltajp1
    real(kind=WP)            :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    integer                  :: overshoot_counter, counter

#include "associate_mesh.h"

    ! --------------------------------------------------------------------------
    ! Vertical advection
    ! --------------------------------------------------------------------------
    ! A piecewise parabolic scheme for uniformly-spaced layers.
    ! See Colella and Woodward, JCP, 1984, 174-201. It can be coded so as to to take 
    ! non-uniformity into account, but this is more cumbersome. This is the version for AB
    ! time stepping
    ! --------------------------------------------------------------------------
    overshoot_counter=0
    counter=0
    do n=1, myDim_nod2D
        !_______________________________________________________________________
        !Interpolate to zbar...depth levels --> all quantities (tracer ...) are 
        ! calculated on mid depth levels 
        ! nzmax ... number of depth levels at node n
        nzmax=nlevels_nod2D(n)
        ! tracer at surface layer
        tv(1)=ttf(1,n)
        ! tracer at surface+1 layer
        !   tv(2)=-ttf(1,n)*min(sign(1.0, Wvel_e(2,n)), 0._WP)+ttf(2,n)*max(sign(1.0, Wvel_e(2,n)), 0._WP)
        tv(2)=0.5*(ttf(1,n)+ttf(2,n))
        ! tacer at bottom-1 layer
        !   tv(nzmax-1)=-ttf(nzmax-2,n)*min(sign(1.0, Wvel_e(nzmax-1,n)), 0._WP)+ttf(nzmax-1,n)*max(sign(1.0, Wvel_e(nzmax-1,n)), 0._WP)
        tv(nzmax-1)=0.5_WP*(ttf(nzmax-2,n)+ttf(nzmax-1,n))
        ! tracer at bottom layer
        tv(nzmax)=ttf(nzmax-1,n)
        
        !_______________________________________________________________________
        ! calc tracer for surface+2 until depth-2 layer
        ! see Colella and Woodward, JCP, 1984, 174-201 --> equation (1.9)
        ! loop over layers (segments)
         do nz=2, nzmax-3
            !___________________________________________________________________
            ! for uniform spaced vertical grids --> piecewise parabolic method (ppm)
            ! equation (1.9)
            ! tv(nz)=(7.0_WP*(ttf(nz-1,n)+ttf(nz,n))-(ttf(nz-2,n)+ttf(nz+1,n)))/12.0_WP
            
            !___________________________________________________________________
            ! for non-uniformity spaced vertical grids --> piecewise parabolic 
            ! method (ppm) see see Colella and Woodward, JCP, 1984, 174-201 
            ! --> full equation (1.6), (1.7) and (1.8)
            dzjm1    = hnode_new(nz-2,n)
            dzj      = hnode_new(nz-1,n)
            dzjp1    = hnode_new(nz,n)
            dzjp2    = hnode_new(nz+1,n)
            ! Be carefull here vertical operation have to be done on NEW vertical mesh !!!
            
            !___________________________________________________________________
            ! equation (1.7)
            ! --> Here deltaj is the average slope in the jth zone of the parabola 
            !     with zone averages a_(j-1) and a_j, a_(j+1)
            ! --> a_j^n
            deltaj   = dzj/(dzjm1+dzj+dzjp1)* &
                      ( &
                       (2._WP*dzjm1+dzj    )/(dzjp1+dzj)*(ttf(nz+1,n)-ttf(nz  ,n)) +  &
                       (dzj    +2._WP*dzjp1)/(dzjm1+dzj)*(ttf(nz  ,n)-ttf(nz-1,n)) &
                      )
            ! --> a_(j+1)^n          
            deltajp1 = dzjp1/(dzj+dzjp1+dzjp2)* &
                      ( &
                       (2._WP*dzj+dzjp1  )/(dzjp2+dzjp1)*(ttf(nz+2,n)-ttf(nz+1,n)) +  &
                       (dzjp1+2._WP*dzjp2)/(dzj  +dzjp1)*(ttf(nz+1,n)-ttf(nz  ,n)) &
                      )
            !___________________________________________________________________
            ! condition (1.8)
            ! --> This modification leads to a somewhat steeper representation of 
            !     discontinuities in the solution. It also guarantees that a_(j+0.5)
            !     lies in the range of values defined by a_j; and a_(j+1);
            if ( (ttf(nz+1,n)-ttf(nz  ,n))*(ttf(nz  ,n)-ttf(nz-1,n)) > 0._WP ) then
                deltaj = min(  abs(deltaj), &
                             2._WP*abs(ttf(nz+1,n)-ttf(nz  ,n)),&
                             2._WP*abs(ttf(nz  ,n)-ttf(nz-1,n)) &
                             )*sign(1.0_WP,deltaj)
            else
                deltaj = 0.0_WP
            endif
            if ( (ttf(nz+2,n)-ttf(nz+1,n))*(ttf(nz+1,n)-ttf(nz  ,n)) > 0._WP ) then
                deltajp1 = min(  abs(deltajp1),&
                               2._WP*abs(ttf(nz+2,n)-ttf(nz+1,n)),&
                               2._WP*abs(ttf(nz+1,n)-ttf(nz,n)) &
                               )*sign(1.0_WP,deltajp1)
            else
                deltajp1 = 0.0_WP
            endif
            !___________________________________________________________________
            ! equation (1.6)
            ! --> calcualte a_(j+0.5)
            ! nz+1 is the interface betweel layers (segments) nz and nz+1
            tv(nz+1)=    ttf(nz,n) &
                        + dzj/(dzj+dzjp1)*(ttf(nz+1,n)-ttf(nz,n)) &
                        + 1._WP/(dzjm1+dzj+dzjp1+dzjp2) * &
                        ( &
                            (2._WP*dzjp1*dzj)/(dzj+dzjp1)* &
                                ((dzjm1+dzj)/(2._WP*dzj+dzjp1) - (dzjp2+dzjp1)/(2._WP*dzjp1+dzj))*(ttf(nz+1,n)-ttf(nz,n)) &
                        - dzj*(dzjm1+dzj)/(2._WP*dzj+dzjp1)*deltajp1 &
                        + dzjp1*(dzjp1+dzjp2)/(dzj+2._WP*dzjp1)*deltaj &
                        )
                       !tv(nz+1)=max(min(ttf(nz, n), ttf(nz+1, n)), min(max(ttf(nz, n), ttf(nz+1, n)), tv(nz+1)))
        end do ! --> do nz=2,nzmax-3
        
        tvert(1:nzmax)=0._WP
        ! loop over layers (segments)
        do nz=1, nzmax-1
            if ((Wvel_e(nz,n)<=0._WP) .AND. (Wvel_e(nz+1,n)>=0._WP)) CYCLE
            counter=counter+1
            aL=tv(nz)
            aR=tv(nz+1)
            if ((aR-ttf(nz, n))*(ttf(nz, n)-aL)<=0._WP) then
                !   write(*,*) aL, ttf(nz, n), aR
                overshoot_counter=overshoot_counter+1
                aL =ttf(nz, n)
                aR =ttf(nz, n)
            end if
            if ((aR-aL)*(ttf(nz, n)-0.5_WP*(aL+aR))> (aR-aL)**2/6._WP) then
                aL =3._WP*ttf(nz, n)-2._WP*aR
            end if
            if ((aR-aL)*(ttf(nz, n)-0.5_WP*(aR+aL))<-(aR-aL)**2/6._WP) then
                aR =3._WP*ttf(nz, n)-2._WP*aL
            end if
            
            dzj   = hnode(nz,n)
            aj=6.0_WP*(ttf(nz, n)-0.5_WP*(aL+aR))
            
            if (Wvel_e(nz,n)>0._WP) then
                x=min(Wvel_e(nz,n)*dt/dzj, 1._WP)
                tvert(nz  )=(-aL-0.5_WP*x*(aR-aL+(1._WP-2._WP/3._WP*x)*aj))
                tvert(nz  )=( tvert(nz)**do_Xmoment ) ! compute 2nd moment for DVD
                tvert(nz  )=tvert(nz)*area(nz,n)*Wvel_e(nz,n)
            end if
            
            if (Wvel_e(nz+1,n)<0._WP) then
                x=min(-Wvel_e(nz+1,n)*dt/dzj, 1._WP)
                tvert(nz+1)=(-aR+0.5_WP*x*(aR-aL-(1._WP-2._WP/3._WP*x)*aj))
                tvert(nz+1)=( tvert(nz+1)**do_Xmoment ) ! compute 2nd moment for DVD
                tvert(nz+1)=tvert(nz+1)*area(nz+1,n)*Wvel_e(nz+1,n)
            end if
        end do
        
        !_______________________________________________________________________
        ! Surface flux
        tvert(1)= -( tv(1)**do_Xmoment )*Wvel_e(1,n)*area(1,n)
        ! Zero bottom flux
        tvert(nzmax)=0.0_WP        
        
        !_______________________________________________________________________
        ! writing vertical ale advection into rhs
        do nz=1, nzmax-1
            ! no division over thickness in ALE !!!
            del_ttf_advvert(nz,n)=del_ttf_advvert(nz,n) + (tvert(nz)-tvert(nz+1))*dt/area(nz,n)
        end do
        
    end do ! --> do n=1, myDim_nod2D
!       if (mype==0) write(*,*) 'PPM overshoot statistics:', real(overshoot_counter)/real(counter)
end subroutine adv_tracers_vert_ppm_ale
!
!
!===============================================================================
! Vertical ALE advection with upwind reconstruction (1st order)
subroutine adv_tracers_vert_upw(ttf, do_Xmoment, mesh)
    use g_config
    use MOD_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_forcing_arrays
    implicit none
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: n, nz, nl1
    integer                  :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP)            :: tvert(mesh%nl), tv
    real(kind=WP)            :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)

#include "associate_mesh.h"

    ! --------------------------------------------------------------------------
    ! Vertical advection
    ! --------------------------------------------------------------------------
    do n=1, myDim_nod2D
        !_______________________________________________________________________
        nl1=nlevels_nod2D(n)-1
        !_______________________________________________________________________
        ! Surface flux
        tvert(1)= -Wvel_e(1,n)*(ttf(1,n)**do_Xmoment)*area(1,n)
        !_______________________________________________________________________
        ! Zero bottom flux
        tvert(nl1+1)=0.0_WP
        !_______________________________________________________________________
        ! Other levels
        do nz=2, nl1
            tv=(ttf(nz-1,n)**do_Xmoment)*min(Wvel_e(nz,n), 0._WP)+ &
               (ttf(nz  ,n)**do_Xmoment)*max(Wvel_e(nz,n), 0._WP)
            tvert(nz)= -tv*area(nz,n)
        end do
        !_______________________________________________________________________
        ! writing vertical ale advection into rhs
        do nz=1, nl1
            ! no division over thickness in ALE !!!
            del_ttf_advvert(nz,n)=del_ttf_advvert(nz,n) + (tvert(nz)-tvert(nz+1))*dt/area(nz,n) 
        end do         
    end do ! --> do n=1, myDim_nod2D
end subroutine adv_tracers_vert_upw
!
!
!===============================================================================
! Vertical ALE advection with central difference reconstruction (2nd order)
subroutine adv_tracers_vert_cdiff(ttf,do_Xmoment, mesh)
    use g_config
    use MOD_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_forcing_arrays
    implicit none
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: n, nz, nl1
    integer                  :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP)            :: tvert(mesh%nl), tv
    real(kind=WP)            :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    
#include "associate_mesh.h"

    ! --------------------------------------------------------------------------
    ! Vertical advection
    ! --------------------------------------------------------------------------
    do n=1, myDim_nod2D
        !_______________________________________________________________________
        nl1=nlevels_nod2D(n)-1
        !_______________________________________________________________________
        ! Surface flux
        tvert(1)= -Wvel_e(1,n)*(ttf(1,n)**do_Xmoment)*area(1,n)        
        !_______________________________________________________________________
        ! Zero bottom flux
        tvert(nl1+1)=0.0_WP        
        !_______________________________________________________________________
        ! Other levels
        do nz=2, nl1
            tv=0.5_WP*(ttf(nz-1,n)+ttf(nz,n))
            tv=tv**do_Xmoment
            tvert(nz)= -tv*Wvel_e(nz,n)*area(nz,n)
        end do
        !_______________________________________________________________________
        ! writing vertical ale advection into rhs
        do nz=1, nl1
            ! no division over thickness in ALE !!!
            del_ttf_advvert(nz,n)=del_ttf_advvert(nz,n) + (tvert(nz)-tvert(nz+1))*dt/area(nz,n) 
        end do         
    end do ! --> do n=1, myDim_nod2D
end subroutine adv_tracers_vert_cdiff
!
!
!===============================================================================
subroutine diff_tracers_ale(tr_num, mesh)
    use mod_mesh
    use g_PARSUP
    use o_arrays
    use o_tracers
    use diff_part_hor_redi_interface
    use diff_ver_part_expl_ale_interface
    use diff_ver_part_redi_expl_interface
    use diff_ver_part_impl_ale_interface
    implicit none
    
    integer, intent(in)      :: tr_num
    integer                  :: n, nzmax
    type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"
    !___________________________________________________________________________
    ! convert tr_arr_old(:,:,tr_num)=ttr_n-0.5   --> prepare to calc ttr_n+0.5
    ! eliminate AB (adams bashfort) interpolates tracer, which is only needed for 
    ! tracer advection. For diffusion only need tracer from previouse time step
    tr_arr_old(:,:,tr_num)=tr_arr(:,:,tr_num) !DS: check that this is the right place!
    
    !___________________________________________________________________________
    ! do horizontal diffusiion
    ! write there also horizontal diffusion rhs to del_ttf which is equal the R_T^n 
    ! in danilovs srcipt
    ! includes Redi diffusivity if Redi=.true.
    call diff_part_hor_redi(mesh) ! seems to be ~9% faster than diff_part_hor
    
    !___________________________________________________________________________
    ! do vertical diffusion: explicite 
    if (.not. i_vert_diff) call diff_ver_part_expl_ale(tr_num, mesh)
    ! A projection of horizontal Redi diffussivity onto vertical. This par contains horizontal
    ! derivatives and has to be computed explicitly!
    if (Redi) call diff_ver_part_redi_expl(mesh)
    
    !___________________________________________________________________________
    ! Update tracers --> calculate T* see Danilov etal "FESOM2 from finite elements
    ! to finite volume" 
    ! T* =  (dt*R_T^n + h^(n-0.5)*T^(n-0.5))/h^(n+0.5)
    do n=1, myDim_nod2D 
        nzmax=nlevels_nod2D(n)-1
        del_ttf(1:nzmax,n)=del_ttf(1:nzmax,n)+tr_arr(1:nzmax,n,tr_num)* &
                                    (hnode(1:nzmax,n)-hnode_new(1:nzmax,n))
        tr_arr(1:nzmax,n,tr_num)=tr_arr(1:nzmax,n,tr_num)+ &
                                    del_ttf(1:nzmax,n)/hnode_new(1:nzmax,n)
        ! WHY NOT ??? --> whats advantage of above --> tested it --> the upper 
        ! equation has a 30% smaller nummerical drift
        !tr_arr(1:nzmax,n,tr_num)=(hnode(1:nzmax,n)*tr_arr(1:nzmax,n,tr_num)+ &
        !                        del_ttf(1:nzmax,n))/hnode_new(1:nzmax,n)
        
    end do
    
    !___________________________________________________________________________
    if (i_vert_diff) then
        ! do vertical diffusion: implicite 
        call diff_ver_part_impl_ale(tr_num, mesh)
        
    end if
    
    !We DO not set del_ttf to zero because it will not be used in this timestep anymore
    !init_tracers will set it to zero for the next timestep
end subroutine diff_tracers_ale
!
!
!===============================================================================
!Vertical diffusive flux(explicit scheme):                                                                            
subroutine diff_ver_part_expl_ale(tr_num, mesh)
    use o_ARRAYS
    use g_forcing_arrays
    use MOD_MESH
    use g_PARSUP
    use g_config,only: dt
    
    implicit none 
    type(t_mesh), intent(in) , target :: mesh    
    real(kind=WP)            :: vd_flux(mesh%nl-1)
    real(kind=WP)            :: rdata,flux,rlx
    integer                  :: nz,nl1,tr_num,n
    real(kind=WP)            :: zinv1,Ty

#include "associate_mesh.h"
    !___________________________________________________________________________    
    do n=1, myDim_nod2D
        nl1=nlevels_nod2D(n)-1
        vd_flux=0._WP
        if (tr_num==1) then
            flux  = -heat_flux(n)/vcpw
            rdata =  Tsurf(n)
            rlx   =  surf_relax_T
        elseif (tr_num==2) then
            flux  =  virtual_salt(n)+relax_salt(n)- real_salt_flux(n)*is_nonlinfs
        else
            flux  = 0._WP
            rdata = 0._WP
            rlx=0._WP
        endif
        
        !_______________________________________________________________________
        !Surface forcing
        vd_flux(1)= flux
        
        !_______________________________________________________________________
        do nz=2,nl1
            !___________________________________________________________________
            zinv1=1.0_WP/(Z_3d_n(nz-1,n)-Z_3d_n(nz,n))
            
            !___________________________________________________________________
!            Ty= Kd(4,nz-1,n)*(Z_3d_n(nz-1,n)-zbar_3d_n(nz,n))*zinv1 *neutral_slope(3,nz-1,n)**2 + &
!                Kd(4,nz,n)*(zbar_3d_n(nz,n)-Z_3d_n(nz,n))*zinv1 *neutral_slope(3,nz,n)**2
            
            vd_flux(nz) = (Kv(nz,n)+Ty)*(tr_arr(nz-1,n,tr_num)-tr_arr(nz,n,tr_num))*zinv1*area(nz,n)
            
        end do
        
        !_______________________________________________________________________
        do nz=1,nl1-1
            del_ttf(nz,n) = del_ttf(nz,n) + (vd_flux(nz) - vd_flux(nz+1))/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))*dt/area(nz,n)
        end do
        del_ttf(nl1,n) = del_ttf(nl1,n) + (vd_flux(nl1)/(zbar_3d_n(nl1,n)-zbar_3d_n(nl1+1,n)))*dt/area(nl1,n)
        
    end do ! --> do n=1, myDim_nod2D
end subroutine diff_ver_part_expl_ale
!
!===============================================================================
! vertical diffusivity augmented with Redi contribution [vertical flux of K(3,3)*d_zT]
subroutine diff_ver_part_impl_ale(tr_num, mesh)
    use MOD_MESH
    use o_PARAM
    use o_ARRAYS
    use i_ARRAYS
    use g_PARSUP
    use g_CONFIG
    use g_forcing_arrays
        use o_mixing_KPP_mod !for ghats _GO_        
        
    implicit none
    type(t_mesh), intent(in) , target :: mesh
    real(kind=WP)            :: bc_surface    
    real(kind=WP)            :: a(mesh%nl), b(mesh%nl), c(mesh%nl), tr(mesh%nl)
    real(kind=WP)            :: cp(mesh%nl), tp(mesh%nl)
    integer                  :: nz, n, nzmax,tr_num
    real(kind=WP)            :: m, zinv, dt_inv, dz
    real(kind=WP)            :: rsss, Ty,Ty1, c1,zinv1,zinv2,v_adv
    real(kind=WP), external  :: TFrez  ! Sea water freeze temperature.
    real(kind=WP)            :: isredi=0._WP

#include "associate_mesh.h"

    !___________________________________________________________________________
    
    if (Redi) isredi=1._WP
    dt_inv=1.0_WP/dt
    Ty    =0.0_WP
    Ty1   =0.0_WP
    
    ! solve equation diffusion equation implicite part: 
    ! --> h^(n+0.5)* (T^(n+0.5)-Tstar) = dt*( K_33*d/dz*(T^(n+0.5)-Tstar) + K_33*d/dz*Tstar )
    ! -->   dTnew = T^(n+0.5)-Tstar
    ! -->   h^(n+0.5)* (dTnew) = dt*(K_33*d/dz*dTnew) + K_33*dt*d/dz*Tstar 
    ! -->   h^(n+0.5)* (dTnew) = dt*(K_33*d/dz*dTnew) + RHS 
    ! -->   solve for dT_new
    !    
    !    ----------- zbar_1, V_1 (Volume eq. to Area)
    ! Z_1 o T_1
    !    ----------- zbar_2, V_2
    ! Z_2 o T_2
    !    ----------- zbar_3, V_3
    ! Z_3 o T_3
    !    ----------- zbar_4
    !        :
    ! --> Difference Quotient at Volume _2:  ddTnew_2/dt + d/dz*K_33 d/dz*dTnew_2 = 0 --> homogene solution 
    ! V2*dTnew_2 *h^(n+0.5) = -dt * [ (dTnew_1-dTnew_2)/(Z_1-Z_2)*V_2 + (dTnew_2-dTnew_3)/(Z_2-Z_3)*V_3 ] + RHS
    !    dTnew_2 *h^(n+0.5) = -dt * [ (dTnew_1-dTnew_2)/(Z_1-Z_2)*V_2 + (dTnew_2-dTnew_3)/(Z_2-Z_3)*V_3/V_2 ] + RHS
    !                                                  |                                 |
    !                                                  v                                 v
    !                                         diffusive flux towards             diffusive flux towards
    !                                         T_2 trough boundary V2             T_2 trough boundary V3 
    !    
    ! --> solve coefficents for homogene part   
    !    dTnew_2 *h^(n+0.5) = -dt * [ a*dTnew_1 + b*dTnew_2 + c*dTnew_3 ] 
    !
    ! --> a = -dt*K_33/(Z_1-Z_2)
    ! 
    ! --> c = -dt*K_33/(Z_2-Z_3)*V_3/V_2
    !
    ! --> b = h^(n+0.5) -[ dt*K_33/(Z_1-Z_2) + dt*K_33/(Z_2-Z_3)*V_3/V_2 ] = -(a+c) + h^(n+0.5)
    
    !___________________________________________________________________________
    ! loop over local nodes
    do n=1,myDim_nod2D  
        
        ! initialise
        a  = 0.0_WP
        b  = 0.0_WP
        c  = 0.0_WP
        tr = 0.0_WP
        tp = 0.0_WP
        cp = 0.0_WP
        
        ! max. number of levels at node n
        nzmax=nlevels_nod2D(n)
        
        !___________________________________________________________________________
        ! Here can not exchange zbar_n & Z_n with zbar_3d_n & Z_3d_n because  
        ! they be calculate from the actualized mesh with hnode_new
        ! calculate new zbar (depth of layers) and Z (mid depths of layers) 
        ! depending on layer thinkness over depth at node n
        ! Be carefull here vertical operation have to be done on NEW vertical mesh !!!
        zbar_n=0.0_WP
        Z_n=0.0_WP
!         zbar_n(nzmax)=zbar(nzmax)
        zbar_n(nzmax)=zbar_n_bot(n)
        Z_n(nzmax-1)=zbar_n(nzmax) + hnode_new(nzmax-1,n)/2.0_WP
        do nz=nzmax-1,2,-1
            zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
            Z_n(nz-1)  = zbar_n(nz) + hnode_new(nz-1,n)/2.0_WP
        end do
        zbar_n(1) = zbar_n(2) + hnode_new(1,n)
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> surface layer 
        nz=1
        
        ! 1/dz(nz)
        zinv2=1.0_WP/(Z_n(nz)-Z_n(nz+1))
        zinv=1.0_WP*dt    ! no .../(zbar(1)-zbar(2)) because of  ALE
        
        ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
        Ty1= (Z_n(nz)     -zbar_n(nz+1))*zinv2 *slope_tapered(3,nz  ,n)**2*Ki(nz  ,n) + &
             (zbar_n(nz+1)-Z_n(   nz+1))*zinv2 *slope_tapered(3,nz+1,n)**2*Ki(nz+1,n)
        Ty1=Ty1*isredi
        ! layer dependent coefficients for for solving dT(1)/dt+d/dz*K_33*d/dz*T(1) = ...
        a(nz)=0.0_WP
        c(nz)=-(Kv(2,n)+Ty1)*zinv2*zinv*area(nz+1,n)/area(nz,n)
        b(nz)=-c(nz)+hnode_new(nz,n)      ! ale
        
        ! update from the vertical advection --> comes from splitting of vert 
        ! velocity into explicite and implicite contribution
        if (tracer_adv/=2) then
            v_adv=zinv*area(2,n)/area(1,n)
            b(1)=b(1)+Wvel_i(1, n)*zinv-min(0._WP, Wvel_i(2, n))*v_adv
            c(1)=c(1)-max(0._WP, Wvel_i(2, n))*v_adv
        end if        
        ! backup zinv2 for next depth level
        zinv1=zinv2
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> 2nd...nl-2 layer
        do nz=2, nzmax-2
        
            ! 1/dz(nz)
            zinv2=1.0_WP/(Z_n(nz)-Z_n(nz+1))
            ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
            Ty = (Z_n(nz-1   )-zbar_n(nz  ))*zinv1 *slope_tapered(3,nz-1,n)**2*Ki(nz-1,n)+ &
                 (zbar_n(nz  )-Z_n(nz     ))*zinv1 *slope_tapered(3,nz  ,n)**2*Ki(nz  ,n)
            Ty1= (Z_n(nz     )-zbar_n(nz+1))*zinv2 *slope_tapered(3,nz  ,n)**2*Ki(nz  ,n)+ &
                 (zbar_n(nz+1)-Z_n(nz+1   ))*zinv2 *slope_tapered(3,nz+1,n)**2*Ki(nz+1,n)
            Ty =Ty *isredi
            Ty1=Ty1*isredi
            ! layer dependent coefficients for for solving dT(nz)/dt+d/dz*K_33*d/dz*T(nz) = ...
            a(nz)=-(Kv(nz,n)  +Ty )*zinv1*zinv
            c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv*area(nz+1,n)/area(nz,n)
            b(nz)=-a(nz)-c(nz)+hnode_new(nz,n)
            
            ! backup zinv2 for next depth level
            zinv1=zinv2
            
            ! update from the vertical advection
            if (tracer_adv/=2) then
                v_adv=zinv
                a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv
                b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv
                v_adv=v_adv*area(nz+1,n)/area(nz,n)
                b(nz)=b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
                c(nz)=c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
            end if
        end do ! --> do nz=2, nzmax-2
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> nl-1 layer
        nz=nzmax-1 
        
        zinv=1.0_WP*dt   ! no ... /(zbar(nzmax-1)-zbar(nzmax)) because of ale
        
        ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
        Ty= (Z_n(nz-1)-zbar_n(nz))   *zinv1 *slope_tapered(3,nz-1,n)**2*Ki(nz-1,n) + &
            (zbar_n(nz)-Z_n(nz)) *zinv1 *slope_tapered(3,nz,n)**2  *Ki(nz,n)
        Ty =Ty *isredi
        ! layer dependent coefficients for for solving dT(nz)/dt+d/dz*K_33*d/dz*T(nz) = ...
        a(nz)=-(Kv(nz,n)+Ty)*zinv1*zinv
        c(nz)=0.0_WP
        b(nz)=-a(nz)+hnode_new(nz,n)
        
        ! update from the vertical advection
        if (tracer_adv/=2) then
            v_adv=zinv
            a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv       
            b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv
        end if
        
        !_______________________________________________________________________
        ! the rhs (inhomogene part): --> rhs = K_33*dt*d/dz*Tstar --> Tstar...tr_arr
        ! solve difference quotient for rhs --> tr
        !  RHS at Volume_2:
        !  
        !  RHS*V_2 = K_33*dt*(T_1-T_2)/(Z_1-Z_2)*V_2 - K_33*dt*(T_2-T_3)/(Z_2-Z_3)*V_3
        !          = -a*T_1 + (a+c)*T_2 - c*T_3
        !
        ! -+--> tr(1) =(a(1)+c(1))*tr_arr(1,n,tr_num)-c(1)*tr_arr(2,n,tr_num)
        !  |--> a(1)=0
        nz=1
        dz=hnode_new(nz,n) ! It would be (zbar(nz)-zbar(nz+1)) if not ALE
        tr(nz)=-(b(nz)-dz)*tr_arr(nz,n,tr_num)-c(nz)*tr_arr(nz+1,n,tr_num)
        !tr(nz)=c(nz)*(tr_arr(nz,n,tr_num) - tr_arr(nz+1,n,tr_num))
        
        
        ! *******************************************************************
        ! nonlocal transport to the rhs (only T and S currently) _GO_
        ! *******************************************************************
        ! rsss will be used later to compute:
        ! 1. the virtual salinity flux 
        ! 2. the contribution from the nonlocal term in KPP for salinity
        if (tr_num==2) then 
            rsss=ref_sss
                if (ref_sss_local) rsss=tr_arr(1,n,2)
        end if
        
        do nz=2,nzmax-2
            dz=hnode_new(nz,n)
            tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num)-(b(nz)-dz)*tr_arr(nz,n,tr_num)-c(nz)*tr_arr(nz+1,n,tr_num)
            !tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num) &
            !       -c(nz)*tr_arr(nz+1,n,tr_num) &
            !       +(a(nz)+c(nz))*tr_arr(nz,n,tr_num)
            
            ! *******************************************************************
            ! nonlocal transport to the rhs (only T and S currently) _GO_
            ! *******************************************************************
!leads to non conservation in 8th digit. needs to be checked!
!            if (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then
!                if (tr_num==1) then ! T
!                    tr(nz)=tr(nz)+(MIN(ghats(nz,n)*Kv(nz,n), 1.0_WP)-MIN(ghats(nz+1,n)*Kv(nz+1,n), 1.0_WP)*area(nz+1,n)/area(nz,n))*heat_flux(n)/vcpw
!                elseif (tr_num==2) then ! S
!                    tr(nz)=tr(nz)-(MIN(ghats(nz,n)*Kv(nz,n), 1.0_WP)-MIN(ghats(nz+1,n)*Kv(nz+1,n), 1.0_WP)*area(nz+1,n)/area(nz,n))*rsss*water_flux(n)
!                end if
!            end if 
        end do
        nz=nzmax-1
        dz=hnode_new(nz,n)
        tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num)-(b(nz)-dz)*tr_arr(nz,n,tr_num)
        !tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num)+a(nz)*tr_arr(nz,n,tr_num)
        
        !_______________________________________________________________________
        ! case of activated shortwave penetration into the ocean, ad 3d contribution
        if (use_sw_pene .and. tr_num==1) then
            do nz=1, nzmax-1
                zinv=1.0_WP*dt  !/(zbar(nz)-zbar(nz+1)) ale!
                tr(nz)=tr(nz)+(sw_3d(nz, n)-sw_3d(nz+1, n)*area(nz+1,n)/area(nz,n))*zinv
            end do
        end if
        
        !_______________________________________________________________________
        !  The first row contains also the boundary condition from heatflux, 
        !  freshwaterflux and relaxation terms
        !  --> tr_arr(1,n,1)*water_flux(n) : latent heatflux contribution due to 
        !      cell volume. If Volume decreases --> temp has to raise, if volume 
        !      expended --> temp has to decrease
        !                           (-)   ^                        (-)   ^ 
        !                            |    |                         |    | 
        !   IN MOMENT: heat_flux ~~~~|~~~~|~~~~   ,  water_flux ~~~~|~~~~|~~~~
        !  (BUT CHECK!)              |    |                         |    |
        !                            v   (+)                        v   (+) 
        !                            
        !
        tr(1)= tr(1)+bc_surface(n, tracer_id(tr_num))        
        
        !_______________________________________________________________________
        ! The forward sweep algorithm to solve the three-diagonal matrix 
        ! problem
        ! 
        !  | b_1 c_1 ...            |   |dTnew_1|
        !  | a_2 b_2 c_2 ...        |   |dTnew_2|
        !  |     a_3 b_3 c_3 ...    | * |dTnew_3| = RHS
        !  |         a_4 b_4 c_4 ...|   |dTnew_3| 
        !  |              :         |   |   :   |
        ! 
        ! 1st: define new coefficents:
        !      --> c'_i = c_i/b_i                               ; i=1
        !          c'_i = c_i/(b_i-a_i*c'_i-1)                  ; i = 2,3,...,n-1
        !      --> rhs'_i = rhs_i/b_i                           ; i=1
        !          rhs'_i = (rhs_i-a_i*d'_i-1)/(b_i-a_i*c'_i-1) ; i = 2,3,...,n-1
        !
        ! 2nd: solution is optained by back substitution
        !      --> dTnew_n = rhs'_n
        !      --> dTnew_i = rhs'_i-c'_i*dTnew_i+1 ; i = n-1,n-2,...,1
        !
        ! initialize c-prime and s,t-prime
        cp(1) = c(1)/b(1)
        tp(1) = tr(1)/b(1)
        
        ! solve for vectors c-prime and t, s-prime
        do nz = 2,nzmax-1
            m = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
        end do
        
        ! start with back substitution 
        tr(nzmax-1) = tp(nzmax-1)
        
        ! solve for x from the vectors c-prime and d-prime
        do nz = nzmax-2, 1, -1
            tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
        end do
        
        !_______________________________________________________________________
        ! update tracer
        ! tr ... dTnew = T^(n+0.5) - T*
        do nz=1,nzmax-1
            ! tr_arr - before ... T*
            tr_arr(nz,n,tr_num)=tr_arr(nz,n,tr_num)+tr(nz)
            ! tr_arr - after ... T^(n+0.5) = dTnew + T* = T^(n+0.5) - T* + T* 
            
        end do
        
    end do ! --> do n=1,myDim_nod2D   
end subroutine diff_ver_part_impl_ale
!
!
!===============================================================================
subroutine diff_ver_part_redi_expl(mesh)
    use o_ARRAYS
    use g_PARSUP
    use MOD_MESH
    USE o_param
    use g_config
    use g_comm_auto
    IMPLICIT NONE
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: elem,k
    integer                  :: n2,nl1,nl2,nz,n
    real(kind=WP)            :: Tx, Ty
    real(kind=WP)            :: tr_xynodes(2,mesh%nl-1,myDim_nod2D+eDim_nod2D), vd_flux(mesh%nl)

#include "associate_mesh.h"

    do n=1, myDim_nod2D
        nl1=nlevels_nod2D(n)-1
        do nz=1, nl1
            Tx=0.0_WP
            Ty=0.0_WP
            do k=1, nod_in_elem2D_num(n)
            elem=nod_in_elem2D(k,n)
            if(nz.LE.(nlevels(elem)-1)) then
                Tx=Tx+tr_xy(1,nz,elem)*elem_area(elem)
                Ty=Ty+tr_xy(2,nz,elem)*elem_area(elem)
            endif
        end do
        tr_xynodes(1,nz,n)=tx/3.0_WP/area(nz,n)
        tr_xynodes(2,nz,n)=ty/3.0_WP/area(nz,n)
        end do
    end do
    
    ! call exchange_nod_begin(tr_xynodes)  !NR the halo is not needed
    
    do n=1, myDim_nod2D
        nl1=nlevels_nod2D(n)-1
        vd_flux=0._WP
        
        !_______________________________________________________________________
        zbar_n=0.0_WP
        Z_n   =0.0_WP
!         zbar_n(nl1+1)=zbar(nl1+1)
        zbar_n(nl1+1)=zbar_n_bot(n)
        Z_n(nl1)=zbar_n(nl1+1) + hnode_new(nl1,n)/2.0_WP
        do nz=nl1, 2, -1
            zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
            Z_n(nz-1)  = zbar_n(nz) + hnode_new(nz-1,n)/2.0_WP
        end do
        zbar_n(1) = zbar_n(2) + hnode_new(1,n)
        
        !_______________________________________________________________________
        do nz=2,nl1
            vd_flux(nz)=(Z_n(nz-1)-zbar_n(nz))*(slope_tapered(1,nz-1,n)*tr_xynodes(1,nz-1,n)+slope_tapered(2,nz-1,n)*tr_xynodes(2,nz-1,n))*Ki(nz-1,n)
                        vd_flux(nz)=vd_flux(nz)+&
                          (zbar_n(nz)-Z_n(nz))  *(slope_tapered(1,nz,n)  *tr_xynodes(1,nz,n)  +slope_tapered(2,nz,n)     *tr_xynodes(2,nz,n))  *Ki(nz,n)
                        vd_flux(nz)=vd_flux(nz)/(Z_n(nz-1)-Z_n(nz))*area(nz,n)
        enddo
        do nz=1,nl1
            del_ttf(nz,n) = del_ttf(nz,n)+(vd_flux(nz) - vd_flux(nz+1))*dt/area(nz,n)
        enddo
    end do
end subroutine diff_ver_part_redi_expl!
!
!
!===============================================================================
subroutine diff_part_hor_redi(mesh)
    use o_ARRAYS
    use g_PARSUP
    use MOD_MESH
    use o_param
    use g_config
    IMPLICIT NONE
    type(t_mesh), intent(in) , target :: mesh
    real(kind=WP)            :: deltaX1,deltaY1,deltaX2,deltaY2
    integer                  :: edge
    integer                  :: n2,nl1,nl2,nz,el(2),elnodes(3),n,enodes(2)
    real(kind=WP)            :: c, Fx, Fy,Tx, Ty, Tx_z, Ty_z, SxTz, SyTz, Tz(2)
    real(kind=WP)            :: rhs1(mesh%nl-1), rhs2(mesh%nl-1), Kh, dz
    real(kind=WP)            :: isredi=0._WP

#include "associate_mesh.h"

    if (Redi) isredi=1._WP
    do edge=1, myDim_edge2D
        rhs1=0.0_WP
        rhs2=0.0_WP
        !_______________________________________________________________________
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        el=edge_tri(:,edge)
        enodes=edges(:,edge)
        nl1=nlevels(el(1))-1
        elnodes=elem2d_nodes(:,el(1))
        !Kh=elem_area(el(1))
        !_______________________________________________________________________
        nl2=0
        if (el(2)>0) then 
            !Kh=0.5_WP*(Kh+elem_area(el(2)))
            nl2=nlevels(el(2))-1
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
        endif
        !Kh=K_hor*Kh/scale_area
        !_______________________________________________________________________
        n2=min(nl1,nl2)
        do nz=1,n2
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=sum(helem(nz, el))/2.0_WP
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=0.5_WP*(tr_xy(1,nz,el(1))+tr_xy(1,nz,el(2)))
            Ty=0.5_WP*(tr_xy(2,nz,el(1))+tr_xy(2,nz,el(2)))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=((deltaX2-deltaX1)*Fy-(deltaY2-deltaY1)*Fx)*dz
            rhs1(nz) = rhs1(nz) + c
            rhs2(nz) = rhs2(nz) - c
        enddo
        
        !_______________________________________________________________________
        do nz=n2+1,nl1
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=helem(nz, el(1))
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=tr_xy(1,nz,el(1))
            Ty=tr_xy(2,nz,el(1))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=(-deltaX1*Fy+deltaY1*Fx)*dz
            rhs1(nz) = rhs1(nz) + c
            rhs2(nz) = rhs2(nz) - c
        end do
        do nz=n2+1,nl2
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=helem(nz, el(2))
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=tr_xy(1,nz,el(2))
            Ty=tr_xy(2,nz,el(2))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=(deltaX2*Fy-deltaY2*Fx)*dz
            rhs1(nz) = rhs1(nz) + c
            rhs2(nz) = rhs2(nz) - c
        end do
        
        !_______________________________________________________________________
        n2=max(nl1,nl2)
        del_ttf(1:n2,enodes(1))=del_ttf(1:n2,enodes(1))+rhs1(1:n2)*dt/area(1:n2,enodes(1))
        del_ttf(1:n2,enodes(2))=del_ttf(1:n2,enodes(2))+rhs2(1:n2)*dt/area(1:n2,enodes(2))
        
    end do
end subroutine diff_part_hor_redi
!
!
!===============================================================================
! this function returns a boundary conditions for a specified thacer ID and surface node
! ID = 0 and 1 are reserved for temperature and salinity
FUNCTION bc_surface(n, id)
  USE o_ARRAYS
  USE g_forcing_arrays
  USE g_PARSUP, only: mype, par_ex
  USE g_config
  implicit none

  REAL(kind=WP)     :: bc_surface
  integer           :: n, id
  character(len=10) :: id_string

  !  --> is_nonlinfs=1.0 for zelvel,zstar ....                            
  !  --> is_nonlinfs=0.0 for linfs
  SELECT CASE (id)
    CASE (0)
        bc_surface=-dt*(heat_flux(n)/vcpw + tr_arr(1,n,1)*water_flux(n)*is_nonlinfs)
    CASE (1)
        ! --> real_salt_flux(:): salt flux due to containment/releasing of salt
        !     by forming/melting of sea ice
        bc_surface= dt*(virtual_salt(n) & !--> is zeros for zlevel/zstar
                    + relax_salt(n) - real_salt_flux(n)*is_nonlinfs)
    CASE (101) ! apply boundary conditions to tracer ID=101
        bc_surface= dt*(prec_rain(n))! - real_salt_flux(n)*is_nonlinfs)
    CASE (301)
        bc_surface=0.0_WP
    CASE (302)
        bc_surface=0.0_WP
    CASE (303)
        bc_surface=0.0_WP
    CASE DEFAULT
      if (mype==0) then
         write (id_string, "(I3)") id
         if (mype==0) write(*,*) 'invalid ID '//trim(id_string)//' specified in boundary conditions'
         if (mype==0) write(*,*) 'the model will stop!'
      end if
      call par_ex
      stop
  END SELECT
  RETURN
END FUNCTION
