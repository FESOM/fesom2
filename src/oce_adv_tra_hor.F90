!===============================================================================================================================
!**************** routines for horizontal tracer advection ***********************
module oce_adv_tra_hor_interfaces
  interface
! (low order upwind)
! returns flux given at edges which contributes with 
! plus sign into 1st. node and with the minus sign into the 2nd node
! IF init_zero=.TRUE.  : flux will be set to zero before computation
! IF init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_hor_upw1(ttf, vel, do_Xmoment, mesh, flux, init_zero)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh), intent(in) , target :: mesh
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, myDim_elem2D+eDim_elem2D)
      integer,       intent(in)         :: do_Xmoment
      real(kind=WP), intent(inout)      :: flux(mesh%nl-1, myDim_edge2D)
      logical, optional                 :: init_zero
    end subroutine
!===============================================================================
! MUSCL
! returns flux given at edges which contributes with 
! plus sign into 1st. node and with the minus sign into the 2nd node
! IF init_zero=.TRUE.  : flux will be set to zero before computation
! IF init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_hor_muscl(ttf, vel, do_Xmoment, mesh, num_ord, flux, init_zero)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh    
      integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
      real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1,  myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, myDim_elem2D+eDim_elem2D)
      real(kind=WP), intent(inout)      :: flux(mesh%nl-1, myDim_edge2D)
      logical, optional                 :: init_zero
    end subroutine
    subroutine adv_tra_hor_muscl_v2(ttf, vel, do_Xmoment, mesh, num_ord, flux, init_zero)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh    
      integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
      real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1,  myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, myDim_elem2D+eDim_elem2D)
      real(kind=WP), intent(inout)      :: flux(mesh%nl-1, myDim_edge2D)
      logical, optional                 :: init_zero
    end subroutine
  end interface
end module
!====================================================================
subroutine adv_tra_hor_upw1(ttf, vel, do_Xmoment, mesh, flux, init_zero)
    use MOD_MESH
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_CONFIG
    use g_comm_auto
    implicit none
    type(t_mesh), intent(in) , target :: mesh    
    integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, myDim_elem2D+eDim_elem2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl-1, myDim_edge2D)
    logical, optional                 :: init_zero
    real(kind=WP)                     :: deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=WP)                     :: a, vflux
    integer                           :: el(2), enodes(2), nz, edge
    integer                           :: n2, nl1, nl2

#include "associate_mesh.h"


    if (present(init_zero))then
       if (init_zero) flux=0.0_WP
    else
       flux=0.0_WP
    end if

    ! The result is the low-order solution horizontal fluxes
    ! They are put into flux
    !___________________________________________________________________________
    do edge=1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        enodes=edges(:,edge)          
        ! local index of element that contribute to edge
        el=edge_tri(:,edge)        
        ! number of layers -1 at elem el(1)
        nl1=nlevels(el(1))-1        
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        a=r_earth*elem_cos(el(1))        
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
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
            ! 1st. low order upwind solution
            ! here already assumed that ed is NOT! a boundary edge so el(2) should exist
            vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) &
                  +(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
            flux(nz, edge)=-0.5_WP*(                                                     &
                         (ttf(nz, enodes(1))**do_Xmoment)*(vflux+abs(vflux))+ &
                         (ttf(nz, enodes(2))**do_Xmoment)*(vflux-abs(vflux)))-flux(nz, edge)
        end do
        !_______________________________________________________________________
        ! remaining segments on the left or on the right
        do nz=n2+1, nl1              
           !_______________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1))                 
           !___________________________________________________________________
           ! 1st. low order upwind solution
           flux(nz, edge)=-0.5_WP*(                                                     &
                         (ttf(nz, enodes(1))**do_Xmoment)*(vflux+abs(vflux))+ &
                         (ttf(nz, enodes(2))**do_Xmoment)*(vflux-abs(vflux))  &
                         )-flux(nz, edge)
        end do
        do nz=n2+1, nl2
                !_______________________________________________________________
                ! volume flux across the segments
                vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
                !_______________________________________________________________
                ! 1st. low order upwind solution
                flux(nz, edge)=-0.5_WP*(                                           &
                             (ttf(nz, enodes(1))**do_Xmoment)*(vflux+abs(vflux))+ &
                             (ttf(nz, enodes(2))**do_Xmoment)*(vflux-abs(vflux)))-flux(nz, edge)
        end do
    end do
end subroutine adv_tra_hor_upw1
!===============================================================================
subroutine adv_tra_hor_muscl(ttf, vel, do_Xmoment, mesh, num_ord, flux, init_zero)
    use MOD_MESH
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_CONFIG
    use g_comm_auto
    implicit none
    type(t_mesh),  intent(in), target :: mesh    
    integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1,  myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, myDim_elem2D+eDim_elem2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl-1, myDim_edge2D)
    logical, optional                 :: init_zero
    real(kind=WP)                     :: deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=WP)                     :: Tmean1, Tmean2, cHO
    real(kind=WP)                     :: c_lo(2)
    real(kind=WP)                     :: a, vflux
    integer                           :: el(2), enodes(2), nz, edge
    integer                           :: n2, nl1, nl2

#include "associate_mesh.h"

    if (present(init_zero))then
       if (init_zero) flux=0.0_WP
    else
       flux=0.0_WP
    end if

    ! The result is the low-order solution horizontal fluxes
    ! They are put into flux
    !___________________________________________________________________________
    do edge=1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        enodes=edges(:,edge)          
        ! local index of element that contribute to edge
        el=edge_tri(:,edge)        
        ! number of layers -1 at elem el(1)
        nl1=nlevels(el(1))-1        
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        a=r_earth*elem_cos(el(1))        
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
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
            c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
            c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
           !___________________________________________________________________
           ! MUSCL-type reconstruction
           ! check if upwind or downwind triagle is necessary
           !
           ! cross product between velocity vector and cross vector edge-elem-center
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
           ! use downwind triangle to interpolate Tracer to edge center with 
           ! fancy scheme --> Linear upwind reconstruction
           ! T_n+0.5 = T_n+1 - 1/2*deltax*GRADIENT
           ! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
           ! T_n+0.5 = T_n+1 - 2/6*(T_n+1-T_n) + 1/6*gradT_down
           ! --> edge_up_dn_grad ... contains already elemental tracer gradient 
           !     of up and dn wind triangle
           ! --> Tmean2 ... edge center interpolated Tracer using tracer
           !     gradient info from upwind triangle
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)
            
            ! use upwind triangle to interpolate Tracer to edge center with 
            ! fancy scheme --> Linear upwind reconstruction
            ! T_n+0.5 = T_n + 1/2*deltax*GRADIENT
            ! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
            ! T_n+0.5 = T_n + 2/6*(T_n+1-T_n) + 1/6*gradT_down
            ! --> Tmean1 ... edge center interpolated Tracer using tracer
            !     gradient info from downwind triangle
            Tmean1=ttf(nz, enodes(1))+ &
                   (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                   edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                   edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)                 
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
            !                         |---> vec_n2 (dy2,-dx2)--> project vec_u2 onto vec_n2 --> -v2*dx2+u2*dy2 -->   sections of dx1,dy1 and dx2,dy2
            !               (dx2,dy2) |                                                                              --> vflux 
            !                         v
            !                         x
            !                         C2 (centroid el(2)) --> (u2,v2)   
            
            ! here already assumed that ed is NOT! a boundary edge so el(2) should exist
            vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) &
                  +(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
                        
            !___________________________________________________________________
            ! (1-num_ord) is done with 3rd order upwind
            cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
            flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
        end do
        !_______________________________________________________________________
        ! remaining segments on the left or on the right
        do nz=n2+1, nl1
           c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
           c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
           !_______________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)
                
           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)
           !_______________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
           cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
        end do
        do nz=n2+1, nl2
           c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
           c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
           !_______________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)
                
           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)  
           !_______________________________________________________________
           ! volume flux across the segments
           vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
           cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
        end do
    end do
end subroutine adv_tra_hor_muscl
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
subroutine adv_tra_hor_muscl_v2(ttf, vel, do_Xmoment, mesh, num_ord, flux, init_zero)
    use MOD_MESH
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_CONFIG
    use g_comm_auto
    implicit none
    type(t_mesh),  intent(in), target :: mesh    
    integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1,  myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, myDim_elem2D+eDim_elem2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl-1, myDim_edge2D)
    logical, optional                 :: init_zero
    real(kind=WP)                     :: deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=WP)                     :: Tmean1, Tmean2, cHO
    real(kind=WP)                     :: c_lo(2),c1
    real(kind=WP)                     :: a, vflux
    integer                           :: el(2), enodes(2), nz, ed
    integer                           :: n2, nl1, nl2
#include "associate_mesh.h"


    if (present(init_zero))then
       if (init_zero) flux=0.0_WP
    else
       flux=0.0_WP
    end if

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
            c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
            c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
            Tmean2=ttf(nz, enodes(2))- &
                    (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                     edge_dxdy(1,ed)*a*edge_up_dn_grad(2,nz,ed)+ &
                     edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(4,nz,ed) &
                    )/6.0_WP*c_lo(2)
            Tmean1=ttf(nz, enodes(1))+ &
                    (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                     edge_dxdy(1,ed)*a*edge_up_dn_grad(1,nz,ed)+ &
                     edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(3,nz,ed) &
                    )/6.0_WP*c_lo(1)
            vflux=(-UV(2,nz,el(1))*deltaX1 + UV(1,nz,el(1))*deltaY1)*helem(nz,el(1)) &
                  +(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
            c1=(vflux+abs(vflux))*( Tmean1**do_Xmoment )+ &
               (vflux-abs(vflux))*( Tmean2**do_Xmoment )
            c1=-0.5_WP*(1.0_WP-num_ord)*c1 - vflux*num_ord*( (0.5_WP*(Tmean1+Tmean2))**do_Xmoment ) 
            flux(nz,ed)=c1-flux(nz,ed)
        end do
        
        do nz=n2+1,nl1
                c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
                c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
                
                Tmean2=ttf(nz, enodes(2))- &
                        (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                        edge_dxdy(1,ed)*a*edge_up_dn_grad(2,nz,ed)+ &
                        edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(4,nz,ed))/6.0_WP*c_lo(2)
                
                Tmean1=ttf(nz, enodes(1))+ &
                        (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                        edge_dxdy(1,ed)*a*edge_up_dn_grad(1,nz,ed)+ &
                        edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(3,nz,ed))/6.0_WP*c_lo(1)
                        
                !_______________________________________________________________
                ! volume flux across the segments
                vflux=(-UV(2,nz,el(1))*deltaX1 + UV(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
                
                !_______________________________________________________________
                ! tracer flux upwind
                c1=(vflux+abs(vflux))*( Tmean1**do_Xmoment )+ &
                   (vflux-abs(vflux))*( Tmean2**do_Xmoment )                
                c1=-0.5_WP*(1.0_WP-num_ord)*c1 - vflux*num_ord*( (0.5_WP*(Tmean1+Tmean2))**do_Xmoment )
                flux(nz,ed)=c1-flux(nz,ed)
         end do ! --> do nz=1+n2,nl1
         do nz=n2+1,nl2
                !_______________________________________________________________
                ! check if upwind or downwind triangle exist, decide if high or 
                ! low order solution is calculated c_lo=1 --> high order, c_lo=0-->low order
                c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
                c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
                
                Tmean2=ttf(nz, enodes(2))- &
                        (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                        edge_dxdy(1,ed)*a*edge_up_dn_grad(2,nz,ed)+ &
                        edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(4,nz,ed))/6.0_WP*c_lo(2)
                
                Tmean1=ttf(nz, enodes(1))+ &
                        (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                        edge_dxdy(1,ed)*a*edge_up_dn_grad(1,nz,ed)+ &
                        edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(3,nz,ed))/6.0_WP*c_lo(1)
                        
                !_______________________________________________________________
                ! volume flux across the segments
                vflux=(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
                
                !_______________________________________________________________
                ! tracer flux upwind
                c1=(vflux+abs(vflux))*( Tmean1**do_Xmoment )+ &
                   (vflux-abs(vflux))*( Tmean2**do_Xmoment )                
                c1=-0.5_WP*(1.0_WP-num_ord)*c1 - vflux*num_ord*( (0.5_WP*(Tmean1+Tmean2))**do_Xmoment )
                flux(nz,ed)=c1-flux(nz,ed)                
        end do ! --> do nz=n2+1,nl2
    end do ! --> do ed=1, myDim_edge2D
end subroutine adv_tra_hor_muscl_v2

