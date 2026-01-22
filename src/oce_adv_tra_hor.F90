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
! a not stable version of MUSCL (reconstruction in the vicinity of bottom topography is not upwind)
! it runs with FCT option only
    subroutine adv_tra_hor_mfct(ttf, vel, do_Xmoment, mesh, num_ord, flux, init_zero)
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
!
!
!===============================================================================
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
    integer                           :: nu12, nl12, nl1, nl2, nu1, nu2

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
        
        ! index off surface layer in case of cavity !=1
        nu1=ulevels(el(1))
        
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        a=r_earth*elem_cos(el(1)) 
        
        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        nu2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
            a=0.5_WP*(a+r_earth*elem_cos(el(2)))
        end if 
        
        !_______________________________________________________________________
        ! nl12 ... minimum number of layers -1 between element el(1) & el(2) that 
        ! contribute to edge ed
        ! nu12 ... upper index of layers between element el(1) & el(2) that 
        ! contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        !                     that means nl1>0, nl2==0, n2=min(nl1,nl2)=0 !!!
        nl12=min(nl1,nl2)
        nu12=max(nu1,nu2)        
        
        !_______________________________________________________________________
        ! (A) goes only into this loop when the edge has only facing element
        ! el(1) --> so the edge is a boundary edge --> this is for ocean 
        ! surface in case of cavity
        do nz=nu1, nu12-1              
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
           
           !____________________________________________________________________
           ! 1st. low order upwind solution
           flux(nz, edge)=-0.5_WP*(                                                     &
                         (ttf(nz, enodes(1))**do_Xmoment)*(vflux+abs(vflux))+ &
                         (ttf(nz, enodes(2))**do_Xmoment)*(vflux-abs(vflux))  &
                         )-flux(nz, edge)
        end do
        
        !_______________________________________________________________________
        ! (B) goes only into this loop when the edge has only facing elemenmt
        ! el(2) --> so the edge is a boundary edge --> this is for ocean 
        ! surface in case of cavity
        if (nu2 > 0) then 
            do nz=nu2, nu12-1
                !___________________________________________________________
                ! volume flux across the segments
                vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
                
                !___________________________________________________________
                ! 1st. low order upwind solution
                flux(nz, edge)=-0.5_WP*(                                           &
                            (ttf(nz, enodes(1))**do_Xmoment)*(vflux+abs(vflux))+ &
                            (ttf(nz, enodes(2))**do_Xmoment)*(vflux-abs(vflux)))-flux(nz, edge)
            end do
        end if     
        
        !_______________________________________________________________________
        ! (C) Both segments
        ! loop over depth layers from top (nu12) to nl12
        ! be carefull !!! --> if ed is a boundary edge, el(2)==0 than nl12=0 so 
        !                     you wont enter in this loop
        do nz=nu12, nl12
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
        ! (D) remaining segments on the left or on the right
        do nz=nl12+1, nl1              
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1))                 
           !____________________________________________________________________
           ! 1st. low order upwind solution
           flux(nz, edge)=-0.5_WP*(                                                     &
                         (ttf(nz, enodes(1))**do_Xmoment)*(vflux+abs(vflux))+ &
                         (ttf(nz, enodes(2))**do_Xmoment)*(vflux-abs(vflux))  &
                         )-flux(nz, edge)
        end do
        
        !_______________________________________________________________________
        ! (E) remaining segments on the left or on the right
        do nz=nl12+1, nl2
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
!
!
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
    integer                           :: nu12, nl12, nl1, nl2, nu1, nu2

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
        
        ! index off surface layer in case of cavity !=1
        nu1=ulevels(el(1))
        
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        a=r_earth*elem_cos(el(1))        
        
        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        nu2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
            a=0.5_WP*(a+r_earth*elem_cos(el(2)))
        end if 
        
        !_______________________________________________________________________
        ! n2 ... minimum number of layers -1 between element el(1) & el(2) that 
        ! contribute to edge ed
        ! nu12 ... upper index of layers between element el(1) & el(2) that 
        ! contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        !                     that means nl1>0, nl2==0, n2=min(nl1,nl2)=0 !!!
        nl12=min(nl1,nl2)        
        nu12=max(nu1,nu2)
        
        !_______________________________________________________________________
        ! (A) goes only into this loop when the edge has only facing element
        ! el(1) --> so the edge is a boundary edge --> this is for ocean 
        ! surface in case of cavity
        do nz=nu1, nu12-1
           c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
           c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
           
           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)
                
           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)
                  
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
           cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
        end do
        
        !_______________________________________________________________________
        ! (B) goes only into this loop when the edge has only facing elemenmt
        ! el(2) --> so the edge is a boundary edge --> this is for ocean 
        ! surface in case of cavity
        if (nu2 > 0) then 
            do nz=nu2, nu12-1
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
        end if 
        
        !_______________________________________________________________________
        ! (C) Both segments
        ! loop over depth layers from top to n2
        ! be carefull !!! --> if ed is a boundary edge, el(2)==0 than n2=0 so 
        !                     you wont enter in this loop
        do nz=nu12, nl12
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
        ! (D) remaining segments on the left or on the right
        do nz=nl12+1, nl1
           c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
           c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
           
           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)
                
           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)
                  
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
           cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
        end do
        
        !_______________________________________________________________________
        ! (E) remaining segments on the left or on the right
        do nz=nl12+1, nl2
           c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
           c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)
           
           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)
                
           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)  
                  
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
           cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
        end do
    end do
end subroutine adv_tra_hor_muscl
!
!
!===============================================================================
subroutine adv_tra_hor_mfct(ttf, vel, do_Xmoment, mesh, num_ord, flux, init_zero)
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
    real(kind=WP)                     :: a, vflux
    integer                           :: el(2), enodes(2), nz, edge
    integer                           :: nu12, nl12, nl1, nl2, nu1, nu2

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
        
        ! index off surface layer in case of cavity !=1
        nu1=ulevels(el(1))
        
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to 
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        a=r_earth*elem_cos(el(1))        
        
        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        nu2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
            a=0.5_WP*(a+r_earth*elem_cos(el(2)))
        end if 
        
        !_______________________________________________________________________
        ! n2 ... minimum number of layers -1 between element el(1) & el(2) that 
        ! contribute to edge ed
        ! nu12 ... upper index of layers between element el(1) & el(2) that 
        ! contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        !                     that means nl1>0, nl2==0, n2=min(nl1,nl2)=0 !!!
        nl12=min(nl1,nl2) 
        nu12=max(nu1,nu2) 
        
        !_______________________________________________________________________
        ! (A) goes only into this loop when the edge has only facing element
        ! el(1) --> so the edge is a boundary edge --> this is for ocean 
        ! surface in case of cavity
        do nz=nu1, nu12-1
           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP
                
           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
                  
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
           cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
        end do
        
        !_______________________________________________________________________
        ! (B) goes only into this loop when the edge has only facing elemenmt
        ! el(2) --> so the edge is a boundary edge --> this is for ocean 
        ! surface in case of cavity
        if (nu2 > 0) then 
            do nz=nu2,nu12-1
            !___________________________________________________________________
            Tmean2=ttf(nz, enodes(2))- &
                    (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                    edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                    edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP
                    
            Tmean1=ttf(nz, enodes(1))+ &
                    (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                    edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                    edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
            !___________________________________________________________________
            ! volume flux across the segments
            vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
            cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
            flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
            end do
        end if
        
        !_______________________________________________________________________
        ! (C) Both segments
        ! loop over depth layers from top to n2
        ! be carefull !!! --> if ed is a boundary edge, el(2)==0 than n2=0 so 
        !                     you wont enter in this loop
        do nz=nu12, nl12
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
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP
            
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
                   edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
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
        ! (D) remaining segments on the left or on the right
        do nz=nl12+1, nl1
           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP
                
           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
                  
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
           cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
        end do
        
        !_______________________________________________________________________
        ! (E) remaining segments on the left or on the right
        do nz=nl12+1, nl2
           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP
                
           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
                  
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
           cHO=(vflux+abs(vflux))*(Tmean1**do_Xmoment) + (vflux-abs(vflux))*(Tmean2**do_Xmoment)
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*( 0.5_WP*(Tmean1+Tmean2))**do_Xmoment-flux(nz,edge)
        end do
    end do
end subroutine adv_tra_hor_mfct

