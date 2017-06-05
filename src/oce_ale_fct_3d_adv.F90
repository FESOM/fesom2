!
!
!===============================================================================
! Caller routine for FCT tracer advection 
subroutine adv_tracer_fct_ale(ttfAB,ttf,num_ord)
	use o_ARRAYS
	use o_MESH
	use o_PARAM
	use g_PARSUP
	implicit none
	
	real(kind=WP)  :: ttfAB(nl-1, myDim_nod2D+eDim_nod2D),ttf(nl-1, myDim_nod2D+eDim_nod2D)
	! real(kind=WP)  :: ttfold(nl-1, myDim_nod2D+eDim_nod2D)
	real(kind=WP)  :: num_ord
	
	! 1st. first calculate Low and High order solution
	call fct_ale_muscl_LH(ttfAB,ttf,num_ord)

!just for using the low order:
!	fct_aec    =0.
!	fct_aec_ver=0.
	
	! 2nd. apply constrained bounds
	call fct_ale(ttf)
	
end subroutine adv_tracer_fct_ale
!
!
!===============================================================================
subroutine fct_ale_muscl_LH(ttfAB,ttf, num_ord)
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_CONFIG
	use g_comm_auto
	implicit none
	
	integer       :: el(2), enodes(2), n, nz, edge
	integer       :: n2, nl1, nl2,tr_num
	real(kind=WP) :: cLO, cHO, deltaX1, deltaY1, deltaX2, deltaY2
	real(kind=WP) :: qc, qu, qd
	real(kind=WP) :: tvert(nl), a, b, c, d, da, db, dg, vflux, Tupw1
	real(kind=WP) :: Tmean, Tmean1, Tmean2, num_ord
	real(kind=WP) :: ttfAB(nl-1,myDim_nod2D+eDim_nod2D), ttf(nl-1,myDim_nod2D+eDim_nod2D)
	
	! --------------------------------------------------------------------------
	! It is assumed that velocity is at n+1/2, hence only tracer field 
	! is AB2 interpolated to n+1/2. 
	! ttf contains tracer at step n
	! ttfAB is AB interpolated (n+1/2)
	! num_ord is the fraction of fourth-order contribution in the HO solution
	! The result is the low-order solution
	! and vertical and horizontal antidiffusive fluxes
	! They are put into fct_LO
	!                   fct_aec
	!                   fct_aec_ver
	! --------------------------------------------------------------------------
	
	!___________________________________________________________________________
	! Clean the low-order solution
	fct_LO=0.0_WP 
	
	!___________________________________________________________________________
	!**************** (A) do horizontal tracer advection ***********************
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
			! MUSCL-type reconstruction
			! check if upwind or downwind triagle is necessary
			!
			! cross product between velocity vector and cross vector edge-elem-center
			! cross product > 0 --> angle vec_v and (dx,dy) --> [0   180] --> upwind triangle
			! cross product < 0 --> angle vec_v and (dx,dy) --> [180 360] --> downwind triangle 
			!
			!	 	o                  o      !	 	o                  o
			!	   / \                / \     !	   / \                / \
			!	  /   \    \ vec_v   /   \    !	  /   \        /     /   \
			!	 /  up \    \       / dn  \   !	 /  up \      /     / dn  \
			!	o-------o----+---->o-------o  !	o-------o----+---->o-------o
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
			Tmean2=ttfAB(nz, enodes(2))- &
					(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
					edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
					edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP
			
			! use upwind triangle to interpolate Tracer to edge center with 
			! fancy scheme --> Linear upwind reconstruction
			! T_n+0.5 = T_n + 1/2*deltax*GRADIENT
			! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
			! T_n+0.5 = T_n + 2/6*(T_n+1-T_n) + 1/6*gradT_down
			! --> Tmean1 ... edge center interpolated Tracer using tracer
			!     gradient info from downwind triangle
			Tmean1=ttfAB(nz, enodes(1))+ &
					(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
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
			vflux=(-UV(2,nz,el(1))*deltaX1 + UV(1,nz,el(1))*deltaY1)*helem(nz,el(1)) &
				  +(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
			
			!___________________________________________________________________
			! 1st. Low order upwind solution
			cLO=-0.5_WP*(ttf(nz, enodes(1))*(vflux+abs(vflux))+ttf(nz, enodes(2))*(vflux-abs(vflux)))
			fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))+cLO     
			fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))-cLO  
			
			!___________________________________________________________________
			! 2nd. High order solution
			! num_ord is the fraction of fourth-order contribution in the HO solution
			! (1-num_ord) is done with 3rd order upwind
			cHO=(vflux+abs(vflux))*Tmean1+(vflux-abs(vflux))*Tmean2
			cHO=-0.5_WP*((1.0_WP-num_ord)*cHO+vflux*num_ord*(Tmean1+Tmean2))
			
			!___________________________________________________________________
			! 3nd. calculate Antidiffusive edge flux: AEF=[HO - LO]
			fct_aec(nz,edge)=cHO-cLO
			
		end do  ! --> do nz=1, n2
		
		!_______________________________________________________________________
		! remaining segments on the left or on the right
		if(nl1>nl2) then
			! be carefull !!! --> if ed is a boundary edge, el(2)==0 than nl1>0 
			!                     and nl2==0, n2=0, so for boundary edges you will 
			!                     skip the previouse do loop and always end up 
			!                     in this part of the if condition
			do nz=n2+1, nl1
				!_______________________________________________________________
				Tmean2=ttfAB(nz, enodes(2))- &
				(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
				edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
				edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP   
				
				Tmean1=ttfAB(nz, enodes(1))+ &
				(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
				edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
				edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
				
				!_______________________________________________________________
				! volume flux across the segments
				vflux=(-UV(2,nz,el(1))*deltaX1 + UV(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
				
				!___________________________________________________________________
				! 1st. Low order upwind solution
				cLO=-0.5_WP*(ttf(nz, enodes(1))*(vflux+abs(vflux))+ttf(nz, enodes(2))*(vflux-abs(vflux)))
				fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))+cLO
				fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))-cLO
				
				!_______________________________________________________________
				! 2nd. High order solution 
				! num_ord is the fraction of fourth-order contribution in the HO solution
				! (1-num_ord) is done with 3rd order upwind
				cHO=(vflux+abs(vflux))*Tmean1+(vflux-abs(vflux))*Tmean2
				cHO=-0.5_WP*((1.0_WP-num_ord)*cHO+vflux*num_ord*(Tmean1+Tmean2))
				
				!_______________________________________________________________
				! 3nd. calculate Antidiffusive edge flux: AEF=-[HO - LO]
				fct_aec(nz,edge)=cHO-cLO
				
			end do ! --> do nz=n2+1, nl1
		else
			do nz=n2+1, nl2
				!_______________________________________________________________
				Tmean2=ttfAB(nz, enodes(2))- &
						(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
						edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
						edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP
				
				Tmean1=ttfAB(nz, enodes(1))+ &
						(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
						edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
						edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP   
						
				!_______________________________________________________________
				! volume flux across the segments
				vflux=(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
				
				!_______________________________________________________________
				! 1st. Low order upwind solution
				cLO=-0.5_WP*(ttf(nz, enodes(1))*(vflux+abs(vflux))+ttf(nz, enodes(2))*(vflux-abs(vflux)))
				fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))+cLO
				fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))-cLO
				
				!_______________________________________________________________
				! 2nd. High order solution 
				! num_ord is the fraction of fourth-order contribution in the HO solution
				! (1-num_ord) is done with 3rd order upwind
				cHO=(vflux+abs(vflux))*Tmean1+(vflux-abs(vflux))*Tmean2
				cHO=-0.5_WP*((1.0_WP-num_ord)*cHO+vflux*num_ord*(Tmean1+Tmean2))
				
				!_______________________________________________________________
				! 3nd. calculate Antidiffusive edge flux: AEF=-[HO - LO]
				fct_aec(nz,edge)=cHO-cLO
			end do ! --> do nz=n2+1, nl2
		end if
	end do
	
	!___________________________________________________________________________
	! ****************** (B) do vertical tracer advection **********************
	!___________________________________________________________________________
	do n=1, myDim_nod2D
		!_______________________________________________________________________
		nl1=nlevels_nod2D(n)
		!_______________________________________________________________________
		! vert. flux at surface layer
		nz=1
		tvert(nz)=-Wvel_e(nz,n)*ttf(nz,n)*area(nz,n)
		fct_aec_ver(nz,n)=0.0_WP
		
		!_______________________________________________________________________
		! vert. flux at surface + 1 layer --> centered differences
		nz=2
		! low order
		cLO=-0.5_WP*(ttf(nz  ,n)*(Wvel_e(nz,n)+abs(Wvel_e(nz,n)))+ &
					 ttf(nz-1,n)*(Wvel_e(nz,n)-abs(Wvel_e(nz,n))))
		tvert(nz)=cLO*area(nz,n)
		
		! high order
		cHO=-0.5_WP*(ttfAB(nz-1,n)+ttfAB(nz,n))*Wvel_e(nz,n)
		
		! antidiffusive flux: (HO-LO)
		fct_aec_ver(nz,n)=(cHO-cLO)*area(nz,n)
		
		!_______________________________________________________________________
		! vert. flux at bottom - 1 layer --> centered differences
		nz=nl1-1
		! low order
		cLO=-0.5_WP*(ttf(nz  ,n)*(Wvel_e(nz,n)+abs(Wvel_e(nz,n)))+ &
					 ttf(nz-1,n)*(Wvel_e(nz,n)-abs(Wvel_e(nz,n))))
		tvert(nz)=cLO*area(nz,n)
		
		! high order
		cHO=-0.5_WP*(ttfAB(nz-1,n)+ttfAB(nz,n))*Wvel_e(nz,n)
		
		! antidiffusive flux: (HO-LO)
		fct_aec_ver(nz,n)=(cHO-cLO)*area(nz,n)
		
		!_______________________________________________________________________
		! vert. flux at bottom layer --> zero bottom flux
		nz=nl1
		tvert(nz)=0.0_WP	
		fct_aec_ver(nz,n)=0.0_WP
		
		!_______________________________________________________________________
		! Be carefull have to do vertical tracer advection here on old vertical grid
		! also horizontal advection is done on old mesh (see helem contains old 
		! mesh information)
		
		!_______________________________________________________________________
		! vert. flux at remaining levels    
		do nz=3,nl1-2
			
			! low order --> First-order upwind estimate
			cLO=-0.5*(ttf(nz  ,n)*(Wvel_e(nz,n)+abs(Wvel_e(nz,n)))+ &
					  ttf(nz-1,n)*(Wvel_e(nz,n)-abs(Wvel_e(nz,n))))
			tvert(nz)=cLO*area(nz,n)
			
			! high order --> centered (4th order)
			qc=(ttfAB(nz-1,n)-ttfAB(nz  ,n))/(Z_3d_n(nz-1,n)-Z_3d_n(nz  ,n))
			qu=(ttfAB(nz  ,n)-ttfAB(nz+1,n))/(Z_3d_n(nz  ,n)-Z_3d_n(nz+1,n))    
			qd=(ttfAB(nz-2,n)-ttfAB(nz-1,n))/(Z_3d_n(nz-2,n)-Z_3d_n(nz-1,n))
			
			Tmean1=ttfAB(nz  ,n)+(2*qc+qu)*(zbar_3d_n(nz,n)-Z_3d_n(nz  ,n))/3.0_WP
			Tmean2=ttfAB(nz-1,n)+(2*qc+qd)*(zbar_3d_n(nz,n)-Z_3d_n(nz-1,n))/3.0_WP
			Tmean =(Wvel_e(nz,n)+abs(Wvel_e(nz,n)))*Tmean1+ &
				   (Wvel_e(nz,n)-abs(Wvel_e(nz,n)))*Tmean2
			cHO=-0.5_WP*(num_ord*(Tmean1+Tmean2)*Wvel_e(nz,n)+(1.0_WP-num_ord)*Tmean)
			
			! antidiffusive flux: (HO-LO)
			fct_aec_ver(nz,n)=(cHO-cLO)*area(nz,n)
			
		end do ! --> do nz=3,nl1-1
		
		!_______________________________________________________________________
		!*************** (C) writing horizontal and vertical *******************
		!***************   low order fct advection into rhs  *******************
		!_______________________________________________________________________
		! writing horizontal and vertical low order fct advection into rhs
		do  nz=1,nlevels_nod2D(n)-1
			fct_LO(nz,n)=ttf(nz,n)+( fct_LO(nz,n)+(tvert(nz)-tvert(nz+1)) )*dt/area(nz,n)
		end do
		
	end do ! --> do n=1, myDim_nod2D
	
	!___________________________________________________________________________
	call exchange_nod(fct_LO) 
	! Summary:   
	! fct_LO contains full low-order solution
	! fct_aec contains antidiffusive component of horizontal flux 
	! fct_aec_ver contains antidiffusive component of vertical fluxes
end subroutine fct_ale_muscl_LH
!
!
!===============================================================================
subroutine fct_ale(ttf)
	!
	! 3D Flux Corrected Transport scheme
	! Limits antidiffusive fluxes==the difference in flux HO-LO
	! LO ==Low-order  (first-order upwind)
	! HO ==High-order (3rd/4th order gradient reconstruction method)
	! Adds limited fluxes to the LO solution   
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_CONFIG
	use g_comm_auto
	implicit none
	integer       :: n, nz, k, elem, enodes(3), num, el(2), nl1, nl2, edge
	real(kind=WP) :: flux, ae,tvert_max(nl-1),tvert_min(nl-1) 
	real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D)
	real*8        :: flux_eps=1e-16
	real*8        :: bignumber=1e3
	integer       :: vlimit=1
	! --------------------------------------------------------------------------
	! ttf is the tracer field on step n
	! del_ttf is the increment 
	! vlimit sets the version of limiting, see below
	! --------------------------------------------------------------------------
	
	!___________________________________________________________________________
	! a1. max, min between old solution and updated low-order solution per node
	do n=1,myDim_nod2D + edim_nod2d
		do nz=1, nlevels_nod2D(n)-1 
			fct_ttf_max(nz,n)=max(fct_LO(nz,n), ttf(nz,n))
			fct_ttf_min(nz,n)=min(fct_LO(nz,n), ttf(nz,n))
		end do
	end do       
	
	!___________________________________________________________________________
	! a2. Admissible increments on elements
	!     (only layers below the first and above the last layer)
	!     look for max, min bounds for each element --> UV_rhs here auxilary array
	do elem=1, myDim_elem2D
		enodes=elem2D_nodes(:,elem)
		do nz=1, nlevels(elem)-1
			UV_rhs(1,nz,elem)=maxval(fct_ttf_max(nz,enodes))
			UV_rhs(2,nz,elem)=minval(fct_ttf_min(nz,enodes))
		end do
		if (nlevels(elem)<=nl-1) then
			do nz=nlevels(elem),nl-1
				UV_rhs(1,nz,elem)=-bignumber
				UV_rhs(2,nz,elem)= bignumber
			end do
		endif
	end do ! --> do elem=1, myDim_elem2D
	
	!___________________________________________________________________________
	! a3. Bounds on clusters and admissible increments
	! Vertical1: In this version we look at the bounds on the clusters
	!            above and below, which leaves wide bounds because typically 
	!            vertical gradients are larger.  
	if(vlimit==1) then
		!Horizontal
		do n=1, myDim_nod2D
			!___________________________________________________________________
			do nz=1,nlevels_nod2D(n)-1
				! max,min horizontal bound in cluster around node n in every 
				! vertical layer
				! nod_in_elem2D     --> elem indices of which node n is surrounded
				! nod_in_elem2D_num --> max number of surrounded elem 
				tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
				tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
			end do
			
			!___________________________________________________________________
			! calc max,min increment of surface layer with respect to low order 
			! solution 
			fct_ttf_max(1,n)=tvert_max(1)-fct_LO(1,n)
			fct_ttf_min(1,n)=tvert_min(1)-fct_LO(1,n)
			
			! calc max,min increment from nz-1:nz+1 with respect to low order 
			! solution at layer nz
			do nz=2,nlevels_nod2D(n)-2  
				fct_ttf_max(nz,n)=maxval(tvert_max(nz-1:nz+1))-fct_LO(nz,n)
				fct_ttf_min(nz,n)=minval(tvert_min(nz-1:nz+1))-fct_LO(nz,n)
			end do
			! calc max,min increment of bottom layer -1 with respect to low order 
			! solution 
			nz=nlevels_nod2D(n)-1
			fct_ttf_max(nz,n)=tvert_max(nz)-fct_LO(nz,n)
			fct_ttf_min(nz,n)=tvert_min(nz)-fct_LO(nz,n)  
		end do
	end if
	
	!___________________________________________________________________________
	! Vertical2: Similar to the version above, but the vertical bounds are more 
	! local  
	if(vlimit==2) then
		do n=1, myDim_nod2D
			do nz=1,nlevels_nod2D(n)-1
				tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
				tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
			end do
			do nz=2, nlevels_nod2D(n)-2
				tvert_max(nz)=max(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
				tvert_min(nz)=min(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
			end do
			do nz=1,nlevels_nod2D(n)-1
				fct_ttf_max(nz,n)=tvert_max(nz)-fct_LO(nz,n)
				fct_ttf_min(nz,n)=tvert_min(nz)-fct_LO(nz,n)  
			end do
		end do
	end if
	
	!___________________________________________________________________________
	! Vertical3: Vertical bounds are taken into account only if they are narrower than the
	!            horizontal ones  
	if(vlimit==3) then
		do n=1, myDim_nod2D
			do nz=1,nlevels_nod2D(n)-1
				tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
				tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
			end do
			do nz=2, nlevels_nod2D(n)-2
				tvert_max(nz)=min(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
				tvert_min(nz)=max(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
			end do
			do nz=1,nlevels_nod2D(n)-1
				fct_ttf_max(nz,n)=tvert_max(nz)-fct_LO(nz,n)
				fct_ttf_min(nz,n)=tvert_min(nz)-fct_LO(nz,n)  
			end do
		end do
	end if
	
	!___________________________________________________________________________
	! b1. Split positive and negative antidiffusive contributions
	! --> sum all positive (fct_plus), negative (fct_minus) antidiffusive 
	!     horizontal element and vertical node contribution to node n and layer nz
	!     see. R. Löhner et al. "finite element flux corrected transport (FEM-FCT)
	!     for the euler and navier stoke equation
	do n=1, myDim_nod2D
		do nz=1,nlevels_nod2D(n)-1
			fct_plus(nz,n)=0.
			fct_minus(nz,n)=0.
		end do
	end do
	
	!Vertical
	do n=1, myDim_nod2D
		do nz=1,nlevels_nod2D(n)-1
! 			fct_plus(nz,n)=fct_plus(nz,n)+ &
! 							(max(0.0_WP,fct_aec_ver(nz,n))+max(0.0_WP,-fct_aec_ver(nz+1,n))) &
! 							/hnode(nz,n)
! 			fct_minus(nz,n)=fct_minus(nz,n)+ &
! 							(min(0.0_WP,fct_aec_ver(nz,n))+min(0.0_WP,-fct_aec_ver(nz+1,n))) &
! 							/hnode(nz,n)
			fct_plus(nz,n) =fct_plus(nz,n) +(max(0.0_WP,fct_aec_ver(nz,n))+max(0.0_WP,-fct_aec_ver(nz+1,n)))
			fct_minus(nz,n)=fct_minus(nz,n)+(min(0.0_WP,fct_aec_ver(nz,n))+min(0.0_WP,-fct_aec_ver(nz+1,n)))
		end do
	end do
	
	!Horizontal
	do edge=1, myDim_edge2D
		enodes(1:2)=edges(:,edge)   
		el=edge_tri(:,edge)
		nl1=nlevels(el(1))-1
		nl2=0
		if(el(2)>0) then
			nl2=nlevels(el(2))-1
		end if   
		do nz=1, max(nl1,nl2)
			fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0, fct_aec(nz,edge))
			fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0, fct_aec(nz,edge))  
			fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0,-fct_aec(nz,edge))
			fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0,-fct_aec(nz,edge)) 
		end do
	end do 
	
	!___________________________________________________________________________
	! b2. Limiting factors
	do n=1,myDim_nod2D
		do nz=1,nlevels_nod2D(n)-1
			flux=fct_plus(nz,n)*dt/area(nz,n)+flux_eps
			fct_plus(nz,n)=min(1.0,fct_ttf_max(nz,n)/flux)
			flux=fct_minus(nz,n)*dt/area(nz,n)-flux_eps
			fct_minus(nz,n)=min(1.0,fct_ttf_min(nz,n)/flux)
		end do
	end do 
	
	! fct_minus and fct_plus must be known to neighbouring PE
	call exchange_nod(fct_plus)
	call exchange_nod(fct_minus)	   
	
	!___________________________________________________________________________
	! b3. Limiting   
	!Vertical
	do n=1, myDim_nod2D
		nz=1
		ae=1.0_8
		flux=fct_aec_ver(nz,n)
		if(flux>=0.0) then 
			ae=min(ae,fct_plus(nz,n))
		else
			ae=min(ae,fct_minus(nz,n))
		end if
		fct_aec_ver(nz,n)=ae*fct_aec_ver(nz,n)    
		do nz=2,nlevels_nod2D(n)-1
			ae=1.0
			flux=fct_aec_ver(nz,n)
			if(flux>=0.) then 
				ae=min(ae,fct_minus(nz-1,n))
				ae=min(ae,fct_plus(nz,n))
			else
				ae=min(ae,fct_plus(nz-1,n))
				ae=min(ae,fct_minus(nz,n))
			end if   
			fct_aec_ver(nz,n)=ae*fct_aec_ver(nz,n)
		end do
	! the bottom flux is always zero 
	end do
	
	!Horizontal
	do edge=1, myDim_edge2D
		enodes(1:2)=edges(:,edge)
		el=edge_tri(:,edge)
		nl1=nlevels(el(1))-1
		nl2=0
		if(el(2)>0) then
			nl2=nlevels(el(2))-1
		end if  
		do nz=1, max(nl1,nl2)
			ae=1.0
			flux=fct_aec(nz,edge)
			if(flux>=0.) then
				ae=min(ae,fct_plus(nz,enodes(1)))
				ae=min(ae,fct_minus(nz,enodes(2)))
			else
				ae=min(ae,fct_minus(nz,enodes(1)))
				ae=min(ae,fct_plus(nz,enodes(2)))
			endif
			fct_aec(nz,edge)=ae*fct_aec(nz,edge)
		end do
	end do
	
	!___________________________________________________________________________
	! c. Update the solution
	! Vertical
	do n=1, myDim_nod2d
		do nz=1,nlevels_nod2D(n)-1  
			del_ttf(nz,n)=del_ttf(nz,n)-ttf(nz,n)+fct_LO(nz,n) + &
					(fct_aec_ver(nz,n)-fct_aec_ver(nz+1,n))*dt/area(nz,n)    
		end do
	end do
	
	! Horizontal
	do edge=1, myDim_edge2D
		enodes(1:2)=edges(:,edge)
		el=edge_tri(:,edge)
		nl1=nlevels(el(1))-1
		nl2=0
		if(el(2)>0) nl2=nlevels(el(2))-1
		do nz=1, max(nl1,nl2)
			del_ttf(nz,enodes(1))=del_ttf(nz,enodes(1))+fct_aec(nz,edge)*dt/area(nz,enodes(1))
			del_ttf(nz,enodes(2))=del_ttf(nz,enodes(2))-fct_aec(nz,edge)*dt/area(nz,enodes(2))
		end do
	end do
	
end subroutine fct_ale
