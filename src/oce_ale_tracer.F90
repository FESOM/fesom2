!
!
!===============================================================================
! Driving routine    Here with ALE changes!!!
subroutine solve_tracers_ale
	use g_parsup
	use o_PARAM, only: tracer_adv,num_tracers
	use o_arrays
	use o_mesh
	use g_comm_auto
	use o_tracers
	
	implicit none
	integer :: tr_num
	real(kind=WP) :: aux_tr(nl-1,myDim_nod2D+eDim_nod2D)
	!___________________________________________________________________________
	! loop over all tracers 
	do tr_num=1,num_tracers
	
		! do tracer AB (Adams-Bashfort) interpolation only for advectiv part 
		! needed
		call init_tracers_AB(tr_num)
		
		! advect tracers
		call adv_tracers_ale(tr_num)
		
		! diffuse tracers 
		call diff_tracers_ale(tr_num)
		
		! relax to salt and temp climatology
		call relax_to_clim(tr_num)
		
		! BRECHSTANGE: prevent temperature from running away to negative values
		! --> should violate tracer conservation
! 		if (tr_num==1) then
! 			where(tr_arr(:,:,1)<-2.1_WP)  tr_arr(:,:,1)=-2.1_WP
! 		end	if
		
		call exchange_nod(tr_arr(:,:,tr_num))
		
	end do
	
end subroutine solve_tracers_ale
!
!
!===============================================================================
subroutine adv_tracers_ale(tr_num)
	use g_parsup
	use g_config
	use o_PARAM, only: tracer_adv
	use o_arrays
	
	implicit none
	integer :: tr_num
	
	! del_ttf ... initialised and setted to zero in call init_tracers_AB(tr_num)
	! --> del_ttf ... equivalent to R_T^n in Danilov etal FESOM2: "from finite element
	!     to finite volume". At the end R_T^n should contain all advection therm and 
	!     the terms due to horizontal duffusion.
	! del_ttf=0d0
	!___________________________________________________________________________
	! horizontal ale tracer advection 
	! here --> add horizontal advection part to del_ttf(nz,n) = del_ttf(nz,n) + ...
	select case (tracer_adv)
		case(1) !MUSCL
			! --> tr_arr_old ... AB interpolated tracer from call init_tracers_AB(tr_num)
			call adv_tracers_muscle_ale(tr_arr_old(:,:,tr_num), 0.85)
			!	call adv_tracers_muscle_ale(tr_arr_old(:,:,tr_num), 0.0) ! use only third order
			!	call adv_tracers_vert_ppm_ale(tr_arr_old(:,:,tr_num))
			!	call adv_tracers_vert_cdiff(tr_arr_old(:,:,tr_num))
			call adv_tracers_vert_upw(tr_arr_old(:,:,tr_num))
			
		case(2) !MUSCL+FCT(3D)
			call adv_tracer_fct_ale(tr_arr_old(:,:,tr_num),tr_arr(:,:,tr_num), 0.85)
! 			call adv_tracer_fct_ale(tr_arr_old(:,:,tr_num),tr_arr(:,:,tr_num), 0.0)
		case default !unknown
			if (mype==0) write(*,*) 'Unknown ALE advection type. Check your namelists.'
			call par_ex(1)
	end select
	
	!___________________________________________________________________________
	! vertical ale tracer advection --> piecewise parabolic method (ppm)
	! here --> add vertical advection part to del_ttf(nz,n) = del_ttf(nz,n) + ...
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
! 			call init_tracers_AB(tr_num)
! 			tr_arr_old(:,:,tr_num)=-(0.5+epsilon)*tr_arr_old(:,:,tr_num)+(1.5+epsilon)*tr_arr(:,:,tr_num)
subroutine adv_tracers_muscle_ale(ttfAB, num_ord)
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_CONFIG
	use g_comm_auto
	implicit none
	integer       :: el(2), enodes(2), n, nz, ed
	integer       :: nl1, nl2, n2
	real(kind=WP) :: c1, deltaX1, deltaY1, deltaX2, deltaY2, vflux=0.0_WP 
	real(kind=WP) :: Tmean1, Tmean2, a
	real(kind=WP) :: ttfAB(nl-1, myDim_nod2D+eDim_nod2D)
	real(kind=WP) :: num_ord
	
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
					 edge_dxdy(1,ed)*a*edge_up_dn_grad(2,nz,ed)+ &
					 edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(4,nz,ed) &
					)/6.0_WP    
			
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
					)/6.0_WP   
					
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
			! tracer flux upwind
			! if vflux (+) --> c1 = 2*vflux*Tmean1
			! if vflux (-) --> c1 = -2*vflux*Tmean2
			! so only use upwind tracer flux !!!!!!!!!!!!!
			c1=(vflux+abs(vflux))*Tmean1+(vflux-abs(vflux))*Tmean2
			
			!___________________________________________________________________
			! combined with centered
			! num_ord is the fraction of fourth-order contribution in the HO solution
			! (1-num_ord) is done with 3rd order upwind
			c1=-0.5_WP*((1.0_WP-num_ord)*c1+vflux*num_ord*(Tmean1+Tmean2))
			
			!___________________________________________________________________
			! write horizontal ale advection into rhs
			del_ttf(nz,enodes(1))=del_ttf(nz,enodes(1))+c1*dt/area(nz,enodes(1))
			del_ttf(nz,enodes(2))=del_ttf(nz,enodes(2))-c1*dt/area(nz,enodes(2))  
			
		end do ! --> do nz=1, n2
		
		!_______________________________________________________________________
		! remaining segments on the left or on the right
		if(nl1>nl2) then 
			! be carefull !!! --> if ed is a boundary edge, el(2)==0 than nl1>0 
			!                     and nl2==0, n2=0, so for boundary edges you will 
			!                     skip the previouse do loop and always end up 
			!                     in this part of the if condition
			do nz=1+n2,nl1
				Tmean2=ttfAB(nz, enodes(2))- &
						(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
						edge_dxdy(1,ed)*a*edge_up_dn_grad(2,nz,ed)+ &
						edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(4,nz,ed))/6.0_WP    
				
				Tmean1=ttfAB(nz, enodes(1))+ &
						(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
						edge_dxdy(1,ed)*a*edge_up_dn_grad(1,nz,ed)+ &
						edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(3,nz,ed))/6.0_WP    
						
				!_______________________________________________________________
				! volume flux across the segments
				vflux=(-UV(2,nz,el(1))*deltaX1 + UV(1,nz,el(1))*deltaY1)*helem(nz,el(1)) 
				
				!_______________________________________________________________
				! tracer flux upwind
				c1=(vflux+abs(vflux))*Tmean1+(vflux-abs(vflux))*Tmean2
				
				!_______________________________________________________________
				! combined with centered
				! num_ord is the fraction of fourth-order contribution in the HO solution
				! (1-num_ord) is done with 3rd order upwind
				c1=-0.5_WP*((1.0_WP-num_ord)*c1+vflux*num_ord*(Tmean1+Tmean2))
				
				!_______________________________________________________________
				! write horizontal ale advection into rhs
				del_ttf(nz,enodes(1))=del_ttf(nz,enodes(1))+c1*dt/area(nz,enodes(1))
				del_ttf(nz,enodes(2))=del_ttf(nz,enodes(2))-c1*dt/area(nz,enodes(2)) 
			end do ! --> do nz=1+n2,nl1
		else
			do nz=n2+1,nl2
				Tmean2=ttfAB(nz, enodes(2))- &
						(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
						edge_dxdy(1,ed)*a*edge_up_dn_grad(2,nz,ed)+ &
						edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(4,nz,ed))/6.0_WP    
				
				Tmean1=ttfAB(nz, enodes(1))+ &
						(2.0_WP*(ttfAB(nz, enodes(2))-ttfAB(nz,enodes(1)))+ &
						edge_dxdy(1,ed)*a*edge_up_dn_grad(1,nz,ed)+ &
						edge_dxdy(2,ed)*r_earth*edge_up_dn_grad(3,nz,ed))/6.0_WP 
						
				!_______________________________________________________________
				! volume flux across the segments
				vflux=(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
				
				!_______________________________________________________________
				! tracer flux upwind
				c1=(vflux+abs(vflux))*Tmean1+(vflux-abs(vflux))*Tmean2
				
				!_______________________________________________________________
				! combined with centered
				! num_ord is the fraction of fourth-order contribution in the HO solution
				! (1-num_ord) is done with 3rd order upwind
				c1=-0.5_WP*((1.0_WP-num_ord)*c1+vflux*num_ord*(Tmean1+Tmean2))
				
				!_______________________________________________________________
				! write horizontal ale advection into rhs
				del_ttf(nz,enodes(1))=del_ttf(nz,enodes(1))+c1*dt/area(nz,enodes(1))
				del_ttf(nz,enodes(2))=del_ttf(nz,enodes(2))-c1*dt/area(nz,enodes(2))  
				
			end do ! --> do nz=n2+1,nl2
		end if ! --> if(nl1>nl2) then 
	end do ! --> do ed=1, myDim_edge2D
end subroutine adv_tracers_muscle_ale
!
!
!===============================================================================
! Vertical ALE advection with PPM reconstruction (5th order)
subroutine adv_tracers_vert_ppm_ale(ttf)
	use g_config
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_forcing_arrays
	implicit none
	integer       :: n, nz, nzmax
	real(kind=WP) :: tvert(nl), tv(nl)
	real(kind=WP) :: dzjm1, dzj, dzjp1, dzjp2, deltaj, deltajp1
	real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D)
	
	! --------------------------------------------------------------------------
	! Vertical advection
	! --------------------------------------------------------------------------
	! A piecewise parabolic scheme for uniformly-spaced layers.
	! See Colella and Woodward, JCP, 1984, 174-201. It can be coded so as to to take 
	! non-uniformity into account, but this is more cumbersome. This is the version for AB
	! time stepping
	! --------------------------------------------------------------------------
	
	do n=1, myDim_nod2D
		!_______________________________________________________________________
		!Interpolate to zbar...depth levels --> all quantities (tracer ...) are 
		! calculated on mid depth levels 
		
		! nzmax ... number of depth levels at node n
		nzmax=nlevels_nod2D(n)
		
		! tracer at surface layer
		tv(1)=ttf(1,n)
		
		! tracer at surface+1 layer
		tv(2)=0.5_WP*(ttf(1,n)+ttf(2,n))
		
		! tacer at bottom-2 layer
		tv(nzmax-2)=0.5_WP*(ttf(nzmax-3,n)+ttf(nzmax-2,n))
		
		! tacer at bottom-1 layer
		tv(nzmax-1)=0.5_WP*(ttf(nzmax-2,n)+ttf(nzmax-1,n))
		
		! tracer at bottom layer
		tv(nzmax)=ttf(nzmax-1,n) 
		
		!_______________________________________________________________________
		! calc tracer for surface+2 until depth-2 layer
		! see Colella and Woodward, JCP, 1984, 174-201 --> equation (1.9)
! 		do nz=2,nzmax-3
		do nz=3,nzmax-3
			!___________________________________________________________________
			! for uniform spaced vertical grids --> piecewise parabolic method (ppm)
			! equation (1.9)
			! tv(nz)=(7.0_8*(ttf(nz-1,n)+ttf(nz,n))-(ttf(nz-2,n)+ttf(nz+1,n)))/12.0_8
			
			!___________________________________________________________________
			! for non-uniformity spaced vertical grids --> piecewise parabolic 
			! method (ppm) see see Colella and Woodward, JCP, 1984, 174-201 
			! --> full equation (1.6), (1.7) and (1.8)
			dzjm1    = Z_3d_n(nz-1,n)
			dzj      = Z_3d_n(nz,n)
			dzjp1    = Z_3d_n(nz+1,n)
			dzjp2    = Z_3d_n(nz+2,n)
			
			!___________________________________________________________________
			! equation (1.7)
			! --> Here deltaj is the average slope in the jth zone of the parabola 
			!     with zone averages a_(j-1) and a_j, a_(j+1)
			! --> a_j^n
			deltaj   = dzj/(dzjm1+dzj+dzjp1)* &
					  ( &
					   (2*dzjm1+dzj    )/(dzjp1+dzj)*(ttf(nz+1,n)-ttf(nz  ,n)) +  &
					   (dzj    +2*dzjp1)/(dzjm1+dzj)*(ttf(nz  ,n)-ttf(nz-1,n)) &
					  )
			! --> a_(j+1)^n		  
			deltajp1 = dzjp1/(dzj+dzjp1+dzjp2)* &
					  ( &
					   (2*dzj+dzjp1  )/(dzjp2+dzjp1)*(ttf(nz+2,n)-ttf(nz+1,n)) +  &
					   (dzjp1+2*dzjp2)/(dzj  +dzjp1)*(ttf(nz+1,n)-ttf(nz  ,n)) &
					  )
			!___________________________________________________________________
			! condition (1.8)
			! --> This modification leads to a somewhat steeper representation of 
			!     discontinuities in the solution. It also guarantees that a_(j+0.5)
			!     lies in the range of values defined by a_j; and a_(j+1);
			if ( (ttf(nz+1,n)-ttf(nz  ,n))*(ttf(nz  ,n)-ttf(nz-1,n))>0 ) then
				deltaj = min(  abs(deltaj), &
							 2*abs(ttf(nz+1,n)-ttf(nz  ,n)),&
							 2*abs(ttf(nz  ,n)-ttf(nz-1,n)) &
							 )*sign(1.0_WP,deltaj)
			else
				deltaj = 0.0_WP
			endif
			if ( (ttf(nz+2,n)-ttf(nz+1,n))*(ttf(nz+1,n)-ttf(nz  ,n))>0 ) then
				deltajp1 = min(  abs(deltajp1),&
							   2*abs(ttf(nz+2,n)-ttf(nz+1,n)),&
							   2*abs(ttf(nz+1,n)-ttf(nz,n)) &
							   )*sign(1.0_WP,deltajp1)
			else
				deltajp1 = 0.0_WP
			endif
			
			!___________________________________________________________________
			! equation (1.6)
			! --> calcualte a_(j+0.5)
			tv(nz)=	ttf(nz,n) &
					+ dzj/(dzj+dzjp1)*(ttf(nz+1,n)-ttf(nz,n)) &
					+ 1/(dzjm1+dzj+dzjp1+dzjp2) * &
					( &
						(2*dzjp1*dzj)/(dzj+dzjp1)* &
							((dzjm1+dzj)/(2*dzj+dzjp1) - (dzjp2+dzjp1)/(2*dzjp1+dzj))*(ttf(nz+1,n)-ttf(nz,n)) &
					   - dzj*(dzjm1+dzj)/(2*dzj+dzjp1)*deltajp1 &
					   + dzjp1*(dzjp1+dzjp2)/(dzj+2*dzjp1)*deltaj &
					)
		end do ! --> do nz=3,nzmax-2
		
		!_______________________________________________________________________
		! Surface flux
		tvert(1)= -tv(1)*Wvel(1,n)*area(1,n)
		
		!_______________________________________________________________________
		! Other levels
		do nz=2, nzmax-1
			tvert(nz)= -tv(nz)*Wvel(nz,n)*area(nz,n)
		end do
		
		!_______________________________________________________________________
		! Zero bottom flux
		tvert(nzmax)=0.0_WP
		
		!_______________________________________________________________________
		! writing vertical ale advection into rhs
		do nz=1, nzmax-1
			! no division over thickness in ALE !!!
			del_ttf(nz,n)=del_ttf(nz,n) + (tvert(nz)-tvert(nz+1))*dt/area(nz,n) 
		end do         
		
	end do ! --> do n=1, myDim_nod2D
	
end subroutine adv_tracers_vert_ppm_ale
!
!
!===============================================================================
! Vertical ALE advection with upwind reconstruction (1st order)
subroutine adv_tracers_vert_upw(ttf)
	use g_config
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_forcing_arrays
	implicit none
	integer       :: n, nz, nl1
	real(kind=WP) :: tvert(nl), tv
	real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D)
	
	! --------------------------------------------------------------------------
	! Vertical advection
	! --------------------------------------------------------------------------
	do n=1, myDim_nod2D
		!_______________________________________________________________________
		nl1=nlevels_nod2D(n)-1
		!_______________________________________________________________________
		! Surface flux
		tvert(1)= -Wvel(1,n)*ttf(1,n)*area(1,n)		
		!_______________________________________________________________________
		! Zero bottom flux
		tvert(nl1+1)=0.0_WP
		!_______________________________________________________________________
		! Other levels
		do nz=2, nl1
			tv=ttf(nz-1,n)*min(Wvel(nz,n), 0._WP)+ttf(nz,n)*max(Wvel(nz,n), 0._WP)
			tvert(nz)= -tv*area(nz,n)
		end do
		!_______________________________________________________________________
		! writing vertical ale advection into rhs
		do nz=1, nl1
			! no division over thickness in ALE !!!
			del_ttf(nz,n)=del_ttf(nz,n) + (tvert(nz)-tvert(nz+1))*dt/area(nz,n) 
			
		end do         
	end do ! --> do n=1, myDim_nod2D
end subroutine adv_tracers_vert_upw
!
!
!===============================================================================
! Vertical ALE advection with central difference reconstruction (2nd order)
subroutine adv_tracers_vert_cdiff(ttf)
	use g_config
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_forcing_arrays
	implicit none
	integer       :: n, nz, nl1
	real(kind=WP) :: tvert(nl), tv
	real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D)
	
	! --------------------------------------------------------------------------
	! Vertical advection
	! --------------------------------------------------------------------------
	do n=1, myDim_nod2D
		!_______________________________________________________________________
		nl1=nlevels_nod2D(n)-1
		!_______________________________________________________________________
		! Surface flux
		tvert(1)= -Wvel(1,n)*ttf(1,n)*area(1,n)		
		!_______________________________________________________________________
		! Zero bottom flux
		tvert(nl1+1)=0.0_WP		
		!_______________________________________________________________________
		! Other levels
		do nz=2, nl1
			tv=0.5_WP*(ttf(nz-1,n)+ttf(nz,n))
			tvert(nz)= -tv*Wvel(nz,n)*area(nz,n)
		end do
		!_______________________________________________________________________
		! writing vertical ale advection into rhs
		do nz=1, nl1
			! no division over thickness in ALE !!!
			del_ttf(nz,n)=del_ttf(nz,n) + (tvert(nz)-tvert(nz+1))*dt/area(nz,n) 
			
		end do         
	end do ! --> do n=1, myDim_nod2D
end subroutine adv_tracers_vert_cdiff
!
!
!===============================================================================
subroutine diff_tracers_ale(tr_num)
	use o_mesh
	use g_PARSUP
	use o_arrays
	use o_tracers
	implicit none
	
	integer, intent(in) :: tr_num
	integer             :: n, nzmax
	
	!___________________________________________________________________________
	! convert tr_arr_old(:,:,tr_num)=ttr_n-0.5   --> prepare to calc ttr_n+0.5
	! eliminate AB (adams bashfort) interpolates tracer, which is only needed for 
	! tracer advection. For diffusion only need tracer from previouse time step
	tr_arr_old(:,:,tr_num)=tr_arr(:,:,tr_num) !DS: check that this is the right place!
	
	!___________________________________________________________________________
	! do horizontal diffusiion
	! write there also horizontal diffusion rhs to del_ttf which is equal the R_T^n 
	! in danilovs srcipt
! 		call diff_part_hor
	call diff_part_hor_sergey ! seems to be ~9% faster than diff_part_hor
	
	!___________________________________________________________________________
	! do vertical diffusion: explicite 
	if (.not. i_vert_diff) call diff_ver_part_expl_ale(tr_num)
		
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
		call diff_ver_part_impl_ale(tr_num)
		
	end if
	
	!We DO not set del_ttf to zero because it will not be used in this timestep anymore
	!init_tracers will set it to zero for the next timestep
	
end subroutine diff_tracers_ale
!
!
!===============================================================================
!Vertical diffusive flux(explicit scheme):                                                                            
subroutine diff_ver_part_expl_ale(tr_num)
	use o_ARRAYS
	use o_MESH
	use g_PARSUP
	use g_config,only: dt
	
	implicit none 
	
	real(kind=WP) :: vd_flux(nl-1)
	real(kind=WP) :: rdata,flux,rlx
	integer       :: nz,nl1,tr_num,n
	real(kind=WP) :: zinv1,Ty
	
	do n=1, myDim_nod2D
		nl1=nlevels_nod2D(n)-1
		vd_flux=0d0
		if (tr_num==1) then
			flux  = -heat_flux(n)/vcpw
			rdata =  Tsurf(n)
			rlx   =  surf_relax_T
		elseif (tr_num==2) then
			flux  =  water_flux(n)*tr_arr(1,n,2)
			rdata =  Ssurf(n)
			rlx   =  surf_relax_S
		else
			flux  = 0d0
			rdata = 0d0
			rlx=0d0
		endif
		
		!_______________________________________________________________________
		!Surface forcing
		vd_flux(1)= flux + rlx*(rdata-tr_arr(1,n,tr_num))
		
		!_______________________________________________________________________
		do nz=2,nl1
			!___________________________________________________________________
			zinv1=1.0_WP/(Z_3d_n(nz-1,n)-Z_3d_n(nz,n))
			
			!___________________________________________________________________
			Ty= Kd(4,nz-1,n)*(Z_3d_n(nz-1,n)-zbar_3d_n(nz,n))*zinv1 *neutral_slope(3,nz-1,n)**2 + &
				Kd(4,nz,n)*(zbar_3d_n(nz,n)-Z_3d_n(nz,n))*zinv1 *neutral_slope(3,nz,n)**2
			
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
!
!===============================================================================
subroutine diff_ver_part_impl_ale(tr_num)
	use o_MESH
	use o_PARAM
	use o_ARRAYS
        use i_ARRAYS
	use g_PARSUP
	use g_CONFIG
	use g_forcing_arrays
		
	implicit none
	
	real(kind=WP)       :: a(nl), b(nl), c(nl), tr(nl)
	real(kind=WP)       :: cp(nl), tp(nl)
	integer             :: nz, n, nzmax,tr_num
	real(kind=WP)       :: m, zinv, dt_inv
	real(kind=WP)       :: rsss, Ty,Ty1,c1,zinv1,zinv2,v_adv
        real(kind=WP), external    :: TFrez  ! Sea water freeze temperature.

	dt_inv=1.0_WP/dt
	Ty=0.0_WP
	Ty1=0.0_WP
	
	! solve equation diffusion equation implicite part: 
	! --> h^(n+0.5)* (T^(n+0.5)-Tstar) = dt*( K_33*d/dz*(T^(n+0.5)-Tstar) + K_33*d/dz*Tstar )
	! -->   dTnew = T^(n+0.5)-Tstar
	! -->   h^(n+0.5)* (dTnew) = dt*(K_33*d/dz*dTnew) + K_33*dt*d/dz*Tstar 
	! -->   h^(n+0.5)* (dTnew) = dt*(K_33*d/dz*dTnew) + RHS 
	! -->   solve for dT_new
	!	
	!	----------- zbar_1, V_1 (Volume eq. to Area)
	! Z_1 o T_1
	!	----------- zbar_2, V_2
	! Z_2 o T_2
	!	----------- zbar_3, V_3
	! Z_3 o T_3
	!	----------- zbar_4
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
		zbar_n(nzmax)=zbar(nzmax)
		Z_n(nzmax-1)=zbar_n(nzmax) + hnode_new(nzmax-1,n)/2.0_WP
		do nz=nzmax-1,2,-1
			zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
			Z_n(nz-1) = zbar_n(nz) + hnode_new(nz-1,n)/2.0_WP
		end do
		zbar_n(1) = zbar_n(2) + hnode_new(1,n)
		
		!_______________________________________________________________________
		! Regular part of coefficients: --> surface layer 
		nz=1
		
		! 1/dz(nz)
		zinv2=1.0_WP/(Z_n(nz)-Z_n(nz+1))
		zinv=1.0_WP*dt    ! no .../(zbar(1)-zbar(2)) because of  ALE
		
		! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
! 		Ty1= Kd(4,nz,n)  *(Z_n(nz)     -zbar_n(nz+1))*zinv2 *neutral_slope(3,nz,n  )**2 + &
! 			 Kd(4,nz+1,n)*(zbar_n(nz+1)-Z_n(nz+1)   )*zinv2 *neutral_slope(3,nz+1,n)**2
		
		! layer dependent coefficients for for solving dT(1)/dt+d/dz*K_33*d/dz*T(1) = ...
		a(nz)=0.0_WP
		c(nz)=-(Kv(2,n)+Ty1)*zinv2*zinv*area(nz+1,n)/area(nz,n)
		b(nz)=-c(nz)+hnode_new(nz,n)      ! ale
		
		! update from the vertical advection --> comes from splitting of vert 
		! velocity into explicite and implicite contribution
! 		v_adv=zinv*area(2,n)/area(1,n)
! 		b(1)=b(1)+Wvel_i(1, n)*zinv-min(0._WP, Wvel_i(2, n))*v_adv
! 		c(1)=c(1)-max(0._WP, Wvel_i(2, n))*v_adv
		
		! backup zinv2 for next depth level
		zinv1=zinv2
		
		!_______________________________________________________________________
		! Regular part of coefficients: --> 2nd...nl-2 layer
		do nz=2, nzmax-2
		
			! 1/dz(nz)
			zinv2=1.0_WP/(Z_n(nz)-Z_n(nz+1))
			
			! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
! 			Ty = Kd(4,nz-1,n)*(Z_n(nz-1)-zbar_n(nz))*zinv1 *neutral_slope(3,nz-1,n)**2 + &
! 				 Kd(4,nz,n)*(zbar_n(nz)-Z_n(nz))*zinv1 *neutral_slope(3,nz,n)**2
! 			Ty1= Kd(4,nz,n)*(Z_n(nz)-zbar_n(nz+1))*zinv2 *neutral_slope(3,nz,n)**2 + &
! 				 Kd(4,nz+1,n)*(zbar_n(nz+1)-Z_n(nz+1))*zinv2 *neutral_slope(3,nz+1,n)**2
			
			! layer dependent coefficients for for solving dT(nz)/dt+d/dz*K_33*d/dz*T(nz) = ...
			a(nz)=-(Kv(nz,n)+Ty)*zinv1*zinv
			c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv*area(nz+1,n)/area(nz,n)
			b(nz)=-a(nz)-c(nz)+hnode_new(nz,n)
			
			! backup zinv2 for next depth level
			zinv1=zinv2
			
			! update from the vertical advection
! 			v_adv=zinv
! 			a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv
! 			b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv
! 			
! 			v_adv=v_adv*area(nz+1,n)/area(nz,n)
! 			b(nz)=b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
! 			c(nz)=c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
		end do ! --> do nz=2, nzmax-2
		
		!_______________________________________________________________________
		! Regular part of coefficients: --> nl-1 layer
		nz=nzmax-1 
		
		zinv=1.0_WP*dt   ! no ... /(zbar(nzmax-1)-zbar(nzmax)) because of ale
		
		! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
! 		Ty= Kd(4,nz-1,n)*(Z_n(nz-1)-zbar_n(nz))*zinv1 *neutral_slope(3,nz-1,n)**2 + &
! 			Kd(4,nz,n)*(zbar_n(nz)-Z_n(nz))*zinv1 *neutral_slope(3,nz,n)**2
		
		! layer dependent coefficients for for solving dT(nz)/dt+d/dz*K_33*d/dz*T(nz) = ...
		a(nz)=-(Kv(nz,n)+Ty)*zinv1*zinv
		c(nz)=0.0_WP
		b(nz)=-a(nz)+hnode_new(nz,n)
		
		! update from the vertical advection
! 		v_adv=zinv
! 		a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv       
! 		b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv
		
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
		tr(nz)=c(nz)*(tr_arr(nz,n,tr_num) - tr_arr(nz+1,n,tr_num))
		do nz=2,nzmax-2
			tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num) &
			       -c(nz)*tr_arr(nz+1,n,tr_num) &
			       +(a(nz)+c(nz))*tr_arr(nz,n,tr_num)
		end do
		nz=nzmax-1
		tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num)+a(nz)*tr_arr(nz,n,tr_num)
		
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
		!      cell volume. If Volume compressed --> temp has tto raise, if volume 
		!      expended --> temp has to decrease
		zinv=1.0_WP*dt    !/(zbar(1)-zbar(2))  ! ale
		if (tr_num==1) then
			tr(1)= tr(1) - zinv*(&
								heat_flux(n)/vcpw & 
  	 						        - tr_arr_old(1,n,1)*water_flux(n) & !
                                                                ! the proper way would be the following (not tested yet)
								- ((1.-a_ice(n))*tr_arr_old(1,n,1)+a_ice(n)*TFrez(tr_arr_old(1,n,2)))*water_flux(n) & !
								- tr_arr(1,n,1)*water_flux(n) & !
								- surf_relax_T*(Tsurf(n)-tr_arr(1,n,1))&
								)
			
		elseif (tr_num==2) then
			!___________________________________________________________________
			! set reference surface salinity if local or global
			rsss=ref_sss
			if(ref_sss_local) rsss = tr_arr(1,n,2)
			
			!___________________________________________________________________
			! on freshwater inflow/outflow or virtual salinity:
			! 1. In zlevel & zstar the freshwater flux is applied in the update of the 
			! ssh matrix when solving the continuity equation of vertically 
			! integrated flow. The salt concentration in the first layer will 
			! be then adjusted according to the change in volume.
			! In this case rsss is forced to be zero by setting ref_sss=0. and ref_sss_local=.false.
			! in routines above.
			! 2. In cases where the volume of the upper layer is fixed (i.e. linfs)  the freshwater flux 
			! 'rsss*water_flux(n)' is applied as a virtual salt boundary condition via the vertical 
			! diffusion operator.
			! --> real_salt_flux(:): salt flux due to containment/releasing of salt
			!     by forming/melting of sea ice
			! --> rsss*water_flux(n) : virtual salt flux 
			tr(1)= tr(1)  +  zinv*( &
									rsss*water_flux(n) &
									+ real_salt_flux(n) &
									+ surf_relax_S*(Ssurf(n)-tr_arr(1,n,2)))
		endif
		
		!_______________________________________________________________________
		! The forward sweep algorithm to solve the three-diagonal matrix 
		! problem
		! 
		!  | b_1 c_1 ...		    |   |dTnew_1|
		!  | a_2 b_2 c_2 ...	    |   |dTnew_2|
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
subroutine diff_part_hor_sergey
	use o_ARRAYS
	use g_PARSUP
	use o_MESH
	USE o_param
	use g_config
	IMPLICIT NONE
	real(kind=WP)   :: deltaX1,deltaY1,deltaX2,deltaY2
	integer         :: edge
	integer         :: n2,nl1,nl2,nz,el(2),elnodes(3),n,enodes(2)
	real(kind=WP)   :: c, Tx, Ty, Kh
	real(kind=WP)   :: rhs1(nl-1), rhs2(nl-1)
	
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
		Kh=elem_area(el(1))
		
		!_______________________________________________________________________
		nl2=0
		if (el(2)>0) then 
			Kh=0.5_WP*(Kh+elem_area(el(2)))
			nl2=nlevels(el(2))-1
			deltaX2=edge_cross_dxdy(3,edge)
			deltaY2=edge_cross_dxdy(4,edge)
		endif
		Kh=K_hor*Kh/scale_area
		
		!_______________________________________________________________________
		n2=min(nl1,nl2)
		do nz=1,n2
			Tx=0.5_WP*Kh*(tr_xy(1,nz,el(1))+tr_xy(1,nz,el(2)))
			Ty=0.5_WP*Kh*(tr_xy(2,nz,el(1))+tr_xy(2,nz,el(2)))
			c=(deltaX2-deltaX1)*Ty-(deltaY2-deltaY1)*Tx
			rhs1(nz) = rhs1(nz) + c
			rhs2(nz) = rhs2(nz) - c
		enddo
		
		!_______________________________________________________________________
		if(nl1>nl2) then
			do nz=n2+1,nl1
				Tx=Kh*tr_xy(1,nz,el(1))
				Ty=Kh*tr_xy(2,nz,el(1))
				c=-deltaX1*Ty+deltaY1*Tx
				rhs1(nz) = rhs1(nz) + c
				rhs2(nz) = rhs2(nz) - c
			end do  
		else
			do nz=n2+1,nl2
				Tx=Kh*tr_xy(1,nz,el(2))
				Ty=Kh*tr_xy(2,nz,el(2))
				c=deltaX2*Ty-deltaY2*Tx
				rhs1(nz) = rhs1(nz) + c
				rhs2(nz) = rhs2(nz) - c
			end do
		endif
		
		!_______________________________________________________________________
		n2=max(nl1,nl2)
		del_ttf(1:n2,enodes(1))=del_ttf(1:n2,enodes(1))+rhs1(1:n2)*dt/area(1:n2,enodes(1))
		del_ttf(1:n2,enodes(2))=del_ttf(1:n2,enodes(2))+rhs2(1:n2)*dt/area(1:n2,enodes(2))
		
	end do
end subroutine diff_part_hor_sergey
