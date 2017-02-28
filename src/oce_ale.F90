! CONTENT:
! ------------
!	subroutine ale_init
!	subroutine init_bottom_elem_thickness
!	subroutine init_thickness_ale
!	subroutine update_thickness_ale
!	subroutine stiff_mat_update
!	subroutine compute_ssh_rhs_ale
!	subroutine compute_hbar_ale
!	subroutine vert_vel_ale
!	subroutine solve_ssh_ale
!	subroutine oce_timestep_ale
!	
!	
!===============================================================================
! allocate & initialise arrays for Arbitrary-Langrangian-Eularian (ALE) method
subroutine ale_init
	USE o_PARAM
	USE o_MESH
	USE g_PARSUP
	USE o_ARRAYS
	USE g_config, only: which_ale
        USE g_forcing_param, only: use_virt_salt
	Implicit NONE
	
	
	!___allocate________________________________________________________________
	! hnode and hnode_new: layer thicknesses at nodes. 
	allocate(hnode(1:nl-1, myDim_nod2D+eDim_nod2D))
	allocate(hnode_new(1:nl-1, myDim_nod2D+eDim_nod2D))
	
	! ssh_rhs_old: auxiliary array to store an intermediate part of the rhs computations.
	allocate(ssh_rhs_old(myDim_nod2D+eDim_nod2D))
	
	! hbar, hbar_old: correspond to the elevation, but on semi-integer time steps.
	allocate(hbar(myDim_nod2D+eDim_nod2D))
	allocate(hbar_old(myDim_nod2D+eDim_nod2D))
	
	! helem: layer thickness at elements. It is interpolated from hnode.
	allocate(helem(1:nl-1, myDim_elem2D))
	
	! dhe: The increment of total fluid depth on elements. It is used to update the matrix
	! of the ssh operator.      
	allocate(dhe(myDim_elem2D))
	
	! zbar_n: depth of layers due to ale thinkness variactions at ervery node n 
	allocate(zbar_n(nl))
	
	! Z_n: mid depth of layers due to ale thinkness variactions at ervery node n 
	allocate(Z_n(nl-1))
	
	! bottom_elem_tickness: changed bottom layer thinkness due to partial cells
	allocate(bottom_elem_thickness(myDim_elem2D))
	
	!___initialize______________________________________________________________
	hbar=0.0_WP
	hbar_old=hbar
	dhe=0.0_WP
	
	! calculate thinkness of partial bottom layer cells
	call init_bottom_elem_thickness
	! if the ale thinkness remain unchanged (like in 'linfs' case) the vitrual salinity flux need to be used
	! otherwise we set the reference salinity to zero
	if ( .not. trim(which_ALE)=='linfs') then
		use_virt_salt=.false.
		! this will force the virtual saltinity flux to be zero
		ref_sss_local=.false.
		ref_sss=0._WP
	else
		use_virt_salt=.true.
	end if
end subroutine ale_init
!
!
!===============================================================================
subroutine init_bottom_elem_thickness
	use o_PARAM
	use o_MESH
	use g_PARSUP
	use o_ARRAYS
	use g_config,only: use_partial_cell
	implicit none
	
	integer :: elem, elnodes(3), nle
	real(kind=WP) :: dd, dd1
	
	! If we use partial cells, the thickness of bottom cell is adjusted.
	! The adjustment is limited. It cannot be more than + (1/2) of the deeper
	! layer, nor -(1/2) of the current layer. 
	if(use_partial_cell) then 
		!Adjust the thickness of bottom cells
		DO elem=1, myDim_elem2D
			elnodes=elem2D_nodes(:,elem) 
			
			! elemental topographic depth
			dd=sum(depth(elnodes))/3.0_WP
			
			! number of full depth levels at elem
			nle=nlevels(elem)    
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
					dd1=zbar(nle-1)-dd   
				else
					! case 1 : max(Z(nle),dd) = dd
					! case 2 : max(Z(nle),dd) = Z(nle)
					dd1=zbar(nle-1)-max(Z(nle),dd)
				end if
			!___________________________________________________________________
			! if topographic depth dd is shallower than depth of deepest full cell 
			! depth level zbar(nle)
			!        : 
			!        : 
			! ______________ zbar(nle-1)--------->+---->+
			!                                     |--dd case1--> dd1=
			!      o dd case1                     |     |
			! -------------- Z(nle-1)(mid-depth)->+     |--dd case 2 --> dd1=
			!      o dd case2 ------------------------->+
			! ______________ zbar(nle) 
			! / / / / / / / 
			!  / / / / / / /
			! / / / / / / / 
			else
				! case 1 : min(Z(nle-1),dd) = Z(nle-1)
				! case 2 : min(Z(nle-1),dd) = dd
				dd1=zbar(nle-1)-min(Z(nle-1),dd) 
			end if	    
			bottom_elem_thickness(elem)=dd1
		END DO
	else
		DO elem=1, myDim_elem2D
			nle=nlevels(elem)
			bottom_elem_thickness(elem)=zbar(nle-1)-zbar(nle)
		end do
	end if 
END subroutine init_bottom_elem_thickness
!
!
!===============================================================================
! initialize thickness arrays based on the current hbar 
subroutine init_thickness_ale
! For z-star case: we stretch scalar thicknesses (nodal) 
! through nlevels_nod2D_min -2 layers. Layer nlevels_nod2D_min-1
! should not be touched if partial cell is implemented (it is).
! In lower layers scalar prisms are modified by the bottom.  
! Important: nlevels_nod2D_min has to be allocated and filled. 
	use g_config,only: dt, which_ale
	use o_PARAM
	use o_MESH
	use g_PARSUP
	use o_ARRAYS
	implicit none
	integer :: n, nz, elem, elnodes(3)
	real(kind=WP) :: dd 
	!___________________________________________________________________________
	! Fill in ssh_rhs_old
	ssh_rhs_old=(hbar-hbar_old)*area(1,:)/dt
	
	! -->see equation (8) FESOM2:from finite elements to finie volume
	eta_n=alpha*hbar_old+(1.0_WP-alpha)*hbar   
	
	if     (trim(which_ale)=='linfs') then
		!_______________________________________________________________________
		! >->->->->->->->->->->->     linear-free-surface    <-<-<-<-<-<-<-<-<-<
		!_______________________________________________________________________
		! no layer thickness variation in any layer
		do n=1,myDim_nod2D+eDim_nod2D
			hnode(1,n)=(zbar(1)-zbar(2))
			do nz=2,nlevels_nod2D(n)-1
				hnode(nz,n)=(zbar(nz)-zbar(nz+1))
			end do      
			do nz=nlevels_nod2D(n),nl-1
				hnode(nz,n)=0.0_WP
			end do
		end do
		do elem=1,myDim_elem2D
			dhe(elem)=0.0_8
			do nz=1,nlevels(elem)-2
				helem(nz,elem)=(zbar(nz)-zbar(nz+1))
			end do
			helem(nlevels(elem)-1,elem)=bottom_elem_thickness(elem)
			Do nz=nlevels(elem),nl-1
				helem(nz,elem)=0.0_WP
			end do
		END DO
	
	elseif (trim(which_ale)=='zlevel') then
		!_______________________________________________________________________
		! >->->->->->->->->->->->->->->     z-level    <-<-<-<-<-<-<-<-<-<-<-<-<
		!_______________________________________________________________________
		! --> include all ssh variations into the top layer 
		do n=1,myDim_nod2D+eDim_nod2D
			
			! put all ssh variation (hbar) into first layer 
			hnode(1,n)=hbar(n)+(zbar(1)-zbar(2))
			
			! leave lower levels untouched
			do nz=2,nlevels_nod2D(n)-1
			hnode(nz,n)=(zbar(nz)-zbar(nz+1))
			end do 
			
			! layer thickness of bottom layer equal 0
			do nz=nlevels_nod2D(n),nl-1
				hnode(nz,n)=0.0_WP
			end do
		end do
		
		do elem=1,myDim_elem2D
			
			elnodes=elem2D_nodes(:,elem) 
			! interpolated ssh variation at element elem
			dhe(elem)=sum(hbar(elnodes))/3.0_WP
			
			! store elemtal ssh varition only in first layer
			helem(1,elem)=dhe(elem)+(zbar(1)-zbar(2))
			
			! lower layers leave untouched 
			do nz=2,nlevels(elem)-2
				helem(nz,elem)=(zbar(nz)-zbar(nz+1))
			end do
			
			! elemental bottom layer thickness
			helem(nlevels(elem)-1,elem)=bottom_elem_thickness(elem)
			
			! fill thickness below bottom layer
			do nz=nlevels(elem),nl-1
				helem(nz,elem)=0.0_WP
			end do
			
		end do
	
	elseif (trim(which_ale)=='zstar' ) then
		!_______________________________________________________________________
		! >->->->->->->->->->->->->->->     z-star     <-<-<-<-<-<-<-<-<-<-<-<-<
		!_______________________________________________________________________
		! --> calcualte layer thinkness at depth layer and node
		do n=1, myDim_nod2D+eDim_nod2D
			! depth anomaly until the last minus one level where the scalar prism is not 
			! intersected with bottom.
			dd=zbar(1)-zbar(nlevels_nod2D_min(n)-1)  
			
			! calc layer thinkness for depth layer nz and node n. distribute hbar surface 
			! elevation linear over verical column
			do nz=1,nlevels_nod2D_min(n)-2
				hnode(nz,n)=(zbar(nz)-zbar(nz+1))*(1.0_WP+hbar(n)/dd)
			end do
			
			! do not distrubute hbar into cells that intersect somehow with bottom layer 
			do nz=nlevels_nod2D_min(n)-1, nlevels_nod2D(n)-1
				hnode(nz,n)=(zbar(nz)-zbar(nz+1))
			end do
			
			! layer thickness of bottom layer equal 0
			do nz=nlevels_nod2D(n),nl-1
				hnode(nz,n)=0.0_WP
			end do
		end do
		
		!_______________________________________________________________________
		! --> calculate mean layer thinkness at element
		do elem=1, myDim_elem2D
			elnodes=elem2D_nodes(:, elem)
			! interpolated ssh variation at element elem
			dhe(elem)=sum(hbar(elnodes))/3.0_WP
			do nz=1,nlevels(elem)-2
				helem(nz,elem)=sum(hnode(nz,elnodes))/3.0_WP
			end do
			
			! elemental bottom layer thickness
			helem(nlevels(elem)-1,elem)=bottom_elem_thickness(elem)
			
			! fill thickness below bottom layer
			do nz=nlevels(elem),nl-1
				helem(nz,elem)=0.0_WP
			end do
		end do
	else
		if (mype==0) then
		write(*,*)
		write(*,*) '____________________________________________________________'
		write(*,*) 'The vertical ALE discretisation ', which_ale,' is currently not supported!!!'
		call par_ex(1)
		end if 
	endif
	
	!___________________________________________________________________________
	hnode_new=hnode  ! Should be initialized, because only variable part is updated.
end subroutine init_thickness_ale
!
!
!===============================================================================
! update thickness arrays based on the current hbar 
subroutine update_thickness_ale
	use o_PARAM
	use o_MESH
	use g_PARSUP
	use o_ARRAYS
	use g_config,only: which_ale
	implicit none
	integer :: n, nz, elem, elnodes(3)
	
	if     (trim(which_ale)=='zlevel') then
		!_______________________________________________________________________
		! >->->->->->->->->->->->->->->     z-level    <-<-<-<-<-<-<-<-<-<-<-<-<
		!_______________________________________________________________________
		! only actualize layer thinkness in first layer 
		do n=1,myDim_nod2D+eDim_nod2D
			hnode(1,n)=hnode_new(1,n)
		end do
		
		do elem=1,myDim_elem2D
			helem(1,elem)=sum(hnode(1,elem2d_nodes(:,elem)))/3.0_WP
		end do
		
	elseif (trim(which_ale)=='zstar' ) then
		!_______________________________________________________________________
		! >->->->->->->->->->->->->->->     z-star     <-<-<-<-<-<-<-<-<-<-<-<-<
		!_______________________________________________________________________
		! --> update layer thinkness at depth layer and node
		do n=1, myDim_nod2D+eDim_nod2D
			do nz=1,nlevels_nod2D_min(n)-2
				hnode(nz,n) = hnode_new(nz,n)
			end do
		end do
		
		!_______________________________________________________________________
		! --> update mean layer thinkness at element
		do elem=1, myDim_elem2D
			elnodes=elem2D_nodes(:, elem)
			do nz=1,nlevels(elem)-2
				helem(nz,elem)=sum(hnode(nz,elnodes))/3.0_WP
			end do
		end do
	endif
end subroutine update_thickness_ale
!
!
!_______________________________________________________________________________
! Update ssh stiffness matrix for a new elevation
subroutine stiff_mat_update
	use g_config,only: dt
	use o_PARAM
	use o_MESH
	use g_PARSUP
	use o_ARRAYS
	!
	implicit none
	integer                             :: n, i, j,  row, ed,n2
	integer                             :: enodes(2), elnodes(3), el(2)
	integer                             :: elem, npos(3), offset, nini, nend
	real(kind=WP)                       :: factor 
	real(kind=WP)                       :: fx(3), fy(3)
	integer, allocatable                :: n_num(:)
	
	allocate(n_num(myDim_nod2D+eDim_nod2D))
	n_num=0
	factor=g*dt*alpha*theta
	!___________________________________________________________________________
	! loop over lcal edges
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
			! loop over number of neghbouring nodes of node-row
			do n=1,SSH_stiff%rowptr(row+1)-SSH_stiff%rowptr(row)
				! nn_pos ... local indice position of neigbouring nodes
				! n2 ... local indice of n-th neighbouring node to node-row
				n2=nn_pos(n,row)
				! n_num(n2) ... global sparse matrix indices of local mesh point n2
				n_num(n2)=offset+n
			end do
			
			!___________________________________________________________________
			dO i=1,2  ! Two elements related to the edge
						! It should be just grad on elements 
				! elem ... local element index to calc grad on that element
				elem=el(i)
				if(elem<1) cycle
				
				! elnodes ... local node indices of nodes that form element elem
				elnodes=elem2D_nodes(:,elem)
				
				! calc gradient on elem
				fx=dhe(elem)*gradient_sca(1:3,elem)
				fy=dhe(elem)*gradient_sca(4:6,elem)
				
				! fx, fy are contribution to -velocity from elem   (-h\nabla\eta)
				fy=fy*(edge_cross_dxdy(2*i-1,ed))- fx*(edge_cross_dxdy(2*i,ed))
				
				if(i==2) fy=-fy
				if(j==2) fy=-fy
				 
				! In the computation above, I've used rules from ssh_rhs (where it is 
				! on the rhs. So the sign is changed in the expression below.
				! npos... sparse matrix indices position of node points elnodes
				npos=n_num(elnodes)
				SSH_stiff%values(npos)=SSH_stiff%values(npos)- fy*factor
			end do
		end do
	end do 
	deallocate(n_num)
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

end subroutine stiff_mat_update
!
!
!===============================================================================
subroutine compute_ssh_rhs_ale
	use g_config,only: which_ALE,dt
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	implicit none
	
	! In the semiimplicit method: 
	! ssh_rhs=-alpha*\nabla\int(U_n+U_rhs)dz-(1-alpha)*...
	! see "FESOM2: from finite elements to finte volumes, S. Danilov..." eq. (11) rhs
	integer       :: ed, el(2), enodes(2),  nz,n
	real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2 
	
	ssh_rhs=0.0_WP
	!___________________________________________________________________________
	! loop over local edges
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
		do nz=1, nlevels(el(1))-1
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
			do nz=1, nlevels(el(2))-1
				c2=c2-alpha*((UV(2,nz,el(2))+UV_rhs(2,nz,el(2)))*deltaX2- &
					(UV(1,nz,el(2))+UV_rhs(1,nz,el(2)))*deltaY2)*helem(nz,el(2))
			end do
		end if
		
		!_______________________________________________________________________
		! calc netto "flux"
		ssh_rhs(enodes(1))=ssh_rhs(enodes(1))+(c1+c2)
		ssh_rhs(enodes(2))=ssh_rhs(enodes(2))-(c1+c2)
	end do
		
	!___________________________________________________________________________
	! take into account water flux
	! at this point ssh_rhs_old= -alpha* nabla* int(u^n + deltau dz) 
	! shown in eq (11) rhs of "FESOM2: from finite elements to finte volumes, S. Danilov..." eq. (11) rhs
	if ( .not. trim(which_ALE)=='linfs') then
		do n=1,myDim_nod2D
			ssh_rhs(n)=ssh_rhs(n)-alpha*water_flux(n)+(1.0_WP-alpha)*ssh_rhs_old(n)
			
! 			ssh_rhs(n)=ssh_rhs(n)-alpha*water_flux(n)*area(1,n)+(1.0_WP-alpha)*ssh_rhs_old(n) !!PS

		end do
	else
		do n=1,myDim_nod2D
			ssh_rhs(n)=ssh_rhs(n)+(1.0_WP-alpha)*ssh_rhs_old(n)
		end do
	end if
end subroutine compute_ssh_rhs_ale
!
!
!===============================================================================
subroutine compute_hbar_ale
	use g_config,only: dt, which_ALE
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_comm_auto
	
	implicit none

	! see "FESOM2: from finite elements to finte volumes, S. Danilov..." 
	! hbar(n+1)-hbar(n)=tau*ssh_rhs_old 
	! ssh_rhs_old=-\nabla\int(U_n)dz-water_flux*area (if free surface)
	! Find new elevation hbar
	
	integer      :: ed, el(2), enodes(2),  nz,n, elnodes(3), elem
	real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2 
	!___________________________________________________________________________
	! compute the rhs
	ssh_rhs_old=0.0_WP
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
		do nz=1, nlevels(el(1))-1
			c1=c1+(UV(2,nz,el(1))*deltaX1-UV(1,nz,el(1))*deltaY1)*helem(nz,el(1))
		end do
		!_______________________________________________________________________
		! if ed is not a boundary edge --> calc depth integral: \nabla\int(U_n)dz 
		! for el(2)
		c2=0.0_WP
		if(el(2)>0) then
			deltaX2=edge_cross_dxdy(3,ed)
			deltaY2=edge_cross_dxdy(4,ed)
			do nz=1, nlevels(el(2))-1
				c2=c2-(UV(2,nz,el(2))*deltaX2-UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
			end do
		end if
		!_______________________________________________________________________
		ssh_rhs_old(enodes(1))=ssh_rhs_old(enodes(1))+(c1+c2)
		ssh_rhs_old(enodes(2))=ssh_rhs_old(enodes(2))-(c1+c2)
		
	end do
	
	!___________________________________________________________________________
	! take into account water flux
	if (.not. trim(which_ALE)=='linfs') then
		ssh_rhs_old(1:myDim_nod2D)=ssh_rhs_old(1:myDim_nod2D)-water_flux(1:myDim_nod2D)*area(1,1:myDim_nod2D)
	end if 
	
	!___________________________________________________________________________
	! update the thickness
	hbar_old=hbar
	hbar(1:myDim_nod2D)=hbar(1:myDim_nod2D)+ssh_rhs_old(1:myDim_nod2D)*dt/area(1,1:myDim_nod2D)
	call exchange_nod(hbar)
		
	!___________________________________________________________________________
	! fill the array for updating the stiffness matrix
	do elem=1,myDim_elem2D
		elnodes=elem2D_nodes(:,elem)
		dhe(elem)=sum(hbar(elnodes)-hbar_old(elnodes))/3.0_WP
	end do

end subroutine compute_hbar_ale
!
!
!===============================================================================
subroutine vert_vel_ale
	use g_config,only: dt, which_ALE
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_comm_auto
	implicit none
	
	integer       :: el(2), enodes(2), n, nz, ed
	real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2, dd, dd1, dddt 
	
	!___________________________________________________________________________
	! Contributions from levels in divergence
	Wvel=0.0_WP
	if (Fer_GM) then
         fer_Wvel(nz,n)=0.0_WP
	end if
	
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
		do nz=nlevels(el(1))-1,1,-1
			! --> h * u_vec * n_vec
			! --> e_vec = (dx,dy), n_vec = (-dy,dx);
			! --> h * u*(-dy) + v*dx
			c1=( UV(2,nz,el(1))*deltaX1 - UV(1,nz,el(1))*deltaY1 )*helem(nz,el(1))
			! inflow(outflow) "flux" to control volume of node enodes1
			Wvel(nz,enodes(1))=Wvel(nz,enodes(1))+c1
			! is equal to outflow(inflow) "flux" to control volume of node enodes2
			Wvel(nz,enodes(2))=Wvel(nz,enodes(2))-c1
			if (Fer_GM) then
				c1=(fer_UV(2,nz,el(1))*deltaX1- &
				fer_UV(1,nz,el(1))*deltaY1)*helem(nz,el(1))     
				fer_Wvel(nz,enodes(1))=fer_Wvel(nz,enodes(1))+c1
				fer_Wvel(nz,enodes(2))=fer_Wvel(nz,enodes(2))-c1
			end if  
		end do
		
		!_______________________________________________________________________
		! if ed is not a boundary edge --> calc div(u_vec*h) for every layer
		! for el(2)
		if(el(2)>0)then
			deltaX2=edge_cross_dxdy(3,ed)
			deltaY2=edge_cross_dxdy(4,ed)
			do nz=nlevels(el(2))-1,1,-1
				c2=-(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
				Wvel(nz,enodes(1))=Wvel(nz,enodes(1))+c2
				Wvel(nz,enodes(2))=Wvel(nz,enodes(2))-c2
				if (Fer_GM) then
					c2=-(fer_UV(2,nz,el(2))*deltaX2- &
					fer_UV(1,nz,el(2))*deltaY2)*helem(nz,el(2))
					fer_Wvel(nz,enodes(1))=fer_Wvel(nz,enodes(1))+c2
					fer_Wvel(nz,enodes(2))=fer_Wvel(nz,enodes(2))-c2
				end if  
			end do
		end if
	end do
	! |
	! |
	! +--> until here Wvel contains the thickness divergence div(u)
	
	!___________________________________________________________________________
	! cumulative summation of div(u_vec*h) vertically
	! W_k = W_k+1 - div(h_k*u_k)
	! W_k ... vertical flux trough 
	do n=1, myDim_nod2D
		do nz=nl-1,1,-1
			Wvel(nz,n)=Wvel(nz,n)+Wvel(nz+1,n)
			if (Fer_GM) then 
				fer_Wvel(nz,n)=fer_Wvel(nz,n)+fer_Wvel(nz+1,n)
			end if
		end do
	end do
	
	!___________________________________________________________________________
	! divide with depth dependent cell area to convert from Vertical flux to 
	! physical vertical velocities in units m/s
	do n=1, myDim_nod2D
		do nz=1,nlevels_nod2D(n)-1
			Wvel(nz,n)=Wvel(nz,n)/area(nz,n)
			if (Fer_GM) then 
				fer_Wvel(nz,n)=fer_Wvel(nz,n)/area(nz,n)          
			end if
		end do
	end do
	! |
	! |--> (A) linear free surface: dh/dt=0 ; W_t-W_b = -div(hu)
	! |
	! |--> (B) Full free surface: k=1 --> W_t-W_b = -div(hu) - dh/dt - Wflux_k1
	! |                           k>1 --> W_t-W_b = -div(hu)
	! |
	! |--> (C) ZSTAR:             W_t-W_b = -div(hu) - dh/dt - Wflux_k1
	!                                                    |
	!                                                    |--> (dh/dt)_k = 1/H*dh/dt
	
	
	!___________________________________________________________________________
	! Correct for free surface (zlevel and zstar)
	if(trim(which_ALE)=='zlevel') then
		! Update the upper level
		! water_flux is positive if out of the ocean
		! Wvel(1,n) should be 0 up to machine precision,
		! this place should be checked.
		do n=1, myDim_nod2D
			Wvel(1,n)=Wvel(1,n)-(hbar(n)-hbar_old(n))/dt-water_flux(n)
			hnode_new(1,n)=hnode(1,n)+hbar(n)-hbar_old(n)
		end do
		
	elseif (trim(which_ALE)=='zstar') then
		! distribute total change in ssh (hbar(n)-hbar_old(n)) over all layers 
		dO n=1, myDim_nod2D
			dd1=zbar(nlevels_nod2D_min(n)-1)
			
			! This is the depth the stretching is applied (area(nz,n)=area(1,n))
			dd=zbar(1)-dd1                      
			
			! how much of (hbar(n)-hbar_old(n)) is distributed into each layer
			! 1/H*dhbar
			dd=(hbar(n)-hbar_old(n))/dd
			
			! 1/H*dhbar/dt
			dddt=dd/dt
			dO nz=1,nlevels_nod2D_min(n)-2
				! why  *(zbar(nz)-dd1) ??? 
				! because here Wvel_k = SUM_k:kmax(div(h_k*v_k))/V_k
				! but Wvel_k = Wvel_k+1 - div(h_k*v_k) - h⁰_k/H*dhbar/dt
				! 				|--> Wvel_k+1 = Wvel_k+2 - div(h_k+1*v_k+1) - h⁰_k+1/H*dhbar/dt
				! 				                  |--> Wvel_k+2 = Wvel_k+3 - div(h_k+2*v_k+2) - h⁰_k+2/H*dhbar/dt
				!
				!     Wvel_k = SUM_i=k:kmax(div(h_i*v_i)) + 1/H*dhbar/dt*SUM_i=k:kmax(h⁰_k)
				! SUM_i=k:kmax(h⁰_k) == (zbar(nz)-dd1)
				Wvel(nz,n)=Wvel(nz,n)-dddt*(zbar(nz)-dd1) 
				hnode_new(nz,n)=(zbar(nz)-zbar(nz+1))*dd
! 				hnode_new(nz,n)=hnode(nz,n) + (zbar(nz)-zbar(nz+1))*dd ! ????
			end do
			Wvel(1,n)=Wvel(1,n)-water_flux(n) 
		end do
		! The implementation here is a bit strange, but this is to avoid 
		! unnecessary multiplications and divisions by area. We use the fact 
		! that we apply stretching only over the part of the column
		! where area(nz,n)=area(1,n)
	endif
	
	!___________________________________________________________________________
	call exchange_nod(Wvel) 
	call exchange_nod(hnode_new)   ! Or extend cycles above  
	if (Fer_GM) call exchange_nod(fer_Wvel)
	
	!___________________________________________________________________________
	! Split implicit vertical velocity onto implicit and explicit components
	if (w_split) then
		do n=1, myDim_nod2D+eDim_nod2D
			do nz=1,nlevels_nod2D(n)-1
				Wvel_e(nz,n)=min(max(Wvel(nz,n), -w_exp_max), w_exp_max)
			end do
		end do
	else
		Wvel_e=Wvel
	end if
	Wvel_i=Wvel-Wvel_e

end subroutine vert_vel_ale
!
!
!===============================================================================
subroutine solve_ssh_ale
use o_PARAM
use o_MESH
use o_ARRAYS
use g_PARSUP
use g_comm_auto
use g_config, only: which_ale
	!
	!
	!___USE PETSC SOLVER________________________________________________________
        ! this is not longer used but is still kept in the code
#ifdef PETSC
implicit none
#include "petscf.h"
integer                         :: myrows
integer                         :: Pmode
real(kind=8)                    :: rinfo(20,20)
integer                         :: maxiter=2000
integer                         :: restarts=15
integer                         :: fillin=3
integer                         :: lutype=2
integer                         :: nrhs=1
real(kind=8)                    :: droptol=1.e-7
real(kind=8)                    :: soltol =1e-10  !1.e-10
logical, save                   :: lfirst=.true.
real(kind=WP), allocatable      :: arr_nod2D(:),arr_nod2D2(:,:),arr_nod2D3(:)
real(kind=8)                    :: cssh1,cssh2,crhs
integer                         :: i
Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB +PET_REPORT + PET_QUIET+ PET_RCM+PET_PCBJ
if (lfirst) then   
   Pmode = Pmode+PET_STRUCT+PET_PMVALS + PET_PCASM+PET_OVL_2 !+PET_PCBJ+PET_ILU
   lfirst=.false.
end if
call PETSC_S(Pmode, 1, ssh_stiff%dim, ssh_stiff%nza, myrows, &
     maxiter,  & 
     restarts, &
     fillin,   &
     droptol,  &  
     soltol,   &
     part, ssh_stiff%rowptr, ssh_stiff%colind, ssh_stiff%values, &
     ssh_rhs, d_eta, &
     rinfo, MPI_COMM_WORLD)
	!
	!
	!___USE PARMS SOLVER (recommended)__________________________________________
#elif defined(PARMS)

  use iso_c_binding, only: C_INT, C_DOUBLE
  implicit none
#include "fparms.h"
logical, save        :: lfirst=.true.
integer(kind=C_INT)  :: ident
integer(kind=C_INT)  :: n3, reuse, new_values
integer(kind=C_INT)  :: maxiter, restart, lutype, fillin
real(kind=C_DOUBLE)  :: droptol, soltol
integer :: n

interface
   subroutine psolver_init(ident, SOL, PCGLOB, PCLOC, lutype, &
        fillin, droptol, maxiter, restart, soltol, &
        part, rowptr, colind, values, reuse, MPI_COMM) bind(C)
     use iso_c_binding, only: C_INT, C_DOUBLE
     integer(kind=C_INT) :: ident, SOL, PCGLOB, PCLOC, lutype, &
                            fillin,  maxiter, restart, &
                            part(*), rowptr(*), colind(*), reuse, MPI_COMM
     real(kind=C_DOUBLE) :: droptol,  soltol, values(*)
   end subroutine psolver_init
end interface
interface
   subroutine psolve(ident, ssh_rhs, values, d_eta, newvalues) bind(C)

     use iso_c_binding, only: C_INT, C_DOUBLE
     integer(kind=C_INT) :: ident, newvalues
     real(kind=C_DOUBLE) :: values(*), ssh_rhs(*), d_eta(*)

   end subroutine psolve
end interface

ident=1
maxiter=2000
restart=15
fillin=3
lutype=2
droptol=1.e-8
soltol=1.e-10

if  (trim(which_ale)=='linfs') then
    reuse=0
    new_values=0
else
    reuse=1      ! For varying coefficients, set reuse=1
    new_values=1 ! and new_values=1, as soon as the coefficients have changed
end if

! reuse=0: matrix remains static
! reuse=1: keeps a copy of the matrix structure to apply scaling of the matrix fast

! new_values=0: matrix coefficients unchanged (compared to the last call of psolve) 
! new_values=1: replaces the matrix values (keeps the structure and the preconditioner) 
! new_values=2: replaces the matrix values and recomputes the preconditioner (keeps the structure)

! new_values>0 requires reuse=1 in psolver_init!

if (lfirst) then
   ! Set SOLCG for CG solver (symmetric, positiv definit matrices only!!)
   !     SOLBICGS for BiCGstab solver (arbitrary matrices)
   ! call psolver_init(ident, SOLCG, PCRAS, PCILUK, lutype, &
   call psolver_init(ident, SOLBICGS, PCRAS, PCILUK, lutype, &
        fillin, droptol, maxiter, restart, soltol, &
        part-1, ssh_stiff%rowptr(:)-ssh_stiff%rowptr(1), &
        ssh_stiff%colind-1, ssh_stiff%values, reuse, MPI_COMM_WORLD)
   lfirst=.false.
end if

   call psolve(ident, ssh_rhs, ssh_stiff%values, d_eta, new_values)

#endif
	!
	!
	!___________________________________________________________________________
call exchange_nod(d_eta) !is this required after calling psolve ?
end subroutine solve_ssh_ale
!
!
!===============================================================================
subroutine impl_vert_visc_ale
USE o_MESH
USE o_PARAM
USE o_ARRAYS
USE g_PARSUP
USE g_CONFIG,only: dt
IMPLICIT NONE

real(kind=WP)              ::  a(nl-1), b(nl-1), c(nl-1), ur(nl-1), vr(nl-1)
real(kind=WP)              ::  cp(nl-1), up(nl-1), vp(nl-1)
integer                    ::  nz, elem, nzmax, elnodes(3)
real(kind=WP)              ::  zinv, m, friction, wu, wd
DO elem=1,myDim_elem2D
	elnodes=elem2D_nodes(:,elem)
	nzmax=nlevels(elem)
	
	!___________________________________________________________________________
	zbar_n=0.0_WP
	Z_n=0.0_WP
	zbar_n(nzmax)=zbar(nzmax)
	Z_n(nzmax-1)=zbar_n(nzmax) + helem(nzmax-1,elem)/2.0_WP
	do nz=nzmax-1,2,-1
		zbar_n(nz) = zbar_n(nz+1) + helem(nz,elem)
		Z_n(nz-1) = zbar_n(nz) + helem(nz-1,elem)/2.0_WP
	end do
	zbar_n(1) = zbar_n(2) + helem(1,elem)
	
	!___________________________________________________________________________
	! Operator
	! Regular part of coefficients:
	do nz=2, nzmax-2
		zinv=1.0_WP*dt/(zbar_n(nz)-zbar_n(nz+1))
		a(nz)=-Av(nz,elem)/(Z_n(nz-1)-Z_n(nz))*zinv
		c(nz)=-Av(nz+1,elem)/(Z_n(nz)-Z_n(nz+1))*zinv
		b(nz)=-a(nz)-c(nz)+1.0_WP
		! Update from the vertical advection
		wu=sum(Wvel_i(nz,   elnodes))/3.
		wd=sum(Wvel_i(nz+1, elnodes))/3.
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
	wu=sum(Wvel_i(nzmax-1, elnodes))/3.
	a(nzmax-1)=a(nzmax-1)+min(0._WP, wu)*zinv
	b(nzmax-1)=b(nzmax-1)+max(0._WP, wu)*zinv
	
	! The first row
	zinv=1.0_WP*dt/(zbar_n(1)-zbar_n(2))
	c(1)=-Av(2,elem)/(Z_n(1)-Z_n(2))*zinv
	a(1)=0.0_WP
	b(1)=-c(1)+1.0_WP
	! Update from the vertical advection
	wu=sum(Wvel_i(1, elnodes))/3.
	wd=sum(Wvel_i(2, elnodes))/3.
	
	b(1)=b(1)+wu*zinv
	b(1)=b(1)-min(0._WP, wd)*zinv
	c(1)=c(1)-max(0._WP, wd)*zinv
	! ===========================
	! The rhs:
	! ===========================
	ur(1:nzmax-1)=UV_rhs(1,1:nzmax-1,elem)
	vr(1:nzmax-1)=UV_rhs(2,1:nzmax-1,elem)
	! The first row contains surface forcing
	ur(1)= ur(1)+zinv*stress_surf(1,elem)/density_0
	vr(1)= vr(1)+zinv*stress_surf(2,elem)/density_0
	! The last row contains bottom friction
	zinv=1.0_WP*dt/(zbar_n(nzmax-1)-zbar_n(nzmax))
	friction=-C_d*sqrt(UV(1,nlevels(elem)-1,elem)**2+ &
				UV(2,nlevels(elem)-1,elem)**2)
	ur(nzmax-1)=ur(nzmax-1)+zinv*friction*UV(1,nzmax-1,elem)
	vr(nzmax-1)=vr(nzmax-1)+zinv*friction*UV(2,nzmax-1,elem)
	! Model solves for the difference to the timestep N and therefore we need to 
	! update the RHS for advective and diffusive contributions at the previous time step
	do nz=2, nzmax-2
		ur(nz)=ur(nz)-a(nz)*UV(1,nz-1,elem)-(b(nz)-1.0_WP)*UV(1,nz,elem)-c(nz)*UV(1,nz+1,elem)
		vr(nz)=vr(nz)-a(nz)*UV(2,nz-1,elem)-(b(nz)-1.0_WP)*UV(2,nz,elem)-c(nz)*UV(2,nz+1,elem)
	end do
	ur(1)=ur(1)-(b(1)-1.0_WP)*UV(1,1,elem)-c(1)*UV(1,2,elem)
	vr(1)=vr(1)-(b(1)-1.0_WP)*UV(2,1,elem)-c(1)*UV(2,2,elem)
	
	ur(nzmax-1)=ur(nzmax-1)-a(nzmax-1)*UV(1,nzmax-2,elem)-(b(nzmax-1)-1.0_WP)*UV(1,nzmax-1,elem)
	vr(nzmax-1)=vr(nzmax-1)-a(nzmax-1)*UV(2,nzmax-2,elem)-(b(nzmax-1)-1.0_WP)*UV(2,nzmax-1,elem)
	! ===========================
	! The sweep algorithm
	! ===========================
	! initialize c-prime and s,t-prime
	cp(1) = c(1)/b(1)
	up(1) = ur(1)/b(1)
	vp(1) = vr(1)/b(1)
	! solve for vectors c-prime and t, s-prime
	do nz = 2,nzmax-1
		m = b(nz)-cp(nz-1)*a(nz)
		cp(nz) = c(nz)/m
		up(nz) = (ur(nz)-up(nz-1)*a(nz))/m
		vp(nz) = (vr(nz)-vp(nz-1)*a(nz))/m
	enddo
	! initialize x
	ur(nzmax-1) = up(nzmax-1)
	vr(nzmax-1) = vp(nzmax-1)
	! solve for x from the vectors c-prime and d-prime
	do nz = nzmax-2, 1, -1
		ur(nz) = up(nz)-cp(nz)*ur(nz+1)
		vr(nz) = vp(nz)-cp(nz)*vr(nz+1)
	end do
	! ===========================
	! RHS update
	! ===========================
	do nz=1,nzmax-1
		UV_rhs(1,nz,elem)=ur(nz)
		UV_rhs(2,nz,elem)=vr(nz)
	end do
	
end do   !!! cycle over elements
end subroutine impl_vert_visc_ale

!
!
!===============================================================================
subroutine oce_timestep_ale(n)
	use g_config, only: logfile_outfreq,dt
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_comm_auto
	
	use ieee_arithmetic        !??????????????????????????????
	
	IMPLICIT NONE
	real(kind=8)      :: t1, t2, t3, t4, t5, t6, t7, t8, t9
	integer           :: n
	integer           :: el,nz,elnodes(3) !??????????????????????????????
	CHARACTER(LEN=80) :: fmt   !??????????????????????????????
	real(kind=8)      :: global_vol_eta ,local_vol_eta ,global_max_eta ,global_min_eta ,local_max_eta ,local_min_eta
	real(kind=8)      :: global_vol_hbar,local_vol_hbar,global_max_hbar,global_min_hbar,local_max_hbar,local_min_hbar
	real(kind=8)      :: global_vol_wflux,local_vol_wflux,global_max_wflux,global_min_wflux,local_max_wflux,local_min_wflux
	real(kind=8)      :: global_vol_hflux,local_vol_hflux,global_max_hflux,global_min_hflux,local_max_hflux,local_min_hflux
!DS 	real(kind=8)      :: global_d_eta, global_wflux, local_d_eta, local_wflux, local_vol, global_vol
	real(kind=8),allocatable  :: aux1(:),aux2(:)

	t1=MPI_Wtime()
	
! 	heat_flux=0.0_WP
! 	water_flux=0.0_WP

	!___________________________________________________________________________
	call pressure_bv               !!!!! HeRE change is made. It is linear EoS now.
	
	!___________________________________________________________________________
	call status_check
	
	!___________________________________________________________________________
	call oce_mixing_PP
	
	!___________________________________________________________________________
	! Current dynamic elevation alpha*hbar(n+1/2)+(1-alpha)*hbar(n-1/2)
	! equation (7) Danlov et.al "the finite volume sea ice ocean model FESOM2
!DS	eta_n=alpha*hbar+(1.0_WP-alpha)*hbar_old
	
	!___________________________________________________________________________
	if(mom_adv/=3) then
		call compute_vel_rhs
	else
		call compute_vel_rhs_vinv
	end if
	
	!___________________________________________________________________________
	if (tau_c > 1.e-12) call viscosity_filt2x
	
	!___________________________________________________________________________
	if(i_vert_visc) call impl_vert_visc_ale
	t2=MPI_Wtime()
	
	!___________________________________________________________________________
	! >->->->->->->->->->->->->     ALE-part starts     <-<-<-<-<-<-<-<-<-<-<-<-
	!___________________________________________________________________________
	! Update stiffness matrix by dhe=hbar(n+1/2)-hbar(n-1/2) on elements
	call stiff_mat_update 
	
	! ssh_rhs=-alpha*\nabla\int(U_n+U_rhs)dz-(1-alpha)*...
	! see "FESOM2: from finite elements to finte volumes, S. Danilov..." eq. (12) rhs
	call compute_ssh_rhs_ale
	
	! Take updated ssh matrix and solve --> new ssh!
	call solve_ssh_ale
	t3=MPI_Wtime() 
	
	! estimate new horizontal velocity u^(n+1)
	! u^(n+1) = u* + [-g * tau * theta * grad(eta^(n+1)-eta^(n)) ]
	call update_vel
	t4=MPI_Wtime() 
	
	! Update to hbar(n+3/2) and compute dhe to be used on the next step
	call compute_hbar_ale
	t5=MPI_Wtime() 
	
	!---------------------------------------------------------------------------
	! Does not belong directly to ALE formalism: Ferrari, Gent, McWiliams parameterisation
	if (Fer_GM) then
		call fer_compute_C_K
		call compute_sigma_xy(tr_arr(:,:,1),tr_arr(:,:,2))
		call fer_solve_Gamma
		call fer_gamma2vel
	end if
	t6=MPI_Wtime() 
	!---------------------------------------------------------------------------
	
	! The main step of ALE procedure --> this is were the magic happens --> here 
	! is decided how change in hbar is distributed over the vertical layers
	call vert_vel_ale 
	t7=MPI_Wtime() 


!if (mstep==3) then
!if (mype==0) then
!write(*,*) 'alpha==', alpha
!do n=1, myDim_nod2D
!   write(*,*) hbar(n)-hbar_old(n), d_eta(n), wvel(1,n)*dt
!end do
!end if
!call par_ex
!stop
!end if
	
	if (mstep==3) then
		if (mype==0) then
			write(*,*) 'alpha==', alpha
			do n=1, myDim_nod2D
				write(*,*) hbar(n)-hbar_old(n), d_eta(n), Wvel(1,n)*dt
			end do
		end if
		call par_ex
		stop
	end if 
	
	! solve tracer equation 
	call solve_tracers_ale
	t8=MPI_Wtime() 
	
	! Update hnode=hnode_new, helem
	call update_thickness_ale  
	t9=MPI_Wtime() 
	
	!___________________________________________________________________________
	! write out execution times for ocean step parts
	if(mod(n,logfile_outfreq)==0 .and. mype==0) then  
		write(*,*) '___ALE OCEAN STEP EXECUTION TIMES______________________'
		write(*,"(A, ES10.3)") '	Oce. Dynamics    :', t2-t1
		write(*,"(A, ES10.3)") '	Oce. Update Vel. :', t4-t3
		write(*,"(A, ES10.3)") '	Oce. Fer-GM.     :', t6-t5
		write(*,*) '   _______________________________'
		write(*,"(A, ES10.3)") '	ALE-Solve SSH    :', t3-t2
		write(*,"(A, ES10.3)") '	ALE-Calc. hbar   :', t5-t4
		write(*,"(A, ES10.3)") '	ALE-Update+W     :', t7-t6
		write(*,"(A, ES10.3)") '	ALE-Solve Tracer :', t8-t7
		write(*,"(A, ES10.3)") '	ALE-Update hnode :', t9-t8
		write(*,*) '   _______________________________'
		write(*,"(A, ES10.3)") '	Oce. TOTAL       :', t9-t1
		write(*,*)
! 		write(*,*) 'average_time (horizontal diffusion) = ', time_sum/n
	end if
	
	!___________________________________________________________________________
	! write out field estimates
! 	if(mod(n,logfile_outfreq)==0 .and. mype==0) then
! 		fmt="(A, ES10.3, A, ES10.3)"
! 		write(*,*) '	___ALE OCEAN STEP: FIELDS______________________________'
! 		write(*,fmt) '	min(hnode) = ',minval(hnode,MASK=(hnode/=0.0)),', max(hnode) = ', maxval(hnode,MASK=(hnode/=0.0))
! 		write(*,fmt) '	min(helem) = ',minval(helem,MASK=(helem/=0.0))          ,', max(helem) = ', maxval(helem,MASK=(helem/=0.0))
! 		write(*,fmt) '	min(hbar)  = ',minval(hbar)           ,', max(hbar)  = ', maxval(hbar)
! 		write(*,fmt) '	min(eta_n) = ',minval(eta_n)          ,', max(eta_n) = ', maxval(eta_n)
! 		write(*,fmt) '	min(U)     = ',minval(UV(1,:,:))      ,', max(U)     = ', maxval(UV(1,:,:))
! 		write(*,fmt) '	min(V)     = ',minval(UV(2,:,:))      ,', max(V)     = ', maxval(UV(2,:,:))
! 		write(*,fmt) '	min(Wvel)  = ',minval(Wvel(:,:))      ,', max(Wvel)  = ', maxval(Wvel(:,:))
! 		write(*,fmt) '	min(temp)  = ',minval(tr_arr(:,:,1))  ,', max(temp)  = ', maxval(tr_arr(:,:,1))
! 		write(*,fmt) '	min(salt)  = ',minval(tr_arr(:,:,2))  ,', max(salt)  = ', maxval(tr_arr(:,:,2))
! 		write(*,fmt) '	min(hflux) = ',minval(heat_flux)      ,'  max(hflux) = ', maxval(heat_flux)
! 		write(*,fmt) '	min(wflux) = ',minval(water_flux)     ,'  max(wflux) = ', maxval(water_flux)
! 		write(*,*)
! 	end if
! #ifdef FALSE	
	!___________________________________________________________________________
	!check for different parameter
	if (mod(n,logfile_outfreq)==0) then
		local_vol_eta=0.0_WP
		local_vol_hbar=0.0_WP
		local_vol_wflux=0.0_WP
		local_vol_hflux=0.0_WP
!DS		local_vol=0.0_WP
!DS		local_d_eta=0.0_WP
!DS		local_wflux=0.0_WP
		do el=1, myDim_elem2D
			elnodes=elem2D_nodes(:, el)
			local_vol_eta=local_vol_eta+sum(eta_n(elnodes))/3.0_WP !*elem_area(el)
			local_vol_hbar=local_vol_hbar+sum(hbar(elnodes))/3.0_WP !*elem_area(el)
			local_vol_wflux=local_vol_wflux+sum(water_flux(elnodes))/3.0_WP !*elem_area(el)
			local_vol_hflux=local_vol_hflux+sum(heat_flux(elnodes))/3.0_WP !*elem_area(el)
!DS 			local_vol=local_vol+elem_area(el)
!DS 			local_d_eta=local_d_eta+sum(d_eta(elnodes))/3.0_WP*elem_area(el)
!DS 			local_wflux=local_wflux+sum(water_flux(elnodes))/3.0_WP*elem_area(el)
		end do
		
		local_vol_eta  = local_vol_eta/myDim_elem2D
		local_max_eta  = maxval(eta_n)
		local_min_eta  = minval(eta_n)
		global_vol_eta=0.0_WP
		global_max_eta=0.0_WP
		global_min_eta=0.0_WP
		call MPI_AllREDUCE(local_vol_eta, global_vol_eta, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
							MPI_COMM_WORLD, MPIerr)
		global_vol_eta = global_vol_eta/npes
		call MPI_AllREDUCE(local_max_eta, global_max_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
							MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(local_min_eta, global_min_eta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
							MPI_COMM_WORLD, MPIerr)
							
		local_vol_hbar = local_vol_hbar/myDim_elem2D
		local_max_hbar = maxval(hbar)
		local_min_hbar = minval(hbar)
		global_vol_hbar=0.0_WP
		global_max_hbar=0.0_WP
		global_min_hbar=0.0_WP
		call MPI_AllREDUCE(local_vol_hbar, global_vol_hbar, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
							MPI_COMM_WORLD, MPIerr)
		global_vol_hbar=global_vol_hbar/npes
		call MPI_AllREDUCE(local_max_hbar, global_max_hbar, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
							MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(local_min_hbar, global_min_hbar, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
							MPI_COMM_WORLD, MPIerr)
							
		local_vol_wflux = local_vol_wflux/myDim_elem2D
		local_max_wflux = maxval(water_flux)
		local_min_wflux = minval(water_flux)
		global_vol_wflux=0.0_WP
		global_max_wflux=0.0_WP
		global_min_wflux=0.0_WP
		call MPI_AllREDUCE(local_vol_wflux, global_vol_wflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
							MPI_COMM_WORLD, MPIerr)
		global_vol_wflux=global_vol_wflux/npes
		call MPI_AllREDUCE(local_max_wflux, global_max_wflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
							MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(local_min_wflux, global_min_wflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
							MPI_COMM_WORLD, MPIerr)
							
		local_vol_hflux = local_vol_hflux/myDim_elem2D
		local_max_hflux = maxval(heat_flux)
		local_min_hflux = minval(heat_flux)
		global_vol_hflux=0.0_WP
		global_max_hflux=0.0_WP
		global_min_hflux=0.0_WP
		call MPI_AllREDUCE(local_vol_hflux, global_vol_hflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
							MPI_COMM_WORLD, MPIerr)
		global_vol_hflux=global_vol_hflux/npes
		call MPI_AllREDUCE(local_max_hflux, global_max_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
							MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(local_min_hflux, global_min_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
							MPI_COMM_WORLD, MPIerr)

!DS             ! integrate d_eta and ssh_rhs over the ocean surface and compute their means
!DS 		global_vol   = 0.0_WP
!DS 		global_d_eta = 0.0_WP
!DS 		global_wflux = 0.0_WP
!DS 		call MPI_AllREDUCE(local_vol, global_vol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
!DS 							MPI_COMM_WORLD, MPIerr)
!DS 		call MPI_AllREDUCE(local_d_eta, global_d_eta, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
!DS 							MPI_COMM_WORLD, MPIerr)
!DS 		call MPI_AllREDUCE(local_wflux, global_wflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
!DS 							MPI_COMM_WORLD, MPIerr)
!DS 		global_wflux=global_wflux/global_vol
!DS 		global_d_eta=global_d_eta/global_vol
!DS             if (mype==0) then
!DS                write(*,*) 'int. d_eta and water_flux =', global_d_eta, global_wflux
!DS             end if
							
		if (mod(n,logfile_outfreq)==0 .and. mype==0) then
			write(*,*) '	___global max/min/mean --> mstep=',mstep,'_________'
			write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") '     eta_n : ',global_max_eta,' | ',global_min_eta,' | ',global_vol_eta
			write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") '     hbar  : ',global_max_hbar,' | ',global_min_hbar,' | ',global_vol_hbar
			write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") '     wflux : ',global_max_wflux,' | ',global_min_wflux,' | ',global_vol_wflux
			write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") '     hflux : ',global_max_hflux,' | ',global_min_hflux,' | ',global_vol_hflux
		endif
		
		
		global_max_hflux=0.0_WP
		global_min_hflux=0.0_WP
		call MPI_AllREDUCE(maxval(tr_arr(:,:,1)), global_max_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(minval(tr_arr(:,:,1)), global_min_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, MPIerr)
		if (mod(n,logfile_outfreq)==0 .and. mype==0) then
			write(*,*) '    temp  : ',global_max_hflux,' | ',global_min_hflux
		end if
		global_max_hflux=0.0_WP
		global_min_hflux=0.0_WP
		call MPI_AllREDUCE(maxval(tr_arr(:,:,2),MASK=(tr_arr(:,:,2)/=0.0)), global_max_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(minval(tr_arr(:,:,2),MASK=(tr_arr(:,:,2)/=0.0)), global_min_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, MPIerr)
		if (mod(n,logfile_outfreq)==0 .and. mype==0) then
			write(*,*) '    salt  : ',global_max_hflux,' | ',global_min_hflux
		end if
		global_max_hflux=0.0_WP
		global_min_hflux=0.0_WP
		call MPI_AllREDUCE(maxval(Wvel(1,:)), global_max_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(minval(Wvel(1,:)), global_min_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, MPIerr)
		if (mod(n,logfile_outfreq)==0 .and. mype==0) then
			write(*,*) '    Wvel(:,1): ',global_max_hflux,' | ',global_min_hflux
		end if
! 		global_max_hflux=0.0_WP
! 		global_min_hflux=0.0_WP
! 		call MPI_AllREDUCE(maxval(hnode,MASK=(hnode/=0.0)), global_max_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, MPIerr)
! 		call MPI_AllREDUCE(minval(hnode,MASK=(hnode/=0.0)), global_min_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, MPIerr)
! 		if (mod(n,logfile_outfreq)==0 .and. mype==0) then
! 			write(*,"(A, ES10.3, A, ES10.3)") '     hnode : ',global_max_hflux,' | ',global_min_hflux
! 		end if
		
		!___________________________________________________________________________
		! check variables for strange values
		do el=1,myDim_nod2d
		
			! check SSH
			if ( (ieee_is_nan(eta_n(el)) .or. &
				eta_n(el)<-10.0 .or. eta_n(el)>10.0)) then
				write(*,*) ' --STOP--> found eta_n become NaN or <-10.0, >10.0'
				write(*,*) 'mype        = ',mype
				write(*,*) 'mstep       = ',n
				write(*,*) 'node        = ',el
				write(*,*) 'eta_n(el)   = ',eta_n(el)
				write(*,*) ' lon,lat    = ',geo_coord_nod2D(:,el)/rad
				call par_ex(1)
			endif
			
			! check HBAR
			if ( (ieee_is_nan(hbar(el)) .or. &
				hbar(el)<-10.0 .or. hbar(el)>10.0)) then
				write(*,*) ' --STOP--> found hbar become NaN or <-10.0, >10.0'
				write(*,*) 'mype        = ',mype
				write(*,*) 'mstep       = ',n
				write(*,*) 'node        = ',el
				write(*,*) 'hbar(el)    = ',hbar(el)
				write(*,*) ' lon,lat    = ',geo_coord_nod2D(:,el)/rad
				call par_ex(1)
			endif
			
			
			
			do nz=1,nlevels_nod2D(el)-1
				! check TEMP
				if ( ieee_is_nan(tr_arr(nz, el,1)) .or. &
					tr_arr(nz, el,1) < -5.0 .or. tr_arr(nz, el,1)>100) then
					write(*,*) ' --STOP--> found temperture become NaN or <-3.0, >100'
					write(*,*) 'mype        = ',mype
					write(*,*) 'mstep       = ',n
					write(*,*) 'node        = ',el
					write(*,*) 'nz          = ',nz
					write(*,*)
					write(*,*) 'temp(nz, el)= ',tr_arr(nz, el,1)
					write(*,*) 'temp(: , el)= ',tr_arr(:, el,1)
					write(*,*)
					write(*,*) 'hflux       = ',heat_flux(el)
					write(*,*)
					write(*,*) 'eta_n       = ',eta_n(el)
					write(*,*)
					write(*,*) 'hnode_new   = ',hnode_new(:,el)
					write(*,*)
					write(*,*) 'Kv          = ',Kv(:,el)
					write(*,*)
					write(*,*) 'Wvel        = ',Wvel(:,el)
					write(*,*)
					write(*,*) ' lon,lat    = ',geo_coord_nod2D(:,el)/rad
					call output (1,n)        ! save (NetCDF)
					call restart(1,n)        ! save (NetCDF)
					call par_ex(1)
				endif
				
				! check SALT
				if ( ieee_is_nan(tr_arr(nz, el,2)) .or.  &
					tr_arr(nz, el,2) < 0 .or. tr_arr(nz, el,2)>50 ) then
					
					write(*,*) ' --STOP--> found salinity become NaN or <0, >50'
					write(*,*) 'mype        = ',mype
					write(*,*) 'mstep       = ',n
					write(*,*) 'node        = ',el
					write(*,*) 'nz          = ',nz
					write(*,*)
					write(*,*) 'salt(nz, el)= ',tr_arr(nz, el,2)
					write(*,*) 'salt(: , el)= ',tr_arr(:, el,2)
					write(*,*)
					write(*,*) 'temp(nz, el)= ',tr_arr(nz, el,1)
					write(*,*) 'temp(: , el)= ',tr_arr(:, el,1)
					write(*,*)
					write(*,*) 'wflux       = ',water_flux(el)
					write(*,*)
					write(*,*) 'eta_n       = ',eta_n(el)
					write(*,*)
					write(*,*) 'hnode_new   = ',hnode_new(:,el)
					write(*,*)
					write(*,*) 'Kv          = ',Kv(:,el)
					write(*,*)
					write(*,*) 'Wvel        = ',Wvel(:,el)
					write(*,*)
					write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,el)/rad
					write(*,*)
! 					call output (0,n)        ! save (NetCDF)
! 					call restart(1,n)        ! save (NetCDF)
! 					pause(30)
					call par_ex(1)
				endif 
			end do
		end do
	end if !-->if (mod(n,logfile_outfreq)==0) then
! #endif	
end subroutine oce_timestep_ale
