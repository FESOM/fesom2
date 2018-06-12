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
	
	integer             :: n, nzmax

	
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
	allocate(zbar_3d_n(nl,myDim_nod2D+eDim_nod2D))
	
	! Z_n: mid depth of layers due to ale thinkness variactions at ervery node n 
	allocate(Z_n(nl-1))
	allocate(Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D)) 
	
	! bottom_elem_tickness: changed bottom layer thinkness due to partial cells
	allocate(bottom_elem_thickness(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
	allocate(zbar_e_bot(myDim_elem2D+eDim_elem2D+eXDim_elem2D)) 
	
	! also change bottom thickness at nodes due to partial cell --> bottom 
	! thickness at nodes is the volume weighted mean of sorounding elemental 
	! thicknesses
	allocate(bottom_node_thickness(myDim_nod2D+eDim_nod2D))
	allocate(zbar_n_bot(myDim_nod2D+eDim_nod2D)) 
	
	!___initialize______________________________________________________________
	hbar=0.0_WP
	hbar_old=hbar
	dhe=0.0_WP
	
	! calculate thickness of partial bottom layer cells
	zbar_n_bot = 0.0
	zbar_e_bot = 0.0
	call init_bottom_elem_thickness
	
	! initialise 3d field of depth levels and mid-depth levels
	zbar_3d_n=0.0_WP
	Z_3d_n   =0.0_WP
	do n=1,myDim_nod2D+eDim_nod2D 
		! max. number of levels at node n
		nzmax=nlevels_nod2D(n)
		
		zbar_3d_n(1:nzmax-1,n)=zbar(1:nzmax-1);
		! in case of partial cells bottom depth is different from zbar(nzmax)
		zbar_3d_n(nzmax,n)=zbar_n_bot(n);
		 
		Z_3d_n(1:nzmax-2,n) =Z(1:nzmax-2);
		! in case of partial cells bottom mid depth is different from Z(nzmax-1)
		Z_3d_n(nzmax-1,n) =zbar_3d_n(nzmax-1,n)+(zbar_n_bot(n)-zbar_3d_n(nzmax-1,n))/2;
		
	end do
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
	use g_comm_auto
	implicit none
	
	integer :: elem, elnodes(3), nle, nln, n, k
	real(kind=WP) :: dd, dd1, hnbot, tvol 
	
	!___________________________________________________________________________
	! If we use partial cells, the thickness of bottom cell is adjusted.
	! The adjustment is limited. It cannot be more than + (1/2) of the deeper
	! layer, nor -(1/2) of the current layer. 
	if(use_partial_cell) then 
		!Adjust the thickness of elemental bottom cells
		do elem=1, myDim_elem2D
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
					zbar_e_bot(elem) = dd
					
				else
					! case 1 : max(Z(nle),dd) = dd
					! case 2 : max(Z(nle),dd) = Z(nle)
					zbar_e_bot(elem) = max(Z(nle),dd)
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
				zbar_e_bot(elem) = min(Z(nle-1),dd)
! 				zbar_e_bot(elem) = zbar(nlevels(elem)) ! partial cells are not allowed to become smaller
				
			end if	    
			dd1=zbar(nle-1)-zbar_e_bot(elem)
			bottom_elem_thickness(elem)=dd1
			
		end do
		
		!_______________________________________________________________________
		! calculate bottom node thickness from weighted mean of sorounding elemental
		! bottom thicknesses
		call exchange_elem(bottom_elem_thickness)
		do n=1,myDim_nod2D+eDim_nod2D
			hnbot= 0.0_WP
			tvol = 0.0_WP
			do k=1, nod_in_elem2D_num(n)
				elem=nod_in_elem2D(k,n)
				tvol=tvol+elem_area(elem)
				hnbot = hnbot + bottom_elem_thickness(elem)*elem_area(elem)
			end do
			bottom_node_thickness(n) = hnbot/tvol
			zbar_n_bot(n)			 = zbar(nlevels_nod2D(n)-1)-bottom_node_thickness(n)
		end do 
		
	else
		do elem=1, myDim_elem2D
			nle=nlevels(elem)
			bottom_elem_thickness(elem)=zbar(nle-1)-zbar(nle)
			zbar_e_bot(elem) = zbar(nle)
		end do
		
		do n=1,myDim_nod2D+eDim_nod2D
			nln = nlevels_nod2D(n)
			bottom_node_thickness(n)=zbar(nln-1)-zbar(nln)
			zbar_n_bot(n) = zbar(nln)
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
	
	if(mype==0) then
		write(*,*) '____________________________________________________________'
		write(*,*) ' --> initialise ALE layerthicknesses, depth levels and middepth levels'
		write(*,*)
	end if
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
! 			do nz=2,nlevels_nod2D(n)-1
			do nz=2,nlevels_nod2D(n)-2
				hnode(nz,n)=(zbar(nz)-zbar(nz+1))
			end do      
			
			! set bottom node thickness
			hnode(nlevels_nod2D(n)-1,n)=bottom_node_thickness(n)
			
			do nz=nlevels_nod2D(n),nl-1
				hnode(nz,n)=0.0_WP
			end do
		end do
		do elem=1,myDim_elem2D
			dhe(elem)=0.0_WP
			do nz=1,nlevels(elem)-2
				helem(nz,elem)=(zbar(nz)-zbar(nz+1))
			end do
			
			! set bottom elem thickness
			helem(nlevels(elem)-1,elem)=bottom_elem_thickness(elem)
			
			Do nz=nlevels(elem),nl-1
				helem(nz,elem)=0.0_WP
			end do
		end do
		
	elseif (trim(which_ale)=='zlevel') then
		!_______________________________________________________________________
		! >->->->->->->->->->->->->->->     z-level    <-<-<-<-<-<-<-<-<-<-<-<-<
		!_______________________________________________________________________
		! --> include all ssh variations into the top layer 
		do n=1,myDim_nod2D+eDim_nod2D
			
			! put all ssh variation (hbar) into first layer 
			hnode(1,n)=hbar(n)+(zbar(1)-zbar(2))
			
			! leave lower levels untouched
! 			do nz=2,nlevels_nod2D(n)-1
			do nz=2,nlevels_nod2D(n)-2
				hnode(nz,n)=(zbar(nz)-zbar(nz+1))
			end do 
			
			! set bottom node thickness
			hnode(nlevels_nod2D(n)-1,n)=bottom_node_thickness(n)
			
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
			
			! do not distribute hbar into cells that intersect somehow with bottom layer 
			do nz=nlevels_nod2D_min(n)-1, nlevels_nod2D(n)-1
				hnode(nz,n)=(zbar(nz)-zbar(nz+1))
			end do
			
			! set bottom node thickness
			hnode(nlevels_nod2D(n)-1,n)=bottom_node_thickness(n)
			
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
	use g_config,only: which_ale,lzstar_lev,min_hnode
	implicit none
	integer :: n, nz, elem, elnodes(3),nzmax
	integer      , dimension(:), allocatable :: idx
	!___________________________________________________________________________
	! >->->->->->->->->->->->->->->       z-level      <-<-<-<-<-<-<-<-<-<-<-<-<
	!___________________________________________________________________________
	if     (trim(which_ale)=='zlevel') then
		
		!_______________________________________________________________________
		! idx is only needed for local star case to estimate over how much 
		! depth layers hnode, depthlevel and mid-depthlevel need to be updated
		allocate(idx(lzstar_lev))
		idx = (/(nz,nz=1,lzstar_lev,1)/)
		
		!_______________________________________________________________________
		do elem=1,myDim_elem2D
			elnodes=elem2D_nodes(:, elem)
			!___________________________________________________________________
			! actualize elemental layer thinkness in first lzstar_lev layers
			if (any(hnode_new(2:lzstar_lev,elnodes(1))-hnode(2:lzstar_lev,elnodes(1))/=0.0_WP)     .or. &
				any(hnode_new(2:lzstar_lev,elnodes(2))-hnode(2:lzstar_lev,elnodes(2))/=0.0_WP) .or. &
				any(hnode_new(2:lzstar_lev,elnodes(3))-hnode(2:lzstar_lev,elnodes(3))/=0.0_WP)      &
				) then
				! --> case local zstar
				! try to limitate over how much layers i realy need to distribute
				! the change in ssh, so that the next loops run only over the 
				! nesseccary levels and not over all lzstar_lev levels
				nz    = max(1 ,maxval(pack(idx,hnode_new(1:lzstar_lev,elnodes(1))-hnode(1:lzstar_lev,elnodes(1))/=0.0_WP)))
				nz    = max(nz,maxval(pack(idx,hnode_new(1:lzstar_lev,elnodes(2))-hnode(1:lzstar_lev,elnodes(2))/=0.0_WP)))
				nz    = max(nz,maxval(pack(idx,hnode_new(1:lzstar_lev,elnodes(3))-hnode(1:lzstar_lev,elnodes(3))/=0.0_WP)))
				nzmax = min(nz,nlevels(elem)-2)
				do nz=1,nzmax
					helem(nz,elem)=sum(hnode_new(nz,elnodes))/3.0_WP
				end do
			!___________________________________________________________________
			! only actualize elemental layer thickness in first layer 
			else
				! --> case normal zlevel
				helem(1,elem)=sum(hnode_new(1,elnodes))/3.0_WP
			end if
		end do
		
		!_______________________________________________________________________
		do n=1,myDim_nod2D+eDim_nod2D
			!___________________________________________________________________
			! actualize layer thinkness in first lzstar_lev layers
			if ( (any(hnode_new(2:lzstar_lev,n)-hnode(2:lzstar_lev,n)/=0.0_WP)) ) then
				! --> case local zstar 
				! try to limitate over how much layers i realy need to distribute
				! the change in ssh, so that the next loops run only over the 
				! nesseccary levels and not over all lzstar_lev levels
				nz = max(1,maxval(pack(idx,hnode_new(1:lzstar_lev,n)-hnode(1:lzstar_lev,n)/=0.0_WP)))
				
				! nlevels_nod2D_min(n)-1 ...would be hnode of partial bottom cell but this
				! one is not allowed to change so go until nlevels_nod2D_min(n)-2
				nzmax = min(nz,nlevels_nod2D_min(n)-2)
				! do not touch zbars_3d_n that are involved in the bottom cell !!!!
				! this ones are set up during initialisation and are not touched afterwards
				! --> nlevels_nod2D_min(n),nlevels_nod2D_min(n)-1
				do nz=nzmax,1,-1
					hnode(nz,n)     = hnode_new(nz,n)
					zbar_3d_n(nz,n) = zbar_3d_n(nz+1,n)+hnode_new(nz,n)
					Z_3d_n(nz,n)    = zbar_3d_n(nz+1,n)+hnode_new(nz,n)/2.0_WP
				end do
			!___________________________________________________________________
			! only actualize layer thinkness in first layer 
			else
				! --> case normal zlevel
				hnode(1,n)    = hnode_new(1,n)
				zbar_3d_n(1,n)= zbar_3d_n(2,n)+hnode_new(1,n)
				Z_3d_n(1,n)   = zbar_3d_n(2,n)+hnode_new(1,n)/2.0_WP
			end if
		end do
		
		!_______________________________________________________________________
		deallocate(idx)
		
	!___________________________________________________________________________
	! >->->->->->->->->->->->->->->       z-star       <-<-<-<-<-<-<-<-<-<-<-<-<
	!___________________________________________________________________________
	elseif (trim(which_ale)=='zstar' ) then
		
		! --> update layer thinkness, depth layer  and mid-depth layer at node
		do n=1, myDim_nod2D+eDim_nod2D
			! actualize 3d depth levels and mid-depth levels from bottom to top
			nzmax = nlevels_nod2D_min(n)-2
			! do not touch zbars_3d_n that are involved in the bottom cell !!!!
			! --> nlevels_nod2D_min(n),nlevels_nod2D_min(n)-1
			do nz=nzmax,1,-1
				hnode(nz,n)     = hnode_new(nz,n)
				zbar_3d_n(nz,n) = zbar_3d_n(nz+1,n) + hnode_new(nz,n)
				Z_3d_n(nz,n)    = zbar_3d_n(nz+1,n) + hnode_new(nz,n)/2.0_WP
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
!===============================================================================
! update thickness arrays based on the current hbar 
subroutine restart_thickness_ale
	use o_PARAM
	use o_MESH
	use g_PARSUP
	use o_ARRAYS
	use g_config,only: which_ale,lzstar_lev,min_hnode
	implicit none
	integer :: n, nz, elem, elnodes(3),nzmax
	integer      , dimension(:), allocatable :: idx
	
	if(mype==0) then
		write(*,*) '____________________________________________________________'
		write(*,*) ' --> restart ALE layerthicknesses, depth levels and middepth levels'
		write(*,*)
	end if
	!___________________________________________________________________________
	! >->->->->->->->->->->->->->->       z-level      <-<-<-<-<-<-<-<-<-<-<-<-<
	!___________________________________________________________________________
	if     (trim(which_ale)=='zlevel') then
		!_______________________________________________________________________
		! idx is only needed for local star case to estimate over how much 
		! depth layers hnode, depthlevel and mid-depthlevel need to be updated
		allocate(idx(lzstar_lev))
		idx = (/(nz,nz=1,lzstar_lev,1)/)
		
		! restart depthlevels (zbar_3d_n) and mitdpethlevels (Z_3d_n)
		do n=1,myDim_nod2D+eDim_nod2D
			if (any(hnode(2:lzstar_lev,n) /=  &
					(zbar(2:lzstar_lev)-zbar(3:lzstar_lev+1))) ) then
				! --> case local zstar 
				! the change in ssh, so that the next loops run only over the 
				! nesseccary levels and not over all lzstar_lev levels
				nz  = max(1,maxval(pack(idx,hnode(1:lzstar_lev,n)/=(zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1)))))
				
				! nlevels_nod2D_min(n)-1 ...would be hnode of partial bottom cell but this
				! one is not allowed to change so go until nlevels_nod2D_min(n)-2
				nzmax = min(nz,nlevels_nod2D_min(n)-2)
				do nz=nzmax,1,-1
					zbar_3d_n(nz,n) = zbar_3d_n(nz+1,n)+hnode(nz,n)
					Z_3d_n(nz,n)    = zbar_3d_n(nz+1,n)+hnode(nz,n)/2.0_WP
				end do
			else
				! --> case normal zlevel
				zbar_3d_n(1,n)= zbar_3d_n(2,n)+hnode(1,n)
				Z_3d_n(1,n)   = zbar_3d_n(2,n)+hnode(1,n)/2.0_WP
			end if
		end do ! --> do n=1,myDim_nod2D+eDim_nod2D
		
		!_______________________________________________________________________
		! restart element layer thinkness (helem) and The increment of total 
		! fluid depth on elements (dhe)
		do elem=1,myDim_elem2D
			elnodes=elem2D_nodes(:,elem)
			!___________________________________________________________________
			if (any(hnode(2:lzstar_lev,elnodes(1))/=(zbar(2:lzstar_lev)-zbar(3:lzstar_lev+1)))     .or. &
				any(hnode(2:lzstar_lev,elnodes(2))/=(zbar(2:lzstar_lev)-zbar(3:lzstar_lev+1))) .or. &
				any(hnode(2:lzstar_lev,elnodes(3))/=(zbar(2:lzstar_lev)-zbar(3:lzstar_lev+1)))      &
				) then
				! --> case local zstar 
				! try to limitate over how much layers i realy need to distribute
				! the change in ssh, so that the next loops run only over the 
				! nesseccary levels and not over all lzstar_lev levels
				nz    = max(1 ,maxval(pack(idx,hnode(1:lzstar_lev,elnodes(1))/=(zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1)))))
				nz    = max(nz,maxval(pack(idx,hnode(1:lzstar_lev,elnodes(2))/=(zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1)))))
				nz    = max(nz,maxval(pack(idx,hnode(1:lzstar_lev,elnodes(3))/=(zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1)))))
				nzmax = min(nz,nlevels_nod2D_min(n)-2)
				do nz=1,nzmax
					helem(nz,elem)=sum(hnode(nz,elnodes))/3.0_WP
				end do
			else
				! --> case normal zlevel
				helem(1,elem)=sum(hnode(1,elnodes))/3.0_WP
			end if
			
			!___________________________________________________________________
			dhe(elem)=sum(hbar(elnodes)-hbar_old(elnodes))/3.0_WP
			
		end do ! --> do elem=1,myDim_elem2D
		
		!_______________________________________________________________________
		deallocate(idx)
		
	!___________________________________________________________________________
	! >->->->->->->->->->->->->->->       z-star       <-<-<-<-<-<-<-<-<-<-<-<-<
	!___________________________________________________________________________
	elseif (trim(which_ale)=='zstar' ) then
		! restart depthlevels (zbar_3d_n) and mitdpethlevels (Z_3d_n)
		! dont forget also at restart zbar_3d_n and Z_3d_n are first initialised 
		! and filled up in ale_init there bottom depth zbar_3d_n(nlevels_nod2d) 
		! ist set according if there are partial cells or not 
		do n=1, myDim_nod2D+eDim_nod2D
			nzmax               =nlevels_nod2D(n)-1
! 			nzmax               =nlevels_nod2D(n)-2
			do nz=nzmax,1,-1
				zbar_3d_n(nz,n) =zbar_3d_n(nz+1,n) + hnode(nz,n)
				Z_3d_n(nz,n)    =zbar_3d_n(nz+1,n) + hnode(nz,n)/2.0_WP
			end do
		end do
		
		!_______________________________________________________________________
		! restart element layer thinkness (helem) and the increment of total 
		! fluid depth on elements (dhe)
		do elem=1, myDim_elem2D
			elnodes=elem2D_nodes(:, elem)
			!___________________________________________________________________
			do nz=1,nlevels(elem)-2
				helem(nz,elem)=sum(hnode(nz,elnodes))/3.0_WP
			end do
			
			!___________________________________________________________________
			dhe(elem)=sum(hbar(elnodes)-hbar_old(elnodes))/3.0_WP
			
		end do
	endif
end subroutine restart_thickness_ale
!
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
subroutine stiff_mat_ale
	use o_PARAM
	use o_MESH
	use g_PARSUP
	use o_ARRAYS, only:zbar_e_bot,bottom_elem_thickness
	use g_CONFIG
	implicit none
	
	integer                             :: n, n1, n2, i, j, row, ed, fileID
	integer                             :: elnodes(3), el(2)
	integer                             :: npos(3), offset, nini, nend
	real(kind=WP)                       :: dmean, ff, factor, a 
	real(kind=WP)                       :: fx(3), fy(3), ax, ay
	integer, allocatable                :: n_num(:), n_pos(:,:), pnza(:), rpnza(:)
	integer, allocatable                :: mapping(:)
	logical                             :: flag
	character*10                        :: mype_string,npes_string
	character*1000                      :: dist_mesh_dir, file_name
	real(kind=WP)                       :: t0, t1
	integer                             :: ierror              ! MPI, return error code
	
	t0=MPI_Wtime()
	!___________________________________________________________________________
	! a)
	ssh_stiff%dim=nod2D
	allocate(ssh_stiff%rowptr(myDim_nod2D+1), ssh_stiff%rowptr_loc(myDim_nod2D+1))
	ssh_stiff%rowptr(1)=1               ! This has to be updated as
										! contiguous numbering is required
	
	allocate(n_num(myDim_nod2D+eDim_nod2D),n_pos(12,myDim_nod2D))
	n_pos=0
	
	!___________________________________________________________________________
	! b) Neighbourhood information
	dO n=1,myDim_nod2d
		! number of neigbouring nodes to node n
		n_num(n)=1
		! position indices of neigbouring nodes
		n_pos(1,n)=n
	end do   
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
	
	!___________________________________________________________________________
	! fill up reduced row vector: indice entry where sparse entrys switch to next 
	! row
	! n_num contains the number of neighbors
	! n_pos -- their indices 
	do n=1,myDim_nod2D
		ssh_stiff%rowptr(n+1) = ssh_stiff%rowptr(n)+n_num(n)
	end do
	
	!___________________________________________________________________________
	! c)
	! how many nonzero entries sparse matrix has
	ssh_stiff%nza = ssh_stiff%rowptr(myDim_nod2D+1)-1								 
	! allocate column ancolindd value array of sparse matrix, have length of nonzero 
	! entrie of sparse matrix 
	allocate(ssh_stiff%colind(ssh_stiff%nza), ssh_stiff%colind_loc(ssh_stiff%nza))
	allocate(ssh_stiff%values(ssh_stiff%nza))
	ssh_stiff%values=0.0_WP  
	
	!___________________________________________________________________________
	! d) 
	! fill up sparse matrix column index 
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

	!!! Thus far everything is in local numbering.	!!!
	!!! We will update it later when the values are !!!
	!!! filled in 									!!!
	
	!___________________________________________________________________________
	! e) fill in 
	!M/dt-alpha*theta*g*dt*\nabla(H\nabla\eta))
	n_num=0 ! is here reinitialised to become auxilary array to switch from local to global node indices
	
	! 1st do secod term of lhs od equation (12) of "FESOM2 from finite element to finite volumes"
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
! 			fy(1:3) = (zbar_e_bot(el(i)))* &
			fy(1:3) = (zbar(nlevels(el(i))-1)-bottom_elem_thickness(el(i)))* &
					  (gradient_sca(1:3,el(i)) * edge_cross_dxdy(2*i  ,ed)   &
					 - gradient_sca(4:6,el(i)) * edge_cross_dxdy(2*i-1,ed))
			
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
	
	! 2nd do first term of lhs od equation (12) of "FESOM2 from finite element to finite volumes"
	! Mass matrix part
	do row=1, myDim_nod2D
		offset = ssh_stiff%rowptr(row)
		SSH_stiff%values(offset) = SSH_stiff%values(offset)+ area(1,row)/dt
	end do
	deallocate(n_pos,n_num)
	
	!___________________________________________________________________________
	! f)
	! =================================
	! Global contiguous numbers:
	! =================================
	! Now we need to exchange between PE to know their 
	! numbers of non-zero entries (nza):
	! ================================= 
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
	
	! =================================
	! replace local nza with a global one
	! =================================
	ssh_stiff%nza=sum(rpnza(1:npes))
	deallocate(rpnza, pnza)
	! ==================================
	! colindices are now local to PE. We need to make them local
	! contiguous
	! ==================================
	! (i) global natural: 
	do n=1,ssh_stiff%rowptr(myDim_nod2D+1)-ssh_stiff%rowptr(1)  
		! myList_nod2D ... contains global node index of every meshpoit that belongs 
		! to a CPU
		! ssh_stiff%colind(n) ... contains local node index (1:myDim_nod2d)
		! myList_nod2D(ssh_stiff%colind(n))  ... converts local index to global index
		ssh_stiff%colind(n)=myList_nod2D(ssh_stiff%colind(n))    
	end do
	
	allocate(mapping(nod2d))
	! 0 proc reads the data in chunks and distributes it between other procs
	write(npes_string,"(I10)") npes
	dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'
	file_name=trim(dist_mesh_dir)//'rpart.out'
	fileID=10
	if (mype==0) then
		write(*,*) 'in stiff_mat_ale, reading ', trim(file_name)
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
		write(*,*) 'Building SSH operator took ', t1-t0, ' seconds'
		write(*,*) '========================='
	endif
end subroutine stiff_mat_ale
!
!
!_______________________________________________________________________________
! Update ssh stiffness matrix for a new elevation
subroutine stiff_mat_ale_update
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

end subroutine stiff_mat_ale_update
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
	real(kind=WP) :: dumc1_1, dumc1_2, dumc2_1, dumc2_2 !!PS
	
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
	! at this point: ssh_rhs     = -alpha * nabla*int(u^n + deltau dz)
	!                ssh_rhs_old = - nabla*int(u^n dz) - water_flux*area
	!                
	! (eta_(n+1)-eta_n)/dt = ssh_rhs - alpha*water_flux*area + (1-alpha)*ssh_rhs_old
	!                      = ssh_rhs
	!                      
	! shown in eq (11) rhs of "FESOM2: from finite elements to finte volumes, S. Danilov..." eq. (11) rhs
	if ( .not. trim(which_ALE)=='linfs') then
		do n=1,myDim_nod2D
			ssh_rhs(n)=ssh_rhs(n)-alpha*water_flux(n)*area(1,n)+(1.0_WP-alpha)*ssh_rhs_old(n)
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
	hbar(1:myDim_nod2D)=hbar_old(1:myDim_nod2D)+ssh_rhs_old(1:myDim_nod2D)*dt/area(1,1:myDim_nod2D)
	
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
	use g_config,only: dt, which_ALE,min_hnode,lzstar_lev
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_comm_auto
	use io_RESTART !!PS
	use i_arrays !!PS
	use g_forcing_arrays !!PS
	implicit none
	
	integer       :: el(2), enodes(2), n, nz, ed
	real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2, dd, dd1, dddt, cflmax
	integer 	  :: nzmax 
	
	!_______________________________
	! --> zlevel with local zstar
	real(kind=WP) :: dhbar_total, dhbar_rest, distrib_dhbar_int  !PS
	real(kind=WP), dimension(:), allocatable :: max_dhbar2distr,cumsum_maxdhbar,distrib_dhbar
	integer      , dimension(:), allocatable :: idx
	
	!___________________________________________________________________________
	! Contributions from levels in divergence
	Wvel=0.0_WP
	if (Fer_GM) then
		fer_Wvel=0.0_WP
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
	end do ! --> do ed=1, myDim_edge2D
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
	! |		Full free surface cases:
	! |		------------------------
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
		allocate(max_dhbar2distr(lzstar_lev),distrib_dhbar(lzstar_lev),idx(lzstar_lev),cumsum_maxdhbar(lzstar_lev))
		idx = (/(nz,nz=1,lzstar_lev,1)/)
		!!PS allocate(max_dhbar2distr(nl-1),distrib_dhbar(nl-1),idx(nl-1),cumsum_maxdhbar(nl-1))
		!!PS idx = (/(nz,nz=1,nl-1,1)/)
		
		do n=1, myDim_nod2D
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
			if (dhbar_total<0.0_WP .and. hnode(1,n)+dhbar_total<=(zbar(1)-zbar(2))*min_hnode ) then 
				! --> do local zstar case 
				!_______________________________________________________________
				! max_dhbar2distr ... how much negative ssh change can be maximal 
				! distributed per layer (must be negativ, if positive or ==0 
				! layer reached already minimum layerthickness)
				max_dhbar2distr = 0.0_WP
				max_dhbar2distr = (zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1))*min_hnode - hnode(1:lzstar_lev,n);
				where (max_dhbar2distr>=0.0_WP) max_dhbar2distr=0.0_WP
				
				!_______________________________________________________________
				! if vertical CFL criteria at a certain node is at its limit 
				! don't take away further layer thickness --> take it than better 
				! from a deeper layer
				where (CFL_z(1:lzstar_lev,n)>=0.95_WP) max_dhbar2distr=0.0_WP
				
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
				
				!_______________________________________________________________
				! calc array for distribution of ssh change over layers
				distrib_dhbar = 0.0_WP
				dhbar_rest    = dhbar_total
				
				! nlevels_nod2D_min(n)-1 ...would be hnode of partial bottom 
				! cell but this one is not allowed to change so go until 
				! nlevels_nod2D_min(n)-2
				nzmax = min(nz,nlevels_nod2D_min(n)-2)
				do nz=1,nzmax
					distrib_dhbar(nz) = max(dhbar_rest,max_dhbar2distr(nz))    
					dhbar_rest        = dhbar_rest - distrib_dhbar(nz)
					dhbar_rest        = min(0.0_WP,dhbar_rest)
				end do
				
				!_______________________________________________________________
				if ( abs(sum(distrib_dhbar)-dhbar_total)>1.0e-10 ) then
					write(*,*) " --> problem <-- with conservation of dhbar distribution over depth"
					write(*,*) "				 there are not enough layers to distribute all "
					write(*,*) "				 all change in ssh "
					write(*,*) "				  > mype        =",mype
					write(*,*) "				  > node        =",n
					write(*,*) "				  > mstep       =",mstep
					write(*,*) "				  > dhbar_total =",dhbar_total
					write(*,*) "				  > dhbar_rest  =",dhbar_rest
				end if 
				
				!_______________________________________________________________
				distrib_dhbar_int = 0.0_WP
				do nz=nzmax,1,-1
					!___________________________________________________________
					! --> integrate ssh distribution from down to up
					distrib_dhbar_int = distrib_dhbar_int + distrib_dhbar(nz)
					
					!___________________________________________________________
					! --> distribute change in ssh over layers in hnode and Wvel
					Wvel(nz,n)        = Wvel(nz,n) - distrib_dhbar_int/dt
					hnode_new(nz,n)   = hnode(nz,n)+ distrib_dhbar(nz)
				end do
				
			!___________________________________________________________________
			! in case local zstar was applied must allow the mesh in case of 
			! positive ssh change to return to the normal zlevel case, that means
			! to first "refill" the subsurface layerthickness and with the rest 
			! than the surface layerthickness
			elseif (dhbar_total>0.0_WP .and. & 
					any(hnode(2:lzstar_lev,n)/=(zbar(2:lzstar_lev)-zbar(3:lzstar_lev+1))) &
					) then
				! --> do return to zlevel
				!_______________________________________________________________
				! max_dhbar2distr ... how much positive ssh change must be 
				! distributed in the subsurface layers to be able to return to 
				! the init layerthickness
				max_dhbar2distr   = 0.0_WP
				max_dhbar2distr   = (zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1)) - hnode(1:lzstar_lev,n);
				! there is no limitation in the surface layer how much positive 
				! ssh change can be put there (1000.0_WP is just an arbitrary 
				! high value that should no be reached by dhbar_total)
				max_dhbar2distr(1)= 1000.0_WP
				
				!_______________________________________________________________
				! try to limitate over how much layers i realy need to distribute
				! the change in ssh, so that the next loops run only over the 
				! nesseccary levels and not over all lzstar_lev levels
				nz = maxval(pack(idx,hnode(1:lzstar_lev,n)/=(zbar(1:lzstar_lev)-zbar(2:lzstar_lev+1))))
				
				! nlevels_nod2D_min(n)-1 ...would be hnode of partial bottom 
				! cell but this one is not allowed to change so go until 
				! nlevels_nod2D_min(n)-2
				nzmax = min(nz,nlevels_nod2D_min(n)-2)
				
				!_______________________________________________________________
				! calc array for distribution of ssh change over layers
				dhbar_rest        = dhbar_total
				distrib_dhbar     = 0.0_WP
				distrib_dhbar_int = 0.0_WP
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
					Wvel(nz,n)        = Wvel(nz,n) - distrib_dhbar_int/dt
					hnode_new(nz,n)   = hnode(nz,n)+ distrib_dhbar(nz)
					
				end do
				
			!___________________________________________________________________
			else
				! --> do normal zlevel case
				! only distribute change in ssh for Wvel and hnode_new into the 
				! surface layer
				Wvel(1,n)      = Wvel(1,n) -dhbar_total/dt
				hnode_new(1,n) = hnode(1,n)+dhbar_total
				
			end if ! --> if (dhbar_total<0 .and. hnode(1,n)+dhbar_total<=... ) then 
			
			!___________________________________________________________________
			! Add surface fresh water flux as upper boundary condition for continutity
			Wvel(1,n) = Wvel(1,n)-water_flux(n)
			
		end do ! --> do n=1, myDim_nod2D
		
		!_______________________________________________________________________
		deallocate(max_dhbar2distr,distrib_dhbar,idx,cumsum_maxdhbar)
		
	!___________________________________________________________________________
	elseif (trim(which_ALE)=='zstar') then
		! distribute total change in ssh (hbar(n)-hbar_old(n)) over all layers 
		do n=1, myDim_nod2D
			!___________________________________________________________________
			! --> be careful Sergey suggest in his paper to use the unperturbed
			!     ocean levels NOT the actual one !!! but spoke with Sergey its not 
			!     so important which to use as long as it is consistent and 
			!     volume is conserved
			dd1=zbar_3d_n(nlevels_nod2D_min(n)-1,n)
			
			! This is the depth the stretching is applied (area(nz,n)=area(1,n))
			dd=zbar_3d_n(1,n)-dd1    
			
			! how much of (hbar(n)-hbar_old(n)) is distributed into each layer
			! 1/H*dhbar
			dd=(hbar(n)-hbar_old(n))/dd
			
			!___________________________________________________________________
			! 1/H*dhbar/dt
			dddt=dd/dt
			
			!___________________________________________________________________
			do nz=1,nlevels_nod2D_min(n)-2
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
			Wvel(1,n)=Wvel(1,n)-water_flux(n) 
			
		end do ! --> do n=1, myDim_nod2D
		! The implementation here is a bit strange, but this is to avoid 
		! unnecessary multiplications and divisions by area. We use the fact 
		! that we apply stretching only over the part of the column
		! where area(nz,n)=area(1,n)
		
	endif ! --> if(trim(which_ALE)=='....') then
	
	if (any(hnode_new<0.0_WP)) then
		write(*,*) ' --> fatal problem <--: layerthickness of a layer became smaller zero'
	endif
	
	!___________________________________________________________________________
	call exchange_nod(Wvel) 
	call exchange_nod(hnode_new)   ! Or extend cycles above  
	if (Fer_GM) call exchange_nod(fer_Wvel)
	
	!___________________________________________________________________________
	! calc vertical CFL criteria for debugging purpose and vertical Wvel splitting
	do n=1, myDim_nod2D+eDim_nod2D
		do nz=1,nlevels_nod2D(n)-1
			! CFL from positiv bottom face prism veloctiy
			c1=max(0.0_WP,Wvel(nz+1,n))*dt/hnode_new(nz,n)
			! CFL from negative top face prism veloctiy
			c2=abs(min(0.0_WP,Wvel(nz  ,n))*dt/hnode_new(nz,n))
			! maximum CFL
			CFL_z(nz,n)=max(c1, c2)
		end do
	end do
    cflmax=maxval(CFL_z(:, 1:myDim_nod2D)) !local CFL maximum is different on each mype
    if (cflmax>1.0) then
       do n=1, myDim_nod2D
          do nz=1,nlevels_nod2D(n)-1
             if (abs(CFL_z(nz,n)-cflmax) < 1.e-12) then
                write(*,*) '***********************************************************'
                write(*,*) 'max. CFL_z = ', cflmax, ' mype = ', mype
	            write(*,*) 'mstep      = ', mstep
                write(*,*) 'glon, glat = ', geo_coord_nod2D(:,n)/rad
                write(*,*) '2D node    = ', myList_nod2D(n)
                write(*,*) 'nz         = ', nz
                write(*,*) '***********************************************************'
             end if
           end do
       end do
    end if
	
	!___________________________________________________________________________
	! Split implicit vertical velocity onto implicit and explicit components
	if (w_split) then
		do n=1, myDim_nod2D+eDim_nod2D
! 			do nz=1,nlevels_nod2D(n)-1
			do nz=1,nlevels_nod2D(n)
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
     rinfo, MPI_COMM_FESOM)
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
	reuse=1     ! For varying coefficients, set reuse=1
	new_values=1 !PS 1 ! and new_values=1, as soon as the coefficients have changed
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
        ssh_stiff%colind-1, ssh_stiff%values, reuse, MPI_COMM_FESOM)
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
	
	!___________________________________________________________________________
	! Here can not exchange zbar_n & Z_n with zbar_3d_n & Z_3d_n because 
	! they run over elements here 
	nzmax =nlevels(elem)
	zbar_n=0.0_WP
	Z_n   =0.0_WP
	! in case of partial cells zbar_n(nzmax) is not any more at zbar(nzmax), 
	! zbar_n(nzmax) is now zbar_e_bot(elem), 
	zbar_n(nzmax)=zbar_e_bot(elem)
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
	use g_config, only: logfile_outfreq
	use o_MESH
	use o_ARRAYS
	use o_PARAM
	use g_PARSUP
	use g_comm_auto
	use io_RESTART !PS
	use i_ARRAYS !PS
	use o_mixing_KPP_mod
	
	IMPLICIT NONE
	real(kind=8)      :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
	integer           :: n
	
	t1=MPI_Wtime()
	!___________________________________________________________________________
	call pressure_bv               !!!!! HeRE change is made. It is linear EoS now.
	call pressure_force
	!___________________________________________________________________________
	! calculate alpha and beta
	! it will be used for KPP, Redi, GM etc. Shall we keep it on in general case?
	call sw_alpha_beta(tr_arr(:,:,1),tr_arr(:,:,2))
	
	! computes the xy gradient of a neutral surface; will be used by Redi, GM etc.
	call compute_sigma_xy(tr_arr(:,:,1),tr_arr(:,:,2))
	
	! compute both: neutral slope and tapered neutral slope. Can be later combined with compute_sigma_xy
	! will be primarily used for computing Redi diffusivities. etc?
	call compute_neutral_slope
	
	!___________________________________________________________________________
	call status_check
	
	!___________________________________________________________________________
	if (trim(mix_scheme)=='KPP') then
		call oce_mixing_KPP(Av, Kv_double)
		Kv=Kv_double(:,:,1)
	else if(trim(mix_scheme)=='PP') then
		call oce_mixing_PP
	else
		stop "!not existing mixing scheme!"
		call par_ex
	end if  
	
	!___________________________________________________________________________
	if(mom_adv/=3) then
		call compute_vel_rhs
	else
		call compute_vel_rhs_vinv
	end if
	
	!___________________________________________________________________________
	call viscosity_filter(2)
	
	!___________________________________________________________________________
	if(i_vert_visc) call impl_vert_visc_ale
	t2=MPI_Wtime()
	
	!___________________________________________________________________________
	! >->->->->->->->->->->->->     ALE-part starts     <-<-<-<-<-<-<-<-<-<-<-<-
	!___________________________________________________________________________
	! Update stiffness matrix by dhe=hbar(n+1/2)-hbar(n-1/2) on elements
	call stiff_mat_ale_update 
	
	! ssh_rhs=-alpha*\nabla\int(U_n+U_rhs)dz-(1-alpha)*...
	! see "FESOM2: from finite elements to finte volumes, S. Danilov..." eq. (12) rhs
	call compute_ssh_rhs_ale
	
	! Take updated ssh matrix and solve --> new ssh!
	call solve_ssh_ale
	t3=MPI_Wtime() 
	
	! estimate new horizontal velocity u^(n+1)
	! u^(n+1) = u* + [-g * tau * theta * grad(eta^(n+1)-eta^(n)) ]
	call update_vel
	! --> eta_(n) --> eta_(n+1) = eta_(n) + deta = eta_(n) + (eta_(n+1) + eta_(n))
	t4=MPI_Wtime() 
	
	! Update to hbar(n+3/2) and compute dhe to be used on the next step
	call compute_hbar_ale
	t5=MPI_Wtime() 
	
	!___________________________________________________________________________
	! Current dynamic elevation alpha*hbar(n+1/2)+(1-alpha)*hbar(n-1/2)
	! equation (7) Danlov et.al "the finite volume sea ice ocean model FESOM2
	! ...if we do it here we don't need to write hbar_old into a restart file...
	eta_n=alpha*hbar+(1.0_WP-alpha)*hbar_old
	! --> eta_(n)
	! call zero_dynamics !DS, zeros several dynamical variables; to be used for testing new implementations!
        
        if (Fer_GM .or. Redi) then
           call init_Redi_GM
        end if
	
	!___________________________________________________________________________
	! Implementation of Gent & McWiliams parameterization after R. Ferrari et al., 2010
	! does not belong directly to ALE formalism
	if (Fer_GM) then
		call fer_solve_Gamma
		call fer_gamma2vel
	end if
	t6=MPI_Wtime() 
	
	!___________________________________________________________________________
	!The main step of ALE procedure --> this is were the magic happens --> here 
	! is decided how change in hbar is distributed over the vertical layers
	call vert_vel_ale 
	t7=MPI_Wtime() 
	
	!___________________________________________________________________________
	! solve tracer equation
	call solve_tracers_ale
	t8=MPI_Wtime() 
	
	!___________________________________________________________________________
	! Update hnode=hnode_new, helem
	call update_thickness_ale  
	t9=MPI_Wtime() 
	
	!___________________________________________________________________________
	! write out global fields for debugging
	call write_step_info(n,logfile_outfreq)
	
	! check model for blowup --> ! write_step_info and check_blowup require 
	! togeather around 2.5% of model runtime
	call check_blowup(n)
	
	!___________________________________________________________________________
	! write out execution times for ocean step parts
	t10=MPI_Wtime() 
	if(mod(n,logfile_outfreq)==0 .and. mype==0) then  
		write(*,*) '___ALE OCEAN STEP EXECUTION TIMES______________________'
		write(*,"(A, ES10.3)") ' 	Oce. Dynamics    :', t2-t1
		write(*,"(A, ES10.3)") ' 	Oce. Update Vel. :', t4-t3
		write(*,"(A, ES10.3)") ' 	Oce. Fer-GM.     :', t6-t5
		write(*,*) '	_______________________________'
		write(*,"(A, ES10.3)") ' 	ALE-Solve SSH    :', t3-t2
		write(*,"(A, ES10.3)") ' 	ALE-Calc. hbar   :', t5-t4
		write(*,"(A, ES10.3)") ' 	ALE-Update+W     :', t7-t6
		write(*,"(A, ES10.3)") ' 	ALE-Solve Tracer :', t8-t7
		write(*,"(A, ES10.3)") ' 	ALE-Update hnode :', t9-t8
		write(*,*) '	_______________________________'
		write(*,"(A, ES10.3)") ' 	check for blowup :', t10-t9
		write(*,*) '	_______________________________'
		write(*,"(A, ES10.3)") ' 	Oce. TOTAL       :', t10-t1
		write(*,*)
		write(*,*)
	end if
	
end subroutine oce_timestep_ale
