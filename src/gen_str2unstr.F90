	module g_str2unstr
	use o_mesh
	implicit none
	save
	real(kind=8)                  :: CplLonMin=-180.
	real(kind=8)                  :: CplLonMax=180.
	real(kind=8)                  :: CplLatMin=-90.
	real(kind=8)                  :: CplLatMax=90.
	real*8, allocatable           :: regtr_arr(:,:,:,:)
	real*8, allocatable           :: regssh(:,:)
	real*8, allocatable           :: regu(:,:,:),regv(:,:,:)
	type(sparse_matrix)           :: oce_2_atm,oce_2_atm_u,oce_2_atm_v
	real*8,allocatable,dimension(:):: aveta_n
	real*8,allocatable,dimension(:,:,:)::avtr_arr
	real*8,allocatable,dimension(:,:)::avUnode
	real*8,allocatable,dimension(:,:)::avVnode
	logical :: setup_sc=.false.,setup_u=.false.,setup_v=.false.
	logical :: setup_done=.false.
	contains
	!=================================================================
	! Determines whether point pt belongs to any triangle in myDim_elem2D
	! point is in lon/lat
	! returns 0 if no triangle found
	!=================================================================
	SUBROUTINE point_in_triangle(el2D,   pt)
	  use g_parsup; use g_comm_auto
	  use o_param
	  use o_mesh
	  use g_rotate_grid

	  implicit none

	  INTEGER, INTENT(OUT)                   :: el2D  
	  REAL(kind=8), DIMENSION(2), INTENT(IN) :: pt
	  
	  integer                                :: elem, elnodes(3), q
	  real(kind=8)                           :: alph, mean_lon, mean_lat, rlon, rlat
	  real(kind=8)                           :: xe(4), ye(4), xt1, xt2, yt1, yt2, x1, x2, y1, y2
	  real(kind=8)                           :: s1, s2, angle
	  
	  el2D=0

	  DO elem=1, myDim_elem2D
	     elnodes(1:3)=elem2D_nodes(:, elem)
	     ! ========
	     ! Very rough criteria to reduce the work
	     ! ========
	     if (rotated_grid) then 
		do q=1,3	
		   rlon=coord_nod2D(1, elnodes(q))
		   rlat=coord_nod2D(2, elnodes(q))
		   call r2g(xe(q), ye(q), rlon, rlat)
		end do
	     else
		write(*,*) 'point_in_triangle, something is wrong'
		call par_ex
		stop	
		xe(1:3)=coord_nod2D(1, elnodes(1:3))
		ye(1:3)=coord_nod2D(2, elnodes(1:3))
	     endif
		    !=====
		    ! Cyclicity:
		    ! Remove 2*pi jumps in x-coordinate
		    ! of nodes in triangle  
		    ! =====
		    xe(2)=xe(2)-xe(1)
		    xe(3)=xe(3)-xe(1)
		    DO q=2,3
		     if(xe(q)> pi) xe(q)=xe(q)-2*pi 
		     if(xe(q)<-pi) xe(q)=xe(q)+2*pi  
		    END DO
		    xe(2)=xe(2)+xe(1)
		    xe(3)=xe(3)+xe(1)
		    ! =====
		    ! Now x(2) and x(3) are measured consistent to x(1) 
		    mean_lon=sum(xe(1:3))/3.0_8
		    mean_lat=sum(ye(1:3))/3.0_8
		    xt1=maxval(abs(xe(1:3)-mean_lon))
		    yt1=maxval(abs(ye(1:3)-mean_lat))
		    ! =====
		    ! Cyclicity
		    ! =====
		    if(xt1> pi) xt1=xt1-2*pi 
		    if(xt1<-pi) xt1=xt1+2*pi  
		    
		    x1=pt(1)-mean_lon
		    if(x1>pi) x1=x1-2*pi 
		    if(x1<-pi) x1=x1+2*pi  
		    
		    if((abs(x1)>xt1).or.(abs(pt(2)-mean_lat)>yt1)) cycle
		    
		    ! ========
		    ! Find if the regular mesh node is within elem
		    ! ========
		    xe(4)=xe(1)
		    ye(4)=ye(1)
		    angle=0.
		    DO q=1,3
		       xt1=xe(q) 
		       xt2=xe(q+1)
		       yt1=ye(q)
		       yt2=ye(q+1)
		       x1=xt1-pt(1)
		       x2=xt2-pt(1)
		       y1=yt1-pt(2)
		       y2=yt2-pt(2)
		       ! =====
		       ! Cyclicity check
		       ! =====
		       if(x1>pi) x1=x1-2*pi 
		       if(x1<-pi) x1=x1+2*pi 
		       if(x2>pi) x2=x2-2*pi 
		       if(x2<-pi) x2=x2+2*pi 

		       alph=x1*y2-x2*y1
		       s1=x1*x1+y1*y1
		       s2=x2*x2+y2*y2
		       alph=alph/sqrt(s1*s2)
		       alph=asin(alph)
		       IF ((xt1-xt2)**2+(yt1-yt2)**2 > max(s1,s2)) THEN
		       if (alph>0) alph=pi-alph
		       if (alph<=0) alph=-pi-alph
		       END IF
		       angle=angle+alph
		     END DO
		     IF (abs(angle)>pi) THEN
		     el2D=elem
		     exit
		     END IF
		 ENDDO
	 ! ============
	 ! A triangle does not belong to a single PE, thus some i,j will be 
	 ! updated on several PEs. It is assumed that elem belongs to a PE 
	 ! that own the first node.
	 ! ============
	  if(el2D>0) then
	    if (elem2D_nodes(1, el2D) > myDim_nod2D)  then
		    el2D=0
	    end if
	  end if 
	END SUBROUTINE point_in_triangle

! =======================================================================
	subroutine build_oce_2_atm(atm_lon,atm_nx,atm_lat,atm_ny,mytype)
	!mytype =  0 -- scalar
	!          1 -- U
	!         -1 -- V
	use g_parsup; use g_comm_auto
	use o_param
	use o_mesh
	use g_rotate_grid

		implicit none
		integer                   :: n, nij, i, j, pos, el, is, ie, case_stat(2)
		real(kind=8)              :: x, y, xr, yr, xx, yy, p(2), rlon, rlat
		real(kind=8)              :: dx, dy
		real(kind=8), allocatable :: xcoord(:), ycoord(:)
		real(kind=8), allocatable :: real_loc(:, :)
		integer, intent(in) :: atm_nx,atm_ny,mytype
		real(kind=8), intent(in) :: atm_lon(atm_nx),atm_lat(atm_ny)
		integer, allocatable      :: do_lin_int(:), my_count(:), gl_count(:), my_mapp2(:)
		real*8 t1,t2,t3
		if (mytype>0) then
!U
		   if (setup_u) return
		else if(mytype<0) then
!V
		   if (setup_v) return
		else
!scalar
		   if (setup_sc) return
		endif
		if (mype==0) write(*,*) 'build_oce_2_atm: START'
	! global arrays of the coordinates are required
		allocate(xcoord(nod2d))
		allocate(ycoord(nod2d))
		allocate(real_loc(2, myDim_nod2d+eDim_nod2d))
		allocate(my_mapp2(nod2D))
		allocate(my_count    (atm_nx*atm_ny))
		allocate(gl_count    (atm_nx*atm_ny))
		allocate(do_lin_int  (atm_nx*atm_ny))
		
	!compute radius dependent on exchange grid	
		dx=360./2./real(atm_nx)
		dy=180./2./real(atm_ny)
		
	!case 1 (normal case)     : the interpolation weight equals my_count/gl_count	
	!case 2 (exceptional case): linear interpolation within the triangle el=do_lin_int(nij)		
		
		do_lin_int  = 0
		my_count    = 0
		gl_count    = 0
		my_mapp2    = 0
		xcoord      = 0.
		ycoord      = 0.
		!oce_2_atm_mask=0
		case_stat   = 0
			
		do n=1, myDim_nod2d+eDim_nod2d
		   x=coord_nod2D(1, n)
		   y=coord_nod2D(2, n)
		   call r2g(real_loc(1, n), real_loc(2, n), x, y)
		end do
		do n=1, myDim_nod2d+eDim_nod2d
		   my_mapp2(myList_nod2D(n))=n
		end do


		call broadcast_nod(real_loc(1, :), xcoord)!broadcast2D(real_loc(1, :), xcoord)
		call broadcast_nod(real_loc(2, :), ycoord)
		t1=MPI_Wtime()
	! start building the atm_2_oce operator
	!1. Dimension=[atm_nx*atm_ny x myDim_nod2d]
	if (mype==0) write(*,*) 'MARK1'
		if (mytype>0) then
		 setup_u=.true.
		 oce_2_atm_U%dim=atm_nx*atm_ny
		 allocate(oce_2_atm_U%rowptr(oce_2_atm_U%dim+1))
		 oce_2_atm_U%rowptr(:)=0
		 oce_2_atm_U%rowptr(1)=1
		elseif (mytype<0) then
		 setup_v=.true.
		 oce_2_atm_V%dim=atm_nx*atm_ny
		 allocate(oce_2_atm_V%rowptr(oce_2_atm_V%dim+1))
		 oce_2_atm_V%rowptr(:)=0
		 oce_2_atm_V%rowptr(1)=1
		else
		 setup_sc=.true.
		 oce_2_atm%dim=atm_nx*atm_ny
		 allocate(oce_2_atm%rowptr(oce_2_atm%dim+1))
		 oce_2_atm%rowptr(:)=0
		 oce_2_atm%rowptr(1)=1
		endif
	!2. rowptr
	if (mype==0) write(*,*) 'MARK2'
		do i=1, atm_nx
		   do j=1, atm_ny
		      nij=(i-1)*atm_ny+j
		      pos=nij+1
		      xr=atm_lon(i)/rad
		      if (xr < -180) xr=xr+360
		      if (xr >  180) xr=xr-360	      
		      yr=atm_lat(j)/rad
		      do n=1, nod2d
			 x=xcoord(n)/rad
			 y=ycoord(n)/rad
			 if (x < -180) x=x+360
			 if (x >  180) x=x-360		 
			 xx=abs(xr-x)
			 xx=min(xx, abs(360-xx))
			 yy=abs(yr-y)
			 if (xx > dx  .OR. yy > dy) CYCLE
			 gl_count(nij)=gl_count(nij)+1
			 if (my_mapp2(n) <= myDim_nod2d .AND. my_mapp2(n)>0) then
			    my_count(nij)=my_count(nij)+1
			 end if
		      end do
	! c > 0 means that the points around [i,j] were found. Otherwise we check 
	! if there is a "rather big" triangle that [i,j] is inside it
		      if (gl_count(nij) <= 3) then
			 my_count(nij)=0
			 gl_count(nij)=0
		      end if
		      if (my_count(nij) > 0) then
			 if (mytype>0) then
			    oce_2_atm_u%rowptr(pos)=oce_2_atm_u%rowptr(pos-1)+my_count(nij)
			 elseif (mytype<0) then
			    oce_2_atm_v%rowptr(pos)=oce_2_atm_v%rowptr(pos-1)+my_count(nij)
			 else
			    oce_2_atm%rowptr(pos)=oce_2_atm%rowptr(pos-1)+my_count(nij)
			 endif
		      else if (my_count(nij)==0 .AND. gl_count(nij)==0) then
			 p(1)=xr*rad
			 p(2)=yr*rad      
			 call point_in_triangle(el,   p)
			 if (mytype>0) then
			   if ((el > 0)) then
	! we will do a linear interpolation to this point
			    oce_2_atm_u%rowptr(pos)=oce_2_atm_u%rowptr(pos-1)+3
			    do_lin_int(nij)      =el		     
			   else
			    oce_2_atm_u%rowptr(pos)=oce_2_atm_u%rowptr(pos-1)+1
			   end if
			 elseif (mytype<0) then
			   if ((el > 0)) then
	! we will do a linear interpolation to this point
			    oce_2_atm_v%rowptr(pos)=oce_2_atm_v%rowptr(pos-1)+3
			    do_lin_int(nij)      =el		     
			   else
			    oce_2_atm_v%rowptr(pos)=oce_2_atm_v%rowptr(pos-1)+1
			   end if
			 else
			   if ((el > 0)) then
	! we will do a linear interpolation to this point
			    oce_2_atm%rowptr(pos)=oce_2_atm%rowptr(pos-1)+3
			    do_lin_int(nij)      =el		     
			   else
			    oce_2_atm%rowptr(pos)=oce_2_atm%rowptr(pos-1)+1
			   end if
			 endif
		       else
			 if (mytype>0) then
			   oce_2_atm_u%rowptr(pos)=oce_2_atm_u%rowptr(pos-1)+1   
			 elseif (mytype<0) then
			   oce_2_atm_v%rowptr(pos)=oce_2_atm_v%rowptr(pos-1)+1   
			 else
			   oce_2_atm%rowptr(pos)=oce_2_atm%rowptr(pos-1)+1   
			 endif
		      end if	      
		    end do
		 end do

		t2=MPI_Wtime()
	if (mype==0) write(*,*) 'MARK3'
	!3. colind & values
		if (mytype>0) then
		n=oce_2_atm_u%rowptr(oce_2_atm_u%dim+1)-1
		allocate(oce_2_atm_u%values(n))
		allocate(oce_2_atm_u%colind(n))
		oce_2_atm_u%colind=1
		oce_2_atm_u%values=0.
		oce_2_atm_u%nza   =n
		do i=1, atm_nx
		   do j=1, atm_ny
		      nij=(i-1)*atm_ny+j
		      if (do_lin_int(nij)==0 .AND. my_count(nij)==0) CYCLE !do nothing
		      pos=oce_2_atm_u%rowptr(nij)	      
		      xr=atm_lon(i)/rad
		      yr=atm_lat(j)/rad	      
		      if (xr < -180) xr=xr+360
		      if (xr >  180) xr=xr-360
	!fill colind and values, CASE 1	      
		      if (my_count(nij) > 0) then
			 case_stat(1)=case_stat(1)+1
	!	         write(*,*) 'oce_2_atm_u: filling the values, CASE 1'	      	      
	!local loop to fill the operator
			 do n=1, myDim_nod2d
			      x=real_loc(1, n)/rad
			      y=real_loc(2, n)/rad
			      if (x < -180) x=x+360
			      if (x >  180) x=x-360		 
			      xx=abs(xr-x)
			      xx=min(xx, abs(360-xx))
			      yy=abs(yr-y)
			      if (xx > dx .OR. yy > dy) CYCLE
			      oce_2_atm_u%colind(pos)=n
			      oce_2_atm_u%values(pos)=1./real(gl_count(nij))
			      pos=pos+1
			 end do		 
		      else if (do_lin_int(nij)>0) then !do linear interpolation to [i, j]
	!fill colind and values, CASE 2
			 case_stat(2)=case_stat(2)+1
	!	         write(*,*) 'oce_2_atm_u: filling the values, CASE 2'
			 el=do_lin_int(nij)
			 call g2r(xr*rad, yr*rad, rlon, rlat)
			 p(1)=rlon/rad
			 p(2)=rlat/rad
			 oce_2_atm_u%colind(pos:pos+2)=elem2D_nodes(:, el)
			 call locbafu_2D(oce_2_atm_u%values(pos:pos+2), el, p)
		      end if
		end do
		end do
		elseif (mytype<0) then
		n=oce_2_atm_v%rowptr(oce_2_atm_v%dim+1)-1
		allocate(oce_2_atm_v%values(n))
		allocate(oce_2_atm_v%colind(n))
		oce_2_atm_v%colind=1
		oce_2_atm_v%values=0.
		oce_2_atm_v%nza   =n
		do i=1, atm_nx
		   do j=1, atm_ny
		      nij=(i-1)*atm_ny+j
		      if (do_lin_int(nij)==0 .AND. my_count(nij)==0) CYCLE !do nothing
		      pos=oce_2_atm_v%rowptr(nij)	      
		      xr=atm_lon(i)/rad
		      yr=atm_lat(j)/rad	      
		      if (xr < -180) xr=xr+360
		      if (xr >  180) xr=xr-360
	!fill colind and values, CASE 1	      
		      if (my_count(nij) > 0) then
			 case_stat(1)=case_stat(1)+1
	!	         write(*,*) 'oce_2_atm_v: filling the values, CASE 1'	      	      
	!local loop to fill the operator
			 do n=1, myDim_nod2d
			      x=real_loc(1, n)/rad
			      y=real_loc(2, n)/rad
			      if (x < -180) x=x+360
			      if (x >  180) x=x-360		 
			      xx=abs(xr-x)
			      xx=min(xx, abs(360-xx))
			      yy=abs(yr-y)
			      if (xx > dx .OR. yy > dy) CYCLE
			      oce_2_atm_v%colind(pos)=n
			      oce_2_atm_v%values(pos)=1./real(gl_count(nij))
			      pos=pos+1
			 end do		 
		      else if (do_lin_int(nij)>0) then !do linear interpolation to [i, j]
	!fill colind and values, CASE 2
			 case_stat(2)=case_stat(2)+1
	!	         write(*,*) 'oce_2_atm_v: filling the values, CASE 2'
			 el=do_lin_int(nij)
			 call g2r(xr*rad, yr*rad, rlon, rlat)
			 p(1)=rlon/rad
			 p(2)=rlat/rad
			 oce_2_atm_v%colind(pos:pos+2)=elem2D_nodes(:, el)
			 call locbafu_2D(oce_2_atm_v%values(pos:pos+2), el, p)
		      end if
		end do
		end do
		else
		n=oce_2_atm%rowptr(oce_2_atm%dim+1)-1
		allocate(oce_2_atm%values(n))
		allocate(oce_2_atm%colind(n))
		oce_2_atm%colind=1
		oce_2_atm%values=0.
		oce_2_atm%nza   =n
		do i=1, atm_nx
		   do j=1, atm_ny
		      nij=(i-1)*atm_ny+j
		      if (do_lin_int(nij)==0 .AND. my_count(nij)==0) CYCLE !do nothing
		      pos=oce_2_atm%rowptr(nij)	      
		      xr=atm_lon(i)/rad
		      yr=atm_lat(j)/rad	      
		      if (xr < -180) xr=xr+360
		      if (xr >  180) xr=xr-360
	!fill colind and values, CASE 1	      
		      if (my_count(nij) > 0) then
			 case_stat(1)=case_stat(1)+1
	!	         write(*,*) 'oce_2_atm: filling the values, CASE 1'	      	      
	!local loop to fill the operator
			 do n=1, myDim_nod2d
			      x=real_loc(1, n)/rad
			      y=real_loc(2, n)/rad
			      if (x < -180) x=x+360
			      if (x >  180) x=x-360		 
			      xx=abs(xr-x)
			      xx=min(xx, abs(360-xx))
			      yy=abs(yr-y)
			      if (xx > dx .OR. yy > dy) CYCLE
			      oce_2_atm%colind(pos)=n
			      oce_2_atm%values(pos)=1./real(gl_count(nij))
			      pos=pos+1
			 end do		 
		      else if (do_lin_int(nij)>0) then !do linear interpolation to [i, j]
	!fill colind and values, CASE 2
			 case_stat(2)=case_stat(2)+1
	!	         write(*,*) 'oce_2_atm: filling the values, CASE 2'
			 el=do_lin_int(nij)
			 call g2r(xr*rad, yr*rad, rlon, rlat)
			 p(1)=rlon/rad
			 p(2)=rlat/rad
			 oce_2_atm%colind(pos:pos+2)=elem2D_nodes(:, el)
			 call locbafu_2D(oce_2_atm%values(pos:pos+2), el, p)
		      end if
		end do
		end do
		endif
	if (mype==0) write(*,*) 'MARK4'
		t3=MPI_Wtime()
		if (mype==0) write(*,*) 'oce_2_atm statistics on PE ', mype, ' : ', case_stat
		if (mype==0) write(*,*) '====================================================='
		if (mype==0) write(*,*) 'rowptr',t2-t1
		if (mype==0) write(*,*) 'colind&val',t3-t2
		if (mype==0) write(*,*) 'Total',t3-t1
		if (mype==0) write(*,*) '====================================================='
		deallocate(xcoord)
		deallocate(ycoord)
		deallocate(real_loc)
		deallocate(my_mapp2)	
		deallocate(do_lin_int)
		deallocate(my_count)
		deallocate(gl_count)
	    end  subroutine build_oce_2_atm
	! =======================================================================
	!=================================================================
	!coords in deg, elem is LOCAL element; 
	!values of the 3 local basisfunctions at the 
	!position 'coords' are returned
	!=================================================================
	SUBROUTINE locbafu_2D(values, elem, coords)
	  USE o_mesh
	  USE o_param
	  IMPLICIT NONE
	  
	  REAL(kind=8), DIMENSION(3), INTENT(OUT) :: values
	  INTEGER,                    INTENT(IN)  :: elem
	  REAL(kind=8), DIMENSION(2), INTENT(IN)  :: coords
	  
	  INTEGER                      :: i,j
	  INTEGER                      :: node
	  INTEGER, DIMENSION(3)        :: local_nodes
	  REAL(kind=8), DIMENSION(2,3) :: local_coords, local_cart_coords
	  REAL(kind=8), DIMENSION(2,2) :: TRANS, TRANS_inv
	  REAL(kind=8)                 :: DET
	  REAL(kind=8), DIMENSION(2)   :: x, x_cart, stdel_coords
	  REAL(kind=8), DIMENSION(2)   :: vec
	  
	  x(1)=coords(1)*rad
	  x(2)=coords(2)*rad
	  do i=1,3
	     node=elem2D_nodes(i,elem)
	     local_nodes(i)=node
	     local_coords(:,i)=coord_nod2D(:,node)
	  end do

	  DO i=1, 2
	     if (local_coords(1,i+1)-local_coords(1,1) > 180.*rad) then 
		 local_coords(1,i+1)=local_coords(1,i+1)-360.*rad
	     end if
	     if (local_coords(1,i+1)-local_coords(1,1) <-180.*rad) then 
		 local_coords(1,i+1)=local_coords(1,i+1)+360.*rad
	     end if
	  END DO  

	  if (x(1)-local_coords(1,1) > 180.*rad) then 
	      x(1)=x(1)-360.*rad
	  end if
	  if (x(1)-local_coords(1,1) <-180.*rad) then 
	      x(1)=x(1)+360.*rad
	  end if
	  !  cartesian coordinates
	  x_cart(1) = r_earth * COS(x(2)) * x(1)
	  x_cart(2) = r_earth * x(2)  

	  do i=1,3
	     local_cart_coords(1,i) = r_earth * COS(local_coords(2,i)) * local_coords(1,i)
	     local_cart_coords(2,i) = r_earth * local_coords(2,i)
	  end do  
	  DO i=1, 2
	     TRANS(:,i) = local_cart_coords(:,i+1)-local_cart_coords(:,1)     
	  END DO  
	  call matrix_inverse_2x2(TRANS, TRANS_inv, DET)  

	  vec=x_cart-local_cart_coords(:,1)
	  stdel_coords = MATMUL(TRANS_inv, vec)
	  call stdbafu_2D(values,  stdel_coords,1)
	END SUBROUTINE locbafu_2D
!========================================================================================
	SUBROUTINE stdbafu_2D(values,x,m) !(stdbafu,x,m)
	  IMPLICIT NONE
	  INTEGER, INTENT(IN)                       :: m
	  REAL(kind=8), DIMENSION(3,m), INTENT(OUT) :: values
	  REAL(kind=8), DIMENSION(2,m), INTENT(IN)  :: x  
	  INTEGER                                   :: k
	  do k=1,m
	     values(1,k)=1.-x(1,k)-x(2,k)
	     values(2,k)=   x(1,k)
	     values(3,k)=          x(2,k)
	  end do
	END SUBROUTINE stdbafu_2D
!========================================================================================
	subroutine  matrix_inverse_2x2 (A, AINV, DET)
	  !
	  !
	  implicit none
	  !
	  real(kind=8), dimension(2,2), intent(IN)  :: A
	  real(kind=8), dimension(2,2), intent(OUT) :: AINV
	  real(kind=8), intent(OUT)                 :: DET
	  !
	  integer                                   :: i,j
	  !
	  DET  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
	  if ( DET .eq. 0.0 )  then
	     do j=1,2
		write(*,*) (A(i,j),i=1,2)
	     end do
	     stop 'SINGULAR 2X2 MATRIX'
	  else
	     AINV(1,1) =  A(2,2)/DET
	     AINV(1,2) = -A(1,2)/DET
	     AINV(2,1) = -A(2,1)/DET
	     AINV(2,2) =  A(1,1)/DET
	  endif
	end subroutine matrix_inverse_2x2
! =======================================================================
        subroutine do_oce_2_atm(oce, atm, atm_lon,atm_nx,atm_lat,atm_ny,nz,mytype)
        use g_parsup; use g_comm_auto
        use o_param
        use o_mesh

        implicit none
!        save
	
	integer, intent(in) :: atm_nx,atm_ny,nz,mytype
	real(kind=8), intent(in) :: atm_lon(atm_nx),atm_lat(atm_ny)
	real(kind=8), intent( IN    )         :: oce(myDim_nod2d+eDim_nod2d)
	real(kind=8), intent( INOUT )         :: atm(atm_nx, atm_ny)
        integer                               :: i, j, k, is, ie, n
        real(kind=8), allocatable             :: row(:)

	allocate(row(atm_nx))
! Apply operator: A*x
	if (mype==0) write(*,*) 'do_oce_2_atm up'
	if (mytype>0) then
	do i=1, atm_nx
	   do j=1, atm_ny
      	      n=(i-1)*atm_ny+j
	      is=oce_2_atm_U%rowptr(n)
	      ie=oce_2_atm_U%rowptr(n+1)-1
		do k=is,ie,1
		  if (nlevels_nod2D(oce_2_atm_U%colind(k))<nz) cycle 
		enddo
	      atm(i, j)=sum(oce_2_atm_U%values(is:ie)*oce(oce_2_atm_U%colind(is:ie)))
	   end do
	end do
	elseif (mytype<0) then
	do i=1, atm_nx
	   do j=1, atm_ny
      	      n=(i-1)*atm_ny+j
	      is=oce_2_atm_V%rowptr(n)
	      ie=oce_2_atm_V%rowptr(n+1)-1
		do k=is,ie,1
		  if (nlevels_nod2D(oce_2_atm_V%colind(k))<nz) cycle 
		enddo
	      atm(i, j)=sum(oce_2_atm_V%values(is:ie)*oce(oce_2_atm_V%colind(is:ie)))
	   end do
	end do
	else
	do i=1, atm_nx
	   do j=1, atm_ny
      	      n=(i-1)*atm_ny+j
	      is=oce_2_atm%rowptr(n)
	      ie=oce_2_atm%rowptr(n+1)-1
		do k=is,ie,1
		  if (nlevels_nod2D(oce_2_atm%colind(k))<nz) cycle 
		enddo
	      atm(i, j)=sum(oce_2_atm%values(is:ie)*oce(oce_2_atm%colind(is:ie)))
	   end do
	end do
	endif
! Broadcast the result (allreduce & sum)
	do j=1, atm_ny
	   row=atm(:, j)
	   atm(:, j)=0.
	   call MPI_Allreduce (row, atm(:, j), atm_nx, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
	end do
	deallocate(row)	
	if (mype==0) write(*,*) 'do_oce_2_atm down'
        end subroutine do_oce_2_atm
! =======================================================================
      subroutine def_exchange_mesh(atm_lon, atm_nx,atm_lat,atm_ny, &
	LonMin,LonMax,LatMin,LatMax)
!       DESCRIPTION
!       This subroutine creates and allocates the intermediate regular grid
!       dependent on the target atmosphere grid (with lon and lat positions)
!       on which the unstructured grid is interpolated to before given to OASIS
        use g_parsup
        use o_param

        implicit none
      
        integer                                :: i,j 
	integer, intent(in) :: atm_nx,atm_ny
	real(kind=8):: atm_lon(atm_nx),atm_lat(atm_ny)
	real*8 :: Dx,Dy
	real*8,intent(in) :: LonMin,LonMax,LatMin,LatMax
        if (LonMax.gt.180.) then
		write(*,*) 'error in def_exchange_mesh'
		call par_ex(1)
	endif


        Dx  = ((LonMax - LonMin) / real(atm_nx) )
        Dy  = ((LatMax - LatMin) / real(atm_ny) )
     
        atm_lon = 0.0
        atm_lat = 0.0

        do i=1,atm_nx
          atm_lon(i)=min(LonMin+(i-1)*Dx,LonMax)
        end do
	
        do i=1,atm_ny
          atm_lat(i)=min(LatMin+(i-1)*Dy,LatMax)
        end do	

        ! =========
        ! Degrees to radians
        ! =========

        atm_lon=atm_lon*rad
        atm_lat=atm_lat*rad

        Dx  =Dx*rad
        Dy  =Dy*rad

        do i=1,atm_nx
        if((atm_lon(i) > pi).OR.(atm_lon(i) < -pi)) then
           write(*,*) 'The regular mesh longitude should be inside &
			 -180 +180 interval ->',atm_lon(i) 
	   call par_ex
	   stop
        end if
        end do

    end subroutine def_exchange_mesh
! =======================================================================
	subroutine clean_av_arrays
	implicit none
	aveta_n=0d0
	avtr_arr=0d0
	avUnode=0d0
	avVnode=0d0
	end subroutine clean_av_arrays
! =======================================================================
	subroutine unstr2str_up(nx,ny)
	use o_mesh
	use o_param
	USE g_PARSUP
	use g_config
	IMPLICIT NONE
	integer     :: elem_size, node_size,nx,ny

	elem_size=myDim_elem2D+eDim_elem2D
	node_size=myDim_nod2D+eDim_nod2D
	!Arrays for averaging unstr data
	!unstr
	allocate(aveta_n(node_size))
	allocate(avtr_arr(nl-1,node_size,num_tracers))
	allocate(avUnode(nl,node_size),avVnode(nl,node_size))
	!regular
	allocate(regtr_arr(nx,ny,nl-1,num_tracers))
	allocate(regu(nx+1,ny,nl-1),regv(nx,ny+1,nl-1))
	allocate(regssh(nx,ny))
	regtr_arr=-999d0
	regssh=-999d0
	regu=-999d0; regv=-999d0
	end subroutine unstr2str_up
! =======================================================================
	subroutine clean_oce_2_atm(mytype)
	use o_mesh
	implicit none
	integer mytype
	deallocate(oce_2_atm%values)
	deallocate(oce_2_atm%colind)
	deallocate(oce_2_atm%rowptr)
	end subroutine clean_oce_2_atm
! =======================================================================
	Subroutine unstr2str(nx,ny)
	use o_mesh
	use o_arrays
	use o_param
	use g_parsup
	implicit none
	integer :: nx,ny,iz,tr
	real*8 :: dx,dy
	real*8 :: LonMin,LonMax,LatMin,LatMax
	real*8,dimension(:),allocatable ::  atm_lon, atm_lat
	real*8 t1,t2,t3,t4

	setup_u=setup_done
	setup_v=setup_done
	setup_sc=setup_done

	dx=(CplLonMax - CplLonMin) / real(nx)
	dy=(CplLatMax - CplLatMin) / real(ny)

!Scalars
        allocate(atm_lon(nx), atm_lat(ny))   ! discrete longitude and latitude

	LonMin=CplLonMin+dx*0.5
	LonMax=CplLonMax-dx*0.5
	LatMin=CplLatMin+dy*0.5
	LatMax=CplLatMax-dy*0.5
	t1=MPI_Wtime()
	call  def_exchange_mesh(atm_lon, nx,atm_lat,ny, &
		LonMin,LonMax,LatMin,LatMax)
	t2=MPI_Wtime()
	do iz=1,nl-1
	t3=MPI_Wtime()
		if (iz==1) then
			call build_oce_2_atm(atm_lon,nx,atm_lat,ny,0)
			call do_oce_2_atm(aveta_n, regssh, atm_lon,nx,atm_lat,ny,iz,0)
			!exit
		t4=MPI_Wtime()
		if (mype==0) write(*,*) '====================================================='
		if (mype==0) write(*,*) 'mesh',t2-t1
		if (mype==0) write(*,*) 'build',t3-t2
		if (mype==0) write(*,*) 'do',t4-t3
		if (mype==0) write(*,*) 'Total',t4-t1
		if (mype==0) write(*,*) '====================================================='
!			write(*,*) 'DBG'
			!call par_ex(1)
		endif
		do tr=1,num_tracers
		!	if (mype==0) write(*,*) 'aaa1',minval(avtr_arr(iz,:,tr)),maxval(avtr_arr(iz,:,tr))
			call do_oce_2_atm(avtr_arr(iz,:,tr), regtr_arr(:,:,iz,tr), atm_lon,nx,atm_lat,ny,iz,0)
		!	if (mype==0) write(*,*) 'aaa2',minval(regtr_arr(:,:,iz,tr)),maxval(regtr_arr(:,:,iz,tr))
		enddo
		!call clean_oce_2_atm
	enddo
	deallocate(atm_lon, atm_lat)
!U h. velocity
        allocate(atm_lon(nx+1), atm_lat(ny))   ! discrete longitude and latitude

	LonMin=CplLonMin
	LonMax=CplLonMax
	LatMin=CplLatMin+dy*0.5
	LatMax=CplLatMax-dy*0.5
	call  def_exchange_mesh(atm_lon, nx+1,atm_lat,ny, &
		LonMin,LonMax,LatMin,LatMax)
	do iz=1,nl-1
		if (iz==1) call build_oce_2_atm(atm_lon,nx+1,atm_lat,ny,1)
		call do_oce_2_atm(avUnode(iz,:), regu(:,:,iz), atm_lon,nx+1,atm_lat,ny,iz,1)
		!call clean_oce_2_atm
	enddo
	deallocate(atm_lon, atm_lat)
!V h. velocity
        allocate(atm_lon(nx), atm_lat(ny+1))   ! discrete longitude and latitude

	LonMin=CplLonMin+dx*0.5
	LonMax=CplLonMax-dx*0.5
	LatMin=CplLatMin
	LatMax=CplLatMax
	call  def_exchange_mesh(atm_lon, nx,atm_lat,ny+1, &
		LonMin,LonMax,LatMin,LatMax)
	do iz=1,nl-1
		if (iz==1) call build_oce_2_atm(atm_lon,nx,atm_lat,ny+1,-1)
		call do_oce_2_atm(avVnode(iz,:), regv(:,:,iz), atm_lon,nx,atm_lat,ny+1,iz,-1)
		!call clean_oce_2_atm
	enddo
	deallocate(atm_lon, atm_lat)
	call fix_boundary(nx,ny)
	end subroutine unstr2str
!==================================================
	subroutine fix_boundary(nx,ny)
! Mesh is supposed to be cyclic only in x-direction
	use o_mesh
	implicit none
	integer, intent(in) :: nx,ny
	integer i,j,k,iz,ii,jj
	integer flag(4)
	do i=1,nx
	do j=1,ny
	flag=0
	!if sum(flag)==4 then node is surrounded by land
	do iz=1,nl-1
		if (regtr_arr(i,j,iz,1)==0) then
		if (sum(flag)==4) then
		do k=1,2
		 regu(i-1+k,j,iz:nl-1)=-999
		 regv(i,j-1+k,iz:nl-1)=-999 
		enddo ! k
		exit
		else !flag
		do k=1,2
		 ii=i-1+2*(k-1)
		 if (ii<1) ii=nx
		 if (ii>nx) ii=1
		 if (regtr_arr(ii,j,iz,1)<-200) then
		 !if (regtr_arr(ii,j,iz,1)==0) then
		  regu(i-1+k,j,iz)=-999
		  flag(k)=1
		 else
		  regu(i-1+k,j,iz)=0d0
		 endif
		 jj=j-1+2*(k-1)
		 if ((jj<1).or.(jj>ny)) then
		  regv(i,j-1+k,iz)=-999
		  flag(2+k)=1
		  cycle
		 endif
		 !if (regtr_arr(i,jj,iz,1)==0) then
		 if (regtr_arr(i,jj,iz,1)<-200) then
		  regv(i,j-1+k,iz)=-999
		  flag(2+k)=1
		 else
		  regv(i,j-1+k,iz)=0d0
		 endif
		enddo! k
		endif !flag
		elseif ((j==1) .or. (j==ny)) then
		 if (j==ny) then
			regv(i,ny+1,iz)=0d0
		 else
			regv(i,1,iz)=0d0
		 endif
		endif !regtr_arr
	enddo
	enddo 
	enddo
	end subroutine fix_boundary
!======================================================================

subroutine get_omask(omask,field)
use g_parsup
implicit none
real*8 :: omask(:,:,:),field(:,:,:)
integer nx,ny,nz
integer i,j,k,l,ii,jj
logical water,land

nx=ubound(field,1)
ny=ubound(field,2)
nz=ubound(field,3)
!write(*,*) 'dfgsdgd',mype,ubound(field,1),ubound(omask,2),ubound(omask,3)
do k=1,nz
do i=1,nx
do j=1,ny
!	if (field(i,j,k)/=0) then
	if (field(i,j,k)==0d0) then
		omask(i,j,k)=0d0
	else
		omask(i,j,k)=1d0
	endif
enddo
enddo
enddo
!filter out islands
do k=1,nz
do i=1,nx
do j=1,ny
  water=.false.
  land=.false.
  do l=1,2 
    ii=i-1+2*(l-1)
    jj=j-1+2*(l-1)
    if (ii>nx) ii=1
    if (ii<1) ii=nx
    if (omask(ii,j,k)/=0d0) then
	water=.true.
    else
	land=.true.
    endif
    if ((jj<1).or.(jj>ny)) cycle
    if (omask(i,jj,k)/=0d0) then
	water=.true.
    else
	land=.true.
    endif
  enddo
  if ((omask(i,j,k)/=0d0) .and. (.not. water) ) omask(i,j,k)=0d0
  if ((omask(i,j,k)==0d0) .and. (.not. land) ) omask(i,j,k)=1d0
enddo
enddo
enddo
end subroutine get_omask
!======================================================================
subroutine apply_omask(omask,field,fill)
implicit none
real*8 :: omask(:,:,:),field(:,:,:),fill
integer nx,ny,nz
integer i,j,k
nx=ubound(field,1)
ny=ubound(field,2)
nz=ubound(field,3)
do k=1,nz
do i=1,nx
do j=1,ny
!       if (field(i,j,k)/=0) then
        if (omask(i,j,k)==0d0) field(i,j,k)=fill
enddo
enddo
enddo
end subroutine apply_omask
!==================================================
subroutine save_projector(mytype)
use g_parsup
use g_config
implicit none
character*10 :: temp
character*100 :: filename
integer i,n,mytype,fid
fid=10+mype
write(temp,'(i10)') mype
if (mytype>0) then
  filename=trim(MeshPath)//'dist/proj/proj_u.'//trim(ADJUSTL(temp))//'.out'
elseif (mytype<0) then
  filename=trim(MeshPath)//'dist/proj/proj_v.'//trim(ADJUSTL(temp))//'.out'
else
  filename=trim(MeshPath)//'dist/proj/proj_sc.'//trim(ADJUSTL(temp))//'.out'
endif
open(unit=fid,file=filename)
if (mytype>0) then
  write(fid,*) oce_2_atm_U%dim
  do i=1,oce_2_atm_U%dim+1
	write(fid,*) oce_2_atm_U%rowptr(i) 
  enddo
  n=oce_2_atm_u%rowptr(oce_2_atm_u%dim+1)-1
  do i=1,n
	write(fid,*) oce_2_atm_u%colind(i)
	write(fid,*) oce_2_atm_u%values(i)
  enddo
elseif (mytype<0) then
  write(fid,*) oce_2_atm_V%dim
  do i=1,oce_2_atm_V%dim+1
	write(fid,*) oce_2_atm_V%rowptr(i) 
  enddo
  n=oce_2_atm_V%rowptr(oce_2_atm_V%dim+1)-1
  do i=1,n
	write(fid,*) oce_2_atm_V%colind(i)
	write(fid,*) oce_2_atm_V%values(i)
  enddo
else
  write(fid,*) oce_2_atm%dim
  do i=1,oce_2_atm%dim+1
	write(fid,*) oce_2_atm%rowptr(i) 
  enddo
  n=oce_2_atm%rowptr(oce_2_atm%dim+1)-1
  do i=1,n
	write(fid,*) oce_2_atm%colind(i)
	write(fid,*) oce_2_atm%values(i)
  enddo
endif
close(fid)
end subroutine save_projector
!==================================================
subroutine load_projector(mytype)
use g_parsup
use g_config
implicit none
character*10 :: temp
character*100 :: filename
integer i,n,mytype,fid
fid=10+mype
write(temp,'(i10)') mype
if (mytype>0) then
  filename=trim(MeshPath)//'dist/proj/proj_u.'//trim(ADJUSTL(temp))//'.out'
elseif (mytype<0) then
  filename=trim(MeshPath)//'dist/proj/proj_v.'//trim(ADJUSTL(temp))//'.out'
else
  filename=trim(MeshPath)//'dist/proj/proj_sc.'//trim(ADJUSTL(temp))//'.out'
endif

open(unit=fid,file=filename)
if (mytype>0) then
  read(fid,*) oce_2_atm_U%dim
  allocate(oce_2_atm_U%rowptr(oce_2_atm_U%dim+1))
  do i=1,oce_2_atm_U%dim+1
        read(fid,*) oce_2_atm_U%rowptr(i)
  enddo
  n=oce_2_atm_u%rowptr(oce_2_atm_u%dim+1)-1
  allocate(oce_2_atm_u%values(n))
  allocate(oce_2_atm_u%colind(n))
  do i=1,n
        read(fid,*) oce_2_atm_u%colind(i)
        read(fid,*) oce_2_atm_u%values(i)
  enddo
elseif (mytype<0) then
  read(fid,*) oce_2_atm_V%dim
  allocate(oce_2_atm_V%rowptr(oce_2_atm_V%dim+1))
  do i=1,oce_2_atm_V%dim+1
        read(fid,*) oce_2_atm_V%rowptr(i)
  enddo
  n=oce_2_atm_V%rowptr(oce_2_atm_V%dim+1)-1
  allocate(oce_2_atm_V%values(n))
  allocate(oce_2_atm_V%colind(n))
  do i=1,n
        read(fid,*) oce_2_atm_V%colind(i)
        read(fid,*) oce_2_atm_V%values(i)
  enddo
else
  read(fid,*) oce_2_atm%dim
  allocate(oce_2_atm%rowptr(oce_2_atm%dim+1))
  do i=1,oce_2_atm%dim+1
        read(fid,*) oce_2_atm%rowptr(i)
  enddo
  n=oce_2_atm%rowptr(oce_2_atm%dim+1)-1
  allocate(oce_2_atm%values(n))
  allocate(oce_2_atm%colind(n))
  do i=1,n
        read(fid,*) oce_2_atm%colind(i)
        read(fid,*) oce_2_atm%values(i)
  enddo
endif
close(fid)
setup_done=.true.
end subroutine load_projector
!==================================================
	end module g_str2unstr
