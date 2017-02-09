! =======================================================================
      subroutine def_exchange_mesh
!       DESCRIPTION
!       This subroutine creates and allocates the intermediate regular grid
!       dependent on the target atmosphere grid (with lon and lat positions)
!       on which the unstructured grid is interpolated to before given to OASIS

        use g_oce_2_reg
        use o_mesh	
        use o_param
        use g_config
        implicit none
      
        integer                                :: i,j
        real(kind=8)                           :: eps=1.e-40

        reg_nx  = int((LonMax - LonMin+eps)/RegDx)
        reg_ny  = int((LatMax - LatMin+eps)/RegDy)

        allocate(reg_lon(reg_nx), reg_lat(reg_ny))   ! discrete longitude and latitude
        reg_lon(:) = 0.0
        reg_lat(:) = 0.0

        do i=1,reg_nx
          reg_lon(i)=min(LonMin+(i-1)*RegDx,LonMax)+RegDx/2.
        end do

        do i=1,reg_ny
          reg_lat(i)=min(LatMin+(i-1)*RegDy,LatMax)+RegDy/2.
        end do	

        ! =========
        ! Degrees to radians
        ! =========

        reg_lon=reg_lon*rad
        reg_lat=reg_lat*rad

        RegDx=RegDx*rad
        RegDy=RegDy*rad
	write(*,*) 'def_exchange_mesh: DONE'
    end subroutine def_exchange_mesh
! =======================================================================
subroutine def_mask
   use o_mesh	
   use g_config
   implicit none

   integer                         :: i, k, n

   allocate(mask_n2(nod2D))
   mask_n2=1.
   if (use_mask) then
      mask_n2=0.
      open(20,file=trim(mask_file),  status='old')
      read(20,*) n
      do i=1, n
         read(20,*) k
         mask_n2(k)=1.
      end do
      close(20)
   end if
end subroutine def_mask
! =======================================================================
    subroutine build_oce_2_reg    
        use g_oce_2_reg
        use o_mesh
        use g_rotate_grid	
        use o_elements	
        use o_param
        use g_config
        implicit none

        integer                   :: n, nij, i, j, pos, el, is, ie
        real(kind=8)              :: x, y, xr, yr, xx, yy, p(2), dx, dy, rlon, rlat
        integer, allocatable      :: do_lin_int(:), my_count(:)
        real(kind=8)              :: eps=1.e-40
	
	allocate(do_lin_int(reg_nx*reg_ny))
	allocate(my_count  (reg_nx*reg_ny))
	
	dx=RegDx/2.
	dy=RegDy/2.

	do_lin_int=0
	my_count  =0
		
! start building the atm_2_oce operator
!1. Dimension=[reg_nx*reg_ny x nod2d]

        oce_2_reg%dim=reg_nx*reg_ny
	write(*,*)  'oce_2_reg%dim=', oce_2_reg%dim
	allocate(oce_2_reg%rowptr(oce_2_reg%dim+1))
	oce_2_reg%rowptr(:)=0
	oce_2_reg%rowptr(1)=1
!2. rowptr
	do i=1, reg_nx
	write(*,*) 'rowptr, i=', i, ' of ', reg_nx	 
	   do j=1, reg_ny
	      nij=(i-1)*reg_ny+j
	      pos=nij+1
              xr=reg_lon(i)
              yr=reg_lat(j)
              if (rotated_grid) then
                 call g2r(reg_lon(i), reg_lat(j), xr, yr)
              endif
	      xr=xr/rad
              yr=yr/rad
              if (xr < -180) xr=xr+360
              if (xr >  180) xr=xr-360

 	      do n=1, nod2d
	         x=coord_nod2D(1,n)/rad
  	         y=coord_nod2D(2,n)/rad
	         if (x < -180) x=x+360
	         if (x >  180) x=x-360		 
		 xx=abs(xr-x)
		 xx=min(xx, abs(360-xx))
		 yy=abs(yr-y)
		 if (xx > dx  .OR. yy > dy) CYCLE
 	         my_count(nij)=my_count(nij)+1
  	      end do
! c > 0 means that the points around [i,j] were found. Otherwise we check 
! if there is a "rather big" triangle that [i,j] is inside it
	      if (my_count(nij) > 0) then
		 oce_2_reg%rowptr(pos)=oce_2_reg%rowptr(pos-1)+my_count(nij)
 	      else		 
	         p(1)=xr*rad
  	         p(2)=yr*rad
 	         call point_in_triangle(el, p)
	         if ((el > 0)) then
! we will do a linear interpolation to this point
		    oce_2_reg%rowptr(pos)=oce_2_reg%rowptr(pos-1)+3
		    do_lin_int(nij)      =el		     
		 else
		    oce_2_reg%rowptr(pos)=oce_2_reg%rowptr(pos-1)+1
		 end if
	      end if	      
	    end do
         end do

	write(*,*) 'build_oce_2_reg/rowptr: DONE'
!3. colind & values
        n=oce_2_reg%rowptr(oce_2_reg%dim+1)-1
	allocate(oce_2_reg%values(n))
	allocate(oce_2_reg%colind(n))
	oce_2_reg%colind=1
	oce_2_reg%values=0.
	oce_2_reg%nza   =n
        do i=1, reg_nx
	write(*,*) 'values, i=', i, ' of ', reg_nx
	   do j=1, reg_ny
   	      nij=(i-1)*reg_ny+j

              if (do_lin_int(nij)==0 .AND. my_count(nij)==0) CYCLE !do nothing
              pos=oce_2_reg%rowptr(nij)	      

              xr=reg_lon(i)
              yr=reg_lat(j)
              if (rotated_grid) then
                 call g2r(reg_lon(i), reg_lat(j), xr, yr)
              endif
	      xr=xr/rad
              yr=yr/rad
              if (xr < -180) xr=xr+360
              if (xr >  180) xr=xr-360

!fill colind and values, CASE 1	      
	      if (my_count(nij) > 0) then
!	         write(*,*) 'oce_2_reg: filling the values, CASE 1'	      	      
!local loop to fill the operator
                 do n=1, nod2d
	              x=coord_nod2D(1,n)/rad
      	              y=coord_nod2D(2,n)/rad
	              if (x < -180) x=x+360
  	              if (x >  180) x=x-360		 
		      xx=abs(xr-x)
		      xx=min(xx, abs(360-xx))
	  	      yy=abs(yr-y)
		      if (xx > dx .OR. yy > dy) CYCLE
		      oce_2_reg%colind(pos)=n
                      !inverse distance weighting will be used
		      oce_2_reg%values(pos)=1./sqrt((xx*cos(yy*rad))**2 + yy**2+eps) !1./real(my_count(nij))		      
		      pos=pos+1
 	         end do
                 is=oce_2_reg%rowptr(nij)
                 ie=oce_2_reg%rowptr(nij+1)-1
                 !normalize inverse distance to sum into one
                 oce_2_reg%values(ie:ie)=oce_2_reg%values(ie:ie)/sum(oce_2_reg%values(ie:ie))
              else if (do_lin_int(nij)>0) then !do linear interpolation to [i, j]
!fill colind and values, CASE 2
!	         write(*,*) 'oce_2_reg: filling the values, CASE 2'
      		 el=do_lin_int(nij)
                 p(1)=xr
                 p(2)=yr
                 oce_2_reg%colind(pos:pos+2)=elem2D_nodes(:, el)
                 call locbafu_2D(oce_2_reg%values(pos:pos+2), el, p)
	      end if
        end do
        end do
write(*,*) 'oce_2_reg:', minval(oce_2_reg%values), maxval(oce_2_reg%values)
	deallocate(do_lin_int)
	deallocate(my_count)
    end  subroutine build_oce_2_reg
! =======================================================================
SUBROUTINE point_in_triangle(el2D,   pt)
  use o_param
  use o_elements
  use o_mesh
  use g_rotate_grid

  implicit none

  INTEGER, INTENT(OUT)                   :: el2D  
  REAL(kind=8), DIMENSION(2), INTENT(IN) :: pt
  
  integer                                :: elem, elnodes(3), q
  real(kind=8)                           :: alpha, mean_lon, mean_lat, rlon, rlat
  real(kind=8)                           :: xe(4), ye(4), xt1, xt2, yt1, yt2, x1, x2, y1, y2
  real(kind=8)                           :: s1, s2, angle
  real(kind=8)                           :: dx(3), dy(3)
  el2D=0

  DO elem=1, elem2D
     elnodes(1:3)=elem2D_nodes(:, elem)

     xe(1:3)=coord_nod2D(1, elnodes(1:3))
     ye(1:3)=coord_nod2D(2, elnodes(1:3))

     dx=pt(1)-xe(1:3)
     dy=pt(2)-ye(1:3)
    
     where (dx> pi) 
           dx=(dx-2.*pi)
     end where

     where (dx<-pi) 
           dx=(dx+2.*pi)
     end where

     if (all(dx*dx(1) > 0.)) cycle
     if (all(dy*dy(1) > 0.)) cycle

	    !=====
            ! Cyclicity:
	    ! Remove 2.*pi jumps in x-coordinate
	    ! of nodes in triangle  
            ! =====
	    xe(2)=xe(2)-xe(1)
	    xe(3)=xe(3)-xe(1)
            DO q=2,3
	     if(xe(q)> pi) xe(q)=xe(q)-2.*pi
             if(xe(q)<-pi) xe(q)=xe(q)+2.*pi
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
	    if(xt1> pi) xt1=xt1-2.*pi 
            if(xt1<-pi) xt1=xt1+2.*pi  
            
	    x1=pt(1)-mean_lon
	    if(x1>pi) x1=x1-2.*pi 
            if(x1<-pi) x1=x1+2.*pi  
	    
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
               if(x1>pi) x1=x1-2.*pi 
               if(x1<-pi) x1=x1+2.*pi 
               if(x2>pi) x2=x2-2.*pi 
               if(x2<-pi) x2=x2+2.*pi 

               alpha=x1*y2-x2*y1
               s1=x1*x1+y1*y1
               s2=x2*x2+y2*y2
               alpha=alpha/sqrt(s1*s2)
               alpha=asin(alpha)
               IF ((xt1-xt2)**2+(yt1-yt2)**2 > max(s1,s2)) THEN
               if (alpha>0) alpha=pi-alpha
               if (alpha<=0) alpha=-pi-alpha
               END IF
               angle=angle+alpha
             END DO
             IF (abs(angle)>pi) THEN
             el2D=elem
	     exit
             END IF
         ENDDO 
END SUBROUTINE point_in_triangle
! =======================================================================
SUBROUTINE locbafu_2D(values, elem, coords)
  USE o_mesh
  USE o_elements
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
! =======================================================================
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
! =======================================================================
subroutine save_oce_2_reg    
  use g_oce_2_reg
  use o_mesh
  use g_rotate_grid	
  use o_elements	
  use o_param
  use g_config
  implicit none
  INTEGER                                   :: i

	
  open(1,file=trim(o2r_filename), form='binary')
  
  write(1)  oce_2_reg%dim, oce_2_reg%nza
  
  do i=1, oce_2_reg%dim+1
     write(1) oce_2_reg%rowptr(i)
  end do

  do i=1, oce_2_reg%nza
     write(1) oce_2_reg%colind(i), oce_2_reg%values(i)
  end do  
  close(1)	
end subroutine save_oce_2_reg
! =======================================================================
subroutine load_oce_2_reg    
  use g_oce_2_reg
  use o_mesh
  use g_rotate_grid	
  use o_elements	
  use o_param
  use g_config
  implicit none
  INTEGER                                   :: i

	
  open(1, file=trim(o2r_filename), form='binary', status='old')
  
  read(1)  oce_2_reg%dim, oce_2_reg%nza
  allocate(oce_2_reg%colind(oce_2_reg%nza), oce_2_reg%values(oce_2_reg%nza))
  allocate(oce_2_reg%rowptr(oce_2_reg%dim+1))
  do i=1, oce_2_reg%dim+1
     read(1) oce_2_reg%rowptr(i)
  end do

  do i=1, oce_2_reg%nza
     read(1) oce_2_reg%colind(i), oce_2_reg%values(i)
  end do  
  close(1)	
write(*,*) 'oce_2_reg statistics:'
write(*,*) 'oce_2_reg%dim        =', oce_2_reg%dim
write(*,*) 'max(oce_2_reg%colind)=', maxval(oce_2_reg%colind)
end subroutine load_oce_2_reg
! =======================================================================
    subroutine do_oce_2_reg(oce, reg, flag)
        use g_oce_2_reg    
        use o_param
        use o_elements
        use o_mesh
        use g_rotate_grid

        implicit none
	
	real(kind=8), intent( IN    )         :: oce(nod2D)
	real(kind=8), intent( INOUT )         :: reg(reg_nx, reg_ny)
	integer,      intent( IN    )         :: flag	
        integer                               :: i, j, is, ie, n

        ! Apply operator: A*x
	do i=1, reg_nx
           !write(*,*) 'i=', i, ' from ', reg_nx
	   do j=1, reg_ny
      	      n=(i-1)*reg_ny+j
	      is=oce_2_reg%rowptr(n)
	      ie=oce_2_reg%rowptr(n+1)-1
	      reg(i, j)=sum(oce_2_reg%values(is:ie)*oce(oce_2_reg%colind(is:ie)))
	   end do
	end do
    end subroutine do_oce_2_reg
