    subroutine build_oce_2_reg_el
        use g_oce_2_reg
        use o_mesh
        use g_rotate_grid	
        use o_elements	
        use o_param
        use g_config
        implicit none

        integer                   :: n, nij, i, j, pos, el, icol
        real(kind=8)              :: p(2)
        real(kind=8), allocatable :: xcoord(:), ycoord(:)
        integer, allocatable      :: do_lin_int(:), my_count(:)
	
! start building the oce_2_reg_el operator (to map elementwise field onto a regular mesh)
!1. Dimension=[reg_nx*reg_ny x elem2d]

        oce_2_reg_el%dim=reg_nx*reg_ny
	write(*,*)  'oce_2_reg_el%dim=', oce_2_reg_el%dim
	allocate(oce_2_reg_el%rowptr(oce_2_reg_el%dim+1))
	allocate(oce_2_reg_el%colind(oce_2_reg_el%dim))
	allocate(oce_2_reg_el%values(oce_2_reg_el%dim))
        oce_2_reg_el%values=0.
       	oce_2_reg_el%colind(:)=0
	oce_2_reg_el%rowptr(:)=0
	oce_2_reg_el%rowptr(1)=1

	do i=1, reg_nx	 
	   do j=1, reg_ny
	      nij=(i-1)*reg_ny+j
	      pos=nij+1
              p(1)=reg_lon(i)
              p(2)=reg_lat(j)
!2. There is only one entry per row (field is constant on elements)
              icol=oce_2_reg_el%rowptr(pos-1)
              oce_2_reg_el%rowptr(pos)=icol+1

              call point_in_triangle(el,   p)
              if (el > 0) then
		 oce_2_reg_el%colind(icol)=el
		 oce_2_reg_el%values=1.
	      else
		 oce_2_reg_el%colind(icol)=1
		 oce_2_reg_el%values=0.            
              end if
	    end do
         end do
	 write(*,*) 'build_oce_2_reg_el: DONE'

    end  subroutine build_oce_2_reg_el
! =======================================================================
subroutine save_oce_2_reg_el
  use g_oce_2_reg
  use o_mesh
  use g_rotate_grid	
  use o_elements	
  use o_param
  use g_config
  implicit none
  INTEGER                                   :: i

	
  open(1,file=trim(outpath)//'oce_2_reg_el.bin', form='binary')
  
  write(1)  oce_2_reg_el%dim, oce_2_reg_el%nza
  
  do i=1, oce_2_reg_el%dim+1
     write(1) oce_2_reg_el%rowptr(i)
  end do

  do i=1, oce_2_reg_el%nza
     write(1) oce_2_reg_el%colind(i), oce_2_reg_el%values(i)
  end do  
  close(1)	
end subroutine save_oce_2_reg_el
! =======================================================================
subroutine load_oce_2_reg_el 
  use g_oce_2_reg
  use o_mesh
  use g_rotate_grid	
  use o_elements	
  use o_param
  use g_config
  implicit none
  INTEGER                                   :: i

	
  open(1, file=trim(outpath)//'oce_2_reg_el.bin', form='binary', status='old')
  
  read(1)  oce_2_reg_el%dim, oce_2_reg_el%nza
  allocate(oce_2_reg_el%colind(oce_2_reg_el%nza), oce_2_reg_el%values(oce_2_reg_el%nza))
  allocate(oce_2_reg_el%rowptr(oce_2_reg_el%dim+1))
  do i=1, oce_2_reg_el%dim+1
     read(1) oce_2_reg_el%rowptr(i)
  end do

  do i=1, oce_2_reg_el%nza
     read(1) oce_2_reg_el%colind(i), oce_2_reg_el%values(i)
  end do  
  close(1)	
end subroutine load_oce_2_reg_el
! =======================================================================
subroutine do_oce_2_reg_el(oce, reg, flag)
     use g_oce_2_reg    
     use o_param
     use o_elements
     use o_mesh
     use g_rotate_grid

     implicit none
     save
	
     real(kind=8), intent( IN    )         :: oce(elem2D)
     real(kind=8), intent( INOUT )         :: reg(reg_nx, reg_ny)
     integer,      intent( IN    )         :: flag	
     integer                               :: i, j, is, ie, n

     ! Apply operator: A*x
     do i=1, reg_nx
	do j=1, reg_ny
      	   n=(i-1)*reg_ny+j
	   is=oce_2_reg_el%rowptr(n)
	   ie=oce_2_reg_el%rowptr(n+1)-1
	   reg(i, j)=sum(oce_2_reg_el%values(is:ie)*oce(oce_2_reg_el%colind(is:ie)))
	end do
     end do
end subroutine do_oce_2_reg_el
