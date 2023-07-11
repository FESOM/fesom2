!
!--------------------------------------------------------------------------
!
subroutine read_2Dmesh
  use g_config
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use g_oce_2_reg
  implicit none

  integer           :: i, n, ind
  real(kind=8)      :: x, y
  open (20,file=trim(MeshPath)//'nod2d.out',  status='old')
  open (21,file=trim(MeshPath)//'elem2d.out', status='old')
  write(*,*) '2D mesh is opened'

  read(20,*) nod2D 
  allocate(coord_nod2D(2, nod2D))
  do i=1, nod2D
     read(20,*) n, x, y, ind
     coord_nod2D(1, i)=x
     coord_nod2D(2, i)=y
  end do
  close(20)
  !
  read(21,*)  elem2D      
  allocate(elem2D_nodes(3,elem2D))  
  do n=1, elem2D
     read(21,*) elem2D_nodes(:,n)     
  end do
  close(21)
  
  
    open (20,file=trim(meshpath)//'nod2d.out',  status='old')
  open (21,file=trim(meshpath)//'elem2d.out', status='old')
  
  write(*,*) 'read_2Dmesh: DONE'
end subroutine read_2Dmesh
!
!---------------------------------------------------------------------------
!
subroutine read_3Dmesh
  use g_config
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use g_oce_2_reg  
  implicit none
  
  integer           :: i, j, n, node, ind
  real(kind=8)      :: x, y, z
  open(10,file=trim(meshpath)//'aux3d.out', status='old')
  open(11,file=trim(meshpath)//'nlvls.out', status='old')
  open(12,file=trim(meshpath)//'elvls.out', status='old')

  allocate(nlvls(nod2D), elvls(elem2D))
  write(*,*) '3D mesh is opened'
  ! Read auxilliary data
  read(10,*) max_num_layers
  nl_1=max_num_layers-1
  allocate(depths(max_num_layers))
  read(10,*) depths
  
  read(11,*) nlvls
  read(12,*) elvls

  close(10)
  close(11)
  close(12)
  write(*,*) 'read_3Dmesh: DONE'  
end subroutine read_3Dmesh
!
!---------------------------------------------------------------------------
!
