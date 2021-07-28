module read_mesh_interface
  interface
    subroutine read_mesh(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module find_levels_interface
  interface
    subroutine find_levels(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module find_levels_cavity_interface
  interface
    subroutine find_levels_cavity(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module test_tri_interface
  interface
    subroutine test_tri(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module load_edges_interface
  interface
    subroutine load_edges(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module find_neighbors_interface
  interface
    subroutine find_neighbors(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module mesh_areas_interface
  interface
    subroutine mesh_areas(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module elem_center_interface
  interface
    subroutine elem_center(elem, x, y, mesh)
      use mod_mesh
      integer      :: elem    
      real(kind=WP) :: x, y
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module edge_center_interface
  interface
    subroutine edge_center(n1, n2, x, y, mesh)
      USE MOD_MESH
      integer                     :: n1, n2
      real(kind=WP)               :: x, y
      type(t_mesh), intent(inout), target :: mesh
    end subroutine
  end interface
end module
module mesh_auxiliary_arrays_interface
  interface
    subroutine mesh_auxiliary_arrays(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module find_levels_min_e2n_interface
  interface
    subroutine find_levels_min_e2n(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module
module check_total_volume_interface
  interface
    subroutine check_total_volume(mesh)
      use mod_mesh
      type(t_mesh), intent(inout)  , target :: mesh
    end subroutine
  end interface
end module

! Driving routine. The distributed mesh information and mesh proper 
! are read from files.
! Auxiliary arrays with mesh information are assembled.
! At the beginning of each routine I list arrays it initializes.
! Array sizes vary (sometimes we need only myDim, yet sometimes more)! 
! S. Danilov, 2012
SUBROUTINE mesh_setup(mesh)
USE MOD_MESH
USE g_parsup
USE g_ROTATE_grid
use read_mesh_interface
use find_levels_interface
use find_levels_cavity_interface
use mesh_auxiliary_arrays_interface
use test_tri_interface
use load_edges_interface
use find_levels_min_e2n_interface
use find_neighbors_interface
use mesh_areas_interface
IMPLICIT NONE

      type(t_mesh), intent(inout) :: mesh

      call set_mesh_transform_matrix  !(rotated grid)
      call read_mesh(mesh)
      call set_par_support(mesh)
!!PS       call find_levels(mesh)
!!PS       
!!PS       if (use_cavity) call find_levels_cavity(mesh)
!!PS         
      call test_tri(mesh)
      call load_edges(mesh)
      call find_neighbors(mesh)
      
      call find_levels(mesh)
      if (use_cavity) call find_levels_cavity(mesh)
      
      call find_levels_min_e2n(mesh)
      call mesh_areas(mesh)
      call mesh_auxiliary_arrays(mesh)
           
END SUBROUTINE mesh_setup
!======================================================================
! Reads distributed mesh
! The mesh will be read only by 0 proc and broadcasted to the others.
SUBROUTINE read_mesh(mesh)
USE o_PARAM
USE g_CONFIG
USE MOD_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_rotate_grid 
IMPLICIT NONE

type(t_mesh), intent(inout), target :: mesh

 integer        :: n, nn, k, m, fileID
 integer        :: error_status !0/1=no error/error
 integer        :: vert_nodes(1000)
 integer        :: nchunk, chunk_size, ipos, iofs, mesh_check, flag_checkisrot, flag_checkmustrot
 real(kind=WP)  :: x, y, rx, ry, xp, yp, dbox, offset
 real(kind=WP)  :: t0, t1
 character*10   :: mype_string,npes_string
 character(len=MAX_PATH)  :: file_name
 character(len=MAX_PATH)  :: dist_mesh_dir
 integer        :: flag_wrongaux3d=0
 integer       :: ierror              ! return error code
 integer, allocatable, dimension(:)        :: mapping
 integer, allocatable, dimension(:,:)      :: ibuff
 real(kind=WP), allocatable, dimension(:,:) :: rbuff
 integer, allocatable, dimension(:,:)      :: auxbuff ! will be used for reading aux3d.out 
 integer fileunit, iostat
 character(32) mesh_checksum

!NR Cannot include the pointers before the targets are allocated...
!NR #include "associate_mesh.h"

  !mesh related files will be read in chunks of chunk_size
  chunk_size=100000
  !==============================
  ! Allocate mapping array (chunk_size)
  ! It will be used for several purposes 
  !==============================
  allocate(mapping(chunk_size))
  allocate(ibuff(chunk_size,4), rbuff(chunk_size,3))

  mapping=0 
  !==============================
  t0=MPI_Wtime()
  write(mype_string,'(i5.5)') mype  
  write(npes_string,"(I10)") npes
  dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'
 
  !=======================
  ! rank partitioning vector
  ! will be read by 0 proc
  !=======================
  if (mype==0) then
     file_name=trim(dist_mesh_dir)//'rpart.out'
     fileID=10
     open(fileID, file=trim(file_name)) 
     allocate(part(npes+1))
     read(fileID,*) n
     error_status=0
     if (n/=npes) error_status=1 !set the error status for consistency in rpart
     part(1)=1
     read(fileID,*) part(2:npes+1)
     DO n=2, npes+1
        part(n)=part(n-1)+part(n)
     END DO
     close(fileID)
  end if
  ! check the error status
  call MPI_BCast(error_status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (error_status/=0) then
     write(*,*) n
     write(*,*) 'error: NPES does not coincide with that of the mesh'
     call par_ex(1)
     STOP
  end if
  ! broadcasting partitioning vector to the other procs
  if (mype/=0) then
     allocate(part(npes+1))
  end if
  call MPI_BCast(part, npes+1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (mype==0) write(*,*) mype,'rpart is read'

  !===========================
  ! Lists of nodes and elements in global indexing. 
  ! every proc reads its file
  !===========================
 
  file_name=trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
  fileID=103+mype  !skip unit range 100--102
    
  open(fileID, file=trim(file_name))
  read(fileID,*) n
 
  read(fileID,*) myDim_nod2D
  read(fileID,*) eDim_nod2D
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D)) 	 
  read(fileID,*) myList_nod2D
	 
  read(fileID,*) myDim_elem2D
  read(fileID,*) eDim_elem2D
  read(fileID,*) eXDim_elem2D
  allocate(myList_elem2D(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  read(fileID,*) myList_elem2D
	
  read(fileID,*) myDim_edge2D
  read(fileID,*) eDim_edge2D
  allocate(myList_edge2D(myDim_edge2D+eDim_edge2D))
  read(fileID,*) myList_edge2D ! m

  close(fileID)
  if (mype==0) write(*,*) 'myLists are read'

  !==============================
  ! read 2d node data
  !==============================
  ! read the nod2D from nod2d.out and check whether it is equal to part(npes+1)-1
  mesh%nod2D=part(npes+1)-1
  allocate(mesh%coord_nod2D(2,myDim_nod2D+eDim_nod2D))
  if (mype==0) then
    file_name=trim(meshpath)//'nod2d.out'
    open(fileID, file=file_name)
    read(fileID,*) n      ! nod2D, we know it already
     error_status=0
     if (n/=mesh%nod2D) error_status=1 !set the error status for consistency between rpart and nod2D
    write(*,*) 'reading '// trim(file_name)   
  end if
  ! check the error status
  call MPI_BCast(error_status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (error_status/=0) then
     write(*,*) n
     write(*,*) 'error: nod2D/=part(npes+1)-1'
     call par_ex(1)
     STOP
  end if

  ! 0 proc reads the data in chunks and distributes it between other procs
  mesh_check=0
  
  ! new rotation pole usually (-40.0°,75.0°)
  call r2g(xp, yp, 0.0_WP*rad, 90.0_WP*rad)
  xp = xp/rad
  yp = yp/rad
  dbox = 2.5_WP
  
  flag_checkisrot=0
  flag_checkmustrot=1
  do nchunk=0, (mesh%nod2D-1)/chunk_size
     !create the mapping for the current chunk
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     !read the chunk into the buffers
     k=min(chunk_size, mesh%nod2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1), rbuff(n,1), rbuff(n,2), ibuff(n,2)
            !___________________________________________________________________
            ! put here offset parameter so that the checkrotation flag also works 
            ! for the pi mesh an others there the points around -40 longitude are 
            ! shifted by 360° --> therefor shift them back with offset parameter
            offset=0.0_WP
            if (rbuff(n,1)>180.0_WP) then 
                offset=-360.0_WP
            elseif (rbuff(n,1)<-180.0_WP) then    
                offset=360.0_WP
            end if    
            !___________________________________________________________________
            ! check if input mesh is already rotated --> force_rotation flag == .False.
            if (force_rotation .and. & 
               (rbuff(n,1)+offset>=xp-dbox .and. rbuff(n,1)+offset<=xp+dbox .and. & 
                rbuff(n,2)       >=yp-dbox .and. rbuff(n,2)       <=yp+dbox)) then
                flag_checkisrot = 1
            !___________________________________________________________________
            ! check if input mesh is already unrotated --> force_rotation flag ==.True.
            elseif ((.not. force_rotation) .and. & 
               (rbuff(n,1)+offset>=xp-dbox .and. rbuff(n,1)+offset<=xp+dbox .and. & 
                rbuff(n,2)       >=yp-dbox .and. rbuff(n,2)       <=yp+dbox)) then
                flag_checkmustrot = 0
            end if    
        end do
     end if
     call MPI_BCast(rbuff(1:k,1), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(rbuff(1:k,2), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     ! fill the local arrays
     do n=1, k
        x=rbuff(n,1)*rad
        y=rbuff(n,2)*rad
        !_______________________________________________________________________
        ! rotate grid points if force_rotation flag .True.
        if (force_rotation) then
           rx=x
           ry=y
           call g2r(rx, ry, x, y)
           
        end if        
        if (mapping(n)>0) then
           mesh_check=mesh_check+1
           mesh%coord_nod2D(1,mapping(n))=x
           mesh%coord_nod2D(2,mapping(n))=y
        end if
     end do
  end do
  if (mype==0) close(fileID)
    
    !___________________________________________________________________________
    ! check if rotation is applied to an already rotated mesh
    if ((mype==0) .and. (force_rotation) .and. (flag_checkisrot==1) .and. (.not. toy_ocean)) then
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: Your input mesh seems to be rotated and you try to' 
        write(*,*) '        rotate it again in FESOM (force_rotation=.true. ) !'
        write(*,*) '        The mesh you loaded suggests that it is already'
        write(*,*) '        rotated because it has ocean points in a box'
        write(*,*) '        around lon,lat = (-40.0, 75.0) (new rotation pole).'
        write(*,*) '        In an unrotated grid should be NO OCEAN points'
        write(*,*) '        over Greenland !!!. So set in the namelist.config '
        write(*,*) '        force_rotation==.False. If I am wrong go to'
        write(*,*) '        oce_mesh.F90:line231 and comment this part'
        write(*,*)
        write(*,*) '        --> check your namelist.config !!!'
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
        write(*,*)
        call par_ex(0)
    !___________________________________________________________________________
    ! check if rotation needs to be applied to an unrotated mesh
    elseif ((mype==0) .and. (.not. force_rotation) .and. (flag_checkmustrot==1) .and. (.not. toy_ocean)) then
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: Your input mesh seems to be unrotated this requires'
        write(*,*) '        that it is rotated in FESOM, but you set force_rotation=.False'
        write(*,*) '        The mesh you loaded suggests that it needs to'
        write(*,*) '        be rotated because it has no ocean points in a box'
        write(*,*) '        around lon,lat = (-40.0, 75.0) (new rotation pole).'
        write(*,*) '        A rotated grid should have some ocean points around'
        write(*,*) '        the rotation pole !!!. So set in namelist.config '
        write(*,*) '        force_rotation=.True. If I am wrong go to'
        write(*,*) '        oce_mesh.F90:line248 and comment this part'
        write(*,*)
        write(*,*) '        --> check your namelist.config !!!'
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
        write(*,*)
        call par_ex(0)
    end if
  
  
  if (mesh_check/=myDim_nod2D+eDim_nod2D) then
     write(*,*) 'ERROR while reading nod2d.out on mype=', mype
     write(*,*) mesh_check, ' values have been read in according to partitioning'
     write(*,*) 'it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
  end if

  !==============================
  ! read 2d elem data
  !==============================
  ! read the elem2D from elem2d.out
  if (mype==0)  then 
     file_name=trim(meshpath)//'elem2d.out'
     open(fileID, file=file_name)
     read(fileID,*) mesh%elem2d
     write(*,*) 'reading '// trim(file_name)   
  end if
  call MPI_BCast(mesh%elem2d, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  allocate(mesh%elem2D_nodes(3, myDim_elem2D+eDim_elem2D+eXDim_elem2D))

  ! 0 proc reads the data in chunks and distributes it between other procs
  do nchunk=0, (mesh%elem2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
        ipos=(myList_elem2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, mesh%elem2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1), ibuff(n,2), ibuff(n,3)
        end do
     end if

     call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,3), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

     do n=1, k
        if (mapping(n)>0) then
           mesh%elem2D_nodes(1,mapping(n))=ibuff(n,1)
           mesh%elem2D_nodes(2,mapping(n))=ibuff(n,2)
           mesh%elem2D_nodes(3,mapping(n))=ibuff(n,3)
        end if
     end do
  end do
  if (mype==0) close(fileID)
  ! nodes in elem2d are in global numbering. convert to local:
  do nchunk=0, (mesh%nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
        do m=1,3
           nn=mesh%elem2D_nodes(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk) then
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              mesh%elem2D_nodes(m,n)=-mapping(iofs) 
           end if
        end do
     end do
  end do
  mesh%elem2D_nodes=-mesh%elem2D_nodes
 if (mype==0) write(*,*) 'elements are read' 
 !==============================
 ! read depth data
 !==============================
 ! 0 proc reads header of aux3d.out and broadcasts it between procs
 allocate(mesh%depth(myDim_nod2D+eDim_nod2D))
 if (mype==0) then !open the file for reading on 0 proc
    file_name=trim(meshpath)//'aux3d.out' 
    open(fileID, file=file_name)
    read(fileID,*) mesh%nl  ! the number of levels 
 end if
 call MPI_BCast(mesh%nl, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
 if (mesh%nl < 3) then
    write(*,*) '!!!Number of levels is less than 3, model will stop!!!'
    call par_ex
    stop
 end if
 allocate(mesh%zbar(mesh%nl))              ! allocate the array for storing the standard depths
 if (mype==0) read(fileID,*) mesh%zbar
 call MPI_BCast(mesh%zbar, mesh%nl, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
 if(mesh%zbar(2)>0) mesh%zbar=-mesh%zbar   ! zbar is negative 
 allocate(mesh%Z(mesh%nl-1))
 mesh%Z=mesh%zbar(1:mesh%nl-1)+mesh%zbar(2:mesh%nl)  ! mid-depths of cells
 mesh%Z=0.5_WP*mesh%Z 

 ! 0 proc reads the data in chunks and distributes it between other procs
 mesh_check=0
 do nchunk=0, (mesh%nod2D-1)/chunk_size
    mapping(1:chunk_size)=0
    do n=1, myDim_nod2D+eDim_nod2D
       ipos=(myList_nod2D(n)-1)/chunk_size
       if (ipos==nchunk) then
          iofs=myList_nod2D(n)-nchunk*chunk_size
          mapping(iofs)=n
       end if
    end do

    k=min(chunk_size, mesh%nod2D-nchunk*chunk_size)
    if (mype==0) then
       do n=1, k
          read(fileID,*) rbuff(n,1)
       end do
       ! check here if aux3d.out contains depth levels (FESOM2.0) or 3d indices
       ! (FESOM1.4) like that check if the proper mesh is loaded. 11000.0 is here 
       ! the maximum depth on earth in marianen trench
       if ( flag_wrongaux3d==0 .and. any(abs(rbuff(1:k,1))>11000.0_WP) ) flag_wrongaux3d=1
    end if
    call MPI_BCast(rbuff(1:k,1), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

    do n=1, k
       x=rbuff(n,1)
       if (x>0) x=-x !deps must be negative!
       if (x>mesh%zbar(5)) x=mesh%zbar(5) !threshold for depth
       if (mapping(n)>0) then
          mesh_check=mesh_check+1
          mesh%depth(mapping(n))=x
        end if
     end do
 end do

 if (mype==0) close(fileID)
 if (mesh_check/=myDim_nod2D+eDim_nod2D) then
    write(*,*) 'ERROR while reading aux3d.out on mype=', mype
    write(*,*) mesh_check, ' values have been read in according to partitioning'
    write(*,*) 'it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
 end if
!_______________________________________________________________________________
! check if the mesh structure of FESOM2.0 and of FESOM1.4 is loaded
if ((mype==0) .and. (flag_wrongaux3d==1)) then
    write(*,*)
    print *, achar(27)//'[33m'
    write(*,*) '____________________________________________________________________'
    write(*,*) ' ERROR: It looks like the mesh you want to use is prepared for ' 
    write(*,*) '        FESOM1.4. Please be aware that the input mesh structure'
    write(*,*) '        between FESOM1.4 and FESOM2.0 is different! The input'
    write(*,*) '        mesh structure of FESOM2.0 contains only the files nod2d.out'
    write(*,*) '        elem2d.out, nlvls.out, elvls.out and aux3d.out, where'
    write(*,*) '        aux3d.out contains the number of depth layers, the depth'
    write(*,*) '        layers and the bottom depth of your mesh vertices.'
    write(*,*) '        So please check your meshpath in namelist.config or your '
    write(*,*) '        mesh directory itself so you use the proper mesh structure'
    write(*,*)
    write(*,*) '____________________________________________________________________'
    print *, achar(27)//'[0m'
    write(*,*)
    call par_ex(0)
end if 

 ! ==============================
 ! Communication information
 ! every proc reads its file
 ! ==============================
 file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out'  
 fileID=103+mype !skip unit range 100--102  
 open(fileID, file=file_name)
 read(fileID,*)  n
 read(fileID,*) com_nod2D%rPEnum
 if (com_nod2D%rPEnum > MAX_NEIGHBOR_PARTITIONS) then
    print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
    stop
 endif
!!$ ALLOCATE(com_nod2D%rPE(com_nod2D%rPEnum))
 read(fileID,*) com_nod2D%rPE(1:com_nod2D%rPEnum)
!!$  ALLOCATE(com_nod2D%rptr(com_nod2D%rPEnum+1))
 read(fileID,*) com_nod2D%rptr(1:com_nod2D%rPEnum+1)
 ALLOCATE(com_nod2D%rlist(eDim_nod2D))
 read(fileID,*) com_nod2D%rlist
	 
 read(fileID,*) com_nod2D%sPEnum
 if (com_nod2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
    print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
    stop
 endif
!!$  ALLOCATE(com_nod2D%sPE(com_nod2D%sPEnum))
 read(fileID,*) com_nod2D%sPE(1:com_nod2D%sPEnum)
!!$  ALLOCATE(com_nod2D%sptr(com_nod2D%sPEnum+1))
 read(fileID,*) com_nod2D%sptr(1:com_nod2D%sPEnum+1)
 n=com_nod2D%sptr(com_nod2D%sPEnum+1)-1
 ALLOCATE(com_nod2D%slist(n))
 read(fileID,*) com_nod2D%slist
	 
 read(fileID,*) com_elem2D%rPEnum
 if (com_elem2D%rPEnum > MAX_NEIGHBOR_PARTITIONS) then
    print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
    stop
 endif
!!$  ALLOCATE(com_elem2D%rPE(com_elem2D%rPEnum))
 read(fileID,*) com_elem2D%rPE(1:com_elem2D%rPEnum)
!!$  ALLOCATE(com_elem2D%rptr(com_elem2D%rPEnum+1))
 read(fileID,*) com_elem2D%rptr(1:com_elem2D%rPEnum+1)
 ALLOCATE(com_elem2D%rlist(eDim_elem2D))
 read(fileID,*) com_elem2D%rlist
	 
 read(fileID,*) com_elem2D%sPEnum
 if (com_elem2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
    print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
    stop
 endif
!!$  ALLOCATE(com_elem2D%sPE(com_elem2D%sPEnum))
 read(fileID,*) com_elem2D%sPE(1:com_elem2D%sPEnum)
!!$  ALLOCATE(com_elem2D%sptr(com_elem2D%sPEnum+1))
 read(fileID,*) com_elem2D%sptr(1:com_elem2D%sPEnum+1)
 n=com_elem2D%sptr(com_elem2D%sPEnum+1)-1
 ALLOCATE(com_elem2D%slist(n))
 read(fileID,*) com_elem2D%slist
	 
 read(fileID,*) com_elem2D_full%rPEnum
 if (com_elem2D_full%rPEnum > MAX_NEIGHBOR_PARTITIONS) then
    print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
    stop
 endif
!!$  ALLOCATE(com_elem2D_full%rPE(com_elem2D_full%rPEnum))
 read(fileID,*) com_elem2D_full%rPE(1:com_elem2D_full%rPEnum)
!!$  ALLOCATE(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1))
 read(fileID,*) com_elem2D_full%rptr(1:com_elem2D_full%rPEnum+1)
 ALLOCATE(com_elem2D_full%rlist(eDim_elem2D+eXDim_elem2D))
 read(fileID,*) com_elem2D_full%rlist
	 
 read(fileID,*) com_elem2D_full%sPEnum
 if (com_elem2D_full%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
    print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
    stop
 endif
!!$  ALLOCATE(com_elem2D_full%sPE(com_elem2D_full%sPEnum))
 read(fileID,*) com_elem2D_full%sPE(1:com_elem2D_full%sPEnum)
!!$  ALLOCATE(com_elem2D_full%sptr(com_elem2D_full%sPEnum+1))
 read(fileID,*) com_elem2D_full%sptr(1:com_elem2D_full%sPEnum+1)
 n=com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)-1
 ALLOCATE(com_elem2D_full%slist(n))
 read(fileID,*) com_elem2D_full%slist

!!$ read(fileID,*) com_edge2D%rPEnum
!!$ ALLOCATE(com_edge2D%rPE(com_edge2D%rPEnum))
!!$ read(fileID,*) com_edge2D%rPE
!!$ ALLOCATE(com_edge2D%rptr(com_edge2D%rPEnum+1))
!!$ read(fileID,*) com_edge2D%rptr
!!$ ALLOCATE(com_edge2D%rlist(eDim_edge2D))
!!$ read(fileID,*) com_edge2D%rlist
!!$	 
!!$ read(fileID,*) com_edge2D%sPEnum
!!$ ALLOCATE(com_edge2D%sPE(com_edge2D%sPEnum))
!!$ read(fileID,*) com_edge2D%sPE
!!$ ALLOCATE(com_edge2D%sptr(com_edge2D%sPEnum+1))
!!$ read(fileID,*) com_edge2D%sptr
!!$ n=com_edge2D%sptr(com_edge2D%sPEnum+1)-1
!!$ ALLOCATE(com_edge2D%slist(n))
!!$ read(fileID,*) com_edge2D%slist
 close(fileID)
 if (mype==0) write(*,*) 'communication arrays are read'
 deallocate(rbuff, ibuff)
 deallocate(mapping)
 
  ! try to calculate checksum and distribute it to every process
  ! the shell command is probably not very portable and might fail, in which case we just do not have a checksum
  mesh%representative_checksum = '                                ' ! we use md5 which is 32 chars long, so set default value to the same length
  if(mype==0) then
    call execute_command_line("md5sum "//trim(MeshPath)//"nod2d.out | cut -d ' ' -f 1 > "//trim(ResultPath)//"mesh_checksum")
    ! we can not check if execute_command_line succeeded (e.g. with cmdstat), as the pipe will swallow any error from the initial command
    ! so we have to thoroughly check if the file exists and if it contains our checksum
    open(newunit=fileunit, file=trim(ResultPath)//"mesh_checksum", action="READ", iostat=iostat)
    if(iostat==0) read(fileunit, *, iostat=iostat) mesh_checksum
    close(fileunit)      
    if(iostat==0 .and. len_trim(mesh_checksum)==32) mesh%representative_checksum = mesh_checksum
  end if
  call MPI_BCAST(mesh%representative_checksum, len(mesh%representative_checksum), MPI_CHAR, 0, MPI_COMM_FESOM, MPIerr)
  mesh%representative_checksum = trim(mesh%representative_checksum) ! if we did not get a checksum, the string is empty

CALL MPI_BARRIER(MPI_COMM_FESOM, MPIerr)
 t1=MPI_Wtime()
 if (mype==0) then
    write(*,*) '========================='
    write(*,*) '2D mesh was read in ', t1-t0, ' seconds'
    write(*,*) '2D mesh info : ', 'nod2D=', mesh%nod2D,' elem2D=', mesh%elem2D,'checksum= ',mesh%representative_checksum
    write(*,*) '========================='
 endif

 END subroutine  read_mesh
!
!
!_______________________________________________________________________________
! load fesom2.0 mesh files: nlvls.out and elvls.out that are created during the 
! partitioning
!_______________________________________________________________________________
subroutine find_levels(mesh)
    use MOD_MESH
    use o_PARAM
    use g_PARSUP
    use g_config
    !
    implicit none
    !
    type(t_mesh), intent(inout), target    :: mesh
    character(len=MAX_PATH)                :: file_name
    integer                                :: ierror   ! MPI return error code
    integer                                :: k, n, fileID
    integer                                :: nchunk, chunk_size, ipos, iofs, mesh_check
    integer, allocatable, dimension(:)     :: mapping
    integer, allocatable, dimension(:)     :: ibuff
    real(kind=WP)                          :: t0, t1

!NR Cannot include the pointers before the targets are allocated...
!NR #include "associate_mesh.h"

    t0=MPI_Wtime()
    !___________________________________________________________________________
    allocate(mesh%nlevels(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
    allocate(mesh%nlevels_nod2D(myDim_nod2D+eDim_nod2D))
    
    !___________________________________________________________________________
    !mesh related files will be read in chunks of chunk_size
    chunk_size=100000
    ! Allocate mapping array (chunk_size), It will be used for several purposes 
    allocate(mapping(chunk_size))
    allocate(ibuff(chunk_size))
 
    !___________________________________________________________________________
    ! Part I: reading levels at elements...
    if (mype==0)  then 
        fileID=10
        file_name=trim(meshpath)//'elvls.out'
        open(fileID, file=file_name)
        write(*,*) 'reading '// trim(file_name)   
    end if
    
    ! 0 proc reads the data in chunks and distributes it between other procs
    mesh_check=0
    do nchunk=0, (mesh%elem2D-1)/chunk_size
        !_______________________________________________________________________
        !create the mapping for the current chunk
        mapping(1:chunk_size)=0
        do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
            ! myList_elem2D(n) contains global element index of the local
            ! element on that CPU
            ! ipos is integer, (myList_elem2D(n)-1)/chunk_size always rounds 
            ! off to integer values
            ! --> ipos is an index to which chunk a global element on a local CPU 
            !     belongs
            ipos=(myList_elem2D(n)-1)/chunk_size
            
            ! if global element chunk index (ipos) lies within the actual chunk
            if (ipos==nchunk) then
                iofs=myList_elem2D(n)-nchunk*chunk_size
                ! connect chunk reduced (iofs) global element index with local
                ! element index n --> mapping(iofs)=n
                mapping(iofs)=n
            end if
        end do
        
        !_______________________________________________________________________
        ! read the chunk piece into the buffer --> done only by one CPU (mype==0)
        ! k ... is actual chunk size, considers also possible change in chunk size
        !       at the end i.e elem2d=130000, nchunk_0 = 100000, nchunk_1=30000
        k=min(chunk_size, mesh%elem2D-nchunk*chunk_size)
        if (mype==0) then
            do n=1, k
                read(fileID,*) ibuff(n)
            end do
        end if
        
        !_______________________________________________________________________
        ! broadcast chunk buffer to all other CPUs (k...size of buffer)
        call MPI_BCast(ibuff(1:k), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
        
        !_______________________________________________________________________
        ! fill the local arrays
        do n=1, k      
            if (mapping(n)>0) then
                mesh_check=mesh_check+1
                mesh%nlevels(mapping(n))=ibuff(n)
            end if
        end do
    end do ! --> do nchunk=0, (mesh%elem2D-1)/chunk_size
    if (mype==0) close(fileID)
    if (mesh_check/=myDim_elem2D+eDim_elem2D+eXDim_elem2D) then
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: while reading elvls.out on mype=', mype
        write(*,*) '        ',mesh_check, ' values have been read in according to partitioning'
        write(*,*) '        it does not equal to myDim_elem2D+eDim_elem2D = ', myDim_elem2D+eDim_elem2D
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
    end if

    !___________________________________________________________________________
    ! Part II: reading levels at nodes...
    if (mype==0)  then 
        file_name=trim(meshpath)//'nlvls.out'
        open(fileID, file=file_name)
        write(*,*) 'reading '// trim(file_name)   
    end if
    
    ! 0 proc reads the data in chunks and distributes it between other procs
    mesh_check=0
    do nchunk=0, (mesh%nod2D-1)/chunk_size
        !_______________________________________________________________________
        !create the mapping for the current chunk
        mapping(1:chunk_size)=0
        do n=1, myDim_nod2D+eDim_nod2D
            ! myList_nod2D(n) contains global vertice index of the local
            ! vertice on that CPU
            ! ipos is integer, (myList_nod2D(n)-1)/chunk_size always rounds 
            ! off to integer values
            ! --> ipos is an index to which chunk a global vertice on a local CPU 
            !     belongs
            ipos=(myList_nod2D(n)-1)/chunk_size
            
            ! if global vertice chunk index (ipos) lies within the actual chunk
            if (ipos==nchunk) then
                iofs=myList_nod2D(n)-nchunk*chunk_size
                ! connect chunk reduced (iofs) global vertice index with local
                ! vertice index n --> mapping(iofs)=n
                mapping(iofs)=n
            end if
        end do
        
        !_______________________________________________________________________
        ! read the chunk piece into the buffer --> done only by one CPU (mype==0)
        ! k ... is actual chunk size, considers also possible change in chunk size
        !       at the end i.e nod2d=130000, nchunk_0 = 100000, nchunk_1=30000
        k=min(chunk_size, mesh%nod2D-nchunk*chunk_size)
        if (mype==0) then
            do n=1, k
                read(fileID,*) ibuff(n)
            end do
        end if
        
        !_______________________________________________________________________
        ! broadcast chunk buffer to all other CPUs (k...size of buffer)
        call MPI_BCast(ibuff(1:k), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
        
        !_______________________________________________________________________
        ! fill the local arrays
        do n=1, k      
            if (mapping(n)>0) then
                mesh_check=mesh_check+1
                mesh%nlevels_nod2D(mapping(n))=ibuff(n)
            end if
        end do
    end do
    if (mype==0) close(fileID)
    if (mesh_check/=myDim_nod2D+eDim_nod2D) then
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: while reading nlvls.out on mype=', mype
        write(*,*) '        ',mesh_check, ' values have been read in according to partitioning'
        write(*,*) '        it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
    end if
    
    !___________________________________________________________________________
    ! deallocate mapping and buffer array
    deallocate(ibuff)
    deallocate(mapping)

    !___________________________________________________________________________
    ! waits until all cpus have reached this points --> all cpus have to be
    ! supplied with cavity info 
    call MPI_BARRIER(MPI_COMM_FESOM, MPIerr)
    
    !___________________________________________________________________________
    ! initializes upper integration boundary index for all vertical vertice and 
    ! element loops, default = 1, but when cavity is used can be different 
    allocate(mesh%ulevels(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
    allocate(mesh%ulevels_nod2D(myDim_nod2D+eDim_nod2D))
    mesh%ulevels=1
    mesh%ulevels_nod2D=1
    
    !___________________________________________________________________________
    t1=MPI_Wtime()
    if (mype==0) then
        write(*,*) '____________________________________________________________________'
        write(*,*) ' --> 3D mesh was read in ', t1-t0, ' seconds'
        write(*,*) ' --> Min/max depth on mype : ', mype, -mesh%zbar(minval(mesh%nlevels)),-mesh%zbar(maxval(mesh%nlevels))
    endif
end subroutine find_levels
!
!
!_______________________________________________________________________________
! load cavity mesh files: cavity_depth, cavity_flag, cavity_nlvls.out and 
! cavity_elvls.out that are created during the partitioning when namelist.config flag
! use_cavity=.True.
!_______________________________________________________________________________
subroutine find_levels_cavity(mesh)
    use MOD_MESH
    use o_PARAM
    use g_PARSUP
    use g_config
    !
    implicit none
    !
    type(t_mesh), intent(inout), target :: mesh
    character(MAX_PATH)                 :: file_name
    integer                             :: ierror   ! MPI return error code
    integer                             :: k, n, fileID
    integer                             :: nchunk, chunk_size, ipos, iofs, mesh_check
    integer, allocatable, dimension(:)  :: mapping
    integer, allocatable, dimension(:)  :: ibuff
    real(kind=WP)                       :: t0, t1
    logical                             :: file_exist=.False.
    integer                             :: elem, elnodes(3), ule,  uln(3), node, j, nz
    integer, allocatable, dimension(:) :: numelemtonode
!NR Cannot include the pointers before the targets are allocated...
!NR #include "associate_mesh.h"
    
    t0=MPI_Wtime()
    !___________________________________________________________________________
    ! allocate arrays, reset pointers
!!PS     allocate(mesh%cavity_flag_e(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
!!PS     allocate(mesh%cavity_flag_n(myDim_nod2D+eDim_nod2D))
    allocate(mesh%cavity_depth(myDim_nod2D+eDim_nod2D))
    
    !___________________________________________________________________________
    ! mesh related files will be read in chunks of chunk_size
    chunk_size=100000
    
    ! Allocate mapping array (chunk_size) --> It will be used for several purposes 
    allocate(mapping(chunk_size))
    allocate(ibuff(chunk_size))
    
    !___________________________________________________________________________
    ! Part I: reading cavity levels at elements...
    if (mype==0)  then 
        fileID=10
        file_name=trim(meshpath)//'cavity_elvls.out'
        file_exist=.False.
        inquire(file=trim(file_name),exist=file_exist)
        if (file_exist) then
            open(fileID, file=file_name)
            write(*,*) ' --> open '// trim(file_name)   
        else
            write(*,*)
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: could not find file: cavity_elvls.out '
            write(*,*) '        Either wrong mesh path, or this file still needs to be'
            write(*,*) '        create during the partitioning with namelist.config flag'
            write(*,*) '        use_cavity=.true.'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            call par_ex
        end if 
    end if
    
    ! 0 proc reads the data in chunks and distributes it between other procs
    mesh_check=0
    do nchunk=0, (mesh%elem2D-1)/chunk_size
        !_______________________________________________________________________
        !create the mapping for the current chunk
        mapping(1:chunk_size)=0
        do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
            ! myList_elem2D(n) contains global element index of the local
            ! element on that CPU
            ! ipos is integer, (myList_elem2D(n)-1)/chunk_size always rounds 
            ! off to integer values
            ! --> ipos is an index to which chunk a global element on a local CPU 
            !     belongs
            ipos=(myList_elem2D(n)-1)/chunk_size
            
            ! if global element chunk index (ipos) lies within the actual chunk
            if (ipos==nchunk) then
                iofs=myList_elem2D(n)-nchunk*chunk_size
                ! connect chunk reduced (iofs) global element index with local
                ! element index n --> mapping(iofs)=n
                mapping(iofs)=n
            end if
        end do
        
        !_______________________________________________________________________
        ! read the chunk piece into the buffer --> done only by one CPU (mype==0)
        ! k ... is actual chunk size, considers also possible change in chunk size
        !       at the end i.e elem2d=130000, nchunk_0 = 100000, nchunk_1=30000
        k=min(chunk_size, mesh%elem2D-nchunk*chunk_size)
        if (mype==0) then
            do n=1, k
                read(fileID,*) ibuff(n)
            end do
        end if
        
        !_______________________________________________________________________
        ! broadcast chunk buffer to all other CPUs (k...size of buffer)
        call MPI_BCast(ibuff(1:k), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
        
        !_______________________________________________________________________
        ! fill the local arrays on each CPU
        do n=1, k
            ! if local mapping == 0, than the global chunk buffer element does not 
            ! belong to local cpu
            if (mapping(n)>0) then
                mesh_check=mesh_check+1
                mesh%ulevels(mapping(n))=ibuff(n)
            end if
        end do
        
    end do ! --> do nchunk=0, (mesh%elem2D-1)/chunk_size
    if (mype==0) close(fileID)
    if (mesh_check/=myDim_elem2D+eDim_elem2D+eXDim_elem2D) then
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: while reading cavity_elvls.out on mype=', mype
        write(*,*) '        ',mesh_check, ' values have been read in according to partitioning'
        write(*,*) '        it does not equal to myDim_elem2D+eDim_elem2D = ', myDim_elem2D+eDim_elem2D
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
        write(*,*)
    end if

    !___________________________________________________________________________
    ! Part II: reading cavity level at nodes
    if (mype==0)  then 
        file_name=trim(meshpath)//'cavity_nlvls.out'
        file_exist=.False.
        inquire(file=trim(file_name),exist=file_exist)
        if (file_exist) then
            open(fileID, file=file_name)
            write(*,*) ' --> open '// trim(file_name)   
        else
            write(*,*)
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: could not find file: cavity_nlvls.out '
            write(*,*) '        Either wrong mesh path, or this file still needs to be'
            write(*,*) '        create during the partitioning with namelist.config flag'
            write(*,*) '        use_cavity=.true.'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            call par_ex
        end if    
    end if
    
    ! 0 proc reads the data in chunks and distributes it between other procs
    mesh_check=0
    do nchunk=0, (mesh%nod2D-1)/chunk_size
        !_______________________________________________________________________
        !create the mapping for the current chunk
        mapping(1:chunk_size)=0
        do n=1, myDim_nod2D+eDim_nod2D
            ! myList_nod2D(n) contains global vertice index of the local
            ! vertice on that CPU
            ! ipos is integer, (myList_nod2D(n)-1)/chunk_size always rounds 
            ! off to integer values
            ! --> ipos is an index to which chunk a global vertice on a local CPU 
            !     belongs
            ipos=(myList_nod2D(n)-1)/chunk_size
            
            ! if global vertice chunk index (ipos) lies within the actual chunk
            if (ipos==nchunk) then
                iofs=myList_nod2D(n)-nchunk*chunk_size
                ! connect chunk reduced (iofs) global vertice index with local
                ! vertice index n --> mapping(iofs)=n
                mapping(iofs)=n
            end if
        end do
        
        !_______________________________________________________________________
        ! read the chunk piece into the buffer --> done only by one CPU (mype==0)
        ! k ... is actual chunk size, considers also possible change in chunk size
        !       at the end i.e nod2d=130000, nchunk_0 = 100000, nchunk_1=30000
        k=min(chunk_size, mesh%nod2D-nchunk*chunk_size)
        if (mype==0) then
            do n=1, k
                read(fileID,*) ibuff(n)
            end do
        end if
        
        !_______________________________________________________________________
        ! broadcast chunk buffer to all other CPUs (k...size of buffer)
        call MPI_BCast(ibuff(1:k), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
        
        !_______________________________________________________________________
        ! fill the local arrays
        do n=1, k      
            if (mapping(n)>0) then
                mesh_check=mesh_check+1
                mesh%ulevels_nod2D(mapping(n))=ibuff(n)
            end if
        end do
    end do ! --> do nchunk=0, (mesh%nod2D-1)/chunk_size
    if (mype==0) close(fileID)
    if (mesh_check/=myDim_nod2D+eDim_nod2D) then
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: while reading cavity_nlvls.out on mype=', mype
        write(*,*) '        ',mesh_check, ' values have been read in according to partitioning'
        write(*,*) '        it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
    end if
    
    !___________________________________________________________________________
    ! Part III: computing cavity flag at nodes and elements
!!PS     mesh%cavity_flag_e = 0
!!PS     do n=1,myDim_elem2D+eDim_elem2D+eXDim_elem2D
!!PS         if (mesh%ulevels(n)>1) mesh%cavity_flag_e(n)=1
!!PS     end do    
!!PS     mesh%cavity_flag_n = 0
!!PS     do n=1,myDim_nod2D+eDim_nod2D
!!PS         if (mesh%ulevels_nod2D(n)>1) mesh%cavity_flag_n(n)=1
!!PS     end do    
!!PS     if (mype==0)  then 
!!PS         file_name=trim(meshpath)//'cavity_flag.out'
!!PS         open(fileID, file=file_name)
!!PS         write(*,*) 'reading '// trim(file_name)   
!!PS     end if
!!PS     
!!PS     ! 0 proc reads the data in chunks and distributes it between other procs
!!PS     mesh_check=0
!!PS     do nchunk=0, (mesh%nod2D-1)/chunk_size
!!PS         !_______________________________________________________________________
!!PS         !create the mapping for the current chunk
!!PS         mapping(1:chunk_size)=0
!!PS         do n=1, myDim_nod2D+eDim_nod2D
!!PS             ! myList_nod2D(n) contains global vertice index of the local
!!PS             ! vertice on that CPU
!!PS             ! ipos is integer, (myList_nod2D(n)-1)/chunk_size always rounds 
!!PS             ! off to integer values
!!PS             ! --> ipos is an index to which chunk a global vertice on a local CPU 
!!PS             !     belongs
!!PS             ipos=(myList_nod2D(n)-1)/chunk_size
!!PS             
!!PS             ! if global vertice chunk index (ipos) lies within the actual chunk
!!PS             if (ipos==nchunk) then
!!PS                 iofs=myList_nod2D(n)-nchunk*chunk_size
!!PS                 ! connect chunk reduced (iofs) global vertice index with local
!!PS                 ! vertice index n --> mapping(iofs)=n
!!PS                 mapping(iofs)=n
!!PS             end if
!!PS         end do
!!PS         
!!PS         !_______________________________________________________________________
!!PS         ! read the chunk piece into the buffer --> done only by one CPU (mype==0)
!!PS         ! k ... is actual chunk size, considers also possible change in chunk size
!!PS         !       at the end i.e nod2d=130000, nchunk_0 = 100000, nchunk_1=30000
!!PS         k=min(chunk_size, mesh%nod2D-nchunk*chunk_size)
!!PS         if (mype==0) then
!!PS             do n=1, k
!!PS                 read(fileID,*) ibuff(n)
!!PS             end do
!!PS         end if
!!PS         
!!PS         !_______________________________________________________________________
!!PS         ! broadcast chunk buffer to all other CPUs (k...size of buffer)
!!PS         call MPI_BCast(ibuff(1:k), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
!!PS         
!!PS         !_______________________________________________________________________
!!PS         ! fill the local arrays
!!PS         do n=1, k      
!!PS             if (mapping(n)>0) then
!!PS                 mesh_check=mesh_check+1
!!PS                 mesh%cavity_flag(mapping(n))=ibuff(n)
!!PS             end if
!!PS         end do
!!PS     end do ! --> do nchunk=0, (mesh%nod2D-1)/chunk_size
!!PS     if (mype==0) close(fileID)
!!PS     if (mesh_check/=myDim_nod2D+eDim_nod2D) then
!!PS         write(*,*)
!!PS         print *, achar(27)//'[33m'
!!PS         write(*,*) '____________________________________________________________________'
!!PS         write(*,*) ' ERROR: while reading cavity_flag.out on mype=', mype
!!PS         write(*,*) '        ',mesh_check, ' values have been read in according to partitioning'
!!PS         write(*,*) '        it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
!!PS         write(*,*) '____________________________________________________________________'
!!PS         print *, achar(27)//'[0m'
!!PS     end if
    
    !___________________________________________________________________________
    ! Part IV: reading cavity depth at nodes
    if (mype==0)  then 
        file_name=trim(meshpath)//'cavity_depth.out'
        file_exist=.False.
        inquire(file=trim(file_name),exist=file_exist)
        if (file_exist) then
            open(fileID, file=file_name)
            write(*,*) ' --> open '// trim(file_name)   
        else
            write(*,*)
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: could not find file: cavity_depth.out '
            write(*,*) '        Wrong mesh path ? This file provides the necessary depth'
            write(*,*) '        information of the cavity.'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            call par_ex
        end if 
    end if
    
    ! 0 proc reads the data in chunks and distributes it between other procs
    mesh_check=0
    do nchunk=0, (mesh%nod2D-1)/chunk_size
        !_______________________________________________________________________
        !create the mapping for the current chunk
        mapping(1:chunk_size)=0
        do n=1, myDim_nod2D+eDim_nod2D
            ! myList_nod2D(n) contains global vertice index of the local
            ! vertice on that CPU
            ! ipos is integer, (myList_nod2D(n)-1)/chunk_size always rounds 
            ! off to integer values
            ! --> ipos is an index to which chunk a global vertice on a local CPU 
            !     belongs
            ipos=(myList_nod2D(n)-1)/chunk_size
            
            ! if global vertice chunk index (ipos) lies within the actual chunk
            if (ipos==nchunk) then
                iofs=myList_nod2D(n)-nchunk*chunk_size
                ! connect chunk reduced (iofs) global vertice index with local
                ! vertice index n --> mapping(iofs)=n
                mapping(iofs)=n
            end if
        end do
        
        !_______________________________________________________________________
        ! read the chunk piece into the buffer --> done only by one CPU (mype==0)
        ! k ... is actual chunk size, considers also possible change in chunk size
        !       at the end i.e nod2d=130000, nchunk_0 = 100000, nchunk_1=30000
        k=min(chunk_size, mesh%nod2D-nchunk*chunk_size)
        if (mype==0) then
            do n=1, k
                read(fileID,*) ibuff(n)
            end do
        end if
        
        !_______________________________________________________________________
        ! broadcast chunk buffer to all other CPUs (k...size of buffer)
        call MPI_BCast(ibuff(1:k), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
        
        !_______________________________________________________________________
        ! fill the local arrays
        do n=1, k      
            if (mapping(n)>0) then
                mesh_check=mesh_check+1
                mesh%cavity_depth(mapping(n))=ibuff(n)
            end if
        end do
    end do ! --> do nchunk=0, (mesh%nod2D-1)/chunk_size
    if (mype==0) close(fileID)
    if (mesh_check/=myDim_nod2D+eDim_nod2D) then
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: while reading cavity_depth.out on mype=', mype
        write(*,*) '        ',mesh_check, ' values have been read in according to partitioning'
        write(*,*) '        it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
    end if
    
    !___________________________________________________________________________
    ! deallocate mapping and buffer array
    deallocate(ibuff)
    deallocate(mapping)

    !___________________________________________________________________________
    ! waits until all cpus have reached this points --> all cpus have to be
    ! supplied with cavity info 
    call MPI_BARRIER(MPI_COMM_FESOM, MPIerr)
    
    !___________________________________________________________________________
    t1=MPI_Wtime()
    if (mype==0) then
        write(*,*) '____________________________________________________________________'
        write(*,*) ' --> cavity info read in ', t1-t0, ' seconds'
    end if
    
    !___________________________________________________________________________
    ! check cavity info 
    do elem=1,myDim_elem2d
        elnodes = mesh%elem2D_nodes(:,elem)
        ule = mesh%ulevels(elem)
        uln = mesh%ulevels_nod2D(elnodes)
        if (ule < maxval(uln)) then 
            write(*,*) ' --> found cavity elem depth shallower than valid cavity node depth, mype=', mype
        end if 
    end do 
    
    
    !___________________________________________________________________________
    allocate(numelemtonode(mesh%nl))
    do node=1, myDim_nod2D+eDim_nod2D
        numelemtonode=0
        !_______________________________________________________________________
        do j=1,mesh%nod_in_elem2D_num(node)
            elem=mesh%nod_in_elem2D(j,node)
            do nz=mesh%ulevels(elem),mesh%nlevels(elem)-1
                numelemtonode(nz) = numelemtonode(nz) + 1
            end do
        end do
        
        !_______________________________________________________________________
        ! check how many triangle elements contribute to every vertice in every layer
        ! every vertice in every layer should be connected to at least two triangle 
        ! elements !
        do nz=mesh%ulevels_nod2D(node),mesh%nlevels_nod2D(node)-1
            if (numelemtonode(nz)== 1) then 
                write(*,*) 'ERROR A: found vertice with just one triangle:', mype, node, nz
            end if 
        end do 
        
    end do
    deallocate(numelemtonode)
    
end subroutine find_levels_cavity
!
!
!_______________________________________________________________________________
! load cavity mesh files: cavity_depth, cavity_flag, cavity_nlvls.out and 
! cavity_elvls.out that are created during the partitioning when namelist.config flag
! use_cavity=.True.
!_______________________________________________________________________________
subroutine find_levels_min_e2n(mesh)
    use MOD_MESH
    use o_PARAM
    use g_PARSUP
    use g_config
    use g_comm_auto
    !
    implicit none
    !
    type(t_mesh), intent(inout), target :: mesh
    integer                             :: node, k
    real(kind=WP)                       :: t0, t1
    
!NR Cannot include the pointers before the targets are allocated...
!NR #include "associate_mesh.h"
    
    t0=MPI_Wtime()
    !___________________________________________________________________________
    allocate(mesh%nlevels_nod2D_min(myDim_nod2D+eDim_nod2D))
    allocate(mesh%ulevels_nod2D_max(myDim_nod2D+eDim_nod2D))
    do node=1, myDim_nod2d
        k=mesh%nod_in_elem2D_num(node)
        ! minimum depth in neigbouring elements around node n
        mesh%nlevels_nod2D_min(node)=minval(mesh%nlevels(mesh%nod_in_elem2D(1:k,node)))
        mesh%ulevels_nod2D_max(node)=maxval(mesh%ulevels(mesh%nod_in_elem2D(1:k,node)))
    end do
    call exchange_nod(mesh%nlevels_nod2D_min)
    call exchange_nod(mesh%ulevels_nod2D_max)
    
    !___________________________________________________________________________
    t1=MPI_Wtime()
    if (mype==0) then
        write(*,*) '____________________________________________________________________'
        write(*,*) ' --> find min/max level e2n in', t1-t0, ' seconds'
    end if
    
end subroutine find_levels_min_e2n
!
!
!
!===========================================================================
SUBROUTINE test_tri(mesh)
USE MOD_MESH
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_rotate_grid
IMPLICIT NONE
! Check the order of nodes in triangles; correct it if necessary to make
! it same sense (clockwise) 
type(t_mesh), intent(inout), target :: mesh
real(kind=WP)               ::  a(2), b(2), c(2),  r
integer                     ::  n, nx, elnodes(3)
real(kind=WP)               :: t0, t1

   t0=MPI_Wtime()
   
   DO n=1, myDim_elem2D
      elnodes=mesh%elem2D_nodes(:,n)
	  
          a=mesh%coord_nod2D(:,elnodes(1))
	  b=mesh%coord_nod2D(:,elnodes(2))-a
	  c=mesh%coord_nod2D(:,elnodes(3))-a
          
	  call trim_cyclic(b(1))
	  call trim_cyclic(c(1))
	  
	    
	  r=b(1)*c(2)-b(2)*c(1)
	  if (r>0._WP) then
	  ! Vector b is to right of c
	  ! Exchange second and third nodes:
	  
	  nx=elnodes(2)
	  elnodes(2)=elnodes(3)
	  elnodes(3)=nx
	  mesh%elem2D_nodes(:,n)=elnodes
      end if
   END DO
   t1=MPI_Wtime()
   if (mype==0) then
      write(*,*) 'test_tri finished in ', t1-t0, ' seconds'
      write(*,*) '========================='
   endif

END SUBROUTINE  test_tri
!=========================================================================
SUBROUTINE load_edges(mesh)
USE MOD_MESH
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
type(t_mesh), intent(inout), target   :: mesh
character(MAX_PATH)                   :: file_name
integer                               :: counter, n, m, nn, k, q, fileID
integer                               :: elems(2), elem
integer                               :: elnodes(3), ed(2), eledges(3)
integer, allocatable                  :: aux(:)         
real(kind=WP)                         :: t0, t1
integer                               :: nchunk, chunk_size, ipos, iofs, mesh_check
integer, allocatable, dimension(:)    :: mapping
integer, allocatable, dimension(:,:)  :: ibuff
integer                               :: ierror              ! return error code

!NR Cannot include the pointers before the targets are allocated...
!NR #include "associate_mesh.h"

t0=MPI_Wtime()

!==============================
! Edge array is already available (we computed it in the init phase)
! (a) Read list of edges and tri containing them from file 
!==============================
!mesh related files will be read in chunks of chunk_size
chunk_size=100000
!==============================
! Allocate mapping array (chunk_size)
! It will be used for several purposes 
!==============================
 allocate(mapping(chunk_size))
 allocate(ibuff(4,chunk_size))
! 0 proc reads the data and sends it to other procs
 fileID=10
 if (mype==0) then
    file_name=trim(meshpath)//'edgenum.out'
    write(*,*) 'reading '// trim(file_name)
    open(fileID, file=trim(file_name))
    read(fileID,*) mesh%edge2D
    read(fileID,*) mesh%edge2D_in
    write(*,*) '2D mesh info : edge2D=', mesh%edge2D
    close(fileID)
  end if
  call MPI_BCast(mesh%edge2D,    1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(mesh%edge2D_in, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

  allocate(mesh%edges(2,myDim_edge2D+eDim_edge2D))
  allocate(mesh%edge_tri(2,myDim_edge2D+eDim_edge2D))

  ! 0 proc reads the data in chunks and distributes it between other procs
  if (mype==0) then
     file_name=trim(meshpath)//'edges.out'
     open(fileID,   file=trim(file_name))
     write(*,*) 'reading '// trim(file_name)

     file_name=trim(meshpath)//'edge_tri.out'
     open(fileID+1, file=trim(file_name))
     write(*,*) 'reading '// trim(file_name)
  end if

 mesh_check=0
 do nchunk=0, (mesh%edge2D-1)/chunk_size
    !create the mapping for the current chunk
    mapping(1:chunk_size)=0
    do n=1, myDim_edge2D+eDim_edge2D
       ipos=(myList_edge2D(n)-1)/chunk_size
       if (ipos==nchunk) then
          iofs=myList_edge2D(n)-nchunk*chunk_size
          mapping(iofs)=n
       end if
    end do
    !read the chunk into the buffers
    k=min(chunk_size, mesh%edge2D-nchunk*chunk_size)
    if (mype==0) then
       do n=1, k
          read(fileID  ,*) ibuff(1:2,n) !edge nodes
          read(fileID+1,*) ibuff(3:4,n) !edge elements
       end do
    end if
    call MPI_BCast(ibuff(1:4,1:k), 4*k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
!    call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
!    call MPI_BCast(ibuff(1:k,3), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
!    call MPI_BCast(ibuff(1:k,4), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
    ! fill the local arrays
    do n=1, k      
       if (mapping(n)>0) then
          mesh_check=mesh_check+1
          mesh%edges   (:, mapping(n))=ibuff(1:2,n)
          mesh%edge_tri(:, mapping(n))=ibuff(3:4,n)
       end if
    end do
 end do

 if (mesh_check/=myDim_edge2D+eDim_edge2D) then
    write(*,*) 'ERROR while reading edges.out/edge_tri.out on mype=', mype
    write(*,*) mesh_check, ' values have been read in according to partitioning'
    write(*,*) 'it does not equal to myDim_edge2D+eDim_edge2D = ', myDim_edge2D+eDim_edge2D
 end if

 if (mype==0) then 
    close(fileID)
    close(fileID+1)
 end if 
 
! =========
! edges are in global numbering, transform it to local indexing
! =========
  mesh_check=0
  do nchunk=0, (mesh%nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_edge2D+eDim_edge2D
        do m=1, 2
           nn=mesh%edges(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk) then
              mesh_check=mesh_check+1
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              mesh%edges(m,n)=-mapping(iofs) 
           end if
        end do
     end do
  end do
  mesh%edges=-mesh%edges
  mesh_check=mesh_check/2
  if (mesh_check/=myDim_edge2D+eDim_edge2D) then
     write(*,*) 'ERROR while transforming edge nodes to local indexing on mype=', mype
     write(*,*) mesh_check, ' edges have been transformed!'
     write(*,*) 'It does not equal to myDim_edge2D+eDim_edge2D = ', myDim_edge2D+eDim_edge2D
  end if  
! =========
! edge_tri are in global numbering, transform it to local indexing
! =========

!*** Vadym pointed to the bug with bilding edge_tri
! if not doing so the land boundary elements will be set to 999 instead of being zeros or negative!
  where (mesh%edge_tri<0)
         mesh%edge_tri=0
  end where
!***
  mesh_check=0
  do nchunk=0, (mesh%elem2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
        ipos=(myList_elem2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_edge2D+eDim_edge2D
        do m=1, 2
           nn=mesh%edge_tri(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk .and. nn > 0) then
              mesh_check=mesh_check+abs(m-2) !only first triangle will contribute to statistic
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              mesh%edge_tri(m,n)=-mapping(iofs) 
           end if
        end do
     end do
  end do
  mesh%edge_tri=-mesh%edge_tri
  if (mesh_check/=myDim_edge2D+eDim_edge2D) then
     write(*,*) 'ERROR while transforming edge elements to local indexing on mype=', mype
     write(*,*) mesh_check, ' edges have been transformed!'
     write(*,*) 'It does not equal to myDim_edge2D+eDim_edge2D = ', myDim_edge2D+eDim_edge2D
  end if  
! =========


! Now the only way to check whether an edge is on boundary is 
! through myList_edge2D(n):  myList_edge2D(n)>edge2D_in == boundary edge

! (b) We need an array inverse to edge_tri listing edges
! of a given triangle 
allocate(mesh%elem_edges(3,myDim_elem2D))
allocate(aux(myDim_elem2D))
aux=0
DO n=1, myDim_edge2D+eDim_edge2D
   DO k=1,2
      q=mesh%edge_tri(k,n)   ! triangle number
	  if((q>0).and.(q<=myDim_elem2D)) then
	  aux(q)=aux(q)+1
	  mesh%elem_edges(aux(q),q)=n
	  end if
   END DO
END DO
deallocate(aux)
! The edges in this list should be ordered so that they
! are listed in the same rotation sense as nodes.
DO elem=1,myDim_elem2D
   elnodes=mesh%elem2D_nodes(:,elem)
   eledges=mesh%elem_edges(:,elem)
   DO q=1,3
      DO k=1,3
         if((mesh%edges(1,eledges(k)).ne.elnodes(q)).and. &
            (mesh%edges(2,eledges(k)).ne.elnodes(q))) then
           mesh%elem_edges(q,elem)=eledges(k)
	   exit
         end if
      END DO
   END DO
END DO
! The edge and elem lists agree in the sense that edge1 does not
! contain node 1 and so on

deallocate(ibuff)
deallocate(mapping)

t1=MPI_Wtime()
if (mype==0) then
   write(*,*) 'load_edges finished in ', t1-t0, ' seconds'
   write(*,*) '========================='
endif

END SUBROUTINE load_edges
!===========================================================================
SUBROUTINE find_neighbors(mesh)
! For each element three its element neighbors are found
! For each node the elements containing it are found
! Allocated are:
! elem_neighbors(3,myDim_elem2D)
! nod_in_elem2D_num(myDim_nod2D)
! nod_in_elem2D(:, myDim_nod2D)
! 

USE o_PARAM
USE MOD_MESH
USE g_PARSUP
USE g_ROTATE_grid
use g_comm_auto
use elem_center_interface
implicit none
type(t_mesh), intent(inout), target :: mesh
integer                     :: elem, eledges(3), elem1, j, n, node, enum,elems(3),count1,count2,exit_flag,i,nz
integer, allocatable        :: temp_i(:)
integer                     :: mymax(npes), rmax(npes)
real(kind=WP)               :: gx,gy,rx,ry
real(kind=WP)               :: t0, t1

!NR Cannot include the pointers before the targets are allocated...
!NR #include "associate_mesh.h"

CALL MPI_BARRIER(MPI_COMM_FESOM, MPIerr)
t0=MPI_Wtime()

 ! =============
 ! elem neighbors == those that share edges
 ! =============
   allocate(mesh%elem_neighbors(3,myDim_elem2D))
   mesh%elem_neighbors=0
  
DO elem=1,myDim_elem2D
   eledges=mesh%elem_edges(:,elem)
   DO j=1,3
   elem1=mesh%edge_tri(1,eledges(j))
   if(elem1==elem) elem1=mesh%edge_tri(2,eledges(j))
   mesh%elem_neighbors(j,elem)=elem1
   END DO
END DO
 ! =============
 ! Node neighbourhood
 ! == elements that contain node n
 ! We need eDim neighborhood too for MUSCL advection. 
 ! And we already have the place allocated for all 
 ! these neighbor elements: it is eDim_elem2D+eXDim_elem2D
 ! =============	 
 allocate(mesh%nod_in_elem2D_num(myDim_nod2D+eDim_nod2D)) 
 mesh%nod_in_elem2D_num=0
 do n=1,myDim_elem2D
    do j=1,3
    node=mesh%elem2D_nodes(j,n)
    if (node>myDim_nod2D) cycle
    mesh%nod_in_elem2D_num(node)=mesh%nod_in_elem2D_num(node)+1
    end do
 end do
CALL MPI_BARRIER(MPI_COMM_FESOM, MPIerr)

 mymax=0
 rmax=0
 mymax(mype+1)=maxval(mesh%nod_in_elem2D_num(1:myDim_nod2D))
 call MPI_AllREDUCE( mymax, rmax, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

 allocate(mesh%nod_in_elem2D(maxval(rmax),myDim_nod2D+eDim_nod2D))
 mesh%nod_in_elem2D=0
 mesh%nod_in_elem2D_num=0
 do n=1,myDim_elem2D   
    do j=1,3
    node=mesh%elem2D_nodes(j,n)
    if (node>myDim_nod2D) cycle 
    mesh%nod_in_elem2D_num(node)=mesh%nod_in_elem2D_num(node)+1
    mesh%nod_in_elem2D(mesh%nod_in_elem2D_num(node),node)=n
    end do
 end do

 call exchange_nod(mesh%nod_in_elem2D_num)
 allocate (temp_i(myDim_nod2D+eDim_nod2D))
 temp_i=0
 DO n=1, maxval(rmax)
       ! Exchange global element numbers
       do j=1,myDim_nod2D
         if (mesh%nod_in_elem2D(n,j)>0) temp_i(j)=myList_elem2D(mesh%nod_in_elem2D(n,j))
       enddo
       call exchange_nod(temp_i)
       mesh%nod_in_elem2D(n,:)=temp_i
 END DO
 deallocate(temp_i)
 ! Substitute back local element numbers
 allocate(temp_i(mesh%elem2D))
 temp_i=0
 Do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
    temp_i(myList_elem2D(n))=n
 END DO
 DO n=1, myDim_nod2D+eDim_nod2D     
    DO j=1, mesh%nod_in_elem2D_num(n)
       mesh%nod_in_elem2D(j,n)=temp_i(mesh%nod_in_elem2D(j,n))
    END DO
 END DO
 deallocate(temp_i)

 ! Among elem_neighbors there can be negative numbers. These correspond to 
 ! boundary elements for which neighbours are absent. However, an element 
 ! should have at least two valid neighbors
 ! ============
 ! Test that there are at least two neighbors at the surface:
 ! ============ 
DO elem=1,myDim_elem2D
   elem1=0
   DO j=1,3
   if(mesh%elem_neighbors(j,elem)>0) elem1=elem1+1
   END DO
   if (elem1<2) then
    write(*,*) 'Insufficient number of neighbors ', myList_elem2D(elem)
    call par_ex(1)
    STOP
   end if
END DO


#if defined (__oasis)
    allocate(mesh%x_corners(myDim_nod2D, maxval(rmax)))
    allocate(mesh%y_corners(myDim_nod2D, maxval(rmax)))
    DO n=1, myDim_nod2D
       DO j=1, mesh%nod_in_elem2D_num(n)
           elem=mesh%nod_in_elem2D(j,n)
           call elem_center(elem, rx, ry, mesh)
           call r2g(gx, gy, rx, ry)
           mesh%x_corners(n, j)=gx/rad
           mesh%y_corners(n, j)=gy/rad
       END DO
       mesh%x_corners(n, mesh%nod_in_elem2D_num(n)+1:maxval(rmax))=mesh%x_corners(n,mesh%nod_in_elem2D_num(n)) !or -999?
       mesh%y_corners(n, mesh%nod_in_elem2D_num(n)+1:maxval(rmax))=mesh%y_corners(n,mesh%nod_in_elem2D_num(n)) !or -999?
       ! to get the max number of corners use size(x_corners, 2)
    END DO
#endif

t1=MPI_Wtime()
if (mype==0) then
   write(*,*) 'find_neighbors finished in ', t1-t0, ' seconds'
   write(*,*) '========================='
endif

END SUBROUTINE find_neighbors
!==========================================================================
subroutine edge_center(n1, n2, x, y, mesh)

USE MOD_MESH
USE o_PARAM
USE g_CONFIG 
implicit none
integer                     :: n1, n2   ! nodes of the edge
real(kind=WP)               :: x, y, a(2), b(2)
type(t_mesh), intent(inout), target :: mesh

a=mesh%coord_nod2D(:,n1)
b=mesh%coord_nod2D(:,n2)
if(a(1)-b(1)> cyclic_length/2.0_WP) a(1)=a(1)-cyclic_length
if(a(1)-b(1)<-cyclic_length/2.0_WP) b(1)=b(1)-cyclic_length
x=0.5_WP*(a(1)+b(1))
y=0.5_WP*(a(2)+b(2))

end subroutine edge_center
!==========================================================================
subroutine elem_center(elem, x, y, mesh)
USE MOD_MESH
USE o_PARAM
USE g_CONFIG  
implicit none
integer      :: elem, elnodes(3), k    
real(kind=WP) :: x, y, ax(3), amin

type(t_mesh), intent(inout), target :: mesh

   elnodes=mesh%elem2D_nodes(:,elem)
   ax=mesh%coord_nod2D(1, elnodes)
   amin=minval(ax)
   DO k=1,3
   if(ax(k)-amin>=cyclic_length/2.0_WP) ax(k)=ax(k)-cyclic_length
   if(ax(k)-amin<-cyclic_length/2.0_WP) ax(k)=ax(k)+cyclic_length
   END DO
   x=sum(ax)/3.0_WP
   y=sum(mesh%coord_nod2D(2,elnodes))/3.0_WP

end subroutine elem_center
!==========================================================================
SUBROUTINE mesh_areas(mesh)
    USE MOD_MESH
    USE o_PARAM
    USE o_arrays, only: dum_3d_n
    USE g_PARSUP
    USE g_ROTATE_GRID
    use g_comm_auto
    IMPLICIT NONE
    ! Collects auxilliary information on the mesh
    ! Allocated and filled in are:
    ! elem_area(myDim_elem2D)
    ! area(nl, myDim_nod2D)

    integer                                   :: n,j,q, elnodes(3), ed(2), elem, nz,nzmin, nzmax
    real(kind=WP)                             :: a(2), b(2), ax, ay, lon, lat, vol, vol2
    real(kind=WP), allocatable,dimension(:)   :: work_array
    integer, allocatable,dimension(:,:)       :: cavity_contribut
    real(kind=WP)                             :: t0, t1
    type(t_mesh), intent(inout), target       :: mesh

    !NR Cannot include the pointers before the targets are allocated...
    !NR #include "associate_mesh.h"

    t0=MPI_Wtime()
    
    ! area of triangles 
    allocate(mesh%elem_area(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
    
    ! area of upper edge and lower edge of scalar cell: size nl x node
    allocate(mesh%area(mesh%nl,myDim_nod2d+eDim_nod2D))
    
    ! "mid" area of scalar cell in case of cavity area \= areasvol, size: nl-1 x node
    allocate(mesh%areasvol(mesh%nl,myDim_nod2d+eDim_nod2D))
    
    ! area inverse
    allocate(mesh%area_inv(mesh%nl,myDim_nod2d+eDim_nod2D))
    allocate(mesh%areasvol_inv(mesh%nl,myDim_nod2d+eDim_nod2D))
    
    ! resolution at nodes 
    allocate(mesh%mesh_resolution(myDim_nod2d+eDim_nod2D))
    
    !___compute triangle areas__________________________________________________
    do n=1, myDim_elem2D
        elnodes=mesh%elem2D_nodes(:,n)
        ay=sum(mesh%coord_nod2D(2,elnodes))/3.0_WP
        ay=cos(ay)
        if (cartesian) ay=1.0_WP
        a = mesh%coord_nod2D(:,elnodes(2))-mesh%coord_nod2D(:,elnodes(1))
        b = mesh%coord_nod2D(:,elnodes(3))-mesh%coord_nod2D(:,elnodes(1))
        call trim_cyclic(a(1))
        call trim_cyclic(b(1))
        a(1)=a(1)*ay
        b(1)=b(1)*ay
        mesh%elem_area(n)=0.5_WP*abs(a(1)*b(2)-b(1)*a(2))
    end do
    call exchange_elem(mesh%elem_area)
    
    !___compute areas of upper/lower scalar cell edge___________________________
    ! areas at different levels (there can be partly land)
    ! --> only areas through which there is exchange are counted
    !
    !-----------------------------~+~~~~~~~+~~~
    ! ############################ |       |   
    ! ############################ |       |   layer k-3
    ! #################### ._______|_______|___area_k-2   
    ! ##  CAVITY  ######## | / / / |       |   
    ! #################### |/ /°/ /|       |   layer k-2 --> Transport:  w_k-2*A_k-1
    ! ############ ._______|_/_/_/_|_______|___area_k-1         -> A_k-1 lower prisma area defines 
    ! ############ |       |       |       |                    scalar area under the cavity
    ! ############ |   °   |       |       |   layer k-1
    !______________|_______|_______|_______|___area_k
    !      |       | / / / |       |       |   
    !      |       |/ /°/ /|       |       |   layer k --> Transport: w_k*A_k
    !______|_______|_/_/_/_|_______|_______|___area_k+1       -> A_k upper prisma face area defines      
    !      |       |       |       |       |                  scalar area of cell   
    !      |       |   °   |       |       |   layer k+1
    !______|_______|_______|_______|_______|___area_k+2
    ! #############|       |       |       |   
    ! #############|   °   |       |       |   layer k+2
    ! #############|_______|_______|_______|___area_k+3
    ! #####################|       |       |   
    ! #####################|       |       |   layer k+3
    ! ##  BOTTOM  #########|_______|_______|___area_k+4
    ! #############################|       |   
    ! #############################|       |   :
    ! #############################|_______|___area_k+5
    ! #########################################
    if (use_cavity) then
        allocate(cavity_contribut(mesh%nl,myDim_nod2d+eDim_nod2D))
        cavity_contribut = 0
    end if 
    
    mesh%area     = 0.0_WP
    do n=1, myDim_nod2D+eDim_nod2D
        do j=1,mesh%nod_in_elem2D_num(n)
            elem=mesh%nod_in_elem2D(j,n)
            
            !___________________________________________________________________
            ! compute scalar area of prisms at different depth layers. In normal 
            ! case without cavity the area of the scalar cell corresponds to the
            ! area of the upper edge of the prism --> if there is cavity its
            ! different. Directly under the cavity the area of scalar cell 
            ! corresponds to the area of the lower edge
            nzmin = mesh%ulevels(elem)
            nzmax = mesh%nlevels(elem)-1
            do nz=nzmin,nzmax
                mesh%area(nz,n)=mesh%area(nz,n)+mesh%elem_area(elem)/3.0_WP
            end do
            
            !___________________________________________________________________
            ! how many ocean-cavity triangles contribute to an upper edge of a 
            ! scalar area 
            if (use_cavity) then
                do nz=1,nzmin-1
                    cavity_contribut(nz,n)=cavity_contribut(nz,n)+1
                end do
            end if 
        end do
    end do
    
    !___compute "mid" scalar cell area__________________________________________
    ! for cavity case: redefine "mid" scalar cell area from upper edge of prism to 
    ! lower edge of prism if a cavity triangle is present at the upper scalar
    ! cell edge 
    mesh%areasvol = 0.0_WP
    if (use_cavity) then
        do n = 1, myDim_nod2D+eDim_nod2D
            nzmin = mesh%ulevels_nod2d(n)
            nzmax = mesh%nlevels_nod2d(n)-1
            do nz=nzmin,nzmax
                if (cavity_contribut(nz,n)>0) then
                    mesh%areasvol(nz,n) = mesh%area(min(nz+1,nzmax),n)
                else
                    mesh%areasvol(nz,n) = mesh%area(nz,n)
                end if
            end do 
        end do
        deallocate(cavity_contribut)
    ! for non cavity case: the "mid" area of the scalar cell always corresponds to 
    ! the area of the upper scalar cell edge    
    else    
        do n = 1, myDim_nod2D+eDim_nod2D
            nzmin = mesh%ulevels_nod2d(n)
            nzmax = mesh%nlevels_nod2d(n)-1
            do nz=nzmin,nzmax
                mesh%areasvol(nz,n) = mesh%area(nz,n)
            end do 
        end do
    end if 
    
    ! update to proper dimension
    ! coordinates are in radians, edge_dxdy are in meters,
    ! and areas are in m^2
    mesh%elem_area = mesh%elem_area*r_earth*r_earth
    mesh%area      = mesh%area     *r_earth*r_earth
    mesh%areasvol  = mesh%areasvol *r_earth*r_earth
 
    call exchange_nod(mesh%area)
    call exchange_nod(mesh%areasvol)
    
    !___compute inverse area____________________________________________________
    mesh%area_inv = 0.0_WP
    do n=1,myDim_nod2d+eDim_nod2D
        nzmin = mesh%ulevels_nod2d(n)
        nzmax = mesh%nlevels_nod2d(n)
        do nz=nzmin,nzmax
!!PS             mesh%area_inv(nz,n) = 1._WP/mesh%area(nz,n)
            if (mesh%area(nz,n) > 0._WP) then
                mesh%area_inv(nz,n) = 1._WP/mesh%area(nz,n)
            else
                mesh%area_inv(nz,n) = 0._WP
            end if
        end do
    end do
    
    if (use_cavity) then
        mesh%areasvol_inv = 0.0_WP
        do n=1,myDim_nod2d+eDim_nod2D
            nzmin = mesh%ulevels_nod2d(n)
            nzmax = mesh%nlevels_nod2d(n)-1
            do nz=nzmin,nzmax
!!PS                 mesh%areasvol_inv(nz,n) = 1._WP/mesh%areasvol(nz,n)
                if (mesh%areasvol(nz,n) > 0._WP) then
                    mesh%areasvol_inv(nz,n) = 1._WP/mesh%areasvol(nz,n)
                else
                    mesh%areasvol_inv(nz,n) = 0._WP
                end if
            end do
        end do
    else
        mesh%areasvol_inv = mesh%area_inv
    endif 
 
    !___compute scalar cell resolution__________________________________________
    allocate(work_array(myDim_nod2D))
    !!PS mesh%mesh_resolution=sqrt(mesh%area(1, :)/pi)*2._WP
    do n=1,myDim_nod2d+eDim_nod2D
        mesh%mesh_resolution(n)=sqrt(mesh%areasvol(mesh%ulevels_nod2d(n),n)/pi)*2._WP
    end do 
 
    ! smooth resolution
    do q=1, 3 !apply mass matrix N times to smooth the field
        do n=1, myDim_nod2D
            vol=0._WP
            work_array(n)=0._WP
            do j=1, mesh%nod_in_elem2D_num(n)
                elem=mesh%nod_in_elem2D(j, n)
                elnodes=mesh%elem2D_nodes(:,elem)
                work_array(n)=work_array(n)+sum(mesh%mesh_resolution(elnodes))/3._WP*mesh%elem_area(elem)
                vol=vol+mesh%elem_area(elem)
            end do
            work_array(n)=work_array(n)/vol
        end do
        do n=1,myDim_nod2D
            mesh%mesh_resolution(n)=work_array(n)
        end do
        call exchange_nod(mesh%mesh_resolution)
    end do
    deallocate(work_array)

    !___compute total ocean areas with/without cavity___________________________
    vol = 0.0_WP
    vol2= 0.0_WP
    do n=1, myDim_nod2D
        vol2=vol2+mesh%area(mesh%ulevels_nod2D(n), n) ! area also under cavity
        if (mesh%ulevels_nod2D(n)>1) cycle
        vol=vol+mesh%area(1, n) ! area only surface  
    end do
    mesh%ocean_area=0.0
    mesh%ocean_areawithcav=0.0
    call MPI_AllREDUCE(vol, mesh%ocean_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
        MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(vol2, mesh%ocean_areawithcav, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
        MPI_COMM_FESOM, MPIerr)
    
    !___write mesh statistics___________________________________________________
    if (mype==0) then
        write(*,*) '____________________________________________________________________'
        write(*,*) ' --> mesh statistics:', mype
        write(*,*)  mype, 'maxArea ',maxval(mesh%elem_area), '   MinArea ', minval(mesh%elem_area)
        write(*,*)  mype, 'maxScArea ',maxval(mesh%area(1,:)), &
                    '   MinScArea ', minval(mesh%area(1,:))
        write(*,*)  mype, 'Edges:    ', mesh%edge2D, ' internal ', mesh%edge2D_in
        if (mype==0) then
            write(*,*) '     > Total ocean surface area is           : ', mesh%ocean_area, ' m^2'
            write(*,*) '     > Total ocean surface area wth cavity is: ', mesh%ocean_areawithcav, ' m^2'
        end if
    endif

    t1=MPI_Wtime()
    if (mype==0) then
        write(*,*) '     > mesh_areas finished in ', t1-t0, ' seconds'
    endif
END SUBROUTINE mesh_areas

!===================================================================

SUBROUTINE mesh_auxiliary_arrays(mesh)
! Collects auxiliary information needed to speed up computations 
! of gradients, div. This also makes implementation of cyclicity 
! much more straightforward
! Allocated and filled in are:
! edge_dxdy(2,myDim_edge2D+eDim_edge2D)
! edge_cross_dxdy(4,myDim_edge2D)
! gradient_node(6,myDim_elem2D)
! gradient_elem(6,myDim_elem2d)
! metric_factor(myDim_elem2D+eDim_elem2D)
! elem_cos(myDim_elem2D+eDim_elem2D)
! coriolis(myDim_elem2D)

USE MOD_MESH
USE o_PARAM
USE i_PARAM
USE g_PARSUP
USE o_ARRAYS
USE g_ROTATE_grid
use g_comm_auto
use elem_center_interface
use edge_center_interface
IMPLICIT NONE

integer              :: n,j,q, elnodes(3), ed(2), elem, el(2), elnodes_(3),node
real(kind=WP)	     :: a(2), b(2), ax, ay, dfactor, lon, lat
real(kind=WP)	     :: deltaX31, deltaX21, deltaY31, deltaY21
real(kind=WP)        :: x(3), y(3), cxx, cxy, cyy, d
real(kind=WP), allocatable :: center_x(:), center_y(:), temp(:) 
real(kind=WP)              :: t0, t1
integer                    :: i, nn, ns
type(t_mesh), intent(inout), target :: mesh

!NR Cannot include the pointers before the targets are allocated...
!NR #include "associate_mesh.h"
t0=MPI_Wtime()

 allocate(mesh%edge_dxdy(2,myDim_edge2D+eDim_edge2D))
 allocate(mesh%edge_cross_dxdy(4,myDim_edge2D+eDim_edge2D))
 allocate(mesh%gradient_sca(6,myDim_elem2D))
 allocate(mesh%gradient_vec(6,myDim_elem2D))
 allocate(mesh%metric_factor(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(mesh%elem_cos(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(coriolis(myDim_elem2D))
 allocate(coriolis_node(myDim_nod2D+eDim_nod2D))
 allocate(mesh%geo_coord_nod2D(2,myDim_nod2D+eDim_nod2D))
 allocate(center_x(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(center_y(myDim_elem2D+eDim_elem2D+eXDim_elem2D)) 


 ! ============
 ! coriolis
 ! ============
 DO n=1,myDim_nod2D+eDim_nod2D 
    call r2g(lon, lat, mesh%coord_nod2D(1,n), mesh%coord_nod2D(2,n))
    coriolis_node(n)=2*omega*sin(lat)	 
 END DO

 DO n=1,myDim_nod2D+eDim_nod2D 
    call r2g(lon, lat, mesh%coord_nod2D(1,n), mesh%coord_nod2D(2,n))
    ! in case of numerical noise at the boundaries
    if (lon > 2._WP*pi) lon=lon-2._WP*pi
    if (lon <-2._WP*pi) lon=lon+2._WP*pi
    mesh%geo_coord_nod2D(1,n)=lon
    mesh%geo_coord_nod2D(2,n)=lat	 
 END DO
 

 DO n=1,myDim_elem2D 
    call elem_center(n, ax, ay, mesh)
    call r2g(lon, lat, ax, ay)
    coriolis(n)=2*omega*sin(lat)	 
 END DO
  
 if(fplane) then 
    coriolis=2*omega*0.71_WP
 end if

 ! ============
 ! cos on elements + metric factor (tan/R_earth) 
 ! ============
 DO n=1,myDim_elem2D
 call elem_center(n, ax, ay, mesh)
 center_x(n)=ax
 center_y(n)=ay
 mesh%elem_cos(n)=cos(ay)
 mesh%metric_factor=tan(ay)/r_earth
 END DO

 call exchange_elem(mesh%metric_factor)
 call exchange_elem(mesh%elem_cos)
 call exchange_elem(center_x)
 call exchange_elem(center_y)  
 if (cartesian) then
    mesh%elem_cos=1.0_WP
    mesh%metric_factor=0.0_WP
 end if

 ! ===========
 ! Distances along the edge
 ! We need them in radian measure!
 ! ===========
 DO n=1, myDim_edge2D+eDim_edge2D
    ed=mesh%edges(:,n)
    a=mesh%coord_nod2D(:,ed(2))-mesh%coord_nod2D(:, ed(1))
    call trim_cyclic(a(1))
      !a(1)=a(1)*aux_cos_edge(n)
      !a=a*r_earth
    mesh%edge_dxdy(:,n)=a
 END DO

 ! ===========
 ! Cross-distances for the edge
 ! They are in physical measure!
 ! ===========
 DO n=1, myDim_edge2D+eDim_edge2D    !!! Difference to FESOM2:
    ! Since we know elem centers on the extended halo of elements
    ! the computations can be carried out for all edges (owned and
    ! the halo ones). 
    ed=mesh%edges(:,n)
    el=mesh%edge_tri(:,n)
    call edge_center(ed(1), ed(2), a(1), a(2), mesh)
    b(1)=center_x(el(1))
    b(2)=center_y(el(1))
    b=b-a
    call trim_cyclic(b(1))
    b(1)=b(1)*mesh%elem_cos(el(1))
    b=b*r_earth
    mesh%edge_cross_dxdy(1:2,n)=b(1:2)
 
   if(el(2)>0) then
      b(1)=center_x(el(2)) 
      b(2)=center_y(el(2))
      b=b-a
      call trim_cyclic(b(1))
      b(1)=b(1)*mesh%elem_cos(el(2))
      b=b*r_earth
      mesh%edge_cross_dxdy(3:4,n)=b(1:2)
   else
      mesh%edge_cross_dxdy(3:4,n)=0.0_WP
   end if
 END DO

! ==========================
! Derivatives of scalar quantities
! ==========================
!_______________________________________________________________________________ 
! calculate gradient/derivative of scalar quantitie via linear shape function  
!                3 (x3,y3)
!                o
!               /'\`
!              /   \
!             /     \
!            /       \
!          '/´        \
!(x1,y1) 1 o---------->o 2 (x2,y2)
!
! f(x,y)      = a0 + a1*x + a2*y --> determine coefficients a0,a1,a2 with nodal 
!                                    values f1,f2,f3 as conditions 
! f(N1,N2,N3) = f_1*N_1 + f_2*N_2 + f_3*N_3 = sum_i=1->3( N_i*f_i)
! N_i... linear 2d triangular basis/shape function 
!        
!_______________________________________________________________________________ 
! N_1(x1,y1) = 1 ; N_1(x2,y2)=N_1(x3,y3)=0
! N_2(x2,y2) = 1 ; N_2(x1,y1)=N_2(x3,y3)=0
! N_3(x3,y3) = 1 ; N_3(x2,y2)=N_3(x1,y1)=0
!
!_______________________________________________________________________________ 
! Coordinate Transform: triangular --> cartesian coord.
! |1|   |  1  1  1 |   |N_1|      1st eq.: N_1    + N_2    + N_3    = 1
! |x| = | x1 x2 x3 | * |N_2|  --> 2nd eq.: N_1*x1 + N_2*x2 + N_3*x3 = x
! |y|   | y1 y2 y3 |   |N_3|      3rd eq.: N_1*y1 + N_2*y2 + N_3*y3 = y
!  |___> cartesian --> triangular
!
! |N_1|          |1|             | x2y3-x3y2 y2-y3 x3-x2 |   |1|
! |N_2| = A^-1 * |x| = 1/2*A_tri*| x3y1-x1y3 y3-y1 x1-x3 | * |x| ;  dy23=y2-y3=-dy32
! |N_3|          |y|             | x1y2-x2y1 y1-y2 x2-x1 |   |y|
!
! dN2/dx = dy31/(2A_tri) , dN3/dx = -dy21/(2A_tri), dN1/dx = -dN2/dx-dN3/dx
!                                                          = (-dy31+dy21)/(2*A_tri)
! dN2/dy = dx31/(2A_tri) , dN3/dy = -dx21/(2A_tri), dN1/dy = -dN2/dy-dN3/dy
!
!_______________________________________________________________________________ 
! gradf =  f1*gradN_1 + f2*gradN_2 + f3*gradN_3
!          --> gradN_1 + gradN_2 + gradN_3 = 0
!       = -f1*(gradN_2+gradN_3) + f2*gradN_2 + f3*gradN_3
!
DO elem=1, myDim_elem2D
   elnodes = mesh%elem2D_nodes(:,elem)
   
   deltaX31 = mesh%coord_nod2D(1,elnodes(3)) - mesh%coord_nod2D(1,elnodes(1))
   call trim_cyclic(deltaX31)
   deltaX31 = mesh%elem_cos(elem)*deltaX31
   
   deltaX21 = mesh%coord_nod2D(1,elnodes(2)) - mesh%coord_nod2D(1,elnodes(1))
   call trim_cyclic(deltaX21)
   deltaX21 = mesh%elem_cos(elem)*deltaX21
   
   deltaY31 = mesh%coord_nod2D(2,elnodes(3)) - mesh%coord_nod2D(2,elnodes(1))
   deltaY21 = mesh%coord_nod2D(2,elnodes(2)) - mesh%coord_nod2D(2,elnodes(1))
   
   dfactor = -0.5_WP*r_earth/mesh%elem_area(elem)
   mesh%gradient_sca(1,elem)=(-deltaY31+deltaY21)*dfactor
   mesh%gradient_sca(2,elem)=deltaY31*dfactor
   mesh%gradient_sca(3,elem)=-deltaY21*dfactor
   
   mesh%gradient_sca(4,elem)=(deltaX31-deltaX21)*dfactor
   mesh%gradient_sca(5,elem)=-deltaX31*dfactor
   mesh%gradient_sca(6,elem)=deltaX21*dfactor
END DO

! ==========================
! Derivatives of vector quantities
! Least squares interpolation is used
! ==========================
!_______________________________________________________________________________ 
! Least square method for gradient reconstruction of elemental velocities
!     o___________o___________o
!      \   P_3   / \   P_1   /          dx_10 = x_1-x_0 ; dy_10 = y_1-y_0
!       \   x   /   \   x   /           dx_20 = x_2-x_0 ; dy_20 = y_2-y_0
!        \     / P_0 \     /            dx_30 = x_3-x_0 ; dy_30 = y_3-y_0   
!         \   /   x   \   /
!          \ /         \ /
!           o-----------o      P_1 = P_0 + dx_10*(dP/dx)_0 + dy_10*(dP_dy)_0
!            \         /
!             \   x   /        P_2 = P_0 + dx_20*(dP/dx)_0 + dy_20*(dP_dy)_0
!              \ P_2 / 
!               \   /          P_3 = P_0 + dx_30*(dP/dx)_0 + dy_30*(dP_dy)_0
!                \ /
!                 o            --> How to calculate dP/dx and dP/dy from P0, P1, 
!                                  P2, P3 ... ?
!
!
! 1st. write as linear equation system ... A*x=f
!
!        x = | dP/dx | ;  f = | P_1 - P_0 | ;  A = | dx_10 dy_10 | 
!            | dP/dy |        | P_2 - P_0 |        | dx_20 dy_20 |
!                             | P_3 - P_0 |        | dx_30 dy_30 |
!   ---> there are 3 equations but only 2 unknowns (dP/dx,dP/dy), no solution to
!        an overestimated problem
!   ---> Least Square Method can be used
!
! 2nd.             A * x =       f  | A^T* , A^T transpose of A
!              A^T*A * x = A^T * f
!                      x = (A^T*A)^-1 * A^T *f 
!
!               
!                A^T = | dx_10 dx_20 dx_30 | 
!                      | dy_10 dy_20 dy_30 |
!
!              A^T*A = | dx_10^2+dx_20^2+dx_30^2            ; dx_10*dy_10+dx_20*dy_20+dx_30*dy_30 |
!                      | dx_10*dy_10+dx_20*dy_20+dx_30*dy_30; dx_10^2+dx_20^2+dx_30^2             |
!
!                    = | sum(dx^2)  ; sum(dx*dy) |
!                      | sum(dx*dy) ; sum(dy^2)  |                    
! 
!          (A^T*A)⁻1 = 1/(sum(dx^2)*sum(dy^2)-sum(dx*dy)^2)* |  sum(dy^2)  -sum(dx*dy) |
!                                      |                     | -sum(dx*dy)  sum(dx^2)  |
!                                      V
!                                     DET
!
!    (A^T*A)⁻1 * A^T = 1/DET * |  sum(dy^2)  -sum(dx*dy) | * | dx_10 dx_20 dx_30 |
!                              | -sum(dx*dy)  sum(dx^2)  |   | dy_10 dy_20 dy_30 |
!
!  dP/dx =  1/DET *(sum(dy^2)*dx_10 - sum(dx*dy)*dy_10) * (P_1-P_0) 
!         + 1/DET *(sum(dy^2)*dx_20 - sum(dx*dy)*dy_20) * (P_2-P_0) 
!         + 1/DET *(sum(dy^2)*dx_30 - sum(dx*dy)*dy_30) * (P_3-P_0) 
!
!  dP/dy =  1/DET *(sum(dx^2)*dy_10 - sum(dx*dy)*dx_10) * (P_1-P_0) 
!         + 1/DET *(sum(dx^2)*dy_20 - sum(dx*dy)*dx_20) * (P_2-P_0) 
!         + 1/DET *(sum(dx^2)*dy_30 - sum(dx*dy)*dx_30) * (P_3-P_0) 
!
DO elem=1,myDim_elem2D
    !elnodes=elem2D_nodes(:,elem)
    a(1)=center_x(elem)
    a(2)=center_y(elem)
    DO j=1,3
        el(1)=mesh%elem_neighbors(j,elem)
        if (el(1)>0) then
                !elnodes_=mesh%elem2D_nodes(:,el(1))
            b(1)=center_x(el(1))
            b(2)=center_y(el(1))
            x(j)=b(1)-a(1)
            call trim_cyclic(x(j))
            y(j)=b(2)-a(2)
        else
            ! Virtual element center is taken
            ed=mesh%edges(:,mesh%elem_edges(j,elem))
            call edge_center(ed(1), ed(2), b(1), b(2), mesh)
            x(j)=(b(1)-a(1))
            call trim_cyclic(x(j))
            x(j)=2*x(j)
            y(j)=2*(b(2)-a(2))
        end if
    END DO
    x=x*mesh%elem_cos(elem)*r_earth
    y=y*r_earth
    cxx=sum(x**2)
    cxy=sum(x*y)
    cyy=sum(y**2)
    d=cxy*cxy-cxx*cyy
    ! coefficients to compute gradients of velocity
    mesh%gradient_vec(1:3,elem)=(cxy*y-cyy*x)/d
    mesh%gradient_vec(4:6,elem)=(cxy*x-cxx*y)/d
END DO
deallocate(center_y, center_x)

    !array of 2D boundary conditions is used in ice_maEVP
    if (whichEVP > 0) then
       allocate(mesh%bc_index_nod2D(myDim_nod2D+eDim_nod2D))
       mesh%bc_index_nod2D=1._WP
       do n=1, myDim_edge2D
          ed=mesh%edges(:, n)
          if (myList_edge2D(n) <= mesh%edge2D_in) cycle
          mesh%bc_index_nod2D(ed)=0._WP
       end do
    end if

#if defined (__oasis)
  nn=0
  ns=0  
  allocate(mesh%lump2d_north(myDim_nod2D), mesh%lump2d_south(myDim_nod2D))
  mesh%lump2d_north=0._WP
  mesh%lump2d_south=0._WP
  do i=1, myDim_nod2D
     if (mesh%geo_coord_nod2D(2, i) > 0) then
        nn=nn+1
        mesh%lump2d_north(i)=mesh%areasvol(mesh%ulevels_nod2d(i), i)
     else
        ns=ns+1     
        mesh%lump2d_south(i)=mesh%area(mesh%ulevels_nod2d(i), i)
     end if	   
  end do   

  if (nn>0) allocate(mesh%ind_north(nn))
  if (ns>0) allocate(mesh%ind_south(ns))
  ns=0
  nn=0
  do i=1, myDim_nod2D
     if (mesh%geo_coord_nod2D(2, i) > 0) then
        nn=nn+1
	mesh%ind_north(nn)=i
     else
        ns=ns+1
	mesh%ind_south(ns)=i	
     end if	     
  end do     
#endif 

    t1=MPI_Wtime()
    if (mype==0) then
       write(*,*) 'mesh_auxiliary_arrays finished in ', t1-t0, ' seconds'
       write(*,*) '========================='
    endif

END SUBROUTINE mesh_auxiliary_arrays
!
!
!_______________________________________________________________________________
SUBROUTINE check_mesh_consistency(mesh)
USE MOD_MESH
USE o_PARAM
USE g_PARSUP
USE g_ROTATE_GRID
  use g_comm_auto
IMPLICIT NONE
! Collects auxilliary information on the mesh
! Allocated and filled in are:
! elem_area(myDim_elem2D)
! area(nl, myDim_nod2D)
type(t_mesh), intent(inout), target :: mesh
integer                     :: nz, n, elem , elnodes(3)
real(kind=WP)	            :: vol_n(mesh%nl), vol_e(mesh%nl), aux(mesh%nl)


   vol_n=0._WP
   vol_e=0._WP

   aux=0._WP
   do n=1, myDim_nod2D
      do nz=mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
         aux(nz)=aux(nz)+mesh%areasvol(nz, n)
      end do
   end do
   call MPI_AllREDUCE(aux, vol_n, mesh%nl, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

   aux=0._WP
   do elem=1, myDim_elem2D
      elnodes=mesh%elem2D_nodes(:, elem)
      if (elnodes(1) > myDim_nod2D) CYCLE
      do nz=mesh%ulevels(elem), mesh%nlevels(elem)-1         
         aux(nz)=aux(nz)+mesh%elem_area(elem)
      end do
   end do
   call MPI_AllREDUCE(aux, vol_e, mesh%nl, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

if (mype==0) then
write(*,*) '***start level area_test***'
do nz=1, mesh%nl
   write(*,*) vol_n(nz), vol_e(nz)
end do
write(*,*) '***end level area_test***'
end if

!call par_ex
!stop
END SUBROUTINE check_mesh_consistency
!
!
!_______________________________________________________________________________
subroutine check_total_volume(mesh)
    USE MOD_MESH
    USE o_PARAM
    USE g_PARSUP
    use g_comm_auto
    use o_ARRAYS
    
    IMPLICIT NONE
    type(t_mesh), intent(inout), target :: mesh
    integer                     :: nz, n, elem , elnodes(3)
    real(kind=WP)	            :: vol_n, vol_e, aux
    
#include "associate_mesh.h"

    !___________________________________________________________________________
    vol_n=0._WP
    vol_e=0._WP
    !___________________________________________________________________________
    ! total ocean volume on nodes
    aux=0._WP
    do n=1, myDim_nod2D
        do nz=ulevels_nod2D(n), nlevels_nod2D(n)-1
            aux=aux+areasvol(nz, n)*hnode(nz,n)
        end do
    end do
    call MPI_AllREDUCE(aux, vol_n, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    !___________________________________________________________________________
    ! total ocean volume on elements
    aux=0._WP
    do elem=1, myDim_elem2D
        elnodes=elem2D_nodes(:, elem)
        if (elnodes(1) > myDim_nod2D) cycle
        do nz=ulevels(elem), nlevels(elem)-1         
            aux=aux+elem_area(elem)*helem(nz,elem)
        end do
    end do
    call MPI_AllREDUCE(aux, vol_e, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)

    !___write mesh statistics___________________________________________________
    if (mype==0) then
        write(*,*) '____________________________________________________________________'
        write(*,*) ' --> ocean volume check:', mype
        write(*,*) '     > Total ocean volume node:', vol_n, ' m^3'
        write(*,*) '     > Total ocean volume elem:', vol_e, ' m^3'
        
    end if

end subroutine check_total_volume
!
!
!_______________________________________________________________________________
