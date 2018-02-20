! Driving routine. The distributed mesh information and mesh proper 
! are read from files.
! Auxiliary arrays with mesh information are assembled.
! At the beginning of each routine I list arrays it initializes.
! Array sizes vary (sometimes we need only myDim, yet sometimes more)! 
! S. Danilov, 2012
SUBROUTINE mesh_setup
USE g_parsup
USE g_ROTATE_grid
IMPLICIT NONE
      call set_mesh_transform_matrix  !(rotated grid)
      call read_mesh
      call set_par_support
      call find_levels
      call test_tri
      call load_edges
      call find_neighbors
      call mesh_areas
      call mesh_auxiliary_arrays
END SUBROUTINE mesh_setup
!======================================================================
! Reads distributed mesh
! The mesh will be read only by 0 proc and broadcasted to the others.
SUBROUTINE read_mesh
USE o_PARAM
USE g_CONFIG
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_rotate_grid 
IMPLICIT NONE

 integer        :: n, nn, k, m, fileID
 integer        :: error_status !0/1=no error/error
 integer        :: vert_nodes(1000)
 integer        :: nchunk, chunk_size, ipos, iofs, mesh_check
 real(kind=WP)  :: x, y, rx, ry
 real(kind=WP)  :: t0, t1
 character*10   :: mype_string,npes_string
 character*500  :: file_name
 character*500  :: dist_mesh_dir
 integer       :: ierror              ! return error code
 integer, allocatable, dimension(:)        :: mapping
 integer, allocatable, dimension(:,:)      :: ibuff
 real(kind=8), allocatable, dimension(:,:) :: rbuff
 integer, allocatable, dimension(:,:)      :: auxbuff ! will be used for reading aux3d.out 


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
  fileID=10+mype  
    
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
  nod2D=part(npes+1)-1
  allocate(coord_nod2D(2,myDim_nod2D+eDim_nod2D))
  if (mype==0) then
    file_name=trim(meshpath)//'nod2d.out'
    open(fileID, file=file_name)
    read(fileID,*) n      ! nod2D, we know it already
     error_status=0
     if (n/=nod2D) error_status=1 !set the error status for consistency between rpart and nod2D
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
  do nchunk=0, (nod2D-1)/chunk_size
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
     k=min(chunk_size, nod2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1), rbuff(n,1), rbuff(n,2), ibuff(n,2)
        end do
     end if
     call MPI_BCast(rbuff(1:k,1), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(rbuff(1:k,2), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     ! fill the local arrays
     do n=1, k
        x=rbuff(n,1)*rad
        y=rbuff(n,2)*rad
        if (force_rotation) then
           rx=x
           ry=y
           call g2r(rx, ry, x, y)
        end if        
        if (mapping(n)>0) then
           mesh_check=mesh_check+1
           coord_nod2D(1,mapping(n))=x
           coord_nod2D(2,mapping(n))=y
        end if
     end do
  end do
  if (mype==0) close(fileID)

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
     read(fileID,*) elem2d
     write(*,*) 'reading '// trim(file_name)   
  end if
  call MPI_BCast(elem2d, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  allocate(elem2D_nodes(3, myDim_elem2D))

  ! 0 proc reads the data in chunks and distributes it between other procs
  do nchunk=0, (elem2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_elem2D
        ipos=(myList_elem2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, elem2D-nchunk*chunk_size)
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
           elem2D_nodes(1,mapping(n))=ibuff(n,1)
           elem2D_nodes(2,mapping(n))=ibuff(n,2)
           elem2D_nodes(3,mapping(n))=ibuff(n,3)
        end if
     end do
  end do
  if (mype==0) close(fileID)
  ! nodes in elem2d are in global numbering. convert to local:
  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_elem2D
        do m=1,3
           nn=elem2D_nodes(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk) then
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              elem2D_nodes(m,n)=-mapping(iofs) 
           end if
        end do
     end do
  end do
  elem2D_nodes=-elem2D_nodes
 if (mype==0) write(*,*) 'elements are read' 
 !==============================
 ! read depth data
 !==============================
 ! 0 proc reads header of aux3d.out and broadcasts it between procs
 allocate(depth(myDim_nod2D+eDim_nod2D))
 if (mype==0) then !open the file for reading on 0 proc
    file_name=trim(meshpath)//'aux3d.out' 
    open(fileID, file=file_name)
    read(fileID,*) nl  ! the number of levels 
 end if
 call MPI_BCast(nl, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
 allocate(zbar(nl))              ! allocate the array for storing the standard depths
 if (mype==0) read(fileID,*) zbar
 call MPI_BCast(zbar, nl, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
 if(zbar(2)>0) zbar=-zbar   ! zbar is negative 
 allocate(Z(nl-1))
 Z=zbar(1:nl-1)+zbar(2:nl)  ! mid-depths of cells
 Z=0.5_WP*Z 

 ! 0 proc reads the data in chunks and distributes it between other procs
 mesh_check=0
 do nchunk=0, (nod2D-1)/chunk_size
    mapping(1:chunk_size)=0
    do n=1, myDim_nod2D+eDim_nod2D
       ipos=(myList_nod2D(n)-1)/chunk_size
       if (ipos==nchunk) then
          iofs=myList_nod2D(n)-nchunk*chunk_size
          mapping(iofs)=n
       end if
    end do

    k=min(chunk_size, nod2D-nchunk*chunk_size)
    if (mype==0) then
       do n=1, k
          read(fileID,*) rbuff(n,1)
       end do
    end if
    call MPI_BCast(rbuff(1:k,1), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

    do n=1, k
       x=rbuff(n,1)
       if (x>0) x=-x !deps must be negative!
       if (x>zbar(5)) x=zbar(5) !threshold for depth
       if (mapping(n)>0) then
          mesh_check=mesh_check+1
          depth(mapping(n))=x
        end if
     end do
 end do

 if (mype==0) close(fileID)
 if (mesh_check/=myDim_nod2D+eDim_nod2D) then
    write(*,*) 'ERROR while reading aux3d.out on mype=', mype
    write(*,*) mesh_check, ' values have been read in according to partitioning'
    write(*,*) 'it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
 end if

 ! ==============================
 ! Communication information
 ! every proc reads its file
 ! ==============================
 file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out'  
 fileID=10+mype  
 open(fileID, file=file_name)
 read(fileID,*)  n
 read(fileID,*) com_nod2D%rPEnum
 ALLOCATE(com_nod2D%rPE(com_nod2D%rPEnum))
 read(fileID,*) com_nod2D%rPE
 ALLOCATE(com_nod2D%rptr(com_nod2D%rPEnum+1))
 read(fileID,*) com_nod2D%rptr
 ALLOCATE(com_nod2D%rlist(eDim_nod2D))
 read(fileID,*) com_nod2D%rlist
	 
 read(fileID,*) com_nod2D%sPEnum
 ALLOCATE(com_nod2D%sPE(com_nod2D%sPEnum))
 read(fileID,*) com_nod2D%sPE
 ALLOCATE(com_nod2D%sptr(com_nod2D%sPEnum+1))
 read(fileID,*) com_nod2D%sptr
 n=com_nod2D%sptr(com_nod2D%sPEnum+1)-1
 ALLOCATE(com_nod2D%slist(n))
 read(fileID,*) com_nod2D%slist
	 
 read(fileID,*) com_elem2D%rPEnum
 ALLOCATE(com_elem2D%rPE(com_elem2D%rPEnum))
 read(fileID,*) com_elem2D%rPE
 ALLOCATE(com_elem2D%rptr(com_elem2D%rPEnum+1))
 read(fileID,*) com_elem2D%rptr
 ALLOCATE(com_elem2D%rlist(eDim_elem2D))
 read(fileID,*) com_elem2D%rlist
	 
 read(fileID,*) com_elem2D%sPEnum
 ALLOCATE(com_elem2D%sPE(com_elem2D%sPEnum))
 read(fileID,*) com_elem2D%sPE
 ALLOCATE(com_elem2D%sptr(com_elem2D%sPEnum+1))
 read(fileID,*) com_elem2D%sptr
 n=com_elem2D%sptr(com_elem2D%sPEnum+1)-1
 ALLOCATE(com_elem2D%slist(n))
 read(fileID,*) com_elem2D%slist
	 
 read(fileID,*) com_elem2D_full%rPEnum
 ALLOCATE(com_elem2D_full%rPE(com_elem2D_full%rPEnum))
 read(fileID,*) com_elem2D_full%rPE
 ALLOCATE(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1))
 read(fileID,*) com_elem2D_full%rptr
 ALLOCATE(com_elem2D_full%rlist(eDim_elem2D+eXDim_elem2D))
 read(fileID,*) com_elem2D_full%rlist
	 
 read(fileID,*) com_elem2D_full%sPEnum
 ALLOCATE(com_elem2D_full%sPE(com_elem2D_full%sPEnum))
 read(fileID,*) com_elem2D_full%sPE
 ALLOCATE(com_elem2D_full%sptr(com_elem2D_full%sPEnum+1))
 read(fileID,*) com_elem2D_full%sptr
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
CALL MPI_BARRIER(MPI_COMM_FESOM, MPIerr)
 t1=MPI_Wtime()
 if (mype==0) then
    write(*,*) '========================='
    write(*,*) '2D mesh was read in ', t1-t0, ' seconds'
    write(*,*) '2D mesh info : ', 'nod2D=', nod2D,' elem2D=', elem2D
    write(*,*) '========================='
 endif

 END subroutine  read_mesh
!============================================================ 
subroutine find_levels
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_config
!
IMPLICIT NONE
!
 character*1000                         :: file_name
 integer                                :: ierror   ! MPI return error code
 integer                                :: k, n, fileID
 integer                                :: nchunk, chunk_size, ipos, iofs, mesh_check
 integer, allocatable, dimension(:)     :: mapping
 integer, allocatable, dimension(:)     :: ibuff
 real(kind=WP)                          :: t0, t1

 t0=MPI_Wtime()
 allocate(nlevels(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(nlevels_nod2D(myDim_nod2D+eDim_nod2D))
 !mesh related files will be read in chunks of chunk_size
 chunk_size=100000
 !==============================
 ! Allocate mapping array (chunk_size)
 ! It will be used for several purposes 
 !==============================
 allocate(mapping(chunk_size))
 allocate(ibuff(chunk_size))
 !==============================
 !Part I: reading levels at elements...
 if (mype==0)  then 
    fileID=10
    file_name=trim(meshpath)//'elvls.out'
    open(fileID, file=file_name)
    write(*,*) 'reading '// trim(file_name)   
 end if
 ! 0 proc reads the data in chunks and distributes it between other procs
 mesh_check=0
 do nchunk=0, (elem2D-1)/chunk_size
    !create the mapping for the current chunk
    mapping(1:chunk_size)=0
    do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
       ipos=(myList_elem2D(n)-1)/chunk_size
       if (ipos==nchunk) then
          iofs=myList_elem2D(n)-nchunk*chunk_size
          mapping(iofs)=n
       end if
    end do
    !read the chunk into the buffers
    k=min(chunk_size, elem2D-nchunk*chunk_size)
    if (mype==0) then
       do n=1, k
          read(fileID,*) ibuff(n)
       end do
    end if
    call MPI_BCast(ibuff(1:k), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
    ! fill the local arrays
    do n=1, k      
       if (mapping(n)>0) then
          mesh_check=mesh_check+1
          nlevels(mapping(n))=ibuff(n)
       end if
    end do
 end do
 if (mype==0) close(fileID)
 if (mesh_check/=myDim_elem2D+eDim_elem2D+eXDim_elem2D) then
    write(*,*) 'ERROR while reading elvls.out on mype=', mype
    write(*,*) mesh_check, ' values have been read in according to partitioning'
    write(*,*) 'it does not equal to myDim_elem2D+eDim_elem2D = ', myDim_elem2D+eDim_elem2D
 end if

 !==============================
 !Part II: reading levels at nodes...
 if (mype==0)  then 
    file_name=trim(meshpath)//'nlvls.out'
    open(fileID, file=file_name)
    write(*,*) 'reading '// trim(file_name)   
 end if
 ! 0 proc reads the data in chunks and distributes it between other procs
 mesh_check=0
 do nchunk=0, (nod2D-1)/chunk_size
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
    k=min(chunk_size, nod2D-nchunk*chunk_size)
    if (mype==0) then
       do n=1, k
          read(fileID,*) ibuff(n)
       end do
    end if
    call MPI_BCast(ibuff(1:k), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
    ! fill the local arrays
    do n=1, k      
       if (mapping(n)>0) then
          mesh_check=mesh_check+1
          nlevels_nod2D(mapping(n))=ibuff(n)
       end if
    end do
 end do
 if (mype==0) close(fileID)
 if (mesh_check/=myDim_nod2D+eDim_nod2D) then
    write(*,*) 'ERROR while reading nelvls.out on mype=', mype
    write(*,*) mesh_check, ' values have been read in according to partitioning'
    write(*,*) 'it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
 end if
 deallocate(ibuff)
 deallocate(mapping)
 !============================== 
CALL MPI_BARRIER(MPI_COMM_FESOM, MPIerr)
 t1=MPI_Wtime()

 if (mype==0) then
    write(*,*) '3D mesh was read in ', t1-t0, ' seconds'
    write(*,*) 'Min/max depth on mype : ', mype, -zbar(minval(nlevels)),-zbar(maxval(nlevels))
    write(*,*) '========================='
 endif
END SUBROUTINE find_levels
!===========================================================================
SUBROUTINE test_tri
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
! Check the order of nodes in triangles; correct it if necessary to make
! it same sense (clockwise) 
real(kind=WP)   ::  a(2), b(2), c(2),  r
integer         ::  n, nx, elnodes(3)
real(kind=WP)   :: t0, t1

   t0=MPI_Wtime()
   
   DO n=1, myDim_elem2D
      elnodes=elem2D_nodes(:,n)
	  
          a=coord_nod2D(:,elnodes(1))
	  b=coord_nod2D(:,elnodes(2))-a
	  c=coord_nod2D(:,elnodes(3))-a
          
	  if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
          if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
	  if(c(1)>cyclic_length/2.) c(1)=c(1)-cyclic_length
          if(c(1)<-cyclic_length/2.) c(1)=c(1)+cyclic_length
	  
	    
	  r=b(1)*c(2)-b(2)*c(1)
	  if (r>0) then
	  ! Vector b is to right of c
	  ! Exchange second and third nodes:
	  
	  nx=elnodes(2)
	  elnodes(2)=elnodes(3)
	  elnodes(3)=nx
	  elem2D_nodes(:,n)=elnodes
      end if
   END DO
   t1=MPI_Wtime()
   if (mype==0) then
      write(*,*) 'test_tri finished in ', t1-t0, ' seconds'
      write(*,*) '========================='
   endif
END SUBROUTINE  test_tri
!=========================================================================
SUBROUTINE load_edges
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
character*1000                        :: file_name
integer                               :: counter, n, m, nn, k, q, fileID
integer                               :: elems(2), elem
integer                               :: elnodes(3), ed(2), eledges(3)
integer, allocatable                  :: aux(:)         
real(kind=WP)                         :: t0, t1
integer                               :: nchunk, chunk_size, ipos, iofs, mesh_check
integer, allocatable, dimension(:)    :: mapping
integer, allocatable, dimension(:,:)  :: ibuff
integer                               :: ierror              ! return error code

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
    read(fileID,*) edge2D
    read(fileID,*) edge2D_in
    write(*,*) '2D mesh info : edge2D=', edge2D
    close(fileID)
  end if
  call MPI_BCast(edge2D,    1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(edge2D_in, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

  allocate(edges(2,myDim_edge2D+eDim_edge2D))
  allocate(edge_tri(2,myDim_edge2D+eDim_edge2D))

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
 do nchunk=0, (edge2D-1)/chunk_size
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
    k=min(chunk_size, edge2D-nchunk*chunk_size)
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
          edges   (:, mapping(n))=ibuff(1:2,n)
          edge_tri(:, mapping(n))=ibuff(3:4,n)
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
  do nchunk=0, (nod2D-1)/chunk_size
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
           nn=edges(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk) then
              mesh_check=mesh_check+1
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              edges(m,n)=-mapping(iofs) 
           end if
        end do
     end do
  end do
  edges=-edges
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
  where (edge_tri<0)
         edge_tri=0
  end where
!***
  mesh_check=0
  do nchunk=0, (elem2D-1)/chunk_size
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
           nn=edge_tri(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk .and. nn > 0) then
              mesh_check=mesh_check+abs(m-2) !only first triangle will contribute to statistic
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              edge_tri(m,n)=-mapping(iofs) 
           end if
        end do
     end do
  end do
  edge_tri=-edge_tri
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
allocate(elem_edges(3,myDim_elem2D))
allocate(aux(myDim_elem2D))
aux=0
DO n=1, myDim_edge2D+eDim_edge2D
   DO k=1,2
      q=edge_tri(k,n)   ! triangle number
	  if((q>0).and.(q<=myDim_elem2D)) then
	  aux(q)=aux(q)+1
	  elem_edges(aux(q),q)=n
	  end if
   END DO
END DO
deallocate(aux)
! The edges in this list should be ordered so that they
! are listed in the same rotation sense as nodes.
DO elem=1,myDim_elem2D
   elnodes=elem2D_nodes(:,elem)
   eledges=elem_edges(:,elem)
   DO q=1,3
      DO k=1,3
         if((edges(1,eledges(k)).ne.elnodes(q)).and. &
            (edges(2,eledges(k)).ne.elnodes(q))) then
           elem_edges(q,elem)=eledges(k)
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
SUBROUTINE find_neighbors
! For each element three its element neighbors are found
! For each node the elements containing it are found
! Allocated are:
! elem_neighbors(3,myDim_elem2D)
! nod_in_elem2D_num(myDim_nod2D)
! nod_in_elem2D(:, myDim_nod2D)
! 

USE o_PARAM
USE o_MESH
USE g_PARSUP
  use g_comm_auto
implicit none
integer               :: elem, eledges(3), elem1, j, n, node, enum,elems(3),count1,count2,exit_flag,i,nz
integer, allocatable  :: temp_i(:)
integer               :: mymax(npes), rmax(npes)
real(kind=WP)         :: t0, t1
CALL MPI_BARRIER(MPI_COMM_FESOM, MPIerr)
t0=MPI_Wtime()

 ! =============
 ! elem neighbors == those that share edges
 ! =============
   allocate(elem_neighbors(3,myDim_elem2D))
   elem_neighbors=0
  
DO elem=1,myDim_elem2D
   eledges=elem_edges(:,elem)
   DO j=1,3
   elem1=edge_tri(1,eledges(j))
   if(elem1==elem) elem1=edge_tri(2,eledges(j))
   elem_neighbors(j,elem)=elem1
   END DO
END DO
 ! =============
 ! Node neighbourhood
 ! == elements that contain node n
 ! We need eDim neighborhood too for MUSCL advection. 
 ! And we already have the place allocated for all 
 ! these neighbor elements: it is eDim_elem2D+eXDim_elem2D
 ! =============	 
 allocate(nod_in_elem2D_num(myDim_nod2D+eDim_nod2D)) 
 nod_in_elem2D_num=0
 do n=1,myDim_elem2D
    do j=1,3
    node=elem2D_nodes(j,n)
    if (node>myDim_nod2D) cycle
    nod_in_elem2D_num(node)=nod_in_elem2D_num(node)+1
    end do
 end do
CALL MPI_BARRIER(MPI_COMM_FESOM, MPIerr)

 mymax=0
 rmax=0
 mymax(mype+1)=maxval(nod_in_elem2D_num(1:myDim_nod2D))
 call MPI_AllREDUCE( mymax, rmax, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
 
 allocate(nod_in_elem2D(maxval(rmax),myDim_nod2D+eDim_nod2D))
 nod_in_elem2D=0
 nod_in_elem2D_num=0
 do n=1,myDim_elem2D   
    do j=1,3
    node=elem2D_nodes(j,n)
    if (node>myDim_nod2D) cycle 
    nod_in_elem2D_num(node)=nod_in_elem2D_num(node)+1
    nod_in_elem2D(nod_in_elem2D_num(node),node)=n
    end do
 end do

 call exchange_nod(nod_in_elem2D_num)
 allocate (temp_i(myDim_nod2D+eDim_nod2D))
 temp_i=0
 DO n=1, maxval(rmax)
       ! Exchange global element numbers
       do j=1,myDim_nod2D
         if (nod_in_elem2D(n,j)>0) temp_i(j)=myList_elem2D(nod_in_elem2D(n,j))
       enddo
       call exchange_nod(temp_i)
       nod_in_elem2D(n,:)=temp_i
 END DO
 deallocate(temp_i)
 ! Substitute back local element numbers
 allocate(temp_i(elem2D))
 temp_i=0
 Do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
    temp_i(myList_elem2D(n))=n
 END DO
 DO n=1, myDim_nod2D+eDim_nod2D     
    DO j=1, nod_in_elem2D_num(n)
       nod_in_elem2D(j,n)=temp_i(nod_in_elem2D(j,n))
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
   if(elem_neighbors(j,elem)>0) elem1=elem1+1
   END DO
   if (elem1<2) then
    write(*,*) 'Insufficient number of neighbors ', myList_elem2D(elem)
    call par_ex(1)
    STOP
   end if
END DO
t1=MPI_Wtime()
if (mype==0) then
   write(*,*) 'find_neighbors finished in ', t1-t0, ' seconds'
   write(*,*) '========================='
endif
END SUBROUTINE find_neighbors
!==========================================================================
subroutine edge_center(n1, n2, x, y)
USE o_MESH
USE o_PARAM
USE g_CONFIG 
implicit none
integer      :: n1, n2   ! nodes of the edge
real(kind=WP) :: x, y, a(2), b(2)

a=coord_nod2D(:,n1)
b=coord_nod2D(:,n2)
if(a(1)-b(1)>cyclic_length/2.0) a(1)=a(1)-cyclic_length
if(a(1)-b(1)<-cyclic_length/2.0) b(1)=b(1)-cyclic_length
x=0.5_WP*(a(1)+b(1))
y=0.5_WP*(a(2)+b(2))
end subroutine edge_center
!==========================================================================
subroutine elem_center(elem, x, y)
USE o_MESH
USE o_PARAM
USE g_CONFIG  
implicit none
integer      :: elem, elnodes(3), k    
real(kind=WP) :: x, y, ax(3), amin

   elnodes=elem2D_nodes(:,elem)
   ax=coord_nod2D(1, elnodes)
   amin=minval(ax)
   DO k=1,3
   if(ax(k)-amin>=cyclic_length/2.0) ax(k)=ax(k)-cyclic_length
   if(ax(k)-amin<-cyclic_length/2.0) ax(k)=ax(k)+cyclic_length
   END DO
   x=sum(ax)/3.0_WP
   y=sum(coord_nod2D(2,elnodes))/3.0_WP
   
end subroutine elem_center
!==========================================================================
SUBROUTINE mesh_areas
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_ROTATE_GRID
use g_comm_auto
IMPLICIT NONE
! Collects auxilliary information on the mesh
! Allocated and filled in are:
! elem_area(myDim_elem2D)
! area(nl, myDim_nod2D)


integer                                   :: n,j,q, elnodes(3), ed(2), elem, nz
real(kind=8)	                          :: a(2), b(2), ax, ay, lon, lat, vol
real(kind=WP), allocatable,dimension(:)   :: work_array
real(kind=WP)                             :: t0, t1

t0=MPI_Wtime()

 allocate(elem_area(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 !allocate(elem_area(myDim_elem2D))
 allocate(area(nl,myDim_nod2d+eDim_nod2D))   !! Extra size just for simplicity
                                             !! in some further routines
 allocate(area_inv(nl,myDim_nod2d+eDim_nod2D)) 
 allocate(mesh_resolution(myDim_nod2d+eDim_nod2D))
 ! ============
 ! The areas of triangles:
 ! ============
 DO n=1, myDim_elem2D
 !DO n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
    elnodes=elem2D_nodes(:,n)
    ay=sum(coord_nod2D(2,elnodes))/3.0_WP
    ay=cos(ay)
    if (cartesian) ay=1.0_WP
    a=coord_nod2D(:,elnodes(2))-coord_nod2D(:,elnodes(1))
    b=coord_nod2D(:,elnodes(3))-coord_nod2D(:,elnodes(1))
    if(a(1)>cyclic_length/2.) a(1)=a(1)-cyclic_length
    if(a(1)<-cyclic_length/2.) a(1)=a(1)+cyclic_length
    if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
    if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
    a(1)=a(1)*ay
    b(1)=b(1)*ay
    elem_area(n)=0.5_WP*abs(a(1)*b(2)-b(1)*a(2))
 END DO
 call exchange_elem(elem_area)
 ! =============
 ! Scalar element 
 ! areas at different levels (there can be partly land)
 ! =============
 
 area=0.0_WP
 DO n=1, myDim_nod2D
    DO j=1,nod_in_elem2D_num(n)
       elem=nod_in_elem2D(j,n)
       DO nz=1,nlevels(elem)-1
       area(nz,n)=area(nz,n)+elem_area(elem)/3.0_WP
       END DO
    END DO
 END DO
 
 ! Only areas through which there is exchange are counted

 ! ===========
 ! Update to proper dimension
 ! ===========
 elem_area=elem_area*r_earth*r_earth
 area=area*r_earth*r_earth
 
 call exchange_nod(area)

do n=1,myDim_nod2d+eDim_nod2D
   do nz=1,nl
      if (area(nz,n) > 0._WP) then
         area_inv(nz,n) = 1._WP/area(nz,n)
      else
         area_inv(nz,n) = 0._WP
      end if
   end do
end do
 ! coordinates are in radians, edge_dxdy are in meters,
 ! and areas are in m^2
 

 allocate(work_array(myDim_nod2D))
 mesh_resolution=sqrt(area(1, :)/pi)*2._WP
 DO q=1, 3 !apply mass matrix N times to smooth the field
    DO n=1, myDim_nod2D
       vol=0._WP
       work_array(n)=0._WP
       DO j=1, nod_in_elem2D_num(n)
          elem=nod_in_elem2D(j, n)
          elnodes=elem2D_nodes(:,elem)
          work_array(n)=work_array(n)+sum(mesh_resolution(elnodes))/3._WP*elem_area(elem)
          vol=vol+elem_area(elem)
       END DO
       work_array(n)=work_array(n)/vol
    END DO
    DO n=1,myDim_nod2D
       mesh_resolution(n)=work_array(n)
    ENDDO
    call exchange_nod(mesh_resolution)
 END DO
 deallocate(work_array)

 vol=0.0_WP
 do n=1, myDim_nod2D
    vol=vol+area(1, n)
 end do
 ocean_area=0.0
 call MPI_AllREDUCE(vol, ocean_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

if (mype==0) then
 write(*,*)  mype, 'Mesh statistics:'
 write(*,*)  mype, 'maxArea ',maxval(elem_area), '   MinArea ', minval(elem_area)
 write(*,*)  mype, 'maxScArea ',maxval(area(1,:)), &
            '   MinScArea ', minval(area(1,:))
 write(*,*)  mype, 'Edges:    ', edge2D, ' internal ', edge2D_in
 if (mype==0) then
    write(*,*) 'Total ocean area is: ', ocean_area, ' m^2'
 end if
endif

t1=MPI_Wtime()
if (mype==0) then
   write(*,*) 'mesh_areas finished in ', t1-t0, ' seconds'
   write(*,*) '========================='
endif
END SUBROUTINE mesh_areas

!===================================================================

SUBROUTINE mesh_auxiliary_arrays
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

USE o_MESH
USE o_PARAM
USE g_PARSUP
USE o_ARRAYS
USE g_ROTATE_grid
use g_comm_auto
IMPLICIT NONE

integer              :: n,j,q, elnodes(3), ed(2), elem, el(2), elnodes_(3)
real(kind=WP)	     :: a(2), b(2), ax, ay, dfactor, lon, lat
real(kind=WP)	     :: deltaX31, deltaX21, deltaY31, deltaY21
real(kind=WP)        :: x(3), y(3), cxx, cxy, cyy, d
real(kind=WP), allocatable :: center_x(:), center_y(:), temp(:) 
real(kind=WP)              :: t0, t1
integer                    :: i, nn, ns

t0=MPI_Wtime()

!real*8,allocatable :: arr2Dglobal(:,:) 
 
 allocate(edge_dxdy(2,myDim_edge2D+eDim_edge2D))
 allocate(edge_cross_dxdy(4,myDim_edge2D+eDim_edge2D))
 allocate(gradient_sca(6,myDim_elem2D))	 
 allocate(gradient_vec(6,myDim_elem2D))
 allocate(metric_factor(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(elem_cos(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(coriolis(myDim_elem2D))
 allocate(coriolis_node(myDim_nod2D+eDim_nod2D))
 allocate(geo_coord_nod2D(2,myDim_nod2D+eDim_nod2D))
 
 allocate(center_x(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(center_y(myDim_elem2D+eDim_elem2D+eXDim_elem2D)) 
 

 ! ============
 ! coriolis
 ! ============
 DO n=1,myDim_nod2D+eDim_nod2D 
 call r2g(lon, lat, coord_nod2D(1,n), coord_nod2D(2,n))
 coriolis_node(n)=2*omega*sin(lat)	 
 END DO
 DO n=1,myDim_nod2D+eDim_nod2D 
 call r2g(lon, lat, coord_nod2D(1,n), coord_nod2D(2,n))
 geo_coord_nod2D(1,n)=lon
 geo_coord_nod2D(2,n)=lat	 
 END DO
 
 DO n=1,myDim_elem2D 
 call elem_center(n, ax, ay)
 call r2g(lon, lat, ax, ay)
 coriolis(n)=2*omega*sin(lat)	 
 END DO
 
 
 if(fplane) then 
 coriolis=2*omega*0.71
 end if
 
 ! ============
 ! cos on elements + metric factor (tan/R_earth) 
 ! ============
 DO n=1,myDim_elem2D
 call elem_center(n, ax, ay)
 center_x(n)=ax
 center_y(n)=ay
 elem_cos(n)=cos(ay)
 metric_factor=tan(ay)/r_earth
 END DO

 call exchange_elem(metric_factor)
 call exchange_elem(elem_cos)
 call exchange_elem(center_x)
 call exchange_elem(center_y)  
 if (cartesian) then
 elem_cos=1.0_WP
 metric_factor=0.0_WP
 end if
 
 ! ===========
 ! Distances along the edge
 ! We need them in radian measure!
 ! ===========
 DO n=1, myDim_edge2D+eDim_edge2D
 ed=edges(:,n)
 a=coord_nod2D(:,ed(2))-coord_nod2D(:, ed(1))
 if(a(1)>cyclic_length/2) a(1)=a(1)-cyclic_length
 if(a(1)<-cyclic_length/2) a(1)=a(1)+cyclic_length
      !a(1)=a(1)*aux_cos_edge(n)
      !a=a*r_earth
 edge_dxdy(:,n)=a
 END DO

 ! ===========
 ! Cross-distances for the edge
 ! They are in physical measure!
 ! ===========
 DO n=1, myDim_edge2D+eDim_edge2D    !!! Difference to FESOM2:
    ! Since we know elem centers on the extended halo of elements
    ! the computations can be carried out for all edges (owned and
    ! the halo ones). 
    ed=edges(:,n)
    el=edge_tri(:,n)
    call edge_center(ed(1), ed(2), a(1), a(2))
    b(1)=center_x(el(1))
    b(2)=center_y(el(1))
    b=b-a
    call trim_cyclic(b(1))
    b(1)=b(1)*elem_cos(el(1))
    b=b*r_earth
    edge_cross_dxdy(1:2,n)=b(1:2)
 
   if(el(2)>0) then
    b(1)=center_x(el(2)) 
    b(2)=center_y(el(2))
    b=b-a
    call trim_cyclic(b(1))
    b(1)=b(1)*elem_cos(el(2))
    b=b*r_earth
    edge_cross_dxdy(3:4,n)=b(1:2)
   else
    edge_cross_dxdy(3:4,n)=0.0_WP
   end if
 END DO

 ! ==========================
 ! Derivatives of scalar quantities
 ! ==========================

DO elem=1, myDim_elem2D
   elnodes=elem2D_nodes(:,elem)
   
   deltaX31=coord_nod2D(1,elnodes(3))-coord_nod2D(1,elnodes(1))
   if(deltaX31>cyclic_length/2) deltaX31=deltaX31-cyclic_length
   if(deltaX31<-cyclic_length/2) deltaX31=deltaX31+cyclic_length
   deltaX31=elem_cos(elem)*deltaX31
   
   deltaX21=coord_nod2D(1,elnodes(2))-coord_nod2D(1,elnodes(1))
   if(deltaX21>cyclic_length/2) deltaX21=deltaX21-cyclic_length
   if(deltaX21<-cyclic_length/2) deltaX21=deltaX21+cyclic_length
   deltaX21=elem_cos(elem)*deltaX21
   
   deltaY31=coord_nod2D(2,elnodes(3))-coord_nod2D(2,elnodes(1))
   deltaY21=coord_nod2D(2,elnodes(2))-coord_nod2D(2,elnodes(1))
   
   dfactor=-0.5_8*r_earth/elem_area(elem)
   gradient_sca(1,elem)=(-deltaY31+deltaY21)*dfactor
   gradient_sca(2,elem)=deltaY31*dfactor
   gradient_sca(3,elem)=-deltaY21*dfactor
   
   gradient_sca(4,elem)=(deltaX31-deltaX21)*dfactor
   gradient_sca(5,elem)=-deltaX31*dfactor
   gradient_sca(6,elem)=deltaX21*dfactor
END DO

 ! ==========================
 ! Derivatives of vector quantities
 ! Least squares interpolation is used
 ! ==========================

   DO elem=1,myDim_elem2D
             !elnodes=elem2D_nodes(:,elem)
      a(1)=center_x(elem)
      a(2)=center_y(elem)
      DO j=1,3
      el(1)=elem_neighbors(j,elem)
      if (el(1)>0) then
             !elnodes_=elem2D_nodes(:,el(1))
      b(1)=center_x(el(1))
      b(2)=center_y(el(1))
      x(j)=b(1)-a(1)
      if(x(j)>cyclic_length/2) x(j)=x(j)-cyclic_length
      if(x(j)<-cyclic_length/2) x(j)=x(j)+cyclic_length
      y(j)=b(2)-a(2)
      else
      ! Virtual element center is taken
      ed=edges(:,elem_edges(j,elem))
      call edge_center(ed(1), ed(2), b(1), b(2))
      x(j)=(b(1)-a(1))
      if(x(j)>cyclic_length/2)   x(j)=x(j)-cyclic_length
      if(x(j)<-cyclic_length/2)  x(j)=x(j)+cyclic_length
      x(j)=2*x(j)
      y(j)=2*(b(2)-a(2))
      end if
      END DO
      x=x*elem_cos(elem)*r_earth
      y=y*r_earth
      cxx=sum(x**2)
      cxy=sum(x*y)
      cyy=sum(y**2)
      d=cxy*cxy-cxx*cyy
	  ! coefficients to compute gradients of velocity
      gradient_vec(1:3,elem)=(cxy*y-cyy*x)/d
      gradient_vec(4:6,elem)=(cxy*x-cxx*y)/d
    END DO
    deallocate(center_y, center_x)


#if defined (__oasis)
  nn=0
  ns=0  
  allocate(lump2d_north(myDim_nod2D), lump2d_south(myDim_nod2D))
  lump2d_north=0.
  lump2d_south=0.
  do i=1, myDim_nod2D
     if (geo_coord_nod2D(2, i) > 0) then
        nn=nn+1
        lump2d_north(i)=area(1, i)
     else
        ns=ns+1     
        lump2d_south(i)=area(1, i)
     end if	   
  end do   

  if (nn>0) allocate(ind_north(nn))
  if (ns>0) allocate(ind_south(ns))
  ns=0
  nn=0
  do i=1, myDim_nod2D
     if (geo_coord_nod2D(2, i) > 0) then
        nn=nn+1
	ind_north(nn)=i
     else
        ns=ns+1
	ind_south(ns)=i	
     end if	     
  end do     
#endif 

    t1=MPI_Wtime()
    if (mype==0) then
       write(*,*) 'mesh_auxiliary_arrays finished in ', t1-t0, ' seconds'
       write(*,*) '========================='
    endif
END SUBROUTINE mesh_auxiliary_arrays

!===================================================================

SUBROUTINE check_mesh_consistency
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_ROTATE_GRID
  use g_comm_auto
IMPLICIT NONE
! Collects auxilliary information on the mesh
! Allocated and filled in are:
! elem_area(myDim_elem2D)
! area(nl, myDim_nod2D)


integer              :: nz, n, elem , elnodes(3)
real(kind=8)	     :: vol_n(nl), vol_e(nl), aux(nl)

   vol_n=0.
   vol_e=0.

   aux=0.
   do n=1, myDim_nod2D
      do nz=1, nlevels_nod2D(n)-1
         aux(nz)=aux(nz)+area(nz, n)
      end do
   end do
   call MPI_AllREDUCE(aux, vol_n, nl, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

   aux=0.
   do elem=1, myDim_elem2D
      elnodes=elem2D_nodes(:, elem)
      if (elnodes(1) > myDim_nod2D) CYCLE
      do nz=1, nlevels(elem)         
         aux(nz)=aux(nz)+elem_area(elem)
      end do
   end do
   call MPI_AllREDUCE(aux, vol_e, nl, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

if (mype==0) then
write(*,*) '***start level area_test***'
do nz=1, nl
   write(*,*) vol_n(nz), vol_e(nz)
end do
write(*,*) '***end level area_test***'
end if

!call par_ex
!stop
END SUBROUTINE check_mesh_consistency
!==================================================================
subroutine trim_cyclic(b)
use o_PARAM
use g_config
implicit none
real(kind=8) :: b
 if(b> cyclic_length/2.0_WP) b=b-cyclic_length
 if(b<-cyclic_length/2.0_WP) b=b+cyclic_length
end subroutine trim_cyclic
!===================================================================

