!=============================================================================
SUBROUTINE com_global2local
USE g_PARSUP
USE o_MESH
IMPLICIT NONE
INTEGER                            :: n, m
INTEGER, ALLOCATABLE, DIMENSION(:) :: temp
allocate(temp(edge2D)) 
! =========
! nodes
! =========
  ! Replace global numbering with a local one
  temp(1:edge2D)=0
  DO n=1, myDim_nod2D+eDim_nod2D
  temp(myList_nod2D(n))=n
  END DO
  DO n=1, com_nod2D%sptr(com_nod2D%sPEnum+1)-1
         m=com_nod2D%slist(n)
	 com_nod2D%slist(n)=temp(m)
  END DO

  DO n=1, com_nod2D%rptr(com_nod2D%rPEnum+1)-1
         m=com_nod2D%rlist(n)
	 com_nod2D%rlist(n)=temp(m)
  END DO
	 ! com_nod2D%rlist should be  
	 ! myDim_nod2D+1:myDim_nod2D+eDim_nod2D

! =========
! edges
! =========
!!$  ! Replace global numbering with a local one
!!$  temp(1:edge2D)=0
!!$  DO n=1, myDim_edge2D+eDim_edge2D
!!$  temp(myList_edge2D(n))=n
!!$  END DO
!!$  DO n=1, com_edge2D%sptr(com_edge2D%sPEnum+1)-1
!!$         m=com_edge2D%slist(n)
!!$	 com_edge2D%slist(n)=temp(m)
!!$  END DO
!!$
!!$  DO n=1, com_edge2D%rptr(com_edge2D%rPEnum+1)-1
!!$         m=com_edge2D%rlist(n)
!!$	 com_edge2D%rlist(n)=temp(m)
!!$  END DO
!!$	 ! com_edge2%rlist should be 
!!$	 ! myDim_edge2D+1:myDim_edge2D+eDim_edge2D
! =========
! elements
! =========
  ! Replace global numbering with a local one
  temp(1:edge2D)=0
  DO n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
  temp(myList_elem2D(n))=n
  END DO
  DO n=1, com_elem2D%sptr(com_elem2D%sPEnum+1)-1
         m=com_elem2D%slist(n)
	 com_elem2D%slist(n)=temp(m)
  END DO

  DO n=1, com_elem2D%rptr(com_elem2D%rPEnum+1)-1
         m=com_elem2D%rlist(n)
	 com_elem2D%rlist(n)=temp(m)
  END DO
! =========
! elements (extra needed for MUSCL advection)
! =========
  
  DO n=1, com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)-1
         m=com_elem2D_full%slist(n)
	 com_elem2D_full%slist(n)=temp(m)
  END DO

  DO n=1, com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1
         m=com_elem2D_full%rlist(n)
	 com_elem2D_full%rlist(n)=temp(m)
  END DO
  
	 ! com_elem2%rlist should be 
	 ! myDim_elem2D+1:myDim_elem2D+eDim_elem2D
deallocate(temp)
END SUBROUTINE com_global2local       	  
!=============================================================================
SUBROUTINE save_dist_mesh
  USE g_CONFIG
  USE o_MESH
  USE o_ARRAYS
  USE g_PARSUP 
  IMPLICIT NONE

  Integer        n, m, fileID, nend, nini,ed(2)
  character*10   mype_string,npes_string
  character*200   file_name
  character*200   dist_mesh_dir
  integer, allocatable, dimension(:)  :: temp, ncount
  integer   n1, n2, flag

  allocate(temp(nod2D))  ! serves for mapping
  allocate(ncount(npes+1))
  write(mype_string,'(i5.5)') mype  
  write(npes_string,"(I10)") npes
  dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'

  ! ==============================
  ! rank partitioning:
  ! ==============================

  if(mype==0) then
     file_name=trim(dist_mesh_dir)//'rpart.out'  
     fileID=103+mype !skip unit range 100--102  
     open(fileID, file=file_name)
     ncount=0;
     DO n=1, nod2D
        m=part(n);
        ncount(m+1)=ncount(m+1)+1
     END DO
     write(fileID,*) npes
     write(fileID,*) ncount(1:npes)
     close(fileID)
  end if


  file_name=trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
  fileID=103+mype !skip unit range 100--102 
  ! =============================   
  ! lists of owned nodes and elements
  ! =============================
  open(fileID, file=file_name)
  write(fileID,*) mype
  write(fileID,*) myDim_nod2D
  write(fileID,*) eDim_nod2D 	 
  write(fileID,*) myList_nod2D(1:myDim_nod2D+eDim_nod2D)

  write(fileID,*) myDim_elem2D
  write(fileID,*) eDim_elem2D
  write(fileID,*) eXDim_elem2D	 	 
  write(fileID,*) myList_elem2D(1:myDim_elem2D +eDim_elem2D +eXDim_elem2D)

  write(fileID,*) myDim_edge2D
  write(fileID,*) eDim_edge2D 	 
  write(fileID,*) myList_edge2D(1:myDim_edge2D +eDim_edge2D)
  close(fileID)       


  ! =========================  
  ! communication information
  ! ========================= 
  call com_global2local   
  file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out' 
  fileID=103+mype  !skip unit range 100--102 
  open(fileID, file=file_name)
  write(fileID,*) mype
  write(fileID,*) com_nod2D%rPEnum
  write(fileID,*) com_nod2D%rPE
  write(fileID,*) com_nod2D%rptr
  write(fileID,*) com_nod2D%rlist
  write(fileID,*) com_nod2D%sPEnum
  write(fileID,*) com_nod2D%sPE
  write(fileID,*) com_nod2D%sptr
  write(fileID,*) com_nod2D%slist
  deallocate(com_nod2D%rPE)
  deallocate(com_nod2D%rptr)
  deallocate(com_nod2D%rlist)
  deallocate(com_nod2D%sPE)
  deallocate(com_nod2D%sptr)
  deallocate(com_nod2D%slist)

  write(fileID,*) com_elem2D%rPEnum
  write(fileID,*) com_elem2D%rPE
  write(fileID,*) com_elem2D%rptr
  write(fileID,*) com_elem2D%rlist
  write(fileID,*) com_elem2D%sPEnum
  write(fileID,*) com_elem2D%sPE
  write(fileID,*) com_elem2D%sptr
  write(fileID,*) com_elem2D%slist
  deallocate(com_elem2D%rPE)
  deallocate(com_elem2D%rptr)
  deallocate(com_elem2D%rlist)
  deallocate(com_elem2D%sPE)
  deallocate(com_elem2D%sptr)
  deallocate(com_elem2D%slist)

  write(fileID,*) com_elem2D_full%rPEnum
  write(fileID,*) com_elem2D_full%rPE
  write(fileID,*) com_elem2D_full%rptr
  write(fileID,*) com_elem2D_full%rlist
  write(fileID,*) com_elem2D_full%sPEnum
  write(fileID,*) com_elem2D_full%sPE
  write(fileID,*) com_elem2D_full%sptr
  write(fileID,*) com_elem2D_full%slist
  deallocate(com_elem2D_full%rPE)
  deallocate(com_elem2D_full%rptr)
  deallocate(com_elem2D_full%rlist)
  deallocate(com_elem2D_full%sPE)
  deallocate(com_elem2D_full%sptr)
  deallocate(com_elem2D_full%slist)
  close(fileID)
  ! ================================
  ! mapping ( PE contiguous 2D numbering) 	 
  ! ================================  
  if(mype==0) then
     file_name=trim(dist_mesh_dir)//'rpart.out'  
     fileID=103+mype !skip unit range 100--102 
     open(fileID, file=file_name)
     ncount=0
     DO n=1, nod2D
        m=part(n)
        ncount(m+1)=ncount(m+1)+1
     END DO
     write(fileID,*) npes
     write(fileID,*) ncount(1:npes)

     ncount=0
     DO n=1, nod2D
        m=part(n)
        ncount(m+2)=ncount(m+2)+1
        temp(n)=ncount(m+2)
     END DO
     ncount(1)=1
     Do n=2,npes+1	  
        ncount(n)=ncount(n)+ncount(n-1)
     end do
     ! Now count == part in range partitioning   

     do n=1,nod2D
        temp(n)=temp(n)+ncount(part(n)+1)-1
        write(fileID,*) temp(n)  
     end do
     close(fileID)
  end if

  deallocate(ncount, temp)
!!$  write(*,*) 'Distributed mesh is saved for rank ', mype
END subroutine  save_dist_mesh


