module oce_local_module
    USE MOD_MESH
    USE MOD_PARTIT
    USE g_CONFIG
    USE o_ARRAYS

    implicit none

    private
    public :: com_global2local, save_dist_mesh

contains

!=============================================================================
SUBROUTINE com_global2local(partit, mesh)
IMPLICIT NONE

type(t_mesh),   intent(in),    target :: mesh
type(t_partit), intent(inout), target :: partit

INTEGER                            :: n, m
INTEGER, ALLOCATABLE, DIMENSION(:) :: temp


allocate(temp(max(mesh%nod2D, mesh%elem2D))) 
! =========
! nodes
! =========
  ! Replace global numbering with a local one
  temp(:)=0
  DO n=1, partit%myDim_nod2D
     temp(partit%myList_nod2D(n))=n
  END DO

  DO n=1, partit%eDim_nod2D
     partit%com_nod2D%rlist(n)=partit%myDim_nod2D+n
  END DO

  DO n=1, partit%com_nod2D%sptr(partit%com_nod2D%sPEnum+1)-1
     m=partit%com_nod2D%slist(n)
     partit%com_nod2D%slist(n)=temp(m)
  END DO

!!$  DO n=1, com_nod2D%rptr(com_nod2D%rPEnum+1)-1
!!$         m=com_nod2D%rlist(n)
!!$	 com_nod2D%rlist(n)=temp(m)
!!$  END DO
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
  temp(:)=0
  DO n=1, partit%myDim_elem2D
     temp(partit%myList_elem2D(n))=n
  END DO

  DO n=1, partit%eDim_elem2D
     temp(partit%com_elem2D%rlist(n))=partit%myDim_elem2D+n
     partit%com_elem2D%rlist(n)=partit%myDim_elem2D+n
  END DO

  DO n=1, partit%com_elem2D%sptr(partit%com_elem2D%sPEnum+1)-1
     m=partit%com_elem2D%slist(n)
     partit%com_elem2D%slist(n)=temp(m)
  END DO

!!$  DO n=1, com_elem2D%rptr(com_elem2D%rPEnum+1)-1
!!$         m=com_elem2D%rlist(n)
!!$	 com_elem2D%rlist(n)=temp(m)
!!$  END DO
! =========
! elements (extra needed for MUSCL advection)
! =========
  
  m = 0
  DO n=1, partit%com_elem2D_full%rptr(partit%com_elem2D_full%rPEnum+1)-1
     if (temp(partit%com_elem2D_full%rlist(n))>0) then
        partit%com_elem2D_full%rlist(n)=temp(partit%com_elem2D_full%rlist(n))
     else
        m = m+1
        partit%com_elem2D_full%rlist(n)=partit%myDim_elem2D+partit%eDim_elem2D+m
     endif
  END DO
  
  DO n=1, partit%com_elem2D_full%sptr(partit%com_elem2D_full%sPEnum+1)-1
     m=partit%com_elem2D_full%slist(n)
     partit%com_elem2D_full%slist(n)=temp(m)
  END DO

!!$  DO n=1, com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1
!!$         m=com_elem2D_full%rlist(n)
!!$	 com_elem2D_full%rlist(n)=temp(m)
!!$  END DO
!!$  
	 ! com_elem2%rlist should be 
	 ! myDim_elem2D+1:myDim_elem2D+eDim_elem2D
deallocate(temp)
END SUBROUTINE com_global2local       	  
!=============================================================================
SUBROUTINE save_dist_mesh(partit, mesh)
  IMPLICIT NONE

  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  Integer        n, m, q, q2, counter, fileID, nend, nini,ed(2)
  character*10   mype_string,npes_string
  character(MAX_PATH) file_name
  character(MAX_PATH) dist_mesh_dir
  integer, allocatable, dimension(:)  :: temp, ncount
  integer   n1, n2, flag, eledges(4)


  write(mype_string,'(i5.5)') partit%mype  
  write(npes_string,"(I10)")  partit%npes
  dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'

  ! ==============================
  ! rank partitioning:
  ! ==============================

!!$  if(mype==0) then
!!$     file_name=trim(dist_mesh_dir)//'rpart.out'  
!!$     fileID=103+mype  !skip unit range 100--102 
!!$     open(fileID, file=file_name)
!!$     ncount=0;
!!$     DO n=1, nod2D
!!$        m=part(n);
!!$        ncount(m+1)=ncount(m+1)+1
!!$     END DO
!!$     write(fileID,*) npes
!!$     write(fileID,*) ncount(1:npes)
!!$     close(fileID)
!!$  end if
!!$

  file_name=trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
  fileID=103+partit%mype !skip unit range 100--102 
  ! =============================   
  ! lists of owned nodes and elements
  ! =============================
  open(fileID, file=file_name)
  write(fileID,*) partit%mype
  write(fileID,*) partit%myDim_nod2D
  write(fileID,*) partit%eDim_nod2D 	 
  write(fileID,*) partit%myList_nod2D(1:partit%myDim_nod2D), partit%com_nod2D%rlist(1:partit%eDim_nod2D)

  ! Create a list for eXDim_elem2D containing elements not in eDim_elem2D
  allocate(temp(partit%com_elem2D_full%rptr(partit%com_elem2D_full%rPEnum+1)-1))
  counter = 0
  DO n=1, partit%com_elem2D_full%rptr(partit%com_elem2D_full%rPEnum+1)-1
     m = partit%com_elem2D_full%rlist(n)
     if (any(partit%com_elem2D%rlist(1:partit%eDim_elem2D) == m)) cycle
     counter = counter + 1
     temp(counter) = m
  end do
  partit%eXDim_elem2D=counter

  write(fileID,*) partit%myDim_elem2D
  write(fileID,*) partit%eDim_elem2D
  write(fileID,*) partit%eXDim_elem2D	 	 
  write(fileID,*) partit%myList_elem2D(1:partit%myDim_elem2D), partit%com_elem2D%rlist(1:partit%eDim_elem2D), temp(1:partit%eXDim_elem2D)
  deallocate(temp)

  allocate(partit%myList_edge2D(4*partit%myDim_elem2D))
  counter = 0
  do n=1, mesh%edge2D
     do q=1,2 
        if (partit%part(mesh%edges(q,n))==partit%mype) then
           counter=counter+1
           partit%myList_edge2D(counter)=n
           exit
        end if
     end do
  end do
  partit%myDim_edge2D=counter   ! It is the total number of my edges 
  Do n=1,partit%myDim_elem2D
     eledges=mesh%elem_edges(:,partit%myList_elem2D(n))
     q2 = merge(3,4,eledges(1) == eledges(4))
     DO q=1,q2
        if((partit%part(mesh%edges(1,eledges(q))).ne.partit%mype).and.(partit%part(mesh%edges(2,eledges(q))).ne.partit%mype) &
             .and. all(partit%myList_edge2D(partit%myDim_edge2D:counter) /= eledges(q))) then
           counter=counter+1 
           partit%myList_edge2D(counter)=eledges(q) 
        end if
     END DO
  END DO
  partit%eDim_edge2D=counter-partit%myDim_edge2D
  if (partit%myDim_edge2D+partit%eDim_edge2D>4*partit%myDim_elem2D) then
     write (*,*) 'Edge array too short in save_dist_mesh!'
     stop
  endif

  write(fileID,*) partit%myDim_edge2D
  write(fileID,*) partit%eDim_edge2D 	 
  write(fileID,*) partit%myList_edge2D(1:partit%myDim_edge2D +partit%eDim_edge2D)
  deallocate(partit%myList_edge2D)
  close(fileID)       


  ! =========================  
  ! communication information
  ! ========================= 
  call com_global2local(partit, mesh)   ! Do not call this subroutine earlier, global numbering is needed!
  file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out' 
  fileID=103+partit%mype  !skip unit range 100--102 
  open(fileID, file=file_name)
  write(fileID,*) partit%mype
  write(fileID,*) partit%com_nod2D%rPEnum
  write(fileID,*) partit%com_nod2D%rPE(1:partit%com_nod2D%rPEnum)
  write(fileID,*) partit%com_nod2D%rptr(1:partit%com_nod2D%rPEnum+1)
  write(fileID,*) partit%com_nod2D%rlist
  write(fileID,*) partit%com_nod2D%sPEnum
  write(fileID,*) partit%com_nod2D%sPE(1:partit%com_nod2D%sPEnum)
  write(fileID,*) partit%com_nod2D%sptr(1:partit%com_nod2D%sPEnum+1)
  write(fileID,*) partit%com_nod2D%slist
  deallocate(partit%myList_nod2D)
!!$  deallocate(com_nod2D%rPE)
!!$  deallocate(com_nod2D%rptr)
  deallocate(partit%com_nod2D%rlist)
!!$  deallocate(com_nod2D%sPE)
!!$  deallocate(com_nod2D%sptr)
  deallocate(partit%com_nod2D%slist)

  write(fileID,*) partit%com_elem2D%rPEnum
  write(fileID,*) partit%com_elem2D%rPE(1:partit%com_elem2D%rPEnum)
  write(fileID,*) partit%com_elem2D%rptr(1:partit%com_elem2D%rPEnum+1)
  write(fileID,*) partit%com_elem2D%rlist
  write(fileID,*) partit%com_elem2D%sPEnum
  write(fileID,*) partit%com_elem2D%sPE(1:partit%com_elem2D%sPEnum)
  write(fileID,*) partit%com_elem2D%sptr(1:partit%com_elem2D%sPEnum+1)
  write(fileID,*) partit%com_elem2D%slist
  deallocate(partit%myList_elem2D)
!!$  deallocate(com_elem2D%rPE)
!!$  deallocate(com_elem2D%rptr)
  deallocate(partit%com_elem2D%rlist)
!!$  deallocate(com_elem2D%sPE)
!!$  deallocate(com_elem2D%sptr)
  deallocate(partit%com_elem2D%slist)

  write(fileID,*) partit%com_elem2D_full%rPEnum
  write(fileID,*) partit%com_elem2D_full%rPE(1:partit%com_elem2D_full%rPEnum)
  write(fileID,*) partit%com_elem2D_full%rptr(1:partit%com_elem2D_full%rPEnum+1)
  write(fileID,*) partit%com_elem2D_full%rlist
  write(fileID,*) partit%com_elem2D_full%sPEnum
  write(fileID,*) partit%com_elem2D_full%sPE(1:partit%com_elem2D_full%sPEnum)
  write(fileID,*) partit%com_elem2D_full%sptr(1:partit%com_elem2D_full%sPEnum+1)
  write(fileID,*) partit%com_elem2D_full%slist
!!$  deallocate(com_elem2D_full%rPE)
!!$  deallocate(com_elem2D_full%rptr)
  deallocate(partit%com_elem2D_full%rlist)
!!$  deallocate(com_elem2D_full%sPE)
!!$  deallocate(com_elem2D_full%sptr)
  deallocate(partit%com_elem2D_full%slist)
  close(fileID)
  ! ================================
  ! mapping ( PE contiguous 2D numbering) 	 
  ! ================================  
  if(partit%mype==0) then
     file_name=trim(dist_mesh_dir)//'rpart.out'  
     fileID=103+partit%mype !skip unit range 100--102 
     open(fileID, file=file_name)

     allocate(temp(mesh%nod2D), ncount(partit%npes+1)) ! mapping 
     ncount=0
     DO n=1, mesh%nod2D
        m=partit%part(n)
        ncount(m+1)=ncount(m+1)+1
     END DO
     write(fileID,*) partit%npes
     write(fileID,*) ncount(1:partit%npes)

     ncount=0
     DO n=1, mesh%nod2D
        m=partit%part(n)
        ncount(m+2)=ncount(m+2)+1
        temp(n)=ncount(m+2)
     END DO
     ncount(1)=1
     Do n=2,partit%npes+1	  
        ncount(n)=ncount(n)+ncount(n-1)
     end do
     ! Now count == part in range partitioning   

     do n=1,mesh%nod2D
        temp(n)=temp(n)+ncount(partit%part(n)+1)-1
        write(fileID,*) temp(n)  
     end do

     deallocate(ncount, temp)
     close(fileID)
  end if
END subroutine  save_dist_mesh

end module oce_local_module
