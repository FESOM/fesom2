module com_global2local_interface
  interface
    subroutine com_global2local(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module

!=============================================================================
SUBROUTINE com_global2local(mesh)
USE g_PARSUP
use MOD_MESH
IMPLICIT NONE

type(t_mesh), intent(in)           , target :: mesh
INTEGER                            :: n, m
INTEGER, ALLOCATABLE, DIMENSION(:) :: temp

#include "associate_mesh_ini.h"

allocate(temp(max(nod2D, elem2D))) 
! =========
! nodes
! =========
  ! Replace global numbering with a local one
  temp(:)=0
  DO n=1, myDim_nod2D
     temp(myList_nod2D(n))=n
  END DO

  DO n=1, eDim_nod2D
     com_nod2D%rlist(n)=myDim_nod2D+n
  END DO

  DO n=1, com_nod2D%sptr(com_nod2D%sPEnum+1)-1
     m=com_nod2D%slist(n)
     com_nod2D%slist(n)=temp(m)
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
  DO n=1, myDim_elem2D
     temp(myList_elem2D(n))=n
  END DO

  DO n=1, eDim_elem2D
     temp(com_elem2D%rlist(n))=myDim_elem2D+n
     com_elem2D%rlist(n)=myDim_elem2D+n
  END DO

  DO n=1, com_elem2D%sptr(com_elem2D%sPEnum+1)-1
     m=com_elem2D%slist(n)
     com_elem2D%slist(n)=temp(m)
  END DO

!!$  DO n=1, com_elem2D%rptr(com_elem2D%rPEnum+1)-1
!!$         m=com_elem2D%rlist(n)
!!$	 com_elem2D%rlist(n)=temp(m)
!!$  END DO
! =========
! elements (extra needed for MUSCL advection)
! =========
  
  m = 0
  DO n=1, com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1
     if (temp(com_elem2D_full%rlist(n))>0) then
        com_elem2D_full%rlist(n)=temp(com_elem2D_full%rlist(n))
     else
        m = m+1
        com_elem2D_full%rlist(n)=myDim_elem2D+eDim_elem2D+m
     endif
  END DO
  
  DO n=1, com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)-1
     m=com_elem2D_full%slist(n)
     com_elem2D_full%slist(n)=temp(m)
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
SUBROUTINE save_dist_mesh(mesh)
  USE g_CONFIG
  USE MOD_MESH
  USE o_ARRAYS
  USE g_PARSUP 
  use com_global2local_interface
  IMPLICIT NONE

  type(t_mesh), intent(in)           , target :: mesh
  Integer        n, m, q, q2, counter, fileID, nend, nini,ed(2)
  character*10   mype_string,npes_string
  character(MAX_PATH) file_name
  character(MAX_PATH) dist_mesh_dir
  integer, allocatable, dimension(:)  :: temp, ncount
  integer   n1, n2, flag, eledges(4)

#include  "associate_mesh_ini.h"

!!$  allocate(temp(nod2D))  ! serves for mapping
!!$  allocate(ncount(npes+1))
  write(mype_string,'(i5.5)') mype  
  write(npes_string,"(I10)") npes
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
  fileID=103+mype !skip unit range 100--102 
  ! =============================   
  ! lists of owned nodes and elements
  ! =============================
  open(fileID, file=file_name)
  write(fileID,*) mype
  write(fileID,*) myDim_nod2D
  write(fileID,*) eDim_nod2D 	 
  write(fileID,*) myList_nod2D(1:myDim_nod2D), com_nod2D%rlist(1:eDim_nod2D)

  ! Create a list for eXDim_elem2D containing elements not in eDim_elem2D
  allocate(temp(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1))
  counter = 0
  DO n=1, com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1
     m = com_elem2D_full%rlist(n)
     if (any(com_elem2D%rlist(1:eDim_elem2D) == m)) cycle
     counter = counter + 1
     temp(counter) = m
  end do
  eXDim_elem2D=counter

  write(fileID,*) myDim_elem2D
  write(fileID,*) eDim_elem2D
  write(fileID,*) eXDim_elem2D	 	 
  write(fileID,*) myList_elem2D(1:myDim_elem2D), com_elem2D%rlist(1:eDim_elem2D), temp(1:eXDim_elem2D)
  deallocate(temp)

  allocate(myList_edge2D(4*myDim_elem2D))
  counter = 0
  do n=1, edge2D
     do q=1,2 
        if (part(edges(q,n))==mype) then
           counter=counter+1
           myList_edge2D(counter)=n
           exit
        end if
     end do
  end do
  myDim_edge2D=counter   ! It is the total number of my edges 
  Do n=1,myDim_elem2D
     eledges=elem_edges(:,myList_elem2D(n))
     q2 = merge(3,4,eledges(1) == eledges(4))
     DO q=1,q2
        if((part(edges(1,eledges(q))).ne.mype).and.(part(edges(2,eledges(q))).ne.mype) &
             .and. all(myList_edge2D(myDim_edge2D:counter) /= eledges(q))) then
           counter=counter+1 
           myList_edge2D(counter)=eledges(q) 
        end if
     END DO
  END DO
  eDim_edge2D=counter-myDim_edge2D
  if (myDim_edge2D+eDim_edge2D>4*myDim_elem2D) then
     write (*,*) 'Edge array too short in save_dist_mesh!'
     stop
  endif

  write(fileID,*) myDim_edge2D
  write(fileID,*) eDim_edge2D 	 
  write(fileID,*) myList_edge2D(1:myDim_edge2D +eDim_edge2D)
  deallocate(myList_edge2D)
  close(fileID)       


  ! =========================  
  ! communication information
  ! ========================= 
  call com_global2local(mesh)   ! Do not call this subroutine earlier, global numbering is needed!
  file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out' 
  fileID=103+mype  !skip unit range 100--102 
  open(fileID, file=file_name)
  write(fileID,*) mype
  write(fileID,*) com_nod2D%rPEnum
  write(fileID,*) com_nod2D%rPE(1:com_nod2D%rPEnum)
  write(fileID,*) com_nod2D%rptr(1:com_nod2D%rPEnum+1)
  write(fileID,*) com_nod2D%rlist
  write(fileID,*) com_nod2D%sPEnum
  write(fileID,*) com_nod2D%sPE(1:com_nod2D%sPEnum)
  write(fileID,*) com_nod2D%sptr(1:com_nod2D%sPEnum+1)
  write(fileID,*) com_nod2D%slist
  deallocate(myList_nod2D)
!!$  deallocate(com_nod2D%rPE)
!!$  deallocate(com_nod2D%rptr)
  deallocate(com_nod2D%rlist)
!!$  deallocate(com_nod2D%sPE)
!!$  deallocate(com_nod2D%sptr)
  deallocate(com_nod2D%slist)

  write(fileID,*) com_elem2D%rPEnum
  write(fileID,*) com_elem2D%rPE(1:com_elem2D%rPEnum)
  write(fileID,*) com_elem2D%rptr(1:com_elem2D%rPEnum+1)
  write(fileID,*) com_elem2D%rlist
  write(fileID,*) com_elem2D%sPEnum
  write(fileID,*) com_elem2D%sPE(1:com_elem2D%sPEnum)
  write(fileID,*) com_elem2D%sptr(1:com_elem2D%sPEnum+1)
  write(fileID,*) com_elem2D%slist
  deallocate(myList_elem2D)
!!$  deallocate(com_elem2D%rPE)
!!$  deallocate(com_elem2D%rptr)
  deallocate(com_elem2D%rlist)
!!$  deallocate(com_elem2D%sPE)
!!$  deallocate(com_elem2D%sptr)
  deallocate(com_elem2D%slist)

  write(fileID,*) com_elem2D_full%rPEnum
  write(fileID,*) com_elem2D_full%rPE(1:com_elem2D_full%rPEnum)
  write(fileID,*) com_elem2D_full%rptr(1:com_elem2D_full%rPEnum+1)
  write(fileID,*) com_elem2D_full%rlist
  write(fileID,*) com_elem2D_full%sPEnum
  write(fileID,*) com_elem2D_full%sPE(1:com_elem2D_full%sPEnum)
  write(fileID,*) com_elem2D_full%sptr(1:com_elem2D_full%sPEnum+1)
  write(fileID,*) com_elem2D_full%slist
!!$  deallocate(com_elem2D_full%rPE)
!!$  deallocate(com_elem2D_full%rptr)
  deallocate(com_elem2D_full%rlist)
!!$  deallocate(com_elem2D_full%sPE)
!!$  deallocate(com_elem2D_full%sptr)
  deallocate(com_elem2D_full%slist)
  close(fileID)
  ! ================================
  ! mapping ( PE contiguous 2D numbering) 	 
  ! ================================  
  if(mype==0) then
     file_name=trim(dist_mesh_dir)//'rpart.out'  
     fileID=103+mype !skip unit range 100--102 
     open(fileID, file=file_name)

     allocate(temp(nod2D), ncount(npes+1)) ! mapping 
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

     deallocate(ncount, temp)
     close(fileID)
  end if
END subroutine  save_dist_mesh
