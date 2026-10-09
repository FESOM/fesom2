module com_global2local_interface
  interface
    subroutine com_global2local(partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_mesh),   intent(in),    target :: mesh
      type(t_partit), intent(inout), target :: partit
    end subroutine com_global2local
  end interface
end module com_global2local_interface

!=============================================================================
! Hybrid meshes: sections appended to dist_N, announced by dist_N/format.out
! (version DIST_FORMAT_VERSION). The FESOM2 part of my_list*.out and com_info*.out
! stays as it is; appended are
!   my_list : eXDim_nod2D and the global ids of the double-halo nodes = vertices of
!             the elements 1..myDim+eDim+eXDim_elem2D that are not among the nodes
!             1..myDim+eDim_nod2D (grouped by owner, ascending within);
!   com_info: com_nod2D_full (every halo node 1..eDim+eXDim_nod2D, received from its
!             owner) and com_edge2D (an edge is sent by its unique owner, the owner of
!             edges(1,edge), to every other PE holding it in 1..myDim+eDim_edge2D,
!             also where it lies in that PE's 1..myDim_edge2D).
! Lists are in local numbering (nodes 1..myDim+eDim+eXDim_nod2D, edges = position in
! myList_edge2D); for one peer both sides are ordered by global id.
module dist_ext
  use MOD_MESH
  use MOD_PARTIT
  implicit none
  integer, parameter :: DIST_FORMAT_VERSION = 1

  type dist_com
     integer              :: rPEnum, sPEnum
     integer, allocatable :: rPE(:), rptr(:), rlist(:), sPE(:), sptr(:), slist(:)
  end type dist_com

contains
!=============================================================================
! eXDim_nod2D list (global ids) and com_nod2D_full. Needs the global numbering
! (call before com_global2local).
subroutine nod_full_lists(partit, mesh, xlist, com)
  type(t_mesh),   intent(in)        :: mesh
  type(t_partit), intent(in)        :: partit
  integer, allocatable, intent(out) :: xlist(:)
  type(dist_com), intent(out)       :: com
  integer, allocatable :: loc(:), pes(:,:), pos(:)
  integer              :: num_recv(0:partit%npes-1), num_send(0:partit%npes-1)
  integer              :: myd, ed, nel, n, m, k, i, j, q, q2, el, el2, p, nx
  logical              :: too_small

  myd = partit%myDim_nod2D
  ed  = partit%eDim_nod2D
  allocate(loc(mesh%nod2D), pos(0:partit%npes-1))
  loc = 0
  do k=1, myd
     loc(partit%myList_nod2D(k)) = k
  end do
  do k=1, ed
     loc(partit%com_nod2D%rlist(k)) = myd+k
  end do
  ! vertices of the local elements (owned + com_elem2D_full receive list) not yet local
  nel = partit%com_elem2D_full%rptr(partit%com_elem2D_full%rPEnum+1)-1
  do k=1, partit%myDim_elem2D+nel
     if (k <= partit%myDim_elem2D) then
        el = partit%myList_elem2D(k)
     else
        el = partit%com_elem2D_full%rlist(k-partit%myDim_elem2D)
     end if
     do q=1,4
        n = mesh%elem2D_nodes(q,el)
        if (loc(n) == 0) loc(n) = -1
     end do
  end do
  ! double-halo nodes grouped by owner, ascending within
  num_recv = 0
  do n=1, mesh%nod2D
     if (loc(n) == -1) num_recv(partit%part(n)) = num_recv(partit%part(n))+1
  end do
  nx = sum(num_recv)
  allocate(xlist(nx))
  call first_pos(num_recv, pos)
  do n=1, mesh%nod2D
     if (loc(n) /= -1) cycle
     p = partit%part(n)
     xlist(pos(p)) = n
     pos(p) = pos(p)+1
  end do
  do k=1, nx
     loc(xlist(k)) = myd+ed+k
  end do

  ! receive: every halo node (eDim + eXDim) from its owner
  num_recv = 0
  do n=1, mesh%nod2D
     if (loc(n) > myd) num_recv(partit%part(n)) = num_recv(partit%part(n))+1
  end do
  call set_pe_ptr(num_recv, com%rPEnum, com%rPE, com%rptr)
  allocate(com%rlist(com%rptr(com%rPEnum+1)-1))
  call first_pos(num_recv, pos)
  do n=1, mesh%nod2D
     if (loc(n) <= myd) cycle
     p = partit%part(n)
     com%rlist(pos(p)) = loc(n)
     pos(p) = pos(p)+1
  end do

  ! send: an owned node n is local to PE p if it is a vertex of an element el that
  ! shares a vertex m with an element el2 having a vertex owned by p
  allocate(pes(MAX_LAENDERECK, myd))
  pes = -1
  num_send = 0
  too_small = .false.
  do k=1, myd
     n = partit%myList_nod2D(k)
     do i=1, mesh%nod_in_elem2D_num(n)
        el = mesh%nod_in_elem2D(i,n)
        do q=1,4
           m = mesh%elem2D_nodes(q,el)
           do j=1, mesh%nod_in_elem2D_num(m)
              el2 = mesh%nod_in_elem2D(j,m)
              do q2=1,4
                 p = partit%part(mesh%elem2D_nodes(q2,el2))
                 if (p /= partit%mype) call add_pe(pes(:,k), p, num_send, too_small)
              end do
           end do
        end do
     end do
  end do
  call set_send(pes, num_send, too_small, com)
  deallocate(loc, pes, pos)
end subroutine nod_full_lists
!=============================================================================
! com_edge2D from myList_edge2D(1:nedges) (global ids, myDim_edge2D + eDim_edge2D).
subroutine edge_lists(partit, mesh, nedges, com)
  type(t_mesh),   intent(in)  :: mesh
  type(t_partit), intent(in)  :: partit
  integer,        intent(in)  :: nedges
  type(dist_com), intent(out) :: com
  integer, allocatable :: loc(:), pes(:,:), pos(:)
  integer              :: num_recv(0:partit%npes-1), num_send(0:partit%npes-1)
  integer              :: e, k, i, q, n, el, p
  logical              :: too_small

  allocate(loc(mesh%edge2D), pos(0:partit%npes-1))
  loc = 0
  do k=1, nedges
     loc(partit%myList_edge2D(k)) = k
  end do

  ! receive: every edge I hold whose unique owner is another PE
  num_recv = 0
  do e=1, mesh%edge2D
     if (loc(e) == 0) cycle
     p = partit%part(mesh%edges(1,e))
     if (p /= partit%mype) num_recv(p) = num_recv(p)+1
  end do
  call set_pe_ptr(num_recv, com%rPEnum, com%rPE, com%rptr)
  allocate(com%rlist(com%rptr(com%rPEnum+1)-1))
  call first_pos(num_recv, pos)
  do e=1, mesh%edge2D
     if (loc(e) == 0) cycle
     p = partit%part(mesh%edges(1,e))
     if (p == partit%mype) cycle
     com%rlist(pos(p)) = loc(e)
     pos(p) = pos(p)+1
  end do

  ! send: a uniquely owned edge goes to the PEs of the vertices of its elements (a PE
  ! holds an edge iff the edge is a side of an element with a vertex owned by it)
  allocate(pes(MAX_LAENDERECK, nedges))
  pes = -1
  num_send = 0
  too_small = .false.
  do k=1, nedges
     e = partit%myList_edge2D(k)
     n = mesh%edges(1,e)
     if (partit%part(n) /= partit%mype) cycle
     do i=1, mesh%nod_in_elem2D_num(n)
        el = mesh%nod_in_elem2D(i,n)
        if (all(mesh%elem_edges(:,el) /= e)) cycle
        do q=1,4
           p = partit%part(mesh%elem2D_nodes(q,el))
           if (p /= partit%mype) call add_pe(pes(:,k), p, num_send, too_small)
        end do
     end do
  end do
  ! uniquely owned edges lie in 1..myDim_edge2D, which is ascending in global id
  call set_send(pes, num_send, too_small, com)
  deallocate(loc, pes, pos)
end subroutine edge_lists
!=============================================================================
! The send side of nod_full_lists / edge_lists: pes(:,k) = the PEs local item k goes to
! (-1 = none), num_send(p) = the number of items for PE p; slist = for each PE
! (ascending) its items in local order
subroutine set_send(pes, num_send, too_small, com)
  integer,        intent(in)    :: pes(:,:), num_send(0:)
  logical,        intent(in)    :: too_small
  type(dist_com), intent(inout) :: com
  integer :: i, k, counter

  if (too_small) then
     print *,'Increase MAX_LAENDERECK in MOD_PARTIT.F90 and recompile'
     stop
  endif
  call set_pe_ptr(num_send, com%sPEnum, com%sPE, com%sptr)
  allocate(com%slist(com%sptr(com%sPEnum+1)-1))
  counter = 0
  do i=1, com%sPEnum
     do k=1, size(pes, 2)
        if (any(pes(:,k) == com%sPE(i))) then
           counter = counter+1
           com%slist(counter) = k
        end if
     end do
  end do
end subroutine set_send
!=============================================================================
! PEs with num(p) > 0 (ascending) and the pointers into their lists; stops on more than
! MAX_NEIGHBOR_PARTITIONS PEs (as gen_comm.F90 for com_nod2D/com_elem2D/com_elem2D_full)
subroutine set_pe_ptr(num, npe, pe, ptr)
  integer,              intent(in)  :: num(0:)
  integer,              intent(out) :: npe
  integer, allocatable, intent(out) :: pe(:), ptr(:)
  integer :: p

  npe = count(num > 0)
  ! the model reads the lists into com_struct arrays of size MAX_NEIGHBOR_PARTITIONS
  if (npe > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in MOD_PARTIT.F90 (and in the model''s src/types/mod_partit.F90) and recompile'
     stop
  endif
  allocate(pe(npe), ptr(npe+1))
  npe = 0
  ptr(1) = 1
  do p=0, size(num)-1
     if (num(p) == 0) cycle
     npe = npe+1
     pe(npe) = p
     ptr(npe+1) = ptr(npe)+num(p)
  end do
end subroutine set_pe_ptr
!=============================================================================
! position of the first entry of PE p in a list grouped by PE (counting sort)
subroutine first_pos(num, pos)
  integer, intent(in)  :: num(0:)
  integer, intent(out) :: pos(0:)
  integer :: p

  pos(0) = 1
  do p=1, size(num)-1
     pos(p) = pos(p-1)+num(p-1)
  end do
end subroutine first_pos
!=============================================================================
! adds PE p to the send list pes of one node/edge (once)
subroutine add_pe(pes, p, num_send, too_small)
  integer, intent(inout) :: pes(:), num_send(0:)
  integer, intent(in)    :: p
  logical, intent(inout) :: too_small
  integer :: k

  do k=1, size(pes)
     if (pes(k) == p) then
        return
     elseif (pes(k) == -1) then
        pes(k) = p
        num_send(p) = num_send(p)+1
        return
     end if
  end do
  too_small = .true.
end subroutine add_pe
!=============================================================================
subroutine write_dist_com(fileID, com)
  integer,        intent(in) :: fileID
  type(dist_com), intent(in) :: com

  write(fileID,*) com%rPEnum
  write(fileID,*) com%rPE
  write(fileID,*) com%rptr
  write(fileID,*) com%rlist
  write(fileID,*) com%sPEnum
  write(fileID,*) com%sPE
  write(fileID,*) com%sptr
  write(fileID,*) com%slist
end subroutine write_dist_com
end module dist_ext

!=============================================================================
SUBROUTINE com_global2local(partit, mesh)
use MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
IMPLICIT NONE

type(t_mesh),   intent(in),    target :: mesh
type(t_partit), intent(inout), target :: partit

INTEGER                            :: n, m
INTEGER, ALLOCATABLE, DIMENSION(:) :: temp

#include "associate_part_def.h"
#include "associate_mesh_ini.h"
#include "associate_part_ass.h"

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
SUBROUTINE save_dist_mesh(partit, mesh)
  USE g_CONFIG
  USE MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE o_ARRAYS 
  use com_global2local_interface
  use dist_ext
  IMPLICIT NONE

  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  Integer        n, m, q, q2, counter, fileID, nend, nini,ed(2)
  character*10   mype_string,npes_string
  character(MAX_PATH) file_name
  character(MAX_PATH) dist_mesh_dir
  integer, allocatable, dimension(:)  :: temp, ncount
  integer   n1, n2, flag, eledges(4)
  integer, allocatable, dimension(:)  :: xlist_nod2D
  type(dist_com) com_nod2D_full, com_edge2D

#include "associate_part_def.h"
#include "associate_mesh_ini.h"
#include "associate_part_ass.h"

  write(mype_string,'(i5.5)') mype  
  write(npes_string,"(I10)")  npes
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

  allocate(partit%myList_edge2D(4*myDim_elem2D))
  myList_edge2D=>partit%myList_edge2D
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

  ! hybrid meshes (dist_N/format.out): double-halo nodes, appended
  call nod_full_lists(partit, mesh, xlist_nod2D, com_nod2D_full)
  call edge_lists(partit, mesh, myDim_edge2D+eDim_edge2D, com_edge2D)
  write(fileID,*) size(xlist_nod2D)
  write(fileID,*) xlist_nod2D
  deallocate(xlist_nod2D)
  deallocate(partit%myList_edge2D)
  close(fileID)       


  ! =========================  
  ! communication information
  ! ========================= 
  call com_global2local(partit, mesh)   ! Do not call this subroutine earlier, global numbering is needed!
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
  deallocate(partit%myList_nod2D)
!!$  deallocate(com_nod2D%rPE)
!!$  deallocate(com_nod2D%rptr)
  deallocate(partit%com_nod2D%rlist)
!!$  deallocate(com_nod2D%sPE)
!!$  deallocate(com_nod2D%sptr)
  deallocate(partit%com_nod2D%slist)

  write(fileID,*) com_elem2D%rPEnum
  write(fileID,*) com_elem2D%rPE(1:com_elem2D%rPEnum)
  write(fileID,*) com_elem2D%rptr(1:com_elem2D%rPEnum+1)
  write(fileID,*) com_elem2D%rlist
  write(fileID,*) com_elem2D%sPEnum
  write(fileID,*) com_elem2D%sPE(1:com_elem2D%sPEnum)
  write(fileID,*) com_elem2D%sptr(1:com_elem2D%sPEnum+1)
  write(fileID,*) com_elem2D%slist
  deallocate(partit%myList_elem2D)
!!$  deallocate(com_elem2D%rPE)
!!$  deallocate(com_elem2D%rptr)
  deallocate(partit%com_elem2D%rlist)
!!$  deallocate(com_elem2D%sPE)
!!$  deallocate(com_elem2D%sptr)
  deallocate(partit%com_elem2D%slist)

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
  deallocate(partit%com_elem2D_full%rlist)
!!$  deallocate(com_elem2D_full%sPE)
!!$  deallocate(com_elem2D_full%sptr)
  deallocate(partit%com_elem2D_full%slist)

  ! hybrid meshes (dist_N/format.out): com_nod2D_full, com_edge2D, appended
  call write_dist_com(fileID, com_nod2D_full)
  call write_dist_com(fileID, com_edge2D)
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

     ! marker: the appended sections of my_list*.out / com_info*.out are present
     file_name=trim(dist_mesh_dir)//'format.out'
     open(fileID, file=file_name)
     write(fileID,*) DIST_FORMAT_VERSION
     close(fileID)
  end if
END subroutine  save_dist_mesh
