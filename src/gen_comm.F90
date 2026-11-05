module mod_gen_comm
    USE MOD_MESH
    USE MOD_PARTIT

    implicit none

    private
    public :: communication_nodn, communication_elemn, mymesh

contains

! Cell-vertex finite-volume version
! Contains: Routines that support parallelization
! set_par_support_ini run in the initialization phase.
! The communication rules are saved. 
! set_par_support in the main phase just allocates memory for buffer 
! arrays, the rest is read together with mesh from saved files.
!=======================================================================
subroutine communication_nodn(partit, mesh)
  implicit none
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer                  :: n, np, prank, el, r_count, s_count, q, i, j, nod, k, l
  integer                  :: num_send(0:partit%npes-1), num_recv(0:partit%npes-1), nd_count
  integer, allocatable     :: recv_from_pe(:), send_to_pes(:,:)
  logical                  :: max_laendereck_too_small=.false.
  integer                  :: IERR
  ! Assume we have 2D partitioning vector in part. Find communication rules
  ! Reduce allocation: find all neighboring PE
  nd_count = count(partit%part(1:mesh%nod2d) == partit%mype)
! write(*,*) nod2d
! write(*,*) MAX_LAENDERECK
! write(*,*) nd_count
! write(*,*) allocated(partit%myList_nod2D)
! write(*,*) partit%mype
  allocate(recv_from_pe(mesh%nod2d), send_to_pes(MAX_LAENDERECK,nd_count), &
           partit%myList_nod2D(nd_count), STAT=IERR)
  if (IERR /= 0) then
     write (*,*) 'Could not allocate arrays in communication_nodn'
     stop
  endif

  nd_count = 0
  do n=1,mesh%nod2D
     ! Checks if element el has nodes that belong to different partitions
     if (partit%part(n) == partit%mype) then
        nd_count = nd_count+1
        partit%myList_nod2D(nd_count)=n
     endif
  end do
  partit%myDim_nod2D=nd_count

  num_send(0:partit%npes-1) = 0
  num_recv(0:partit%npes-1) = 0
  recv_from_pe(1:mesh%nod2d) = -1
  send_to_pes(1:MAX_LAENDERECK,1:nd_count) = -1

  ! For the local nodes, run through the patch and collect all nodes in the patch
  ! (It would be simpler to re-use the adjacency matrix, but it is not a global
  !  variable... and I am lazy and want to keep the data structure simple)

  do l=1,nd_count
     n = partit%myList_nod2D(l)
     do i = 1, mesh%nod_in_elem2D_num(n)
        ! Over all elements el that the node n is part of
        el = mesh%nod_in_elem2D(i,n)

        ! Checks, if elements are quads or triangles
        q = 4    ! quads as default
        if (mesh%elem2D_nodes(1,el) == mesh%elem2D_nodes(4,el)) q=3  ! triangle

        do j = 1, q
           ! Over all nodes in every element el
           nod = mesh%elem2D_nodes(j,el)
           
           ! Checks, if node j is not in another partitionen
           if (partit%part(nod) /= partit%mype) then
              ! Checks, if not already considered to be received from this
              ! node from partition part(nod)
              if (recv_from_pe(nod) == -1) then  ! nod already collected to be received?
                 ! We have to receive this node from part(nod)
                 ! Add plus one node to the total number of
                 num_recv(partit%part(nod)) = num_recv(partit%part(nod)) + 1
                 ! ???
                 recv_from_pe(nod) = partit%part(nod)   ! recv_from_pe(recv_count) = nod  ! no new information, just handy
              endif
              ! Checks, if all possible connected partition
              ! And we have to send n to part(nod). Do we know this already?
              do k=1,MAX_LAENDERECK    !???
                 if (send_to_pes(k,l) == partit%part(nod)) then
                    exit  ! already collected
                 elseif (send_to_pes(k,l) == -1) then
                    send_to_pes(k,l) = partit%part(nod)
                    num_send(partit%part(nod)) = num_send(partit%part(nod)) + 1
                    exit
                 elseif (k== MAX_LAENDERECK) then
                    max_laendereck_too_small = .true.  ! Problem
                 endif
              enddo
           endif
        enddo
     enddo
  enddo

  if (max_laendereck_too_small) then
     print *,'Increase MAX_LAENDERECK in gen_modules_partitioning.F90 and recompile'
     stop
  endif

  ! Now, build the send and recv communication data structure
  ! To how many PE needs to be send and which PEs
!!$  com_nod2D%rPEnum = count(num_recv(0:npes-1) > 0)
!!$  com_nod2D%sPEnum = count(num_send(0:npes-1) > 0)

!!$  if (com_nod2D%rPEnum > MAX_NEIGHBOR_PARTITIONS .or.  &
!!$       com_nod2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
!!$     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
!!$     stop
!!$  endif
!!$  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
!!$  allocate(com_nod2D%sPE(com_nod2D%sPEnum))

  r_count = 0
  s_count = 0
  partit%com_nod2D%rptr(1) = 1
  partit%com_nod2D%sptr(1) = 1

  do np = 0, partit%npes-1
     if(num_recv(np) /= 0) then
        r_count = r_count+1
        partit%com_nod2D%rPE(r_count) = np
        partit%com_nod2D%rptr(r_count+1) =  partit%com_nod2D%rptr(r_count)+ num_recv(np)
     end if
     if(num_send(np) /= 0) then
        s_count = s_count+1
        partit%com_nod2D%sPE(s_count) = np
        partit%com_nod2D%sptr(s_count+1) =  partit%com_nod2D%sptr(s_count)+ num_send(np)
     end if
  enddo
  partit%com_nod2D%rPEnum = r_count
  partit%com_nod2D%sPEnum = s_count
  if (partit%com_nod2D%rPEnum > MAX_NEIGHBOR_PARTITIONS .or.  &
       partit%com_nod2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif

  ! Counts the number of node for each partition PE mype has to send/receive
  ! In ascending order of PE number
!!$  r_count = 0
!!$  s_count = 0
!!$  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1)) 
!!$  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1))

!!$  com_nod2D%rptr(1) = 1
!!$  com_nod2D%sptr(1) = 1
!!$
!!$  do r_count = 1, com_nod2D%rPEnum
!!$     np = com_nod2D%rPE(r_count)
!!$     com_nod2D%rptr(r_count+1) =  com_nod2D%rptr(r_count)+ num_recv(np)
!!$  enddo
!!$  do s_count = 1, com_nod2D%sPEnum
!!$     np = com_nod2D%sPE(s_count)
!!$     com_nod2D%sptr(s_count+1) =  com_nod2D%sptr(s_count)+ num_send(np)
!!$  enddo

  ! Lists themselves

  r_count = 0
  partit%eDim_nod2D=partit%com_nod2D%rptr(partit%com_nod2D%rPEnum+1)-1   
  allocate(partit%com_nod2D%rlist(partit%eDim_nod2D), &
           partit%com_nod2D%slist(partit%com_nod2D%sptr(partit%com_nod2D%sPEnum+1)-1), STAT=IERR) 
  if (IERR /= 0) then
     write (*,*) 'Could not allocate arrays in communication_nodn'
     stop
  endif

  do np = 1,partit%com_nod2D%rPEnum
     prank = partit%com_nod2D%rPE(np)
     do n = 1, mesh%nod2D
        if (recv_from_pe(n) == prank) then
           r_count = r_count+1
           partit%com_nod2D%rlist(r_count) = n
        end if
     end do
  end do

  s_count = 0
!!$  allocate(com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1)) 
  do np = 1,partit%com_nod2D%sPEnum
     prank = partit%com_nod2D%sPE(np)
     do l = 1, nd_count
        n = partit%myList_nod2D(l)
        if(any(send_to_pes(:,l) == prank)) then 
           s_count = s_count+1
           partit%com_nod2D%slist(s_count) = n
        end if
     end do
  end do

  ! Summary of this piece: mype receives
  ! information on external 2D nodes from
  ! comm_nod2D%rPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%rPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%rptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)
  ! Putting everything into structure takes many operations, but
  ! without the structure we will need to many names and arrays
  ! Do not forget that we need also send part.

  ! mype sends its data to
  ! comm_nod2D%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%sptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)

  deallocate(recv_from_pe, send_to_pes)
end subroutine communication_nodn

!==========================================================================
subroutine communication_elemn(partit, mesh)
  implicit none

  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer, allocatable     :: recv_from_pe(:), send_to_pes(:,:)
  logical                  :: max_laendereck_too_small=.false.
  integer                  :: n, k, ep, np, prank, el, nod
  integer                  :: p, q, j, elem, i, l, r_count, s_count, el_count
  integer                  :: num_send(0:partit%npes-1), num_recv(0:partit%npes-1)
  integer                  :: IERR
  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules. An elem is external to element n if neither of its nodes 
  ! belongs to PE, but it is among the neighbors. Element n belongs to PE if 
  ! any of its nodes does. 
  
  ! This routine takes into account 
  ! com_elem2D_full: all  neighbors  
  ! com_elem2D:      only those sharing an edge 


  !===========================================
  !  com_elem2D
  !===========================================

  allocate(recv_from_pe(mesh%elem2D), STAT=IERR)
  if (IERR /= 0) then
     write (*,*) 'Could not allocate arrays in communication_elemn'
     stop
  endif

  el_count = 0
  do el=1,mesh%elem2D
     ! Checks if element el has nodes that belong to different partitions
     if (any(partit%part(mesh%elem2D_nodes(1:4,el)) == partit%mype)) then
        el_count = el_count+1
        recv_from_pe(el_count) = el
     endif
  end do
  partit%myDim_elem2D=el_count

  allocate(partit%myList_elem2D(el_count), send_to_pes(MAX_LAENDERECK,el_count), STAT=IERR)
  if (IERR /= 0) then
     write (*,*) 'Could not allocate arrays in communication_elemn'
     stop
  endif

  partit%myList_elem2D(1:el_count) = recv_from_pe(1:el_count)
  num_send(0:partit%npes-1) = 0
  num_recv(0:partit%npes-1) = 0
  recv_from_pe(1:mesh%elem2D) = -1
  send_to_pes(1:MAX_LAENDERECK,1:el_count) = -1
  
  ! For the local elements, collect all adjacent (sharing an edge) elements
  ! that belong to other PEs
  
  do l=1,el_count
     el = partit%myList_elem2D(l)

     ! Checks if triangles or quads
     q = merge(3, 4, mesh%elem2D_nodes(1,el) == mesh%elem2D_nodes(4,el))

     do n = 1, q                       ! cycle through all nodes

        ! Neighboring element elem that shares an edge with element el
        elem = mesh%elem_neighbors(n,el)

        if (elem < 1) cycle  ! boundary, "ghost element"

        ! Check for elements to be received
        if (all(partit%part(mesh%elem2D_nodes(1:4,elem)) /= partit%mype) .and. recv_from_pe(elem)==-1) then
           ! elem to be received already collected?
           ! We have to receive elem from PE ep:
           ep = partit%part(mesh%elem2D_nodes(1,elem))  ! PE of first node is "main" owner
           num_recv(ep) = num_recv(ep) + 1
           recv_from_pe(elem) = ep  
        endif

        ! Check for elements to be sent to
        ! And maybe, we have to send el to the owners of elem
        ! 1. Is partition mype the main owner of element el?
        if (partit%part(mesh%elem2D_nodes(1,el)) == partit%mype .and. any(partit%part(mesh%elem2D_nodes(1:4,elem)) /= partit%mype)) then

           ! 2. Who owns element elem and needs to get element el? We must check all nodes!
           p=merge(3, 4, mesh%elem2D_nodes(1,elem) == mesh%elem2D_nodes(4,elem))

           do i=1,p
              ep = partit%part(mesh%elem2D_nodes(i,elem))

              ! 3. Is ep also an owner of el and no send is needed? This excludes also mype==ep.
              if (any(partit%part(mesh%elem2D_nodes(1:q,el)) == ep)) cycle

              ! 4. Ok, for the owner ep, check if sending el is already collected
              do k=1,MAX_LAENDERECK
                 if (send_to_pes(k,l) == ep) then
                    exit  ! already collected
                 elseif (send_to_pes(k,l) == -1) then
                    send_to_pes(k,l) = ep
                    num_send(ep) = num_send(ep) + 1
                    exit
                 elseif (k== MAX_LAENDERECK) then
                    max_laendereck_too_small = .true.  ! Problem
                 endif
              enddo
           enddo
        end if
     end do
  enddo

  if (max_laendereck_too_small) then
     print *,'Increase MAX_LAENDERECK in gen_modules_partitioning.F90 and recompile'
     stop
  endif
    
! Now, build the send and recv communication data structure
  r_count = 0
  s_count = 0
  partit%com_elem2D%rptr(1) = 1
  partit%com_elem2D%sptr(1) = 1

  do np = 0, partit%npes-1
     if(num_recv(np) /= 0) then
        r_count = r_count+1
        partit%com_elem2D%rPE(r_count) = np
        partit%com_elem2D%rptr(r_count+1) =  partit%com_elem2D%rptr(r_count)+ num_recv(np)
     end if
     if(num_send(np) /= 0) then
        s_count = s_count+1
        partit%com_elem2D%sPE(s_count) = np
        partit%com_elem2D%sptr(s_count+1) =  partit%com_elem2D%sptr(s_count)+ num_send(np)
     end if
  enddo

  partit%com_elem2D%rPEnum = r_count
  partit%com_elem2D%sPEnum = s_count
  if (partit%com_elem2D%rPEnum > MAX_NEIGHBOR_PARTITIONS .or.  &
       partit%com_elem2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif

  ! Lists themselves

  r_count = 0
  partit%eDim_elem2D=partit%com_elem2D%rptr(partit%com_elem2D%rPEnum+1)-1   
  allocate(partit%com_elem2D%rlist(partit%eDim_elem2D))
  do np = 1,partit%com_elem2D%rPEnum
     prank = partit%com_elem2D%rPE(np)
     do el = 1, mesh%elem2D
        if (recv_from_pe(el) == prank) then
           r_count = r_count+1
           partit%com_elem2D%rlist(r_count) = el
        end if
     end do
  end do
  
  s_count = 0
  allocate(partit%com_elem2D%slist(partit%com_elem2D%sptr(partit%com_elem2D%sPEnum+1)-1)) 
  do np = 1,partit%com_elem2D%sPEnum
     prank = partit%com_elem2D%sPE(np)
     do l = 1, el_count
        el = partit%myList_elem2D(l)
        if( any(send_to_pes(:,l) == prank)) then 
           s_count = s_count+1
           partit%com_elem2D%slist(s_count) = el
        end if
     end do
  end do
  
  !===========================================
  !  com_elem2D_full
  !===========================================

  ! The halo relations that are already determined can be kept.
  ! Just add the elements connected only via nodes.
  ! num_send(0:npes-1) = 0
  ! num_recv(0:npes-1) = 0
  ! recv_from_pe(1:elem2D) = -1
  ! send_to_pes(1:MAX_LAENDERECK,1:elem2D) = -1
  
  ! For the local elements, collect all adjacent (sharing a node) elements
  ! that belong to other PEs

  do l=1,el_count
     el = partit%myList_elem2D(l)

        ! Checks if triangles or quads
        q = merge(3, 4, mesh%elem2D_nodes(1,el) == mesh%elem2D_nodes(4,el))

        do n = 1, q                       ! cycle through all nodes

           nod = mesh%elem2D_nodes(n,el)

           ! Loop over all elements that belong to node nod 
           do j = 1, mesh%nod_in_elem2D_num(nod)  ! and for each node, through its patch
              elem = mesh%nod_in_elem2D(j,nod)

              ! Check for elements to be received
              if (all(partit%part(mesh%elem2D_nodes(1:4,elem)) /= partit%mype) .and. recv_from_pe(elem)==-1) then
                 ! elem to be received already collected?
                 ! We have to receive elem from PE ep:
                 ep = partit%part(mesh%elem2D_nodes(1,elem))  ! PE of first node is "main" owner
                 num_recv(ep) = num_recv(ep) + 1
                 recv_from_pe(elem) = ep  
              endif

              ! Check for elements to be sent to
              ! And maybe, we have to send el to the owners of elem
              ! This gets more complicated:
              ! 1. Is partition mype the main owner of element el?
              if (partit%part(mesh%elem2D_nodes(1,el)) == partit%mype .and. any(partit%part(mesh%elem2D_nodes(1:4,elem)) /= partit%mype)) then

                 ! 2. Who owns element elem and needs to get element el? We must check all nodes!
                 p=merge(3, 4, mesh%elem2D_nodes(1,elem) == mesh%elem2D_nodes(4,elem))

                 do i=1,p
                    ep = partit%part(mesh%elem2D_nodes(i,elem))

                    ! 3. Is ep also an owner of el and no send is needed? This excludes also mype==ep.
                    if (any(partit%part(mesh%elem2D_nodes(1:q,el)) == ep)) cycle

                    ! 4. Ok, for the owner ep, check if sending el is already collected
                    do k=1,MAX_LAENDERECK
                       if (send_to_pes(k,l) == ep) then
                          exit  ! already collected
                       elseif (send_to_pes(k,l) == -1) then
                          send_to_pes(k,l) = ep
                          num_send(ep) = num_send(ep) + 1
                          exit
                       elseif (k== MAX_LAENDERECK) then
                          max_laendereck_too_small = .true.  ! Problem
                       endif
                    enddo
                 end do
              endif
           end do
        end do
  enddo
  
  if (max_laendereck_too_small) then
     print *,'Increase MAX_LAENDERECK in gen_modules_partitioning.F90 and recompile'
     stop
  endif

  ! Now, build the send and recv communication data structure

  r_count = 0
  s_count = 0
  partit%com_elem2D_full%rptr(1) = 1
  partit%com_elem2D_full%sptr(1) = 1

  do np = 0, partit%npes-1
     if(num_recv(np) /= 0) then
        r_count = r_count+1
        partit%com_elem2D_full%rPE(r_count) = np
        partit%com_elem2D_full%rptr(r_count+1) =  partit%com_elem2D_full%rptr(r_count)+ num_recv(np)
     end if
     if(num_send(np) /= 0) then
        s_count = s_count+1
        partit%com_elem2D_full%sPE(s_count) = np
        partit%com_elem2D_full%sptr(s_count+1) =  partit%com_elem2D_full%sptr(s_count)+ num_send(np)
     end if
  enddo

  partit%com_elem2D_full%rPEnum = r_count
  partit%com_elem2D_full%sPEnum = s_count

  ! Lists themselves

  r_count = 0
  allocate(partit%com_elem2D_full%rlist(partit%com_elem2D_full%rptr(partit%com_elem2D_full%rPEnum+1)-1)) 
  do np = 1,partit%com_elem2D_full%rPEnum
     prank = partit%com_elem2D_full%rPE(np)
     do el = 1, mesh%elem2D
        if (recv_from_pe(el) == prank) then
           r_count = r_count+1
           partit%com_elem2D_full%rlist(r_count) = el
        end if
     end do
  end do

  s_count = 0
  allocate(partit%com_elem2D_full%slist(partit%com_elem2D_full%sptr(partit%com_elem2D_full%sPEnum+1)-1)) 
  do np = 1,partit%com_elem2D_full%sPEnum
     prank = partit%com_elem2D_full%sPE(np)
     do l = 1, el_count
        el = partit%myList_elem2D(l)
        if( any(send_to_pes(:,l) == prank)) then 
           s_count = s_count+1
           partit%com_elem2D_full%slist(s_count) = el
        end if
     end do
  end do

  deallocate(recv_from_pe, send_to_pes)
end subroutine communication_elemn
!==========================================================================
subroutine mymesh(partit, mesh)
  implicit none

  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer                  :: n, counter, q, k, elem, q2, eledges(4)
  integer, allocatable     :: aux(:)
  !======= NODES 

  ! Owned nodes + external nodes which I need:
!!$  myDim_nod2D=count(part(:) == mype)
!!$  eDim_nod2D=com_nod2D%rptr(com_nod2D%rPEnum+1)-1
  ! Check if the length of myList_nod2D is sufficient

!!$  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D))
!!$  counter=0   
!!$  do n=1, nod2D
!!$     if (part(n)==mype) then
!!$        counter=counter+1
!!$        myList_nod2D(counter)=n
!!$     end if
!!$  end do
!!$  myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)=&
!!$       com_nod2D%rlist(1:eDim_nod2D)
  ! Summary:  	     
  ! myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)
  ! contains external nodes which mype needs;    
  ! myList_nod2D(1:myDim_nod2D) contains owned nodes

  !======= ELEMENTS
  ! 2D elements 
  ! Element belongs to PE if any of its nodes is owned by PE
  ! Element is external if it is a neighbor and does not contain owned nodes
  ! The external part is needed for FVCOM type of discretization.
  
!!$  counter=0
!!$  do n=1, elem2D
!!$     q2 = merge(3,4,elem2D_nodes(1,n) == elem2D_nodes(4,n))
!!$     do q=1,q2
!!$        if(part(elem2D_nodes(q,n))==mype) then
!!$           counter=counter+1
!!$           myList_elem2D(counter)=n
!!$           exit
!!$        end if
!!$     end do
!!$  end do
!!$  myDim_elem2D=counter
!!$  eDim_elem2D=com_elem2D%rptr(com_elem2D%rPEnum+1)-1   
  ! =======
  ! full element neighbourhood requires 
  ! a longer list     
  ! ======= 
!!$  allocate(aux(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1))    
!!$  aux=0
!!$  do n=1,com_elem2D%rptr(com_elem2D%rPEnum+1)-1         
!!$     k=com_elem2D%rlist(n)
!!$     do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
!!$       if(com_elem2D_full%rlist(q)==k) aux(q)=1
!!$     end do
!!$  end do
!!$  ! Test: 
!!$  if(sum(aux).ne.eDim_elem2D) write(*,*) 'mymesh problem'
!!$  
!!$  counter=0
!!$  do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
!!$     if(aux(q)==0) counter=counter+1   
!!$  end do
!!$  eXDim_elem2D=counter
!!$
!!$  myList_elem2D(myDim_elem2D+1:myDim_elem2D+eDim_elem2D)=&
!!$       com_elem2D%rlist(1:eDim_elem2D)
!!$  counter=myDim_elem2D+eDim_elem2D
!!$  do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
!!$     if(aux(q)==0) then
!!$     counter=counter+1
!!$     myList_elem2D(counter)=com_elem2D_full%rlist(q)
!!$     end if   
!!$  end do
!!$  deallocate(aux)
              
  ! Summary: 
  ! myList_elem2D(1:myDim_elem2D) contains elements with at least single owned node.
  ! myList_elem2D(myDim_elem2D+1:myDim_elem2D+eDim_elem2D) contains elements-neighbours.
  ! They are needed when interpolation is done in FV type of code.
  ! The final piece from myDim_elem2D+eDim_elem2D to 
  ! myDim_elem2D+eDim_elem2D+eXDim_elem2D is needed for MUSCL-type advection

! ======== EDGES 
  ! Owned edges (both nodes are mine)+ shared edges I do computations 
  ! at (only one node is mine; some other PE updates them
  ! simultaneously with me but I do not care ) + 
  ! external edges which I need (neither of nodes is mine, but they
  ! belong to elements in myList_elem:
  counter=0
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
     elem=partit%myList_elem2D(n)
     eledges=mesh%elem_edges(:,elem)
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
  ! Summary:  	     
  ! myList_edge2D(myDim_edge2D+1:myDim_edge2D+eDim_edge2D)
  ! contains external edges which mype needs;    
  ! myList_edge2D(1:myDim_edge2D) contains owned edges +
  ! shared edges which mype updates
end subroutine mymesh
!=================================================================

end module mod_gen_comm
