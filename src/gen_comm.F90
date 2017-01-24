! Cell-vertex finite-volume version
! Contains: Routines that support parallelization
! set_par_support_ini run in the initialization phase.
! The communication rules are saved. 
! set_par_support in the main phase just allocates memory for buffer 
! arrays, the rest is read together with mesh from saved files.

!KK: moved par_ex,set_par_support,set_par_support_ini to module g_PARSUP

! ===============================================================
!=======================================================================
subroutine communication_nodn
  use o_MESH
  use g_PARSUP
  implicit none

  integer n,np, nz, prank, elem, elnodes(3), epe(3), counter, nini, nend
  integer, allocatable :: aux(:,:), pnum(:,:),pmap(:)

  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules
  ! allocate(aux(npes,nod2D))
  ! Reduce allocation: find all neighboring PE
  allocate(pnum(npes,npes))

  do elem=1,elem2D
     elnodes=elem2D_nodes(:,elem)
     epe=part(elnodes)+1
     if(epe(1).ne.epe(2)) then
        pnum(epe(1), epe(2))=1
        pnum(epe(2), epe(1))=1
     end if

     if(epe(2).ne.epe(3)) then
        pnum(epe(3), epe(2))=1
        pnum(epe(2), epe(3))=1
     end if
     if(epe(1).ne.epe(3)) then
        pnum(epe(1), epe(3))=1
        pnum(epe(3), epe(1))=1
     end if
  end do
  ! mype interacts with a limited number of PEs
  allocate(pmap(npes))
  pmap=0
  counter=0
  DO n=1,npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        pmap(n)=counter
     end if
     if(n==mype+1) then
        counter=counter+1
        pmap(n)=counter
     end if
     
  END DO    
  allocate(aux(counter,nod2D)) 
  ! This has a much smaller size if npes is large
  ! but we need pmap to address it
  aux=0
  pnum=0
  
  
  do elem=1,elem2D
     elnodes=elem2D_nodes(:,elem)
     epe=part(elnodes)+1
     if(epe(1).ne.epe(2)) then
        if(pmap(epe(1)).ne.0) aux(pmap(epe(1)), elnodes(2))=1
        if(pmap(epe(2)).ne.0) aux(pmap(epe(2)), elnodes(1))=1
     end if

     if(epe(2).ne.epe(3)) then
        if(pmap(epe(3)).ne.0) aux(pmap(epe(3)), elnodes(2))=1
        if(pmap(epe(2)).ne.0) aux(pmap(epe(2)), elnodes(3))=1
     end if
     if(epe(1).ne.epe(3)) then
        if(pmap(epe(1)).ne.0) aux(pmap(epe(1)), elnodes(3))=1
        if(pmap(epe(3)).ne.0) aux(pmap(epe(3)), elnodes(1))=1
     end if
  end do

  do n=1, nod2D
     do np=1, npes
        if(pmap(np).ne.0) then
        if(aux(pmap(np),n).ne.0) then 
           pnum(np,part(n)+1)=pnum(np,part(n)+1)+1
        end if
	end if
     end do
  end do


  ! We know how many external nodes each PE needs
  ! This is the 'receive' list   
  ! com_nod2D for 2D nodes
  ! 

  ! The number of external PE I receive information from
  com_nod2D%rPEnum=0
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        com_nod2D%rPEnum=com_nod2D%rPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_nod2D%rPE(counter)=n-1
     end if
  end do


  ! Ptr to list of external nodes ordered by external PE ranks

  counter=0
  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1)) 
  com_nod2D%rptr(1)=1
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_nod2D%rptr(counter+1)=com_nod2D%rptr(counter)+ pnum(mype+1,n)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_nod2D%rlist(com_nod2D%rptr(com_nod2D%rPEnum+1)-1)) 
  do np=1,com_nod2D%rPEnum
     prank=com_nod2D%rPE(np)
     do n=1, nod2D
        if((aux(pmap(mype+1),n)==1).and.(part(n)==prank)) then
           counter=counter+1
           com_nod2D%rlist(counter)=n
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
  ! Do not forget that we need also send part, and we need analogous
  ! information for 3D nodes and edges.

  ! SENDING PART
  com_nod2D%sPEnum=0
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        com_nod2D%sPEnum=com_nod2D%sPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_nod2D%sPE(com_nod2D%sPEnum))
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_nod2D%sPE(counter)=n-1
     end if
  end do

  ! Ptr to list of external nodes ordered by external PE ranks
  counter=0
  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1)) 
  com_nod2D%sptr(1)=1
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_nod2D%sptr(counter+1)=com_nod2D%sptr(counter)+ pnum(n,mype+1)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1)) 
  do np=1,com_nod2D%sPEnum
     prank=com_nod2D%sPE(np)
     do n=1, nod2D
        if(pmap(prank+1).ne.0) then 
        if((aux(pmap(prank+1),n)==1).and.(part(n)==mype)) then
           counter=counter+1
           com_nod2D%slist(counter)=n
        end if
	end if
     end do
  end do

  ! mype sends its data to
  ! comm_nod2D%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%sptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)

  deallocate(pmap, pnum,aux)

end subroutine communication_nodn

!==========================================================================
subroutine communication_elemn
  use o_MESH
  use g_PARSUP 
  implicit none
  integer n, k, ep, np1, np2, np3, prank, elem, elemn(3), counter
  integer, allocatable :: aux(:,:), pnum(:,:), pmap(:)

  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules. An elem is external to element n if neither of its nodes 
  ! belongs to PE, but it is among the neighbors. Element n belongs to PE if 
  ! any of its nodes does. 

  !allocate(aux(npes,elem2D))
  ! This version uses a reduced memory allocation 
  
  
  allocate(pnum(npes,npes))
  ! Find PEs I must communicate with
  pnum=0
  do elem=1,elem2D
     elemn=elem_neighbors(:,elem)
    do k=1,3 
     ep=part(elem2D_nodes(k,elem))
     do n=1,3
        if(elemn(n)<1) cycle
	np1=part(elem2D_nodes(1,elemn(n)))
	np2=part(elem2D_nodes(2,elemn(n)))
	np3=part(elem2D_nodes(3,elemn(n)))
	if((np1.ne.ep).and.(np2.ne.ep).and.(np3.ne.ep)) then
	   pnum(ep+1,np1+1)=1
	   pnum(ep+1,np2+1)=1
	   pnum(ep+1,np3+1)=1
        end if
     end do
    end do 
  end do
  ! pnum =1 for all PEs that can in principle communicate with each other 
  ! but: mype communicates with a subset of npes
  allocate(pmap(npes))
  pmap=0
  counter=0
  DO n=1,npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        pmap(n)=counter
     end if
     if(n==mype+1) then
        counter=counter+1
        pmap(n)=counter
     end if
     
  END DO    
  allocate(aux(counter,elem2D)) 
  ! This has a much smaller size if npes is large
  ! but we need pmap to address it
  
  aux=0
  pnum=0
  do elem=1,elem2D
     elemn=elem_neighbors(:,elem)
    do k=1,3 
     ep=part(elem2D_nodes(k,elem))
     if(pmap(ep+1).ne.0) then
     do n=1,3
        if(elemn(n)<1) cycle
	np1=part(elem2D_nodes(1,elemn(n)))
	np2=part(elem2D_nodes(2,elemn(n)))
	np3=part(elem2D_nodes(3,elemn(n)))
	if((np1.ne.ep).and.(np2.ne.ep).and.(np3.ne.ep)) then
	   aux(pmap(ep+1),elemn(n))=1
        end if
     end do
     end if
    end do 
  end do
  
  do n=1, elem2D
     do ep=1, npes
        if(pmap(ep).ne.0) then
        if(aux(pmap(ep),n).ne.0) then 
           pnum(ep,part(elem2D_nodes(1,n))+1)=pnum(ep,part(elem2D_nodes(1,n))+1)+1
        end if
	end if
     end do
  end do
  ! I take information on the external element from the PE 
  ! of the first node of this element.
  
  ! We know how many external elements each PE needs
  ! This is the 'receive' list   
  ! com_elem2D for 2D elements
  ! 

  ! The number of external PEs I receive information from
  com_elem2D%rPEnum=0
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        com_elem2D%rPEnum=com_elem2D%rPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_elem2D%rPE(com_elem2D%rPEnum))
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_elem2D%rPE(counter)=n-1
     end if
  end do


  ! Ptr to list of external elements ordered by external PE ranks

  counter=0
  allocate(com_elem2D%rptr(com_elem2D%rPEnum+1)) 
  com_elem2D%rptr=1
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_elem2D%rptr(counter+1)=com_elem2D%rptr(counter)+ pnum(mype+1,n)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_elem2D%rlist(com_elem2D%rptr(com_elem2D%rPEnum+1)-1)) 
  do ep=1,com_elem2D%rPEnum
     prank=com_elem2D%rPE(ep)
     do n=1, elem2D
        if((aux(pmap(mype+1),n)==1).and.(part(elem2D_nodes(1,n))==prank)) then
           counter=counter+1
           com_elem2D%rlist(counter)=n
        end if
     end do
  end do

  ! Summary of this piece: mype receives
  ! information on external 2D elements from
  ! comm_elem2D%rPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_elem2D%rPE(:)
  ! Pointers to external node numbers are in
  ! comm_elem2D%rptr(:)
  ! The node numbers are in 
  ! comm_elem2D%list(:)
  
  ! SENDING PART
  com_elem2D%sPEnum=0
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        com_elem2D%sPEnum=com_elem2D%sPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_elem2D%sPE(com_elem2D%sPEnum))
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_elem2D%sPE(counter)=n-1
     end if
  end do

  ! Ptr to list of external nodes ordered by external PE ranks
  counter=0
  allocate(com_elem2D%sptr(com_elem2D%sPEnum+1)) 
  com_elem2D%sptr=1
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_elem2D%sptr(counter+1)=com_elem2D%sptr(counter)+ pnum(n,mype+1)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_elem2D%slist(com_elem2D%sptr(com_elem2D%sPEnum+1)-1)) 
  do ep=1,com_elem2D%sPEnum
     prank=com_elem2D%sPE(ep)
     if(pmap(prank+1).ne.0) then
     do n=1, elem2D
        if((aux(pmap(prank+1),n)==1).and.(part(elem2D_nodes(1,n))==mype)) then
           counter=counter+1
           com_elem2D%slist(counter)=n
        end if
     end do
     end if
  end do

  ! mype sends its data to
  ! comm_elem2D%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_elem2D%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_elem2D%sptr(:)
  ! The node numbers are in 
  ! comm_elem2D%list(:)

  deallocate(pnum)
  deallocate(pmap)
  deallocate(aux)
end subroutine communication_elemn
!========================================================================
subroutine communication_elem_fulln
  use o_MESH
  use g_PARSUP 
  implicit none
  integer n, k, ep, np1, np2, np3, prank, elem, elemn(30), counter, enodes(3),q
  integer, allocatable, dimension(:,:) :: aux, pnum, aux_e_list
  integer, allocatable, dimension(:) :: aux_e, aux_e_num, pmap 

  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules. An elem is external to element n if neither of its nodes 
  ! belongs to PE. Element n belongs to PE if any of its nodes does.
  
  ! This routine takes into account all posible neighbors, not only 
  ! those sharing an edge as in communication_elem

  !allocate(aux(npes,elem2D))
  allocate(pnum(npes,npes))
  pnum=0
  
  ! =========
  ! Find elements that have common nodes with a given element
  ! =========
  
  allocate(aux_e(elem2D), aux_e_num(elem2D), aux_e_list(30,elem2D))
  aux_e=0
  DO n=1, elem2D
     enodes=elem2D_nodes(:,n)
     ep=1
     aux_e_list(ep,n)=n
     DO q=1,3
       DO k=1, nod_in_elem2D_num(enodes(q))
        elem=nod_in_elem2D(k,enodes(q))
        if(aux_e(elem)==0) then
	   aux_e(elem)=1 
	   if(ep>30) then
	   write(*,*) 'ERROR in communication_elem_fulln'
	   call par_ex(1)
	   STOP
	   end if
           ep=ep+1
	   aux_e_list(ep,n)=elem   
        end if
       END DO
     END DO
     aux_e_num(n)=ep
     aux_e(aux_e_list(1:ep,n))=0
  END DO   
  
  
  ! ========
  ! How many PEs are involved?   
  ! ========     	 
  do elem=1,elem2D
     elemn= aux_e_list(:,elem) !elem_neighbors(:,elem)
    do k=1,3 
     ep=part(elem2D_nodes(k,elem))
     do n=1, aux_e_num(elem) !3
        if(elemn(n)<1) cycle
	np1=part(elem2D_nodes(1,elemn(n)))
	np2=part(elem2D_nodes(2,elemn(n)))
	np3=part(elem2D_nodes(3,elemn(n)))
	if((np1.ne.ep).and.(np2.ne.ep).and.(np3.ne.ep)) then
	   pnum(ep+1,np1+1)=1
	   pnum(ep+1,np2+1)=1
	   pnum(ep+1,np3+1)=1
        end if
     end do
    end do 
  end do
   ! pnum =1 for all PEs that can in principle communicate with each other 
  ! but: mype communicates with a subset of npes
  allocate(pmap(npes))
  pmap=0
  counter=0
  DO n=1,npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        pmap(n)=counter
     end if
     if(n==mype+1) then
        counter=counter+1
        pmap(n)=counter
     end if
     
  END DO    
  allocate(aux(counter,elem2D)) 
  ! This has a much smaller size if npes is large
  ! but we need pmap to address it
  aux=0
  pnum=0
  
  
  do elem=1,elem2D
     elemn= aux_e_list(:,elem) !elem_neighbors(:,elem)
    do k=1,3 
     ep=part(elem2D_nodes(k,elem))
     if(pmap(ep+1).ne.0) then
     do n=1, aux_e_num(elem) !3
        if(elemn(n)<1) cycle
	np1=part(elem2D_nodes(1,elemn(n)))
	np2=part(elem2D_nodes(2,elemn(n)))
	np3=part(elem2D_nodes(3,elemn(n)))
	if((np1.ne.ep).and.(np2.ne.ep).and.(np3.ne.ep)) then
	   aux(pmap(ep+1),elemn(n))=1
        end if
     end do
     end if
    end do 
  end do
  
  deallocate(aux_e_list, aux_e_num, aux_e)
  
  do n=1, elem2D
     do ep=1, npes
        if(pmap(ep).ne.0) then
	if(aux(pmap(ep),n).ne.0) then 
           pnum(ep,part(elem2D_nodes(1,n))+1)=pnum(ep,part(elem2D_nodes(1,n))+1)+1
        end if
	end if
     end do
  end do
  ! I take information on the external element from the PE 
  ! of the first node of this element.
  
  ! We know how many external elements each PE needs
  ! This is the 'receive' list   
  ! com_elem2D for 2D elements
  ! 

  ! The number of external PEs I receive information from
  com_elem2D_full%rPEnum=0
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        com_elem2D_full%rPEnum=com_elem2D_full%rPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_elem2D_full%rPE(com_elem2D_full%rPEnum))
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_elem2D_full%rPE(counter)=n-1
     end if
  end do


  ! Ptr to list of external nodes ordered by external PE ranks

  counter=0
  allocate(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)) 
  com_elem2D_full%rptr=1
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_elem2D_full%rptr(counter+1)=com_elem2D_full%rptr(counter)+ pnum(mype+1,n)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_elem2D_full%rlist(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1)) 
  do ep=1,com_elem2D_full%rPEnum
     prank=com_elem2D_full%rPE(ep)
     do n=1, elem2D
        if((aux(pmap(mype+1),n)==1).and.(part(elem2D_nodes(1,n))==prank)) then
           counter=counter+1
           com_elem2D_full%rlist(counter)=n
        end if
     end do
  end do

  ! Summary of this piece: mype receives
  ! information on external 2D elements from
  ! comm_elem2D_full%rPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_elem2D_full%rPE(:)
  ! Pointers to external node numbers are in
  ! comm_elem2D_full%rptr(:)
  ! The node numbers are in 
  ! comm_elem2D_full%list(:)
  
  ! SENDING PART
  com_elem2D_full%sPEnum=0
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        com_elem2D_full%sPEnum=com_elem2D_full%sPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_elem2D_full%sPE(com_elem2D_full%sPEnum))
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_elem2D_full%sPE(counter)=n-1
     end if
  end do

  ! Ptr to list of external nodes ordered by external PE ranks
  counter=0
  allocate(com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)) 
  com_elem2D_full%sptr=1
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_elem2D_full%sptr(counter+1)=com_elem2D_full%sptr(counter)+ pnum(n,mype+1)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_elem2D_full%slist(com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)-1)) 
  do ep=1,com_elem2D_full%sPEnum
     prank=com_elem2D_full%sPE(ep)
     if(pmap(prank+1).ne.0) then
     do n=1, elem2D
        if((aux(pmap(prank+1),n)==1).and.(part(elem2D_nodes(1,n))==mype)) then
           counter=counter+1
           com_elem2D_full%slist(counter)=n
        end if
     end do
     end if
  end do

  ! mype sends its data to
  ! comm_elem2D_full%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_elem2D_full%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_elem2D_full%sptr(:)
  ! The node numbers are in 
  ! comm_elem2D_full%list(:)

  deallocate(pnum)
  deallocate(pmap)
  deallocate(aux)
end subroutine communication_elem_fulln
!==========================================================================
subroutine communication_edgen
  use o_MESH
  use g_PARSUP 
  implicit none
  integer n,np, nz, nu,prank, elem, eledges(3), counter
  integer, allocatable, dimension(:,:) :: aux, pnum
  integer, allocatable, dimension(:) :: pmap

  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules. An edge is external to a node n if the triple
  ! n, edges(:,edge) forms an element while both PEs of edge 
  ! are different from PE of n

  ! an edge is owned by PE if any its nodes belongs to PE 

  !allocate(aux(npes,edge2D))
  allocate(pnum(npes,npes))
  pnum=0
  allocate(pmap(npes))
  pmap=0
  
  do elem=1,elem2D
     eledges=elem_edges(:,elem)
      ! if both PEs of edge opposing a node
      ! are different from node's PE then the edge should be communicated
     do n=1,3
        np=elem2D_nodes(n,elem)
        ! find edge that does not contain np
        counter=0
        do nz=1,3
	   if((edges(1,eledges(nz)).ne.np).and. &
                (edges(2,eledges(nz)).ne.np)) then
	      counter=nz
	      exit
	   end if
        end do
        if((part(np).ne.part(edges(1,eledges(counter)))).and.&   
	     (part(np).ne.part(edges(2,eledges(counter))))) then 
	   pnum(part(np)+1, part(edges(1,eledges(counter)))+1)=1
	   pnum(part(np)+1, part(edges(2,eledges(counter)))+1)=1
        end if
     end do
  end do
  counter=0
  DO n=1,npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        pmap(n)=counter
     end if
     if(n==mype+1) then
        counter=counter+1
        pmap(n)=counter
     end if
  END DO    
  allocate(aux(counter,edge2D)) 
  ! This has a much smaller size if npes is large
  ! but we need pmap to address it
  aux=0
  pnum=0
  do elem=1,elem2D
     eledges=elem_edges(:,elem)
      do n=1,3
        np=elem2D_nodes(n,elem)
        ! find edge that does not contain np
        counter=0
	if(pmap(part(np)+1).ne.0) then
        do nz=1,3
	   if((edges(1,eledges(nz)).ne.np).and. &
                (edges(2,eledges(nz)).ne.np)) then
	      counter=nz
	      exit
	   end if
        end do
        if((part(np).ne.part(edges(1,eledges(counter)))).and.&   
	     (part(np).ne.part(edges(2,eledges(counter))))) then 
	   aux(pmap(part(np)+1),eledges(counter))=1
        end if
	end if
      end do
  end do

  do n=1, edge2D
     do np=1, npes
        if(pmap(np).ne.0) then
        if(aux(pmap(np),n).ne.0) then 
           pnum(np,part(edges(1,n))+1)=pnum(np,part(edges(1,n))+1)+1
        end if
	end if
     end do
  end do

  ! We know how many external edges each PE needs
  ! This is the 'receive' list   
  ! com_edge2D for 2D nodes
  ! 

  ! The number of external PEs I receive information from
  com_edge2D%rPEnum=0
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        com_edge2D%rPEnum=com_edge2D%rPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_edge2D%rPE(com_edge2D%rPEnum))
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_edge2D%rPE(counter)=n-1
     end if
  end do


  ! Ptr to list of external nodes ordered by external PE ranks

  counter=0
  allocate(com_edge2D%rptr(com_edge2D%rPEnum+1)) 
  com_edge2D%rptr=1
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_edge2D%rptr(counter+1)=com_edge2D%rptr(counter)+ pnum(mype+1,n)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_edge2D%rlist(com_edge2D%rptr(com_edge2D%rPEnum+1)-1)) 
  do np=1,com_edge2D%rPEnum
     prank=com_edge2D%rPE(np)
     do n=1, edge2D
        if((aux(pmap(mype+1),n)==1).and.(part(edges(1,n))==prank)) then
           counter=counter+1
           com_edge2D%rlist(counter)=n
        end if
     end do
  end do

  ! Summary of this piece: mype receives
  ! information on external 2D edges from
  ! comm_edge2D%rPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_edge2D%rPE(:)
  ! Pointers to external node numbers are in
  ! comm_edge2D%rptr(:)
  ! The node numbers are in 
  ! comm_edge2D%list(:)
  ! Putting everything into structure takes many operations, but
  ! without the structure we will need to many names and arrays
  ! Do not forget that we need also send part, and we need analogous
  ! information for 3D nodes and edges.

  ! SENDING PART
  com_edge2D%sPEnum=0
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        com_edge2D%sPEnum=com_edge2D%sPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_edge2D%sPE(com_edge2D%sPEnum))
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_edge2D%sPE(counter)=n-1
     end if
  end do

  ! Ptr to list of external nodes ordered by external PE ranks
  counter=0
  allocate(com_edge2D%sptr(com_edge2D%sPEnum+1)) 
  com_edge2D%sptr=1
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_edge2D%sptr(counter+1)=com_edge2D%sptr(counter)+ pnum(n,mype+1)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_edge2D%slist(com_edge2D%sptr(com_edge2D%sPEnum+1)-1)) 
  do np=1,com_edge2D%sPEnum
     prank=com_edge2D%sPE(np)
     if(pmap(prank+1).ne.0) then
     do n=1, edge2D
        if((aux(pmap(prank+1),n)==1).and.(part(edges(1,n))==mype)) then
           counter=counter+1
           com_edge2D%slist(counter)=n
        end if
     end do
     end if
  end do

  ! mype sends its data to
  ! comm_edge2D%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_edge2D%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_edge2D%sptr(:)
  ! The node numbers are in 
  ! comm_edge2D%list(:)

  deallocate(pnum)
  deallocate(pmap)
  deallocate(aux)

end subroutine communication_edgen
!========================================================================
subroutine mymesh
  use o_MESH
  use g_PARSUP 
  implicit none
  integer                :: n, counter, q, k, elem
  integer, allocatable   :: aux(:)

  !======= NODES 

  ! Owned nodes + external nodes which I need:
  counter=0
  do n=1, nod2D
     if (part(n)==mype) counter=counter+1
  end do
  myDim_nod2D=counter
  eDim_nod2D=com_nod2D%rptr(com_nod2D%rPEnum+1)-1   
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D))
  counter=0   
  do n=1, nod2D
     if (part(n)==mype) then
        counter=counter+1
        myList_nod2D(counter)=n
     end if
  end do
  myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)=&
       com_nod2D%rlist
  ! Summary:  	     
  ! myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)
  ! contains external nodes which mype needs;    
  ! myList_nod2D(1:myDim_nod2D) contains owned nodes

  !======= ELEMENTS
  ! 2D elements 
  ! Element belongs to PE if any of its nodes is owned by PE
  ! Element is external if it is a neighbor and does not contain owned nodes
  ! The external part is needed for FVCOM type of discretization.
  
  counter=0
  do n=1, elem2D
     do q=1,3
        if(part(elem2D_nodes(q,n))==mype) then
           counter=counter+1
           exit
        end if
     end do
  end do
  myDim_elem2D=counter
  eDim_elem2D=com_elem2D%rptr(com_elem2D%rPEnum+1)-1   
  ! =======
  ! full element neighbourhood requires 
  ! a longer list     
  ! ======= 
  allocate(aux(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1))    
  aux=0
  do n=1,com_elem2D%rptr(com_elem2D%rPEnum+1)-1         
     k=com_elem2D%rlist(n)
     do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
       if(com_elem2D_full%rlist(q)==k) aux(q)=1
     end do
  end do
  ! Test: 
  if(sum(aux).ne.eDim_elem2D) write(*,*) 'mymesh problem'
  
  counter=0
  do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
     if(aux(q)==0) counter=counter+1   
  end do
     eXDim_elem2D=counter
  allocate(myList_elem2D(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  counter=0
  do n=1, elem2D
     do q=1,3
        if(part(elem2D_nodes(q,n))==mype) then
           counter=counter+1
           myList_elem2D(counter)=n
           exit
        end if
     end do
  end do
  myList_elem2D(myDim_elem2D+1:myDim_elem2D+eDim_elem2D)=&
       com_elem2D%rlist
  counter=myDim_elem2D+eDim_elem2D
  do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
     if(aux(q)==0) then
     counter=counter+1
     myList_elem2D(counter)=com_elem2D_full%rlist(q)
     end if   
  end do
  deallocate(aux)
              
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
  do n=1, edge2D
     do q=1,2 
        if (part(edges(q,n))==mype) then
           counter=counter+1
           exit
        end if
     end do
  end do
  myDim_edge2D=counter
  eDim_edge2D=com_edge2D%rptr(com_edge2D%rPEnum+1)-1   
  allocate(myList_edge2D(myDim_edge2D+eDim_edge2D))
  counter=0   
  do n=1, edge2D
     do q=1,2
        if (part(edges(q,n))==mype) then
           counter=counter+1
           myList_edge2D(counter)=n
           exit
        end if
     end do
  end do
  myList_edge2D(myDim_edge2D+1:myDim_edge2D+eDim_edge2D)=&
       com_edge2D%rlist

  ! Summary:  	     
  ! myList_edge2D(myDim_edge2D+1:myDim_edge2D+eDim_edge2D)
  ! contains external edges which mype needs;    
  ! myList_edge2D(1:myDim_edge2D) contains owned edges +
  ! shared edges which mype updates

end subroutine mymesh
!=================================================================
#ifndef FVOM_INIT
subroutine status_check
use g_config
use g_parsup
implicit none
integer :: res
res=0
call MPI_Allreduce (pe_status, res, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, MPIerr)
if (res /= 0 ) then
    if (mype==0) write(*,*) 'Something Broke. Flushing and stopping...'
!!! a restart file must be written here !!!
    call par_ex(1)
endif
end subroutine status_check
#endif
