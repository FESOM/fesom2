!==========================================================
module g_PARSUP
USE o_PARAM
! Variables to organize parallel work  
implicit none
save

#ifdef PETSC
#include "finclude/petsc.h"
#else
 include 'mpif.h'
#endif

 integer                                :: MPI_COMM_FESOM
 integer, parameter   :: MAX_LAENDERECK=8
 integer, parameter   :: MAX_NEIGHBOR_PARTITIONS=32
  type com_struct
     integer                                       :: rPEnum ! the number of PE I receive info from 
     integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: rPE    ! their list
     integer, dimension(MAX_NEIGHBOR_PARTITIONS+1) :: rptr   ! allocatables to the list of nodes
     integer, dimension(:), allocatable            :: rlist  ! the list of nodes
     integer                                       :: sPEnum ! send part 
     integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: sPE
     integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: sptr
     integer, dimension(:), allocatable            :: slist
     integer, dimension(:), allocatable            :: req    ! request for MPI_Wait
     integer                                       :: nreq   ! number of requests for MPI_Wait
                                                             ! (to combine halo exchange of several fields)
  end type com_struct

  type(com_struct)   :: com_nod2D
!!$  type(com_struct)   :: com_edge2D
  type(com_struct), target :: com_elem2D
  type(com_struct), target :: com_elem2D_full
 
  ! MPI Datatypes for interface exchange

  ! Edge fields (2D)
  integer, allocatable       :: s_mpitype_edge2D(:),         r_mpitype_edge2D(:)   

  ! Element fields (2D; 2D integer; 3D with nl-1 or nl levels, 1 - 4 values)
  !                 small halo and / or full halo
  integer, allocatable, target :: s_mpitype_elem2D(:,:),       r_mpitype_elem2D(:,:) 
  integer, allocatable         :: s_mpitype_elem2D_full_i(:),  r_mpitype_elem2D_full_i(:) 
  integer, allocatable, target :: s_mpitype_elem2D_full(:,:),  r_mpitype_elem2D_full(:,:) 
  integer, allocatable, target :: s_mpitype_elem3D(:,:,:),     r_mpitype_elem3D(:,:,:) 
  integer, allocatable, target :: s_mpitype_elem3D_full(:,:,:),r_mpitype_elem3D_full(:,:,:) 

  ! Nodal fields (2D; 2D integer; 3D with nl-1 or nl levels, one, two, or three values)
  integer, allocatable       :: s_mpitype_nod2D(:),     r_mpitype_nod2D(:) 
  integer, allocatable       :: s_mpitype_nod2D_i(:),   r_mpitype_nod2D_i(:)
  integer, allocatable       :: s_mpitype_nod3D(:,:,:), r_mpitype_nod3D(:,:,:) 

  ! general MPI part
  integer            :: MPIERR
  integer            :: npes
  integer            :: mype
  integer            :: maxPEnum=100
  integer, allocatable, dimension(:)  :: part

  ! Mesh partition
  integer                             :: myDim_nod2D, eDim_nod2D
  integer, allocatable, dimension(:)  :: myList_nod2D
  integer                             :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
  integer, allocatable, dimension(:)  :: myList_elem2D
  integer                             :: myDim_edge2D, eDim_edge2D
  integer, allocatable, dimension(:)  :: myList_edge2D

  integer :: pe_status = 0 ! if /=0 then something is wrong 

  integer, allocatable ::  remPtr_nod2D(:),  remList_nod2D(:)
  integer, allocatable ::  remPtr_elem2D(:), remList_elem2D(:)

  logical :: elem_full_flag  
!$OMP threadprivate(com_nod2D,com_elem2D,com_elem2D_full)
!$OMP threadprivate(mype)
!$OMP threadprivate(myDim_nod2D, eDim_nod2D, myList_nod2D)
!$OMP threadprivate(myDim_elem2D, eDim_elem2D, eXDim_elem2D, myList_elem2D)
!$OMP threadprivate(myDim_edge2D, eDim_edge2D, myList_edge2D)
  

contains
subroutine par_init    ! initializes MPI
  implicit none


  integer :: i

#ifndef __oasis
  call MPI_Comm_Size(MPI_COMM_WORLD,npes,i)
  call MPI_Comm_Rank(MPI_COMM_WORLD,mype,i)
  MPI_COMM_FESOM=MPI_COMM_WORLD
#else
  call MPI_Comm_Size(MPI_COMM_FESOM,npes,i)
  call MPI_Comm_Rank(MPI_COMM_FESOM,mype,i)
#endif  

  if(mype==0) then
  write(*,*) 'MPI has been initialized'
  write(*, *) 'Running on ', npes, ' PEs'
  end if
end subroutine par_init
!=================================================================
subroutine par_ex(abort)       ! finalizes MPI
#if defined (__oasis)
  use mod_prism 
#endif
  implicit none
  integer,optional :: abort

#ifndef __oasis
  if (present(abort)) then
     if (mype==0) write(*,*) 'Run finished unexpectedly!'
     call MPI_ABORT( MPI_COMM_FESOM, 1 )
  else
     call  MPI_Barrier(MPI_COMM_FESOM,MPIerr)
     call  MPI_Finalize(MPIerr)
  endif
#else
#ifndef __oifs
!for coupling with ECHAM only
  if (.not. present(abort)) then
     if (mype==0) print *, 'FESOM calls MPI_Barrier before calling prism_terminate'
     call  MPI_Barrier(MPI_COMM_WORLD, MPIerr)
  end if
#endif
  call prism_terminate_proto(MPIerr)
  if (mype==0) print *, 'FESOM calls MPI_Barrier before calling MPI_Finalize'
  call  MPI_Barrier(MPI_COMM_WORLD, MPIerr)
  
  if (mype==0) print *, 'FESOM calls MPI_Finalize'
  call MPI_Finalize(MPIerr)
#endif
   if (mype==0) print *, 'fesom should stop with exit status = 0'
end subroutine par_ex
!=================================================================
subroutine set_par_support_ini

  use iso_c_binding, only: idx_t=>C_INT32_T
  use o_MESH
  use g_config
  implicit none

  integer         :: n, j, k, nini, nend, ierr
  integer(idx_t)  :: np(10)

  interface 
     subroutine partit(n,ptr,adj,wgt,np,part) bind(C)
       use iso_c_binding, only: idx_t=>C_INT32_T
       integer(idx_t), intent(in)  :: n, ptr(*), adj(*), wgt(*), np(*)
       integer(idx_t), intent(out) :: part(*)
     end subroutine partit
  end interface

  ! Construct partitioning vector
  if (n_levels<1 .OR. n_levels>10) then
     print *,'Number of hierarchic partition levels is out of range [1-10]! Aborting...'
     call MPI_ABORT( MPI_COMM_FESOM, 1 )
  end if

  np(:) = n_part(:)               ! Number of partitions on each hierarchy level
  if (n_part(1) == 0) then        ! Backward compatibility case: Take the number of 
     np(1) = npes                 ! partitions from the number of MPI processes
     n_levels = 1
  end if
  if (n_levels < 10) then         ! 0 is an indicator of the last hierarchy level
     np(n_levels+1) = 0 
  end if

  allocate(part(nod2D))
  part=0

  npes = PRODUCT(np(1:n_levels))
  if(npes<2) then
     print *,'Total number of parallel partitions is less than one! Aborting...'
     stop
  end if
  
  write(*,*) 'Calling partit'
  call partit(ssh_stiff%dim, ssh_stiff%rowptr, ssh_stiff%colind, &
       nlevels_nod2D, np, part)

!!$  call check_partitioning

  write(*,*) 'Partitioning is done.'

! The stiffness matrix is no longer needed. 
  deallocate(ssh_stiff%rowptr)
  deallocate(ssh_stiff%colind)
        
  !NR No longer needed - last use was as weight for partitioning
  deallocate(nlevels_nod2D)
end subroutine set_par_support_ini

!=======================================================================

subroutine check_partitioning

  ! In general, METIS 5 has several advantages compared to METIS 4, e.g.,
  !   * neighbouring tasks get neighbouring partitions (important for multicore computers!)
  !   * lower maximum of weights per partition (better load balancing)
  !   * lower memory demand
  !
  ! BUT: there might be outliers, single nodes connected to their partition by
  !      only one edge or even completely isolated. This spoils everything :-(
  !
  ! This routine checks for isolated nodes and moves them to an adjacent partition,
  ! trying not to spoil the load balance.

  use o_MESH
  integer :: i, j, k, n, n_iso, n_iter, is, ie, kmax, np
  integer :: nod_per_partition(2,0:npes-1)
  integer :: max_nod_per_part(2), min_nod_per_part(2)
  integer :: average_nod_per_part(2), node_neighb_part(100)
  logical :: already_counted, found_part

  integer :: max_adjacent_nodes
  integer, allocatable :: ne_part(:), ne_part_num(:), ne_part_load(:,:)

  ! Check load balancing
  do i=0,npes-1
     nod_per_partition(1,i) = count(part(:) == i)
     nod_per_partition(2,i) = sum(nlevels_nod2D,part(:) == i)
  enddo

  min_nod_per_part(1) = minval( nod_per_partition(1,:))
  min_nod_per_part(2) = minval( nod_per_partition(2,:))

  max_nod_per_part(1) = maxval( nod_per_partition(1,:))
  max_nod_per_part(2) = maxval( nod_per_partition(2,:))

  average_nod_per_part(1) = nod2D / npes
  average_nod_per_part(2) = sum(nlevels_nod2D(:)) / npes

  ! Now check for isolated nodes (connect by one or even no edge to other
  ! nodes of its partition) and repair, if possible

  max_adjacent_nodes = maxval(ssh_stiff%rowptr(2:nod2D+1) - ssh_stiff%rowptr(1:nod2D))
  allocate(ne_part(max_adjacent_nodes), ne_part_num(max_adjacent_nodes), &
       ne_part_load(2,max_adjacent_nodes))

  isolated_nodes_check: do n_iter = 1, 10
     print *,' '
     print *,'Check for isolated nodes, new iteration ========'
     n_iso = 0
     do n=1,nod2D
        is = ssh_stiff%rowptr(n)
        ie = ssh_stiff%rowptr(n+1) -1

        node_neighb_part(1:ie-is) = part(ssh_stiff%colind(is:ie))
        if (count(node_neighb_part(1:ie-is) == part(n)) <= 1) then

           n_iso = n_iso+1
           print *,'Isolated node',n, 'in partition', part(n)
           print *,'Neighbouring nodes are in partitions',  node_neighb_part(1:ie-is)

           ! count the adjacent nodes of the other PEs

           np=1
           ne_part(1) = node_neighb_part(1)
           ne_part_num(1) = 1
           ne_part_load(1,1) = nod_per_partition(1,ne_part(1)) + 1
           ne_part_load(2,1) = nod_per_partition(2,ne_part(1)) + nlevels_nod2D(n)

           do i=1,ie-is
              if (node_neighb_part(i)==part(n)) cycle
              already_counted = .false.
              do k=1,np
                 if (node_neighb_part(i) == ne_part(k)) then
                    ne_part_num(k) = ne_part_num(k) + 1
                    already_counted = .true.
                    exit
                 endif
              enddo
              if (.not. already_counted) then
                 np = np+1
                 ne_part(np) = node_neighb_part(i)
                 ne_part_num(np) = 1
                 ne_part_load(1,np) = nod_per_partition(1,ne_part(np)) + 1
                 ne_part_load(2,np) = nod_per_partition(2,ne_part(np)) + nlevels_nod2D(n)
              endif
           enddo

           ! Now, check for two things: The load balance, and if 
           ! there is more than one node of that partition.
           ! Otherwise, it would become isolated again.

           ! Best choice would be the partition with most adjacent nodes (edgecut!)
           ! Choose, if it does not decrease the load balance. 
           !        (There might be two partitions with the same number of adjacent
           !         nodes. Don't care about this here)

           kmax = maxloc(ne_part_num(1:np),1)

           if (ne_part_num(kmax) <= 1) then
              print *,'Sorry, no chance to solve an isolated node problem'
              exit isolated_nodes_check
           endif

           if  (ne_part_load(1,kmax) <= max_nod_per_part(1) .and. &
                ne_part_load(2,kmax) <= max_nod_per_part(2) ) then
              k = kmax
           else
              ! Don't make it too compicated. Reject partitions that have only one
              ! adjacent node. Take the next not violating the load balance.
              found_part = .false.
              do k=1,np
                 if (ne_part_num(k)==1 .or. k==kmax) cycle

                 if  (ne_part_load(1,k) <= max_nod_per_part(1) .and. &
                      ne_part_load(2,k) <= max_nod_per_part(2) ) then

                    found_part = .true.
                    exit
                 endif
              enddo

              if (.not. found_part) then
                 ! Ok, don't think to much. Simply go for minimized edge cut.
                 k = kmax
              endif
           endif

           ! Adjust the load balancing

           nod_per_partition(1,ne_part(k)) = nod_per_partition(1,ne_part(k)) + 1 
           nod_per_partition(2,ne_part(k)) = nod_per_partition(2,ne_part(k)) + nlevels_nod2D(n)
           nod_per_partition(1,part(n))    = nod_per_partition(1,part(n)) - 1
           nod_per_partition(2,part(n))    = nod_per_partition(2,part(n)) - nlevels_nod2D(n)

           ! And, finally, move nod n to other partition        
           part(n) = ne_part(k)
           print *,'Node',n,'is moved to part',part(n)
        endif
     enddo

     if (n_iso==0) then
        print *,'No isolated nodes found'
        exit isolated_nodes_check 
     endif
     ! Check for isolated nodes again
  end do isolated_nodes_check

  deallocate(ne_part, ne_part_num, ne_part_load)

  if (n_iso > 0) then
     print *,'++++++++++++++++++++++++++++++++++++++++++++'
     print *,'+  WARNING: PARTITIONING LOOKS REALLY BAD  +'
     print *,'+  It was not possible to remove all       +'
     print *,'+  isolated nodes. Consider to rerun with  +'
     print *,'+  new metis random seed (see Makefile.in) +'
     print *,'++++++++++++++++++++++++++++++++++++++++++++'
  endif

  print *,'=== LOAD BALANCING ==='
  print *,'2D nodes: min, aver, max per part',min_nod_per_part(1), &
       average_nod_per_part(1),max_nod_per_part(1)

  write(*,"('2D nodes: percent min, aver, max ',f8.3,'%, 100%, ',f8.3,'%')") &
       100.*real(min_nod_per_part(1)) / real(average_nod_per_part(1)), &
       100.*real(max_nod_per_part(1)) / real(average_nod_per_part(1))

  print *,'3D nodes: Min, aver, max per part',min_nod_per_part(2), &
       average_nod_per_part(2),max_nod_per_part(2)
  write(*,"('3D nodes: percent min, aver, max ',f8.3,'%, 100%, ',f8.3,'%')") &
       100.*real(min_nod_per_part(2)) / real(average_nod_per_part(2)), &
       100.*real(max_nod_per_part(2)) / real(average_nod_per_part(2))

end subroutine check_partitioning
!=======================================================================

subroutine set_par_support
  use o_MESH
  implicit none

  integer   n, offset
  integer :: i, max_nb, nb, nini, nend, nl1, n_val
  integer, allocatable :: blocklen(:),     displace(:)
  integer, allocatable :: blocklen_tmp(:), displace_tmp(:)

  !
  ! In the distributed memory version, most of the job is already done 
  ! at the initialization phase and is taken into account in read_mesh
  ! routine. Here, MPI datatypes are built and buffers for MPI wait requests
  ! are allocated. 

   if (npes > 1) then

!================================================
! MPI REQUEST BUFFERS
!================================================
!!$      allocate(com_edge2D%req(          3*com_edge2D%rPEnum +      3*com_edge2D%sPEnum))
      allocate(com_nod2D%req(            3*com_nod2D%rPEnum +       3*com_nod2D%sPEnum))
      allocate(com_elem2D%req(          3*com_elem2D%rPEnum +      3*com_elem2D%sPEnum))
      allocate(com_elem2D_full%req(3*com_elem2D_full%rPEnum + 3*com_elem2D_full%sPEnum))

!================================================
! MPI DATATYPES
!================================================
      ! Build MPI Data types for halo exchange: Edges
!!$      allocate(r_mpitype_edge2D(com_edge2D%rPEnum))  ! 2D
!!$      allocate(s_mpitype_edge2D(com_edge2D%sPEnum))  

      ! Upper limit for the length of the local interface between the neighbor PEs 
!!$      max_nb = max(maxval(com_edge2D%rptr(2:com_edge2D%rPEnum+1) - com_edge2D%rptr(1:com_edge2D%rPEnum)), &
!!$                   maxval(com_edge2D%sptr(2:com_edge2D%sPEnum+1) - com_edge2D%sptr(1:com_edge2D%sPEnum)))

!!$      allocate(displace(max_nb),     blocklen(max_nb))
!!$
!!$      do n=1,com_edge2D%rPEnum
!!$         nb = 1
!!$         nini = com_edge2D%rptr(n)
!!$         nend = com_edge2D%rptr(n+1) - 1
!!$         displace(:) = 0
!!$         displace(1) = com_edge2D%rlist(nini) -1  ! C counting, start at 0
!!$         blocklen(:) = 1
!!$         do i=nini+1, nend
!!$            if (com_edge2D%rlist(i) /= com_edge2D%rlist(i-1) + 1) then  
!!$               ! New block
!!$               nb = nb+1
!!$               displace(nb) = com_edge2D%rlist(i) -1
!!$            else
!!$               blocklen(nb) = blocklen(nb)+1
!!$            endif
!!$         enddo
!!$         
!!$         call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_DOUBLE_PRECISION, r_mpitype_edge2D(n), MPIerr)
!!$
!!$         call MPI_TYPE_COMMIT(r_mpitype_edge2D(n),   MPIerr) 
!!$      enddo
!!$
!!$      do n=1,com_edge2D%sPEnum
!!$         nb = 1
!!$         nini = com_edge2D%sptr(n)
!!$         nend = com_edge2D%sptr(n+1) - 1
!!$         displace(:) = 0
!!$         displace(1) = com_edge2D%slist(nini) -1  ! C counting, start at 0
!!$         blocklen(:) = 1
!!$         do i=nini+1, nend
!!$            if (com_edge2D%slist(i) /= com_edge2D%slist(i-1) + 1) then  
!!$               ! New block
!!$               nb = nb+1
!!$               displace(nb) = com_edge2D%slist(i) -1
!!$            else
!!$               blocklen(nb) = blocklen(nb)+1
!!$            endif
!!$         enddo
!!$         
!!$         call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_DOUBLE_PRECISION, s_mpitype_edge2D(n), MPIerr)
!!$
!!$         call MPI_TYPE_COMMIT(s_mpitype_edge2D(n),   MPIerr) 
!!$
!!$      enddo
!!$
!!$      deallocate(displace, blocklen)


      ! Build MPI Data types for halo exchange: Elements
      allocate(r_mpitype_elem2D(com_elem2D%rPEnum,4))     ! 2D, small halo
      allocate(s_mpitype_elem2D(com_elem2D%sPEnum,4))
      allocate(r_mpitype_elem2D_full_i(com_elem2D_full%rPEnum))   ! 2D, wide halo, integer
      allocate(s_mpitype_elem2D_full_i(com_elem2D_full%sPEnum))

      allocate(r_mpitype_elem2D_full(com_elem2D_full%rPEnum,4))     ! 2D, wide halo 
      allocate(s_mpitype_elem2D_full(com_elem2D_full%sPEnum,4))

      allocate(r_mpitype_elem3D(com_elem2D%rPEnum, nl-1:nl,4))     ! 3D, small halo 
      allocate(s_mpitype_elem3D(com_elem2D%sPEnum, nl-1:nl,4))

      allocate(r_mpitype_elem3D_full(com_elem2D_full%rPEnum, nl-1:nl,4))     ! 3D, wide halo
      allocate(s_mpitype_elem3D_full(com_elem2D_full%sPEnum, nl-1:nl,4))
      
      
      ! Upper limit for the length of the local interface between the neighbor PEs 
      max_nb = max(  &
           maxval(com_elem2D%rptr(2:com_elem2D%rPEnum+1) - com_elem2D%rptr(1:com_elem2D%rPEnum)), &
           maxval(com_elem2D%sptr(2:com_elem2D%sPEnum+1) - com_elem2D%sptr(1:com_elem2D%sPEnum)), &
           maxval(com_elem2D_full%rptr(2:com_elem2D_full%rPEnum+1) - com_elem2D_full%rptr(1:com_elem2D_full%rPEnum)), &
           maxval(com_elem2D_full%sptr(2:com_elem2D_full%sPEnum+1) - com_elem2D_full%sptr(1:com_elem2D_full%sPEnum)))
      
      allocate(displace(max_nb),     blocklen(max_nb))
      allocate(displace_tmp(max_nb), blocklen_tmp(max_nb))

      
      do n=1,com_elem2D%rPEnum
         nb = 1
         nini = com_elem2D%rptr(n)
         nend = com_elem2D%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_elem2D%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_elem2D%rlist(i) /= com_elem2D%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_elem2D%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         DO n_val=1,4

            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val 
            displace_tmp(1:nb) = displace(1:nb)*n_val 

            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_DOUBLE_PRECISION, &
                 r_mpitype_elem2D(n,n_val), MPIerr)

            call MPI_TYPE_COMMIT(r_mpitype_elem2D(n,n_val), MPIerr) 

            DO nl1=nl-1, nl

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_DOUBLE_PRECISION, & 
                    r_mpitype_elem3D(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(r_mpitype_elem3D(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO  
      enddo

      do n=1,com_elem2D%sPEnum
         nb = 1
         nini = com_elem2D%sptr(n)
         nend = com_elem2D%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_elem2D%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_elem2D%slist(i) /= com_elem2D%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_elem2D%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
                  
         DO n_val=1,4

            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val 
            displace_tmp(1:nb) = displace(1:nb)*n_val 

            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_DOUBLE_PRECISION, &
                 s_mpitype_elem2D(n, n_val), MPIerr)

            call MPI_TYPE_COMMIT(s_mpitype_elem2D(n, n_val),   MPIerr) 
 
            DO nl1=nl-1, nl

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_DOUBLE_PRECISION, & 
                    s_mpitype_elem3D(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(s_mpitype_elem3D(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO  
      enddo
      
      do n=1,com_elem2D_full%rPEnum
         nb = 1
         nini = com_elem2D_full%rptr(n)
         nend = com_elem2D_full%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_elem2D_full%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_elem2D_full%rlist(i) /= com_elem2D_full%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_elem2D_full%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         call MPI_TYPE_INDEXED(nb, blocklen,displace,MPI_INTEGER, r_mpitype_elem2D_full_i(n),MPIerr)

         call MPI_TYPE_COMMIT(r_mpitype_elem2D_full_i(n), MPIerr)

         DO n_val=1,4

            call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_DOUBLE_PRECISION, &
                 r_mpitype_elem2D_full(n,n_val), MPIerr)
            call MPI_TYPE_COMMIT(r_mpitype_elem2D_full(n, n_val),   MPIerr)

            DO nl1=nl-1, nl

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_DOUBLE_PRECISION, & 
                    r_mpitype_elem3D_full(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(r_mpitype_elem3D_full(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO   
      enddo

      do n=1,com_elem2D_full%sPEnum
         nb = 1
         nini = com_elem2D_full%sptr(n)
         nend = com_elem2D_full%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_elem2D_full%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_elem2D_full%slist(i) /= com_elem2D_full%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_elem2D_full%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         call MPI_TYPE_INDEXED(nb, blocklen,displace,MPI_INTEGER, s_mpitype_elem2D_full_i(n), MPIerr)

         call MPI_TYPE_COMMIT(s_mpitype_elem2D_full_i(n), MPIerr)  
 
         DO n_val=1,4
            call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_DOUBLE_PRECISION, &
                 s_mpitype_elem2D_full(n,n_val),MPIerr)
            call MPI_TYPE_COMMIT(s_mpitype_elem2D_full(n,n_val),   MPIerr)
  
            DO nl1=nl-1, nl

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_DOUBLE_PRECISION, & 
                    s_mpitype_elem3D_full(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(s_mpitype_elem3D_full(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO
      enddo

      deallocate(displace,     blocklen)
      deallocate(displace_tmp, blocklen_tmp)


   ! Build MPI Data types for halo exchange: Nodes

      allocate(r_mpitype_nod2D(com_nod2D%rPEnum))     ! 2D
      allocate(s_mpitype_nod2D(com_nod2D%sPEnum))
      allocate(r_mpitype_nod2D_i(com_nod2D%rPEnum))   ! 2D integer
      allocate(s_mpitype_nod2D_i(com_nod2D%sPEnum))   

      allocate(r_mpitype_nod3D(com_nod2D%rPEnum,nl-1:nl,3))  ! 3D with nl-1 or nl layers, 1-3 values 
      allocate(s_mpitype_nod3D(com_nod2D%sPEnum,nl-1:nl,3))
  
      ! Upper limit for the length of the local interface between the neighbor PEs 
      max_nb = max(maxval(com_nod2D%rptr(2:com_nod2D%rPEnum+1) - com_nod2D%rptr(1:com_nod2D%rPEnum)), &
                   maxval(com_nod2D%sptr(2:com_nod2D%sPEnum+1) - com_nod2D%sptr(1:com_nod2D%sPEnum)))

      allocate(displace(max_nb),     blocklen(max_nb))
      allocate(displace_tmp(max_nb), blocklen_tmp(max_nb))

      do n=1,com_nod2D%rPEnum
         nb = 1
         nini = com_nod2D%rptr(n)
         nend = com_nod2D%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_nod2D%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_nod2D%rlist(i) /= com_nod2D%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_nod2D%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_DOUBLE_PRECISION, & 
              r_mpitype_nod2D(n),     MPIerr)

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_INTEGER, & 
              r_mpitype_nod2D_i(n),   MPIerr)

         call MPI_TYPE_COMMIT(r_mpitype_nod2D(n),     MPIerr)    
         call MPI_TYPE_COMMIT(r_mpitype_nod2D_i(n),   MPIerr)     

         DO nl1=nl-1, nl
            DO n_val=1,3

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_DOUBLE_PRECISION, & 
                    r_mpitype_nod3D(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(r_mpitype_nod3D(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO
      enddo

      do n=1,com_nod2D%sPEnum
         nb = 1
         nini = com_nod2D%sptr(n)
         nend = com_nod2D%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_nod2D%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_nod2D%slist(i) /= com_nod2D%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_nod2D%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_DOUBLE_PRECISION, & 
              s_mpitype_nod2D(n),     MPIerr)

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_INTEGER, & 
              s_mpitype_nod2D_i(n),   MPIerr)

         call MPI_TYPE_COMMIT(s_mpitype_nod2D(n),     MPIerr)    
         call MPI_TYPE_COMMIT(s_mpitype_nod2D_i(n),   MPIerr)     

         DO nl1=nl-1, nl
            DO n_val=1,3

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_DOUBLE_PRECISION, & 
                    s_mpitype_nod3D(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(s_mpitype_nod3D(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO
      enddo

      deallocate(blocklen,     displace)
      deallocate(blocklen_tmp, displace_tmp)

   endif

   call init_gatherLists

  if(mype==0) write(*,*) 'Communication arrays are set' 

end subroutine set_par_support


!===================================================================
subroutine init_gatherLists

  use o_MESH
  implicit none
  
  integer :: n2D, e2D, sum_loc_elem2D
  integer :: n, estart, nstart

  if (mype==0) then

     if (npes > 1) then

        allocate(remPtr_nod2D(npes))
        allocate(remPtr_elem2D(npes))

        remPtr_nod2D(1) = 1
        remPtr_elem2D(1) = 1
        
        do n=1, npes-1
           call MPI_RECV(n2D, 1, MPI_INTEGER, n, 0, MPI_COMM_FESOM, MPI_STATUS_IGNORE, MPIerr )
           call MPI_RECV(e2D, 1, MPI_INTEGER, n, 1, MPI_COMM_FESOM, MPI_STATUS_IGNORE, MPIerr )

           remPtr_nod2D(n+1)  = remPtr_nod2D(n)  + n2D
           remPtr_elem2D(n+1) = remPtr_elem2D(n) + e2D 
        enddo



        allocate(remList_nod2D(remPtr_nod2D(npes)))   ! this should be nod2D - myDim_nod2D
        allocate(remList_elem2D(remPtr_elem2D(npes))) ! this is > elem2D, because the elements overlap.
                                                      ! Consider optimization: avoid multiple communication
                                                      ! of the same elem from different PEs.

        do n=1, npes-1
           nstart = remPtr_nod2D(n)
           n2D    = remPtr_nod2D(n+1) - remPtr_nod2D(n) 
           call MPI_RECV(remList_nod2D(nstart), n2D, MPI_INTEGER, n, 2, MPI_COMM_FESOM, &
                                                           MPI_STATUS_IGNORE, MPIerr ) 
           estart = remPtr_elem2D(n)
           e2D    = remPtr_elem2D(n+1) - remPtr_elem2D(n)
           call MPI_RECV(remList_elem2D(estart),e2D, MPI_INTEGER, n, 3, MPI_COMM_FESOM, &
                                                           MPI_STATUS_IGNORE, MPIerr ) 

        enddo
     end if
  else

     call MPI_SEND(myDim_nod2D,   1,            MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
     call MPI_SEND(myDim_elem2D,  1,            MPI_INTEGER, 0, 1, MPI_COMM_FESOM, MPIerr )
     call MPI_SEND(myList_nod2D,  myDim_nod2D,  MPI_INTEGER, 0, 2, MPI_COMM_FESOM, MPIerr )
     call MPI_SEND(myList_elem2D, myDim_elem2D, MPI_INTEGER, 0, 3, MPI_COMM_FESOM, MPIerr )

  endif

end subroutine init_gatherLists


end module g_PARSUP
