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
 integer, parameter   :: MAX_LAENDERECK=16
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
  integer provided_mpi_thread_support_level
  character(:), allocatable :: provided_mpi_thread_support_level_name

#ifndef __oasis
  call MPI_Comm_Size(MPI_COMM_WORLD,npes,i)
  call MPI_Comm_Rank(MPI_COMM_WORLD,mype,i)
  MPI_COMM_FESOM=MPI_COMM_WORLD
#else
  call MPI_Comm_Size(MPI_COMM_FESOM,npes,i)
  call MPI_Comm_Rank(MPI_COMM_FESOM,mype,i)
#endif  

  if(mype==0) then
    call MPI_Query_thread(provided_mpi_thread_support_level, i)
    if(provided_mpi_thread_support_level == MPI_THREAD_SINGLE) then
      provided_mpi_thread_support_level_name = "MPI_THREAD_SINGLE"
    else if(provided_mpi_thread_support_level == MPI_THREAD_FUNNELED) then
      provided_mpi_thread_support_level_name = "MPI_THREAD_FUNNELED"
    else if(provided_mpi_thread_support_level == MPI_THREAD_SERIALIZED) then
      provided_mpi_thread_support_level_name = "MPI_THREAD_SERIALIZED"
    else if(provided_mpi_thread_support_level == MPI_THREAD_MULTIPLE) then
      provided_mpi_thread_support_level_name = "MPI_THREAD_MULTIPLE"
    else
      provided_mpi_thread_support_level_name = "unknown"
    end if
    write(*,*) 'MPI has been initialized, provided MPI thread support level: ', &
         provided_mpi_thread_support_level_name,provided_mpi_thread_support_level
    write(*, *) 'Running on ', npes, ' PEs'
  end if
end subroutine par_init
!=================================================================
subroutine par_ex(abort)       ! finalizes MPI
#ifndef __oifs
!For standalone and coupled ECHAM runs
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
  if (.not. present(abort)) then
     if (mype==0) print *, 'FESOM calls MPI_Barrier before calling prism_terminate'
     call  MPI_Barrier(MPI_COMM_WORLD, MPIerr)
  end if
  call prism_terminate_proto(MPIerr)
  if (mype==0) print *, 'FESOM calls MPI_Barrier before calling MPI_Finalize'
  call  MPI_Barrier(MPI_COMM_WORLD, MPIerr)
  
  if (mype==0) print *, 'FESOM calls MPI_Finalize'
  call MPI_Finalize(MPIerr)
#endif
  if (mype==0) print *, 'fesom should stop with exit status = 0'
#endif
#if defined (__oifs)
!OIFS coupling doesnt call prism_terminate_proto and uses MPI_COMM_FESOM
  implicit none
  integer,optional :: abort
  if (present(abort)) then
	if (mype==0) write(*,*) 'Run finished unexpectedly!'
	call MPI_ABORT( MPI_COMM_FESOM, 1 )
  else
	call  MPI_Barrier(MPI_COMM_FESOM,MPIerr)
	call  MPI_Finalize(MPIerr)
  endif
#endif

end subroutine par_ex
!=======================================================================
subroutine set_par_support(mesh)
  use MOD_MESH
  implicit none

  type(t_mesh), intent(in)         , target :: mesh
  integer                          :: n, offset
  integer                          :: i, max_nb, nb, nini, nend, nl1, n_val
  integer, allocatable             :: blocklen(:),     displace(:)
  integer, allocatable             :: blocklen_tmp(:), displace_tmp(:)

#include "associate_mesh.h"  
  !
  ! In the distributed memory version, most of the job is already done 
  ! at the initialization phase and is taken into account in read_mesh
  ! routine. Here, MPI datatypes are built and buffers for MPI wait requests
  ! are allocated. 

   if (npes > 1) then

!================================================
! MPI REQUEST BUFFERS
!================================================
      allocate(com_nod2D%req(            3*com_nod2D%rPEnum +       3*com_nod2D%sPEnum))
      allocate(com_elem2D%req(          3*com_elem2D%rPEnum +      3*com_elem2D%sPEnum))
      allocate(com_elem2D_full%req(3*com_elem2D_full%rPEnum + 3*com_elem2D_full%sPEnum))

!================================================
! MPI DATATYPES
!================================================

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
