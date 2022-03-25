module mod_parsup
  interface
  subroutine par_ex(COMM, mype, abort)
     USE MOD_PARTIT
     implicit none
     integer,           intent(in)   :: COMM
     integer,           intent(in)   :: mype
     integer, optional, intent(in)   :: abort
  end subroutine
  end interface
end module mod_parsup

module par_support_interfaces
  interface
  subroutine par_init(partit)
     USE o_PARAM
     USE MOD_PARTIT
     USE MOD_PARSUP
     implicit none
     type(t_partit), intent(inout), target :: partit
  end subroutine

  subroutine init_mpi_types(partit, mesh)
     use MOD_MESH
     USE MOD_PARTIT
     USE MOD_PARSUP
     implicit none
     type(t_partit), intent(inout), target :: partit
     type(t_mesh),   intent(in), target :: mesh
  end subroutine

  subroutine init_gatherLists(partit)
     USE MOD_PARTIT
     USE MOD_PARSUP
     implicit none
     type(t_partit), intent(inout), target :: partit    
  end subroutine
  end interface
end module

subroutine par_init(partit)    ! initializes MPI
  USE o_PARAM
  USE MOD_PARTIT
  USE MOD_PARSUP
  implicit none
  type(t_partit), intent(inout), target :: partit
  integer                               :: i
  integer                               :: provided_mpi_thread_support_level
  character(:), allocatable             :: provided_mpi_thread_support_level_name

#if defined __oasis || defined  __ifsinterface
  ! use comm from coupler or ifs
#else
  partit%MPI_COMM_FESOM=MPI_COMM_WORLD ! use global comm if not coupled (e.g. no __oasis or __ifsinterface)
#endif  
  call MPI_Comm_Size(partit%MPI_COMM_FESOM,partit%npes,i)
  call MPI_Comm_Rank(partit%MPI_COMM_FESOM,partit%mype,i)
 

  if(partit%mype==0) then
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
    write(*, *) 'Running on                   ', partit%npes, ' PEs'
#if defined(_OPENMP)
    write(*, *) 'This is MPI/OpenMP run, with ', OMP_GET_MAX_THREADS(), ' threads per PE'
#endif
  end if
end subroutine par_init
!=================================================================
subroutine par_ex(COMM, mype, abort)       ! finalizes MPI
USE MOD_PARTIT
#ifndef __oifs
!For standalone and coupled ECHAM runs
#if defined (__oasis)
  use mod_prism 
#endif
  implicit none
  integer,           intent(in)   :: COMM
  integer,           intent(in)   :: mype
  integer, optional, intent(in)   :: abort
  integer                         :: error

#ifndef __oasis
  if (present(abort)) then
     if (partit%mype==0) write(*,*) 'Run finished unexpectedly!'
     call MPI_ABORT(COMM, 1 )
  else
     call  MPI_Barrier(COMM, error)
     call  MPI_Finalize(error)
  endif
#else
  if (.not. present(abort)) then
     if (partit%mype==0) print *, 'FESOM calls MPI_Barrier before calling prism_terminate'
     call  MPI_Barrier(MPI_COMM_WORLD, error)
  end if
  call prism_terminate_proto(error)
  if (partit%mype==0) print *, 'FESOM calls MPI_Barrier before calling MPI_Finalize'
  call  MPI_Barrier(MPI_COMM_WORLD, error)
  
  if (partit%mype==0) print *, 'FESOM calls MPI_Finalize'
  call MPI_Finalize(error)
#endif
  if (partit%mype==0) print *, 'fesom should stop with exit status = 0'
#endif
#if defined (__oifs)
!OIFS coupling doesnt call prism_terminate_proto and uses COMM instead of MPI_COMM_WORLD
  implicit none
  integer,           intent(in)   :: COMM
  integer,           intent(in)   :: mype
  integer, optional, intent(in)   :: abort
  integer                         :: error
  if (present(abort)) then
    if (mype==0) write(*,*) 'Run finished unexpectedly!'
    call MPI_ABORT(COMM, 1 )
  else
    call  MPI_Barrier(COMM, error)
    call  MPI_Finalize(error)
  endif
#endif

end subroutine par_ex
!=======================================================================
subroutine init_mpi_types(partit, mesh)
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  implicit none

  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  integer                          :: n, offset
  integer                          :: i, max_nb, nb, nini, nend, nl1, n_val
  integer, allocatable             :: blocklen(:),     displace(:)
  integer, allocatable             :: blocklen_tmp(:), displace_tmp(:)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  !
  ! In the distributed memory version, most of the job is already done 
  ! at the initialization phase and is taken into account in read_mesh
  ! routine. Here, MPI datatypes are built and buffers for MPI wait requests
  ! are allocated. 

   if (npes > 1) then

!================================================
! MPI REQUEST BUFFERS
!================================================
      if (.not. allocated(com_nod2D%req))        allocate(com_nod2D%req(            3*com_nod2D%rPEnum  + 3*com_nod2D%sPEnum))
      if (.not. allocated(com_elem2D%req))       allocate(com_elem2D%req(           3*com_elem2D%rPEnum + 3*com_elem2D%sPEnum))
      if (.not. allocated(com_elem2D_full%req))  allocate(com_elem2D_full%req(3*com_elem2D_full%rPEnum  + 3*com_elem2D_full%sPEnum))
!================================================
! MPI DATATYPES
!================================================
      ! Build MPI Data types for halo exchange: Elements
      allocate(partit%r_mpitype_elem2D(com_elem2D%rPEnum,4))     ! 2D, small halo
      allocate(partit%s_mpitype_elem2D(com_elem2D%sPEnum,4))
      allocate(partit%r_mpitype_elem2D_full_i(com_elem2D_full%rPEnum))   ! 2D, wide halo, integer
      allocate(partit%s_mpitype_elem2D_full_i(com_elem2D_full%sPEnum))
      allocate(partit%r_mpitype_elem2D_full(com_elem2D_full%rPEnum,4))     ! 2D, wide halo 
      allocate(partit%s_mpitype_elem2D_full(com_elem2D_full%sPEnum,4))
      allocate(partit%r_mpitype_elem3D(com_elem2D%rPEnum, nl-1:nl,4))     ! 3D, small halo 
      allocate(partit%s_mpitype_elem3D(com_elem2D%sPEnum, nl-1:nl,4))
      allocate(partit%r_mpitype_elem3D_full(com_elem2D_full%rPEnum, nl-1:nl,4))     ! 3D, wide halo
      allocate(partit%s_mpitype_elem3D_full(com_elem2D_full%sPEnum, nl-1:nl,4))
!after the allocation we just reassotiate ALL pointers again here
#include "associate_part_ass.h"
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

               ! r_mpitype_elem3D shape is e.g. 7,2,4 and its bounds 1:7,1:2,1:4 but the args n,nl1,n_val are 1,47,1 and thus OUT OF BOUNDS
               ! the second dimension of r_mpitype_elem3D is probably always 2 (from nl-1 to nl)
             if(.not. (all(lbound(r_mpitype_elem3D) .le. [n,nl1,n_val]) .and. all(ubound(r_mpitype_elem3D) .ge. [n,nl1,n_val])) ) then
               print *,"out of bounds error, lbound:",lbound(r_mpitype_elem3D), "indices:", n,nl1,n_val, "ubound:", ubound(r_mpitype_elem3D), __FILE__,__LINE__
               stop 1
             end if
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

               ! r_mpitype_elem3D shape is e.g. 7,2,4 and its bounds 1:7,1:2,1:4 but the args n,nl1,n_val are 1,47,1 and thus OUT OF BOUNDS
               ! the second dimension of r_mpitype_elem3D is probably always 2 (from nl-1 to nl)
             if(.not. (all(lbound(r_mpitype_elem3D_full) .le. [n,nl1,n_val]) .and. all(ubound(r_mpitype_elem3D_full) .ge. [n,nl1,n_val])) ) then
               print *,"out of bounds error, lbound:",lbound(r_mpitype_elem3D_full), "indices:", n,nl1,n_val, "ubound:", ubound(r_mpitype_elem3D_full), __FILE__,__LINE__
               stop 1
             end if
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

      allocate(partit%r_mpitype_nod2D(com_nod2D%rPEnum))     ! 2D
      allocate(partit%s_mpitype_nod2D(com_nod2D%sPEnum))
      allocate(partit%r_mpitype_nod2D_i(com_nod2D%rPEnum))   ! 2D integer
      allocate(partit%s_mpitype_nod2D_i(com_nod2D%sPEnum))   

      allocate(partit%r_mpitype_nod3D(com_nod2D%rPEnum,nl-1:nl,3))  ! 3D with nl-1 or nl layers, 1-3 values 
      allocate(partit%s_mpitype_nod3D(com_nod2D%sPEnum,nl-1:nl,3))
!after the allocation we just reassotiate ALL pointers again here
#include "associate_part_ass.h"
  
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
end subroutine init_mpi_types
!===================================================================
subroutine init_gatherLists(partit)
  USE MOD_PARTIT
  USE MOD_PARSUP
  implicit none
  type(t_partit), intent(inout), target :: partit    
  integer                               :: n2D, e2D, sum_loc_elem2D
  integer                               :: n, estart, nstart
#include "associate_part_def.h"
#include "associate_part_ass.h"
  if (mype==0) then

     if (npes > 1) then

        allocate(partit%remPtr_nod2D(npes))
        allocate(partit%remPtr_elem2D(npes))
!reassociate the pointers to the just allocated arrays
#include "associate_part_ass.h"
        remPtr_nod2D(1) = 1
        remPtr_elem2D(1) = 1
        
        do n=1, npes-1
           call MPI_RECV(n2D, 1, MPI_INTEGER, n, 0, MPI_COMM_FESOM, MPI_STATUS_IGNORE, MPIerr )
           call MPI_RECV(e2D, 1, MPI_INTEGER, n, 1, MPI_COMM_FESOM, MPI_STATUS_IGNORE, MPIerr )

           remPtr_nod2D(n+1)  = remPtr_nod2D(n)  + n2D
           remPtr_elem2D(n+1) = remPtr_elem2D(n) + e2D 
        enddo

        allocate(partit%remList_nod2D(remPtr_nod2D(npes)))   ! this should be nod2D - myDim_nod2D
        allocate(partit%remList_elem2D(remPtr_elem2D(npes))) ! this is > elem2D, because the elements overlap.
                                                      ! Consider optimization: avoid multiple communication
                                                      ! of the same elem from different PEs.
!reassociate the pointers to the just allocated arrays
#include "associate_part_ass.h"

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
!$OMP MASTER
#if defined(_OPENMP)
    allocate(partit%plock(myDim_elem2D+eDim_elem2D)) !allocate with maximum dimention (nELEM> nNODE) 
    do n=1, myDim_elem2D+eDim_elem2D
!experiments showd that OPENMP5 implementation of the lock (201811) is >10% more efficient
!make sure you use OPENMP v. 5.0
#if _OPENMP >= 201811
       call omp_init_lock_with_hint(partit%plock(n),omp_sync_hint_speculative+omp_sync_hint_uncontended)
#else
       call omp_init_lock(partit%plock(n))
#endif
    enddo
#endif
!$OMP END MASTER
end subroutine init_gatherLists
!===================================================================
subroutine status_check(partit)
USE MOD_PARTIT
USE MOD_PARSUP
implicit none
type(t_partit), intent(inout), target :: partit
integer                               :: res
res=0
call MPI_Allreduce (partit%pe_status, res, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
if (res /= 0 ) then
    if (partit%mype==0) write(*,*) 'Something Broke. Flushing and stopping...'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
endif
end subroutine status_check

