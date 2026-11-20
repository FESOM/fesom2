module par_support_module
    USE o_PARAM
    USE MOD_PARTIT
    USE iso_fortran_env, only: output_unit, error_unit
    USE MOD_MESH
#ifdef __MULTIO
  USE iom
  USE mpp_io
#endif

    implicit none

    private
    public :: par_init, par_abort, par_finalize, init_mpi_types, init_mpi_types_wet, init_mpi_types_fbin, &
              init_gatherLists, status_check

contains

subroutine par_init(partit)    ! initializes MPI

  implicit none
  type(t_partit), intent(inout), target :: partit
  integer                               :: i
  integer                               :: provided_mpi_thread_support_level
  character(:), allocatable             :: provided_mpi_thread_support_level_name
#if defined __oasis || defined  __ifsinterface || defined  __yac
  ! use comm from coupler or ifs
#else
  partit%MPI_COMM_FESOM=MPI_COMM_WORLD ! use global comm if not coupled (e.g. no __coupled or __ifsinterface)
#endif
  call MPI_Comm_Size(partit%MPI_COMM_FESOM,partit%npes,i)
  call MPI_Comm_Rank(partit%MPI_COMM_FESOM,partit%mype,i) 

  if(partit%mype==0) then
#if !defined(__PGI)
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
#endif
    write(*, *) 'Running on                   ', partit%npes, ' PEs'
#if defined(_OPENMP)
    write(*, *) 'This is MPI/OpenMP run, with ', OMP_GET_MAX_THREADS(), ' threads per PE'
#endif
  end if
end subroutine par_init
!=================================================================
subroutine par_abort(COMM, mype)           ! aborts the run (error exit)

  implicit none
  integer,           intent(in)   :: COMM
  integer,           intent(in)   :: mype
  integer                         :: error

! par_abort is the error exit: it always calls MPI_ABORT with a non-zero error
! code. The normal end of a run is par_finalize.
  write(*,*) 'Run finished unexpectedly! FESOM aborts on rank ', mype

  ! flush stdout and stderr so that error messages are visible before MPI_ABORT
  flush(output_unit)
  flush(error_unit)

#if !defined(__oasis) && !defined(__yac)
  ! standalone and ifsinterface
  call MPI_ABORT(COMM, 1, error)
#else
  ! coupled (oasis or yac): abort MPI_COMM_WORLD so that the coupled partners stop too
  call MPI_ABORT(MPI_COMM_WORLD, 1, error)
#endif

end subroutine par_abort
!=================================================================
subroutine par_finalize(COMM, mype)        ! finalizes MPI (normal end of the run)

! In case we are letting oasis orchestrate MPI, we need to shut down through
! oasis as well, thus we are including it here.
#if defined (__oasis)
#if defined (__oifs)
  !For OpenIFS coupled runs we use the new OASIS nameing scheme (oasis)
  use mod_oasis
#else
  !For ECHAM coupled runs we use the old OASIS nameing scheme (prism / prism_proto)
  use mod_prism
#endif
         ! oifs/echam
#endif
         ! oasis

  implicit none
  integer,           intent(in)   :: COMM
  integer,           intent(in)   :: mype
  integer                         :: error

  flush(output_unit)
  flush(error_unit)

! The model is finalized here: MPI_Barrier and MPI_Finalize, or the termination
! of the coupler. Error exits must call par_abort instead.
!---------------------------------------------------------------
#if !defined(__oasis) && !defined(__yac)
          ! TODO: this is where fesom standalone, and ifsinterface get to in case of normal exit.
          !1. barrier may hang so have to be careful
          !2. when using fesom as interface using finalize is bad here as there may 
          !   be other MPI tasks running in calling library like IFS, better 
          !   better practice in that case would be to free the communicator.
  call  MPI_Barrier(COMM, error)
  call  MPI_Finalize(error)

#else 
! From here on the two coupled options
!-------------------------------------
#if defined (__oifs)
  !For OpenIFS coupled runs we use the new OASIS nameing scheme (oasis)
#ifdef __oasis
  call oasis_terminate
#elif __yac
  CALL yac_ffinalize()
  CALL MPI_Finalize(error)
#endif
#elif defined(__yac)
 CALL yac_ffinalize()
 CALL MPI_Finalize(error)
#else
  !For ECHAM coupled runs we use the old OASIS nameing scheme (prism / prism_proto)
  if (mype==0) print *, 'FESOM calls MPI_Barrier before calling prism_terminate'
  call  MPI_Barrier(MPI_COMM_WORLD, error)
  call prism_terminate_proto(error)
  if (mype==0) print *, 'FESOM calls MPI_Barrier before calling MPI_Finalize'
  call  MPI_Barrier(MPI_COMM_WORLD, error)
  
  if (mype==0) print *, 'FESOM calls MPI_Finalize'
  call MPI_Finalize(error)
#endif /* oifs or echam */
#endif /* coupled: oasis&yac */

! Regardless of standalone, OpenIFS oder ECHAM coupling, if we reach to this point
! we should be fine shutting the whole model down
if (mype==0) print *, 'fesom should stop with exit status = 0'

end subroutine par_finalize
!=======================================================================
subroutine init_mpi_types(partit, mesh)
  implicit none

  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  integer                          :: n, offset
  integer                          :: i, max_nb, nb, nini, nend, nl1, n_val
  integer, allocatable             :: blocklen(:),     displace(:)
  integer, allocatable             :: blocklen_tmp(:), displace_tmp(:)

!PS npes            => partit%npes
!PS com_nod2D       => partit%com_nod2D
!PS com_elem2D      => partit%com_elem2D
!PS com_elem2D_full => partit%com_elem2D_full

  !
  ! In the distributed memory version, most of the job is already done 
  ! at the initialization phase and is taken into account in read_mesh
  ! routine. Here, MPI datatypes are built and buffers for MPI wait requests
  ! are allocated. 
  if (partit%npes > 1) then

!================================================
! MPI REQUEST BUFFERS
!================================================
      if (.not. allocated(partit%com_nod2D%req))        allocate(partit%com_nod2D%req(            3*partit%com_nod2D%rPEnum  + 3*partit%com_nod2D%sPEnum))
      if (.not. allocated(partit%com_elem2D%req))       allocate(partit%com_elem2D%req(           3*partit%com_elem2D%rPEnum + 3*partit%com_elem2D%sPEnum))
      if (.not. allocated(partit%com_elem2D_full%req))  allocate(partit%com_elem2D_full%req(3*partit%com_elem2D_full%rPEnum  + 3*partit%com_elem2D_full%sPEnum))

!================================================
! MPI DATATYPES
!================================================
      ! Build MPI Data types for halo exchange: Elements
      allocate(partit%r_mpitype_elem2D(partit%com_elem2D%rPEnum,4))     ! 2D, small halo
      allocate(partit%s_mpitype_elem2D(partit%com_elem2D%sPEnum,4))
      allocate(partit%r_mpitype_elem2D_full_i(partit%com_elem2D_full%rPEnum))   ! 2D, wide halo, integer
      allocate(partit%s_mpitype_elem2D_full_i(partit%com_elem2D_full%sPEnum))
      allocate(partit%r_mpitype_elem2D_full(partit%com_elem2D_full%rPEnum,4))     ! 2D, wide halo 
      allocate(partit%s_mpitype_elem2D_full(partit%com_elem2D_full%sPEnum,4))
      allocate(partit%r_mpitype_elem3D(partit%com_elem2D%rPEnum, mesh%nl-1:mesh%nl,4))     ! 3D, small halo 
      allocate(partit%s_mpitype_elem3D(partit%com_elem2D%sPEnum, mesh%nl-1:mesh%nl,4))
      allocate(partit%r_mpitype_elem3D_full(partit%com_elem2D_full%rPEnum, mesh%nl-1:mesh%nl,4))     ! 3D, wide halo
      allocate(partit%s_mpitype_elem3D_full(partit%com_elem2D_full%sPEnum, mesh%nl-1:mesh%nl,4))

!after the allocation we just reassotiate ALL pointers again here
      ! Upper limit for the length of the local interface between the neighbor PEs 
      max_nb = max(  &
           maxval(partit%com_elem2D%rptr(2:partit%com_elem2D%rPEnum+1) - partit%com_elem2D%rptr(1:partit%com_elem2D%rPEnum)), &
           maxval(partit%com_elem2D%sptr(2:partit%com_elem2D%sPEnum+1) - partit%com_elem2D%sptr(1:partit%com_elem2D%sPEnum)), &
           maxval(partit%com_elem2D_full%rptr(2:partit%com_elem2D_full%rPEnum+1) - partit%com_elem2D_full%rptr(1:partit%com_elem2D_full%rPEnum)), &
           maxval(partit%com_elem2D_full%sptr(2:partit%com_elem2D_full%sPEnum+1) - partit%com_elem2D_full%sptr(1:partit%com_elem2D_full%sPEnum)))
      
      allocate(displace(max_nb),     blocklen(max_nb))
      allocate(displace_tmp(max_nb), blocklen_tmp(max_nb))

      do n=1,partit%com_elem2D%rPEnum
         nb = 1
         nini = partit%com_elem2D%rptr(n)
         nend = partit%com_elem2D%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = partit%com_elem2D%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (partit%com_elem2D%rlist(i) /= partit%com_elem2D%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = partit%com_elem2D%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         DO n_val=1,4

            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val 
            displace_tmp(1:nb) = displace(1:nb)*n_val 

            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
                 partit%r_mpitype_elem2D(n,n_val), partit%MPIerr)

            call MPI_TYPE_COMMIT(partit%r_mpitype_elem2D(n,n_val), partit%MPIerr) 

            DO nl1=mesh%nl-1, mesh%nl

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               ! r_mpitype_elem3D shape is e.g. 7,2,4 and its bounds 1:7,1:2,1:4 but the args n,nl1,n_val are 1,47,1 and thus OUT OF BOUNDS
               ! the second dimension of r_mpitype_elem3D is probably always 2 (from nl-1 to nl)
             if(.not. (all(lbound(partit%r_mpitype_elem3D) .le. [n,nl1,n_val]) .and. all(ubound(partit%r_mpitype_elem3D) .ge. [n,nl1,n_val])) ) then
               print *,"out of bounds error, lbound:",lbound(partit%r_mpitype_elem3D), "indices:", n,nl1,n_val, "ubound:",ubound(partit%r_mpitype_elem3D), __FILE__,__LINE__
               stop 1
             end if
               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, & 
                    partit%r_mpitype_elem3D(n,nl1,n_val),  partit%MPIerr)

               call MPI_TYPE_COMMIT(partit%r_mpitype_elem3D(n,nl1,n_val),  partit%MPIerr)  
            ENDDO
         ENDDO  
      enddo

      do n=1,partit%com_elem2D%sPEnum
         nb = 1
         nini = partit%com_elem2D%sptr(n)
         nend = partit%com_elem2D%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = partit%com_elem2D%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (partit%com_elem2D%slist(i) /= partit%com_elem2D%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = partit%com_elem2D%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
                  
         DO n_val=1,4

            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val 
            displace_tmp(1:nb) = displace(1:nb)*n_val 

            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
                 partit%s_mpitype_elem2D(n, n_val), partit%MPIerr)

            call MPI_TYPE_COMMIT(partit%s_mpitype_elem2D(n, n_val),   partit%MPIerr) 
 
            DO nl1=mesh%nl-1, mesh%nl

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, & 
                    partit%s_mpitype_elem3D(n,nl1,n_val),  partit%MPIerr)

               call MPI_TYPE_COMMIT(partit%s_mpitype_elem3D(n,nl1,n_val),  partit%MPIerr)  
            ENDDO
         ENDDO  
      enddo
      
      do n=1,partit%com_elem2D_full%rPEnum
         nb = 1
         nini = partit%com_elem2D_full%rptr(n)
         nend = partit%com_elem2D_full%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = partit%com_elem2D_full%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (partit%com_elem2D_full%rlist(i) /= partit%com_elem2D_full%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = partit%com_elem2D_full%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         call MPI_TYPE_INDEXED(nb, blocklen,displace,MPI_INTEGER, partit%r_mpitype_elem2D_full_i(n),partit%MPIerr)

         call MPI_TYPE_COMMIT(partit%r_mpitype_elem2D_full_i(n), partit%MPIerr)

         DO n_val=1,4

            call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_WP, &
                 partit%r_mpitype_elem2D_full(n,n_val), partit%MPIerr)
            call MPI_TYPE_COMMIT(partit%r_mpitype_elem2D_full(n, n_val),   partit%MPIerr)

            DO nl1=mesh%nl-1, mesh%nl

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               ! r_mpitype_elem3D shape is e.g. 7,2,4 and its bounds 1:7,1:2,1:4 but the args n,nl1,n_val are 1,47,1 and thus OUT OF BOUNDS
               ! the second dimension of r_mpitype_elem3D is probably always 2 (from nl-1 to nl)
             if(.not. (all(lbound(partit%r_mpitype_elem3D_full) .le. [n,nl1,n_val]) .and. all(ubound(partit%r_mpitype_elem3D_full) .ge. [n,nl1,n_val])) ) then
               print *,"out of bounds error, lbound:",lbound(partit%r_mpitype_elem3D_full), "indices:", n,nl1,n_val, "ubound:", ubound(partit%r_mpitype_elem3D_full), __FILE__,__LINE__
               stop 1
             end if
               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, & 
                    partit%r_mpitype_elem3D_full(n,nl1,n_val),  partit%MPIerr)

               call MPI_TYPE_COMMIT(partit%r_mpitype_elem3D_full(n,nl1,n_val),  partit%MPIerr)  
            ENDDO
         ENDDO   
      enddo

      do n=1,partit%com_elem2D_full%sPEnum
         nb = 1
         nini = partit%com_elem2D_full%sptr(n)
         nend = partit%com_elem2D_full%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = partit%com_elem2D_full%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (partit%com_elem2D_full%slist(i) /= partit%com_elem2D_full%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = partit%com_elem2D_full%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         call MPI_TYPE_INDEXED(nb, blocklen,displace,MPI_INTEGER, partit%s_mpitype_elem2D_full_i(n), partit%MPIerr)

         call MPI_TYPE_COMMIT(partit%s_mpitype_elem2D_full_i(n), partit%MPIerr)  
 
         DO n_val=1,4
            call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_WP, &
                 partit%s_mpitype_elem2D_full(n,n_val), partit%MPIerr)
            call MPI_TYPE_COMMIT(partit%s_mpitype_elem2D_full(n,n_val),   partit%MPIerr)
  
            DO nl1=mesh%nl-1, mesh%nl

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, & 
                    partit%s_mpitype_elem3D_full(n,nl1,n_val),  partit%MPIerr)

               call MPI_TYPE_COMMIT(partit%s_mpitype_elem3D_full(n,nl1,n_val),  partit%MPIerr)  
            ENDDO
         ENDDO
      enddo

      deallocate(displace,     blocklen)
      deallocate(displace_tmp, blocklen_tmp)


   ! Build MPI Data types for halo exchange: Nodes

      allocate(partit%r_mpitype_nod2D(partit%com_nod2D%rPEnum))     ! 2D
      allocate(partit%s_mpitype_nod2D(partit%com_nod2D%sPEnum))
      allocate(partit%r_mpitype_nod2D_i(partit%com_nod2D%rPEnum))   ! 2D integer
      allocate(partit%s_mpitype_nod2D_i(partit%com_nod2D%sPEnum))   

      allocate(partit%r_mpitype_nod3D(partit%com_nod2D%rPEnum,mesh%nl-1:mesh%nl,3))  ! 3D with nl-1 or nl layers, 1-3 values 
      allocate(partit%s_mpitype_nod3D(partit%com_nod2D%sPEnum,mesh%nl-1:mesh%nl,3))
#if defined(USE_SINGLE_PRECISION)
      allocate(partit%r_mpitype_nod2D_r8(partit%com_nod2D%rPEnum))
      allocate(partit%s_mpitype_nod2D_r8(partit%com_nod2D%sPEnum))
      allocate(partit%r_mpitype_nod3D_r8(partit%com_nod2D%rPEnum,mesh%nl-1:mesh%nl,3))
      allocate(partit%s_mpitype_nod3D_r8(partit%com_nod2D%sPEnum,mesh%nl-1:mesh%nl,3))
#endif

!after the allocation we just reassotiate ALL pointers again here
  
      ! Upper limit for the length of the local interface between the neighbor PEs 
      max_nb = max(maxval(partit%com_nod2D%rptr(2:partit%com_nod2D%rPEnum+1) - partit%com_nod2D%rptr(1:partit%com_nod2D%rPEnum)), &
                   maxval(partit%com_nod2D%sptr(2:partit%com_nod2D%sPEnum+1) - partit%com_nod2D%sptr(1:partit%com_nod2D%sPEnum)))

      allocate(displace(max_nb),     blocklen(max_nb))
      allocate(displace_tmp(max_nb), blocklen_tmp(max_nb))

      do n=1,partit%com_nod2D%rPEnum
         nb = 1
         nini = partit%com_nod2D%rptr(n)
         nend = partit%com_nod2D%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = partit%com_nod2D%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (partit%com_nod2D%rlist(i) /= partit%com_nod2D%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = partit%com_nod2D%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_WP, & 
              partit%r_mpitype_nod2D(n),     partit%MPIerr)

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_INTEGER, & 
              partit%r_mpitype_nod2D_i(n),   partit%MPIerr)

         call MPI_TYPE_COMMIT(partit%r_mpitype_nod2D(n),     partit%MPIerr)    
         call MPI_TYPE_COMMIT(partit%r_mpitype_nod2D_i(n),   partit%MPIerr)     
#if defined(USE_SINGLE_PRECISION)
         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_WP_FULL, &
              partit%r_mpitype_nod2D_r8(n),  partit%MPIerr)
         call MPI_TYPE_COMMIT(partit%r_mpitype_nod2D_r8(n),  partit%MPIerr)
#endif

         DO nl1=mesh%nl-1, mesh%nl
            DO n_val=1,3

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, & 
                    partit%r_mpitype_nod3D(n,nl1,n_val),  partit%MPIerr)

               call MPI_TYPE_COMMIT(partit%r_mpitype_nod3D(n,nl1,n_val),  partit%MPIerr)  
#if defined(USE_SINGLE_PRECISION)
               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP_FULL, &
                    partit%r_mpitype_nod3D_r8(n,nl1,n_val),  partit%MPIerr)
               call MPI_TYPE_COMMIT(partit%r_mpitype_nod3D_r8(n,nl1,n_val),  partit%MPIerr)
#endif
            ENDDO
         ENDDO
      enddo

      do n=1,partit%com_nod2D%sPEnum
         nb = 1
         nini = partit%com_nod2D%sptr(n)
         nend = partit%com_nod2D%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = partit%com_nod2D%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (partit%com_nod2D%slist(i) /= partit%com_nod2D%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = partit%com_nod2D%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_WP, & 
              partit%s_mpitype_nod2D(n),     partit%MPIerr)

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_INTEGER, & 
              partit%s_mpitype_nod2D_i(n),   partit%MPIerr)

         call MPI_TYPE_COMMIT(partit%s_mpitype_nod2D(n),     partit%MPIerr)    
         call MPI_TYPE_COMMIT(partit%s_mpitype_nod2D_i(n),   partit%MPIerr)     
#if defined(USE_SINGLE_PRECISION)
         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_WP_FULL, &
              partit%s_mpitype_nod2D_r8(n),  partit%MPIerr)
         call MPI_TYPE_COMMIT(partit%s_mpitype_nod2D_r8(n),  partit%MPIerr)
#endif

         DO nl1=mesh%nl-1, mesh%nl
            DO n_val=1,3

               blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nl1 
               displace_tmp(1:nb) = displace(1:nb)*n_val*nl1 

               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, & 
                    partit%s_mpitype_nod3D(n,nl1,n_val),  partit%MPIerr)

               call MPI_TYPE_COMMIT(partit%s_mpitype_nod3D(n,nl1,n_val),  partit%MPIerr)  
#if defined(USE_SINGLE_PRECISION)
               call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP_FULL, &
                    partit%s_mpitype_nod3D_r8(n,nl1,n_val),  partit%MPIerr)
               call MPI_TYPE_COMMIT(partit%s_mpitype_nod3D_r8(n,nl1,n_val),  partit%MPIerr)
#endif
            ENDDO
         ENDDO
      enddo

      deallocate(blocklen,     displace)
      deallocate(blocklen_tmp, displace_tmp)
   endif
end subroutine init_mpi_types

!=======================================================================
! Rebuild the 3D halo datatypes so that each halo point carries only its wet
! levels (ulevels..nlevels-1 for nl-1 layer arrays, ulevels..nlevels for nl
! interface arrays) instead of whole columns, and keep the level bounds of all
! local points for the packed device exchange. Called from mesh_setup after the
! level arrays are known (halo_wet_levels in namelist.config).
subroutine init_mpi_types_wet(partit, mesh)
  implicit none
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  integer          :: n, nl1, n_val, sz, ierr
  integer(8)       :: bytes(4), tot(4)

  partit%halo_wet = .true.
  partit%wet_nl   = mesh%nl
  allocate(partit%wet_ulev_nod(size(mesh%ulevels_nod2D)), partit%wet_nlev_nod(size(mesh%nlevels_nod2D)))
  allocate(partit%wet_ulev_elem(size(mesh%ulevels)),      partit%wet_nlev_elem(size(mesh%nlevels)))
  partit%wet_ulev_nod  = mesh%ulevels_nod2D
  partit%wet_nlev_nod  = mesh%nlevels_nod2D
  partit%wet_ulev_elem = mesh%ulevels
  partit%wet_nlev_elem = mesh%nlevels
  if (partit%npes <= 1) return

  bytes = 0
  do n = 1, partit%com_nod2D%rPEnum
     call MPI_TYPE_SIZE(partit%r_mpitype_nod3D(n,mesh%nl-1,1), sz, ierr); bytes(1) = bytes(1) + sz
  end do
  do n = 1, partit%com_elem2D%rPEnum
     call MPI_TYPE_SIZE(partit%r_mpitype_elem3D(n,mesh%nl-1,1), sz, ierr); bytes(3) = bytes(3) + sz
  end do
  do nl1 = mesh%nl-1, mesh%nl
     do n_val = 1, 4
        do n = 1, partit%com_elem2D%rPEnum
           call wet_type(partit%com_elem2D%rlist, partit%com_elem2D%rptr(n), partit%com_elem2D%rptr(n+1)-1, &
                mesh%ulevels, mesh%nlevels, nl1, n_val, MPI_WP, partit%r_mpitype_elem3D(n,nl1,n_val))
        end do
        do n = 1, partit%com_elem2D%sPEnum
           call wet_type(partit%com_elem2D%slist, partit%com_elem2D%sptr(n), partit%com_elem2D%sptr(n+1)-1, &
                mesh%ulevels, mesh%nlevels, nl1, n_val, MPI_WP, partit%s_mpitype_elem3D(n,nl1,n_val))
        end do
        do n = 1, partit%com_elem2D_full%rPEnum
           call wet_type(partit%com_elem2D_full%rlist, partit%com_elem2D_full%rptr(n), partit%com_elem2D_full%rptr(n+1)-1, &
                mesh%ulevels, mesh%nlevels, nl1, n_val, MPI_WP, partit%r_mpitype_elem3D_full(n,nl1,n_val))
        end do
        do n = 1, partit%com_elem2D_full%sPEnum
           call wet_type(partit%com_elem2D_full%slist, partit%com_elem2D_full%sptr(n), partit%com_elem2D_full%sptr(n+1)-1, &
                mesh%ulevels, mesh%nlevels, nl1, n_val, MPI_WP, partit%s_mpitype_elem3D_full(n,nl1,n_val))
        end do
        if (n_val > 3) cycle
        do n = 1, partit%com_nod2D%rPEnum
           call wet_type(partit%com_nod2D%rlist, partit%com_nod2D%rptr(n), partit%com_nod2D%rptr(n+1)-1, &
                mesh%ulevels_nod2D, mesh%nlevels_nod2D, nl1, n_val, MPI_WP, partit%r_mpitype_nod3D(n,nl1,n_val))
#if defined(USE_SINGLE_PRECISION)
           call wet_type(partit%com_nod2D%rlist, partit%com_nod2D%rptr(n), partit%com_nod2D%rptr(n+1)-1, &
                mesh%ulevels_nod2D, mesh%nlevels_nod2D, nl1, n_val, MPI_WP_FULL, partit%r_mpitype_nod3D_r8(n,nl1,n_val))
#endif
        end do
        do n = 1, partit%com_nod2D%sPEnum
           call wet_type(partit%com_nod2D%slist, partit%com_nod2D%sptr(n), partit%com_nod2D%sptr(n+1)-1, &
                mesh%ulevels_nod2D, mesh%nlevels_nod2D, nl1, n_val, MPI_WP, partit%s_mpitype_nod3D(n,nl1,n_val))
#if defined(USE_SINGLE_PRECISION)
           call wet_type(partit%com_nod2D%slist, partit%com_nod2D%sptr(n), partit%com_nod2D%sptr(n+1)-1, &
                mesh%ulevels_nod2D, mesh%nlevels_nod2D, nl1, n_val, MPI_WP_FULL, partit%s_mpitype_nod3D_r8(n,nl1,n_val))
#endif
        end do
     end do
  end do
  do n = 1, partit%com_nod2D%rPEnum
     call MPI_TYPE_SIZE(partit%r_mpitype_nod3D(n,mesh%nl-1,1), sz, ierr); bytes(2) = bytes(2) + sz
  end do
  do n = 1, partit%com_elem2D%rPEnum
     call MPI_TYPE_SIZE(partit%r_mpitype_elem3D(n,mesh%nl-1,1), sz, ierr); bytes(4) = bytes(4) + sz
  end do
  call MPI_REDUCE(bytes, tot, 4, MPI_INTEGER8, MPI_SUM, 0, partit%MPI_COMM_FESOM, ierr)
  if (partit%mype == 0) write(*,'(a,2(f7.3,a))') ' 3D halo exchange of wet levels only, bytes per field wet/full: nodes ', &
       real(tot(2))/max(1_8,tot(1)), ', elements ', real(tot(4))/max(1_8,tot(3))

contains

  subroutine wet_type(list, nini, nend, ulev, nlev, nl1, n_val, base, mtype)
    integer, intent(in)    :: list(:), nini, nend, ulev(:), nlev(:), nl1, n_val, base
    integer, intent(inout) :: mtype
    integer, allocatable   :: bl(:), dp(:)
    integer                :: i, p, top, bot, ierr

    allocate(bl(max(1, nend-nini+1)), dp(max(1, nend-nini+1)))
    do i = nini, nend
       p   = list(i)
       top = ulev(p)
       bot = min(nlev(p) - (mesh%nl - nl1), nl1)
       bl(i-nini+1) = max(0, bot-top+1)*n_val
       dp(i-nini+1) = ((p-1)*nl1 + top-1)*n_val
    end do
    call MPI_TYPE_FREE(mtype, ierr)
    call MPI_TYPE_INDEXED(nend-nini+1, bl, dp, base, mtype, ierr)
    call MPI_TYPE_COMMIT(mtype, ierr)
    deallocate(bl, dp)
  end subroutine wet_type
end subroutine init_mpi_types_wet



!===================================================================
! Build MPI Data types for halo exchange of spectral bin (fbin) fields.
! These are used by IDEMIX2 where the first array dimension is nfbin
! (number of spectral frequency bins) instead of nl (number of vertical levels).
! Must be called after init_mpi_types and before first spectral bin exchange.
subroutine init_mpi_types_fbin(nfbin, partit)
    implicit none

    integer,        intent(in)           :: nfbin
    type(t_partit), intent(inout), target :: partit
    integer                          :: n, n_val
    integer                          :: i, max_nb, nb, nini, nend
    integer, allocatable             :: blocklen(:),     displace(:)
    integer, allocatable             :: blocklen_tmp(:), displace_tmp(:)


    if (partit%npes <= 1) return

    ! Store nfbin for runtime checks in exchange routines
    partit%nfbin_mpi = nfbin

    !================================================================
    ! Allocate MPI datatype arrays for spectral bin exchange
    ! Full spectral column exchanged at once (mirroring standard nod3D/elem3D pattern)
    !================================================================
    allocate(partit%r_mpitype_nod3D_fbin(      partit%com_nod2D%rPEnum      , nfbin:nfbin, 3)) ! nodes, nfbin bins, 1-3 values
    allocate(partit%s_mpitype_nod3D_fbin(      partit%com_nod2D%sPEnum      , nfbin:nfbin, 3))
    allocate(partit%r_mpitype_elem2D_fbin(     partit%com_elem2D%rPEnum     , nfbin:nfbin   )) ! elems 2D small halo
    allocate(partit%s_mpitype_elem2D_fbin(     partit%com_elem2D%sPEnum     , nfbin:nfbin   ))
    allocate(partit%r_mpitype_elem2D_full_fbin(partit%com_elem2D_full%rPEnum, nfbin:nfbin   )) ! elems 2D full halo
    allocate(partit%s_mpitype_elem2D_full_fbin(partit%com_elem2D_full%sPEnum, nfbin:nfbin   ))
    allocate(partit%r_mpitype_elem3D_fbin(     partit%com_elem2D%rPEnum     , nfbin:nfbin, 4)) ! elems 3D small halo, 1-4 values
    allocate(partit%s_mpitype_elem3D_fbin(     partit%com_elem2D%sPEnum     , nfbin:nfbin, 4))
    allocate(partit%r_mpitype_elem3D_full_fbin(partit%com_elem2D_full%rPEnum, nfbin:nfbin, 4)) ! elems 3D full halo
    allocate(partit%s_mpitype_elem3D_full_fbin(partit%com_elem2D_full%sPEnum, nfbin:nfbin, 4))

    !================================================================
    ! Build MPI Data types for element fields (small halo)
    !================================================================
    max_nb = max(  &
                maxval(partit%com_elem2D%rptr(     2:partit%com_elem2D%rPEnum+1     ) - partit%com_elem2D%rptr(     1:partit%com_elem2D%rPEnum     )), &
                maxval(partit%com_elem2D%sptr(     2:partit%com_elem2D%sPEnum+1     ) - partit%com_elem2D%sptr(     1:partit%com_elem2D%sPEnum     )), &
                maxval(partit%com_elem2D_full%rptr(2:partit%com_elem2D_full%rPEnum+1) - partit%com_elem2D_full%rptr(1:partit%com_elem2D_full%rPEnum)), &
                maxval(partit%com_elem2D_full%sptr(2:partit%com_elem2D_full%sPEnum+1) - partit%com_elem2D_full%sptr(1:partit%com_elem2D_full%sPEnum)))

    allocate(displace(max_nb),     blocklen(max_nb))
    allocate(displace_tmp(max_nb), blocklen_tmp(max_nb))

    ! --- Receive types for elements (small halo) ---
    do n=1,partit%com_elem2D%rPEnum
        nb = 1
        nini = partit%com_elem2D%rptr(n)
        nend = partit%com_elem2D%rptr(n+1) - 1
        displace(:) = 0
        displace(1) = partit%com_elem2D%rlist(nini) -1
        blocklen(:) = 1
        do i=nini+1, nend
            if (partit%com_elem2D%rlist(i) /= partit%com_elem2D%rlist(i-1) + 1) then
            nb = nb+1
            displace(nb) = partit%com_elem2D%rlist(i) -1
            else
            blocklen(nb) = blocklen(nb)+1
            endif
        enddo
        ! 2D fbin type: full spectral column at once (nfbin per element)
        blocklen_tmp(1:nb) = blocklen(1:nb)*nfbin
        displace_tmp(1:nb) = displace(1:nb)*nfbin
        call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
            partit%r_mpitype_elem2D_fbin(n,nfbin), partit%MPIerr)
        call MPI_TYPE_COMMIT(partit%r_mpitype_elem2D_fbin(n,nfbin), partit%MPIerr)
        ! 3D fbin types: full spectral column, n_val values per bin
        DO n_val=1,4
            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nfbin
            displace_tmp(1:nb) = displace(1:nb)*n_val*nfbin
            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
                partit%r_mpitype_elem3D_fbin(n,nfbin,n_val), partit%MPIerr)
            call MPI_TYPE_COMMIT(partit%r_mpitype_elem3D_fbin(n,nfbin,n_val), partit%MPIerr)
        ENDDO
    enddo

    ! --- Send types for elements (small halo) ---
    do n=1,partit%com_elem2D%sPEnum
        nb = 1
        nini = partit%com_elem2D%sptr(n)
        nend = partit%com_elem2D%sptr(n+1) - 1
        displace(:) = 0
        displace(1) = partit%com_elem2D%slist(nini) -1
        blocklen(:) = 1
        do i=nini+1, nend
            if (partit%com_elem2D%slist(i) /= partit%com_elem2D%slist(i-1) + 1) then
            nb = nb+1
            displace(nb) = partit%com_elem2D%slist(i) -1
            else
            blocklen(nb) = blocklen(nb)+1
            endif
        enddo
        ! 2D fbin type: full spectral column at once
        blocklen_tmp(1:nb) = blocklen(1:nb)*nfbin
        displace_tmp(1:nb) = displace(1:nb)*nfbin
        call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
            partit%s_mpitype_elem2D_fbin(n,nfbin), partit%MPIerr)
        call MPI_TYPE_COMMIT(partit%s_mpitype_elem2D_fbin(n,nfbin), partit%MPIerr)
        ! 3D fbin types: full spectral column, n_val values per bin
        DO n_val=1,4
            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nfbin
            displace_tmp(1:nb) = displace(1:nb)*n_val*nfbin
            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
                partit%s_mpitype_elem3D_fbin(n,nfbin,n_val), partit%MPIerr)
            call MPI_TYPE_COMMIT(partit%s_mpitype_elem3D_fbin(n,nfbin,n_val), partit%MPIerr)
        ENDDO
    enddo

    !================================================================
    ! Build MPI Data types for element fields (full halo)
    !================================================================
    ! --- Receive types for elements (full halo) ---
    do n=1,partit%com_elem2D_full%rPEnum
        nb = 1
        nini = partit%com_elem2D_full%rptr(n)
        nend = partit%com_elem2D_full%rptr(n+1) - 1
        displace(:) = 0
        displace(1) = partit%com_elem2D_full%rlist(nini) -1
        blocklen(:) = 1
        do i=nini+1, nend
            if (partit%com_elem2D_full%rlist(i) /= partit%com_elem2D_full%rlist(i-1) + 1) then
            nb = nb+1
            displace(nb) = partit%com_elem2D_full%rlist(i) -1
            else
            blocklen(nb) = blocklen(nb)+1
            endif
        enddo
        ! 2D fbin type: full spectral column at once
        blocklen_tmp(1:nb) = blocklen(1:nb)*nfbin
        displace_tmp(1:nb) = displace(1:nb)*nfbin
        call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
            partit%r_mpitype_elem2D_full_fbin(n,nfbin), partit%MPIerr)
        call MPI_TYPE_COMMIT(partit%r_mpitype_elem2D_full_fbin(n,nfbin), partit%MPIerr)
        ! 3D fbin types: full spectral column, n_val values per bin
        DO n_val=1,4
            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nfbin
            displace_tmp(1:nb) = displace(1:nb)*n_val*nfbin
            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
                partit%r_mpitype_elem3D_full_fbin(n,nfbin,n_val), partit%MPIerr)
            call MPI_TYPE_COMMIT(partit%r_mpitype_elem3D_full_fbin(n,nfbin,n_val), partit%MPIerr)
        ENDDO
    enddo

    ! --- Send types for elements (full halo) ---
    do n=1,partit%com_elem2D_full%sPEnum
        nb = 1
        nini = partit%com_elem2D_full%sptr(n)
        nend = partit%com_elem2D_full%sptr(n+1) - 1
        displace(:) = 0
        displace(1) = partit%com_elem2D_full%slist(nini) -1
        blocklen(:) = 1
        do i=nini+1, nend
            if (partit%com_elem2D_full%slist(i) /= partit%com_elem2D_full%slist(i-1) + 1) then
            nb = nb+1
            displace(nb) = partit%com_elem2D_full%slist(i) -1
            else
            blocklen(nb) = blocklen(nb)+1
            endif
        enddo
        ! 2D fbin type: full spectral column at once
        blocklen_tmp(1:nb) = blocklen(1:nb)*nfbin
        displace_tmp(1:nb) = displace(1:nb)*nfbin
        call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
            partit%s_mpitype_elem2D_full_fbin(n,nfbin), partit%MPIerr)
        call MPI_TYPE_COMMIT(partit%s_mpitype_elem2D_full_fbin(n,nfbin), partit%MPIerr)
        ! 3D fbin types: full spectral column, n_val values per bin
        DO n_val=1,4
            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nfbin
            displace_tmp(1:nb) = displace(1:nb)*n_val*nfbin
            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
                partit%s_mpitype_elem3D_full_fbin(n,nfbin,n_val), partit%MPIerr)
            call MPI_TYPE_COMMIT(partit%s_mpitype_elem3D_full_fbin(n,nfbin,n_val), partit%MPIerr)
        ENDDO
    enddo

    deallocate(displace,     blocklen)
    deallocate(displace_tmp, blocklen_tmp)

    !================================================================
    ! Build MPI Data types for nodal fields
    !================================================================
    max_nb = max(maxval(partit%com_nod2D%rptr(2:partit%com_nod2D%rPEnum+1) - partit%com_nod2D%rptr(1:partit%com_nod2D%rPEnum)), &
                 maxval(partit%com_nod2D%sptr(2:partit%com_nod2D%sPEnum+1) - partit%com_nod2D%sptr(1:partit%com_nod2D%sPEnum)))

    allocate(displace(max_nb),     blocklen(max_nb))
    allocate(displace_tmp(max_nb), blocklen_tmp(max_nb))

    ! --- Receive types for nodes ---
    do n=1,partit%com_nod2D%rPEnum
        nb = 1
        nini = partit%com_nod2D%rptr(n)
        nend = partit%com_nod2D%rptr(n+1) - 1
        displace(:) = 0
        displace(1) = partit%com_nod2D%rlist(nini) -1
        blocklen(:) = 1
        do i=nini+1, nend
            if (partit%com_nod2D%rlist(i) /= partit%com_nod2D%rlist(i-1) + 1) then
            nb = nb+1
            displace(nb) = partit%com_nod2D%rlist(i) -1
            else
            blocklen(nb) = blocklen(nb)+1
            endif
        enddo
        ! Full spectral column: n_val*nfbin contiguous elements per node
        DO n_val=1,3
            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nfbin
            displace_tmp(1:nb) = displace(1:nb)*n_val*nfbin
            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
                partit%r_mpitype_nod3D_fbin(n,nfbin,n_val), partit%MPIerr)
            call MPI_TYPE_COMMIT(partit%r_mpitype_nod3D_fbin(n,nfbin,n_val), partit%MPIerr)
        ENDDO
    enddo

    ! --- Send types for nodes ---
    do n=1,partit%com_nod2D%sPEnum
        nb = 1
        nini = partit%com_nod2D%sptr(n)
        nend = partit%com_nod2D%sptr(n+1) - 1
        displace(:) = 0
        displace(1) = partit%com_nod2D%slist(nini) -1
        blocklen(:) = 1
        do i=nini+1, nend
            if (partit%com_nod2D%slist(i) /= partit%com_nod2D%slist(i-1) + 1) then
            nb = nb+1
            displace(nb) = partit%com_nod2D%slist(i) -1
            else
            blocklen(nb) = blocklen(nb)+1
            endif
        enddo
        ! Full spectral column: n_val*nfbin contiguous elements per node
        DO n_val=1,3
            blocklen_tmp(1:nb) = blocklen(1:nb)*n_val*nfbin
            displace_tmp(1:nb) = displace(1:nb)*n_val*nfbin
            call MPI_TYPE_INDEXED(nb, blocklen_tmp, displace_tmp, MPI_WP, &
                partit%s_mpitype_nod3D_fbin(n,nfbin,n_val), partit%MPIerr)
            call MPI_TYPE_COMMIT(partit%s_mpitype_nod3D_fbin(n,nfbin,n_val), partit%MPIerr)
        ENDDO
    enddo

    deallocate(blocklen,     displace)
    deallocate(blocklen_tmp, displace_tmp)

end subroutine init_mpi_types_fbin
!===================================================================
subroutine init_gatherLists(partit)
  implicit none
  type(t_partit), intent(inout), target :: partit    
  integer                               :: n2D, e2D, sum_loc_elem2D
  integer                               :: n, estart, nstart
  if (partit%mype==0) then

     ! also on a single rank: the I/O gathers index remPtr_nod2D(npes) unconditionally
     if (partit%npes >= 1) then

        allocate(partit%remPtr_nod2D(partit%npes))
        allocate(partit%remPtr_elem2D(partit%npes))
!reassociate the pointers to the just allocated arrays
        partit%remPtr_nod2D(1) = 1
        partit%remPtr_elem2D(1) = 1
        
        do n=1, partit%npes-1
           call MPI_RECV(n2D, 1, MPI_INTEGER, n, 0, partit%MPI_COMM_FESOM, MPI_STATUS_IGNORE, partit%MPIerr )
           call MPI_RECV(e2D, 1, MPI_INTEGER, n, 1, partit%MPI_COMM_FESOM, MPI_STATUS_IGNORE, partit%MPIerr )

           partit%remPtr_nod2D(n+1)  = partit%remPtr_nod2D(n)  + n2D
           partit%remPtr_elem2D(n+1) = partit%remPtr_elem2D(n) + e2D 
        enddo

        allocate(partit%remList_nod2D(partit%remPtr_nod2D(partit%npes)))   ! this should be nod2D - myDim_nod2D
        allocate(partit%remList_elem2D(partit%remPtr_elem2D(partit%npes))) ! this is > elem2D, because the elements overlap.
                                                      ! Consider optimization: avoid multiple communication
                                                      ! of the same elem from different PEs.
!reassociate the pointers to the just allocated arrays

        do n=1, partit%npes-1
           nstart = partit%remPtr_nod2D(n)
           n2D    = partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n) 
           call MPI_RECV(partit%remList_nod2D(nstart), n2D, MPI_INTEGER, n, 2, partit%MPI_COMM_FESOM, &
                                                           MPI_STATUS_IGNORE, partit%MPIerr ) 
           estart = partit%remPtr_elem2D(n)
           e2D    = partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n)
           call MPI_RECV(partit%remList_elem2D(estart),e2D, MPI_INTEGER, n, 3, partit%MPI_COMM_FESOM, &
                                                           MPI_STATUS_IGNORE, partit%MPIerr ) 

        enddo
     end if
  else

     call MPI_SEND(partit%myDim_nod2D,   1,            MPI_INTEGER, 0, 0, partit%MPI_COMM_FESOM, partit%MPIerr )
     call MPI_SEND(partit%myDim_elem2D,  1,            MPI_INTEGER, 0, 1, partit%MPI_COMM_FESOM, partit%MPIerr )
     call MPI_SEND(partit%myList_nod2D,  partit%myDim_nod2D,  MPI_INTEGER, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )
     call MPI_SEND(partit%myList_elem2D, partit%myDim_elem2D, MPI_INTEGER, 0, 3, partit%MPI_COMM_FESOM, partit%MPIerr )

  endif
end subroutine init_gatherLists
!===================================================================
subroutine status_check(partit)
implicit none
type(t_partit), intent(inout), target :: partit
integer                               :: res
res=0
call MPI_Allreduce (partit%pe_status, res, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
if (res /= 0 ) then
    if (partit%mype==0) write(*,*) 'Something Broke. Flushing and stopping...'
    call par_abort(partit%MPI_COMM_FESOM, partit%mype)
endif
end subroutine status_check

end module par_support_module
