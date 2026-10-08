program mpi_knem_stale_repro
    ! Standalone reproducer for LESSONS L61: Open MPI 4.1.2 (Levante spack build), pml ob1 +
    ! btl vader with the knem single-copy mechanism delivers the payload of an EARLIER message
    ! when the send buffer's heap pages were returned to the kernel (glibc trims the heap top
    ! with brk) and re-grown between two messages at the same virtual address.
    !
    ! Correct MPI only: contiguous allocatable buffers passed by their first element, both
    ! requests completed by MPI_Waitall before the buffers are freed, one message per pair and
    ! call. Each call allocates its buffers (as mod_halo's core_* do), fills the send buffer
    ! with a call-unique value and checks every received word.
    !
    ! Run by tools/run_mpi_knem_stale_repro.sh, which forces the glibc trim with
    ! MALLOC_TRIM_THRESHOLD_=0 MALLOC_TOP_PAD_=0 (free() then shrinks the heap on every call).
    ! Exit code 1 if any rank received a wrong word.
    use mpi
    implicit none
    integer, parameter :: NCALLS = 400, NL = 69       ! 69 x 60..90 doubles = 33..50 KB
    integer :: ierr, mype, npes, peer, icall, ncol, nbad, nbad_all, req(2)
    integer :: ncols(4) = (/ 60, 75, 60, 90 /)
    real(8), allocatable :: sbuf(:,:), rbuf(:,:)
    real(8) :: expect

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)
    if (npes /= 2) then
        if (mype == 0) write(*,'(a)') 'mpi_knem_stale_repro: run with exactly 2 ranks'
        call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    end if
    peer = 1 - mype
    nbad = 0
    do icall = 1, NCALLS
        ncol = ncols(mod(icall-1, 4) + 1)
        allocate(rbuf(NL, ncol), sbuf(NL, ncol))
        call MPI_Irecv(rbuf(1,1), NL*ncol, MPI_DOUBLE_PRECISION, peer, 1, MPI_COMM_WORLD, req(1), ierr)
        sbuf = real(icall, 8)*1000.0d0 + real(mype, 8)
        call MPI_Isend(sbuf(1,1), NL*ncol, MPI_DOUBLE_PRECISION, peer, 1, MPI_COMM_WORLD, req(2), ierr)
        call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
        expect = real(icall, 8)*1000.0d0 + real(peer, 8)
        if (any(rbuf /= expect)) then
            nbad = nbad + 1
            if (nbad <= 3) write(*,'(a,i0,a,i0,a,f10.1,a,f10.1,a,f10.1)') '  rank ', mype, &
                ' call ', icall, ': first word ', rbuf(1,1), ', last word ', rbuf(NL,ncol), &
                ', expected ', expect
        end if
        deallocate(rbuf, sbuf)
    end do
    call MPI_Allreduce(nbad, nbad_all, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (mype == 0) write(*,'(a,i0,a,i0,a)') 'mpi_knem_stale_repro: ', nbad_all, &
        ' wrong receives in ', 2*NCALLS, merge(' -> STALE PAYLOAD', ' -> ok           ', nbad_all > 0)
    call MPI_Finalize(ierr)
    if (nbad_all > 0) stop 1
end program mpi_knem_stale_repro
