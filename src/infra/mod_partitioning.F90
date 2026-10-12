module mod_partitioning
    ! Parallel decomposition setup (decision D7). Structure + dist_<NP>/ reading
    ! transcribed from FESOM2 v2.7.3 src/oce_mesh.F90:250-880 and
    ! src/gen_modules_partitioning.F90 (par_init/par_ex). Adds the in-memory
    ! 1-rank partition synthesis (D7): METIS cannot produce dist_1, so when
    ! npes==1 we build the identity local<->global map with no file read, keeping
    ! "serial == 1-rank MPI" as the bit-identity anchor.
    !
    ! Files are ASCII free-field (read(*,*)); ids on disk are 1-based (no shift).
    ! init_mpi_types (precompiled MPI_TYPE_INDEXED datatypes) is intentionally NOT
    ! built: the M0.5 halo packs manually via com_struct (broadcast-only), so the
    ! s/r_mpitype_* fields stay unallocated for v1.
    !
    ! Hybrid meshes: dist_N/format.out (version DIST_FORMAT_VERSION, written by fesom_ini)
    ! announces sections appended to the FESOM2 files -- my_list: eXDim_nod2D + the
    ! global ids of the double-halo nodes (appended to myList_nod2D); com_info:
    ! com_nod2D_full and com_edge2D in the com_struct layout. Without format.out the
    ! FESOM2 layout is read as before and dist_format stays 0.
    use mpi
    use mod_partit, only: t_partit, com_struct, MAX_NEIGHBOR_PARTITIONS
    implicit none
    private
    public :: par_init, par_ex, set_partition, read_mesh_dims, set_empty_comms, die

    integer, parameter :: DIST_FORMAT_VERSION = 1

contains

    subroutine par_init(partit)
        ! Initialise MPI (if needed) and fill the communicator / rank / size.
        type(t_partit), intent(inout) :: partit
        integer :: ierr
        logical :: inited
        call MPI_Initialized(inited, ierr)
        if (.not. inited) call MPI_Init(ierr)
        partit%MPI_COMM_FESOM = MPI_COMM_WORLD
        call MPI_Comm_size(partit%MPI_COMM_FESOM, partit%npes, ierr)
        call MPI_Comm_rank(partit%MPI_COMM_FESOM, partit%mype, ierr)
        if (partit%mype == 0) write(*,'(a,i0,a)') 'par_init: running on ', partit%npes, ' PE(s)'
    end subroutine par_init

    subroutine par_ex(comm, mype, abort)
        ! Finalise MPI (or abort). Transcribed from gen_modules_partitioning.F90:87.
        integer,           intent(in) :: comm
        integer,           intent(in) :: mype
        integer, optional, intent(in) :: abort
        integer :: error
        if (present(abort)) then
            if (mype == 0) write(*,*) 'Run finished unexpectedly!'
            call MPI_Abort(comm, 1, error)
        else
            call MPI_Barrier(comm, error)
            call MPI_Finalize(error)
        end if
    end subroutine par_ex

    subroutine read_mesh_dims(mesh_dir, nod2D, elem2D, edge2D, edge2D_in)
        ! Read the global counts from the mesh headers (nod2d.out, elem2d.out line 1;
        ! edgenum.out lines 1,2). Needed by the 1-rank synthesis.
        character(len=*), intent(in)  :: mesh_dir
        integer,          intent(out) :: nod2D, elem2D, edge2D, edge2D_in
        integer :: u, ios
        open(newunit=u, file=trim(mesh_dir)//'/nod2d.out', status='old', action='read', iostat=ios)
        read(u,*) nod2D; close(u)
        open(newunit=u, file=trim(mesh_dir)//'/elem2d.out', status='old', action='read', iostat=ios)
        read(u,*) elem2D; close(u)
        open(newunit=u, file=trim(mesh_dir)//'/edgenum.out', status='old', action='read', iostat=ios)
        read(u,*) edge2D
        read(u,*, iostat=ios) edge2D_in
        if (ios /= 0) edge2D_in = edge2D
        close(u)
    end subroutine read_mesh_dims

    subroutine set_partition(partit, mesh_dir)
        ! Build partit: synthesize the trivial partition for npes==1, else read
        ! dist_<npes>/ for this PE.
        type(t_partit),   intent(inout) :: partit
        character(len=*), intent(in)    :: mesh_dir
        if (partit%npes == 1) then
            call synthesize_1rank(partit, mesh_dir)
        else
            call read_dist_partition(partit, mesh_dir)
        end if
    end subroutine set_partition

    subroutine synthesize_1rank(partit, mesh_dir)
        ! Identity local<->global partition (D7): no dist file, no neighbors.
        type(t_partit),   intent(inout) :: partit
        character(len=*), intent(in)    :: mesh_dir
        integer :: nod2D, elem2D, edge2D, edge2D_in, i
        call read_mesh_dims(mesh_dir, nod2D, elem2D, edge2D, edge2D_in)

        partit%myDim_nod2D = nod2D; partit%eDim_nod2D = 0; partit%eXDim_nod2D = 0
        allocate(partit%myList_nod2D(nod2D))
        partit%myList_nod2D = [(i, i=1, nod2D)]

        partit%myDim_elem2D = elem2D; partit%eDim_elem2D = 0; partit%eXDim_elem2D = 0
        allocate(partit%myList_elem2D(elem2D))
        partit%myList_elem2D = [(i, i=1, elem2D)]

        partit%myDim_edge2D = edge2D; partit%eDim_edge2D = 0
        allocate(partit%myList_edge2D(edge2D))
        partit%myList_edge2D = [(i, i=1, edge2D)]

        allocate(partit%part(2)); partit%part = [1, nod2D+1]   ! cumulative, npes+1

        call set_empty_comms(partit)
        partit%pe_status = 0
    end subroutine synthesize_1rank

    subroutine set_empty_comms(partit)
        ! No neighbours: empty com_structs (rPEnum=sPEnum=0; halo is a no-op), the
        ! appended lists (com_nod2D_full, com_edge2D) present and empty. Also used by the
        ! analytic patches' 1-rank partition (mod_mesh_analytic synth_partit).
        type(t_partit), intent(inout) :: partit
        call clear_com(partit%com_nod2D)
        call clear_com(partit%com_elem2D)
        call clear_com(partit%com_elem2D_full)
        call clear_com(partit%com_nod2D_full)
        call clear_com(partit%com_edge2D)
        call clear_com(partit%com_elem2D_shared)
        if (allocated(partit%com_nod2D_full%rlist)) deallocate(partit%com_nod2D_full%rlist)
        if (allocated(partit%com_nod2D_full%slist)) deallocate(partit%com_nod2D_full%slist)
        if (allocated(partit%com_edge2D%rlist))     deallocate(partit%com_edge2D%rlist)
        if (allocated(partit%com_edge2D%slist))     deallocate(partit%com_edge2D%slist)
        allocate(partit%com_nod2D_full%rlist(0), partit%com_nod2D_full%slist(0))
        allocate(partit%com_edge2D%rlist(0), partit%com_edge2D%slist(0))
        partit%dist_format = DIST_FORMAT_VERSION          ! the new lists exist (empty)
    end subroutine set_empty_comms

    subroutine clear_com(com)
        type(com_struct), intent(inout) :: com
        com%rPEnum = 0; com%sPEnum = 0; com%nreq = 0
        com%rptr(1) = 1; com%sptr(1) = 1
    end subroutine clear_com

    subroutine read_dist_partition(partit, mesh_dir)
        ! Read dist_<npes>/ for this PE: rpart.out (counts), my_list, com_info.
        ! Transcribed from FESOM2 v2.7.3 oce_mesh.F90:262-879.
        type(t_partit),   intent(inout) :: partit
        character(len=*), intent(in)    :: mesh_dir
        character(len=:), allocatable :: dist_dir
        character(len=16) :: npes_string, mype_string
        character(len=512) :: fname
        integer :: u, ios, n, npes, mype
        integer, allocatable :: xlist(:)
        logical :: has_format

        npes = partit%npes; mype = partit%mype
        write(npes_string, '(I10)') npes
        write(mype_string, '(I5.5)') mype
        dist_dir = trim(mesh_dir)//'/dist_'//trim(adjustl(npes_string))//'/'

        ! ---- rpart.out: npes check + cumulative node-count boundaries (part) ----
        fname = dist_dir//'rpart.out'
        open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
        if (ios /= 0) call die('cannot open '//trim(fname), partit)
        read(u,*) n
        if (n /= npes) call die('rpart npes mismatch', partit)
        allocate(partit%part(npes+1))
        partit%part(1) = 1
        read(u,*) partit%part(2:npes+1)
        do n = 2, npes+1
            partit%part(n) = partit%part(n-1) + partit%part(n)
        end do
        close(u)

        ! ---- format.out: appended sections present (hybrid meshes) ----
        fname = dist_dir//'format.out'
        inquire(file=trim(fname), exist=has_format)
        partit%dist_format = 0
        if (has_format) then
            open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
            if (ios /= 0) call die('cannot open '//trim(fname), partit)
            read(u,*, iostat=ios) partit%dist_format
            close(u)
            if (ios /= 0 .or. partit%dist_format /= DIST_FORMAT_VERSION) &
                call die('unsupported version in '//trim(fname), partit)
        end if

        ! ---- my_list<mype>.out: local dims + global-id lists ----
        fname = dist_dir//'my_list'//trim(mype_string)//'.out'
        open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
        if (ios /= 0) call die('cannot open '//trim(fname), partit)
        read(u,*) n
        read(u,*) partit%myDim_nod2D
        read(u,*) partit%eDim_nod2D
        allocate(partit%myList_nod2D(partit%myDim_nod2D + partit%eDim_nod2D))
        read(u,*) partit%myList_nod2D
        read(u,*) partit%myDim_elem2D
        read(u,*) partit%eDim_elem2D
        read(u,*) partit%eXDim_elem2D
        allocate(partit%myList_elem2D(partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D))
        read(u,*) partit%myList_elem2D
        read(u,*) partit%myDim_edge2D
        read(u,*) partit%eDim_edge2D
        allocate(partit%myList_edge2D(partit%myDim_edge2D + partit%eDim_edge2D))
        read(u,*) partit%myList_edge2D
        partit%eXDim_nod2D = 0
        if (partit%dist_format >= 1) then
            ! double-halo nodes: local numbers myDim+eDim+1 .. myDim+eDim+eXDim_nod2D
            read(u,*, iostat=ios) partit%eXDim_nod2D
            if (ios /= 0) call die('no eXDim_nod2D in '//trim(fname), partit)
            allocate(xlist(partit%eXDim_nod2D))
            read(u,*, iostat=ios) xlist
            if (ios /= 0) call die('short eXDim_nod2D list in '//trim(fname), partit)
            partit%myList_nod2D = [partit%myList_nod2D, xlist]
            deallocate(xlist)
        end if
        close(u)

        ! ---- com_info<mype>.out: the three com_structs ----
        fname = dist_dir//'com_info'//trim(mype_string)//'.out'
        open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
        if (ios /= 0) call die('cannot open '//trim(fname), partit)
        read(u,*) n   ! header
        call read_com(u, partit%com_nod2D,       partit%eDim_nod2D, 'com_nod2D', fname, partit)
        call read_com(u, partit%com_elem2D,      partit%eDim_elem2D, 'com_elem2D', fname, partit)
        call read_com(u, partit%com_elem2D_full, partit%eDim_elem2D + partit%eXDim_elem2D, &
                      'com_elem2D_full', fname, partit)
        if (partit%dist_format >= 1) then
            call read_com(u, partit%com_nod2D_full, partit%eDim_nod2D + partit%eXDim_nod2D, &
                          'com_nod2D_full', fname, partit)
            call read_com(u, partit%com_edge2D, -1, 'com_edge2D', fname, partit)
            if (partit%com_nod2D_full%rptr(partit%com_nod2D_full%rPEnum+1) - 1 /= &
                partit%eDim_nod2D + partit%eXDim_nod2D) &
                call die('com_nod2D_full does not cover the eDim+eXDim halo nodes in '// &
                         trim(fname), partit)
        end if
        close(u)
        partit%pe_status = 0
    end subroutine read_dist_partition

    subroutine read_com(u, com, rlist_size, name, fname, partit)
        ! One com_struct from an open com_info unit (oce_mesh.F90:804-828 pattern).
        ! rlist_size < 0: the receive-list length is taken from rptr (com_edge2D).
        ! rPE/rptr/sPE/sptr are fixed-size: stop on more than MAX_NEIGHBOR_PARTITIONS PEs.
        integer,          intent(in)    :: u
        type(com_struct), intent(inout) :: com
        integer,          intent(in)    :: rlist_size
        character(len=*), intent(in)    :: name, fname
        type(t_partit),   intent(in)    :: partit
        integer :: n
        read(u,*) com%rPEnum
        call check_pe_count(com%rPEnum, 'receive')
        read(u,*) com%rPE(1:com%rPEnum)
        read(u,*) com%rptr(1:com%rPEnum+1)
        n = rlist_size
        if (n < 0) n = com%rptr(com%rPEnum+1) - 1
        allocate(com%rlist(n))
        read(u,*) com%rlist
        read(u,*) com%sPEnum
        call check_pe_count(com%sPEnum, 'send')
        read(u,*) com%sPE(1:com%sPEnum)
        read(u,*) com%sptr(1:com%sPEnum+1)
        n = com%sptr(com%sPEnum+1) - 1
        allocate(com%slist(n))
        read(u,*) com%slist
    contains
        subroutine check_pe_count(npe, side)
            integer,          intent(in) :: npe
            character(len=*), intent(in) :: side
            character(len=16) :: cnum, cmax
            if (npe >= 0 .and. npe <= MAX_NEIGHBOR_PARTITIONS) return
            write(cnum,'(i0)') npe
            write(cmax,'(i0)') MAX_NEIGHBOR_PARTITIONS
            call die(name//': '//trim(cnum)//' '//side//' PEs, outside 0..MAX_NEIGHBOR_PARTITIONS='// &
                     trim(cmax)//' in '//trim(fname)//' (increase MAX_NEIGHBOR_PARTITIONS in '// &
                     'src/types/mod_partit.F90 and in fesom_ini, then recompile)', partit)
        end subroutine check_pe_count
    end subroutine read_com

    subroutine die(msg, partit)
        ! Print msg and abort every rank (MPI_Abort): a failing rank must not leave the
        ! others waiting in a collective. Public for the setup-time list builders
        ! (mod_shared_elem).
        character(len=*), intent(in) :: msg
        type(t_partit),   intent(in) :: partit
        write(*,'(a)') 'set_partition: '//msg
        call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
    end subroutine die

end module mod_partitioning
