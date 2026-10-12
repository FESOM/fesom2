module mod_partit
    ! Parallel-decomposition type (decision D5: STATIC). Structure from the
    ! tracer_dwarf MOD_PARTIT (== FESOM2 v2.7.3), trimmed to v1 scope: iceberg
    ! async communicators and OpenMP locks are dropped (YAGNI). com_struct is the
    ! language-agnostic halo contract (D7).
    use mod_binary_arrays, only: write_bin_array, read_bin_array, &
                                 write1d_int_static, read1d_int_static
    use mpi
    implicit none
    save

    integer, parameter :: MAX_NEIGHBOR_PARTITIONS = 32

    type com_struct
        integer :: rPEnum = 0                                    ! # PEs I receive from
        integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: rPE     ! their list
        integer, dimension(MAX_NEIGHBOR_PARTITIONS+1) :: rptr    ! ptrs into rlist (1-based cumulative)
        integer, dimension(:), allocatable            :: rlist   ! received node list
        integer :: sPEnum = 0                                    ! send part
        integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: sPE
        integer, dimension(MAX_NEIGHBOR_PARTITIONS+1) :: sptr
        integer, dimension(:), allocatable            :: slist
        integer, dimension(:), allocatable            :: req     ! MPI_Wait requests (runtime; not serialized)
        integer :: nreq = 0
    contains
        procedure :: write_com => write_com_struct
        procedure :: read_com  => read_com_struct
    end type com_struct

    type t_partit
        integer :: MPI_COMM_FESOM           ! FESOM communicator
        integer :: npes = 1                 ! number of PEs
        integer :: mype = 0                 ! this PE rank
        integer :: MPIERR

        type(com_struct) :: com_nod2D
        type(com_struct) :: com_elem2D
        type(com_struct) :: com_elem2D_full
        ! Hybrid meshes, read only when dist_N/format.out exists (dist_format >= 1):
        ! com_nod2D_full receives every halo node 1..eDim+eXDim_nod2D from its owner,
        ! com_edge2D every held edge from its unique owner (owner of edges(1,edge)).
        type(com_struct) :: com_nod2D_full
        type(com_struct) :: com_edge2D
        integer :: dist_format = 0          ! dist_N/format.out version (0: FESOM2 layout only)
        ! Shared elements (src/infra/mod_shared_elem.F90; runtime, not serialized): an element
        ! is owned by every rank that owns one of its vertices, so a border element has 2-3
        ! owned copies. com_elem2D_shared sends the copy of the designated owner (the highest
        ! vertex-owner rank) to the other owners: receive blocks = this rank's non-designated
        ! copies, send blocks = the elements it is the designated owner of. Built at the end
        ! of read_mesh_local; sync_shared_elem = FESOM3_SHARED_SYNC (update_vel syncs uv).
        type(com_struct) :: com_elem2D_shared
        logical :: sync_shared_elem = .true.

        integer, allocatable, dimension(:) :: part   ! global owner map (nod2D)

        ! Local dimensions and global-id lists (1-based local -> global).
        ! myList_nod2D holds myDim+eDim+eXDim_nod2D ids; eXDim_nod2D (double halo: the
        ! vertices of the local elements 1..myDim+eDim+eXDim_elem2D that are no owned or
        ! halo node) is 0 without dist_N/format.out.
        integer :: myDim_nod2D = 0, eDim_nod2D = 0, eXDim_nod2D = 0
        integer, allocatable, dimension(:) :: myList_nod2D
        integer :: myDim_elem2D = 0, eDim_elem2D = 0, eXDim_elem2D = 0
        integer, allocatable, dimension(:) :: myList_elem2D
        integer :: myDim_edge2D = 0, eDim_edge2D = 0
        integer, allocatable, dimension(:) :: myList_edge2D
        integer :: pe_status = 0

        ! Precompiled MPI indexed datatypes for interface exchange (runtime; not
        ! serialized). Built by init_mpi_types (M0.4).
        integer, allocatable :: s_mpitype_nod2D(:),      r_mpitype_nod2D(:)
        integer, allocatable :: s_mpitype_nod2D_i(:),    r_mpitype_nod2D_i(:)
        integer, allocatable :: s_mpitype_nod3D(:,:,:),  r_mpitype_nod3D(:,:,:)
        integer, allocatable :: s_mpitype_elem2D(:,:),       r_mpitype_elem2D(:,:)
        integer, allocatable :: s_mpitype_elem2D_full_i(:),  r_mpitype_elem2D_full_i(:)
        integer, allocatable :: s_mpitype_elem2D_full(:,:),  r_mpitype_elem2D_full(:,:)
        integer, allocatable :: s_mpitype_elem3D(:,:,:),     r_mpitype_elem3D(:,:,:)
        integer, allocatable :: s_mpitype_elem3D_full(:,:,:),r_mpitype_elem3D_full(:,:,:)
    contains
        procedure :: write_unformatted => write_t_partit
        procedure :: read_unformatted  => read_t_partit
        generic   :: write(unformatted) => write_unformatted
        generic   :: read(unformatted)  => read_unformatted
    end type t_partit

contains

    subroutine write_com_struct(tstruct, unit)
        class(com_struct), intent(in) :: tstruct
        integer,           intent(in) :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        write(unit, iostat=iostat, iomsg=iomsg) tstruct%rPEnum
        call write1d_int_static(tstruct%rPE,  unit, iostat, iomsg)
        call write1d_int_static(tstruct%rptr, unit, iostat, iomsg)
        call write_bin_array(tstruct%rlist,   unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) tstruct%sPEnum
        call write1d_int_static(tstruct%sPE,  unit, iostat, iomsg)
        call write1d_int_static(tstruct%sptr, unit, iostat, iomsg)
        call write_bin_array(tstruct%slist,   unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) tstruct%nreq
    end subroutine write_com_struct

    subroutine read_com_struct(tstruct, unit)
        class(com_struct), intent(inout) :: tstruct
        integer,           intent(in)    :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        read(unit, iostat=iostat, iomsg=iomsg) tstruct%rPEnum
        call read1d_int_static(tstruct%rPE,  unit, iostat, iomsg)
        call read1d_int_static(tstruct%rptr, unit, iostat, iomsg)
        call read_bin_array(tstruct%rlist,   unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) tstruct%sPEnum
        call read1d_int_static(tstruct%sPE,  unit, iostat, iomsg)
        call read1d_int_static(tstruct%sptr, unit, iostat, iomsg)
        call read_bin_array(tstruct%slist,   unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) tstruct%nreq
    end subroutine read_com_struct

    subroutine write_t_partit(partit, unit, iostat, iomsg)
        class(t_partit), intent(in)    :: partit
        integer,         intent(in)    :: unit
        integer,         intent(out)   :: iostat
        character(*),    intent(inout) :: iomsg
        call partit%com_nod2D%write_com(unit)
        call partit%com_elem2D%write_com(unit)
        call partit%com_elem2D_full%write_com(unit)
        call partit%com_nod2D_full%write_com(unit)
        call partit%com_edge2D%write_com(unit)
        write(unit, iostat=iostat, iomsg=iomsg) partit%dist_format
        write(unit, iostat=iostat, iomsg=iomsg) partit%npes
        write(unit, iostat=iostat, iomsg=iomsg) partit%mype
        call write_bin_array(partit%part, unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) partit%myDim_nod2D
        write(unit, iostat=iostat, iomsg=iomsg) partit%eDim_nod2D
        write(unit, iostat=iostat, iomsg=iomsg) partit%eXDim_nod2D
        call write_bin_array(partit%myList_nod2D, unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) partit%myDim_elem2D
        write(unit, iostat=iostat, iomsg=iomsg) partit%eDim_elem2D
        write(unit, iostat=iostat, iomsg=iomsg) partit%eXDim_elem2D
        call write_bin_array(partit%myList_elem2D, unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) partit%myDim_edge2D
        write(unit, iostat=iostat, iomsg=iomsg) partit%eDim_edge2D
        call write_bin_array(partit%myList_edge2D, unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) partit%pe_status
    end subroutine write_t_partit

    subroutine read_t_partit(partit, unit, iostat, iomsg)
        class(t_partit), intent(inout) :: partit
        integer,         intent(in)    :: unit
        integer,         intent(out)   :: iostat
        character(*),    intent(inout) :: iomsg
        call partit%com_nod2D%read_com(unit)
        call partit%com_elem2D%read_com(unit)
        call partit%com_elem2D_full%read_com(unit)
        call partit%com_nod2D_full%read_com(unit)
        call partit%com_edge2D%read_com(unit)
        read(unit, iostat=iostat, iomsg=iomsg) partit%dist_format
        read(unit, iostat=iostat, iomsg=iomsg) partit%npes
        read(unit, iostat=iostat, iomsg=iomsg) partit%mype
        call read_bin_array(partit%part, unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) partit%myDim_nod2D
        read(unit, iostat=iostat, iomsg=iomsg) partit%eDim_nod2D
        read(unit, iostat=iostat, iomsg=iomsg) partit%eXDim_nod2D
        call read_bin_array(partit%myList_nod2D, unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) partit%myDim_elem2D
        read(unit, iostat=iostat, iomsg=iomsg) partit%eDim_elem2D
        read(unit, iostat=iostat, iomsg=iomsg) partit%eXDim_elem2D
        call read_bin_array(partit%myList_elem2D, unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) partit%myDim_edge2D
        read(unit, iostat=iostat, iomsg=iomsg) partit%eDim_edge2D
        call read_bin_array(partit%myList_edge2D, unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) partit%pe_status
    end subroutine read_t_partit

end module mod_partit
