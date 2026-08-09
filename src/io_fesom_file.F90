 ! synopsis: generic implementation to asynchronously read/write FESOM mesh variable(s) with distributed cell or element data in 2D or 3D to/from a NetCDF file
module io_fesom_file_module
  use io_netcdf_file_module
  use async_threads_module
  use io_redistribute
  use MOD_PARTIT
  implicit none
  public fesom_file_type
  public set_parallel_write, parallel_write_enabled
  private
  
  
  type var_info
    integer var_index
    real(kind=8), pointer :: external_local_data_ptr(:,:) => null()
    real(kind=8), allocatable, dimension(:,:) :: local_data_copy
    real(kind=8), allocatable :: global_level_data(:)
    integer :: global_level_data_size = 0
    logical is_elem_based
    logical :: is_icepackvar2=.false.
    character(:), allocatable :: varname ! todo: maybe use a getter in netcdf_file_type to get the name
  end type var_info
  
  
  type dim_info
    integer idx
    integer len ! better query the len from the netcdf_file_type ?
  end type dim_info

  
  type, extends(netcdf_file_type) :: fesom_file_type ! todo maybe: do not inherit but use composition to have different implementations for the iorank and non-io ranks
    private
    integer time_dimidx
    integer time_varidx
    type(var_info) var_infos(20); integer :: nvar_infos = 0 ! todo: allow dynamically allocated size without messing with shallow copied pointers
    type(dim_info), allocatable :: used_mesh_dims(:) ! the dims we add for our variables, we need to identify them when adding our mesh related variables
    integer :: rec_cnt = -1
    integer :: iorank = 0
    integer :: fesom_file_index
    type(thread_type) thread
    logical :: thread_running = .false.
    integer :: comm
    type(t_partit), pointer :: partit
    logical gather_and_write
  contains
    procedure, public :: async_read_and_scatter_variables, async_gather_and_write_variables, join, init, is_iorank, rec_count, time_varindex, time_dimindex
    procedure, public :: read_variables_raw, write_variables_raw
    procedure, public :: close_file ! inherited procedures we overwrite
    procedure, public :: read_and_scatter_variables
    generic, public :: specify_node_var => specify_node_var_2d, specify_node_var_2dicepack, specify_node_var_3d 
    generic, public :: specify_elem_var => specify_elem_var_2d, specify_elem_var_3d
    procedure, private :: specify_node_var_2d, specify_node_var_2dicepack, specify_node_var_3d
    procedure, private :: specify_elem_var_2d, specify_elem_var_3d
    procedure, private :: gather_and_write_variables
    procedure, private :: redist_and_write_variables
    procedure, public  :: is_writer, is_lead_writer, writer_comm
  end type fesom_file_type
  
  
  integer, save :: m_nod2d
  integer, save :: m_elem2d
  integer, save :: m_nl
  integer, save :: m_ncat

  ! --- parallel (collective) write path -------------------------------------
  ! Off by default, so the gather path stays the shipped behaviour and remains
  ! available for A/B comparison. Enable with set_parallel_write before any
  ! fesom_file_type is initialised.
  !
  ! The schedules depend only on the partition and the writer count, not on the
  ! file, so all files share one nodal and one element schedule, built once.
  logical, save :: m_parallel_write = .false.
  integer, save :: m_n_writers = 0
  type(redist_type), save, target :: m_sched_nod, m_sched_elem
  logical, save :: m_sched_ready = .false.
  integer, save :: m_writer_comm = -1   !< communicator holding exactly the writers
  

  type fesom_file_type_ptr
    class(fesom_file_type), pointer :: ptr
  end type fesom_file_type_ptr
  type(fesom_file_type_ptr), allocatable, save :: all_fesom_files(:)


contains


  function is_iorank(this) result(x)
    class(fesom_file_type), intent(in) :: this
    logical x
    x = (this%partit%mype == this%iorank)
  end function is_iorank


  ! return the number of timesteps of the file if a file is attached or return the default value of -1
  function rec_count(this) result(x)
    class(fesom_file_type), intent(inout) :: this
    integer x
    ! EO parameters
    integer, allocatable :: time_shape(:)
    
    if(this%rec_cnt == -1 .and. this%is_attached()) then
      ! update from file if rec_cnt has never been used before
      call this%read_var_shape(this%time_varidx, time_shape)
      this%rec_cnt = time_shape(1)
    end if
    
    x = this%rec_cnt
  end function rec_count


  function time_varindex(this) result(x)
    class(fesom_file_type), intent(in) :: this
    integer x
    x = this%time_varidx
  end function time_varindex


  function time_dimindex(this) result(x)
    class(fesom_file_type), intent(in) :: this
    integer x
    x = this%time_dimidx
  end function time_dimindex
  
  
  subroutine init(this, mesh_nod2d, mesh_elem2d, mesh_nl, partit, mesh_ncat) ! todo: would like to call it initialize but Fortran is rather cluncky with overwriting base type procedures
    use io_netcdf_workaround_module
    use io_gather_module
    use MOD_PARTIT
    class(fesom_file_type), target, intent(inout) :: this
    integer mesh_nod2d
    integer mesh_elem2d
    integer mesh_nl
    integer, optional :: mesh_ncat
    type(t_partit), target :: partit
    ! EO parameters
    type(fesom_file_type_ptr), allocatable :: tmparr(:)
    logical async_netcdf_allowed
    integer err
    integer provided_mpi_thread_support_level

    call init_io_gather(partit)

    ! get hold of our mesh data for later use (assume the mesh instance will not change)
    m_nod2d  = mesh_nod2d
    m_elem2d = mesh_elem2d
    m_nl     = mesh_nl
    !PS mesh_ncat ... icepack number of ice thickness classes,
    if (present(mesh_ncat)) then 
        m_ncat   = mesh_ncat
    else    
        m_ncat   = 0
    end if 
    
    call this%netcdf_file_type%initialize()

    allocate(this%used_mesh_dims(0))

    this%time_dimidx = this%add_dim_unlimited('time')

    this%time_varidx = this%add_var_double('time', [this%time_dimidx])

    ! add this instance to global array
    ! the array is being used to identify the instance in an async call
    if( .not. allocated(all_fesom_files)) then
      allocate(all_fesom_files(1))
    else
      allocate( tmparr(size(all_fesom_files)+1) )
      tmparr(1:size(all_fesom_files)) = all_fesom_files
      deallocate(all_fesom_files)
      call move_alloc(tmparr, all_fesom_files)
    end if
    all_fesom_files(size(all_fesom_files))%ptr => this
    this%fesom_file_index = size(all_fesom_files)

    this%partit => partit
    ! set up async output
    
    this%iorank = next_io_rank(partit%MPI_COMM_FESOM, async_netcdf_allowed, partit)

    call MPI_Comm_dup(partit%MPI_COMM_FESOM, this%comm, err)

    call this%thread%initialize(async_worker, this%fesom_file_index)
    if(.not. async_netcdf_allowed) call this%thread%disable_async()
  
    ! check if we have multi thread support available in the MPI library
    ! tough MPI_THREAD_FUNNELED should be enough here, at least on cray-mpich 7.5.3 async mpi calls fail if we do not have support level 'MPI_THREAD_MULTIPLE'
    ! on cray-mpich we only get level 'MPI_THREAD_MULTIPLE' if 'MPICH_MAX_THREAD_SAFETY=multiple' is set in the environment
    call MPI_Query_thread(provided_mpi_thread_support_level, err)
    if(provided_mpi_thread_support_level < MPI_THREAD_MULTIPLE) call this%thread%disable_async()

    ! The collective path MUST NOT run on a worker thread. Two independent
    ! reasons: (1) Levante's HDF5 1.12.1 reports Threadsafety: OFF, and with a
    ! shared writer set several files' workers would land on the same rank and
    ! call HDF5 concurrently; (2) a collective put_var issued from a worker while
    ! the main thread proceeds deadlocks at the next matched collective.
    if(m_parallel_write) then
      call this%thread%disable_async()
      ! Schedules and the writer communicator must exist before the first file
      ! is created, because the create itself is collective over the writers.
      call ensure_schedules(partit, this%comm)
    end if
  end subroutine init
  
  
  subroutine read_and_scatter_variables(this)
#ifdef ENABLE_NVHPC_WORKAROUNDS
    use nvfortran_subarray_workaround_module
#endif
    use io_scatter_module
    class(fesom_file_type), target :: this
    ! EO parameters
    integer i,lvl, nlvl
    logical is_2d
    integer last_rec_idx
    type(var_info), pointer :: var
    real(kind=8), allocatable :: laux(:)
    integer mpierr
    real(kind=8) :: total_netcdf_time, total_scatter_time, total_barrier_time
    real(kind=8) :: t_start, t_end

    ! Initialize timing variables
    total_netcdf_time = 0.0d0
    total_scatter_time = 0.0d0
    total_barrier_time = 0.0d0

    last_rec_idx = this%rec_count()
    
    do i=1, this%nvar_infos
      var => this%var_infos(i)
    
    
      if (var%is_icepackvar2) then 
        !PS for icepack external_local_data_ptr has still dimension [nod2, ncat]
        !PS but usual fesom output would have [nlev, nod2] so we need to change here
        !PS dim= parameter to get the proper dimension 
        nlvl = size(var%external_local_data_ptr,dim=2)
        allocate(laux( size(var%external_local_data_ptr,dim=1) )) ! i.e. myDim_elem2D+eDim_elem2D or myDim_nod2D+eDim_nod2D
      else
        nlvl = size(var%external_local_data_ptr,dim=1)
        allocate(laux( size(var%external_local_data_ptr,dim=2) )) ! i.e. myDim_elem2D+eDim_elem2D or myDim_nod2D+eDim_nod2D
      end if 
      is_2d = (nlvl == 1)

      if(this%is_iorank()) then
        ! todo: choose how many levels we read at once
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( var%global_level_data_size ))
      else
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( 0 ))
      end if

      do lvl=1, nlvl
        if(this%is_iorank()) then
          t_start = MPI_Wtime()
          if(is_2d) then
            call this%read_var(var%var_index, [1,last_rec_idx], [size(var%global_level_data),1], var%global_level_data)
          else
            call this%read_var(var%var_index, [1,lvl,last_rec_idx], [size(var%global_level_data),1,1], var%global_level_data)
          end if
          t_end = MPI_Wtime()
          total_netcdf_time = total_netcdf_time + (t_end - t_start)
        end if

        t_start = MPI_Wtime()
        if(var%is_elem_based) then
          call scatter_elem2D(var%global_level_data, laux, this%iorank, this%comm, this%partit)
        else
          call scatter_nod2D(var%global_level_data, laux, this%iorank, this%comm, this%partit)
        end if
        t_end = MPI_Wtime()
        total_scatter_time = total_scatter_time + (t_end - t_start)
#ifdef ENABLE_NVHPC_WORKAROUNDS
  if(var%varname=='u') then
    dynamics_workaround%uv(1,lvl,:) = laux
  else if(var%varname=='v') then
    dynamics_workaround%uv(2,lvl,:) = laux
  else if(var%varname=='urhs_AB') then
    dynamics_workaround%uv_rhsAB(1,lvl,:) = laux
  else if(var%varname=='vrhs_AB') then
    dynamics_workaround%uv_rhsAB(2,lvl,:) = laux
  else
#endif
       if (var%is_icepackvar2) then 
            ! the data from our pointer is not contiguous (if it is 3D data), so we can not pass the pointer directly to MPI
            ! PS Again icepack variable demension is [nod2, ncat] thats why we need 
            ! PS to switch here the index where we sort in the 2d slices
            var%external_local_data_ptr(:,lvl) = laux ! todo: remove this buffer and pass the data directly to MPI (change order of data layout to be levelwise or do not gather levelwise but by columns)
       else
            ! the data from our pointer is not contiguous (if it is 3D data), so we can not pass the pointer directly to MPI
            var%external_local_data_ptr(lvl,:) = laux ! todo: remove this buffer and pass the data directly to MPI (change order of data layout to be levelwise or do not gather levelwise but by columns)
       end if 
#ifdef ENABLE_NVHPC_WORKAROUNDS
  end if
#endif
      end do
      deallocate(laux)
    end do
    
    ! Print timing breakdown for restart reading (rank 0 only)
    if(this%partit%mype == 0 .and. this%nvar_infos > 0) then
      write(*,'(A)') '=================================='
      if(allocated(this%var_infos(1)%varname)) then
        write(*,'(A,A)') 'RESTART READ TIMING for file: ', trim(this%var_infos(1)%varname)
      else
        write(*,'(A)') 'RESTART READ TIMING BREAKDOWN'
      end if
      write(*,'(A,F10.3,A)') '  NetCDF read time:  ', total_netcdf_time, ' seconds'
      write(*,'(A,F10.3,A)') '  MPI scatter time:  ', total_scatter_time, ' seconds'
      if(total_barrier_time > 0.0d0) then
        write(*,'(A,F10.3,A)') '  MPI barrier time:  ', total_barrier_time, ' seconds'
      end if
      write(*,'(A,F10.3,A)') '  TOTAL time:        ', total_netcdf_time + total_scatter_time + total_barrier_time, ' seconds'
      write(*,'(A)') '=================================='
    end if
  end subroutine


  !> Select the collective write path. Call before any file is initialised.
  !> n_writers <= 0 means "as many as the block-size guard allows".
  subroutine set_parallel_write(enable, n_writers)
    logical, intent(in) :: enable
    integer, intent(in) :: n_writers
    m_parallel_write = enable
    m_n_writers = n_writers
  end subroutine set_parallel_write


  logical function parallel_write_enabled()
    parallel_write_enabled = m_parallel_write
  end function parallel_write_enabled


  !> Build the nodal and element schedules once, from the partition.
  !> Elements use myInd_elem2D_shrinked, NOT myList_elem2D: an element belongs
  !> to a PE if ANY of its nodes does (gen_comm.F90:265), so sum(myDim_elem2D)
  !> exceeds elem2D and there is no bijection. The shrinked set assigns each
  !> element to the rank owning its first node (oce_mesh.F90:888-898).
  subroutine ensure_schedules(partit, comm)
    type(t_partit), intent(in) :: partit
    integer, intent(in) :: comm
    integer :: ierr, k
    integer, allocatable :: gelem(:)

    if(m_sched_ready) return

    call redist_build(m_sched_nod, partit%myList_nod2D(1:partit%myDim_nod2D), &
                      partit%myDim_nod2D, m_nod2d, m_n_writers, comm, &
                      partit%mype==0, ierr)
    if(ierr /= 0) then
      if(partit%mype==0) write(*,*) 'io_fesom_file: nodal redist_build failed, ierr=', ierr
      stop 1
    end if

    allocate(gelem(max(1, partit%myDim_elem2D_shrinked)))
    do k = 1, partit%myDim_elem2D_shrinked
      gelem(k) = partit%myList_elem2D( partit%myInd_elem2D_shrinked(k) )
    end do
    call redist_build(m_sched_elem, gelem, partit%myDim_elem2D_shrinked, &
                      m_elem2d, m_n_writers, comm, .false., ierr)
    if(ierr /= 0) then
      if(partit%mype==0) write(*,*) 'io_fesom_file: element redist_build failed, ierr=', ierr
      stop 1
    end if
    deallocate(gelem)

    ! Writer sub-communicator. The nodal and element schedules select the same
    ! ranks -- both use redist_writer_rank with the same npes and effective
    ! writer count -- so one communicator serves both. Assert that rather than
    ! rely on it silently.
    if(m_sched_nod%is_writer .neqv. m_sched_elem%is_writer) then
      write(*,*) 'io_fesom_file: nodal and element writer sets disagree on rank ', partit%mype
      stop 1
    end if
    call MPI_Comm_split(comm, merge(1,0,m_sched_nod%is_writer), partit%mype, m_writer_comm, ierr)

    m_sched_ready = .true.
  end subroutine ensure_schedules


  !> Communicator containing exactly the writer ranks. Valid after init.
  integer function writer_comm(this)
    class(fesom_file_type), intent(in) :: this
    writer_comm = m_writer_comm
  end function writer_comm


  !> One designated writer, for per-record scalars like iter and time that must
  !> be written once rather than once per writer.
  logical function is_lead_writer(this)
    class(fesom_file_type), intent(in) :: this
    is_lead_writer = m_parallel_write .and. m_sched_ready .and. (m_sched_nod%widx == 0)
  end function is_lead_writer


  !> Collective write: one MPI_Alltoallv moves a whole field into block order,
  !> then every writer issues one put_var per level with identical call counts.
  !> Writers whose block is empty still call, with count 0, so the collective
  !> stays matched.
  subroutine redist_and_write_variables(this)
    class(fesom_file_type), target :: this
    integer :: i, lvl, nlvl, k, ierr, nloc
    type(var_info), pointer :: var
    type(redist_type), pointer :: sched
    real(kind=8), allocatable :: sbuf(:,:), blk(:,:), lvlbuf(:), flat(:,:)
    real(kind=8) :: empty(0)

    call ensure_schedules(this%partit, this%comm)

    ! Every writer has the same file open, so they all read the same record
    ! count and agree without a broadcast.
    if(this%is_writer()) this%rec_cnt = this%rec_count()+1

    do i=1, this%nvar_infos
      var => this%var_infos(i)
      nlvl = size(var%local_data_copy, dim=1)

      if(var%is_elem_based) then
        sched => m_sched_elem
        nloc = this%partit%myDim_elem2D_shrinked
        allocate(sbuf(nlvl, max(1,nloc)))
        do k = 1, nloc
          sbuf(:,k) = var%local_data_copy(:, this%partit%myInd_elem2D_shrinked(k))
        end do
      else
        sched => m_sched_nod
        nloc = this%partit%myDim_nod2D
        allocate(sbuf(nlvl, max(1,nloc)))
        do k = 1, nloc
          sbuf(:,k) = var%local_data_copy(:, k)
        end do
      end if

      allocate(blk(nlvl, max(1, sched%count)))
      call redist_exchange(sched, nlvl, sbuf, blk, this%comm, ierr)
      if(ierr /= 0) then
        write(*,*) 'io_fesom_file: redist_exchange failed for ', trim(var%varname)
        stop 1
      end if

      if(this%is_writer()) then
        ! ONE write per variable, not one per level. On the output path the same
        ! change was worth 1.35x on NG5/1024 purely by collapsing 69 collective
        ! calls into one; with a chunk spanning several levels it also stops
        ! HDF5 rewriting a partially-filled chunk once per level.
        !
        ! write_var takes `values` as TARGET and calls c_loc on it, so the array
        ! handed over must be contiguous -- passing a strided section blk(lvl,:)
        ! silently wrote permuted data for every 3D field while 2D ones, where
        ! blk(1,:) happens to be contiguous, looked correct. Hence the explicit
        ! transpose into (count, nlvl), which is also the file's own order.
        if(nlvl == 1) then
          if(sched%count > 0) then
            allocate(lvlbuf(sched%count))
            lvlbuf(1:sched%count) = blk(1, 1:sched%count)
            call this%write_var(var%var_index, [sched%first, this%rec_cnt], &
                                [sched%count, 1], lvlbuf)
            deallocate(lvlbuf)
          else
            call this%write_var(var%var_index, [1, this%rec_cnt], [0, 1], empty)
          end if
        else
          if(sched%count > 0) then
            allocate(flat(sched%count, nlvl))
            do lvl=1, nlvl
              flat(1:sched%count, lvl) = blk(lvl, 1:sched%count)
            end do
            call this%write_var(var%var_index, [sched%first, 1, this%rec_cnt], &
                                [sched%count, nlvl, 1], reshape(flat, [sched%count*nlvl]))
            deallocate(flat)
          else
            call this%write_var(var%var_index, [1, 1, this%rec_cnt], [0, nlvl, 1], empty)
          end if
        end if
      end if

      deallocate(sbuf, blk)
    end do
  end subroutine redist_and_write_variables


  logical function is_writer(this)
    class(fesom_file_type), intent(in) :: this
    is_writer = m_parallel_write .and. m_sched_ready .and. m_sched_nod%is_writer
  end function is_writer


  subroutine gather_and_write_variables(this)
    use io_gather_module
    class(fesom_file_type), target :: this
    ! EO parameters
    integer i,lvl, nlvl
    logical is_2d
    real(kind=8), allocatable :: laux(:)
    type(var_info), pointer :: var
    integer mpierr

    ! Collective path: every rank participates, the record index is decided by
    ! the caller (io_restart) and broadcast, and there is no gather at all.
    if(m_parallel_write) then
      call this%redist_and_write_variables()
      return
    end if

    if(this%is_iorank()) this%rec_cnt = this%rec_count()+1

    do i=1, this%nvar_infos
      var => this%var_infos(i)

      nlvl = size(var%local_data_copy,dim=1)
      is_2d = (nlvl == 1)
      allocate(laux( size(var%local_data_copy,dim=2) )) ! i.e. myDim_elem2D+eDim_elem2D or myDim_nod2D+eDim_nod2D

      if(this%is_iorank()) then
        ! todo: choose how many levels we write at once
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( var%global_level_data_size ))
      else
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( 0 ))
      end if

      do lvl=1, nlvl
        ! the data from our pointer is not contiguous (if it is 3D data), so we can not pass the pointer directly to MPI
        laux = var%local_data_copy(lvl,:) ! todo: remove this buffer and pass the data directly to MPI (change order of data layout to be levelwise or do not gather levelwise but by columns)

        if(var%is_elem_based) then
          call gather_elem2D(laux, var%global_level_data, this%iorank, 42, this%comm, this%partit)
        else
          call gather_nod2D (laux, var%global_level_data, this%iorank, 42, this%comm, this%partit)
        end if

        if(this%is_iorank()) then
          if(is_2d) then
            call this%write_var(var%var_index, [1,this%rec_cnt], [size(var%global_level_data),1], var%global_level_data)
          else
           call this%write_var(var%var_index, [1,lvl,this%rec_cnt], [size(var%global_level_data),1,1], var%global_level_data)
          end if
        end if
      end do
      deallocate(laux)
    end do
    
    if(this%is_iorank()) call this%flush_file() ! flush the file to disk after each write
  end subroutine


  subroutine read_variables_raw(this, fileunit)
#ifdef ENABLE_NVHPC_WORKAROUNDS
    use nvfortran_subarray_workaround_module
#endif
    class(fesom_file_type), target :: this
    integer, intent(in) :: fileunit
    ! EO parameters
    integer i
    type(var_info), pointer :: var
    integer status
    
    do i=1, this%nvar_infos
      var => this%var_infos(i)
#ifdef ENABLE_NVHPC_WORKAROUNDS
      if(var%varname=='u') then
        read(fileunit) dynamics_workaround%uv(1,:,:)
      else if(var%varname=='v') then
        read(fileunit) dynamics_workaround%uv(2,:,:)
      else if(var%varname=='urhs_AB') then
        read(fileunit) dynamics_workaround%uv_rhsAB(1,:,:)
      else if(var%varname=='vrhs_AB') then
        read(fileunit) dynamics_workaround%uv_rhsAB(2,:,:)
      else
#endif
      read(fileunit) var%external_local_data_ptr ! directly use external_local_data_ptr, use the local_data_copy only when called asynchronously
#ifdef ENABLE_NVHPC_WORKAROUNDS
      end if
#endif
    end do
  end subroutine


  subroutine write_variables_raw(this, fileunit)
#ifdef ENABLE_NVHPC_WORKAROUNDS
    use nvfortran_subarray_workaround_module
#endif
    class(fesom_file_type), target :: this
    integer, intent(in) :: fileunit
    ! EO parameters
    integer i
    type(var_info), pointer :: var
    
    do i=1, this%nvar_infos
      var => this%var_infos(i)
#ifdef ENABLE_NVHPC_WORKAROUNDS
      if(var%varname=='u') then
        write(fileunit) dynamics_workaround%uv(1,:,:)
      else if(var%varname=='v') then
        write(fileunit) dynamics_workaround%uv(2,:,:)
      else if(var%varname=='urhs_AB') then
        write(fileunit) dynamics_workaround%uv_rhsAB(1,:,:)
      else if(var%varname=='vrhs_AB') then
        write(fileunit) dynamics_workaround%uv_rhsAB(2,:,:)
      else
#endif
      write(fileunit) var%external_local_data_ptr ! directly use external_local_data_ptr, use the local_data_copy only when called asynchronously
#ifdef ENABLE_NVHPC_WORKAROUNDS
      end if
#endif
    end do
  end subroutine


  subroutine join(this)
    class(fesom_file_type) this
    ! EO parameters
    
    if(this%thread_running) call this%thread%join()
    this%thread_running = .false.    
  end subroutine


  subroutine async_read_and_scatter_variables(this)
    class(fesom_file_type), target :: this

    call assert(.not. this%thread_running, __LINE__)

    this%gather_and_write = .false.
    call this%thread%run()
    this%thread_running = .true.
  end subroutine


  subroutine async_gather_and_write_variables(this)
#ifdef ENABLE_NVHPC_WORKAROUNDS
use nvfortran_subarray_workaround_module
#endif
    class(fesom_file_type), target :: this
    ! EO parameters
    integer i
    type(var_info), pointer :: var
    
    call assert(.not. this%thread_running, __LINE__)

    ! copy data so we can write the current values asynchronously
    do i=1, this%nvar_infos
      var => this%var_infos(i)
      if(.not. allocated(var%local_data_copy)) allocate( var%local_data_copy(size(var%external_local_data_ptr,dim=1), size(var%external_local_data_ptr,dim=2)) )
#ifdef ENABLE_NVHPC_WORKAROUNDS
      if(var%varname=='u') then
        var%local_data_copy = dynamics_workaround%uv(1,:,:)
      else if(var%varname=='v') then
        var%local_data_copy = dynamics_workaround%uv(2,:,:)
      else if(var%varname=='urhs_AB') then
        var%local_data_copy = dynamics_workaround%uv_rhsAB(1,:,:)
      else if(var%varname=='vrhs_AB') then
        var%local_data_copy = dynamics_workaround%uv_rhsAB(2,:,:)
      else
#endif
        if (var%is_icepackvar2) then 
            var%local_data_copy = transpose(var%external_local_data_ptr)
        else
            var%local_data_copy = var%external_local_data_ptr
        end if   
#ifdef ENABLE_NVHPC_WORKAROUNDS
      end if
#endif
    end do
    
    this%gather_and_write = .true.
    call this%thread%run()
    this%thread_running = .true.
  end subroutine


  subroutine async_worker(fesom_file_index)
    integer, intent(in) :: fesom_file_index
    ! EO parameters
    type(fesom_file_type), pointer :: f

    f => all_fesom_files(fesom_file_index)%ptr

    if(f%gather_and_write) then
      call f%gather_and_write_variables()
    else
      call f%read_and_scatter_variables()
    end if
  end subroutine


  ! use separate procedures to specify node based or element based variables
  ! if we would otherwise specify the vars only via the sizes of their dimensions,
  ! we have to assign the corresponding dimindx somewhere else, which would be error prone
  subroutine specify_node_var_2d(this, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:) ! todo: be able to set precision
    ! EO parameters
    real(8), pointer :: external_local_data_ptr(:,:)
    type(dim_info) level_diminfo

!PS     write(*,*) "--> specify_node_var_2d:", __LINE__, __FILE__
    level_diminfo = obtain_diminfo(this, m_nod2d)
   
    external_local_data_ptr(1:1,1:size(local_data)) => local_data(:)
    call specify_variable(this, name, [level_diminfo%idx, this%time_dimidx], level_diminfo%len, external_local_data_ptr, .false., longname, units)    
  end subroutine
  
  subroutine specify_node_var_2dicepack(this, name, longname, units, local_data, ncat)
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:)
    integer, intent(in) :: ncat! todo: be able to set precision
    ! EO parameters
    type(dim_info) level_diminfo, ncat_diminfo
    
!PS     write(*,*) "--> specify_node_var_2dicepack:", __LINE__, __FILE__, size(local_data)
    level_diminfo = obtain_diminfo(this, m_nod2d)    
    ncat_diminfo = obtain_diminfo(this, size(local_data, dim=2))

    call specify_variable(this, name, [level_diminfo%idx, ncat_diminfo%idx, this%time_dimidx], level_diminfo%len, local_data, .false., longname, units, ncat)
    
  end subroutine


  subroutine specify_node_var_3d(this, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision
    ! EO parameters
    type(dim_info) level_diminfo, depth_diminfo
    
!PS     write(*,*) "--> specify_node_var_3d:", __LINE__, __FILE__
    level_diminfo = obtain_diminfo(this, m_nod2d)    
    depth_diminfo = obtain_diminfo(this, size(local_data, dim=1))

    call specify_variable(this, name, [level_diminfo%idx, depth_diminfo%idx, this%time_dimidx], level_diminfo%len, local_data, .false., longname, units)
  end subroutine


  subroutine specify_elem_var_2d(this, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:) ! todo: be able to set precision
    ! EO parameters
    real(8), pointer :: external_local_data_ptr(:,:)
    type(dim_info) level_diminfo

!PS     write(*,*) "--> specify_elem_var_2d:", __LINE__, __FILE__
    level_diminfo = obtain_diminfo(this, m_elem2d)

    external_local_data_ptr(1:1,1:size(local_data)) => local_data(:)
    call specify_variable(this, name, [level_diminfo%idx, this%time_dimidx], level_diminfo%len, external_local_data_ptr, .true., longname, units)    
  end subroutine


  subroutine specify_elem_var_3d(this, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision
    ! EO parameters
    type(dim_info) level_diminfo, depth_diminfo

!PS     write(*,*) "--> specify_elem_var_3d:", __LINE__, __FILE__
    level_diminfo = obtain_diminfo(this, m_elem2d)
    depth_diminfo = obtain_diminfo(this, size(local_data, dim=1))
    
   call specify_variable(this, name, [level_diminfo%idx, depth_diminfo%idx, this%time_dimidx], level_diminfo%len, local_data, .true., longname, units)
  end subroutine
  
  
  function obtain_diminfo(this, len) result(info)
    type(fesom_file_type), intent(inout) :: this
    type(dim_info) info
    integer len
    ! EO parameters
    integer i
    type(dim_info), allocatable :: tmparr(:)
    
    do i=1, size(this%used_mesh_dims)
      if(this%used_mesh_dims(i)%len == len) then
        info = this%used_mesh_dims(i)
        return
      end if
    end do
    
    ! the dim has not been added yet, see if it is one of our allowed mesh related dims
    if     (len == m_nod2d)  then
      info = dim_info( idx=this%add_dim('node', len), len=len)
    else if(len == m_elem2d) then
      info = dim_info( idx=this%add_dim('elem', len), len=len)
    else if(len == m_nl-1  ) then
      info = dim_info( idx=this%add_dim('nz_1', len), len=len)
    else if(len == m_nl    ) then
      info = dim_info( idx=this%add_dim('nz'  , len), len=len)
    else if(len == m_ncat  ) then
      info = dim_info( idx=this%add_dim('ncat', len), len=len)  
    else
      print *, "error in line ",__LINE__, __FILE__," can not find dimension with size",len
      stop 1
    end if
    
    ! append the new dim to our list of used dims, i.e. the dims we use for the mesh based variables created via #specify_variable
    ! assume the used_mesh_dims array is allocated
    allocate( tmparr(size(this%used_mesh_dims)+1) )
    tmparr(1:size(this%used_mesh_dims)) = this%used_mesh_dims
    deallocate(this%used_mesh_dims)
    call move_alloc(tmparr, this%used_mesh_dims)
    this%used_mesh_dims( size(this%used_mesh_dims) ) = info
  end function


  subroutine specify_variable(this, name, dim_indices, global_level_data_size, local_data, is_elem_based, longname, units, ncat)
    type(fesom_file_type), intent(inout)            :: this
    character(len=*)     , intent(in)               :: name
    integer              , intent(in)               :: dim_indices(:)
    integer              , intent(in)               :: global_level_data_size
    real(kind=8)         , intent(inout), target    :: local_data(:,:) ! todo: be able to set precision?
    logical              , intent(in)               :: is_elem_based
    character(len=*)     , intent(in)               :: units, longname
    integer              , intent(in)   , optional  :: ncat
    
    ! EO parameters
    integer var_index

    var_index = this%add_var_double(name, dim_indices)
    call this%add_var_att(var_index, "units", units)
    call this%add_var_att(var_index, "long_name", longname)
    
    call assert(this%nvar_infos < size(this%var_infos), __LINE__)
    this%nvar_infos = this%nvar_infos+1
    this%var_infos(this%nvar_infos)%var_index = var_index
    this%var_infos(this%nvar_infos)%external_local_data_ptr => local_data
    this%var_infos(this%nvar_infos)%global_level_data_size = global_level_data_size
    this%var_infos(this%nvar_infos)%is_elem_based = is_elem_based
    this%var_infos(this%nvar_infos)%varname = name
    if (present(ncat)) this%var_infos(this%nvar_infos)%is_icepackvar2=.true.
    
  end subroutine
  
  
  subroutine close_file(this)
    class(fesom_file_type), intent(inout) :: this

    if(this%thread_running) call this%thread%join()
    this%thread_running = .false.    
   
    this%rec_cnt = -1 ! reset state (should probably be done in all the open_ procedures, not here)
    call this%netcdf_file_type%close_file()
  end subroutine  


  subroutine assert(val, line)
    logical, intent(in) :: val
    integer, intent(in) :: line
    ! EO parameters
    if(.not. val) then
      print *, "error in line ",line, __FILE__
      stop 1
    end if
  end subroutine


  subroutine assert_nc(status, line)
    integer, intent(in) :: status
    integer, intent(in) :: line
    ! EO parameters
    include "netcdf.inc"
    if(status /= nf_noerr) then
      print *, "error in line ",line, __FILE__, ' ', trim(nf_strerror(status))
      stop 1
    endif   
  end subroutine

end module
