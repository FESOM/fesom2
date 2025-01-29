  integer,          pointer     :: MPI_COMM_FESOM       ! FESOM communicator (for ocean only runs if often a copy of MPI_COMM_WORLD)
  integer,          pointer     :: MPI_COMM_FESOM_IB    ! FESOM communicator copy for icebergs LA: 2023-05-22
  type(com_struct), pointer     :: com_nod2D
  type(com_struct), pointer     :: com_elem2D
  type(com_struct), pointer     :: com_elem2D_full
  integer                       :: ub, lb ! to work with r(s)_mpitype_elem3D(nod3D)

  integer, dimension(:),     pointer :: s_mpitype_edge2D, r_mpitype_edge2D
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D, r_mpitype_elem2D
  integer, dimension(:),     pointer :: s_mpitype_elem2D_full_i, r_mpitype_elem2D_full_i
  integer, dimension(:,:),   pointer :: s_mpitype_elem2D_full,   r_mpitype_elem2D_full
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D,        r_mpitype_elem3D
  integer, dimension(:,:,:), pointer :: s_mpitype_elem3D_full,   r_mpitype_elem3D_full

  integer, dimension(:),     pointer :: s_mpitype_nod2D, r_mpitype_nod2D
  integer, dimension(:),     pointer :: s_mpitype_nod2D_i, r_mpitype_nod2D_i
  integer, dimension(:,:,:), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

  integer, pointer   :: MPIERR
  integer, pointer   :: MPIERR_IB       ! copy for icebergs LA: 2023-05-22
  integer, pointer   :: npes
  integer, pointer   :: mype
  integer, pointer   :: maxPEnum

  integer, dimension(:), pointer     :: part

  ! Mesh partition
  integer, pointer                :: myDim_nod2D, eDim_nod2D
  integer, dimension(:), pointer  :: myList_nod2D
  integer, pointer                :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
  integer, dimension(:), pointer  :: myList_elem2D
  integer, pointer                :: myDim_edge2D, eDim_edge2D
  integer, dimension(:), pointer  :: myList_edge2D

  integer, pointer                :: pe_status

  integer, dimension(:), pointer  ::  remPtr_nod2D(:),  remList_nod2D(:)
  integer, dimension(:), pointer  ::  remPtr_elem2D(:), remList_elem2D(:)

  logical, pointer                :: elem_full_flag
