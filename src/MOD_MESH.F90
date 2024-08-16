!==========================================================
MODULE MOD_MESH
USE O_PARAM
USE MOD_WRITE_BINARY_ARRAYS
USE MOD_READ_BINARY_ARRAYS
USE,     intrinsic    :: ISO_FORTRAN_ENV, only : int32
IMPLICIT NONE
SAVE
integer, parameter    :: MAX_ADJACENT=32 ! Max allowed number of adjacent nodes

TYPE SPARSE_MATRIX
     integer :: nza
     integer :: dim
     real(kind=WP),  allocatable,   dimension(:) :: values
     integer(int32), allocatable,   dimension(:) :: colind
     integer(int32), allocatable,   dimension(:) :: rowptr
     integer(int32), allocatable,   dimension(:) :: colind_loc
     integer(int32), allocatable,   dimension(:) :: rowptr_loc
     real(kind=WP),  allocatable,   dimension(:) :: pr_values !preconditioner values
END TYPE SPARSE_MATRIX

TYPE T_MESH
integer                                     :: nod2D           ! the number of 2D nodes
real(kind=WP)                               :: ocean_area
real(kind=WP)                               :: ocean_areawithcav
real(kind=WP), allocatable, dimension(:,:)  :: coord_nod2D
real(kind=WP), allocatable, dimension(:,:)  :: geo_coord_nod2D
integer                                     :: edge2D       ! the number of 2D edges
integer                                     :: edge2D_in    ! the number of internal 2D edges
integer                                     :: elem2D       ! the number of 2D elements
integer, allocatable, dimension(:,:)        :: elem2D_nodes ! elem2D_nodes(:,n) lists; 3 nodes of element n
integer, allocatable, dimension(:,:)        :: edges        ! edge(:,n) lists 2 nodes; edge n
integer, allocatable, dimension(:,:)        :: edge_tri     ! edge_tri(:,n) lists 2
                                                            ! elements containing edge n: the first one is to left
                                                            ! of the line directed to the second node
integer,       allocatable, dimension(:,:)  :: elem_edges   ! elem_edges(:,n) are edges of element n.
real(kind=WP), allocatable, dimension(:)    :: elem_area
real(kind=WP), allocatable, dimension(:,:)  :: edge_dxdy
real(kind=WP), allocatable, dimension(:,:)  :: edge_cross_dxdy
real(kind=WP), allocatable, dimension(:)    :: elem_cos
real(kind=WP), allocatable, dimension(:)    :: metric_factor
integer,       allocatable, dimension(:,:)  :: elem_neighbors
integer,       allocatable, dimension(:,:)  :: nod_in_elem2D
real(kind=WP), allocatable, dimension(:,:)  :: x_corners
real(kind=WP), allocatable, dimension(:,:)  :: y_corners ! cornes for the scalar points
integer,       allocatable, dimension(:)    :: nod_in_elem2D_num
real(kind=WP), allocatable, dimension(:)    :: depth                ! depth(n) is the depths at node n
real(kind=WP), allocatable, dimension(:,:)  :: gradient_vec
                                                           ! coefficients of linear reconstruction
                                                           ! of velocities on elements
real(kind=WP), allocatable, dimension(:,:)  :: gradient_sca ! Coefficients to compute gradient of scalars
                                                           ! on elements
INTEGER,       ALLOCATABLE, DIMENSION(:)    :: bc_index_nod2D(:)
                                                           ! vertical structure
!
!
!___vertical mesh info__________________________________________________________
! total number of layers
integer                                     :: nl

! initial layer, mid-depth layer and element depth
real(kind=WP), allocatable, dimension(:)    :: zbar
real(kind=WP), allocatable, dimension(:)    :: Z
real(kind=WP), allocatable, dimension(:)    :: elem_depth

! upper boudnary index of all vertical vertice/element loops, default==1 but when
! cavity is used becomes index of cavity-ocean boundary at vertices and elements
integer,       allocatable, dimension(:)    :: ulevels
integer,       allocatable, dimension(:)    :: ulevels_nod2D
integer,       allocatable, dimension(:)    :: ulevels_nod2D_max

! number of levels at elem and vertices considering bottom topography
integer,       allocatable, dimension(:)    :: nlevels
integer,       allocatable, dimension(:)    :: nlevels_nod2D
integer,       allocatable, dimension(:)    :: nlevels_nod2D_min

!
!
!___horizontal mesh info________________________________________________________
real(kind=WP), allocatable, dimension(:,:)  :: area
real(kind=WP), allocatable, dimension(:,:)  :: area_inv
real(kind=WP), allocatable, dimension(:,:)  :: areasvol
real(kind=WP), allocatable, dimension(:,:)  :: areasvol_inv
real(kind=WP), allocatable, dimension(:)    :: mesh_resolution

!
!
!___cavity mesh info____________________________________________________________
! level index of cavity-ocean boundary at vertices and elements
! --> see: ulevels, ulevels_nod2D (fvom_main)

! vertice/element yes=1/no=0 flag if cavity exists
integer,       allocatable, dimension(:)    :: cavity_flag_n
integer,       allocatable, dimension(:)    :: cavity_flag_e

! depth of cavity-ocean interface
real(kind=WP), allocatable, dimension(:)    :: cavity_depth


real(kind=WP), allocatable, dimension(:,:)  :: cavity_nrst_cavlpnt_xyz

!
!
!___Elevation stiffness matrix__________________________________________________
type(sparse_matrix)                         :: ssh_stiff

!#if defined (__oasis)
real(kind=WP), allocatable, dimension(:)    :: lump2d_south
real(kind=WP), allocatable, dimension(:)    :: lump2d_north
integer,       allocatable, dimension(:)    :: ind_south
integer,       allocatable, dimension(:)    :: ind_north
!#endif

integer                                       :: nn_size
integer, allocatable, dimension(:)            :: nn_num
integer, allocatable, dimension(:,:)          :: nn_pos

!_______________________________________________________________________________
! Arrays added for ALE implementation:
! --> layer thinkness at node and depthlayer for t=n and t=n+1
!------------------------------
real(kind=WP), allocatable,dimension(:,:)   :: hnode
real(kind=WP), allocatable,dimension(:,:)   :: hnode_new
real(kind=WP), allocatable,dimension(:,:)   :: zbar_3d_n
real(kind=WP), allocatable,dimension(:,:)   :: Z_3d_n
! LA 2023-01-31 add icebergs
real(kind=WP), allocatable,dimension(:,:)   :: Z_3d_n_ib

!------------------------------
! --> layer thinkness at elements, interpolated from hnode
real(kind=WP), allocatable,dimension(:,:)   :: helem

! --> thinkness of bottom elem (important for partial cells)
real(kind=WP), allocatable,dimension(:)     :: bottom_elem_thickness
real(kind=WP), allocatable,dimension(:)     :: bottom_node_thickness

! --> The increment of total fluid depth on elements. It is used to update the matrix
real(kind=WP), allocatable,dimension(:)     :: dhe

! --> hbar, hbar_old: correspond to the elevation, but on semi-integer time steps.
real(kind=WP), allocatable,dimension(:)     :: hbar
real(kind=WP), allocatable,dimension(:)     :: hbar_old

! --> auxiliary array to store depth of layers and depth of mid level due to changing
!     layer thinkness at every node
!real(kind=WP), allocatable,dimension(:)     :: zbar_n
!real(kind=WP), allocatable,dimension(:)     :: Z_n

! new bottom depth at node and element due to partial cells
real(kind=WP), allocatable,dimension(:)     :: zbar_n_bot
real(kind=WP), allocatable,dimension(:)     :: zbar_e_bot

! new depth of cavity-ocean interface at node and element due to partial cells
real(kind=WP), allocatable,dimension(:)     :: zbar_n_srf
real(kind=WP), allocatable,dimension(:)     :: zbar_e_srf

character(:), allocatable :: representative_checksum

!
!
!___coriolis force______________________________________________________________
real(kind=WP), allocatable, dimension(:)    :: coriolis
real(kind=WP), allocatable, dimension(:)    :: coriolis_node

contains
#if defined(__PGI)
  procedure, private write_t_mesh
  procedure, private read_t_mesh
#else
  procedure write_t_mesh
  procedure read_t_mesh
#endif
  generic :: write(unformatted) => write_t_mesh
  generic :: read(unformatted)  => read_t_mesh
END TYPE T_MESH

contains

! Unformatted writing for t_mesh
subroutine write_t_mesh(mesh, unit, iostat, iomsg)
    IMPLICIT NONE
    class(t_mesh), intent(in)         :: mesh
    integer,      intent(in)          :: unit
    integer,      intent(out)         :: iostat
    character(*), intent(inout)       :: iomsg
    integer                           :: i, j, k
    integer                           :: s1, s2, s3
    ! write records (giving sizes for the allocation for arrays)
    write(unit, iostat=iostat, iomsg=iomsg) mesh%nod2D
    write(unit, iostat=iostat, iomsg=iomsg) mesh%ocean_area
    write(unit, iostat=iostat, iomsg=iomsg) mesh%ocean_areawithcav
    write(unit, iostat=iostat, iomsg=iomsg) mesh%edge2D
    write(unit, iostat=iostat, iomsg=iomsg) mesh%edge2D_in
    write(unit, iostat=iostat, iomsg=iomsg) mesh%elem2D
    call write_bin_array(mesh%elem2D_nodes, unit, iostat, iomsg)
    call write_bin_array(mesh%edges,        unit, iostat, iomsg)
    call write_bin_array(mesh%edge_tri,     unit, iostat, iomsg)
    call write_bin_array(mesh%elem_edges,   unit, iostat, iomsg)
    call write_bin_array(mesh%elem_area,    unit, iostat, iomsg)
    call write_bin_array(mesh%edge_dxdy,    unit, iostat, iomsg)

    call write_bin_array(mesh%edge_cross_dxdy,     unit, iostat, iomsg)
    call write_bin_array(mesh%elem_cos,            unit, iostat, iomsg)
    call write_bin_array(mesh%metric_factor,       unit, iostat, iomsg)
    call write_bin_array(mesh%elem_neighbors,      unit, iostat, iomsg)
    call write_bin_array(mesh%nod_in_elem2D,       unit, iostat, iomsg)
    call write_bin_array(mesh%x_corners,           unit, iostat, iomsg)
    call write_bin_array(mesh%y_corners,           unit, iostat, iomsg)
    call write_bin_array(mesh%nod_in_elem2D_num,   unit, iostat, iomsg)
    call write_bin_array(mesh%depth,               unit, iostat, iomsg)
    call write_bin_array(mesh%gradient_vec,        unit, iostat, iomsg)
    call write_bin_array(mesh%gradient_sca,        unit, iostat, iomsg)
    call write_bin_array(mesh%bc_index_nod2D,      unit, iostat, iomsg)

    write(unit, iostat=iostat, iomsg=iomsg) mesh%nl

    call write_bin_array(mesh%zbar,              unit, iostat, iomsg)
    call write_bin_array(mesh%Z,                 unit, iostat, iomsg)
    call write_bin_array(mesh%elem_depth,        unit, iostat, iomsg)
    call write_bin_array(mesh%ulevels,           unit, iostat, iomsg)
    call write_bin_array(mesh%ulevels_nod2D,     unit, iostat, iomsg)
    call write_bin_array(mesh%ulevels_nod2D_max, unit, iostat, iomsg)
    call write_bin_array(mesh%nlevels,           unit, iostat, iomsg)
    call write_bin_array(mesh%nlevels_nod2D,     unit, iostat, iomsg)
    call write_bin_array(mesh%nlevels_nod2D_min, unit, iostat, iomsg)
    call write_bin_array(mesh%area,              unit, iostat, iomsg)
    call write_bin_array(mesh%area_inv,          unit, iostat, iomsg)
    call write_bin_array(mesh%areasvol,          unit, iostat, iomsg)
    call write_bin_array(mesh%areasvol_inv,      unit, iostat, iomsg)
    call write_bin_array(mesh%mesh_resolution,   unit, iostat, iomsg)

    call write_bin_array(mesh%cavity_flag_n,     unit, iostat, iomsg)
    call write_bin_array(mesh%cavity_flag_e,     unit, iostat, iomsg)
    call write_bin_array(mesh%cavity_depth,      unit, iostat, iomsg)
    call write_bin_array(mesh%cavity_nrst_cavlpnt_xyz,   unit, iostat, iomsg)

    write(unit, iostat=iostat, iomsg=iomsg) mesh%ssh_stiff%dim
    write(unit, iostat=iostat, iomsg=iomsg) mesh%ssh_stiff%nza

    call write_bin_array(mesh%ssh_stiff%rowptr,     unit, iostat, iomsg)
    call write_bin_array(mesh%ssh_stiff%colind,     unit, iostat, iomsg)
    call write_bin_array(mesh%ssh_stiff%values,     unit, iostat, iomsg)
    call write_bin_array(mesh%ssh_stiff%colind_loc, unit, iostat, iomsg)
    call write_bin_array(mesh%ssh_stiff%rowptr_loc, unit, iostat, iomsg)

    call write_bin_array(mesh%lump2d_south,            unit, iostat, iomsg)
    call write_bin_array(mesh%lump2d_north,            unit, iostat, iomsg)
    call write_bin_array(mesh%ind_south,               unit, iostat, iomsg)
    call write_bin_array(mesh%ind_north,               unit, iostat, iomsg)
    write(unit, iostat=iostat, iomsg=iomsg) mesh%nn_size
    call write_bin_array(mesh%nn_num,                  unit, iostat, iomsg)
    call write_bin_array(mesh%nn_pos,                  unit, iostat, iomsg)
    call write_bin_array(mesh%hnode,                   unit, iostat, iomsg)
    call write_bin_array(mesh%hnode_new,               unit, iostat, iomsg)
    call write_bin_array(mesh%zbar_3d_n,               unit, iostat, iomsg)
    call write_bin_array(mesh%Z_3d_n,                  unit, iostat, iomsg)
    call write_bin_array(mesh%Z_3d_n_ib,               unit, iostat, iomsg)
    call write_bin_array(mesh%helem,                   unit, iostat, iomsg)
    call write_bin_array(mesh%bottom_elem_thickness,   unit, iostat, iomsg)
    call write_bin_array(mesh%bottom_node_thickness,   unit, iostat, iomsg)
    call write_bin_array(mesh%dhe,                     unit, iostat, iomsg)
    call write_bin_array(mesh%hbar,                    unit, iostat, iomsg)
    call write_bin_array(mesh%hbar_old,                unit, iostat, iomsg)
!   call write_bin_array(mesh%zbar_n,                  unit, iostat, iomsg)
!   call write_bin_array(mesh%Z_n,                     unit, iostat, iomsg)
    call write_bin_array(mesh%zbar_n_bot,              unit, iostat, iomsg)
    call write_bin_array(mesh%zbar_e_bot,              unit, iostat, iomsg)
    call write_bin_array(mesh%zbar_n_srf,              unit, iostat, iomsg)
    call write_bin_array(mesh%zbar_e_srf,              unit, iostat, iomsg)
!   call write_bin_array(mesh%representative_checksum, unit, iostat, iomsg)
    call write_bin_array(mesh%coriolis,                unit, iostat, iomsg)
    call write_bin_array(mesh%coriolis_node,           unit, iostat, iomsg)

end subroutine write_t_mesh

! Unformatted reading for t_mesh
subroutine read_t_mesh(mesh, unit, iostat, iomsg)
    IMPLICIT NONE
    class(t_mesh), intent(inout)       :: mesh
    integer,       intent(in)          :: unit
    integer,       intent(out)         :: iostat
    character(*),  intent(inout)       :: iomsg

    ! write records (giving sizes for the allocation for arrays)
    read(unit, iostat=iostat, iomsg=iomsg) mesh%nod2D
    read(unit, iostat=iostat, iomsg=iomsg) mesh%ocean_area
    read(unit, iostat=iostat, iomsg=iomsg) mesh%ocean_areawithcav
    read(unit, iostat=iostat, iomsg=iomsg) mesh%edge2D
    read(unit, iostat=iostat, iomsg=iomsg) mesh%edge2D_in
    read(unit, iostat=iostat, iomsg=iomsg) mesh%elem2D

    call read_bin_array(mesh%elem2D_nodes, unit, iostat, iomsg)
    call read_bin_array(mesh%edges,        unit, iostat, iomsg)
    call read_bin_array(mesh%edge_tri,     unit, iostat, iomsg)
    call read_bin_array(mesh%elem_edges,   unit, iostat, iomsg)
    call read_bin_array(mesh%elem_area,    unit, iostat, iomsg)
    call read_bin_array(mesh%edge_dxdy,    unit, iostat, iomsg)

    call read_bin_array(mesh%edge_cross_dxdy,     unit, iostat, iomsg)
    call read_bin_array(mesh%elem_cos,            unit, iostat, iomsg)
    call read_bin_array(mesh%metric_factor,       unit, iostat, iomsg)
    call read_bin_array(mesh%elem_neighbors,      unit, iostat, iomsg)
    call read_bin_array(mesh%nod_in_elem2D,       unit, iostat, iomsg)
    call read_bin_array(mesh%x_corners,           unit, iostat, iomsg)
    call read_bin_array(mesh%y_corners,           unit, iostat, iomsg)
    call read_bin_array(mesh%nod_in_elem2D_num,   unit, iostat, iomsg)
    call read_bin_array(mesh%depth,               unit, iostat, iomsg)
    call read_bin_array(mesh%gradient_vec,        unit, iostat, iomsg)
    call read_bin_array(mesh%gradient_sca,        unit, iostat, iomsg)
    call read_bin_array(mesh%bc_index_nod2D,      unit, iostat, iomsg)

    read(unit, iostat=iostat, iomsg=iomsg) mesh%nl

    call read_bin_array(mesh%zbar,              unit, iostat, iomsg)
    call read_bin_array(mesh%Z,                 unit, iostat, iomsg)
    call read_bin_array(mesh%elem_depth,        unit, iostat, iomsg)
    call read_bin_array(mesh%ulevels,           unit, iostat, iomsg)
    call read_bin_array(mesh%ulevels_nod2D,     unit, iostat, iomsg)
    call read_bin_array(mesh%ulevels_nod2D_max, unit, iostat, iomsg)
    call read_bin_array(mesh%nlevels,           unit, iostat, iomsg)
    call read_bin_array(mesh%nlevels_nod2D,     unit, iostat, iomsg)
    call read_bin_array(mesh%nlevels_nod2D_min, unit, iostat, iomsg)
    call read_bin_array(mesh%area,              unit, iostat, iomsg)
    call read_bin_array(mesh%area_inv,          unit, iostat, iomsg)
    call read_bin_array(mesh%areasvol,          unit, iostat, iomsg)
    call read_bin_array(mesh%areasvol_inv,      unit, iostat, iomsg)
    call read_bin_array(mesh%mesh_resolution,   unit, iostat, iomsg)

    call read_bin_array(mesh%cavity_flag_n,     unit, iostat, iomsg)
    call read_bin_array(mesh%cavity_flag_e,     unit, iostat, iomsg)
    call read_bin_array(mesh%cavity_depth,      unit, iostat, iomsg)
    call read_bin_array(mesh%cavity_nrst_cavlpnt_xyz,   unit, iostat, iomsg)

    read(unit, iostat=iostat, iomsg=iomsg) mesh%ssh_stiff%dim
    read(unit, iostat=iostat, iomsg=iomsg) mesh%ssh_stiff%nza

    call read_bin_array(mesh%ssh_stiff%rowptr,     unit, iostat, iomsg)
    call read_bin_array(mesh%ssh_stiff%colind,     unit, iostat, iomsg)
    call read_bin_array(mesh%ssh_stiff%values,     unit, iostat, iomsg)
    call read_bin_array(mesh%ssh_stiff%colind_loc, unit, iostat, iomsg)
    call read_bin_array(mesh%ssh_stiff%rowptr_loc, unit, iostat, iomsg)

    call read_bin_array(mesh%lump2d_south,            unit, iostat, iomsg)
    call read_bin_array(mesh%lump2d_north,            unit, iostat, iomsg)
    call read_bin_array(mesh%ind_south,               unit, iostat, iomsg)
    call read_bin_array(mesh%ind_north,               unit, iostat, iomsg)
    read(unit, iostat=iostat, iomsg=iomsg) mesh%nn_size
    call read_bin_array(mesh%nn_num,                  unit, iostat, iomsg)
    call read_bin_array(mesh%nn_pos,                  unit, iostat, iomsg)
    call read_bin_array(mesh%hnode,                   unit, iostat, iomsg)
    call read_bin_array(mesh%hnode_new,               unit, iostat, iomsg)
    call read_bin_array(mesh%zbar_3d_n,               unit, iostat, iomsg)
    call read_bin_array(mesh%Z_3d_n,                  unit, iostat, iomsg)
    call read_bin_array(mesh%Z_3d_n_ib,               unit, iostat, iomsg)
    call read_bin_array(mesh%helem,                   unit, iostat, iomsg)
    call read_bin_array(mesh%bottom_elem_thickness,   unit, iostat, iomsg)
    call read_bin_array(mesh%bottom_node_thickness,   unit, iostat, iomsg)
    call read_bin_array(mesh%dhe,                     unit, iostat, iomsg)
    call read_bin_array(mesh%hbar,                    unit, iostat, iomsg)
    call read_bin_array(mesh%hbar_old,                unit, iostat, iomsg)
!   call read_bin_array(mesh%zbar_n,                  unit, iostat, iomsg)
!   call read_bin_array(mesh%Z_n,                     unit, iostat, iomsg)
    call read_bin_array(mesh%zbar_n_bot,              unit, iostat, iomsg)
    call read_bin_array(mesh%zbar_e_bot,              unit, iostat, iomsg)
    call read_bin_array(mesh%zbar_n_srf,              unit, iostat, iomsg)
    call read_bin_array(mesh%zbar_e_srf,              unit, iostat, iomsg)
!   call read_bin_array(mesh%representative_checksum, unit, iostat, iomsg)
    call read_bin_array(mesh%coriolis,                unit, iostat, iomsg)
    call read_bin_array(mesh%coriolis_node,           unit, iostat, iomsg)

end subroutine read_t_mesh
end module MOD_MESH
!==========================================================
