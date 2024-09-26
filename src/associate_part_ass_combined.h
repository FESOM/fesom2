
! normally in associate_part_def.h
integer                       :: ub, lb ! to work with r(s)_mpitype_elem3D(nod3D)

! Leave these as pointers because 1. they are not used inside kernels and more importantly 2. the bounds remapping
! used for the 3d mpi datatypes gets broken if switched to the associate block.
integer, dimension(:), pointer  ::  remPtr_nod2D(:),  remList_nod2D(:)
integer, dimension(:), pointer  ::  remPtr_elem2D(:), remList_elem2D(:)

integer, dimension(:),     pointer :: s_mpitype_edge2D, r_mpitype_edge2D
integer, dimension(:,:),   pointer :: s_mpitype_elem2D, r_mpitype_elem2D
integer, dimension(:),     pointer :: s_mpitype_elem2D_full_i, r_mpitype_elem2D_full_i
integer, dimension(:,:),   pointer :: s_mpitype_elem2D_full,   r_mpitype_elem2D_full
integer, dimension(:,:,:), pointer :: s_mpitype_elem3D,        r_mpitype_elem3D
integer, dimension(:,:,:), pointer :: s_mpitype_elem3D_full,   r_mpitype_elem3D_full

integer, dimension(:),     pointer :: s_mpitype_nod2D, r_mpitype_nod2D
integer, dimension(:),     pointer :: s_mpitype_nod2D_i, r_mpitype_nod2D_i
integer, dimension(:,:,:), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

lb=lbound(partit%s_mpitype_elem3D, 2)
ub=ubound(partit%s_mpitype_elem3D, 2)


associate ( &
MPI_COMM_FESOM  => partit%MPI_COMM_FESOM, &
com_nod2D       => partit%com_nod2D, &
com_elem2D      => partit%com_elem2D, &
com_elem2D_full => partit%com_elem2D_full, &
myDim_nod2D     => partit%myDim_nod2D, &
eDim_nod2D      => partit%eDim_nod2D, &
myDim_elem2D    => partit%myDim_elem2D, &
eDim_elem2D     => partit%eDim_elem2D, &
eXDim_elem2D    => partit%eXDim_elem2D, &
myDim_edge2D    => partit%myDim_edge2D, &
eDim_edge2D     => partit%eDim_edge2D, &
pe_status       => partit%pe_status, &
elem_full_flag  => partit%elem_full_flag, &
MPIERR          => partit%MPIERR, &
npes            => partit%npes, &
mype            => partit%mype, &
maxPEnum        => partit%maxPEnum, &
part            => partit%part,     &

myList_nod2D                 => partit%myList_nod2D(1:partit%myDim_nod2D + partit%eDim_nod2D), &
myList_elem2D  => partit%myList_elem2D(1:partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D), &
myList_edge2D               => partit%myList_edge2D(1:partit%myDim_edge2D + partit%eDim_edge2D) )


if (allocated(partit%remPtr_nod2D)) then
   remPtr_nod2D  (1:npes)                => partit%remPtr_nod2D(:)
   remList_nod2D (1:remPtr_nod2D(npes))  => partit%remList_nod2D(:)
end if

if (allocated(partit%remPtr_elem2D)) then
    remPtr_elem2D (1:npes)                => partit%remPtr_elem2D(:)
    remList_elem2D(1:remPtr_elem2D(npes)) => partit%remList_elem2D(:)
end if


s_mpitype_elem2D(1:com_elem2D%sPEnum, 1:4) => partit%s_mpitype_elem2D(:,:)
r_mpitype_elem2D(1:com_elem2D%rPEnum, 1:4) => partit%r_mpitype_elem2D(:,:)

s_mpitype_elem2D_full_i(1:com_elem2D_full%sPEnum) => partit%s_mpitype_elem2D_full_i(:)
r_mpitype_elem2D_full_i(1:com_elem2D_full%rPEnum) => partit%r_mpitype_elem2D_full_i(:)

s_mpitype_elem2D_full(1:com_elem2D_full%sPEnum, 1:4) => partit%s_mpitype_elem2D_full(:,:)
r_mpitype_elem2D_full(1:com_elem2D_full%rPEnum, 1:4) => partit%r_mpitype_elem2D_full(:,:)

s_mpitype_elem3D(1:com_elem2D%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D(:,:,:)
r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D(:,:,:)

s_mpitype_elem3D_full(1:com_elem2D_full%sPEnum, lb:ub, 1:4) => partit%s_mpitype_elem3D_full(:,:,:)
r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)

s_mpitype_nod2D(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D(:)
r_mpitype_nod2D(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D(:)

s_mpitype_nod2D_i(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D_i(:)
r_mpitype_nod2D_i(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D_i(:)

s_mpitype_nod3D(1:com_nod2D%sPEnum, lb:ub, 1:3) => partit%s_mpitype_nod3D(:,:,:)
r_mpitype_nod3D(1:com_nod2D%rPEnum, lb:ub, 1:3) => partit%r_mpitype_nod3D(:,:,:)

