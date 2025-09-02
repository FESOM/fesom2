MPI_COMM_FESOM_WORLD                => partit%MPI_COMM_FESOM_WORLD
MPI_COMM_FESOM_SAME_RANK_IN_GROUPS  => partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS
MPI_COMM_FESOM  => partit%MPI_COMM_FESOM
MPI_COMM_FESOM_IB       => partit%MPI_COMM_FESOM_IB
com_nod2D       => partit%com_nod2D
com_elem2D      => partit%com_elem2D
com_elem2D_full => partit%com_elem2D_full
myDim_nod2D     => partit%myDim_nod2D
eDim_nod2D      => partit%eDim_nod2D
myDim_elem2D    => partit%myDim_elem2D
eDim_elem2D     => partit%eDim_elem2D
eXDim_elem2D    => partit%eXDim_elem2D
myDim_edge2D    => partit%myDim_edge2D
eDim_edge2D     => partit%eDim_edge2D
pe_status       => partit%pe_status
elem_full_flag  => partit%elem_full_flag
MPIERR          => partit%MPIERR
MPIERR_IB       => partit%MPIERR_IB
npes            => partit%npes
mype            => partit%mype
my_fesom_group  => partit%my_fesom_group
maxPEnum        => partit%maxPEnum
part            => partit%part

lb=lbound(partit%s_mpitype_elem3D, 2)
ub=ubound(partit%s_mpitype_elem3D, 2)

myList_nod2D (1:myDim_nod2D +eDim_nod2D)                => partit%myList_nod2D(:)
myList_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)  => partit%myList_elem2D(:)
myList_edge2D(1:myDim_edge2D+eDim_edge2D)               => partit%myList_edge2D(:)

if (allocated(partit%remPtr_nod2D)) then
   remPtr_nod2D  (1:npes)                => partit%remPtr_nod2D(:)
end if
if (allocated(partit%remList_nod2D)) then
   remList_nod2D (1:remPtr_nod2D(npes))  => partit%remList_nod2D(:)
end if 

if (allocated(partit%remPtr_elem2D)) then
    remPtr_elem2D (1:npes)                => partit%remPtr_elem2D(:)
end if
if (allocated(partit%remList_elem2D)) then
    remList_elem2D(1:remPtr_elem2D(npes)) => partit%remList_elem2D(:)
end if 

if (allocated(partit%s_mpitype_elem2D)) then
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

   r_mpitype_elem3D(1:com_elem2D%rPEnum, lb:ub, 1:4)           => partit%r_mpitype_elem3D(:,:,:)
   r_mpitype_elem3D_full(1:com_elem2D_full%rPEnum, lb:ub, 1:4) => partit%r_mpitype_elem3D_full(:,:,:)
end if 

if (allocated(partit%s_mpitype_nod2D)) then
   s_mpitype_nod2D(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D(:)
   r_mpitype_nod2D(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D(:)

   s_mpitype_nod2D_i(1:com_nod2D%sPEnum) => partit%s_mpitype_nod2D_i(:)
   r_mpitype_nod2D_i(1:com_nod2D%rPEnum) => partit%r_mpitype_nod2D_i(:)

   s_mpitype_nod3D(1:com_nod2D%sPEnum, lb:ub, 1:3) => partit%s_mpitype_nod3D(:,:,:)
   r_mpitype_nod3D(1:com_nod2D%rPEnum, lb:ub, 1:3) => partit%r_mpitype_nod3D(:,:,:)
end if 