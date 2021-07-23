module io_mesh_info
use g_PARSUP
use MOD_MESH
use g_config
use g_comm_auto
use o_ARRAYS

implicit none
#include "netcdf.inc"
private
public :: write_mesh_info
INTERFACE my_put_vara
            MODULE PROCEDURE my_put_vara_double_1D
            MODULE PROCEDURE my_put_vara_double_2D
            MODULE PROCEDURE my_put_vara_int_1D
            MODULE PROCEDURE my_put_vara_int_2D
END INTERFACE
!INTERFACE
!subroutine my_def_var(ncid, short_name, vtype, dsize, dids, id, att_text)
!           integer, intent(in)        :: ncid, dsize, dids(dsize), vtype
!           character(*), intent(in)   :: short_name, att_text
!           integer,      intent(inout):: id
!end subroutine my_def_var
!
!subroutine my_def_dim(ncid, short_name, value, id)
!           integer,      intent(in)   :: ncid, value
!           character(*), intent(in)   :: short_name
!           integer,      intent(inout):: id
!end subroutine my_def_dim
!END INTERFACE
contains
!
!-------------------------------------------------------------------------
! this routine stores most of metadata used in FESOM. Shall be called at the cold start once during the simulation. 
! info: fesom.mesh.diag.nc is 77MB for the CORE II mesh with 47 vertical levels
subroutine write_mesh_info(mesh)
implicit none
  type(t_mesh), intent(in)  , target :: mesh
  integer                    :: status, ncid, j
  integer                    :: nod_n_id, elem_n_id, edge_n_id, nod_part_id, elem_part_id
  integer                    :: nl_id, nl1_id
  integer                    :: id_2, id_3, id_4, id_N
  integer                    :: i, k, N_max
  integer                    :: nod_area_id, elem_area_id
  integer, target            :: z_id, zbar_id
  integer                    :: edges_id, edge_tri_id, edge_cross_dxdy_id
  integer                    :: nlevels_nod2D_id, nlevels_id
  integer                    :: nod_in_elem2D_num_id, nod_in_elem2D_id
  integer                    :: gradient_sca_x_id, gradient_sca_y_id
  integer                    :: zbar_e_bot_id,zbar_n_bot_id
  integer                    :: elem_id
  integer                    :: nod_id
  character(100)             :: longname
  character(2000)            :: filename
  real(kind=WP), allocatable :: rbuffer(:), lrbuffer(:)
  integer, allocatable       :: ibuffer(:), lbuffer(:)
  character(2000)            :: att_text, short_name
  integer                    :: vtype
  integer, pointer           :: pid

#include "associate_mesh.h"

  
  call MPI_AllREDUCE(maxval(nod_in_elem2D_num), N_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_FESOM, MPIerr)

  filename=trim(ResultPath)//runid//'.mesh.diag.nc'
  call my_create(filename, IOR(NF_CLOBBER,IOR(NF_NETCDF4,NF_CLASSIC_MODEL)), ncid)

  !Define the dimensions
  call my_def_dim(ncid, 'nod2',    nod2D,  nod_n_id)
  call my_def_dim(ncid, 'edg_n',   edge2d, edge_n_id)
  call my_def_dim(ncid, 'elem',    elem2d, elem_n_id)
  call my_def_dim(ncid, 'nz',      nl,     nl_id)
  call my_def_dim(ncid, 'nz1',     nl-1,   nl1_id)
  call my_def_dim(ncid, 'n2',      2,      id_2)
  call my_def_dim(ncid, 'n3',      3,      id_3)
  call my_def_dim(ncid, 'n4',      4,      id_4)
  call my_def_dim(ncid, 'N',       N_max,  id_N)

  !Define the variables
  ! 1D
  call my_def_var(ncid, 'nz',                NF_DOUBLE, 1, (/nl_id /),    zbar_id,              'depth of levels'                       )
  call my_def_var(ncid, 'nz1',               NF_DOUBLE, 1, (/nl1_id/),    z_id,                 'depth of layers'                       )
  call my_def_var(ncid, 'elem_area',         NF_DOUBLE, 1, (/elem_n_id/), elem_area_id,         'element areas'                         )
  call my_def_var(ncid, 'nlevels_nod2D',     NF_INT,    1, (/nod_n_id/),  nlevels_nod2D_id,     'number of levels below nodes'          )
  call my_def_var(ncid, 'nlevels',           NF_INT,    1, (/elem_n_id/), nlevels_id,           'number of levels below elements'       )
  call my_def_var(ncid, 'nod_in_elem2D_num', NF_INT,    1, (/nod_n_id/),  nod_in_elem2D_num_id, 'number of elements containing the node')
  call my_def_var(ncid, 'nod_part',          NF_INT,    1, (/nod_n_id/),  nod_part_id,          'nodal partitioning at the cold start'  )
  call my_def_var(ncid, 'elem_part',         NF_INT,    1, (/elem_n_id/), elem_part_id,         'element partitioning at the cold start')
  call my_def_var(ncid, 'zbar_e_bottom',     NF_DOUBLE, 1, (/elem_n_id/), zbar_e_bot_id,        'element bottom depth')
  call my_def_var(ncid, 'zbar_n_bottom',     NF_DOUBLE, 1, (/nod_n_id/) , zbar_n_bot_id,        'nodal bottom depth')
  ! 2D
  call my_def_var(ncid, 'nod_area',          NF_DOUBLE, 2, (/nod_n_id, nl_id/), nod_area_id,        'nodal areas'                 )
  call my_def_var(ncid, 'elements',          NF_INT,    2, (/elem_n_id, id_3/), elem_id,            'elements'                    )
  call my_def_var(ncid, 'nodes',             NF_DOUBLE, 2, (/nod_n_id,  id_2/), nod_id,             'nodal geo. coordinates'      )
  call my_def_var(ncid, 'nod_in_elem2D',     NF_INT,    2, (/nod_n_id, id_N/),  nod_in_elem2D_id,   'elements containing the node')
  call my_def_var(ncid, 'edges',             NF_INT,    2, (/edge_n_id, id_2/), edges_id,           'edges'                       )
  call my_def_var(ncid, 'edge_tri',          NF_INT,    2, (/edge_n_id, id_2/), edge_tri_id,        'edge triangles'              )
  call my_def_var(ncid, 'edge_cross_dxdy',   NF_DOUBLE, 2, (/edge_n_id, id_4/), edge_cross_dxdy_id, 'edge cross distancess'       )
  call my_def_var(ncid, 'gradient_sca_x',    NF_DOUBLE, 2, (/id_3, elem_n_id/), gradient_sca_x_id,  'x component of a gradient at nodes of an element')
  call my_def_var(ncid, 'gradient_sca_y',    NF_DOUBLE, 2, (/id_3, elem_n_id/), gradient_sca_y_id,  'y component of a gradient at nodes of an element')
  call my_nf_enddef(ncid)

  ! vercical levels/layers
  call my_put_vara(ncid, zbar_id, 1, nl, zbar)
  call my_put_vara(ncid, z_id, 1, nl-1, Z)

  ! nodal areas
  allocate(rbuffer(nod2D))
  do k=1, nl
     call gather_nod(area(k, :), rbuffer)
     call my_put_vara(ncid, nod_area_id, (/1, k/), (/nod2D, 1/), rbuffer)
  end do
  deallocate(rbuffer)

  ! element areas
  allocate(rbuffer(elem2D))
  call gather_elem(elem_area(1:myDim_elem2D), rbuffer)
  call my_put_vara(ncid, elem_area_id, 1, elem2D, rbuffer)
  deallocate(rbuffer)

  ! elements
  allocate(ibuffer(elem2D))
  allocate(lbuffer(myDim_elem2D))  
  do i=1, 3
     do k=1, myDim_elem2D
        lbuffer(k)=myList_nod2D(elem2d_nodes(i, k))
     end do
     call gather_elem(lbuffer, ibuffer)
     call my_put_vara(ncid, elem_id, (/1, i/), (/elem2D, 1/), ibuffer)
  end do
  deallocate(lbuffer, ibuffer)

  ! number of levels below elements
  allocate(ibuffer(elem2D))
  call gather_elem(nlevels(1:myDim_elem2D), ibuffer)
  call my_put_vara(ncid, nlevels_id, 1, elem2D, ibuffer)
  deallocate(ibuffer)

  ! number of levels below nodes
  allocate(ibuffer(nod2D))
  call gather_nod(nlevels_nod2D(1:myDim_nod2D), ibuffer)
  call my_put_vara(ncid, nlevels_nod2D_id, 1, nod2D, ibuffer)
  deallocate(ibuffer)
  
  ! number of elements containing the node
  allocate(ibuffer(nod2D))
  call gather_nod(nod_in_elem2D_num(1:myDim_nod2D), ibuffer)
  call my_put_vara(ncid, nod_in_elem2D_num_id, 1, nod2D, ibuffer)
  deallocate(ibuffer)

  ! elements containing the node
  allocate(ibuffer(nod2D))
  allocate(lbuffer(myDim_nod2D))
  DO i=1, N_max
     lbuffer=0
        do k=1, myDim_nod2D
           if ((nod_in_elem2D_num(k)>=i)) then
              lbuffer(k)=myList_elem2D(nod_in_elem2D(i, k))
           end if
        end do
     call gather_nod(lbuffer, ibuffer)
     call my_put_vara(ncid, nod_in_elem2D_id, (/1, i/), (/nod2D, 1/), ibuffer)
  END DO
  deallocate(lbuffer, ibuffer)

  ! nodal partitioning
  allocate(ibuffer(nod2D))
  allocate(lbuffer(myDim_nod2D))  
  lbuffer=mype
  call gather_nod(lbuffer, ibuffer)
  call my_put_vara(ncid, nod_part_id, 1, nod2D, ibuffer)
  deallocate(lbuffer, ibuffer)

  ! element partitioning
  allocate(ibuffer(elem2D))
  allocate(lbuffer(myDim_elem2D))  
  lbuffer=mype
  call gather_elem(lbuffer, ibuffer)
  call my_put_vara(ncid, elem_part_id, 1, elem2D, ibuffer)
  deallocate(lbuffer, ibuffer)


  ! nodes (GEO coordinates)
  allocate(rbuffer(nod2D))
  do i=1, 2
     call gather_nod(geo_coord_nod2D(i, 1:myDim_nod2D), rbuffer)
     call my_put_vara(ncid, nod_id, (/1, i/), (/nod2D, 1/), rbuffer)
  end do
  deallocate(rbuffer)

  ! edges
  allocate(ibuffer(edge2D))
  allocate(lbuffer(myDim_edge2D))
  do i=1, 2
     do k=1, myDim_edge2D
        lbuffer(k)=myList_nod2D(edges(i, k))
     end do
     call gather_edge(lbuffer, ibuffer)
     call my_put_vara(ncid, edges_id, (/1, i/), (/edge2D, 1/), ibuffer)
  end do
  deallocate(lbuffer, ibuffer)

  ! edge triangles
  allocate(ibuffer(edge2D))
  allocate(lbuffer(myDim_edge2D))
  do i=1, 2
     do k=1, myDim_edge2D
        if (edge_tri(i,k) > 0) then
           lbuffer(k) = myList_elem2D(edge_tri(i,k))
        else
           lbuffer(k) = 0
        endif
     end do
     call gather_edge(lbuffer, ibuffer)
     call my_put_vara(ncid, edge_tri_id, (/1, i/), (/edge2D, 1/), ibuffer)
  end do
  deallocate(lbuffer, ibuffer)

  ! edge cross distances
  allocate(rbuffer(edge2D))
  allocate(lrbuffer(myDim_edge2D))
  do i=1, 4
     lrbuffer=edge_cross_dxdy(i, 1:myDim_edge2D)
     call gather_edge(lrbuffer, rbuffer)
     call my_put_vara(ncid, edge_cross_dxdy_id, (/1, i/), (/edge2D, 1/), rbuffer)
  end do
  deallocate(rbuffer, lrbuffer)


  ! X component of gadient at elements
  allocate(rbuffer(elem2D))
  do i=1, 3
     call gather_elem(gradient_sca(i, 1:myDim_elem2D), rbuffer)
     call my_put_vara(ncid, gradient_sca_x_id, (/4-i, 1/), (/1, elem2D/), rbuffer) ! (4-i), NETCDF will permute otherwise
  end do
  deallocate(rbuffer)

  ! Y component of gadient at elements
  allocate(rbuffer(elem2D))
  do i=1, 3
     call gather_elem(gradient_sca(i+3, 1:myDim_elem2D), rbuffer)
     call my_put_vara(ncid, gradient_sca_y_id, (/4-i, 1/), (/1, elem2D/), rbuffer) ! (4-i), NETCDF will permute otherwise
  end do
  deallocate(rbuffer)
  
  ! nodal bottom depth (take into account partial cells if used)
  allocate(rbuffer(nod2D))
  call gather_nod(zbar_n_bot(1:myDim_nod2D), rbuffer)
  call my_put_vara(ncid, zbar_n_bot_id, 1, nod2D, rbuffer)
  deallocate(rbuffer)

  ! element bottom depth (take into account partial cells if used)
  allocate(rbuffer(elem2D))
  call gather_elem(zbar_e_bot(1:myDim_elem2D), rbuffer)
  call my_put_vara(ncid, zbar_e_bot_id, 1, elem2D, rbuffer)
  deallocate(rbuffer)
  
  call my_close(ncid)
  
end subroutine write_mesh_info
!
!============================================================================
!
subroutine my_def_dim(ncid, short_name, value, id)
IMPLICIT NONE

integer,      intent(in)   :: ncid, value
character(*), intent(in)   :: short_name
integer,      intent(inout):: id
integer                    :: ierror, status

if (mype==0) then
   status =  nf_def_dim(ncid, trim(short_name), value, id)
end if

call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status)

end subroutine my_def_dim
!
!============================================================================
!
subroutine my_def_var(ncid, short_name, vtype, dsize, dids, id, att_text)
IMPLICIT NONE

integer, intent(in)        :: ncid, dsize, dids(dsize), vtype
character(*), intent(in)   :: short_name, att_text
integer,      intent(inout):: id
integer                    :: ierror, status

if (mype==0) then
   status = nf_def_var(ncid, trim(short_name), vtype, dsize, dids, id)
end if

call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status)

if (mype==0) then
   status = nf_put_att_text(ncid, id, 'long_name', len_trim(att_text), trim(att_text));
end if

call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status)

end subroutine my_def_var
!
!============================================================================
!
subroutine my_nf_enddef(ncid)
IMPLICIT NONE
integer, intent(in)        :: ncid
integer                    :: ierror, status

if (mype==0) then
   status = nf_enddef(ncid)
end if

call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status)
end subroutine my_nf_enddef
!
!============================================================================
!
subroutine my_put_vara_double_1D(ncid, varid, start, N, var)
IMPLICIT NONE

integer, intent(in)  :: ncid, varid, start, N
real(kind=WP)        :: var(:)
integer              :: ierror, status


  if (mype==0) status=nf_put_vara_double(ncid, varid, start, N, var)
  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status)

end subroutine my_put_vara_double_1D
!
!============================================================================
!
subroutine my_put_vara_double_2D(ncid, varid, start, N, var)
IMPLICIT NONE

integer, intent(in)  :: ncid, varid, start(:), N(:)
real(kind=WP)        :: var(:)
integer              :: ierror, status

  if (mype==0) status=nf_put_vara_double(ncid, varid, start, N, var)
  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status)

end subroutine my_put_vara_double_2D
!
!============================================================================
!
subroutine my_put_vara_int_1D(ncid, varid, start, N, var)
IMPLICIT NONE

integer, intent(in)  :: ncid, varid, start, N
integer              :: var(:)
integer              :: ierror, status


  if (mype==0) status=nf_put_vara_int(ncid, varid, start, N, var)
  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status)

end subroutine my_put_vara_int_1D
!
!============================================================================
!
subroutine my_put_vara_int_2D(ncid, varid, start, N, var)
IMPLICIT NONE

integer, intent(in)  :: ncid, varid, start(:), N(:)
integer              :: var(:)
integer              :: ierror, status

  if (mype==0) status=nf_put_vara_int(ncid, varid, start, N, var)
  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status)

end subroutine my_put_vara_int_2D
!
!============================================================================
!
subroutine my_create(filename, opt, ncid)
IMPLICIT NONE
integer, intent(in)        :: opt, ncid
character(*), intent(in)   :: filename
integer                    :: ierror, status
  if (mype==0) then  ! create a file
     ! create a file
     status = nf_create(filename, opt, ncid)
     if (status.ne.nf_noerr) call handle_err(status)
  end if
  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status)
end subroutine my_create
!
!============================================================================
!
subroutine my_close(ncid)
IMPLICIT NONE
integer, intent(in)        :: ncid
integer                    :: ierror, status

if (mype==0) status = nf_close(ncid)

call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status)
end subroutine my_close

end module io_mesh_info
