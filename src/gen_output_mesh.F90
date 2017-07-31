!
!-------------------------------------------------------------------------
! this routine stores most of metadata used in FESOM. Shall be called at the cold start once during the simulation. 
! info: fesom.mesh.diag.nc is 77MB for the CORE II mesh with 47 vertical levels
subroutine write_mesh_diag

use g_PARSUP
use o_MESH
use g_config
use g_comm_auto
use o_ARRAYS
implicit none

#include "netcdf.inc" 

  integer                    :: status, ncid, j
  integer                    :: nod_n_id, elem_n_id, edge_n_id
  integer                    :: nl_id, nl1_id
  integer                    :: id_2, id_3, id_4, id_N
  integer                    :: i, k, N_max
  integer                    :: nod_area_id, elem_area_id
  integer, target            :: z_id, zbar_id
  integer                    :: edges_id, edge_tri_id, edge_cross_dxdy_id
  integer                    :: nlevels_nod2D_id, nlevels_id
  integer                    :: nod_in_elem2D_num_id, nod_in_elem2D_id
  integer                    :: gradient_sca_x_id, gradient_sca_y_id
  integer                    :: elem_id
  integer                    :: nod_id
  character(100)             :: longname
  character(2000)            :: filename
  real(kind=WP), allocatable :: rbuffer(:)
  integer, allocatable       :: ibuffer(:), lbuffer(:)
  character(2000)            :: att_text, short_name
  integer                    :: vtype
  integer, pointer           :: pid

  
  call MPI_AllREDUCE(maxval(nod_in_elem2D_num), N_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, MPIerr)


  if (mype==0) then  ! create a file

     filename=trim(ResultPath)//runid//'.mesh.diag.nc'

     ! create a file
     status = nf_create(filename, IOR(NF_CLOBBER,IOR(NF_NETCDF4,NF_CLASSIC_MODEL)), ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nod_n', nod2d, nod_n_id)
     if (status .ne. nf_noerr) call handle_err(status)
 
     status = nf_def_dim(ncid, 'edg_n', edge2d, edge_n_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_dim(ncid, 'elem_n', elem2d, elem_n_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_dim(ncid, 'nl',  nl,   nl_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_dim(ncid, 'nl1', nl-1, nl1_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_dim(ncid, 'n2',  2,   id_2)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_dim(ncid, 'n3',  3,   id_3)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_dim(ncid, 'n4',  4,   id_4)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_dim(ncid, 'N',   N_max, id_N)
     if (status .ne. nf_noerr) call handle_err(status)

     pid=>zbar_id; short_name='zbar'; att_text='depth of levels'; vtype=NF_DOUBLE

     status = nf_def_var(ncid, 'Z', NF_DOUBLE, 1, nl1_id, z_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='depth of layers'
     status = nf_put_att_text(ncid, z_id, 'long_name', len_trim(att_text), trim(att_text));
    
     status = nf_def_var(ncid, 'nod_area', NF_DOUBLE, 2, (/nod_n_id, nl_id/), nod_area_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='nodal areas'
     status = nf_put_att_text(ncid, nod_area_id, 'long_name', len_trim(att_text), trim(att_text));
     
     status = nf_def_var(ncid, 'elem_area', NF_DOUBLE, 1, elem_n_id, elem_area_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='element areas'
     status = nf_put_att_text(ncid, elem_area_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'elem', NF_INT, 2, (/elem_n_id, id_3/), elem_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='elements'
     status = nf_put_att_text(ncid, elem_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'nodes', NF_DOUBLE, 2, (/nod_n_id, id_2/), nod_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='nodal geo. coordinates'
     status = nf_put_att_text(ncid, nod_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'nlevels_nod2D', NF_INT, 1, nod_n_id, nlevels_nod2D_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='number of levels below nodes'
     status = nf_put_att_text(ncid, nlevels_nod2D_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'nlevels', NF_INT, 1, elem_n_id, nlevels_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='number of levels below elements'
     status = nf_put_att_text(ncid, nlevels_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'nod_in_elem2D_num', NF_INT, 1, nod_n_id, nod_in_elem2D_num_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='number of elements containing the node'
     status = nf_put_att_text(ncid, nod_in_elem2D_num_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'nod_in_elem2D', NF_INT, 2, (/nod_n_id, id_N/), nod_in_elem2D_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='elements containing the node'
     status = nf_put_att_text(ncid, nod_in_elem2D_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'edges', NF_INT, 2, (/edge_n_id, id_2/), edges_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='edges'
     status = nf_put_att_text(ncid, edges_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'edge_tri', NF_INT, 2, (/edge_n_id, id_2/), edge_tri_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='edge triangles'
     status = nf_put_att_text(ncid, edge_tri_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'edge_cross_dxdy', NF_DOUBLE, 2, (/edge_n_id, id_4/), edge_cross_dxdy_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='edge triangles'
     status = nf_put_att_text(ncid, edge_cross_dxdy_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'gradient_sca_x', NF_DOUBLE, 2, (/elem_n_id, id_3/), gradient_sca_x_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='x component of a gradient at nodes of an element'
     status = nf_put_att_text(ncid, gradient_sca_x_id, 'long_name', len_trim(att_text), trim(att_text));

     status = nf_def_var(ncid, 'gradient_sca_y', NF_DOUBLE, 2, (/elem_n_id, id_3/), gradient_sca_y_id)
     if (status .ne. nf_noerr) call handle_err(status)
     att_text='y component of a gradient at nodes of an element'
     status = nf_put_att_text(ncid, gradient_sca_y_id, 'long_name', len_trim(att_text), trim(att_text));


     status = nf_enddef(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  ! vercical levels/layers
  if (mype==0) then
     status=nf_put_vara_double(ncid, zbar_id, 1, nl, zbar)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, z_id, 1, nl-1, Z)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  ! nodal areas
  allocate(rbuffer(nod2D))
  do k=1, nl
     call gather_nod(area(k, :), rbuffer)
     if (mype==0) then
        status=nf_put_vara_double(ncid, nod_area_id, (/1, k/), (/nod2D, 1/), rbuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do
  deallocate(rbuffer)

  ! element areas
  allocate(rbuffer(elem2D))
  call gather_elem(elem_area(1:myDim_elem2D), rbuffer)
  if (mype==0) then
     status=nf_put_vara_double(ncid, elem_area_id, 1, elem2D, rbuffer)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  deallocate(rbuffer)

  ! elements
  allocate(ibuffer(elem2D))
  allocate(lbuffer(myDim_elem2D))  
  do i=1, 3
     do k=1, myDim_elem2D
        lbuffer(k)=myList_nod2D(elem2d_nodes(i, k))
     end do
     call gather_elem(lbuffer, ibuffer)
     if (mype==0) then
        status=nf_put_vara_int(ncid, elem_id, (/1, i/), (/elem2D, 1/), ibuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do
  deallocate(lbuffer, ibuffer)

  ! number of levels below elements
  allocate(ibuffer(elem2D))
  call gather_elem(nlevels(1:myDim_elem2D), ibuffer)
  if (mype==0) then
     status=nf_put_vara_int(ncid, nlevels_id, 1, elem2D, ibuffer)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  deallocate(ibuffer)

  ! number of levels below nodes
  allocate(ibuffer(nod2D))
  call gather_nod(nlevels_nod2D(1:myDim_nod2D), ibuffer)
  if (mype==0) then
     status=nf_put_vara_int(ncid, nlevels_nod2D_id, 1, nod2D, ibuffer)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  deallocate(ibuffer)
  
  ! number of elements containing the node
  allocate(ibuffer(nod2D))
  call gather_nod(nod_in_elem2D_num(1:myDim_nod2D), ibuffer)
  if (mype==0) then
     status=nf_put_vara_int(ncid, nod_in_elem2D_num_id, 1, nod2D, ibuffer)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  deallocate(ibuffer)

  ! elements containing the node
  allocate(ibuffer(nod2D))
  allocate(lbuffer(myDim_nod2D))
  DO i=1, N_max
     lbuffer=0
        do k=1, myDim_nod2D
           if ((nod_in_elem2D(i, k) > 0) .and. (N_max<=nod_in_elem2D_num(k))) then
              lbuffer(k)=myList_elem2D(nod_in_elem2D(i, k))
           end if
        end do
     call gather_nod(lbuffer, ibuffer)
     if (mype==0) then
        status=nf_put_vara_int(ncid, nod_in_elem2D_id, (/1, i/), (/nod2D, 1/), ibuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  END DO
  deallocate(lbuffer, ibuffer)

  ! nodes (GEO coordinates)
  allocate(rbuffer(nod2D))
  do i=1, 2
     call gather_nod(geo_coord_nod2D(i, 1:myDim_nod2D), rbuffer)
     if (mype==0) then
        status=nf_put_vara_double(ncid, nod_id, (/1, i/), (/nod2D, 1/), rbuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do
  deallocate(rbuffer)

  ! edges
  allocate(ibuffer(edge2D))
  allocate(lbuffer(myDim_edge2D))
  do i=1, 2
     do k=1, myDim_edge2D
        lbuffer(k)=myList_edge2D(edges(i, k))
     end do
     call gather_edge(lbuffer, ibuffer)
     if (mype==0) then
        status=nf_put_vara_int(ncid, edges_id, (/1, i/), (/edge2D, 1/), ibuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do
  deallocate(lbuffer, ibuffer)

  ! edge triangles
  allocate(ibuffer(edge2D))
  allocate(lbuffer(myDim_edge2D))
  do i=1, 2
     do k=1, myDim_edge2D
        lbuffer(k)=myList_elem2D(edge_tri(i, k))
     end do
     call gather_edge(lbuffer, ibuffer)
     if (mype==0) then
        status=nf_put_vara_int(ncid, edge_tri_id, (/1, i/), (/edge2D, 1/), ibuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do
  deallocate(lbuffer, ibuffer)

  ! edge cross distances
  allocate(rbuffer(edge2D))
  do i=1, 4
     call gather_edge(edge_cross_dxdy(i, :), rbuffer)
     if (mype==0) then
        status=nf_put_vara_double(ncid, edge_cross_dxdy_id, (/1, i/), (/edge2D, 1/), rbuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do
  deallocate(rbuffer)


  ! X component of gadient at elements
  allocate(rbuffer(elem2D))
  do i=1, 3
     call gather_elem(gradient_sca(i, 1:myDim_elem2D), rbuffer)
     if (mype==0) then
        status=nf_put_vara_double(ncid, gradient_sca_x_id, (/1, i/), (/elem2D, 1/), rbuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do
  deallocate(rbuffer)

  ! Y component of gadient at elements
  allocate(rbuffer(elem2D))
  do i=1, 3
     call gather_elem(gradient_sca(i+3, 1:myDim_elem2D), rbuffer)
     if (mype==0) then
        status=nf_put_vara_double(ncid, gradient_sca_y_id, (/1, i/), (/elem2D, 1/), rbuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do
  deallocate(rbuffer)

  if (mype==0) then  
     status = nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
end subroutine write_mesh_diag
!
!============================================================================
!
subroutine my_def_var(ncid, pid, short_name, att_text, vtype)
use g_PARSUP
USE o_MESH
IMPLICIT NONE

integer, intent(in)        :: ncid, pid
character(:), intent(in)   :: short_name, att_text
integer                    :: status, lstatus

if (mype==0) then
   status = nf_def_var(ncid, trim(short_name), vtype, 1, nl_id, pid)
   if (status .ne. nf_noerr) call handle_err(status)
   lstatus = nf_put_att_text(ncid, pid, 'long_name', len_trim(att_text), trim(att_text));
end if



end subroutine my_def_var
