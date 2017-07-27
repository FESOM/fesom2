!
!-------------------------------------------------------------------------
!
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
  integer                    :: id_2, id_3
  integer                    :: i, k
  integer                    :: nod_area_id, elem_area_id
  integer                    :: elem_id
  character(100)             :: longname
  character(2000)            :: filename
  real(kind=WP), allocatable :: rbuffer(:)
  integer, allocatable       :: ibuffer(:)

  if (mype==0) then  ! create a file

     filename=trim(ResultPath)//runid//'.mesh.diag.nc'

     ! create a file
     status = nf_create(filename, nf_clobber, ncid)
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
    
     status = nf_def_var(ncid, 'nod_area', NF_DOUBLE, 2, (/nod_n_id, nl_id/), nod_area_id)
     if (status .ne. nf_noerr) call handle_err(status)
     
     status = nf_def_var(ncid, 'elem_area', NF_DOUBLE, 1, elem_n_id, elem_area_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_var(ncid, 'elem', NF_INT, 2, (/elem_n_id, id_3/), elem_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_enddef(ncid)
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

  allocate(ibuffer(elem2D))  
  ! elements
  do i=1, 3
     call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
     call gather_elem(elem2d_nodes(i, 1:myDim_elem2D), ibuffer)
     if (mype==0) then
        status=nf_put_vara_int(ncid, elem_id, (/1, i/), (/elem2D, 1/), ibuffer)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do
  deallocate(ibuffer)
  

  if (mype==0) then  
     status = nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
end subroutine write_mesh_diag

