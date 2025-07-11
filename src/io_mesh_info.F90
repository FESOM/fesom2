module io_mesh_info
USE MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
use g_config
use g_comm_auto
use o_PARAM

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
contains
!
! -------------------------------------------------------------------------
subroutine write_mesh_info(partit, mesh)
! this routine stores most of metadata used in FESOM. Shall be called at the cold start once during the simulation.
!
! .. info::
!
!    fesom.mesh.diag.nc is 77MB for the CORE II mesh with 47 vertical levels
implicit none
  type(t_mesh),   intent(in)   , target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer                    :: status, ncid, j
  integer                    :: nod_n_id, elem_n_id, edge_n_id, nod_part_id, elem_part_id
  integer                    :: nl_id, nl1_id
  integer                    :: id_2, id_3, id_4, id_N
  integer                    :: i, k, N_max
  integer                    :: nod_area_id, elem_area_id
  integer, target            :: z_id, zbar_id
  integer                    :: edges_id, edge_tri_id, edge_cross_dxdy_id
  integer                    :: nlevels_nod2D_id, nlevels_id, ulevels_nod2D_id, ulevels_id
  integer                    :: nod_in_elem2D_num_id, nod_in_elem2D_id
  integer                    :: gradient_sca_x_id, gradient_sca_y_id
  integer                    :: gradient_vec_x_id, gradient_vec_y_id
  integer                    :: zbar_e_bot_id, zbar_n_bot_id, zbar_e_srf_id, zbar_n_srf_id
  integer                    :: elem_id
  integer                    :: nod_id
  integer                    :: lon_id, lat_id
  ! UGRID NAMES START HERE
  integer                    :: face_node_id, edge_nodes_id, face_edges_id, face_links_id, edge_face_links_id
  character(100)             :: longname
  character(2000)            :: filename
  real(kind=WP), allocatable :: rbuffer(:), lrbuffer(:)
  integer, allocatable       :: ibuffer(:), lbuffer(:)
  character(2000)            :: att_text, short_name
  integer                    :: vtype
  integer, pointer           :: pid

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"


  call MPI_AllREDUCE(maxval(nod_in_elem2D_num), N_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_FESOM, MPIerr)

  filename=trim(ResultPath)//trim(runid)//'.mesh.diag.nc'
  call my_create(filename, IOR(NF_CLOBBER,IOR(NF_NETCDF4,NF_CLASSIC_MODEL)), ncid, partit)

  ! NOTE(PG): For naming conventions of node, edge, face, please see here:
  ! https://tinyurl.com/43z4m6cd
  !
  ! Define the dimensions
  call my_def_dim(ncid, 'nod2', nod2D,  nod_n_id,  partit)
  call my_def_dim(ncid, 'edg_n', edge2D, edge_n_id, partit)
  call my_def_dim(ncid, 'elem', elem2D, elem_n_id, partit)
  call my_def_dim(ncid, 'nz',    nl,     nl_id,     partit)
  call my_def_dim(ncid, 'nz1',   nl-1,   nl1_id,    partit)
  call my_def_dim(ncid, 'n2',    2,      id_2,      partit)
  call my_def_dim(ncid, 'n3',    3,      id_3,      partit)
  call my_def_dim(ncid, 'n4',    4,      id_4,      partit)
  call my_def_dim(ncid, 'N',     N_max,  id_N,      partit)

  ! Add ugrid information to the fesom.mesh.diag.nc file:
  call add_fesom_ugrid_info(ncid, partit)
  !Define the variables
  ! 1D
  call my_def_var(ncid,                                &  ! NetCDF Variable File Handle ID
    'nz',                                              &  ! Short Name
    NF_DOUBLE,                                         &  ! Variable Type
    1,                                                 &  ! Variable Dimensionality
    (/nl_id /),                                        &  ! Dimension IDs
    zbar_id,                                           &  ! ID
    'depth of levels',                                 &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit,                                            &  ! Partitioning
    "",                                                &  ! Standard Name
    "meters"                                           &  ! Units
  )
  ! PG: Add `positive` attribute according to CF-Conventions for depth:
  ! https://cfconventions.org/cf-conventions/cf-conventions.html#vertical-coordinate
  !
  ! Note that the sign (negative in normal FESOM code) is inverted during the saving process.
  if (partit%mype==0) then
    status = nf_put_att_text(ncid, zbar_id, 'positive', len_trim("down"), trim("down"));
  end if

  call my_def_var(ncid,                                &  ! NetCDF Variable File Handle ID
    'nz1',                                             &  ! Short Name
    NF_DOUBLE,                                         &  ! Variable Type
    1,                                                 &  ! Variable Dimensionality
    (/nl1_id/),                                        &  ! Dimension IDs
    z_id,                                              &  ! ID
    'depth of layers',                                 &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit                                             &  ! Partitioning
  )
  ! PG: Add `positive` attribute according to CF-Conventions for depth:
  ! https://cfconventions.org/cf-conventions/cf-conventions.html#vertical-coordinate
  !
  ! Note that the sign (negative in normal FESOM code) is inverted during the saving process.
  if (partit%mype==0) then
    status = nf_put_att_text(ncid, z_id, 'positive', len_trim("down"), trim("down"));
  end if

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'elem_area',                                              &  ! Short Name
    NF_DOUBLE,                                                &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/elem_n_id/),                                            &  ! Dimension IDs
    elem_area_id,                                             &  ! ID
    'element areas',                                          &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit                                                    &  ! Partitioning
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'nlevels_nod2D',                                          &  ! Short Name
    NF_INT,                                                   &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/nod_n_id/),                                             &  ! Dimension IDs
    nlevels_nod2D_id,                                         &  ! ID
    'number of levels below nodes',                           &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit                                                    &  ! Partitioning
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'nlevels',                                                &  ! Short Name
    NF_INT,                                                   &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/elem_n_id/),                                            &  ! Dimension IDs
    nlevels_id,                                               &  ! ID
    'number of levels below elements',                        &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit                                                    &  ! Partitioning
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'nod_in_elem2D_num',                                      &  ! Short Name
    NF_INT,                                                   &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/nod_n_id/),                                             &  ! Dimension IDs
    nod_in_elem2D_num_id,                                     &  ! ID
    'number of elements containing the node',                 &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit                                                    &  ! Partitioning
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'nod_part',                                               &  ! Short Name
    NF_INT,                                                   &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/nod_n_id/),                                             &  ! Dimension IDs
    nod_part_id,                                              &  ! ID
    'nodal partitioning at the cold start',                   &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit                                                    &  ! Partitioning
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'elem_part',                                              &  ! Short Name
    NF_INT,                                                   &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/elem_n_id/),                                            &  ! Dimension IDs
    elem_part_id,                                             &  ! ID
    'element partitioning at the cold start',                 &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit                                                    &  ! Partitioning
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'zbar_e_bottom',                                          &  ! Short Name
    NF_DOUBLE,                                                &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/elem_n_id/),                                            &  ! Dimension IDs
    zbar_e_bot_id,                                            &  ! ID
    'element bottom dep',                                     &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topologyth',
    partit                                                    &  ! Partitioning
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'zbar_n_bottom',                                          &  ! Short Name
    NF_DOUBLE,                                                &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/nod_n_id/),                                             &  ! Dimension IDs
    zbar_n_bot_id,                                            &  ! ID
    'nodal bottom depth',                                     &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology',
    partit                                                    &  ! Partitioning
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'lon',                                                    &  ! Short Name
    NF_DOUBLE,                                                &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/nod_n_id/),                                             &  ! Dimension IDs
    lon_id,                                                   &  ! ID
    'longitude',                                              &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit,                                                   &  ! Partitioning
    ! Everything after partit is an optional argument
    "longitude",                                              &
    "degrees_east"                                            &
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'lat',                                                    &  ! Short Name
    NF_DOUBLE,                                                &  ! Variable Type
    1,                                                        &  ! Variable Dimensionality
    (/nod_n_id/),                                             &  ! Dimension IDs
    lat_id,                                                   &  ! ID
    'latitude',                                               &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit,                                                   &  ! Partitioning
    ! Everything after partit is an optional argument
    "latitude",                                               &
    "degrees_north"                                           &
  )
  if (use_cavity) then
    call my_def_var(ncid, 'ulevels_nod2D',     NF_INT,    1, (/nod_n_id/),  ulevels_nod2D_id,     'number of levels above nodes',           partit)
    call my_def_var(ncid, 'ulevels',           NF_INT,    1, (/elem_n_id/), ulevels_id,           'number of levels above elements',        partit)
    call my_def_var(ncid, 'zbar_e_surface',    NF_DOUBLE, 1, (/elem_n_id/), zbar_e_srf_id,        'element surface depth',                  partit)
    call my_def_var(ncid, 'zbar_n_surface',    NF_DOUBLE, 1, (/nod_n_id/) , zbar_n_srf_id,        'nodal surface depth',                    partit)
  endif

  ! 2D
  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'nod_area',                                               &  ! Short Name
    NF_DOUBLE,                                                &  ! Variable Type
    2,                                                        &  ! Variable Dimensionality
    (/nod_n_id, nl_id/),                                      &  ! Dimension IDs
    nod_area_id,                                              &  ! ID
    'nodal areas',                                            &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit                                                    &  ! Partitioning
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
    'face_nodes',                                             &  ! Short Name
    NF_INT,                                                   &  ! Variable Type
    2,                                                        &  ! Variable Dimensionality
    (/elem_n_id, id_3/),                                      &  ! Dimension IDs
    face_node_id,                                             &  ! ID
    "Maps every triangular face to its three corner nodes.",  &  ! Long Name. NOTE(PG): Taken from https://ugrid-conventions.github.io/ugrid-conventions/#2d-triangular-mesh-topology
    partit,                                                   &  ! Partitioning
    ! Everything after partit is an optional argument
    "face",                                                   &  ! Standard Name
    "",                                                       &  ! Units
    "face_node_connectivity",                                 &  ! CF Role
    start_index = 1                                           &  ! Start index
  )
  if (partit%mype==0) then
    status = nf_put_att_text(ncid, face_node_id, 'location', len_trim("face"), trim("face"));
  end if

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
  "edge_nodes",                                               &  ! Short Name
  NF_INT,                                                     &  ! Variable Type
  2,                                                          &  ! Variable Dimensionality
  (/edge_n_id, id_2/),                                        &  ! Dimension IDs
  edge_nodes_id,                                              &  ! ID
  "Maps every edge to the two nodes that it connects",        &  ! Long Name
  partit,                                                     &  ! Partitioning
  ! Everything after partit is an optional argument
  "edge",                                                     &  ! Standard Name
  "",                                                         &  ! Units
  "edge_node_connectivity",                                   &  ! CF Role
  start_index = 1                                             &  ! Start index
  )

  call my_def_var(ncid,                                       &  ! NetCDF Variable File Handle ID
  'face_edges',                                               &  ! Short Name
  NF_INT,                                                     &  ! Variable Type
  2,                                                          &  ! Variable Dimensionality
  (/elem_n_id, id_3/),                                        &  ! Dimension IDs
  face_edges_id,                                              &  ! ID
  'Maps every triangular face to its three edges.',           &  ! Long Name
  partit,                                                     &  ! Partitioning
  ! Everything after partit is an optional argument
  "face_edges",                                               &  ! Standard Name
  "",                                                         &  ! Units
  "face_edge_connectivity"                                    &  ! CF Role
  )

call my_def_var(ncid,                                         &  ! NetCDF Variable File Handle ID
  "face_links",                                               &  ! Short Name
  NF_INT,                                                     &  ! Variable Type
  2,                                                          &  ! Variable Dimensionality
  (/elem_n_id, id_3/),                                        &  ! Dimension IDs
  face_links_id,                                              &  ! ID
  "neighbor faces for faces",                                 &  ! Long Name
  partit,                                                     &  ! Partitioning
  ! Everything after partit is an optional argument
  "face_links",                                               &  ! Standard Name
  "",                                                         &  ! Units
  "face_face_connectivity",                                   &  ! CF Role
  "missing neighbor faces are indicated using _FillValue",    &  ! Comment
  -999                                                        &  ! Missing Value
  )

call my_def_var(ncid,                                         &  ! NetCDF Variable File Handle ID
  "edge_face_links",                                          &  ! Short Name
  NF_INT,                                                     &  ! Variable Type
  2,                                                          &  ! Variable Dimensionality
  (/edge_n_id,  id_2/),                                       &  ! Dimension IDs
  edge_face_links_id,                                         &  ! ID
  "neighbor faces for edges",                                 &  ! Long Name
  partit,                                                     &  ! Partitioning
  "edge_face_links",                                          &  ! Standard Name
  "",                                                         &  ! Units
  "edge_face_connectivity",                                   &  ! CF Role
  "missing neighbor faces are indicated using _FillValue",    &  ! Comment
  -999,                                                       &  ! Missing Value
  1                                                           &  ! Start index
  )

  call my_def_var(ncid, 'nod_in_elem2D',     NF_INT,    2, (/nod_n_id, id_N/),  nod_in_elem2D_id,   'elements containing the node', partit)
  call my_def_var(ncid, 'edge_cross_dxdy',   NF_DOUBLE, 2, (/edge_n_id, id_4/), edge_cross_dxdy_id, 'edge cross distancess',        partit)
  call my_def_var(ncid, 'gradient_sca_x',    NF_DOUBLE, 2, (/elem_n_id, id_3/), gradient_sca_x_id,  'x component of a gradient at nodes of an element', partit)
  call my_def_var(ncid, 'gradient_sca_y',    NF_DOUBLE, 2, (/elem_n_id, id_3/), gradient_sca_y_id,  'y component of a gradient at nodes of an element', partit)
  call my_def_var(ncid, 'gradient_vec_x',    NF_DOUBLE, 2, (/elem_n_id, id_3/), gradient_vec_x_id,  'x component of a gradient at elements of an element', partit)
  call my_def_var(ncid, 'gradient_vec_y',    NF_DOUBLE, 2, (/elem_n_id, id_3/), gradient_vec_y_id,  'y component of a gradient at elements of an element', partit)
  call my_nf_enddef(ncid, partit)

  ! NOTE(PG): Same order as definition!
  ! NOTE(PG): We invert the sign of ``zbar`` here to conform with CF-Conventions
  ! https://cfconventions.org/cf-conventions/cf-conventions.html#vertical-coordinate
  ! Mid and full levels
  call my_put_vara(ncid, zbar_id, 1, nl, -zbar, partit)
  call my_put_vara(ncid, z_id, 1, nl-1, -Z, partit)

  ! nodal areas
  allocate(rbuffer(nod2D))
  do k=1, nl
     ! area     ...area of prism top face
     ! areasvol ...area of scalar volume cell
     ! better write areasvol to meshdiag, without cavity area==areasvol but with
     ! cavity area!=areasvol, at least where is cavity
     !!PS      call gather_nod(area(k, :), rbuffer, partit)
     call gather_nod(areasvol(k, :), rbuffer, partit)
     call my_put_vara(ncid, nod_area_id, (/1, k/), (/nod2D, 1/), rbuffer, partit)
  end do
  deallocate(rbuffer)

  ! WRITE(*,*) "elem_area"
  allocate(rbuffer(elem2D))
  call gather_elem(elem_area(1:myDim_elem2D), rbuffer, partit)
  call my_put_vara(ncid, elem_area_id, 1, elem2D, rbuffer, partit)
  deallocate(rbuffer)

  ! WRITE(*,*) "nlevels_nod2D"
  allocate(ibuffer(nod2D))
  call gather_nod(nlevels_nod2D(1:myDim_nod2D), ibuffer, partit)
  call my_put_vara(ncid, nlevels_nod2D_id, 1, nod2D, ibuffer, partit)
  deallocate(ibuffer)

  ! WRITE(*,*) "nlevels"
  allocate(ibuffer(elem2D))
  call gather_elem(nlevels(1:myDim_elem2D), ibuffer, partit)
  call my_put_vara(ncid, nlevels_id, 1, elem2D, ibuffer, partit)
  deallocate(ibuffer)

  ! WRITE(*,*) "nod_in_elem2D_num"
  allocate(ibuffer(nod2D))
  call gather_nod(nod_in_elem2D_num(1:myDim_nod2D), ibuffer, partit)
  call my_put_vara(ncid, nod_in_elem2D_num_id, 1, nod2D, ibuffer, partit)
  deallocate(ibuffer)

  ! WRITE(*,*) "nod_part"
  allocate(ibuffer(nod2D))
  allocate(lbuffer(myDim_nod2D))
  lbuffer=mype
  call gather_nod(lbuffer, ibuffer, partit)
  call my_put_vara(ncid, nod_part_id, 1, nod2D, ibuffer, partit)
  deallocate(lbuffer, ibuffer)

  ! WRITE(*,*) "elem_part"
  allocate(ibuffer(elem2D))
  allocate(lbuffer(myDim_elem2D))
  lbuffer=mype
  call gather_elem(lbuffer, ibuffer, partit)
  call my_put_vara(ncid, elem_part_id, 1, elem2D, ibuffer, partit)
  deallocate(lbuffer, ibuffer)

  ! WRITE(*,*) "zbar_e_bottom"
  allocate(rbuffer(elem2D))
  call gather_elem(zbar_e_bot(1:myDim_elem2D), rbuffer, partit)
  call my_put_vara(ncid, zbar_e_bot_id, 1, elem2D, rbuffer, partit)
  deallocate(rbuffer)

  ! WRITE(*,*) "zbar_n_bottom"
  allocate(rbuffer(nod2D))
  call gather_nod(zbar_n_bot(1:myDim_nod2D), rbuffer, partit)
  call my_put_vara(ncid, zbar_n_bot_id, 1, nod2D, rbuffer, partit)
  deallocate(rbuffer)

  ! WRITE(*,*) "lon and lat"
  ! nodes (GEO coordinates)
  allocate(rbuffer(nod2D))
  do i=1, 2
     call gather_nod(geo_coord_nod2D(i, 1:myDim_nod2D), rbuffer, partit)
     rbuffer = rbuffer/rad
     if (i == 1) then
        call my_put_vara(ncid, lon_id, 1, nod2D, rbuffer, partit)
     else
        call my_put_vara(ncid, lat_id, 1, nod2D, rbuffer, partit)
     endif
  end do
  deallocate(rbuffer)

  ! WRITE(*,*) "nod_area"
  allocate(rbuffer(nod2D))
  do k=1, nl
     call gather_nod(area(k, :), rbuffer, partit)
     call my_put_vara(ncid, nod_area_id, (/1, k/), (/nod2D, 1/), rbuffer, partit)
  end do
  deallocate(rbuffer)

  ! WRITE(*,*) "face_nodes"
  ! used to be called: elements
  allocate(ibuffer(elem2D))
  allocate(lbuffer(myDim_elem2D))
  do i=1, 3
     do k=1, myDim_elem2D
        lbuffer(k)=myList_nod2D(elem2d_nodes(i, k))
     end do
     call gather_elem(lbuffer, ibuffer, partit)
     call my_put_vara(ncid, face_node_id, (/1, i/), (/elem2D, 1/), ibuffer, partit)
  end do
  deallocate(lbuffer, ibuffer)

  ! WRITE(*,*) "edge_nodes"
  ! This used to be called "edges" before UGRID Change
  allocate(ibuffer(edge2D))
  allocate(lbuffer(myDim_edge2D))
  do i=1, 2
     do k=1, myDim_edge2D
        lbuffer(k)=myList_nod2D(edges(i, k))
     end do
     call gather_edge(lbuffer, ibuffer, partit)
     call my_put_vara(ncid, edge_nodes_id, (/1, i/), (/edge2D, 1/), ibuffer, partit)
  end do
  deallocate(lbuffer, ibuffer)

  ! WRITE(*,*) "face_edges"
  ! This used to be called elem_edge
  allocate(ibuffer(elem2D))
  do i=1, 3
     call gather_elem(elem_edges(i,1:myDim_elem2D), ibuffer, partit)
     call my_put_vara(ncid, face_edges_id, (/1, i/), (/elem2D, 1/), ibuffer, partit)
  end do
  deallocate(ibuffer)

  ! WRITE(*,*) "face_links"
  ! This used to be called elem_neighbors
  allocate(lbuffer(myDim_elem2D), ibuffer(elem2D))
  do i=1, 3
     lbuffer(1:myDim_elem2D) = elem_neighbors(i, 1:myDim_elem2D)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
     do k=1, myDim_elem2D
        if (elem_neighbors(i, k) > 0) then
          lbuffer(k) = elem_neighbors(i, k)
        else
          lbuffer(k) = -999
        end if
     end do
!$OMP END PARALLEL DO
     call gather_elem(lbuffer, ibuffer, partit)
     call my_put_vara(ncid, face_links_id, (/1, i/), (/elem2D, 1/), ibuffer, partit)
  end do
  deallocate(lbuffer, ibuffer)

  ! WRITE(*,*) "edge_face_links"
  ! Used to be called edge_tri
  allocate(ibuffer(edge2D))
  allocate(lbuffer(myDim_edge2D))
  do i=1, 2
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
     do k=1, myDim_edge2D
        if (edge_tri(i,k) > 0) then
           lbuffer(k) = myList_elem2D(edge_tri(i,k))
        else
           lbuffer(k) = -999  ! Changed from missing value of 0 before
        endif
     end do
!$OMP END PARALLEL DO
     call gather_edge(lbuffer, ibuffer, partit)
     call my_put_vara(ncid, edge_face_links_id, (/1, i/), (/edge2D, 1/), ibuffer, partit)
  end do
  deallocate(lbuffer, ibuffer)

  ! WRITE(*,*) "nod_in_elem2D"
  ! elements containing the node
  allocate(ibuffer(nod2D))
  allocate(lbuffer(myDim_nod2D))
  DO i=1, N_max
     lbuffer=0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
     do k=1, myDim_nod2D
        if ((nod_in_elem2D_num(k)>=i)) then
           lbuffer(k)=myList_elem2D(nod_in_elem2D(i, k))
        end if
     end do
!$OMP END PARALLEL DO
     call gather_nod(lbuffer, ibuffer, partit)
     call my_put_vara(ncid, nod_in_elem2D_id, (/1, i/), (/nod2D, 1/), ibuffer, partit)
  END DO
  deallocate(lbuffer, ibuffer)


  ! WRITE(*,*) "edge_cross_dxdy"
  allocate(rbuffer(edge2D))
  allocate(lrbuffer(myDim_edge2D))
  do i=1, 4
     lrbuffer=edge_cross_dxdy(i, 1:myDim_edge2D)
     call gather_edge(lrbuffer, rbuffer, partit)
     call my_put_vara(ncid, edge_cross_dxdy_id, (/1, i/), (/edge2D, 1/), rbuffer, partit)
  end do
  deallocate(rbuffer, lrbuffer)

  ! WRITE(*,*) "gradient_sca_x"
    ! X component of gadient at elements
  allocate(rbuffer(elem2D))
  do i=1, 3
     call gather_elem(gradient_sca(i, 1:myDim_elem2D), rbuffer, partit)
     call my_put_vara(ncid, gradient_sca_x_id, (/1, 4-i/), (/elem2D, 1/), rbuffer, partit) ! (4-i), NETCDF will permute otherwise
  end do
  deallocate(rbuffer)

  ! WRITE(*,*) "gradient_sca_y"
    ! Y component of gadient at elements
  allocate(rbuffer(elem2D))
  do i=1, 3
     call gather_elem(gradient_sca(i+3, 1:myDim_elem2D), rbuffer, partit)
     call my_put_vara(ncid, gradient_sca_y_id, (/1, 4-i/), (/elem2D, 1/), rbuffer, partit)! (4-i), NETCDF will permute otherwise
  end do
  deallocate(rbuffer)

  ! WRITE(*,*) "gradient_vec_x"
    ! X component of gadient at elements
  allocate(rbuffer(elem2D))
  do i=1, 3
     call gather_elem(gradient_vec(i, 1:myDim_elem2D), rbuffer, partit)
     call my_put_vara(ncid, gradient_vec_x_id, (/1, 4-i/), (/elem2D, 1/), rbuffer, partit) ! (4-i), NETCDF will permute otherwise
  end do
  deallocate(rbuffer)

  ! WRITE(*,*) "gradient_vec_y"
    ! Y component of gadient at elements
  allocate(rbuffer(elem2D))
  do i=1, 3
     call gather_elem(gradient_vec(i+3, 1:myDim_elem2D), rbuffer, partit)
     call my_put_vara(ncid, gradient_vec_y_id, (/1, 4-i/), (/elem2D, 1/), rbuffer, partit)! (4-i), NETCDF will permute otherwise
  end do
  deallocate(rbuffer)

  ! element bottom depth (take into account partial cells if used)
  allocate(rbuffer(elem2D))
  call gather_elem(zbar_e_bot(1:myDim_elem2D), rbuffer, partit)
  call my_put_vara(ncid, zbar_e_bot_id, 1, elem2D, rbuffer, partit)
  deallocate(rbuffer)

  if (use_cavity) then
    ! number of levels above elements
    allocate(ibuffer(elem2D))
    call gather_elem(ulevels(1:myDim_elem2D), ibuffer, partit)
    call my_put_vara(ncid, ulevels_id, 1, elem2D, ibuffer, partit)
    deallocate(ibuffer)

    ! number of levels above nodes
    allocate(ibuffer(nod2D))
    call gather_nod(ulevels_nod2D(1:myDim_nod2D), ibuffer, partit)
    call my_put_vara(ncid, ulevels_nod2D_id, 1, nod2D, ibuffer, partit)
    deallocate(ibuffer)

    ! nodal surface depth (take into account partial cells if used)
    allocate(rbuffer(nod2D))
    call gather_nod(zbar_n_srf(1:myDim_nod2D), rbuffer, partit)
    call my_put_vara(ncid, zbar_n_srf_id, 1, nod2D, rbuffer, partit)
    deallocate(rbuffer)

    ! element surface depth (take into account partial cells if used)
    allocate(rbuffer(elem2D))
    call gather_elem(zbar_e_srf(1:myDim_elem2D), rbuffer, partit)
    call my_put_vara(ncid, zbar_e_srf_id, 1, elem2D, rbuffer, partit)
    deallocate(rbuffer)
  endif

  call my_close(ncid, partit)

end subroutine write_mesh_info
!
!============================================================================
!
subroutine my_def_dim(ncid, short_name, value, id, partit)
IMPLICIT NONE

type(t_partit), intent(inout) :: partit
integer,        intent(in)    :: ncid, value
character(*),   intent(in)    :: short_name
integer,        intent(inout) :: id
integer                       :: ierror, status

if (partit%mype==0) then
   status =  nf_def_dim(ncid, trim(short_name), value, id)
end if

call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

end subroutine my_def_dim
!
!============================================================================
!
!
subroutine add_fesom_ugrid_info(ncid, partit)
! Writes UGRID information to the mesh diag files to describe the mesh topology.
!
!
IMPLICIT NONE
type(t_partit), intent(inout) :: partit
integer, intent(in)           :: ncid
integer                       :: info_id
integer                       :: fesom_mesh_id
integer                       :: ierror, status

! NOTE(PG): Adds global Conventions attribute
if (partit%mype==0) then
   status = nf_put_att_text(ncid, NF_GLOBAL, 'Conventions', len_trim('UGRID-1.0'), 'UGRID-1.0')
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)



! NOTE(PG): Order of calls should mirror the UGRID website; otherwise I cannot keep track:

call my_def_dim(ncid, "info", 1, info_id, partit)
if (partit%mype==0) then
   status = nf_def_var(ncid, 'fesom_mesh', NF_INT, 0, info_id, fesom_mesh_id)
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'cf_role', len_trim('mesh_topology'), trim('mesh_topology'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'long_name', len_trim('Topology data of 2D unstructured mesh'), trim('Topology data of 2D unstructured mesh'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_int(ncid, fesom_mesh_id, 'topology_dimension', NF_INT, 1, 2);
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'node_coordinates', len_trim('lon lat'), trim('lon lat'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'face_node_connectivity', len_trim('face_nodes'), trim('face_nodes'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'face_dimension', len_trim('elem') , trim('elem'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'edge_node_connectivity', len_trim('edge_nodes'), trim('edge_nodes'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'edge_dimension', len_trim('edg_n'), trim('edg_n'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

!FIXME(PG): edge coordinates missing
!Mesh2:edge_coordinates = "Mesh2_edge_x Mesh2_edge_y" ; // optional attribute (requires edge_node_connectivity)
!FIXME(PG): face coordinates missing
!Mesh2:face_coordinates = "Mesh2_face_x Mesh2_face_y" ; // optional attribute
if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'face_edge_connectivity', len_trim('face_edges'), trim('face_edges'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'face_face_connectivity', len_trim('face_links'), trim('face_links'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

if (partit%mype==0) then
  status = nf_put_att_text(ncid, fesom_mesh_id, 'edge_face_connectivity', len_trim('edge_face_links'), trim('edge_face_links'));
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

end subroutine add_fesom_ugrid_info
!============================================================================
!
!
!============================================================================
!
! NOTE/IDEA(CR, PG): This could be an INTERFACE subroutine (overloading) to give either missing value as int or real.
subroutine my_def_var(ncid, short_name, vtype, dsize, dids, id, att_text, partit, standard_name, units, cf_role, comment, missing_value, start_index)
IMPLICIT NONE

type(t_partit), intent(inout):: partit
integer,        intent(in)   :: ncid, dsize, dids(dsize), vtype
character(*),   intent(in)   :: short_name, att_text
integer,        intent(inout):: id
integer                      :: ierror, status
character(*), intent(in), optional :: standard_name, units, cf_role
character(*), intent(in), optional :: comment
integer,      intent(in), optional :: missing_value
integer,      intent(in), optional :: start_index


if (partit%mype==0) then
   status = nf_def_var(ncid, trim(short_name), vtype, dsize, dids, id)
end if

if (partit%mype==0) then
   status = nf_put_att_text(ncid, id, 'long_name', len_trim(att_text), trim(att_text));
end if

if (partit%mype==0) then
  if (present(standard_name)) then
      status = nf_put_att_text(ncid, id, 'standard_name', len_trim(standard_name), trim(standard_name))
    endif
endif

if (partit%mype==0) then
  if (present(units)) then
     status = nf_put_att_text(ncid, id, 'units', len_trim(units), trim(units));
  endif
endif

if (partit%mype==0) then
  if (present(cf_role)) then
     status = nf_put_att_text(ncid, id, 'cf_role', len_trim(cf_role), trim(cf_role));
  endif
endif

if (partit%mype==0) then
  if (present(comment)) then
      status = nf_put_att_text(ncid, id, 'comment', len_trim(comment), trim(comment));
  endif
endif


if (partit%mype==0) then
  if (present(missing_value)) then
      status = nf_put_att_int(ncid, id, '_FillValue', NF_INT, 1, missing_value);
  endif
endif

if (partit%mype==0) then
  if (present(start_index)) then
    status = nf_put_att_int(ncid, id, 'start_index', NF_INT, 1, start_index);
  endif
end if
call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)

call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)





end subroutine my_def_var
!
!============================================================================
!
subroutine my_nf_enddef(ncid, partit)
IMPLICIT NONE
type(t_partit), intent(inout) :: partit
integer, intent(in)           :: ncid
integer                       :: ierror, status

if (partit%mype==0) then
   status = nf_enddef(ncid)
end if

call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)
end subroutine my_nf_enddef
!
!============================================================================
!
subroutine my_put_vara_double_1D(ncid, varid, start, N, var, partit)
IMPLICIT NONE
type(t_partit), intent(inout) :: partit
integer, intent(in)           :: ncid, varid, start, N
real(kind=WP)                 :: var(:)
integer                       :: ierror, status

  if (partit%mype==0) status=nf_put_vara_double(ncid, varid, start, N, var)
  call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status, partit)

end subroutine my_put_vara_double_1D
!
!============================================================================
!
subroutine my_put_vara_double_2D(ncid, varid, start, N, var, partit)
IMPLICIT NONE
type(t_partit), intent(inout) :: partit
integer, intent(in)           :: ncid, varid, start(:), N(:)
real(kind=WP)                 :: var(:)
integer                       :: ierror, status

  if (partit%mype==0) status=nf_put_vara_double(ncid, varid, start, N, var)
  call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status, partit)

end subroutine my_put_vara_double_2D
!
!============================================================================
!
subroutine my_put_vara_int_1D(ncid, varid, start, N, var, partit)
IMPLICIT NONE
type(t_partit), intent(inout) :: partit
integer,        intent(in)    :: ncid, varid, start, N
integer                       :: var(:)
integer                       :: ierror, status

  if (partit%mype==0) status=nf_put_vara_int(ncid, varid, start, N, var)
  call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status, partit)

end subroutine my_put_vara_int_1D
!
!============================================================================
!
subroutine my_put_vara_int_2D(ncid, varid, start, N, var, partit)
IMPLICIT NONE

type(t_partit), intent(inout) :: partit
integer,        intent(in)    :: ncid, varid, start(:), N(:)
integer                       :: var(:)
integer                       :: ierror, status

  if (partit%mype==0) status=nf_put_vara_int(ncid, varid, start, N, var)
  call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status, partit)

end subroutine my_put_vara_int_2D
!
!============================================================================
!
subroutine my_create(filename, opt, ncid, partit)
IMPLICIT NONE
type(t_partit), intent(inout):: partit
integer,        intent(in)   :: opt, ncid
character(*),   intent(in)   :: filename
integer                      :: ierror, status
  if (partit%mype==0) then  ! create a file
     ! create a file
     status = nf_create(filename, opt, ncid)
     if (status.ne.nf_noerr) call handle_err(status, partit)
  end if
  call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
  if (status .ne. nf_noerr) call handle_err(status, partit)
end subroutine my_create
!
!============================================================================
!
subroutine my_close(ncid, partit)
IMPLICIT NONE
type(t_partit), intent(inout) :: partit
integer,        intent(in)    :: ncid
integer                       :: ierror, status

if (partit%mype==0) status = nf_close(ncid)

call MPI_BCast(status, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
if (status .ne. nf_noerr) call handle_err(status, partit)
end subroutine my_close
end module io_mesh_info
