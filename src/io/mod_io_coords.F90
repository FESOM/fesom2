module mod_io_coords
    ! Shared ushow/xarray coordinate-embedding helpers for FESOM3 Zarr stores (restart Stage 2,
    ! Task 2.1). Factored VERBATIM out of mod_io_means so the restart checkpoint stores
    ! (mod_io_restart) are byte-identical in their lon/lat coordinate + _ARRAY_DIMENSIONS +
    ! UGRID/CF-attr embedding to the M9 field-output stores. Pure refactor: same variable names
    ! written, same values, same attribute strings, same float precision, same chunking — the
    ! on-disk bytes are unchanged.
    !
    ! What "makes a Zarr store ushow/xarray-viewable" (the embedding extracted here):
    !   - the 1-D lon/lat coordinate variables (node coords AND the element-centroid coords with
    !     the r2g rotation) — io_coords_compute (owned values) + io_coords_init_lonlat (array
    !     handles) + io_coords_put (redistribute + write the chunk data);
    !   - the xarray `_ARRAY_DIMENSIONS` attribute + the CF long_name/units/standard_name on the
    !     lon/lat coordinate variables — io_coords_define_lonlat.
    !
    ! The horizontal entity is selected by `entity` (DECOMP_NODE / DECOMP_ELEM): nodes embed the
    ! geographic node coords (geo_coord_nod2D) and write data on the node decomp; elements embed
    ! the element-CENTROID coords (mod_mesh_areas elem_center: the area-weighted centroid of the
    ! sub-triangles (1,3,4)+(1,2,3) of the rotated vertices, cyclic-wrapped, r2g-rotated to
    ! geographic; hybrid-mesh plan Task 14) and write on the element decomp. FESOM2 io_r2g:3004
    ! took the plain mean of the 3 rotated vertices WITHOUT the cyclic wrap; for triangles the two
    ! agree to roundoff except on elements across the rotated dateline, where the plain mean lies
    ! on the wrong side of the globe (wrong coordinate and wrong vector rotation). The owned-entity count
    ! `nO` is the caller's partit-derived local count (means_init/mod_io_restart from local_dims).
    !
    ! The data-variable `_ARRAY_DIMENSIONS` (which differ by store shape: a growing time dim for
    ! means output, a snapshot for restart) stay with each writer — they are set via the existing
    ! zattr_str_arr primitive, not here.
    use mod_precision,   only: WP
    use mod_constants,   only: rad
    use mod_mesh,        only: t_mesh
    use mod_mesh_rotate, only: r2g
    use mod_mesh_areas,  only: elem_center
    use mod_io_zarr
    use mod_io_decomp
    implicit none
    private

    public :: io_coords_compute, io_coords_init_lonlat, io_coords_define_lonlat, io_coords_put

contains

    ! Compute the owned-entity coordinates embedded into each store, for entity = DECOMP_NODE or
    ! DECOMP_ELEM (lon/lat = GEOGRAPHIC degrees for the embed; rlon/rlat = ROTATED radians for a
    ! later vector_r2g; nlev = the 3-D below-bottom mask count). Formerly means_init:
    !   NODE: lon/lat = geo_coord_nod2D/rad; rlon/rlat = coord_nod2D; nlev = nlevels_nod2D.
    !   ELEM: rotated centroid = elem_center (rlon/rlat; the longitude may lie outside
    !         [-pi, pi] after the wrap, r2g does not care), its r2g image /rad (lon/lat),
    !         nlev = nlevels(e).
    ! The caller pre-allocates the output arrays to (at least) nO and only 1..nO is written.
    subroutine io_coords_compute(entity, mesh, nO, lon, lat, rlon, rlat, nlev)
        integer,      intent(in)  :: entity, nO
        type(t_mesh), intent(in)  :: mesh
        real(WP),     intent(out) :: lon(:), lat(:), rlon(:), rlat(:)
        integer,      intent(out) :: nlev(:)
        integer  :: i, e
        real(WP) :: rcx, rcy, gcx, gcy
        if (entity == DECOMP_ELEM) then
            do e = 1, nO
                call elem_center(mesh, e, rcx, rcy)
                rlon(e) = rcx; rlat(e) = rcy
                call r2g(gcx, gcy, rcx, rcy)               ! rotated centroid -> geographic (rad)
                lon(e) = real(gcx/rad, WP); lat(e) = real(gcy/rad, WP)
                nlev(e) = mesh%nlevels(e)
            end do
        else
            do i = 1, nO
                lon(i)  = real(mesh%geo_coord_nod2D(1,i)/rad, WP)
                lat(i)  = real(mesh%geo_coord_nod2D(2,i)/rad, WP)
                rlon(i) = real(mesh%coord_nod2D(1,i), WP)
                rlat(i) = real(mesh%coord_nod2D(2,i), WP)
                nlev(i) = mesh%nlevels_nod2D(i)
            end do
        end if
    end subroutine io_coords_compute

    ! Initialize the lon/lat coordinate-array handles (1-D, dim N, chunk C, <f8, no fill). All
    ! ranks call this (the handles are used by io_coords_put on every writer). VERBATIM.
    subroutine io_coords_init_lonlat(a_lon, a_lat, N, C)
        type(t_zarr_array), intent(out) :: a_lon, a_lat
        integer,            intent(in)  :: N, C
        call zarr_array_init(a_lon, 'lon', [N], [C], '<f8', has_fill=.false.)
        call zarr_array_init(a_lat, 'lat', [N], [C], '<f8', has_fill=.false.)
    end subroutine io_coords_init_lonlat

    ! Define the lon/lat coordinate arrays + their _ARRAY_DIMENSIONS + CF attrs in the store (the
    ! metadata-writing rank only — store-create ordering). hdim = the horizontal dim name ('nod2'
    ! node / 'elem' element-centroid). VERBATIM. Restart and means both call this.
    subroutine io_coords_define_lonlat(store, a_lon, a_lat, hdim)
        type(t_zarr_store), intent(inout) :: store
        type(t_zarr_array), intent(in)    :: a_lon, a_lat
        character(len=*),   intent(in)    :: hdim
        type(t_zarr_attrs) :: at
        ! lon / lat entity coords (node or element-centroid, dim = hdim)
        call zarr_attrs_init(at)
        call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: hdim])
        call zattr_str(at, 'long_name', 'longitude'); call zattr_str(at, 'units', 'degrees_east')
        call zattr_str(at, 'standard_name', 'longitude')
        call zarr_define_array(store, a_lon, at)
        call zarr_attrs_init(at)
        call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: hdim])
        call zattr_str(at, 'long_name', 'latitude'); call zattr_str(at, 'units', 'degrees_north')
        call zattr_str(at, 'standard_name', 'latitude')
        call zarr_define_array(store, a_lat, at)
    end subroutine io_coords_define_lonlat

    ! Redistribute a static 1-D owned entity field (node or element coord) to canonical order and
    ! let each writer write its chunks [c] (the lon/lat coord-data embed). D = the field's entity
    ! decomp. VERBATIM (formerly mod_io_means%put_static — io-independent).
    subroutine io_coords_put(D, store, arr, owned)
        type(t_io_decomp),  intent(in) :: D
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        real(WP),           intent(in) :: owned(:)
        real(WP), allocatable :: buf(:)
        integer :: c, lo
        allocate(buf(max(1, D%w_nbuf)))
        call decomp_redistribute(D, owned, buf, 0.0_WP)
        do c = D%w_first_chunk, D%w_last_chunk
            lo = (c - D%w_first_chunk)*D%C + 1
            call zarr_write_chunk(store, arr, [c], buf(lo:lo + D%C - 1))
        end do
    end subroutine io_coords_put

end module mod_io_coords
