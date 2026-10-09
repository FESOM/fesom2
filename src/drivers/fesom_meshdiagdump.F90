program fesom_meshdiagdump
    ! M9 Stage 1 driver: load a mesh (1-rank or dist_<npes>), compute geometry, and write
    ! fesom.mesh.diag.zarr via mod_io_meshdiag. Mirrors fesom_geomdump's standalone mesh setup
    ! (same pi rotation 50/15/-90, cyclic 360, force_rotation) so geo_coord_nod2D matches the
    ! FESOM2 oracle. tools/zarr_diff.py --meshdiag then compares vs FESOM2 fesom.mesh.diag.nc.
    !
    ! Hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 14): REMAPPED: the
    ! store has the hybrid layout (4-column face_nodes/face_edges/face_links/gradient_sca_x/y,
    ! a triangle's 4th = _FillValue). tools/check_meshdiag.py checks it (ctest meshdiag_*);
    ! the zarr_diff.py --meshdiag FESOM2 gate is retired.
    !
    !   FESOM3_MESH_DIR      mesh dir (default: pi)
    !   FESOM3_MESHDIAG_OUT  output store path (default: ./fesom.mesh.diag.zarr)
    !   FESOM3_CHUNK_HORIZ   horizontal chunk size (default 500000; clamp to N)
    !   FESOM3_N_WRITERS     writer-rank subset (default auto = min(npes, nchunks))
    use mpi
    use mod_precision,    only: WP, MP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_part_bounds,  only: local_dims
    use mod_io_meshdiag,  only: meshdiag_write
    implicit none
    character(len=512) :: mesh_dir, out_path
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    integer :: nsw, e, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_MESHDIAG_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'fesom.mesh.diag.zarr'

    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)

    ! zbar_e_bot is set by the lifecycle drivers (not compute_geometry); replicate the full-cell
    ! value zbar(nlevels(e)) here so the standalone mesh.diag has it (mirror fesom_lifecycle_mr).
    call local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)
    if (.not. allocated(mesh%zbar_e_bot)) then
        allocate(mesh%zbar_e_bot(nElemF)); mesh%zbar_e_bot = 0.0_MP
        do e = 1, nElemF
            mesh%zbar_e_bot(e) = mesh%zbar(mesh%nlevels(e))
        end do
    end if

    if (partit%mype == 0) &
        write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_meshdiagdump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', mesh%nl

    call meshdiag_write(trim(out_path), mesh, partit)

    if (partit%mype == 0) write(*,'(a)') 'MESHDIAGDUMP OK'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_meshdiagdump
