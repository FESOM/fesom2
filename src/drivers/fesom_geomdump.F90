program fesom_geomdump
    ! M1 geometry byte-gate driver (closes deferred M0.7). Loads pi 1-rank with the
    ! SAME rotation as the FESOM2 pi run (alpha/beta/gamma = 50/15/-90, cyclic 360,
    ! force_rotation), computes geometry, and writes a binary dump in the FESOM2
    ! oracle's format. tools/geom_diff.py then compares vs the FESOM2 geom dump
    ! (tools/run_geomdump_pi.sh) for max|delta|=0.
    !
    ! Hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 14): the byte gate is
    ! retired; mod_geom_dump is REMAPPED to the FESOM2 layout (triangle meshes only, stops on
    ! quads) so geom_diff.py --signal=<tol> remains a roundoff diagnostic.
    !
    !   FESOM3_MESH_DIR  mesh dir   (default: pi)
    !   FESOM3_GEOM_OUT  out path   (default: geom_f3.bin)
    use mpi
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_geom_dump,    only: geom_dump_write
    implicit none
    character(len=512) :: mesh_dir, out_path
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    integer :: nsw

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_GEOM_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'geom_f3.bin'

    call par_init(partit)
    ! Populate the partition (myList_* + com-structs): identity synthesis at npes==1,
    ! dist_<npes>/ read at npes>1. read_mesh then builds the LOCAL mesh (M2.12a).
    call set_partition(partit, trim(mesh_dir))

    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    if (partit%mype == 0) &
        write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_geomdump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', mesh%nl
    write(*,'(a,i0,a,i0)') 'fesom_geomdump: rank ', partit%mype, ' CW swaps = ', nsw
    call geom_dump_write(mesh, partit, trim(out_path))
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_geomdump
