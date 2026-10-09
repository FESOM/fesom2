module test_util
    ! Scaffolding shared by the hybrid-mesh self-tests (test_coupling_mixed, test_geom_mixed,
    ! test_ice_dyn_mixed, test_ice_mass, test_pgf_mixed, test_ssh_stiff, test_updown_mixed):
    ! the mesh setup from FESOM3_MESH_DIR / FESOM3_TEST_ROTATE, the failure counter and the
    ! final report, cross-rank reductions (mod_halo) that are identities at 1 rank -- an
    ! analytic patch's own 1-rank partition has no communicator --, node values hashed from
    ! the global id, and the full-cell depths the SSH matrix needs.
    use mod_precision,    only: WP, MP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_part_bounds,  only: owned_bounds
    use mod_halo,         only: allreduce_sum, allreduce_max
    implicit none
    private
    public :: setup_test_mesh, check, finish, max_all, sum_all, hval, set_depths

    integer :: nfail = 0      ! failed checks on this rank

    interface sum_all
        module procedure sum_all_r0, sum_all_r1, sum_all_i0
    end interface

contains

    subroutine setup_test_mesh(mesh, partit, default_dir, cartesian, partition_np1)
        ! mesh = FESOM3_MESH_DIR (default default_dir) read unrotated or, with
        ! FESOM3_TEST_ROTATE=1, with the standard rotation 50/15/-90 (pi, CORE2 -- the
        ! unrotated pi has elements around the geographic pole, which have no planar lon/lat
        ! embedding); compute_geometry(cartesian). The partition is read at N ranks and, if
        ! partition_np1, also at 1 rank (the 1-rank synthesis).
        type(t_mesh),     intent(inout) :: mesh
        type(t_partit),   intent(inout) :: partit
        character(len=*), intent(in)    :: default_dir
        logical,          intent(in)    :: cartesian, partition_np1
        character(len=512) :: mesh_dir, cbuf
        integer :: nsw
        mesh_dir = ''
        call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
        if (len_trim(mesh_dir) == 0) mesh_dir = default_dir
        if (partit%npes > 1 .or. partition_np1) call set_partition(partit, trim(mesh_dir))
        cbuf = ''
        call get_environment_variable('FESOM3_TEST_ROTATE', cbuf)
        if (trim(cbuf) == '1') then
            call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                           force_rotation=.true., n_cw_swaps=nsw)
        else
            call read_mesh(mesh, partit, trim(mesh_dir), 0.0_WP, 0.0_WP, 0.0_WP, 360.0_WP, &
                           force_rotation=.false., n_cw_swaps=nsw)
        end if
        call compute_geometry(mesh, partit, cartesian=cartesian)
    end subroutine setup_test_mesh

    subroutine check(cond, name, mype)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: name
        integer,          intent(in) :: mype
        if (.not. cond) then
            nfail = nfail + 1
            write(*,'(a,i0,a)') '  FAIL [rank ', mype, ']: '//name
        end if
    end subroutine check

    subroutine finish(partit, test_name)
        ! the failures of all ranks; "<test_name> (npes=N): ALL PASS" or the count; MPI
        ! finalised; error stop 1 on a failure
        type(t_partit),   intent(inout) :: partit
        character(len=*), intent(in)    :: test_name
        integer :: nfail_g
        nfail_g = nfail
        call sum_all(nfail_g, partit)
        if (partit%mype == 0) then
            if (nfail_g == 0) then
                write(*,'(a,i0,a)') test_name//' (npes=', partit%npes, '): ALL PASS'
            else
                write(*,'(a,i0,a,i0,a)') test_name//' (npes=', partit%npes, '): ', nfail_g, ' FAILURE(S)'
            end if
        end if
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        if (nfail_g /= 0) error stop 1
    end subroutine finish

    subroutine max_all(x, p)
        real(kind=WP),  intent(inout) :: x
        type(t_partit), intent(in)    :: p
        if (p%npes > 1) call allreduce_max(x, p)
    end subroutine max_all

    subroutine sum_all_r0(x, p)
        real(kind=WP),  intent(inout) :: x
        type(t_partit), intent(in)    :: p
        if (p%npes > 1) call allreduce_sum(x, p)
    end subroutine sum_all_r0

    subroutine sum_all_r1(x, p)
        real(kind=WP),  intent(inout) :: x(:)
        type(t_partit), intent(in)    :: p
        if (p%npes > 1) call allreduce_sum(x, p)
    end subroutine sum_all_r1

    subroutine sum_all_i0(k, p)
        integer,        intent(inout) :: k
        type(t_partit), intent(in)    :: p
        if (p%npes > 1) call allreduce_sum(k, p)
    end subroutine sum_all_i0

    real(kind=WP) function hval(gid, k)
        ! a node value hashed from the global id (identical on owned and halo copies)
        integer, intent(in) :: gid, k
        hval = sin(1.7_WP*real(gid, WP) + 0.9_WP*real(k, WP))
    end function hval

    subroutine set_depths(m, p)
        ! full-cell bottom depth and a zero SSH change, as the drivers set them (the
        ! init_stiff_mat_ale inputs)
        type(t_mesh),   intent(inout) :: m
        type(t_partit), intent(in)    :: p
        integer :: nNodO, nNodL, nEdgeO, nElemO, e
        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        allocate(m%zbar_e_bot(nElemO), m%dhe(nElemO))
        do e = 1, nElemO
            m%zbar_e_bot(e) = m%zbar(m%nlevels(e))
        end do
        m%dhe = 0.0_MP
    end subroutine set_depths

end module test_util
