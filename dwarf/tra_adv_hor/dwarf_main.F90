! Horizontal tracer advection dwarf. Each MPI rank reads its dump (dump_dir/tra_adv_<rank>.bin) and
! times the per-tracer horizontal sequence: low-order flux, low-order tendency gather, high-order
! (MFCT) flux, horizontal flux-to-tendency gather. No communication: ranks only share the node's
! memory bandwidth, which is what the 1 vs 8 rank comparison measures.
! Variants: ref (FESOM code, arrays passed through pointers as in do_oce_adv_tra), ref_c (same,
! contiguous pointers), merged (one level loop per edge, mesh-only edge data precomputed),
! fused (low- and high-order flux in one edge pass, then the gathers).
! Usage: mpirun -n N dwarf.x <dump_dir> [repeats]
program dwarf_tra_adv_hor
    use mpi
    use o_PARAM
    use MOD_MESH
    use MOD_PARTIT
    use dwarf_data
    use kernels_ref
    use kernels_opt
    implicit none
    type(t_case), target :: c
    type(t_edge_meta) :: m
    character(len=512) :: dir, arg
    character(len=600) :: fname
    integer :: rank, npes, ierr, nrep, it, nl, ne, nn, nel, k, nv
    real(kind=WP), allocatable, target :: flux_lo(:,:), flux(:,:)
    real(kind=WP_full), allocatable, target :: lo_tend(:,:), dttf_h(:,:)
    real(kind=WP), allocatable :: scale(:,:)
    real(kind=WP), allocatable :: chk_lo(:,:), chk_tend(:,:), chk_ho(:,:), chk_dt(:,:)
    real(kind=WP), pointer :: pttf(:,:), pttfab(:,:), pflux(:,:)
    real(kind=WP), pointer, contiguous :: cttf(:,:), cttfab(:,:), cflux(:,:)
    real(real64) :: t(4), tmin(4, 4), tt
    character(len=8) :: vname(4) = ['ref     ', 'ref_c   ', 'merged  ', 'fused   ']
    real(real64) :: work_lay, work_edges

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)
    call get_command_argument(1, dir)
    nrep = 20
    if (command_argument_count() > 1) then
        call get_command_argument(2, arg); read(arg, *) nrep
    end if
    write(fname, '(a,a,i5.5,a)') trim(dir), '/tra_adv_', rank, '.bin'
    call read_dump(trim(fname), c)
    nl = c%mesh%nl; ne = c%partit%myDim_edge2D; nn = c%partit%myDim_nod2D + c%partit%eDim_nod2D
    nel = c%partit%myDim_elem2D + c%partit%eDim_elem2D
    allocate(flux_lo(nl-1, ne), flux(nl-1, ne), lo_tend(nl-1, nn), dttf_h(nl-1, nn), scale(nl-1, nn))
    scale = 0
    do k = 1, nn
        where (c%mesh%areasvol(1:nl-1, k) > 0) scale(:, k) = c%dt/c%mesh%areasvol(1:nl-1, k)
    end do
    call edge_meta_init(m, c%partit, c%mesh)
    work_edges = ne
    work_lay = sum(real(m%nzb - m%nzt + 1, real64))
    pttf => c%ttf; pttfab => c%ttfAB; pflux => flux
    cttf => c%ttf; cttfab => c%ttfAB; cflux => flux

    tmin = huge(1.0_real64)
    do nv = 1, 4
        do it = 0, nrep
            t = 0
            select case (nv)
            case (1, 2)
                call MPI_Barrier(MPI_COMM_WORLD, ierr); tt = MPI_Wtime()
                if (nv == 1) then
                    call adv_tra_hor_upw1(c%vel, pttf, c%partit, c%mesh, flux_lo, o_init_zero=.true.)
                else
                    call adv_tra_hor_upw1(c%vel, cttf, c%partit, c%mesh, flux_lo, o_init_zero=.true.)
                end if
                t(1) = MPI_Wtime() - tt; tt = MPI_Wtime()
                call lo_gather(flux_lo, lo_tend, c%partit, c%mesh)
                t(2) = MPI_Wtime() - tt
                flux = flux_lo
                tt = MPI_Wtime()
                if (nv == 1) then
                    call adv_tra_hor_mfct(c%vel, pttfab, c%partit, c%mesh, c%opth, pflux, c%grad, o_init_zero=.false.)
                else
                    call adv_tra_hor_mfct(c%vel, cttfab, c%partit, c%mesh, c%opth, cflux, c%grad, o_init_zero=.false.)
                end if
                t(3) = MPI_Wtime() - tt
                dttf_h = 0; tt = MPI_Wtime()
                call hor_gather(c%dt, dttf_h, flux, c%partit, c%mesh)
                t(4) = MPI_Wtime() - tt
            case (3)
                call MPI_Barrier(MPI_COMM_WORLD, ierr); tt = MPI_Wtime()
                call upw1_merged(m, c%vel, c%ttf, c%mesh%helem, flux_lo, nl, ne, nn, nel)
                t(1) = MPI_Wtime() - tt
                lo_tend = 0; tt = MPI_Wtime()
                call gather_meta(m, flux_lo, c%mesh%nod_in_edge2D_sgn, c%mesh%nod_in_edge2D_num, c%mesh%nod_in_edge2D, &
                                 lo_tend, nl=nl, ne=ne, nn=nn)
                t(2) = MPI_Wtime() - tt
                flux = flux_lo
                tt = MPI_Wtime()
                call mfct_merged(m, c%opth, c%vel, c%ttfAB, c%mesh%helem, c%grad, flux, nl, ne, nn, nel)
                t(3) = MPI_Wtime() - tt
                dttf_h = 0; tt = MPI_Wtime()
                call gather_meta(m, flux, c%mesh%nod_in_edge2D_sgn, c%mesh%nod_in_edge2D_num, c%mesh%nod_in_edge2D, &
                                 dttf_h, scale, nl, ne, nn)
                t(4) = MPI_Wtime() - tt
            case (4)
                call MPI_Barrier(MPI_COMM_WORLD, ierr); tt = MPI_Wtime()
                call fused_upw1_mfct(m, c%opth, c%vel, c%ttf, c%ttfAB, c%mesh%helem, c%grad, flux_lo, flux, nl, ne, nn, nel)
                t(1) = MPI_Wtime() - tt
                lo_tend = 0; tt = MPI_Wtime()
                call gather_meta(m, flux_lo, c%mesh%nod_in_edge2D_sgn, c%mesh%nod_in_edge2D_num, c%mesh%nod_in_edge2D, &
                                 lo_tend, nl=nl, ne=ne, nn=nn)
                t(2) = MPI_Wtime() - tt
                dttf_h = 0; tt = MPI_Wtime()
                call gather_meta(m, flux, c%mesh%nod_in_edge2D_sgn, c%mesh%nod_in_edge2D_num, c%mesh%nod_in_edge2D, &
                                 dttf_h, scale, nl, ne, nn)
                t(4) = MPI_Wtime() - tt
            end select
            if (it > 0) tmin(:, nv) = min(tmin(:, nv), t)
        end do
        ! results of this variant
        if (nv == 1) then
            call check('ref flux_lo vs FESOM', flux_lo, c%flux_lo_ref, 0)
            call check('ref lo_tend vs FESOM', real(lo_tend, WP), c%lo_tend_ref, 0)
            call check('ref flux_ho vs FESOM', flux, c%flux_ho_ref, 0)
            chk_lo = flux_lo; chk_tend = real(lo_tend, WP); chk_ho = flux; chk_dt = real(dttf_h, WP)
        else
            call check(trim(vname(nv))//' flux_lo', flux_lo, chk_lo, 1)
            call check(trim(vname(nv))//' lo_tend', real(lo_tend, WP), chk_tend, 1)
            call check(trim(vname(nv))//' flux_ho', flux, chk_ho, 1)
            call check(trim(vname(nv))//' dttf_h ', real(dttf_h, WP), chk_dt, 1)
        end if
    end do
    call report()
    call rank_file()
    call MPI_Finalize(ierr)

contains

    ! mode 0: must be bit-identical; 1: report the largest difference relative to the field's maximum
    subroutine check(what, a, b, mode)
        character(len=*), intent(in) :: what
        real(kind=WP), intent(in) :: a(:,:), b(:,:)
        integer, intent(in) :: mode
        real(real64) :: d, s, dg, sg
        integer :: n
        d = maxval(abs(a - b)); s = maxval(abs(b))
        n = count(a /= b)
        call MPI_Allreduce(d, dg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(s, sg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, n, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (rank == 0) then
            if (n == 0) then
                write(*,'(a,a)') '  check ', what//': bit-identical'
            else
                write(*,'(a,a,i10,a,es10.2)') '  check ', what//': values differing', n, ', max rel diff', dg/max(sg, tiny(1.0_real64))
            end if
        end if
        if (mode == 0 .and. n /= 0 .and. rank == 0) write(*,*) '  ** reference kernel does not reproduce FESOM **'
    end subroutine check

    subroutine report()
        real(real64) :: tmax(4,4), tsum(4,4), wl(2), wls(2), wlmax(2)
        call MPI_Reduce(tmin, tmax, 16, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce(tmin, tsum, 16, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        wl = [work_edges, work_lay]
        call MPI_Reduce(wl, wls, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce(wl, wlmax, 2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
        if (rank /= 0) return
        write(*,'(a,i4,a,i4,a,f12.0,a,f12.0)') ' ranks', npes, '  repeats', nrep, '  edges (sum)', wls(1), &
             '  edge-levels (sum)', wls(2)
        write(*,'(a,2f8.3)') ' work imbalance max/mean: edges, edge-levels', wlmax(1)*npes/wls(1), wlmax(2)*npes/wls(2)
        write(*,'(a)') ' best time per call in ms, slowest rank (mean over ranks in brackets)'
        write(*,'(a8,5a20)') 'variant', 'LO flux', 'LO gather', 'HO flux', 'flux gather', 'total'
        do nv = 1, 4
            write(*,'(a8,5(f10.3,a,f7.3,a))') vname(nv), (1e3*tmax(k,nv), ' (', 1e3*tsum(k,nv)/npes, ')', k = 1, 4), &
                 1e3*sum(tmax(:,nv)), ' (', 1e3*sum(tsum(:,nv))/npes, ')'
        end do
    end subroutine report

    ! per-rank work and time of the reference sequence (dwarf_ranks.txt), for a cost model
    subroutine rank_file()
        real(real64) :: row(7)
        real(real64), allocatable :: tab(:,:)
        integer :: r, u
        row = [work_edges, work_lay, real(c%partit%myDim_nod2D, real64), tmin(1,1), tmin(2,1), tmin(3,1), tmin(4,1)]
        allocate(tab(7, npes))
        call MPI_Gather(row, 7, MPI_DOUBLE_PRECISION, tab, 7, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (rank /= 0) return
        open(newunit=u, file='dwarf_ranks.txt', status='replace')
        write(u, '(a)') 'rank edges edge_levels own_nodes t_lo t_lo_gather t_ho t_gather'
        do r = 1, npes
            write(u, '(i6, 3f12.0, 4es14.6)') r-1, tab(:, r)
        end do
        close(u)
    end subroutine rank_file
end program dwarf_tra_adv_hor
