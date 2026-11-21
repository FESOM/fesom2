! Sea ice on its own partition (ice_dist_name in namelist.config). The ocean keeps the canonical ice
! state on its partition, which all ocean-side code, output and restarts use. Around each ice step
! the state moves to a second ice object on the ice partition and back:
! - nodes: every local node (owned and halo) of the receiving partition gets the value of the rank
!   that owns the node on the sending partition;
! - elements: every owned element (one with an owned node) of the receiving partition gets the
!   value from a rank that owns one of its nodes on the sending partition, itself if it does.
! The plans are built once from global node and element numbers; each move is one MPI_Alltoallv
! for all fields of a kind.
module ice_partition_module
    use mpi
    use o_param, only: WP
    use MOD_PARTIT
    use MOD_MESH
    use MOD_ICE
    implicit none
    private
    public :: ice_part_setup_plans, ice_part_alloc_forc, ice_part_to_ice, ice_part_to_ocean, ice_part_mesh_extras

    type t_plan
        integer, allocatable :: scnt(:), sdsp(:), rcnt(:), rdsp(:)
        integer, allocatable :: sidx(:), ridx(:)        ! local indices: packed on send, filled on receive
        integer              :: ns = 0, nr = 0
    end type t_plan

    type t_pair
        real(kind=WP), pointer :: o(:) => null(), i(:) => null()   ! ocean side, ice side
        integer                :: kind = 1                         ! 1 nodes, 2 elements
    end type t_pair

    type(t_plan), save, target :: pn_o2i, pn_i2o, pe_o2i, pe_i2o
    type(t_pair), allocatable, save :: pairs(:)
    integer, save :: npairs = 0

contains

    ! Global owner of every node on partition p
    subroutine node_owner(p, nod2D, owner)
        type(t_partit), intent(in) :: p
        integer,        intent(in) :: nod2D
        integer,        intent(out):: owner(nod2D)
        integer, allocatable :: cnt(:), dsp(:), glob(:)
        integer :: r, k, ierr
        allocate(cnt(p%npes), dsp(p%npes))
        call MPI_Allgather(p%myDim_nod2D, 1, MPI_INTEGER, cnt, 1, MPI_INTEGER, p%MPI_COMM_FESOM, ierr)
        dsp(1) = 0
        do r = 2, p%npes
            dsp(r) = dsp(r-1) + cnt(r-1)
        end do
        allocate(glob(sum(cnt)))
        call MPI_Allgatherv(p%myList_nod2D, p%myDim_nod2D, MPI_INTEGER, glob, cnt, dsp, MPI_INTEGER, &
                            p%MPI_COMM_FESOM, ierr)
        owner = -1
        do r = 1, p%npes
            do k = 1, cnt(r)
                owner(glob(dsp(r)+k)) = r-1
            end do
        end do
    end subroutine node_owner

    ! Plan to fill the receiving entities (global numbers rglob, sending rank of each in rsrc) from
    ! the sending partition, whose local entities have global numbers sglob(1:nsl)
    subroutine build_plan(comm, npes, nrl, rglob, rsrc, nsl, sglob, nglob, pl)
        integer,      intent(in)    :: comm, npes, nrl, rglob(:), rsrc(:), nsl, sglob(:), nglob
        type(t_plan), intent(inout) :: pl
        integer, allocatable :: g2l(:), req(:), pos(:), got(:)
        integer :: k, r, ierr
        allocate(pl%rcnt(npes), pl%rdsp(npes), pl%scnt(npes), pl%sdsp(npes))
        pl%rcnt = 0
        do k = 1, nrl
            pl%rcnt(rsrc(k)+1) = pl%rcnt(rsrc(k)+1) + 1
        end do
        pl%rdsp(1) = 0
        do r = 2, npes
            pl%rdsp(r) = pl%rdsp(r-1) + pl%rcnt(r-1)
        end do
        pl%nr = nrl
        allocate(pl%ridx(max(1,nrl)), req(max(1,nrl)), pos(npes))
        pos = pl%rdsp
        do k = 1, nrl
            r = rsrc(k)+1
            pos(r) = pos(r) + 1
            pl%ridx(pos(r)) = k
            req(pos(r)) = rglob(k)
        end do
        ! the requests go to the senders
        call MPI_Alltoall(pl%rcnt, 1, MPI_INTEGER, pl%scnt, 1, MPI_INTEGER, comm, ierr)
        pl%sdsp(1) = 0
        do r = 2, npes
            pl%sdsp(r) = pl%sdsp(r-1) + pl%scnt(r-1)
        end do
        pl%ns = sum(pl%scnt)
        allocate(got(max(1,pl%ns)), pl%sidx(max(1,pl%ns)))
        call MPI_Alltoallv(req, pl%rcnt, pl%rdsp, MPI_INTEGER, got, pl%scnt, pl%sdsp, &
                           MPI_INTEGER, comm, ierr)
        allocate(g2l(nglob))
        g2l = 0
        do k = nsl, 1, -1
            g2l(sglob(k)) = k
        end do
        do k = 1, pl%ns
            pl%sidx(k) = g2l(got(k))
            if (pl%sidx(k) == 0) then
                write(*,*) 'ice_partition: entity', got(k), 'requested but not local on the sender'
                call MPI_Abort(comm, 1, ierr)
            end if
        end do
    end subroutine build_plan

    ! Plans for both directions, nodes and elements (po, mo: ocean; pi, mi: sea ice)
    subroutine ice_part_setup_plans(po, mo, pi, mi)
        type(t_partit), intent(in) :: po, pi
        type(t_mesh),   intent(in) :: mo, mi
        call plans_one_way(po, mo, pi, mi, pn_o2i, pe_o2i)
        call plans_one_way(pi, mi, po, mo, pn_i2o, pe_i2o)
        if (po%mype == 0) write(*,'(a,2i9)') ' sea ice partition: nodes / elements moved per field on rank 0: ', &
             pn_o2i%nr, pe_o2i%nr
    end subroutine ice_part_setup_plans

    ! ps, ms: sending partition; pr, mr: receiving partition. The send part of a plan lives on the
    ! rank of the sending partition, the receive part on the same rank of the receiving one.
    subroutine plans_one_way(ps, ms, pr, mr, pln, ple)
        type(t_partit), intent(in) :: ps, pr
        type(t_mesh),   intent(in) :: ms, mr
        type(t_plan), intent(inout) :: pln, ple
        integer, allocatable :: own(:), rglob(:), rsrc(:)
        integer :: n, e, k, nl, nodes(3)
        allocate(own(ms%nod2D))
        call node_owner(ps, ms%nod2D, own)
        ! nodes
        nl = pr%myDim_nod2D + pr%eDim_nod2D
        allocate(rglob(nl), rsrc(nl))
        do n = 1, nl
            rglob(n) = pr%myList_nod2D(n)
            rsrc(n)  = own(rglob(n))
        end do
        call build_plan(pr%MPI_COMM_FESOM, pr%npes, nl, rglob, rsrc, ps%myDim_nod2D + ps%eDim_nod2D, &
                        ps%myList_nod2D, ms%nod2D, pln)
        deallocate(rglob, rsrc)
        ! elements
        nl = pr%myDim_elem2D
        allocate(rglob(nl), rsrc(nl))
        do e = 1, nl
            rglob(e) = pr%myList_elem2D(e)
            nodes = pr%myList_nod2D(mr%elem2D_nodes(1:3, e))
            rsrc(e) = own(nodes(1))
            do k = 1, 3
                if (own(nodes(k)) == pr%mype) rsrc(e) = pr%mype
            end do
        end do
        call build_plan(pr%MPI_COMM_FESOM, pr%npes, nl, rglob, rsrc, ps%myDim_elem2D, ps%myList_elem2D, &
                        ms%elem2D, ple)
    end subroutine plans_one_way

    ! Mesh arrays the ice step uses that the ocean setup builds outside mesh_setup: the row pointer
    ! of the SSH stiffness matrix, which the ice FCT mass matrix shares (one entry for the node and
    ! one per incident edge, as init_stiff_mat_ale), and the node neighbour lists (as muscl_adv_init)
    subroutine ice_part_mesh_extras(partit, mesh)
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(inout), target :: mesh
        integer :: n, k, n1, n2
        allocate(mesh%ssh_stiff%rowptr(partit%myDim_nod2D+1))
        mesh%ssh_stiff%rowptr(1) = 1
        do n = 1, partit%myDim_nod2D
            mesh%ssh_stiff%rowptr(n+1) = mesh%ssh_stiff%rowptr(n) + mesh%nod_in_edge2D_num(n) + 1
        end do
        mesh%nn_size = 0
        do n = 1, partit%myDim_nod2D
            mesh%nn_size = max(mesh%nn_size, mesh%ssh_stiff%rowptr(n+1) - mesh%ssh_stiff%rowptr(n))
        end do
        allocate(mesh%nn_num(partit%myDim_nod2D), mesh%nn_pos(mesh%nn_size, partit%myDim_nod2D))
        do n = 1, partit%myDim_nod2D
            mesh%nn_num(n)   = 1
            mesh%nn_pos(1,n) = n
        end do
        do n = 1, partit%myDim_nod2D
            do k = 1, mesh%nod_in_edge2D_num(n)
                n1 = mesh%nod_in_edge2D(k,n)
                n2 = mesh%edges(1,n1)
                if (n2 == n) n2 = mesh%edges(2,n1)
                mesh%nn_pos(mesh%nn_num(n)+1,n) = n2
                mesh%nn_num(n) = mesh%nn_num(n)+1
            end do
        end do
    end subroutine ice_part_mesh_extras

    ! ice%forc on the ice partition (on the ocean side ice_import allocates it)
    subroutine ice_part_alloc_forc(ice, pi)
        type(t_ice),    intent(inout) :: ice
        type(t_partit), intent(in)    :: pi
        integer :: n
        n = pi%myDim_nod2D + pi%eDim_nod2D
        associate(f => ice%forc)
        allocate(f%shortwave(n), f%longwave(n), f%tair(n), f%shum(n), f%u_wind(n), f%v_wind(n), f%runoff(n))
        allocate(f%ch_atm_oce_arr(n), f%ce_atm_oce_arr(n), f%prec_rain(n), f%prec_snow(n), f%evaporation(n))
        allocate(f%ice_sublimation(n), f%flice(n), f%fw_ice(n), f%fw_snw(n), f%real_salt_flux(n))
        allocate(f%hf_qlat(n), f%hf_qsen(n), f%hf_qradtot(n), f%hf_qswr(n), f%hf_qlwr(n), f%hf_qlwrout(n))
        allocate(f%sidm_dyn(n), f%sidm_cut(n), f%sidm_mpre(n), f%sisn_dyn(n), f%sisn_cut(n), f%sisn_spre(n))
        allocate(f%sidc_dyn(n), f%sidc_cut(n), f%sidc_apre(n))
        f%shortwave = 0; f%longwave = 0; f%tair = 0; f%shum = 0; f%u_wind = 0; f%v_wind = 0; f%runoff = 0
        f%ch_atm_oce_arr = 0; f%ce_atm_oce_arr = 0; f%prec_rain = 0; f%prec_snow = 0; f%evaporation = 0
        f%ice_sublimation = 0; f%flice = 0; f%fw_ice = 0; f%fw_snw = 0; f%real_salt_flux = 0
        f%hf_qlat = 0; f%hf_qsen = 0; f%hf_qradtot = 0; f%hf_qswr = 0; f%hf_qlwr = 0; f%hf_qlwrout = 0
        f%sidm_dyn = 0; f%sidm_cut = 0; f%sidm_mpre = 0; f%sisn_dyn = 0; f%sisn_cut = 0; f%sisn_spre = 0
        f%sidc_dyn = 0; f%sidc_cut = 0; f%sidc_apre = 0
        end associate
    end subroutine ice_part_alloc_forc

    subroutine add(o, i, kind)
        real(kind=WP), allocatable, target, intent(inout) :: o(:), i(:)
        integer, intent(in) :: kind
        type(t_pair), allocatable :: tmp(:)
        if (.not. (allocated(o) .and. allocated(i))) return
        if (.not. allocated(pairs)) allocate(pairs(64))
        if (npairs == size(pairs)) then
            allocate(tmp(2*npairs)); tmp(1:npairs) = pairs; call move_alloc(tmp, pairs)
        end if
        npairs = npairs + 1
        pairs(npairs)%o => o
        pairs(npairs)%i => i
        pairs(npairs)%kind = kind
    end subroutine add

    ! All ice fields that carry state or coupling between steps. Not moved: geometry-derived
    ! constants of each partition (inv_mass, inv_areamass, fct_massmatrix) and FCT flux scratch.
    subroutine build_pairs(io, ii)
        type(t_ice), intent(inout), target :: io, ii
        integer :: k
        npairs = 0
        do k = 1, min(size(io%data), size(ii%data))
            call add(io%data(k)%values,         ii%data(k)%values,         1)
            call add(io%data(k)%values_old,     ii%data(k)%values_old,     1)
            call add(io%data(k)%values_rhs,     ii%data(k)%values_rhs,     1)
            call add(io%data(k)%values_div_rhs, ii%data(k)%values_div_rhs, 1)
            call add(io%data(k)%dvalues,        ii%data(k)%dvalues,        1)
            call add(io%data(k)%valuesl,        ii%data(k)%valuesl,        1)
        end do
        call add(io%uice, ii%uice, 1);         call add(io%vice, ii%vice, 1)
        call add(io%uice_rhs, ii%uice_rhs, 1); call add(io%vice_rhs, ii%vice_rhs, 1)
        call add(io%uice_old, ii%uice_old, 1); call add(io%vice_old, ii%vice_old, 1)
        call add(io%uice_aux, ii%uice_aux, 1); call add(io%vice_aux, ii%vice_aux, 1)
        call add(io%uice_ib, ii%uice_ib, 1);   call add(io%vice_ib, ii%vice_ib, 1)
        call add(io%stress_atmice_x, ii%stress_atmice_x, 1); call add(io%stress_atmice_y, ii%stress_atmice_y, 1)
        call add(io%stress_iceoce_x, ii%stress_iceoce_x, 1); call add(io%stress_iceoce_y, ii%stress_iceoce_y, 1)
        call add(io%srfoce_temp, ii%srfoce_temp, 1); call add(io%srfoce_salt, ii%srfoce_salt, 1)
        call add(io%srfoce_ssh, ii%srfoce_ssh, 1)
        call add(io%srfoce_u, ii%srfoce_u, 1);       call add(io%srfoce_v, ii%srfoce_v, 1)
        call add(io%flx_fw, ii%flx_fw, 1);           call add(io%flx_h, ii%flx_h, 1)
        call add(io%h_ice, ii%h_ice, 1);             call add(io%h_snow, ii%h_snow, 1)
        call add(io%work%fct_tmax, ii%work%fct_tmax, 1);   call add(io%work%fct_tmin, ii%work%fct_tmin, 1)
        call add(io%work%fct_plus, ii%work%fct_plus, 1);   call add(io%work%fct_minus, ii%work%fct_minus, 1)
        call add(io%work%ice_strength, ii%work%ice_strength, 1)
        call add(io%thermo%t_skin, ii%thermo%t_skin, 1);   call add(io%thermo%thdgr, ii%thermo%thdgr, 1)
        call add(io%thermo%thdgrsn, ii%thermo%thdgrsn, 1); call add(io%thermo%thdgra, ii%thermo%thdgra, 1)
        call add(io%thermo%thdgr_old, ii%thermo%thdgr_old, 1); call add(io%thermo%ustar, ii%thermo%ustar, 1)
        call add(io%thermo%dyngr, ii%thermo%dyngr, 1);     call add(io%thermo%dyngrsn, ii%thermo%dyngrsn, 1)
        call add(io%thermo%dyngra, ii%thermo%dyngra, 1)
        call add(io%thermo%apnd, ii%thermo%apnd, 1); call add(io%thermo%hpnd, ii%thermo%hpnd, 1)
        call add(io%thermo%ipnd, ii%thermo%ipnd, 1)
        associate(a => io%forc, b => ii%forc)
        call add(a%shortwave, b%shortwave, 1); call add(a%longwave, b%longwave, 1)
        call add(a%tair, b%tair, 1);           call add(a%shum, b%shum, 1)
        call add(a%u_wind, b%u_wind, 1);       call add(a%v_wind, b%v_wind, 1)
        call add(a%runoff, b%runoff, 1)
        call add(a%ch_atm_oce_arr, b%ch_atm_oce_arr, 1); call add(a%ce_atm_oce_arr, b%ce_atm_oce_arr, 1)
        call add(a%prec_rain, b%prec_rain, 1); call add(a%prec_snow, b%prec_snow, 1)
        call add(a%evaporation, b%evaporation, 1); call add(a%ice_sublimation, b%ice_sublimation, 1)
        call add(a%flice, b%flice, 1);         call add(a%fw_ice, b%fw_ice, 1)
        call add(a%fw_snw, b%fw_snw, 1);       call add(a%real_salt_flux, b%real_salt_flux, 1)
        call add(a%hf_qlat, b%hf_qlat, 1);     call add(a%hf_qsen, b%hf_qsen, 1)
        call add(a%hf_qradtot, b%hf_qradtot, 1); call add(a%hf_qswr, b%hf_qswr, 1)
        call add(a%hf_qlwr, b%hf_qlwr, 1);     call add(a%hf_qlwrout, b%hf_qlwrout, 1)
        call add(a%sidm_dyn, b%sidm_dyn, 1);   call add(a%sidm_cut, b%sidm_cut, 1)
        call add(a%sidm_mpre, b%sidm_mpre, 1); call add(a%sisn_dyn, b%sisn_dyn, 1)
        call add(a%sisn_cut, b%sisn_cut, 1);   call add(a%sisn_spre, b%sisn_spre, 1)
        call add(a%sidc_dyn, b%sidc_dyn, 1);   call add(a%sidc_cut, b%sidc_cut, 1)
        call add(a%sidc_apre, b%sidc_apre, 1)
        end associate
        call add(io%work%sigma11, ii%work%sigma11, 2); call add(io%work%sigma12, ii%work%sigma12, 2)
        call add(io%work%sigma22, ii%work%sigma22, 2)
        call add(io%work%eps11, ii%work%eps11, 2);     call add(io%work%eps12, ii%work%eps12, 2)
        call add(io%work%eps22, ii%work%eps22, 2)
        call add(io%alpha_evp_array, ii%alpha_evp_array, 2); call add(io%beta_evp_array, ii%beta_evp_array, 2)
    end subroutine build_pairs

    ! Move all fields of one kind: to_ice selects the direction
    subroutine move(kind, to_ice, comm, npes)
        integer, intent(in) :: kind, comm, npes
        logical, intent(in) :: to_ice
        type(t_plan), pointer :: ps, pr
        real(kind=WP), allocatable :: sb(:), rb(:)
        integer, allocatable :: sc(:), sd(:), rc(:), rd(:)
        integer :: p, nf, f, k, r, ierr
        real(kind=WP), pointer :: src(:), dst(:)
        if (kind == 1) then
            if (to_ice) then; ps => pn_o2i; else; ps => pn_i2o; end if
        else
            if (to_ice) then; ps => pe_o2i; else; ps => pe_i2o; end if
        end if
        pr => ps
        nf = count(pairs(1:npairs)%kind == kind)
        allocate(sb(max(1, nf*ps%ns)), rb(max(1, nf*pr%nr)), sc(npes), sd(npes), rc(npes), rd(npes))
        sc = nf*ps%scnt; sd = nf*ps%sdsp; rc = nf*pr%rcnt; rd = nf*pr%rdsp
        f = 0
        do p = 1, npairs
            if (pairs(p)%kind /= kind) cycle
            f = f + 1
            if (to_ice) then; src => pairs(p)%o; else; src => pairs(p)%i; end if
            do r = 1, npes
                do k = 1, ps%scnt(r)
                    sb(sd(r) + (f-1)*ps%scnt(r) + k) = src(ps%sidx(ps%sdsp(r) + k))
                end do
            end do
        end do
        call MPI_Alltoallv(sb, sc, sd, MPI_WP, rb, rc, rd, MPI_WP, comm, ierr)
        f = 0
        do p = 1, npairs
            if (pairs(p)%kind /= kind) cycle
            f = f + 1
            if (to_ice) then; dst => pairs(p)%i; else; dst => pairs(p)%o; end if
            do r = 1, npes
                do k = 1, pr%rcnt(r)
                    dst(pr%ridx(pr%rdsp(r) + k)) = rb(rd(r) + (f-1)*pr%rcnt(r) + k)
                end do
            end do
        end do
    end subroutine move

    subroutine ice_part_to_ice(io, ii, po)
        type(t_ice),    intent(inout), target :: io, ii
        type(t_partit), intent(in)            :: po
        if (npairs == 0) call build_pairs(io, ii)
        ii%ice_update          = io%ice_update
        ii%ice_steps_since_upd = io%ice_steps_since_upd
        call move(1, .true., po%MPI_COMM_FESOM, po%npes)
        call move(2, .true., po%MPI_COMM_FESOM, po%npes)
    end subroutine ice_part_to_ice

    subroutine ice_part_to_ocean(io, ii, po)
        type(t_ice),    intent(inout), target :: io, ii
        type(t_partit), intent(in)            :: po
        call move(1, .false., po%MPI_COMM_FESOM, po%npes)
        call move(2, .false., po%MPI_COMM_FESOM, po%npes)
    end subroutine ice_part_to_ocean
end module ice_partition_module
