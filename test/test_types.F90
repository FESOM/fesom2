program test_types
    ! M0.3 gate: allocate every derived type, serialize via write(unformatted),
    ! deserialize into a fresh instance via read(unformatted), assert max|delta|=0.
    use mpi
    use mod_precision, only: WP, MP
    use mod_mesh,   only: t_mesh
    use mod_partit, only: t_partit
    use mod_dyn,    only: t_dyn
    use mod_tracer, only: t_tracer
    use mod_ice,    only: t_ice
    implicit none

    integer :: ierr, nfail
    nfail = 0
    call MPI_Init(ierr)

    call test_mesh()
    call test_partit()
    call test_dyn()
    call test_tracer()
    call test_ice()

    if (nfail == 0) then
        write(*,'(a)') 'test_types: ALL PASS'
    else
        write(*,'(a,i0,a)') 'test_types: ', nfail, ' FAILURE(S)'
    end if
    call MPI_Finalize(ierr)
    if (nfail /= 0) error stop 1

contains

    subroutine fail(msg)
        character(len=*), intent(in) :: msg
        nfail = nfail + 1
        write(*,'(a)') '  FAIL: '//msg
    end subroutine

    integer function scratch_write_open() result(u)
        open(newunit=u, file='test_types.bin', form='unformatted', &
             access='stream', status='replace', action='write')
    end function
    integer function scratch_read_open() result(u)
        open(newunit=u, file='test_types.bin', form='unformatted', &
             access='stream', status='old', action='read')
    end function

    subroutine test_mesh()
        type(t_mesh) :: a, b
        integer :: u, ios, i, j
        a%nod2D = 7; a%elem2D = 5; a%edge2D = 11; a%edge2D_in = 9; a%nl = 4
        a%ocean_area = 1234.5_MP
        allocate(a%coord_nod2D(2, a%nod2D), a%elem2D_nodes(4, a%elem2D))
        allocate(a%elem2D_nnodes(a%elem2D), a%nlevels(a%elem2D))
        allocate(a%area(a%nod2D), a%hnode(a%nl, a%nod2D))   ! area is 1-D (bottom at vertices)
        allocate(a%ssh_stiff%rowptr(a%nod2D+1), a%ssh_stiff%values(13))
        a%ssh_stiff%dim = a%nod2D; a%ssh_stiff%nza = 13
        do j = 1, a%nod2D
            a%coord_nod2D(1, j) = real(j, MP) * 0.5_MP
            a%coord_nod2D(2, j) = real(j, MP) * (-0.25_MP)
        end do
        do j = 1, a%elem2D
            a%elem2D_nnodes(j) = 3
            a%nlevels(j) = mod(j, a%nl) + 1
            do i = 1, 4
                a%elem2D_nodes(i, j) = i + j
            end do
        end do
        a%area  = [(real(i, MP)*1.5_MP, i=1, a%nod2D)]
        a%hnode = reshape([(real(i, MP)*3.0_MP, i=1, a%nl*a%nod2D)], [a%nl, a%nod2D])
        a%ssh_stiff%rowptr = [(i, i=1, a%nod2D+1)]
        a%ssh_stiff%values = [(real(i, MP)*3.0_MP, i=1, 13)]

        u = scratch_write_open(); write(u, iostat=ios) a; close(u)
        u = scratch_read_open();  read(u,  iostat=ios) b; close(u)

        if (b%nod2D /= a%nod2D .or. b%elem2D /= a%elem2D .or. b%nl /= a%nl) call fail('mesh scalars')
        if (b%ocean_area /= a%ocean_area) call fail('mesh ocean_area')
        if (maxval(abs(b%coord_nod2D - a%coord_nod2D)) /= 0.0_MP) call fail('mesh coord_nod2D')
        if (any(b%elem2D_nodes /= a%elem2D_nodes)) call fail('mesh elem2D_nodes')
        if (any(b%elem2D_nnodes /= a%elem2D_nnodes)) call fail('mesh elem2D_nnodes')
        if (any(b%nlevels /= a%nlevels)) call fail('mesh nlevels')
        if (maxval(abs(b%area - a%area)) /= 0.0_MP) call fail('mesh area')
        if (maxval(abs(b%hnode - a%hnode)) /= 0.0_MP) call fail('mesh hnode')
        if (b%ssh_stiff%nza /= a%ssh_stiff%nza) call fail('mesh ssh_stiff%nza')
        if (maxval(abs(b%ssh_stiff%values - a%ssh_stiff%values)) /= 0.0_MP) call fail('mesh ssh_stiff%values')
        if (any(b%ssh_stiff%rowptr /= a%ssh_stiff%rowptr)) call fail('mesh ssh_stiff%rowptr')
    end subroutine

    subroutine test_partit()
        type(t_partit) :: a, b
        integer :: u, ios, i
        a%npes = 4; a%mype = 2
        a%myDim_nod2D = 6; a%eDim_nod2D = 2
        a%myDim_elem2D = 5; a%eDim_elem2D = 1; a%eXDim_elem2D = 1
        a%myDim_edge2D = 9; a%eDim_edge2D = 3; a%pe_status = 0
        allocate(a%part(20), a%myList_nod2D(a%myDim_nod2D + a%eDim_nod2D))
        allocate(a%myList_elem2D(6), a%myList_edge2D(12))
        a%part = [(mod(i,4), i=1,20)]
        a%myList_nod2D = [(i*3, i=1, size(a%myList_nod2D))]
        a%myList_elem2D = [(i*7, i=1,6)]
        a%myList_edge2D = [(i*2, i=1,12)]
        a%com_nod2D%rPEnum = 2; a%com_nod2D%sPEnum = 1; a%com_nod2D%nreq = 3
        a%com_nod2D%rPE(1:2) = [1, 3]
        a%com_nod2D%rptr(1:3) = [1, 3, 5]
        allocate(a%com_nod2D%rlist(4), a%com_nod2D%slist(2))
        a%com_nod2D%rlist = [10, 11, 12, 13]
        a%com_nod2D%slist = [20, 21]
        a%eXDim_nod2D = 1; a%dist_format = 1
        a%com_edge2D%rPEnum = 1; a%com_edge2D%rptr(1:2) = [1, 3]
        allocate(a%com_edge2D%rlist(2), a%com_edge2D%slist(0))
        a%com_edge2D%rlist = [5, 7]

        u = scratch_write_open(); write(u, iostat=ios) a; close(u)
        u = scratch_read_open();  read(u,  iostat=ios) b; close(u)

        if (b%npes /= a%npes .or. b%mype /= a%mype) call fail('partit npes/mype')
        if (b%myDim_nod2D /= a%myDim_nod2D .or. b%eDim_nod2D /= a%eDim_nod2D) call fail('partit nod dims')
        if (b%eXDim_elem2D /= a%eXDim_elem2D) call fail('partit eXDim_elem2D')
        if (any(b%part /= a%part)) call fail('partit part')
        if (any(b%myList_nod2D /= a%myList_nod2D)) call fail('partit myList_nod2D')
        if (any(b%myList_elem2D /= a%myList_elem2D)) call fail('partit myList_elem2D')
        if (b%com_nod2D%rPEnum /= a%com_nod2D%rPEnum) call fail('partit com rPEnum')
        if (any(b%com_nod2D%rptr(1:3) /= a%com_nod2D%rptr(1:3))) call fail('partit com rptr')
        if (any(b%com_nod2D%rlist /= a%com_nod2D%rlist)) call fail('partit com rlist')
        if (any(b%com_nod2D%slist /= a%com_nod2D%slist)) call fail('partit com slist')
        if (b%eXDim_nod2D /= a%eXDim_nod2D .or. b%dist_format /= a%dist_format) &
            call fail('partit eXDim_nod2D/dist_format')
        if (any(b%com_edge2D%rlist /= a%com_edge2D%rlist)) call fail('partit com_edge2D rlist')
    end subroutine

    subroutine test_dyn()
        type(t_dyn) :: a, b
        integer :: u, ios, i
        integer, parameter :: ne = 5, nn = 7, nl = 4
        a%opt_visc = 7; a%momadv_opt = 2; a%AB_order = 2
        a%visc_gamma0 = 0.025_WP
        ! the w split parameters, all NON-default (type defaults .false. / 0.5 / 1.0)
        a%use_wsplit = .true.; a%wsplit_mincfl = 0.3_WP; a%wsplit_maxcfl = 0.8_WP
        allocate(a%uv(2, nl-1, ne), a%uv_rhsAB(2, nl-1, ne, 2))
        allocate(a%uvnode(2, nl-1, nn), a%w(nl, nn))
        allocate(a%eta_n(nn), a%ssh_rhs(nn))
        allocate(a%solverinfo%rr(nn))
        allocate(a%work%uvnode_rhs(2, nl-1, nn))
        a%uv = reshape([(real(i, WP)*0.1_WP, i=1, 2*(nl-1)*ne)], [2, nl-1, ne])
        a%uv_rhsAB = reshape([(real(i, WP)*0.01_WP, i=1, 2*(nl-1)*ne*2)], [2, nl-1, ne, 2])
        a%uvnode = reshape([(real(i, WP)*0.2_WP, i=1, 2*(nl-1)*nn)], [2, nl-1, nn])
        a%w = reshape([(real(i, WP)*0.3_WP, i=1, nl*nn)], [nl, nn])
        a%eta_n = [(real(i, WP)*0.5_WP, i=1, nn)]
        a%ssh_rhs = [(real(i, WP)*0.7_WP, i=1, nn)]
        a%solverinfo%rr = [(real(i, WP)*1.1_WP, i=1, nn)]
        a%work%uvnode_rhs = a%uvnode * 5.0_WP

        u = scratch_write_open(); write(u, iostat=ios) a; close(u)
        u = scratch_read_open();  read(u,  iostat=ios) b; close(u)

        if (b%opt_visc /= a%opt_visc .or. b%momadv_opt /= a%momadv_opt) call fail('dyn config ints')
        if (b%visc_gamma0 /= a%visc_gamma0) call fail('dyn visc_gamma0')
        if (b%use_wsplit .neqv. a%use_wsplit) call fail('dyn use_wsplit')
        if (b%wsplit_mincfl /= a%wsplit_mincfl) call fail('dyn wsplit_mincfl')
        if (b%wsplit_maxcfl /= a%wsplit_maxcfl) call fail('dyn wsplit_maxcfl')
        if (maxval(abs(b%uv - a%uv)) /= 0.0_WP) call fail('dyn uv')
        if (maxval(abs(b%uv_rhsAB - a%uv_rhsAB)) /= 0.0_WP) call fail('dyn uv_rhsAB (4D)')
        if (maxval(abs(b%uvnode - a%uvnode)) /= 0.0_WP) call fail('dyn uvnode')
        if (maxval(abs(b%w - a%w)) /= 0.0_WP) call fail('dyn w')
        if (maxval(abs(b%eta_n - a%eta_n)) /= 0.0_WP) call fail('dyn eta_n')
        if (maxval(abs(b%solverinfo%rr - a%solverinfo%rr)) /= 0.0_WP) call fail('dyn solverinfo%rr')
        if (maxval(abs(b%work%uvnode_rhs - a%work%uvnode_rhs)) /= 0.0_WP) call fail('dyn work%uvnode_rhs')
    end subroutine

    subroutine test_tracer()
        type(t_tracer) :: a, b
        integer :: u, ios, i
        integer, parameter :: nn = 7, nl = 4
        a%num_tracers = 2
        allocate(a%data(2))
        allocate(a%data(1)%values(nl-1, nn), a%data(2)%values(nl-1, nn))
        allocate(a%data(2)%valuesold(nl-1, nn, 2))
        allocate(a%work%del_ttf(nl-1, nn))
        a%data(1)%ID = 1; a%data(1)%tra_adv_hor = 'MUSCL'; a%data(1)%tra_adv_lim = 'FCT'
        a%data(2)%ID = 2; a%data(2)%tra_adv_ver = 'QR4C'
        a%data(1)%values = reshape([(real(i, WP)*0.4_WP, i=1, (nl-1)*nn)], [nl-1, nn])
        a%data(2)%values = a%data(1)%values * (-1.0_WP)
        a%data(2)%valuesold = reshape([(real(i, WP)*0.05_WP, i=1, (nl-1)*nn*2)], [nl-1, nn, 2])
        a%work%del_ttf = a%data(1)%values + 100.0_WP

        u = scratch_write_open(); write(u, iostat=ios) a; close(u)
        u = scratch_read_open();  read(u,  iostat=ios) b; close(u)

        if (b%num_tracers /= a%num_tracers) call fail('tracer num_tracers')
        if (b%data(1)%ID /= 1 .or. b%data(2)%ID /= 2) call fail('tracer IDs')
        if (trim(b%data(1)%tra_adv_hor) /= 'MUSCL') call fail('tracer tra_adv_hor')
        if (trim(b%data(1)%tra_adv_lim) /= 'FCT')   call fail('tracer tra_adv_lim')
        if (trim(b%data(2)%tra_adv_ver) /= 'QR4C')  call fail('tracer tra_adv_ver')
        if (maxval(abs(b%data(1)%values - a%data(1)%values)) /= 0.0_WP) call fail('tracer data(1)%values')
        if (maxval(abs(b%data(2)%valuesold - a%data(2)%valuesold)) /= 0.0_WP) call fail('tracer data(2)%valuesold (3D)')
        if (maxval(abs(b%work%del_ttf - a%work%del_ttf)) /= 0.0_MP) call fail('tracer work%del_ttf')
    end subroutine

    subroutine test_ice()
        ! M3a: smoke-test the restructured t_ice (data(1:3) ice tracers a_ice/m_ice/m_snow
        ! + work EVP stress/strain on elements). Restart binary I/O for ice is deferred
        ! (cold-start lifecycle only through M3; the M3 byte-gates use mod_dump/FADVHDMP,
        ! not the type's write(unformatted)) — so this is a field-access smoke test.
        type(t_ice) :: a
        integer :: i
        integer, parameter :: nn = 7, ne = 5
        a%ice_dt = 1800.0_WP
        a%num_itracers = 3
        allocate(a%data(3))
        do i = 1, 3
            allocate(a%data(i)%values(nn)); a%data(i)%ID = i
        end do
        allocate(a%uice(nn))
        allocate(a%work%sigma11(ne), a%work%sigma12(ne), a%work%eps11(ne))
        a%data(1)%values = [(real(i, WP)*0.1_WP, i=1, nn)]     ! a_ice
        a%data(2)%values = a%data(1)%values * 2.0_WP            ! m_ice
        a%uice  = [(real(i, WP)*0.02_WP, i=1, nn)]
        a%work%sigma11 = [(real(i, WP)*3.0_WP, i=1, ne)]
        a%work%eps11 = a%work%sigma11 * 0.5_WP

        if (a%ice_dt /= 1800.0_WP) call fail('ice ice_dt')
        if (a%num_itracers /= 3) call fail('ice num_itracers')
        if (size(a%data) /= 3) call fail('ice data size')
        if (a%data(2)%ID /= 2) call fail('ice data ID')
        if (maxval(abs(a%data(2)%values - 2.0_WP*a%data(1)%values)) /= 0.0_WP) call fail('ice m_ice')
        if (maxval(abs(a%work%eps11 - 0.5_WP*a%work%sigma11)) /= 0.0_WP) call fail('ice eps11')
    end subroutine

end program test_types
