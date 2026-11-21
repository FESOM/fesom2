! Minimal FESOM types for the horizontal tracer advection dwarf, and the reader of the per-rank dump
! written by FESOM with FESOM_DWARF_DUMP (src/dwarf_dump.F90). Component names match FESOM's
! t_mesh / t_partit so that the reference kernels are verbatim copies.
module o_PARAM
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    integer, parameter :: WP = real64, WP_full = real64
    real(kind=WP), parameter :: r_earth = 6367500.0_WP
    integer :: acc_vl = 64
end module o_PARAM

module MOD_MESH
    use o_PARAM
    implicit none
    type t_mesh
        integer :: nl
        integer,       allocatable :: edges(:,:), edge_tri(:,:), nlevels(:), ulevels(:)
        integer,       allocatable :: nlevels_nod2D(:), ulevels_nod2D(:)
        integer,       allocatable :: nod_in_edge2D_num(:), nod_in_edge2D(:,:), nod_in_edge2D_sgn(:,:)
        real(kind=WP), allocatable :: edge_cross_dxdy(:,:), edge_dxdy(:,:), elem_cos(:), helem(:,:), areasvol(:,:)
    end type t_mesh
end module MOD_MESH

module MOD_PARTIT
    implicit none
    type t_partit
        integer :: myDim_nod2D, eDim_nod2D, myDim_elem2D, eDim_elem2D, eXDim_elem2D, myDim_edge2D, eDim_edge2D
        integer :: mype = 0, npes = 1
    end type t_partit
end module MOD_PARTIT

module dwarf_data
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use o_PARAM
    use MOD_MESH
    use MOD_PARTIT
    implicit none
    type t_rec
        character(len=32) :: name
        integer :: typ, rnk, shp(3)
        integer(int32), allocatable :: iv(:)
        real(real64),   allocatable :: rv(:)
    end type t_rec
    type t_case
        type(t_mesh)   :: mesh
        type(t_partit) :: partit
        real(kind=WP)  :: dt, opth
        real(kind=WP), allocatable :: ttf(:,:), ttfAB(:,:), vel(:,:,:), grad(:,:,:)
        real(kind=WP), allocatable :: flux_lo_ref(:,:), lo_tend_ref(:,:), flux_ho_ref(:,:)
    end type t_case
contains
    subroutine read_dump(fname, c)
        character(len=*), intent(in) :: fname
        type(t_case), intent(out) :: c
        type(t_rec), allocatable :: r(:)
        integer :: u, n, ios, d(8)
        type(t_rec) :: t
        allocate(r(64)); n = 0
        open(newunit=u, file=fname, access='stream', form='unformatted', status='old', action='read')
        do
            read(u, iostat=ios) t%name, t%typ, t%rnk
            if (ios /= 0) exit
            t%shp = 1
            if (t%rnk > 0) read(u) t%shp(1:t%rnk)
            if (t%typ == 1) then
                allocate(t%iv(product(t%shp(1:max(1,t%rnk))))); read(u) t%iv
            else
                if (t%rnk == 0) then; allocate(t%rv(1)); else; allocate(t%rv(product(t%shp(1:t%rnk)))); end if
                read(u) t%rv
            end if
            n = n + 1; r(n) = t
            if (allocated(t%iv)) deallocate(t%iv)
            if (allocated(t%rv)) deallocate(t%rv)
        end do
        close(u)
        d = geti1('dims')
        c%mesh%nl = d(1)
        c%partit%myDim_nod2D = d(2); c%partit%eDim_nod2D = d(3); c%partit%myDim_elem2D = d(4)
        c%partit%eDim_elem2D = d(5); c%partit%eXDim_elem2D = d(6); c%partit%myDim_edge2D = d(7)
        c%partit%eDim_edge2D = d(8)
        if (abs(getr0('r_earth') - r_earth) > 0) stop 'r_earth differs'
        c%dt = getr0('dt'); c%opth = getr0('opth')
        c%mesh%edges = geti2('edges');           c%mesh%edge_tri = geti2('edge_tri')
        c%mesh%nlevels = geti1('nlevels');       c%mesh%ulevels = geti1('ulevels')
        c%mesh%nlevels_nod2D = geti1('nlevels_nod2D'); c%mesh%ulevels_nod2D = geti1('ulevels_nod2D')
        c%mesh%nod_in_edge2D_num = geti1('nod_in_edge2D_num')
        c%mesh%nod_in_edge2D = geti2('nod_in_edge2D'); c%mesh%nod_in_edge2D_sgn = geti2('nod_in_edge2D_sgn')
        c%mesh%edge_cross_dxdy = getr2('edge_cross_dxdy'); c%mesh%edge_dxdy = getr2('edge_dxdy')
        c%mesh%elem_cos = getr1('elem_cos'); c%mesh%helem = getr2('helem'); c%mesh%areasvol = getr2('areasvol')
        c%ttf = getr2('ttf'); c%ttfAB = getr2('ttfAB'); c%vel = getr3('vel'); c%grad = getr3('edge_up_dn_grad')
        c%flux_lo_ref = getr2('flux_lo'); c%lo_tend_ref = getr2('lo_tend'); c%flux_ho_ref = getr2('flux_ho')
    contains
        integer function find(nm)
            character(len=*), intent(in) :: nm
            do find = 1, n
                if (trim(r(find)%name) == nm) return
            end do
            write(*,*) 'missing record ', nm; stop 1
        end function find
        function geti1(nm) result(a)
            character(len=*), intent(in) :: nm
            integer, allocatable :: a(:)
            integer :: k
            k = find(nm); a = r(k)%iv
        end function geti1
        function geti2(nm) result(a)
            character(len=*), intent(in) :: nm
            integer, allocatable :: a(:,:)
            integer :: k
            k = find(nm); a = reshape(r(k)%iv, r(k)%shp(1:2))
        end function geti2
        real(real64) function getr0(nm)
            character(len=*), intent(in) :: nm
            integer :: k
            k = find(nm); getr0 = r(k)%rv(1)
        end function getr0
        function getr1(nm) result(a)
            character(len=*), intent(in) :: nm
            real(real64), allocatable :: a(:)
            integer :: k
            k = find(nm); a = r(k)%rv
        end function getr1
        function getr2(nm) result(a)
            character(len=*), intent(in) :: nm
            real(real64), allocatable :: a(:,:)
            integer :: k
            k = find(nm); a = reshape(r(k)%rv, r(k)%shp(1:2))
        end function getr2
        function getr3(nm) result(a)
            character(len=*), intent(in) :: nm
            real(real64), allocatable :: a(:,:,:)
            integer :: k
            k = find(nm); a = reshape(r(k)%rv, r(k)%shp(1:3))
        end function getr3
    end subroutine read_dump
end module dwarf_data
