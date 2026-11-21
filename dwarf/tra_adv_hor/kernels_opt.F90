! Variants of the horizontal advection kernels.
! - edge_meta: per-edge data that only depends on the mesh (computed once, as it could be at model
!   start): level range of each adjacent element, the cross vectors, and a = r_earth*mean(cos).
! - *_merged: one level loop per edge from the top of either element to the bottom of either,
!   with 0/1 weights for the element contributions, instead of up to five loops (A-E). Longer,
!   fewer loops; the weights remove the branches. Boundary edges keep the one-sided loop.
! - fused_upw1_mfct: one edge pass computes the low-order flux (ttf) and the high-order flux
!   (ttfAB, as flux_ho - flux_lo, the antidiffusive flux FESOM's MFCT returns with FCT on),
!   so velocity, layer thickness and geometry are read once instead of twice.
module kernels_opt
    use o_PARAM
    use MOD_MESH
    use MOD_PARTIT
    implicit none
    type t_edge_meta
        integer,       allocatable :: n1(:), n2(:), e1(:), e2(:)      ! nodes, elements (e2 = e1 if boundary)
        integer,       allocatable :: nu1(:), nl1(:), nu2(:), nl2(:)  ! level ranges (nl2 = 0 if boundary)
        integer,       allocatable :: nzt(:), nzb(:)                  ! top and bottom over both elements
        real(kind=WP), allocatable :: dx1(:), dy1(:), dx2(:), dy2(:), a(:), ex(:), ey(:)
    end type t_edge_meta
contains

    subroutine edge_meta_init(m, partit, mesh)
        type(t_edge_meta), intent(out) :: m
        type(t_partit), intent(in) :: partit
        type(t_mesh),   intent(in) :: mesh
        integer :: e, ne
        ne = partit%myDim_edge2D
        allocate(m%n1(ne), m%n2(ne), m%e1(ne), m%e2(ne), m%nu1(ne), m%nl1(ne), m%nu2(ne), m%nl2(ne), &
                 m%nzt(ne), m%nzb(ne), m%dx1(ne), m%dy1(ne), m%dx2(ne), m%dy2(ne), m%a(ne), m%ex(ne), m%ey(ne))
        do e = 1, ne
            m%n1(e) = mesh%edges(1,e); m%n2(e) = mesh%edges(2,e)
            m%e1(e) = mesh%edge_tri(1,e)
            m%nu1(e) = mesh%ulevels(m%e1(e)); m%nl1(e) = mesh%nlevels(m%e1(e))-1
            m%dx1(e) = mesh%edge_cross_dxdy(1,e); m%dy1(e) = mesh%edge_cross_dxdy(2,e)
            m%a(e) = r_earth*mesh%elem_cos(m%e1(e))
            if (mesh%edge_tri(2,e) > 0) then
                m%e2(e) = mesh%edge_tri(2,e)
                m%nu2(e) = mesh%ulevels(m%e2(e)); m%nl2(e) = mesh%nlevels(m%e2(e))-1
                m%dx2(e) = mesh%edge_cross_dxdy(3,e); m%dy2(e) = mesh%edge_cross_dxdy(4,e)
                m%a(e) = 0.5_WP*(m%a(e)+r_earth*mesh%elem_cos(m%e2(e)))
            else
                m%e2(e) = m%e1(e); m%nu2(e) = 0; m%nl2(e) = 0
                m%dx2(e) = 0; m%dy2(e) = 0
            end if
            m%nzt(e) = m%nu1(e); if (m%nu2(e) > 0) m%nzt(e) = min(m%nu1(e), m%nu2(e))
            m%nzb(e) = max(m%nl1(e), m%nl2(e))
            m%ex(e) = mesh%edge_dxdy(1,e)*m%a(e)
            m%ey(e) = mesh%edge_dxdy(2,e)*r_earth
        end do
    end subroutine edge_meta_init

    ! low-order flux, one merged level loop per edge (o_init_zero = .true.)
    subroutine upw1_merged(m, vel, ttf, helem, flux, nl, ne, nn, nel)
        type(t_edge_meta), intent(in) :: m
        integer, intent(in) :: nl, ne, nn, nel
        real(kind=WP), intent(in)    :: vel(2, nl-1, nel), ttf(nl-1, nn), helem(nl-1, nel)
        real(kind=WP), intent(inout) :: flux(nl-1, ne)
        integer :: e, nz, e1, e2, n1, n2, nu1, nl1, nu2, nl2, nzt, nzb
        real(kind=WP) :: w1, w2, vflux, dx1, dy1, dx2, dy2
        do e = 1, ne
            e1 = m%e1(e); e2 = m%e2(e); n1 = m%n1(e); n2 = m%n2(e)
            nu1 = m%nu1(e); nl1 = m%nl1(e); nu2 = m%nu2(e); nl2 = m%nl2(e)
            if (nu2 == 0) nu2 = huge(nu2)
            nzt = m%nzt(e); nzb = m%nzb(e)
            dx1 = m%dx1(e); dy1 = m%dy1(e); dx2 = m%dx2(e); dy2 = m%dy2(e)
            flux(:, e) = 0.0_WP
            do nz = nzt, nzb
                w1 = merge(1.0_WP, 0.0_WP, nz >= nu1 .and. nz <= nl1)
                w2 = merge(1.0_WP, 0.0_WP, nz >= nu2 .and. nz <= nl2)
                vflux = w1*(-vel(2,nz,e1)*dx1 + vel(1,nz,e1)*dy1)*helem(nz,e1) &
                      + w2*( vel(2,nz,e2)*dx2 - vel(1,nz,e2)*dy2)*helem(nz,e2)
                flux(nz, e) = -0.5_WP*(ttf(nz,n1)*(vflux+abs(vflux)) + ttf(nz,n2)*(vflux-abs(vflux)))
            end do
        end do
    end subroutine upw1_merged

    ! high-order flux minus the input flux (o_init_zero = .false.), one merged level loop per edge
    subroutine mfct_merged(m, num_ord, vel, ttf, helem, grad, flux, nl, ne, nn, nel)
        type(t_edge_meta), intent(in) :: m
        integer, intent(in) :: nl, ne, nn, nel
        real(kind=WP), intent(in)    :: num_ord
        real(kind=WP), intent(in)    :: vel(2, nl-1, nel), ttf(nl-1, nn), helem(nl-1, nel), grad(4, nl-1, ne)
        real(kind=WP), intent(inout) :: flux(nl-1, ne)
        integer :: e, nz, e1, e2, n1, n2, nu1, nl1, nu2, nl2, nzt, nzb
        real(kind=WP) :: w1, w2, vflux, t1, t2, tm1, tm2, cHO, dx1, dy1, dx2, dy2, ex, ey
        do e = 1, ne
            e1 = m%e1(e); e2 = m%e2(e); n1 = m%n1(e); n2 = m%n2(e)
            nu1 = m%nu1(e); nl1 = m%nl1(e); nu2 = m%nu2(e); nl2 = m%nl2(e)
            if (nu2 == 0) nu2 = huge(nu2)
            nzt = m%nzt(e); nzb = m%nzb(e)
            dx1 = m%dx1(e); dy1 = m%dy1(e); dx2 = m%dx2(e); dy2 = m%dy2(e)
            ex = m%ex(e); ey = m%ey(e)
            do nz = nzt, nzb
                w1 = merge(1.0_WP, 0.0_WP, nz >= nu1 .and. nz <= nl1)
                w2 = merge(1.0_WP, 0.0_WP, nz >= nu2 .and. nz <= nl2)
                t1 = ttf(nz,n1); t2 = ttf(nz,n2)
                tm2 = t2 - (2.0_WP*(t2-t1) + ex*grad(2,nz,e) + ey*grad(4,nz,e))/6.0_WP
                tm1 = t1 + (2.0_WP*(t2-t1) + ex*grad(1,nz,e) + ey*grad(3,nz,e))/6.0_WP
                vflux = w1*(-vel(2,nz,e1)*dx1 + vel(1,nz,e1)*dy1)*helem(nz,e1) &
                      + w2*( vel(2,nz,e2)*dx2 - vel(1,nz,e2)*dy2)*helem(nz,e2)
                cHO = (vflux+abs(vflux))*tm1 + (vflux-abs(vflux))*tm2
                flux(nz,e) = -0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(tm1+tm2) - flux(nz,e)
            end do
        end do
    end subroutine mfct_merged

    ! both fluxes in one edge pass: flux_lo (from ttf) and flux_ho - flux_lo (from ttfAB)
    subroutine fused_upw1_mfct(m, num_ord, vel, ttf, ttfab, helem, grad, flux_lo, flux, nl, ne, nn, nel)
        type(t_edge_meta), intent(in) :: m
        integer, intent(in) :: nl, ne, nn, nel
        real(kind=WP), intent(in)  :: num_ord
        real(kind=WP), intent(in)  :: vel(2, nl-1, nel), ttf(nl-1, nn), ttfab(nl-1, nn), helem(nl-1, nel)
        real(kind=WP), intent(in)  :: grad(4, nl-1, ne)
        real(kind=WP), intent(out) :: flux_lo(nl-1, ne), flux(nl-1, ne)
        integer :: e, nz, e1, e2, n1, n2, nu1, nl1, nu2, nl2, nzt, nzb
        real(kind=WP) :: w1, w2, vflux, t1, t2, tm1, tm2, cHO, flo, dx1, dy1, dx2, dy2, ex, ey
        do e = 1, ne
            e1 = m%e1(e); e2 = m%e2(e); n1 = m%n1(e); n2 = m%n2(e)
            nu1 = m%nu1(e); nl1 = m%nl1(e); nu2 = m%nu2(e); nl2 = m%nl2(e)
            if (nu2 == 0) nu2 = huge(nu2)
            nzt = m%nzt(e); nzb = m%nzb(e)
            dx1 = m%dx1(e); dy1 = m%dy1(e); dx2 = m%dx2(e); dy2 = m%dy2(e)
            ex = m%ex(e); ey = m%ey(e)
            flux_lo(:, e) = 0.0_WP
            flux(:, e) = 0.0_WP
            do nz = nzt, nzb
                w1 = merge(1.0_WP, 0.0_WP, nz >= nu1 .and. nz <= nl1)
                w2 = merge(1.0_WP, 0.0_WP, nz >= nu2 .and. nz <= nl2)
                vflux = w1*(-vel(2,nz,e1)*dx1 + vel(1,nz,e1)*dy1)*helem(nz,e1) &
                      + w2*( vel(2,nz,e2)*dx2 - vel(1,nz,e2)*dy2)*helem(nz,e2)
                flo = -0.5_WP*(ttf(nz,n1)*(vflux+abs(vflux)) + ttf(nz,n2)*(vflux-abs(vflux)))
                flux_lo(nz, e) = flo
                t1 = ttfab(nz,n1); t2 = ttfab(nz,n2)
                tm2 = t2 - (2.0_WP*(t2-t1) + ex*grad(2,nz,e) + ey*grad(4,nz,e))/6.0_WP
                tm1 = t1 + (2.0_WP*(t2-t1) + ex*grad(1,nz,e) + ey*grad(3,nz,e))/6.0_WP
                cHO = (vflux+abs(vflux))*tm1 + (vflux-abs(vflux))*tm2
                flux(nz,e) = -0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(tm1+tm2) - flo
            end do
        end do
    end subroutine fused_upw1_mfct

    ! node gather of edge fluxes with precomputed per-node edge level ranges
    subroutine gather_meta(m, flux, sgn_k, nin_num, nin, out, scale, nl, ne, nn)
        type(t_edge_meta), intent(in) :: m
        integer, intent(in) :: nl, ne, nn, nin_num(nn), nin(:,:), sgn_k(:,:)
        real(kind=WP), intent(in)    :: flux(nl-1, ne)
        real(kind=WP), intent(in), optional :: scale(:,:)
        real(kind=WP_full), intent(inout) :: out(nl-1, nn)
        integer :: n, k, e, nz, nzt, nzb
        if (present(scale)) then
            do n = 1, nn
                do k = 1, nin_num(n)
                    e = nin(k,n); nzt = m%nzt(e); nzb = m%nzb(e)
                    do nz = nzt, nzb
                        out(nz,n) = out(nz,n) + sgn_k(k,n)*real(flux(nz,e), WP_full)*scale(nz,n)
                    end do
                end do
            end do
        else
            do n = 1, nn
                do k = 1, nin_num(n)
                    e = nin(k,n); nzt = m%nzt(e); nzb = m%nzb(e)
                    do nz = nzt, nzb
                        out(nz,n) = out(nz,n) + sgn_k(k,n)*real(flux(nz,e), WP_full)
                    end do
                end do
            end do
        end if
    end subroutine gather_meta
end module kernels_opt
