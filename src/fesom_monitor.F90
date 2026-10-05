!===============================================================================
! Monitor record: per component, a set of named global scalars (area integral,
! minimum, maximum, count of non-finite values) of the model state after a time
! step. It is filled once per step in the post-step phase of the driver and read
! by the step log and the blow-up check.
!===============================================================================
module fesom_monitor_module
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use o_PARAM, only: WP, WP_full
    use MOD_MESH
    use MOD_PARTIT
    use MOD_TRACER
    use MOD_DYN
    use MOD_ICE
    use g_config, only: use_ice
    use o_ARRAYS, only: water_flux, heat_flux, pgf_x, pgf_y, Av, Kv, density_dmoc
    use diagnostics, only: ldiag_dMOC

    implicit none
    private
    public :: t_monitor, monitor_fill_ocean, monitor_value, monitor_integral, monitor_nonfinite

    integer, parameter :: MON_NAME_LEN = 16
    integer, parameter :: MON_MAX_ENTRIES = 32

    ! A statistic that is not computed for an entry stays 0.
    type t_mon_entry
        character(len=MON_NAME_LEN) :: name     = ''
        real(kind=WP_full)          :: integral = 0.0_WP_full
        real(kind=WP_full)          :: vmin     = 0.0_WP_full
        real(kind=WP_full)          :: vmax     = 0.0_WP_full
        integer(kind=int64)         :: nonfinite = 0
    end type t_mon_entry

    type t_monitor
        character(len=8)                              :: component = ''
        integer                                       :: step      = -1
        integer                                       :: n         = 0
        type(t_mon_entry), dimension(MON_MAX_ENTRIES) :: e
    end type t_monitor

contains

!
!
!===============================================================================
! Value of statistic `what` ('min' or 'max') of entry `name`, in WP.
function monitor_value(mon, name, what) result(v)
    type(t_monitor), intent(in) :: mon
    character(len=*), intent(in) :: name, what
    real(kind=WP)                :: v
    integer                      :: i

    i = entry_index(mon, name)
    if (i == 0) then
        v = 0.0_WP
        return
    end if
    select case (what)
    case ('min')
        v = real(mon%e(i)%vmin, WP)
    case ('max')
        v = real(mon%e(i)%vmax, WP)
    case default
        v = 0.0_WP
    end select
end function monitor_value
!
!
!===============================================================================
! Integral of entry `name` in WP_full (0 if it has none).
function monitor_integral(mon, name) result(v)
    type(t_monitor), intent(in)  :: mon
    character(len=*), intent(in) :: name
    real(kind=WP_full)           :: v
    integer                      :: i

    i = entry_index(mon, name)
    v = 0.0_WP_full
    if (i > 0) v = mon%e(i)%integral
end function monitor_integral
!
!
!===============================================================================
! Number of non-finite values of entry `name` (0 if it has none recorded).
function monitor_nonfinite(mon, name) result(k)
    type(t_monitor), intent(in)  :: mon
    character(len=*), intent(in) :: name
    integer(kind=int64)          :: k
    integer                      :: i

    i = entry_index(mon, name)
    k = 0
    if (i > 0) k = mon%e(i)%nonfinite
end function monitor_nonfinite
!
!
!===============================================================================
integer function entry_index(mon, name)
    type(t_monitor), intent(in) :: mon
    character(len=*), intent(in) :: name
    integer                      :: i

    entry_index = 0
    do i = 1, mon%n
        if (trim(mon%e(i)%name) == trim(name)) then
            entry_index = i
            return
        end if
    end do
end function entry_index
!
!
!===============================================================================
subroutine set_entry(mon, name, integral, vmin, vmax, nonfinite)
    type(t_monitor), intent(inout) :: mon
    character(len=*), intent(in)   :: name
    real(kind=WP_full), intent(in), optional :: integral
    real(kind=WP), intent(in), optional :: vmin, vmax
    integer(kind=int64), intent(in), optional :: nonfinite
    integer                        :: i

    i = entry_index(mon, name)
    if (i == 0) then
        mon%n = mon%n + 1
        i     = mon%n
        mon%e(i) = t_mon_entry()
        mon%e(i)%name = name
    end if
    if (present(integral)) mon%e(i)%integral = integral
    if (present(vmin))     mon%e(i)%vmin     = real(vmin, WP_full)
    if (present(vmax))     mon%e(i)%vmax     = real(vmax, WP_full)
    if (present(nonfinite)) mon%e(i)%nonfinite = nonfinite
end subroutine set_entry
!
!
!===============================================================================
! Fill the ocean record of step istep.
! Every step, in one pass over the owned nodes: min, max (of finite values) and
! the number of non-finite values of eta, deta, wvel1 and hnode1, and of temp
! and salt on the wet levels. With `full` also the rest, in one more pass:
!   integrals: area-weighted surface means over owned nodes (eta, hbar, dhbar,
!              wflux); 3D tracer fields have none yet.
!   min/max:   over owned nodes (elements for Av); dens over levels 1..nl-1
!              with exact zeros, i.e. dry cells, excluded.
subroutine monitor_fill_ocean(mon, istep, full, ice, dynamics, tracers, partit, mesh)
    type(t_monitor), intent(inout)        :: mon
    integer        , intent(in)           :: istep
    logical        , intent(in)           :: full
    type(t_ice)    , intent(in)   , target :: ice
    type(t_dyn)    , intent(in)   , target :: dynamics
    type(t_tracer) , intent(in)   , target :: tracers
    type(t_partit) , intent(inout), target :: partit
    type(t_mesh)   , intent(in)   , target :: mesh

    integer        :: n, nz
    real(kind=WP)  :: x
    ! every-step entries: 1 eta, 2 deta, 3 wvel1, 4 hnode1, 5 temp, 6 salt
    integer, parameter  :: NS = 6
    character(len=MON_NAME_LEN), parameter :: sname(NS) = &
        [character(len=MON_NAME_LEN) :: 'eta', 'deta', 'wvel1', 'hnode1', 'temp', 'salt']
    real(kind=WP)       :: smin(NS), smax(NS)
    integer(kind=int64) :: snf(NS), gnf(NS)
    real(kind=WP_full)  :: sbuf(3*NS), gbuf(3*NS)
    real(kind=WP), dimension(:,:), pointer :: temp, salt
    real(kind=WP_full) :: wgt
    ! log-only entries (see below)
    integer, parameter  :: NLOG = 16
    real(kind=WP)       :: lmin(NLOG), lmax(NLOG)
    real(kind=WP_full)  :: lbuf(2*NLOG), gbuf_l(2*NLOG), lsum(4), gsum(4)
    real(kind=WP), dimension(:,:,:), pointer :: UVnode
    real(kind=WP), dimension(:,:)  , pointer :: Wvel, CFL_z
    real(kind=WP), dimension(:)    , pointer :: eta_n, d_eta, m_ice
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UVnode => dynamics%uvnode(:,:,:)
    Wvel   => dynamics%w(:,:)
    CFL_z  => dynamics%cfl_z(:,:)
    eta_n  => dynamics%eta_n(:)
    if ( .not. dynamics%use_ssh_se_subcycl) d_eta => dynamics%d_eta(:)
    m_ice  => ice%data(2)%values(:)
    temp   => tracers%data(1)%values(:,:)
    salt   => tracers%data(2)%values(:,:)

    mon%component = 'ocean'
    mon%step      = istep
    mon%n         = 0

    !___________________________________________________________________________
    ! every step: min/max of finite values and non-finite count
    smin = huge(smin)
    smax = -huge(smax)
    snf  = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, x) REDUCTION(min:smin) REDUCTION(max:smax) REDUCTION(+:snf)
    do n=1, myDim_nod2D
        x = eta_n(n)
        if (ieee_is_finite(x)) then
            smin(1) = min(smin(1), x); smax(1) = max(smax(1), x)
        else
            snf(1) = snf(1) + 1
        end if
        if ( .not. dynamics%use_ssh_se_subcycl) then
            x = d_eta(n)
        else
            x = hbar(n)-hbar_old(n)
        end if
        if (ieee_is_finite(x)) then
            smin(2) = min(smin(2), x); smax(2) = max(smax(2), x)
        else
            snf(2) = snf(2) + 1
        end if
        x = Wvel(1,n)
        if (ieee_is_finite(x)) then
            smin(3) = min(smin(3), x); smax(3) = max(smax(3), x)
        else
            snf(3) = snf(3) + 1
        end if
        x = hnode(1,n)
        if (ieee_is_finite(x)) then
            smin(4) = min(smin(4), x); smax(4) = max(smax(4), x)
        else
            snf(4) = snf(4) + 1
        end if
        do nz=ulevels_nod2D(n), nlevels_nod2D(n)-1
            x = temp(nz,n)
            if (ieee_is_finite(x)) then
                smin(5) = min(smin(5), x); smax(5) = max(smax(5), x)
            else
                snf(5) = snf(5) + 1
            end if
            x = salt(nz,n)
            if (ieee_is_finite(x)) then
                smin(6) = min(smin(6), x); smax(6) = max(smax(6), x)
            else
                snf(6) = snf(6) + 1
            end if
        end do
    end do
!$OMP END PARALLEL DO
    ! One reduction per step: min as max of the negated value, and the largest
    ! per-rank non-finite count, all exact in WP_full. The global counts are
    ! summed only when some rank has a non-finite value, which every rank sees.
    sbuf(1:NS)        = -real(smin, WP_full)
    sbuf(NS+1:2*NS)   =  real(smax, WP_full)
    sbuf(2*NS+1:3*NS) =  real(snf,  WP_full)
    call MPI_AllREDUCE(sbuf, gbuf, 3*NS, MPI_WP_FULL, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    if (any(gbuf(2*NS+1:3*NS) > 0.0_WP_full)) then
        call MPI_AllREDUCE(snf, gnf, NS, MPI_INTEGER8, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    else
        gnf = 0
    end if
    do n=1, NS
        call set_entry(mon, trim(sname(n)), vmin=real(-gbuf(n), WP), vmax=real(gbuf(NS+n), WP), &
                       nonfinite=gnf(n))
    end do
    if (.not. full) return

    !___________________________________________________________________________
    ! log-only fields, in one pass over the owned nodes and one over the owned
    ! elements:
    !   min/max: 1 hbar, 2 wflux, 3 hflux, 4 wvel2, 5 uvel1, 6 uvel2, 7 vvel1,
    !            8 vvel2, 9 hnode2, 16 dens (levels 1..nl-1, exact zeros excluded)
    !   max:     10 m_ice, 11 cfl_z, 12 pgf_x, 13 pgf_y (levels 1..nl-1),
    !            14 Kv (levels 1..nl), 15 Av (elements, levels 1..nl)
    lmin = huge(lmin)
    lmax = -huge(lmax)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, x) REDUCTION(min:lmin) REDUCTION(max:lmax)
    do n=1, myDim_nod2D
        lmin(1) = min(lmin(1), hbar(n));          lmax(1) = max(lmax(1), hbar(n))
        lmin(2) = min(lmin(2), water_flux(n));    lmax(2) = max(lmax(2), water_flux(n))
        lmin(3) = min(lmin(3), heat_flux(n));     lmax(3) = max(lmax(3), heat_flux(n))
        lmin(4) = min(lmin(4), Wvel(2,n));        lmax(4) = max(lmax(4), Wvel(2,n))
        lmin(5) = min(lmin(5), UVnode(1,1,n));    lmax(5) = max(lmax(5), UVnode(1,1,n))
        lmin(6) = min(lmin(6), UVnode(1,2,n));    lmax(6) = max(lmax(6), UVnode(1,2,n))
        lmin(7) = min(lmin(7), UVnode(2,1,n));    lmax(7) = max(lmax(7), UVnode(2,1,n))
        lmin(8) = min(lmin(8), UVnode(2,2,n));    lmax(8) = max(lmax(8), UVnode(2,2,n))
        lmin(9) = min(lmin(9), hnode(2,n));       lmax(9) = max(lmax(9), hnode(2,n))
        if (use_ice) lmax(10) = max(lmax(10), m_ice(n))
        do nz=1, nl-1
            lmax(11) = max(lmax(11), CFL_z(nz,n))
            lmax(12) = max(lmax(12), pgf_x(nz,n))
            lmax(13) = max(lmax(13), pgf_y(nz,n))
        end do
        do nz=1, nl
            lmax(14) = max(lmax(14), Kv(nz,n))
        end do
        if (ldiag_dMOC) then
            do nz=1, nl-1
                x = density_dmoc(nz,n)
                if (x /= 0.0_WP) then
                    lmin(16) = min(lmin(16), x); lmax(16) = max(lmax(16), x)
                end if
            end do
        end if
    end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz) REDUCTION(max:lmax)
    do n=1, myDim_elem2D
        do nz=1, nl
            lmax(15) = max(lmax(15), Av(nz,n))
        end do
    end do
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    ! area integrals of the surface fields: in WP_full and in node order, so
    ! they neither degrade with WP nor depend on the number of threads
    lsum = 0.0_WP_full
    do n=1, myDim_nod2D
       wgt     = real(areasvol(ulevels_nod2D(n), n), WP_full)
       lsum(1) = lsum(1) + wgt*real(eta_n(n), WP_full)
       lsum(2) = lsum(2) + wgt*real(hbar(n), WP_full)
       lsum(3) = lsum(3) + wgt*(real(hbar(n), WP_full)-real(hbar_old(n), WP_full))
       lsum(4) = lsum(4) + wgt*real(water_flux(n), WP_full)
    end do

    !___________________________________________________________________________
    ! two reductions: all min/max (min as max of the negated value) and the sums
    lbuf(1:NLOG)        = -real(lmin, WP_full)
    lbuf(NLOG+1:2*NLOG) =  real(lmax, WP_full)
    call MPI_AllREDUCE(lbuf, gbuf_l, 2*NLOG, MPI_WP_FULL, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(lsum, gsum, 4, MPI_WP_FULL, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    gsum = gsum/real(ocean_areawithcav, WP_full)

    call set_entry(mon, 'eta',    integral=gsum(1))
    call set_entry(mon, 'hbar',   gsum(2), lv(1, 'min'), lv(1, 'max'))
    call set_entry(mon, 'dhbar',  integral=gsum(3))
    call set_entry(mon, 'wflux',  gsum(4), lv(2, 'min'), lv(2, 'max'))
    call set_entry(mon, 'hflux',  vmin=lv(3, 'min'), vmax=lv(3, 'max'))
    call set_entry(mon, 'wvel2',  vmin=lv(4, 'min'), vmax=lv(4, 'max'))
    call set_entry(mon, 'uvel1',  vmin=lv(5, 'min'), vmax=lv(5, 'max'))
    call set_entry(mon, 'uvel2',  vmin=lv(6, 'min'), vmax=lv(6, 'max'))
    call set_entry(mon, 'vvel1',  vmin=lv(7, 'min'), vmax=lv(7, 'max'))
    call set_entry(mon, 'vvel2',  vmin=lv(8, 'min'), vmax=lv(8, 'max'))
    call set_entry(mon, 'hnode2', vmin=lv(9, 'min'), vmax=lv(9, 'max'))
    if (use_ice)    call set_entry(mon, 'm_ice', vmax=lv(10, 'max'))
    if (ldiag_dMOC) call set_entry(mon, 'dens',  vmin=lv(16, 'min'), vmax=lv(16, 'max'))
    call set_entry(mon, 'cfl_z',  vmax=lv(11, 'max'))
    call set_entry(mon, 'pgf_x',  vmax=lv(12, 'max'))
    call set_entry(mon, 'pgf_y',  vmax=lv(13, 'max'))
    call set_entry(mon, 'Av',     vmax=lv(15, 'max'))
    call set_entry(mon, 'Kv',     vmax=lv(14, 'max'))

contains

    ! reduced min or max of log entry i; 0 if no value entered it (dens with
    ! every level masked)
    real(kind=WP) function lv(i, what)
        integer,          intent(in) :: i
        character(len=*), intent(in) :: what
        real(kind=WP_full)           :: v
        if (what == 'min') then
            v = -gbuf_l(i)
        else
            v =  gbuf_l(NLOG+i)
        end if
        if (abs(v) == real(huge(1.0_WP), WP_full)) v = 0.0_WP_full
        lv = real(v, WP)
    end function lv

end subroutine monitor_fill_ocean

end module fesom_monitor_module
