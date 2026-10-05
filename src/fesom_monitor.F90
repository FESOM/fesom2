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
    use g_support, only: omp_min_max_sum1, omp_min_max_sum2

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
! and salt on the wet levels. With `full` also the rest:
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
    real(kind=WP)  :: loc, glo, x
    ! every-step entries: 1 eta, 2 deta, 3 wvel1, 4 hnode1, 5 temp, 6 salt
    integer, parameter  :: NS = 6
    character(len=MON_NAME_LEN), parameter :: sname(NS) = &
        [character(len=MON_NAME_LEN) :: 'eta', 'deta', 'wvel1', 'hnode1', 'temp', 'salt']
    real(kind=WP)       :: smin(NS), smax(NS)
    integer(kind=int64) :: snf(NS), gnf(NS)
    real(kind=WP_full)  :: sbuf(3*NS), gbuf(3*NS)
    real(kind=WP), dimension(:,:), pointer :: temp, salt
    real(kind=WP_full) :: wgt
    real(kind=WP_full) :: loc_eta, loc_hbar, loc_dhbar, loc_wflux
    real(kind=WP_full) :: int_eta, int_hbar, int_dhbar, int_wflux
    real(kind=WP)  :: vmin, vmax
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
    ! area integrals of the surface fields: in WP_full and in node order, so
    ! they neither degrade with WP nor depend on the number of threads
    loc_eta   = 0.0_WP_full
    loc_hbar  = 0.0_WP_full
    loc_dhbar = 0.0_WP_full
    loc_wflux = 0.0_WP_full
    do n=1, myDim_nod2D
       wgt       = real(areasvol(ulevels_nod2D(n), n), WP_full)
       loc_eta   = loc_eta   + wgt*real(eta_n(n), WP_full)
       loc_hbar  = loc_hbar  + wgt*real(hbar(n), WP_full)
       loc_dhbar = loc_dhbar + wgt*(real(hbar(n), WP_full)-real(hbar_old(n), WP_full))
       loc_wflux = loc_wflux + wgt*real(water_flux(n), WP_full)
    end do
    call MPI_AllREDUCE(loc_eta  , int_eta  , 1, MPI_WP_FULL, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(loc_hbar , int_hbar , 1, MPI_WP_FULL, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(loc_dhbar, int_dhbar, 1, MPI_WP_FULL, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(loc_wflux, int_wflux, 1, MPI_WP_FULL, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    wgt      = real(ocean_areawithcav, WP_full)
    int_eta  = int_eta  /wgt
    int_hbar = int_hbar /wgt
    int_dhbar= int_dhbar/wgt
    int_wflux= int_wflux/wgt

    !___________________________________________________________________________
    ! 2D node fields
    call set_entry(mon, 'eta', integral=int_eta)
    call minmax1(hbar,       vmin, vmax);  call set_entry(mon, 'hbar',  int_hbar,  vmin, vmax)
    call set_entry(mon, 'dhbar', integral=int_dhbar)
    call minmax1(water_flux, vmin, vmax);  call set_entry(mon, 'wflux', int_wflux, vmin, vmax)
    call minmax1(heat_flux,  vmin, vmax);  call set_entry(mon, 'hflux', vmin=vmin, vmax=vmax)
    call minmax1(Wvel(2,:),     vmin, vmax);  call set_entry(mon, 'wvel2', vmin=vmin, vmax=vmax)
    call minmax1(UVnode(1,1,:), vmin, vmax);  call set_entry(mon, 'uvel1', vmin=vmin, vmax=vmax)
    call minmax1(UVnode(1,2,:), vmin, vmax);  call set_entry(mon, 'uvel2', vmin=vmin, vmax=vmax)
    call minmax1(UVnode(2,1,:), vmin, vmax);  call set_entry(mon, 'vvel1', vmin=vmin, vmax=vmax)
    call minmax1(UVnode(2,2,:), vmin, vmax);  call set_entry(mon, 'vvel2', vmin=vmin, vmax=vmax)
    call minmax1(hnode(2,1:myDim_nod2D), vmin, vmax);  call set_entry(mon, 'hnode2', vmin=vmin, vmax=vmax)
    if (use_ice) then
        loc=omp_min_max_sum1(m_ice, 1, myDim_nod2D, 'max', partit)
        call MPI_AllREDUCE(loc, glo, 1, MPI_WP, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        call set_entry(mon, 'm_ice', vmax=glo)
    end if

    !___________________________________________________________________________
    ! 3D node fields
    if (ldiag_dMOC) then
        call minmax2_masked(density_dmoc, vmin, vmax);        call set_entry(mon, 'dens', vmin=vmin, vmax=vmax)
    end if
    call max2(CFL_z, nl-1, myDim_nod2D,  vmax);  call set_entry(mon, 'cfl_z', vmax=vmax)
    call max2(pgf_x, nl-1, myDim_nod2D,  vmax);  call set_entry(mon, 'pgf_x', vmax=vmax)
    call max2(pgf_y, nl-1, myDim_nod2D,  vmax);  call set_entry(mon, 'pgf_y', vmax=vmax)
    call max2(Av,    nl,   myDim_elem2D, vmax);  call set_entry(mon, 'Av',    vmax=vmax)
    call max2(Kv,    nl,   myDim_nod2D,  vmax);  call set_entry(mon, 'Kv',    vmax=vmax)

contains

    subroutine minmax1(arr, gmin, gmax)
        real(kind=WP), intent(in)  :: arr(:)
        real(kind=WP), intent(out) :: gmin, gmax
        real(kind=WP)              :: l
        l=omp_min_max_sum1(arr, 1, myDim_nod2D, 'min', partit)
        call MPI_AllREDUCE(l, gmin, 1, MPI_WP, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        l=omp_min_max_sum1(arr, 1, myDim_nod2D, 'max', partit)
        call MPI_AllREDUCE(l, gmax, 1, MPI_WP, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    end subroutine minmax1

    subroutine minmax2_masked(arr, gmin, gmax)
        real(kind=WP), intent(in)  :: arr(:,:)
        real(kind=WP), intent(out) :: gmin, gmax
        real(kind=WP)              :: l
        l=omp_min_max_sum2(arr, 1, nl-1, 1, myDim_nod2D, 'min', partit, 0.0_WP)
        call MPI_AllREDUCE(l, gmin, 1, MPI_WP, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        l=omp_min_max_sum2(arr, 1, nl-1, 1, myDim_nod2D, 'max', partit, 0.0_WP)
        call MPI_AllREDUCE(l, gmax, 1, MPI_WP, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    end subroutine minmax2_masked

    subroutine max2(arr, nlev, ncol, gmax)
        real(kind=WP), intent(in)  :: arr(:,:)
        integer,       intent(in)  :: nlev, ncol
        real(kind=WP), intent(out) :: gmax
        real(kind=WP)              :: l
        l=omp_min_max_sum2(arr, 1, nlev, 1, ncol, 'max', partit)
        call MPI_AllREDUCE(l, gmax, 1, MPI_WP, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    end subroutine max2

end subroutine monitor_fill_ocean

end module fesom_monitor_module
