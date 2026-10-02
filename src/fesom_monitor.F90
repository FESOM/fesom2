!===============================================================================
! Monitor record: per component, a set of named global scalars (area integral,
! minimum, maximum) of the model state after a time step. It is filled once in
! the post-step phase of the driver and read by the step log.
!===============================================================================
module fesom_monitor_module
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
    public :: t_monitor, monitor_fill_ocean, monitor_value

    integer, parameter :: MON_NAME_LEN = 16
    integer, parameter :: MON_MAX_ENTRIES = 32

    ! A statistic that is not computed for an entry stays 0.
    type t_mon_entry
        character(len=MON_NAME_LEN) :: name     = ''
        real(kind=WP_full)          :: integral = 0.0_WP_full
        real(kind=WP_full)          :: vmin     = 0.0_WP_full
        real(kind=WP_full)          :: vmax     = 0.0_WP_full
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
! Value of statistic `what` ('int', 'min' or 'max') of entry `name`, in WP.
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
    case ('int')
        v = real(mon%e(i)%integral, WP)
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
subroutine set_entry(mon, name, integral, vmin, vmax)
    type(t_monitor), intent(inout) :: mon
    character(len=*), intent(in)   :: name
    real(kind=WP), intent(in), optional :: integral, vmin, vmax
    integer                        :: i

    i = entry_index(mon, name)
    if (i == 0) then
        mon%n = mon%n + 1
        i     = mon%n
        mon%e(i) = t_mon_entry()
        mon%e(i)%name = name
    end if
    if (present(integral)) mon%e(i)%integral = real(integral, WP_full)
    if (present(vmin))     mon%e(i)%vmin     = real(vmin, WP_full)
    if (present(vmax))     mon%e(i)%vmax     = real(vmax, WP_full)
end subroutine set_entry
!
!
!===============================================================================
! Fill the ocean record of step istep.
!   integrals: area-weighted surface means over owned nodes (eta, hbar, dhbar,
!              wflux); 3D tracer fields have none yet.
!   min/max:   over owned nodes (elements for Av); 3D node fields over levels
!              1..nl-1 with exact zeros, i.e. dry cells, excluded.
subroutine monitor_fill_ocean(mon, istep, ice, dynamics, tracers, partit, mesh)
    type(t_monitor), intent(inout)        :: mon
    integer        , intent(in)           :: istep
    type(t_ice)    , intent(in)   , target :: ice
    type(t_dyn)    , intent(in)   , target :: dynamics
    type(t_tracer) , intent(in)   , target :: tracers
    type(t_partit) , intent(inout), target :: partit
    type(t_mesh)   , intent(in)   , target :: mesh

    integer        :: n
    real(kind=WP)  :: loc, glo
    real(kind=WP)  :: loc_eta, loc_hbar, loc_dhbar, loc_wflux
    real(kind=WP)  :: int_eta, int_hbar, int_dhbar, int_wflux
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

    mon%component = 'ocean'
    mon%step      = istep
    mon%n         = 0

    !___________________________________________________________________________
    ! area integrals of the surface fields
    loc_eta   = 0.
    loc_hbar  = 0.
    loc_dhbar = 0.
    loc_wflux = 0.
#if !defined(__openmp_reproducible)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n) REDUCTION(+:loc_eta, loc_hbar, loc_dhbar, loc_wflux)
#endif
    do n=1, myDim_nod2D
       loc_eta   = loc_eta   + areasvol(ulevels_nod2D(n), n)*eta_n(n)
       loc_hbar  = loc_hbar  + areasvol(ulevels_nod2D(n), n)*hbar(n)
       loc_dhbar = loc_dhbar + areasvol(ulevels_nod2D(n), n)*(hbar(n)-hbar_old(n))
       loc_wflux = loc_wflux + areasvol(ulevels_nod2D(n), n)*water_flux(n)
    end do
#if !defined(__openmp_reproducible)
!$OMP END PARALLEL DO
#endif
    call MPI_AllREDUCE(loc_eta  , int_eta  , 1, MPI_WP, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(loc_hbar , int_hbar , 1, MPI_WP, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(loc_dhbar, int_dhbar, 1, MPI_WP, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(loc_wflux, int_wflux, 1, MPI_WP, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    int_eta  = int_eta  /ocean_areawithcav
    int_hbar = int_hbar /ocean_areawithcav
    int_dhbar= int_dhbar/ocean_areawithcav
    int_wflux= int_wflux/ocean_areawithcav

    !___________________________________________________________________________
    ! 2D node fields
    call minmax1(eta_n,      vmin, vmax);  call set_entry(mon, 'eta',   int_eta,   vmin, vmax)
    call minmax1(hbar,       vmin, vmax);  call set_entry(mon, 'hbar',  int_hbar,  vmin, vmax)
    call set_entry(mon, 'dhbar', integral=int_dhbar)
    call minmax1(water_flux, vmin, vmax);  call set_entry(mon, 'wflux', int_wflux, vmin, vmax)
    call minmax1(heat_flux,  vmin, vmax);  call set_entry(mon, 'hflux', vmin=vmin, vmax=vmax)
    if ( .not. dynamics%use_ssh_se_subcycl) then
        call minmax1(d_eta, vmin, vmax)
    else
        call minmax1(hbar-hbar_old, vmin, vmax)
    end if
    call set_entry(mon, 'deta', vmin=vmin, vmax=vmax)
    call minmax1(Wvel(1,:),     vmin, vmax);  call set_entry(mon, 'wvel1', vmin=vmin, vmax=vmax)
    call minmax1(Wvel(2,:),     vmin, vmax);  call set_entry(mon, 'wvel2', vmin=vmin, vmax=vmax)
    call minmax1(UVnode(1,1,:), vmin, vmax);  call set_entry(mon, 'uvel1', vmin=vmin, vmax=vmax)
    call minmax1(UVnode(1,2,:), vmin, vmax);  call set_entry(mon, 'uvel2', vmin=vmin, vmax=vmax)
    call minmax1(UVnode(2,1,:), vmin, vmax);  call set_entry(mon, 'vvel1', vmin=vmin, vmax=vmax)
    call minmax1(UVnode(2,2,:), vmin, vmax);  call set_entry(mon, 'vvel2', vmin=vmin, vmax=vmax)
    call minmax1(hnode(1,1:myDim_nod2D), vmin, vmax);  call set_entry(mon, 'hnode1', vmin=vmin, vmax=vmax)
    call minmax1(hnode(2,1:myDim_nod2D), vmin, vmax);  call set_entry(mon, 'hnode2', vmin=vmin, vmax=vmax)
    if (use_ice) then
        loc=omp_min_max_sum1(m_ice, 1, myDim_nod2D, 'max', partit)
        call MPI_AllREDUCE(loc, glo, 1, MPI_WP, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        call set_entry(mon, 'm_ice', vmax=glo)
    end if

    !___________________________________________________________________________
    ! 3D node fields, dry cells (exact zeros) excluded
    call minmax2_masked(tracers%data(1)%values, vmin, vmax);  call set_entry(mon, 'temp', vmin=vmin, vmax=vmax)
    call minmax2_masked(tracers%data(2)%values, vmin, vmax);  call set_entry(mon, 'salt', vmin=vmin, vmax=vmax)
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
