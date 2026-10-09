module oce_initial_state
    ! 3-D climatology initial conditions for T/S (M2.11b do_ic3d). Transcribed from
    ! FESOM2 v2.7.3 gen_ic3d.F90 (MODULE g_ic3d) + the extrapolation helper
    ! gen_support.F90:400 (extrap_nod3D):
    !   nc_readGrid   ( 68-217)  read lon/lat/depth axes; +2 periodic-lon halo; cyclic
    !   nc_ic3d_ini   (220-301)  per-node bilinear source indices (non-cavity branch)
    !   getcoeffld    (303-491)  read the 3-D double cube; NaN->dummy; spatial bilinear
    !                            -> data1d; vertical LINEAR interp onto Z_3d_n
    !   extrap_nod3D  (gs:400)   fill dummy: Gauss-Seidel horizontal neighbour-average
    !                            sweep (to convergence) + downward vertical fill
    !   do_ic3d       (493-644)  orchestrate per tracer; zero land/bottom; Kelvin guard;
    !                            insitu2pot (t_insitu=.true.)
    !
    ! CONFIG (CORE2 work_core/namelist.tra): n_ic3d=2, idlist=2,1, filelist=2x
    ! phc3.0_winter.nc, varlist='salt','temp', t_insitu=.true. So salt (ID 2) is read
    ! FIRST into tracers%data(2)%values, temp (ID 1) SECOND into data(1)%values; then
    ! data(1) is converted in-situ->potential. phc3.0_winter.nc: depth=33 (0..5500 m,
    ! positive down), lat=180 (-89.5..89.5 ascending), lon=360 (0.5..359.5).
    !
    ! BIT-IDENTITY NOTES (gate target = post-init data(1)/data(2)%values vs the FESOM2
    ! oracle's Tclim/Sclim, max|delta|=0 on CORE2 1-rank):
    ! * The IC data is read as real(8) (FESOM2 nf_get_vara_double), promoted from the
    !   on-disk double by netCDF identically on both sides; phc3.0 land points are NaN
    !   (no _FillValue, data in [-2.1,41.4]) so the missing-value mask reduces to
    !   ieee_is_nan(v) .or. v<-0.99*dummy .or. v>dummy (the FESOM2 ==FILL_VALUE branch
    !   is subsumed by >dummy for NF_FILL_DOUBLE, vacuous here).
    ! * The bilinear weights / vertical-interp coefficients are per-node-independent WP
    !   arithmetic on byte-identical operands (geo_coord_nod2D/rad geometry-proven on
    !   CORE2 by M2.11a; Z_3d_n built linfs full-cell = reference mid-depth at init;
    !   nc_lon/nc_lat/nc_depth read identically). binarysearch is reused verbatim from
    !   mod_forcing_read (forcing_binarysearch, identical d=1e-9 bisection).
    ! * extrap_nod3D is THE order-dependent step: its Gauss-Seidel sweep reads/writes a
    !   work_array in NODE order, accumulating neighbours via nod_in_elem2D/elem2D_nodes
    !   in their stored order. At 1-rank (exchange_nod a no-op, eDim=0) the node order +
    !   nod_in_elem2D order match FESOM2's global order (geometry-gate proven, L9), so
    !   the sweep is deterministic and byte-identical.
    !
    ! SCOPE: no cavity (use_cavity=.false. on CORE2), no RECOM, T/S only. M2.12-MVP added
    ! the MULTI-RANK path via an OPTIONAL partit (absent => the proven 1-rank path verbatim;
    ! present+npes>1 => owned / owned+halo loop bounds + the FESOM2 extrap_nod3D exchange_nod
    ! + allreduce_max(glob_max), gated per-rank vs same-partition FESOM2, L8). The cavity
    ! bilinear/vertical branches stay v1-dropped (deferred to a cavity re-gate). do_ic3d takes
    ! a t_ic3d_config the driver fills from the namelist (no namelist reader in v1).
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_precision,    only: WP
    use mod_constants,    only: rad
    use mod_mesh,         only: t_mesh, MAX_NV
    use mod_partit,       only: t_partit
    use mod_tracer,       only: t_tracer
    use mod_io_netcdf,    only: nc_open_read, nc_close, nc_dimlen, nc_get_axis_dp, &
                                nc_get_var3d_dp
    use mod_forcing_read, only: forcing_binarysearch
    use mod_part_bounds,  only: owned_bounds, is_multirank
    use mod_halo,         only: exchange_nod, allreduce_max
    use oce_pressure_bv,  only: insitu2pot
    implicit none
    private
    public :: t_ic3d_config, do_ic3d
    public :: extrap_nod3D   ! also exercised by test_pgf_mixed (E1)

    integer, parameter :: IC_MAX = 10

    type :: t_ic3d_config
        integer            :: n_ic3d   = 0
        integer            :: idlist(IC_MAX)   = 0
        character(len=512) :: filelist(IC_MAX) = ''   ! full paths (ClimateDataPath prepended)
        character(len=64)  :: varlist(IC_MAX)  = ''
        logical            :: t_insitu  = .true.
        logical            :: ic_cyclic = .true.
        real(WP)           :: dummy     = 1.e10_WP    ! FESOM2 g_config default
    end type t_ic3d_config

    ! Module-saved netCDF grid + per-node bilinear indices (mirrors g_ic3d module
    ! state): allocated in nc_readGrid/do_ic3d, freed in nc_end/do_ic3d.
    real(WP), allocatable, save :: nc_lon(:), nc_lat(:), nc_depth(:)
    integer,  save :: nc_Nlon, nc_Nlat, nc_Ndepth
    integer,  allocatable, save :: bilin_indx_i(:), bilin_indx_j(:)

contains

    !===========================================================================
    ! Read lon/lat/depth axes; +2 periodic-lon halo; cyclic +-360 shift. FESOM2
    ! nc_readGrid (68-217). The lon halo mirrors the wrap columns; depth is forced
    ! positive (nc_depth(2)<0 -> negate); lat is NOT flipped (the IC reader, unlike
    ! the surface forcing reader, keeps the file order).
    subroutine nc_readGrid(filename, ic)
        character(len=*),     intent(in) :: filename
        type(t_ic3d_config),  intent(in) :: ic
        integer :: ncid, nlon0
        real(real64), allocatable :: tmp(:)
        ncid   = nc_open_read(filename)
        nlon0  = nc_dimlen(ncid, ['LON      ','lon      ','longitude'])
        nc_Nlat   = nc_dimlen(ncid, ['LAT      ','lat      ','latitude '])
        nc_Ndepth = nc_dimlen(ncid, ['depth'])
        nc_Nlon   = nlon0 + 2                       ! +2 periodic-lon halo
        allocate(nc_lon(nc_Nlon), nc_lat(nc_Nlat), nc_depth(nc_Ndepth))
        ! lat (no halo)
        call nc_get_axis_dp(ncid, ['LAT      ','lat      ','latitude '], nc_lat)
        ! lon into interior 2:Nlon-1, then periodic halo
        allocate(tmp(nlon0))
        call nc_get_axis_dp(ncid, ['LON      ','lon      ','longitude'], tmp)
        nc_lon(2:nc_Nlon-1) = tmp
        nc_lon(1)       = nc_lon(nc_Nlon-1)
        nc_lon(nc_Nlon) = nc_lon(2)
        deallocate(tmp)
        ! depth (positive down)
        call nc_get_axis_dp(ncid, ['depth'], nc_depth)
        if (nc_depth(2) < 0.0_WP) nc_depth = -nc_depth
        call nc_close(ncid)
        ! cyclic lon: shift the two halo columns by +-360
        if (ic%ic_cyclic) then
            nc_lon(1)       = nc_lon(1)       - 360._WP
            nc_lon(nc_Nlon) = nc_lon(nc_Nlon) + 360._WP
        end if
    end subroutine nc_readGrid

    !===========================================================================
    ! Per-node bilinear source indices (FESOM2 nc_ic3d_ini, non-cavity branch
    ! 281-300). Reads the grid first (nc_readGrid is called from here, as in FESOM2).
    subroutine nc_ic3d_ini(filename, ic, mesh, partit)
        character(len=*),    intent(in) :: filename
        type(t_ic3d_config), intent(in) :: ic
        type(t_mesh),        intent(in) :: mesh
        type(t_partit),      intent(in), optional :: partit
        integer  :: i, nNodO, nNodL, nEdgeO, nElemO
        real(WP) :: x, y
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        call nc_readGrid(filename, ic)
        do i = 1, nNodO            ! FESOM2 nc_ic3d_ini: do i=1, myDim_nod2D (owned)
            x = mesh%geo_coord_nod2D(1, i)/rad
            y = mesh%geo_coord_nod2D(2, i)/rad
            if (x < 0._WP)   x = x + 360._WP
            if (x > 360._WP) x = x - 360._WP
            if (x <= nc_lon(nc_Nlon) .and. x >= nc_lon(1)) then
                call forcing_binarysearch(nc_Nlon, nc_lon, x, bilin_indx_i(i))
            else                                    ! NO extrapolation in space
                bilin_indx_i(i) = -1
            end if
            if (y <= nc_lat(nc_Nlat) .and. y >= nc_lat(1)) then
                call forcing_binarysearch(nc_Nlat, nc_lat, y, bilin_indx_j(i))
            else
                bilin_indx_j(i) = -1
            end if
        end do
    end subroutine nc_ic3d_ini

    !===========================================================================
    ! Read the 3-D double cube, mask missing -> dummy, spatial bilinear + vertical
    ! linear interp onto Z_3d_n. FESOM2 getcoeffld (303-491, non-cavity branch).
    subroutine getcoeffld(filename, varname, ic, values, mesh, partit)
        character(len=*),    intent(in)    :: filename, varname
        type(t_ic3d_config), intent(in)    :: ic
        type(t_mesh),        intent(in)    :: mesh
        real(WP),            intent(inout) :: values(:,:)   ! (nl-1, nNodL) — assumed-shape so
                                                            ! values(:,:)=dummy spans the LOCAL size
        type(t_partit),      intent(in), optional :: partit
        integer  :: ncid, i, j, ii, ip1, jp1, k, d_indx, d_indx_p1, nl1, ul1
        integer  :: nNodO, nNodL, nEdgeO, nElemO
        real(WP) :: cf_a, cf_b, delta_d, denom, x1, x2, y1, y2, x, y, d1, d2
        real(WP) :: dummy
        real(real64), allocatable :: raw(:,:,:)
        real(WP),     allocatable :: ncdata(:,:,:), data1d(:)

        dummy = ic%dummy
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        allocate(ncdata(nc_Nlon, nc_Nlat, nc_Ndepth), data1d(nc_Ndepth))
        allocate(raw(nc_Nlon-2, nc_Nlat, nc_Ndepth))
        ncdata = 0.0_WP
        data1d = 0.0_WP
        values(:,:) = dummy
        ncid = nc_open_read(filename)
        call nc_get_var3d_dp(ncid, [varname], raw)
        call nc_close(ncid)
        ncdata(2:nc_Nlon-1, :, :) = raw
        ncdata(1,       :, :) = ncdata(nc_Nlon-1, :, :)        ! periodic lon halo
        ncdata(nc_Nlon, :, :) = ncdata(2,         :, :)
        ! replace NaN / out-of-range by dummy (phc3.0 land = NaN; no _FillValue)
        do k = 1, nc_Ndepth
            do j = 1, nc_Nlat
                do i = 1, nc_Nlon
                    if (ieee_is_nan(ncdata(i,j,k)) .or. &
                        ncdata(i,j,k) < -0.99_WP*dummy .or. ncdata(i,j,k) > dummy) then
                        ncdata(i,j,k) = dummy
                    end if
                end do
            end do
        end do
        deallocate(raw)

        ! bilinear space interp + vertical linear interp (data on a regular grid).
        ! Owned nodes only (FESOM2 getcoeffld: do ii=1, myDim_nod2D); the halo is filled
        ! by extrap_nod3D's exchange_nod (per-node independent => halo == owner's value).
        do ii = 1, nNodO
            nl1 = mesh%nlevels_nod2D(ii) - 1
            ul1 = mesh%ulevels_nod2D(ii)
            i   = bilin_indx_i(ii)
            j   = bilin_indx_j(ii)
            ip1 = i + 1
            jp1 = j + 1
            x = mesh%geo_coord_nod2D(1, ii)/rad
            y = mesh%geo_coord_nod2D(2, ii)/rad
            if (x < 0._WP)   x = x + 360._WP
            if (x > 360._WP) x = x - 360._WP
            if (min(i,j) > 0) then
                if (any(ncdata(i:ip1, j:jp1, 1) > dummy*0.99_WP)) cycle
                x1 = nc_lon(i);  x2 = nc_lon(ip1)
                y1 = nc_lat(j);  y2 = nc_lat(jp1)
                ! if point inside forcing domain
                denom = (x2 - x1)*(y2 - y1)
                data1d(:) = ( ncdata(i,j,:)   * (x2-x)*(y2-y) + ncdata(ip1,j,:)    * (x-x1)*(y2-y) + &
                              ncdata(i,jp1,:) * (x2-x)*(y-y1) + ncdata(ip1,jp1,:)  * (x-x1)*(y-y1) ) / denom
                where (ncdata(i,j,:)   > 0.99_WP*dummy .OR. ncdata(ip1,j,:)   > 0.99_WP*dummy .OR. &
                       ncdata(i,jp1,:) > 0.99_WP*dummy .OR. ncdata(ip1,jp1,:) > 0.99_WP*dummy)
                    data1d(:) = dummy
                end where
                ! vertical interp (non-cavity)
                do k = ul1, nl1
                    call forcing_binarysearch(nc_Ndepth, nc_depth, -mesh%Z_3d_n(k,ii), d_indx)
                    if (d_indx < nc_Ndepth .and. d_indx > 0) then
                        d_indx_p1 = d_indx + 1
                        delta_d   = nc_depth(d_indx+1) - nc_depth(d_indx)
                        d1 = data1d(d_indx)
                        d2 = data1d(d_indx_p1)
                        if ((d1 < 0.99_WP*dummy) .and. (d2 < 0.99_WP*dummy)) then
                            cf_a = (d2 - d1)/delta_d
                            cf_b = d1 - cf_a*nc_depth(d_indx)
                            values(k,ii) = -cf_a*mesh%Z_3d_n(k,ii) + cf_b
                        end if
                    elseif (d_indx == 0) then
                        values(k,ii) = data1d(1)
                    end if
                end do
            end if
        end do
        deallocate(ncdata, data1d)
    end subroutine getcoeffld

    !===========================================================================
    ! Fill remaining dummy values: horizontal Gauss-Seidel neighbour-average sweep
    ! (to convergence over the surface) + downward vertical fill. FESOM2
    ! extrap_nod3D (gen_support.F90:400-507). 1-rank: exchange_nod is a no-op and the
    ! MPI_AllREDUCE collapses to the local maxval, so they are dropped.
    subroutine extrap_nod3D(arr, ic, mesh, partit)
        type(t_ic3d_config), intent(in)    :: ic
        type(t_mesh),        intent(in)    :: mesh
        real(WP),            intent(inout) :: arr(:,:)   ! (nl-1, nNodL) assumed-shape (FESOM2 form)
        type(t_partit),      intent(in), optional :: partit
        integer  :: n, nl1, nz, k, j, el, cnt
        integer  :: nNodO, nNodL, nEdgeO, nElemO
        integer  :: enodes(MAX_NV)
        logical  :: success
        real(WP) :: val, glob_max, dummy
        real(WP), allocatable :: work_array(:)
        dummy = ic%dummy
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        allocate(work_array(nNodL))
        ! FESOM2 extrap_nod3D (gen_support.F90:418-424): initial halo exchange so the
        ! owned-node Gauss-Seidel sees valid halo neighbours; glob_max is the GLOBAL
        ! surface max (it only drives the loop COUNT). 1-rank: both are no-ops. An owned
        ! node's nod_in_elem2D references only OWNED elements (build_nod_in_elem_local),
        ! so elem2D_nodes(:,el) is in-bounds; halo node VALUES come from work_array(nNodL).
        if (is_multirank(partit)) call exchange_nod(arr, partit)
        glob_max = maxval(arr(1,:))
        if (is_multirank(partit)) call allreduce_max(glob_max, partit)
        do while (glob_max > 0.99_WP*dummy)
            ! horizontal extrapolation
            do nz = 1, mesh%nl-1
                work_array = arr(nz,:)
                success = .true.
                do while (success)               ! runs as long as success==.true.
                    success = .false.
                    do n = 1, nNodO            ! OWNED nodes (FESOM2: do n=1, myDim_nod2D)
                        if ((work_array(n) > 0.99_WP*dummy) .and. (mesh%nlevels_nod2D(n) > nz)) then
                            cnt = 0
                            val = 0._WP
                            do k = 1, mesh%nod_in_elem2D_num(n)
                                el = mesh%nod_in_elem2D(k, n)
                                if (nz > mesh%nlevels(el)) cycle
                                enodes = mesh%elem2D_nodes(:, el)
                                do j = 1, mesh%elem2D_nnodes(el)   ! UNIQUE vertices (a triangle's slot 4 repeats slot 1)
                                    if (enodes(j) == 0) cycle
                                    if ((work_array(enodes(j)) < 0.99_WP*dummy) .and. &
                                        (mesh%nlevels_nod2D(enodes(j)) > nz)) then
                                        val = val + work_array(enodes(j))
                                        cnt = cnt + 1
                                    end if
                                end do
                            end do
                            if (cnt > 0) then
                                work_array(n) = val/real(cnt, WP)
                                success = .true.
                            end if
                        end if
                    end do
                end do
                arr(nz,:) = work_array
            end do
            ! propagate the filled owned values into neighbours' halos, then re-check the
            ! global surface max (FESOM2 gen_support.F90:484-488).
            if (is_multirank(partit)) call exchange_nod(arr, partit)
            glob_max = maxval(arr(1,:))
            if (is_multirank(partit)) call allreduce_max(glob_max, partit)
        end do
        ! vertical extrapolation (owned nodes; FESOM2: do n=1, myDim_nod2D)
        do n = 1, nNodO
            nl1 = mesh%nlevels_nod2D(n) - 1
            do nz = 2, nl1
                if (arr(nz,n) > 0.99_WP*dummy) arr(nz,n) = arr(nz-1,n)
            end do
        end do
        if (is_multirank(partit)) call exchange_nod(arr, partit)   ! FESOM2 final exchange
        deallocate(work_array)
    end subroutine extrap_nod3D

    !===========================================================================
    ! Orchestrate: per (file,var) read+interp+extrap into the matching tracer; zero
    ! land/bottom; Kelvin guard; in-situ->potential. FESOM2 do_ic3d (493-644).
    subroutine do_ic3d(tracers, ic, mesh, partit)
        type(t_tracer),      intent(inout) :: tracers
        type(t_ic3d_config), intent(in)    :: ic
        type(t_mesh),        intent(in)    :: mesh
        type(t_partit),      intent(in), optional :: partit
        integer :: n, ct, nNodO, nNodL, nEdgeO, nElemO
        real(WP) :: dummy
        dummy = ic%dummy
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        allocate(bilin_indx_i(nNodL), bilin_indx_j(nNodL))
        do n = 1, ic%n_ic3d
            do ct = 1, tracers%num_tracers
                if (tracers%data(ct)%ID == ic%idlist(n)) then
                    call nc_ic3d_ini(trim(ic%filelist(n)), ic, mesh, partit)
                    call getcoeffld(trim(ic%filelist(n)), trim(ic%varlist(n)), ic, &
                                    tracers%data(ct)%values, mesh, partit)
                    call nc_end()
                    call extrap_nod3D(tracers%data(ct)%values, ic, mesh, partit)
                    exit
                elseif (ct == tracers%num_tracers) then
                    write(*,*) 'do_ic3d: idlist contains tracer ID not in tracer list: ', ic%idlist(n)
                    error stop 1
                end if
            end do
        end do
        deallocate(bilin_indx_i, bilin_indx_j)

        do ct = 1, tracers%num_tracers
            ! set remaining dummy values (bottom topography / unfilled) to 0
            where (tracers%data(ct)%values > 0.9_WP*dummy)
                tracers%data(ct)%values = 0.0_WP
            end where
            ! ensure bottom is zero (no cavity -> no surface zeroing); owned+halo
            ! (FESOM2 do_ic3d: do n=1, myDim_nod2D+eDim_nod2D).
            do n = 1, nNodL
                tracers%data(ct)%values(mesh%nlevels_nod2D(n):mesh%nl-1, n) = 0.0_WP
            end do
        end do
        ! convert temperature from Kelvin -> degC (vacuous for phc3.0 Celsius data)
        where (tracers%data(1)%values(:,:) > 100._WP)
            tracers%data(1)%values(:,:) = tracers%data(1)%values(:,:) - 273.15_WP
        end where

        if (ic%t_insitu) then
            call insitu2pot(tracers%data(1)%values, tracers%data(2)%values, mesh, partit)
        end if
    end subroutine do_ic3d

    !===========================================================================
    subroutine nc_end()
        deallocate(nc_lon, nc_lat, nc_depth)
    end subroutine nc_end

end module oce_initial_state
