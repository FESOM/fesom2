module mod_tracks_geometry
  ! Pure-Fortran port of tripyview's transect geometry pipeline:
  !   _do_calc_csect_vec
  !   _do_find_intersected_edges_fastnew
  !   _do_build_path (+ __add_upsection_elem2path / __add_downsection_elem2path)
  !   _do_compute_distance_from_startpoint
  !
  ! Reference: /work/ab0246/a270092/software/tripyview/tripyview/sub_transect.py
  !
  ! All inputs are GLOBAL mesh arrays — caller must gather/broadcast them.
  ! Module has no MPI, no XIOS, no FESOM dependencies; safe to unit-test
  ! in isolation against tripyview output.
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private

  integer, parameter, public :: TG_WP = selected_real_kind(13)

  real(TG_WP), parameter :: PI       = 3.14159265358979323846_TG_WP
  real(TG_WP), parameter :: DEG2RAD  = PI / 180.0_TG_WP
  real(TG_WP), parameter :: R_EARTH_KM = 6367.5_TG_WP

  type, public :: transect_t
    character(len=64)         :: name = ''
    integer                   :: M    = 0   ! # crossed edges
    integer                   :: P    = 0   ! # path triangle entries
    integer                   :: Pxy  = 0   ! # path_xy vertices (differs from P near land)

    integer,     allocatable  :: edge_cut_i   (:)        ! (M) global edge id (-1 = inserted land slot)
    integer,     allocatable  :: edge_cut_ni  (:, :)     ! (2, M) global node ids
    real(TG_WP), allocatable  :: edge_cut_lint(:)        ! (M)
    real(TG_WP), allocatable  :: edge_cut_P   (:, :)     ! (2, M) lon, lat (degrees)
    real(TG_WP), allocatable  :: edge_cut_midP(:, :)     ! (2, M) lon, lat (degrees)
    real(TG_WP), allocatable  :: edge_cut_evec(:, :)     ! (2, M) edge along-vector (degrees)
    real(TG_WP), allocatable  :: edge_cut_dist(:)        ! (M) km from start

    integer,     allocatable  :: path_ei         (:)     ! (P) global elem id, -1 = land slot
    integer,     allocatable  :: path_ni         (:, :)  ! (3, P) global node ids
    integer,     allocatable  :: path_cut_ni     (:, :)  ! (2, P) global node ids of entering edge
    real(TG_WP), allocatable  :: path_dx         (:)     ! (P) midpoint→centroid x (metres)
    real(TG_WP), allocatable  :: path_dy         (:)     ! (P) midpoint→centroid y (metres)
    real(TG_WP), allocatable  :: path_nvec_cs    (:, :)  ! (2, P) section-normal direction
    real(TG_WP), allocatable  :: path_centroid_xy(:, :)  ! (2, P) lon, lat of path centroid

    real(TG_WP), allocatable  :: path_xy   (:, :)        ! (2, Pxy) alternating midP/centroid
    real(TG_WP), allocatable  :: path_dist (:)           ! (Pxy) km from start
  end type transect_t

  public :: analyse_transect
  public :: free_transect

  ! -- private helpers --
  type :: subseg_t
    integer                   :: ncut = 0
    integer,     allocatable  :: edge_cut_i   (:)
    integer,     allocatable  :: edge_cut_ni  (:, :)
    real(TG_WP), allocatable  :: edge_cut_lint(:)
    real(TG_WP), allocatable  :: edge_cut_P   (:, :)
    real(TG_WP), allocatable  :: edge_cut_midP(:, :)
    real(TG_WP), allocatable  :: edge_cut_evec(:, :)
    real(TG_WP)               :: e_vec(2)   = 0.0_TG_WP
    real(TG_WP)               :: n_vec(2)   = 0.0_TG_WP
    real(TG_WP)               :: alpha      = 0.0_TG_WP

    integer                   :: P    = 0
    integer                   :: Pxy  = 0
    integer,     allocatable  :: path_ei         (:)
    integer,     allocatable  :: path_ni         (:, :)
    integer,     allocatable  :: path_cut_ni     (:, :)
    real(TG_WP), allocatable  :: path_dx         (:)
    real(TG_WP), allocatable  :: path_dy         (:)
    real(TG_WP), allocatable  :: path_nvec_cs    (:, :)   ! (2, P) per-sub-segment n_vec
    real(TG_WP), allocatable  :: path_centroid_xy(:, :)
    real(TG_WP), allocatable  :: path_xy         (:, :)
  end type subseg_t

contains

! -----------------------------------------------------------------------------
!> Full geometry pipeline: take a CSV polyline and a global FESOM mesh,
!> return a populated transect_t holding the list of crossed edges, the
!> alternating triangle path, the cross-vectors and the cumulative
!> great-circle distance.
!>
!> in:  n_nod, n_edge, n_elem   array sizes for the global mesh
!>      lon_nod, lat_nod        global node coordinates [deg, geographic]
!>      edge                    global edge -> 2 node ids (1-based)
!>      edge_tri                global edge -> 2 element ids (-1 = boundary)
!>      edge_cross_dxdy         midpoint-to-centroid (4, n_edge) [m]:
!>                                rows 1:2 left tri, rows 3:4 right tri
!>      elem_nodes              global element -> 3 node ids (1-based)
!>      n_csi, csi_lon/csi_lat  CSV polyline (degrees)
!>      name                    transect label
!> out: transect                fully populated transect_t
! -----------------------------------------------------------------------------
  subroutine analyse_transect(n_nod, lon_nod, lat_nod, &
                              n_edge, edge, edge_tri, edge_cross_dxdy, &
                              n_elem, elem_nodes, &
                              n_csi, csi_lon, csi_lat, name, &
                              transect)
    integer,          intent(in)    :: n_nod, n_edge, n_elem, n_csi
    real(TG_WP),      intent(in)    :: lon_nod(n_nod), lat_nod(n_nod)        ! degrees, geographic
    integer,          intent(in)    :: edge(2, n_edge)                       ! global node ids (1-based)
    integer,          intent(in)    :: edge_tri(2, n_edge)                   ! global elem ids (1-based, -1 = boundary)
    real(TG_WP),      intent(in)    :: edge_cross_dxdy(4, n_edge)            ! 1:2 LEFT, 3:4 RIGHT (metres)
    integer,          intent(in)    :: elem_nodes(3, n_elem)                 ! global node ids (1-based)
    real(TG_WP),      intent(in)    :: csi_lon(n_csi), csi_lat(n_csi)        ! degrees
    character(len=*), intent(in)    :: name
    type(transect_t), intent(out)   :: transect

    type(subseg_t), allocatable :: segs(:)
    real(TG_WP)                 :: poly_lon(n_csi), poly_lat(n_csi)
    real(TG_WP)                 :: auxx, auxy, auxn, alpha_tot
    real(TG_WP)                 :: n_vec_tot(2)
    integer                     :: ii, n_segs
    logical                     :: crosses_dl
    real(TG_WP)                 :: px0, px1

    transect%name = name

    ! ------- whole-transect orientation (Px[-1]-Px[0], Py[-1]-Py[0]) -------
    poly_lon = csi_lon
    poly_lat = csi_lat
    auxx = poly_lon(n_csi) - poly_lon(1)
    auxy = poly_lat(n_csi) - poly_lat(1)
    if (auxx >  180.0_TG_WP) auxx = auxx - 360.0_TG_WP
    if (auxx < -180.0_TG_WP) auxx = auxx + 360.0_TG_WP
    auxn = sqrt(auxx*auxx + auxy*auxy)
    if (auxn > 0.0_TG_WP) then
      auxx = auxx / auxn
      auxy = auxy / auxn
    end if
    alpha_tot = atan2(auxy, auxx) * 180.0_TG_WP / PI

    ! tripyview flips the polyline if total bearing falls in [-180, -90]
    if (alpha_tot >= -180.0_TG_WP .and. alpha_tot <= -90.0_TG_WP) then
      poly_lon = poly_lon(n_csi:1:-1)
      poly_lat = poly_lat(n_csi:1:-1)
      auxx = -auxx
      auxy = -auxy
    end if
    n_vec_tot = (/ -auxy, auxx /)

    ! ------- per-subsegment work -------
    allocate(segs(n_csi - 1))
    n_segs = 0
    do ii = 1, n_csi - 1
      px0 = poly_lon(ii)
      px1 = poly_lon(ii + 1)
      crosses_dl = abs(px1 - px0) > 180.0_TG_WP
      if (crosses_dl) then
        if (px0 < 0.0_TG_WP) px0 = px0 + 360.0_TG_WP
        if (px1 < 0.0_TG_WP) px1 = px1 + 360.0_TG_WP
      end if

      call calc_csect_vec(px0, px1, poly_lat(ii), poly_lat(ii + 1), segs(ii))

      call find_intersected_edges(n_nod, lon_nod, lat_nod, n_edge, edge, &
                                  px0, px1, poly_lat(ii), poly_lat(ii + 1), &
                                  crosses_dl, segs(ii))

      ! shift edge_cut_P back to [-180,180] for downstream concat/plotting
      if (crosses_dl .and. segs(ii)%ncut > 0) then
        where (segs(ii)%edge_cut_P   (1, :) > 180.0_TG_WP) &
          segs(ii)%edge_cut_P   (1, :) = segs(ii)%edge_cut_P   (1, :) - 360.0_TG_WP
        where (segs(ii)%edge_cut_midP(1, :) > 180.0_TG_WP) &
          segs(ii)%edge_cut_midP(1, :) = segs(ii)%edge_cut_midP(1, :) - 360.0_TG_WP
      end if

      if (segs(ii)%ncut == 0) cycle    ! no crossings — skip this sub-segment

      call build_path(n_nod, lon_nod, lat_nod, n_elem, elem_nodes, &
                      edge_tri, edge_cross_dxdy, &
                      ii, n_csi - 1, segs(ii))
      n_segs = n_segs + 1
    end do

    ! ------- concat sub-segments into the public transect_t -------
    call concat_subtransects(segs, n_csi - 1, transect)
    deallocate(segs)

    ! ------- flip sign if total normal points SW/S (positive transport S->N, W->E) -------
    ! Matches tripyview's _do_build_path:875 which gates on n_vec_tot
    ! and flips path_dx/path_dy uniformly across all sub-segments.
    ! path_nvec_cs is left as-is from concat_subtransects: each path step
    ! carries its sub-segment's own n_vec (tripyview's _do_build_path:855).
    if (transect%P > 0) then
      if (n_vec_tot(1) < 0.0_TG_WP .or. n_vec_tot(2) < 0.0_TG_WP) then
        transect%path_dx = -transect%path_dx
        transect%path_dy = -transect%path_dy
      end if
    end if

    ! Note: tripyview's _do_insert_landpts adds synthetic NaN slots
    ! between consecutive boundary edges (polyline crossing a continent).
    ! Not ported -- the within-segment boundary case is handled by
    ! push_ghost at build_path time. Add later if a use case appears.

    call compute_distance_from_startpoint(transect)
  end subroutine analyse_transect

  !> Deallocate every allocatable array in t and reset the counters.
  !> Safe to call on a partially-built or already-freed transect.
  subroutine free_transect(t)
    type(transect_t), intent(inout) :: t
    if (allocated(t%edge_cut_i))     deallocate(t%edge_cut_i)
    if (allocated(t%edge_cut_ni))    deallocate(t%edge_cut_ni)
    if (allocated(t%edge_cut_lint))  deallocate(t%edge_cut_lint)
    if (allocated(t%edge_cut_P))     deallocate(t%edge_cut_P)
    if (allocated(t%edge_cut_midP))  deallocate(t%edge_cut_midP)
    if (allocated(t%edge_cut_evec))  deallocate(t%edge_cut_evec)
    if (allocated(t%edge_cut_dist))  deallocate(t%edge_cut_dist)
    if (allocated(t%path_ei))        deallocate(t%path_ei)
    if (allocated(t%path_ni))        deallocate(t%path_ni)
    if (allocated(t%path_cut_ni))    deallocate(t%path_cut_ni)
    if (allocated(t%path_dx))        deallocate(t%path_dx)
    if (allocated(t%path_dy))        deallocate(t%path_dy)
    if (allocated(t%path_nvec_cs))   deallocate(t%path_nvec_cs)
    if (allocated(t%path_centroid_xy)) deallocate(t%path_centroid_xy)
    if (allocated(t%path_xy))        deallocate(t%path_xy)
    if (allocated(t%path_dist))      deallocate(t%path_dist)
    t%M   = 0
    t%P   = 0
    t%Pxy = 0
  end subroutine free_transect

! -----------------------------------------------------------------------------
!> Fill seg%e_vec, seg%n_vec and seg%alpha (radians) for the lat-lon
!> straight sub-segment from (x0, y0) to (x1, y1) in degrees. Port of
!> tripyview's _do_calc_csect_vec.
! -----------------------------------------------------------------------------
  subroutine calc_csect_vec(x0, x1, y0, y1, seg)
    real(TG_WP),     intent(in)    :: x0, x1, y0, y1
    type(subseg_t),  intent(inout) :: seg
    real(TG_WP) :: dx, dy, dn

    dx = x1 - x0
    dy = y1 - y0
    dn = sqrt(dx*dx + dy*dy)
    if (dn > 0.0_TG_WP) then
      seg%e_vec = (/  dx/dn,  dy/dn /)
      ! Per-sub-segment normal with tripyview's conditional flip:
      ! going N (dy > 0) or W (dx < 0) -> clockwise rotation;
      ! otherwise counter-clockwise. Matches sub_transect.py:437-440.
      if (dy > 0.0_TG_WP .or. dx < 0.0_TG_WP) then
        seg%n_vec = (/  dy/dn, -dx/dn /)
      else
        seg%n_vec = (/ -dy/dn,  dx/dn /)
      end if
    else
      seg%e_vec = 0.0_TG_WP
      seg%n_vec = 0.0_TG_WP
    end if
    seg%alpha = atan2(seg%e_vec(2), seg%e_vec(1))   ! radians
  end subroutine calc_csect_vec

! -----------------------------------------------------------------------------
!> For one sub-segment, find every mesh edge it crosses and record the
!> per-edge geometry needed downstream. Bounding-box pre-filter widens
!> to 180 deg in longitude at high latitudes, then a signed-distance
!> test on the line equation a*x + b*y + c = 0 picks crossings and
!> sorts them by progression along the sub-segment. Port of tripyview's
!> _do_find_intersected_edges_fastnew.
!>
!> When crosses_dl is true the sub-segment wraps the dateline; both the
!> input endpoints (x0, x1) and the candidate edge endpoints are
!> shifted to a common [0, 360] frame for the intersection test.
! -----------------------------------------------------------------------------
  subroutine find_intersected_edges(n_nod, lon_nod, lat_nod, n_edge, edge, &
                                    x0, x1, y0, y1, crosses_dl, seg)
    integer,         intent(in)    :: n_nod, n_edge
    real(TG_WP),     intent(in)    :: lon_nod(n_nod), lat_nod(n_nod)
    integer,         intent(in)    :: edge(2, n_edge)
    real(TG_WP),     intent(in)    :: x0, x1, y0, y1
    logical,         intent(in)    :: crosses_dl
    type(subseg_t),  intent(inout) :: seg

    ! Bounding box (with polar widening)
    real(TG_WP) :: Pxmin, Pxmax, Pymin, Pymax, Pdx, Pdy
    real(TG_WP) :: a, b, c, dx10, dy10, dd10
    integer     :: n_cand, n_hit
    integer,     allocatable :: cand_idx(:), hit_idx(:), srt(:)
    real(TG_WP), allocatable :: x0e(:), y0e(:), x1e(:), y1e(:), dxe(:), dye(:)
    real(TG_WP), allocatable :: s0(:), s1(:), tparam(:), xi(:), yi(:), fac(:)
    real(TG_WP), allocatable :: emin_lon(:), emax_lon(:), emin_lat(:), emax_lat(:)
    real(TG_WP), allocatable :: fac_hit(:)
    integer     :: e, j, k
    real(TG_WP) :: lon_a, lon_b

    Pdx = 10.0_TG_WP   ! tripyview default
    Pdy = 10.0_TG_WP
    Pxmin = min(x0, x1)
    Pxmax = max(x0, x1)
    Pymin = min(y0, y1)
    Pymax = max(y0, y1)
    ! polar latitude widening (tripyview: Pymin<-70 or Pymax>80 -> Pdx=180)
    if (Pymin < -70.0_TG_WP .or. Pymax > 80.0_TG_WP) Pdx = 180.0_TG_WP
    Pxmin = Pxmin - Pdx;  Pxmax = Pxmax + Pdx
    Pymin = Pymin - Pdy;  Pymax = Pymax + Pdy

    ! Build per-edge bbox in the right longitude frame
    allocate(emin_lon(n_edge), emax_lon(n_edge), &
             emin_lat(n_edge), emax_lat(n_edge))
    do e = 1, n_edge
      lon_a = lon_nod(edge(1, e))
      lon_b = lon_nod(edge(2, e))
      if (crosses_dl) then
        if (lon_a < 0.0_TG_WP) lon_a = lon_a + 360.0_TG_WP
        if (lon_b < 0.0_TG_WP) lon_b = lon_b + 360.0_TG_WP
      end if
      emin_lon(e) = min(lon_a, lon_b)
      emax_lon(e) = max(lon_a, lon_b)
      emin_lat(e) = min(lat_nod(edge(1, e)), lat_nod(edge(2, e)))
      emax_lat(e) = max(lat_nod(edge(1, e)), lat_nod(edge(2, e)))
    end do

    ! Candidate list: edges whose bbox lies inside expanded section bbox
    ! Tripyview uses `min >= Pxmin & max <= Pxmax & ...` (entirely inside),
    ! which is a tight pre-filter that still includes any edge a finite
    ! section could touch given the +Pdx margin.
    n_cand = 0
    allocate(cand_idx(n_edge))
    do e = 1, n_edge
      if (emin_lon(e) >= Pxmin .and. emax_lon(e) <= Pxmax .and. &
          emin_lat(e) >= Pymin .and. emax_lat(e) <= Pymax) then
        n_cand = n_cand + 1
        cand_idx(n_cand) = e
      end if
    end do
    deallocate(emin_lon, emax_lon, emin_lat, emax_lat)

    if (n_cand == 0) then
      deallocate(cand_idx)
      seg%ncut = 0
      return
    end if

    ! Section line coefficients L(x,y) = a*x + b*y + c = 0
    a = y0    - y1
    b = x1    - x0
    c = x0*y1 - x1*y0
    dx10 = x1 - x0
    dy10 = y1 - y0
    dd10 = dx10*dx10 + dy10*dy10

    allocate(x0e(n_cand), y0e(n_cand), x1e(n_cand), y1e(n_cand), &
             dxe(n_cand), dye(n_cand), s0(n_cand), s1(n_cand), &
             tparam(n_cand), xi(n_cand), yi(n_cand), fac(n_cand))
    do j = 1, n_cand
      e = cand_idx(j)
      x0e(j) = lon_nod(edge(1, e))
      y0e(j) = lat_nod(edge(1, e))
      x1e(j) = lon_nod(edge(2, e))
      y1e(j) = lat_nod(edge(2, e))
      if (crosses_dl) then
        if (x0e(j) < 0.0_TG_WP) x0e(j) = x0e(j) + 360.0_TG_WP
        if (x1e(j) < 0.0_TG_WP) x1e(j) = x1e(j) + 360.0_TG_WP
      end if
      dxe(j) = x1e(j) - x0e(j)
      dye(j) = y1e(j) - y0e(j)
      s0(j)  = a*x0e(j) + b*y0e(j) + c
      s1(j)  = a*x1e(j) + b*y1e(j) + c
    end do

    n_hit = 0
    allocate(hit_idx(n_cand), fac_hit(n_cand))
    do j = 1, n_cand
      if ((s1(j) - s0(j)) == 0.0_TG_WP) cycle    ! parallel — no isolated crossing
      tparam(j) = s0(j) / (s0(j) - s1(j))
      xi(j) = x0e(j) + tparam(j) * dxe(j)
      yi(j) = y0e(j) + tparam(j) * dye(j)
      fac(j) = ( (xi(j) - x0)*dx10 + (yi(j) - y0)*dy10 ) / dd10
      if (s0(j)*s1(j) < 0.0_TG_WP .and. fac(j) >= 0.0_TG_WP .and. fac(j) <= 1.0_TG_WP) then
        n_hit = n_hit + 1
        hit_idx(n_hit) = j
        fac_hit(n_hit) = fac(j)
      end if
    end do

    if (n_hit == 0) then
      seg%ncut = 0
      deallocate(cand_idx, x0e, y0e, x1e, y1e, dxe, dye, s0, s1, &
                 tparam, xi, yi, fac, hit_idx, fac_hit)
      return
    end if

    ! Sort by progression along the segment
    allocate(srt(n_hit))
    call argsort_real(fac_hit(1:n_hit), srt)

    ! Allocate seg storage
    allocate(seg%edge_cut_i   (n_hit))
    allocate(seg%edge_cut_ni  (2, n_hit))
    allocate(seg%edge_cut_lint(n_hit))
    allocate(seg%edge_cut_P   (2, n_hit))
    allocate(seg%edge_cut_midP(2, n_hit))
    allocate(seg%edge_cut_evec(2, n_hit))

    do k = 1, n_hit
      j = hit_idx(srt(k))
      e = cand_idx(j)
      seg%edge_cut_i   (k)    = e
      seg%edge_cut_ni  (1, k) = edge(1, e)
      seg%edge_cut_ni  (2, k) = edge(2, e)
      seg%edge_cut_lint(k)    = tparam(j)
      seg%edge_cut_P   (1, k) = xi(j)
      seg%edge_cut_P   (2, k) = yi(j)
      seg%edge_cut_midP(1, k) = x0e(j) + 0.5_TG_WP * dxe(j)
      seg%edge_cut_midP(2, k) = y0e(j) + 0.5_TG_WP * dye(j)
      seg%edge_cut_evec(1, k) = dxe(j)
      seg%edge_cut_evec(2, k) = dye(j)
    end do
    seg%ncut = n_hit

    deallocate(cand_idx, x0e, y0e, x1e, y1e, dxe, dye, s0, s1, &
               tparam, xi, yi, fac, hit_idx, fac_hit, srt)
  end subroutine find_intersected_edges

! -----------------------------------------------------------------------------
!> Build the alternating midpoint -> centroid -> midpoint path that the
!> sub-segment makes through the mesh triangles, with per-step cross-
!> vectors (path_dx, path_dy). Each crossed edge contributes an
!> up-section step (incoming triangle, signed -edge_cross_dxdy) and a
!> down-section step (outgoing triangle, +edge_cross_dxdy). Boundary
!> edges inject ghost slots with NaN dx/dy and path_ei = -1. Port of
!> tripyview's _do_build_path + __add_*section_elem2path.
!>
!> ncsi is this sub-segment's 1-based index, ncs the total sub-segment
!> count; both influence boundary-edge ghost handling at sub-segment
!> ends. edge_cross_dxdy is (4, n_edge) with rows 1:2 = left tri,
!> rows 3:4 = right tri.
! -----------------------------------------------------------------------------
  subroutine build_path(n_nod, lon_nod, lat_nod, n_elem, elem_nodes, &
                        edge_tri, edge_cross_dxdy, ncsi, ncs, seg)
    integer,        intent(in)    :: n_nod, n_elem
    real(TG_WP),    intent(in)    :: lon_nod(n_nod), lat_nod(n_nod)
    integer,        intent(in)    :: elem_nodes(3, n_elem)
    integer,        intent(in)    :: edge_tri(2, *)
    real(TG_WP),    intent(in)    :: edge_cross_dxdy(4, *)
    integer,        intent(in)    :: ncsi, ncs              ! 1-based sub-seg id, # sub-segs
    type(subseg_t), intent(inout) :: seg

    integer  :: edi, nced, eg, el(2)
    real(TG_WP) :: alpha, ca, sa, theta, auxx, auxy

    integer,     allocatable :: bei(:), bni(:, :), bcni(:, :)
    real(TG_WP), allocatable :: bdx(:), bdy(:), bcen(:, :), bxy(:, :)
    integer :: pcap, pxycap, pp, pxy_p

    nced = seg%ncut
    alpha = seg%alpha
    ca = cos(-alpha); sa = sin(-alpha)

    ! Worst-case sizing: each edi adds <= 3 path entries and <= 2 path_xy entries,
    ! plus the start-of-section centroid. Grow geometric if we ever exceed.
    pcap   = max(8, 3*nced + 4)
    pxycap = max(8, 2*nced + 4)
    allocate(bei (pcap))
    allocate(bni (3, pcap))
    allocate(bcni(2, pcap))
    allocate(bdx (pcap))
    allocate(bdy (pcap))
    allocate(bcen(2, pcap))
    allocate(bxy (2, pxycap))
    pp = 0; pxy_p = 0

    do edi = 1, nced
      eg = seg%edge_cut_i(edi)
      el(1) = edge_tri(1, eg)
      el(2) = edge_tri(2, eg)

      ! Rotate the edge unit-vector by -alpha; theta tells us which side
      ! of the section the edge points to (left if theta > 0, right if < 0).
      auxx = seg%edge_cut_evec(1, edi)*ca - seg%edge_cut_evec(2, edi)*sa
      auxy = seg%edge_cut_evec(1, edi)*sa + seg%edge_cut_evec(2, edi)*ca
      theta = atan2(auxy, auxx)

      ! ------------------ upsection element ------------------
      call add_upsection(edi, nced, ncsi, ncs, theta, el, eg, &
                         n_nod, lon_nod, lat_nod, n_elem, elem_nodes, &
                         edge_cross_dxdy, seg, &
                         pp, pxy_p, pcap, pxycap, &
                         bei, bni, bcni, bdx, bdy, bcen, bxy)

      ! midpoint of crossed edge (always)
      call grow_xy(bxy, pxycap, pxy_p + 1)
      pxy_p = pxy_p + 1
      bxy(1, pxy_p) = seg%edge_cut_midP(1, edi)
      bxy(2, pxy_p) = seg%edge_cut_midP(2, edi)

      ! ------------------ downsection element ------------------
      call add_downsection(edi, nced, theta, el, eg, &
                           n_nod, lon_nod, lat_nod, n_elem, elem_nodes, &
                           edge_cross_dxdy, seg, &
                           pp, pxy_p, pcap, pxycap, &
                           bei, bni, bcni, bdx, bdy, bcen, bxy)
    end do

    ! Trim and move to seg
    seg%P   = pp
    seg%Pxy = pxy_p
    allocate(seg%path_ei         (pp))
    allocate(seg%path_ni         (3, pp))
    allocate(seg%path_cut_ni     (2, pp))
    allocate(seg%path_dx         (pp))
    allocate(seg%path_dy         (pp))
    allocate(seg%path_nvec_cs    (2, pp))
    allocate(seg%path_centroid_xy(2, pp))
    allocate(seg%path_xy         (2, pxy_p))

    seg%path_ei         = bei(1:pp)
    seg%path_ni         = bni(:, 1:pp)
    seg%path_cut_ni     = bcni(:, 1:pp)
    seg%path_dx         = bdx(1:pp)
    seg%path_dy         = bdy(1:pp)
    seg%path_centroid_xy = bcen(:, 1:pp)
    seg%path_xy         = bxy(:, 1:pxy_p)

    ! Tag every path step inside this sub-segment with its own n_vec,
    ! matching tripyview's _do_build_path:855-857.
    seg%path_nvec_cs(1, :) = seg%n_vec(1)
    seg%path_nvec_cs(2, :) = seg%n_vec(2)

    deallocate(bei, bni, bcni, bdx, bdy, bcen, bxy)
  end subroutine build_path

  !> Append the up-section triangle for the current crossed edge to the
  !> growing path. theta's sign decides whether the up-section side is
  !> the LEFT (>=0) or the RIGHT (<0) mesh triangle. el(1:2) are the two
  !> elem ids around the crossed edge; eg is the global edge id. Boundary
  !> edges (el(.)==-1) push a NaN ghost slot; the first edge of a sub-
  !> segment also adds a centroid path_xy point.
  subroutine add_upsection(edi, nced, ncsi, ncs, theta, el, eg, &
                           n_nod, lon_nod, lat_nod, n_elem, elem_nodes, &
                           edge_cross_dxdy, seg, &
                           pp, pxy_p, pcap, pxycap, &
                           bei, bni, bcni, bdx, bdy, bcen, bxy)
    integer,        intent(in)    :: edi, nced, ncsi, ncs
    real(TG_WP),    intent(in)    :: theta
    integer,        intent(in)    :: el(2), eg
    integer,        intent(in)    :: n_nod, n_elem
    real(TG_WP),    intent(in)    :: lon_nod(n_nod), lat_nod(n_nod)
    integer,        intent(in)    :: elem_nodes(3, n_elem)
    real(TG_WP),    intent(in)    :: edge_cross_dxdy(4, *)
    type(subseg_t), intent(in)    :: seg
    integer,        intent(inout) :: pp, pxy_p, pcap, pxycap
    integer,        allocatable, intent(inout) :: bei(:), bni(:, :), bcni(:, :)
    real(TG_WP),    allocatable, intent(inout) :: bdx(:), bdy(:), bcen(:, :), bxy(:, :)
    real(TG_WP) :: cx, cy

    if (theta >= 0.0_TG_WP) then
      ! upsection = LEFT triangle (el(1)), always interior because el(1) > 0
      if (edi == 1) then
        call elem_centroid(el(1), n_nod, lon_nod, lat_nod, n_elem, elem_nodes, cx, cy)
        call grow_xy(bxy, pxycap, pxy_p + 1)
        pxy_p = pxy_p + 1
        bxy(1, pxy_p) = cx; bxy(2, pxy_p) = cy
      end if
      call elem_centroid(el(1), n_nod, lon_nod, lat_nod, n_elem, elem_nodes, cx, cy)
      call grow_p(bei, bni, bcni, bdx, bdy, bcen, pcap, pp + 1)
      pp = pp + 1
      bei(pp)        = el(1)
      bni(:, pp)     = elem_nodes(:, el(1))
      bcni(:, pp)    = seg%edge_cut_ni(:, edi)
      bdx(pp)        = -edge_cross_dxdy(1, eg)
      bdy(pp)        = -edge_cross_dxdy(2, eg)
      bcen(1, pp)    = cx; bcen(2, pp) = cy
    else
      ! upsection = RIGHT triangle (el(2)); may be -1 for boundary
      if (el(2) > 0) then
        if (edi == 1) then
          call elem_centroid(el(2), n_nod, lon_nod, lat_nod, n_elem, elem_nodes, cx, cy)
          call grow_xy(bxy, pxycap, pxy_p + 1)
          pxy_p = pxy_p + 1
          bxy(1, pxy_p) = cx; bxy(2, pxy_p) = cy
        end if
        call elem_centroid(el(2), n_nod, lon_nod, lat_nod, n_elem, elem_nodes, cx, cy)
        call grow_p(bei, bni, bcni, bdx, bdy, bcen, pcap, pp + 1)
        pp = pp + 1
        bei(pp)     = el(2)
        bni(:, pp)  = elem_nodes(:, el(2))
        bcni(:, pp) = seg%edge_cut_ni(:, edi)
        bdx(pp)     = -edge_cross_dxdy(3, eg)
        bdy(pp)     = -edge_cross_dxdy(4, eg)
        bcen(1, pp) = cx; bcen(2, pp) = cy
      else
        ! Boundary right triangle missing → ghost entry, with the bookkeeping
        ! tripyview applies to keep path_xy/path_ei coherent across sub-segments.
        call grow_xy(bxy, pxycap, pxy_p + 1)
        pxy_p = pxy_p + 1
        bxy(1, pxy_p) = seg%edge_cut_midP(1, edi)
        bxy(2, pxy_p) = seg%edge_cut_midP(2, edi)
        call push_ghost(bei, bni, bcni, bdx, bdy, bcen, pcap, pp)

        if (edi /= 1 .and. edi /= nced) then
          call push_ghost(bei, bni, bcni, bdx, bdy, bcen, pcap, pp)
        else if (ncsi /= 1 .and. edi == 1) then
          call grow_xy(bxy, pxycap, pxy_p + 1)
          pxy_p = pxy_p + 1
          bxy(1, pxy_p) = seg%edge_cut_midP(1, edi)
          bxy(2, pxy_p) = seg%edge_cut_midP(2, edi)
          call push_ghost(bei, bni, bcni, bdx, bdy, bcen, pcap, pp)
        else if (ncsi /= ncs .and. edi == nced) then
          call grow_xy(bxy, pxycap, pxy_p + 1)
          pxy_p = pxy_p + 1
          bxy(1, pxy_p) = seg%edge_cut_midP(1, edi)
          bxy(2, pxy_p) = seg%edge_cut_midP(2, edi)
          call push_ghost(bei, bni, bcni, bdx, bdy, bcen, pcap, pp)
        end if
      end if
    end if
  end subroutine add_upsection

  !> Symmetric counterpart to add_upsection: same arguments, opposite
  !> triangle choice (theta>=0 -> RIGHT, theta<0 -> LEFT). Always inserts
  !> a centroid path_xy point and an elem/node entry (or a NaN ghost slot
  !> for boundary edges).
  subroutine add_downsection(edi, nced, theta, el, eg, &
                             n_nod, lon_nod, lat_nod, n_elem, elem_nodes, &
                             edge_cross_dxdy, seg, &
                             pp, pxy_p, pcap, pxycap, &
                             bei, bni, bcni, bdx, bdy, bcen, bxy)
    integer,        intent(in)    :: edi, nced
    real(TG_WP),    intent(in)    :: theta
    integer,        intent(in)    :: el(2), eg
    integer,        intent(in)    :: n_nod, n_elem
    real(TG_WP),    intent(in)    :: lon_nod(n_nod), lat_nod(n_nod)
    integer,        intent(in)    :: elem_nodes(3, n_elem)
    real(TG_WP),    intent(in)    :: edge_cross_dxdy(4, *)
    type(subseg_t), intent(in)    :: seg
    integer,        intent(inout) :: pp, pxy_p, pcap, pxycap
    integer,        allocatable, intent(inout) :: bei(:), bni(:, :), bcni(:, :)
    real(TG_WP),    allocatable, intent(inout) :: bdx(:), bdy(:), bcen(:, :), bxy(:, :)
    real(TG_WP) :: cx, cy

    if (theta >= 0.0_TG_WP) then
      ! downsection = RIGHT triangle (el(2)); may be -1 for boundary
      if (el(2) > 0) then
        call elem_centroid(el(2), n_nod, lon_nod, lat_nod, n_elem, elem_nodes, cx, cy)
        call grow_xy(bxy, pxycap, pxy_p + 1)
        pxy_p = pxy_p + 1
        bxy(1, pxy_p) = cx; bxy(2, pxy_p) = cy

        call grow_p(bei, bni, bcni, bdx, bdy, bcen, pcap, pp + 1)
        pp = pp + 1
        bei(pp)     = el(2)
        bni(:, pp)  = elem_nodes(:, el(2))
        bcni(:, pp) = seg%edge_cut_ni(:, edi)
        bdx(pp)     = edge_cross_dxdy(3, eg)
        bdy(pp)     = edge_cross_dxdy(4, eg)
        bcen(1, pp) = cx; bcen(2, pp) = cy
      else
        call grow_xy(bxy, pxycap, pxy_p + 1)
        pxy_p = pxy_p + 1
        bxy(1, pxy_p) = seg%edge_cut_midP(1, edi)
        bxy(2, pxy_p) = seg%edge_cut_midP(2, edi)
        call push_ghost(bei, bni, bcni, bdx, bdy, bcen, pcap, pp)
      end if
    else
      ! downsection = LEFT triangle (el(1)), always interior
      call elem_centroid(el(1), n_nod, lon_nod, lat_nod, n_elem, elem_nodes, cx, cy)
      call grow_xy(bxy, pxycap, pxy_p + 1)
      pxy_p = pxy_p + 1
      bxy(1, pxy_p) = cx; bxy(2, pxy_p) = cy

      call grow_p(bei, bni, bcni, bdx, bdy, bcen, pcap, pp + 1)
      pp = pp + 1
      bei(pp)     = el(1)
      bni(:, pp)  = elem_nodes(:, el(1))
      bcni(:, pp) = seg%edge_cut_ni(:, edi)
      bdx(pp)     = edge_cross_dxdy(1, eg)
      bdy(pp)     = edge_cross_dxdy(2, eg)
      bcen(1, pp) = cx; bcen(2, pp) = cy
    end if
  end subroutine add_downsection

  !> Append one synthetic NaN path step (-1 element id, NaN dx/dy, NaN
  !> centroid). Keeps the path index in sync with the edge index when a
  !> boundary edge has no far-side triangle.
  subroutine push_ghost(bei, bni, bcni, bdx, bdy, bcen, pcap, pp)
    integer,     allocatable, intent(inout) :: bei(:), bni(:, :), bcni(:, :)
    real(TG_WP), allocatable, intent(inout) :: bdx(:), bdy(:), bcen(:, :)
    integer,                  intent(inout) :: pcap, pp
    real(TG_WP) :: nan
    nan = ieee_value(0.0_TG_WP, ieee_quiet_nan)
    call grow_p(bei, bni, bcni, bdx, bdy, bcen, pcap, pp + 1)
    pp = pp + 1
    bei(pp)     = -1
    bni(:, pp)  = -1
    bcni(:, pp) = -1
    bdx(pp)     = nan
    bdy(pp)     = nan
    bcen(:, pp) = nan
  end subroutine push_ghost

! -----------------------------------------------------------------------------
!> Concatenate the per-sub-segment scratch into a single transect_t,
!> allocated at the total size (sum of ncut and P over all sub-segments).
!> Copies edge_cut_* and path_* blocks in CSV-vertex order; zero-
!> initialises path_nvec_cs; skips sub-segments with ncut == 0.
! -----------------------------------------------------------------------------
  subroutine concat_subtransects(segs, n_subs, t)
    type(subseg_t),   intent(in)    :: segs(:)
    integer,          intent(in)    :: n_subs
    type(transect_t), intent(inout) :: t
    integer :: ii, tot_M, tot_P, tot_Pxy
    integer :: off_M, off_P, off_Pxy

    tot_M = 0; tot_P = 0; tot_Pxy = 0
    do ii = 1, n_subs
      tot_M   = tot_M   + segs(ii)%ncut
      tot_P   = tot_P   + segs(ii)%P
      tot_Pxy = tot_Pxy + segs(ii)%Pxy
    end do

    t%M   = tot_M
    t%P   = tot_P
    t%Pxy = tot_Pxy
    if (tot_M == 0) return

    allocate(t%edge_cut_i      (tot_M))
    allocate(t%edge_cut_ni     (2, tot_M))
    allocate(t%edge_cut_lint   (tot_M))
    allocate(t%edge_cut_P      (2, tot_M))
    allocate(t%edge_cut_midP   (2, tot_M))
    allocate(t%edge_cut_evec   (2, tot_M))

    allocate(t%path_ei         (tot_P))
    allocate(t%path_ni         (3, tot_P))
    allocate(t%path_cut_ni     (2, tot_P))
    allocate(t%path_dx         (tot_P))
    allocate(t%path_dy         (tot_P))
    allocate(t%path_nvec_cs    (2, tot_P))
    allocate(t%path_centroid_xy(2, tot_P))

    allocate(t%path_xy         (2, tot_Pxy))

    off_M = 0; off_P = 0; off_Pxy = 0
    do ii = 1, n_subs
      if (segs(ii)%ncut == 0) cycle
      t%edge_cut_i   (off_M+1 : off_M+segs(ii)%ncut)        = segs(ii)%edge_cut_i
      t%edge_cut_ni  (:, off_M+1 : off_M+segs(ii)%ncut)     = segs(ii)%edge_cut_ni
      t%edge_cut_lint(off_M+1 : off_M+segs(ii)%ncut)        = segs(ii)%edge_cut_lint
      t%edge_cut_P   (:, off_M+1 : off_M+segs(ii)%ncut)     = segs(ii)%edge_cut_P
      t%edge_cut_midP(:, off_M+1 : off_M+segs(ii)%ncut)     = segs(ii)%edge_cut_midP
      t%edge_cut_evec(:, off_M+1 : off_M+segs(ii)%ncut)     = segs(ii)%edge_cut_evec
      off_M = off_M + segs(ii)%ncut

      t%path_ei         (off_P+1 : off_P+segs(ii)%P)        = segs(ii)%path_ei
      t%path_ni         (:, off_P+1 : off_P+segs(ii)%P)     = segs(ii)%path_ni
      t%path_cut_ni     (:, off_P+1 : off_P+segs(ii)%P)     = segs(ii)%path_cut_ni
      t%path_dx         (off_P+1 : off_P+segs(ii)%P)        = segs(ii)%path_dx
      t%path_dy         (off_P+1 : off_P+segs(ii)%P)        = segs(ii)%path_dy
      t%path_nvec_cs    (:, off_P+1 : off_P+segs(ii)%P)     = segs(ii)%path_nvec_cs
      t%path_centroid_xy(:, off_P+1 : off_P+segs(ii)%P)     = segs(ii)%path_centroid_xy
      off_P = off_P + segs(ii)%P

      t%path_xy(:, off_Pxy+1 : off_Pxy+segs(ii)%Pxy)        = segs(ii)%path_xy
      off_Pxy = off_Pxy + segs(ii)%Pxy
    end do
  end subroutine concat_subtransects

! -----------------------------------------------------------------------------
!> Fill t%path_dist (length Pxy) and t%edge_cut_dist (length M) with
!> cumulative great-circle distance from the section start, in km.
!> Uses cartesian unit-sphere coordinates so the arc length reduces to
!> arccos(p_i . p_{i+1}) * R_earth.
! -----------------------------------------------------------------------------
  subroutine compute_distance_from_startpoint(t)
    type(transect_t), intent(inout) :: t
    real(TG_WP), allocatable :: xx(:), yy(:), zz(:), seg(:)
    real(TG_WP) :: dot, cum
    integer :: i

    if (t%Pxy > 0) then
      if (.not. allocated(t%path_dist)) allocate(t%path_dist(t%Pxy))
      allocate(xx(t%Pxy), yy(t%Pxy), zz(t%Pxy))
      do i = 1, t%Pxy
        call lonlat_to_cart3d(t%path_xy(1, i), t%path_xy(2, i), xx(i), yy(i), zz(i))
      end do
      allocate(seg(t%Pxy - 1))
      do i = 1, t%Pxy - 1
        dot = xx(i)*xx(i+1) + yy(i)*yy(i+1) + zz(i)*zz(i+1)
        if (dot > 1.0_TG_WP) dot = 1.0_TG_WP
        if (dot < -1.0_TG_WP) dot = -1.0_TG_WP
        seg(i) = acos(dot) * R_EARTH_KM
      end do
      cum = 0.0_TG_WP
      t%path_dist(1) = 0.0_TG_WP
      do i = 1, t%Pxy - 1
        cum = cum + seg(i)
        t%path_dist(i + 1) = cum
      end do
      deallocate(xx, yy, zz, seg)
    end if

    if (t%M > 0) then
      if (.not. allocated(t%edge_cut_dist)) allocate(t%edge_cut_dist(t%M))
      allocate(xx(t%M), yy(t%M), zz(t%M))
      do i = 1, t%M
        call lonlat_to_cart3d(t%edge_cut_midP(1, i), t%edge_cut_midP(2, i), &
                              xx(i), yy(i), zz(i))
      end do
      t%edge_cut_dist(1) = 0.0_TG_WP
      cum = 0.0_TG_WP
      do i = 1, t%M - 1
        dot = xx(i)*xx(i+1) + yy(i)*yy(i+1) + zz(i)*zz(i+1)
        if (dot > 1.0_TG_WP) dot = 1.0_TG_WP
        if (dot < -1.0_TG_WP) dot = -1.0_TG_WP
        cum = cum + acos(dot) * R_EARTH_KM
        t%edge_cut_dist(i + 1) = cum
      end do
      deallocate(xx, yy, zz)
    end if
  end subroutine compute_distance_from_startpoint

! -----------------------------------------------------------------------------
! small helpers
! -----------------------------------------------------------------------------

  !> Project (lon, lat) in degrees onto the unit sphere.
  subroutine lonlat_to_cart3d(lon_deg, lat_deg, x, y, z)
    real(TG_WP), intent(in)  :: lon_deg, lat_deg
    real(TG_WP), intent(out) :: x, y, z
    real(TG_WP) :: lon, lat, clat
    lon = lon_deg * DEG2RAD
    lat = lat_deg * DEG2RAD
    clat = cos(lat)
    x = clat * cos(lon)
    y = clat * sin(lon)
    z = sin(lat)
  end subroutine lonlat_to_cart3d

  !> Centroid (cx, cy) of element eid's three vertex nodes in degrees,
  !> with a [-180, 180] wrap when the element straddles the dateline.
  subroutine elem_centroid(eid, n_nod, lon_nod, lat_nod, n_elem, elem_nodes, cx, cy)
    integer,     intent(in)  :: eid, n_nod, n_elem
    real(TG_WP), intent(in)  :: lon_nod(n_nod), lat_nod(n_nod)
    integer,     intent(in)  :: elem_nodes(3, n_elem)
    real(TG_WP), intent(out) :: cx, cy
    real(TG_WP) :: xv(3), yv(3)
    integer :: i
    do i = 1, 3
      xv(i) = lon_nod(elem_nodes(i, eid))
      yv(i) = lat_nod(elem_nodes(i, eid))
    end do
    ! periodic boundary shift if the triangle straddles the dateline
    if (maxval(xv) - minval(xv) > 180.0_TG_WP) then
      if (count(xv > 0.0_TG_WP) > count(xv < 0.0_TG_WP)) then
        where (xv < 0.0_TG_WP) xv = xv + 360.0_TG_WP
      else
        where (xv > 0.0_TG_WP) xv = xv - 360.0_TG_WP
      end if
    end if
    cx = sum(xv) / 3.0_TG_WP
    cy = sum(yv) / 3.0_TG_WP
  end subroutine elem_centroid

  !> Insertion-sort permutation: idx is filled such that arr(idx(i))
  !> is non-decreasing. O(n^2); fine for the few-hundred-crossings size.
  subroutine argsort_real(arr, idx)
    real(TG_WP), intent(in)  :: arr(:)
    integer,     intent(out) :: idx(size(arr))
    integer :: i, j, tmp
    do i = 1, size(arr)
      idx(i) = i
    end do
    do i = 2, size(arr)
      tmp = idx(i)
      j = i - 1
      do while (j >= 1)
        if (arr(idx(j)) <= arr(tmp)) exit
        idx(j + 1) = idx(j)
        j = j - 1
      end do
      idx(j + 1) = tmp
    end do
  end subroutine argsort_real

  !> Geometric-growth resize of a 2 x cap buffer so it can hold at least
  !> `want` entries along its second axis. No-op when cap is large enough.
  subroutine grow_xy(buf, cap, want)
    real(TG_WP), allocatable, intent(inout) :: buf(:, :)
    integer,                  intent(inout) :: cap
    integer,                  intent(in)    :: want
    real(TG_WP), allocatable :: tmp(:, :)
    integer :: ncap
    if (want <= cap) return
    ncap = max(want, 2*cap)
    allocate(tmp(2, ncap))
    tmp(:, 1:cap) = buf
    call move_alloc(tmp, buf)
    cap = ncap
  end subroutine grow_xy

  !> Geometric-growth resize of the bundle of path-step scratch buffers
  !> so they can all hold at least `want` entries. No-op when cap is
  !> large enough.
  subroutine grow_p(bei, bni, bcni, bdx, bdy, bcen, cap, want)
    integer,     allocatable, intent(inout) :: bei(:), bni(:, :), bcni(:, :)
    real(TG_WP), allocatable, intent(inout) :: bdx(:), bdy(:), bcen(:, :)
    integer,                  intent(inout) :: cap
    integer,                  intent(in)    :: want
    integer,     allocatable :: tei(:), tni(:, :), tcni(:, :)
    real(TG_WP), allocatable :: tdx(:), tdy(:), tcen(:, :)
    integer :: ncap
    if (want <= cap) return
    ncap = max(want, 2*cap)
    allocate(tei(ncap), tni(3, ncap), tcni(2, ncap), &
             tdx(ncap), tdy(ncap), tcen(2, ncap))
    tei  (1:cap)    = bei
    tni  (:, 1:cap) = bni
    tcni (:, 1:cap) = bcni
    tdx  (1:cap)    = bdx
    tdy  (1:cap)    = bdy
    tcen (:, 1:cap) = bcen
    call move_alloc(tei,  bei)
    call move_alloc(tni,  bni)
    call move_alloc(tcni, bcni)
    call move_alloc(tdx,  bdx)
    call move_alloc(tdy,  bdy)
    call move_alloc(tcen, bcen)
    cap = ncap
  end subroutine grow_p

end module mod_tracks_geometry
