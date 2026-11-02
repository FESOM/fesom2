!===============================================================================
! real64 views of FESOM fields for REcoM
!
! REcoM keeps all its fields in real64 (its own WP = 8, src/recom). In a
! double-precision FESOM build (WP = real64) the views below are pointers to
! FESOM's own arrays: nothing is copied and nothing changes. In a single-precision
! build (WP = real32) they point to real64 copies made here, and the tracer state
! REcoM changes is rounded back into the WP tracers by recom_tracers_back, with the
! rounding residual going into the tracer carry (tracers%work%carry, see
! diff_tracers_ale). The carry keeps REcoM's sources and sinks that are smaller
! than the real32 resolution of the tracer: a deep remineralisation step of DIC is
! about 1e-3 of the real32 spacing at 2000 mmol m-3 and would otherwise be lost.
!===============================================================================
module recom_views
#if defined(__recom)
    use, intrinsic :: ISO_FORTRAN_ENV, only: real64
    use o_PARAM, only: WP, WP_full
    use MOD_MESH
    use MOD_PARTIT
    use MOD_TRACER
    use recom_glovar, only: tracers_info_type
    implicit none
    private

    public :: recom_views_mesh, recom_tracer_pointers, recom_tracer_r8, recom_tracers_back, recom_mpitypes

    ! geometry as REcoM reads it; refreshed by recom_views_mesh
    real(real64), pointer, public :: z_3d_n_r8(:,:)          => null()
    real(real64), pointer, public :: zbar_3d_n_r8(:,:)       => null()
    real(real64), pointer, public :: area_r8(:,:)            => null()
    real(real64), pointer, public :: areasvol_r8(:,:)        => null()
    real(real64), pointer, public :: elem_area_r8(:)         => null()
    real(real64), pointer, public :: geo_coord_nod2D_r8(:,:) => null()

#if defined(USE_SINGLE_PRECISION)
    real(real64), allocatable, target :: z_3d_n_c(:,:), zbar_3d_n_c(:,:), area_c(:,:), areasvol_c(:,:)
    real(real64), allocatable, target :: elem_area_c(:), geo_coord_nod2D_c(:,:)
    ! tracers as REcoM reads and writes them, (level, node, tracer)
    real(real64), allocatable, target :: tracers_c(:,:,:)
#endif

contains

!-------------------------------------------------------------------------------
! point the geometry views to the current mesh; the vertical coordinates change
! every step (ALE), so call this before REcoM uses them in a step
subroutine recom_views_mesh(mesh)
    type(t_mesh), intent(in), target :: mesh
#if defined(USE_SINGLE_PRECISION)
    call copy2(mesh%z_3d_n,          z_3d_n_c)
    call copy2(mesh%zbar_3d_n,       zbar_3d_n_c)
    call copy2(mesh%area,            area_c)
    call copy2(mesh%areasvol,        areasvol_c)
    call copy2(mesh%geo_coord_nod2D, geo_coord_nod2D_c)
    if (.not. allocated(elem_area_c)) allocate(elem_area_c(size(mesh%elem_area)))
    elem_area_c = real(mesh%elem_area, real64)
    z_3d_n_r8          => z_3d_n_c
    zbar_3d_n_r8       => zbar_3d_n_c
    area_r8            => area_c
    areasvol_r8        => areasvol_c
    elem_area_r8       => elem_area_c
    geo_coord_nod2D_r8 => geo_coord_nod2D_c
#else
    z_3d_n_r8          => mesh%z_3d_n
    zbar_3d_n_r8       => mesh%zbar_3d_n
    area_r8            => mesh%area
    areasvol_r8        => mesh%areasvol
    elem_area_r8       => mesh%elem_area
    geo_coord_nod2D_r8 => mesh%geo_coord_nod2D
#endif
end subroutine recom_views_mesh

#if defined(USE_SINGLE_PRECISION)
subroutine copy2(src, dst)
    real(kind=WP),             intent(in)    :: src(:,:)
    real(real64), allocatable, intent(inout) :: dst(:,:)
    if (.not. allocated(dst)) allocate(dst(size(src,1), size(src,2)))
    dst = real(src, real64)
end subroutine copy2
#endif

!-------------------------------------------------------------------------------
! point info%data_pointers (allocated by the caller) to the tracers as REcoM sees
! them; in single precision this copies the current WP tracers into real64. With
! only, just the tracers with these IDs are copied and pointed to (the others stay
! null), for callers that read a few tracers
subroutine recom_tracer_pointers(tracers, info, only)
    type(t_tracer),          intent(inout), target :: tracers
    type(tracers_info_type), intent(inout)      :: info
    integer, optional,       intent(in)         :: only(:)
    integer :: n
#if defined(USE_SINGLE_PRECISION)
    call alloc_tracers_c(tracers)
    do n = 1, tracers%num_tracers
        if (present(only)) then
            if (.not. any(only == tracers%data(n)%ID)) cycle
        end if
        tracers_c(:,:,n) = real(tracers%data(n)%values, real64)
        info%data_pointers(n)%tracer_data => tracers_c(:,:,n)
    end do
#else
    do n = 1, tracers%num_tracers
        info%data_pointers(n)%tracer_data => tracers%data(n)%values
    end do
#endif
end subroutine recom_tracer_pointers

!-------------------------------------------------------------------------------
! tracer n as REcoM reads it (read only), refreshed from the WP tracer
function recom_tracer_r8(tracers, n) result(p)
    type(t_tracer), intent(in), target :: tracers
    integer,        intent(in)         :: n
    real(real64),   pointer            :: p(:,:)
#if defined(USE_SINGLE_PRECISION)
    call alloc_tracers_c(tracers)
    tracers_c(:,:,n) = real(tracers%data(n)%values, real64)
    p => tracers_c(:,:,n)
#else
    p => tracers%data(n)%values
#endif
end function recom_tracer_r8

#if defined(USE_SINGLE_PRECISION)
subroutine alloc_tracers_c(tracers)
    type(t_tracer), intent(in) :: tracers
    if (.not. allocated(tracers_c)) &
        allocate(tracers_c(size(tracers%data(1)%values,1), size(tracers%data(1)%values,2), tracers%num_tracers))
end subroutine alloc_tracers_c
#endif

!-------------------------------------------------------------------------------
! after REcoM has changed the tracers through recom_tracer_pointers: round them
! back into the WP tracers (every node, so that the halo REcoM exchanged in real64
! stays consistent between ranks). With carry, the content REcoM added below the
! WP resolution, (X - WP(X))*hnode on the wet levels of the own nodes, goes into
! the tracer carry, where the next tracer update (diff_tracers_ale) picks it up.
! Nothing to do in double precision, where REcoM wrote into the tracers directly.
subroutine recom_tracers_back(tracers, mesh, partit, carry)
    type(t_tracer), intent(inout), target :: tracers
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(in),    target :: partit
    logical,        intent(in)            :: carry
#if defined(USE_SINGLE_PRECISION)
    integer       :: tr, node, nz
    real(kind=WP) :: v

    do tr = 1, tracers%num_tracers
        if (carry) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, v)
            do node = 1, partit%myDim_nod2D
                do nz = mesh%ulevels_nod2D(node), mesh%nlevels_nod2D(node)-1
                    v = real(tracers_c(nz,node,tr), WP)
                    tracers%work%carry(nz,node,tr) = tracers%work%carry(nz,node,tr) &
                        + real((tracers_c(nz,node,tr) - real(v, real64))*mesh%hnode(nz,node), WP)
                end do
            end do
!$OMP END PARALLEL DO
        end if
        tracers%data(tr)%values = real(tracers_c(:,:,tr), WP)
    end do
#endif
end subroutine recom_tracers_back

!-------------------------------------------------------------------------------
! the node halo datatypes for REcoM's real64 fields
subroutine recom_mpitypes(partit, s2, r2, s3, r3)
    type(t_partit), intent(in), target :: partit
    integer, pointer :: s2(:), r2(:), s3(:,:,:), r3(:,:,:)
#if defined(USE_SINGLE_PRECISION)
    s2 => partit%s_mpitype_nod2D_r8
    r2 => partit%r_mpitype_nod2D_r8
    s3 => partit%s_mpitype_nod3D_r8
    r3 => partit%r_mpitype_nod3D_r8
#else
    s2 => partit%s_mpitype_nod2D
    r2 => partit%r_mpitype_nod2D
    s3 => partit%s_mpitype_nod3D
    r3 => partit%r_mpitype_nod3D
#endif
end subroutine recom_mpitypes
#endif
end module recom_views
