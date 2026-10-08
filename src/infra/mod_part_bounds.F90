module mod_part_bounds
    ! M2.12b: local owned/halo loop bounds for the partitioned advection kernels.
    !
    ! The advection kernels keep their proven explicit-shape array dummies and per-
    ! element arithmetic UNCHANGED (so the codegen — and thus the byte-match vs FESOM2,
    ! incl. any vectorised divide, LESSONS L29 — is preserved); the only multi-rank
    ! change is the LOOP BOUNDS (and added halo exchanges in the callers). owned_bounds
    ! returns those bounds from an OPTIONAL partit:
    !   - partit absent OR npes==1  -> the GLOBAL mesh counts (the proven 1-rank path is
    !     byte-for-byte unchanged; 1-rank callers simply omit partit).
    !   - npes>1                     -> the partition's owned (myDim) and owned+halo
    !     (myDim+eDim) dims.
    ! At 1-rank nNodL==nNodO==mesh%nod2D etc., so an "owned+halo" loop collapses to the
    ! global loop exactly as before.
    use mod_mesh,   only: t_mesh
    use mod_partit, only: t_partit
    implicit none
    private
    public :: owned_bounds, is_multirank, local_dims

contains

    pure subroutine local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, &
                               nElemO, nElemL, nElemF)
        ! M2.12c: the full set of LOCAL loop/allocation bounds the dynamics kernels need.
        !   nNodO  = owned nodes (myDim)            nNodL  = owned+halo nodes (myDim+eDim)
        !   nEdgeO = owned edges (myDim)            nEdgeL = owned+halo edges (myDim+eDim)
        !   nElemO = owned elements (myDim)         nElemL = owned+eDim elements
        !   nElemF = owned+eDim+eXDim elements (the full element halo).
        ! partit absent OR npes==1 -> all = the mesh GLOBAL counts (the proven 1-rank path
        ! is byte-for-byte unchanged); npes>1 -> the partition dims. Mirrors
        ! mod_mesh_areas::local_bounds but adds nElemL (visc U_c init = myDim+eDim).
        type(t_mesh),   intent(in)  :: mesh
        integer,        intent(out) :: nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF
        type(t_partit), intent(in), optional :: partit
        if (present(partit)) then
            if (partit%npes > 1) then
                nNodO  = partit%myDim_nod2D
                nNodL  = partit%myDim_nod2D  + partit%eDim_nod2D
                nEdgeO = partit%myDim_edge2D
                nEdgeL = partit%myDim_edge2D + partit%eDim_edge2D
                nElemO = partit%myDim_elem2D
                nElemL = partit%myDim_elem2D + partit%eDim_elem2D
                nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
                return
            end if
        end if
        nNodO  = mesh%nod2D;  nNodL  = mesh%nod2D
        nEdgeO = mesh%edge2D; nEdgeL = mesh%edge2D
        nElemO = mesh%elem2D; nElemL = mesh%elem2D; nElemF = mesh%elem2D
    end subroutine local_dims

    pure subroutine owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit, nNodX)
        ! nNodO  = owned nodes ; nNodL = owned+halo nodes
        ! nEdgeO = owned edges ; nElemO = owned elements
        ! nNodX  = owned+halo+double-halo nodes (myDim+eDim+eXDim_nod2D): the size of the
        !          fields refreshed by exchange_nod_full (= nNodL without dist_N/format.out)
        type(t_mesh),   intent(in)  :: mesh
        integer,        intent(out) :: nNodO, nNodL, nEdgeO, nElemO
        type(t_partit), intent(in), optional :: partit
        integer,        intent(out), optional :: nNodX
        if (present(partit)) then
            if (partit%npes > 1) then
                nNodO  = partit%myDim_nod2D
                nNodL  = partit%myDim_nod2D + partit%eDim_nod2D
                nEdgeO = partit%myDim_edge2D
                nElemO = partit%myDim_elem2D
                if (present(nNodX)) nNodX = nNodL + partit%eXDim_nod2D
                return
            end if
        end if
        nNodO  = mesh%nod2D;  nNodL  = mesh%nod2D
        nEdgeO = mesh%edge2D; nElemO = mesh%elem2D
        if (present(nNodX)) nNodX = mesh%nod2D
    end subroutine owned_bounds

    pure logical function is_multirank(partit)
        ! .true. only when a partit is present AND npes>1 (safe — no short-circuit
        ! reliance). Guards the halo exchanges so 1-rank callers never exchange.
        type(t_partit), intent(in), optional :: partit
        is_multirank = .false.
        if (present(partit)) is_multirank = (partit%npes > 1)
    end function is_multirank

end module mod_part_bounds
