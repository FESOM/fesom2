module oce_tracer_grad
    ! Elemental tracer gradient tr_xy, transcribed from FESOM2 v2.7.3
    ! oce_tracer_mod.F90:149-190 (SUBROUTINE tracer_gradient_elements).
    !
    ! tr_xy(1:2, nz, elem) is the horizontal gradient of the tracer field ttf,
    ! reconstructed per element from the Green-Gauss coefficients gradient_sca over the 4
    ! vertex slots (x 1:4, y 5:8; a triangle's slot 4 has coefficient 0; built in
    ! mod_mesh_areas). It feeds the MUSCL kernels (on the fly) and the horizontal diffusion.
    !
    ! Clean-architecture change vs FESOM2 (D7, plan "USE-globals -> explicit type
    ! arguments"): FESOM2 writes the module-global tr_xy (o_ARRAYS); here tr_xy is
    ! an explicit intent(out) argument. 1-rank only (myDim_elem2D == elem2D); the
    ! multi-rank loop bounds + halo exchange of tr_xy are deferred to M1.5/M2.12.
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh, MAX_NV
    use mod_partit,       only: t_partit
    use mod_part_bounds,  only: owned_bounds
    implicit none
    private
    public :: tracer_gradient_elements, tracer_gradient_z

contains

    subroutine tracer_gradient_z(ttf, tr_z, mesh, partit)
        ! Vertical tracer gradient tr_z(nz,n) = (ttf(nz-1,n)-ttf(nz,n))/dz, dz the mean of the
        ! two adjacent layer thicknesses (FESOM2 oce_tracer_mod.F90:211-251). Feeds the Redi
        ! K13/K23 isopycnal terms in diff_part_hor_redi (M4d). Loop owned+halo (nNodL) so tr_z
        ! is halo-valid (FESOM2 loops myDim+eDim then exchanges; here the owned+halo loop reads
        ! the already-valid ttf halo). tr_z=0 at the top/bottom interface.
        type(t_mesh),  intent(in)  :: mesh
        real(kind=WP), intent(in)  :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(out) :: tr_z(mesh%nl, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit
        integer :: n, nz, nzmin, nzmax
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: dz

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do n = 1, nNodL
            nzmax = mesh%nlevels_nod2D(n)
            nzmin = mesh%ulevels_nod2D(n)
            do nz = nzmin+1, nzmax-1
                dz = 0.5_WP*(mesh%hnode(nz-1,n)+mesh%hnode(nz,n))
                tr_z(nz, n) = (ttf(nz-1,n)-ttf(nz,n))/dz
            end do
            tr_z(nzmin, n) = 0.0_WP
            tr_z(nzmax, n) = 0.0_WP
        end do
    end subroutine tracer_gradient_z

    subroutine tracer_gradient_elements(ttf, tr_xy, mesh, partit)
        ! computes elemental gradient of tracer  (oce_tracer_mod.F90:149-190).
        ! M2.12b: optional partit -> loop over OWNED elements (the caller exchanges
        ! tr_xy to the full halo afterwards). 1-rank (partit absent) loops all elements.
        type(t_mesh),  intent(in)  :: mesh
        real(kind=WP), intent(in)  :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(out) :: tr_xy(2, mesh%nl-1, mesh%elem2D)
        type(t_partit), intent(in), optional :: partit
        integer :: elem, elnodes(MAX_NV), nz, nzmin, nzmax
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:, elem)
            nzmin = mesh%ulevels(elem)
            nzmax = mesh%nlevels(elem)
            do nz = nzmin, nzmax-1
                tr_xy(1, nz, elem) = sum(mesh%gradient_sca(1:4, elem) * ttf(nz, elnodes))
                tr_xy(2, nz, elem) = sum(mesh%gradient_sca(5:8, elem) * ttf(nz, elnodes))
            end do
        end do
    end subroutine tracer_gradient_elements

end module oce_tracer_grad
