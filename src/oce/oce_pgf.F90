module oce_pgf
    ! Hydrostatic pressure gradient force (PGF) for which_ALE='linfs' with full
    ! cells. Transcribed from FESOM2 v2.7.3 oce_ale_pressure_bv.F90:
    !   pressure_force_4_linfs          (507-569)  dispatcher: branches on
    !                                              use_partial_cell/use_cavity/which_pgf
    !   pressure_force_4_linfs_fullcell (575-614)  the full-cell PGF kernel (ported)
    !
    ! M2.2 SCOPE — consumes the M2.1 hpressure (top-down hydrostatic integral) and
    ! contracts it with gradient_sca (the per-element linear shape-function
    ! derivatives, geometry-gated) into the element PGF pgf_x/pgf_y — the horizontal
    ! pressure term of the momentum RHS. NO hydrostatic integration here (that lives
    ! in oce_pressure_bv / M2.1). The dispatcher's non-fullcell branches (nemo,
    ! shchepetkin, cubicspline, easypgf, cavity, partial cells) are v1-dropped and
    ! UNGATED on pi (use_partial_cell=.false., no cavity ⇒ the fullcell branch); the
    ! dispatcher itself is pure namelist branching (no numerics) so it is not ported.
    !
    !   pgf_x(nz,elem) = Σ_k gradient_sca(k,  elem) * hpressure(nz, elnodes(k)) / density_0
    !   pgf_y(nz,elem) = Σ_k gradient_sca(4+k,elem) * hpressure(nz, elnodes(k)) / density_0
    !   (k = 1..4; a triangle's slot 4 repeats slot 1 with gradient 0)
    !
    ! BIT-IDENTITY NOTES:
    !  - density_0 is a RUNTIME divisor, but byte-identical on both sides (= 1030),
    !    so the -no-prec-div reciprocal matches (L7/L10). Transcribe the expression
    !    VERBATIM: the /density_0 is INSIDE the sum (each of the 3 products divided,
    !    then summed) — NOT sum(...)/density_0.
    !  - gradient_sca and the elem2D_nodes order are geometry-gate-proven (L9); the
    !    same gradient_sca contraction already byte-matched in oce_tracer_grad (M1.1).
    !  - Hybrid meshes: unified 4-slot formulas (gradient_sca 1:4 / 5:8, element means
    !    with the kite weights elem_w), triangles equal to the former 3-slot code to
    !    roundoff, not bitwise.
    !
    ! 1-rank only: the element loop runs to elem2D (== myDim_elem2D at 1-rank);
    ! pgf_x/pgf_y are element-local (no halo dependency, FESOM2 does not exchange
    ! them). The outputs are intent(inout) and NOT zeroed inside (FESOM2 leaves
    ! below-bottom entries at their allocation-zero) — the CALLER pre-zeros them.
    use mod_precision, only: WP
    use mod_mesh,      only: t_mesh, MAX_NV
    use mod_constants, only: density_0, g
    use mod_partit,    only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: pressure_force_4_linfs_fullcell
    public :: pressure_force_4_zxxxx_shchepetkin

contains

    !===========================================================================
    subroutine pressure_force_4_linfs_fullcell(hpressure, mesh, pgf_x, pgf_y, partit)
        ! hpressure: (nl, nod2D) in (M2.1 output; only levels 1..nlevels-1 read).
        ! pgf_x/pgf_y: (nl-1, elem2D) inout; the caller pre-zeros them (see header).
        ! M2.12c: optional partit -> owned element loop (FESOM2 :593 do elem=1,
        ! myDim_elem2D). pgf is element-local (no halo dependency, FESOM2 does not
        ! exchange it); compute_vel_rhs reads it at owned elements only.
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: hpressure(mesh%nl, mesh%nod2D)
        real(kind=WP), intent(inout) :: pgf_x(mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: pgf_y(mesh%nl-1, mesh%elem2D)
        type(t_partit), intent(in), optional :: partit
        integer :: elem, elnodes(MAX_NV), nle, ule, nlz
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        !_______________________________________________________________________
        ! loop over elements
        do elem=1, nElemO
            !___________________________________________________________________
            ! number of levels at elem
            nle = mesh%nlevels(elem)-1
            ule = mesh%ulevels(elem)

            !___________________________________________________________________
            ! node indices of elem
            elnodes = mesh%elem2D_nodes(:, elem)

            !___________________________________________________________________
            ! loop over mid-depth levels to calculate the pressure gradient force
            ! (pgf) --> from top to bottom
            do nlz=ule,nle
                pgf_x(nlz,elem) = sum(mesh%gradient_sca(1:4,elem)*hpressure(nlz,elnodes)/density_0)
                pgf_y(nlz,elem) = sum(mesh%gradient_sca(5:8,elem)*hpressure(nlz,elnodes)/density_0)
            end do
        end do
    end subroutine pressure_force_4_linfs_fullcell

    !===========================================================================
    subroutine pressure_force_4_zxxxx_shchepetkin(density_m_rho0, mesh, pgf_x, pgf_y, partit)
        ! Shchepetkin & McWilliams (2003) density-Jacobian PGF for the full free
        ! surface (which_ALE='zlevel'/'zstar', which_pgf='shchepetkin'). Transcribed
        ! VERBATIM from FESOM2 v2.7.3 oce_ale_pressure_bv.F90:2104-2339
        ! (subroutine pressure_force_4_zxxxx_shchepetkin). M6a-1.
        !
        ! Self-contained: integrates the density_m_rho0 gradient over the (evolving)
        ! Z_3d_n/helem column directly (eqn 1.7) -> NO hpressure dependency (unlike
        ! the linfs PGF). Computes the central second-order vertical density gradient
        ! drho_dz via a Newton interpolation polynomial on the non-equidistant Z_3d_n.
        !
        !   density_m_rho0: (nl-1, nod2D) in (the M2.1 pressure_bv output).
        !   pgf_x/pgf_y:    (nl-1, elem2D) inout; caller pre-zeros (kernel writes ule..nle only).
        !
        ! BIT-IDENTITY NOTES:
        !  - g/density_0 are runtime divisors, byte-identical on both sides (L7/L10).
        !  - Hybrid meshes: per-vertex arrays over the 4 slots (a triangle's slot 4 repeats
        !    slot 1), gradients over gradient_sca 1:4 / 5:8, and the element mean of the
        !    vertex drho_dz with the kite weights elem_w (FESOM2: sum(drho_dz)/3). FESOM2's
        !    `(/1.0,1.0,1.0/)*z_n(nlz)` is the scalar z_n(nlz) (1.0 exact).
        !  - L29 watch: the subsurface bulk loop does length-4 array divides on LOCAL
        !    arrays.
        !  - M6a-1: optional partit -> owned element loop; pgf is element-local (no
        !    halo dependency, FESOM2 does not exchange it).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: density_m_rho0(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: pgf_x(mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: pgf_y(mesh%nl-1, mesh%elem2D)
        type(t_partit), intent(in), optional :: partit
        integer :: elem, elnodes(MAX_NV), nle, ule, nlz, ni, idx(MAX_NV)
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: int_dp_dx(2), drho_dx, drho_dy, drho_dz(MAX_NV), dz_dx, dz_dy, aux_sum
        real(kind=WP) :: dx10(MAX_NV), dx20(MAX_NV), dx21(MAX_NV), df10(MAX_NV), df21(MAX_NV)
        real(kind=WP) :: zbar_n(mesh%nl), z_n(mesh%nl-1)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        !_______________________________________________________________________
        ! loop over elements
        do elem=1, nElemO
            ! nle...number of mid-depth levels at elem
            nle           = mesh%nlevels(elem)-1
            ule           = mesh%ulevels(elem)
            elnodes       = mesh%elem2D_nodes(:, elem)

            !___________________________________________________________________
            ! mid-depth element levels Z_n (from helem + zbar_e_bot)
            zbar_n        = 0.0_WP
            z_n           = 0.0_WP
            zbar_n(nle+1) = mesh%zbar_e_bot(elem)
            z_n(nle)      = zbar_n(nle+1) + mesh%helem(nle,elem)*0.5_WP
            do nlz=nle,ule+1,-1
                zbar_n(nlz) = zbar_n(nlz+1) + mesh%helem(nlz,elem)
                z_n(nlz-1)  = zbar_n(nlz)   + mesh%helem(nlz-1,elem)*0.5_WP
            end do
            zbar_n(ule)   = zbar_n(ule+1) + mesh%helem(ule,elem)

            !___________________________________________________________________
            ! surface layer pressure gradient
            nlz = ule
            idx = nlz - mesh%ulevels_nod2D(elnodes)
            do ni=1,MAX_NV
                if (idx(ni)==0) then
                    ! elemental surface index nlz is also surface index of vertice ni
                    dx10(ni)   = mesh%Z_3d_n(nlz+1,elnodes(ni))-mesh%Z_3d_n(nlz  ,elnodes(ni))
                    dx21(ni)   = mesh%Z_3d_n(nlz+2,elnodes(ni))-mesh%Z_3d_n(nlz+1,elnodes(ni))
                    dx20(ni)   = mesh%Z_3d_n(nlz+2,elnodes(ni))-mesh%Z_3d_n(nlz  ,elnodes(ni))
                    df10(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz+2,elnodes(ni))-density_m_rho0(nlz+1,elnodes(ni))
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                 *((z_n(nlz)-mesh%Z_3d_n(nlz+1,elnodes(ni))) + &
                                   (z_n(nlz)-mesh%Z_3d_n(nlz  ,elnodes(ni))))
                else
                    ! elemental surface index nlz is deeper than surface index of vertice ni
                    dx10(ni)   = mesh%Z_3d_n(nlz  ,elnodes(ni))-mesh%Z_3d_n(nlz-1,elnodes(ni))
                    dx21(ni)   = mesh%Z_3d_n(nlz+1,elnodes(ni))-mesh%Z_3d_n(nlz  ,elnodes(ni))
                    dx20(ni)   = mesh%Z_3d_n(nlz+1,elnodes(ni))-mesh%Z_3d_n(nlz-1,elnodes(ni))
                    df10(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                 *((z_n(nlz)-mesh%Z_3d_n(nlz  ,elnodes(ni))) + &
                                   (z_n(nlz)-mesh%Z_3d_n(nlz-1,elnodes(ni))))
                end if
            end do

            ! zonal surface pressure gradient
            drho_dx         = sum(mesh%gradient_sca(1:4,elem)*density_m_rho0(nlz,elnodes))
            dz_dx           = sum(mesh%gradient_sca(1:4,elem)*mesh%Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dx-sum(mesh%elem_w(:,elem)*drho_dz)*dz_dx)*mesh%helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(1)    = aux_sum

            ! meridional surface pressure gradient
            drho_dy         = sum(mesh%gradient_sca(5:8,elem)*density_m_rho0(nlz,elnodes))
            dz_dy           = sum(mesh%gradient_sca(5:8,elem)*mesh%Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dy-sum(mesh%elem_w(:,elem)*drho_dz)*dz_dy)*mesh%helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(2)    = aux_sum

            !___________________________________________________________________
            ! subsurface bulk layers (until one layer above the bottom)
            do nlz=ule+1,nle-1
                dx10            = mesh%Z_3d_n(nlz  ,elnodes)-mesh%Z_3d_n(nlz-1,elnodes)
                dx21            = mesh%Z_3d_n(nlz+1,elnodes)-mesh%Z_3d_n(nlz  ,elnodes)
                dx20            = mesh%Z_3d_n(nlz+1,elnodes)-mesh%Z_3d_n(nlz-1,elnodes)
                df10            = density_m_rho0(nlz  ,elnodes)-density_m_rho0(nlz-1,elnodes)
                df21            = density_m_rho0(nlz+1,elnodes)-density_m_rho0(nlz  ,elnodes)
                drho_dz         = df10/dx10 + (dx10*df21-dx21*df10)/(dx20*dx21*dx10)&
                                  *((z_n(nlz)-mesh%Z_3d_n(nlz,elnodes)) + &
                                    (z_n(nlz)-mesh%Z_3d_n(nlz-1,elnodes)))

                ! zonal bulk pressure gradient
                drho_dx         = sum(mesh%gradient_sca(1:4,elem)*density_m_rho0(nlz,elnodes))
                dz_dx           = sum(mesh%gradient_sca(1:4,elem)*mesh%Z_3d_n(nlz,elnodes))
                aux_sum         = (drho_dx-sum(mesh%elem_w(:,elem)*drho_dz)*dz_dx)*mesh%helem(nlz,elem)*g/density_0
                pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
                int_dp_dx(1)    = int_dp_dx(1) + aux_sum

                ! meridional bulk pressure gradient
                drho_dy         = sum(mesh%gradient_sca(5:8,elem)*density_m_rho0(nlz,elnodes))
                dz_dy           = sum(mesh%gradient_sca(5:8,elem)*mesh%Z_3d_n(nlz,elnodes))
                aux_sum         = (drho_dy-sum(mesh%elem_w(:,elem)*drho_dz)*dz_dy)*mesh%helem(nlz,elem)*g/density_0
                pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
                int_dp_dx(2)    = int_dp_dx(2) + aux_sum
            end do

            !___________________________________________________________________
            ! bottom layer pressure gradient
            nlz = nle
            idx = mesh%nlevels_nod2D(elnodes)-1 - nlz
            do ni=1,MAX_NV
                if (idx(ni)==0) then
                    ! elemental bottom index nlz is also bottom index of vertice ni
                    dx10(ni)   = mesh%Z_3d_n(nlz-1,elnodes(ni))-mesh%Z_3d_n(nlz-2,elnodes(ni))
                    dx21(ni)   = mesh%Z_3d_n(nlz  ,elnodes(ni))-mesh%Z_3d_n(nlz-1,elnodes(ni))
                    dx20(ni)   = mesh%Z_3d_n(nlz  ,elnodes(ni))-mesh%Z_3d_n(nlz-2,elnodes(ni))
                    df10(ni)   = density_m_rho0(nlz-1,elnodes(ni))-density_m_rho0(nlz-2,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                 *((z_n(nlz)-mesh%Z_3d_n(nlz-1,elnodes(ni))) + &
                                   (z_n(nlz)-mesh%Z_3d_n(nlz-2,elnodes(ni))))
                else
                    ! elemental bottom index nlz is shallower than bottom index of vertice ni
                    dx10(ni)   = mesh%Z_3d_n(nlz  ,elnodes(ni))-mesh%Z_3d_n(nlz-1,elnodes(ni))
                    dx21(ni)   = mesh%Z_3d_n(nlz+1,elnodes(ni))-mesh%Z_3d_n(nlz  ,elnodes(ni))
                    dx20(ni)   = mesh%Z_3d_n(nlz+1,elnodes(ni))-mesh%Z_3d_n(nlz-1,elnodes(ni))
                    df10(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                 *((z_n(nlz)-mesh%Z_3d_n(nlz  ,elnodes(ni))) + &
                                   (z_n(nlz)-mesh%Z_3d_n(nlz-1,elnodes(ni))))
                end if
            end do

            ! zonal bottom pressure gradient
            drho_dx         = sum(mesh%gradient_sca(1:4,elem)*density_m_rho0(nlz,elnodes))
            dz_dx           = sum(mesh%gradient_sca(1:4,elem)*(mesh%Z_3d_n(nlz,elnodes)))
            aux_sum         = (drho_dx-sum(mesh%elem_w(:,elem)*drho_dz)*dz_dx)*mesh%helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP

            ! meridional bottom pressure gradient
            drho_dy         = sum(mesh%gradient_sca(5:8,elem)*density_m_rho0(nlz,elnodes))
            dz_dy           = sum(mesh%gradient_sca(5:8,elem)*(mesh%Z_3d_n(nlz,elnodes)))
            aux_sum         = (drho_dy-sum(mesh%elem_w(:,elem)*drho_dz)*dz_dy)*mesh%helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
        end do
    end subroutine pressure_force_4_zxxxx_shchepetkin

end module oce_pgf
