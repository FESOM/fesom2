!---------------------------------------------------------------------------
!Implementation of mixed layer eddy restratification after Fox-Kemper et al., 2008, 2011
!Contains:
!  mle_add_gamma
!===========================================================================
module mle_interface
    interface
        subroutine mle_add_gamma(partit, mesh)
            use mod_mesh
            USE MOD_PARTIT
            USE MOD_PARSUP
            type(t_partit), intent(inout), target :: partit
            type(t_mesh),   intent(in),    target :: mesh
        end subroutine mle_add_gamma
    end interface
end module mle_interface
!
!===============================================================================
!
subroutine mle_add_gamma(partit, mesh)
    !  Adds the mixed layer eddy streamfunction to the GM one, before
    !  fer_gamma2vel converts both to a bolus velocity.
    !
    !  Gamma_mle = Ce * (ds/Lf) * H**2 * mu(z) * grad_h(b) / sqrt(f**2 + tau**-2)
    !  mu(z)     = (1-xi**2) * (1 + 5/21*xi**2),  xi = 2*(z-ztop)/H - 1
    !
    !  H       is the mixed layer depth, limited by mle_hmax
    !  grad_h(b) is the buoyancy gradient averaged over the mixed layer
    !  ds/Lf   is the grid scale over the frontal width, limited by mle_resscale_max.
    !          The grid-scale buoyancy gradient underestimates the frontal one by
    !          about Lf/ds, so this is a field on a variable resolution mesh and
    !          mesh_resolution(n) is used per node.
    !  Lf      = max(N*H/sqrt(f**2+tau**-2), mle_Lf_min), the mixed layer
    !          deformation radius with a floor. Using the regularised f makes the
    !          factor vanish towards the equator rather than diverge.
    !  H       is optionally passed through a one-sided running mean
    !          (mle_mld_decay_time): it deepens at once but decays slowly, so
    !          restratification continues after a convective event. Not carried
    !          in the restart; it respins within a few decay times.
    !  Ce, Lf, tau and the limits are namelist-tunable (oce_dyn)
    !
    !  GM and MLE do not double count: the ODM95 taper drives GM to zero wherever
    !  the neutral slope exceeds ODM95_Scr, which is the whole mixed layer front.
    !
    !  Gamma is parallel to grad_h(b), as in fer_solve_Gamma, whose interior limit
    !  is Gamma = K*grad_h(b)/N**2. mu vanishes at both ends of the range, so the
    !  boundary conditions match those of the GM Gamma.
    !
    !  Levels are restricted to ulevels_nod2D_max..nlevels_nod2D_min as in
    !  fer_solve_Gamma. This is required, not cosmetic: fer_gamma2vel forms
    !  u* = [Gamma(nz)-Gamma(nz+1)]/h, so the column integral of the bolus
    !  transport telescopes to Gamma(top)-Gamma(bottom) and vanishes only if
    !  Gamma is zero at both ends of the range the ELEMENT has. Where a node is
    !  deeper than its surrounding elements the residual forces the barotropic mode.
    !
    !  The column is then scaled so that max|d(Gamma)/dz| <= mle_ustar_max. Psi
    !  grows as H**2 and as 1/f, so without this it is unbounded at the equator.
    !
    ! Output: fer_gamma(1:2,:,:) incremented, mle_psi(:,:) diagnostic
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE o_PARAM
    USE o_ARRAYS, ONLY: sigma_xy, fer_gamma, MLD2, bvfreq, mle_psi, mle_hbar
    USE g_CONFIG
    use g_comm_auto
    IMPLICIT NONE

    type(t_partit), intent(inout), target  :: partit
    type(t_mesh),   intent(in),    target  :: mesh

    integer       :: n, nz, nzmin, nzmax, nzmin_s, nzmax_s, nlev_ml
    real(kind=WP) :: hml, dxlf, fdenom, prefac, r, xi, mu, zlev, wsum, thick, umax
    real(kind=WP) :: n2bar, nml, lfront
    real(kind=WP) :: ztop, zspan
    real(kind=WP) :: sigbar(2)
    real(kind=WP) :: zbar_n(mesh%nl), gcol(2, mesh%nl)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    if (.not. use_mle) return
    ! With GM off nothing else fills fer_gamma, so it has to start clean each step.
    if (.not. Fer_GM) fer_gamma = 0.0_WP
    ! Diagnostic: the MLE streamfunction on its own, since fer_uv holds the sum of
    ! the MLE and GM bolus velocities.
    mle_psi = 0.0_WP
    r = g / density_0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax, nlev_ml, hml, dxlf, &
!$OMP                                     fdenom, prefac, xi, mu, zlev, wsum, thick, &
!$OMP                                     umax, ztop, zspan, nzmin_s, nzmax_s, &
!$OMP                                     n2bar, nml, lfront, &
!$OMP                                     sigbar, zbar_n, gcol)
    DO n=1, myDim_nod2D
        nzmin = ulevels_nod2D(n)
        nzmax = nlevels_nod2D(n)
        if (nzmax - nzmin < 3) cycle

        !_______________________________________________________________________
        ! level depths below the local surface, positive downward, rebuilt from
        ! hnode_new exactly as fer_solve_Gamma does so the two stay consistent
        zbar_n(nzmin) = 0.0_WP
        do nz = nzmin, nzmax-1
            zbar_n(nz+1) = zbar_n(nz) + hnode_new(nz,n)
        end do

        !_______________________________________________________________________
        ! Level range every element containing this node has, as in fer_solve_Gamma.
        ! Required for the bolus transport to integrate to zero; see the header.
        nzmin_s = ulevels_nod2D_max(n)
        nzmax_s = nlevels_nod2D_min(n)
        if (nzmax_s - nzmin_s < 3) cycle

        !_______________________________________________________________________
        ! mixed layer depth, limited by mle_hmax: Psi grows as H**2 and a column
        ! convecting to the bottom is outside the regime this was derived for
        ztop = zbar_n(nzmin_s)
        hml  = min(abs(MLD2(n)), mle_hmax)
        hml  = min(hml, zbar_n(nzmax_s))
        if (mle_mld_decay_time > 0.0_WP) then
            mle_hbar(n) = max(hml, (dt*hml + mle_mld_decay_time*mle_hbar(n)) &
                                   / (dt + mle_mld_decay_time))
            hml = min(mle_hbar(n), zbar_n(nzmax_s))
        end if
        zspan = hml - ztop
        if (zspan <= zbar_n(nzmin_s+1) - ztop) cycle   ! thinner than one layer

        !_______________________________________________________________________
        ! mixed-layer average of the horizontal density gradient
        sigbar = 0.0_WP
        n2bar  = 0.0_WP
        wsum   = 0.0_WP
        nlev_ml = nzmin_s
        do nz = nzmin_s, nzmax_s-1
            if (zbar_n(nz) >= hml) exit
            thick    = min(zbar_n(nz+1), hml) - zbar_n(nz)
            sigbar   = sigbar + sigma_xy(1:2,nz,n) * thick
            n2bar    = n2bar + max(bvfreq(nz,n), 0.0_WP) * thick
            wsum     = wsum + thick
            nlev_ml  = nz+1
        end do
        if (wsum <= 0.0_WP) cycle
        sigbar = sigbar / wsum
        n2bar  = n2bar / wsum

        !_______________________________________________________________________
        ! resolution factor and the equatorial regularisation
        fdenom = sqrt(mesh%coriolis_node(n)**2 + 1.0_WP/(mle_tau*mle_tau))
        nml    = sqrt(max(n2bar, 0.0_WP))
        lfront = max(nml*zspan/fdenom, mle_Lf_min)
        dxlf   = min(mesh_resolution(n)/lfront, mle_resscale_max)
        prefac = mle_Ce * dxlf * zspan*zspan / fdenom

        !_______________________________________________________________________
        ! vertical structure, and add to the GM streamfunction.  Gamma is parallel
        ! to grad_h(b) = -(g/rho0)*sigma_xy, matching the GM sign (see header).
        gcol = 0.0_WP
        do nz = nzmin_s, min(nlev_ml, nzmax_s)
            zlev = zbar_n(nz)
            if (zlev > hml) exit
            ! xi runs -1 at the top of the range to +1 at the mixed layer base
            xi = 2.0_WP*(zlev - ztop)/zspan - 1.0_WP
            mu = (1.0_WP - xi*xi) * (1.0_WP + (5.0_WP/21.0_WP)*xi*xi)
            if (mu <= 0.0_WP) cycle
            gcol(1:2,nz) = prefac * mu * ( -r * sigbar(1:2) )
        end do

        !_______________________________________________________________________
        ! Limit the induced bolus velocity u* = d(Gamma)/dz. One factor for the
        ! whole column, so mu(z) is not distorted.
        umax = 0.0_WP
        do nz = nzmin_s, nzmax_s-1
            thick = max(hnode_new(nz,n), 1.0_WP)
            umax  = max(umax, sqrt( (gcol(1,nz)-gcol(1,nz+1))**2 &
                                  + (gcol(2,nz)-gcol(2,nz+1))**2 ) / thick)
        end do
        if (umax > mle_ustar_max) gcol = gcol * (mle_ustar_max/umax)

        do nz = nzmin_s, nzmax_s
            fer_gamma(1:2,nz,n) = fer_gamma(1:2,nz,n) + gcol(1:2,nz)
            mle_psi(nz,n) = sqrt(gcol(1,nz)*gcol(1,nz) + gcol(2,nz)*gcol(2,nz))
        end do
    END DO
!$OMP END PARALLEL DO

    call exchange_nod(fer_gamma, partit)
end subroutine mle_add_gamma
