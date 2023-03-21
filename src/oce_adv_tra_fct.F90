module oce_adv_tra_fct_interfaces
  interface
    subroutine oce_adv_tra_fct_init(twork, partit, mesh)
      use MOD_MESH
      use MOD_TRACER
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_mesh),  intent(in),    target        :: mesh
      type(t_partit),intent(inout), target        :: partit
      type(t_tracer_work), intent(inout), target  :: twork
    end subroutine

    subroutine oce_tra_adv_fct(dt, ttf, lo, adf_h, adf_v, fct_ttf_min, fct_ttf_max, fct_plus, fct_minus, AUX, partit, mesh)
      use MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      real(kind=WP), intent(in),    target :: dt
      type(t_partit),intent(inout), target :: partit
      type(t_mesh),  intent(in),    target :: mesh
      real(kind=WP), intent(inout)      :: fct_ttf_min(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(inout)      :: fct_ttf_max(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: lo (mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(inout)      :: adf_h(mesh%nl-1, partit%myDim_edge2D)
      real(kind=WP), intent(inout)      :: adf_v(mesh%nl,   partit%myDim_nod2D)
      real(kind=WP), intent(inout)      :: fct_plus(mesh%nl-1, partit%myDim_nod2D)
      real(kind=WP), intent(inout)      :: fct_minus(mesh%nl,  partit%myDim_nod2D)
      real(kind=WP), intent(inout)      :: AUX(:,:,:) !a large auxuary array
    end subroutine
 end interface
end module
!
!
!===============================================================================
subroutine oce_adv_tra_fct_init(twork, partit, mesh)
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    implicit none
    integer                                    :: my_size
    type(t_mesh),        intent(in) ,   target :: mesh
    type(t_partit),      intent(inout), target :: partit
    type(t_tracer_work), intent(inout), target :: twork
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    my_size=myDim_nod2D+eDim_nod2D
    allocate(twork%fct_LO(nl-1, my_size))        ! Low-order solution
    allocate(twork%adv_flux_hor(nl-1,partit%myDim_edge2D)) ! antidiffusive hor. contributions / from edges
    allocate(twork%adv_flux_ver(nl, partit%myDim_nod2D))   ! antidiffusive ver. fluxes / from nodes

    allocate(twork%fct_ttf_max(nl-1, my_size),twork%fct_ttf_min(nl-1, my_size))
    allocate(twork%fct_plus(nl-1, my_size),   twork%fct_minus(nl-1, my_size))
    ! Initialize with zeros:
    twork%fct_LO=0.0_WP
    twork%adv_flux_hor=0.0_WP
    twork%adv_flux_ver=0.0_WP
    twork%fct_ttf_max=0.0_WP
    twork%fct_ttf_min=0.0_WP
    twork%fct_plus=0.0_WP
    twork%fct_minus=0.0_WP

    if (mype==0) write(*,*) 'FCT is initialized'
end subroutine oce_adv_tra_fct_init

!
!
!===============================================================================
subroutine oce_tra_adv_fct(dt, ttf, lo, adf_h, adf_v, fct_ttf_min, fct_ttf_max, fct_plus, fct_minus, AUX, partit, mesh)
    !
    ! 3D Flux Corrected Transport scheme
    ! Limits antidiffusive fluxes==the difference in flux HO-LO
    ! LO ==Low-order  (first-order upwind)
    ! HO ==High-order (3rd/4th order gradient reconstruction method)
    ! Adds limited fluxes to the LO solution
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE g_comm_auto
    implicit none
    real(kind=WP), intent(in),    target :: dt
    type(t_mesh),  intent(in),    target :: mesh
    type(t_partit),intent(inout), target :: partit
    real(kind=WP), intent(inout)      :: fct_ttf_min(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: fct_ttf_max(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: lo (mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: adf_h(mesh%nl-1, partit%myDim_edge2D)
    real(kind=WP), intent(inout)      :: adf_v(mesh%nl,   partit%myDim_nod2D)
    real(kind=WP), intent(inout)      :: fct_plus (mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: fct_minus(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(inout)      :: AUX(:,:,:) !a large auxuary array, let us use twork%edge_up_dn_grad(1:4, 1:NL-2, 1:partit%myDim_edge2D) to save space
    integer                           :: n, nz, k, elem, enodes(3), num, el(2), nl1, nl2, nu1, nu2, nl12, nu12, edge
    real(kind=WP)                     :: flux, ae, tvert_max(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D), tvert_min(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP)                     :: flux_eps=1e-16
    real(kind=WP)                     :: bignumber=1e3
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, k, elem, enodes, num, el, nl1, nl2, nu1, nu2, nl12, nu12, edge, &
!$OMP                          flux, ae,tvert_max, tvert_min)
    ! --------------------------------------------------------------------------
    ! ttf is the tracer field on step n
    ! del_ttf is the increment
    ! vlimit sets the version of limiting, see below
    ! --------------------------------------------------------------------------
    !___________________________________________________________________________
    ! a1. max, min between old solution and updated low-order solution per node

    !$ACC DATA CREATE(tvert_max, tvert_min)

!$OMP DO
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1,myDim_nod2D + edim_nod2d
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        !$ACC LOOP VECTOR
        do nz=nu1, nl1-1
            fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
            fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
        end do
        !$ACC END LOOP
    end do
    !$ACC END PARALLEL LOOP
!$OMP END DO
    !___________________________________________________________________________
    ! a2. Admissible increments on elements
    !     (only layers below the first and above the last layer)
    !     look for max, min bounds for each element --> AUX here auxilary array
!$OMP DO
    !$ACC PARALLEL LOOP GANG PRIVATE(enodes) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do elem=1, myDim_elem2D
        enodes=elem2D_nodes(:,elem)
        nu1 = ulevels(elem)
        nl1 = nlevels(elem)
        !$ACC LOOP VECTOR
        do nz=nu1, nl1-1
            AUX(1,nz,elem)=max(fct_ttf_max(nz,enodes(1)), fct_ttf_max(nz,enodes(2)), fct_ttf_max(nz,enodes(3)))
            AUX(2,nz,elem)=min(fct_ttf_min(nz,enodes(1)), fct_ttf_min(nz,enodes(2)), fct_ttf_min(nz,enodes(3)))
        end do
        !$ACC END LOOP
        if (nl1<=nl-1) then
            !$ACC LOOP VECTOR
            do nz=nl1,nl-1
                AUX(1,nz,elem)=-bignumber
                AUX(2,nz,elem)= bignumber
            end do
            !$ACC END LOOP
        endif
    end do ! --> do elem=1, myDim_elem2D
    !$ACC END PARALLEL LOOP
!$OMP END DO
    !___________________________________________________________________________
    ! a3. Bounds on clusters and admissible increments
    ! Vertical1: In this version we look at the bounds on the clusters
    !            above and below, which leaves wide bounds because typically
    !            vertical gradients are larger.
        !Horizontal
!$OMP DO
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1, myDim_nod2D
       nu1 = ulevels_nod2D(n)
       nl1 = nlevels_nod2D(n)
       !___________________________________________________________________
       !$ACC LOOP VECTOR
       do nz=nu1,nl1-1
          ! max,min horizontal bound in cluster around node n in every
          ! vertical layer
          ! nod_in_elem2D     --> elem indices of which node n is surrounded
          ! nod_in_elem2D_num --> max number of surrounded elem
          tvert_max(nz, n) = AUX(1,nz, nod_in_elem2D(1, n))
          tvert_min(nz, n) = AUX(2,nz, nod_in_elem2D(1, n))
          !$ACC LOOP SEQ
          do elem=2,nod_in_elem2D_num(n)
              tvert_max(nz, n) = dmax1(tvert_max(nz, n), AUX(1,nz, nod_in_elem2D(elem,n)))
              tvert_min(nz, n) = dmin1(tvert_min(nz, n), AUX(2,nz, nod_in_elem2D(elem,n)))
          end do
          !$ACC END LOOP
       end do
       !$ACC END LOOP
    end do
    !$ACC END PARALLEL LOOP
!$OMP END DO

!$OMP DO
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1, myDim_nod2D
       nu1 = ulevels_nod2D(n)
       nl1 = nlevels_nod2D(n)
       !___________________________________________________________________
       ! calc max,min increment of surface layer with respect to low order
       ! solution
       fct_ttf_max(nu1,n)=tvert_max(nu1, n)-LO(nu1,n)
       fct_ttf_min(nu1,n)=tvert_min(nu1, n)-LO(nu1,n)
       ! calc max,min increment from nz-1:nz+1 with respect to low order
       ! solution at layer nz
       !$ACC LOOP VECTOR
       do nz=nu1+1,nl1-2
          fct_ttf_max(nz,n)=dmax1(tvert_max(nz-1, n), tvert_max(nz, n), tvert_max(nz+1, n))-LO(nz,n)
          fct_ttf_min(nz,n)=dmin1(tvert_min(nz-1, n), tvert_min(nz, n), tvert_min(nz+1, n))-LO(nz,n)
       end do
       !$ACC END LOOP
       ! calc max,min increment of bottom layer -1 with respect to low order
       ! solution
       nz=nl1-1
       fct_ttf_max(nz,n)=tvert_max(nz, n)-LO(nz,n)
       fct_ttf_min(nz,n)=tvert_min(nz, n)-LO(nz,n)
    end do
   !$ACC END PARALLEL LOOP
!$OMP END DO
    !___________________________________________________________________________
    ! b1. Split positive and negative antidiffusive contributions
    ! --> sum all positive (fct_plus), negative (fct_minus) antidiffusive
    !     horizontal element and vertical node contribution to node n and layer nz
    !     see. R. Löhner et al. "finite element flux corrected transport (FEM-FCT)
    !     for the euler and navier stoke equation
!$OMP DO
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1, myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
       !$ACC LOOP VECTOR
        do nz=nu1,nl1-1
            fct_plus(nz,n)=0._WP
            fct_minus(nz,n)=0._WP
        end do
       !$ACC END LOOP
    end do
    !$ACC END PARALLEL LOOP
!$OMP END DO

    !Vertical
!$OMP DO
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1, myDim_nod2D
       nu1 = ulevels_nod2D(n)
       nl1 = nlevels_nod2D(n)
       !$ACC LOOP VECTOR
       do nz=nu1,nl1-1
          fct_plus(nz,n) =fct_plus(nz,n) +(max(0.0_WP,adf_v(nz,n))+max(0.0_WP,-adf_v(nz+1,n)))
          fct_minus(nz,n)=fct_minus(nz,n)+(min(0.0_WP,adf_v(nz,n))+min(0.0_WP,-adf_v(nz+1,n)))
       end do
       !$ACC END LOOP
    end do
    !$ACC END PARALLEL LOOP
!$OMP END DO

!$OMP DO
    !Horizontal
#if !defined(DISABLE_OPENACC_ATOMICS)
    !$ACC PARALLEL LOOP GANG PRIVATE(enodes, el) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#else
    !$ACC UPDATE SELF(fct_plus, fct_minus, adf_h)
#endif
    do edge=1, myDim_edge2D
       enodes(1:2)=edges(:,edge)
       el=edge_tri(:,edge)
       nl1=nlevels(el(1))-1
       nu1=ulevels(el(1))
       nl2=0
       nu2=0
       if (el(2)>0) then
          nl2=nlevels(el(2))-1
          nu2=ulevels(el(2))
       end if

       nl12 = max(nl1,nl2)
       nu12 = nu1
       if (nu2>0) nu12 = min(nu1,nu2)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_set_lock(partit%plock(enodes(1)))
#else
!$OMP ORDERED
#endif
#if !defined(DISABLE_OPENACC_ATOMICS)
       !$ACC LOOP VECTOR
#endif
       do nz=nu12, nl12
#if !defined(DISABLE_OPENACC_ATOMICS)
          !$ACC ATOMIC UPDATE
#endif
          fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0_WP, adf_h(nz,edge))
#if !defined(DISABLE_OPENACC_ATOMICS)
          !$ACC ATOMIC UPDATE
#endif
          fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0_WP, adf_h(nz,edge))
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       end do
       call omp_unset_lock(partit%plock(enodes(1)))
       call omp_set_lock  (partit%plock(enodes(2)))
       do nz=nu12, nl12
#endif
#if !defined(DISABLE_OPENACC_ATOMICS)
          !$ACC ATOMIC UPDATE
#endif
          fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0_WP,-adf_h(nz,edge))
#if !defined(DISABLE_OPENACC_ATOMICS)
          !$ACC ATOMIC UPDATE
#endif
          fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0_WP,-adf_h(nz,edge))
       end do
#if !defined(DISABLE_OPENACC_ATOMICS)
       !$ACC END LOOP
#endif
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_unset_lock(partit%plock(enodes(2)))
#else
!$OMP END ORDERED
#endif
    end do
#if !defined(DISABLE_OPENACC_ATOMICS)
    !$ACC END PARALLEL LOOP
#else
    !$ACC UPDATE DEVICE(fct_plus, fct_minus)
#endif
!$OMP END DO

    !___________________________________________________________________________
    ! b2. Limiting factors
!$OMP DO
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1,myDim_nod2D
        nu1=ulevels_nod2D(n)
        nl1=nlevels_nod2D(n)
        !$ACC LOOP VECTOR
        do nz=nu1,nl1-1
            flux=fct_plus(nz,n)*dt/areasvol(nz,n)+flux_eps
            fct_plus(nz,n)=min(1.0_WP,fct_ttf_max(nz,n)/flux)
            flux=fct_minus(nz,n)*dt/areasvol(nz,n)-flux_eps
            fct_minus(nz,n)=min(1.0_WP,fct_ttf_min(nz,n)/flux)
        end do
        !$ACC END LOOP
    end do
    !$ACC END PARALLEL LOOP
!$OMP END DO

    ! fct_minus and fct_plus must be known to neighbouring PE
!$OMP MASTER
    call exchange_nod(fct_plus, fct_minus, partit, luse_g2g = .true.)
!$OMP END MASTER
!$OMP BARRIER

    !___________________________________________________________________________
    ! b3. Limiting
    !Vertical
!$OMP DO
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1, myDim_nod2D
        nu1=ulevels_nod2D(n)
        nl1=nlevels_nod2D(n)

        !_______________________________________________________________________
        nz=nu1
        ae=1.0_WP
        flux=adf_v(nz,n)
        if(flux>=0.0_WP) then
            ae=min(ae,fct_plus(nz,n))
        else
            ae=min(ae,fct_minus(nz,n))
        end if
        adf_v(nz,n)=ae*adf_v(nz,n)

        !_______________________________________________________________________
        !$ACC LOOP VECTOR
        do nz=nu1+1,nl1-1
            ae=1.0_WP
            flux=adf_v(nz,n)
            if(flux>=0._WP) then
                ae=min(ae,fct_minus(nz-1,n))
                ae=min(ae,fct_plus(nz,n))
            else
                ae=min(ae,fct_plus(nz-1,n))
                ae=min(ae,fct_minus(nz,n))
            end if
            adf_v(nz,n)=ae*adf_v(nz,n)
        end do
        !$ACC END LOOP
    ! the bottom flux is always zero
    end do
    !$ACC END PARALLEL LOOP
!$OMP END DO

    !Horizontal
!$OMP DO

    !$ACC PARALLEL LOOP GANG PRIVATE(enodes, el) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do edge=1, myDim_edge2D
        enodes(1:2)=edges(:,edge)
        el=edge_tri(:,edge)
        nu1=ulevels(el(1))
        nl1=nlevels(el(1))-1
        nl2=0
        nu2=0
        if(el(2)>0) then
            nu2=ulevels(el(2))
            nl2=nlevels(el(2))-1
        end if

        nl12 = max(nl1,nl2)
        nu12 = nu1
        if (nu2>0) nu12 = min(nu1,nu2)

        !$ACC LOOP VECTOR
        do nz=nu12, nl12
            ae=1.0_WP
            flux=adf_h(nz,edge)

            if(flux>=0._WP) then
                ae=min(ae,fct_plus(nz,enodes(1)))
                ae=min(ae,fct_minus(nz,enodes(2)))
            else
                ae=min(ae,fct_minus(nz,enodes(1)))
                ae=min(ae,fct_plus(nz,enodes(2)))
            endif

            adf_h(nz,edge)=ae*adf_h(nz,edge)
        end do
        !$ACC END LOOP
    end do
    !$ACC END PARALLEL LOOP
!$OMP END DO

!$ACC END DATA
!$OMP END PARALLEL
end subroutine oce_tra_adv_fct
