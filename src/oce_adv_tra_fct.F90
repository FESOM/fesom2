module oce_adv_tra_fct_interfaces
  interface
    subroutine oce_adv_tra_fct_init(mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh), intent(in), target  :: mesh
    end subroutine

    subroutine oce_tra_adv_fct(dttf_h, dttf_v, ttf, lo, adf_h, adf_v, mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP), intent(inout)      :: dttf_h(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: dttf_v(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: lo (mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: adf_h(mesh%nl-1, myDim_edge2D)
      real(kind=WP), intent(inout)      :: adf_v(mesh%nl,  myDim_nod2D)
    end subroutine
 end interface
end module
!
!
!===============================================================================
subroutine oce_adv_tra_fct_init(mesh)
    use MOD_MESH
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use profiling_nvtx
    implicit none
    integer                  :: my_size
    type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

    call nvtxStartRange("oce_adv_tra_fct_init")
! Copy mesh information to gpu once at initialization
!$acc enter data copyin(nlevels_nod2D,ulevels_nod2D,nl,nlevels,ulevels,elem2D_nodes,nod_in_elem2D,nod_in_elem2D_num,&
!$acc& edges,edge_tri,area,edge_cross_dxdy,edge_dxdy,elem_cos,elem_area,nlevels_nod2D_min,ulevels_nod2D_max)

    my_size=myDim_nod2D+eDim_nod2D
    allocate(fct_LO(nl-1, my_size))        ! Low-order solution
    allocate(adv_flux_hor(nl-1,myDim_edge2D)) ! antidiffusive hor. contributions / from edges
    allocate(adv_flux_ver(nl, myDim_nod2D))   ! antidiffusive ver. fluxes / from nodes

    allocate(fct_ttf_max(nl-1, my_size),fct_ttf_min(nl-1, my_size))
    allocate(fct_plus(nl-1, my_size),fct_minus(nl-1, my_size))
    ! Initialize with zeros:
    fct_LO=0.0_WP
    adv_flux_hor=0.0_WP
    adv_flux_ver=0.0_WP
    fct_ttf_max=0.0_WP
    fct_ttf_min=0.0_WP
    fct_plus=0.0_WP
    fct_minus=0.0_WP
! Allocate gpu arrays for fct kernels.
!$acc enter data create(fct_ttf_max,fct_ttf_min,UV_rhs,fct_plus,fct_minus, fct_LO, adv_flux_hor, adv_flux_ver)
    
    if (mype==0) write(*,*) 'FCT is initialized'

    call nvtxEndRange
end subroutine oce_adv_tra_fct_init
!
!
!===============================================================================
subroutine oce_tra_adv_fct(dttf_h, dttf_v, ttf, lo, adf_h, adf_v, mesh)
    !
    ! 3D Flux Corrected Transport scheme
    ! Limits antidiffusive fluxes==the difference in flux HO-LO
    ! LO ==Low-order  (first-order upwind)
    ! HO ==High-order (3rd/4th order gradient reconstruction method)
    ! Adds limited fluxes to the LO solution
    use MOD_MESH
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_CONFIG
    use g_comm_auto
    use openacc_params
    use profiling_nvtx
    implicit none
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(inout)      :: dttf_h(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: dttf_v(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: lo (mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: adf_h(mesh%nl-1, myDim_edge2D)
    real(kind=WP), intent(inout)      :: adf_v(mesh%nl,  myDim_nod2D)
    integer                           :: n, nz, k, elem, enodes(3), num, el(2), nl1, nl2, nu1, nu2, nl12, nu12, edge
    real(kind=WP)                     :: flux, ae,tvert_max(mesh%nl-1),tvert_min(mesh%nl-1)
    real(kind=WP)                     :: flux_eps=1e-16
    real(kind=WP)                     :: bignumber=1e3
    integer                           :: vlimit=1

#include "associate_mesh.h"

    call nvtxStartRange("oce_adv_tra_fct")
    ! ------------      --------------------------------------------------------------
    ! ttf is the tracer field on step n
    ! del_ttf is the increment
    ! vlimit sets the version of limiting, see below
    ! --------------------------------------------------------------------------
    !___________________________________________________________________________
    ! a1. max, min between old solution and updated low-order solution per node
    ! Double loop over blocks then threads. Copy LO array, rest is either static mesh info or copied in 
    ! calling routine do_oce_adv_tra
    call nvtxStartRange("loop a1")
    !$acc parallel loop gang present(fct_ttf_max,fct_ttf_min,LO,ttf,nlevels_nod2D,ulevels_nod2D)&
    !$acc& private(nl1,nu1) vector_length(z_vector_length)
    do n=1,myDim_nod2D + edim_nod2d
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        !$acc loop vector
        do nz=nu1, nl1-1
            fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
            fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
        end do
    end do
    call nvtxEndRange

    !___________________________________________________________________________
    ! a2. Admissible increments on elements
    !     (only layers below the first and above the last layer)
    !     look for max, min bounds for each element --> UV_rhs here auxiliary array

    ! Double loop over blocks then threads. All arrays are either static mesh info or copied in
    ! calling routine do_oce_adv_tra, enodes are private array
    call nvtxStartRange("loop a2")
    !$acc parallel loop gang present(UV_rhs,nl,ulevels,nlevels,elem2D_nodes,fct_ttf_min,fct_ttf_max) &
    !$acc& private(enodes,nu1,nl1) vector_length(z_vector_length)
    do elem=1, myDim_elem2D
        enodes=elem2D_nodes(:,elem)
        nu1 = ulevels(elem)
        nl1 = nlevels(elem)
        !$acc loop vector
        do nz=nu1, nl1-1
            UV_rhs(1,nz,elem)=maxval(fct_ttf_max(nz,enodes))
            UV_rhs(2,nz,elem)=minval(fct_ttf_min(nz,enodes))
        end do
        if (nl1<=nl-1) then
            !$acc loop vector
            do nz=nl1,nl-1
                UV_rhs(1,nz,elem)=-bignumber
                UV_rhs(2,nz,elem)= bignumber
            end do
        endif
    end do ! --> do elem=1, myDim_elem2D
    call nvtxEndRange

    !___________________________________________________________________________
    ! a3. Bounds on clusters and admissible increments
    ! Vertical1: In this version we look at the bounds on the clusters
    !            above and below, which leaves wide bounds because typically
    !            vertical gradients are larger.

    ! Double loop over blocks then threads. All arrays are either static mesh info or copied in
    ! calling routine do_oce_adv_tra, enodes are private array
    if(vlimit==1) then
        !Horizontal
        call nvtxStartRange("loop a1 vlimit=1")
        !$acc parallel loop gang present(ulevels_nod2D,nlevels_nod2D,nod_in_elem2D,nod_in_elem2D_num,UV_rhs,&
        !$acc& fct_ttf_min,fct_ttf_max,LO) private(tvert_min,tvert_max)&
        !$acc& private(nu1,nl1) vector_length(z_vector_length)
        do n=1, myDim_nod2D
            nu1 = ulevels_nod2D(n)
            nl1 = nlevels_nod2D(n)

            !___________________________________________________________________
            !$acc loop vector
            do nz=nu1,nl1-1
                ! max,min horizontal bound in cluster around node n in every
                ! vertical layer
                ! nod_in_elem2D     --> elem indices of which node n is surrounded
                ! nod_in_elem2D_num --> max number of surrounded elem
                tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
                tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
            end do

            !___________________________________________________________________
            ! calc max,min increment of surface layer with respect to low order
            ! solution
            fct_ttf_max(nu1,n)=tvert_max(nu1)-LO(nu1,n)
            fct_ttf_min(nu1,n)=tvert_min(nu1)-LO(nu1,n)

            ! calc max,min increment from nz-1:nz+1 with respect to low order
            ! solution at layer nz
            !$acc loop vector
            do nz=nu1+1,nl1-2
                fct_ttf_max(nz,n)=maxval(tvert_max(nz-1:nz+1))-LO(nz,n)
                fct_ttf_min(nz,n)=minval(tvert_min(nz-1:nz+1))-LO(nz,n)
            end do
            ! calc max,min increment of bottom layer -1 with respect to low order
            ! solution
            nz=nl1-1
            fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
            fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)
        end do
        call nvtxEndRange
    end if

    !___________________________________________________________________________
    ! Vertical2: Similar to the version above, but the vertical bounds are more
    ! local
    if(vlimit==2) then
        call nvtxStartRange("loop a2 vlimit=2")
        !$acc parallel loop gang present(ulevels_nod2D,nlevels_nod2D,nod_in_elem2D,nod_in_elem2D_num,UV_rhs,&
        !$acc& fct_ttf_min,fct_ttf_max,LO) private(nu1,nl1,tvert_min,tvert_max)&
        !$acc& vector_length(z_vector_length)
        do n=1, myDim_nod2D
            nu1 = ulevels_nod2D(n)
            nl1 = nlevels_nod2D(n)
            !$acc loop vector
            do nz=nu1,nl1-1
                tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
                tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
            end do
            !$acc loop vector
            do nz=nu1+1, nl1-2
                tvert_max(nz)=max(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
                tvert_min(nz)=min(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
            end do
            !$acc loop vector
            do nz=nu1,nl1-1
                fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
                fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)
            end do
        end do
        call nvtxEndRange
    end if

    !___________________________________________________________________________
    ! Vertical3: Vertical bounds are taken into account only if they are narrower than the
    !            horizontal ones
    if(vlimit==3) then
        call nvtxStartRange("loop a3 vlimit=3")
        !$acc parallel loop gang present(ulevels_nod2D,nlevels_nod2D,nod_in_elem2D,nod_in_elem2D_num,UV_rhs,&
        !$acc& fct_ttf_min,fct_ttf_max,LO) private(nu1,nl1,tvert_min,tvert_max)&
        !$acc& vector_length(z_vector_length)
        do n=1, myDim_nod2D
            nu1 = ulevels_nod2D(n)
            nl1 = nlevels_nod2D(n)
            !$acc loop vector
            do nz=nu1, nl1-1
                tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
                tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
            end do
            !$acc loop vector
            do nz=nu1+1, nl1-2
                tvert_max(nz)=min(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
                tvert_min(nz)=max(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
            end do
            !$acc loop vector
            do nz=nu1, nl1-1
                fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
                fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)
            end do
        end do
        call nvtxEndRange
    end if

    !___________________________________________________________________________
    ! b1. Split positive and negative antidiffusive contributions
    ! --> sum all positive (fct_plus), negative (fct_minus) antidiffusive
    !     horizontal element and vertical node contribution to node n and layer nz
    !     see. R. Löhner et al. "finite element flux corrected transport (FEM-FCT)
    !     for the euler and navier stoke equation
    call nvtxStartRange("loop b1")

    call nvtxStartRange("double loop")
    ! Double loop over blocks then threads. All arrays are either static mesh info or copied in
    ! calling routine do_oce_adv_tra, enodes are private array
    !$acc parallel loop gang present(ulevels_nod2D,nlevels_nod2D,fct_plus,fct_minus)&
    !$acc& private(nu1,nl1) vector_length(z_vector_length)
    do n=1, myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        !$acc loop vector
        do nz=nu1,nl1-1
            fct_plus(nz,n)=0._WP
            fct_minus(nz,n)=0._WP
        end do
    end do
    call nvtxEndRange

    !Vertical
    call nvtxStartRange("vertical")

    ! Wait for stream copying adf_v onto GPU, then perform kernel
    !$acc parallel loop gang present(ulevels_nod2D,nlevels_nod2D,fct_plus,fct_minus,adf_v)&
    !$acc& private(nu1,nl1) vector_length(z_vector_length) async(stream_ver_adv_tra)
    do n=1, myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        !$acc loop vector
        do nz=nu1,nl1-1
!             fct_plus(nz,n)=fct_plus(nz,n)+ &
!                             (max(0.0_WP,adf_v(nz,n))+max(0.0_WP,-adf_v(nz+1,n))) &
!                             /hnode(nz,n)
!             fct_minus(nz,n)=fct_minus(nz,n)+ &
!                             (min(0.0_WP,adf_v(nz,n))+min(0.0_WP,-adf_v(nz+1,n))) &
!                             /hnode(nz,n)
            !$acc atomic
            fct_plus(nz,n) =fct_plus(nz,n) +(max(0.0_WP,adf_v(nz,n))+max(0.0_WP,-adf_v(nz+1,n)))
            !$acc atomic
            fct_minus(nz,n)=fct_minus(nz,n)+(min(0.0_WP,adf_v(nz,n))+min(0.0_WP,-adf_v(nz+1,n)))
        end do
    end do
    call nvtxEndRange

    !Horizontal
    call nvtxStartRange("horizontal")

    ! Wait for stream copying adf_h onto GPU, then perform kernel
    !$acc parallel loop gang present(nlevels,ulevels,edges,edge_tri,fct_plus,fct_minus,adf_h)&
    !$acc& private(enodes,el,nl1,nu1,nl2,nu2,nl12,nu12)&
    !$acc& vector_length(z_vector_length) async(stream_hor_adv_tra)
    do edge=1, myDim_edge2D
        enodes(1:2)=edges(:,edge)
        el=edge_tri(:,edge)
        nl1=nlevels(el(1))-1
        nu1=ulevels(el(1))
        nl2=0
        nu2=0
        if(el(2)>0) then
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
        end if

        nl12 = max(nl1,nl2)
        nu12 = nu1
        if (nu2>0) nu12 = min(nu1,nu2)
        
        !$acc loop vector
        do nz=nu12, nl12
            !$acc atomic
            fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0_WP, adf_h(nz,edge))
            !$acc atomic
            fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0_WP, adf_h(nz,edge))
            !$acc atomic
            fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0_WP,-adf_h(nz,edge))
            !$acc atomic
            fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0_WP,-adf_h(nz,edge))
        end do
    end do
    call nvtxEndRange
    
    !$acc wait(stream_ver_adv_tra)
    !$acc wait(stream_hor_adv_tra)
    call nvtxEndRange
    !___________________________________________________________________________
    ! b2. Limiting factors

    call nvtxStartRange("loop b2")
    ! Double loop over blocks then threads. Transfer fct_plus, fct_minus to cpu for halo exchange
    !$acc parallel loop gang present(nlevels_nod2D,ulevels_nod2D,fct_plus,fct_minus,fct_ttf_min,fct_ttf_max,area)&
    !$acc& private(nu1,nl1,flux)&
    !$acc& vector_length(z_vector_length)
    do n=1,myDim_nod2D
        nu1=ulevels_nod2D(n)
        nl1=nlevels_nod2D(n)
        !$acc loop vector
        do nz=nu1,nl1-1
            flux=fct_plus(nz,n)*dt/area(nz,n)+flux_eps
            fct_plus(nz,n)=min(1.0_WP,fct_ttf_max(nz,n)/flux)
            flux=fct_minus(nz,n)*dt/area(nz,n)-flux_eps
            fct_minus(nz,n)=min(1.0_WP,fct_ttf_min(nz,n)/flux)
        end do
    end do
    
    !$acc update self(fct_plus,fct_minus)
    call nvtxEndRange

    ! fct_minus and fct_plus must be known to neighbouring PE
    call exchange_nod(fct_plus, fct_minus)

    !___________________________________________________________________________
    ! b3. Limiting
    call nvtxStartRange("loop b3")
    !Vertical
    call nvtxStartRange("vertival")

    ! Double loop over blocks then threads. Transfer result adf_v back to cpu.
    ! Schedule kernel asynchronously on stream 3, because next one is independent

    !$acc parallel loop gang present(ulevels_nod2D,nlevels_nod2D,fct_plus,fct_minus,adf_v)&
    !$acc& private(nz,ae,flux,nu1,nl1)&
    !$acc& vector_length(z_vector_length) async(stream_ver_adv_tra)
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
        !$acc loop vector
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
    ! the bottom flux is always zero
    end do

    call exchange_nod_end  ! fct_plus, fct_minus

    !$acc update device(fct_plus,fct_minus)
    call nvtxEndRange

    !Horizontal
    call nvtxStartRange("horizontal")

    ! Double loop over blocks then threads. Transfer result adf_h back to cpu.
    ! Schedule kernel asynchronously on stream 4, because next one is independent

    !$acc parallel loop gang present(nlevels,ulevels,edges,edge_tri,fct_plus,fct_minus,adf_h)&
    !$acc& private(enodes,el,nl1,nl2,nu1,nu2,ae,flux,nl12,nu12) vector_length(z_vector_length) async(stream_hor_adv_tra)
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

        !$acc loop vector
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
    end do
    call nvtxEndRange

    call nvtxEndRange
    call nvtxEndRange
end subroutine oce_tra_adv_fct
