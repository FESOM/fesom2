!=======================================================================
!
! This submodule contains the subroutines 
! for the advection of sea ice tracers
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!
!=======================================================================

module icedrv_advection

    use icedrv_kinds
    use icedrv_constants
    use icedrv_system,      only: icedrv_system_abort

    implicit none

    public :: fct_init_icepack, tracer_advection_icepack

    private

    real(kind=dbl_kind), allocatable, dimension(:) :: &
         d_tr,        trl,                            &
         rhs_tr,      rhs_trdiv,                      &
         icepplus,    iceppminus,                     &
         icefluxes,   mass_matrix                   

    ! Variables needed for advection

    contains

    subroutine tg_rhs_icepack(mesh, trc)
    
        use mod_mesh
        use i_arrays
        use i_param
        use g_parsup
        use o_param
        use g_config
    
        implicit none
    
        ! Input - output
    
        type(t_mesh),        target,       intent(in)     :: mesh
        real(kind=dbl_kind), dimension(:), intent(inout)  :: trc
    
        ! Local variables
    
        real(kind=dbl_kind)              :: diff,       um,    vm,    vol, & 
                                            entries(3), dx(3), dy(3)
        integer(kind=int_kind)           :: n,          q,     row,        &
                                            elem,       elnodes(3)
    
    #include "associate_mesh.h"
    
        ! Taylor-Galerkin (Lax-Wendroff) rhs
      
        do row = 1, nx_nh
           rhs_tr(row)=c0
        enddo
    
        ! Velocities at nodes
    
        do elem = 1, nx_elem_nh  !assembling rhs over elements
    
           elnodes = elem2D_nodes(:,elem)
    
           ! Derivatives
           dx  = gradient_sca(1:3,elem)
           dy  = gradient_sca(4:6,elem)
           vol = elem_area(elem)
           um  = sum(uvel_elem(elnodes))
           vm  = sum(vvel_elem(elnodes))
    
           ! Diffusivity
    
           diff = ice_diff * sqrt( elem_area(elem) / scale_area )
           do n = 1, 3
              row = elnodes(n)
              do q = 1, 3
                 entries(q) = vol*ice_dt*((dx(n)*(um+uvel_elem(elnodes(q))) +      &
                              dy(n)*(vm+vvel_elem(elnodes(q))))/12.0_WP -          &
                              diff*(dx(n)*dx(q)+ dy(n)*dy(q)) -                &
                              0.5_WP*ice_dt*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0_WP)
              enddo
              rhs_tr(row)=rhs_tr(row)+sum(entries*trc(elnodes))
           enddo
        enddo
    
    end subroutine tg_rhs_icepack
    
    !=======================================================================
    
    subroutine fct_init_icepack(mesh)
    
        use o_param
        use o_mesh
        use i_arrays
        use g_parsup
    
        type(t_mesh), intent(in), target :: mesh
          
        implicit none
      
        ! Initialization of arrays necessary to implement FCT algorithm
        allocate(trl(nx))   ! low-order solutions
        allocate(d_tr(nx))  ! increments of high
                            ! order solutions
        allocate(icefluxes(nx_elem_nh, 3))
        allocate(icepplus(nx), icepminus(nx))
        allocate(rhs_tr(n_nx),  rhs_trdiv(nx))
        allocate(mass_matrix(sum(nn_num(1:nx_nh))))

      
        trl(:)    = c0
        d_tr(:)   = c0
        rhs_tr(:) = c0
        rhs_trdiv(:)   = c0
        icefluxes(:,:) = c0
        icepplus(:)    = c0
        icepminus(:)   = c0
        mass_matrix(:,:) = c0
      
        ! Fill in  the mass matrix
        call fill_mass_matrix_icepack(mesh)
      
        if (mype==0) write(*,*) 'Icepack FCT is initialized'
    
    end subroutine fct_init_icepack
    
    !=======================================================================
    
    subroutine fill_mass_matrix_icepack(mesh)
    
        use mod_mesh
        use o_mesh
        use i_param
        use i_arrays
        use g_parsup
      
        implicit none
      
        integer(kind=int_kind)                 :: n, n1, n2, row
        integer(kind=int_kind)                 :: elem, elnodes(3), q, offset, col, ipos
        integer(kind=int_kind), allocatable    :: col_pos(:)
        real(kind=dbl_kind)                    :: aa
        integer(kind=int_kind)                 :: flag=0 ,iflag=0
        type(t_mesh), intent(in), target       :: mesh
      
    #include "associate_mesh.h"
      
        allocate(col_pos(nx))
          
        do elem=1,nx_elem_nh
           elnodes=elem2D_nodes(:,elem)
           do n = 1, 3
              row = elnodes(n)
              if ( row > nx_nh ) cycle
              ! Global-to-local neighbourhood correspondence
              do q = 1, nn_num(row)
                 col_pos(nn_pos(q,row))=q
              enddo
              offset = ssh_stiff%rowptr(row) - ssh_stiff%rowptr(1)
              do q = 1, 3
                 col  = elnodes(q)
                 ipos = offset+col_pos(col)
                 mass_matrix(ipos) = mass_matrix(ipos) + elem_area(elem) / 12.0_WP
                 if ( q == n ) then
                 mass_matrix(ipos) = mass_matrix(ipos) + elem_area(elem) / 12.0_WP
                 end if
              enddo
           enddo
        enddo
      
        ! TEST: area == sum of row entries in mass_matrix:
        do q = 1, nx_nh
           offset = ssh_stiff%rowptr(q)   - ssh_stiff%rowptr(1) + 1
           n      = ssh_stiff%rowptr(q+1) - ssh_stiff%rowptr(1)
           aa     = sum(mass_matrix(offset:n))
           if ( abs(area(1,q)-aa) > p1) then
              iflag = q
              flag  = 1
           endif
        enddo
      
        if ( flag > 0 ) then
           offset = ssh_stiff%rowptr(iflag)   - ssh_stiff%rowptr(1)+1
           n      = ssh_stiff%rowptr(iflag+1) - ssh_stiff%rowptr(1)
           aa     = sum(mass_matrix(offset:n))
      
           write(*,*) '#### MASS MATRIX PROBLEM', mype, iflag, aa, area(1,iflag)
        endif
        deallocate(col_pos)
    
    end subroutine fill_mass_matrix_icepack
    
    !=======================================================================
    
    subroutine solve_low_order_icepack(mesh, trc)
    
        !============================
        ! Low-order solution
        !============================
        !
        ! It is assumed that m_ice, a_ice and m_snow from the previous time step
        ! are known at 1:myDim_nod2D+eDim_nod2D.
        ! We add diffusive contribution to the rhs. The diffusion operator
        ! is implemented as the difference between the consistent and lumped mass
        ! matrices acting on the field from the previous time step. The consistent
        ! mass matrix on the lhs is replaced with the lumped one.
    
        use mod_mesh
        use o_mesh
        use i_arrays
        use i_parm
        use g_parsup
        use g_comm_auto
      
        implicit none
      
        integer(kind=int_kind)           :: row, clo, clo2, cn, location(100)
        real   (kind=dbl_kind)           :: gamma
        type(t_mesh),        target,       intent(in)     :: mesh
        real(kind=dbl_kind), dimension(:), intent(inout)  :: trc
    
    #include "associate_mesh.h"
      
        gamma = ice_gamma_fct       ! Added diffusivity parameter
                                    ! Adjust it to ensure posivity of solution
      
        do row = 1, nx_nh
           clo  = ssh_stiff%rowptr(row)   - ssh_stiff%rowptr(1) + 1
           clo2 = ssh_stiff%rowptr(row+1) - ssh_stiff%rowptr(1)
           cn   = clo2 - clo + 1
           location(1:cn) = nn_pos(1:cn, row)
           trl(row) = (rhs_tr(row) + gamma * sum(mass_matrix(clo:clo2) * &
                      trc(location(1:cn)))) / area(1,row) +              &
                      (1.0_WP-gamma) * trc(row)
        enddo
      
        call exchange_nod(trl)
      
        ! Low-order solution must be known to neighbours
    
    end subroutine solve_low_order_icepack
    
    !=======================================================================
    
    subroutine solve_high_order_icepack(mesh, trc)
    
        use mod_mesh
        use o_mesh
        use i_arrays
        use i_parm
        use g_parsup
        use g_comm_auto
    
        implicit none
      
        integer(kind=int_kind)             :: n,i,clo,clo2,cn,location(100),row
        real   (kind=double_kind)          :: rhs_new
        integer(kind=int_kind), parameter  :: num_iter_solve = 3
        type(t_mesh),        target,       intent(in)     :: mesh
        real(kind=dbl_kind), dimension(:), intent(inout)  :: trc
      
    #include "associate_mesh.h"
      
        ! Taylor-Galerkin solution
       
        ! First guess solution
        do row = 1, nx_nh
           d_tr(row) = rhs_tr(row) / area(1,row)
        end do
      
        call exchange_nod(d_tr)
      
        ! Iterate
        do n = 1, num_iter_solve - 1
           do row = 1, nx_nh
              clo  = ssh_stiff%rowptr(row) - ssh_stiff%rowptr(1) + 1
              clo2 = ssh_stiff%rowptr(row+1) - ssh_stiff%rowptr(1)
              cn   = clo2 - clo + 1
              location(1:cn) = nn_pos(1:cn,row)
              rhs_new  = rhs_tr(row) - sum(mass_matrix(clo:clo2) * d_tr(location(1:cn)))
              trl(row) = d_tr(row) + rhs_new / area(1,row)
           enddo
           do row = 1, nx_nh
              d_tr(row) = trl(row)
           enddo
           call exchange_nod(d_tr)
        enddo
      
    end subroutine solve_high_order_icepack
    
    !=======================================================================
    
    subroutine fem_fct_icepack(mesh, trc)
    
        !============================
        ! Flux corrected transport algorithm for tracer advection
        !============================
        !
        ! It is based on Loehner et al. (Finite-element flux-corrected
        ! transport (FEM-FCT) for the Euler and Navier-Stokes equation,
        ! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
        ! Turek. (kuzmin@math.uni-dortmund.de)
    
        use mod_mesh
        use o_mesh
        use i_arrays
        use o_param
        use i_parm
        use g_parsup
        use g_comm_auto
    
        integer(kind=int_kind)                            :: icoef(3,3), n, q, elem, elnodes(3), row
        real   (kind=dbl_kind), allocatable, dimension(:) :: tmax, tmin
        real   (kind=dbl_kind)                            :: vol, flux, ae, gamma
        type(t_mesh),        target,       intent(in)     :: mesh  
        real(kind=dbl_kind), dimension(:), intent(inout)  :: trc
    
    #include "associate_mesh.h"
      
        gamma = ice_gamma_fct        ! It should coinside with gamma in
                                     ! ts_solve_low_order
    
        !==========================
        ! Compute elemental antidiffusive fluxes to nodes
        !==========================
        ! This is the most unpleasant part ---
        ! it takes memory and time. For every element
        ! we need its antidiffusive contribution to
        ! each of its 3 nodes
       
        allocate(tmax(nx_nh), tmin(nx_nh))
       
        ! Auxiliary elemental operator (mass matrix- lumped mass matrix)
       
        icoef = 1
       
        do n = 1, 3   ! three upper nodes
                      ! Cycle over rows  row=elnodes(n)
                 icoef(n,n) = -2
        enddo
       
        do elem = 1, nx_elem_nh
           elnodes = elem2D_nodes(:,elem)
           vol     = elem_area(elem)
           do q = 1, 3
              icefluxes(elem,q) = -sum(icoef(:,q) * (gamma * trc(elnodes) +       &
                                  d_tr(elnodes))) * (vol / area(1,elnodes(q))) / 12.0_WP
           enddo
        enddo
       
        !==========================
        ! Screening the low-order solution
        !==========================
        ! TO BE ADDED IF FOUND NECESSARY
        ! Screening means comparing low-order solutions with the
        ! solution on the previous time step and using whichever
        ! is greater/smaller in computations of max/min below
       
        !==========================
        ! Cluster min/max
        !==========================
       
        do row = 1, nx_nh
           n = nn_num(row)
           tmax(row) = maxval(trl(nn_pos(1:n,row)))
           tmin(row) = minval(trl(nn_pos(1:n,row)))
           ! Admissible increments
           tmax(row) = tmax(row) - trl(row)
           tmin(row) = tmin(row) - trl(row)
        enddo
       
        !=========================
        ! Sums of positive/negative fluxes to node row
        !=========================
       
        icepplus = c0
        icepminus = c0
        do elem = 1, nx_elem_nh
           elnodes = elem2D_nodes(:,elem)
           do q = 1, 3
              n    = elnodes(q)
              flux = icefluxes(elem,q)
              if ( flux > 0 ) then
                 icepplus(n) = icepplus(n) + flux
              else
                 icepminus(n) = icepminus(n) + flux
              endif
           enddo
        enddo
       
        !========================
        ! The least upper bound for the correction factors
        !========================
    
        do n = 1, nx_nh
           flux = icepplus(n)
           if ( abs(flux) > 0 ) then
              icepplus(n) = min(1.0,tmax(n) / flux)
           else
              icepplus(n) = c0
           endif
       
           flux = icepminus(n)
           if ( abs(flux) > 0 ) then
              icepminus(n) = min(1.0,tmin(n) / flux)
           else
              icepminus(n)=c0
           endif
         enddo
       
         ! pminus and pplus are to be known to neighbouting PE
         call exchange_nod(icepminus, icepplus)
       
         !========================
         ! Limiting
         !========================
       
         do elem = 1, nx_elem_nh
            elnodes = elem2D_nodes(:,elem)
            ae = c1
            do q = 1, 3
               n    = elnodes(q)
               flux = icefluxes(elem,q)
               if ( flux >=c0 ) ae = min(ae, icepplus(n))
               if ( flux < c0 ) ae = min(ae, icepminus(n))
            enddo
            icefluxes(elem,:) = ae * icefluxes(elem,:)
         enddo
       
       
         !==========================
         ! Update the solution
         !==========================
       
          do n = 1, nx_nh
             trc(n) = trl(n)
          end do
          do elem = 1, nx_elem_nh
             elnodes = elem2D_nodes(:,elem)
             do q = 1, 3
                n = elnodes(q)
                trc(n) = trc(n) + icefluxes(elem,q)
             enddo
          enddo
       
          call exchange_nod(trc)
       
          deallocate(tmin, tmax)
    
    end subroutine fem_fct_icepack
    
    !=======================================================================
    
    subroutine tg_rhs_div_icepack(mesh, trc)
    
        use mod_mesh
        use o_mesh
        use i_arrays
        use o_param
        use i_parm
        use g_parsup
        use g_comm_auto
    
        implicit none
    
        real   (kind=dbl_kind)    :: diff, entries(3), um, vm, vol, dx(3), dy(3)
        integer(kind=int_kind)    :: n, q, row, elem, elnodes(3)
        real   (kind=dbl_kind)    :: c_1, c_2, c_3, c_4, c_x, entries2(3)
        type(t_mesh),        target,       intent(in)     :: mesh
        real(kind=dbl_kind), dimension(:), intent(inout)  :: trc    

    
    #include "associate_mesh.h"
    
        ! Computes the rhs in a Taylor-Galerkin way (with urrayspwind 
        ! type of correction for the advection operator).
        ! In this version I tr to split divergent term off, 
        ! so that FCT works without it.
    
        do row = 1, nx_nh
           ! row=myList_nod2D(m)
           rhs_tr(row)    = c0
           rhs_trdiv(row) = c0
        enddo
      
        do elem = 1, nx_elem_nh   !! assembling rhs over elements
                                  !! elem=myList_elem2D(m)
           elnodes = elem2D_nodes(:,elem)
      
           ! Derivatives
           dx  = gradient_sca(1:3,elem)
           dy  = gradient_sca(4:6,elem)
           vol = elem_area(elem)
           um  = sum(uvel_elem(elnodes))
           vm  = sum(vvel_elem(elnodes))
      
           ! This is exact computation (no assumption of u=const 
           ! on elements used in the standard version)
           c_1 = (um*um+sum(uvel_elem(elnodes)*uvel_elem(elnodes))) / 12.0_dbl_kind
           c_2 = (vm*vm+sum(vvel_elem(elnodes)*vvel_elem(elnodes))) / 12.0_dbl_kind
           c_3 = (um*vm+sum(vvel_elem(elnodes)*uvel_elem(elnodes))) / 12.0_dbl_kind
           c_4 = sum(dx*uvel_elem(elnodes)+dy*vvel_elem(elnodes))
      
           do n = 1, 3
              row = elnodes(n)
      
              do q = 1, 3
                 entries(q)  = vol*ice_dt*((c1-p5*ice_dt*c_4)*(dx(n)*(um+uvel_elem(elnodes(q)))+ &
                               dy(n)*(vm+vvel_elem(elnodes(q))))/12.0_dbl_kind                 - &
                               p5*ice_dt*(c_1*dx(n)*dx(q)+c_2*dy(n)*dy(q)+c_3*(dx(n)*dy(q)+dx(q)*dy(n))))
                 entries2(q) = p5*ice_dt*(dx(n)*(um+uvel_elem(elnodes(q)))                     + &
                               dy(n)*(vm+vvel_elem(elnodes(q)))-dx(q)*(um+uvel_elem(row))      - &
                               dy(q)*(vm+vvel_elem(row)))
              enddo
              c_x = vol*ice_dt*c_4*(sum(trc(elnodes))+trc(elnodes(n))+sum(entries2*trc(elnodes))) / 12.0_dbl_kind
              rhs_tr(row)    = rhs_tr(row) + sum(entries * trc(elnodes)) + c_x
              rhs_trdiv(row) = rhs_trdiv(row) - c_x
           enddo
        enddo
    
    end subroutine tg_rhs_div_icepack
    
    !=======================================================================
    
    subroutine update_for_div_icepack(mesh, trc)
    
        use mod_mesh
        use o_mesh
        use i_arrays
        use o_param
        use i_parm
        use g_parsup
        use g_comm_auto
    
        implicit none
    
        integer(kind=int_kind)            :: n, i, clo, clo2, cn, &
                                             location(100), row
        real   (kind=dbl_kind)            :: rhs_new
        integer(kind=int_kind), parameter :: num_iter_solve=3
        type(t_mesh),        target,       intent(in)     :: mesh
        real(kind=dbl_kind), dimension(:), intent(inout)  :: trc
    
    #include "associate_mesh.h"
    
        ! Computes Taylor-Galerkin solution
        ! first approximation
      
        do row = 1, nx_nh
           d_tr(row) = rhs_trdiv(row) / area(1,row)
        enddo
      
        call exchange_nod(d_tr)
      
        ! Iterate
      
        do n = 1, num_iter_solve-1
           do row = 1, nx_nh
              clo            = ssh_stiff%rowptr(row)   - ssh_stiff%rowptr(1) + 1
              clo2           = ssh_stiff%rowptr(row+1) - ssh_stiff%rowptr(1)
              cn             = clo2 - clo + 1
              location(1:cn) = nn_pos(1:cn, row)
              rhs_new        = rhs_trdiv(row) - sum(mass_matrix(clo:clo2) * d_tr(location(1:cn)))
              trl(row)       = d_tr(row) + rhs_new / area(1,row)
           enddo
           do row = 1, nx_nh
              d_tr(row) = trl(row)
           enddo
           call exchange_nod(d_tr)
        enddo
      
        trc = trc + d_tr
    
    end subroutine update_for_div_icepack 
    
    !=======================================================================
    
    subroutine fct_solve_icepack(mesh, trc)
          
        implicit none

        real(kind=dbl_kind), dimension(:), intent(inout) :: trc      
        type(t_mesh),        target,       intent(in)    :: mesh
      
        ! Driving sequence
        call ice_TG_rhs_div(mesh, trc)
        call solve_high_order_icepack(mesh, trc) ! uses arrays of low-order solutions as temp
                                                 ! storage. It should preceed the call of low
                                                 ! order solution.
        call solve_low_order_icepack(mesh, trc)
        call fem_fct_icepack(mesh, trc)
        call ice_update_for_div(mesh, trc)

    end subroutine fct_solve_icepack

    !=======================================================================

    subroutine tracer_advection_icepack(mesh, trc)

        implicit none

        real(kind=dbl_kind), dimension(:), intent(inout) :: trc
        type(t_mesh),        target,       intent(in)    :: mesh


    end subroutine tracer_advection_icepack

    !=======================================================================

    subroutine work_to_state (nx,                        &
                              ntrcr,                     &
                              narr,     trcr_depend,     &
                              aicen,    trcrn,           &
                              vicen,    vsnon,           &
                              aice0,    works)

      use icedrv_main,        only: ncat, nslyr, nilyr, salinz
      use icepack_intfc
      use icedrv_system,      only: icedrv_system_abort
      use icedrv_flux,        only: salinz

      integer (kind=int_kind), intent(in) ::     &
         nx      , & ! block dimensions
         ntrcr   , & ! number of tracers in use
         narr        ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::     &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx,ncat), intent(out) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume per unit area of ice          (m)
         vsnon       ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx,ntrcr,ncat), intent(out) ::     &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx), intent(out) :: &
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
         works     ! work array

      ! local variables

      integer (kind=int_kind) :: &
         nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, ktherm

      logical (kind=log_kind) :: &
         tr_pond_cesm, tr_pond_lvl, tr_pond_topo, heat_capacity

      integer (kind=int_kind) ::      &
         k, i, n, it   , & ! counting indices
         narrays       , & ! counter for number of state variable arrays
         nt_qsno       , &
         nt_qice       , &
         nt_sice

      real (kind=dbl_kind) :: &
         rhos       , &
         rhoi       , &
         Lfresh     , &
         Tsmelt

      real (kind=dbl_kind), dimension(ncat) :: &
         tmp, exc, puny

      real (kind=dbl_kind), parameter :: &
         small = 0.000001_dbl_kind

      character(len=*), parameter :: subname = '(state_to_work)'

      call icepack_query_tracer_flags(tr_pond_cesm_out=tr_pond_cesm,              &
           tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
           nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno,                              &
           nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_Tsfc_out=nt_Tsfc) 
      call icepack_query_parameters(rhoi_out=rhoi,     rhos_out=rhoi,                   &
                                    Lfresh_out=Lfresh, heat_capacity_out=heat_capacity, &
                                    Tsmelt_out=Tsmelt, ktherm_out=ktherm,               &
                                    puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)

      ! Open water fraction

      trcrn(:,:,:) = c0
      aicen(:,:)   = c0
      vicen(:,:)   = c0
      vsnon(:,:)   = c0

      do i = 1, nx
         if (works(i,1) <= puny) then
             aice0(i) = c0
         else if (works(i,1) >= c1) then
             aice0(i) = c1
         else
             aice0(i) = works(i,1)
         end if
      enddo
      narrays = 1

      ! Sea ice area and volume per unit area of ice and snow

      do n=1,ncat
         do i = 1, nx
            if (works(i,narrays+1) > c1) then
                works(i,narrays+1) = c1
            end if
            if (works(i,narrays+1) <= small .or. works(i,narrays+2) <= small) then
                works(i,narrays+1) = c0
                works(i,narrays+2) = c0
                works(i,narrays+3) = c0
            end if
            if (works(i,narrays+3) <= small) then
               works(i,narrays+3) = c0
            end if
            aicen(i,n) = works(i,narrays+1)
            vicen(i,n) = works(i,narrays+2)
            vsnon(i,n) = works(i,narrays+3)
         end do

      narrays = narrays + 3 + ntrcr
      end do

      do i = 1, nx  ! For each grid cell
         if (sum(aicen(i,:)) > c1) then
            tmp(:) = c0
            exc(:) = c0
            do n = 1, ncat
               if (aicen(i,n) > puny) tmp(n) = c1
            end do
            do n = 1, ncat
                exc(n) = max(c0,(sum(aicen(i,:)) - c1))  &
                         * aicen(i,n) / sum(aicen(i,:))
            end do
            do n = 1, ncat
               aicen(i,n) = max(c0,aicen(i,n) - exc(n))
               aice0      = max(c0,sum(aicen(i,:)))
            end do
         end if
      end do

      narrays = 1

      do n=1, ncat

         narrays = narrays + 3

         do it = 1, ntrcr

            if (trcr_depend(it) == 0) then
               do i = 1, nx
                  if (aicen(i,n) > c0) then
                     if (it == nt_Tsfc) then
                         trcrn(i,it,n) = min(c0,works(i,narrays+it)/aicen(i,n))
                     else
                         trcrn(i,it,n) = works(i,narrays+it) / aicen(i,n)
                     end if
                  end if
               enddo
            elseif (trcr_depend(it) == 1) then
               do i = 1, nx
                  if (vicen(i,n) > c0) then
                      if (it >= nt_qice .and. it < nt_qice+nilyr) then
                          trcrn(i,it,n) = min(c0,works(i,narrays+it)/vicen(i,n))
                          if (.not. heat_capacity) trcrn(i,it,n) = -rhoi * Lfresh
                      else if (it >= nt_sice .and. it < nt_sice+nilyr) then
                          trcrn(i,it,n) = max(c0,works(i,narrays+it)/vicen(i,n))
                      end if
                  end if
               enddo
            elseif (trcr_depend(it) == 2) then
               do i = 1, nx
                   if (vsnon(i,n) > c0) then
                       if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                           trcrn(i,it,n) = min(c0,works(i,narrays+it)/vsnon(i,n)) - rhos*Lfresh
                           if (.not. heat_capacity) trcrn(i,it,n) = -rhos * Lfresh
                       end if
                   end if
               enddo
            ! Tracers not yet checked or implemented
            !elseif (trcr_depend(it) == 2+nt_alvl) then
            !   do i = 1, nx
            !      works(i,narrays+it) = aicen(i,n) &
            !                          * trcrn(i,nt_alvl,n) &
            !                          * trcrn(i,it,n)
            !   enddo
            !elseif (trcr_depend(it) == 2+nt_apnd .and. &
            !        tr_pond_cesm .or. tr_pond_topo) then
            !   do i = 1, nx
            !      works(i,narrays+it) = aicen(i,n) &
            !                          * trcrn(i,nt_apnd,n) &
            !                          * trcrn(i,it,n)
            !   enddo
            !elseif (trcr_depend(it) == 2+nt_apnd .and. &
            !        tr_pond_lvl) then
            !   do i = 1, nx
            !      works(i,narrays+it) = aicen(i,n) &
            !                          * trcrn(i,nt_alvl,n) &
            !                          * trcrn(i,nt_apnd,n) &
            !                          * trcrn(i,it,n)
            !   enddo
            !elseif (trcr_depend(it) == 2+nt_fbri) then
            !   do i = 1, nx
            !      works(i,narrays+it) = vicen(i,n) &
            !                          * trcrn(i,nt_fbri,n) &
            !                          * trcrn(i,it,n)
            !   enddo
            endif
         enddo

         narrays = narrays + ntrcr

      enddo                     ! number of categories 


end module icedrv_advection
