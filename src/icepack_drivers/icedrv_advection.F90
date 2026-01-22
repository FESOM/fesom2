!=======================================================================
!
! This submodule contains the subroutines 
! for the advection of sea ice tracers
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!
!=======================================================================

submodule (icedrv_main) icedrv_advection

    use icedrv_kinds
    use icedrv_constants
    use icedrv_system,      only: icedrv_system_abort
    use g_comm_auto,        only: exchange_nod
    use icepack_intfc,      only: icepack_warnings_flush,         &
                                  icepack_warnings_aborted,       &
                                  icepack_query_tracer_indices,   &
                                  icepack_query_tracer_flags,     &
                                  icepack_query_parameters,       &
                                  icepack_query_tracer_sizes

    implicit none

    real(kind=dbl_kind), allocatable, dimension(:)   :: &
         d_tr,        trl,                              &
         rhs_tr,      rhs_trdiv,                        &
         icepplus,    icepminus,                        &
         mass_matrix                   

    real(kind=dbl_kind), allocatable, dimension(:,:) :: &
         icefluxes                   

    ! Variables needed for advection

    contains

    subroutine tg_rhs_icepack(mesh, trc)
    
        use mod_mesh
        use i_param
        use g_parsup
        use o_param
        use g_config
    
        implicit none
    
        ! Input - output
    
        type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc
    
        ! Local variables
    
        real(kind=dbl_kind)              :: diff,       um,    vm,    vol, & 
                                            entries(3), dx(3), dy(3)
        integer(kind=int_kind)           :: n,          q,     row,        &
                                            elem,       elnodes(3)
    
#include "../associate_mesh.h"
    
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
           um  = sum(uvel(elnodes))
           vm  = sum(vvel(elnodes))
    
           ! Diffusivity
    
           diff = ice_diff * sqrt( elem_area(elem) / scale_area )
           do n = 1, 3
              row = elnodes(n)
              do q = 1, 3
                 entries(q) = vol*ice_dt*((dx(n)*(um+uvel(elnodes(q))) +      &
                              dy(n)*(vm+vvel(elnodes(q))))/12.0_WP -          &
                              diff*(dx(n)*dx(q)+ dy(n)*dy(q)) -               &
                              0.5_WP*ice_dt*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0_WP)
              enddo
              rhs_tr(row)=rhs_tr(row)+sum(entries*trc(elnodes))
           enddo
        enddo
    
    end subroutine tg_rhs_icepack
    
    !=======================================================================
    
    module subroutine init_advection_icepack(mesh)
    
        use o_param
        use o_mesh
        use g_parsup
        use mod_mesh
    
        implicit none
      
        type(t_mesh), intent(in), target :: mesh
          
        ! Initialization of arrays necessary to implement FCT algorithm
        allocate(trl(nx))   ! low-order solutions
        allocate(d_tr(nx))  ! increments of high
                            ! order solutions
        allocate(icefluxes(nx_elem_nh, 3))
        allocate(icepplus(nx), icepminus(nx))
        allocate(rhs_tr(nx),  rhs_trdiv(nx))
        allocate(mass_matrix(sum(nn_num(1:nx_nh))))

      
        trl(:)    = c0
        d_tr(:)   = c0
        rhs_tr(:) = c0
        rhs_trdiv(:)   = c0
        icefluxes(:,:) = c0
        icepplus(:)    = c0
        icepminus(:)   = c0
        mass_matrix(:) = c0
      
        ! Fill in  the mass matrix
        call fill_mass_matrix_icepack(mesh)
      
        if (mype==0) write(*,*) 'Icepack FCT is initialized'
    
    end subroutine init_advection_icepack
    
    !=======================================================================
    
    subroutine fill_mass_matrix_icepack(mesh)
    
        use mod_mesh
        use o_mesh
        use i_param
        use g_parsup
      
        implicit none
      
        integer(kind=int_kind)                 :: n, n1, n2, row
        integer(kind=int_kind)                 :: elem, elnodes(3), q, offset, col, ipos
        integer(kind=int_kind), allocatable    :: col_pos(:)
        real(kind=dbl_kind)                    :: aa
        integer(kind=int_kind)                 :: flag=0 ,iflag=0
        type(t_mesh), intent(in), target       :: mesh
      
#include "../associate_mesh.h"
      
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
        use i_param
        use g_parsup
  
      
        implicit none
      
        integer(kind=int_kind)           :: row, clo, clo2, cn, location(100)
        real   (kind=dbl_kind)           :: gamma
        type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc
    
#include "../associate_mesh.h"
      
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
        use i_param
        use g_parsup
   
    
        implicit none
      
        integer(kind=int_kind)             :: n,i,clo,clo2,cn,location(100),row
        real   (kind=dbl_kind)          :: rhs_new
        integer(kind=int_kind), parameter  :: num_iter_solve = 3
        type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc
      
#include "../associate_mesh.h"
      
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
        use o_param
        use i_param
        use g_parsup
   
    
        integer(kind=int_kind)                            :: icoef(3,3), n, q, elem, elnodes(3), row
        real   (kind=dbl_kind), allocatable, dimension(:) :: tmax, tmin
        real   (kind=dbl_kind)                            :: vol, flux, ae, gamma
        type(t_mesh),        target,        intent(in)     :: mesh  
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc
    
#include "../associate_mesh.h"
      
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
        use o_param
        use i_param
        use g_parsup
  
    
        implicit none
    
        real   (kind=dbl_kind)    :: diff, entries(3), um, vm, vol, dx(3), dy(3)
        integer(kind=int_kind)    :: n, q, row, elem, elnodes(3)
        real   (kind=dbl_kind)    :: c_1, c_2, c_3, c_4, c_x, entries2(3)
        type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc    

#include "../associate_mesh.h"
    
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
           um  = sum(uvel(elnodes))
           vm  = sum(vvel(elnodes))
      
           ! This is exact computation (no assumption of u=const 
           ! on elements used in the standard version)
           c_1 = (um*um+sum(uvel(elnodes)*uvel(elnodes))) / 12.0_dbl_kind
           c_2 = (vm*vm+sum(vvel(elnodes)*vvel(elnodes))) / 12.0_dbl_kind
           c_3 = (um*vm+sum(vvel(elnodes)*uvel(elnodes))) / 12.0_dbl_kind
           c_4 = sum(dx*uvel(elnodes)+dy*vvel(elnodes))
      
           do n = 1, 3
              row = elnodes(n)
      
              do q = 1, 3
                 entries(q)  = vol*ice_dt*((c1-p5*ice_dt*c_4)*(dx(n)*(um+uvel(elnodes(q)))+ &
                               dy(n)*(vm+vvel(elnodes(q))))/12.0_dbl_kind                 - &
                               p5*ice_dt*(c_1*dx(n)*dx(q)+c_2*dy(n)*dy(q)+c_3*(dx(n)*dy(q)+dx(q)*dy(n))))
                 entries2(q) = p5*ice_dt*(dx(n)*(um+uvel(elnodes(q)))                     + &
                               dy(n)*(vm+vvel(elnodes(q)))-dx(q)*(um+uvel(row))      - &
                               dy(q)*(vm+vvel(row)))
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
        use o_param
        use i_param
        use g_parsup
    
    
        implicit none
    
        integer(kind=int_kind)            :: n, i, clo, clo2, cn, &
                                             location(100), row
        real   (kind=dbl_kind)            :: rhs_new
        integer(kind=int_kind), parameter :: num_iter_solve=3
        type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc
    
#include "../associate_mesh.h"
    
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
        
        use mod_mesh     
     
        implicit none

        real(kind=dbl_kind), dimension(nx), intent(inout) :: trc      
        type(t_mesh),        target,        intent(in)    :: mesh
      
        ! Driving sequence
        call tg_rhs_div_icepack(mesh, trc)
        call solve_high_order_icepack(mesh, trc) ! uses arrays of low-order solutions as temp
                                                 ! storage. It should preceed the call of low
                                                 ! order solution.
        call solve_low_order_icepack(mesh, trc)
        call fem_fct_icepack(mesh, trc)
        call update_for_div_icepack(mesh, trc)

    end subroutine fct_solve_icepack

    !=======================================================================

    module subroutine tracer_advection_icepack(mesh)

        use mod_mesh
        use icepack_intfc,        only: icepack_aggregate
        use icepack_itd,          only: cleanup_itd
        use g_config,             only: dt

        implicit none
            
        integer (kind=int_kind) :: ntrcr, ntrace, narr, nbtrcr, i,     &
                                   nt,   nt1,    k
        integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno,                          & 
                                   nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
                                   nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S
        logical (kind=log_kind) :: tr_pond_topo, tr_pond_lvl, tr_pond_cesm,            &
                                   tr_pond,      tr_aero,     tr_FY,                   &
                                   tr_iage,      heat_capacity
        real    (kind=dbl_kind) :: puny
      
        ! Tracer dependencies and additional arrays
      
        integer (kind=int_kind), dimension(:),    allocatable    ::    &
                tracer_type    , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
                depend             ! tracer dependencies (see below)
      
        logical (kind=log_kind), dimension (:),   allocatable   ::     &
                has_dependents    ! true if a tracer has dependent tracers
      
        real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
               works
      
        type(t_mesh),        target,       intent(in)    :: mesh

        call icepack_query_parameters(heat_capacity_out=heat_capacity,    &
                                      puny_out=puny)
        call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
        call icepack_query_tracer_flags(                                  &
               tr_iage_out=tr_iage, tr_FY_out=tr_FY,                      &
               tr_aero_out=tr_aero, tr_pond_out=tr_pond,                  &
               tr_pond_cesm_out=tr_pond_cesm,                             &
               tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
               nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,       &
               nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_alvl_out=nt_alvl,           &
               nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd,       &
               nt_ipnd_out=nt_ipnd, nt_bgc_Nit_out=nt_bgc_Nit, nt_bgc_S_out=nt_bgc_S)
      
        narr   = 1 + ncat * (3 + ntrcr) ! max number of state variable arrays
      
        ! Allocate works array
      
        if (allocated(works)) deallocate(works)
        allocate ( works(nx,narr) )

        works(:,:) = c0
      
        call state_to_work (ntrcr, narr, works(:,:))
      
        ! Advect each tracer
      
        do nt = 1, narr    
              call fct_solve_icepack ( mesh, works(:,nt) )
        end do
      
        call work_to_state (ntrcr, narr, works(:,:))
          
        ! cut off icepack
      
        call cut_off_icepack
     
        do i=1,nx
           if (ncat < 0) then ! Do we really need this?
     
              call cleanup_itd  (dt,                     ntrcr,                &
                                 nilyr,                  nslyr,                &
                                 ncat,                   hin_max(:),           &
                                 aicen(i,:),             trcrn(i,1:ntrcr,:),   &
                                 vicen(i,:),             vsnon(i,:),           &
                                 aice0(i),               aice(i),              &
                                 n_aero,                                       &
                                 nbtrcr,                 nblyr,                &
                                 tr_aero,                                      &
                                 tr_pond_topo,                                 &
                                 heat_capacity,                                &
                                 first_ice(i,:),                               &
                                 trcr_depend(1:ntrcr),   trcr_base(1:ntrcr,:), &
                                 n_trcr_strata(1:ntrcr), nt_strata(1:ntrcr,:), &
                                 fpond(i),               fresh(i),             &
                                 fsalt(i),               fhocn(i),             &
                                 faero_ocn(i,:),         fzsal(i),             &
                                 flux_bio(i,1:nbtrcr))
      
              call icepack_aggregate (ncat,                    &
                                     aicen(i,:),               &
                                     trcrn(i,1:ntrcr,:),       &
                                     vicen(i,:),               &
                                     vsnon(i,:),               &
                                     aice (i),                 &
                                     trcr (i,1:ntrcr),         &
                                     vice (i),                 &
                                     vsno (i),                 &
                                     aice0(i),                 &
                                     ntrcr,                    &
                                     trcr_depend  (1:ntrcr),   &
                                     trcr_base    (1:ntrcr,:), &
                                     n_trcr_strata(1:ntrcr),   &
                                     nt_strata    (1:ntrcr,:))
           end if
        end do

        deallocate(works)

    end subroutine tracer_advection_icepack

    !=======================================================================

    subroutine work_to_state (ntrcr, narr, works)

        use icepack_intfc,    only: icepack_compute_tracers

        integer (kind=int_kind), intent(in) :: ntrcr, narr

        real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, ktherm
  
        logical (kind=log_kind) :: &
           tr_lvl, tr_pond_cesm, tr_pond_lvl, tr_pond_topo, heat_capacity
  
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
           tmp, exc

        real (kind=dbl_kind) :: puny
  
        real (kind=dbl_kind), parameter :: &
           small = 0.000001_dbl_kind
  
        character(len=*), parameter :: subname = '(state_to_work)'
  
        call icepack_query_tracer_flags(tr_pond_cesm_out=tr_pond_cesm,              &
             tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo,            &
             tr_lvl_out=tr_lvl)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno,                              &
             nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_Tsfc_out=nt_Tsfc) 
        call icepack_query_parameters(rhoi_out=rhoi,     rhos_out=rhos,                   &
                                      Lfresh_out=Lfresh, heat_capacity_out=heat_capacity, &
                                      Tsmelt_out=Tsmelt, ktherm_out=ktherm,               &
                                      puny_out=puny)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)
   
        trcrn(:,:,:) = c0
        aicen(:,:)   = c0
        vicen(:,:)   = c0
        vsnon(:,:)   = c0
        aice0(:)     = c0

        ! Open water fraction
  
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

           do i = 1, nx
              call icepack_compute_tracers(ntrcr=ntrcr,trcr_depend=trcr_depend(:),    &
                                           atrcrn = works(i,narrays+1:narrays+ntrcr), &
                                           aicen  = aicen(i,n),                       &
                                           vicen  = vicen(i,n),                       &
                                           vsnon  = vsnon(i,n),                       &
                                           trcr_base     = trcr_base(:,:),            &
                                           n_trcr_strata = n_trcr_strata(:),          &
                                           nt_strata     = nt_strata(:,:),            &
                                           trcrn  = trcrn(i,:,n)) 
           enddo
           
           narrays = narrays + ntrcr
        enddo
 
!        do n=1, ncat
!  
!           narrays = narrays + 3
!  
!           do it = 1, ntrcr
!  
!              if (trcr_depend(it) == 0) then
!                 do i = 1, nx
!                    if (aicen(i,n) > c0) then
!                       if (it == nt_Tsfc) then
!                           trcrn(i,it,n) = min(c0,works(i,narrays+it)/aicen(i,n))
!                       else if (it == nt_alvl .or. it == nt_apnd) then
!                           trcrn(i,it,n) = max(c0,min(c1,works(i,narrays+it) / aicen(i,n)))
!                       endif
!                    end if
!                 enddo
!              elseif (trcr_depend(it) == 1) then
!                 do i = 1, nx
!                    if (vicen(i,n) > c0) then
!                        if (it >= nt_qice .and. it < nt_qice+nilyr) then
!                            trcrn(i,it,n) = min(c0,works(i,narrays+it)/vicen(i,n))
!                            if (.not. heat_capacity) trcrn(i,it,n) = -rhoi * Lfresh
!                        else if (it >= nt_sice .and. it < nt_sice+nilyr) then
!                            trcrn(i,it,n) = max(c0,works(i,narrays+it)/vicen(i,n))
!                        else
!                            trcrn(i,it,n) = max(c0,works(i,narrays+it)/vicen(i,n))
!                        end if
!                    end if
!                 enddo
!              elseif (trcr_depend(it) == 2) then
!                 do i = 1, nx
!                     if (vsnon(i,n) > c0) then
!                         if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
!                             trcrn(i,it,n) = min(c0,works(i,narrays+it)/vsnon(i,n)) - rhos*Lfresh
!                             if (.not. heat_capacity) trcrn(i,it,n) = -rhos * Lfresh
!                         else
!                             trcrn(i,it,n) = min(c0,works(i,narrays+it)/vsnon(i,n)) - rhos*Lfresh
!                         end if
!                     end if
!                 enddo
!              elseif (trcr_depend(it) == 2+nt_alvl) then
!                 do i = 1, nx
!                    if (trcrn(i,nt_alvl,n) > small) then
!                       trcrn(i,it,n) = max( c0, works(i,narrays+it)  &
!                                           / trcrn(i,nt_alvl,n) )     
!                    endif
!                 enddo
!              elseif (trcr_depend(it) == 2+nt_apnd .and. &
!                      tr_pond_cesm .or. tr_pond_topo) then
!                 do i = 1, nx
!                    if (trcrn(i,nt_apnd,n) > small) then
!                       trcrn(i,it,n) = max( c0, works(i,narrays+it)  &
!                                           / trcrn(i,nt_apnd,n) )      
!                    endif
!                 enddo
!              elseif (trcr_depend(it) == 2+nt_apnd .and. &
!                      tr_pond_lvl) then
!                 do i = 1, nx
!                    if (trcrn(i,nt_apnd,n) > small .and. trcrn(i,nt_alvl,n) > small) then
!                       trcrn(i,it,n) = max( c0, works(i,narrays+it)  &
!                                           / trcrn(i,nt_apnd,n)      &
!                                           / trcrn(i,nt_alvl,n)      &
!                                           / aicen(i,n) )
!                    endif
!                 enddo
!              elseif (trcr_depend(it) == 2+nt_fbri) then
!                 do i = 1, nx
!                        works(i,narrays+it) = vicen(i,n) &
!                                            / trcrn(i,nt_fbri,n) &
!                                            / trcrn(i,it,n)
!                 enddo
!              endif
!           enddo
!  
!           narrays = narrays + ntrcr
!  
!        enddo                 ! number of categories 

!        if (mype == 0) write(*,*) 'Tracer salinity: ', nt_sice, ' - ', (nt_sice + nilyr - 1)
!        if (mype == 0) write(*,*) 'ktherm: ', ktherm

        if (ktherm == 1) then ! For bl99 themodynamics
                              ! always ridefine salinity
                              ! after advection
           do i = 1, nx       ! For each grid cell
              do k = 1, nilyr
                 trcrn(i,nt_sice+k-1,:) = salinz(i,k)
              end do          ! nilyr
           end do
        end if                ! ktherm==1

    end subroutine work_to_state

    !=======================================================================

    subroutine state_to_work (ntrcr, narr, works)

        integer (kind=int_kind), intent(in) :: ntrcr, narr
  
        real (kind=dbl_kind), dimension(nx,narr), intent (out) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc
  
        logical (kind=log_kind) :: &
           tr_pond_cesm, tr_pond_lvl, tr_pond_topo
  
        integer (kind=int_kind) ::      &
           i, n, it   , & ! counting indices
           narrays    , & ! counter for number of state variable arrays
           nt_qsno
  
        real (kind=dbl_kind) :: &
           rhos       , &
           Lfresh
  
        character(len=*), parameter :: subname = '(state_to_work)'
  
        call icepack_query_tracer_flags(tr_pond_cesm_out=tr_pond_cesm, &
             tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc)
        call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)
  
        do i = 1, nx
           works(i,1) = aice0(i)
        enddo
        narrays = 1
  
        do n=1, ncat
  
           do i = 1, nx
              works(i,narrays+1) = aicen(i,n)
              works(i,narrays+2) = vicen(i,n)
              works(i,narrays+3) = vsnon(i,n)
           enddo                  ! i
           narrays = narrays + 3
  
           do it = 1, ntrcr
              if (trcr_depend(it) == 0) then
                 do i = 1, nx
                    works(i,narrays+it) = aicen(i,n)*trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 1) then
                 do i = 1, nx
                    works(i,narrays+it) = vicen(i,n)*trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 2) then
                 do i = 1, nx
                    if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                        works(i,narrays+it) = vsnon(i,n)*trcrn(i,it,n) 
                    else
                        works(i,narrays+it) = vsnon(i,n)*trcrn(i,it,n)
                    end if
                 enddo
              elseif (trcr_depend(it) == 2+nt_alvl) then
                 do i = 1, nx
                    works(i,narrays+it) = aicen(i,n) &
                                        * trcrn(i,nt_alvl,n) &
                                        * trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_cesm .or. tr_pond_topo) then
                 do i = 1, nx
                    works(i,narrays+it) = aicen(i,n) &
                                        * trcrn(i,nt_apnd,n) &
                                        * trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_lvl) then
                 do i = 1, nx
                    works(i,narrays+it) = aicen(i,n) &
                                        * trcrn(i,nt_alvl,n) &
                                        * trcrn(i,nt_apnd,n) &
                                        * trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 2+nt_fbri) then
                 do i = 1, nx
                    works(i,narrays+it) = vicen(i,n) &
                                        * trcrn(i,nt_fbri,n) &
                                        * trcrn(i,it,n)
                 enddo
              endif
           enddo
           narrays = narrays + ntrcr
  
        enddo                     ! n
  
        if (narr /= narrays .and. mype == 0 ) write(ice_stderr,*)      &
            "Wrong number of arrays in transport bound call"

    end subroutine state_to_work

    !=======================================================================

    module subroutine cut_off_icepack

        use icepack_intfc,         only: icepack_compute_tracers
        use icepack_intfc,         only: icepack_aggregate
        use icepack_intfc,         only: icepack_init_trcr
        use icepack_intfc,         only: icepack_sea_freezing_temperature
        use icepack_therm_shared,  only: calculate_Tin_from_qin
        use icepack_mushy_physics, only: icepack_mushy_temperature_mush
  
        ! local variables
  
        real (kind=dbl_kind), dimension(nilyr) :: &
           qin            , & ! ice enthalpy (J/m3)
           qin_max        , & ! maximum ice enthalpy (J/m3)
           zTin               ! initial ice temperature
  
        real (kind=dbl_kind), dimension(nslyr) :: &
           qsn            , & ! snow enthalpy (J/m3)
           zTsn               ! initial snow temperature
        integer (kind=int_kind) ::      &
           i, n, k, it    , & ! counting indices
           narrays        , & ! counter for number of state variable arrays
           icells         , & ! number of ocean/ice cells
           ktherm         , &
           ntrcr
  
         real (kind=dbl_kind), dimension(ncat) :: &
           aicecat
  
        real (kind=dbl_kind) ::         &
           rhos,       Lfresh,          &
           cp_ice,     cp_ocn,          &
           qrd_snow,   qrd_ice,         &
           Tsfc,       exc,             &
           depressT,   Tmin,            &
           T_air_C,    hice,            &
           puny,       Tsmelt,          &
           small,      rhoi,            &
           hpnd_max          
           
           
  
        logical (kind=log_kind) :: tr_brine, tr_lvl, flag_snow, flag_cold_ice, flag_warm_ice, &
                                   tr_pond_cesm, tr_pond_topo, tr_pond_lvl, tr_FY, tr_iage,   &
                                   heat_capacity
        integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno, nt_sice
        integer (kind=int_kind) :: nt_fbri, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, nt_FY, nt_iage
  
        character(len=*), parameter :: subname = '(cut_off_icepack)'
  
        call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
        call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_lvl_out=tr_lvl,      &
          tr_pond_cesm_out=tr_pond_cesm, tr_pond_topo_out=tr_pond_topo,                &
          tr_pond_lvl_out=tr_pond_lvl, tr_FY_out=tr_FY, tr_iage_out=tr_iage            )
        call icepack_query_tracer_indices( nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice,   &
          nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_FY_out=nt_FY,                   &
          nt_fbri_out=nt_fbri, nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl,               &         
          nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd,               &
          nt_iage_out=nt_iage                                                          )
        call icepack_query_parameters(rhos_out=rhos, rhoi_out=rhoi, Lfresh_out=Lfresh, &
                                      cp_ice_out=cp_ice, cp_ocn_out=cp_ocn)
        call icepack_query_parameters(depressT_out=depressT, puny_out=puny, &
          Tsmelt_out=Tsmelt, ktherm_out=ktherm, heat_capacity_out=heat_capacity)
        call icepack_warnings_flush(ice_stderr)
  
        small = puny        
        Tmin  = -100.0_dbl_kind
  
        if (.not. heat_capacity) then ! for 0 layer thermodynamics
           do n = 1, ncat
              do i = 1, nx
                 if (trcrn(i,nt_Tsfc,n) > Tf(i) .or. trcrn(i,nt_Tsfc,n)< Tmin) then
                     trcrn(i,nt_Tsfc,n) = min(Tf(i), (T_air(i) + 273.15_dbl_kind))
                 endif
              enddo
           enddo
        endif
  
        if (heat_capacity) then    ! only for bl99 and mushy thermodynamics
  
        ! Here we should implement some conditions to check the tracers
        ! when ice is present, particularly enthalpy, surface temperature
        ! and salinity.
  
        do n = 1, ncat                      ! For each thickness cathegory
             do i = 1, nx                   ! For each grid point
  
                  call icepack_init_trcr(T_air(i),    Tf(i),         &
                                         salinz(i,:), Tmltz(i,:),    &
                                         Tsfc,                       &
                                         nilyr,        nslyr,        &
                                         qin   (:),    qsn  (:))
                  call icepack_warnings_flush(ice_stderr)
                  if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
                  file=__FILE__, line=__LINE__)  

                  ! Correct qin profile for melting temperatures

                  ! Maximum ice enthalpy 

                  qin_max(:) = rhoi * cp_ocn * (Tmltz(i,:) - 0.1_dbl_kind)  

                  if (vicen(i,n) > small .and. aicen(i,n) > small) then
  
                      ! Condition on surface temperature
                      if (trcrn(i,nt_Tsfc,n) > Tsmelt .or. trcrn(i,nt_Tsfc,n) < Tmin) then
                          trcrn(i,nt_Tsfc,n) = Tsfc
                      endif
  
                      ! Condition on ice enthalpy
  
                      flag_warm_ice = .false.
                      flag_cold_ice = .false.
                      flag_snow     = .false.
  
                      do k = 1, nilyr  ! Check for problems
  
                         if (ktherm == 2) then
                             zTin(k) = icepack_mushy_temperature_mush(trcrn(i,nt_qice+k-1,n),trcrn(i,nt_sice+k-1,n))
                         else
                             zTin(k) = calculate_Tin_from_qin(trcrn(i,nt_qice+k-1,n),Tmltz(i,k))
                         endif
  
                         if (zTin(k) <  Tmin      ) flag_cold_ice = .true.
                         if (zTin(k) >= Tmltz(i,k)) flag_warm_ice = .true.
  
                      enddo !nilyr
  
                      if (flag_cold_ice) then
  
                          trcrn(i,nt_Tsfc,n) = Tsfc
  
                          do k = 1, nilyr
                               trcrn(i,nt_qice+k-1,n) = min(qin_max(k), qin(k))
                          enddo        ! nilyr
  
                          if (vsnon(i,n) > small) then ! Only if there is snow
                                                       ! on top of the sea ice
                              do k = 1, nslyr
                                    trcrn(i,nt_qsno+k-1,n) = qsn(k)
                              enddo   
                          else                        ! No snow 
                              trcrn(i,nt_qsno:nt_qsno+nslyr-1,n) = c0
                          endif
  
                      end if            ! flag cold ice
  
                      if (flag_warm_ice) then         ! This sea ice should have melted already 
  
                          aicen(i,n)   = c0
                          vicen(i,n)   = c0
                          vsnon(i,n)   = c0
                          trcrn(i,:,n) = c0
                          trcrn(i,nt_Tsfc,n) = Tf(i)
  
                      end if
  
                      if (vsnon(i,n) > small) then
  
                          flag_snow = .false.
  
                          do k = 1, nslyr  
                             if (trcrn(i,nt_qsno+k-1,n) >= -rhos*Lfresh) flag_snow = .true.
                             zTsn(k) = (Lfresh + trcrn(i,nt_qsno+k-1,n)/rhos)/cp_ice
                             if (zTsn(k) < Tmin) flag_snow = .true.
                          end do 
  
                          if (flag_snow) then
                              trcrn(i,nt_Tsfc,n) = Tsfc
                              do k = 1, nslyr
                                   trcrn(i,nt_qsno+k-1,n) = qsn(k)
                              enddo        ! nslyr
                          endif            ! flag snow
                      endif ! vsnon(i,n) > c0
  
                  else
  
                      aicen(i,n)   = c0
                      vicen(i,n)   = c0
                      vsnon(i,n)   = c0
                      trcrn(i,:,n) = c0
                      trcrn(i,nt_Tsfc,n) = Tf(i)
  
                  endif
  
             enddo   ! nx
        enddo        ! ncat

        ! Melt ponds and level ice cutoff

        do n = 1, ncat
           do i = 1, nx
              if (aicen(i,n) > c0) then
                 hpnd_max = 0.9_dbl_kind * vicen(i,n) / aicen(i,n)
                 ! Sea ice age
                 if (tr_iage) then
                     if (trcrn(i,nt_iage,n) < 0.000001_dbl_kind) trcrn(i,nt_iage,n) = c0
                 end if
                 ! First year ice fraction
                 if (tr_FY) then
                     if (trcrn(i,nt_FY,n) < 0.000001_dbl_kind) trcrn(i,nt_FY,n) = c0
                     if (trcrn(i,nt_FY,n) > c1) trcrn(i,nt_FY,n) = c1
                 end if
                 ! Level ice
                 if (tr_lvl) then
                     if (trcrn(i,nt_alvl,n) > aicen(i,n)) then
                         trcrn(i,nt_alvl,n) = aicen(i,n)
                     elseif (trcrn(i,nt_alvl,n) < 0.000001_dbl_kind) then
                         trcrn(i,nt_alvl,n) = c0
                     endif
                     if (trcrn(i,nt_vlvl,n) < 0.000001_dbl_kind .or. trcrn(i,nt_alvl,n) < 0.000001_dbl_kind) trcrn(i,nt_vlvl,n) = c0
                     if (trcrn(i,nt_vlvl,n) > vicen(i,n)) trcrn(i,nt_vlvl,n) = vicen(i,n)
                 end if
                 ! CESM melt pond parameterization
                 if (tr_pond_cesm) then
                     if (trcrn(i,nt_apnd,n) > c1) then
                         trcrn(i,nt_apnd,n) = c1
                     elseif (trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) then
                         trcrn(i,nt_apnd,n) = c0
                     endif
                     if (trcrn(i,nt_hpnd,n) < 0.000001_dbl_kind .or. trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) trcrn(i,nt_hpnd,n) = c0
                     if (trcrn(i,nt_hpnd,n) > hpnd_max) trcrn(i,nt_hpnd,n) = hpnd_max
                 end if
                 ! Topo and level melt pond parameterization
                 if (tr_pond_topo .or. tr_pond_lvl) then
                     if (trcrn(i,nt_apnd,n) > c1) then
                         trcrn(i,nt_apnd,n) = c1
                     elseif (trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) then
                         trcrn(i,nt_apnd,n) = c0
                     endif
                     if (trcrn(i,nt_hpnd,n) < 0.000001_dbl_kind .or. trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) trcrn(i,nt_hpnd,n) = c0
                     if (trcrn(i,nt_hpnd,n) > hpnd_max) trcrn(i,nt_hpnd,n) = hpnd_max
                     if (trcrn(i,nt_ipnd,n) < 0.000001_dbl_kind .or. trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) trcrn(i,nt_ipnd,n) = c0
                 end if
                 ! Dynamic salt
                 if (tr_brine) then
                     if (trcrn(i,nt_fbri,n) < 0.000001_dbl_kind) trcrn(i,nt_fbri,n) = c0
                 endif
              endif
           enddo
        enddo

        do i = 1, nx
           aice(i) = c0
           vice(i) = c0
           vsno(i) = c0
           do it = 1, ntrcr
              trcr(i,it) = c0
           enddo
           call icepack_aggregate (ncat,                    &
                                  aicen(i,:),               &
                                  trcrn(i,1:ntrcr,:),       &
                                  vicen(i,:),               &
                                  vsnon(i,:),               &
                                  aice (i),                 &
                                  trcr (i,1:ntrcr),         &
                                  vice (i),                 &
                                  vsno (i),                 &
                                  aice0(i),                 &
                                  ntrcr,                    &
                                  trcr_depend  (1:ntrcr),   &
                                  trcr_base    (1:ntrcr,:), &
                                  n_trcr_strata(1:ntrcr),   &
                                  nt_strata    (1:ntrcr,:))
        end do
  
        end if ! heat_capacity
  
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
            file=__FILE__, line=__LINE__)

    end subroutine cut_off_icepack

end submodule icedrv_advection


