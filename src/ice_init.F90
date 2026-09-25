!> @brief
!> ice_init.F90
!! Initialisation of the sea-ice derived type.
!!
!! Separate from MOD_ICE.F90, which defines the types themselves: these routines
!! need MOD_MESH, MOD_PARTIT and ice_meltponds, and MOD_ICE is used by most of the
!! model, so those dependencies stay out of it.

module ice_init_module
    USE MOD_ICE
    USE MOD_PARTIT
    use par_support_module, only: par_ex
    USE MOD_MESH
    USE o_param, only: WP
    USE ice_meltponds, only: init_meltponds

    implicit none

    private
    public :: ice_init, ice_init_toyocean_dummy

contains

!
!
!_______________________________________________________________________________
! interface to initialise derived type for sea ice
!
!
!_______________________________________________________________________________
! initialise derived type for sea ice
subroutine ice_init(ice, partit, mesh)
    IMPLICIT NONE
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer        :: elem_size, node_size, n, ed(2)
    integer, save  :: nm_unit  = 105       ! unit to open namelist file, skip 100-102 for cray
    integer        :: iost
    !___________________________________________________________________________
    ! define ice namelist parameter
    integer        :: whichEVP, evp_rheol_steps, ice_ave_steps
    real(kind=WP)  :: Pstar, ellipse, c_pressure, delta_min, ice_gamma_fct, &
                      ice_diff, theta_io, alpha_evp, beta_evp, c_aevp, Cd_oce_ice
    namelist /ice_dyn/ whichEVP, Pstar, ellipse, c_pressure, delta_min, evp_rheol_steps, &
                       Cd_oce_ice, ice_gamma_fct, ice_diff, theta_io, ice_ave_steps, &
                       alpha_evp, beta_evp, c_aevp
    logical        :: snowdist, new_iclasses, use_meltponds, snowmelt_tgate
    integer        :: open_water_albedo, iclasses
    real(kind=WP)  :: Sice, h0, h0_s, emiss_ice, emiss_wat, albsn, albsnm, albi, &
                      albim, albw, con, consn, hmin, armin, c_melt, h_cutoff, h_ml, h_snowscale, alb_tramp
    namelist /ice_therm/ Sice, iclasses, h0, h0_s, hmin, armin,  emiss_ice, emiss_wat, albsn, albsnm, albi, &
                         albim, albw, con, consn,  snowdist, new_iclasses, open_water_albedo, c_melt, h_cutoff, h_ml, use_meltponds, &
                         h_snowscale, snowmelt_tgate, alb_tramp

    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    
    !___________________________________________________________________________
    ! pre initialise namelist parameters with the defaults, for the case they are not 
    ! mentioned in the namelist.ice, because otherwise they would be overwritten 
    ! with garbage
    whichEVP         = ice%whichEVP
    Pstar            = ice%pstar
    ellipse          = ice%ellipse
    c_pressure       = ice%c_pressure
    delta_min        = ice%delta_min
    evp_rheol_steps  = ice%evp_rheol_steps
    Cd_oce_ice       = ice%cd_oce_ice
    ice_gamma_fct    = ice%ice_gamma_fct
    ice_diff         = ice%ice_diff
    theta_io         = ice%theta_io
    ice_ave_steps    = ice%ice_ave_steps
    alpha_evp        = ice%alpha_evp
    beta_evp         = ice%beta_evp
    c_aevp           = ice%c_aevp
    con              = ice%thermo%con
    consn            = ice%thermo%consn
    Sice             = ice%thermo%Sice
    iclasses         = ice%thermo%iclasses
    h0               = ice%thermo%h0
    h0_s             = ice%thermo%h0_s
    hmin             = ice%thermo%hmin
    armin            = ice%thermo%armin
    emiss_ice        = ice%thermo%emiss_ice
    emiss_wat        = ice%thermo%emiss_wat
    albsn            = ice%thermo%albsn
    albsnm           = ice%thermo%albsnm
    albi             = ice%thermo%albi
    albim            = ice%thermo%albim
    albw             = ice%thermo%albw
    h_ml             = ice%thermo%h_ml
    h_snowscale      = ice%thermo%h_snowscale
    snowmelt_tgate   = ice%thermo%snowmelt_tgate
    alb_tramp        = ice%thermo%alb_tramp
    snowdist         = ice%thermo%snowdist
    new_iclasses     = ice%thermo%new_iclasses
    open_water_albedo= ice%thermo%open_water_albedo
    c_melt           = ice%thermo%c_melt
    h_cutoff         = ice%thermo%h_cutoff
    ! cc and cl are computed internally, not read from namelist:
    ! ice%thermo%cc = ice%thermo%rhowat*4190.0
    ! ice%thermo%cl = ice%thermo%rhoice*3.34e5

    !___________________________________________________________________________
    ! now parameters from derived types are overwriiten by the parameters
    ! in namelist.ice
    ! open and read namelist.ice for I/O
        open(unit=nm_unit, file='namelist.ice', form='formatted', access='sequential', status='old', iostat=iost )
    if (iost == 0) then
        if (mype==0) write(*,*) '     file   : ', 'namelist.ice',' open ok'
    else
        if (mype==0) write(*,*) 'ERROR: --> bad opening file   : ', 'namelist.ice',' ; iostat=',iost
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    end if
    read(nm_unit, nml=ice_dyn  , iostat=iost)
    read(nm_unit, nml=ice_therm, iostat=iost)
    close(nm_unit)

    !___________________________________________________________________________
    ! set parameters in ice derived type from namelist.ice --> namelist /ice_dyn/
    ! now they get written back into the derived types
    ice%whichEVP        = whichEVP
    ice%pstar           = Pstar
    ice%ellipse         = ellipse
    ice%c_pressure      = c_pressure
    ice%delta_min       = delta_min
    ice%evp_rheol_steps = evp_rheol_steps
    ice%cd_oce_ice      = Cd_oce_ice
    ice%ice_gamma_fct   = ice_gamma_fct
    ice%ice_diff        = ice_diff
    ice%theta_io        = theta_io
    ice%ice_ave_steps   = ice_ave_steps
    ice%alpha_evp       = alpha_evp
    ice%beta_evp        = beta_evp
    ice%c_aevp          = c_aevp

    ! set parameters in ice derived type from namelist.ice --> namelist /ice_therm/
    ice%thermo%con      = con
    ice%thermo%consn    = consn
    ice%thermo%Sice     = Sice
    ice%thermo%iclasses = iclasses
    ice%thermo%h0       = h0
    ice%thermo%h0_s     = h0_s
    ice%thermo%hmin     = hmin
    ice%thermo%armin    = armin
    ice%thermo%emiss_ice= emiss_ice
    ice%thermo%emiss_wat= emiss_wat
    ice%thermo%albsn    = albsn
    ice%thermo%albsnm   = albsnm
    ice%thermo%albi     = albi
    ice%thermo%albim    = albim
    ice%thermo%albw     = albw
    ice%thermo%h_ml     = h_ml
    ice%thermo%h_snowscale = h_snowscale
    ice%thermo%snowmelt_tgate = snowmelt_tgate
    ice%thermo%alb_tramp = alb_tramp
    ice%thermo%snowdist = snowdist
    ice%thermo%new_iclasses=new_iclasses
    ice%thermo%open_water_albedo=open_water_albedo
    ice%thermo%use_meltponds = use_meltponds
    if (use_meltponds) call init_meltponds('namelist.ice', partit%mype)
    ice%thermo%c_melt   = c_melt
    ice%thermo%h_cutoff = h_cutoff    
    ice%thermo%cc       =ice%thermo%rhowat*4190.0  ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
    ice%thermo%cl       =ice%thermo%rhoice*3.34e5  ! Volumetr. latent heat of ice fusion [J/m**3](cl=rhoice*Lf)

    !___________________________________________________________________________
    call set_ice_tracer_layout(ice)

    !___________________________________________________________________________
    ! define local vertice & elem array size
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D +eDim_nod2D

    !___________________________________________________________________________
    ! allocate/initialise arrays in ice derived type
    ! initialise velocity and stress related arrays in ice derived type
    allocate(ice%uice(                 node_size))
    allocate(ice%uice_rhs(             node_size))
    allocate(ice%uice_old(             node_size))
    allocate(ice%vice(                 node_size))
    allocate(ice%vice_rhs(             node_size))
    allocate(ice%vice_old(             node_size))
    allocate(ice%stress_atmice_x(      node_size))
    allocate(ice%stress_iceoce_x(      node_size))
    allocate(ice%stress_atmice_y(      node_size))
    allocate(ice%stress_iceoce_y(      node_size))
    allocate(ice%h_ice          (      node_size))
    allocate(ice%h_snow         (      node_size))    
    ice%uice            = 0.0_WP
    ice%uice_rhs        = 0.0_WP
    ice%uice_old        = 0.0_WP
    ice%stress_atmice_x = 0.0_WP
    ice%stress_iceoce_x = 0.0_WP
    ice%vice            = 0.0_WP
    ice%vice_rhs        = 0.0_WP
    ice%vice_old        = 0.0_WP
    ice%stress_atmice_y = 0.0_WP
    ice%stress_iceoce_y = 0.0_WP
    ice%h_ice           = 0.0_WP
    ice%h_snow          = 0.0_WP    
    if (ice%whichEVP /= 0) then
        allocate(ice%uice_aux(         node_size))
        allocate(ice%vice_aux(         node_size))
        ice%uice_aux    = 0.0_WP
        ice%vice_aux    = 0.0_WP
    end if
    if (ice%whichEVP == 2) then
        allocate(ice%alpha_evp_array(  elem_size))   ! element-indexed in stress_tensor_a/find_alpha_field_a
        allocate(ice%beta_evp_array(   node_size))   ! node-indexed in ice_momentum_step_a
        ice%alpha_evp_array = ice%alpha_evp
        ice%beta_evp_array  = ice%alpha_evp
    end if

    !___________________________________________________________________________
    ! initialise surface ocean arrays in ice derived type
    allocate(ice%srfoce_u(             node_size))
    allocate(ice%srfoce_v(             node_size))
    allocate(ice%srfoce_temp(          node_size))
    allocate(ice%srfoce_salt(          node_size))
    allocate(ice%srfoce_ssh(           node_size))
    ice%srfoce_u         = 0.0_WP
    ice%srfoce_v         = 0.0_WP
    ice%srfoce_temp      = 0.0_WP
    ice%srfoce_salt      = 0.0_WP
    ice%srfoce_ssh       = 0.0_WP

    allocate(ice%flx_fw(node_size))
    allocate(ice%flx_h( node_size))
    ice%flx_fw           = 0.0_WP
    ice%flx_h            = 0.0_WP
    
    !___________________________________________________________________________
    ! initialse data array of ice derived type containing "ice tracer" that have
    ! to be advected: a_ice (index=1), m_ice (index=2), m_snow (index=3),
    ! ice_temp (index=4, only when coupled)
    allocate(ice%data(ice%num_itracers))
    do n = 1, ice%num_itracers
        allocate(ice%data(n)%values(    node_size))
        allocate(ice%data(n)%values_old(node_size))
        allocate(ice%data(n)%values_rhs(node_size))
        allocate(ice%data(n)%values_div_rhs(node_size))
        allocate(ice%data(n)%dvalues(   node_size))
        allocate(ice%data(n)%valuesl(   node_size))
        ice%data(n)%ID             = n
        ice%data(n)%values         = 0.0_WP
        ice%data(n)%values_old     = 0.0_WP
        ice%data(n)%values_rhs     = 0.0_WP
        ice%data(n)%values_div_rhs = 0.0_WP
        ice%data(n)%dvalues        = 0.0_WP
        ice%data(n)%valuesl        = 0.0_WP
        if (n==4) ice%data(n)%values = 265.15_WP
    end do

    !___________________________________________________________________________
    ! initialse work array of ice derived type
    allocate(ice%work%fct_tmax(        node_size))
    allocate(ice%work%fct_tmin(        node_size))
    allocate(ice%work%fct_plus(        node_size))
    allocate(ice%work%fct_minus(       node_size))
    allocate(ice%work%fct_fluxes(      elem_size, 3))
    ice%work%fct_tmax    = 0.0_WP
    ice%work%fct_tmin    = 0.0_WP
    ice%work%fct_plus    = 0.0_WP
    ice%work%fct_minus   = 0.0_WP
    ice%work%fct_fluxes  = 0.0_WP

    allocate(ice%work%fct_massmatrix(sum(nn_num(1:myDim_nod2D))))
    ice%work%fct_massmatrix = 0.0_WP

    allocate(ice%work%sigma11(         elem_size))
    allocate(ice%work%sigma12(         elem_size))
    allocate(ice%work%sigma22(         elem_size))
    allocate(ice%work%eps11(           elem_size))
    allocate(ice%work%eps12(           elem_size))
    allocate(ice%work%eps22(           elem_size))
    ice%work%sigma11     = 0.0_WP
    ice%work%sigma12     = 0.0_WP
    ice%work%sigma22     = 0.0_WP
    ice%work%eps11       = 0.0_WP
    ice%work%eps12       = 0.0_WP
    ice%work%eps22       = 0.0_WP

    allocate(ice%work%inv_areamass(node_size))
    allocate(ice%work%inv_mass(    node_size))
    allocate(ice%work%ice_strength(node_size))
    ice%work%inv_areamass = 0.0_WP
    ice%work%inv_mass     = 0.0_WP
    ice%work%ice_strength = 0.0_WP

    !___________________________________________________________________________
    ! initialse thermo array of ice derived type
    allocate(ice%thermo%ustar(         node_size))
    allocate(ice%thermo%t_skin(        node_size))
    allocate(ice%thermo%thdgr(         node_size))
    allocate(ice%thermo%thdgrsn(       node_size))
    allocate(ice%thermo%thdgra(        node_size))
    allocate(ice%thermo%thdgr_old(     node_size))
    ! melt pond arrays
    allocate(ice%thermo%apnd(          node_size))
    allocate(ice%thermo%hpnd(          node_size))
    allocate(ice%thermo%ipnd(          node_size))
    ice%thermo%ustar     = 0.0_WP
    ice%thermo%t_skin    = 0.0_WP
    ice%thermo%thdgr     = 0.0_WP
    ice%thermo%thdgrsn   = 0.0_WP
    ice%thermo%thdgra    = 0.0_WP
    ice%thermo%thdgr_old = 0.0_WP
    ! initialize melt pond arrays
    ice%thermo%apnd      = 0.0_WP
    ice%thermo%hpnd      = 0.0_WP
    ice%thermo%ipnd      = 0.0_WP

    !___________________________________________________________________________
    ! initialse coupling array of ice derived type 
#if defined (__cpl_enabled)
    allocate(ice%atmcoupl%oce_flx_h(     node_size))
    allocate(ice%atmcoupl%ice_flx_h(     node_size))
    allocate(ice%atmcoupl%tmpoce_flx_h(  node_size))
    allocate(ice%atmcoupl%tmpice_flx_h(  node_size))
    ice%atmcoupl%oce_flx_h     = 0.0_WP
    ice%atmcoupl%ice_flx_h     = 0.0_WP
    ice%atmcoupl%tmpoce_flx_h  = 0.0_WP
    ice%atmcoupl%tmpice_flx_h  = 0.0_WP
    allocate(ice%atmcoupl%ice_alb(       node_size))
    allocate(ice%atmcoupl%enthalpyoffuse(node_size))
    allocate(ice%atmcoupl%runoff_liquid(node_size))
    allocate(ice%atmcoupl%runoff_solid(node_size))
    ice%atmcoupl%ice_alb       = 0.6_WP
    ice%atmcoupl%enthalpyoffuse= 0.0_WP
    ice%atmcoupl%runoff_liquid= 0.0_WP
    ice%atmcoupl%runoff_solid= 0.0_WP
    allocate(ice%atmcoupl%flx_qres(node_size))
    allocate(ice%atmcoupl%flx_qcon(node_size))
    ice%atmcoupl%flx_qres      = 0.0_WP
    ice%atmcoupl%flx_qcon      = 0.0_WP
    ! 0 = "no send yet": ice_surftemp falls back to the local t as anchor
    ! until the first actual OASIS transmission populates it.
    allocate(ice%atmcoupl%ist_ref(node_size))
    ice%atmcoupl%ist_ref       = 0.0_WP
#endif /* (__cpl_enabled) */

    !___________________________________________________________________________
    ! --> took from oce_mesh.F90 --> subroutine mesh_auxiliary_arrays(partit, mesh)
    ! to here since namelist.ice is now read in ice_init where whichEVP is not available
    ! when  mesh_auxiliary_arrays is called
    !array of 2D boundary conditions is used in ice_maEVP
    
    ! LA 2023-05-24 initiate bc_index_nod2D also for whichEVP==0
        allocate(mesh%bc_index_nod2D(myDim_nod2D+eDim_nod2D))
        mesh%bc_index_nod2D=1._WP
        do n=1, myDim_edge2D
            ed=mesh%edges(:, n)
            if (myList_edge2D(n) <= mesh%edge2D_in) cycle
            mesh%bc_index_nod2D(ed)=0._WP
        end do
        
    !___________________________________________________________________________
    ! initialise Lettis sea-ice cmip6 paramters! Idealy i would add them to the 
    ! ice derived type. but when i do this now i will mess up previous created raw/bin
    ! restart files. We should add this at a latter point 
    allocate(ice%thermo%dyngr(node_size), ice%thermo%dyngrsn(node_size), ice%thermo%dyngra(node_size))
    ice%thermo%dyngr   = 0.0_WP
    ice%thermo%dyngrsn = 0.0_WP
    ice%thermo%dyngra  = 0.0_WP
    
end subroutine ice_init  
!
!
!
!
!
!_______________________________________________________________________________
! initialise derived type for sea ice
subroutine ice_init_toyocean_dummy(ice, partit, mesh)
    IMPLICIT NONE
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer        :: node_size, n
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    call set_ice_tracer_layout(ice)

    !___________________________________________________________________________
    ! define local vertice & elem array size
    node_size=myDim_nod2D+eDim_nod2D

    !___________________________________________________________________________
    ! allocate/initialise arrays in ice derived type
    ! initialise velocity and stress related arrays in ice derived type
    allocate(ice%uice(                  node_size))
    allocate(ice%vice(                  node_size))
    ice%uice               = 0.0_WP
    ice%vice               = 0.0_WP
    allocate(ice%data(ice%num_itracers))
    do n = 1, ice%num_itracers
        allocate(ice%data(n)%values(    node_size))
        allocate(ice%data(n)%values_old(node_size))
        ice%data(n)%ID          = n
        ice%data(n)%values      = 0.0_WP
        ice%data(n)%values_old  = 0.0_WP
    end do

    allocate(ice%srfoce_temp(           node_size))
    allocate(ice%srfoce_salt(           node_size))
    allocate(ice%srfoce_ssh(            node_size))
    allocate(ice%srfoce_u(              node_size))
    allocate(ice%srfoce_v(              node_size))
    allocate(ice%stress_iceoce_x(       node_size))
    allocate(ice%stress_iceoce_y(       node_size))
    allocate(ice%stress_atmice_x(       node_size))
    allocate(ice%stress_atmice_y(       node_size))
    allocate(ice%flx_h(                 node_size))
    allocate(ice%flx_fw(                node_size))
    allocate(ice%thermo%thdgr(          node_size))
    allocate(ice%thermo%thdgrsn(        node_size))
    ice%srfoce_temp        = 0.0_WP
    ice%srfoce_salt        = 0.0_WP
    ice%srfoce_ssh         = 0.0_WP
    ice%srfoce_u           = 0.0_WP
    ice%srfoce_v           = 0.0_WP
    ice%stress_iceoce_x    = 0.0_WP
    ice%stress_iceoce_y    = 0.0_WP
    ice%stress_atmice_x    = 0.0_WP
    ice%stress_atmice_y    = 0.0_WP
    ice%flx_h              = 0.0_WP
    ice%flx_fw             = 0.0_WP
    ice%thermo%thdgr       = 0.0_WP
    ice%thermo%thdgrsn     = 0.0_WP
end subroutine ice_init_toyocean_dummy

end module ice_init_module
