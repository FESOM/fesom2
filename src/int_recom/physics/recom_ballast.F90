!===============================================================================
! MODULE: recom_ballast  (interface modules) + ballast subroutines
!
! Purpose:
!   Ballast-driven particle sinking parameterisation for REcoM.
!   Extracted from recom_sinking.F90 (Step 4 restructuring).
!
!   These three routines compute the density- and viscosity-based sinking
!   velocity scaling factors that are passed back into the FESOM tracer
!   transport layer (oce_ale_tracer.F90).
!
! Subroutines:
!   get_seawater_viscosity  — Sharqawy (2010) viscosity from T/S
!   get_particle_density    — Cram et al. (2018) bulk particle density
!   ballast                 — assemble scaling_density/viscosity arrays
!
! Usage (in callers):
!   use ballast_interface
!   use get_particle_density_interface
!   use get_seawater_viscosity_interface
!
! Origin:
!   Extracted from src/int_recom/recom_sinking.F90 lines 709-989.
!   Original implementation: O. Gürses, AWI Bremerhaven.
!===============================================================================
 
! ---------------------------------------------------------------------------
! Interface modules — keep the same pattern as the sinking interface modules
! already used in oce_ale_tracer.F90
! ---------------------------------------------------------------------------
 
module get_seawater_viscosity_interface
  interface
    subroutine get_seawater_viscosity(tr_num, tracers, partit, mesh)
      use mod_mesh
      use MOD_PARTIT
      use MOD_PARSUP
      use mod_tracer
      integer       , intent(in)   , target :: tr_num
      type(t_tracer), intent(inout), target :: tracers
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine
  end interface
end module
 
module get_particle_density_interface
  interface
    subroutine get_particle_density(tracers, partit, mesh)
      use mod_mesh
      use MOD_PARTIT
      use MOD_PARSUP
      use mod_tracer
      type(t_tracer), intent(inout), target :: tracers
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine
  end interface
end module
 
module ballast_interface
  interface
    subroutine ballast(tr_num, tracers, partit, mesh)
      use mod_mesh
      use MOD_PARTIT
      use MOD_PARSUP
      use mod_tracer
      integer       , intent(in)   , target :: tr_num
      type(t_tracer), intent(inout), target :: tracers
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine
  end interface
end module
 
! ---------------------------------------------------------------------------
! Implementations
! ---------------------------------------------------------------------------
 
!-------------------------------------------------------------------------------
! Subroutine to approximate seawater viscosity with current temperature
! based on Cram et al. (2018)
!
! Neglecting salinity effects, which are much smaller than those of temperature.
! https://bitbucket.org/ohnoplus/ballasted-sinking/src/master/tools/waterviscosity.m
!-------------------------------------------------------------------------------
subroutine get_seawater_viscosity(tr_num, tracers, partit, mesh)
 
  use recom_config
  use recom_glovar
    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_TRACER
    use o_PARAM
    use o_ARRAYS
    use g_CONFIG
    use g_forcing_arrays
    use g_comm_auto
    use g_clock
    use g_rotate_grid
 
    implicit none
 
!!  temp [degrees C] Ocean temperature
!!  salt [g/kg or n.d.] Ocean salinity
!!  seawater_visc_3D [kg m-1 s-1] Ocean viscosity
 
    real(kind=8),dimension(1)             :: A, B, mu_w
    integer                               :: row, k, nzmin, nzmax
 
    integer       , intent(in)   , target :: tr_num
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
 
#include "../../associate_part_def.h"
#include "../../associate_mesh_def.h"
#include "../../associate_part_ass.h"
#include "../../associate_mesh_ass.h"
 
    seawater_visc_3D(:,:) = 0.0
    do row=1,myDim_nod2d
     !if (ulevels_nod2D(row)>1) cycle
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)
 
        do k=nzmin, nzmax
     ! Eq from Sharaway 2010
     ! validity:
     !  0<temp<180 degC
     !  0<salt<0.15 kg/kg
     ! Note: because salinity is expected to be in kg/kg, use conversion factor 0.001 below!
            A(1) = 1.541 + 1.998*0.01*tracers%data(1)%values(k,row) - 9.52*1e-5*tracers%data(1)%values(k,row)*tracers%data(1)%values(k,row)
            B(1) = 7.974 - 7.561*0.01*tracers%data(1)%values(k,row) + 4.724*1e-4*tracers%data(1)%values(k,row)*tracers%data(1)%values(k,row)
            mu_w(1) = 4.2844*1.0e-5 + (1.0/(0.157*(tracers%data(1)%values(k,row)+64.993)*(tracers%data(1)%values(k,row)+64.993)-91.296))
            seawater_visc_3D(k,row) = mu_w(1) * (1.0 + A(1)*tracers%data(2)%values(k,row)*0.001 + B(1)*tracers%data(2)%values(k,row)*0.001*tracers%data(2)%values(k,row)*0.001)
        enddo
    end do
 
end subroutine get_seawater_viscosity
 
!-------------------------------------------------------------------------------
! Subroutine: calculate density of particle
! depending on composition (detC, detOpal, detCaCO3) based on Cram et al. (2018)
!-------------------------------------------------------------------------------
subroutine get_particle_density(tracers, partit, mesh)
 
    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_TRACER
 
    use recom_config
    use recom_glovar
    USE o_PARAM
    USE o_ARRAYS
    USE g_CONFIG
    use g_forcing_arrays
    use g_comm_auto
    use g_clock
    use g_rotate_grid
 
    implicit none
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
 
    integer                               :: row, k, nzmin, nzmax, tr_num, num_tracers
 
    real(kind=8)                          :: a1(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D) ! [n.d.] fraction of carbon in detritus class
    real(kind=8)                          :: a2(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D) ! [n.d.] fraction of nitrogen in detritus class
    real(kind=8)                          :: a3(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D) ! [n.d.] fraction of Opal in detritus class
    real(kind=8)                          :: a4(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D) ! [n.d.] fraction of CaCO3 in detritus class
    real(kind=8)                          :: b1(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=8)                          :: b2(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=8)                          :: b3(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=8)                          :: b4(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=8)                          :: aux(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
 
#include "../../associate_part_def.h"
#include "../../associate_mesh_def.h"
#include "../../associate_part_ass.h"
#include "../../associate_mesh_ass.h"
 
    num_tracers=tracers%num_tracers
 
    rho_particle1 = 0.0
    b1 = 0.0
    b2 = 0.0
    b3 = 0.0
    b4 = 0.0
    aux = 0.0
 
! Below guarantees non-negative tracer field
    do tr_num=1,num_tracers
        if (tracers%data(tr_num)%ID==1008)  b1 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetc      ! [mmol m-3] detritus carbon
        if (tracers%data(tr_num)%ID==1007)  b2 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetn      ! [mmol m-3] detritus nitrogen
        if (tracers%data(tr_num)%ID==1017)  b3 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetsi     ! [mmol m-3] detritus Si
        if (tracers%data(tr_num)%ID==1021)  b4 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetcal    ! [mmol m-3] detritus CaCO3
    end do
 
    do row=1,myDim_nod2d
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)
        aux(nzmin:nzmax,row) = b1(nzmin:nzmax,row)+b2(nzmin:nzmax,row)+b3(nzmin:nzmax,row)+b4(nzmin:nzmax,row)
        a1(nzmin:nzmax,row)  = b1(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a2(nzmin:nzmax,row)  = b2(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a3(nzmin:nzmax,row)  = b3(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a4(nzmin:nzmax,row)  = b4(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        rho_particle1(nzmin:nzmax,row) = rho_CaCO3*a4(nzmin:nzmax,row) + rho_opal*a3(nzmin:nzmax,row) + rho_POC*a1(nzmin:nzmax,row) + rho_PON*a2(nzmin:nzmax,row)
    end do
 
    if (enable_3zoo2det) then
    rho_particle2 = 0.0
    b1 = 0.0
    b2 = 0.0
    b3 = 0.0
    b4 = 0.0
    aux = 0.0
    do tr_num=1,num_tracers
        if (tracers%data(tr_num)%ID==1026)  b1 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetz2c
        if (tracers%data(tr_num)%ID==1025)  b2 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetz2n
        if (tracers%data(tr_num)%ID==1027)  b3 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetz2si
        if (tracers%data(tr_num)%ID==1028)  b4 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetz2calc
    end do
 
    do row=1,myDim_nod2d+eDim_nod2D   ! myDim is sufficient
        !if (ulevels_nod2D(row)>1) cycle
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)
        aux(nzmin:nzmax,row) = b1(nzmin:nzmax,row)+b2(nzmin:nzmax,row)+b3(nzmin:nzmax,row)+b4(nzmin:nzmax,row)
        a1(nzmin:nzmax,row)  = b1(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a2(nzmin:nzmax,row)  = b2(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a3(nzmin:nzmax,row)  = b3(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a4(nzmin:nzmax,row)  = b4(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        rho_particle2(nzmin:nzmax,row) = rho_CaCO3*a4(nzmin:nzmax,row) + rho_opal*a3(nzmin:nzmax,row) + rho_POC*a1(nzmin:nzmax,row) + rho_PON*a2(nzmin:nzmax,row)
    end do
    endif
 
end subroutine get_particle_density
 
!-------------------------------------------------------------------------------
! Subroutine: calculate ballasting
!
! Assembles the seawater-density- and viscosity-based sinking velocity scaling
! factors (scaling_density1/2_3D, scaling_visc_3D) that FESOM uses when
! computing tracer transport with use_ballasting = .true.
!-------------------------------------------------------------------------------
subroutine ballast(tr_num, tracers, partit, mesh)
 
    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_TRACER
 
    use recom_config
    use recom_glovar
 
    use o_PARAM
    use o_ARRAYS
    use g_CONFIG
    use g_forcing_arrays
    use g_comm_auto
    use g_clock
    use g_rotate_grid
    use mvars
    use mdepth2press
    use gsw_mod_toolbox, only: gsw_sa_from_sp,gsw_ct_from_pt,gsw_rho
 
    implicit none
    integer       , intent(in)   , target :: tr_num
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    integer                               :: row, k, nzmin, nzmax
    real(kind=8)                          :: depth_pos(1)
    real(kind=8)                          :: pres(1)
    real(kind=8)                          :: sa(1)
    real(kind=8)                          :: ct(1)
    real(kind=8)                          :: rho_seawater(1)
    real(kind=8)                          :: Lon_degree(1)
    real(kind=8)                          :: Lat_degree(1)
 
#include "../../associate_part_def.h"
#include "../../associate_mesh_def.h"
#include "../../associate_part_ass.h"
#include "../../associate_mesh_ass.h"
 
  ! For ballasting, calculate scaling factors here and pass them to FESOM,
  ! where sinking velocities are calculated.
  ! -----
  ! If ballasting is used, sinking velocities are a function of:
  !   a) particle composition (= density)
  !   b) sea water viscosity
  !   c) depth (currently for small detritus only)
  !   d) a constant reference sinking speed
  ! -----
 
     !___________________________________________________________________________
     ! loop over local nodes
     do row=1,myDim_nod2D
         ! max. number of levels at node n
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)
         !! lon
        Lon_degree(1)=geo_coord_nod2D(1,row)/rad !! convert from rad to degree
         !! lat
        Lat_degree(1)=geo_coord_nod2D(2,row)/rad !! convert from rad to degree
 
        ! get scaling vectors -> these need to be passed to FESOM to get sinking velocities
        ! get local seawater density
        do k=nzmin, nzmax
 
           !! level depth
           depth_pos(1) = abs(Z_3d_n(k,row))  ! take depth of tracers instead of levels abs(zbar_3d_n(k,row))
 
           call depth2press(depth_pos(1), Lat_degree(1), pres, 1)  ! pres is output of function,1=number of records
           sa           = gsw_sa_from_sp(tracers%data(2)%values(k,row), pres, Lon_degree(1), Lat_degree(1))
           ct           = gsw_ct_from_pt(sa,tracers%data(1)%values(k,row))
           rho_seawater = gsw_rho(sa, ct, pres)
 
           ! (i.e. no density scaling)
           scaling_density1_3D(k,row)=1.0
           scaling_density2_3D(k,row)=1.0
 
           if (use_density_scaling) then
              scaling_density1_3D(k,row) = (rho_particle1(k,row)-rho_seawater(1))/(rho_ref_part-rho_ref_water)
              scaling_density2_3D(k,row) = (rho_particle2(k,row)-rho_seawater(1))/(rho_ref_part-rho_ref_water)
           endif
 
            scaling_visc_3D(k,row)=1.0
 
            if (use_viscosity_scaling) then
                if (seawater_visc_3D(k,row)==0) then
                    scaling_visc_3D(k,row)=1.0
                else
                    scaling_visc_3D(k,row)= visc_ref_water/seawater_visc_3D(k,row)
                endif
            endif
 
        end do
        rho_particle1(nzmax+1,row) = rho_particle1(nzmax,row)
        rho_particle2(nzmax+1,row) = rho_particle2(nzmax,row)
        scaling_visc_3D(nzmax+1,row) = scaling_visc_3D(nzmax,row)
     end do
    ! in the unlikely case that rho_particle(k)-rho_seawater(1)<0,
    ! prevent the scaling factor from being negative
    if (any(scaling_density1_3D(:,:) <= tiny)) scaling_density1_3D(:,:) = 1.0_WP
 
    if (enable_3zoo2det) then
        if (any(scaling_density2_3D(:,:) <= tiny)) scaling_density2_3D(:,:) = 1.0_WP
    endif
 
end subroutine ballast
