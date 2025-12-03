!==========================================================
! CMOR Variables Diagnostics Module for CMIP6/CMIP7 output
! Ported from FESOM 1.4, adapted for FESOM 2 level-based grid
!==========================================================
module cmor_variables_diag
  use o_PARAM
  use MOD_MESH
  use MOD_PARTIT
  use MOD_PARSUP
  use MOD_TRACER
  use MOD_DYN
  use MOD_ICE
  use o_ARRAYS
  use g_config
  use g_comm_auto
  
  implicit none
  
  public :: init_cmor_diag, compute_cmor_diag, lcmor_diag
  public :: volo, opottemptend, pbo, soga, thetaoga, tos, sos
  public :: siarean, siareas, siextentn, siextents, sivoln, sivols
  
  private
  
  ! Control flag for CMOR diagnostics
  logical :: lcmor_diag = .false.
  
  ! CMOR diagnostic variables
  logical, save :: initialized = .false.
  logical, save :: first_call = .true.
  real(kind=WP), save :: volo                            ! Ocean volume [m^3]
  real(kind=WP), save, allocatable :: opottemptend(:)    ! Ocean potential temperature tendency [W/m^2] - 2D field
  real(kind=WP), save, allocatable :: pbo(:)             ! Sea water pressure at sea floor [Pa] - 2D field
  real(kind=WP), save, allocatable :: tos(:)             ! Sea surface temperature [degC] - 2D field
  real(kind=WP), save, allocatable :: sos(:)             ! Sea surface salinity [psu] - 2D field
  real(kind=WP), save :: soga, thetaoga                  ! Global mean salinity and temperature
  real(kind=WP), save :: siarean, siareas                ! Sea ice area north/south [10^12 m^2]
  real(kind=WP), save :: siextentn, siextents            ! Sea ice extent north/south [10^12 m^2]
  real(kind=WP), save :: sivoln, sivols                  ! Sea ice volume north/south [10^9 m^3]
  
  ! Auxiliary arrays for temperature tendency
  ! FESOM 2 uses (level, nod2D) indexing
  real(kind=WP), save, allocatable :: previous_temp(:,:)  ! (nl-1, myDim_nod2D)
  
  namelist /cmor_diag/ lcmor_diag

contains

  !=================================================================
  ! Initialize CMOR diagnostics
  !=================================================================
  subroutine init_cmor_diag(partit, mesh)
    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    integer :: n2, k, ierr
    real(kind=WP) :: volo_local
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    if (mype == 0) then
        write(*,*) '____________________________________________________________'
        write(*,*) ' --> Initializing CMOR diagnostics for CMIP6/CMIP7'
    end if

    ! Compute ocean volume by summing over all levels and 2D nodes
    ! FESOM 2: loop over (level, nod2D), not nod3D
    volo_local = 0.0_WP
    do n2 = 1, myDim_nod2D
       do k = ulevels_nod2D(n2), nlevels_nod2D(n2)-1
          volo_local = volo_local + area(k, n2) * hnode(k, n2)
       end do
    end do
    call MPI_AllREDUCE(volo_local, volo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
  
    ! Allocate 2D arrays (indexed by nod2D)
    allocate(opottemptend(myDim_nod2D))
    allocate(pbo(myDim_nod2D))
    allocate(tos(myDim_nod2D))
    allocate(sos(myDim_nod2D))
    
    ! Allocate previous temperature array (level, nod2D)
    allocate(previous_temp(nl-1, myDim_nod2D))
    
    opottemptend = 0.0_WP
    pbo = 0.0_WP
    tos = 0.0_WP
    sos = 0.0_WP
    previous_temp = 0.0_WP

    initialized = .true.
    first_call = .true.
    
    if (mype == 0) then
        write(*,*) '    Global ocean volume = ', volo, ' m^3'
        write(*,*) '    CMOR diagnostics initialized'
    end if
    
  end subroutine init_cmor_diag


  !=================================================================
  ! Compute CMOR diagnostics
  !=================================================================
  subroutine compute_cmor_diag(tracers, ice, dynamics, partit, mesh)
    implicit none
    type(t_tracer), intent(in),    target :: tracers
    type(t_ice),    intent(in),    target :: ice
    type(t_dyn),    intent(in),    target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    
    integer :: n2, k, ierr, ku, kl
    real(kind=WP) :: z_up, z_lo, dens_val, dz
    real(kind=WP), dimension(:,:), pointer :: temp, salt
    real(kind=WP), dimension(:),   pointer :: eta_n
    
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    if (.not. initialized) call init_cmor_diag(partit, mesh)
    
    ! Point to tracer data: FESOM 2 uses (level, nod2D) indexing
    temp => tracers%data(1)%values(:,:)
    salt => tracers%data(2)%values(:,:)
    eta_n => dynamics%eta_n(:)
    
    ! Initialize previous_temp on first call
    if (first_call) then
        do n2 = 1, myDim_nod2D
            do k = ulevels_nod2D(n2), nlevels_nod2D(n2)-1
                previous_temp(k, n2) = temp(k, n2)
            end do
        end do
        first_call = .false.
    end if
    
    ! Initialize arrays
    opottemptend = 0.0_WP
    pbo = 0.0_WP
    
    !=================================================================
    ! Compute bottom pressure and temperature tendency
    ! FESOM 2: loop over 2D nodes, then levels
    !=================================================================
    do n2 = 1, myDim_nod2D
        ku = ulevels_nod2D(n2)  ! Upper level (surface or cavity)
        kl = nlevels_nod2D(n2)  ! Bottom level + 1
        
        ! Bottom pressure: SSH contribution
        pbo(n2) = eta_n(n2) * density_0 * g
        
        ! Loop through water column levels
        z_up = zbar_3d_n(ku, n2)  ! Top of first layer
        
        do k = ku, kl-1
            z_lo = zbar_3d_n(k+1, n2)  ! Bottom of current layer
            dz = abs(z_up - z_lo)
            
            ! Use density_m_rho0 if available, else use density_0
            if (allocated(density_m_rho0)) then
                dens_val = density_0 + density_m_rho0(k, n2)
            else
                dens_val = density_0
            end if
            
            ! Bottom pressure: steric contribution (accumulate through column)
            pbo(n2) = pbo(n2) + g * dz * dens_val
            
            ! Temperature tendency: dT/dt * rho * Cp * dz
            ! Tendency is expressed as rate of change of heat content per unit area
            opottemptend(n2) = opottemptend(n2) + &
                (temp(k, n2) - previous_temp(k, n2)) / dt * vcpw * hnode(k, n2)
            
            z_up = z_lo
        end do
    end do
    
    !=================================================================
    ! Compute global mean salinity and temperature
    ! Volume-weighted average over all levels
    !=================================================================
    soga = 0.0_WP
    thetaoga = 0.0_WP
    do n2 = 1, myDim_nod2D
        do k = ulevels_nod2D(n2), nlevels_nod2D(n2)-1
            soga = soga + salt(k, n2) * area(k, n2) * hnode(k, n2)
            thetaoga = thetaoga + temp(k, n2) * area(k, n2) * hnode(k, n2)
        end do
    end do
    call MPI_AllREDUCE(MPI_IN_PLACE, soga, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, thetaoga, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    soga = soga / volo
    thetaoga = thetaoga / volo
  
    !=================================================================
    ! Compute sea surface temperature and salinity
    ! Use uppermost wet level for each node
    !=================================================================
    do n2 = 1, myDim_nod2D
        k = ulevels_nod2D(n2)  ! Surface level (or first wet level in cavity)
        tos(n2) = temp(k, n2)
        sos(n2) = salt(k, n2)
    end do
    
    !=================================================================
    ! Compute sea ice diagnostics
    !=================================================================
    siarean = 0.0_WP
    siareas = 0.0_WP
    siextentn = 0.0_WP
    siextents = 0.0_WP
    sivoln = 0.0_WP
    sivols = 0.0_WP
    
    do n2 = 1, myDim_nod2D
        if (coord_nod2D(2, n2) > 0.0_WP) then
            ! Northern hemisphere
            siarean = siarean + ice%data(1)%values(n2) * area(1, n2)
            if (ice%data(1)%values(n2) >= 0.15_WP) then
                siextentn = siextentn + area(1, n2)
            end if
            sivoln = sivoln + ice%data(2)%values(n2) * area(1, n2)
        else
            ! Southern hemisphere
            siareas = siareas + ice%data(1)%values(n2) * area(1, n2)
            if (ice%data(1)%values(n2) >= 0.15_WP) then
                siextents = siextents + area(1, n2)
            end if
            sivols = sivols + ice%data(2)%values(n2) * area(1, n2)
        end if
    end do
    
    call MPI_AllREDUCE(MPI_IN_PLACE, siarean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, siareas, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, siextentn, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, siextents, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, sivoln, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, sivols, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    
    ! Convert to proper units
    siarean = siarean / 1.0e12_WP      ! to 10^12 m^2
    siareas = siareas / 1.0e12_WP      ! to 10^12 m^2
    siextentn = siextentn / 1.0e12_WP  ! to 10^12 m^2
    siextents = siextents / 1.0e12_WP  ! to 10^12 m^2
    sivoln = sivoln / 1.0e9_WP         ! to 10^9 m^3
    sivols = sivols / 1.0e9_WP         ! to 10^9 m^3
    
    ! Update previous tracer values for next time step
    do n2 = 1, myDim_nod2D
        do k = ulevels_nod2D(n2), nlevels_nod2D(n2)-1
            previous_temp(k, n2) = temp(k, n2)
        end do
    end do
    
  end subroutine compute_cmor_diag

end module cmor_variables_diag
