module runoff_scaling_interface
  interface
    subroutine runoff_scaling_init()
    end subroutine runoff_scaling_init

    subroutine runoff_scaling(runoff_in, partit, mesh)
      use MOD_PARTIT
      use MOD_MESH
      use o_PARAM, only: WP
      real(kind=WP), intent(inout) :: runoff_in(:)
      type(t_partit), intent(inout), target :: partit
      type(t_mesh),   intent(in), target :: mesh
    end subroutine runoff_scaling
  end interface
end module runoff_scaling_interface

module runoff_scaling_state
  use o_PARAM , only: WP
  implicit none
  save

  logical :: loaded = .false.
  real(kind=WP), allocatable :: runoff_ref(:)
end module runoff_scaling_state

! read reference file and apply values to callable variable
! For now just monthly reference
! runoff reference filename required to be 'runoff_ref_YYYY.dat'
! file structure: 12 lines with monthly ref value in Sv

subroutine runoff_scaling_init()
  use g_config
  use g_clock
  use runoff_scaling_state

  implicit none
  
  character(len=300)        :: filename
  integer                   :: iunit, ios, m
  !real(kind=WP),allocatable :: runoff_ref(:)
  !logical                   :: loaded = .false.

  select case (trim(runoff_scaling_method))

  case ('ref')
    allocate(runoff_ref(12))
    write(filename, '(A,"/runoff_ref_",I4.4,".dat")') trim(runoff_dir), yearnew

    iunit = 99
    open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR: Could not open runoff reference file:", trim(filename)
      stop
    end if
    
    ! if (trim(runoff_scaling_time) == 'y') then
    ! write cases for y/m/d -> how do leap years work?
    ! as reference and model years should be the same, i think i just need to find the leap year flag and set
    ! one do to 365 and one to 366

    do m = 1, 12
      read(iunit, *, iostat=ios) runoff_ref(m)
      if (ios /= 0) then
        write(*,*) "ERROR reading month ", m, " in ", trim(filename)
        stop
      end if
    end do
    close(iunit)

    loaded = .true.
    write(*,*) "Runoff scaling: Loaded ", trim(filename)
    write (*,*) "Reference runoff: ", runoff_ref

  case ('const')

    write(*,*) "Case const"
    ! runoff_ref = runoff_const_ref

  case ('mult')
  
    write(*,*) "Case mult"
    ! runoff_ref = runoff_mult_factor

  case default
  
    write(*,*) "Case default"

  end select

end subroutine runoff_scaling_init

! take the runoff for each timestep, calc the spatial integral for the southern ocean, compare with
! the reference for the appropriate time, and finally apply the scaling factor from the integral
! difference to each SO node
!
! For now just monthly reference

subroutine runoff_scaling(runoff_in, partit, mesh)
  use MOD_PARTIT
  use MOD_MESH
  use g_support
  use g_forcing_arrays
  use g_clock
  use runoff_scaling_state

  implicit none

  real(kind=WP), intent(inout)          :: runoff_in(:)
  integer                               :: i
  !real(kind=WP), allocatable            :: runoff_ref(:)
  real(kind=WP), allocatable            :: runoff_mask(:)
  real(kind=WP)                         :: runoff_south
  real(kind=WP)                         :: runoff_factor
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in), target    :: mesh
  !logical                               :: loaded
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  if (mype==0) write (*,*) "In runoff_scaling()"
  if (mype==0) write (*,*) "size of runoff_in:   ", size(runoff_in)

  allocate(runoff_mask(size(runoff_in)))
  
  if (mype==0) write (*,*) "allocated runoff_mask"
  if (mype==0) write (*,*) "size of runoff_mask: ", size(runoff_mask)

  runoff_mask = 0.0_WP
  if (mype==0) write (*,*) "size of runoff_mask: ", size(runoff_mask)

  if (mype==0) write (*,*) "set runoff_mask to 0"

  where (geo_coord_nod2D(2, :) < -60.0_WP * rad)
    runoff_mask = runoff_in
  end where

  if (mype==0) write (*,*) "size of geo_coord_nod2D: ", size(geo_coord_nod2D(2, :))
  if (mype==0) write (*,*) "size of runoff_mask: ", size(runoff_mask)
  if (mype==0) write (*,*) "size of runoff_in:   ", size(runoff_in)
  if (mype==0) write (*,*) "set runoff_mask to runoff_in fuer < -60"

  call integrate_nod(runoff_mask, runoff_south, partit, mesh)
 
  if (mype==0) write (*,*) "runoff south: ", runoff_south

  select case (trim(runoff_scaling_method))

  case ('ref')

    if (.not. loaded) then
      write(*,*) "ERROR: runoff_scaling called before initialization."
      stop
    end if

    if (runoff_south == 0.0_WP) then
      runoff_factor = 1.0_WP
      if (mype==0) write(*,*) "Warning: runoff_south = 0, skipping scaling, setting factor to 1.0"
    else
      runoff_factor = runoff_ref(month) * 1.0e6_WP / runoff_south
      if (mype==0) write(*,*) "runoff factor:", runoff_factor
    end if
  
    !$OMP PARALLEL DO
    do i = 1, myDim_nod2D + eDim_nod2D
      if (geo_coord_nod2D(2, i) < -60.0_WP * rad) then
        runoff_in(i) = runoff_in(i) * runoff_factor
      end if
    end do
    !$OMP END PARALLEL DO

  case ('const')

    if (runoff_south == 0.0_WP) then
      runoff_factor = 1.0_WP
      if (mype==0) write(*,*) "Warning: runoff_south = 0, skipping scaling, setting factor to 1.0"
    else
      runoff_factor = runoff_const_ref * 1.0e6_WP / runoff_south
      if (mype==0) write(*,*) "runoff factor:", runoff_factor
    end if

    !$OMP PARALLEL DO
    do i = 1, myDim_nod2D + eDim_nod2D
      if (geo_coord_nod2D(2, i) < -60.0_WP * rad) then
        runoff_in(i) = runoff_in(i) * runoff_factor
      end if
    end do
    !$OMP END PARALLEL DO

  case ('mult')

    !$OMP PARALLEL DO
    do i = 1, myDim_nod2D + eDim_nod2D
      if (geo_coord_nod2D(2, i) < -60.0_WP * rad) then
        runoff_in(i) = runoff_in(i) * runoff_mult_factor
      end if
    end do
    !$OMP END PARALLEL DO

  case default
  
    write(*,*) "Ended up in default case in runoff_scaling"

  end select
    
  deallocate(runoff_mask)
end subroutine runoff_scaling

! To Do:
! - how does file read?
! - check how month from g_clock actually works
! - check runoff scaling
! - 
!
! Additional ideas/later to dos:
! - set the scalling to daily/monthly/yearly
! - allow for constant fwf addition -> set hosing value to be applied only on runoff nodes
!

