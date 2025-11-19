module runoff_scaling_interface
  interface
    subroutine runoff_scaling_init(runoff_dir)
      use g_clock
      character(*), intent(in) :: runoff_dir
    end subroutine runoff_scaling_init

    subroutine runoff_scaling(runoff, month, partit, mesh)
      use g_glock
      real(WP), intent(inout) :: runoff(:)
      integer, intent(in) :: month
      type(t_partit), intent(in), target :: partit
      type(t_mesh), intent(in), target :: mesh
    end subroutine runoff_scaling
  end interface
end module runoff_scaling_interface

! read reference file and apply values to callable variable
! For now just monthly reference
! runoff reference filename required to be 'runoff_ref_YYYY.dat'
! file structure: 12 lines with monthly ref value in Sv

subroutine runoff_scaling_init()
  use g_clock

  implicit none

  character(len=*)   :: runoff_dir
  character(len=300) :: file_name
  integer            :: inunit, ios, m

  write(filename, '(A,"/runoff_ref_",I4.4,".dat")') trim(runoff_dir), yearnew

  iunit = 99
  open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(*,*) "ERROR: Could not open runoff reference file:", trim(filename)
    stop
  end if

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

end subroutine runoff_scaling_init

! take the runoff for each timestep, calc the spatial integral for the southern ocean, compare with
! the reference for the appropriate time, and finally apply the scaling factor from the integral
! difference to each SO node
!
! For now just monthly reference

subroutine runoff_scaling(runoff)
  use MOD_PARTIT
  use MOD_MESH
  use g_support
  use g_forcing_arrays
  use g_clock

  implicit none

  real(kind=WP), intent(inout) :: runoff(:)
  integer                      :: i
  real(kind=WP), allocatable :: runoff_ref
  real(kind=WP), allocatable :: runoff_mask(:)
  real(kind=WP) :: runoff_south
  real(kind=WP) :: runoff_factor
  type(t_partit), intent(in), target :: partit
  type(t_mesh),   intent(in), target :: mesh

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  allocate(runoff_mask(size(runoff)))

  if (.not. loaded) then
    write(*,*) "ERROR: runoff_scaling called before initialization."
    stop
  end if
  
  runoff_mask = 0.0_WP

  where (geo_coord_nod2D(2, :) < -60.0_WP * rad)
    runoff_masked = runoff
  end where
  call integrate_nod(runoff_mask, runoff_south, partit, mesh)

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
      runoff(i) = runoff(i) * runoff_factor
    end if
  end do
  !$OMP END PARALLEL DO

  deallocate(runoff_mask)
end subroutine runoff_scaling()

! To Do:
! - where do I need to add oce_runoff_scaling to have it compile with the model?
! - how does file read?
! - check how month from g_clock actually works
! - clean other files from my previous attempts
! - check runoff scaling
! - 
!
! Additional ideas/later to dos:
! - set the scalling to daily/monthly/yearly
! - set different flags for ref or factor scaling
! - allow for constant fwf addition -> set hosing value to be applied only on runoff nodes
!

