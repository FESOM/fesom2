! Load river biogeochemical variables from NetCDF file R2OMIP
! Variables: DIC, DIN, DOCl, DOCs, POC

! Subroutine to read and broadcast a single variable

subroutine load_river_variable(ncid, vari, model_2Darray, partit, mesh)  

  use, intrinsic :: ISO_FORTRAN_ENV, only: real64

  use g_config
  use o_param
  USE MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_rotate_grid
  use netcdf
  implicit none

  type(t_mesh),   intent(in)   , target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer                       :: n
  integer                       :: status, ncid, varid
  integer                       :: RIstart, RIcount

  real(real64), allocatable :: ncdata(:)

  real(real64),  intent(out)	:: model_2Darray(partit%myDim_nod2D+partit%eDim_nod2D)

  character(*),  intent(in)     :: vari
  integer                       :: ierror

#include "../../associate_part_def.h"
#include "../../associate_mesh_def.h"
#include "../../associate_part_ass.h"
#include "../../associate_mesh_ass.h"

    ! Allocate temporary array for full domain
    allocate(ncdata(mesh%nod2D))

    if (mype == 0) then
        status = nf90_inq_varid(ncid, trim(vari), varid)
        RIstart = 1
        RIcount = nod2D
        status = nf90_get_var(ncid, varid, ncdata, start=(/RIstart/), count=(/RIcount/))
    endif

    ! Broadcast to all processes
    call MPI_BCast(ncdata, nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

    ! Extract local partition
    model_2Darray=ncdata(myList_nod2D)

    ! Clean up
    deallocate(ncdata)

end subroutine load_river_variable

!-----------------------------------------------------------------------------
! Subroutine: river_sanity_check
!
! Purpose: Compute and print global min/max of a 2D riverine input array
!          across all MPI ranks as a quick plausibility check after reading.
!
! Arguments:
!   river_var  - local 2D array of riverine input values
!   varname    - short label printed in the diagnostic output
!-----------------------------------------------------------------------------
!subroutine river_sanity_check(river_var, varname, partit, mesh)
subroutine river_sanity_check(locmax, locmin, varname, partit, mesh)

    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_MESH
    !use mpi
 
    implicit none

    ! Input parameters
    type(t_partit), intent(inout),   target          :: partit
    type(t_mesh), intent(inout),     target          :: mesh

    !real(8),          intent(in) :: river_var(:)
    character(len=*), intent(in) :: varname
 
    real(8) :: locmax, locmin, glo
    !integer :: MPIerr


#include "../../associate_part_def.h"
#include "../../associate_mesh_def.h"
#include "../../associate_part_ass.h"
#include "../../associate_mesh_ass.h"

    !locmax = maxval(river_var)
    !locmin = minval(river_var)
 
    call MPI_AllReduce(locmax, glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    if (mype == 0) write(*,'(A,A,A,ES14.6)') '  |-> global max riverine ', trim(varname), '. = ', glo
 
    call MPI_AllReduce(locmin, glo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
    if (mype == 0) write(*,'(A,A,A,ES14.6)') '  |-> global min riverine ', trim(varname), '. = ', glo
 
end subroutine river_sanity_check


