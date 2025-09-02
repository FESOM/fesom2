module landice_water_init_interface
  interface
    subroutine landice_water_init(partit, mesh)
      use mod_partit
      use mod_mesh
      type(t_partit), intent(in), target :: partit
      type(t_mesh),   intent(in), target :: mesh
    end subroutine landice_water_init
  end interface
end module landice_water_init_interface

subroutine landice_water_init(partit, mesh)
  ! init land ice melting rate
  use MOD_PARTIT
  use MOD_MESH
  use o_PARAM , only: WP
  use g_comm
  use g_forcing_param
  use g_forcing_arrays

  implicit none

  integer                     :: n, i, j, num_reg, num_nod, num_days
  integer                     :: n_loc, fileID
  integer, allocatable        :: temp_arr2d(:), nodes_in_region(:)
  real(kind=WP)               :: vol, vol_sum, aux
  real(kind=WP), allocatable  :: totalflux(:)
  character*300               :: file_name
  character                   :: c_reg_ind
  type(t_partit), intent(in), target :: partit
  type(t_mesh),   intent(in), target :: mesh
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  ! read the number of region and the total yearly discharge in each region:
  file_name=trim(fwf_path)//'landice_yearly_mass_loss.out'
  fileID=160
  open(fileID, file=file_name)
  read(fileID,*) num_reg
  allocate(totalflux(num_reg))
  read(fileID,*) totalflux      !Surface Flux (m/s) = water flux (Sv)*1e6/area_sum
  close(fileID)


  allocate(temp_arr2d(nod2D))
  temp_arr2d=0
  do n=1, myDim_nod2D  !note: eDim_nod2D should not be included in this case
     temp_arr2d(myList_nod2D(n))=n
  end do

  do i=1,num_reg

     !read in nodes in the region
     write(c_reg_ind,'(i1)') i   !assume that num_reg is less than 10
     file_name=trim(fwf_path)//'landice_nodes_in_region_'//c_reg_ind//'.out'
     fileID=160
     open(fileID, file=file_name)
     read(fileID,*) num_nod
     allocate(nodes_in_region(num_nod))
     read(fileID,*) nodes_in_region
     close(fileID)


     aux=totalflux(i)  !m/s

     do n=1,num_nod
        n_loc=temp_arr2d(nodes_in_region(n))
        if(n_loc>0) then
           runoff_landice(n_loc)=aux
        end if
     end do

     deallocate(nodes_in_region)
  end do

  landice_season=0.0
  landice_season(landice_start_mon:landice_end_mon)=1.0

  deallocate(temp_arr2d, totalflux)

  if(mype==0) write(*,*) 'Land-ice melt water fluxes prepared.'
  if(mype==0) write(*,*) 'freshwater flux (m/s)=', aux
  if(mype==0) write(*,*) 'landice_season=', landice_season

end subroutine landice_water_init

