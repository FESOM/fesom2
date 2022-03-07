module iceberg_water_init_interface
  interface
    subroutine iceberg_water_init(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module

subroutine iceberg_water_init(mesh)
  ! init iceberg melting rate
  use MOD_MESH
  use o_MESH
  USE g_CONFIG
  use o_ARRAYS
  use i_ARRAYS
  use g_comm_auto
  use g_PARSUP
  use o_param
  use g_support
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_clock

  implicit none

  integer                     :: n, i, j, num_nod1, num_nod, num_days
  integer                     :: n_loc, fileID
  integer, allocatable        :: temp_arr2d(:), nodes_in_region(:)
  real(kind=8), allocatable   :: iceberg_flux(:)
  character*300               :: file_name
  type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

  allocate(temp_arr2d(nod2D))
  temp_arr2d=0
  do n=1, myDim_nod2D  !note: eDim_nod2D should not be included in this case
     temp_arr2d(myList_nod2D(n))=n
  end do

  ! read the iceberg discharge (m/s) in each region:
  file_name=trim(iceberg_path)//'iceberg_flux.out'
  fileID=160
  open(fileID, file=file_name)
  read(fileID,*) num_nod1
  allocate(iceberg_flux(num_nod1))
  read(fileID,*) iceberg_flux    !Flux(m/s)=iceberg_dischanges(pism)/86400/365/1000 
  close(fileID)




     !read in nodes in the region
     file_name=trim(iceberg_path)//'iceberg_nodes.out'
     fileID=160
     open(fileID, file=file_name)
     read(fileID,*) num_nod
     allocate(nodes_in_region(num_nod))
     read(fileID,*) nodes_in_region
     close(fileID)


     do n=1,num_nod
        n_loc=temp_arr2d(nodes_in_region(n))
        if(n_loc>0) then
           runoff_iceberg(n_loc)=iceberg_flux(n)
        end if
     end do

     deallocate(nodes_in_region)
     deallocate(iceberg_flux)

  iceberg_season=0.0
  iceberg_season(iceberg_start_mon:iceberg_end_mon)=1.0

  deallocate(temp_arr2d)

  if(mype==0) write(*,*) 'iceberg melt water fluxes prepared.'
  if(mype==0) write(*,*) 'iceberg_season=', iceberg_season

end subroutine iceberg_water_init

