module age_tracer_init_interface
  interface
    subroutine age_tracer_init(partit, mesh)
      use mod_partit
      use mod_mesh
      type(t_partit), intent(in), target :: partit
      type(t_mesh),   intent(in), target :: mesh
    end subroutine age_tracer_init
  end interface
end module age_tracer_init_interface

subroutine age_tracer_init(partit, mesh)
  use MOD_PARTIT
  use MOD_MESH
  use g_comm_auto
  use g_forcing_param
  use g_forcing_arrays

  implicit none

  integer                     :: n, i, j, num_reg, num_nod, num_days
  integer                     :: n_loc, fileID
  integer, allocatable        :: temp_arr2d(:), nodes_release(:)
  character*300               :: file_name
  type(t_partit), intent(in), target :: partit
  type(t_mesh), intent(in) , target :: mesh

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  if (use_age_mask) then
     allocate(temp_arr2d(nod2D))
     temp_arr2d=0
     do n=1, myDim_nod2D  !note: eDim_nod2D should not be included in this case
        temp_arr2d(myList_nod2D(n))=n
     end do

     file_name=trim(age_tracer_path)//'age_tracer_release_nodes.out'
     fileID=160
     open(fileID, file=file_name)
     read(fileID,*) num_nod
     allocate(nodes_release(num_nod))
     read(fileID,*) nodes_release
     close(fileID)

     do n=1,num_nod
        n_loc=temp_arr2d(nodes_release(n))
        if(n_loc>0) then
           age_tracer_loc_index(n_loc)=1
        end if
     end do

     deallocate(nodes_release)
     deallocate(temp_arr2d)

     if(mype==0) write(*,*) 'Age tracer file prepared: ', file_name
     if(mype==0) write(*,*) 'index_age_tracer: ', index_age_tracer 
  end if

end subroutine age_tracer_init

