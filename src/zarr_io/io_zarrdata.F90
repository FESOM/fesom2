module io_zarrdata

  use o_PARAM, only : WP
  use, intrinsic :: iso_fortran_env, only: real64, real32
  use io_data_strategy_module
  use async_threads_module
  implicit none
  public:: output_zarr
  private
contains

subroutine output_zarr(istep)
  use g_clock
  use mod_mesh
  use g_PARSUP
  use io_gather_module
  use o_arrays
#if defined (__icepack)
  use icedrv_main,    only: init_io_icepack
#endif
  use zarr_util

  implicit none

  integer       :: istep
  logical, save :: lfirst=.true.
  integer       :: n, k, i, nsteps
  logical       :: do_output
  !type(t_mesh), intent(in) , target :: mesh
  character(:), allocatable :: filepath
  real(real64)                  :: rtime !timestamp of the record
  character(len=9999) :: root, filename, dirnames, text, zgroup, zarray, tmp
  character(len=200) :: variable_path
  character(len=5) :: mype_key
  integer  :: fileunit
  real(real32) :: rdata(myDim_nod2d)
  integer, save :: counter=0

!  type(Meandata), pointer :: entry
!
  filename = "test2.zarr"
  dirnames = trim(root) // "/sst"  
  if (lfirst) then
     if (mype==0) then
      root = create_zarr_group(trim(filename))
     else
      root = filename 
     end if 
     ! not necessary
     !variable_path = create_zarr_group('sst',root) ! also adds dirs for each pe to hold its arrays
     ! testing
     !call write_zarr_array(istep, 'sst', reshape(tr_arr(1,1:myDim_nod2d,1),[size(tr_arr(1,1:myDim_nod2d,1))]), filename)
     lfirst=.false.
  end if

  call write_zarr_array(istep, 'sst', reshape(tr_arr(1,1:myDim_nod2d,1),[size(tr_arr(1,1:myDim_nod2d,1))]), filename)

end subroutine output_zarr

function create_zarr_group(group_name, root) result (group)
  use zarr_util 
  use g_PARSUP
  implicit none
  character(len=*), intent(in)               :: group_name
  character(len=*), intent(in), optional     :: root
  character(200)               :: tmp, group ! in future group is a type?
  character(9999)               :: text
  integer        :: fileunit
  logical        :: file_exists
  
  ! need to handle restart run
  ! is inquire to check if group exists better then mkdir 
  if (present(root)) then
     group = trim(root) // '/' //trim(group_name)
  else
     group = trim(group_name) 
  end if
  call execute_command_line ('mkdir -p ' // adjustl(trim(group)))

  tmp=trim(group) // "/.zgroup"
  INQUIRE(FILE=trim(tmp), EXIST=file_exists)  
  if (mype == 0) then
    if (.not.file_exists) then  
      open(newunit=fileunit, file=trim(tmp))
      call get_zgroup_json(text)
      write(fileunit, '(a)') trim(text)
      flush(fileunit)
      close(fileunit)

      tmp=trim(group)//"/.zattrs"
      open(newunit=fileunit, file=trim(tmp))
      call get_zattrs_json(text)
      write(fileunit, '(a)') trim(text)
      flush(fileunit)
      close(fileunit)
    end if
  end if


end function create_zarr_group



!writes and appends in time
subroutine write_zarr_array(istep, variable, rarray, root) !, mensions, add_time_dim)
  use zarr_util 
  use g_PARSUP
  implicit none
  character(len=*), intent(in)               :: variable, root
  integer, intent(in)           :: istep
  character(1000)               :: variable_group, group
  character(9999)               :: text 
  real(kind=WP), intent(in) :: rarray(:) 
  character(2000)               ::  tmp
  integer          :: data_kind, fileunit, nsteps
  character(1), parameter :: endianess="<" ! later find from func
  logical,save :: lfirst=.true.
  character(5), save :: mype_string
  character(10) :: dtype
  !logical                 :: add_time_dim=.true.
  !character(len=*), intent(in) :: dimensions(:)
   

  if (lfirst) then
     variable_group = create_zarr_group(variable, root) 
     
     write(mype_string,'(I0)') mype
     tmp=trim(variable_group) // "/" // trim(mype_string) 
     !should we inquire? before?
     call execute_command_line ('mkdir -p ' // adjustl(trim(tmp)))

     data_kind = sizeof(rarray)/size(rarray) !use storage_size intrensic?
     write(tmp, '(I0)') data_kind
     dtype = byteorder() // "f" // tmp
     
     call get_run_steps(nsteps)

     !should we inquire? before or replace?
     tmp=trim(variable_group) // "/" // trim(mype_string) 
     tmp=trim(tmp) // "/.zarray"
     open(newunit=fileunit, file=trim(tmp)) 
     !call get_zarray_json(text, shape(rarray), shape(rarray), dtype) !dtype)
     !call get_zarray_json(text, [1,shape(rarray)], [nsteps, shape(rarray)], dtype)
     call get_zarray_json(text, [shape(rarray),1], [shape(rarray),nsteps], dtype)
     write(fileunit, '(a)') trim(text)
     flush(fileunit)
     close(fileunit)

     tmp=trim(variable_group) // "/" // trim(mype_string) 
     tmp=trim(tmp) // "/.zattrs"
     open(newunit=fileunit, file=trim(tmp))
     !text = '{"_ARRAY_DIMENSIONS":["nod2d"]}'
     !text = '{"_ARRAY_DIMENSIONS":["time","nod2d"]}'
     text = '{"_ARRAY_DIMENSIONS":["nod2d", "time"]}' ! works for dask.array.concatinate
     write(fileunit, '(a)') trim(text)
     flush(fileunit)
     close(fileunit)
     lfirst = .false. 
  end if  
  
  write(tmp, '(I0)') istep-1  
  !tmp=trim(root) // '/' //trim(variable) // "/" // trim(mype_string) // '/0'
  !tmp=trim(root) // '/' //trim(variable) // "/" // trim(mype_string) // '/' // trim(tmp)// '.0'
  tmp=trim(root) // '/' //trim(variable) // "/" // trim(mype_string) // '/0.' // trim(tmp)
  call write_real_array(tmp, rarray)

end subroutine write_zarr_array



subroutine write_real_array(path, rarray)
      use o_PARAM 
      character(len=*), intent(in) :: path
      !integer, intent(in) :: array_kind !can we auto decipher with sizeof(arr)/size(arr)      !and use king 
      !real(kind=selected_real_kind(array_kind)), intent(in) :: rarray 
      real(kind=WP), intent(in) :: rarray(:) 
      integer  :: fileunit
      ! later when zarr and fesom are seperated fesom's wrapper should handle default kind and pe info
      open(newunit=fileunit, file=trim(path), form="unformatted", access='stream', action='write') !, recl=rec_array)
      write(fileunit) rarray
      flush(fileunit)
      close(fileunit)
           
end subroutine write_real_array

end module

