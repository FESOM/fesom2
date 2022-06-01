module io_zarrdata

  use o_PARAM, only : WP
  use, intrinsic :: iso_fortran_env, only: real64, real32
  use io_data_strategy_module
  use async_threads_module
  implicit none
  public:: output_zarr
  private
  integer, save:: nsteps
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
  integer       :: n, k, i!, nsteps
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
  !dirnames = trim(root) // "/sst"  
  if (lfirst) then
     if (mype==0) then
      root = create_zarr_group(trim(filename))
     else
      root = filename 
     end if
     
     call get_run_steps(nsteps)
     ! not necessary
     !variable_path = create_zarr_group('sst',root) ! also adds dirs for each pe to hold its arrays
     ! testing
     !call write_zarr_array(istep, 'sst', reshape(tr_arr(1,1:myDim_nod2d,1),[size(tr_arr(1,1:myDim_nod2d,1))]), filename)
     lfirst=.false.
  end if
  ! reshape, inefficient but because we dont have interface for diffrernt dims
  call write_zarr_array(istep, 'sst','scalar', reshape(tr_arr(1,1:myDim_nod2d,1),[size(tr_arr(1,1:myDim_nod2d,1))]), filename)
  call write_zarr_array(istep, 'ssh','scalar', eta_n(1:myDim_nod2d), filename)
  

  call write_zarr_array(istep, 'temp','vector',reshape(tr_arr(:,1:myDim_nod2d,1),[size(tr_arr(:,1:myDim_nod2d,1))]), filename)
  call write_zarr_array(istep, 'salt','vector',reshape(tr_arr(:,1:myDim_nod2d,2),[size(tr_arr(:,1:myDim_nod2d,2))]), filename)
  call write_zarr_array(istep, 'unod','vector',reshape(Unode(1,:,1:myDim_nod2d),[size(Unode(1,:,1:myDim_nod2d))]), filename)
  call write_zarr_array(istep, 'vnod','vector',reshape(Unode(2,:,1:myDim_nod2d),[size(Unode(2,:,1:myDim_nod2d))]), filename)

  
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
  !call execute_command_line ('mkdir -p ' // adjustl(trim(group)))
  call mkdir(adjustl(trim(group))) ! fix for failing execute_command_line at times

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
subroutine write_zarr_array(istep, variable, vartype, rarray, root) !, dimensions, add_time_dim)
  use zarr_util 
  use mod_mesh
  use g_PARSUP
  use o_arrays
  use diagnostics

  implicit none
  character(len=*), intent(in)               :: variable, root
  integer, intent(in)           :: istep
  character(1000)               :: variable_group, group
  character(9999)               :: text
  character(6)                 :: vartype ! scalar(1d per time) or 3d
  real(kind=WP), intent(in) :: rarray(:) 
  character(2000)               ::  tmp
  integer          :: data_kind, fileunit !, nsteps
  character(1), parameter :: endianess="<" ! later find from func
  logical,save :: lfirst=.true.
  character(5), save :: mype_string
  character(10) :: dtype
  !logical                 :: add_time_dim=.true.
  !character(len=*), intent(in) :: dimensions(:)
   
  ! this redundant block will run on every call should be only done once for a variable, 
  ! for that, we need to track variable IO lists
  variable_group = create_zarr_group(variable, root) 
     
  write(mype_string,'(I0)') mype
  tmp=trim(variable_group) // "/" // trim(mype_string) 
  !should we inquire? before?
  call execute_command_line ('mkdir -p ' // adjustl(trim(tmp)))
  !call mkdir(adjustl(trim(tmp))) ! fix for failing execute_command_line at times

  data_kind = 8!sizeof(rarray)/size(rarray) !use storage_size intrensic?
  write(tmp, '(I0)') data_kind
  dtype = byteorder() // "f" // tmp
     
  !call get_run_steps(nsteps)

  !should we inquire? before or replace?
  tmp=trim(variable_group) // "/" // trim(mype_string) 
  tmp=trim(tmp) // "/.zarray"
  open(newunit=fileunit, file=trim(tmp)) 
  !call get_zarray_json(text, shape(rarray), shape(rarray), dtype) !dtype)
  !call get_zarray_json(text, [1,shape(rarray)], [nsteps, shape(rarray)], dtype)
  if (vartype == "scalar") then
     !call get_zarray_json(text, [shape(rarray),1], [shape(rarray),nsteps], dtype)
     call get_zarray_json(text, [myDim_nod2d,1], [myDim_nod2d,nsteps], dtype)
  else
     call get_zarray_json(text, [size(rarray)/myDim_nod2d,myDim_nod2d,1], [size(rarray)/myDim_nod2d,myDim_nod2d,nsteps], dtype)
  end if
             
  write(fileunit, '(a)') trim(text)
  flush(fileunit)
  close(fileunit)

  tmp=trim(variable_group) // "/" // trim(mype_string) 
  tmp=trim(tmp) // "/.zattrs"
  open(newunit=fileunit, file=trim(tmp))
  !text = '{"_ARRAY_DIMENSIONS":["nod2d"]}'
  !text = '{"_ARRAY_DIMENSIONS":["time","nod2d"]}' 
  if (vartype == "scalar") then
     text = '{"_ARRAY_DIMENSIONS":["nod2d", "time"]}' ! works for dask.array.concatinate
  else
     text = '{"_ARRAY_DIMENSIONS":["lev", "nod2d", "time"]}' ! works for dask.array.concatinate
  end if ! missing if nothing?

  write(fileunit, '(a)') trim(text)
  flush(fileunit)
  close(fileunit)
  ! end block redundant 
  write(tmp, '(I0)') istep-1  
  !tmp=trim(root) // '/' //trim(variable) // "/" // trim(mype_string) // '/0'
  !tmp=trim(root) // '/' //trim(variable) // "/" // trim(mype_string) // '/' // trim(tmp)// '.0'
  if (vartype == "scalar") then 
        tmp=trim(root) // '/' //trim(variable) // "/" // trim(mype_string) // '/0.' // trim(tmp)
  else
        tmp=trim(root) // '/' //trim(variable) // "/" // trim(mype_string) // '/0.0.' // trim(tmp)
  end if
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

! from fesom's upstream branch https://github.com/FESOM/fesom2/commit/64029aa9a1a2fef52eca5a34090d5baf52f71821 by @hegish
! using EXECUTE_COMMAND_LINE to call mkdir sometimes fail (EXECUTE_COMMAND_LINE is forked as an new process, which may be the problem)
! try to use the C mkdir as an alternative
subroutine mkdir(path)
    use iso_c_binding
    character(len=*), intent(in) :: path
    ! EO parameters
    integer result
    character(:,kind=C_CHAR), allocatable :: pathcopy

    interface
      function mkdir_c(path, mode) bind(c,name="mkdir")
        use iso_c_binding
        integer(c_int) mkdir_c
        character(kind=c_char,len=1) path(*)
        integer(c_int), value :: mode
      end function
    end interface

    pathcopy = path ! we need to pass an array of c_char to the C funcktion (this is not a correct type conversion, but Fortran characters seem to be of the same kind as c_char)
    ! result is 0 if the dir has been created from this call, otherwise -1
    ! the mode will not exactly be what we pass here, as it is subtracted by the umask bits (and possibly more)
    result = mkdir_c(pathcopy//C_NULL_CHAR, int(o'777', c_int))
  end subroutine mkdir

end module

