module forcing_provider_module
  implicit none
  public forcing_provider
  private
  
  type forcing_provider_type
    private
    contains
    procedure, public :: read_forcingdata
  end type
  type(forcing_provider_type), save :: forcing_provider ! handle to singleton instance
  
  
  contains
  subroutine read_forcingdata(this, filepath, varname, t_index1, forcingdata1)
    class(forcing_provider_type), intent(in) :: this
    character(len=*), intent(in) :: filepath
    character(len=*), intent(in) :: varname
    integer, intent(in) :: t_index1
    real(8), intent(out) :: forcingdata1(:,:)
    ! EO args
    real(4), allocatable :: values(:,:,:)
    
    call read_netcdf_timesteps(filepath, varname, t_index1, t_index1, values)
    
    forcingdata1 = values(:,:,1)
  end subroutine
  
  
  subroutine read_netcdf_timesteps(filepath, varname, timeindex_first, timeindex_last, values)
    character(len=*), intent(in) :: filepath
    character(len=*), intent(in) :: varname
    integer, intent(in) :: timeindex_first, timeindex_last
    real(4), allocatable, intent(out) :: values(:,:,:)
    ! EO args
    integer, parameter :: timedim_index = 3
    include "netcdf.inc" ! old netcdf fortran interface required?
    integer status
    integer fileid
    integer varid
    integer dim_size
    integer, allocatable, dimension(:) :: dimids
    integer, allocatable, dimension(:) :: dim_sizes
    integer i
    integer, allocatable, dimension(:) :: starts, sizes
   
    ! assume netcdf variable like: float q(time, lat, lon)
    status = nf_open(filepath,NF_NOWRITE,fileid)
    status = nf_inq_varid(fileid, varname, varid)
    status = nf_inq_varndims(fileid, varid, dim_size)
    allocate(dimids(dim_size))
    status = nf_inq_vardimid(fileid, varid, dimids)
    
    allocate(dim_sizes(dim_size))
    do i=1, dim_size
      status = nf_inq_dimlen(fileid, dimids(i), dim_sizes(i))
    end do
    
    ! todo: check if variable datatype is single precision (f77 real)
    
    allocate(starts(dim_size))
    allocate(sizes(dim_size))
    ! assert timedim_index == 3 && dim_size == 3
    
    ! todo: make this work if we have more than 3 dimensions and also if timedim_index != 3
    ! todo: check if values is already allocated
    starts(1) = 1
    starts(2) = 1
    starts(3) = timeindex_first
    sizes(1) = dim_sizes(1)
    sizes(2) = dim_sizes(2)
    sizes(3) = timeindex_last-timeindex_first+1
    allocate(values(sizes(1),sizes(2),sizes(3)))
    
    status = nf_get_vara_real(fileid, varid, starts, sizes, values)
    
    status = nf_close(fileid)
  end subroutine
  
end module
