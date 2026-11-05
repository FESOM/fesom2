#ifdef ENABLE_NVHPC_WORKAROUNDS
module nvfortran_subarray_workaround_module
  ! nvfortran passes a section of a derived-type component (e.g. dynamics%uv(1,:,:)) to the def_*_var
  ! routines as a temporary copy, so the stored pointer dangles: the restart I/O re-points it by name.
  use o_PARAM, only: WP
  use MOD_DYN
  use MOD_TRACER
  implicit none

  type(t_dyn),    pointer, save :: dynamics_workaround => null()
  type(t_tracer), pointer, save :: tracers_workaround  => null()
  ! restart name of every tracer (temp, salt, ...), set in ini_ocean_io
  character(len=500), allocatable, save :: tracer_restart_name(:)

contains

  ! the array a restart variable really refers to: the section of the derived-type
  ! component for the variables registered as such sections, default_ptr otherwise
  function nvhpc_restart_data(varname, default_ptr) result(p)
    character(len=*), intent(in)          :: varname
    real(kind=WP), pointer, intent(in)    :: default_ptr(:,:)
    real(kind=WP), pointer                :: p(:,:)
    integer                               :: j

    p => default_ptr
    if (associated(dynamics_workaround)) then
      select case (varname)
      case ('u');        p => dynamics_workaround%uv(1,:,:)
      case ('v');        p => dynamics_workaround%uv(2,:,:)
      case ('urhs_AB');  p => dynamics_workaround%uv_rhsAB(1,1,:,:)
      case ('vrhs_AB');  p => dynamics_workaround%uv_rhsAB(1,2,:,:)
      case ('urhs_AB3'); p => dynamics_workaround%uv_rhsAB(2,1,:,:)
      case ('vrhs_AB3'); p => dynamics_workaround%uv_rhsAB(2,2,:,:)
      end select
    end if
    if (associated(tracers_workaround) .and. allocated(tracer_restart_name)) then
      do j=1, size(tracer_restart_name)
        if (varname == trim(tracer_restart_name(j))//'_M1') p => tracers_workaround%data(j)%valuesold(1,:,:)
        if (varname == trim(tracer_restart_name(j))//'_M2') p => tracers_workaround%data(j)%valuesold(2,:,:)
      end do
    end if
  end function nvhpc_restart_data
end module nvfortran_subarray_workaround_module
#endif
