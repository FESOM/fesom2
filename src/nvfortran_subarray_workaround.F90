#ifdef ENABLE_NVHPC_WORKAROUNDS
module nvfortran_subarray_workaround_module
  use MOD_DYN
  implicit none

  type(t_dyn), pointer, save :: dynamics_workaround
end module nvfortran_subarray_workaround_module
#endif
