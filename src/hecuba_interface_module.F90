module hecuba_interface_module
  implicit none

contains

    subroutine hecuba_init_datamodel(thread_idx)
      use iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: thread_idx
      ! EO args
      call init_datamodel_c(thread_idx)
    end subroutine
    
    subroutine hecuba_output(ts, name, description, units, freq, freq_unit, accuracy )
      use iso_c_binding
      use g_clock
      use g_PARSUP
      use mod_mesh
      use fesom_version_info_module
      use g_config
      use i_PARAM
      use o_PARAM
    
      implicit none
      ! variable declarations 
      integer(c_int), intent(in), value :: ts
      character(kind=c_char), dimension(100), intent(in) :: name
      character(kind=c_char), dimension(500), intent(in) :: description
      character(kind=c_char), dimension(100), intent(in) :: units
      integer(c_int), intent(in), value :: freq
      character(kind=c_char), dimension(*), intent(in) :: freq_unit
      integer(c_int), intent(in), value :: accuracy
      type(t_mesh) mesh
      
      !TODO: check that name is coming the wrong way with spaces causing error on creating table
      call hecuba_output_c(ts, name, description, units, freq, freq_unit, accuracy )
      
    end subroutine
    

end module
