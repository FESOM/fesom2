module mod_hecuba
  implicit none
  interface
     subroutine hecuba_start_session(expname) &
                        bind(C, name="initHecubaSession")
              use iso_c_binding, only: c_char 
              character(kind=c_char) :: expname(*)
      end subroutine hecuba_start_session

      subroutine write_array(expname, varname, timestep, chunk, cdata, arr_size) &
                    bind(c,name="sendMetricsToHecuba")
              use iso_c_binding, only: c_char, c_int, c_ptr, c_float, c_double 
              character(kind=c_char) :: expname(*)
              character(kind=c_char) :: varname(*)
              integer(kind=c_int),value :: timestep
              integer, value :: chunk
              real(kind=c_double) :: cdata(*)
              integer(kind=c_int),value :: arr_size
      end subroutine write_array
      
      subroutine write_scalar(expname, varname, timestep, chunk, cdata, arr_size) &
                    bind(c,name="sendMetricsToHecuba")
              use iso_c_binding, only: c_char, c_int, c_ptr, c_float, c_double 
              character(kind=c_char) :: expname(*)
              character(kind=c_char) :: varname(*)
              integer(kind=c_int),value :: timestep
              integer, value :: chunk
              real(kind=c_double) :: cdata
              integer(kind=c_int),value :: arr_size
      end subroutine write_scalar
  
      subroutine hecuba_set_metadata(expname, key, value_str) &
                      bind(c,name="setExpMetaData")
              use iso_c_binding, only: c_char 
              character(kind=c_char) :: expname(*)
              character(kind=c_char) :: key(*)
              character(kind=c_char) :: value_str(*)
      end subroutine hecuba_set_metadata
      
      subroutine get_prune_state(expname, state) bind(C, name="getPruneState")
            use iso_c_binding, only: c_char, c_bool
            character(kind=c_char),intent(in) :: expname(*)
            logical(kind=c_bool),intent(out)   :: state
      end subroutine get_prune_state

  end interface

contains 
  subroutine write_time(expname, timestep, cdata) 
     use iso_c_binding, only: c_char, c_int, c_ptr, c_float, c_double, c_null_char
     character(kind=c_char) :: expname(*)
     integer(kind=c_int),value :: timestep
     real(kind=c_double) :: cdata
     
     call write_array(expname, "time"//c_null_char, timestep, -1, (/cdata/), 1)
     !call write_scalar(expname, varname, timestep, -1, cdata, arr_size)
  end subroutine write_time
 
  ! for coordinates that dont change in time
  ! coordinates that are chunked
  subroutine write_chunked_coordinate(expname, varname, timestep, chunk, cdata, arr_size) 
     use iso_c_binding, only: c_char, c_int, c_ptr, c_float, c_double 
     character(kind=c_char) :: expname(*) 
     character(kind=c_char) :: varname(*)
     integer(kind=c_int)    :: timestep
     integer(kind=c_int)    :: chunk
     real(kind=c_double)    :: cdata(*)
     integer(kind=c_int)    :: arr_size
     if (timestep<0) then 
        call write_array(expname, varname, -1, chunk, cdata, arr_size)
     else
        write(*,*) "warning skipping writing dimension to hecuba as both timestep and chunk are present use write_array instead" 
     end if
  end subroutine write_chunked_coordinate
  ! coordinates that are global not chunked
  subroutine write_global_coordinate(expname, varname, cdata, arr_size) 
     use iso_c_binding, only: c_char, c_int, c_double 
     character(kind=c_char) :: expname(*)
     character(kind=c_char) :: varname(*)
     real(kind=c_double) :: cdata(*)
     integer(kind=c_int),value :: arr_size
     
     call write_array(expname, varname, -1, -1, cdata, arr_size)
  end subroutine write_global_coordinate
end module
