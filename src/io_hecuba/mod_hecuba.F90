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
end interface

end module
