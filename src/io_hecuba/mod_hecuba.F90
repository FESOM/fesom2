module mod_hecuba
implicit none
interface
        subroutine fhecuba_start_session() &
                        bind(C, name="start_hecuba_session")
        end subroutine fhecuba_start_session
        subroutine fload_datamodel() &
                        bind(C, name="load_datamodel")
        end subroutine
        subroutine putrealval(key, valueC, arr_size) &
                    bind(c,name="hecuba_put_array_val_C")
              use iso_c_binding, only: c_char, c_int, c_ptr, c_float, c_double 
              character(kind=c_char) :: key(*);
              real(kind=c_double) :: valueC(*);
              integer(kind=c_int) :: arr_size;

      end subroutine putrealval
      !below wecan import kind from param or
      !create an interface with multi subroutines dealing with types
      !ex:  INTERFACE write_array
      !      MODULE PROCEDURE put_int, real... defined in routines
      ! END INTERFACE
      !or using parametrized types

!      subroutine putrealval2(varname, ctime, chunk, valueC, arr_size) &
!                    bind(c,name="hecuba_put_array_val_C2")
!              use iso_c_binding, only: c_char, c_int, c_ptr, c_float, c_double 
!              character(kind=c_char) :: varname(*);
!              real(kind=8),value :: ctime;
!              integer, value :: chunk;
!              real(kind=c_double) :: valueC(*);
!              integer(kind=c_int) :: arr_size;
!
!      end subroutine putrealval2
      ! temp fix subroutines till options are explored that supports more types
      subroutine write_array(varname, ctime, chunk, valueC, arr_size) &
                    bind(c,name="hecuba_put_array_val_C2")
              use iso_c_binding, only: c_char, c_int, c_ptr, c_float, c_double 
              character(kind=c_char) :: varname(*);
              real(kind=8),value :: ctime;
              integer, value :: chunk;
              real(kind=c_double) :: valueC(*);
              integer(kind=c_int) :: arr_size;

      end subroutine write_array
!      interface write_array
!             module procedure write_array_r4, write_array_r8
!      end interface  
end interface

!below works but error that program cant see, make it public?
!interface write_array
!      module procedure write_array_r8
!end interface 

!contains 
!      subroutine write_array_r8(varname, ctime, chunk, valueC, arr_size) &
!                    bind(c,name="hecuba_put_array_val_C2")
!              use iso_c_binding, only: c_char, c_int, c_ptr, c_float, c_double 
!              character(kind=c_char) :: varname(*);
!              real(kind=8),value :: ctime;
!              integer, value :: chunk;
!              real(kind=c_double) :: valueC(*);
!              integer(kind=c_int) :: arr_size;
!
!      end subroutine write_array_r8
!
!      subroutine write_array_r4(varname, ctime, chunk, valueC, arr_size) &
!                    bind(c,name="hecuba_put_array_val_C2")
!              use iso_c_binding, only: c_char, c_int, c_ptr, c_float, c_double 
!              character(kind=c_char) :: varname(*);
!              real(kind=4),value :: ctime;
!              integer, value :: chunk;
!              real(kind=c_double) :: valueC(*);
!              integer(kind=c_int) :: arr_size;
!
!      end subroutine write_array_r4

end module