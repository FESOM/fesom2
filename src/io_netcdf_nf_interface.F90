module io_netcdf_nf_interface
implicit none

  interface
    function nf_get_vara_double(ncid, varid, start, counts, dvals) result(status)
      integer, intent(in) :: ncid, varid
      integer, intent(in) :: start(*), counts(*)
      real(8), intent(out) :: dvals(*)
      integer status
    end function nf_get_vara_double

    function nf_get_vara_real(ncid, varid, start, counts, dvals) result(status)
      integer, intent(in) :: ncid, varid
      integer, intent(in) :: start(*), counts(*)
      real(4), intent(out) :: dvals(*)
      integer status
    end function nf_get_vara_real
  end interface


  interface nf_get_vara_x
    procedure nf_get_vara_real, nf_get_vara_double
  end interface


  interface
    function nf_put_vara_double(ncid, varid, start, counts, dvals) result(status)
      integer, intent(in) :: ncid, varid
      integer, intent(in) :: start(*), counts(*)
      real(8), intent(in) :: dvals(*)
      integer status
    end function nf_put_vara_double

    function nf_put_vara_real(ncid, varid, start, counts, dvals) result(status)
      integer, intent(in) :: ncid, varid
      integer, intent(in) :: start(*), counts(*)
      real(4), intent(in) :: dvals(*)
      integer status
    end function nf_put_vara_real
  end interface


  interface nf_put_vara_x
    procedure nf_put_vara_real, nf_put_vara_double
  end interface

end module io_netcdf_nf_interface
