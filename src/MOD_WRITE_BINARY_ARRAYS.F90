!==========================================================
!
!------------------------------------------------------------------------------------------
! useful interface (write_bin_array) for writing arbitary binary arrays into an opened file
MODULE MOD_WRITE_BINARY_ARRAYS
use o_PARAM
private
public :: write_bin_array, write1d_int_static
INTERFACE write_bin_array
          MODULE PROCEDURE write1d_real, write1d_int, write1d_char, write2d_real, write2d_int, write3d_real, write3d_int, write4d_real, write4d_int
END INTERFACE
contains

subroutine write1d_real(arr, unit, iostat, iomsg)
    real(kind=WP), intent(in), allocatable  :: arr(:)
    integer,       intent(in)               :: unit
    integer,       intent(out)              :: iostat
    character(*),  intent(inout)            :: iomsg
    integer                                 :: s1

    if (allocated(arr)) then
       s1=size(arr, 1)
       write(unit, iostat=iostat, iomsg=iomsg) s1
       write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1)
    else
       s1=0
       write(unit, iostat=iostat, iomsg=iomsg) s1
    end if
end subroutine write1d_real

subroutine write1d_int(arr, unit, iostat, iomsg)
    integer,       intent(in), allocatable  :: arr(:)
    integer,       intent(in)               :: unit
    integer,       intent(out)              :: iostat
    character(*),  intent(inout)            :: iomsg
    integer                                 :: s1

    if (allocated(arr)) then
       s1=size(arr, 1)
       write(unit, iostat=iostat, iomsg=iomsg) s1
       write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1)
    else
       s1=0
       write(unit, iostat=iostat, iomsg=iomsg) s1
    end if
end subroutine write1d_int

subroutine write1d_char(arr, unit, iostat, iomsg)
    character,     intent(in), allocatable  :: arr(:)
    integer,       intent(in)               :: unit
    integer,       intent(out)              :: iostat
    character(*),  intent(inout)            :: iomsg
    integer                                 :: s1

    if (allocated(arr)) then
       s1=size(arr, 1)
       write(unit, iostat=iostat, iomsg=iomsg) s1
       write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1)
    else
       s1=0
       write(unit, iostat=iostat, iomsg=iomsg) s1
    end if
end subroutine write1d_char

subroutine write1d_int_static(arr, unit, iostat, iomsg)
    IMPLICIT NONE
    integer,       intent(in)               :: arr(:)
    integer,       intent(in)               :: unit
    integer,       intent(out)              :: iostat
    character(*),  intent(inout)            :: iomsg
    integer                                 :: s1

    s1=size(arr, 1)
    write(unit, iostat=iostat, iomsg=iomsg) s1
    write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1)
end subroutine write1d_int_static

subroutine write2d_real(arr, unit, iostat, iomsg)
    real(kind=WP), intent(in), allocatable :: arr(:,:)
    integer,       intent(in)              :: unit
    integer,       intent(out)             :: iostat
    character(*),  intent(inout)           :: iomsg
    integer                                :: s1, s2

    if (allocated(arr)) then
       s1=size(arr, 1)
       s2=size(arr, 2)
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2
       write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2)
    else
       s1=0
       s2=0
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2
    end if
end subroutine write2d_real

subroutine write2d_int(arr, unit, iostat, iomsg)
    integer,       intent(in), allocatable :: arr(:,:)
    integer,       intent(in)              :: unit
    integer,       intent(out)             :: iostat
    character(*),  intent(inout)           :: iomsg
    integer                                :: s1, s2

    if (allocated(arr)) then
       s1=size(arr, 1)
       s2=size(arr, 2)
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2
       write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2)
    else
       s1=0
       s2=0
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2
    end if
end subroutine write2d_int


subroutine write3d_real(arr, unit, iostat, iomsg)
    real(kind=WP), intent(in), allocatable :: arr(:,:,:)
    integer,       intent(in)              :: unit
    integer,       intent(out)             :: iostat
    character(*),  intent(inout)           :: iomsg
    integer                                :: s1, s2, s3

    if (allocated(arr)) then
       s1=size(arr, 1)
       s2=size(arr, 2)
       s3=size(arr, 3)
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3
       write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2, 1:s3)
    else
       s1=0
       s2=0
       s3=0
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3
    end if
end subroutine write3d_real

subroutine write3d_int(arr, unit, iostat, iomsg)
    integer,       intent(in), allocatable :: arr(:,:,:)
    integer,       intent(in)              :: unit
    integer,       intent(out)             :: iostat
    character(*),  intent(inout)           :: iomsg
    integer                                :: s1, s2, s3

    if (allocated(arr)) then
       s1=size(arr, 1)
       s2=size(arr, 2)
       s3=size(arr, 3)
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3
       write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2, 1:s3)
    else
       s1=0
       s2=0
       s3=0
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3
    end if
end subroutine write3d_int
subroutine write4d_real(arr, unit, iostat, iomsg)
    real(kind=WP), intent(in), allocatable :: arr(:,:,:,:)
    integer,       intent(in)              :: unit
    integer,       intent(out)             :: iostat
    character(*),  intent(inout)           :: iomsg
    integer                                :: s1, s2, s3, s4

    if (allocated(arr)) then
       s1=size(arr, 1)
       s2=size(arr, 2)
       s3=size(arr, 3)
       s4=size(arr, 4)
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3, s4
       write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2, 1:s3, 1:s4)
    else
       s1=0
       s2=0
       s3=0
       s4=0
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3, s4
    end if
end subroutine write4d_real

subroutine write4d_int(arr, unit, iostat, iomsg)
    integer,       intent(in), allocatable :: arr(:,:,:,:)
    integer,       intent(in)              :: unit
    integer,       intent(out)             :: iostat
    character(*),  intent(inout)           :: iomsg
    integer                                :: s1, s2, s3, s4

    if (allocated(arr)) then
       s1=size(arr, 1)
       s2=size(arr, 2)
       s3=size(arr, 3)
       s4=size(arr, 4)
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3, s4
       write(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2, 1:s3, 1:s4)
    else
       s1=0
       s2=0
       s3=0
       s4=0
       write(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3, s4
    end if
end subroutine write4d_int

end module MOD_WRITE_BINARY_ARRAYS
!==========================================================
