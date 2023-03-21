!==========================================================
!
!------------------------------------------------------------------------------------------
! useful interface (read_bin_array) for reading arbitary binary arrays into an opened file
MODULE MOD_READ_BINARY_ARRAYS
use o_PARAM
private
public :: read_bin_array, read1d_int_static
INTERFACE read_bin_array
          MODULE PROCEDURE read1d_real, read1d_int, read1d_char, read2d_real, read2d_int, read3d_real, read3d_int
END INTERFACE
contains
subroutine read1d_real(arr, unit, iostat, iomsg)
    real(kind=WP), intent(inout), allocatable :: arr(:)
    integer,       intent(in)                 :: unit
    integer,       intent(out)                :: iostat
    character(*),  intent(inout)              :: iomsg
    integer                                   :: s1

    read(unit, iostat=iostat, iomsg=iomsg) s1
    if (s1==0) return
    if (.not. allocated(arr)) allocate(arr(s1))
    read(unit, iostat=iostat, iomsg=iomsg) arr(1:s1)
end subroutine read1d_real

subroutine read1d_int(arr, unit, iostat, iomsg)
    integer,       intent(inout), allocatable :: arr(:)
    integer,       intent(in)                 :: unit
    integer,       intent(out)                :: iostat
    character(*),  intent(inout)              :: iomsg
    integer                                   :: s1

    read(unit, iostat=iostat, iomsg=iomsg) s1
    if (s1==0) return
    if (.not. allocated(arr)) allocate(arr(s1))
    read(unit, iostat=iostat, iomsg=iomsg) arr(1:s1)
end subroutine read1d_int

subroutine read1d_char(arr, unit, iostat, iomsg)
    character,     intent(inout), allocatable :: arr(:)
    integer,       intent(in)                 :: unit
    integer,       intent(out)                :: iostat
    character(*),  intent(inout)              :: iomsg
    integer                                   :: s1

    read(unit, iostat=iostat, iomsg=iomsg) s1
    if (s1==0) return
    if (.not. allocated(arr)) allocate(arr(s1))
    read(unit, iostat=iostat, iomsg=iomsg) arr(1:s1)
end subroutine read1d_char

subroutine read1d_int_static(arr, unit, iostat, iomsg)
    IMPLICIT NONE
    integer,       intent(inout)              :: arr(:)
    integer,       intent(in)                 :: unit
    integer,       intent(out)                :: iostat
    character(*),  intent(inout)              :: iomsg
    integer                                   :: s1

    read(unit, iostat=iostat, iomsg=iomsg) s1
    if (s1==0) return
    read(unit, iostat=iostat, iomsg=iomsg) arr(1:s1)
end subroutine read1d_int_static

subroutine read2d_real(arr, unit, iostat, iomsg)
    real(kind=WP), intent(inout), allocatable :: arr(:,:)
    integer,       intent(in)                 :: unit
    integer,       intent(out)                :: iostat
    character(*),  intent(inout)              :: iomsg
    integer                                   :: s1, s2

    read(unit, iostat=iostat, iomsg=iomsg) s1, s2
    if ((s1==0) .or. (s2==0)) return
    if (.not. allocated(arr)) allocate(arr(s1, s2))
    read(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2)
end subroutine read2d_real

subroutine read2d_int(arr, unit, iostat, iomsg)
    integer, intent(inout), allocatable :: arr(:,:)
    integer,       intent(in)           :: unit
    integer,       intent(out)          :: iostat
    character(*),  intent(inout)        :: iomsg
    integer                             :: s1, s2

    read(unit, iostat=iostat, iomsg=iomsg) s1, s2
    if ((s1==0) .or. (s2==0)) return
    if (.not. allocated(arr)) allocate(arr(s1, s2))
    read(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2)
end subroutine read2d_int

subroutine read3d_real(arr, unit, iostat, iomsg)
    real(kind=WP), intent(inout), allocatable :: arr(:,:,:)
    integer,       intent(in)                 :: unit
    integer,       intent(out)                :: iostat
    character(*),  intent(inout)              :: iomsg
    integer                                   :: s1, s2, s3

    read(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3
    if ((s1==0) .or. (s2==0) .or. (s3==0)) return
    if (.not. allocated(arr)) allocate(arr(s1,s2,s3))
    read(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2, 1:s3)
end subroutine read3d_real

subroutine read3d_int(arr, unit, iostat, iomsg)
    integer, intent(inout), allocatable :: arr(:,:,:)
    integer,       intent(in)           :: unit
    integer,       intent(out)          :: iostat
    character(*),  intent(inout)        :: iomsg
    integer                             :: s1, s2, s3

    read(unit, iostat=iostat, iomsg=iomsg) s1, s2, s3
    if ((s1==0) .or. (s2==0) .or. (s3==0)) return
    if (.not. allocated(arr)) allocate(arr(s1,s2,s3))
    read(unit, iostat=iostat, iomsg=iomsg) arr(1:s1, 1:s2, 1:s3)
end subroutine read3d_int
end module MOD_READ_BINARY_ARRAYS
!==========================================================
