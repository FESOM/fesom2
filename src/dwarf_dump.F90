! Dump of the inputs and outputs of the horizontal tracer advection of one call, per rank, for the
! standalone dwarf in dwarf/tra_adv_hor. Enabled by the environment variable FESOM_DWARF_DUMP=<k>:
! the k-th call of do_oce_adv_tra for tracer 1 is dumped to FESOM_DWARF_DIR (default dwarf_dump)
! as tra_adv_<rank>.bin. Each record: name (32 characters), type (1 integer, 2 real64), rank,
! shape, values; all as unformatted stream.
module dwarf_dump_module
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    private
    public :: dwarf_dump_active, dwarf_dump_open, dwarf_dump_close, dput

    integer, save :: unit = -1, ncall = 0, target_call = -1
    logical, save :: checked = .false.

    interface dput
        module procedure put_i1, put_i2, put_r0, put_r1, put_r2, put_r3
    end interface

contains

    ! true for the call to dump; counts calls (tracer 1 only)
    logical function dwarf_dump_active(tr_num)
        integer, intent(in) :: tr_num
        character(len=32) :: v
        integer :: l, st
        if (.not. checked) then
            checked = .true.
            call get_environment_variable('FESOM_DWARF_DUMP', v, l, st)
            if (st == 0 .and. l > 0) read(v, *) target_call
        end if
        dwarf_dump_active = .false.
        if (target_call < 0 .or. tr_num /= 1) return
        ncall = ncall + 1
        dwarf_dump_active = (ncall == target_call)
    end function dwarf_dump_active

    subroutine dwarf_dump_open(mype)
        integer, intent(in) :: mype
        character(len=512) :: dir
        character(len=600) :: fname
        integer :: l, st
        call get_environment_variable('FESOM_DWARF_DIR', dir, l, st)
        if (st /= 0 .or. l == 0) dir = 'dwarf_dump'
        call execute_command_line('mkdir -p '//trim(dir))
        write(fname, '(a,a,i5.5,a)') trim(dir), '/tra_adv_', mype, '.bin'
        open(newunit=unit, file=trim(fname), access='stream', form='unformatted', status='replace')
    end subroutine dwarf_dump_open

    subroutine dwarf_dump_close()
        if (unit > 0) close(unit)
        unit = -1
    end subroutine dwarf_dump_close

    subroutine head(name, typ, shp)
        character(len=*), intent(in) :: name
        integer,          intent(in) :: typ, shp(:)
        character(len=32) :: nm
        nm = name
        write(unit) nm, int(typ, int32), int(size(shp), int32), int(shp, int32)
    end subroutine head

    subroutine put_i1(name, a)
        character(len=*), intent(in) :: name
        integer, intent(in) :: a(:)
        call head(name, 1, shape(a)); write(unit) int(a, int32)
    end subroutine put_i1
    subroutine put_i2(name, a)
        character(len=*), intent(in) :: name
        integer, intent(in) :: a(:,:)
        call head(name, 1, shape(a)); write(unit) int(a, int32)
    end subroutine put_i2
    subroutine put_r0(name, a)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: a
        integer :: s(0)
        call head(name, 2, s); write(unit) a
    end subroutine put_r0
    subroutine put_r1(name, a)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: a(:)
        call head(name, 2, shape(a)); write(unit) a
    end subroutine put_r1
    subroutine put_r2(name, a)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: a(:,:)
        call head(name, 2, shape(a)); write(unit) a
    end subroutine put_r2
    subroutine put_r3(name, a)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: a(:,:,:)
        call head(name, 2, shape(a)); write(unit) a
    end subroutine put_r3
end module dwarf_dump_module
