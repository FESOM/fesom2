!==========================================================
MODULE MOD_TRACER
USE O_PARAM
USE, intrinsic :: ISO_FORTRAN_ENV
USE MOD_WRITE_BINARY_ARRAYS
USE MOD_READ_BINARY_ARRAYS
IMPLICIT NONE
SAVE

TYPE T_TRACER_DATA
real(kind=WP), allocatable, dimension(:,:)  :: values, valuesAB  ! instant values & Adams-Bashfort interpolation
logical                                     :: smooth_bh_tra=.false.
real(kind=WP)                               :: gamma0_tra, gamma1_tra, gamma2_tra
logical                                     :: i_vert_diff =.false.
character(20)                               :: tra_adv_hor, tra_adv_ver, tra_adv_lim ! type of the advection scheme for this tracer
real(kind=WP)                               :: tra_adv_ph  = 1.  ! a parameter to be used in horizontal advection (for MUSCL it is the fraction of fourth-order contribution in the solution)
real(kind=WP)                               :: tra_adv_pv  = 1.  ! a parameter to be used in horizontal advection (for QR4C  it is the fraction of fourth-order contribution in the solution)
integer                                     :: ID

contains
  procedure WRITE_T_TRACER_DATA
  procedure READ_T_TRACER_DATA
  generic :: write(unformatted) => WRITE_T_TRACER_DATA
  generic :: read(unformatted)  => READ_T_TRACER_DATA
END TYPE T_TRACER_DATA


TYPE T_TRACER_WORK
!auxuary arrays to work with tracers:
real(kind=WP), allocatable         :: del_ttf(:,:)
real(kind=WP), allocatable         :: del_ttf_advhoriz(:,:),del_ttf_advvert(:,:)
!_______________________________________________________________________________
! in case ldiag_DVD=.true. --> calculate discrete variance decay (DVD)
real(kind=WP), allocatable                    :: tr_dvd_horiz(:,:,:), tr_dvd_vert(:,:,:)
! The fct part
real(kind=WP),allocatable,dimension(:,:)      :: fct_LO          ! Low-order solution
real(kind=WP),allocatable,dimension(:,:)      :: adv_flux_hor    ! Antidif. horiz. contrib. from edges / backup for iterafive fct scheme
real(kind=WP),allocatable,dimension(:,:)      :: adv_flux_ver    ! Antidif. vert. fluxes from nodes    / backup for iterafive fct scheme

real(kind=WP),allocatable,dimension(:,:)      :: fct_ttf_max,fct_ttf_min
real(kind=WP),allocatable,dimension(:,:)      :: fct_plus,fct_minus
! MUSCL type reconstruction
integer,allocatable,dimension(:)              :: nboundary_lay
integer,allocatable,dimension(:,:)            :: edge_up_dn_tri
real(kind=WP),allocatable,dimension(:,:,:)    :: edge_up_dn_grad

contains
  procedure WRITE_T_TRACER_WORK
  procedure READ_T_TRACER_WORK
  generic :: write(unformatted) => WRITE_T_TRACER_WORK
  generic :: read(unformatted)  => READ_T_TRACER_WORK
END TYPE T_TRACER_WORK

! auxury type for reading namelist.tra
TYPE NML_TRACER_LIST_TYPE
        INTEGER                 :: ID         =-1
        CHARACTER(len=4)        :: adv_hor    ='NONE'
        CHARACTER(len=4)        :: adv_ver    ='NONE'
        CHARACTER(len=4)        :: adv_lim    ='NONE'
        REAL(kind=WP)           :: adv_ph     =1.
        REAL(kind=WP)           :: adv_pv     =1.
END TYPE NML_TRACER_LIST_TYPE

TYPE T_TRACER
! total number of tracers:
integer                                     :: num_tracers=2
type(t_tracer_data), allocatable            :: data(:)
type(t_tracer_work)                         :: work
! general options for all tracers (can be moved to T_TRACER is needed)
! bharmonic diffusion for tracers. We recommend to use this option in very high resolution runs (Redi is generally off there).
logical                       :: smooth_bh_tra = .false.
real(kind=WP)                 :: gamma0_tra    = 0.0005
real(kind=WP)                 :: gamma1_tra    = 0.0125
real(kind=WP)                 :: gamma2_tra    = 0.
logical                       :: i_vert_diff   = .true.

contains
procedure WRITE_T_TRACER
procedure READ_T_TRACER
generic :: write(unformatted) => WRITE_T_TRACER
generic :: read(unformatted)  => READ_T_TRACER
END TYPE T_TRACER

contains

! Unformatted writing for T_TRACER_DATA
subroutine WRITE_T_TRACER_DATA(tdata, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_TRACER_DATA), intent(in)     :: tdata
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg

    call write_bin_array(tdata%values,   unit, iostat, iomsg)
    call write_bin_array(tdata%valuesAB, unit, iostat, iomsg)
    write(unit, iostat=iostat, iomsg=iomsg) tdata%smooth_bh_tra
    write(unit, iostat=iostat, iomsg=iomsg) tdata%gamma0_tra
    write(unit, iostat=iostat, iomsg=iomsg) tdata%gamma1_tra
    write(unit, iostat=iostat, iomsg=iomsg) tdata%gamma2_tra
    write(unit, iostat=iostat, iomsg=iomsg) tdata%i_vert_diff
    write(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_hor
    write(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_ver
    write(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_lim
    write(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_ph
    write(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_pv
    write(unit, iostat=iostat, iomsg=iomsg) tdata%ID
end subroutine WRITE_T_TRACER_DATA

! Unformatted reading for T_TRACER_DATA
subroutine READ_T_TRACER_DATA(tdata, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_TRACER_DATA), intent(inout)  :: tdata
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg

    call read_bin_array(tdata%values,   unit, iostat, iomsg)
    call read_bin_array(tdata%valuesAB, unit, iostat, iomsg)
    read(unit, iostat=iostat, iomsg=iomsg) tdata%smooth_bh_tra
    read(unit, iostat=iostat, iomsg=iomsg) tdata%gamma0_tra
    read(unit, iostat=iostat, iomsg=iomsg) tdata%gamma1_tra
    read(unit, iostat=iostat, iomsg=iomsg) tdata%gamma2_tra
    read(unit, iostat=iostat, iomsg=iomsg) tdata%i_vert_diff
    read(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_hor
    read(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_ver
    read(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_lim
    read(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_ph
    read(unit, iostat=iostat, iomsg=iomsg) tdata%tra_adv_pv
    read(unit, iostat=iostat, iomsg=iomsg) tdata%ID
end subroutine READ_T_TRACER_DATA

! Unformatted writing for T_TRACER_WORK
subroutine WRITE_T_TRACER_WORK(twork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_TRACER_WORK), intent(in)     :: twork
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg

    call write_bin_array(twork%del_ttf,          unit, iostat, iomsg)
    call write_bin_array(twork%del_ttf_advhoriz, unit, iostat, iomsg)
    call write_bin_array(twork%del_ttf_advvert,  unit, iostat, iomsg)
    call write_bin_array(twork%tr_dvd_horiz,     unit, iostat, iomsg)
    call write_bin_array(twork%tr_dvd_vert,      unit, iostat, iomsg)
    call write_bin_array(twork%fct_LO,           unit, iostat, iomsg)
    call write_bin_array(twork%adv_flux_hor,     unit, iostat, iomsg)
    call write_bin_array(twork%adv_flux_ver,     unit, iostat, iomsg)
    call write_bin_array(twork%fct_ttf_max,      unit, iostat, iomsg)
    call write_bin_array(twork%fct_ttf_min,      unit, iostat, iomsg)
    call write_bin_array(twork%fct_plus,         unit, iostat, iomsg)
    call write_bin_array(twork%fct_minus,        unit, iostat, iomsg)
    call write_bin_array(twork%nboundary_lay,    unit, iostat, iomsg)
    call write_bin_array(twork%edge_up_dn_tri,   unit, iostat, iomsg)
    call write_bin_array(twork%edge_up_dn_grad,  unit, iostat, iomsg)
end subroutine WRITE_T_TRACER_WORK

! Unformatted reading for T_TRACER_WORK
subroutine READ_T_TRACER_WORK(twork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_TRACER_WORK), intent(inout)  :: twork
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg

    call read_bin_array(twork%del_ttf,          unit, iostat, iomsg)
    call read_bin_array(twork%del_ttf_advhoriz, unit, iostat, iomsg)
    call read_bin_array(twork%del_ttf_advvert,  unit, iostat, iomsg)
    call read_bin_array(twork%tr_dvd_horiz,     unit, iostat, iomsg)
    call read_bin_array(twork%tr_dvd_vert,      unit, iostat, iomsg)
    call read_bin_array(twork%fct_LO,           unit, iostat, iomsg)
    call read_bin_array(twork%adv_flux_hor,     unit, iostat, iomsg)
    call read_bin_array(twork%adv_flux_ver,     unit, iostat, iomsg)
    call read_bin_array(twork%fct_ttf_max,      unit, iostat, iomsg)
    call read_bin_array(twork%fct_ttf_min,      unit, iostat, iomsg)
    call read_bin_array(twork%fct_plus,         unit, iostat, iomsg)
    call read_bin_array(twork%fct_minus,        unit, iostat, iomsg)
    call read_bin_array(twork%nboundary_lay,    unit, iostat, iomsg)
    call read_bin_array(twork%edge_up_dn_tri,   unit, iostat, iomsg)
    call read_bin_array(twork%edge_up_dn_grad,  unit, iostat, iomsg)
end subroutine READ_T_TRACER_WORK

! Unformatted writing for T_TRACER
subroutine WRITE_T_TRACER(tracer, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_TRACER),      intent(in)     :: tracer
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    integer                              :: i

    write(unit, iostat=iostat, iomsg=iomsg) tracer%num_tracers
    do i=1, tracer%num_tracers
       write(unit, iostat=iostat, iomsg=iomsg) tracer%data(i)
    end do
    write(unit, iostat=iostat, iomsg=iomsg)    tracer%work
    write(unit, iostat=iostat, iomsg=iomsg)    tracer%smooth_bh_tra
    write(unit, iostat=iostat, iomsg=iomsg)    tracer%gamma0_tra
    write(unit, iostat=iostat, iomsg=iomsg)    tracer%gamma1_tra
    write(unit, iostat=iostat, iomsg=iomsg)    tracer%gamma2_tra
    write(unit, iostat=iostat, iomsg=iomsg)    tracer%i_vert_diff
end subroutine WRITE_T_TRACER

! Unformatted reading for T_TRACER
subroutine READ_T_TRACER(tracer, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_TRACER),      intent(inout)  :: tracer
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    integer                              :: i

    read(unit, iostat=iostat, iomsg=iomsg)    tracer%num_tracers
!   write(*,*) 'number of tracers to read: ', tracer%num_tracers
    allocate(tracer%data(tracer%num_tracers))
    do i=1, tracer%num_tracers
       read(unit, iostat=iostat, iomsg=iomsg) tracer%data(i)
!      write(*,*) 'tracer info:', tracer%data(i)%ID, TRIM(tracer%data(i)%tra_adv_hor), TRIM(tracer%data(i)%tra_adv_ver), TRIM(tracer%data(i)%tra_adv_lim)
    end do
    read(unit, iostat=iostat, iomsg=iomsg)    tracer%work
    read(unit, iostat=iostat, iomsg=iomsg)    tracer%smooth_bh_tra
    read(unit, iostat=iostat, iomsg=iomsg)    tracer%gamma0_tra
    read(unit, iostat=iostat, iomsg=iomsg)    tracer%gamma1_tra
    read(unit, iostat=iostat, iomsg=iomsg)    tracer%gamma2_tra
    read(unit, iostat=iostat, iomsg=iomsg)    tracer%i_vert_diff
end subroutine READ_T_TRACER
end module MOD_TRACER
!==========================================================

