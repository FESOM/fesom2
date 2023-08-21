!==========================================================
MODULE MOD_DYN
USE O_PARAM
USE, intrinsic :: ISO_FORTRAN_ENV
USE MOD_WRITE_BINARY_ARRAYS
USE MOD_READ_BINARY_ARRAYS
IMPLICIT NONE
SAVE

!
!
!_______________________________________________________________________________
TYPE T_SOLVERINFO
    integer       :: ident   = 1
    integer       :: maxiter = 2000
    integer       :: restart = 15
    integer       :: fillin  = 3
    integer       :: lutype  = 2
    real(kind=WP) :: droptol = 1.e-8
!!! PARMS Solver
    real(kind=WP) :: soltol  = 1e-10  ! default for PARMS
    logical       :: use_parms = .TRUE.
!!!
!!! Sergey's Solver
!   real(kind=WP)  :: soltol  = 1e-5  ! default for PARMS
!   logical        :: use_parms = .FALSE.
!!!
    real(kind=WP), allocatable   :: rr(:), zz(:), pp(:), App(:)
    contains
    procedure WRITE_T_SOLVERINFO
    procedure READ_T_SOLVERINFO
END TYPE T_SOLVERINFO
!
!
!_______________________________________________________________________________
TYPE T_DYN_WORK
    real(kind=WP), allocatable, dimension(:,:,:) :: uvnode_rhs
    real(kind=WP), allocatable, dimension(:,:)   :: u_c, v_c

    ! easy backscatter contribution
    real(kind=WP), allocatable, dimension(:,:)   :: u_b, v_b
    contains
    procedure WRITE_T_DYN_WORK
    procedure READ_T_DYN_WORK
END TYPE T_DYN_WORK
!
!
!_______________________________________________________________________________
! set main structure for dynamicss, contains viscosity options and parameters +
! option for momentum advection
TYPE T_DYN
!___________________________________________________________________________
    ! instant zonal merdional velocity & Adams-Bashfort rhs
    real(kind=WP), allocatable, dimension(:,:,:):: uv, uv_rhs, uv_rhsAB, fer_uv

    ! horizontal velocities at nodes
    real(kind=WP), allocatable, dimension(:,:,:):: uvnode

    ! instant vertical vel arrays
    real(kind=WP), allocatable, dimension(:,:)  :: w, w_e, w_i, w_old, cfl_z, fer_w

    ! sea surface height arrays
    real(kind=WP), allocatable, dimension(:)    :: eta_n, d_eta, ssh_rhs, ssh_rhs_old
    
    ! LA: 2023-05-17 iceberg arrays
    real(kind=WP), allocatable, dimension(:)    :: eta_n_ib ! kh 18.03.21 additional array for asynchronous iceberg computations
    real(kind=WP), allocatable, dimension(:,:,:):: uv_ib    ! kh 18.03.21 additional array for asynchronous iceberg computations
    
    !___________________________________________________________________________
    ! summarizes solver input parameter
    type(t_solverinfo)                          :: solverinfo

    !___________________________________________________________________________
    ! put dynmiacs working arrays
    type(t_dyn_work)                            :: work

    !___________________________________________________________________________
    ! opt_visc=...
    ! 5=Kinematic (easy) Backscatter
    ! 6=Biharmonic flow aware (viscosity depends on velocity Laplacian)
    ! 7=Biharmonic flow aware (viscosity depends on velocity differences)
    ! 8=Dynamic Backscatter
    logical                                     :: check_opt_visc= .true.
    integer                                     :: opt_visc      = 5

    ! gamma0 [m/s],   backgroung viscosity= gamma0*len, it should be as small
    !                 as possible (keep it < 0.01 m/s).
    ! gamma1 [nodim], for computation of the flow aware viscosity
    ! gamma2 [s/m],   is only used in easy backscatter option
    real(kind=WP)                               :: visc_gamma0   = 0.03
    real(kind=WP)                               :: visc_gamma1   = 0.1
    real(kind=WP)                               :: visc_gamma2   = 0.285

    ! coefficient for returned sub-gridscale energy, to be used with opt_visc=5
    ! (easy backscatter)
    real(kind=WP)                               :: visc_easybsreturn = 1.5

    logical                                     :: use_ivertvisc = .true.
    integer                                     :: momadv_opt    = 2

    ! Switch on free slip
    logical                                     :: use_freeslip  = .false.

    ! do implicite, explicite spliting of vertical velocity
    logical                                     :: use_wsplit    = .false.
    ! maximum allowed CFL criteria in vertical (0.5 < w_max_cfl < 1.)
    ! in older FESOM it used to be w_exp_max=1.e-3
    real(kind=WP)                               :: wsplit_maxcfl= 1.0
    ! energy diagnostic part: will be computed inside the model ("hard integration"):
    logical                                      :: ldiag_ke       = .true.
    ! different contributions to velocity change. will be computed inside the code.
    real(kind=WP), allocatable, dimension(:,:,:) :: ke_adv, ke_cor, ke_pre, ke_hvis, ke_vvis, ke_du2, ke_umean, ke_u2mean
    real(kind=WP), allocatable, dimension(:,:)   :: ke_wind, ke_drag
    ! same as above but multiplied by velocity. we need both for later computation of turbulent fluxes
    real(kind=WP), allocatable, dimension(:,:,:) :: ke_adv_xVEL, ke_cor_xVEL, ke_pre_xVEL, ke_hvis_xVEL, ke_vvis_xVEL
    real(kind=WP), allocatable, dimension(:,:)   :: ke_wind_xVEL, ke_drag_xVEL
    real(kind=WP), allocatable, dimension(:,:)   :: ke_wrho         !we use pressure to compute (W*dens) as it appeares much easier to compute (P*dW) instead of (dP*w)
    real(kind=WP), allocatable, dimension(:,:)   :: ke_dW, ke_Pfull !for later computation of turbulent fluxes from the term above
    real(kind=WP), allocatable, dimension(:,:,:) :: ke_adv_AB, ke_cor_AB
    real(kind=WP), allocatable, dimension(:,:,:) :: ke_rhs_bak
    ! surface fields to compute APE generation
    real(kind=WP), allocatable, dimension(:)     :: ke_J, ke_D, ke_G, ke_D2, ke_n0, ke_JD, ke_GD, ke_swA, ke_swB

    !___________________________________________________________________________
    contains
#if defined(__PGI)
     procedure, private WRITE_T_DYN
     procedure, private READ_T_DYN
#else
     procedure WRITE_T_DYN
     procedure READ_T_DYN
#endif
     generic :: write(unformatted) => WRITE_T_DYN
     generic :: read(unformatted)  => READ_T_DYN
END TYPE T_DYN

contains

!
!
!_______________________________________________________________________________
! set unformatted writing and reading for T_DYN_WORK
subroutine WRITE_T_SOLVERINFO(tsolverinfo, unit)
    IMPLICIT NONE
    class(T_SOLVERINFO),  intent(in)     :: tsolverinfo
    integer,              intent(in)     :: unit
    integer                              :: iostat
    character(len=1024)                  :: iomsg
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%ident
    write(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%maxiter
    write(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%restart
    write(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%fillin
    write(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%lutype
    write(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%droptol
    write(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%soltol
    call write_bin_array(tsolverinfo%rr,  unit, iostat, iomsg)
    call write_bin_array(tsolverinfo%zz,  unit, iostat, iomsg)
    call write_bin_array(tsolverinfo%pp,  unit, iostat, iomsg)
    call write_bin_array(tsolverinfo%App, unit, iostat, iomsg)
end subroutine WRITE_T_SOLVERINFO

subroutine READ_T_SOLVERINFO(tsolverinfo, unit)
    IMPLICIT NONE
    class(T_SOLVERINFO),  intent(inout)     :: tsolverinfo
    integer,              intent(in)     :: unit
    integer                              :: iostat
    character(len=1024)                  :: iomsg
    read(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%ident
    read(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%maxiter
    read(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%restart
    read(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%fillin
    read(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%lutype
    read(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%droptol
    read(unit, iostat=iostat, iomsg=iomsg) tsolverinfo%soltol
    call read_bin_array(tsolverinfo%rr,  unit, iostat, iomsg)
    call read_bin_array(tsolverinfo%zz,  unit, iostat, iomsg)
    call read_bin_array(tsolverinfo%pp,  unit, iostat, iomsg)
    call read_bin_array(tsolverinfo%App, unit, iostat, iomsg)
end subroutine READ_T_SOLVERINFO

!
!
!_______________________________________________________________________________
! set unformatted writing and reading for T_DYN_WORK
subroutine WRITE_T_DYN_WORK(twork, unit)
    IMPLICIT NONE
    class(T_DYN_WORK), intent(in)        :: twork
    integer,              intent(in)     :: unit
    integer                              :: iostat
    character(len=1024)                  :: iomsg
    call write_bin_array(twork%uvnode_rhs, unit, iostat, iomsg)
    call write_bin_array(twork%u_c,        unit, iostat, iomsg)
    call write_bin_array(twork%v_c,        unit, iostat, iomsg)
    call write_bin_array(twork%u_b,        unit, iostat, iomsg)
    call write_bin_array(twork%v_b,        unit, iostat, iomsg)
end subroutine WRITE_T_DYN_WORK

subroutine READ_T_DYN_WORK(twork, unit)
    IMPLICIT NONE
    class(T_DYN_WORK), intent(inout)        :: twork
    integer,              intent(in)     :: unit
    integer                              :: iostat
    character(len=1024)                  :: iomsg
    call read_bin_array(twork%uvnode_rhs, unit, iostat, iomsg)
    call read_bin_array(twork%u_c,        unit, iostat, iomsg)
    call read_bin_array(twork%v_c,        unit, iostat, iomsg)
    call read_bin_array(twork%u_b,        unit, iostat, iomsg)
    call read_bin_array(twork%v_b,        unit, iostat, iomsg)
end subroutine READ_T_DYN_WORK

!
!
!_______________________________________________________________________________
! set unformatted writing and reading for T_DYN
subroutine WRITE_T_DYN(dynamics, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN),         intent(in)     :: dynamics
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg

    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%opt_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_gamma0
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_gamma1
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_gamma2
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_easybsreturn

    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%use_ivertvisc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%momadv_opt

    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%use_freeslip
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%use_wsplit
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%wsplit_maxcfl

    !___________________________________________________________________________
    call dynamics%solverinfo%WRITE_T_SOLVERINFO(unit)

    !___________________________________________________________________________
    call dynamics%work%WRITE_T_DYN_WORK(unit)

    !___________________________________________________________________________
    call write_bin_array(dynamics%uv        , unit, iostat, iomsg)
    call write_bin_array(dynamics%uv_rhs    , unit, iostat, iomsg)
    call write_bin_array(dynamics%uv_rhsAB  , unit, iostat, iomsg)
    call write_bin_array(dynamics%uvnode    , unit, iostat, iomsg)
    call write_bin_array(dynamics%w         , unit, iostat, iomsg)
    call write_bin_array(dynamics%w_e       , unit, iostat, iomsg)
    call write_bin_array(dynamics%w_i       , unit, iostat, iomsg)
    call write_bin_array(dynamics%cfl_z     , unit, iostat, iomsg)
    if (Fer_GM) then
        call write_bin_array(dynamics%fer_w , unit, iostat, iomsg)
        call write_bin_array(dynamics%fer_uv, unit, iostat, iomsg)
    end if


end subroutine WRITE_T_DYN

subroutine READ_T_DYN(dynamics, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN),         intent(inout)  :: dynamics
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg

    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%opt_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_gamma0
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_gamma1
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_gamma2
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_easybsreturn

    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%use_ivertvisc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%momadv_opt

    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%use_freeslip
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%use_wsplit
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%wsplit_maxcfl

    !___________________________________________________________________________
    call dynamics%solverinfo%READ_T_SOLVERINFO(unit)

    !___________________________________________________________________________
    call dynamics%work%READ_T_DYN_WORK(unit)

    !___________________________________________________________________________
    call read_bin_array(dynamics%uv        , unit, iostat, iomsg)
    call read_bin_array(dynamics%uv_rhs    , unit, iostat, iomsg)
    call read_bin_array(dynamics%uv_rhsAB  , unit, iostat, iomsg)
    call read_bin_array(dynamics%uvnode    , unit, iostat, iomsg)
    call read_bin_array(dynamics%w         , unit, iostat, iomsg)
    call read_bin_array(dynamics%w_e       , unit, iostat, iomsg)
    call read_bin_array(dynamics%w_i       , unit, iostat, iomsg)
    call read_bin_array(dynamics%cfl_z     , unit, iostat, iomsg)
    if (Fer_GM) then
        call read_bin_array(dynamics%fer_w     , unit, iostat, iomsg)
        call read_bin_array(dynamics%fer_uv    , unit, iostat, iomsg)
    end if

end subroutine READ_T_DYN

END MODULE MOD_DYN
