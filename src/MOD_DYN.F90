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
        generic :: write(unformatted) => WRITE_T_SOLVERINFO
        generic :: read(unformatted)  => READ_T_SOLVERINFO
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
        generic :: write(unformatted) => WRITE_T_DYN_WORK
        generic :: read(unformatted)  => READ_T_DYN_WORK
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
    real(kind=WP), allocatable, dimension(:,:)  :: w, w_e, w_i, cfl_z, fer_w
    
    ! sea surface height arrays
    real(kind=WP), allocatable, dimension(:)    :: eta_n, d_eta, ssh_rhs, ssh_rhs_old
    
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

    !___________________________________________________________________________
    contains
        procedure WRITE_T_DYN
        procedure READ_T_DYN
        generic :: write(unformatted) => WRITE_T_DYN
        generic :: read(unformatted)  => READ_T_DYN
END TYPE T_DYN

contains

!
!
!_______________________________________________________________________________
! set unformatted writing and reading for T_DYN_WORK
subroutine WRITE_T_SOLVERINFO(tsolverinfo, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_SOLVERINFO),  intent(in)     :: tsolverinfo
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
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

subroutine READ_T_SOLVERINFO(tsolverinfo, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_SOLVERINFO),  intent(inout)     :: tsolverinfo
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
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
subroutine WRITE_T_DYN_WORK(twork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN_WORK), intent(in)        :: twork
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    call write_bin_array(twork%uvnode_rhs, unit, iostat, iomsg)
    call write_bin_array(twork%u_c,        unit, iostat, iomsg)
    call write_bin_array(twork%v_c,        unit, iostat, iomsg)
    call write_bin_array(twork%u_b,        unit, iostat, iomsg)
    call write_bin_array(twork%v_b,        unit, iostat, iomsg)
end subroutine WRITE_T_DYN_WORK

subroutine READ_T_DYN_WORK(twork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN_WORK), intent(inout)        :: twork
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
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
    
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%work
    
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%solverinfo
    
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
    
end subroutine WRITE_T_DYN

subroutine READ_T_DYN(dynamics, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_DYN),         intent(inout)  :: dynamics
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg
    
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
    
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%work
    
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
    
end subroutine READ_T_DYN

END MODULE MOD_DYN
