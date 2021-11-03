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
TYPE T_solverinfo
    integer       :: ident   = 1
    integer       :: maxiter = 2000
    integer       :: restart = 15
    integer       :: fillin  = 3
    integer       :: lutype  = 2
    real(kind=WP) :: droptol=1.e-8
    real(kind=WP) :: soltol =1e-10  !1.e-10
    
END TYPE T_solverinfo

!
!
!_______________________________________________________________________________
! set main structure for dynamicss, contains viscosity options and parameters + 
! option for momentum advection 
TYPE T_DYN
    ! instant zonal merdional velocity & Adams-Bashfort rhs
    real(kind=WP), allocatable, dimension(:,:,:):: uv, uv_rhs, uv_rhsAB, fer_uv  

    ! horizontal velocities at nodes
    real(kind=WP), allocatable, dimension(:,:,:):: uvnode, uvnode_rhs
    
    ! instant vertical vel arrays
    real(kind=WP), allocatable, dimension(:,:)  :: w, w_e, w_i, cfl_z, fer_w
    
    ! sea surface height arrays
    real(kind=WP), allocatable, dimension(:)    :: eta_n, d_eta, ssh_rhs, ssh_rhs_old
    
    ! summarizes solver input parameter
    type(t_solverinfo)                          :: solverinfo
    
    
    ! visc_option=...
    ! 1=Harmonic Leith parameterization;
    ! 2=Laplacian+Leith+biharmonic background
    ! 3=Biharmonic Leith parameterization
    ! 4=Biharmonic flow aware
    ! 5=Kinematic (easy) Backscatter
    ! 6=Biharmonic flow aware (viscosity depends on velocity Laplacian)
    ! 7=Biharmonic flow aware (viscosity depends on velocity differences)
    ! 8=Dynamic Backscatter
    integer                                     :: visc_opt      = 5      

    ! gamma0 [m/s],   backgroung viscosity= gamma0*len, it should be as small 
    !                 as possible (keep it < 0.01 m/s).
    ! gamma1 [nodim], for computation of the flow aware viscosity
    ! gamma2 [s/m],   is only used in easy backscatter option
    real(kind=WP)                               :: gamma0_visc   = 0.03
    real(kind=WP)                               :: gamma1_visc   = 0.1
    real(kind=WP)                               :: gamma2_visc   = 0.285

    ! div_c the strength of the modified Leith viscosity, nondimensional, 0.3 -- 1.0
    ! leith the strength of the Leith viscosity
    real(kind=WP)                               :: div_c_visc    = 0.5
    real(kind=WP)                               :: leith_c_visc  = 0.05
    
    ! coefficient for returned sub-gridscale energy, to be used with visc_option=5 
    ! (easy backscatter)
    real(kind=WP)                               :: easy_bs_return= 1.5    

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
    call write_bin_array(dynamics%uvnode_rhs, unit, iostat, iomsg)
    
    call write_bin_array(dynamics%w         , unit, iostat, iomsg)
    call write_bin_array(dynamics%w_e       , unit, iostat, iomsg)
    call write_bin_array(dynamics%w_i       , unit, iostat, iomsg)
    call write_bin_array(dynamics%cfl_z     , unit, iostat, iomsg)
    
    if (Fer_GM) then
        call write_bin_array(dynamics%fer_w     , unit, iostat, iomsg)
        call write_bin_array(dynamics%fer_uv    , unit, iostat, iomsg)
    end if 
    
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_opt
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma0_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma1_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma2_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%div_c_visc
    write(unit, iostat=iostat, iomsg=iomsg) dynamics%leith_c_visc
    
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
    call read_bin_array(dynamics%uvnode_rhs, unit, iostat, iomsg)
    
    call read_bin_array(dynamics%w         , unit, iostat, iomsg)
    call read_bin_array(dynamics%w_e       , unit, iostat, iomsg)
    call read_bin_array(dynamics%w_i       , unit, iostat, iomsg)
    call read_bin_array(dynamics%cfl_z     , unit, iostat, iomsg)
    
    if (Fer_GM) then
        call read_bin_array(dynamics%fer_w     , unit, iostat, iomsg)
        call read_bin_array(dynamics%fer_uv    , unit, iostat, iomsg)
    end if
    
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%visc_opt
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma0_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma1_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%gamma2_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%div_c_visc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%leith_c_visc
    
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%use_ivertvisc
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%momadv_opt
    
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%use_freeslip
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%use_wsplit
    read(unit, iostat=iostat, iomsg=iomsg) dynamics%wsplit_maxcfl
    
end subroutine READ_T_DYN

END MODULE MOD_DYN