module cvmix_utils

!BOP
!\newpage
! !MODULE: cvmix_utils
!
! !AUTHOR:
!  Michael N. Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains routines that are called by multiple modules but don't
!  specifically compute anything mixing related.
!\\
!\\

! !USES:

   use cvmix_kinds_and_types, only : cvmix_r8,                                &
                                     cvmix_strlen,                            &
                                     CVMIX_SUM_OLD_AND_NEW_VALS,              &
                                     CVMIX_MAX_OLD_AND_NEW_VALS,              &
                                     CVMIX_OVERWRITE_OLD_VAL

!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:
  public :: cvmix_update_wrap
  public :: cvmix_att_name
  public :: cvmix_update_tke
  public :: solve_tridiag

!EOP

contains

!BOP

! !IROUTINE: cvmix_update_wrap
! !INTERFACE:

  subroutine cvmix_update_wrap(old_vals, nlev, Mdiff_out, new_Mdiff,          &
                               Tdiff_out, new_Tdiff, Sdiff_out, new_Sdiff)

! !DESCRIPTION:
!  Update diffusivity values based on \verb|old_vals| (either overwrite, sum, or find
!  the level-by-level max)
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    integer, intent(in) :: old_vals, nlev
    real(cvmix_r8), dimension(nlev+1), optional, intent(in) :: new_Mdiff,     &
                                                               new_Tdiff,     &
                                                               new_Sdiff

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(nlev+1), optional, intent(inout) :: Mdiff_out,  &
                                                                  Tdiff_out,  &
                                                                  Sdiff_out

!EOP
!BOC

    integer :: kw

    select case (old_vals)
      case (CVMIX_SUM_OLD_AND_NEW_VALS)
        if ((present(Mdiff_out)).and.(present(new_Mdiff))) &
          Mdiff_out = Mdiff_out + new_Mdiff
        if ((present(Tdiff_out)).and.(present(new_Tdiff))) &
          Tdiff_out = Tdiff_out + new_Tdiff
        if ((present(Sdiff_out)).and.(present(new_Sdiff))) &
          Sdiff_out = Sdiff_out + new_Sdiff
      case (CVMIX_MAX_OLD_AND_NEW_VALS)
        do kw=1,nlev+1
          if ((present(Mdiff_out)).and.(present(new_Mdiff))) &
            Mdiff_out(kw) = max(Mdiff_out(kw), new_Mdiff(kw))
          if ((present(Tdiff_out)).and.(present(new_Tdiff))) &
            Tdiff_out(kw) = max(Tdiff_out(kw), new_Tdiff(kw))
          if ((present(Sdiff_out)).and.(present(new_Sdiff))) &
            Sdiff_out(kw) = max(Sdiff_out(kw), new_Sdiff(kw))
        end do
      case (CVMIX_OVERWRITE_OLD_VAL)
        if ((present(Mdiff_out)).and.(present(new_Mdiff))) &
          Mdiff_out = new_Mdiff
        if ((present(Tdiff_out)).and.(present(new_Tdiff))) &
          Tdiff_out = new_Tdiff
        if ((present(Sdiff_out)).and.(present(new_Sdiff))) &
          Sdiff_out = new_Sdiff
      case DEFAULT
        print*, "ERROR: do not know how to handle old values!"
        stop 1
    end select

!EOC

  end subroutine cvmix_update_wrap

!BOP

! !IROUTINE: cvmix_att_name
! !INTERFACE:

  function cvmix_att_name(varname)

! !DESCRIPTION:
!  Given a variable short name, returns the precise name of the desired
!  attribute in the cvmix\_data\_type structure.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname

! !OUTPUT PARAMETERS:
    character(len=cvmix_strlen) :: cvmix_att_name

!EOP
!BOC

    select case(trim(varname))
      ! Scalars
      case ("nlev", "NumberLevels", "NumberOfLevels")
        cvmix_att_name = "nlev"
      case ("max_nlev", "MaxNumberLevels", "MaxNumberOfLevels")
        cvmix_att_name = "max_nlev"
      case ("depth", "ocn_depth", "OceanDepth", "DepthOfOcean")
        cvmix_att_name = "OceanDepth"
      case ('BoundaryLayerDepth','OBL_depth')
        cvmix_att_name = "BoundaryLayerDepth"
      case ("SSH", "surf_hgt", "SeaSurfaceHeight", "SurfaceHeight", "height")
        cvmix_att_name = "SeaSurfaceHeight"
      case ("surf_fric", "SurfaceFriction")
        cvmix_att_name = "SurfaceFriction"
      case ("surf_buoy", "SurfaceBuoyancy", "SurfaceBuoyancyForcing")
        cvmix_att_name = "SurfaceBuoyancyForcing"
      case ("lat", "latitude", "Latitude")
        cvmix_att_name = "lat"
      case ("lon", "longitude", "Longitude")
        cvmix_att_name = "lon"
      case ("coriolis", "Coriolis", "CoriolisFreq", "CoriolisFrequency")
        cvmix_att_name = "Coriolis"
      case ("kOBL_depth", "BoundaryLayerDepthIndex")
        cvmix_att_name = "kOBL_depth"
      ! QL, 150610
      case ("LangmuirEnhancementFactor", "EnhancementFactor",   &
              "langmuir_Efactor")
        cvmix_att_name = "LangmuirEnhancementFactor"
      case ("SurfaceStokesDrift", "stokes_drift")
        cvmix_att_name = "SurfaceStokesDrift"
! IDEMIX/TKE
      case("forc_iw_bottom")
        cvmix_att_name = "forc_iw_bottom"
      case("forc_iw_surface")
        cvmix_att_name = "forc_iw_surface"
      case("forc_tke_surf")
        cvmix_att_name = "forc_tke_surf"
      case("forc_rho_surf")
        cvmix_att_name = "forc_rho_surf"
      case("dtime")
        cvmix_att_name = "dtime"
! Tidal (Simmons et al)
!      case("forc_tidal_bottom")
!        cvmix_att_name = "forc_tidal_bottom"
 

      ! Variables on level interfaces
      case ("zw", "zw_iface")
        cvmix_att_name = "zw_iface"
      case ("dzw", "dzw_iface")
        cvmix_att_name = "dzw"
      case ("Mdiff", "Udiff", "MomentumDiff", "MomentumDiffusivity")
        cvmix_att_name = "Mdiff_iface"
      case ("Tdiff", "TempDiff", "TemperatureDiff", "TemperatureDiffusivity")
        cvmix_att_name = "Tdiff_iface"
      case ("Sdiff", "SaltDiff", "SalinityDiff", "SalinityDiffusivity")
        cvmix_att_name = "Sdiff_iface"
      case ("Ri", "Ri_iface", "Richardson", "ShearRichardson",                &
            "RichardsonNumber", "ShearRichardsonNumber",                      &
            "ShearRichardson_iface")
        cvmix_att_name = "ShearRichardson_iface"
      case ("buoy", "buoy_iface", "N", "Nsqr", "BuoyancyFreq", "SqrBuoyancy", &
            "SqrBuoyancyFreq", "SqrBuoyancyFreq_iface")
        cvmix_att_name = "SqrBuoyancyFreq_iface"
      case ("kpp_transport", "kpp_nonlocal", "nonlocal_transport",            &
            "nonlocal", "kpp_nonlocal_iface")
        ! Note: this isn't an attribute in the data type, but put / get
        !       uses this as short hand for "both Tnonlocal and Snonlocal"
        cvmix_att_name = "kpp_nonlocal_iface"
      case ("Tnonlocal", "KPP_T_Nonlocal", "kpp_Tnonlocal", "kpp_Ttransport", &
            "kpp_Tnonlocal_iface")
        cvmix_att_name = "kpp_Tnonlocal_iface"
      case ("Snonlocal", "KPP_S_Nonlocal", "kpp_Snonlocal", "kpp_Stransport", &
            "kpp_Snonlocal_iface")
        cvmix_att_name = "kpp_Snonlocal_iface"

! IDEMIX / TKE
      case("KappaM_iface")
        cvmix_att_name = "KappaM_iface"
      case("KappaH_iface")
        cvmix_att_name = "KappaH_iface"
      case("Ssqr_iface")
        cvmix_att_name = "Ssqr_iface"
      case("Nsqr_iface")
        cvmix_att_name = "Nsqr_iface"
      case("TKE")
        cvmix_att_name = "TKE"
      case("E_iw")
        cvmix_att_name = "E_iw"
      case("iw_diss")
        cvmix_att_name = "iw_diss"
      case("alpha_c")
        cvmix_att_name = "alpha_c"

      ! Variables on level centers
      case ("z","zt","zt_cntr")
        cvmix_att_name = "zt_cntr"
      case ("dz", "dzt", "CellThickness")
        cvmix_att_name = "dzt"
      case ("rho", "dens", "WaterDensity", "WaterDensity_cntr")
        cvmix_att_name = "WaterDensity_cntr"
      case ("rho_lwr", "dens_lwr", "AdiabWaterDensity",                       &
            "AdiabWaterDensity_cntr")
        cvmix_att_name = "AdiabWaterDensity_cntr"
      case ("Rib", "Ri_bulk", "BulkRichardson", "BulkRichardsonNumber",       &
            "BulkRichardson_cntr")
        cvmix_att_name = "BulkRichardson_cntr"
      case ("Rrho", "strat_param")
        ! Note: this isn't an attribute in the data type, but the I/O routines
        !       use it to denote strat_param_num / strat_param_denom
        cvmix_att_name = "strat_param"
      case ("Rrho_num", "strat_param_num")
        cvmix_att_name = "strat_param_num"
      case ("Rrho_denom", "strat_param_denom")
        cvmix_att_name = "strat_param_denom"
      case ("Buoyancy","buoyancy","buoyancy_cntr")
        cvmix_att_name = "buoyancy_cntr"
      case ("U", "Vx", "Vx_cntr")
        cvmix_att_name = "Vx_cntr"
      case ("V", "Vy", "Vy_cntr")
        cvmix_att_name = "Vy_cntr"
      case ("SimmonsCoeff", "TidalCoeff")
        cvmix_att_name = "SimmonsCoeff"
      case ("VertDep", "VertDep_iface", "vert_dep")
        cvmix_att_name = "VertDep_iface"
      case DEFAULT
        print*, "ERROR: ", trim(varname), " is not tied to an attribute of ", &
                "the cvmix_data_type structure."
        stop 1
    end select

!EOC
  end function cvmix_att_name


!=================================================================================
  subroutine cvmix_update_tke(old_vals,nlev, tke_diss_out,new_tke_diss, tke_out,new_tke, KappaM_out, new_KappaM,           &
                              KappaH_out, new_KappaH,E_iw_out, new_E_iw,&
                              iw_diss_out, new_iw_diss)
!
!! !DESCRIPTION:
!!  Update diffusivity values based on old_vals (either overwrite, sum, or find
!!  the level-by-level max)
!!  TKE & IDEMIX both use only the overwrite option
!
!! !INPUT PARAMETERS:
    integer, intent(in) :: old_vals, nlev
    real(cvmix_r8), dimension(nlev+1), optional, intent(in) ::     &
      new_KappaM                                                  ,& !
      new_KappaH                                                  ,& !
      new_E_iw                                                    ,& !
      new_iw_diss                                                 ,& !
      new_tke                                                     ,& !
      new_tke_diss                                                   !

!! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(nlev+1), optional, intent(inout) ::  & !
      KappaM_out                                                  ,& !
      KappaH_out                                                  ,& !
      E_iw_out                                                    ,& !
      iw_diss_out                                                 ,& !
      tke_out                                                     ,& !
      tke_diss_out                                                   !  


    integer :: kw


    select case (old_vals)
      case (CVMIX_SUM_OLD_AND_NEW_VALS)
        if ((present(KappaM_out)).and.(present(new_KappaM))) &
          KappaM_out = KappaM_out + new_KappaM
        if ((present(KappaH_out)).and.(present(new_KappaH))) &
          KappaH_out = KappaH_out + new_KappaH

      case (CVMIX_MAX_OLD_AND_NEW_VALS)
        do kw=1,nlev+1
          if ((present(KappaM_out)).and.(present(new_KappaM))) &
            KappaM_out(kw) = max(KappaM_out(kw), new_KappaM(kw))
          if ((present(KappaH_out)).and.(present(new_KappaH))) &
            KappaH_out(kw) = max(KappaH_out(kw), new_KappaH(kw))
        end do
!
!      !This is the only option IDEMIX&TKE use to update to new values  
      case (CVMIX_OVERWRITE_OLD_VAL)
        if ((present(tke_diss_out)).and.(present(new_tke_diss))) then
          tke_diss_out = new_tke_diss
          print*, "updating tke_diss=",tke_diss_out
        end if
        if ((present(KappaH_out)).and.(present(new_KappaH))) &
          KappaH_out = new_KappaH

        if ((present(tke_out)).and.(present(new_tke))) then
          tke_out = new_tke
          print*, "updating tke"
         end if
        if ((present(KappaM_out)).and.(present(new_KappaM))) &
          KappaM_out = new_KappaM

        if ((present(E_iw_out)).and.(present(new_E_iw))) then
          E_iw_out = new_E_iw
        print*, "updating Internal wave energy"
        end if
        if ((present(iw_diss_out)).and.(present(new_iw_diss))) &
          iw_diss_out = new_iw_diss
      case DEFAULT
        print*, "ERROR: do not know how to handle old values!"
        stop 1
    end select

  end subroutine cvmix_update_tke


!=================================================================================
  subroutine solve_tridiag(a,b,c,d,x,n)
        implicit none
  !---------------------------------------------------------------------------------
  !        a - sub-diagonal (means it is the diagonal below the main diagonal)
  !        b - the main diagonal
  !        c - sup-diagonal (means it is the diagonal above the main diagonal)
  !        d - right part
  !        x - the answer
  !        n - number of equations
  !---------------------------------------------------------------------------------
          integer,intent(in) :: n
          real*8,dimension(n),intent(in) :: a,b,c,d
          real*8,dimension(n),intent(out) :: x
          real*8,dimension(n) :: cp,dp
          real*8 :: m,fxa
          integer i
  
  ! initialize c-prime and d-prime
          cp(1) = c(1)/b(1)
          dp(1) = d(1)/b(1)
  ! solve for vectors c-prime and d-prime
           do i = 2,n
             m = b(i)-cp(i-1)*a(i)
             fxa = 1D0/m
             cp(i) = c(i)*fxa
             dp(i) = (d(i)-dp(i-1)*a(i))*fxa
           enddo
  ! initialize x
           x(n) = dp(n)
  ! solve for x from the vectors c-prime and d-prime
          do i = n-1, 1, -1
            x(i) = dp(i)-cp(i)*x(i+1)
          end do
  end subroutine solve_tridiag

!!!subroutine solve_tridiag(a,b,c,d,x,n)
!!!
!!!! This subroutine solves a tri-diagonal matrix
!!!!---------------------------------------------------------------------------------
!!!!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!!!!        b - the main diagonal
!!!!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!!!!        d - right part
!!!!        x - the answer
!!!!        n - number of equations
!!!!---------------------------------------------------------------------------------
!!!
!!!  implicit none
!!!  integer,intent(in) :: n
!!!  real(cvmix_r8),dimension(n),intent(in) :: a,b,c,d
!!!  real(cvmix_r8),dimension(n),intent(out) :: x
!!!  real(cvmix_r8),dimension(n) :: cp,dp
!!!  real(cvmix_r8) :: m,fxa
!!!  integer i
!!!
!!!  ! initialize c-prime and d-prime
!!!  cp(n) = c(n)/b(n)
!!!  dp(n) = d(n)/b(n)
!!!
!!!  ! solve for vectors c-prime and d-prime
!!!  do i = n-1,1,-1
!!!    m = b(i)-cp(i+1)*a(i)
!!!    fxa = 1.0/m
!!!    cp(i) = c(i)*fxa
!!!    dp(i) = (d(i)-dp(i+1)*a(i))*fxa
!!!  enddo
!!!
!!!  ! initialize x 
!!!  x(1) = dp(1)
!!!
!!!  ! solve for x from the vectors c-prime and d-prime
!!!  do i =  2, n
!!!   x(i) = dp(i)-cp(i)*x(i-1)
!!!  end do
!!!
!!!end subroutine solve_tridiag




end module cvmix_utils
