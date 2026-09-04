module cvmix_utils_addon

!BOP
!\newpage
! !MODULE: cvmix_utils_addon
!
! !DESCRIPTION:
!  addon subroutines for cvmix_tke and cvmix_idemix
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
  public :: cvmix_update_tke
  public :: solve_tridiag

!EOP

contains


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




end module cvmix_utils_addon