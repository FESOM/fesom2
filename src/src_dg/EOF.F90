!
!--------------------------------------------------------------------------------------------
 ! generate random vector
!--------------------------------------------------------------------------------------------
!

SUBROUTINE init_random_seed()
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
          
CALL SYSTEM_CLOCK(COUNT=clock)
          
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
          
DEALLOCATE(seed)
END SUBROUTINE




subroutine myrandom(random)

IMPLICIT NONE
INTEGER, PARAMETER :: R=selected_int_kind(9)
Real(kind = 8), Intent(OUT)      :: random
INTEGER(R), PARAMETER :: multiplier=16807, modulus=2147483647, quotient=127773, rest=2836
REAL(kind = 8), SAVE :: below
INTEGER(R), SAVE :: marsaglia=-1, park_miller=-1
INTEGER(R) :: i,k,j,L
!INTEGER :: m, n, clock
!INTEGER, ALLOCATABLE :: seed
INTEGER :: m, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
          
CALL SYSTEM_CLOCK(COUNT=clock)
          
seed = clock + 37 * (/ (m - 1, m = 1, n) /)
CALL RANDOM_SEED(PUT = seed)

!CALL SYSTEM_CLOCK(COUNT=clock)
!write(*,*), 'clock', clock
!seed = MOD(clock,100000)




!if (seed /= 0) then
below = nearest(1.0,-1.0)
park_miller=ior(ieor(888889999, abs(seed(1))),1)
marsaglia=ieor(777755555,abs(seed(1)))


!Marsaglia random number generator
marsaglia=ieor(marsaglia,ishft(marsaglia,13))
marsaglia=ieor(marsaglia,ishft(marsaglia,-17))
marsaglia=ieor(marsaglia,ishft(marsaglia,5))

!Park-Miller generator
i=park_miller/quotient
park_miller=multiplier*(park_miller-i*quotient)-rest*i
if (park_miller < 0) park_miller=park_miller+modulus
random = below*(ior(iand(modulus,ieor(marsaglia,park_miller)),1))/modulus

DEALLOCATE(seed)
end subroutine myrandom

!--------------------------------------------------------------------------------------------
 ! Box-Muller
!--------------------------------------------------------------------------------------------

subroutine rnd_seed
    implicit none
    integer,     parameter :: SP = kind(1.0)
    integer(SP) :: seed_size, clock
    integer(SP), allocatable :: seed(:)

    call system_clock(clock)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = clock
    call random_seed(put=seed)
    deallocate(seed)
end subroutine rnd_seed

subroutine rnd(r)
    integer,     parameter :: SP = kind(1.0)
    integer(SP), parameter :: DP = selected_real_kind(2 * precision(1.0_SP))                 
    real(DP),    parameter :: PI = 4.0_DP * atan(1.0_DP)  
    real(kind = 8), intent(out) :: r(2)
    real(kind = 8) :: r_u(2)  

    call random_number(r_u(1))
    call random_number(r_u(2))
    r(1) = sqrt(-2 * log(r_u(1))) * cos(2 * PI * r_u(2))
    r(2) = sqrt(-2 * log(r_u(1))) * sin(2 * PI * r_u(2))


end subroutine rnd

!--------------------------------------------------------------------------------------------
 ! EOF * PC
!--------------------------------------------------------------------------------------------


subroutine PCA
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto



IMPLICIT NONE


REAL(kind = 8) :: r

!call myrandom(ran)

call rnd_seed
call rnd(r)


SB_PC1 = SB_PC1 * (SB_l**(dt/86400)) + sqrt(dt/86400)* r !* SB_s 


call rnd_seed
call rnd(r)
SB_PC2 = SB_PC2 * (SB_l2**(dt/86400)) + sqrt(dt/86400)* r ! * SB_s2


call rnd_seed
call rnd(r)
SB_PC3 = SB_PC3 * (SB_l3**(dt/86400)) + sqrt(dt/86400)* r !* SB_s3


call rnd_seed
call rnd(r)
SB_PC4 = SB_PC4 * (SB_l4**(dt/86400)) + sqrt(dt/86400)* r! * SB_s4


call rnd_seed
call rnd(r)
SB_PC5 = SB_PC5 * (SB_l5**(dt/86400)) + sqrt(dt/86400)* r! * SB_s5


call rnd_seed
call rnd(r)
SB_PC6 = SB_PC6 * (SB_l6**(dt/86400)) + sqrt(dt/86400)* r! * SB_s6


call rnd_seed
call rnd(r)
SB_PC7 = SB_PC7 * (SB_l7**(dt/86400)) + sqrt(dt/86400)* r! * SB_s7


!call rnd_seed
!call rnd(r)
!SB_PC8 = SB_PC8 * (SB_l8**(dt/86400)) + sqrt(dt/86400)* r! * SB_s8

!call rnd_seed
!call rnd(r)
!SB_PC9 = SB_PC9 * (SB_l9**(dt/86400)) + sqrt(dt/86400)* r! * SB_s9

!call rnd_seed
!call rnd(r)
!SB_PC10 = SB_PC10 * (SB_l10**(dt/86400)) + sqrt(dt/86400)* r! * SB_s10

!call rnd_seed
!call rnd(r)
!SB_PC11 = SB_PC11 * (SB_l11**(dt/86400)) + sqrt(dt/86400)* r! * SB_s11

!call rnd_seed
!call rnd(r)
!SB_PC12 = SB_PC12 * (SB_l12**(dt/86400)) + sqrt(dt/86400)* r! * SB_s12

!call rnd_seed
!call rnd(r)
!SB_PC13 = SB_PC13 * (SB_l13**(dt/86400)) + sqrt(dt/86400)* r! * SB_s13

!call rnd_seed
!call rnd(r)
!SB_PC14 = SB_PC14 * (SB_l14**(dt/86400)) + sqrt(dt/86400)* r! * SB_s14

!call rnd_seed
!call rnd(r)
!SB_PC15 = SB_PC15 * (SB_l15**(dt/86400)) + sqrt(dt/86400)* r! * SB_s15

!call rnd_seed
!call rnd(r)
!SB_PC16 = SB_PC16 * (SB_l16**(dt/86400)) + sqrt(dt/86400)* r! * SB_s16

!call rnd_seed
!call rnd(r)
!SB_PC17 = SB_PC17 * (SB_l17**(dt/86400)) + sqrt(dt/86400)* r! * SB_s17


end subroutine PCA






