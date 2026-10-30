! ========================================================================
! Halo exchange routines + broadcast routines that collect information
! on the entire field (needed for output)
! The routines here are very similar, difference is the data type and
! exchange pattern.
! exchange_nod2D_i(arr(myDim_nod2D+eDim_nod2D))    INTEGER
! exchange_nod2D(arr(myDim_nod2D+eDim_nod2D))      WP
! exchange_nod3D(arr(nl-1,myDim_nod2D+eDim_nod2D)) WP
! exchange_nod3D_full(arr(nl,myDim_nod2D+eDim_nod2D)) WP
! exchange_edge2D(edge_array2D)     WP  not used currently  !!! no buffer!!!
! exchange_edge3D(edge_array3D)     WP  not used currently  !!! no buffer!!!
! exchange_elem3D(elem_array3D)     WP
! exchange_elem2d_full
! exchange_elem2d_full_i
! ========================================================================

module g_comm

  use, intrinsic :: ISO_FORTRAN_ENV, only: int16, int32, real32, real64
  use o_PARAM, only: WP
#ifdef ENABLE_OPENACC
  use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
#endif

  implicit none

#ifdef ENABLE_OPENACC
  ! Device-to-device exchange of one or two nodal 2D fields (exchange_nod2D_begin, exchange_nod2D_2fields_begin
  ! with luse_g2g: EVP subcycles, SSH solver): halo values packed into contiguous device buffers in slist/rlist
  ! order, unpacked in exchange_nod_end.
  real(kind=WP), allocatable, target :: g2g_sbuf(:), g2g_rbuf(:)
  integer,       allocatable         :: g2g_slist(:), g2g_rlist(:)
  real(kind=WP), pointer             :: g2g_f1(:) => null(), g2g_f2(:) => null()
  logical                            :: g2g_unpack_pending = .false.
  integer                            :: g2g_nfields = 0
  ! The same for whole columns (3D fields, several values per point) of nodal (kind 1) and elemental
  ! (kind 2) fields, one message per neighbour; unpacked in exchange_nod_end / exchange_elem_end.
  type g2g_cols_state
     real(kind=WP), allocatable :: sbuf(:), rbuf(:)
     integer,       allocatable :: slist(:), rlist(:)
     type(c_ptr)                :: f(3)
     integer                    :: nf = 0, inner = 0, ntot = 0
     logical                    :: pending = .false.
  end type g2g_cols_state
  type(g2g_cols_state), target :: g2g_cols(2)
#endif

contains

#ifdef DEBUG
! General version of the communication routine for 2D nodal fields
! Only needed in debug mode
subroutine check_mpi_comm(rn, sn, r_mpitype, s_mpitype, rPE, sPE, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
integer,        intent(in)         :: sn, rn, r_mpitype(:), s_mpitype(:), rPE(:), sPE(:)
integer                            :: n, sdebug, rdebug, status(MPI_STATUS_SIZE), request
DO n=1,rn
   CALL MPI_TYPE_SIZE(r_mpitype(n), rdebug, partit%MPIerr)
   CALL MPI_ISEND(rdebug, 1, MPI_INTEGER, rPE(n), 10, partit%MPI_COMM_FESOM, request, partit%MPIerr)
END DO
DO n=1, sn
   call MPI_RECV(sdebug, 1, MPI_INTEGER, sPE(n), 10, partit%MPI_COMM_FESOM,    &
        status, partit%MPIerr)
   call MPI_TYPE_SIZE(s_mpitype(n), rdebug, partit%MPIerr)
   if (sdebug /= rdebug) then
      print *, "Mismatching MPI send/recieve message lengths."
      print *,"Send/receive process numbers: ", partit%mype, '/', sPE(n)
      print *,"Number of send/receive bytes: ", sdebug, '/', rdebug
      call MPI_ABORT( partit%MPI_COMM_FESOM, 1 )
      end if
END DO
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)
END SUBROUTINE check_mpi_comm
#endif


subroutine exchange_nod2D_i(nod_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
integer,        intent(inout)         :: nod_array2D(:)
logical,        intent(in),optional   :: luse_g2g
if (partit%npes > 1) then
    call exchange_nod2D_i_begin(nod_array2D, partit, luse_g2g)
    call exchange_nod_end(partit)
endif
END SUBROUTINE exchange_nod2D_i

!=============================================================================
! General version of the communication routine for 2D nodal fields
subroutine exchange_nod2D_i_begin(nod_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
integer,        intent(inout)         :: nod_array2D(:)
integer                               :: n, sn, rn
logical,        intent(in),optional   :: luse_g2g
logical                               :: lg2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

  if (partit%npes > 1) then

     sn=partit%com_nod2D%sPEnum
     rn=partit%com_nod2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, partit%r_mpitype_nod2D_i, partit%s_mpitype_nod2D_i,         &
          partit%com_nod2D%rPE, partit%com_nod2D%sPE)
#endif
     if (lg2g) then


     !$ACC HOST_DATA USE_DEVICE(nod_array2D)
         DO n=1,rn

            call MPI_IRECV(nod_array2D, 1, partit%r_mpitype_nod2D_i(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
         END DO

         DO n=1, sn

            call MPI_ISEND(nod_array2D, 1, partit%s_mpitype_nod2D_i(n), partit%com_nod2D%sPE(n), &
                 partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
         END DO
         !$ACC END HOST_DATA
     else

         DO n=1,rn
    
            call MPI_IRECV(nod_array2D, 1, partit%r_mpitype_nod2D_i(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
         END DO
    
         DO n=1, sn
    
            call MPI_ISEND(nod_array2D, 1, partit%s_mpitype_nod2D_i(n), partit%com_nod2D%sPE(n), &
                 partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
         END DO

     endif
     partit%com_nod2D%nreq = rn+sn

  endif
END SUBROUTINE exchange_nod2D_i_begin

! ========================================================================
! General version of the communication routine for 2D nodal fields
subroutine exchange_nod2D(nod_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod_array2D(:)
logical,        intent(in),optional   :: luse_g2g

 if (partit%npes > 1) then
    call exchange_nod2D_begin(nod_array2D, partit, luse_g2g)
    call exchange_nod_end(partit)
 end if

END SUBROUTINE exchange_nod2D

! ========================================================================
! General version of the communication routine for 2D nodal fields
subroutine exchange_nod2D_begin(nod_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout), target :: nod_array2D(:)
integer                               :: n, sn, rn
logical,        intent(in),optional   :: luse_g2g
logical                               :: lg2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

  if (partit%npes > 1) then

     sn=partit%com_nod2D%sPEnum
     rn=partit%com_nod2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, partit%r_mpitype_nod2D, partit%s_mpitype_nod2D,           &
          partit%com_nod2D%rPE, partit%com_nod2D%sPE)
#endif

     if(lg2g) then
#ifdef ENABLE_OPENACC
         call g2g_exchange_1field_begin(nod_array2D, partit)
#else
         !$ACC HOST_DATA USE_DEVICE(nod_array2D)

         DO n=1,rn
            call MPI_IRECV(nod_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
         END DO
         DO n=1, sn
            call MPI_ISEND(nod_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                 partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
         END DO

        !$ACC END HOST_DATA
#endif

     else 
         DO n=1,rn
            call MPI_IRECV(nod_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
         END DO
         DO n=1, sn
            call MPI_ISEND(nod_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                 partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
         END DO
     endif
     partit%com_nod2D%nreq = rn+sn

  end if

END SUBROUTINE exchange_nod2D_begin
!===============================================
! General version of the communication routine for 2D nodal fields
subroutine exchange_nod2D_2fields(nod1_array2D, nod2_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod1_array2D(:)
real(kind=WP),   intent(inout)         :: nod2_array2D(:)
logical,        intent(in),optional   :: luse_g2g

 if (partit%npes > 1) then
    call exchange_nod2D_2fields_begin(nod1_array2D, nod2_array2D, partit, luse_g2g)
    call exchange_nod_end(partit)
 end if

END SUBROUTINE exchange_nod2D_2fields

! ========================================================================
! General version of the communication routine for 2D nodal fields
subroutine exchange_nod2D_2fields_begin(nod1_array2D, nod2_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout), target :: nod1_array2D(:)
real(kind=WP),   intent(inout), target :: nod2_array2D(:)
integer                               :: n, sn, rn
logical,        intent(in),optional   :: luse_g2g
logical                               :: lg2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

if (partit%npes > 1) then

  sn=partit%com_nod2D%sPEnum
  rn=partit%com_nod2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, partit%r_mpitype_nod2D, partit%s_mpitype_nod2D,           &
          partit%com_nod2D%rPE, partit%com_nod2D%sPE)
#endif
  if (lg2g) then 
#ifdef ENABLE_OPENACC
    call g2g_exchange_2fields_begin(nod1_array2D, nod2_array2D, partit)
#else
    !$ACC HOST_DATA USE_DEVICE(nod1_array2D, nod2_array2D) 
  
    DO n=1,rn
       call MPI_IRECV(nod1_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n),      partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n-1), partit%MPIerr)
  
       call MPI_IRECV(nod2_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n)+partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n),   partit%MPIerr)
    END DO
    DO n=1, sn
       call MPI_ISEND(nod1_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype,      partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n-1), partit%MPIerr)
  
       call MPI_ISEND(nod2_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype+partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n),   partit%MPIerr)
    END DO
  
    !$ACC END HOST_DATA
#endif
  else  
    DO n=1,rn
       call MPI_IRECV(nod1_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n),      partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n-1), partit%MPIerr)
  
       call MPI_IRECV(nod2_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n)+partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n),   partit%MPIerr)
    END DO
    DO n=1, sn
       call MPI_ISEND(nod1_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype,      partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n-1), partit%MPIerr)
  
       call MPI_ISEND(nod2_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype+partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n),   partit%MPIerr)
    END DO
  endif
  partit%com_nod2D%nreq = 2*(rn+sn)


end if

END SUBROUTINE exchange_nod2D_2fields_begin

!===============================================
subroutine exchange_nod2D_3fields(nod1_array2D, nod2_array2D, nod3_array2D, partit, luse_g2g)
! General version of the communication routine for 2D nodal fields
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod1_array2D(:)
real(kind=WP),   intent(inout)         :: nod2_array2D(:)
real(kind=WP),   intent(inout)         :: nod3_array2D(:)
logical,        intent(in),optional   :: luse_g2g


 if (partit%npes > 1) then
    call exchange_nod2D_3fields_begin(nod1_array2D, nod2_array2D, nod3_array2D, partit, luse_g2g)
    call exchange_nod_end(partit)
 end if

END SUBROUTINE exchange_nod2D_3fields

! ========================================================================
subroutine exchange_nod2D_3fields_begin(nod1_array2D, nod2_array2D, nod3_array2D, partit, luse_g2g)
! General version of the communication routine for 2D nodal fields
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod1_array2D(:)
real(kind=WP),   intent(inout)         :: nod2_array2D(:)
real(kind=WP),   intent(inout)         :: nod3_array2D(:)
integer                               :: n, sn, rn
logical,        intent(in),optional   :: luse_g2g
logical                               :: lg2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

 if (partit%npes > 1) then

  sn=partit%com_nod2D%sPEnum
  rn=partit%com_nod2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, partit%r_mpitype_nod2D, partit%s_mpitype_nod2D,           &
          partit%com_nod2D%rPE, partit%com_nod2D%sPE)
#endif

  if (lg2g) then
#ifdef ENABLE_OPENACC
     call g2g_cols_begin(1, 1, size(nod1_array2D), partit, nod1_array2D, nod2_array2D, nod3_array2D)
     return
#endif

   !$ACC HOST_DATA USE_DEVICE(nod1_array2D, nod2_array2D, nod3_array2D)

    DO n=1,rn
       call MPI_IRECV(nod1_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n),        partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*n-2), partit%MPIerr)
  
       call MPI_IRECV(nod2_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n)+partit%npes,   partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*n-1), partit%MPIerr)
  
       call MPI_IRECV(nod3_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n)+2*partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*n),   partit%MPIerr)
    END DO
    DO n=1, sn
       call MPI_ISEND(nod1_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype,        partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*rn+3*n-2), partit%MPIerr)
  
       call MPI_ISEND(nod2_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype+partit%npes,   partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*rn+3*n-1), partit%MPIerr)
  
       call MPI_ISEND(nod3_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype+2*partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*rn+3*n),   partit%MPIerr)
    END DO
  
   !$ACC END HOST_DATA
  else
    DO n=1,rn
       call MPI_IRECV(nod1_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n),        partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*n-2), partit%MPIerr)
  
       call MPI_IRECV(nod2_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n)+partit%npes,   partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*n-1), partit%MPIerr)
  
       call MPI_IRECV(nod3_array2D, 1, partit%r_mpitype_nod2D(n), partit%com_nod2D%rPE(n), &
                 partit%com_nod2D%rPE(n)+2*partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*n),   partit%MPIerr)
    END DO
    DO n=1, sn
       call MPI_ISEND(nod1_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype,        partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*rn+3*n-2), partit%MPIerr)
  
       call MPI_ISEND(nod2_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype+partit%npes,   partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*rn+3*n-1), partit%MPIerr)
  
       call MPI_ISEND(nod3_array2D, 1, partit%s_mpitype_nod2D(n), partit%com_nod2D%sPE(n), &
                      partit%mype+2*partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(3*rn+3*n),   partit%MPIerr)
    END DO
  
  endif
  partit%com_nod2D%nreq = 3*(rn+sn)

end if

END SUBROUTINE exchange_nod2D_3fields_begin

! ========================================================================
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
subroutine exchange_nod3D(nod_array3D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod_array3D(:,:)
logical,        intent(in),optional   :: luse_g2g

if (partit%npes > 1) then
   call exchange_nod3D_begin(nod_array3D, partit, luse_g2g)
   call exchange_nod_end(partit)
endif

END SUBROUTINE exchange_nod3D

! ========================================================================
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
subroutine exchange_nod3D_begin(nod_array3D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod_array3D(:,:)
integer                               :: n, sn, rn
integer                               :: nz, nl1
logical,        intent(in),optional   :: luse_g2g
logical                               :: lg2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

 if (partit%npes > 1) then
    sn=partit%com_nod2D%sPEnum
    rn=partit%com_nod2D%rPEnum

    nl1=ubound(nod_array3D,1)

    if ((nl1<ubound(partit%r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(partit%r_mpitype_nod3D, 2))) then
       if (partit%mype==0) then
          print *,'Subroutine exchange_nod3D not implemented for',nl1,'layers.'
          print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
       endif
       call par_abort(partit%MPI_COMM_FESOM, partit%mype)
    endif

    ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
    call check_mpi_comm(rn, sn, partit%r_mpitype_nod3D(:,nl1,1), partit%s_mpitype_nod3D(:,nl1,1), &
         partit%com_nod2D%rPE, partit%com_nod2D%sPE)
#endif
    if(lg2g) then
#ifdef ENABLE_OPENACC
       call g2g_cols_begin(1, size(nod_array3D,1), size(nod_array3D,2), partit, nod_array3D)
       return
#endif
      !$ACC HOST_DATA USE_DEVICE(nod_array3D) 

      DO n=1,rn
         call MPI_IRECV(nod_array3D, 1, partit%r_mpitype_nod3D(n,nl1,1), partit%com_nod2D%rPE(n), &
              partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
      END DO
      DO n=1, sn
         call MPI_ISEND(nod_array3D, 1, partit%s_mpitype_nod3D(n,nl1,1), partit%com_nod2D%sPE(n), &
              partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
      END DO

  
      !$ACC END HOST_DATA
    else
      DO n=1,rn
         call MPI_IRECV(nod_array3D, 1, partit%r_mpitype_nod3D(n,nl1,1), partit%com_nod2D%rPE(n), &
              partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
      END DO
      DO n=1, sn
         call MPI_ISEND(nod_array3D, 1, partit%s_mpitype_nod3D(n,nl1,1), partit%com_nod2D%sPE(n), &
              partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
      END DO
    endif
    partit%com_nod2D%nreq = rn+sn

 endif
END SUBROUTINE exchange_nod3D_begin

! ========================================================================
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
subroutine exchange_nod3D_2fields(nod1_array3D,nod2_array3D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod1_array3D(:,:)
real(kind=WP),   intent(inout)         :: nod2_array3D(:,:)
logical,        intent(in),optional   :: luse_g2g

if (partit%npes > 1) then
   call exchange_nod3D_2fields_begin(nod1_array3D,nod2_array3D, partit, luse_g2g)
   call exchange_nod_end(partit)
endif
END SUBROUTINE exchange_nod3D_2fields

! ========================================================================
subroutine exchange_nod3D_2fields_begin(nod1_array3D,nod2_array3D, partit, luse_g2g)
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod1_array3D(:,:)
real(kind=WP),   intent(inout)         :: nod2_array3D(:,:)
integer                               :: n, sn, rn
integer                               :: nz, nl1, nl2
logical,        intent(in),optional   :: luse_g2g
logical                               :: lg2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

 if (partit%npes > 1) then
    sn=partit%com_nod2D%sPEnum
    rn=partit%com_nod2D%rPEnum

    nl1 = ubound(nod1_array3D,1)

    if ((nl1<ubound(partit%r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(partit%r_mpitype_nod3D, 2))) then
       if (partit%mype==0) then
          print *,'Subroutine exchange_nod3D not implemented for',nl1,'layers.'
          print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
       endif
       call par_abort(partit%MPI_COMM_FESOM, partit%mype)
    endif

    nl2 = ubound(nod2_array3D,1)
    if ((nl2<ubound(partit%r_mpitype_nod3D, 2)-1) .or. (nl2>ubound(partit%r_mpitype_nod3D, 2))) then
       if (partit%mype==0) then
          print *,'Subroutine exchange_nod3D not implemented for',nl2,'layers.'
          print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
       endif
       call par_abort(partit%MPI_COMM_FESOM, partit%mype)
    endif

#ifdef DEBUG
    call check_mpi_comm(rn, sn, partit%r_mpitype_nod3D(:,nl1,1), partit%s_mpitype_nod3D(:,nl1,1), &
         partit%com_nod2D%rPE, partit%com_nod2D%sPE)
#endif

  if(lg2g) then 
#ifdef ENABLE_OPENACC
     call g2g_cols_begin(1, size(nod1_array3D,1), size(nod1_array3D,2), partit, nod1_array3D, nod2_array3D)
     return
#endif
     !$ACC HOST_DATA USE_DEVICE(nod1_array3D, nod2_array3D) 

    DO n=1,rn
       call MPI_IRECV(nod1_array3D, 1, partit%r_mpitype_nod3D(n,nl1,1), partit%com_nod2D%rPE(n), &
            partit%com_nod2D%rPE(n),      partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n-1), partit%MPIerr)

       call MPI_IRECV(nod2_array3D, 1, partit%r_mpitype_nod3D(n,nl2,1), partit%com_nod2D%rPE(n), &
            partit%com_nod2D%rPE(n)+partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n  ), partit%MPIerr)
    END DO

    DO n=1, sn
       call MPI_ISEND(nod1_array3D, 1, partit%s_mpitype_nod3D(n,nl1,1), partit%com_nod2D%sPE(n), &
            partit%mype,      partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n-1), partit%MPIerr)

       call MPI_ISEND(nod2_array3D, 1, partit%s_mpitype_nod3D(n,nl2,1), partit%com_nod2D%sPE(n), &
            partit%mype+partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n), partit%MPIerr)
    END DO

     !$ACC END HOST_DATA
  else
    DO n=1,rn
       call MPI_IRECV(nod1_array3D, 1, partit%r_mpitype_nod3D(n,nl1,1), partit%com_nod2D%rPE(n), &
            partit%com_nod2D%rPE(n),      partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n-1), partit%MPIerr)

       call MPI_IRECV(nod2_array3D, 1, partit%r_mpitype_nod3D(n,nl2,1), partit%com_nod2D%rPE(n), &
            partit%com_nod2D%rPE(n)+partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n  ), partit%MPIerr)
    END DO

    DO n=1, sn
       call MPI_ISEND(nod1_array3D, 1, partit%s_mpitype_nod3D(n,nl1,1), partit%com_nod2D%sPE(n), &
            partit%mype,      partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n-1), partit%MPIerr)

       call MPI_ISEND(nod2_array3D, 1, partit%s_mpitype_nod3D(n,nl2,1), partit%com_nod2D%sPE(n), &
            partit%mype+partit%npes, partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n), partit%MPIerr)
    END DO

  endif
  partit%com_nod2D%nreq = 2*(rn+sn)

 endif
END SUBROUTINE exchange_nod3D_2fields_begin
! ========================================================================
subroutine exchange_nod3D_n(nod_array3D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)          :: nod_array3D(:,:,:)
logical,        intent(in),optional    :: luse_g2g
if (partit%npes>1) then
   call exchange_nod3D_n_begin(nod_array3D, partit, luse_g2g)
   call exchange_nod_end(partit)
endif

END SUBROUTINE exchange_nod3D_n

!=================================================
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
subroutine exchange_nod3D_n_begin(nod_array3D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)          :: nod_array3D(:,:,:)
integer                                :: n, sn, rn
integer                                :: nz, nl1, n_val
logical,        intent(in),optional    :: luse_g2g
logical                                :: lg2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

if (partit%npes>1) then
 ! nod_array3D(n_val,nl1,nod2D_size)
  nl1   = ubound(nod_array3D,2)
  n_val = ubound(nod_array3D,1)
  if ((nl1<ubound(partit%r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(partit%r_mpitype_nod3D, 2)) .or. (n_val > 3)) then
     ! This routine also works for swapped dimensions nod_array3D(nl1,n_val, nod2D_size)
     nl1   = ubound(nod_array3D,1)
     n_val = ubound(nod_array3D,2)

     if ((nl1<ubound(partit%r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(partit%r_mpitype_nod3D, 2)) .or. (n_val > 3)) then
        if (partit%mype==0) then
           print *,'Subroutine exchange_nod3D_n not implemented for'
           print *,nl1,'layers and / or ',n_val,'values per element.'
           print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
        endif
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
     endif
  endif
  sn=partit%com_nod2D%sPEnum
  rn=partit%com_nod2D%rPEnum

  ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
  call check_mpi_comm(rn, sn, partit%r_mpitype_nod3D(:,nl1,n_val),                 &
       partit%s_mpitype_nod3D(:,nl1,n_val), partit%com_nod2D%rPE, partit%com_nod2D%sPE)
#endif

  if(lg2g) then
#ifdef ENABLE_OPENACC
     call g2g_cols_begin(1, size(nod_array3D,1)*size(nod_array3D,2), size(nod_array3D,3), partit, nod_array3D)
     return
#endif
    !$ACC HOST_DATA USE_DEVICE(nod_array3D)
  
    DO n=1,rn
       call MPI_IRECV(nod_array3D, 1, partit%r_mpitype_nod3D(n,nl1,n_val), partit%com_nod2D%rPE(n), &
            partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
    END DO
  
    DO n=1, sn
       call MPI_ISEND(nod_array3D, 1, partit%s_mpitype_nod3D(n,nl1,n_val), partit%com_nod2D%sPE(n), &
            partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
    END DO
  
  !$ACC END HOST_DATA

  else
    DO n=1,rn
       call MPI_IRECV(nod_array3D, 1, partit%r_mpitype_nod3D(n,nl1,n_val), partit%com_nod2D%rPE(n), &
            partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
    END DO
  
    DO n=1, sn
       call MPI_ISEND(nod_array3D, 1, partit%s_mpitype_nod3D(n,nl1,n_val), partit%com_nod2D%sPE(n), &
            partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
    END DO
  endif
  partit%com_nod2D%nreq = rn+sn

 endif


END SUBROUTINE exchange_nod3D_n_begin

!=======================================
! AND WAITING
!=======================================

SUBROUTINE exchange_nod_end(partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit

if (partit%npes > 1) &
  call MPI_WAITALL(partit%com_nod2D%nreq, partit%com_nod2D%req, MPI_STATUSES_IGNORE, partit%MPIerr)
#ifdef ENABLE_OPENACC
if (g2g_unpack_pending) call g2g_unpack_2fields(partit)
if (g2g_cols(1)%pending) call g2g_cols_unpack(1, partit)
#endif

END SUBROUTINE exchange_nod_end

#ifdef ENABLE_OPENACC
! Pack the send halo of two device fields, then post contiguous device-buffer messages
! (same tags as the derived-datatype version). Unpacked by g2g_unpack_2fields.
SUBROUTINE g2g_exchange_2fields_begin(f1, f2, partit)
USE MOD_PARTIT
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),  intent(inout), target :: f1(:), f2(:)
integer :: i, n, sn, rn, ns, nr, a, l

sn=partit%com_nod2D%sPEnum
rn=partit%com_nod2D%rPEnum
ns=partit%com_nod2D%sptr(sn+1)-1
nr=partit%com_nod2D%rptr(rn+1)-1
call g2g_alloc(partit)

!$ACC PARALLEL LOOP DEFAULT(PRESENT)
do i=1, ns
    g2g_sbuf(i)   =f1(g2g_slist(i))
    g2g_sbuf(ns+i)=f2(g2g_slist(i))
end do
!$ACC END PARALLEL LOOP

!$ACC HOST_DATA USE_DEVICE(g2g_sbuf, g2g_rbuf)
do n=1, rn
    a=partit%com_nod2D%rptr(n)
    l=partit%com_nod2D%rptr(n+1)-a
    call MPI_IRECV(g2g_rbuf(a),    l, MPI_WP, partit%com_nod2D%rPE(n), partit%com_nod2D%rPE(n), &
                   partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n-1), partit%MPIerr)
    call MPI_IRECV(g2g_rbuf(nr+a), l, MPI_WP, partit%com_nod2D%rPE(n), partit%com_nod2D%rPE(n)+partit%npes, &
                   partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*n),   partit%MPIerr)
end do
do n=1, sn
    a=partit%com_nod2D%sptr(n)
    l=partit%com_nod2D%sptr(n+1)-a
    call MPI_ISEND(g2g_sbuf(a),    l, MPI_WP, partit%com_nod2D%sPE(n), partit%mype, &
                   partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n-1), partit%MPIerr)
    call MPI_ISEND(g2g_sbuf(ns+a), l, MPI_WP, partit%com_nod2D%sPE(n), partit%mype+partit%npes, &
                   partit%MPI_COMM_FESOM, partit%com_nod2D%req(2*rn+2*n),   partit%MPIerr)
end do
!$ACC END HOST_DATA

g2g_f1 => f1
g2g_f2 => f2
g2g_nfields = 2
g2g_unpack_pending = .true.
END SUBROUTINE g2g_exchange_2fields_begin

! The same for one field (requests and tags as in exchange_nod2D_begin).
SUBROUTINE g2g_exchange_1field_begin(f1, partit)
USE MOD_PARTIT
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),  intent(inout), target :: f1(:)
integer :: i, n, sn, rn, ns, nr, a, l

sn=partit%com_nod2D%sPEnum
rn=partit%com_nod2D%rPEnum
ns=partit%com_nod2D%sptr(sn+1)-1
nr=partit%com_nod2D%rptr(rn+1)-1
call g2g_alloc(partit)

!$ACC PARALLEL LOOP DEFAULT(PRESENT)
do i=1, ns
    g2g_sbuf(i)=f1(g2g_slist(i))
end do
!$ACC END PARALLEL LOOP

!$ACC HOST_DATA USE_DEVICE(g2g_sbuf, g2g_rbuf)
do n=1, rn
    a=partit%com_nod2D%rptr(n)
    l=partit%com_nod2D%rptr(n+1)-a
    call MPI_IRECV(g2g_rbuf(a), l, MPI_WP, partit%com_nod2D%rPE(n), partit%com_nod2D%rPE(n), &
                   partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
end do
do n=1, sn
    a=partit%com_nod2D%sptr(n)
    l=partit%com_nod2D%sptr(n+1)-a
    call MPI_ISEND(g2g_sbuf(a), l, MPI_WP, partit%com_nod2D%sPE(n), partit%mype, &
                   partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
end do
!$ACC END HOST_DATA

g2g_f1 => f1
g2g_nfields = 1
g2g_unpack_pending = .true.
END SUBROUTINE g2g_exchange_1field_begin

! Buffers for two fields and the send/receive node lists on the device, once.
SUBROUTINE g2g_alloc(partit)
USE MOD_PARTIT
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
integer :: ns, nr
if (allocated(g2g_sbuf)) return
ns=partit%com_nod2D%sptr(partit%com_nod2D%sPEnum+1)-1
nr=partit%com_nod2D%rptr(partit%com_nod2D%rPEnum+1)-1
allocate(g2g_sbuf(max(1,2*ns)), g2g_rbuf(max(1,2*nr)), g2g_slist(max(1,ns)), g2g_rlist(max(1,nr)))
g2g_slist(1:ns)=partit%com_nod2D%slist(1:ns)
g2g_rlist(1:nr)=partit%com_nod2D%rlist(1:nr)
!$ACC ENTER DATA CREATE(g2g_sbuf, g2g_rbuf) COPYIN(g2g_slist, g2g_rlist)
END SUBROUTINE g2g_alloc

SUBROUTINE g2g_unpack_2fields(partit)
USE MOD_PARTIT
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),  pointer               :: f1(:), f2(:)
integer :: i, nr

nr=partit%com_nod2D%rptr(partit%com_nod2D%rPEnum+1)-1
f1 => g2g_f1
if (g2g_nfields == 2) then
    f2 => g2g_f2
!$ACC PARALLEL LOOP DEFAULT(PRESENT)
    do i=1, nr
        f1(g2g_rlist(i))=g2g_rbuf(i)
        f2(g2g_rlist(i))=g2g_rbuf(nr+i)
    end do
!$ACC END PARALLEL LOOP
else
!$ACC PARALLEL LOOP DEFAULT(PRESENT)
    do i=1, nr
        f1(g2g_rlist(i))=g2g_rbuf(i)
    end do
!$ACC END PARALLEL LOOP
end if
g2g_unpack_pending = .false.
g2g_nfields = 0
nullify(g2g_f1, g2g_f2)
END SUBROUTINE g2g_unpack_2fields

! Packed device exchange of nf (1-3) contiguous device fields over the nodal (kind=1) or elemental (kind=2)
! halo: whole columns of the listed points, same values as the MPI indexed datatypes, one message per neighbour.
SUBROUTINE g2g_cols_begin(kind, inner, ntot, partit, a1, a2, a3)
USE MOD_PARTIT
IMPLICIT NONE
integer,        intent(in)                     :: kind, inner, ntot
type(t_partit), intent(inout), target          :: partit
real(kind=WP),  intent(inout), target          :: a1(*)
real(kind=WP),  intent(inout), target, optional:: a2(*), a3(*)
type(g2g_cols_state), pointer :: st
type(com_struct),     pointer :: com
real(kind=WP), pointer :: f1(:), f2(:), f3(:), sb(:), rb(:)
integer, pointer       :: sl(:)
integer :: i, j, n, sn, rn, ns, nr, nf, a, l, blk

st => g2g_cols(kind)
if (kind==1) then
   com => partit%com_nod2D
else
   com => partit%com_elem2D
end if
sn=com%sPEnum
rn=com%rPEnum
ns=com%sptr(sn+1)-1
nr=com%rptr(rn+1)-1
nf=1
if (present(a2)) nf=2
if (present(a3)) nf=3
blk=nf*inner
call g2g_cols_alloc(st, com, blk*ns, blk*nr)

st%f(1)=c_loc(a1(1)); call c_f_pointer(st%f(1), f1, [inner*ntot])
if (nf>=2) then; st%f(2)=c_loc(a2(1)); call c_f_pointer(st%f(2), f2, [inner*ntot]); end if
if (nf==3) then; st%f(3)=c_loc(a3(1)); call c_f_pointer(st%f(3), f3, [inner*ntot]); end if
sb => st%sbuf
rb => st%rbuf
sl => st%slist
!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(PRESENT)
do i=1, ns
   do j=1, inner
      sb((i-1)*blk+j)=f1((sl(i)-1)*inner+j)
   end do
end do
!$ACC END PARALLEL LOOP
if (nf>=2) then
!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(PRESENT)
do i=1, ns
   do j=1, inner
      sb((i-1)*blk+inner+j)=f2((sl(i)-1)*inner+j)
   end do
end do
!$ACC END PARALLEL LOOP
end if
if (nf==3) then
!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(PRESENT)
do i=1, ns
   do j=1, inner
      sb((i-1)*blk+2*inner+j)=f3((sl(i)-1)*inner+j)
   end do
end do
!$ACC END PARALLEL LOOP
end if

!$ACC HOST_DATA USE_DEVICE(sb, rb)
do n=1, rn
   a=(com%rptr(n)-1)*blk+1
   l=(com%rptr(n+1)-com%rptr(n))*blk
   call MPI_IRECV(rb(a), l, MPI_WP, com%rPE(n), com%rPE(n), partit%MPI_COMM_FESOM, com%req(n), partit%MPIerr)
end do
do n=1, sn
   a=(com%sptr(n)-1)*blk+1
   l=(com%sptr(n+1)-com%sptr(n))*blk
   call MPI_ISEND(sb(a), l, MPI_WP, com%sPE(n), partit%mype, partit%MPI_COMM_FESOM, com%req(rn+n), partit%MPIerr)
end do
!$ACC END HOST_DATA
com%nreq=rn+sn
st%nf=nf
st%inner=inner
st%ntot=ntot
st%pending=.true.
END SUBROUTINE g2g_cols_begin

SUBROUTINE g2g_cols_unpack(kind, partit)
USE MOD_PARTIT
IMPLICIT NONE
integer,        intent(in)            :: kind
type(t_partit), intent(inout), target :: partit
type(g2g_cols_state), pointer :: st
real(kind=WP), pointer :: f(:), rb(:)
integer, pointer       :: rl(:)
integer :: i, j, k, nr, inner, blk

st => g2g_cols(kind)
if (kind==1) then
   nr=partit%com_nod2D%rptr(partit%com_nod2D%rPEnum+1)-1
else
   nr=partit%com_elem2D%rptr(partit%com_elem2D%rPEnum+1)-1
end if
inner=st%inner
blk=st%nf*inner
rb => st%rbuf
rl => st%rlist
do k=1, st%nf
   call c_f_pointer(st%f(k), f, [inner*st%ntot])
!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(PRESENT)
   do i=1, nr
      do j=1, inner
         f((rl(i)-1)*inner+j)=rb((i-1)*blk+(k-1)*inner+j)
      end do
   end do
!$ACC END PARALLEL LOOP
end do
st%pending=.false.
st%nf=0
END SUBROUTINE g2g_cols_unpack

! Device buffers (grown on demand) and lists of one halo kind.
SUBROUTINE g2g_cols_alloc(st, com, nsb, nrb)
USE MOD_PARTIT
IMPLICIT NONE
type(g2g_cols_state), intent(inout) :: st
type(com_struct),     intent(in)    :: com
integer,              intent(in)    :: nsb, nrb
integer :: ns, nr
if (.not. allocated(st%slist)) then
   ns=com%sptr(com%sPEnum+1)-1
   nr=com%rptr(com%rPEnum+1)-1
   allocate(st%slist(max(1,ns)), st%rlist(max(1,nr)))
   st%slist(1:ns)=com%slist(1:ns)
   st%rlist(1:nr)=com%rlist(1:nr)
   !$ACC ENTER DATA COPYIN(st%slist, st%rlist)
end if
if (allocated(st%sbuf)) then
   if (size(st%sbuf)>=nsb .and. size(st%rbuf)>=nrb) return
   !$ACC EXIT DATA DELETE(st%sbuf, st%rbuf)
   deallocate(st%sbuf, st%rbuf)
end if
allocate(st%sbuf(max(1,nsb)), st%rbuf(max(1,nrb)))
!$ACC ENTER DATA CREATE(st%sbuf, st%rbuf)
END SUBROUTINE g2g_cols_alloc
#endif

SUBROUTINE exchange_elem_end(partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit

  if (partit%npes > 1) then
     if (partit%elem_full_flag) then
        call MPI_WAITALL(partit%com_elem2D_full%nreq, &
             partit%com_elem2D_full%req, MPI_STATUSES_IGNORE, partit%MPIerr)
     else
        call MPI_WAITALL(partit%com_elem2D%nreq, &
             partit%com_elem2D%req, MPI_STATUSES_IGNORE, partit%MPIerr)
     endif
  end if
#ifdef ENABLE_OPENACC
  if (g2g_cols(2)%pending) call g2g_cols_unpack(2, partit)
#endif
END SUBROUTINE exchange_elem_end
!=============================================================================
subroutine exchange_elem3D(elem_array3D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)          :: elem_array3D(:,:)
logical,        intent(in),optional    :: luse_g2g

call exchange_elem3D_begin(elem_array3D, partit, luse_g2g)
call exchange_elem_end(partit)

END SUBROUTINE exchange_elem3D
!===========================================
! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
subroutine exchange_elem3D_begin(elem_array3D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)       :: elem_array3D(:,:)
integer                             :: n, sn, rn, nl1
logical                                :: lg2g
logical,        intent(in),optional    :: luse_g2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

if (partit%npes> 1) then

   nl1=ubound(elem_array3D,1)

   if (ubound(elem_array3D,2)<=partit%myDim_elem2D+partit%eDim_elem2D) then

      partit%elem_full_flag = .false.

      sn=partit%com_elem2D%sPEnum
      rn=partit%com_elem2D%rPEnum

      if (nl1==ubound(partit%r_mpitype_elem3D, 2) .or. nl1==ubound(partit%r_mpitype_elem3D, 2)-1) then

         ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
         call check_mpi_comm(rn, sn, partit%r_mpitype_elem3D(:,nl1,1), partit%s_mpitype_elem3D(:,nl1,1), &
              partit%com_elem2D%rPE, partit%com_elem2D%sPE)
#endif

        if(lg2g) then
#ifdef ENABLE_OPENACC
           call g2g_cols_begin(2, size(elem_array3D,1), size(elem_array3D,2), partit, elem_array3D)
           return
#endif
         !$ACC HOST_DATA USE_DEVICE(elem_array3D)
           DO n=1,rn
              call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem3D(n,nl1,1), partit%com_elem2D%rPE(n), &
                   partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, &
                   partit%com_elem2D%req(n), partit%MPIerr)
           END DO
           DO n=1, sn
              call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem3D(n,nl1,1), partit%com_elem2D%sPE(n), &
                   partit%mype,    partit%MPI_COMM_FESOM, &
                   partit%com_elem2D%req(rn+n), partit%MPIerr)
           END DO
         !$ACC END HOST_DATA
        else 
           DO n=1,rn
              call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem3D(n,nl1,1), partit%com_elem2D%rPE(n), &
                   partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, &
                   partit%com_elem2D%req(n), partit%MPIerr)
           END DO
           DO n=1, sn
              call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem3D(n,nl1,1), partit%com_elem2D%sPE(n), &
                   partit%mype,    partit%MPI_COMM_FESOM, &
                   partit%com_elem2D%req(rn+n), partit%MPIerr)
           END DO
        endif 
      elseif (nl1 <= 4) then
         ! In fact, this is a 2D-array with up to 4 values, e.g. derivatives

         ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
         call check_mpi_comm(rn, sn, partit%r_mpitype_elem2D(:,nl1), partit%s_mpitype_elem2D(:,nl1), &
              partit%com_elem2D%rPE, partit%com_elem2D%sPE)
#endif
        if(lg2g) then
#ifdef ENABLE_OPENACC
           call g2g_cols_begin(2, size(elem_array3D,1), size(elem_array3D,2), partit, elem_array3D)
           return
#endif
          !$ACC HOST_DATA USE_DEVICE(elem_array3D)
           DO n=1,rn
              call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem2D(n,nl1), partit%com_elem2D%rPE(n), &
                   partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, &
                   partit%com_elem2D%req(n), partit%MPIerr)
           END DO
           DO n=1, sn
              call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem2D(n,nl1), partit%com_elem2D%sPE(n), &
                   partit%mype,    partit%MPI_COMM_FESOM, &
                   partit%com_elem2D%req(rn+n), partit%MPIerr)
           END DO
         !$ACC END HOST_DATA
  
        else 
           DO n=1,rn
              call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem2D(n,nl1), partit%com_elem2D%rPE(n), &
                   partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, &
                   partit%com_elem2D%req(n), partit%MPIerr)
           END DO
           DO n=1, sn
              call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem2D(n,nl1), partit%com_elem2D%sPE(n), &
                   partit%mype,    partit%MPI_COMM_FESOM, &
                   partit%com_elem2D%req(rn+n), partit%MPIerr)
           END DO

        endif
      else
         if (partit%mype==0) print *,'Sorry, no MPI datatype prepared for',nl1,'values per element (exchange_elem3D)'
         call par_abort(partit%MPI_COMM_FESOM, partit%mype)
      endif

      partit%com_elem2D%nreq = rn+sn

   else

      partit%elem_full_flag = .true.

      sn=partit%com_elem2D_full%sPEnum
      rn=partit%com_elem2D_full%rPEnum

      if (nl1==ubound(partit%r_mpitype_elem3D_full, 2) .or. nl1==ubound(partit%r_mpitype_elem3D_full, 2)-1) then
         ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
         call check_mpi_comm(rn, sn, partit%r_mpitype_elem3D_full(:,nl1,1), &
              partit%s_mpitype_elem3D_full(:,nl1,1), partit%com_elem2D_full%rPE, partit%com_elem2D_full%sPE)
#endif

         if(lg2g) then
           !$ACC HOST_DATA USE_DEVICE(elem_array3D)
           DO n=1,rn
              call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem3D_full(n,nl1,1), &
                   partit%com_elem2D_full%rPE(n), &
                   partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, &
                   partit%com_elem2D_full%req(n), partit%MPIerr)
           END DO
           DO n=1, sn
              call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem3D_full(n,nl1,1), &
                   partit%com_elem2D_full%sPE(n), &
                   partit%mype,    partit%MPI_COMM_FESOM, &
                   partit%com_elem2D_full%req(rn+n), partit%MPIerr)
           END DO
           !$ACC END HOST_DATA
         else
           DO n=1,rn
              call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem3D_full(n,nl1,1), &
                   partit%com_elem2D_full%rPE(n), &
                   partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, &
                   partit%com_elem2D_full%req(n), partit%MPIerr)
           END DO
           DO n=1, sn
              call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem3D_full(n,nl1,1), &
                   partit%com_elem2D_full%sPE(n), &
                   partit%mype,    partit%MPI_COMM_FESOM, &
                   partit%com_elem2D_full%req(rn+n), partit%MPIerr)
           END DO
         endif
      elseif (nl1 <= 4) then
         ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
         call check_mpi_comm(rn, sn, partit%r_mpitype_elem2D_full(:,nl1), &
              partit%s_mpitype_elem2D_full(:,nl1), partit%com_elem2D_full%rPE, partit%com_elem2D_full%sPE)
#endif

         ! In fact, this is a 2D-array with up to 4 values, e.g. derivatives
         if(lg2g) then
           !$ACC HOST_DATA USE_DEVICE(elem_array3D)
           DO n=1,rn
              call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem2D_full(n,nl1), &
                   partit%com_elem2D_full%rPE(n), &
                   partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, &
                   partit%com_elem2D_full%req(n), partit%MPIerr)
           END DO
           DO n=1, sn
              call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem2D_full(n,nl1), &
                   partit%com_elem2D_full%sPE(n), &
                   partit%mype,    partit%MPI_COMM_FESOM, &
                   partit%com_elem2D_full%req(rn+n), partit%MPIerr)
           END DO
           !$ACC END HOST_DATA
         else
           DO n=1,rn
              call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem2D_full(n,nl1), &
                   partit%com_elem2D_full%rPE(n), &
                   partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, &
                   partit%com_elem2D_full%req(n), partit%MPIerr)
           END DO
           DO n=1, sn
              call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem2D_full(n,nl1), &
                   partit%com_elem2D_full%sPE(n), &
                   partit%mype,    partit%MPI_COMM_FESOM, &
                   partit%com_elem2D_full%req(rn+n), partit%MPIerr)
           END DO
         endif
      else
         if (partit%mype==0) print *,'Sorry, no MPI datatype prepared for',nl1,'values per element (exchange_elem3D)'
         call par_abort(partit%MPI_COMM_FESOM, partit%mype)
      endif

      partit%com_elem2D_full%nreq = rn+sn

   endif

endif

END SUBROUTINE exchange_elem3D_begin

!=============================================================================
! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
subroutine exchange_elem3D_n(elem_array3D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)       :: elem_array3D(:,:,:)
logical,        intent(in),optional    :: luse_g2g

 if (partit%npes> 1) then
    call exchange_elem3D_n_begin(elem_array3D, partit, luse_g2g)
    call exchange_elem_end(partit)
 endif
END SUBROUTINE exchange_elem3D_n
!=============================================================================
subroutine exchange_elem3D_n_begin(elem_array3D, partit, luse_g2g)
! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)       :: elem_array3D(:,:,:)
integer                             :: n, sn, rn, n_val, nl1
logical                                :: lg2g
logical,        intent(in),optional    :: luse_g2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

 if (partit%npes> 1) then
 nl1   = ubound(elem_array3D,2)
 n_val = ubound(elem_array3D,1)

  if ((nl1<ubound(partit%r_mpitype_elem3D, 2)-1) .or. (nl1>ubound(partit%r_mpitype_elem3D, 2)) .or. (n_val > 4)) then

     ! This routine also works for swapped dimensions elem_array3D(nl1,n_val, elem2D_size)
     nl1= ubound(elem_array3D,1)
     n_val = ubound(elem_array3D,2)

     if ((nl1<ubound(partit%r_mpitype_elem3D, 2)-1) .or. (nl1>ubound(partit%r_mpitype_elem3D, 2)) .or. (n_val > 4)) then
        if (partit%mype==0) then
           print *,'Subroutine exchange_elem3D_n not implemented for'
           print *,nl1,'layers and / or ',n_val,'values per element.'
           print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
        endif
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
     endif
  endif

 if (ubound(elem_array3D,3)<=partit%myDim_elem2D+partit%eDim_elem2D) then

    partit%elem_full_flag = .false.

     sn=partit%com_elem2D%sPEnum
     rn=partit%com_elem2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, partit%r_mpitype_elem3D(:,nl1,n_val), &
          partit%s_mpitype_elem3D(:,nl1,n_val), partit%com_elem2D%rPE, partit%com_elem2D%sPE)
#endif
     if(lg2g) then
#ifdef ENABLE_OPENACC
        call g2g_cols_begin(2, size(elem_array3D,1)*size(elem_array3D,2), size(elem_array3D,3), partit, elem_array3D)
        return
#endif
       !$ACC HOST_DATA USE_DEVICE(elem_array3D)
       DO n=1,rn
          call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem3D(n,nl1,n_val), partit%com_elem2D%rPE(n), &
                         partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D%req(n), partit%MPIerr)
       END DO
       DO n=1, sn
          call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem3D(n,nl1,n_val), partit%com_elem2D%sPE(n), &
                         partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D%req(rn+n), partit%MPIerr)
       END DO
       !$ACC END HOST_DATA
     else
       DO n=1,rn
          call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem3D(n,nl1,n_val), partit%com_elem2D%rPE(n), &
                         partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D%req(n), partit%MPIerr)
       END DO
       DO n=1, sn
          call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem3D(n,nl1,n_val), partit%com_elem2D%sPE(n), &
                         partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D%req(rn+n), partit%MPIerr)
       END DO
     endif
     partit%com_elem2D%nreq = rn+sn

  else

     partit%elem_full_flag = .true.

     sn=partit%com_elem2D_full%sPEnum
     rn=partit%com_elem2D_full%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, partit%r_mpitype_elem3D_full(:,nl1,n_val), &
          partit%s_mpitype_elem3D_full(:,nl1,n_val), partit%com_elem2D_full%rPE, partit%com_elem2D_full%sPE)
#endif
     if(lg2g) then
       !$ACC HOST_DATA USE_DEVICE(elem_array3D)

       DO n=1,rn
          call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem3D_full(n,nl1,n_val), partit%com_elem2D_full%rPE(n), &
                         partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(n), partit%MPIerr)
       END DO
       DO n=1, sn
          call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem3D_full(n,nl1,n_val), partit%com_elem2D_full%sPE(n), &
                         partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(rn+n), partit%MPIerr)
       END DO
       !$ACC END HOST_DATA
     else
       DO n=1,rn
          call MPI_IRECV(elem_array3D, 1, partit%r_mpitype_elem3D_full(n,nl1,n_val), partit%com_elem2D_full%rPE(n), &
                         partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(n), partit%MPIerr)
       END DO
       DO n=1, sn
          call MPI_ISEND(elem_array3D, 1, partit%s_mpitype_elem3D_full(n,nl1,n_val), partit%com_elem2D_full%sPE(n), &
                         partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(rn+n), partit%MPIerr)
       END DO
     endif

     partit%com_elem2D_full%nreq = rn+sn

  end if


endif
END SUBROUTINE exchange_elem3D_n_begin
!========================================================================
! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
subroutine exchange_elem2D(elem_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)       :: elem_array2D(:)
logical,        intent(in),optional    :: luse_g2g

 if (partit%npes> 1) then
    call exchange_elem2D_begin(elem_array2D, partit, luse_g2g)
    call exchange_elem_end(partit)
 end if

END SUBROUTINE exchange_elem2D
!========================================================================
! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
subroutine exchange_elem2D_begin(elem_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)       :: elem_array2D(:)
integer                             :: n, sn, rn
logical                                :: lg2g
logical,        intent(in),optional    :: luse_g2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

if (partit%npes> 1) then

  if (ubound(elem_array2D,1)<=partit%myDim_elem2D+partit%eDim_elem2D) then

     partit%elem_full_flag = .false.

     sn=partit%com_elem2D%sPEnum
     rn=partit%com_elem2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, partit%r_mpitype_elem2D(:,1), partit%s_mpitype_elem2D(:,1), &
          partit%com_elem2D%rPE, partit%com_elem2D%sPE)
#endif

     if(lg2g) then
#ifdef ENABLE_OPENACC
        call g2g_cols_begin(2, 1, size(elem_array2D), partit, elem_array2D)
        return
#endif
       !$ACC HOST_DATA USE_DEVICE(elem_array2D)
       DO n=1,rn
          call MPI_IRECV(elem_array2D, 1, partit%r_mpitype_elem2D(n,1), partit%com_elem2D%rPE(n), &
                         partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D%req(n), partit%MPIerr)
       END DO
       DO n=1, sn
          call MPI_ISEND(elem_array2D, 1, partit%s_mpitype_elem2D(n,1), partit%com_elem2D%sPE(n), &
                         partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D%req(rn+n), partit%MPIerr)
       END DO
       !$ACC END HOST_DATA
     else
       DO n=1,rn
          call MPI_IRECV(elem_array2D, 1, partit%r_mpitype_elem2D(n,1), partit%com_elem2D%rPE(n), &
                         partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D%req(n), partit%MPIerr)
       END DO
       DO n=1, sn
          call MPI_ISEND(elem_array2D, 1, partit%s_mpitype_elem2D(n,1), partit%com_elem2D%sPE(n), &
                         partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D%req(rn+n), partit%MPIerr)
       END DO
     endif
     partit%com_elem2D%nreq = rn+sn

  else
     partit%elem_full_flag = .true.

     sn=partit%com_elem2D_full%sPEnum
     rn=partit%com_elem2D_full%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, partit%r_mpitype_elem2D_full(:,1), partit%s_mpitype_elem2D_full(:,1), &
          partit%com_elem2D_full%rPE, partit%com_elem2D_full%sPE)
#endif

     if(lg2g) then
       !$ACC HOST_DATA USE_DEVICE(elem_array2D)
       DO n=1,rn
          call MPI_IRECV(elem_array2D, 1, partit%r_mpitype_elem2D_full(n,1), partit%com_elem2D_full%rPE(n), &
                         partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(n), partit%MPIerr)
       END DO
       DO n=1, sn
          call MPI_ISEND(elem_array2D, 1, partit%s_mpitype_elem2D_full(n,1), partit%com_elem2D_full%sPE(n), &
                         partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(rn+n), partit%MPIerr)
       END DO
       !$ACC END HOST_DATA
     else
       DO n=1,rn
          call MPI_IRECV(elem_array2D, 1, partit%r_mpitype_elem2D_full(n,1), partit%com_elem2D_full%rPE(n), &
                         partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(n), partit%MPIerr)
       END DO
       DO n=1, sn
          call MPI_ISEND(elem_array2D, 1, partit%s_mpitype_elem2D_full(n,1), partit%com_elem2D_full%sPE(n), &
                         partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(rn+n), partit%MPIerr)
       END DO
     endif

     partit%com_elem2D_full%nreq = rn+sn

  end if

end if

END SUBROUTINE exchange_elem2D_begin
! ========================================================================
!Exchange with ALL(!) the neighbours
subroutine exchange_elem2D_i(elem_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer,        intent(inout)       :: elem_array2D(:)
integer                             :: n, sn, rn
logical,        intent(in),optional    :: luse_g2g

 if (partit%npes> 1) then
    call exchange_elem2D_i_begin(elem_array2D, partit, luse_g2g)
    call exchange_elem_end(partit)
end if

END SUBROUTINE exchange_elem2D_i
!=============================================================================
!Exchange with ALL(!) the neighbours
subroutine exchange_elem2D_i_begin(elem_array2D, partit, luse_g2g)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer,        intent(inout)       :: elem_array2D(:)
integer                             :: n, sn, rn
logical                                :: lg2g
logical,        intent(in),optional    :: luse_g2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

 if (partit%npes> 1) then

    partit%elem_full_flag = .true.

    sn=partit%com_elem2D_full%sPEnum
    rn=partit%com_elem2D_full%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
    call check_mpi_comm(rn, sn, partit%r_mpitype_elem2D_full_i, partit%s_mpitype_elem2D_full_i,        &
          partit%com_elem2D_full%rPE, partit%com_elem2D_full%sPE)
#endif
    if(lg2g) then
      !$ACC HOST_DATA USE_DEVICE(elem_array2D)
      DO n=1,rn
         call MPI_IRECV(elem_array2D, 1, partit%r_mpitype_elem2D_full_i(n), partit%com_elem2D_full%rPE(n), &
              partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(n), partit%MPIerr)
      END DO
  
      DO n=1, sn
  
         call MPI_ISEND(elem_array2D, 1, partit%s_mpitype_elem2D_full_i(n), partit%com_elem2D_full%sPE(n), &
              partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(rn+n), partit%MPIerr)
      END DO
      !$ACC END HOST_DATA
    else
      DO n=1,rn
         call MPI_IRECV(elem_array2D, 1, partit%r_mpitype_elem2D_full_i(n), partit%com_elem2D_full%rPE(n), &
              partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(n), partit%MPIerr)
      END DO
  
      DO n=1, sn
  
         call MPI_ISEND(elem_array2D, 1, partit%s_mpitype_elem2D_full_i(n), partit%com_elem2D_full%sPE(n), &
              partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(rn+n), partit%MPIerr)
      END DO
    endif

    partit%com_elem2D_full%nreq = rn+sn

end if

END SUBROUTINE exchange_elem2D_i_begin
! ========================================================================
! Broadcast routines
! Many because of different sizes.
! ========================================================================
subroutine broadcast_nod3D(arr3D, arr3Dglobal, partit)
! Distribute the nodal information available on 0 PE to other PEs
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
INTEGER                             :: nz, counter,nl1
integer                             ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
real(kind=WP)                        ::  arr3D(:,:)
real(kind=WP)                        ::  arr3Dglobal(:,:)
integer                             :: node_size
INTEGER,      ALLOCATABLE, DIMENSION(:) ::  irecvbuf
real(kind=WP), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf

node_size=partit%myDim_nod2D+partit%eDim_nod2D
nl1=ubound(arr3D,1)
IF ( partit%mype == 0 ) THEN
    if (partit%npes>1) then
    arr3D(:,1:node_size)=arr3Dglobal(:,partit%myList_nod2D(1:node_size))
    end if
    DO  n = 1, partit%npes-1
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, partit%MPI_COMM_FESOM, status, partit%MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(sendbuf(nTS*nl1), irecvbuf(nTS))

       CALL MPI_RECV(irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, partit%MPI_COMM_FESOM, status, partit%MPIerr )
       counter=0
       DO i = 1, nTS
          DO nz=1, nl1
             counter=counter+1
             sendbuf(counter) = arr3Dglobal(nz,irecvbuf(i))
          ENDDO
       ENDDO

       CALL MPI_SEND(sendbuf(1), nTS*nl1, MPI_WP, &
                   sender, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

       DEALLOCATE(irecvbuf, sendbuf)
    ENDDO
ELSE
    CALL MPI_SEND( node_size, 1, MPI_INTEGER, 0, 0, partit%MPI_COMM_FESOM, partit%MPIerr )
    CALL MPI_SEND( partit%myList_nod2D(1), node_size, MPI_INTEGER, 0, 1, &
                   partit%MPI_COMM_FESOM, partit%MPIerr )

    ALLOCATE(recvbuf(node_size*nl1))
    CALL MPI_RECV( recvbuf(1), node_size*nl1, MPI_WP, 0, &
                      2, partit%MPI_COMM_FESOM, status, partit%MPIerr )
    counter=0
    DO n = 1, node_size
       DO nz=1, nl1
         counter=counter+1
         arr3D(nz,n)=recvbuf(counter)
       ENDDO
    ENDDO

    DEALLOCATE(recvbuf)
ENDIF
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)
end subroutine broadcast_nod3D
!
!============================================================================
!
subroutine broadcast_nod2D(arr2D, arr2Dglobal, partit)
! A 2D version of the previous routine
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(in), target       :: partit
integer :: MPIerr
real(kind=WP)                             :: arr2D(:)
real(kind=WP)                             :: arr2Dglobal(:)
integer                                  :: i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER,      ALLOCATABLE, DIMENSION(:)  :: irecvbuf
real(kind=WP), ALLOCATABLE, DIMENSION(:)  :: sendbuf
integer                                  :: node_size

node_size=partit%myDim_nod2D+partit%eDim_nod2D

IF ( partit%mype == 0 ) THEN
    if (partit%npes>1) then
    arr2D(1:node_size)=arr2Dglobal(partit%myList_nod2D(1:node_size))
    end if
    DO  n = 1, partit%npes-1
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, partit%MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(sendbuf(nTS), irecvbuf(nTS))

       CALL MPI_RECV(irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, partit%MPI_COMM_FESOM, status, MPIerr )
       DO i = 1, nTS
             sendbuf(i) = arr2Dglobal(irecvbuf(i))
       ENDDO

       CALL MPI_SEND(sendbuf(1), nTS, MPI_WP, &
                   sender, 2, partit%MPI_COMM_FESOM, MPIerr )

       DEALLOCATE(irecvbuf, sendbuf)
    ENDDO
ELSE
    CALL MPI_SEND( node_size, 1, MPI_INTEGER, 0, 0, partit%MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( partit%myList_nod2D(1), node_size, MPI_INTEGER, 0, 1, &
                   partit%MPI_COMM_FESOM, MPIerr )
    CALL MPI_RECV( arr2D(1), node_size, MPI_WP, 0, &
                      2, partit%MPI_COMM_FESOM, status, MPIerr )
ENDIF
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,MPIerr)
end subroutine broadcast_nod2D
!
!============================================================================
!
subroutine broadcast_elem3D(arr3D, arr3Dglobal, partit)
! Distribute the elemental information available on 0 PE to other PEs
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(in), target      :: partit
integer :: MPIerr
INTEGER                                 :: nz, counter,nl1
integer                                 :: i, n, nTS, sender, status(MPI_STATUS_SIZE)
real(kind=WP)                            :: arr3D(:,:)
real(kind=WP)                            :: arr3Dglobal(:,:)
integer                                 :: elem_size

INTEGER,      ALLOCATABLE, DIMENSION(:) ::  irecvbuf
real(kind=WP), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf

elem_size=partit%myDim_elem2D+partit%eDim_elem2D

nl1=ubound(arr3D,1)
IF ( partit%mype == 0 ) THEN
    if (partit%npes>1) then
    arr3D(:,1:elem_size)=arr3Dglobal(:,partit%myList_elem2D(1:elem_size))
    end if
    DO  n = 1, partit%npes-1
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, partit%MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(sendbuf(nTS*nl1), irecvbuf(nTS))

       CALL MPI_RECV(irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, partit%MPI_COMM_FESOM, status, MPIerr )
       counter=0
       DO i = 1, nTS
          DO nz=1, nl1
             counter=counter+1
             sendbuf(counter) = arr3Dglobal(nz,irecvbuf(i))
          ENDDO
       ENDDO

       CALL MPI_SEND(sendbuf(1), nTS*nl1, MPI_WP, &
                   sender, 2, partit%MPI_COMM_FESOM, MPIerr )

       DEALLOCATE(irecvbuf, sendbuf)
    ENDDO
ELSE
    CALL MPI_SEND( elem_size, 1, MPI_INTEGER, 0, 0, partit%MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( partit%myList_elem2D(1), elem_size, MPI_INTEGER, 0, 1, &
                   partit%MPI_COMM_FESOM, MPIerr )

    ALLOCATE(recvbuf(elem_size*nl1))
    CALL MPI_RECV( recvbuf(1), elem_size*nl1, MPI_WP, 0, &
                      2, partit%MPI_COMM_FESOM, status, MPIerr )
    counter=0
    DO n = 1, elem_size
       DO nz=1, nl1
       counter=counter+1
       arr3D(nz,n)=recvbuf(counter)
       ENDDO
    ENDDO

    DEALLOCATE(recvbuf)
ENDIF
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,MPIerr)
end subroutine broadcast_elem3D
!
!============================================================================
!
subroutine broadcast_elem2D(arr2D, arr2Dglobal, partit)
! A 2D version of the previous routine
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(in), target      :: partit
integer :: MPIerr
integer                                 :: i, n, nTS, sender, status(MPI_STATUS_SIZE)
real(kind=WP)                            :: arr2D(:)
real(kind=WP)                            :: arr2Dglobal(:)
integer                                 :: elem_size
INTEGER,      ALLOCATABLE, DIMENSION(:) :: irecvbuf
real(kind=WP), ALLOCATABLE, DIMENSION(:) :: sendbuf

elem_size=partit%myDim_elem2D+partit%eDim_elem2D

IF ( partit%mype == 0 ) THEN
    if (partit%npes>1) then
    arr2D(1:elem_size)=arr2Dglobal(partit%myList_elem2D(1:elem_size))
    end if
    DO  n = 1, partit%npes-1
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, partit%MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(sendbuf(1:nTS), irecvbuf(nTS))

       CALL MPI_RECV(irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, partit%MPI_COMM_FESOM, status, MPIerr )
       DO i = 1, nTS
             sendbuf(i) = arr2Dglobal(irecvbuf(i))
       ENDDO

       CALL MPI_SEND(sendbuf(1), nTS, MPI_WP, &
                   sender, 2, partit%MPI_COMM_FESOM, MPIerr )

       DEALLOCATE(irecvbuf, sendbuf)
    ENDDO
ELSE
    CALL MPI_SEND( elem_size, 1, MPI_INTEGER, 0, 0, partit%MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( partit%myList_elem2D(1), elem_size, MPI_INTEGER, 0, 1, &
                   partit%MPI_COMM_FESOM, MPIerr )
    CALL MPI_RECV( arr2D(1), elem_size, MPI_WP, 0, &
                      2, partit%MPI_COMM_FESOM, status, MPIerr )
ENDIF
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,MPIerr)
end subroutine broadcast_elem2D
!
!============================================================================
! Make nodal information available to master PE
! Use only with 3D arrays stored in (vertical, horizontal) way
subroutine gather_nod3D(arr3D, arr3D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
INTEGER                                :: nl1
integer                                :: n
real(real64)                           :: arr3D(:,:)
real(real64)                           :: arr3D_global(:,:)
real(real64), allocatable              :: recvbuf(:,:)
integer                                :: req(partit%npes-1)
integer                                :: start, n3D

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(nl1,ubound(arr3D_global,2)))

      do  n = 1, partit%npes-1
         n3D = (partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n))*nl1
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_DOUBLE_PRECISION, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr3D_global(1:nl1,partit%myList_nod2D(1:partit%myDim_nod2D)) = arr3D(1:nl1,1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr3D_global(1:nl1, partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1:nl1, 1 : partit%remPtr_nod2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   call MPI_SEND( arr3D, partit%myDim_nod2D*nl1, MPI_DOUBLE_PRECISION, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

end if
end subroutine gather_nod3D
!
!============================================================================
!
subroutine gather_real4_nod3D(arr3D, arr3D_global, partit)

! Make nodal information available to master PE
!
! Use only with 3D arrays stored in (vertical, horizontal) way
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
INTEGER                             :: nl1
integer                             :: n
real(real32)                        :: arr3D(:,:)
real(real32)                        :: arr3D_global(:,:)
real(real32), allocatable           :: recvbuf(:,:)
integer                             :: req(partit%npes-1)
integer                             :: start, n3D

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(nl1,ubound(arr3D_global,2)))

      do  n = 1, partit%npes-1
         n3D = (partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n))*nl1
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr3D_global(1:nl1,partit%myList_nod2D(1:partit%myDim_nod2D)) = arr3D(1:nl1,1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr3D_global(1:nl1, partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1:nl1, 1 : partit%remPtr_nod2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   call MPI_SEND( arr3D, partit%myDim_nod2D*nl1, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

end if
end subroutine gather_real4_nod3D
!=======================================================

subroutine gather_int2_nod3D(arr3D, arr3D_global, partit)

! Make nodal information available to master PE
!
! Use only with 3D arrays stored in (vertical, horizontal) way
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
INTEGER                             :: nl1
integer                             :: n
integer(int16)                      :: arr3D(:,:)
integer(int16)                      :: arr3D_global(:,:)
integer(int16), allocatable         :: recvbuf(:,:)
integer                             :: req(partit%npes-1)
integer                             :: start, n3D


if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(nl1,ubound(arr3D_global,2)))

      do  n = 1, partit%npes-1
         n3D = (partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n))*nl1
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_SHORT, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr3D_global(1:nl1,partit%myList_nod2D(1:partit%myDim_nod2D)) = arr3D(1:nl1,1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr3D_global(1:nl1, partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1:nl1, 1 : partit%remPtr_nod2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   call MPI_SEND( arr3D, partit%myDim_nod2D*nl1, MPI_SHORT, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

end if
end subroutine gather_int2_nod3D
!==============================================
subroutine gather_nod2D(arr2D, arr2D_global, partit)
! Make nodal information available to master PE
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                             :: n
real(real64)                        :: arr2D(:)
real(real64)                        :: arr2D_global(:)
real(real64), allocatable           :: recvbuf(:)
integer                             :: req(partit%npes-1)
integer                             :: start, n2D

if (partit%npes> 1) then

CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, partit%npes-1
         n2D   = partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n)
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_DOUBLE_PRECISION, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(partit%myList_nod2D(1:partit%myDim_nod2D)) = arr2D(1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1 : partit%remPtr_nod2D(partit%npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)

   endif

ELSE

   call MPI_SEND( arr2D(1:partit%myDim_nod2D), partit%myDim_nod2D, MPI_DOUBLE_PRECISION, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

endif
end subroutine gather_nod2D
!==============================================
subroutine gather_real4_nod2D(arr2D, arr2D_global, partit)
! Make nodal information available to master PE
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                             :: n
real(real32)                        :: arr2D(:)
real(real32)                        :: arr2D_global(:)
real(real32),   allocatable         :: recvbuf(:)
integer                             :: req(partit%npes-1)
integer                             :: start, n2D

if (partit%npes> 1) then

CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, partit%npes-1
         n2D   = partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n)
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(partit%myList_nod2D(1:partit%myDim_nod2D)) = arr2D(1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1 : partit%remPtr_nod2D(partit%npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)

   endif

ELSE

   call MPI_SEND( arr2D(1:partit%myDim_nod2D), partit%myDim_nod2D, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

endif
end subroutine gather_real4_nod2D

!==============================================
! Make nodal information available to master PE
subroutine gather_int2_nod2D(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                             :: n
integer(int16)                      :: arr2D(:)
integer(int16)                      :: arr2D_global(:)
integer(int16), allocatable         :: recvbuf(:)
integer                             :: req(partit%npes-1)
integer                             :: start, n2D

if (partit%npes> 1) then

CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, partit%npes-1
         n2D   = partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n)
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_SHORT, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(partit%myList_nod2D(1:partit%myDim_nod2D)) = arr2D(1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1 : partit%remPtr_nod2D(partit%npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)

   endif

ELSE

   call MPI_SEND( arr2D(1:partit%myDim_nod2D), partit%myDim_nod2D, MPI_SHORT, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

endif
end subroutine gather_int2_nod2D

!============================================================================
subroutine gather_elem3D(arr3D, arr3D_global, partit)
! Make element information available to master PE
!
! Use only with 3D arrays stored in (vertical, horizontal) way
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
INTEGER                             :: nl1
integer                             :: n
real(real64)                        :: arr3D(:,:)
real(real64)                        :: arr3D_global(:,:)
real(real64),   allocatable         :: recvbuf(:,:)
integer                             :: req(partit%npes-1)
integer                             :: start, e3D, ende, err_alloc
integer                             :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
!
      allocate(recvbuf(nl1,partit%remPtr_elem2D(partit%npes)))

      do  n = 1, partit%npes-1
         e3D = (partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n))*nl1
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_DOUBLE_PRECISION, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr3D_global(1:nl1,partit%myList_elem2D(1:partit%myDim_elem2D)) = arr3D(1:nl1,1:partit%myDim_elem2D)


      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr3D_global(1:nl1, partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                        = recvbuf(1:nl1, 1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   call MPI_SEND( arr3D, partit%myDim_elem2D*nl1, MPI_DOUBLE_PRECISION, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

endif
end subroutine gather_elem3D

!===================================================================
! Make element information available to master PE
! Use only with 3D arrays stored in (vertical, horizontal) way
subroutine gather_real4_elem3D(arr3D, arr3D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
INTEGER                                :: nl1
integer                                :: n
real(real32)                           :: arr3D(:,:)
real(real32)                           :: arr3D_global(:,:)
real(real32),   allocatable            :: recvbuf(:,:)
integer                                :: req(partit%npes-1)
integer                                :: start, e3D, ende, err_alloc
integer                                :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
!
      allocate(recvbuf(nl1,partit%remPtr_elem2D(partit%npes)))

      do  n = 1, partit%npes-1
         e3D = (partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n))*nl1
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr3D_global(1:nl1,partit%myList_elem2D(1:partit%myDim_elem2D)) = arr3D(1:nl1,1:partit%myDim_elem2D)


      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr3D_global(1:nl1, partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                        = recvbuf(1:nl1, 1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   call MPI_SEND( arr3D, partit%myDim_elem2D*nl1, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

endif
end subroutine gather_real4_elem3D


!===================================================================
! Make element information available to master PE
! Use only with 3D arrays stored in (vertical, horizontal) way
subroutine gather_int2_elem3D(arr3D, arr3D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
INTEGER                             :: nl1
integer                             :: n
integer(int16)                      :: arr3D(:,:)
integer(int16)                      :: arr3D_global(:,:)
integer(int16), allocatable         :: recvbuf(:,:)
integer                             :: req(partit%npes-1)
integer                             :: start, e3D, ende, err_alloc
integer                             :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
!
      allocate(recvbuf(nl1,partit%remPtr_elem2D(partit%npes)))

      do  n = 1, partit%npes-1
         e3D = (partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n))*nl1
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_SHORT, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr3D_global(1:nl1,partit%myList_elem2D(1:partit%myDim_elem2D)) = arr3D(1:nl1,1:partit%myDim_elem2D)


      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr3D_global(1:nl1, partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                        = recvbuf(1:nl1, 1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   call MPI_SEND( arr3D, partit%myDim_elem2D*nl1, MPI_SHORT, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

endif
end subroutine gather_int2_elem3D


!==============================================
! Make element information available to master PE
subroutine gather_elem2D(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                             :: n
real(real64)                        :: arr2D(:)
real(real64)                        :: arr2D_global(:)
real(real64), allocatable           :: recvbuf(:)
integer                             :: req(partit%npes-1)
integer                             :: start, e2D

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then

      allocate(recvbuf(partit%remPtr_elem2D(partit%npes)))

      do  n = 1, partit%npes-1
         e2D   = partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n)
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_DOUBLE_PRECISION, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(partit%myList_elem2D(1:partit%myDim_elem2D)) = arr2D(1:partit%myDim_elem2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                       = recvbuf(1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)

   endif

ELSE

   call MPI_SEND( arr2D(1:partit%myDim_elem2D), partit%myDim_elem2D, MPI_DOUBLE_PRECISION, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF
end if

end subroutine gather_elem2D

!================================================
! Make element information available to master PE
subroutine gather_real4_elem2D(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                             :: n
real(real32)                        :: arr2D(:)
real(real32)                        :: arr2D_global(:)
real(real32), allocatable           :: recvbuf(:)
integer                             :: req(partit%npes-1)
integer                             :: start, e2D


if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then

      allocate(recvbuf(partit%remPtr_elem2D(partit%npes)))

      do  n = 1, partit%npes-1
         e2D   = partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n)
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(partit%myList_elem2D(1:partit%myDim_elem2D)) = arr2D(1:partit%myDim_elem2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                       = recvbuf(1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)

   endif

ELSE

   call MPI_SEND( arr2D(1:partit%myDim_elem2D), partit%myDim_elem2D, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF
end if

end subroutine gather_real4_elem2D

!================================================
! Make element information available to master PE
subroutine gather_int2_elem2D(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                                :: n
integer(int16)                         :: arr2D(:)
integer(int16)                         :: arr2D_global(:)
integer(int16), allocatable            :: recvbuf(:)
integer                                :: req(partit%npes-1)
integer                                :: start, e2D

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then

      allocate(recvbuf(partit%remPtr_elem2D(partit%npes)))

      do  n = 1, partit%npes-1
         e2D   = partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n)
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_SHORT, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(partit%myList_elem2D(1:partit%myDim_elem2D)) = arr2D(1:partit%myDim_elem2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                       = recvbuf(1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)

   endif

ELSE

   call MPI_SEND( arr2D(1:partit%myDim_elem2D), partit%myDim_elem2D, MPI_SHORT, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF
end if

end subroutine gather_int2_elem2D


!============================================================================
! Make nodal information available to master PE
! Use only with 3D arrays stored in (vertical, horizontal) way
subroutine gather_real8to4_nod3D(arr3D, arr3D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
INTEGER                             :: nl1
integer                             :: n
real(real64)                        :: arr3D(:,:)
real(real32)                        :: arr3D_global(:,:)
integer                             :: req(partit%npes-1)
integer                             :: start, n3D, ierr
real(real32), allocatable           :: recvbuf(:,:)
real(real32), allocatable           :: sendbuf(:,:)

if (partit%npes> 1) then

CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)
nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(nl1, ubound(arr3D_global,2)))

      do  n = 1, partit%npes-1
         n3D = (partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n))*nl1
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr3D_global(1:nl1,partit%myList_nod2D(1:partit%myDim_nod2D)) = arr3D(1:nl1,1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr3D_global(1:nl1, partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                     = recvbuf(1:nl1, 1 : partit%remPtr_nod2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   allocate(sendbuf(nl1,partit%myDim_nod2D))
   sendbuf(1:nl1,1:partit%myDim_nod2D) = arr3D(1:nl1,1:partit%myDim_nod2D)

   call MPI_SEND(sendbuf, partit%myDim_nod2D*nl1, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )
   deallocate(sendbuf)

ENDIF

end if

end subroutine gather_real8to4_nod3D
!==============================================
! Make nodal information available to master PE
subroutine gather_real8to4_nod2D(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                             :: n
real(real64)                        :: arr2D(:)
real(real32)                        :: arr2D_global(:)
real(real32)                        :: sendbuf(partit%myDim_nod2D)
real(real64), allocatable           :: recvbuf(:)
integer                             :: req(partit%npes-1)
integer                             :: start, n2D

! Consider MPI-datatypes to recv directly into arr2D_global!

 if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)
IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, partit%npes-1
         n2D   = partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n)
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(partit%myList_nod2D(1:partit%myDim_nod2D)) = arr2D(1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1 : partit%remPtr_nod2D(partit%npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)

   endif

ELSE
   sendbuf(1:partit%myDim_nod2D) = real(arr2D(1:partit%myDim_nod2D),real32)

   call MPI_SEND(sendbuf, partit%myDim_nod2D, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

end if
end subroutine gather_real8to4_nod2D
!============================================================================
subroutine gather_real8to4_elem3D(arr3D, arr3D_global, partit)
! Make element information available to master PE
! Use only with 3D arrays stored in (vertical, horizontal) way
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
INTEGER                             :: nl1
integer                             :: n
real(real64)                        :: arr3D(:,:)
real(real32)                        :: arr3D_global(:,:)
integer                             :: req(partit%npes-1)
integer                             :: start, e3D
real(real32), allocatable           :: recvbuf(:,:)
real(real32), allocatable           :: sendbuf(:,:)


if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)
nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(nl1,partit%remPtr_elem2D(partit%npes)))

      do  n = 1, partit%npes-1
         e3D = (partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n))*nl1
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr3D_global(1:nl1,partit%myList_elem2D(1:partit%myDim_elem2D)) = arr3D(1:nl1,1:partit%myDim_elem2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr3D_global(1:nl1, partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                     = recvbuf(1:nl1, 1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   allocate(sendbuf(nl1,partit%myDim_elem2D))
   sendbuf(1:nl1,1:partit%myDim_elem2D) = arr3D(1:nl1,1:partit%myDim_elem2D)

   call MPI_SEND(sendbuf, partit%myDim_elem2D*nl1, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )
   deallocate(sendbuf)
ENDIF

end if
end subroutine gather_real8to4_elem3D
!================================================
! Make element information available to master PE
subroutine gather_real8to4_elem2D(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                             :: n
real(real64)                        :: arr2D(:)
real(real32)                        :: arr2D_global(:)
real(real32), allocatable           :: recvbuf(:)
real(real32)                        :: sendbuf(partit%myDim_elem2D)
integer                             :: req(partit%npes-1)
integer                             :: start, e2D

if (partit%npes> 1) then

CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)
! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(partit%remPtr_elem2D(partit%npes)))

      do  n = 1, partit%npes-1
         e2D   = partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n)
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(partit%myList_elem2D(1:partit%myDim_elem2D)) = arr2D(1:partit%myDim_elem2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                       = recvbuf(1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)

   endif

ELSE

   sendbuf(1:partit%myDim_elem2D) = real(arr2D(1:partit%myDim_elem2D),real32)
   call MPI_SEND(sendbuf, partit%myDim_elem2D, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

end if
end subroutine gather_real8to4_elem2D
!==============================================
subroutine gather_elem2D_i(arr2D, arr2D_global, partit)
! Make element information available to master PE
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
integer                       :: n
integer                       :: arr2D(:)
integer                       :: arr2D_global(:)
integer, allocatable          :: recvbuf(:)
integer                       :: req(partit%npes-1)
integer                       :: start, e2D


CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)
  ! Consider MPI-datatypes to recv directly into arr2D_global!
  IF ( partit%mype == 0 ) THEN
     if (partit%npes > 1) then
        allocate(recvbuf(partit%remPtr_elem2D(partit%npes)))
        do  n = 1, partit%npes-1
            e2D   = partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n)
            start = partit%remPtr_elem2D(n)
            call MPI_IRECV(recvbuf(start), e2D, MPI_INTEGER, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
        enddo
        arr2D_global(partit%myList_elem2D(1:partit%myDim_elem2D)) = arr2D(1:partit%myDim_elem2D)
        call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)
        arr2D_global(partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                          = recvbuf(1 : partit%remPtr_elem2D(partit%npes)-1)
        deallocate(recvbuf)
     else
        arr2D_global(:) = arr2D(:)
     endif
  ELSE
     ! Modify arr2D to be an explicitly contiguous array with arr2D --> arr2D(1:myDim_elem2D)
     call MPI_SEND(arr2D(1:partit%myDim_elem2D), partit%myDim_elem2D, MPI_INTEGER, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )
  ENDIF
end subroutine gather_elem2D_i
!==============================================
! Make nodal information available to master PE
subroutine gather_nod2D_i(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
integer              :: n
integer              :: arr2D(:)
integer              :: arr2D_global(:)
integer, allocatable :: recvbuf(:)
integer              :: req(partit%npes-1)
integer              :: start, n2D

if (partit%npes> 1) then

CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(ubound(arr2D_global, 1)))
      do  n = 1, partit%npes-1
         n2D   = partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n)
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_INTEGER, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(partit%myList_nod2D(1:partit%myDim_nod2D)) = arr2D(1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1 : partit%remPtr_nod2D(partit%npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)

   endif

ELSE

   call MPI_SEND( arr2D(1:partit%myDim_nod2D), partit%myDim_nod2D, MPI_INTEGER, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

endif
end subroutine gather_nod2D_i
!============================================================================
! A 2D version of the previous routine
subroutine gather_edg2D(arr2D, arr2Dglobal, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(in), target      :: partit
integer :: MPIerr
real(kind=WP)                            :: arr2D(:)
real(kind=WP)                            :: arr2Dglobal(:)
integer                                 ::  i, n, buf_size, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:)      ::  ibuf
real(kind=WP), ALLOCATABLE, DIMENSION(:) ::  rbuf

IF ( partit%mype == 0 ) THEN
    arr2Dglobal(partit%myList_edge2D(1:partit%myDim_edge2D))=arr2D(1:partit%myDim_edge2D)
    DO  n = 1, partit%npes-1
       CALL MPI_RECV( buf_size, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                      0, partit%MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(rbuf(buf_size), ibuf(buf_size))

!PS        CALL MPI_RECV(ibuf(1), buf_size, MPI_INTEGER, sender, &
!PS                       1, MPI_COMM_FESOM, status, MPIerr )
!PS        CALL MPI_RECV(rbuf(1), buf_size, MPI_WP, sender, &
!PS                       2, MPI_COMM_FESOM, status, MPIerr )
                      
       CALL MPI_RECV(ibuf, buf_size, MPI_INTEGER, sender, &
                      1, partit%MPI_COMM_FESOM, status, MPIerr )
       CALL MPI_RECV(rbuf, buf_size, MPI_WP, sender, &
                      2, partit%MPI_COMM_FESOM, status, MPIerr )               
       arr2Dglobal(ibuf)=rbuf
       DEALLOCATE(ibuf, rbuf)
    ENDDO
ELSE
!PS     CALL MPI_SEND( myDim_edge2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
!PS     CALL MPI_SEND( myList_edge2D(1), myDim_edge2D, MPI_INTEGER, 0, 1, &
!PS                    MPI_COMM_FESOM, MPIerr )
!PS     CALL MPI_SEND( arr2D(1), myDim_edge2D, MPI_WP, 0, 2,&
!PS                    MPI_COMM_FESOM, MPIerr )
                   
    CALL MPI_SEND( partit%myDim_edge2D, 1, MPI_INTEGER, 0, 0, partit%MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( partit%myList_edge2D(1:partit%myDim_edge2D), partit%myDim_edge2D, MPI_INTEGER, 0, 1, &
                   partit%MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( arr2D(1:partit%myDim_edge2D), partit%myDim_edge2D, MPI_WP, 0, 2,&
                   partit%MPI_COMM_FESOM, MPIerr )
                   
ENDIF
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,MPIerr)
end subroutine gather_edg2D
!
!============================================================================
! A 2D version of the previous routine
subroutine gather_edg2D_i(arr2D, arr2Dglobal, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                                :: arr2D(:)
integer                                :: arr2Dglobal(:)
integer                                :: i, n, buf_size, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:)     :: ibuf, vbuf

IF ( partit%mype == 0 ) THEN
    arr2Dglobal(partit%myList_edge2D(1:partit%myDim_edge2D))=arr2D(1:partit%myDim_edge2D)
    DO  n = 1, partit%npes-1
       CALL MPI_RECV( buf_size, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                      0, partit%MPI_COMM_FESOM, status, partit%MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(ibuf(buf_size), vbuf(buf_size))

!PS        CALL MPI_RECV(ibuf(1), buf_size, MPI_INTEGER, sender, &
!PS                       1, MPI_COMM_FESOM, status, MPIerr )
!PS        CALL MPI_RECV(vbuf(1), buf_size, MPI_INTEGER, sender, &
!PS                       2, MPI_COMM_FESOM, status, MPIerr )

       CALL MPI_RECV(ibuf, buf_size, MPI_INTEGER, sender, &
                      1, partit%MPI_COMM_FESOM, status, partit%MPIerr )
       CALL MPI_RECV(vbuf, buf_size, MPI_INTEGER, sender, &
                      2, partit%MPI_COMM_FESOM, status, partit%MPIerr )               
       arr2Dglobal(ibuf)=vbuf
       DEALLOCATE(ibuf, vbuf)
    ENDDO
ELSE
!PS     CALL MPI_SEND( myDim_edge2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
!PS     CALL MPI_SEND( myList_edge2D(1), myDim_edge2D, MPI_INTEGER, 0, 1, &
!PS                    MPI_COMM_FESOM, MPIerr )
!PS     CALL MPI_SEND( arr2D(1), myDim_edge2D, MPI_INTEGER, 0, 2,&
!PS                    MPI_COMM_FESOM, MPIerr )
                   
    CALL MPI_SEND( partit%myDim_edge2D, 1, MPI_INTEGER, 0, 0, partit%MPI_COMM_FESOM, partit%MPIerr )
    CALL MPI_SEND( partit%myList_edge2D(1:partit%myDim_edge2D), partit%myDim_edge2D, MPI_INTEGER, 0, 1, &
                   partit%MPI_COMM_FESOM, partit%MPIerr )
    CALL MPI_SEND( arr2D(1:partit%myDim_edge2D), partit%myDim_edge2D, MPI_INTEGER, 0, 2,&
                   partit%MPI_COMM_FESOM, partit%MPIerr )
                   
                   
ENDIF
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)
end subroutine gather_edg2D_i
!==============================================

! ========================================================================
! Gather routines for spectral bin (fbin) fields used by IDEMIX2.
! These handle arrays with dimensions (nfbin, horizontal) instead of (nl, horizontal).
! ========================================================================

!============================================================================
! Make nodal spectral bin information available to master PE
! arr3D(nfbin, nod2D)
subroutine gather_nod2D_fbin(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                                :: nfbin
integer                                :: n
real(real64)                           :: arr2D(:,:)
real(real64)                           :: arr2D_global(:,:)
real(real64), allocatable              :: recvbuf(:,:)
integer                                :: req(partit%npes-1)
integer                                :: start, n2D

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nfbin=ubound(arr2D,1)

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(nfbin,ubound(arr2D_global,2)))

      do  n = 1, partit%npes-1
         n2D = (partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n))*nfbin
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n2D, MPI_DOUBLE_PRECISION, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(1:nfbin,partit%myList_nod2D(1:partit%myDim_nod2D)) = arr2D(1:nfbin,1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(1:nfbin, partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1:nfbin, 1 : partit%remPtr_nod2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr2D_global(:,:) = arr2D(:,:)
   endif

ELSE

   call MPI_SEND( arr2D, partit%myDim_nod2D*nfbin, MPI_DOUBLE_PRECISION, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

end if
end subroutine gather_nod2D_fbin

!============================================================================
! Make element spectral bin information available to master PE
! arr2D(nfbin, elem2D)
subroutine gather_elem2D_fbin(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                                :: nfbin
integer                                :: n
real(real64)                           :: arr2D(:,:)
real(real64)                           :: arr2D_global(:,:)
real(real64), allocatable              :: recvbuf(:,:)
integer                                :: req(partit%npes-1)
integer                                :: start, n2D

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nfbin=ubound(arr2D,1)

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(nfbin,ubound(arr2D_global,2)))

      do  n = 1, partit%npes-1
         n2D = (partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n))*nfbin
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), n2D, MPI_DOUBLE_PRECISION, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(1:nfbin,partit%myList_elem2D(1:partit%myDim_elem2D)) = arr2D(1:nfbin,1:partit%myDim_elem2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(1:nfbin, partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                       = recvbuf(1:nfbin, 1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr2D_global(:,:) = arr2D(:,:)
   endif

ELSE

   call MPI_SEND( arr2D, partit%myDim_elem2D*nfbin, MPI_DOUBLE_PRECISION, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

end if
end subroutine gather_elem2D_fbin

!============================================================================
! Make nodal spectral bin information available to master PE (real4 version)
! arr2D(nfbin, nod2D)
subroutine gather_real4_nod2D_fbin(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                                :: nfbin
integer                                :: n
real(real32)                           :: arr2D(:,:)
real(real32)                           :: arr2D_global(:,:)
real(real32), allocatable              :: recvbuf(:,:)
integer                                :: req(partit%npes-1)
integer                                :: start, n2D

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nfbin=ubound(arr2D,1)

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(nfbin,ubound(arr2D_global,2)))

      do  n = 1, partit%npes-1
         n2D = (partit%remPtr_nod2D(n+1) - partit%remPtr_nod2D(n))*nfbin
         start = partit%remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n2D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(1:nfbin,partit%myList_nod2D(1:partit%myDim_nod2D)) = arr2D(1:nfbin,1:partit%myDim_nod2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(1:nfbin, partit%remList_nod2D(1 : partit%remPtr_nod2D(partit%npes)-1)) &
                       = recvbuf(1:nfbin, 1 : partit%remPtr_nod2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr2D_global(:,:) = arr2D(:,:)
   endif

ELSE

   call MPI_SEND( arr2D, partit%myDim_nod2D*nfbin, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

end if
end subroutine gather_real4_nod2D_fbin

!============================================================================
! Make element spectral bin information available to master PE (real4 version)
! arr2D(nfbin, elem2D)
subroutine gather_real4_elem2D_fbin(arr2D, arr2D_global, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
integer                                :: nfbin
integer                                :: n
real(real32)                           :: arr2D(:,:)
real(real32)                           :: arr2D_global(:,:)
real(real32), allocatable              :: recvbuf(:,:)
integer                                :: req(partit%npes-1)
integer                                :: start, n2D

if (partit%npes> 1) then
CALL MPI_BARRIER(partit%MPI_COMM_FESOM,partit%MPIerr)

nfbin=ubound(arr2D,1)

IF ( partit%mype == 0 ) THEN

   if (partit%npes>1) then
      allocate(recvbuf(nfbin,ubound(arr2D_global,2)))

      do  n = 1, partit%npes-1
         n2D = (partit%remPtr_elem2D(n+1) - partit%remPtr_elem2D(n))*nfbin
         start = partit%remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), n2D, MPI_REAL, n, 2, partit%MPI_COMM_FESOM, req(n), partit%MPIerr)
      enddo

      arr2D_global(1:nfbin,partit%myList_elem2D(1:partit%myDim_elem2D)) = arr2D(1:nfbin,1:partit%myDim_elem2D)

      call MPI_WAITALL(partit%npes-1, req, MPI_STATUSES_IGNORE, partit%MPIerr)

      arr2D_global(1:nfbin, partit%remList_elem2D(1 : partit%remPtr_elem2D(partit%npes)-1)) &
                       = recvbuf(1:nfbin, 1 : partit%remPtr_elem2D(partit%npes)-1)

      deallocate(recvbuf)

   else
      arr2D_global(:,:) = arr2D(:,:)
   endif

ELSE

   call MPI_SEND( arr2D, partit%myDim_elem2D*nfbin, MPI_REAL, 0, 2, partit%MPI_COMM_FESOM, partit%MPIerr )

ENDIF

end if
end subroutine gather_real4_elem2D_fbin

!==============================================

! ========================================================================
! Halo exchange routines for spectral bin (fbin) fields used by IDEMIX2.
! The first array dimension is nfbin (number of spectral frequency bins)
! instead of nl (number of vertical levels).
! Naming convention:
!   exchange_nod2D_fbin:  nodal 2D spectral field  arr(nfbin, nod2D)
!   exchange_nod3D_fbin:  nodal 3D spectral field  arr(n_val, nfbin, nod2D)
!   exchange_elem2D_fbin: elem  2D spectral field  arr(nfbin, elem2D)
!   exchange_elem3D_fbin: elem  3D spectral field  arr(n_val, nfbin, elem2D)
! ========================================================================

! ========================================================================
! 2D nodal spectral field:  arr(nfbin, myDim_nod2D+eDim_nod2D)
! ========================================================================
subroutine exchange_nod2D_fbin(nod_array, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod_array(:,:)
if (partit%npes > 1) then
   call exchange_nod2D_fbin_begin(nod_array, partit)
   call exchange_nod_end(partit)
endif
END SUBROUTINE exchange_nod2D_fbin

subroutine exchange_nod2D_fbin_begin(nod_array, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod_array(:,:)
integer                               :: n, sn, rn

if (partit%npes > 1) then
   if (partit%nfbin_mpi == 0) then
      if (partit%mype==0) print *,'ERROR: exchange_nod2D_fbin called but init_mpi_types_fbin was not called.'
      call par_abort(partit%MPI_COMM_FESOM, partit%mype)
   endif

   sn=partit%com_nod2D%sPEnum
   rn=partit%com_nod2D%rPEnum

   ! Exchange full spectral column at once (n_val=1 for 2D array)
   DO n=1,rn
      call MPI_IRECV(nod_array, 1, partit%r_mpitype_nod3D_fbin(n,partit%nfbin_mpi,1), partit%com_nod2D%rPE(n), &
           partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
   END DO
   DO n=1, sn
      call MPI_ISEND(nod_array, 1, partit%s_mpitype_nod3D_fbin(n,partit%nfbin_mpi,1), partit%com_nod2D%sPE(n), &
           partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
   END DO
   partit%com_nod2D%nreq = rn+sn
   call MPI_WAITALL(partit%com_nod2D%nreq, partit%com_nod2D%req, MPI_STATUSES_IGNORE, partit%MPIerr)
endif
END SUBROUTINE exchange_nod2D_fbin_begin

! ========================================================================
! 3D nodal spectral field:  arr(n_val, nfbin, myDim_nod2D+eDim_nod2D)
! ========================================================================
subroutine exchange_nod3D_fbin(nod_array, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod_array(:,:,:)
if (partit%npes > 1) then
   call exchange_nod3D_fbin_begin(nod_array, partit)
   call exchange_nod_end(partit)
endif
END SUBROUTINE exchange_nod3D_fbin

subroutine exchange_nod3D_fbin_begin(nod_array, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target :: partit
real(kind=WP),   intent(inout)         :: nod_array(:,:,:)
integer                               :: n, sn, rn, n_val

if (partit%npes > 1) then
   if (partit%nfbin_mpi == 0) then
      if (partit%mype==0) print *,'ERROR: exchange_nod3D_fbin called but init_mpi_types_fbin was not called.'
      call par_abort(partit%MPI_COMM_FESOM, partit%mype)
   endif

   ! nod_array(n_val, nfbin, nod2D)
   n_val = ubound(nod_array,1)
   if (n_val < 1 .or. n_val > 3) then
      if (partit%mype==0) print *,'exchange_nod3D_fbin: n_val must be 1..3, got', n_val
      call par_abort(partit%MPI_COMM_FESOM, partit%mype)
   endif

   sn=partit%com_nod2D%sPEnum
   rn=partit%com_nod2D%rPEnum

   ! Exchange full spectral column at once
   DO n=1,rn
      call MPI_IRECV(nod_array, 1, partit%r_mpitype_nod3D_fbin(n,partit%nfbin_mpi,n_val), partit%com_nod2D%rPE(n), &
           partit%com_nod2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_nod2D%req(n), partit%MPIerr)
   END DO
   DO n=1, sn
      call MPI_ISEND(nod_array, 1, partit%s_mpitype_nod3D_fbin(n,partit%nfbin_mpi,n_val), partit%com_nod2D%sPE(n), &
           partit%mype, partit%MPI_COMM_FESOM, partit%com_nod2D%req(rn+n), partit%MPIerr)
   END DO
   partit%com_nod2D%nreq = rn+sn
   call MPI_WAITALL(partit%com_nod2D%nreq, partit%com_nod2D%req, MPI_STATUSES_IGNORE, partit%MPIerr)
endif
END SUBROUTINE exchange_nod3D_fbin_begin

! ========================================================================
! 2D elemental spectral field:  arr(nfbin, myDim_elem2D+eDim_elem2D)
! ========================================================================
subroutine exchange_elem2D_fbin(elem_array, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)          :: elem_array(:,:)

if (partit%npes > 1) then
   call exchange_elem2D_fbin_begin(elem_array, partit)
   call exchange_elem_end(partit)
endif
END SUBROUTINE exchange_elem2D_fbin

subroutine exchange_elem2D_fbin_begin(elem_array, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)          :: elem_array(:,:)
integer                                :: n, sn, rn

if (partit%npes > 1) then
   if (partit%nfbin_mpi == 0) then
      if (partit%mype==0) print *,'ERROR: exchange_elem2D_fbin called but init_mpi_types_fbin was not called.'
      call par_abort(partit%MPI_COMM_FESOM, partit%mype)
   endif

   if (ubound(elem_array,2)<=partit%myDim_elem2D+partit%eDim_elem2D) then
      ! small halo
      partit%elem_full_flag = .false.
      sn=partit%com_elem2D%sPEnum
      rn=partit%com_elem2D%rPEnum

      ! Exchange full spectral column at once
      DO n=1,rn
         call MPI_IRECV(elem_array, 1, partit%r_mpitype_elem2D_fbin(n,partit%nfbin_mpi), partit%com_elem2D%rPE(n), &
              partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D%req(n), partit%MPIerr)
      END DO
      DO n=1, sn
         call MPI_ISEND(elem_array, 1, partit%s_mpitype_elem2D_fbin(n,partit%nfbin_mpi), partit%com_elem2D%sPE(n), &
              partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D%req(rn+n), partit%MPIerr)
      END DO
      partit%com_elem2D%nreq = rn+sn
      call MPI_WAITALL(partit%com_elem2D%nreq, partit%com_elem2D%req, MPI_STATUSES_IGNORE, partit%MPIerr)
   else
      ! full halo
      partit%elem_full_flag = .true.
      sn=partit%com_elem2D_full%sPEnum
      rn=partit%com_elem2D_full%rPEnum

      ! Exchange full spectral column at once
      DO n=1,rn
         call MPI_IRECV(elem_array, 1, partit%r_mpitype_elem2D_full_fbin(n,partit%nfbin_mpi), partit%com_elem2D_full%rPE(n), &
              partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(n), partit%MPIerr)
      END DO
      DO n=1, sn
         call MPI_ISEND(elem_array, 1, partit%s_mpitype_elem2D_full_fbin(n,partit%nfbin_mpi), partit%com_elem2D_full%sPE(n), &
              partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(rn+n), partit%MPIerr)
      END DO
      partit%com_elem2D_full%nreq = rn+sn
      call MPI_WAITALL(partit%com_elem2D_full%nreq, partit%com_elem2D_full%req, MPI_STATUSES_IGNORE, partit%MPIerr)
   endif
endif
END SUBROUTINE exchange_elem2D_fbin_begin

! ========================================================================
! 3D elemental spectral field:  arr(n_val, nfbin, myDim_elem2D+eDim_elem2D)
! ========================================================================
subroutine exchange_elem3D_fbin(elem_array, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)          :: elem_array(:,:,:)

if (partit%npes > 1) then
   call exchange_elem3D_fbin_begin(elem_array, partit)
   call exchange_elem_end(partit)
endif
END SUBROUTINE exchange_elem3D_fbin

subroutine exchange_elem3D_fbin_begin(elem_array, partit)
use MOD_MESH
USE MOD_PARTIT
use par_support_module, only: par_abort
IMPLICIT NONE
type(t_partit), intent(inout), target  :: partit
real(kind=WP),   intent(inout)          :: elem_array(:,:,:)
integer                                :: n, sn, rn, n_val

if (partit%npes > 1) then
   if (partit%nfbin_mpi == 0) then
      if (partit%mype==0) print *,'ERROR: exchange_elem3D_fbin called but init_mpi_types_fbin was not called.'
      call par_abort(partit%MPI_COMM_FESOM, partit%mype)
   endif

   ! elem_array(n_val, nfbin, elem2D)
   n_val = ubound(elem_array,1)
   if (n_val < 1 .or. n_val > 4) then
      if (partit%mype==0) print *,'exchange_elem3D_fbin: n_val must be 1..4, got', n_val
      call par_abort(partit%MPI_COMM_FESOM, partit%mype)
   endif

   if (ubound(elem_array,3)<=partit%myDim_elem2D+partit%eDim_elem2D) then
      ! small halo
      partit%elem_full_flag = .false.
      sn=partit%com_elem2D%sPEnum
      rn=partit%com_elem2D%rPEnum

      ! Exchange full spectral column at once
      DO n=1,rn
         call MPI_IRECV(elem_array, 1, partit%r_mpitype_elem3D_fbin(n,partit%nfbin_mpi,n_val), partit%com_elem2D%rPE(n), &
              partit%com_elem2D%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D%req(n), partit%MPIerr)
      END DO
      DO n=1, sn
         call MPI_ISEND(elem_array, 1, partit%s_mpitype_elem3D_fbin(n,partit%nfbin_mpi,n_val), partit%com_elem2D%sPE(n), &
              partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D%req(rn+n), partit%MPIerr)
      END DO
      partit%com_elem2D%nreq = rn+sn
      call MPI_WAITALL(partit%com_elem2D%nreq, partit%com_elem2D%req, MPI_STATUSES_IGNORE, partit%MPIerr)
   else
      ! full halo
      partit%elem_full_flag = .true.
      sn=partit%com_elem2D_full%sPEnum
      rn=partit%com_elem2D_full%rPEnum

      ! Exchange full spectral column at once
      DO n=1,rn
         call MPI_IRECV(elem_array, 1, partit%r_mpitype_elem3D_full_fbin(n,partit%nfbin_mpi,n_val), partit%com_elem2D_full%rPE(n), &
              partit%com_elem2D_full%rPE(n), partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(n), partit%MPIerr)
      END DO
      DO n=1, sn
         call MPI_ISEND(elem_array, 1, partit%s_mpitype_elem3D_full_fbin(n,partit%nfbin_mpi,n_val), partit%com_elem2D_full%sPE(n), &
              partit%mype, partit%MPI_COMM_FESOM, partit%com_elem2D_full%req(rn+n), partit%MPIerr)
      END DO
      partit%com_elem2D_full%nreq = rn+sn
      call MPI_WAITALL(partit%com_elem2D_full%nreq, partit%com_elem2D_full%req, MPI_STATUSES_IGNORE, partit%MPIerr)
   endif
endif
END SUBROUTINE exchange_elem3D_fbin_begin

! ========================================================================

end module g_comm



module g_comm_auto
use g_comm
implicit none
interface exchange_nod
      module procedure exchange_nod2D
      module procedure exchange_nod2D_i
      module procedure exchange_nod2D_2fields
      module procedure exchange_nod2D_3fields
      module procedure exchange_nod3D
      module procedure exchange_nod3D_2fields
      module procedure exchange_nod3D_n
end interface exchange_nod

interface exchange_nod_begin
      module procedure exchange_nod2D_begin
      module procedure exchange_nod2D_i_begin
      module procedure exchange_nod2D_2fields_begin
      module procedure exchange_nod2D_3fields_begin
      module procedure exchange_nod3D_begin
      module procedure exchange_nod3D_2fields_begin
      module procedure exchange_nod3D_n_begin
end interface exchange_nod_begin

!!$interface exchange_edge
!!$      module procedure exchange_edge2D
!!$!      module procedure exchange_edge3D  ! not available, not used
!!$end interface exchange_edge

interface exchange_elem
      module procedure exchange_elem3D
      module procedure exchange_elem3D_n
      module procedure exchange_elem2d
      module procedure exchange_elem2d_i
end interface exchange_elem

interface exchange_elem_begin
      module procedure exchange_elem3D_begin
      module procedure exchange_elem3D_n_begin
      module procedure exchange_elem2d_begin
      module procedure exchange_elem2d_i_begin
end interface exchange_elem_begin

interface exchange_nod_fbin
      module procedure exchange_nod2D_fbin
      module procedure exchange_nod3D_fbin
end interface exchange_nod_fbin

interface exchange_nod_fbin_begin
      module procedure exchange_nod2D_fbin_begin
      module procedure exchange_nod3D_fbin_begin
end interface exchange_nod_fbin_begin

interface exchange_elem_fbin
      module procedure exchange_elem2D_fbin
      module procedure exchange_elem3D_fbin
end interface exchange_elem_fbin

interface exchange_elem_fbin_begin
      module procedure exchange_elem2D_fbin_begin
      module procedure exchange_elem3D_fbin_begin
end interface exchange_elem_fbin_begin


interface broadcast_nod
      module procedure broadcast_nod3D
      module procedure broadcast_nod2D
end interface broadcast_nod

interface broadcast_elem
      module procedure broadcast_elem3D
      module procedure broadcast_elem2D
end interface broadcast_elem

interface gather_nod
      module procedure gather_nod3D
      module procedure gather_nod2D
      module procedure gather_real4_nod3D
      module procedure gather_real4_nod2D
      module procedure gather_int2_nod3D
      module procedure gather_int2_nod2D
      module procedure gather_real8to4_nod3D
      module procedure gather_real8to4_nod2D
      module procedure gather_nod2D_i
end interface gather_nod

interface gather_nod_fbin
      module procedure gather_nod2D_fbin
      module procedure gather_real4_nod2D_fbin
end interface gather_nod_fbin

interface gather_elem
      module procedure gather_elem3D
      module procedure gather_elem2D
      module procedure gather_real4_elem3D
      module procedure gather_real4_elem2D
      module procedure gather_int2_elem3D
      module procedure gather_int2_elem2D
      module procedure gather_real8to4_elem3D
      module procedure gather_real8to4_elem2D
      module procedure gather_elem2D_i
end interface gather_elem

interface gather_elem_fbin
      module procedure gather_elem2D_fbin
      module procedure gather_real4_elem2D_fbin
end interface gather_elem_fbin

interface gather_edge
      module procedure gather_edg2D
      module procedure gather_edg2D_i
end interface gather_edge


private  ! hides items not listed on public statement
public :: exchange_nod,exchange_elem,broadcast_nod,broadcast_elem, &
          gather_nod, gather_elem, exchange_nod_begin, exchange_nod_end, exchange_elem_begin, &
          exchange_elem_end, gather_edge, &
          exchange_nod_fbin, exchange_elem_fbin, &
          exchange_nod_fbin_begin, exchange_elem_fbin_begin, &
          gather_nod_fbin, gather_elem_fbin
end module g_comm_auto
