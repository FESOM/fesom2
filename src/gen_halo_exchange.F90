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

  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none

contains

#ifdef DEBUG
  ! Only needed in debug mode
  subroutine check_mpi_comm(rn, sn, r_mpitype, s_mpitype, rPE, sPE)
    USE g_PARSUP
    IMPLICIT NONE

    ! General version of the communication routine for 2D nodal fields

    integer, intent(in) :: sn, rn, r_mpitype(:), s_mpitype(:), rPE(:), sPE(:)
    integer  :: n, sdebug, rdebug, status(MPI_STATUS_SIZE), request

    DO n=1,rn    
       call MPI_TYPE_SIZE(r_mpitype(n), rdebug, MPIerr)
       CALL MPI_ISEND(rdebug, 1, MPI_INTEGER, rPE(n), 10, MPI_COMM_FESOM, request, MPIerr)
    END DO

    DO n=1, sn
       call MPI_RECV(sdebug, 1, MPI_INTEGER, sPE(n), 10, MPI_COMM_FESOM,    &
            status, MPIerr)
       call MPI_TYPE_SIZE(s_mpitype(n), rdebug, MPIerr)
       if (sdebug /= rdebug) then
          print *, "Mismatching MPI send/recieve message lengths."
          print *,"Send/receive process numbers: ", mype, '/', sPE(n)
          print *,"Number of send/receive bytes: ", sdebug, '/', rdebug
          call MPI_ABORT( MPI_COMM_FESOM, 1 )
       end if
    END DO
    CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

  END SUBROUTINE check_mpi_comm
#endif


subroutine exchange_nod2D_i(nod_array2D)

USE g_PARSUP
IMPLICIT NONE
 
 integer, intent(inout)  :: nod_array2D(:)

 if (npes > 1) then
    call exchange_nod2D_i_begin(nod_array2D)
    call exchange_nod_end  
endif
END SUBROUTINE exchange_nod2D_i

!=============================================================================

subroutine exchange_nod2D_i_begin(nod_array2D)
  USE o_MESH
  USE g_PARSUP
  IMPLICIT NONE

  ! General version of the communication routine for 2D nodal fields

  integer, intent(inout)  :: nod_array2D(:)
  integer  :: n, sn, rn

  if (npes > 1) then

     sn=com_nod2D%sPEnum
     rn=com_nod2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, r_mpitype_nod2D_i, s_mpitype_nod2D_i,         &
          com_nod2D%rPE, com_nod2D%sPE)
#endif

     DO n=1,rn    

        call MPI_IRECV(nod_array2D, 1, r_mpitype_nod2D_i(n), com_nod2D%rPE(n), &
             com_nod2D%rPE(n), MPI_COMM_FESOM, com_nod2D%req(n), MPIerr) 
     END DO

     DO n=1, sn

        call MPI_ISEND(nod_array2D, 1, s_mpitype_nod2D_i(n), com_nod2D%sPE(n), &
             mype, MPI_COMM_FESOM, com_nod2D%req(rn+n), MPIerr)
     END DO

     com_nod2D%nreq = rn+sn

  endif
END SUBROUTINE exchange_nod2D_i_begin

! ========================================================================
subroutine exchange_nod2D(nod_array2D)

USE g_PARSUP
IMPLICIT NONE

! General version of the communication routine for 2D nodal fields
 
 real(real64), intent(inout)  :: nod_array2D(:)

 if (npes > 1) then
    call exchange_nod2D_begin(nod_array2D)  
    call exchange_nod_end
 end if
 
END SUBROUTINE exchange_nod2D

! ========================================================================
subroutine exchange_nod2D_begin(nod_array2D)
  USE o_MESH
  USE g_PARSUP
  IMPLICIT NONE

  ! General version of the communication routine for 2D nodal fields

  real(real64), intent(inout)  :: nod_array2D(:)

  integer  :: n, sn, rn

  if (npes > 1) then

     sn=com_nod2D%sPEnum
     rn=com_nod2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, r_mpitype_nod2D, s_mpitype_nod2D,           &
          com_nod2D%rPE, com_nod2D%sPE)
#endif

     DO n=1,rn         
        call MPI_IRECV(nod_array2D, 1, r_mpitype_nod2D(n), com_nod2D%rPE(n), &
             com_nod2D%rPE(n), MPI_COMM_FESOM, com_nod2D%req(n), MPIerr) 
     END DO
     DO n=1, sn
        call MPI_ISEND(nod_array2D, 1, s_mpitype_nod2D(n), com_nod2D%sPE(n), &
             mype, MPI_COMM_FESOM, com_nod2D%req(rn+n), MPIerr)
     END DO

     com_nod2D%nreq = rn+sn

  end if

END SUBROUTINE exchange_nod2D_begin
!===============================================
subroutine exchange_nod2D_2fields(nod1_array2D, nod2_array2D)

USE g_PARSUP
IMPLICIT NONE

! General version of the communication routine for 2D nodal fields
 
 real(real64), intent(inout)  :: nod1_array2D(:)
 real(real64), intent(inout)  :: nod2_array2D(:)

 if (npes > 1) then
    call exchange_nod2D_2fields_begin(nod1_array2D, nod2_array2D)  
    call exchange_nod_end
 end if
 
END SUBROUTINE exchange_nod2D_2fields

! ========================================================================
subroutine exchange_nod2D_2fields_begin(nod1_array2D, nod2_array2D)
USE o_MESH
USE g_PARSUP
IMPLICIT NONE

! General version of the communication routine for 2D nodal fields
 
 real(real64), intent(inout)  :: nod1_array2D(:)
 real(real64), intent(inout)  :: nod2_array2D(:)

 integer  :: n, sn, rn

 if (npes > 1) then

  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, r_mpitype_nod2D, s_mpitype_nod2D,           &
          com_nod2D%rPE, com_nod2D%sPE)
#endif

  DO n=1,rn         
     call MPI_IRECV(nod1_array2D, 1, r_mpitype_nod2D(n), com_nod2D%rPE(n), &
               com_nod2D%rPE(n),      MPI_COMM_FESOM, com_nod2D%req(2*n-1), MPIerr) 
 
     call MPI_IRECV(nod2_array2D, 1, r_mpitype_nod2D(n), com_nod2D%rPE(n), &
               com_nod2D%rPE(n)+npes, MPI_COMM_FESOM, com_nod2D%req(2*n),   MPIerr) 
  END DO  
  DO n=1, sn
     call MPI_ISEND(nod1_array2D, 1, s_mpitype_nod2D(n), com_nod2D%sPE(n), &
                    mype,      MPI_COMM_FESOM, com_nod2D%req(2*rn+2*n-1), MPIerr)

     call MPI_ISEND(nod2_array2D, 1, s_mpitype_nod2D(n), com_nod2D%sPE(n), &
                    mype+npes, MPI_COMM_FESOM, com_nod2D%req(2*rn+2*n),   MPIerr)
  END DO

   com_nod2D%nreq = 2*(rn+sn)

end if
 
END SUBROUTINE exchange_nod2D_2fields_begin

!===============================================
subroutine exchange_nod2D_3fields(nod1_array2D, nod2_array2D, nod3_array2D)

USE g_PARSUP
IMPLICIT NONE

! General version of the communication routine for 2D nodal fields
 
 real(real64), intent(inout)  :: nod1_array2D(:)
 real(real64), intent(inout)  :: nod2_array2D(:)
 real(real64), intent(inout)  :: nod3_array2D(:)

 if (npes > 1) then
    call exchange_nod2D_3fields_begin(nod1_array2D, nod2_array2D, nod3_array2D)  
    call exchange_nod_end
 end if
 
END SUBROUTINE exchange_nod2D_3fields

! ========================================================================
subroutine exchange_nod2D_3fields_begin(nod1_array2D, nod2_array2D, nod3_array2D)
USE o_MESH
USE g_PARSUP
IMPLICIT NONE

! General version of the communication routine for 2D nodal fields
 
 real(real64), intent(inout)  :: nod1_array2D(:)
 real(real64), intent(inout)  :: nod2_array2D(:)
 real(real64), intent(inout)  :: nod3_array2D(:)


 integer  :: n, sn, rn

 if (npes > 1) then

  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, r_mpitype_nod2D, s_mpitype_nod2D,           &
          com_nod2D%rPE, com_nod2D%sPE)
#endif

  DO n=1,rn         
     call MPI_IRECV(nod1_array2D, 1, r_mpitype_nod2D(n), com_nod2D%rPE(n), &
               com_nod2D%rPE(n),        MPI_COMM_FESOM, com_nod2D%req(3*n-2), MPIerr) 
 
     call MPI_IRECV(nod2_array2D, 1, r_mpitype_nod2D(n), com_nod2D%rPE(n), &
               com_nod2D%rPE(n)+npes,   MPI_COMM_FESOM, com_nod2D%req(3*n-1), MPIerr) 

     call MPI_IRECV(nod3_array2D, 1, r_mpitype_nod2D(n), com_nod2D%rPE(n), &
               com_nod2D%rPE(n)+2*npes, MPI_COMM_FESOM, com_nod2D%req(3*n),   MPIerr) 
  END DO  
  DO n=1, sn
     call MPI_ISEND(nod1_array2D, 1, s_mpitype_nod2D(n), com_nod2D%sPE(n), &
                    mype,        MPI_COMM_FESOM, com_nod2D%req(3*rn+3*n-2), MPIerr)

     call MPI_ISEND(nod2_array2D, 1, s_mpitype_nod2D(n), com_nod2D%sPE(n), &
                    mype+npes,   MPI_COMM_FESOM, com_nod2D%req(3*rn+3*n-1), MPIerr)

     call MPI_ISEND(nod3_array2D, 1, s_mpitype_nod2D(n), com_nod2D%sPE(n), &
                    mype+2*npes, MPI_COMM_FESOM, com_nod2D%req(3*rn+3*n),   MPIerr)
  END DO

   com_nod2D%nreq = 3*(rn+sn)

end if
 
END SUBROUTINE exchange_nod2D_3fields_begin

! ========================================================================
subroutine exchange_nod3D(nod_array3D)

USE g_PARSUP
IMPLICIT NONE

real(real64), intent(inout) :: nod_array3D(:,:) 
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
 
if (npes > 1) then
   call exchange_nod3D_begin(nod_array3D)
   call exchange_nod_end
endif
END SUBROUTINE exchange_nod3D

! ========================================================================
subroutine exchange_nod3D_begin(nod_array3D)
USE o_MESH
USE g_PARSUP
IMPLICIT NONE


real(real64), intent(inout) :: nod_array3D(:,:) 
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
 
 integer  :: n, sn, rn
 integer  :: nz, nl1

 if (npes > 1) then
    sn=com_nod2D%sPEnum
    rn=com_nod2D%rPEnum

    nl1=ubound(nod_array3D,1)

    if ((nl1<ubound(r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(r_mpitype_nod3D, 2))) then
       if (mype==0) then
          print *,'Subroutine exchange_nod3D not implemented for',nl1,'layers.'
          print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
       endif
       call par_ex(1)
    endif

    ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
    call check_mpi_comm(rn, sn, r_mpitype_nod3D(:,nl1,1), s_mpitype_nod3D(:,nl1,1), &
         com_nod2D%rPE, com_nod2D%sPE)
#endif

    DO n=1,rn    
       call MPI_IRECV(nod_array3D, 1, r_mpitype_nod3D(n,nl1,1), com_nod2D%rPE(n), &
            com_nod2D%rPE(n), MPI_COMM_FESOM, com_nod2D%req(n), MPIerr) 
    END DO

    DO n=1, sn
       call MPI_ISEND(nod_array3D, 1, s_mpitype_nod3D(n,nl1,1), com_nod2D%sPE(n), &
            mype, MPI_COMM_FESOM, com_nod2D%req(rn+n), MPIerr)
    END DO

    com_nod2D%nreq = rn+sn

 endif
END SUBROUTINE exchange_nod3D_begin

! ========================================================================
subroutine exchange_nod3D_2fields(nod1_array3D,nod2_array3D)

USE g_PARSUP
IMPLICIT NONE

real(real64), intent(inout) :: nod1_array3D(:,:) 
real(real64), intent(inout) :: nod2_array3D(:,:) 
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
 
if (npes > 1) then
   call exchange_nod3D_2fields_begin(nod1_array3D,nod2_array3D)
   call exchange_nod_end
endif
END SUBROUTINE exchange_nod3D_2fields

! ========================================================================
subroutine exchange_nod3D_2fields_begin(nod1_array3D,nod2_array3D)
USE o_MESH
USE g_PARSUP
IMPLICIT NONE


real(real64), intent(inout) :: nod1_array3D(:,:) 
real(real64), intent(inout) :: nod2_array3D(:,:) 
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
 
 integer  :: n, sn, rn
 integer  :: nz, nl1, nl2

 if (npes > 1) then
    sn=com_nod2D%sPEnum
    rn=com_nod2D%rPEnum

    nl1 = ubound(nod1_array3D,1)

    if ((nl1<ubound(r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(r_mpitype_nod3D, 2))) then
       if (mype==0) then
          print *,'Subroutine exchange_nod3D not implemented for',nl1,'layers.'
          print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
       endif
       call par_ex(1)
    endif

    nl2 = ubound(nod2_array3D,1)
    if ((nl2<ubound(r_mpitype_nod3D, 2)-1) .or. (nl2>ubound(r_mpitype_nod3D, 2))) then
       if (mype==0) then
          print *,'Subroutine exchange_nod3D not implemented for',nl2,'layers.'
          print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
       endif
       call par_ex(1)
    endif

#ifdef DEBUG
    call check_mpi_comm(rn, sn, r_mpitype_nod3D(:,nl1,1), s_mpitype_nod3D(:,nl1,1), &
         com_nod2D%rPE, com_nod2D%sPE)
#endif

    DO n=1,rn    
       call MPI_IRECV(nod1_array3D, 1, r_mpitype_nod3D(n,nl1,1), com_nod2D%rPE(n), &
            com_nod2D%rPE(n),      MPI_COMM_FESOM, com_nod2D%req(2*n-1), MPIerr)  

       call MPI_IRECV(nod2_array3D, 1, r_mpitype_nod3D(n,nl2,1), com_nod2D%rPE(n), &
            com_nod2D%rPE(n)+npes, MPI_COMM_FESOM, com_nod2D%req(2*n  ), MPIerr) 
    END DO

    DO n=1, sn
       call MPI_ISEND(nod1_array3D, 1, s_mpitype_nod3D(n,nl1,1), com_nod2D%sPE(n), &
            mype,      MPI_COMM_FESOM, com_nod2D%req(2*rn+2*n-1), MPIerr)

       call MPI_ISEND(nod2_array3D, 1, s_mpitype_nod3D(n,nl2,1), com_nod2D%sPE(n), &
            mype+npes, MPI_COMM_FESOM, com_nod2D%req(2*rn+2*n), MPIerr)
    END DO

    com_nod2D%nreq = 2*(rn+sn)

 endif
END SUBROUTINE exchange_nod3D_2fields_begin
! ========================================================================
subroutine exchange_nod3D_n(nod_array3D)
USE o_MESH
USE g_PARSUP
IMPLICIT NONE

real(real64), intent(inout) :: nod_array3D(:,:,:) 
 
if (npes>1) then
   call exchange_nod3D_n_begin(nod_array3D)
   call exchange_nod_end
endif

END SUBROUTINE exchange_nod3D_n

!=================================================

subroutine exchange_nod3D_n_begin(nod_array3D)
USE o_MESH
USE g_PARSUP
IMPLICIT NONE

real(real64), intent(inout) :: nod_array3D(:,:,:) 
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
 
 integer  :: n, sn, rn
 integer  :: nz, nl1, n_val

if (npes>1) then
 ! nod_array3D(n_val,nl1,nod2D_size)
  nl1= ubound(nod_array3D,2)
  n_val = ubound(nod_array3D,1)

  if ((nl1<ubound(r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(r_mpitype_nod3D, 2)) .or. (n_val > 3)) then

     ! This routine also works for swapped dimensions nod_array3D(nl1,n_val, nod2D_size) 
     nl1   = ubound(nod_array3D,1)
     n_val = ubound(nod_array3D,2)

     if ((nl1<ubound(r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(r_mpitype_nod3D, 2)) .or. (n_val > 3)) then
        if (mype==0) then
           print *,'Subroutine exchange_nod3D_n not implemented for'
           print *,nl1,'layers and / or ',n_val,'values per element.'
           print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
        endif
        call par_ex(1)
     endif
  endif
  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum

  ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
  call check_mpi_comm(rn, sn, r_mpitype_nod3D(:,nl1,n_val),                 &
       s_mpitype_nod3D(:,nl1,n_val), com_nod2D%rPE, com_nod2D%sPE)
#endif

  DO n=1,rn    
     call MPI_IRECV(nod_array3D, 1, r_mpitype_nod3D(n,nl1,n_val), com_nod2D%rPE(n), &
          com_nod2D%rPE(n), MPI_COMM_FESOM, com_nod2D%req(n), MPIerr) 
  END DO
 
  DO n=1, sn
     call MPI_ISEND(nod_array3D, 1, s_mpitype_nod3D(n,nl1,n_val), com_nod2D%sPE(n), &
          mype, MPI_COMM_FESOM, com_nod2D%req(rn+n), MPIerr)
  END DO

  com_nod2D%nreq = rn+sn

 endif


END SUBROUTINE exchange_nod3D_n_begin

!=======================================
! AND WAITING
!=======================================
 
SUBROUTINE exchange_nod_end
  USE g_PARSUP

if (npes > 1) &
  call MPI_WAITALL(com_nod2D%nreq, com_nod2D%req, MPI_STATUSES_IGNORE, MPIerr)

END SUBROUTINE exchange_nod_end

SUBROUTINE exchange_elem_end

  USE g_PARSUP

  if (npes > 1) then
     if (elem_full_flag) then
        call MPI_WAITALL(com_elem2D_full%nreq, &
             com_elem2D_full%req, MPI_STATUSES_IGNORE, MPIerr)     
     else
        call MPI_WAITALL(com_elem2D%nreq, &
             com_elem2D%req, MPI_STATUSES_IGNORE, MPIerr)     
     endif
  end if
END SUBROUTINE exchange_elem_end
! ========================================================================

!nr  Not used, no MPI datatype built (yet)
!
!!$subroutine exchange_edge3D(edge_array3D)
!!$  use o_MESH
!!$  use g_PARSUP 
!!$  implicit none
!!$  
!!$ ! Communication of edge based data stored in (vertical, horizontal) format 
!!$
!!$ INTEGER  :: sreq(maxPEnum)
!!$ INTEGER  :: rreq(maxPEnum)
!!$ INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
!!$ INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
!!$ integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
!!$ integer  :: nz, nh, nc
!!$ real(real64) :: edge_array3D(nl-1,edge2D) 
!!$  
!!$  sn=com_edge2D%sPEnum
!!$  rn=com_edge2D%rPEnum
!!$  ! Put data to be communicated into send buffer 
!!$
!!$
!!$  do n=1, sn
!!$     nini=com_edge2D%sptr(n)
!!$     nend=com_edge2D%sptr(n+1) - 1
!!$      nc=0
!!$      DO nh=nini, nend
!!$         DO nz=1, nl-1
!!$	 nc=nc+1
!!$         s_buff_edge3D(n)%array(nc)=edge_array3D(nz,com_edge2D%slist(nh))
!!$	 END DO
!!$      END DO	 
!!$  end do
!!$
!!$
!!$  do n=1, sn
!!$     dest=com_edge2D%sPE(n)
!!$     nini=com_edge2D%sptr(n)
!!$     offset=(com_edge2D%sptr(n+1) - nini)*(nl-1)
!!$
!!$     call MPI_ISEND(s_buff_edge3D(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
!!$          MPI_COMM_FESOM, sreq(n), MPIerr)
!!$  end do
!!$  do n=1, rn
!!$     source=com_edge2D%rPE(n)
!!$     nini=com_edge2D%rptr(n)
!!$     offset=(com_edge2D%rptr(n+1) - nini)*(nl-1)
!!$
!!$     call MPI_IRECV(r_buff_edge3D(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
!!$          source, MPI_COMM_FESOM, rreq(n), MPIerr) 
!!$  end do
!!$
!!$  call MPI_WAITALL(sn,sreq,sstat, MPIerr)
!!$  call MPI_WAITALL(rn,rreq,rstat, MPIerr)
!!$
!!$  ! Put received data to their destination
!!$
!!$  do n=1, rn
!!$     nini=com_edge2D%rptr(n)
!!$     nend=com_edge2D%rptr(n+1) - 1
!!$      nc=0
!!$      DO nh=nini, nend
!!$         DO nz=1, nl-1
!!$	 nc=nc+1
!!$	 edge_array3D(nz,com_edge2D%rlist(nh))=r_buff_edge3D(n)%array(nc)
!!$         END DO
!!$      END DO	  
!!$  end do
!!$
!!$end subroutine exchange_edge3D
!==========================================================================

!!$subroutine exchange_edge2D(edge_array2D)
!!$  use o_MESH
!!$  use g_PARSUP 
!!$  implicit none
!!$
!!$! General version of the communication routine for 2D edge fields
!!$! This routine is not split, it is used only once during setup. 
!!$ real(real64), intent(inout)  :: edge_array2D(:)
!!$
!!$ integer  :: n, sn, rn
!!$
!!$ if (npes> 1) then
!!$  sn=com_edge2D%sPEnum
!!$  rn=com_edge2D%rPEnum
!!$
!!$  DO n=1,rn         
!!$     call MPI_IRECV(edge_array2D, 1, r_mpitype_edge2D(n), com_edge2D%rPE(n), &
!!$               com_edge2D%rPE(n), MPI_COMM_FESOM, com_edge2D%req(n), MPIerr) 
!!$  END DO  
!!$  DO n=1, sn
!!$     call MPI_ISEND(edge_array2D, 1, s_mpitype_edge2D(n), com_edge2D%sPE(n), &
!!$                    mype, MPI_COMM_FESOM, com_edge2D%req(rn+n), MPIerr)
!!$  END DO
!!$  
!!$  call MPI_WAITALL(rn+sn,com_edge2D%req,MPI_STATUSES_IGNORE, MPIerr)
!!$
!!$  endif
!!$
!!$end subroutine exchange_edge2D
!=============================================================================
subroutine exchange_elem3D(elem_array3D)

USE g_PARSUP 
IMPLICIT NONE

 real(real64), intent(inout) :: elem_array3D(:,:) 

 call exchange_elem3D_begin(elem_array3D)
 call exchange_elem_end

END SUBROUTINE exchange_elem3D
!===========================================
subroutine exchange_elem3D_begin(elem_array3D)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format

real(real64), intent(inout) :: elem_array3D(:,:) 
integer  :: n, sn, rn, nl1

if (npes> 1) then

   nl1=ubound(elem_array3D,1)

   if (ubound(elem_array3D,2)<=myDim_elem2D+eDim_elem2D) then

      elem_full_flag = .false.

      sn=com_elem2D%sPEnum
      rn=com_elem2D%rPEnum

      if (nl1==ubound(r_mpitype_elem3D, 2) .or. nl1==ubound(r_mpitype_elem3D, 2)-1) then

         ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
         call check_mpi_comm(rn, sn, r_mpitype_elem3D(:,nl1,1), s_mpitype_elem3D(:,nl1,1), &
              com_elem2D%rPE, com_elem2D%sPE)
#endif

         DO n=1,rn         
            call MPI_IRECV(elem_array3D, 1, r_mpitype_elem3D(n,nl1,1), com_elem2D%rPE(n), &
                 com_elem2D%rPE(n), MPI_COMM_FESOM, &
                 com_elem2D%req(n), MPIerr)
         END DO
         DO n=1, sn
            call MPI_ISEND(elem_array3D, 1, s_mpitype_elem3D(n,nl1,1), com_elem2D%sPE(n), &
                 mype,    MPI_COMM_FESOM, &
                 com_elem2D%req(rn+n), MPIerr)
         END DO

      elseif (nl1 <= 4) then
         ! In fact, this is a 2D-array with up to 4 values, e.g. derivatives

         ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
         call check_mpi_comm(rn, sn, r_mpitype_elem2D(:,nl1), s_mpitype_elem2D(:,nl1), &
              com_elem2D%rPE, com_elem2D%sPE)
#endif

         DO n=1,rn         
            call MPI_IRECV(elem_array3D, 1, r_mpitype_elem2D(n,nl1), com_elem2D%rPE(n), &
                 com_elem2D%rPE(n), MPI_COMM_FESOM, &
                 com_elem2D%req(n), MPIerr) 
         END DO
         DO n=1, sn
            call MPI_ISEND(elem_array3D, 1, s_mpitype_elem2D(n,nl1), com_elem2D%sPE(n), &
                 mype,    MPI_COMM_FESOM, &
                 com_elem2D%req(rn+n), MPIerr)
         END DO
      else
         if (mype==0) print *,'Sorry, no MPI datatype prepared for',nl1,'values per element (exchange_elem3D)'
         call par_ex(1)
      endif

      com_elem2D%nreq = rn+sn

   else

      elem_full_flag = .true.

      sn=com_elem2D_full%sPEnum
      rn=com_elem2D_full%rPEnum

      if (nl1==ubound(r_mpitype_elem3D_full, 2) .or. nl1==ubound(r_mpitype_elem3D_full, 2)-1) then
         ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
         call check_mpi_comm(rn, sn, r_mpitype_elem3D_full(:,nl1,1), &
              s_mpitype_elem3D_full(:,nl1,1), com_elem2D_full%rPE, com_elem2D_full%sPE)
#endif

         DO n=1,rn         
            call MPI_IRECV(elem_array3D, 1, r_mpitype_elem3D_full(n,nl1,1), &
                 com_elem2D_full%rPE(n), &
                 com_elem2D_full%rPE(n), MPI_COMM_FESOM, &
                 com_elem2D_full%req(n), MPIerr) 
         END DO
         DO n=1, sn
            call MPI_ISEND(elem_array3D, 1, s_mpitype_elem3D_full(n,nl1,1), &
                 com_elem2D_full%sPE(n), & 
                 mype,    MPI_COMM_FESOM, &
                 com_elem2D_full%req(rn+n), MPIerr)
         END DO
      elseif (nl1 <= 4) then
         ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
         call check_mpi_comm(rn, sn, r_mpitype_elem2D_full(:,nl1), &
              s_mpitype_elem2D_full(:,nl1), com_elem2D_full%rPE, com_elem2D_full%sPE)
#endif

         ! In fact, this is a 2D-array with up to 4 values, e.g. derivatives
         DO n=1,rn         
            call MPI_IRECV(elem_array3D, 1, r_mpitype_elem2D_full(n,nl1), &
                 com_elem2D_full%rPE(n), &
                 com_elem2D_full%rPE(n), MPI_COMM_FESOM, &
                 com_elem2D_full%req(n), MPIerr) 
         END DO
         DO n=1, sn
            call MPI_ISEND(elem_array3D, 1, s_mpitype_elem2D_full(n,nl1), &
                 com_elem2D_full%sPE(n), & 
                 mype,    MPI_COMM_FESOM, &
                 com_elem2D_full%req(rn+n), MPIerr)
         END DO
      else
         if (mype==0) print *,'Sorry, no MPI datatype prepared for',nl1,'values per element (exchange_elem3D)'
         call par_ex(1)
      endif

      com_elem2D_full%nreq = rn+sn

   endif

endif

END SUBROUTINE exchange_elem3D_begin

!=============================================================================
subroutine exchange_elem3D_n(elem_array3D)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
 
 real(real64), intent(inout) :: elem_array3D(:,:,:) 

 if (npes> 1) then
    call exchange_elem3D_n_begin(elem_array3D)
    call exchange_elem_end
 endif
END SUBROUTINE exchange_elem3D_n
!=============================================================================
subroutine exchange_elem3D_n_begin(elem_array3D)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
 
 real(real64), intent(inout) :: elem_array3D(:,:,:) 
 integer  :: n, sn, rn, n_val, nl1

 if (npes> 1) then
 nl1   = ubound(elem_array3D,2)
 n_val = ubound(elem_array3D,1)
 
  if ((nl1<ubound(r_mpitype_elem3D, 2)-1) .or. (nl1>ubound(r_mpitype_elem3D, 2)) .or. (n_val > 4)) then

     ! This routine also works for swapped dimensions elem_array3D(nl1,n_val, elem2D_size) 
     nl1= ubound(elem_array3D,1)
     n_val = ubound(elem_array3D,2)

     if ((nl1<ubound(r_mpitype_elem3D, 2)-1) .or. (nl1>ubound(r_mpitype_elem3D, 2)) .or. (n_val > 4)) then
        if (mype==0) then
           print *,'Subroutine exchange_elem3D_n not implemented for'
           print *,nl1,'layers and / or ',n_val,'values per element.'
           print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
        endif
        call par_ex(1)
     endif
  endif

 if (ubound(elem_array3D,3)<=myDim_elem2D+eDim_elem2D) then

    elem_full_flag = .false.

     sn=com_elem2D%sPEnum
     rn=com_elem2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, r_mpitype_elem3D(:,nl1,n_val), &
          s_mpitype_elem3D(:,nl1,n_val), com_elem2D%rPE, com_elem2D%sPE)
#endif

     DO n=1,rn         
        call MPI_IRECV(elem_array3D, 1, r_mpitype_elem3D(n,nl1,n_val), com_elem2D%rPE(n), &
                       com_elem2D%rPE(n), MPI_COMM_FESOM, com_elem2D%req(n), MPIerr) 
     END DO
     DO n=1, sn
        call MPI_ISEND(elem_array3D, 1, s_mpitype_elem3D(n,nl1,n_val), com_elem2D%sPE(n), &
                       mype, MPI_COMM_FESOM, com_elem2D%req(rn+n), MPIerr)
     END DO

     com_elem2D%nreq = rn+sn

  else
     
     elem_full_flag = .true.

     sn=com_elem2D_full%sPEnum
     rn=com_elem2D_full%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, r_mpitype_elem3D_full(:,nl1,n_val), &
          s_mpitype_elem3D_full(:,nl1,n_val), com_elem2D_full%rPE, com_elem2D_full%sPE)
#endif

     DO n=1,rn         
        call MPI_IRECV(elem_array3D, 1, r_mpitype_elem3D_full(n,nl1,n_val), com_elem2D_full%rPE(n), &
                       com_elem2D_full%rPE(n), MPI_COMM_FESOM, com_elem2D_full%req(n), MPIerr) 
     END DO
     DO n=1, sn
        call MPI_ISEND(elem_array3D, 1, s_mpitype_elem3D_full(n,nl1,n_val), com_elem2D_full%sPE(n), &
                       mype, MPI_COMM_FESOM, com_elem2D_full%req(rn+n), MPIerr)
     END DO

     com_elem2D_full%nreq = rn+sn

  end if     
  

endif  
END SUBROUTINE exchange_elem3D_n_begin
!========================================================================
subroutine exchange_elem2D(elem_array2D)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
 
 real(real64), intent(inout) :: elem_array2D(:) 

 if (npes> 1) then
    call exchange_elem2D_begin(elem_array2D)
    call exchange_elem_end
 end if
  
END SUBROUTINE exchange_elem2D
!========================================================================
subroutine exchange_elem2D_begin(elem_array2D)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
 
 real(real64), intent(inout) :: elem_array2D(:) 
 integer  :: n, sn, rn

 if (npes> 1) then

  if (ubound(elem_array2D,1)<=myDim_elem2D+eDim_elem2D) then

     elem_full_flag = .false.

     sn=com_elem2D%sPEnum
     rn=com_elem2D%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, r_mpitype_elem2D(:,1), s_mpitype_elem2D(:,1), &
          com_elem2D%rPE, com_elem2D%sPE)
#endif

     DO n=1,rn         
        call MPI_IRECV(elem_array2D, 1, r_mpitype_elem2D(n,1), com_elem2D%rPE(n), &
                       com_elem2D%rPE(n), MPI_COMM_FESOM, com_elem2D%req(n), MPIerr) 
     END DO
     DO n=1, sn
        call MPI_ISEND(elem_array2D, 1, s_mpitype_elem2D(n,1), com_elem2D%sPE(n), &
                       mype, MPI_COMM_FESOM, com_elem2D%req(rn+n), MPIerr)
     END DO

     com_elem2D%nreq = rn+sn

  else
     elem_full_flag = .true.

     sn=com_elem2D_full%sPEnum
     rn=com_elem2D_full%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, r_mpitype_elem2D_full(:,1), s_mpitype_elem2D_full(:,1), &
          com_elem2D_full%rPE, com_elem2D_full%sPE)
#endif

     DO n=1,rn         
        call MPI_IRECV(elem_array2D, 1, r_mpitype_elem2D_full(n,1), com_elem2D_full%rPE(n), &
                       com_elem2D_full%rPE(n), MPI_COMM_FESOM, com_elem2D_full%req(n), MPIerr) 
     END DO
     DO n=1, sn
        call MPI_ISEND(elem_array2D, 1, s_mpitype_elem2D_full(n,1), com_elem2D_full%sPE(n), &
                       mype, MPI_COMM_FESOM, com_elem2D_full%req(rn+n), MPIerr)
     END DO

     com_elem2D_full%nreq = rn+sn

  end if     

end if
  
END SUBROUTINE exchange_elem2D_begin
! ========================================================================
subroutine exchange_elem2D_i(elem_array2D)
!Exchange with ALL(!) the neighbours
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

 integer, intent(inout)  :: elem_array2D(:)

 integer  :: n, sn, rn

 if (npes> 1) then
    call exchange_elem2D_i_begin(elem_array2D)
    call exchange_elem_end
end if

END SUBROUTINE exchange_elem2D_i
!=============================================================================
subroutine exchange_elem2D_i_begin(elem_array2D)
!Exchange with ALL(!) the neighbours
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

 integer, intent(inout)  :: elem_array2D(:)

 integer  :: n, sn, rn

 if (npes> 1) then

    elem_full_flag = .true.

    sn=com_elem2D_full%sPEnum
    rn=com_elem2D_full%rPEnum

     ! Check MPI point-to-point communication for consistency
#ifdef DEBUG
     call check_mpi_comm(rn, sn, r_mpitype_elem2D_full_i, s_mpitype_elem2D_full_i,        &
          com_elem2D_full%rPE, com_elem2D_full%sPE)
#endif

    DO n=1,rn       
       call MPI_IRECV(elem_array2D, 1, r_mpitype_elem2D_full_i(n), com_elem2D_full%rPE(n), &
            com_elem2D_full%rPE(n), MPI_COMM_FESOM, com_elem2D_full%req(n), MPIerr) 
    END DO
 
    DO n=1, sn
     
       call MPI_ISEND(elem_array2D, 1, s_mpitype_elem2D_full_i(n), com_elem2D_full%sPE(n), &
            mype, MPI_COMM_FESOM, com_elem2D_full%req(rn+n), MPIerr)
    END DO

    com_elem2D_full%nreq = rn+sn
  
end if

END SUBROUTINE exchange_elem2D_i_begin
!=============================================================================



! ========================================================================
! Broadcast routines
! Many because of different sizes.
! ========================================================================
subroutine broadcast_nod3D(arr3D, arr3Dglobal)
! Distribute the nodal information available on 0 PE to other PEs
use g_PARSUP
USE o_MESH

IMPLICIT NONE

INTEGER      :: nz, counter,nl1
integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:) ::  irecvbuf
real(real64) ::  arr3D(:,:)
real(real64) ::  arr3Dglobal(:,:)
real(real64), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf
integer       :: node_size

node_size=myDim_nod2D+eDim_nod2D
nl1=ubound(arr3D,1)
IF ( mype == 0 ) THEN
    if (npes>1) then
    arr3D(:,1:node_size)=arr3Dglobal(:,myList_nod2D(1:node_size))
    end if
    DO  n = 1, npes-1
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(sendbuf(nTS*nl1), irecvbuf(nTS))

       CALL MPI_RECV(irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, MPI_COMM_FESOM, status, MPIerr )
       counter=0
       DO i = 1, nTS
          DO nz=1, nl1
             counter=counter+1
             sendbuf(counter) = arr3Dglobal(nz,irecvbuf(i))
          ENDDO
       ENDDO

       CALL MPI_SEND(sendbuf(1), nTS*nl1, MPI_DOUBLE_PRECISION, &
                   sender, 2, MPI_COMM_FESOM, MPIerr )

       DEALLOCATE(irecvbuf, sendbuf)
    ENDDO
ELSE
    CALL MPI_SEND( node_size, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( myList_nod2D(1), node_size, MPI_INTEGER, 0, 1, &
                   MPI_COMM_FESOM, MPIerr )

    ALLOCATE(recvbuf(node_size*nl1))
    CALL MPI_RECV( recvbuf(1), node_size*nl1, MPI_DOUBLE_PRECISION, 0, &
                      2, MPI_COMM_FESOM, status, MPIerr )
    counter=0
    DO n = 1, node_size
       DO nz=1, nl1
         counter=counter+1
         arr3D(nz,n)=recvbuf(counter)
       ENDDO
    ENDDO

    DEALLOCATE(recvbuf)
ENDIF
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
end subroutine broadcast_nod3D
!
!============================================================================
!
subroutine broadcast_nod2D(arr2D, arr2Dglobal)
! A 2D version of the previous routine
use g_PARSUP
USE o_MESH
IMPLICIT NONE

real(real64) ::  arr2D(:)
real(real64) ::  arr2Dglobal(:)

integer                                  ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:)       ::  irecvbuf
real(real64), ALLOCATABLE, DIMENSION(:)  ::  sendbuf
integer       :: node_size

node_size=myDim_nod2D+eDim_nod2D

IF ( mype == 0 ) THEN
    if (npes>1) then
    arr2D(1:node_size)=arr2Dglobal(myList_nod2D(1:node_size))
    end if
    DO  n = 1, npes-1
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(sendbuf(nTS), irecvbuf(nTS))

       CALL MPI_RECV(irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, MPI_COMM_FESOM, status, MPIerr )
       DO i = 1, nTS
             sendbuf(i) = arr2Dglobal(irecvbuf(i))
       ENDDO

       CALL MPI_SEND(sendbuf(1), nTS, MPI_DOUBLE_PRECISION, &
                   sender, 2, MPI_COMM_FESOM, MPIerr )

       DEALLOCATE(irecvbuf, sendbuf)
    ENDDO
ELSE
    CALL MPI_SEND( node_size, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( myList_nod2D(1), node_size, MPI_INTEGER, 0, 1, &
                   MPI_COMM_FESOM, MPIerr )
    CALL MPI_RECV( arr2D(1), node_size, MPI_DOUBLE_PRECISION, 0, &
                      2, MPI_COMM_FESOM, status, MPIerr )
ENDIF
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
end subroutine broadcast_nod2D
!
!============================================================================
!
subroutine broadcast_elem3D(arr3D, arr3Dglobal)
! Distribute the elemental information available on 0 PE to other PEs
use g_PARSUP
USE o_MESH

IMPLICIT NONE

INTEGER      :: nz, counter,nl1
integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:) ::  irecvbuf
real(real64) ::  arr3D(:,:)
real(real64) ::  arr3Dglobal(:,:)
real(real64), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf
integer       :: elem_size

elem_size=myDim_elem2D+eDim_elem2D

nl1=ubound(arr3D,1)
IF ( mype == 0 ) THEN
    if (npes>1) then
    arr3D(:,1:elem_size)=arr3Dglobal(:,myList_elem2D(1:elem_size))
    end if
    DO  n = 1, npes-1
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(sendbuf(nTS*nl1), irecvbuf(nTS))

       CALL MPI_RECV(irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, MPI_COMM_FESOM, status, MPIerr )
       counter=0
       DO i = 1, nTS
          DO nz=1, nl1
             counter=counter+1
             sendbuf(counter) = arr3Dglobal(nz,irecvbuf(i))
          ENDDO
       ENDDO

       CALL MPI_SEND(sendbuf(1), nTS*nl1, MPI_DOUBLE_PRECISION, &
                   sender, 2, MPI_COMM_FESOM, MPIerr )

       DEALLOCATE(irecvbuf, sendbuf)
    ENDDO
ELSE
    CALL MPI_SEND( elem_size, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( myList_elem2D(1), elem_size, MPI_INTEGER, 0, 1, &
                   MPI_COMM_FESOM, MPIerr )

    ALLOCATE(recvbuf(elem_size*nl1))
    CALL MPI_RECV( recvbuf(1), elem_size*nl1, MPI_DOUBLE_PRECISION, 0, &
                      2, MPI_COMM_FESOM, status, MPIerr )
    counter=0
    DO n = 1, elem_size
       DO nz=1, nl1
       counter=counter+1
       arr3D(nz,n)=recvbuf(counter)
       ENDDO
    ENDDO

    DEALLOCATE(recvbuf)
ENDIF
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
end subroutine broadcast_elem3D
!
!============================================================================
!
subroutine broadcast_elem2D(arr2D, arr2Dglobal)
! A 2D version of the previous routine
use g_PARSUP
USE o_MESH
IMPLICIT NONE

integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:)       ::  irecvbuf

real(real64) ::  arr2D(:)
real(real64) ::  arr2Dglobal(:)
real(real64), ALLOCATABLE, DIMENSION(:) ::  sendbuf
integer       :: elem_size

elem_size=myDim_elem2D+eDim_elem2D



IF ( mype == 0 ) THEN
    if (npes>1) then
    arr2D(1:elem_size)=arr2Dglobal(myList_elem2D(1:elem_size))
    end if
    DO  n = 1, npes-1
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(sendbuf(1:nTS), irecvbuf(nTS))

       CALL MPI_RECV(irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, MPI_COMM_FESOM, status, MPIerr )
       DO i = 1, nTS
             sendbuf(i) = arr2Dglobal(irecvbuf(i))
       ENDDO

       CALL MPI_SEND(sendbuf(1), nTS, MPI_DOUBLE_PRECISION, &
                   sender, 2, MPI_COMM_FESOM, MPIerr )

       DEALLOCATE(irecvbuf, sendbuf)
    ENDDO
ELSE
    CALL MPI_SEND( elem_size, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( myList_elem2D(1), elem_size, MPI_INTEGER, 0, 1, &
                   MPI_COMM_FESOM, MPIerr )
    CALL MPI_RECV( arr2D(1), elem_size, MPI_DOUBLE_PRECISION, 0, &
                      2, MPI_COMM_FESOM, status, MPIerr )
ENDIF
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
end subroutine broadcast_elem2D
!
!============================================================================
!
subroutine gather_nod3D(arr3D, arr3D_global)

! Make nodal information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(real64) ::  arr3D(:,:)
real(real64) ::  arr3D_global(:,:)
real(real64), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D

 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,ubound(arr3D_global,2)))

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl1
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_nod2D*nl1, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

end if
end subroutine gather_nod3D
!
!============================================================================
!
subroutine gather_real4_nod3D(arr3D, arr3D_global)

! Make nodal information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(real32)  ::  arr3D(:,:)
real(real32)  ::  arr3D_global(:,:)
real(real32), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D

 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,ubound(arr3D_global,2)))

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl1
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_nod2D*nl1, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

end if
end subroutine gather_real4_nod3D
!=======================================================

subroutine gather_int2_nod3D(arr3D, arr3D_global)

! Make nodal information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH

IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

integer(int16) ::  arr3D(:,:)
integer(int16) ::  arr3D_global(:,:)
integer(int16), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D

 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,ubound(arr3D_global,2)))

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl1
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_SHORT, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_nod2D*nl1, MPI_SHORT, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

end if
end subroutine gather_int2_nod3D
!==============================================
subroutine gather_nod2D(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real64) ::  arr2D(:)
real(real64) ::  arr2D_global(:)
real(real64), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, n2D

 if (npes> 1) then

CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

endif
end subroutine gather_nod2D
!==============================================
subroutine gather_real4_nod2D(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real32)  ::  arr2D(:)
real(real32)  ::  arr2D_global(:)
real(real32), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, n2D

 if (npes> 1) then

CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

endif
end subroutine gather_real4_nod2D

!==============================================
subroutine gather_int2_nod2D(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

integer(int16) ::  arr2D(:)
integer(int16) ::  arr2D_global(:)
integer(int16), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, n2D

 if (npes> 1) then

CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_SHORT, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_SHORT, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

endif
end subroutine gather_int2_nod2D

!============================================================================
subroutine gather_elem3D(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(real64) ::  arr3D(:,:)
real(real64) ::  arr3D_global(:,:)
real(real64), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D, ende, err_alloc
integer        :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( mype == 0 ) THEN
   
   if (npes>1) then
! 
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl1
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   

      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                        = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_elem2D*nl1, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

endif
end subroutine gather_elem3D

!===================================================================

subroutine gather_real4_elem3D(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER        :: nl1
integer        :: n

real(real32)   ::  arr3D(:,:)
real(real32)   ::  arr3D_global(:,:)
real(real32), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D, ende, err_alloc
integer        :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( mype == 0 ) THEN
   
   if (npes>1) then
! 
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl1
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   

      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                        = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_elem2D*nl1, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

endif
end subroutine gather_real4_elem3D


!===================================================================

subroutine gather_int2_elem3D(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

integer(int16) ::  arr3D(:,:)
integer(int16) ::  arr3D_global(:,:)
integer(int16), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D, ende, err_alloc
integer        :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( mype == 0 ) THEN
   
   if (npes>1) then
! 
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl1
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_SHORT, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   

      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                        = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_elem2D*nl1, MPI_SHORT, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

endif
end subroutine gather_int2_elem3D


!==============================================
subroutine gather_elem2D(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real64) ::  arr2D(:)
real(real64) ::  arr2D_global(:)
real(real64), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, e2D


 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then

      allocate(recvbuf(remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_elem2D, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF
end if

end subroutine gather_elem2D

!==============================================
subroutine gather_real4_elem2D(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real32)  ::  arr2D(:)
real(real32)  ::  arr2D_global(:)
real(real32), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, e2D


 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then

      allocate(recvbuf(remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_elem2D, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF
end if

end subroutine gather_real4_elem2D

!==============================================
subroutine gather_int2_elem2D(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

integer(int16) ::  arr2D(:)
integer(int16) ::  arr2D_global(:)
integer(int16), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, e2D


 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then

      allocate(recvbuf(remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_SHORT, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_elem2D, MPI_SHORT, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF
end if

end subroutine gather_int2_elem2D


!============================================================================
subroutine gather_real8to4_nod3D(arr3D, arr3D_global)

! Make nodal information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(real64)  ::  arr3D(:,:)
real(real32)   ::  arr3D_global(:,:)
real(real32), allocatable :: recvbuf(:,:)
real(real32), allocatable :: sendbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D, ierr

 if (npes> 1) then

CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1, ubound(arr3D_global,2)))

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl1
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                     = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   allocate(sendbuf(nl1,myDim_nod2D))
   sendbuf(1:nl1,1:myDim_nod2D) = arr3D(1:nl1,1:myDim_nod2D)
   
   call MPI_SEND(sendbuf, myDim_nod2D*nl1, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   deallocate(sendbuf)
   
ENDIF

end if

end subroutine gather_real8to4_nod3D
!==============================================
subroutine gather_real8to4_nod2D(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real64)  ::  arr2D(:)
real(real32)   :: arr2D_global(:)
real(real32)   :: sendbuf(myDim_nod2D)
real(real64), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, n2D

! Consider MPI-datatypes to recv directly into arr2D_global!

 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   sendbuf(1:myDim_nod2D) = real(arr2D(1:myDim_nod2D),real32)

   call MPI_SEND(sendbuf, myDim_nod2D, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

end if
end subroutine gather_real8to4_nod2D
!==============================================
!============================================================================
subroutine gather_real8to4_elem3D(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(real64) ::  arr3D(:,:)
real(real32)  ::  arr3D_global(:,:)
real(real32), allocatable :: recvbuf(:,:)
real(real32), allocatable :: sendbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D


 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl1
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                     = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   allocate(sendbuf(nl1,myDim_elem2D))
   sendbuf(1:nl1,1:myDim_elem2D) = arr3D(1:nl1,1:myDim_elem2D)
   
   call MPI_SEND(sendbuf, myDim_elem2D*nl1, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   deallocate(sendbuf)
ENDIF

end if
end subroutine gather_real8to4_elem3D
!==============================================
subroutine gather_real8to4_elem2D(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real64) ::  arr2D(:)
real(real32)  ::  arr2D_global(:)
real(real32), allocatable :: recvbuf(:)
real(real32)  :: sendbuf(myDim_elem2D)
integer        :: req(npes-1)
integer        :: start, e2D


 if (npes> 1) then

CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   sendbuf(1:myDim_elem2D) = real(arr2D(1:myDim_elem2D),real32)
   call MPI_SEND(sendbuf, myDim_elem2D, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

end if
end subroutine gather_real8to4_elem2D
!==============================================
subroutine gather_elem2D_i(arr2D, arr2D_global)
! Make element information available to master PE 
  use g_PARSUP
  use o_MESH
  IMPLICIT NONE

  integer                       :: n
  integer                       :: arr2D(:)
  integer                       :: arr2D_global(:)
  integer, allocatable          :: recvbuf(:)
  integer                       :: req(npes-1)
  integer                       :: start, e2D
  CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
  ! Consider MPI-datatypes to recv directly into arr2D_global!
  IF ( mype == 0 ) THEN
     if (npes > 1) then
        allocate(recvbuf(remPtr_elem2D(npes)))
        do  n = 1, npes-1
            e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
            start = remPtr_elem2D(n)
            call MPI_IRECV(recvbuf(start), e2D, MPI_INTEGER, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
        enddo      
        arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
        call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
        arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                          = recvbuf(1 : remPtr_elem2D(npes)-1)
        deallocate(recvbuf)
     else
        arr2D_global(:) = arr2D(:)
     endif
  ELSE
     call MPI_SEND(arr2D, myDim_elem2D, MPI_INTEGER, 0, 2, MPI_COMM_FESOM, MPIerr )
  ENDIF
end subroutine gather_elem2D_i
!============================================================================
subroutine gather_nod2D_i(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer              :: n
integer              :: arr2D(:)
integer              :: arr2D_global(:)
integer, allocatable :: recvbuf(:)
integer              :: req(npes-1)
integer              :: start, n2D

if (npes> 1) then

CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global, 1)))
      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_INTEGER, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_INTEGER, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

endif
end subroutine gather_nod2D_i
!============================================================================
!
subroutine gather_edg2D(arr2D, arr2Dglobal)
! A 2D version of the previous routine
use g_PARSUP
USE o_MESH
IMPLICIT NONE

real(real64) ::  arr2D(:)
real(real64) ::  arr2Dglobal(:)

integer                                  ::  i, n, buf_size, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:)       ::  ibuf
REAL(real64), ALLOCATABLE, DIMENSION(:) ::  rbuf

IF ( mype == 0 ) THEN
    arr2Dglobal(myList_edge2D(1:myDim_edge2D))=arr2D(1:myDim_edge2D)
    DO  n = 1, npes-1
       CALL MPI_RECV( buf_size, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                      0, MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(rbuf(buf_size), ibuf(buf_size))

       CALL MPI_RECV(ibuf(1), buf_size, MPI_INTEGER, sender, &
                      1, MPI_COMM_FESOM, status, MPIerr )

       CALL MPI_RECV(rbuf(1), buf_size, MPI_DOUBLE_PRECISION, sender, &
                      2, MPI_COMM_FESOM, status, MPIerr )
       arr2Dglobal(ibuf)=rbuf
       DEALLOCATE(ibuf, rbuf)
    ENDDO
ELSE
    CALL MPI_SEND( myDim_edge2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( myList_edge2D(1), myDim_edge2D, MPI_INTEGER, 0, 1, &
                   MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( arr2D(1), myDim_edge2D, MPI_DOUBLE_PRECISION, 0, 2,&
                   MPI_COMM_FESOM, MPIerr )
ENDIF
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
end subroutine gather_edg2D
!
!============================================================================
!
subroutine gather_edg2D_i(arr2D, arr2Dglobal)
! A 2D version of the previous routine
use g_PARSUP
USE o_MESH
IMPLICIT NONE

integer  ::  arr2D(:)
integer  ::  arr2Dglobal(:)

integer                                  ::  i, n, buf_size, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:)       ::  ibuf, vbuf

IF ( mype == 0 ) THEN
    arr2Dglobal(myList_edge2D(1:myDim_edge2D))=arr2D(1:myDim_edge2D)
    DO  n = 1, npes-1
       CALL MPI_RECV( buf_size, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                      0, MPI_COMM_FESOM, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(ibuf(buf_size), vbuf(buf_size))

       CALL MPI_RECV(ibuf(1), buf_size, MPI_INTEGER, sender, &
                      1, MPI_COMM_FESOM, status, MPIerr )

       CALL MPI_RECV(vbuf(1), buf_size, MPI_INTEGER, sender, &
                      2, MPI_COMM_FESOM, status, MPIerr )
       arr2Dglobal(ibuf)=vbuf
       DEALLOCATE(ibuf, vbuf)
    ENDDO
ELSE
    CALL MPI_SEND( myDim_edge2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( myList_edge2D(1), myDim_edge2D, MPI_INTEGER, 0, 1, &
                   MPI_COMM_FESOM, MPIerr )
    CALL MPI_SEND( arr2D(1), myDim_edge2D, MPI_INTEGER, 0, 2,&
                   MPI_COMM_FESOM, MPIerr )
ENDIF
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)
end subroutine gather_edg2D_i
!==============================================

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

interface gather_edge
      module procedure gather_edg2D
      module procedure gather_edg2D_i
end interface gather_edge


private  ! hides items not listed on public statement 
public :: exchange_nod,exchange_elem,broadcast_nod,broadcast_elem, &
          gather_nod, gather_elem, exchange_nod_begin, exchange_nod_end, exchange_elem_begin, &
          exchange_elem_end, gather_edge
end module g_comm_auto
