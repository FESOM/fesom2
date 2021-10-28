module g_picocpl
  USE mod_mesh
  use g_parsup
  use g_comm_auto
  use o_arrays
  implicit none
  private
  public :: init_picocpl, fesom2pico
  logical,       dimension(:),  allocatable  :: coast_Antc
  integer                                    :: pico_N
  real(kind=WP), dimension(:),  allocatable  :: work_array
  real(kind=WP), dimension(:),  allocatable  :: f2pico_data_t, f2pico_data_s, f2pico_count
  real(kind=WP), dimension(:),  allocatable  :: pico_data_t, pico_data_s, pico_z_up, pico_z_lo, pico_lon
  real(kind=WP)                              :: lvol, gvol
!
!--------------------------------------------------------------------------------------------
!
contains
!
!--------------------------------------------------------------------------------------------
!
! 1. look up the nodes at coast of Antarctic (only vertical walls will be compared with PICO)
! 2. allocate arrays
! 3. read PICO data
subroutine init_picocpl(mesh)
  IMPLICIT NONE
  type(t_mesh), intent(in)          , target :: mesh
  integer                                    :: i, ed, nd, el, k
  integer                                    :: fileID, ierror

#include "associate_mesh.h"
  ! this will contain nodes around Antarctic
  allocate(coast_Antc(myDim_nod2D+eDim_nod2D))
  coast_Antc=.FALSE.

  DO ed=1, myDim_edge2D
     if (sum(geo_coord_nod2D(2, edges(:,ed)))/2.>-60.*rad) CYCLE
     if (edge_tri(2,ed)<=0) then
         coast_Antc(edges(:,ed))=.TRUE.
     end if
  END DO

  ! do not read but create the synthetic PICO data
  ! 360 points (longitudes)
  !!
  ! read PICO data
  !!
  !360
  !0.   10. 500. -1. 34.
  !...
  !359. 10. 500. -1. 34.

  if (mype==0) then
     open(fileID, file='test.dat')
     read(fileID,*) pico_N
     write (*,*) 'start reading ', pico_N, ' data points!'
  end if

  call MPI_BCast(pico_N, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)  
  allocate(pico_lon(pico_N), pico_z_up(pico_N), pico_z_lo(pico_N), pico_data_t(pico_N), pico_data_s(pico_N))

  if (mype==0) then
     do i=1, pico_N
        read(fileID,*) pico_lon(i), pico_data_s(i), pico_data_t(i), pico_z_up(i), pico_z_lo(i)
     end do
     pico_data_t=pico_data_t-273.
  end if

  call MPI_BCast(pico_lon,    pico_N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(pico_data_s, pico_N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(pico_data_t, pico_N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(pico_z_up,   pico_N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(pico_z_lo,   pico_N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

  !!!
  ! modelled data interpolated onto PICo points
  allocate(f2pico_data_t(pico_N), f2pico_data_s(pico_N), f2pico_count(pico_N))
  ! some temporary array will be required
  allocate(work_array(pico_N))
end subroutine init_picocpl
!
!--------------------------------------------------------------------------------------------
!
! FESOM 2 PICO
! here we interpolate FESOM data onto PICO mesh
! one can do restoring to PICO as well (if required)
subroutine fesom2pico(mesh)
  IMPLICIT NONE
  type(t_mesh), intent(in), target :: mesh
  integer                          :: n, k, ipico
  real(kind=WP)                    :: x
  real(kind=WP)                    :: gammaT=1.e-4, gammaS=1.e-4
#include "associate_mesh.h"
  f2pico_data_t=0.
  f2pico_data_s=0.
  f2pico_count =0.
  lvol         =0.
  DO n=1, myDim_nod2D
     IF (coast_Antc(n)) then
        x=geo_coord_nod2D(1,n)
        work_array=(pico_lon-x/rad)
        where (work_array>=180.)
               work_array=work_array-360.
        end where
        where (work_array<=-180.)
               work_array=work_array+360.
        end where
        ipico=minloc(abs(work_array), 1)
        DO k=1, nlevels_nod2D(n)-1
           if (((zbar_3d_n(k,   n) <= pico_z_up(ipico)) .AND. & !FESOM segment within the PICO segment
               (zbar_3d_n(k+1, n) >= pico_z_lo(ipico))) .OR.  & !PICO segment within the FESOM segment
               ((zbar_3d_n(k,   n) <= pico_z_up(ipico)) .AND. &
               (zbar_3d_n(k+1, n) <= pico_z_lo(ipico)))) then
               if (pico_data_s(ipico) < 1.) CYCLE !SKIP 
               tr_arr(k, n, 1)=tr_arr(k, n, 1)+gammaT*(pico_data_t(ipico)-tr_arr(k, n, 1)) !apply nudging
               tr_arr(k, n, 2)=tr_arr(k, n, 2)+gammaS*(pico_data_s(ipico)-tr_arr(k, n, 2)) !apply nudging
               lvol=lvol+area(k,n)*hnode_new(k,n) !will the the total volume of 3D points where nuging is applied

               f2pico_data_t(ipico)=f2pico_data_t(ipico)+tr_arr(k, n, 1)
               f2pico_data_s(ipico)=f2pico_data_s(ipico)+tr_arr(k, n, 2)
               f2pico_count (ipico)=f2pico_count (ipico)+1.
           end if
        END DO
     END IF
  END DO

  !in case of nudging
  call exchange_nod(tr_arr(:,:,1))
  call exchange_nod(tr_arr(:,:,2))

  work_array=0.
  call MPI_AllREDUCE(f2pico_data_t , work_array , pico_N, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  f2pico_data_t=work_array

  work_array=0.
  call MPI_AllREDUCE(f2pico_data_s , work_array , pico_N, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  f2pico_data_s=work_array

  work_array=0.
  call MPI_AllREDUCE(f2pico_count , work_array , pico_N, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  f2pico_count=work_array

  call MPI_AllREDUCE(lvol, gvol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)

  where (f2pico_count>0.1)
         f2pico_data_t=f2pico_data_t/f2pico_count
         f2pico_data_s=f2pico_data_s/f2pico_count
  end where

  ! f2pico_data_t<->pico_data_t
  ! f2pico_data_s<->pico_data_s
  ! do comparison with PICO & output
  if (mype==0) then
     if (mstep==1) write(*,*) 'the 3D volume of PICO nudged area [m3] is : ', gvol
     if (mod(mstep, 10)==1) then
     OPEN(UNIT=2010, FILE='fesom.dat', POSITION='APPEND') 
     write(2010,*)'****************** ', mstep, ' ******************'
     do n=1, pico_N
        if (f2pico_count(n) > 0.1) then
           write(2010,*) INT(pico_lon(n)), nint(f2pico_count(n)), pico_data_t(n), f2pico_data_t(n), pico_data_s(n), f2pico_data_s(n)
        end if     
     end do
     CLOSE(2010)
     end if
  end if
end subroutine fesom2pico

end module g_picocpl


