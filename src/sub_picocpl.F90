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
!
!--------------------------------------------------------------------------------------------
!
contains
!
!--------------------------------------------------------------------------------------------
!
! 1. find_coast_Antc shall be called once to look up the nodes at coast of Antarctic
! 2. allocate arrays
! 3. read PICO data
subroutine init_picocpl(mesh)
  IMPLICIT NONE
  type(t_mesh), intent(in)          , target :: mesh
  integer                                    :: i, ed, nd, el, k

#include "associate_mesh.h"
  allocate(coast_Antc(myDim_nod2D+eDim_nod2D))
  coast_Antc=.FALSE.

  DO ed=1, myDim_edge2D
     if (sum(geo_coord_nod2D(2, edges(:,ed)))/2.>-60.*rad) CYCLE
     if (edge_tri(2,ed)<=0) then
         el=edge_tri(1, ed)
         coast_Antc(edges(:,ed))=.TRUE.         
     end if
  END DO

  pico_N=360
  allocate(pico_lon(pico_N), pico_z_up(pico_N), pico_z_lo(pico_N), pico_data_t(pico_N), pico_data_s(pico_N))
  allocate(f2pico_data_t(pico_N), f2pico_data_s(pico_N), f2pico_count(pico_N))
  allocate(work_array(pico_N))
  pico_lon= (/(i, i=0, 359)/)
  pico_z_up=   0.
  pico_z_lo=-500.
  pico_data_t=-1.0
  pico_data_s=34.0

end subroutine init_picocpl
!
!--------------------------------------------------------------------------------------------
!
! FESOM 2 PICO
subroutine fesom2pico(mesh)
  IMPLICIT NONE
  type(t_mesh), intent(in), target :: mesh
  integer                          :: n, k, ipico
  real(kind=WP)                    :: x

#include "associate_mesh.h"
  f2pico_data_t=0.
  f2pico_data_s=0.
  f2pico_count =0.
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
           if ((zbar_3d_n(k,   n) <= pico_z_up(ipico)) .AND. &
               (zbar_3d_n(k+1, n) >= pico_z_lo(ipico))) then
               f2pico_data_t(ipico)=f2pico_data_t(ipico)+tr_arr(k, n, 1)
               f2pico_data_s(ipico)=f2pico_data_s(ipico)+tr_arr(k, n, 2)
               f2pico_count (ipico)=f2pico_count (ipico)+1.
           end if
        END DO
     END IF
  END DO
  where (f2pico_count>0.1)
         f2pico_data_t=f2pico_data_t/f2pico_count
         f2pico_data_s=f2pico_data_s/f2pico_count
  end where
  work_array=0.
  call MPI_AllREDUCE(f2pico_data_t , work_array , pico_N, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  f2pico_data_t=work_array

  work_array=0.
  call MPI_AllREDUCE(f2pico_data_s , work_array , pico_N, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  f2pico_data_s=work_array

  work_array=0.
  call MPI_AllREDUCE(f2pico_count , work_array , pico_N, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  f2pico_count=work_array
  if (mype==0) then
     do n=1, pico_N
        write(*,*) pico_lon(n), nint(f2pico_count(n)), f2pico_data_t(n), f2pico_data_s(n)
     end do
  end if
end subroutine fesom2pico

end module g_picocpl


