module mod_healpix
use o_param
contains

subroutine regrid_to_healpix(arr_in, a_healpix, w_healpix, nside, partit, mesh)
    use g_config
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    implicit none

    type(t_mesh),   intent(in),    target        :: mesh
    type(t_partit), intent(inout), target        :: partit
    real(kind=WP), intent(in)                    :: arr_in(partit%myDim_nod2D)
    real(kind=WP), intent(inout), allocatable    :: a_healpix(:), w_healpix(:)
    real(kind=WP)                                :: w_loc(partit%myDim_nod2D), a_loc(partit%myDim_nod2D)
    integer                                      :: i_loc(partit%myDim_nod2D)
    real(kind=WP), allocatable                   :: recvbuf_a(:)
    real(kind=WP), allocatable                   :: recvbuf_w(:)
    integer,       allocatable                   :: recvbuf_i(:)
    real(kind=WP), allocatable                   :: w_glo(:), a_glo(:)
    integer, intent(in)                          :: nside
    integer                                      :: req(3)
    integer                                      :: indices(partit%myDim_nod2D)
    integer                                      :: i, ii, n, nn, myDim_nod2D_max
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"


    call MPI_AllREDUCE(partit%myDim_nod2D, myDim_nod2D_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    if (mype==0) then
        if (.NOT. allocated(a_healpix) ) allocate(a_healpix(12*nside**2))
        if (.NOT. allocated(w_healpix) ) allocate(w_healpix(12*nside**2))
        allocate(recvbuf_i(myDim_nod2D_max), recvbuf_a(myDim_nod2D_max), recvbuf_w(myDim_nod2D_max))
        a_healpix=0.0_WP
        w_healpix=0.0_WP
    end if

    w_loc =0.0_WP
    a_loc =0.0_WP
    i_loc =0

    do n=1, myDim_nod2D
       call ang2pix_ring(nside, 0.5_WP*pi - geo_coord_nod2D(2,n), geo_coord_nod2D(1,n), i_loc(n))
       w_loc(n)=areasvol(ulevels_nod2D(n),n)
       a_loc(n)=arr_in(n)*w_loc(n)
    end do

    IF (mype == 0 ) THEN
        do n = 1, npes-1
            nn   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
           call MPI_IRECV(recvbuf_i, nn, MPI_INTEGER, n, 1, MPI_COMM_FESOM, req(1), MPIerr)
           call MPI_IRECV(recvbuf_a, nn, MPI_REAL,    n, 2, MPI_COMM_FESOM, req(2), MPIerr)
           call MPI_IRECV(recvbuf_w, nn, MPI_REAL,    n, 3, MPI_COMM_FESOM, req(3), MPIerr)
           call MPI_WAITALL(3, req, MPI_STATUSES_IGNORE, MPIerr)
           do i=1, nn
              ii=recvbuf_i(i)
              a_healpix(ii)=a_healpix(ii)+recvbuf_a(i)
              w_healpix(ii)=w_healpix(ii)+recvbuf_w(i)
           end do
        enddo
        do i=1, myDim_nod2D
           ii=i_loc(i)
           a_healpix(ii)=a_healpix(ii)+a_loc(i)
           w_healpix(ii)=w_healpix(ii)+w_loc(i)
        end do
        deallocate(recvbuf_i, recvbuf_a, recvbuf_w)
        where (w_healpix > 0.0_WP)
              a_healpix=a_healpix/w_healpix
        end where
    ELSE
        call MPI_SEND(i_loc, myDim_nod2D, MPI_INTEGER, 0, 1, MPI_COMM_FESOM, MPIerr )
        call MPI_SEND(a_loc, myDim_nod2D, MPI_REAL,    0, 2, MPI_COMM_FESOM, MPIerr )
        call MPI_SEND(w_loc, myDim_nod2D, MPI_REAL,    0, 3, MPI_COMM_FESOM, MPIerr )
    ENDIF
end subroutine regrid_to_healpix

subroutine ang2pix_ring(nside, theta, phi, ipix)
    implicit none
    integer, intent(in) :: nside
    real(8), intent(in) :: theta, phi
    integer, intent(out) :: ipix
    real(8) :: z, z0, z1, z2, az, za, tt
    integer :: jp, jm, ifp, ifm, ip, jm1, jm2, irn, kshift

    ! Check if nside is a power of 2
    if (mod(nside, 2) /= 0 .or. nside <= 0) then
      print *, "Invalid value of nside. It must be a power of 2."
      ipix = -1
      return
    endif

    ! Calculate z coordinate
    z = cos(theta)
  
    ! Compute pixel number in the phi direction
    az = mod(phi, 2 * pi)
    za = acos(z)
    tt = az / (2 * pi) * real(nside)
  
    ! Determine the pixel number
    jp = nside * (0.5 + tt - 0.5)
    jm = nside * (0.5 + tt + 0.5)
  
    ! Compute the shift for the mod function
    ifp = int(jp)
    ifm = int(jm)
    ip = mod(ifp, nside)
    jm1 = mod(ifm, nside)
    jm2 = mod(ifm - 1, nside)
  
    ! Find the ring number
    irn = nside + min(jp - ifp, jm - ifm)
  
    ! Calculate the pixel index
    kshift = mod(irn, 2)
    ipix = (irn * (irn - 1)) / 2 + ip + kshift * nside * (nside - 1) / 2 + 1
  
  end subroutine ang2pix_ring

end module mod_healpix