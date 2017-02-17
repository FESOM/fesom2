! Routines for reading inital or re-start files
! Adapted from FESOM, changes involve: (i) vertical-horizontal storage
! (ii) horizontal velocities on elements
! (iii) T, S at midlevels, w at full levels
! SD, 18.04.2012
module g_input
contains
!oce_input
!age_tracer_input
!passive_tracer_input
!ice_input
!read_prepared_initial_ice
!read_init_ts
!ice_input_unformatted
!oce_input_unformatted
subroutine oce_input
  ! read restart fields for ocean dynamics and active tracer variables
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------
  
  use o_param
  use o_mesh
  use o_arrays
  use g_clock
  use g_config
  use g_PARSUP
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j, dimid_rec, nrec
  integer                   :: ssh_varid
  integer                   :: tra_varid(num_tracer), tra_varid_ab(num_tracer)
  integer                   :: u_varid, v_varid
  integer                   :: w_varid, we_varid, wi_varid
  integer                   :: urhs_varid_AB, vrhs_varid_AB
  integer                   :: istart(2), icount(2), n2
  integer                   :: istart3(3), icount3(3)
  character(100)            :: filename
  character(100)            :: trname
  real(kind=8), allocatable :: aux2(:), aux3(:,:) 
  real(kind=8)              :: lval, gval
  allocate(aux2(nod2D), aux3(nl-1,elem2D)) 
  ! Sizes: u, v (nl-1, elem_size)
  !        w (nl, nod_size)
  !        T, S, ... (nl-1, nod_size)
  !        u_ice, v_ice nod_size
  
  ! open files
  filename=trim(ResultPath)//runid//'.'//cyearold//'.oce.restart.nc'

  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! inquire variable id
  status=nf_inq_varid(ncid, 'ssh', ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'u', u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'v', v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'urhs_AB', urhs_varid_AB)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'vrhs_AB', vrhs_varid_AB)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'w', w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'w_expl', we_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'w_impl', wi_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  do j=1,num_tracers
     SELECT CASE (j)
       CASE(1)
         trname='temp'
       CASE(2)
         trname='salt'
       CASE DEFAULT
         write(trname,'(A3,i1)') 'ptr', j
     END SELECT
     status=nf_inq_varid(ncid, trim(trname), tra_varid(j))
     if (status .ne. nf_noerr) call handle_err(status)
     ! Adams–Bashforth part
     trname=trim(trname)//'_AB'
     status=nf_inq_varid(ncid, trim(trname), tra_varid_ab(j))
     if (status .ne. nf_noerr) call handle_err(status)
  end do

  ! read variables
  ! which record to read
  if(restartflag=='last') then
     status = nf_inq_dimid(ncid, 'T', dimid_rec)
     if(status .ne. nf_noerr) call handle_err(status)
     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
     if(status .ne. nf_noerr) call handle_err(status)
  else
     read(restartflag,'(i4)') nrec
  end if

  ! 2d fields
  istart=(/1,nrec/)
  icount=(/nod2d, 1/)
  status=nf_get_vara_double(ncid, ssh_varid, istart, icount, aux2)
  if (status .ne. nf_noerr) call handle_err(status)
  eta_n=aux2(myList_nod2D)
  hbar=eta_n
  hbar_old=eta_n

  ! 3d fields, velocities
  istart3=(/1,1,nrec/)
  icount3=(/nl-1, elem2D, 1/)
  n2=myDim_elem2D+eDim_elem2D

  status=nf_get_vara_double(ncid, u_varid, istart3, icount3, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  UV(1,:,1:n2)=aux3(:, myList_elem2D(1:n2))

  status=nf_get_vara_double(ncid, v_varid, istart3, icount3, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  UV(2,:,1:n2)=aux3(:, myList_elem2D(1:n2))
  
  ! 3d fields UV_rhs_AB
  status=nf_get_vara_double(ncid, urhs_varid_AB, istart3, icount3, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  UV_rhsAB(1,:,1:n2)=aux3(:, myList_elem2D(1:n2))

  status=nf_get_vara_double(ncid, vrhs_varid_AB, istart3, icount3, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  UV_rhsAB(2,:,1:n2)=aux3(:, myList_elem2D(1:n2))

  deallocate(aux3)
  allocate(aux3(nl,nod2D))
  icount3=(/nl, nod2D, 1/)
  n2=myDim_nod2D+eDim_nod2D

  status=nf_get_vara_double(ncid, w_varid, istart3, icount3, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  Wvel=aux3(:,myList_nod2D)

  status=nf_get_vara_double(ncid, we_varid, istart3, icount3, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  Wvel_e=aux3(:,myList_nod2D)

  status=nf_get_vara_double(ncid, wi_varid, istart3, icount3, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  Wvel_i=aux3(:,myList_nod2D)

  deallocate(aux3)
  allocate(aux3(nl-1,nod2D)) 
  icount3=(/nl-1, nod2D, 1/)

  !T,S and passive tracers
  do j=1,num_tracers
     status=nf_get_vara_double(ncid, tra_varid(j), istart3, icount3, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
     tr_arr(:,:,j)=aux3(:,myList_nod2D)
     ! Adams–Bashforth part
     status=nf_get_vara_double(ncid, tra_varid_ab(j), istart3, icount3, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
     tr_arr_old(:,:,j)=aux3(:,myList_nod2D)
  end do

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  deallocate(aux3, aux2)
 lval=0.0_WP
 gval=0.0_WP
 do n2=1, myDim_nod2D
    lval=lval+area(1, n2)*eta_n(n2)
 end do
 call MPI_AllREDUCE(lval, gval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
if (mype==0) write(*,*) 'The global integral of SSH after restart is:', gval
end subroutine oce_input
!
!-------------------------------------------------------------------------
!
subroutine ice_input
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------

  use o_mesh
  use i_arrays
  use g_config
  use g_clock
  use g_PARSUP
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, nrec
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: istart(2), icount(2),n2
  character(100)            :: filename
  real(kind=8), allocatable :: aux2(:)

  allocate(aux2(nod2D))  

  ! open files
  filename=trim(ResultPath)//runid//'.'//cyearold//'.ice.restart.nc'
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! inquire variable id
  status=nf_inq_varid(ncid, 'area', area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hice', hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'uice', uice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'vice', vice_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! read variables

  ! which record to read
  if(restartflag=='last') then
     status = nf_inq_dimid(ncid, 'T', dimid_rec)
     if(status .ne. nf_noerr) call handle_err(status)
     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
     if(status .ne. nf_noerr) call handle_err(status)
  else
     read(restartflag,'(i4)') nrec
  end if

  istart=(/1,nrec/)
  icount=(/nod2d, 1/)
  status=nf_get_vara_double(ncid, area_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  a_ice=aux2(myList_nod2D)     
  status=nf_get_vara_double(ncid, hice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_ice=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, hsnow_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_snow=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, uice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  u_ice=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, vice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  v_ice=aux2(myList_nod2D) 
  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  deallocate(aux2)   
end subroutine ice_input
!
!-------------------------------------------------------------------------
!
subroutine read_prepared_initial_ice
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------
  
  use o_mesh
  use i_arrays
  use g_config
  use g_clock
  use g_PARSUP
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, nrec
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: istart(2), icount(2)
  character(100)            :: filename
  real(kind=8), allocatable :: aux2(:)

  allocate(aux2(nod2D))  

  ! open files
  filename=trim(ResultPath)//runid//'.'//'initial_ice.nc'
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) then
     print*,'ERROR: CANNOT READ initial ice FILE CORRECTLY !'
     print*,'Error in opening netcdf file'//filename
     call par_ex(1) 
     stop
  endif
  
  ! inquire variable id
  status=nf_inq_varid(ncid, 'area', area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hice', hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  
  ! read variables

  ! which record to read
  if(restartflag=='last') then
     status = nf_inq_dimid(ncid, 'T', dimid_rec)
     if(status .ne. nf_noerr) call handle_err(status)
     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
     if(status .ne. nf_noerr) call handle_err(status)
  else
     read(restartflag,'(i4)') nrec
  end if

  istart=(/1,nrec/)
  icount=(/nod2d, 1/)
  status=nf_get_vara_double(ncid, area_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  a_ice=aux2(myList_nod2D)     
  status=nf_get_vara_double(ncid, hice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_ice=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, hsnow_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_snow=aux2(myList_nod2D)      
  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  deallocate(aux2)   

  !cutoff
  
  where(a_ice>1.0)
     a_ice=1.0
  end where

  where(a_ice<0.0)
     a_ice=0.0
  end where
  
  where(m_ice<0.0)
     m_ice=0.0
  end where

  where(m_snow<0.0)
     m_snow=0.0
  end where
  
end subroutine read_prepared_initial_ice
!
!-----------------------------------------------------------------------------
!
subroutine read_init_ts
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------

  use o_mesh
  use o_param
  use o_arrays
  use g_config
  use g_PARSUP
  implicit none
  !
  integer                     :: i, j, n, nz
  integer                     :: num_lat_reg, num_lon_reg, num_lay_reg
  real(kind=8)                :: pp, pr, tt, ss, lon, lat, tbott, sbott
  real(kind=8), external      :: ptheta
  real(kind=8), allocatable   :: lon_reg(:), lat_reg(:), lay_reg(:)
  real(kind=8), allocatable   :: raw_data(:,:,:)
  real(kind=8), allocatable   :: temp_x(:), temp_y(:)

  ! open global T/S data files
  open(19,file=trim(ClimateDataPath)//trim(OceClimaDataName), status='old')

  ! read reg. grid
  read(19,*) num_lon_reg, num_lat_reg, num_lay_reg
  allocate(lon_reg(num_lon_reg))
  allocate(lat_reg(num_lat_reg))
  allocate(lay_reg(num_lay_reg))
  read(19,*) lon_reg
  read(19,*) lat_reg
  read(19,*) lay_reg
  allocate(raw_data(num_lon_reg,num_lat_reg,num_lay_reg))

  ! model grid coordinates
  allocate(temp_x(myDim_nod2d+eDim_nod2D), temp_y(myDim_nod2d+eDim_nod2D)) 
  do n=1, myDim_nod2d+eDim_nod2D                    
     temp_x(n)=geo_coord_nod2D(1,n)/rad
     temp_y(n)=geo_coord_nod2D(2,n)/rad
     ! change lon range to [0 360]
     if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0_WP  
  end do

  ! read raw data and do interpolation
  do i=1, num_lay_reg
     do j=1, num_lat_reg
        read(19, *) raw_data(:,j,i)
     end do
  end do  
  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
       lon_reg, lat_reg, lay_reg, raw_data, nl-1, myDim_nod2D+eDim_nod2D, &
       temp_x, temp_y, Z, tr_arr(:,:,1))

  do i=1, num_lay_reg
     do j=1, num_lat_reg
        read(19, *) raw_data(:,j,i)       
     end do
  end do 

  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
       lon_reg, lat_reg, lay_reg,raw_data, nl-1, myDim_nod2d+eDim_nod2D, &
       temp_x, temp_y, Z, tr_arr(:,:,2))

  close(19) 

! This is recommended cutoff for the cold start initialization. Uncomment if the model explodes !!!
! where (tr_arr(:,:,2)<=28._WP)
!       tr_arr(:,:,2)=28._WP
! end where

  ! Convert in situ temperature into potential temperature
  pr=0.0_WP
  do n=1,myDim_nod2d+eDim_nod2D
     do nz=1, nlevels_nod2D(n)-1    
        tt=tr_arr(nz,n,1)
        ss=tr_arr(nz,n,2)
        pp=abs(Z(nz))
        tr_arr(nz,n,1)=ptheta(ss, tt, pp, pr)
     end do	
  end do
  deallocate(temp_y, temp_x, raw_data, lay_reg, lat_reg, lon_reg)
end subroutine read_init_ts
end module g_input
