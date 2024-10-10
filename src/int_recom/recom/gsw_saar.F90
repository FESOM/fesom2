!==========================================================================
elemental function gsw_saar (p, long, lat)
!==========================================================================
!
! Calculates the Absolute Salinity Anomaly Ratio, SAAR at a geographic point.
!
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_saar : Absolute Salinity Anomaly Ratio               [unitless]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_add_barrier, gsw_add_mean, gsw_util_indx

use gsw_mod_error_functions, only : gsw_error_code

use gsw_mod_saar_data

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: p, long, lat

real (r8) :: gsw_saar

integer :: indx0, indy0, indz0, i, j, k, nmean

real (r8), dimension(4) :: saar, saar_old
real (r8) :: p0_original, lon0_in, sa_upper, sa_lower, dlong, dlat
real (r8) :: r1, s1, t1, saar_mean, ndepth_max, p_tmp, long_tmp

character (*), parameter :: func_name = "gsw_saar"

if (lat .lt. -86.0_r8 .or. lat .gt. 90.0_r8) then
   gsw_saar = gsw_error_code(2,func_name)
   return
end if

long_tmp = long
if (long_tmp.lt.0_r8) long_tmp = long_tmp + 360.0_r8

dlong = longs_ref(2) - longs_ref(1)
dlat = lats_ref(2) - lats_ref(1)

indx0 = floor(1.0_r8 + (nx-1.0_r8)*(long_tmp-longs_ref(1)) / &
		(longs_ref(nx)-longs_ref(1)))
if (indx0.eq.nx) indx0 = nx-1

indy0 = floor(1.0_r8 + (ny-1.0_r8)*(lat-lats_ref(1))/(lats_ref(ny)-lats_ref(1)))
if (indy0.eq.ny) indy0 = ny-1

! Look for the maximum valid "ndepth_ref" value around our point.
! Note: invalid "ndepth_ref" values are NaNs (a hangover from the codes
! Matlab origins). These Nans will not satisfy the ".gt. 0.0" test below,
! so will not be recognised as valid values.
ndepth_max = -1.0_r8
do k = 1,4
   if (ndepth_ref(indy0+delj(k),indx0+deli(k)).gt.0.0_r8) &
      ndepth_max = max(ndepth_max,DBLE(ndepth_ref(indy0+delj(k),indx0+deli(k))))
end do

! If we are a long way from the ocean then there will be no valid "ndepth_ref"
! values near the point (ie. surrounded by NaNs) - so just return SAAR = 0.0
if (ndepth_max.eq.-1.0_r8) then
   gsw_saar = 0.0_r8 
   return
end if 

p0_original = p
p_tmp = p
if (p_tmp.gt.p_ref(int(ndepth_max))) p_tmp = p_ref(int(ndepth_max))
call gsw_util_indx(p_ref,nz,p_tmp,indz0)
    
r1 = (long_tmp-longs_ref(indx0))/(longs_ref(indx0+1)-longs_ref(indx0))
s1 = (lat-lats_ref(indy0))/(lats_ref(indy0+1)-lats_ref(indy0))
t1 = (p_tmp-p_ref(indz0))/(p_ref(indz0+1)-p_ref(indz0))

do k = 1,4
   saar(k) = saar_ref(indz0,indy0+delj(k),indx0+deli(k))
end do

if (  longs_pan(1).le.long_tmp .and. long_tmp.le.longs_pan(npan)-0.001_r8 .and.&
    lats_pan(npan).le.lat      .and.      lat.le.lats_pan(1)) then
  saar_old = saar
  call gsw_add_barrier(saar_old,long_tmp,lat,longs_ref(indx0), &
  		       lats_ref(indy0),dlong,dlat,saar)
else if (abs(sum(saar)).ge.1e10_r8) then
  saar_old = saar
  call gsw_add_mean(saar_old,saar)
end if

sa_upper = (1.0_r8-s1)*(saar(1) + r1*(saar(2)-saar(1))) + s1*(saar(4) + &
		r1*(saar(3)-saar(4)))

do k = 1,4
   saar(k) = saar_ref(indz0+1,indy0+delj(k),indx0+deli(k))
end do

if (  longs_pan(1).le.long_tmp .and. long_tmp.le.longs_pan(npan)-0.001_r8 .and.&
    lats_pan(npan).le.lat      .and.      lat.le.lats_pan(1)) then
   saar_old = saar
   call gsw_add_barrier(saar_old,long_tmp,lat,longs_ref(indx0), &
   			lats_ref(indy0),dlong,dlat,saar)
else if (abs(sum(saar)).ge.1e10_r8) then 
   saar_old = saar
   call gsw_add_mean(saar_old,saar)
end if

sa_lower = (1.0_r8-s1)*(saar(1) + r1*(saar(2)-saar(1))) + s1*(saar(4) + &
		r1*(saar(3)-saar(4)))
if (abs(sa_lower).ge.1e10_r8) sa_lower = sa_upper
gsw_saar = sa_upper + t1*(sa_lower-sa_upper)

if (abs(gsw_saar).ge.1e10_r8) gsw_saar = gsw_error_code(3,func_name)

return
end function

!--------------------------------------------------------------------------
