!==========================================================================
pure subroutine gsw_util_indx (x, n, z, k)
!==========================================================================
!
!  Finds the index of the value in a monotonically increasing array
!
!  x	 :  array of monotonically increasing values
!  n     :  length of the array
!  z     :  value to be indexed
!
!  k      : index K :- if x(k) <= z < x(k+1), or
!               n-1 :- if z = x(n)
!--------------------------------------------------------------------------

use gsw_mod_kinds

integer, intent(in) :: n
integer, intent(out) :: k
real (r8), intent(in), dimension(n) :: x
real (r8), intent(in) :: z

integer :: ku, kl, km

if(z.gt.x(1).and.z.lt.x(n)) then
   kl=1
   ku=n
   do while (ku-kl.gt.1)
      km=0.5_r8*(ku+kl)
      if(z.gt.x(km)) then
         kl=km
      else
         ku=km
      endif
   end do
   k=kl
   if(z.eq.x(k+1)) then 
     k = k+1
   end if
elseif (z.le.x(1)) then
      k = 1
else	!if (z.ge.x(n)) then - removed (GBH 3/6/2015) so z=NaN has somewhere to go (otherwise k is undefined and gives segmentation fault)
      k = n-1
end if

return
end subroutine

!--------------------------------------------------------------------------

