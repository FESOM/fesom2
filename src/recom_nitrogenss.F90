subroutine recom_nitogenss(mesh)

  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use recom_config
  use REcoM_ciso

   use g_clock
  use o_PARAM
  use g_PARSUP
  use g_rotate_grid
  use g_config
  use mod_MESH
  use i_arrays 		! a_ice, m_ice 
  use o_param           ! num_tracers
  use i_param
  use o_arrays
  use g_forcing_arrays  ! press_air
  use g_comm_auto
  use i_therm_param
  use g_comm
  use g_support
  implicit none

  type(t_mesh), intent(in) , target :: mesh

  Integer                           :: node, nz,  id, n,  nl1, ul1, k, nlevels_nod2D_minimum
  Real(kind=8)                      :: wflux(mesh%nl)
  Real(kind=WP)                     :: tv, Fc
  real(kind=WP)                     :: Vben(mesh%nl),  aux(mesh%nl), LocDenit(myDim_nod2D)
!

#include "../associate_mesh.h"

if (.not. NitrogenSS) DenitBen = 0.0_WP
if (.not. NitrogenSS) return


if (mype==0 .and. (mod(mstep,100)==0))  print *, achar(27)//'[37m'//'         --> recom_nitogenss'//achar(27)//'[0m'
   aux=0._WP
   Vben=0._WP
   LocDenit=0._WP

   do n=1, myDim_nod2D ! needs exchange_nod in the end
      nl1=nlevels_nod2D(n)-1
      ul1=ulevels_nod2D(n)

      Vben = VDet
      if (allow_var_sinking) Vben =( Vdet_a * abs(zbar_3d_n(:,n)) + VDet) !--- [m/d] ---

      k=nod_in_elem2D_num(n)
      ! Screening minimum depth in neigbouring nodes around node n
      nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

      do nz=nlevels_nod2D_minimum, nl1
         tv = tr_arr(nz,n,10)*Vben(nz) !idetc [mmolC m-2 d-1]
         Fc = max(tiny, 0.1d0 * tv)                                            ! Fc= Liable carbon flux, Conversion of detC from [mmolC/m2/d] => [umolC/cm2/d]        
         aux(nz) = -0.9543 + 0.7662 * log(Fc) - 0.2350 * log(Fc) * log(Fc)     ! Denitrification factor following Middelburg et al. (1996)
         aux(nz) = exp(aux(nz) * log(10.d0))                                   ! = 10^Denit, following Middelburg et al. (1996)
         aux(nz) = - q_NC_Denit * aux(nz) * (area(nz,n)-area(nz+1,n)) * 10.d0  ! Conversion of detC from [umolC/cm2/d] => [mmolC/m2/d]
      end do

      nz=nl1
      tv = tr_arr(nz,n,10)*Vben(nz) !idetc
      Fc = max(tiny,0.1d0 * tv)    
      aux(nz) = -0.9543 + 0.7662 * log(Fc) - 0.2350 * log(Fc) * log(Fc) ! Denitrification factor following Middelburg et al. (1996)
      aux(nz) = exp(aux(nz) * log(10.d0))                              ! = 10^Denit, following Middelburg et al. (1996) [umolC/cm2/d] 
      aux(nz) = - q_NC_Denit * aux(nz) *(area(nz+1,n)) * 10.d0  ![mmolN/d]    

      ! nss array will be applied to tracer_id(tr_num)==1001 !idin
      do nz=ul1,nl1
         nss(nz,n) = nss(nz,n) + (aux(nz)/SecondsPerDay)*dt/area(nz,n)/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n)) ! we remove it from idin
         LocDenit(n) = LocDenit(n) - (aux(nz)/SecondsPerDay)/area(nz,n) !(aux(nz))*dt
      end do

      DenitBen(n)= DenitBen(n) + LocDenit(n)


! Below is the code from recom
!    wFluxDet(2)    = Vben_det * state(Nn,idetc)
!    Fc             = 0.1d0 * wFluxDet(2)                                     ! Conversion of detC from [mmolC/m2/day] => [umol/cm2/day]
!    if (NitrogenSS) then
      !LocDenit       = -0.9543 + 0.7662 * log(Fc) - 0.2350 * log(Fc) * log(Fc) ! Denitrification factor following Middelburg et al. (1996)
      !LocDenit       = exp(LocDenit * log(10.d0))                                 ! = 10^Denit, following Middelburg et al. (1996)
      !LocDenit       = q_NC_Denit * LocDenit * wFluxDet(2)
      !sms(Nn,idin)   = sms(Nn,idin)    -   LocDenit * dt * recipdzF(Nn)   ! Part of benthos din is removed and lost due to denitrification
!      LocDenit = zero
!    else
!      LocDenit = zero
!    end if                     

   end do

   call exchange_nod(DenitBen)

end subroutine recom_nitogenss



