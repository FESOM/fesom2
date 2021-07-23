! cal_rejected_salt:
! calculates the total salt released via brine rejection during ice formation [unit: psu m3]
! app_rejected_salt:
! apply the subgrid scale salt plume parameterization by
! distributing salt released through brine rejection inside the mixed layer
! accorcing to the specified vertical distribution function.
! Currently only support the linear free surface option!
! Ref: Duffy1997, Duffy1999, Nguyen2009
! Originaly coded by Qiang Wang in FESOM 1.4 
!--------------------------------------------------------
subroutine cal_rejected_salt(mesh)
use g_parsup
use o_arrays
use mod_mesh
use g_comm_auto
use o_tracers
use g_forcing_arrays, only: thdgr
use i_ARRAYS,         only: S_oc_array
use i_therm_param,    only: rhoice, rhowat, Sice
use g_config,         only: dt
implicit none

integer         :: row
real(kind=WP)   :: aux
type(t_mesh), intent(in), target :: mesh

#include  "associate_mesh.h"

aux=rhoice/rhowat*dt
do row=1, myDim_nod2d +eDim_nod2D! myDim is sufficient
   if (thdgr(row)>0.0_WP .and. ulevels_nod2D(row)==1) then
      ice_rejected_salt(row)= &
      (S_oc_array(row)-Sice)*thdgr(row)*aux*area(1, row)
      !unit: psu m3
   else
      ice_rejected_salt(row)=0.0_WP
   end if
end do

end subroutine cal_rejected_salt
!
!----------------------------------------------------------------------------
!
subroutine app_rejected_salt(mesh)
  use g_parsup
  use o_arrays
  use mod_mesh
  use g_comm_auto
  use o_tracers
  implicit none

  integer         :: row, k, nod, nup, nlo, kml, nzmin, nzmax
  real(kind=WP)    :: zsurf, rhosurf, drhodz, spar(100)

  integer         :: n_distr
  real(kind=WP)    :: drhodz_cri, rho_cri
  data drhodz_cri /0.01_WP/  !kg/m3/m  !NH   !Nguyen2011
  data n_distr /5/
  data rho_cri /0.4_WP/      !kg/m3    !SH   !Duffy1999

  type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

  do row=1,myDim_nod2d+eDim_nod2D   ! myDim is sufficient
     if (ulevels_nod2D(row)>1) cycle
     if (ice_rejected_salt(row)<=0.0_WP) cycle
     ! do not parameterize brine rejection in regions with low salinity
     ! 1. it leads to further decrease of SSS
     ! 2. in case of non zero salinity of ice (the well accepted value is 5psu) the SSS might become negative
     nzmin = ulevels_nod2D(row)
     nzmax = nlevels_nod2D(row)
     !!PS if (tr_arr(1,row,2) < 10.0_WP) cycle
     if (tr_arr(nzmin,row,2) < 10.0_WP) cycle
     if (geo_coord_nod2D(2,row)>0.0_WP) then  !NH
        kml=1
        !!PS spar(1)=0.0_WP
        spar(nzmin)=0.0_WP
        
        !!PS do k=1, nlevels_nod2D(row)
        do k=nzmin, nzmax
           drhodz=bvfreq(k, row)*density_0/g
           if (drhodz>=drhodz_cri .or. Z_3d_n(k,row)<-50.0_WP) exit
           kml=kml+1
           spar(k+1)=area(k+1,row)*hnode(k+1,row)*(Z_3d_n(1,row)-Z_3d_n(k+1,row))**n_distr
        end do

        !!PS if (kml>1) then
        !!PS    tr_arr(1,row,2)=tr_arr(1,row,2)-ice_rejected_salt(row)/area(1,row)/hnode(1,row)
        !!PS    spar(2:kml)=spar(2:kml)/sum(spar(2:kml))
        !!PS    do k=2,kml
        !!PS       tr_arr(k,row,2)=tr_arr(k,row,2)+ice_rejected_salt(row)*spar(k)/area(k,row)/hnode(k,row)
        !!PS    end do
        !!PS endif
        if (kml>nzmin) then
           tr_arr(nzmin,row,2)=tr_arr(nzmin,row,2)-ice_rejected_salt(row)/areasvol(1,row)/hnode(1,row)
           spar(nzmin+1:kml)=spar(nzmin+1:kml)/sum(spar(nzmin+1:kml))
           do k=nzmin+1,kml
              tr_arr(k,row,2)=tr_arr(k,row,2)+ice_rejected_salt(row)*spar(k)/areasvol(k,row)/hnode(k,row)
           end do
        endif
      endif
  end do

end subroutine app_rejected_salt
