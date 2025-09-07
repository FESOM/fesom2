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
subroutine cal_rejected_salt(ice, partit, mesh)
use o_arrays
USE MOD_ICE
use mod_mesh
USE MOD_PARTIT
USE MOD_PARSUP
use g_comm_auto
use o_tracers
use g_config,         only: dt
implicit none

integer         :: row
real(kind=WP)   :: aux
type(t_ice)   , intent(in), target :: ice
type(t_mesh)  , intent(in), target :: mesh
type(t_partit), intent(in), target :: partit
real(kind=WP), dimension(:)  , pointer :: thdgr, S_oc_array
real(kind=WP)                , pointer :: rhoice, rhowat, Sice
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
thdgr      => ice%thermo%thdgr
S_oc_array => ice%srfoce_salt
rhoice     => ice%thermo%rhoice
rhowat     => ice%thermo%rhowat
Sice       => ice%thermo%Sice

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
subroutine app_rejected_salt(ttf, partit, mesh)
  use o_arrays
  use mod_mesh
  USE MOD_PARTIT
  USE MOD_PARSUP
  use o_tracers
  use g_comm_auto
  implicit none

  integer         :: row, k, nod, nup, nlo, kml, nzmin, nzmax
  real(kind=WP)    :: zsurf, rhosurf, drhodz, spar(100)=0.0_WP

  integer         :: n_distr
  real(kind=WP)    :: drhodz_cri, rho_cri
  data drhodz_cri /0.01_WP/  !kg/m3/m  !NH   !Nguyen2011
  data n_distr /5/
  data rho_cri /0.4_WP/      !kg/m3    !SH   !Duffy1999

  type(t_mesh),   intent(in), target :: mesh
  type(t_partit), intent(in), target :: partit
  real(kind=WP),  intent (inout)     :: ttf(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  do row=1,myDim_nod2d+eDim_nod2D   ! myDim is sufficient
     if (ulevels_nod2D(row)>1) cycle
     if (ice_rejected_salt(row)<=0.0_WP) cycle
     ! do not parameterize brine rejection in regions with low salinity
     ! 1. it leads to further decrease of SSS
     ! 2. in case of non zero salinity of ice (the well accepted value is 5psu) the SSS might become negative
     nzmin = ulevels_nod2D(row)
     nzmax = nlevels_nod2D(row)
     if (ttf(nzmin,row) < 10.0_WP) cycle
     !if (geo_coord_nod2D(2,row)>0.0_WP) then  !NH
        kml=1
        spar(nzmin)=0.0_WP
        do k=nzmin, nzmax-1
           drhodz=bvfreq(k, row)*density_0/g
           if (drhodz>=drhodz_cri .or. Z_3d_n(k,row)<-80.0_WP) exit
           kml=kml+1
           spar(k+1)=area(k+1,row)*hnode(k+1,row)*(Z_3d_n(1,row)-Z_3d_n(k+1,row))**n_distr
        end do
        
        !_______________________________________________________________________
        !PS Here make sure that kml<=nzmax-1, otherwise you will write 
        !PS ttf(k,row)=ttf(k,row)+... into the bottom topography, which will cause 
        !PS a division by zero, due to areasvol(kml,row)=0.0_WP wherever is bottom
        !PS and can trigger the Nan Checker
        kml = min(kml, nzmax-1)
        
        !_______________________________________________________________________
        if (kml>nzmin) then
            ttf(nzmin,row)=ttf(nzmin,row)-ice_rejected_salt(row)/areasvol(1,row)/hnode(1,row)
            spar(nzmin+1:kml)=spar(nzmin+1:kml)/sum(spar(nzmin+1:kml))
            do k=nzmin+1,kml
                ttf(k,row)=ttf(k,row)+ice_rejected_salt(row)*spar(k)/areasvol(k,row)/hnode(k,row)
            end do
        endif
      !endif
  end do

end subroutine app_rejected_salt
