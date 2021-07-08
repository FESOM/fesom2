!==============================================================================
!
! Simple initialization, forcing and output, just for tests
! for ocean and ice.
! ============================================================================ 
! ============================================================================
subroutine initial_state_test(mesh)
  use MOD_MESH
  use o_ARRAYS
  use o_PARAM
  use g_PARSUP
  !
  implicit none
  integer                            :: elem, n, nz, elnodes(3)
  integer                            :: elevation, strat, wind, cooling, tperturb
  real(kind=WP)                      :: lon, lat, a, dst  
  real(kind=WP)                      :: minlat,maxlat,tt,rwidth 
  type(t_mesh), intent(in)           , target :: mesh

#include "associate_mesh.h"
 

! Now updated for the box mesh, it was originally designed for hex mesh.
! In that case, the southern boundary is 40, the northern 48.83, and 0:18 the
! longitudinal extent.



! Default values
 stress_surf=0.0
 tr_arr(:,:,1)=20.0_WP
 Tsurf=tr_arr(1,:,1)
 heat_flux=0.0_WP
 tr_arr(:,:,2)=35.0_WP
 Ssurf=tr_arr(1,:,2)
 water_flux=0.0_WP
 relax2clim=0.0
 
 elevation=0
 strat=1
 wind=1
 cooling=0
 tperturb=0
 surf_relax_T=0 !10.0/10.0/24.0/3600.
 surf_relax_S=0.
 
 
 ! Stratification
  if(strat==1) then
  DO n=1, myDim_nod2D+eDim_nod2D
     DO nz=1, nlevels_nod2D(n)-1 
   !  tr_arr(nz,n,1)=tr_arr(nz,n,1)- 8.2e-3*abs(Z(nz))
     tr_arr(nz,n,1)=tr_arr(nz,n,1)-0.95_WP*20*tanh(abs(Z(nz))/300)-abs(Z(nz))/2400.0_WP

     END DO
  END DO   
  end if

Tsurf=tr_arr(1,:,1)

  if (tperturb==0) then
  ! Temperature perturbation
  do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)
     lon=coord_nod2D(1,n)
     dst=sqrt((lat-37.5*rad)**2+(lon-4.5*rad)**2)
     if (dst>1.5*rad) cycle
     do nz=1, nlevels_nod2D(n)-1 
        tr_arr(nz,n,1)=tr_arr(nz,n,1)+0.1*exp(-(dst/(1.5*rad))**2)*sin(pi*abs(Z(nz))/1600)
     end do
  end do
  end if


  if (cooling==1) then
  ! Surface cooling
  do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)
     lon=coord_nod2D(1,n)
     dst=sqrt((lat-37.5*rad)**2+(lon-4.5*rad)**2)
     if (dst>3.7*rad) cycle
        Tsurf(n)=Tsurf(n)-1*exp(-(dst/(2.2*rad))**2)
  end do
  end if

#ifdef false  
  if (wind==1) then
  DO elem=1, myDim_elem2D
     elnodes=elem2d_nodes(:,elem)
     lat=sum(coord_nod2D(2,elnodes))/3.0_WP
     lon=sum(coord_nod2D(1,elnodes))/3.0_WP
     stress_surf(1,elem)=-0.2 *cos(pi*(lat-30.0*rad)/(15.0*rad)) !(8.83*rad))  
          ! 40 is the south boundary of the hex box
  END DO
  end if
#endif

  if (wind==1) then
  DO elem=1, myDim_elem2D
     elnodes=elem2d_nodes(:,elem)
     lat=sum(coord_nod2D(2,elnodes))/3.0_WP
     lon=sum(coord_nod2D(1,elnodes))/3.0_WP
      !stress_surf(1,elem)=0.1 *cos(pi*(lat-40.0*rad)/(1500000.0/r_earth))* &
      !     exp(-((lat-40.0*rad)/(1500000.0/r_earth))**2)   
          ! 40 is the center of domain
     stress_surf(1,elem)=0.1 *cos(pi*(lat-35.0*rad)/(1250000.0/r_earth))* &
           exp(-((lat-35.0*rad)/(1250000.0/r_earth))**2)* &
           (1.0_WP-0.5_WP*((lat-35.0*rad)/(1250000.0/r_earth)))
          ! 35 is the center of domain    
   END DO
  end if

  ! Fix for too low salinity
  where (tr_arr(:,:,2)<20.4) tr_arr(:,:,2)=20.4
end subroutine initial_state_test
! ====================================================================

subroutine initial_state_channel_test(mesh)
  use MOD_MESH
  use o_ARRAYS
  use o_PARAM
  use g_PARSUP
  use g_CONFIG
  !
  implicit none
  integer                            :: elem, n, nz, elnodes(3)
  integer                            :: strat, wind, elevation 
  real(kind=WP)                      :: lon, lat, a, dst 
  type(t_mesh), intent(in)           , target :: mesh

#include "associate_mesh.h"

  ! Default values
 stress_surf=0.0
 tr_arr(:,:,1)=20.0_WP
 Tsurf=tr_arr(1,:,1)
 heat_flux=0.0_WP
 tr_arr(:,:,2)=35.0_WP
 Ssurf=tr_arr(1,:,2)
 water_flux=0.0_WP
 
  strat=1
  wind=0
  elevation=0
  
  lat=30.0*rad
  if (strat==1) then 
  do n=1, myDim_nod2D+eDim_nod2D
     dst=coord_nod2D(2, n)-lat
     do nz=1, nlevels(n)-1
        tr_arr(nz,n,1)=25.-0.5e-5*r_earth*dst- 8.2e-3*abs(Z(nz))
!	tr_arr(nz,n,1)=(25.-0.5e-5*r_earth*dst)*exp(Z(nz)/800)
     end do
  end do
  end if
 
  if (wind==1) then
   DO elem=1, myDim_elem2D
     call elem_center(elem, lon, lat, mesh)
     stress_surf(1,elem)=-0.2 *cos(pi*(lat-30.0*rad)/(15*rad))  
      ! 40 is the south boundary of the box
  END DO
  end if
  
  Tsurf=tr_arr(1,:,1)
  Ssurf=tr_arr(1,:,2)
  Tclim=tr_arr(:,:,1)
  Sclim=tr_arr(:,:,2)
  
  ! small perturbation:
  if (strat==1) then 
  do n=1, myDim_nod2D+eDim_nod2D
     dst=coord_nod2D(2, n)-30.0*rad
     do nz=1, nlevels(n)-1
        tr_arr(nz,n,1)=tr_arr(nz,n,1)-0.2*sin(2*pi*dst/(15.0*rad))*sin(pi*Z(nz)/1500.0) &
	       *(sin(8*pi*coord_nod2D(1,n)/(20.0*rad))+ &
	       0.5*sin(3*pi*coord_nod2D(1,n)/(20.0*rad)))
     end do
  end do
  end if
  
  if(elevation==1) then 
  eta_n=0.01*(coord_nod2D(2,:)-30.0*rad)/(15.0*rad)
  end if
  
  ! relaxation to climatology:
  Do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)
     if(lat>43.5*rad) relax2clim(n)=clim_relax*(1.0-(45*rad-lat)/(1.5*rad))
     if(lat<31.5*rad) relax2clim(n)=clim_relax*(1.0+(30*rad-lat)/(1.5*rad))
  END DO
 return
 ! advection scheme tests 
   
   dst=45.0*rad-30.0*rad;
   DO n=1, myDim_nod2D+eDim_nod2D
   lat=coord_nod2D(2,n)-30.0*rad
   lon=coord_nod2D(1,n)
      eta_n(n)=(1000000./pi)*sin(pi*lat/dst)*sin(2*pi*lon/(20*rad))
      !eta_n(n)=(1000000./pi)*sin(pi*lat/dst)*sin(pi*lon/(10*rad))
   end do 
   
     
   Do n=1, myDim_elem2D
     UV(1,:,n)=-sum(gradient_sca(4:6,n)*eta_n(elem2D_nodes(:,n)))
     UV(2,:,n)=sum(gradient_sca(1:3,n)*eta_n(elem2D_nodes(:,n)))
   END DO
  
   
   !Do n=1, elem2D
   !call elem_center(n, lon, lat, mesh)
   !lat=lat-30.0*rad
   !UV(1,:,n)=-(20*rad/dst)*0.1*cos(pi*lat/dst)*sin(2*pi*lon/(20*rad))
   !UV(2,:,n)= 0.2*sin(pi*lat/dst)*cos(2*pi*lon/(20*rad))    
   !end do
   relax2clim=0.
   tr_arr(:,:,1)=20.0
   Tsurf=tr_arr(1,:,1)
   surf_relax_T=0.
   surf_relax_S=0.
   !U_n=-0.3
   !V_n=0.
! Temperature perturbation
  do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)-32.5*rad
     lon=coord_nod2D(1,n)-5.0*rad
     if (lon>cyclic_length/2) lon=lon-cyclic_length
     if (lon<-cyclic_length/2) lon=lon+cyclic_length
     dst=sqrt((lat)**2+(lon)**2)
     if (dst>1.5*rad) cycle
     do nz=1, nlevels_nod2D(n)-1 
        !if(abs(Z(nz)+500)<300) then
        tr_arr(nz,n,1)=tr_arr(nz,n,1)+1.0*cos(pi*dst/2.0/1.5/rad)       !exp(-(dst/(1.5*rad))**2)
	!end if
     end do
  end do 
end subroutine initial_state_channel_test
! ====================================================================
subroutine initial_state_channel_narrow_test(mesh)
  use MOD_MESH
  use o_ARRAYS
  use o_PARAM
  use g_PARSUP
  use g_CONFIG
  !
  implicit none
  integer                            :: elem, n, nz, elnodes(3)
  integer                            :: strat, wind, elevation 
  real(kind=WP)                      :: lon, lat, a, dst 
  type(t_mesh), intent(in)           , target :: mesh

#include "associate_mesh.h"

  ! Default values
 stress_surf=0.0
 tr_arr(:,:,1)=20.0_WP
 Tsurf=tr_arr(1,:,1)
 heat_flux=0.0_WP
 tr_arr(:,:,2)=35.0_WP
 Ssurf=tr_arr(1,:,2)
 water_flux=0.0_WP
 
  strat=1
  wind=0
  elevation=0
  
  lat=30.0*rad
  if (strat==1) then 
  do n=1, myDim_nod2D+eDim_nod2D
     dst=coord_nod2D(2, n)-lat
     do nz=1, nlevels_nod2D(n)-1
        tr_arr(nz,n,1)=25.-0.5e-5*r_earth*dst- 8.2e-3*abs(Z(nz))
!	tr_arr(nz,n,1)=(25.-0.5e-5*r_earth*dst)*exp(Z(nz)/800)
     end do
  end do
  end if
 
  if (wind==1) then
   DO elem=1, myDim_elem2D
     call elem_center(elem, lon, lat, mesh)
     stress_surf(1,elem)=-0.2 *cos(pi*(lat-30.0*rad)/(10*rad))  
      ! 40 is the south boundary of the box
  END DO
  end if
  
  Tsurf=tr_arr(1,:,1)
  Ssurf=tr_arr(1,:,2)
  Tclim=tr_arr(:,:,1)
  Sclim=tr_arr(:,:,2)
  
  ! small perturbation:
  if (strat==1) then 
  do n=1, myDim_nod2D+eDim_nod2D
     dst=coord_nod2D(2, n)-30.0*rad
     do nz=1, nlevels(n)-1
        tr_arr(nz,n,1)=tr_arr(nz,n,1)-0.1*sin(pi*dst/(10.0*rad))*sin(pi*Z(nz)/1600.0) &
	       *(sin(4*pi*coord_nod2D(1,n)/(10.0*rad))+0.5*sin(3*pi*coord_nod2D(1,n)/(10.0*rad)))
     end do
  end do
  end if
  
  if(elevation==1) then 
  eta_n=0.01*(coord_nod2D(2,:)-30.0*rad)/(10.0*rad)
  end if
  
  ! relaxation to climatology:
  Do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)
     if(lat>38.5*rad) relax2clim(n)=clim_relax*(1.0-(40*rad-lat)/(1.5*rad))
     if(lat<31.5*rad) relax2clim(n)=clim_relax*(1.0+(30*rad-lat)/(1.5*rad))
  END DO
!T_rhsAB=tr_arr(:,:,1)   in case upwind1
!S_rhsAB=tr_arr(:,:,2)
! Advection experiments:
return
   UV(1,:,:)=-0.3
   UV(2,:,:)=0.
   
   dst=maxval(coord_nod2D(2,:))-30.0*rad;
   DO n=1, myDim_nod2D+eDim_nod2D
   lat=coord_nod2D(2,n)-30.0*rad
   lon=coord_nod2D(1,n)
      eta_n(n)=(1000000./pi)*sin(pi*lat/dst)*sin(2*pi*lon/(10*rad))
      !eta_n(n)=(1000000./pi)*sin(pi*lat/dst)*sin(pi*lon/(10*rad))
   end do 
   
     
   Do n=1, myDim_elem2D
     UV(1,:,n)=-sum(gradient_sca(4:6,n)*eta_n(elem2D_nodes(:,n)))
     UV(2,:,n)=sum(gradient_sca(1:3,n)*eta_n(elem2D_nodes(:,n)))
   END DO
  
   
   Do n=1, myDim_elem2D
   call elem_center(n, lon, lat, mesh)
   lat=lat-30.0*rad
   UV(1,:,n)=-0.1*(dst/10.0/rad)*cos(pi*lat/dst)*sin(2*pi*lon/(10*rad))
   UV(2,:,n)= 0.2*sin(pi*lat/dst)*cos(2*pi*lon/(10*rad))    
   end do
   
   
   
   relax2clim=0.
   tr_arr(:,:,1)=20.0
  
! Temperature perturbation
  do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)-32.5*rad
     lon=coord_nod2D(1,n)-5.0*rad
     if (lon>cyclic_length/2) lon=lon-cyclic_length
     if (lon<-cyclic_length/2) lon=lon+cyclic_length
     dst=sqrt((lat)**2+(lon)**2)
     if (dst>1.5*rad) cycle
     do nz=1, nlevels_nod2D(n)-1 
        tr_arr(nz,n,1)=tr_arr(nz,n,1)+1.0*cos(pi*dst/2.0/1.5/rad)       !exp(-(dst/(1.5*rad))**2)
     end do
  end do 
end subroutine initial_state_channel_narrow_test
! ================================================================ 
subroutine init_fields_na_test(mesh)
  use MOD_MESH
  use o_PARAM
  use o_ARRAYS
  use g_PARSUP
  !
  implicit none
  integer                            :: n, nz, nd
  real(kind=WP)                      :: maxlat, minlat, rwidth, lat,lon 
  logical                            :: c_status
  real(kind=WP)                      :: p0, ss, tt,pr
  type(t_mesh), intent(in)           , target :: mesh

#include "associate_mesh.h"

  c_status = .false.

  ! ===================
  ! Fill the model fields with dummy values
  ! ===================
  
  ! Default values
 stress_surf=0.0
 tr_arr(:,:,1)=20.0_WP
 Tsurf=tr_arr(1,:,1)
 heat_flux=0.0_WP
 tr_arr(:,:,2)=35.0_WP
 Ssurf=tr_arr(1,:,2)
 water_flux=0.0_WP
 
  ! ===================
  ! Initialize T, S from files
  ! ===================

  !call get_TS_mean('gur', c_status)

  ! ===================
  ! If database contains in situ
  ! temperature, transform it to 
  ! potential temperature
  ! ===================
  if(c_status) then
     pr=0.
     do n=1,myDim_nod2D+eDim_nod2D 
        DO nz=1,nlevels_nod2D(n)-1
        tt=tr_arr(nz,n,1)
        ss=tr_arr(nz,n,2)
        p0=abs(Z(nz))
        call ptheta(ss, tt, p0, pr, tr_arr(nz,n,1))
	END DO
     end do
     write(*,*) 'In situ temperature is converted to potential temperature'
  end if
  Tclim=tr_arr(:,:,1)
  Sclim=tr_arr(:,:,2)
  do n=1, myDim_nod2D+eDim_nod2D
     Tsurf(n)=tr_arr(1,n,1)
     Ssurf(n)=tr_arr(1,n,2)
  end do

  ! ====================
  ! Specify where restoring to 
  ! climatology is applied
  ! ====================
  ! relaxation to climatology:
  maxlat=80.0*rad
  minlat=-28.0*rad
  rwidth=10.0*rad
  Do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)
     if(lat>maxlat-rwidth) then
     relax2clim(n)=clim_relax*(cos(pi*0.5*(maxlat-lat)/rwidth))**2
     end if
     if(lat<minlat+rwidth) then
     relax2clim(n)=clim_relax*(cos(pi*0.5*(lat-minlat)/rwidth))**2
     end if
  END DO
  ! These were northern and southern zones
  ! Mediterranean salinity tongue:
  maxlat=36*rad    ! lat. of center
  minlat=-4*rad    ! lon. of center
  rwidth= 6*rad
  Do n=1, myDim_nod2D+eDim_nod2D
     lon=coord_nod2D(1,n)
     lat=coord_nod2D(2,n)
     tt=sqrt((lon-minlat)**2+(lat-maxlat)**2)
     if(tt<rwidth) relax2clim(n)=clim_relax*(cos(pi*0.5*tt/rwidth))**2
  END DO
  ! N-E  corner
  maxlat=73*rad    ! lat. of center
  minlat=18*rad    ! lon. of center
  rwidth= 6*rad
  Do n=1, myDim_nod2D+eDim_nod2D
     lon=coord_nod2D(1,n)
     lat=coord_nod2D(2,n)
     tt=sqrt((lon-minlat)**2+(lat-maxlat)**2)
     if(tt<rwidth) relax2clim(n)=clim_relax*(cos(pi*0.5*tt/rwidth))**2
  END DO
   
   
     
  !  
  ! Fix for too low salinity
  where (tr_arr(:,:,2)<20.4)
     tr_arr(:,:,2)=20.4
  end where
end subroutine init_fields_na_test  
! ================================================================  
subroutine init_fields_global_test(mesh)
  use MOD_MESH
  use o_PARAM
  use o_ARRAYS
  use g_PARSUP
  !
  implicit none
  integer                            :: n, nz, nd
  real(kind=WP)                      :: maxlat, minlat, rwidth, lat, lon 
  logical                            :: c_status
  real(kind=WP)                      :: p0, ss, tt,pr
  type(t_mesh), intent(in)           , target :: mesh

#include "associate_mesh.h"

  c_status = .false.

  ! ===================
  ! Fill the model fields with dummy values
  ! ===================
  
  ! Default values
 stress_surf=0.0
 tr_arr(:,:,1)=20.0_WP
 Tsurf=tr_arr(1,:,1)
 heat_flux=0.0_WP
 tr_arr(:,:,2)=35.0_WP
 Ssurf=tr_arr(1,:,2)
 water_flux=0.0_WP
 
  ! ===================
  ! Initialize T, S from files
  ! ===================

  !call get_TS_mean('woa', c_status)

  ! ===================
  ! If database contains in situ
  ! temperature, transform it to 
  ! potential temperature
  ! ===================
  if(c_status) then
     pr=0.
     do n=1,myDim_nod2D+eDim_nod2D
        DO nz=1,nlevels_nod2D(n)-1
        tt=tr_arr(nz,n,1)
        ss=tr_arr(nz,n,2)
        p0=abs(Z(nz))
        call ptheta(ss, tt, p0, pr, tr_arr(nz,n,1))
	END DO
     end do
     write(*,*) 'In situ temperature is converted to potential temperature'
  end if
  Tclim=tr_arr(:,:,1)
  Sclim=tr_arr(:,:,2)
  do n=1, myDim_nod2D+eDim_nod2D
     Tsurf(n)=tr_arr(1,n,1)
     Ssurf(n)=tr_arr(1,n,2)
  end do

  ! ====================
  ! Specify where restoring to 
  ! climatology is applied
  ! ====================
  ! relaxation to climatology:
  maxlat=80.0*rad
  !minlat=minval(coord_nod2D(2,:))
  rwidth=20.0*rad
  Do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)
     if(lat>maxlat-rwidth) relax2clim(n)=clim_relax*(1.0-(maxlat-lat)/rwidth)
      !if(lat<minlat+rwidth) relax2clim(n)=clim_relax*(1.0+(minlat-lat)/rwidth)
  END DO
     
  !  
  ! Fix for too low salinity
  where (tr_arr(:,:,2)<32.)
     tr_arr(:,:,2)=32.
  end where
  ! Fix for Mediterranean Sea
  Do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)/rad
     lon=coord_nod2D(1,n)/rad
     if((lon>15.0).and.(lon<40.0).and.(lat>30.0).and.(lat<40.0)) then  
     DO nz=1,nlevels_nod2D(n)-1
        if(tr_arr(nz,n,2)<38.0)  tr_arr(nz,n,2)=38.0 
     END DO
     end if
  end do
end subroutine init_fields_global_test
! ================================================================ 
! ====================================================================

subroutine initial_state_channel_dima_test(mesh)
  use MOD_MESH
  use o_ARRAYS
  use o_PARAM
  use g_PARSUP
  !
  implicit none
  integer                            :: elem, n, nz, elnodes(3)
  integer                            :: strat, wind, elevation 
  real(kind=WP)                      :: lon, lat, a, dst 
  type(t_mesh), intent(in)           , target :: mesh

#include "associate_mesh.h"

! Default values
 stress_surf=0.0
 tr_arr(:,:,1)=20.0_WP
 Tsurf=tr_arr(1,:,1)
 heat_flux=0.0_WP
 tr_arr(:,:,2)=35.0_WP
 Ssurf=tr_arr(1,:,2)
 water_flux=0.0_WP
 
  strat=1
  wind=0
  elevation=0
  
  lat=30.0*rad
  if (strat==1) then 
  do n=1, myDim_nod2D+eDim_nod2D
     dst=coord_nod2D(2, n)-lat
     do nz=1, nlevels(n)-1
        tr_arr(nz,n,1)=25.-0.5e-5*r_earth*dst- 8.2e-3*abs(Z(nz))
!	tr_arr(nz,n,1)=(25.-0.5e-5*r_earth*dst)*exp(Z(nz)/800)
     end do
  end do
  end if
 
  if (wind==1) then
   DO elem=1, myDim_elem2D
     call elem_center(elem, lon, lat, mesh)
     stress_surf(1,elem)=-0.2 *cos(pi*(lat-30.0*rad)/(15*rad))  
      ! 40 is the south boundary of the box
  END DO
  end if
  
  Tsurf=tr_arr(1,:,1)
  Ssurf=tr_arr(1,:,2)
  Tclim=tr_arr(:,:,1)
  Sclim=tr_arr(:,:,2)
  
  ! small perturbation:
  if (strat==1) then 
  do n=1, myDim_nod2D+eDim_nod2D
     dst=coord_nod2D(2, n)-30.0*rad
     do nz=1, nlevels(n)-1
        tr_arr(nz,n,1)=tr_arr(nz,n,1)-0.1*sin(pi*dst/(15.0*rad))*sin(pi*Z(nz)/1500.0) &
	       *(sin(8*pi*coord_nod2D(1,n)/(40.0*rad))+sin(5*pi*coord_nod2D(1,n)/(40.0*rad)))
     end do
  end do
  end if
  
  if(elevation==1) then 
  eta_n=0.01*(coord_nod2D(2,:)-30.0*rad)/(15.0*rad)
  end if
  
  ! relaxation to climatology:
  Do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)
     if(lat>43.5*rad) relax2clim(n)=clim_relax*(1.0-(45*rad-lat)/(1.5*rad))
     if(lat<31.5*rad) relax2clim(n)=clim_relax*(1.0+(30*rad-lat)/(1.5*rad))
  END DO
end subroutine initial_state_channel_dima_test
! ====================================================================
subroutine ice_init_fields_test(mesh)
!
! Simple initialization for a box model to test the dynamical part.
! No thermodinamics is initialized here 
!
use mod_mesh
use i_arrays
use i_param
use o_param
use g_PARSUP
use o_ARRAYS
use g_CONFIG
use g_comm_auto

IMPLICIT NONE
real(kind=WP)              :: xmin, xmax, ymin, ymax, Lx, Ly, meanf
integer                    :: n, elnodes(3)   
type(t_mesh), intent(in)   , target :: mesh

#include "associate_mesh.h"

   
  coriolis=1.4e-4  ! redefines Coriolis
  coriolis_node=1.4e-4  
   ! Set initial thickness and area coverage:
   m_ice=2.0
   m_snow=0.0
   u_ice=0.0
   v_ice=0.0
   stress_atmice_x=0.0
   stress_atmice_y=0.0
   ! a_ice is defined later
   
   
   ! Set ocean velocity (stationary in time):
   xmin=0.0_WP*rad
   xmax=20.0_WP*rad    !10.0_WP*rad
   ymin=30._WP*rad     !30._WP*rad
   ymax=45.0_WP*rad    !40.0_WP*rad
   Lx=xmax-xmin
   Ly=ymax-ymin
   
   DO n=1, myDim_nod2D+eDim_nod2D
   a_ice(n)=(coord_nod2d(1,n)-xmin)/Lx      
   END DO
   
   DO n=1, myDim_nod2D+eDim_nod2D
   U_w(n)=0.1*(2*(coord_nod2d(2,n)-ymin)-Ly)/Ly
   V_w(n)=-0.1*(2*(coord_nod2d(1,n)-xmin)-Lx)/Lx
   END DO
   m_ice=m_ice*a_ice
  
   ! Elevation computed approximately, from the geostrophy:
   meanf= 1.4e-4*r_earth   !2*omega*sin(yc)*r_earth
   DO n=1, myDim_nod2d+eDim_nod2D 
      elevation(n)=-0.1*meanf/g *((coord_nod2d(2,n)-ymin)**2/Ly- &
                     (coord_nod2d(2,n)-ymin)+ &  
                    (coord_nod2d(1,n)-xmin)**2/Lx -&
		    (coord_nod2d(1,n)-xmin))
   END DO
end subroutine ice_init_fields_test
! =============================================================================
Subroutine ice_update_forcing_test(step, mesh)
!
! Here only simple wind variability is introduced
!
use mod_mesh
use i_arrays
use i_param
use o_param
use i_therm_param
use g_PARSUP
use g_forcing_arrays
USE g_CONFIG
IMPLICIT NONE
real(kind=WP)             :: xmin, xm, ym, ymin, Lx, Ly, td, cdwin
integer                   :: step, n, elnodes(3)    
type(t_mesh), intent(in)  , target :: mesh

#include "associate_mesh.h"

   cdwin=0.00225_WP
   ! Set wind velocity (stationary in time):
   xmin=0.0_WP*rad
   Lx=20.0_WP*rad-xmin
   ymin=30.0_WP*rad
   Ly=45.0_WP*rad-ymin
   td=4*3600*24.0_WP
   
   DO n=1, myDim_nod2D+eDim_nod2D
   xm=coord_nod2d(1,n)
   ym=coord_nod2d(2,n)
   u_wind(n)=5.0+(sin(2*pi*step*dt/td)-3.0)*sin(2*pi*(xm-xmin)/Lx) &
                *sin(pi*(ym-ymin)/Ly) 
		  
   v_wind(n)=5.0+(sin(2*pi*step*dt/td)-3.0)*sin(2*pi*(ym-ymin)/Ly) &
                *sin(pi*(xm-xmin)/Lx) 
   END DO
   ! wind to stress:
   
   stress_atmice_x = rhoair*cdwin*sqrt(u_wind**2+v_wind**2)*u_wind
   stress_atmice_y = rhoair*cdwin*sqrt(u_wind**2+v_wind**2)*v_wind
end subroutine ice_update_forcing_test
!
!==============================================================================
! Simple initialization for tests for GM with the real geometry
! ============================================================================ 
subroutine ini_global_ocean(mesh)
  use MOD_MESH
  use o_ARRAYS
  use o_PARAM
  use g_PARSUP
  USE g_ROTATE_grid
  !
  implicit none
  integer                            :: n, nz
  real(kind=WP)                      :: minlat,maxlat, lon, lat, val
  type(t_mesh), intent(in)           , target :: mesh

#include "associate_mesh.h"

 tr_arr(:,:,1)=20.0_WP
 tr_arr(:,:,2)=34.0_WP


 call r2g(lon, maxlat, coord_nod2D(1,1), coord_nod2D(2,1))
 call r2g(lon, minlat, coord_nod2D(1,1), coord_nod2D(2,1))
 DO n=2,myDim_nod2D+eDim_nod2D 
    call r2g(lon, lat, coord_nod2D(1,n), coord_nod2D(2,n))
    maxlat=max(maxlat, lat)
    minlat=min(minlat, lat)
 END DO

 call MPI_AllREDUCE(minlat, val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
 minlat=val
 call MPI_AllREDUCE(maxlat, val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
 maxlat=val

 ! Stratification
  DO n=1, myDim_nod2D+eDim_nod2D
     call r2g(lon, lat, coord_nod2D(1,n), coord_nod2D(2,n))
     DO nz=1, nlevels_nod2D(n)-1 
        tr_arr(nz,n,1)=tr_arr(nz,n,1)-(lat-minlat)/(maxlat-minlat)*2.0_WP
     END DO
  END DO
end subroutine ini_global_ocean
! ====================================================================
!
!==============================================================================
! Zero the dynamicsl variables and forcing to allow for debugging of new implementations
! ============================================================================ 
subroutine zero_dynamics
   use g_parsup
   use o_arrays
   use g_comm_auto
   use o_tracers
   use g_forcing_arrays
   implicit none

  water_flux    =0._WP
  real_salt_flux=0._WP
  surf_relax_S  =0._WP
  heat_flux     =0._WP
  UV            =0._WP
  Wvel          =0._WP
end subroutine zero_dynamics
! ====================================================================

