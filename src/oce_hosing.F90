!
! Freshwater hosing experiments: impose a freshwater anomaly around Antarctica,
! either as a surface virtual salinity flux or distributed over a depth range.
!
! Both routines are opt-in and do nothing unless use_hosing is set in
! namelist.config. See &run_config: use_hosing, hosing_mode, hosing_hSv.
!
module hosing_interface
    interface
        subroutine fw_surf_anomaly(hSv, tracers, partit, mesh)
        use o_PARAM
        USE MOD_TRACER
        USE MOD_PARTIT
        use MOD_MESH
        type(t_tracer), intent(in),    target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        real(kind=WP) , intent(in)            :: hSv
        end subroutine fw_surf_anomaly

        subroutine fw_depth_anomaly(tts, ttt, hSv, partit, mesh)
        use o_PARAM
        USE MOD_PARTIT
        use MOD_MESH
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(inout), target :: partit
        real(kind=WP),  intent(in)            :: hSv
        real(kind=WP),  intent(inout)         :: tts(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP),  intent(inout)         :: ttt(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
        end subroutine fw_depth_anomaly
    end interface
end module hosing_interface

subroutine fw_surf_anomaly(hSv, tracers, partit, mesh)  !!!! subroutine to apply freshwater at the surface

    use o_PARAM
    use o_ARRAYS
    USE MOD_TRACER
    USE MOD_PARTIT
    use MOD_MESH
    use MOD_PARSUP
    USE g_CONFIG
    use g_comm_auto
    use g_support
    implicit none

    type(t_partit), intent(inout), target  :: partit
    type(t_mesh),   intent(in),    target  :: mesh
    type(t_tracer), intent(in),    target  :: tracers
    real(kind=WP),  intent(in)             :: hSv
    real(kind=WP)                          :: x, y, net
    integer                                :: n, ed(2)
    real(kind=WP), dimension(:,:), pointer :: temp
    real(kind=WP), parameter               :: lhf = 334000. !Latent heat of fusion for sea ice (J/Kg)
    real(kind=WP), parameter               :: rhofwt = 1000. !Fresh water density (Kg/m3)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
      temp          => tracers%data(1)%values(:,:)

      hosing_flux=0._WP
      hosing_heat_flux=0._WP

     !!freshwater anomaly is applied on all surface nodes south of 60S (uncomment the next 4 lines)
     !do n=1, myDim_nod2D+eDim_nod2D
     !   y = geo_coord_nod2D(2,n)/rad
     !   if (y<-60._WP) hosing_flux(n)=1._WP
     !end do

     !!freshwater anomaly is applied only on surface coastal nodes south of 60S (uncomment the next 6 lines)
     !!the mask is created looping on all the mesh edges and considering only the nodes belonging to non internal edges 
     do n=1, myDim_edge2D
        ed=mesh%edges(:, n)
        if (myList_edge2D(n) <= mesh%edge2D_in) cycle
        y = sum(geo_coord_nod2D(2,ed))/2._WP/rad
        if (y<-60._WP) hosing_flux(ed)=1._WP
     end do

     !!this is to increase the number of adjacent nodes on which the anomaly is applied, starting from the coastal nodes. 
     !!We use first a smoothing function, choosing the number of nodes (here 8) to smooth the freshwater anomaly and then equalize the value for all nodes
     !call smooth_nod (hosing_flux, 8, partit, mesh) 
     !do n=1,myDim_nod2d+eDim_nod2D
     !   if (hosing_flux(n)>0.0_WP) hosing_flux(n)=1.0_WP
     !end do
     
     !!computing the area-integrated value of hosing_flux (net is the spatial integral) and then normalizing it so that its area integral equals 1
     !!and then multiplying for the actual freshwater flux to impose (hSv)
     call integrate_nod(hosing_flux, net, partit, mesh)

     if (abs(net)>1.e-6) then
        hosing_flux=hosing_flux/net*hSv*1.e6 !the freshwater anomaly is transformed from Sv in m/s (Sv*1.e6 =  m/s)
     end if
     
     water_flux=water_flux-hosing_flux !freshwater flux is negative going into the ocean, so the sign is negative 
     
     !!Computing the heat used to melt the ice corresponding to the freshwater flux
     !!mass_wat = hosing_flux*rhofwt !mass flux of fresh water (Kg/m2*s) 
     !!the heat flux is in W/m2 so no multiplication for the node area is needed 
     !!mass_ice = hosing_flux_ice*rhoice =  hosing_flux*rhofwt !mass of melted ice (Kg/m2*s)
     !! in therms of mass they are equal
     !!hosing_latent_heat_flux = mass_ice*Lf !(J/m2*s=W/m2) 
     
     !!to remove latent heat of fusion in the freshwater experiment uncomment the next line
     !hosing_heat_flux = hosing_flux*rhofwt*lhf !(W/m2)
     
     !!Computing the heat flux change due to the freshwater temperature 
     !!by default freshwater temperature is assumed equal to surrounding temperature so heat flux should be removed if we want the freshwater to have a different temperature
     !!since we are adding super cool freshwater (water at freezing temp) we have a heat flux going out of the ocean, i.e. positive)   
     !!with t_freezing=-1.8, (temp(1,:)-t_freezing) is always positive (temp is almost always >-1.8, so hosing_heat_flux>0)
   
     !!to modify latent heat to take into account different freshwater temperature uncomment the next 2 lines
     !hosing_heat_flux=(temp(1,:)+1.8)*hosing_flux*vcpw  
     !heat_flux=heat_flux+hosing_heat_flux !(W/m2) !heat flux is negative going into the ocean, the hosing heat flux goes out of the ocean, so the sign is positive 

!!print the freshwater amount in the log file (it is printed for every processor so it is a lot or printouts)     
!write(*,*) mype, 'hSv=', hSv, 'net=', net, 'hf=', minval(hosing_flux, 1), maxval(hosing_flux, 1)

end subroutine fw_surf_anomaly

subroutine fw_depth_anomaly(tts, ttt, hSv, partit, mesh) !!!! subroutine to apply freshwater at a chosen depth
    use o_PARAM
    use o_ARRAYS
    USE MOD_PARTIT
    use MOD_MESH
    use MOD_PARSUP
    use g_comm_auto
    use g_support
    use o_tracers
    use g_config,         only: dt
    implicit none

    integer                               :: n, ed(2), row, k, nz, nlev1, nlev2, nzmin, nzmax, elem1, elem2
    real(kind=WP)                         :: x, y, net, spar(100)
    real(kind=WP),  intent(in)            :: hSv
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP),  intent (inout)        :: tts(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP),  intent (inout)        :: ttt(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), save, allocatable      :: mask(:), mask3D(:,:)
    real(kind=WP), parameter              :: lhf = 334000. !Latent heat of fusion for sea ice (J/Kg)
    real(kind=WP), parameter              :: rhofwt = 1000. !Fresh water density (Kg/m3)
    logical, save                         :: lfirst=.true.
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    hosing_flux3D=0._WP
    hosing_heat_flux3D=0._WP
    
    !!generating the 3D mask to apply the anomaly (only for the first iteration)
    if (lfirst) then
       allocate(mask(myDim_nod2D+eDim_nod2D), mask3D(mesh%nl-1, myDim_nod2D+eDim_nod2D))
       mask =0._WP
       mask3D=0._WP
       spar=0._WP
       !!the 2D mask is created in the same way as for the surface application first selecting nodes belonging to coastal edges
       !!and then using the smoothing function to increase it including 5 rows of adjacent nodes
       do n=1, myDim_edge2D
          ed=mesh%edges(:, n)
          if (myList_edge2D(n) <= mesh%edge2D_in) cycle
             y = sum(geo_coord_nod2D(2,ed))/2._WP/rad
          if (y<-60._WP) mask(ed)=1._WP
       end do
       
       call smooth_nod (mask, 5, partit, mesh)
       where (mask>1.e-12)
             mask=1.0_WP
       end where

       !!create a 3D mask to apply the freshwater anomaly between the chosen levels nzmin and nzmax 
       do n=1, myDim_edge2D
          ed = mesh%edges(:, n)
          !!select all the elements adjacent to the edges with both nodes in the 2D mask
          if (any(mask(ed)<1.e-12)) cycle 
          elem1 = mesh%edge_tri(1, n)
          elem2 = mesh%edge_tri(2, n)
          !!determine how many vertical levels exist for each element
          nlev1 = mesh%nlevels(elem1)
          nlev2 = 1
          if (elem2 > 0) nlev2 = mesh%nlevels(elem2)
          !!choose a vertical range of levels where the freshwater anomaly will be applied
          nzmin = 16
          nzmax = 19
          !!select all nodes belonging to edges fully inside the mask2D in the coosen vertical levels
          do k = nzmin, nzmax
             if (k < nlev1 .or. k < nlev2) then !!Apply mask only to levels that exist in at least one of the two elements
             mask3D(k, ed) = 1.0_WP
             end if
          end do
          
       end do
       
       !!collapse the 3D mask back to a 2D mask (a node is kept in mask only if it has at least one active depth level in mask3D)
       !!to compute the area integral and normalize the mask horizontally
       do row=1,myDim_nod2d+eDim_nod2D
          if (sum(mask3D(:, row))<1.e-12) mask(row)=0._WP
       end do
       
       call integrate_nod(mask, net, partit, mesh)
       if (abs(net)>1.e-6) then
          mask=mask/net*1.e6 !hSv will be applied later as it changes in time
       end if
       
       !!Building the 3D distribution
       do row=1,myDim_nod2d+eDim_nod2D 
          spar(:mesh%nl)=0._WP
          !!build a volume-weighted vertical shape inside the chosen levels (to distribute the anomaly proportionally to the volume)
          do k=ulevels_nod2D(row), nlevels_nod2D(row)-1
             spar(k)=areasvol(k,row)*hnode(k,row)*mask3D(k, row)                           
          end do
          !!normalize the vertical distribution so it sums to 1
          spar(:nlevels_nod2D(row)-1)=spar(:nlevels_nod2D(row)-1)/max(sum(spar(:nlevels_nod2D(row)-1)), 1.e-12)
          !!convert the vertically normalized weights so that ts vertical integral at each node matches the horizontal forcing defined by mask(row) 
          do k=ulevels_nod2D(row), nlevels_nod2D(row)-1
             spar(k)=spar(k)/areasvol(k,row)/hnode(k,row)*mask(row)*area(1,row)*dt
          end do
          mask3D(:nlevels_nod2D(row)-1,row)=spar(:nlevels_nod2D(row)-1)
       end do

       lfirst=.false. !set to false so that the mask generation is not repeated
    !!printout the area/volume-integrated flux to check it is correct
    !call integrate_nod(mask3D, net, partit, mesh)
    !write(*,*) 'integrated mask3D=', net
    !call integrate_nod(mask, net, partit, mesh)
    !write(*,*) 'integrated mask2D=', net
    end if

    !!updating the salinity field 
    !!(Snew = Sold (1-dS) where dS=mask3D*hSv is the fraction of the gridcell volume where saltwater is replaced by freshwater during the timestep)
    tts=tts-mask3D*hSv*tts
    !call integrate_nod(tts, net, partit, mesh)
    !write(*,*) 'integrated salinity=', net

    !hosing_heat_flux3D = mask3D*hsv*rhofwt*lhf    
    !ttt = ttt-hosing_heat_flux3D/vcpw !temperature field modification to remove latent heat of fusion
    !ttt=ttt-mask3D*hSv*(ttt+1.8) !temperature field modification to have meltwaterwater added at -1.8 deg.
    !call integrate_nod(ttt, net, partit, mesh)
    !write(*,*) 'integrated temperature=', net

    hosing_flux3D=mask3D*hSv

end subroutine fw_depth_anomaly
