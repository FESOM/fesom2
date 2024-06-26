subroutine recom_sinking_new(tr_num,mesh)

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

  Integer                           :: node, nz,  id, nzmin, nzmax, n,  tr_num, k, nlevels_nod2D_minimum
  Real(kind=8)                      :: wflux(mesh%nl)      ! BALL
  Real(kind=WP)                     :: vd_flux(mesh%nl)
  Real(kind=8)                      :: dz_trr(mesh%nl), aux
  Real(kind=8)                      :: wLoc,wM,wPs
  Real(kind=8)                      :: Rjp,Rj,Rjm

  Real(kind=8)                      :: cfl, d0, d1, thetaP, thetaM, psiP, psiM
  Real(kind=8)                      :: onesixth	= 	1.d0/6.d0
  Real(kind=8)                      :: dt_sink, c1, c2     ! BALL added dt_sink
  Real(kind=8)                      :: Vsink, tv
  Real(kind=8),dimension(mesh%nl)   :: Wvel_flux

#include "../associate_mesh.h"

!< Constant sinking velocities (we prescribe them under namelist recom)
!< This hardcoded part is temporary 
!< .OG. 07.07.2021

    Vsink=0.0_WP

    if (tracer_id(tr_num)==1007 .or.    &  !idetn
        tracer_id(tr_num)==1008 .or.    &  !idetc
        tracer_id(tr_num)==1017 .or.    &  !idetsi
        tracer_id(tr_num)==1021 ) then     !idetcal
	   
            Vsink = VDet 
          
    elseif(tracer_id(tr_num)==1004 .or. &  !iphyn
        tracer_id(tr_num)==1005 .or.    &  !iphyc
        !tracer_id(tr_num)==1020 .or.    &  !iphycal
        tracer_id(tr_num)==1006 ) then     !ipchl

            Vsink = VPhy

    elseif(tracer_id(tr_num)==1013 .or. &  !idian
        tracer_id(tr_num)==1014 .or.    &  !idiac
        tracer_id(tr_num)==1016 .or.    &  !idiasi
        tracer_id(tr_num)==1015 ) then     !idchl

            Vsink = VDia

#if defined (__coccos)
    elseif(tracer_id(tr_num)==1029 .or. &  !icocn
        tracer_id(tr_num)==1030 .or.    &  !icocc
        tracer_id(tr_num)==1031 ) then     !icchl

            Vsink = VCocco
#endif
    elseif(tracer_id(tr_num)==1020) then   !iphycal
       
#if defined (__coccos)
               Vsink = VCocco
#else
               Vsink = VPhy
#endif
            
#if defined (__3Zoo2Det)
    elseif(tracer_id(tr_num)==1025 .or. &  !idetz2n
         tracer_id(tr_num)==1026 .or. &  !idetz2c
         tracer_id(tr_num)==1027 .or. &  !idetz2si
         tracer_id(tr_num)==1028 ) then  !idetz2calc 
            
            Vsink = VDet_zoo2            
#endif
    end if

!! No sinking if background sinking velocity is less than 0.1 m/day
if (Vsink .gt. 0.1) then

   do n = 1,myDim_nod2D
      if (ulevels_nod2D(n)>1) cycle 
      nzmin = ulevels_nod2D(n)
      nzmax = nlevels_nod2D(n)-1

      ! distance between tracer points, surface and bottom dz_trr is half 
      ! the layer thickness
      dz_trr                = 0.0d0
      dz_trr(nzmin+1:nzmax) = abs(Z_3d_n(nzmin:nzmax-1,n)-Z_3d_n(nzmin+1:nzmax,n))
      dz_trr(nzmin)         = hnode(nzmin,n)/2.0d0
      dz_trr(nzmax+1)       = hnode(nzmax,n)/2.0d0

      ! Vertical sinking velocity for BCG tracers
      Wvel_flux(nzmin:nzmax+1)= 0.d0  
                                                                     
      do nz=nzmin,nzmax+1

         Wvel_flux(nz) = -Vsink/SecondsPerDay ! allow_var_sinking = .false.

         if (allow_var_sinking) then 
            if (use_ballasting) then 

! Apply ballasting on slow sinking detritus
!if (any(recom_sinking_tracer_id == tracer_id(tr_num))) then 

                if (tracer_id(tr_num)==1007 .or.    &  !idetn
                    tracer_id(tr_num)==1008 .or.    &  !idetc
                    tracer_id(tr_num)==1017 .or.    &  !idetsi
                    tracer_id(tr_num)==1021 ) then     !idetcal  

                     Wvel_flux(nz) = w_ref1 * scaling_density1_3D(nz,n) * scaling_visc_3D(nz,n)

                     if (depth_scaling1.gt.0.0) Wvel_flux(nz) = Wvel_flux(nz) + (depth_scaling1 * abs(zbar_3d_n(nz,n)))

                     if (abs(Wvel_flux(nz)) .gt. max_sinking_velocity) Wvel_flux(nz) = max_sinking_velocity

                     ! sinking velocity [m d-1] surface --> bottom (negative)
                     Wvel_flux(nz) = -1.0d0 * Wvel_flux(nz)/SecondsPerDay ! now in [m s-1]
                 end if
            else ! use_ballasting = .false.
               Wvel_flux(nz) = -((Vdet_a * abs(zbar_3d_n(nz,n))/SecondsPerDay) + Vsink/SecondsPerDay)
            endif
         end if

#if defined (__3Zoo2Det)

! Apply ballasting on fast sinking detritus

            ! We assume constant sinking for second detritus
            if(tracer_id(tr_num)==1025 .or. &  !idetz2n
               tracer_id(tr_num)==1026 .or. &  !idetz2c
               tracer_id(tr_num)==1027 .or. &  !idetz2si
               tracer_id(tr_num)==1028 ) then  !idetz2calc  
 
               if (use_ballasting) then    ! NEW BALL

                  Wvel_flux(nz) = w_ref2*scaling_density2_3D(nz,n)*scaling_visc_3D(nz,n)

                  if (depth_scaling2.gt.0.0) Wvel_flux(nz) = Wvel_flux(nz) + (depth_scaling2 * abs(zbar_3d_n(nz,n)))

                  if (abs(Wvel_flux(nz)) .gt. max_sinking_velocity) Wvel_flux(nz) = max_sinking_velocity

                ! sinking velocity [m d-1] surface --> bottom (negative)
                  Wvel_flux(nz) = -1.0d0 * Wvel_flux(nz)/SecondsPerDay ! now in [m s-1]

               else ! use_ballasting = .false.
                  Wvel_flux(nz) = -VDet_zoo2/SecondsPerDay ! --> VDet_zoo2 ! NEW BALL changed -Vsink to -VDet_zoo2
               end if
            endif ! second detritus tracers
#endif

         if (tracer_id(tr_num)==1021) Sinkvel1_tr(nz,n,tr_num) = Wvel_flux(nz) !-1.0d0/SecondsPerDay  !idetcal  
#if defined (__3Zoo2Det) 
         if (tracer_id(tr_num)==1028) Sinkvel2_tr(nz,n,tr_num) = Wvel_flux(nz)  !idetz2calc
#endif
      end do

      wflux   = 0.d0
      dt_sink = dt
      vd_flux = 0.0d0

if (1) then ! 3rd Order DST Sceheme with flux limiting. This code comes from old recom

      k=nod_in_elem2D_num(n)
      ! Screening minimum depth in neigbouring nodes around node n
      nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

      vd_flux(nzmin:nzmax+1)= 0.0_WP

      do nz=nzmax, nzmin+1,-1

         Rjp = tr_arr(nz,n,tr_num)              - tr_arr(min(nz+1,nzmax),n,tr_num)
         Rj  = tr_arr(max(nzmin,nz-1),n,tr_num) - tr_arr(nz,n,tr_num) 
         Rjm = tr_arr(max(nzmin,nz-2),n,tr_num) - tr_arr(max(nzmin,nz-1),n,tr_num)

         cfl = abs(Wvel_flux(nz) * dt_sink / dz_trr(nz)) !(Z_n(nz-1)-Z_n(nz)))       ! [m/day] * [day] * [1/m]  ! NEW BALL changed dt to dt_sink

         wPs = Wvel_flux(nz) + abs(Wvel_flux(nz)) ! --> Positive vertical velocity
         wM  = Wvel_flux(nz) - abs(Wvel_flux(nz)) ! --> Negative vertical velocity

         d0 = (2.d0 - cfl)*(1.d0 - cfl)*onesixth
         d1 = (1.d0 - cfl*cfl)*onesixth
	
         thetaP = Rjm/(1.d-20+Rj)
         psiP = d0 + d1*thetaP
         psiP = max(0.d0, min(min(1.d0,psiP), &
            (1.d0-cfl)/(1.d-20+cfl)*thetaP))

         thetaM = Rjp/(1.d-20 + Rj)	
         psiM = d0 + d1*thetaM
         psiM = max(0.d0, min(min(1.d0,psiM), &
            (1.d0-cfl)/(1.d-20-cfl)*thetaM))

         tv= (0.5 * wPs * (tr_arr(nz,n,tr_num)              + psiM * Rj)+ &
	      0.5 * wM  * (tr_arr(max(nzmin,nz-1),n,tr_num) + psiP * Rj))
         vd_flux(nz)= - tv*area(nz,n)
      end do
end if ! 3rd Order DST Sceheme with flux limiting

if (0) then ! simple upwind

      ! Surface flux
      vd_flux(nzmin)= 0.0_WP

      ! Bottom flux
      vd_flux(nzmax+1)= 0.0_WP

      k=nod_in_elem2D_num(n)
      ! Screening minimum depth in neigbouring nodes around node n
      nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

      do nz=nzmin+1,nzmax !nlevels_nod2D_minimum-1
!         tv = tr_arr(nz,n,tr_num)                                ! simple scheme        - test1
!         tv = 0.5_WP*(tr_arr(nz-1,n,tr_num)+tr_arr(nz,n,tr_num)) ! consider both layers - test2  
!         tv = tv*Wvel_flux(nz) ! Wvel_flux is negative
         tv = - 0.5* & ! - test3
            (tr_arr(nz-1,n,tr_num)*(Wvel_flux(nz)-abs(Wvel_flux(nz))) + &
             tr_arr(nz  ,n,tr_num)*(Wvel_flux(nz)+abs(Wvel_flux(nz))))
         vd_flux(nz)= tv*area(nz,n)

      end do
end if ! simple upwind
      do nz=nzmin,nzmax
         vert_sink(nz,n) = vert_sink(nz,n) + (vd_flux(nz)-vd_flux(nz+1))*dt/areasvol(nz,n)/hnode_new(nz,n) !/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
      end do
   end do
end if ! Vsink .gt. 0.1

end subroutine recom_sinking_new
!-------------------------------------------------------------------------------  
! Subroutine calculate ballasting
!-------------------------------------------------------------------------------  
!subroutine ballast(tr_num,mesh)
subroutine ballast(mesh)
  use recom_config
  use recom_GloVar
  USE mod_MESH
  USE o_PARAM
  USE o_ARRAYS
  USE g_PARSUP
  USE g_CONFIG
  use g_forcing_arrays
  use g_comm_auto
  use i_param
  use i_arrays
  use i_therm_param
  use g_clock
  use g_rotate_grid
  use g_comm
  use mvars                                                                ! NEW MOCSY
  use mdepth2press                                                         ! NEW BALL
  use gsw_mod_toolbox, only: gsw_sa_from_sp,gsw_ct_from_pt,gsw_rho         ! NEW BALL

  implicit none

  integer                                                 :: row, k, nzmin, nzmax!, tr_num
  Real(kind=8)                                            :: depth_pos(1)         ! NEW BALL pos. depth -> for ballasting
  Real(kind=8)                                            :: pres(1)              ! NEW BALL local pressure
  Real(kind=8)                                            :: sa(1)                ! NEW BALL local absolute salinity
  Real(kind=8)                                            :: ct(1)                ! NEW BALL local conservative temperature
  Real(kind=8)                                            :: rho_seawater(1)      ! NEW BALL local seawater density
  Real(kind=8)                                            :: Lon_degree(1)        ! NEW BALL local seawater density
  Real(kind=8)                                            :: Lat_degree(1)        ! NEW BALL local seawater density
  type(t_mesh), intent(in), target :: mesh

#include  "../associate_mesh.h"

  ! For ballasting, calculate scaling factors here and pass them to FESOM, where sinking velocities are calculated
     ! -----
     ! If ballasting is used, sinking velocities are a function of a) particle composition (=density),
     ! b) sea water viscosity, c) depth (currently for small detritus only), and d) a constant reference sinking speed
     ! -----

! check oce_ale_tracer.F90
!     call get_seawater_viscosity(mesh) ! seawater_visc_3D
!     call get_particle_density(mesh) ! rho_particle = density of particle class 1 and 2 

     !___________________________________________________________________________
     ! loop over local nodes
     do row=1,myDim_nod2D 
         ! max. number of levels at node n
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)
         !! lon
        Lon_degree(1)=geo_coord_nod2D(1,row)/rad !! convert from rad to degree
         !! lat
        Lat_degree(1)=geo_coord_nod2D(2,row)/rad !! convert from rad to degree

        ! get scaling vectors -> these need to be passed to FESOM to get sinking velocities
        ! get local seawater density
        do k=nzmin, nzmax

           !! level depth 
           depth_pos(1) = abs(Z_3d_n(k,row))  ! take depth of tracers instead of levels !abs(zbar_3d_n(k,row))

           call depth2press(depth_pos(1), Lat_degree(1), pres, 1)  ! pres is output of function,1=number of records
           sa           = gsw_sa_from_sp(tr_arr(k,row,2), pres, Lon_degree(1), Lat_degree(1))
           ct           = gsw_ct_from_pt(sa, tr_arr(k,row,1))
           rho_seawater = gsw_rho(sa, ct, pres)

           ! (i.e. no density scaling)
           scaling_density1_3D(k,row)=1.0
           scaling_density2_3D(k,row)=1.0

              if (use_density_scaling) then
                    !if (tr_arr(k,row,10)>0.001) then ! idetc only apply ballasting above a certain biomass
                       scaling_density1_3D(k,row) = (rho_particle1(k,row)-rho_seawater(1))/(rho_ref_part-rho_ref_water)
                    !endif 

#if defined (__3Zoo2Det)
                    !if (tr_arr(k,row,28)>0.001) then ! idetz2c only apply ballasting above a certain biomass
                       scaling_density2_3D(k,row) = (rho_particle2(k,row)-rho_seawater(1))/(rho_ref_part-rho_ref_water)
                    !endif 
#endif
              endif

            scaling_visc_3D(k,row)=1.0

            if (use_viscosity_scaling) then
                if (seawater_visc_3D(k,row)==0) then 
                    scaling_visc_3D(k,row)=1.0
                else
                    scaling_visc_3D(k,row)= visc_ref_water/seawater_visc_3D(k,row)
                endif
            endif
        end do
        scaling_visc_3D(nzmax+1,row) = scaling_visc_3D(nzmax,row)
     end do

    ! in the unlikely (if possible at all...) case that rho_particle(k)-rho_seawater(1)<0, prevent the scaling factor from being negative

    if (any(scaling_density1_3D(:,:) <= tiny)) scaling_density1_3D(:,:) = 1.0_WP      ! tiny = 2.23D-16
#if defined (__3Zoo2Det)
    if (any(scaling_density2_3D(:,:) <= tiny)) scaling_density2_3D(:,:) = 1.0_WP      ! tiny = 2.23D-16    
#endif
end subroutine ballast
!-------------------------------------------------------------------------------  
! Subroutine calculate density of particle 
! depending on composition (detC, detOpal, detCaCO3) based on Cram et al. (2018)
!-------------------------------------------------------------------------------  
subroutine get_particle_density(mesh)          ! NEW BALL developed by Cara and Onur  
  use recom_config
  use recom_GloVar
  USE mod_MESH
  USE o_PARAM
  USE o_ARRAYS
  USE g_PARSUP
  USE g_CONFIG
  use g_forcing_arrays
  use g_comm_auto
  use i_param
  use i_arrays
  use i_therm_param
  use g_clock
  use g_rotate_grid
  use g_comm

  implicit none
 
  integer                                                 :: row, k, nzmin, nzmax, tr_num
  type(t_mesh), intent(in), target                        :: mesh

  real(kind=8)                                            :: a1(mesh%nl-1, myDim_nod2D+eDim_nod2D) ! [n.d.] fraction of carbon in detritus class
  real(kind=8)                                            :: a2(mesh%nl-1, myDim_nod2D+eDim_nod2D) ! [n.d.] fraction of nitrogen in detritus class
  real(kind=8)                                            :: a3(mesh%nl-1, myDim_nod2D+eDim_nod2D) ! [n.d.] fraction of Opal in detritus class
  real(kind=8)                                            :: a4(mesh%nl-1, myDim_nod2D+eDim_nod2D) ! [n.d.] fraction of CaCO3 in detritus class
  real(kind=8)                                            :: b1(mesh%nl-1, myDim_nod2D+eDim_nod2D)
  real(kind=8)                                            :: b2(mesh%nl-1, myDim_nod2D+eDim_nod2D)
  real(kind=8)                                            :: b3(mesh%nl-1, myDim_nod2D+eDim_nod2D)
  real(kind=8)                                            :: b4(mesh%nl-1, myDim_nod2D+eDim_nod2D)
  real(kind=8)                                            :: aux(mesh%nl-1, myDim_nod2D+eDim_nod2D)

#include "../associate_mesh.h"

  rho_particle1 = 0.0
  b1 = 0.0
  b2 = 0.0
  b3 = 0.0
  b4 = 0.0
  aux = 0.0

! OG below guarantees non-negative tracer field
  do tr_num=1,num_tracers
     if (tracer_id(tr_num)==1008)  b1 = max(tiny,tr_arr(:,:,tr_num)) !idetc      ! [mmol m-3] detritus carbon
     if (tracer_id(tr_num)==1007)  b2 = max(tiny,tr_arr(:,:,tr_num)) !idetn      ! [mmol m-3] detritus nitrogen
     if (tracer_id(tr_num)==1017)  b3 = max(tiny,tr_arr(:,:,tr_num)) !idetsi     ! [mmol m-3] detritus Si
     if (tracer_id(tr_num)==1021)  b4 = max(tiny,tr_arr(:,:,tr_num)) !idetcal    ! [mmol m-3] detritus CaCO3
  end do

  do row=1,myDim_nod2d
     nzmin = ulevels_nod2D(row)
     nzmax = nlevels_nod2D(row)
     aux(nzmin:nzmax,row) = b1(nzmin:nzmax,row)+b2(nzmin:nzmax,row)+b3(nzmin:nzmax,row)+b4(nzmin:nzmax,row)
     a1(nzmin:nzmax,row)  = b1(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
     a2(nzmin:nzmax,row)  = b2(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
     a3(nzmin:nzmax,row)  = b3(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
     a4(nzmin:nzmax,row)  = b4(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
     rho_particle1(nzmin:nzmax,row) = rho_CaCO3*a4(nzmin:nzmax,row) + rho_opal*a3(nzmin:nzmax,row) + rho_POC*a1(nzmin:nzmax,row) + rho_PON*a2(nzmin:nzmax,row)
     rho_particle1(nzmax+1,row) = rho_particle1(nzmax,row)
  end do

#if defined (__3Zoo2Det)
     rho_particle2 = 0.0
     b1 = 0.0
     b2 = 0.0
     b3 = 0.0
     b4 = 0.0
     aux = 0.0
     do tr_num=1,num_tracers
        if (tracer_id(tr_num)==1026)  b1 = max(tiny,tr_arr(:,:,tr_num)) !idetz2c
        if (tracer_id(tr_num)==1025)  b2 = max(tiny,tr_arr(:,:,tr_num)) !idetz2n
        if (tracer_id(tr_num)==1027)  b3 = max(tiny,tr_arr(:,:,tr_num)) !idetz2si
        if (tracer_id(tr_num)==1028)  b4 = max(tiny,tr_arr(:,:,tr_num)) !idetz2calc 
     end do

     do row=1,myDim_nod2d+eDim_nod2D   ! myDim is sufficient
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)
        aux(nzmin:nzmax,row) = b1(nzmin:nzmax,row)+b2(nzmin:nzmax,row)+b3(nzmin:nzmax,row)+b4(nzmin:nzmax,row)
        a1(nzmin:nzmax,row)  = b1(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a2(nzmin:nzmax,row)  = b2(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a3(nzmin:nzmax,row)  = b3(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a4(nzmin:nzmax,row)  = b4(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        rho_particle2(nzmin:nzmax,row) = rho_CaCO3*a4(nzmin:nzmax,row) + rho_opal*a3(nzmin:nzmax,row) + rho_POC*a1(nzmin:nzmax,row) + rho_PON*a2(nzmin:nzmax,row)
        rho_particle2(nzmax+1,row) = rho_particle2(nzmax,row)
     end do
#endif

end subroutine get_particle_density
!-------------------------------------------------------------------------------   
! Subroutine to approximate seawater viscosity with current temperature 
! based on Cram et al. (2018)
!------------------------------------------------------------------------------- 

! neglecting salinity effects, which are much smaller than those of temperature
! https://bitbucket.org/ohnoplus/ballasted-sinking/src/master/tools/waterviscosity.m

subroutine get_seawater_viscosity(mesh)!,Nn,Temp,Salt,seawater_visc,mesh) ! NEW BALL developed by Cara and Onur    
  
  use recom_config
  use recom_GloVar
  USE mod_MESH
  USE o_PARAM
  USE o_ARRAYS
  USE g_PARSUP
  USE g_CONFIG
  use g_forcing_arrays
  use g_comm_auto
  use i_param
  use i_arrays
  use i_therm_param
  use g_clock
  use g_rotate_grid
  use g_comm

  implicit none

!  integer, intent(in)                                     :: Nn   !< Total number of nodes in the vertical
!  real(kind=8),dimension(mesh%nl-1)        ,intent(in)    :: Temp !< [degrees C] Ocean temperature
!  real(kind=8),dimension(mesh%nl-1)        ,intent(in)    :: Salt !< [g/kg or n.d.] Ocean salinity

!  real(kind=8),dimension(mesh%nl-1)                       :: seawater_visc !<[kg m-1 s-1] Ocean viscosity
  real(kind=8),dimension(1)                               :: A,B,mu_w
  integer                                                 :: row, k, nzmin, nzmax

  type(t_mesh), intent(in) , target :: mesh

#include "../associate_mesh.h"

  seawater_visc_3D(:,:) = 0.0
  do row=1,myDim_nod2d
     !if (ulevels_nod2D(row)>1) cycle
! Do we need a kind of constriction here?
! i.e., if (seawater_visc_3D(row)<=0.0_WP) cycle
     nzmin = ulevels_nod2D(row)
     nzmax = nlevels_nod2D(row)

     do k=nzmin, nzmax
     ! Eq from Sharaway 2010
     ! validity: 
     !  0<temp<180°C
     !  0<salt<0.15 kg/kg
     ! Note: because salinity is expected to be in kg/kg, use conversion factor 0.001 below!
        A(1)             = 1.541 + 1.998*0.01*tr_arr(k,row,1) - 9.52*1e-5*tr_arr(k,row,1)*tr_arr(k,row,1)
        B(1)             = 7.974 - 7.561*0.01*tr_arr(k,row,1) + 4.724*1e-4*tr_arr(k,row,1)*tr_arr(k,row,1)
        mu_w(1)          = 4.2844*1.0e-5 + (1.0/(0.157*(tr_arr(k,row,1)+64.993)*(tr_arr(k,row,1)+64.993)-91.296))
        seawater_visc_3D(k,row) = mu_w(1) * (1.0 + A(1)*tr_arr(k,row,2)*0.001 + B(1)*tr_arr(k,row,2)*0.001*tr_arr(k,row,2)*0.001)
     enddo
  end do

end subroutine get_seawater_viscosity
