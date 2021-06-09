! CONTENT:
! ------------
!    subroutine recom_init
!    subroutine recom_forcing
!    subroutine recom_sms
!
! initially written by REcoM group, adapted by  O. Gurses 02.03.2020

subroutine recom(mesh)

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
  implicit none
  type(t_mesh), intent(in) , target :: mesh
! ======================================================================================
!! Depth information

!! zbar (depth of layers) and Z (mid depths of layers)
!! zbar is negative 
!! zbar(nl) allocate the array for storing the standard depths
!! Z(nl-1)  mid-depths of cells

!! max. number of levels at node n
!! nzmax = nlevels_nod2D(n)
!! u_ice and v_ice are at nodes
!! u_w, v_w are at nodes (interpolated from elements)
!! u_wind and v_wind are always at nodes
! ======================================================================================

  real(kind=8)               :: SW, Loc_slp
  integer                    :: tr_num
  integer                    :: nz, n, nzmin, nzmax  
  integer                    :: idiags

  real(kind=8)               :: Sali
  real (kind=8), allocatable :: Temp(:),  zr(:), PAR(:)
  real(kind=8),  allocatable :: C(:,:)
  character(len=2)           :: tr_num_name
#include "../associate_mesh.h"

  allocate(Temp(nl-1), zr(nl-1) , PAR(nl-1))
  allocate(C(nl-1,bgc_num))

  if (.not. use_REcoM) return

! ======================================================================================
!************************* READ SURFACE BOUNDARY FILES *********************************			

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> Atm_input'//achar(27)//'[0m'

  call Atm_input(mesh)        !<  read surface atmospheric deposition for Fe, N, CO2
  call River_input(mesh)      !<  read riverine input
  call Erosion_input(mesh)    !<  read erosion input

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> bio_fluxes'//achar(27)//'[0m'

  call bio_fluxes(mesh)       !<  alkalinity restoring/ virtual flux is possible

! ======================================================================================
!********************************* LOOP STARTS *****************************************			

  do n=1, myDim_nod2D  ! needs exchange_nod in the end
     if (ulevels_nod2D(n)>1) cycle 
            nzmin = ulevels_nod2D(n)
!            nzmax = nlevels_nod2D(n)-1

     !!---- Number of vertical layers
     nzmax = nlevels_nod2D(n)-1

     !!---- This is needed for piston velocity 
     Loc_ice_conc = a_ice(n) 

     !!---- Mean sea level pressure 
     Loc_slp = press_air(n)

     !!---- Benthic layers
     LocBenthos(1:benthos_num) = Benthos(n,1:benthos_num)

     !!---- Local conc of [H+]-ions from last time time step. Stored in LocVar
     !!---- used as first guess for H+ conc.in subroutine CO2flux (provided by recom_init)
     Hplus = GloHplus(n)                                  

     !!---- Interpolated wind from atmospheric forcing 
     !!---- temporarily stored in module LocVar
     ULoc = sqrt(u_wind(n)**2+v_wind(n)**2)

     !!---- Atmospheric CO2 in LocVar                                                                        
     LocAtmCO2         = AtmCO2(month)   
     if (ciso) then
        LocAtmCO2_13 = AtmCO2_13(month)
        LocAtmCO2_14 = AtmCO2_14(month)
        r_atm_13 = LocAtmCO2_13(1) / LocAtmCO2(1)
        r_atm_14 = LocAtmCO2_14(1) / LocAtmCO2(1)
     end if

     !!---- Shortwave penetration
     SW = parFrac * shortwave(n)
     SW = SW * (1.d0 - a_ice(n))

     !!---- Temperature in water column
     Temp(1:nzmax) = tr_arr(1:nzmax, n, 1)

     !!---- Surface salinity
     Sali = tr_arr(1,       n, 2)

     !!---- Biogeochemical tracers
     C(1:nzmax,1:bgc_num) = tr_arr(1:nzmax, n, 3:num_tracers)             

     !!---- Depth of the nodes in the water column 
     zr(1:nzmax) = Z_3d_n(1:nzmax, n)                          

     !!---- The PAR in the local water column is initialized
     PAR(1:nzmax) = 0.d0                                        

     !!---- a_ice(row): Ice concentration in the local node
     FeDust = GloFeDust(n) * (1 - a_ice(n)) * dust_sol    
     NDust = GloNDust(n)  * (1 - a_ice(n))

     allocate(Diags3Dloc(nzmax,8))
     Diags3Dloc(:,:) = 0.d0

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> REcoM_Forcing'//achar(27)//'[0m'


! ======================================================================================
!******************************** RECOM FORCING ****************************************
			
     call REcoM_Forcing(zr, n, nzmax, C, SW, Loc_slp, Temp, Sali, PAR, mesh)

     tr_arr(1:nzmax, n, 3:num_tracers)       = C(1:nzmax, 1:bgc_num)
   
     !!---- Local variables that have been changed during the time-step are stored so they can be saved
     Benthos(n,1:benthos_num)     = LocBenthos(1:benthos_num)                                ! Updating Benthos values

     Diags2D(n,1:8)               = LocDiags2D(1:8)                                ! Updating diagnostics
     GloPCO2surf(n)               = pco2surf(1)
     GlodPCO2surf(n)              = dpco2surf(1)

     GloCO2flux(n)                = dflux(1)
     GloCO2flux_seaicemask(n)     = co2flux_seaicemask(1)                 !  [mmol/m2/s]
     GloO2flux_seaicemask(n)      = o2flux_seaicemask(1)                  !  [mmol/m2/s]
     if (ciso) then
!        tr_arr(1:nzmax, n, 25:40)    = C(1:nzmax,23:38)
        GloCO2flux_seaicemask_13(n)     = co2flux_seaicemask_13(1)        !  [mmol/m2/s]
        GloCO2flux_seaicemask_14(n)     = co2flux_seaicemask_14(1)        !  [mmol/m2/s]
     end if

     GloHplus(n)                  = hplus
     AtmFeInput(n)                = FeDust
     AtmNInput(n)                 = NDust 
     DenitBen(n)                  = LocDenit

!if (mype==0) write (*,*) "decayBenthos_Si mmol/m2/day= " , decayBenthos(3)
     GlodecayBenthos(n, 1:benthos_num) = decayBenthos(1:benthos_num)/SecondsPerDay  ! convert from [mmol/m2/d] to [mmol/m2/s]  
!     GloWflux(n, 1) = wFluxDia(3)/SecondsPerDay  ! convert from [mmol/m2/d] to [mmol/m2/s]  
!     GloWflux(n, 2) = wFluxDet(3)/SecondsPerDay  ! convert from [mmol/m2/d] to [mmol/m2/s]
!if (mype==0) write (*,*) "decayBenthos_Si mmol/m2/s=   " , GlodecayBenthos(n, 3)

!    GlowFluxDet(n,1:benthos_num)=wFluxDet(1:benthos_num)/SecondsPerDay  ! convert from [mmol/m2/d] to [mmol/m2/s]
!if   (NitrogenSS) LocDenit !!!!!
!     GlowFluxPhy(n,1:benthos_num)=wFluxPhy(1:benthos_num)/SecondsPerDay  ! convert from [mmol/m2/d] to [mmol/m2/s]
!     GlowFluxDia(n,1:benthos_num)=wFluxDia(1:benthos_num)/SecondsPerDay  ! convert from [mmol/m2/d] to [mmol/m2/s]

     PAR3D(1:nzmax,n)             = PAR(1:nzmax) !     PAR3D(inds(1:nn))   = PAR(1:nn)
   
     do idiags = 1,diags3d_num
       Diags3D(1:nzmax,n,idiags)  = Diags3Dloc(1:nzmax,idiags) ! 1=NPPnano, 2=NPPdia
     end do

     deallocate(Diags3Dloc)

  end do

! ======================================================================================
!************************** EXCHANGE NODAL INFORMATION *********************************			

  do tr_num=3, bgc_num+2 
    call exchange_nod(tr_arr(:,:,tr_num))
  end do

  do n=1, benthos_num
    call exchange_nod(Benthos(:,n))
  end do
  
  do n=1, 8
    call exchange_nod(Diags2D(:,n))	
  end do

  call exchange_nod(GloPCO2surf)	
  call exchange_nod(GloCO2flux)	
  call exchange_nod(GloCO2flux_seaicemask)

  do n=1, 4
    call exchange_nod(GlodecayBenthos(:,n))
  end do
  if (ciso) then
    call exchange_nod(GloPCO2surf_13)
    call exchange_nod(GloPCO2surf_14)
    call exchange_nod(GloCO2flux_13)
    call exchange_nod(GloCO2flux_14)
    call exchange_nod(GloCO2flux_seaicemask_13)
    call exchange_nod(GloCO2flux_seaicemask_14)  
  end if
  call exchange_nod(GloO2flux_seaicemask)	
  call exchange_nod(GloHplus)	
  call exchange_nod(AtmFeInput)	
  call exchange_nod(AtmNInput)	
  call exchange_nod(DenitBen)	
  call exchange_nod(PAR3D)	
!  do n=1, 2
!     call exchange_nod(Diags3D(:,:,n))	
!  end do
end subroutine recom
! ======================================================================================
! Alkalinity restoring to climatology                                 			*
! =========================================================================== bio_fluxes 
subroutine bio_fluxes(mesh)

  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use recom_config

  use mod_MESH
  USE g_CONFIG
  use o_ARRAYS
  use i_ARRAYS
  use g_comm_auto
  use g_forcing_param, only: use_virt_salt
  use g_forcing_arrays
  use g_PARSUP
  use g_support
  use i_therm_param

  implicit none
  integer                           :: n, elem, elnodes(3),n1
  real(kind=WP)                     :: ralk, net
  type(t_mesh), intent(in) , target :: mesh
#include "../associate_mesh.h"
  
!___________________________________________________________________
! on freshwater inflow/outflow or virtual alkalinity:
  ! 1. In zlevel & zstar the freshwater flux is applied in the update of the 
  ! ssh matrix when solving the continuity equation of vertically 
  ! integrated flow. The alcalinity concentration in the first layer will 
  ! be then adjusted according to the change in volume.

  ! In this case ralk is forced to be zero by setting ref_alk=0. and ref_alk_local=.false.

  ! 2. In cases where the volume of the upper layer is fixed (i.e. linfs)  the freshwater flux 
  ! 'ralk*water_flux(n)' is applied as a virtual alkalinity boundary condition via the vertical 
  ! diffusion operator.

  ! --> ralk*water_flux(n) : virtual alkalinity flux 
  ! virtual alkalinity flux

!  if (use_virt_alk) then ! OG in case of virtual alkalinity flux
!     ralk=ref_alk
!     do n=1, myDim_nod2D+eDim_nod2D
!        if (ref_alk_local) ralk = tr_arr(1,n,5)
!        virtual_alk(n)=ralk*water_flux(n) 
!     end do
!  end if

  ! Alkalinity restoring to climatology

  do n=1, myDim_nod2D+eDim_nod2D
     relax_alk(n)=surf_relax_Alk*(Alk_surf(n)-tr_arr(1,n,2+ialk)) ! 1 temp, 2 salt
  end do

!  if (mype==0) then
!      if (mype==0)  write(*,*) 'Alk_surf  = ', Alk_surf
!      if (mype==0)  write(*,*) 'tr_arr  = ', tr_arr(1,:,5)
!      if (mype==0)  write(*,*) 'relax_alk  = ', relax_alk
!  endif


  ! 2. virtual alkalinity flux
!  if (use_virt_alk) then ! is already zero otherwise
!     call integrate_nod(virtual_alk, net)
!     virtual_alk=virtual_alk-net/ocean_area
!  end if


  ! 3. restoring to Alkalinity climatology
  call integrate_nod(relax_alk, net, mesh)

!  if (mype==0) then
!      if (mype==0)  write(*,*) 'ocean_area  = ', ocean_area
!      if (mype==0)  write(*,*) 'net  = ', net
!      if (mype==0)  write(*,*) 'net/ocean_area  = ', net/ocean_area
!  endif

  relax_alk=relax_alk-net/ocean_area  ! at ocean surface layer


!  if (mype==0) then
!     write(*,*) '____________________________________________________________'
!     write(*,*) ' --> relax_alk,  = ', relax_alk
!  endif

end subroutine bio_fluxes
