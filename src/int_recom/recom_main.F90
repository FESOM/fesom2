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
  use g_support
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

  real(kind=8)               :: Sali, net, net1, net2
  real(kind=8), allocatable  :: Temp(:), Sali_depth(:), zr(:), PAR(:)
  real(kind=8),  allocatable :: C(:,:)
  real(kind=8),  allocatable :: CO2_watercolumn(:)                                        ! NEW MOCSY
  real(kind=8),  allocatable :: pH_watercolumn(:)                                         ! NEW MOCSY
  real(kind=8),  allocatable :: pCO2_watercolumn(:)                                       ! NEW MOCSY
  real(kind=8),  allocatable :: HCO3_watercolumn(:)                                       ! NEW MOCSY
  real(kind=8),  allocatable :: CO3_watercolumn(:)                                        ! NEW DISS
  real(kind=8),  allocatable :: OmegaC_watercolumn(:)                                     ! NEW DISS
  real(kind=8),  allocatable :: kspc_watercolumn(:)                                       ! NEW DISS
  real(kind=8),  allocatable :: rhoSW_watercolumn(:)                                      ! NEW DISS

  character(len=2)           :: tr_num_name
#include "../associate_mesh.h"

  allocate(Temp(nl-1), Sali_depth(nl-1), zr(nl-1) , PAR(nl-1))
  allocate(CO2_watercolumn(nl-1), pH_watercolumn(nl-1), pCO2_watercolumn(nl-1) , HCO3_watercolumn(nl-1))
  allocate(CO3_watercolumn(nl-1), OmegaC_watercolumn(nl-1), kspc_watercolumn(nl-1) , rhoSW_watercolumn(nl-1))
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
!     if (ulevels_nod2D(n)>1) cycle 
!            nzmin = ulevels_nod2D(n)

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
     Sali_depth(1:nzmax)= tr_arr(1:nzmax, n, 2)                                    ! NEW MOCSY

     !!---- CO2 in the watercolumn                                                 ! NEW MOCSY
     CO2_watercolumn(1:nzmax)    = CO23D(1:nzmax, n)                               ! NEW MOCSY
     pH_watercolumn(1:nzmax)     = pH3D(1:nzmax, n)                                ! NEW MOCSY
     pCO2_watercolumn(1:nzmax)   = pCO23D(1:nzmax, n)                              ! NEW MOCSY
     HCO3_watercolumn(1:nzmax)   = HCO33D(1:nzmax, n)                              ! NEW MOCSY
     CO3_watercolumn(1:nzmax)    = CO33D(1:nzmax, n)                               ! NEW DISS
     OmegaC_watercolumn(1:nzmax) = OmegaC3D(1:nzmax, n)                            ! NEW DISS
     kspc_watercolumn(1:nzmax)   = kspc3D(1:nzmax, n)                              ! NEW DISS
     rhoSW_watercolumn(1:nzmax)  = rhoSW3D(1:nzmax, n)                             ! NEW DISS

     !!---- Biogeochemical tracers
     C(1:nzmax,1:bgc_num) = tr_arr(1:nzmax, n, num_tracers-bgc_num+1:num_tracers)             

     !!---- Depth of the nodes in the water column 
     zr(1:nzmax) = Z_3d_n(1:nzmax, n)                          

     !!---- The PAR in the local water column is initialized
     PAR(1:nzmax) = 0.d0                                        

     !!---- a_ice(row): Ice concentration in the local node
     FeDust = GloFeDust(n) * (1 - a_ice(n)) * dust_sol    
     NDust = GloNDust(n)  * (1 - a_ice(n))

if (Diags) then
!     allocate(Diags3Dloc(nzmax,diags3d_num))
!     Diags3Dloc(:,:) = 0.d0
     !!---- Allocate 3D diagnostics
     allocate(vertgrazmeso_tot(nl-1), vertgrazmeso_n(nl-1), vertgrazmeso_d(nl-1), vertgrazmeso_c(nl-1))
     vertgrazmeso_tot = 0.d0
     vertgrazmeso_n   = 0.d0
     vertgrazmeso_d   = 0.d0
     vertgrazmeso_c   = 0.d0

     allocate(vertrespmeso(nl-1), vertrespmacro(nl-1), vertrespmicro(nl-1))
     vertrespmeso  = 0.d0
     vertrespmacro = 0.d0
     vertrespmicro = 0.d0

     allocate(vertcalcdiss(nl-1), vertcalcif(nl-1))
     vertcalcdiss = 0.d0
     vertcalcif   = 0.d0

     allocate(vertaggn(nl-1), vertaggd(nl-1), vertaggc(nl-1))
     vertaggn = 0.d0
     vertaggd = 0.d0
     vertaggc = 0.d0

     allocate(vertdocexn(nl-1), vertdocexd(nl-1), vertdocexc(nl-1))
     vertdocexn = 0.d0
     vertdocexd = 0.d0
     vertdocexc = 0.d0

     allocate(vertrespn(nl-1), vertrespd(nl-1), vertrespc(nl-1))
     vertrespn = 0.d0
     vertrespd = 0.d0
     vertrespc = 0.d0

     !!---- Allocate 2D diagnostics
     allocate(vertNPPn(nl-1), vertGPPn(nl-1), vertNNAn(nl-1), vertChldegn(nl-1)) 
     vertNPPn = 0.d0
     vertGPPn = 0.d0
     vertNNAn = 0.d0
     Chldegn = 0.d0

     allocate(vertNPPd(nl-1), vertGPPd(nl-1), vertNNAd(nl-1), vertChldegd(nl-1)) 
     vertNPPd = 0.d0
     vertGPPd = 0.d0
     vertNNAd = 0.d0
     Chldegd = 0.d0

     allocate(vertNPPc(nl-1), vertGPPc(nl-1), vertNNAc(nl-1), vertChldegc(nl-1)) 
     vertNPPc = 0.d0
     vertGPPc = 0.d0
     vertNNAc = 0.d0
     Chldegc = 0.d0
end if

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> REcoM_Forcing'//achar(27)//'[0m'

! ======================================================================================
!******************************** RECOM FORCING ****************************************
     call REcoM_Forcing(zr, n, nzmax, C, SW, Loc_slp, Temp, Sali, Sali_depth &
           , CO2_watercolumn                                     & ! NEW MOCSY CO2 for the whole watercolumn
           , pH_watercolumn                                      & ! NEW MOCSY pH for the whole watercolumn
           , pCO2_watercolumn                                    & ! NEW MOCSY pCO2 for the whole watercolumn
           , HCO3_watercolumn                                    & ! NEW MOCSY HCO3 for the whole watercolumn
           , CO3_watercolumn                                     & ! NEW DISS CO3 for the whole watercolumn
           , OmegaC_watercolumn                                  & ! NEW DISS OmegaC for the whole watercolumn
           , kspc_watercolumn                                    & ! NEW DISS stoichiometric solubility product for calcite [mol^2/kg^2]
           , rhoSW_watercolumn                                   & ! NEW DISS in-situ density of seawater [mol/m^3]
           , PAR, mesh)

     tr_arr(1:nzmax, n, num_tracers-bgc_num+1:num_tracers)       = C(1:nzmax, 1:bgc_num)

     !!---- Local variables that have been changed during the time-step are stored so they can be saved
     Benthos(n,1:benthos_num)     = LocBenthos(1:benthos_num)                                ! Updating Benthos values

!     Diags2D(n,1:12)               = LocDiags2D(1:12)                                ! Updating diagnostics

if (Diags) then
     !!---- Updating 2D diagnostics
     NPPn(n) = locNPPn
     NPPd(n) = locNPPd
     GPPn(n) = locGPPn
     GPPd(n) = locGPPd
     NNAn(n) = locNNAn
     NNAd(n) = locNNAd
     Chldegn(n) = locChldegn
     Chldegd(n) = locChldegd
     NPPc(n) = locNPPc
     GPPc(n) = locGPPc
     NNAc(n) = locNNAc
     Chldegc(n) = locChldegc

     !!---- Updating 3D diagnostics
     grazmeso_tot(1:nzmax,n) = vertgrazmeso_tot(1:nzmax)
     grazmeso_n(1:nzmax,n)   = vertgrazmeso_n(1:nzmax)
     grazmeso_d(1:nzmax,n)   = vertgrazmeso_d(1:nzmax)
     grazmeso_c(1:nzmax,n)   = vertgrazmeso_c(1:nzmax)
     respmeso(1:nzmax,n)     = vertrespmeso(1:nzmax)
     respmacro(1:nzmax,n)    = vertrespmacro(1:nzmax)
     respmicro(1:nzmax,n)    = vertrespmicro(1:nzmax)
     calcdiss(1:nzmax,n)     = vertcalcdiss(1:nzmax)
     calcif(1:nzmax,n)       = vertcalcif(1:nzmax)
     aggn(1:nzmax,n)         = vertaggn(1:nzmax)
     aggd(1:nzmax,n)         = vertaggd(1:nzmax)
     aggc(1:nzmax,n)         = vertaggc(1:nzmax)
     docexn(1:nzmax,n)       = vertdocexn(1:nzmax)
     docexd(1:nzmax,n)       = vertdocexd(1:nzmax)
     docexc(1:nzmax,n)       = vertdocexc(1:nzmax)
     respn(1:nzmax,n)        = vertrespn(1:nzmax)
     respd(1:nzmax,n)        = vertrespd(1:nzmax)
     respc(1:nzmax,n)        = vertrespc(1:nzmax)
     NPPn3D(1:nzmax,n)       = vertNPPn(1:nzmax)
     NPPd3D(1:nzmax,n)       = vertNPPd(1:nzmax)
     NPPc3D(1:nzmax,n)       = vertNPPc(1:nzmax)
   
!     do idiags = 1,diags3d_num
!       Diags3D(1:nzmax,n,idiags)  = Diags3Dloc(1:nzmax,idiags) ! 1=NPPnano, 2=NPPdia
!     end do

!     deallocate(Diags3Dloc)

     !!---- Deallocating 2D diagnostics
     deallocate(vertNPPn,vertGPPn,vertNNAn,vertChldegn) 
     deallocate(vertNPPd,vertGPPd,vertNNAd,vertChldegd) 
     deallocate(vertNPPc,vertGPPc,vertNNAc,vertChldegc) 

     !!---- Deallocating 3D Diagnistics
     deallocate(vertgrazmeso_tot, vertgrazmeso_n, vertgrazmeso_d, vertgrazmeso_c)
     deallocate(vertrespmeso, vertrespmacro, vertrespmicro)
     deallocate(vertcalcdiss, vertcalcif)
     deallocate(vertaggn, vertaggd, vertaggc)
     deallocate(vertdocexn, vertdocexd, vertdocexc)
     deallocate(vertrespn, vertrespd, vertrespc)
endif

     GloPCO2surf(n)               = pco2surf(1)
     GlodPCO2surf(n)              = dpco2surf(1)

     GloCO2flux(n)                = dflux(1)
     GloCO2flux_seaicemask(n)     = co2flux_seaicemask(1)                 !  [mmol/m2/s]
     GloO2flux_seaicemask(n)      = o2flux_seaicemask(1)                  !  [mmol/m2/s]
     if (ciso) then
        GloCO2flux_seaicemask_13(n)     = co2flux_seaicemask_13(1)        !  [mmol/m2/s]
        GloCO2flux_seaicemask_14(n)     = co2flux_seaicemask_14(1)        !  [mmol/m2/s]
     end if

     GloHplus(n)                  = ph(1) !hplus
     AtmFeInput(n)                = FeDust
     AtmNInput(n)                 = NDust 
!     DenitBen(n)                  = LocDenit

     GlodecayBenthos(n, 1:benthos_num) = decayBenthos(1:benthos_num)/SecondsPerDay ! convert from [mmol/m2/d] to [mmol/m2/s]  

     PAR3D(1:nzmax,n)             = PAR(1:nzmax) !     PAR3D(inds(1:nn))   = PAR(1:nn)
     CO23D(1:nzmax,n)             = CO2_watercolumn(1:nzmax)       ! NEW MOCSY
     pH3D(1:nzmax,n)              = pH_watercolumn(1:nzmax)        ! NEW MOCSY
     pCO23D(1:nzmax,n)            = pCO2_watercolumn(1:nzmax)      ! NEW MOCSY 
     HCO33D(1:nzmax,n)            = HCO3_watercolumn(1:nzmax)      ! NEW MOCSY
     CO33D(1:nzmax,n)             = CO3_watercolumn(1:nzmax)       ! NEW MOCSY
     OmegaC3D(1:nzmax,n)          = OmegaC_watercolumn(1:nzmax)    ! NEW DISS
     kspc3D(1:nzmax,n)            = kspc_watercolumn(1:nzmax)      ! NEW DISS
     rhoSW3D(1:nzmax,n)           = rhoSW_watercolumn(1:nzmax)     ! NEW DISS

  end do

! ======================================================================================
!************************** EXCHANGE NODAL INFORMATION *********************************			

!if (mype==0) print *, "num_tracers = ", num_tracers
!if (mype==0) print *, "bgc_num = ", bgc_num
!if (mype==0) print *, "num_tracers - bgc_num = ", num_tracers-bgc_num
  do tr_num=num_tracers-bgc_num+1, num_tracers !bgc_num+2 
    call exchange_nod(tr_arr(:,:,tr_num))
  end do

  do n=1, benthos_num
    call exchange_nod(Benthos(:,n))
  end do
  
!  do n=1, 12
!    call exchange_nod(Diags2D(:,n))
!  end do

  if (Diags) then
    call exchange_nod(NPPn)
    call exchange_nod(NPPd)
    call exchange_nod(GPPn)
    call exchange_nod(GPPd)
    call exchange_nod(NNAn)
    call exchange_nod(NNAd)
    call exchange_nod(Chldegn)
    call exchange_nod(Chldegd)
    call exchange_nod(NPPc)
    call exchange_nod(GPPc)
    call exchange_nod(NNAc)
    call exchange_nod(Chldegc)
  endif

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
!  call exchange_nod(DenitBen)	

  call exchange_nod(PAR3D)
  call exchange_nod(CO23D)       ! NEW ms
  call exchange_nod(pH3D)        ! NEW ms
  call exchange_nod(pCO23D)      ! NEW ms
  call exchange_nod(HCO33D)      ! NEW ms
  call exchange_nod(CO33D)       ! NEW ms
  call exchange_nod(OmegaC3D)    ! NEW ms
  call exchange_nod(kspc3D)      ! NEW ms
  call exchange_nod(rhoSW3D)	 ! NEW ms	

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
  if (.not. restore_alkalinity) return
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
