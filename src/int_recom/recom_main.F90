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

  logical :: do_update = .false. 

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

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'   --> River_input'//achar(27)//'[0m'
  call River_input(mesh)      !<  read riverine input

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//' --> Erosion_input'//achar(27)//'[0m'
  call Erosion_input(mesh)    !<  read erosion input

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> bio_fluxes'//achar(27)//'[0m'
  call bio_fluxes(mesh)       !<  alkalinity restoring/ virtual flux is possible

  if (use_atbox) then    ! MERGE
! Prognostic atmospheric isoCO2
    call recom_atbox(mesh)
!   optional I/O of isoCO2 and inferred cosmogenic 14C production; this may cost some CPU time
    if (ciso .and. ciso_14) then
      call annual_event(do_update)
      if (do_update .and. mype==0) write (*, fmt = '(a50,2x,i6,4(2x,f6.2))') &
                                         'Year, xCO2 (ppm), cosmic 14C flux (at / cm² / s):', &
                                          yearold, x_co2atm(1), x_co2atm_13(1), x_co2atm_14(1), cosmic_14(1) * production_rate_to_flux_14
    end if
  end if

! MERGE
! ======================================================================================
!************************* READ BOTTOM BOUNDARY FILES*********************************
if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> Sed_input'//achar(27)//'[0m'
if (use_MEDUSA .and. (sedflx_num .ne. 0)) then
   call Sed_input(mesh) !! --> sedimentary input from MEDUSA
end if ! use_MEDUSA and sedflx_num not 0

! ======================================================================================
!********************************* LOOP STARTS *****************************************			

  do n=1, myDim_nod2D  ! needs exchange_nod in the end
!     if (ulevels_nod2D(n)>1) cycle 
!            nzmin = ulevels_nod2D(n)

     !!---- Number of vertical layers
     nzmax = nlevels_nod2D(n)-1

     !!---- This is needed for piston velocity 
     Loc_ice_conc = a_ice(n) 

     !!---- Mean sea level pressure ! MERGE
#if defined (__oasis) 
!!   MB: This is an ad-hoc patch for AWIESM-2.1 and needs to be improved:
!!   We should consider air pressure provided by ECHAM.
     Loc_slp           = pa2atm
#else
     Loc_slp = press_air(n)
#endif 

     !!---- Benthic layers
     LocBenthos(1:benthos_num) = Benthos(n,1:benthos_num)
!CV: It is not clear to me whether this is still needed

     !!---- Local conc of [H+]-ions from last time time step. Stored in LocVar
     !!---- used as first guess for H+ conc.in subroutine CO2flux (provided by recom_init)
     Hplus = GloHplus(n)                                  

     !!---- Interpolated wind from atmospheric forcing 
     !!---- temporarily stored in module LocVar
#if defined (__oasis)
!    Derive 10m-wind speed from wind stress fields, see module recom_ciso.
!    This is an ad-hoc solution as long as 10m-winds are not handled from OASIS.
     Uloc              = wind_10(stress_atmoce_x(n), stress_atmoce_y(n))
#else
     ULoc = sqrt(u_wind(n)**2+v_wind(n)**2)
#endif

     !!---- Atmospheric CO2 in LocVar                                                                        
     LocAtmCO2         = AtmCO2(month)

! Update of prognostic atmospheric CO2 values
     if (use_atbox) then
       LocAtmCO2                   = x_co2atm(1)
       if (ciso) then
         LocAtmCO2_13              = x_co2atm_13(1)
         if (ciso_14) LocAtmCO2_14 = x_co2atm_14(1)
       end if
     else
! Consider prescribed atmospheric CO2 values
       if (ciso) then
         LocAtmCO2_13              = AtmCO2_13(month)
         if (ciso_14) then
!          Latitude of nodal point n 
           lat_val = geo_coord_nod2D(2,n) / rad
!          Zonally binned NH / SH / TZ 14CO2 input values
           LocAtmCO2_14 = AtmCO2_14(lat_zone(lat_val), month)
         end if
       end if
     end if  ! use_atbox

     if (ciso) then
       r_atm_13                    = LocAtmCO2_13(1) / LocAtmCO2(1)
       if (ciso_14) r_atm_14       = LocAtmCO2_14(1) / LocAtmCO2(1)
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

     !!---- Allocate 3D diagnostics    ! Comment Miriam (2/2024): changed grazing from a 3D to a 2D diagnostic while completing all grazing fluxes
!     allocate(vertgrazmeso_tot(nl-1), vertgrazmeso_n(nl-1), vertgrazmeso_d(nl-1))
!     vertgrazmeso_tot = 0.d0
!     vertgrazmeso_n   = 0.d0
!     vertgrazmeso_d   = 0.d0

!#if defined (__coccos)
!     allocate(vertgrazmeso_c(nl-1))
!     vertgrazmeso_c   = 0.d0
!#endif

     allocate(vertrespmeso(nl-1))
     vertrespmeso  = 0.d0

#if defined (__3Zoo2Det)
     allocate(vertrespmacro(nl-1), vertrespmicro(nl-1))
     vertrespmacro = 0.d0
     vertrespmicro = 0.d0
#endif

     allocate(vertcalcdiss(nl-1), vertcalcif(nl-1))
     vertcalcdiss = 0.d0
     vertcalcif   = 0.d0

     allocate(vertaggn(nl-1), vertaggd(nl-1))
     vertaggn = 0.d0
     vertaggd = 0.d0

     allocate(vertdocexn(nl-1), vertdocexd(nl-1))
     vertdocexn = 0.d0
     vertdocexd = 0.d0

     allocate(vertrespn(nl-1), vertrespd(nl-1))
     vertrespn = 0.d0
     vertrespd = 0.d0

#if defined (__coccos)
     allocate(vertaggc(nl-1), vertdocexc(nl-1), vertrespc(nl-1))
     vertaggc = 0.d0
     vertdocexc = 0.d0
     vertrespc = 0.d0
#endif

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

#if defined (__coccos)
     allocate(vertNPPc(nl-1), vertGPPc(nl-1), vertNNAc(nl-1), vertChldegc(nl-1)) 
     vertNPPc = 0.d0
     vertGPPc = 0.d0
     vertNNAc = 0.d0
     Chldegc = 0.d0
#endif

     if (Grazing_detritus) then
        allocate(vertgrazmeso_tot(nl-1), vertgrazmeso_n(nl-1), vertgrazmeso_d(nl-1))
        vertgrazmeso_tot = 0.d0
        vertgrazmeso_n   = 0.d0
        vertgrazmeso_d   = 0.d0
#if defined (__coccos)
        allocate(vertgrazmeso_c(nl-1))
        vertgrazmeso_c   = 0.d0
#endif
        allocate(vertgrazmeso_det(nl-1))
        vertgrazmeso_det = 0.d0
#if defined (__3Zoo2Det)
        allocate(vertgrazmeso_mic(nl-1), vertgrazmeso_det2(nl-1))
        vertgrazmeso_mic  = 0.d0
        vertgrazmeso_det2 = 0.d0
        allocate(vertgrazmacro_tot(nl-1), vertgrazmacro_n(nl-1), vertgrazmacro_d(nl-1))
        vertgrazmacro_tot = 0.d0
        vertgrazmacro_n   = 0.d0
        vertgrazmacro_d   = 0.d0
#if defined (__coccos)
        allocate(vertgrazmacro_c(nl-1))
        vertgrazmacro_c   = 0.d0
#endif
        allocate(vertgrazmacro_mes(nl-1), vertgrazmacro_det(nl-1), vertgrazmacro_mic(nl-1), vertgrazmacro_det2(nl-1))
        vertgrazmacro_mes = 0.d0
        vertgrazmacro_det = 0.d0
        vertgrazmacro_mic = 0.d0
        vertgrazmacro_det2= 0.d0
        allocate(vertgrazmicro_tot(nl-1), vertgrazmicro_n(nl-1), vertgrazmicro_d(nl-1))
        vertgrazmicro_tot = 0.d0
        vertgrazmicro_n   = 0.d0
        vertgrazmicro_d   = 0.d0
#if defined (__coccos)
        allocate(vertgrazmicro_c(nl-1))
        vertgrazmicro_c   = 0.d0
#endif
#endif
     endif !Grazing_detritus
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
     GPPn(n) = locGPPn
     NNAn(n) = locNNAn
     Chldegn(n) = locChldegn

     NPPd(n) = locNPPd
     GPPd(n) = locGPPd
     NNAd(n) = locNNAd
     Chldegd(n) = locChldegd

#if defined (__coccos)
     NPPc(n) = locNPPc
     GPPc(n) = locGPPc
     NNAc(n) = locNNAc
     Chldegc(n) = locChldegc
#endif

     if (Grazing_detritus) then
     ! Mesozooplankton
     grazmeso_tot(n) = locgrazmeso_tot
     grazmeso_n(n)   = locgrazmeso_n
     grazmeso_d(n)   = locgrazmeso_d
#if defined (__coccos)
     grazmeso_c(n)   = locgrazmeso_c
#endif
     grazmeso_det(n) = locgrazmeso_det
#if defined (__3Zoo2Det)
     grazmeso_mic(n) = locgrazmeso_mic
     grazmeso_det2(n)= locgrazmeso_det2
#endif

#if defined (__3Zoo2Det)
     ! Macrozooplankton
     grazmacro_tot(n) = locgrazmacro_tot
     grazmacro_n(n)   = locgrazmacro_n
     grazmacro_d(n)   = locgrazmacro_d
#if defined (__coccos)
     grazmacro_c(n)   = locgrazmacro_c
#endif
     grazmacro_mes(n) = locgrazmacro_mes
     grazmacro_det(n) = locgrazmacro_det
     grazmacro_mic(n) = locgrazmacro_mic
     grazmacro_det2(n)= locgrazmacro_det2

     ! Microzooplankton
     grazmicro_tot(n) = locgrazmicro_tot
     grazmicro_n(n)   = locgrazmicro_n
     grazmicro_d(n)   = locgrazmicro_d
#if defined (__coccos)
     grazmicro_c(n)   = locgrazmicro_c
#endif
     
#endif
     endif
     
     !!---- Updating 3D diagnostics
!     grazmeso_tot(1:nzmax,n) = vertgrazmeso_tot(1:nzmax)
!     grazmeso_n(1:nzmax,n)   = vertgrazmeso_n(1:nzmax)
!     grazmeso_d(1:nzmax,n)   = vertgrazmeso_d(1:nzmax)
!#if defined (__coccos)
!     grazmeso_c(1:nzmax,n)   = vertgrazmeso_c(1:nzmax)
!#endif

     respmeso(1:nzmax,n)     = vertrespmeso(1:nzmax)
#if defined (__3Zoo2Det)
     respmacro(1:nzmax,n)    = vertrespmacro(1:nzmax)
     respmicro(1:nzmax,n)    = vertrespmicro(1:nzmax)
#endif
     calcdiss(1:nzmax,n)     = vertcalcdiss(1:nzmax)
     calcif(1:nzmax,n)       = vertcalcif(1:nzmax)

     aggn(1:nzmax,n)         = vertaggn(1:nzmax)
     docexn(1:nzmax,n)       = vertdocexn(1:nzmax)
     respn(1:nzmax,n)        = vertrespn(1:nzmax)
     NPPn3D(1:nzmax,n)       = vertNPPn(1:nzmax)

     aggd(1:nzmax,n)         = vertaggd(1:nzmax)
     respd(1:nzmax,n)        = vertrespd(1:nzmax)
     docexd(1:nzmax,n)       = vertdocexd(1:nzmax)
     NPPd3D(1:nzmax,n)       = vertNPPd(1:nzmax)

#if defined (__coccos)
     aggc(1:nzmax,n)         = vertaggc(1:nzmax)
     docexc(1:nzmax,n)       = vertdocexc(1:nzmax)
     respc(1:nzmax,n)        = vertrespc(1:nzmax)
     NPPc3D(1:nzmax,n)       = vertNPPc(1:nzmax)
#endif

     !!---- Deallocating 2D diagnostics
     deallocate(vertNPPn,vertGPPn,vertNNAn,vertChldegn) 
     deallocate(vertNPPd,vertGPPd,vertNNAd,vertChldegd) 
#if defined (__coccos)
     deallocate(vertNPPc,vertGPPc,vertNNAc,vertChldegc) 
#endif

     if (Grazing_detritus) then
        deallocate(vertgrazmeso_tot, vertgrazmeso_n,  vertgrazmeso_d)
#if defined (__coccos)
        deallocate(vertgrazmeso_c)
#endif
        deallocate(vertgrazmeso_det)
#if defined (__3Zoo2Det)
        deallocate(vertgrazmeso_mic, vertgrazmeso_det2)
        deallocate(vertgrazmacro_tot, vertgrazmacro_n, vertgrazmacro_d)
#if defined (__coccos)
        deallocate(vertgrazmacro_c)
#endif
        deallocate(vertgrazmacro_mes, vertgrazmacro_det, vertgrazmacro_mic, vertgrazmacro_det2)
        deallocate(vertgrazmicro_tot, vertgrazmicro_n, vertgrazmicro_d)
#if defined (__coccos)
        deallocate(vertgrazmicro_c)
#endif        
#endif
     endif ! Grazing_detritus
        

     !!---- Deallocating 3D Diagnistics
     !     deallocate(vertgrazmeso_tot, vertgrazmeso_n, vertgrazmeso_d, vertrespmeso)
     deallocate(vertrespmeso)
#if defined (__3Zoo2Det)
     deallocate(vertrespmacro, vertrespmicro)
#endif
     deallocate(vertcalcdiss, vertcalcif)
     deallocate(vertaggn, vertdocexn, vertrespn)
     deallocate(vertaggd, vertdocexd, vertrespd)
#if defined (__coccos)
!    deallocate(vertgrazmeso_c)
     deallocate(vertaggc, vertdocexc, vertrespc)
#endif
endif

     GloPCO2surf(n)               = pco2surf(1)
     GlodPCO2surf(n)              = dpco2surf(1)

     GloCO2flux(n)                = dflux(1)
     GloCO2flux_seaicemask(n)     = co2flux_seaicemask(1)                 !  [mmol/m2/s]
     GloO2flux_seaicemask(n)      = o2flux_seaicemask(1)                  !  [mmol/m2/s]
     if (ciso) then
        GloCO2flux_seaicemask_13(n)     = co2flux_seaicemask_13(1)        !  [mmol/m2/s]
        if (ciso_14) then
            GloCO2flux_seaicemask_14(n) = co2flux_seaicemask_14(1)        !  [mmol/m2/s]
        end if
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
#if defined (__coccos)
    call exchange_nod(NPPc)
    call exchange_nod(GPPc)
    call exchange_nod(NNAc)
    call exchange_nod(Chldegc)
#endif
    call exchange_nod(grazmeso_tot)
    call exchange_nod(grazmeso_n)
    call exchange_nod(grazmeso_d)
#if defined (__coccos)
    call exchange_nod(grazmeso_c)
#endif
    call exchange_nod(grazmeso_det)
#if defined (__3Zoo2Det)
    call exchange_nod(grazmeso_mic)
    call exchange_nod(grazmeso_det2)
    call exchange_nod(grazmacro_tot)
    call exchange_nod(grazmacro_n)
    call exchange_nod(grazmacro_d)
#if defined (__coccos)
    call exchange_nod(grazmacro_c)
#endif
    call exchange_nod(grazmacro_mes)
    call exchange_nod(grazmacro_det)
    call exchange_nod(grazmacro_mic)
    call exchange_nod(grazmacro_det2)
    call exchange_nod(grazmicro_tot)
    call exchange_nod(grazmicro_n)
    call exchange_nod(grazmicro_d)
#if defined (__coccos)
    call exchange_nod(grazmicro_c)
#endif
#endif
  endif

  call exchange_nod(GloPCO2surf)	
  call exchange_nod(GloCO2flux)	
  call exchange_nod(GloCO2flux_seaicemask)
  if (ciso) then
    call exchange_nod(GloPCO2surf_13)
    call exchange_nod(GloCO2flux_13)
    call exchange_nod(GloCO2flux_seaicemask_13)
    if (ciso_14) then
      call exchange_nod(GloPCO2surf_14)
      call exchange_nod(GloCO2flux_14)
      call exchange_nod(GloCO2flux_seaicemask_14)
    end if 
  end if

  do n=1, benthos_num  !4
    call exchange_nod(GlodecayBenthos(:,n))
  end do
 
  call exchange_nod(GloO2flux_seaicemask)	
  call exchange_nod(GloHplus)	
  call exchange_nod(AtmFeInput)	
  call exchange_nod(AtmNInput)	
!  call exchange_nod(DenitBen)	

  call exchange_nod(PAR3D)
  call exchange_nod(CO23D)
  call exchange_nod(pH3D)
  call exchange_nod(pCO23D)
  call exchange_nod(HCO33D)
  call exchange_nod(CO33D)
  call exchange_nod(OmegaC3D)
  call exchange_nod(kspc3D)
  call exchange_nod(rhoSW3D)

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
