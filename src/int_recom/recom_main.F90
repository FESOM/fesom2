! 24.03.2023
! OG
!===============================================================================
! Main REcoM 
module recom_interface
    interface
        subroutine recom(ice, dynamics, tracers, partit, mesh)
            use mod_mesh
            USE MOD_PARTIT
            USE MOD_PARSUP
            use mod_tracer
            use MOD_DYN
            use MOD_ICE
            type(t_dyn)   , intent(inout), target :: dynamics
            type(t_ice)   , intent(inout), target :: ice
            type(t_tracer), intent(inout), target :: tracers
            type(t_partit), intent(inout), target :: partit
            type(t_mesh)  , intent(inout), target :: mesh
        end subroutine
    end interface
end module

module bio_fluxes_interface
    interface
        subroutine bio_fluxes(tracers, partit, mesh)
            use recom_declarations
            use recom_locvar
            use recom_glovar
            use recom_config

            use mod_mesh
            USE MOD_PARTIT
            USE MOD_PARSUP
            use mod_tracer

            use g_config
            use o_arrays
            use g_comm_auto
            use g_forcing_arrays
            use g_support

            type(t_tracer), intent(inout), target :: tracers
            type(t_partit), intent(inout), target :: partit
            type(t_mesh), intent(inout), target :: mesh

      end subroutine
    end interface
end module

subroutine recom(ice, dynamics, tracers, partit, mesh)
    use g_config
    use MOD_MESH
    use MOD_TRACER
    use MOD_DYN
    USE MOD_ICE
    use o_ARRAYS
    use o_PARAM
    USE MOD_PARTIT
    USE MOD_PARSUP

    use recom_declarations
    use bio_fluxes_interface
    use recom_locvar
    use recom_glovar
    use recom_config
    use recom_ciso
    use g_clock
    use g_forcing_arrays, only: press_air, u_wind, v_wind, shortwave
    use g_comm_auto
    IMPLICIT NONE

    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    type(t_ice)   , intent(inout), target :: ice
    !___________________________________________________________________________

    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice

! ======================================================================================
!! Depth information

!! zbar(nl) allocate the array for storing the standard depths, it is negativ
!! Z(nl-1)  mid-depths of cells

!! max. number of levels at node n
!! nzmax = nlevels_nod2D(n)
!! u_ice and v_ice are at nodes
!! u_w, v_w are at nodes (interpolated from elements)
!! u_wind and v_wind are always at nodes
! ======================================================================================

    real(kind=8)               :: SW, Loc_slp
    integer                    :: tr_num, num_tracers
    integer                    :: nz, n, nzmin, nzmax, nu1, nl1
    integer                    :: idiags

    real(kind=8)               :: Sali
    logical                    :: do_update = .false. 

    real(kind=8),  allocatable :: Temp(:), Sali_depth(:), zr(:), PAR(:)
    real(kind=8),  allocatable :: C(:,:)

    !! * Mocsy *
    real(kind=8),  allocatable :: CO2_watercolumn(:)
    real(kind=8),  allocatable :: pH_watercolumn(:)
    real(kind=8),  allocatable :: pCO2_watercolumn(:)
    real(kind=8),  allocatable :: HCO3_watercolumn(:)

    !! * Diss *
    real(kind=8),  allocatable :: CO3_watercolumn(:)
    real(kind=8),  allocatable :: OmegaC_watercolumn(:)
    real(kind=8),  allocatable :: kspc_watercolumn(:)
    real(kind=8),  allocatable :: rhoSW_watercolumn(:)
    real(kind=WP)              :: ttf_rhs_bak (mesh%nl-1, tracers%num_tracers) ! local variable ! OG - tra_diag

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    allocate(Temp(nl-1), Sali_depth(nl-1), zr(nl-1) , PAR(nl-1))
    allocate(C(nl-1, bgc_num))
    allocate(CO2_watercolumn(nl-1), pH_watercolumn(nl-1), pCO2_watercolumn(nl-1) , HCO3_watercolumn(nl-1))
    allocate(CO3_watercolumn(nl-1), OmegaC_watercolumn(nl-1), kspc_watercolumn(nl-1) , rhoSW_watercolumn(nl-1))

    !< ice concentration [0 to 1]

    a_ice       => ice%data(1)%values(:)
    num_tracers = tracers%num_tracers

    !< alkalinity restoring to climatology
    !< virtual flux is possible

    if (restore_alkalinity) call bio_fluxes(tracers, partit, mesh)
    if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> bio_fluxes'//achar(27)//'[0m'

  if (use_atbox) then    ! MERGE
! Prognostic atmospheric isoCO2
    call recom_atbox(partit,mesh)
!   optional I/O of isoCO2 and inferred cosmogenic 14C production; this may cost some CPU time
    if (ciso .and. ciso_14) then
      call annual_event(do_update)
      if (do_update .and. mype==0) write (*, fmt = '(a50,2x,i6,4(2x,f6.2))') &
                                         'Year, xCO2 (ppm), cosmic 14C flux (at / cm² / s):', &
                                          yearold, x_co2atm(1), x_co2atm_13(1), x_co2atm_14(1), cosmic_14(1) * production_rate_to_flux_14
    end if
  end if
! ======================================================================================
!********************************* LOOP STARTS *****************************************

    do n=1, myDim_nod2D  ! needs exchange_nod in the end
!     if (ulevels_nod2D(n)>1) cycle
!       nzmin = ulevels_nod2D(n)

        !!---- Number of vertical layers
        nzmax = nlevels_nod2D(n)-1

        !!---- This is needed for piston velocity
        Loc_ice_conc = a_ice(n)

        !!---- Mean sea level pressure
#if defined (__oasis)
!!      MB: This is an ad-hoc patch for AWIESM-2.1 and needs to be improved:
!!      We should consider air pressure provided by ECHAM.
        Loc_slp = pa2atm
#else
        Loc_slp = press_air(n)
#endif

        !!---- Benthic layers
        LocBenthos(1:benthos_num) = Benthos(n,1:benthos_num)

        !!---- Local conc of [H+]-ions from last time step. Decleared and saved in LocVar.
        !!---- used as first guess for H+ conc. in subroutine CO2flux (provided by recom_init)
        Hplus = GloHplus(n)

        !!---- Interpolated wind from atmospheric forcing
#if defined (__oasis)
        !! Derive 10m-wind speed from wind stress fields, see module recom_ciso.
        !! This is an ad-hoc solution as long as 10m-winds are not handled from OASIS.
        Uloc = wind_10(stress_atmoce_x(n), stress_atmoce_y(n))
#else
        ULoc = sqrt(u_wind(n)**2+v_wind(n)**2)
#endif

        !!---- Atmospheric CO2 in LocVar
        LocAtmCO2 = AtmCO2(month)

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
        Temp(1:nzmax) = tracers%data(1)%values(1:nzmax, n)

        !!---- Surface salinity
        Sali                = tracers%data(2)%values(1, n)
        Sali_depth(1:nzmax) = tracers%data(2)%values(1:nzmax, n)


        !!---- CO2 in the watercolumn

        !! * Mocsy *
        CO2_watercolumn(1:nzmax)    = CO23D(1:nzmax, n)
        pH_watercolumn(1:nzmax)     = pH3D(1:nzmax, n)
        pCO2_watercolumn(1:nzmax)   = pCO23D(1:nzmax, n)
        HCO3_watercolumn(1:nzmax)   = HCO33D(1:nzmax, n)
        !! * Diss *
        CO3_watercolumn(1:nzmax)    = CO33D(1:nzmax, n)
        OmegaC_watercolumn(1:nzmax) = OmegaC3D(1:nzmax, n)
        kspc_watercolumn(1:nzmax)   = kspc3D(1:nzmax, n)
        rhoSW_watercolumn(1:nzmax)  = rhoSW3D(1:nzmax, n)

        !!---- Biogeochemical tracers
        do tr_num = num_tracers-bgc_num+1, num_tracers
            C(1:nzmax, tr_num-2) = tracers%data(tr_num)%values(1:nzmax, n)
        end do

        ttf_rhs_bak = 0.0 ! OG - tra_diag


           do tr_num=1, num_tracers
        if (tracers%data(tr_num)%ltra_diag) then ! OG - tra_diag
              ttf_rhs_bak(1:nzmax,tr_num) = tracers%data(tr_num)%values(1:nzmax, n)
        end if
           end do

        !!---- Depth of the nodes in the water column
        zr(1:nzmax) = Z_3d_n(1:nzmax, n)

        !!---- The PAR in the local water column is initialized
        PAR(1:nzmax) = 0.d0

        !!---- a_ice(row): Ice concentration in the local node
        FeDust = GloFeDust(n) * (1.d0 - a_ice(n)) * dust_sol
        NDust  = GloNDust(n)  * (1.d0 - a_ice(n))

        if (Diags) then

        !! * Allocate 3D diagnostics *
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

            !! * Allocate 2D diagnostics *
            allocate(vertNPPn(nl-1), vertGPPn(nl-1), vertNNAn(nl-1), vertChldegn(nl-1)) 
            vertNPPn = 0.d0
            vertGPPn = 0.d0
            vertNNAn = 0.d0
            vertChldegn  = 0.d0

            allocate(vertNPPd(nl-1), vertGPPd(nl-1), vertNNAd(nl-1), vertChldegd(nl-1)) 
            vertNPPd = 0.d0
            vertGPPd = 0.d0
            vertNNAd = 0.d0
            vertChldegd  = 0.d0

#if defined (__coccos)
            allocate(vertNPPc(nl-1), vertGPPc(nl-1), vertNNAc(nl-1), vertChldegc(nl-1)) 
            vertNPPc = 0.d0
            vertGPPc = 0.d0
            vertNNAc = 0.d0
            vertChldegc  = 0.d0
#endif
        end if

        if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> REcoM_Forcing'//achar(27)//'[0m'

! ======================================================================================
!******************************** RECOM FORCING ****************************************
        call REcoM_Forcing(zr, n, nzmax, C, SW, Loc_slp , Temp, Sali, Sali_depth &
           , CO2_watercolumn                                     & ! NEW MOCSY CO2 for the whole watercolumn
           , pH_watercolumn                                      & ! NEW MOCSY pH for the whole watercolumn
           , pCO2_watercolumn                                    & ! NEW MOCSY pCO2 for the whole watercolumn
           , HCO3_watercolumn                                    & ! NEW MOCSY HCO3 for the whole watercolumn
           , CO3_watercolumn                                     & ! NEW DISS CO3 for the whole watercolumn
           , OmegaC_watercolumn                                  & ! NEW DISS OmegaC for the whole watercolumn
           , kspc_watercolumn                                    & ! NEW DISS stoichiometric solubility product for calcite [mol^2/kg^2]
           , rhoSW_watercolumn                                   & ! NEW DISS in-situ density of seawater [mol/m^3]
                           , PAR, ice, dynamics, tracers, partit, mesh)

        do tr_num = num_tracers-bgc_num+1, num_tracers !bgc_num+2
            tracers%data(tr_num)%values(1:nzmax, n) = C(1:nzmax, tr_num-2)
        end do

        ! recom_sms

           do tr_num=1, num_tracers
        if (tracers%data(tr_num)%ltra_diag) then ! OG - tra_diag
             tracers%work%tra_recom_sms(1:nzmax,n,tr_num) = tracers%data(tr_num)%values(1:nzmax, n) - ttf_rhs_bak(1:nzmax,tr_num)
             !if (mype==0)  print *,  tra_recom_sms(:,:,tr_num)
        end if

           end do

        !!---- Local variables that have been changed during the time-step are stored so they can be saved
        Benthos(n,1:benthos_num) = LocBenthos(1:benthos_num)
        GlodecayBenthos(n, 1:benthos_num) = decayBenthos(1:benthos_num)/SecondsPerDay ! convert from [mmol/m2/d] to [mmol/m2/s]

        if (Diags) then

            !! * Update 2D diagnostics *
            NPPn(n) = locNPPn
            NPPd(n) = locNPPd
            NPPc(n) = locNPPc
            GPPn(n) = locGPPn
            GPPd(n) = locGPPd
            GPPc(n) = locGPPc
            NNAn(n) = locNNAn
            NNAd(n) = locNNAd
            NNAc(n) = locNNAc
            Chldegn(n) = locChldegn
            Chldegd(n) = locChldegd
            Chldegc(n) = locChldegc

            !! * Update 3D diagnostics *
            respmeso     (1:nzmax,n) = vertrespmeso     (1:nzmax)
#if defined (__3Zoo2Det)
            respmacro    (1:nzmax,n) = vertrespmacro    (1:nzmax)
            respmicro    (1:nzmax,n) = vertrespmicro    (1:nzmax)
#endif
            calcdiss     (1:nzmax,n) = vertcalcdiss     (1:nzmax)
            calcif       (1:nzmax,n) = vertcalcif       (1:nzmax)

            aggn         (1:nzmax,n) = vertaggn         (1:nzmax)
            docexn       (1:nzmax,n) = vertdocexn       (1:nzmax)
            respn        (1:nzmax,n) = vertrespn        (1:nzmax)
            NPPn3D       (1:nzmax,n) = vertNPPn         (1:nzmax)

            aggd         (1:nzmax,n) = vertaggd         (1:nzmax)
            docexd       (1:nzmax,n) = vertdocexd       (1:nzmax)
            respd        (1:nzmax,n) = vertrespd        (1:nzmax)
            NPPd3D       (1:nzmax,n) = vertNPPd         (1:nzmax)

#if defined (__coccos)
            aggc         (1:nzmax,n) = vertaggc         (1:nzmax)
            docexc       (1:nzmax,n) = vertdocexc       (1:nzmax)
            respc        (1:nzmax,n) = vertrespc        (1:nzmax)
            NPPc3D       (1:nzmax,n) = vertNPPc         (1:nzmax)
#endif

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> ciso after REcoM_Forcing'//achar(27)//'[0m'

            !! * Deallocating 2D diagnostics *
            deallocate(vertNPPn, vertGPPn, vertNNAn, vertChldegn) 
            deallocate(vertNPPd, vertGPPd, vertNNAd, vertChldegd)
#if defined (__coccos)
            deallocate(vertNPPc, vertGPPc, vertNNAc, vertChldegc) 
#endif 

            !! * Deallocating 3D Diagnostics *
            deallocate(vertrespmeso)
#if defined (__3Zoo2Det)
            deallocate(vertrespmacro, vertrespmicro)
#endif
            deallocate(vertcalcdiss, vertcalcif)
            deallocate(vertaggn, vertdocexn, vertrespn)
            deallocate(vertaggd, vertdocexd, vertrespd)
#if defined (__coccos)
!           deallocate(vertgrazmeso_c)
            deallocate(vertaggc, vertdocexc, vertrespc)
#endif

        end if 

        AtmFeInput(n)            = FeDust
        AtmNInput(n)             = NDust
        GloHplus(n)              = ph(1)

        GloPCO2surf(n)           = pco2surf(1)
        GlodPCO2surf(n)          = dpco2surf(1)
        GloCO2flux(n)            = dflux(1)                   !  [mmol/m2/d]
        GloCO2flux_seaicemask(n) = co2flux_seaicemask(1)      !  [mmol/m2/s]
        GloO2flux_seaicemask(n)  = o2flux_seaicemask(1)       !  [mmol/m2/s]
    if (ciso) then
        GloCO2flux_seaicemask_13(n)     = co2flux_seaicemask_13(1)        !  [mmol/m2/s]
        if (ciso_14) then
            GloCO2flux_seaicemask_14(n) = co2flux_seaicemask_14(1)        !  [mmol/m2/s]
        end if
     end if
        GloO2flux(n)             = oflux(1)                   !  [mmol/m2/d]

        PAR3D(1:nzmax,n)         = PAR(1:nzmax)

        !! * Mocsy *
        CO23D(1:nzmax,n)         = CO2_watercolumn(1:nzmax)
        pH3D(1:nzmax,n)          = pH_watercolumn(1:nzmax)
        pCO23D(1:nzmax,n)        = pCO2_watercolumn(1:nzmax)
        HCO33D(1:nzmax,n)        = HCO3_watercolumn(1:nzmax)

        !! * Diss *
        CO33D(1:nzmax,n)         = CO3_watercolumn(1:nzmax)
        OmegaC3D(1:nzmax,n)      = OmegaC_watercolumn(1:nzmax)
        kspc3D(1:nzmax,n)        = kspc_watercolumn(1:nzmax)
        rhoSW3D(1:nzmax,n)       = rhoSW_watercolumn(1:nzmax)
    end do

! ======================================================================================
!************************** EXCHANGE NODAL INFORMATION *********************************

    do tr_num=num_tracers-bgc_num+1, num_tracers
        call exchange_nod(tracers%data(tr_num)%values(:,:), partit)
    end do

    call exchange_nod(GloPCO2surf, partit)
    call exchange_nod(GlodPCO2surf, partit)
    call exchange_nod(GloCO2flux, partit)
    call exchange_nod(GloCO2flux_seaicemask, partit)

    call exchange_nod(GloO2flux_seaicemask, partit)
    if (ciso) then
        call exchange_nod(GloPCO2surf_13, partit)
        call exchange_nod(GloCO2flux_13, partit)
        call exchange_nod(GloCO2flux_seaicemask_13, partit)
        if (ciso_14) then
            call exchange_nod(GloPCO2surf_14, partit)
            call exchange_nod(GloCO2flux_14, partit)
            call exchange_nod(GloCO2flux_seaicemask_14, partit)
        end if 
    end if
    do n=1, benthos_num
        call exchange_nod(Benthos(:,n), partit)
    end do

    if (Diags) then
        call exchange_nod(NPPn, partit)
        call exchange_nod(NPPd, partit)
        call exchange_nod(GPPn, partit)
        call exchange_nod(GPPd, partit)
        call exchange_nod(NNAn, partit)
        call exchange_nod(NNAd, partit)
        call exchange_nod(Chldegn, partit)
        call exchange_nod(Chldegd, partit)
#if defined (__coccos)
        call exchange_nod(NPPc, partit)
        call exchange_nod(GPPc, partit)
        call exchange_nod(NNAc, partit)
        call exchange_nod(Chldegc, partit)
#endif
    endif

    do n=1, benthos_num
        call exchange_nod(GlodecayBenthos(:,n), partit)
    end do

    call exchange_nod(GloHplus, partit)
    call exchange_nod(AtmFeInput, partit)
    call exchange_nod(AtmNInput, partit)

    call exchange_nod(PAR3D, partit)

    call exchange_nod(CO23D, partit)
    call exchange_nod(pH3D, partit)
    call exchange_nod(pCO23D, partit)
    call exchange_nod(HCO33D, partit)
    call exchange_nod(CO33D, partit)
    call exchange_nod(OmegaC3D, partit)
    call exchange_nod(kspc3D, partit)
    call exchange_nod(rhoSW3D, partit)

end subroutine recom

! ======================================================================================
! Alkalinity restoring to climatology                                 	     
! ======================================================================================
subroutine bio_fluxes(tracers, partit, mesh)

    use recom_declarations
    use recom_locvar
    use recom_glovar
    use recom_config

    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use mod_tracer

    use g_config
    use o_arrays
    use g_comm_auto
    use g_forcing_arrays
    use g_support

    implicit none
    integer                               :: n, elem, elnodes(3),n1
    real(kind=WP)                         :: ralk, net

    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh

    !___________________________________________________________________________
    real(kind=WP), dimension(:,:), pointer :: alkalinity

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    alkalinity => tracers%data(2+ialk)%values(:,:) ! 1 temp, 2 salt, 3 din, 4 dic, 5 alk
!___________________________________________________________________
! on freshwater inflow/outflow or virtual alkalinity:
  ! 1. In zlevel & zstar the freshwater flux is applied in the update of the 
  ! ssh matrix when solving the continuity equation of vertically 
  ! integrated flow. The alkalinity concentration in the first layer will 
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
!        if (ref_alk_local) ralk = tracers%data(2+ialk)%values(1, n)
!        virtual_alk(n)=ralk*water_flux(n) 
!     end do
!  end if

    !___________________________________________________________________________
    ! Balance alkalinity restoring to climatology
    do n=1, myDim_nod2D+eDim_nod2D
!        relax_alk(n)=surf_relax_Alk * (Alk_surf(n) - tracers%data(2+ialk)%values(1, n)) 
!        relax_alk(n)=surf_relax_Alk * (Alk_surf(n) - alkalinity(ulevels_nod2d(n),n)
        relax_alk(n)=surf_relax_Alk * (Alk_surf(n) - alkalinity(1, n))
    end do

  ! 2. virtual alkalinity flux
!  if (use_virt_alk) then ! is already zero otherwise
!     call integrate_nod(virtual_alk, net, partit, mesh)
!     virtual_alk=virtual_alk-net/ocean_area
!  end if


  ! 3. restoring to Alkalinity climatology
    call integrate_nod(relax_alk, net, partit, mesh)

    relax_alk=relax_alk-net/ocean_area  ! at ocean surface layer

end subroutine bio_fluxes
