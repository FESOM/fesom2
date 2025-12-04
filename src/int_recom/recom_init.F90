! ------------
! 23.03.2023
! OG   
!===============================================================================
! allocate & initialise arrays for REcoM
module recom_init_interface
    interface
        subroutine recom_init(tracers, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in)  ,  target :: mesh
        end subroutine
    end interface
end module
!
!
!_______________________________________________________________________________
subroutine recom_init(tracers, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    use g_clock
    use REcoM_declarations
    use REcoM_GloVar
    use REcoM_locVar
    use recom_config
    use REcoM_ciso
    implicit none
#include "netcdf.inc"
    !___________________________________________________________________________
    ! pointer on necessary derived types 
    integer                                 :: n, k, row, nzmin, nzmax, i, id
    integer                                 :: elem_size, node_size, num_tracers

    real(kind=WP)                           :: locDINmax, locDINmin, locDICmax, locDICmin, locAlkmax, glo
    real(kind=WP)                           :: locAlkmin, locDSimax, locDSimin, locDFemax, locDFemin
    real(kind=WP)                           :: locO2max, locO2min
    real(kind=WP)                           :: locDICremin, locDICremax ! initialization of DIC remin (added by Sina)

    type(t_tracer), intent(inout), target   :: tracers
    type(t_partit), intent(inout), target   :: partit
    type(t_mesh),   intent(in) ,   target   :: mesh

    ! After reading tracer namelist - validate actual IDs
    integer, dimension(tracers%num_tracers) :: tracer_id_array

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    elem_size   = myDim_elem2D+eDim_elem2D
    node_size   = myDim_nod2D+eDim_nod2D
    num_tracers = tracers%num_tracers


!! *** Allocate and initialize ***

    !! * Fe and N deposition as surface boundary condition *
    allocate(GloFeDust             ( node_size ))
    allocate(AtmFeInput            ( node_size ))
    allocate(GloNDust              ( node_size ))
    allocate(AtmNInput             ( node_size ))

    !! * River nutrients as surface boundary condition *
    allocate(RiverDIN2D            ( node_size ))
    allocate(RiverDON2D            ( node_size ))
    allocate(RiverDOC2D            ( node_size ))
    allocate(RiverDSi2D            ( node_size ))
    allocate(RiverDIC2D            ( node_size ))
    allocate(RiverAlk2D            ( node_size ))
    allocate(RiverFe               ( node_size ))

    !! * Erosion nutrients as surface boundary condition *
    allocate(ErosionTON2D          ( node_size ))
    allocate(ErosionTOC2D          ( node_size ))
    allocate(ErosionTSi2D          ( node_size ))

    !! * Alkalinity restoring to climatology *
    allocate(relax_alk             ( node_size ))
    allocate(virtual_alk           ( node_size ))

    allocate(cosAI                 ( node_size ))
    allocate(GloPCO2surf           ( node_size ))
    allocate(GloCO2flux            ( node_size ))
    allocate(GloO2flux             ( node_size ))
    allocate(GloCO2flux_seaicemask ( node_size ))
    allocate(GloO2flux_seaicemask  ( node_size ))
    allocate(GlodPCO2surf          ( node_size ))
    allocate(GlodecayBenthos       ( node_size, benthos_num ))
    allocate(Benthos               ( node_size, benthos_num ))
    allocate(Benthos_tr            ( node_size, benthos_num, num_tracers )) ! kh 25.03.22 buffer per tracer index
    allocate(GloHplus              ( node_size ))
    allocate(DenitBen              ( node_size ))
    allocate(PistonVelocity        ( node_size ))
    allocate(alphaCO2              ( node_size ))


    allocate(LocBenthos            ( benthos_num ))
    allocate(decayBenthos          ( benthos_num ))     ! [1/day] Decay rate of detritus in the benthic layer
    allocate(PAR3D                 ( nl-1, node_size ))


    GloFeDust             = 0.d0
    AtmFeInput            = 0.d0
    GloNDust              = 0.d0
    AtmNInput             = 0.d0

    RiverDIN2D            = 0.d0
    RiverDON2D            = 0.d0
    RiverDOC2D            = 0.d0
    RiverDSi2D            = 0.d0
    RiverDIC2D            = 0.d0
    RiverAlk2D            = 0.d0
    RiverFe               = 0.d0

    ErosionTON2D          = 0.d0
    ErosionTON2D          = 0.d0
    ErosionTSi2D          = 0.d0

    relax_alk             = 0.d0
    virtual_alk           = 0.d0

    cosAI                 = 0.d0
    GloPCO2surf           = 0.d0
    GloCO2flux            = 0.d0
    GloCO2flux_seaicemask = 0.d0
    GloO2flux_seaicemask  = 0.d0
    GlodPCO2surf          = 0.d0
    GlodecayBenthos       = 0.d0
    Benthos               = 0.d0
    Benthos_tr(:,:,:)     = 0.0d0 ! kh 25.03.22
    GloHplus              = exp(-8.d0 * log(10.d0)) ! = 10**(-8)
    DenitBen              = 0.d0
    PistonVelocity        = 0.d0
    alphaCO2              = 0.d0

    LocBenthos            = 0.d0
    decayBenthos          = 0.d0
    PAR3D                 = 0.d0

!    pco2surf           = 0.d0
!    dflux              = 0.d0
!    oflux              = 0.d0
!    co2flux_seaicemask = 0.d0
!    o2flux_seaicemask  = 0.d0
!    dpco2surf          = 0.d0
!    co2                = 0.d0

    if (Diags) then

!! *** Allocate 2D diagnostics ***
    allocate(NPPn    ( node_size ))
    allocate(NPPd    ( node_size ))
    allocate(NPPc    ( node_size ))
    allocate(NPPp    ( node_size ))
    allocate(GPPn    ( node_size ))
    allocate(GPPd    ( node_size ))
    allocate(GPPc    ( node_size ))
    allocate(GPPp    ( node_size ))
    allocate(NNAn    ( node_size ))
    allocate(NNAd    ( node_size ))
    allocate(NNAc    ( node_size ))
    allocate(NNAp    ( node_size ))
    allocate(Chldegn ( node_size ))
    allocate(Chldegd ( node_size ))
    allocate(Chldegc ( node_size ))
    allocate(Chldegp ( node_size ))

    NPPn    = 0.d0
    NPPd    = 0.d0
    NPPc    = 0.d0
    NPPp    = 0.d0
    GPPn    = 0.d0
    GPPd    = 0.d0
    GPPc    = 0.d0
    GPPp    = 0.d0
    NNAn    = 0.d0
    NNAd    = 0.d0
    NNAc    = 0.d0
    NNAp    = 0.d0
    Chldegn = 0.d0
    Chldegd = 0.d0
    Chldegc = 0.d0
    Chldegp = 0.d0

  allocate(grazmeso_tot(node_size))
  allocate(grazmeso_n(node_size))
  allocate(grazmeso_d(node_size))
  allocate(grazmeso_c(node_size))
  allocate(grazmeso_p(node_size))
  allocate(grazmeso_det(node_size))
  allocate(grazmeso_mic(node_size))
  allocate(grazmeso_det2(node_size))

  grazmeso_tot = 0.d0
  grazmeso_n   = 0.d0
  grazmeso_d   = 0.d0
  grazmeso_c   = 0.d0
  grazmeso_p   = 0.d0
  grazmeso_det = 0.d0
  grazmeso_mic = 0.d0
  grazmeso_det2= 0.d0

  allocate(grazmacro_tot(node_size))
  allocate(grazmacro_n(node_size))
  allocate(grazmacro_d(node_size))
  allocate(grazmacro_c(node_size))
  allocate(grazmacro_p(node_size))
  allocate(grazmacro_mes(node_size))
  allocate(grazmacro_det(node_size))
  allocate(grazmacro_mic(node_size))
  allocate(grazmacro_det2(node_size))

  grazmacro_tot = 0.d0
  grazmacro_n = 0.d0
  grazmacro_d = 0.d0
  grazmacro_c = 0.d0
  grazmacro_p = 0.d0
  grazmacro_mes = 0.d0
  grazmacro_det = 0.d0
  grazmacro_mic = 0.d0
  grazmacro_det2= 0.d0

  allocate(grazmicro_tot(node_size))
  allocate(grazmicro_n(node_size))
  allocate(grazmicro_d(node_size))
  allocate(grazmicro_c(node_size))
  allocate(grazmicro_p(node_size))

  grazmicro_tot = 0.d0
  grazmicro_n = 0.d0
  grazmicro_d = 0.d0
  grazmicro_c = 0.d0
  grazmicro_p = 0.d0

!! *** Allocate 3D diagnostics ***
    allocate(respmeso     ( nl-1, node_size ))
    allocate(respmacro    ( nl-1, node_size ))
    allocate(respmicro    ( nl-1, node_size ))
    allocate(calcdiss     ( nl-1, node_size ))
    allocate(calcif       ( nl-1, node_size ))
    allocate(aggn         ( nl-1, node_size ))
    allocate(aggd         ( nl-1, node_size ))
    allocate(aggc         ( nl-1, node_size ))
    allocate(aggp         ( nl-1, node_size ))
    allocate(docexn       ( nl-1, node_size ))
    allocate(docexd       ( nl-1, node_size ))
    allocate(docexc       ( nl-1, node_size ))
    allocate(docexp       ( nl-1, node_size ))
    allocate(respn        ( nl-1, node_size ))
    allocate(respd        ( nl-1, node_size ))
    allocate(respc        ( nl-1, node_size ))
    allocate(respp        ( nl-1, node_size ))
    allocate(NPPn3D       ( nl-1, node_size ))
    allocate(NPPd3D       ( nl-1, node_size ))
    allocate(NPPc3D       ( nl-1, node_size ))
    allocate(NPPp3D       ( nl-1, node_size ))

    respmeso     = 0.d0
    respmacro    = 0.d0
    respmicro    = 0.d0
    calcdiss     = 0.d0
    calcif       = 0.d0
    aggn         = 0.d0
    aggd         = 0.d0
    aggc         = 0.d0
    aggp         = 0.d0
    docexn       = 0.d0
    docexd       = 0.d0
    docexc       = 0.d0
    docexp       = 0.d0
    respn        = 0.d0
    respd        = 0.d0
    respc        = 0.d0
    respp        = 0.d0
    NPPn3D       = 0.d0
    NPPd3D       = 0.d0
    NPPc3D       = 0.d0
    NPPp3D       = 0.d0

!! From Hannahs new temperature function (not sure if needed as diagnostic):

  allocate(TTemp_diatoms  (nl-1,node_size))
  allocate(TTemp_phyto    (nl-1,node_size))
  allocate(TTemp_cocco    (nl-1,node_size))
  allocate(TTemp_phaeo    (nl-1,node_size))

  TTemp_diatoms  (:,:) = 0.d0
  TTemp_phyto    (:,:) = 0.d0
  TTemp_cocco    (:,:) = 0.d0
  TTemp_phaeo    (:,:) = 0.d0

  allocate(TPhyCO2        (nl-1,node_size))
  allocate(TDiaCO2        (nl-1,node_size))
  allocate(TCoccoCO2      (nl-1,node_size))
  allocate(TPhaeoCO2      (nl-1,node_size))

  TPhyCO2        (:,:) = 0.d0
  TDiaCO2        (:,:) = 0.d0
  TCoccoCO2      (:,:) = 0.d0
  TPhaeoCO2      (:,:) = 0.d0

  allocate(TqlimitFac_phyto     (nl-1,node_size))
  allocate(TqlimitFac_diatoms   (nl-1,node_size))
  allocate(TqlimitFac_cocco     (nl-1,node_size))
  allocate(TqlimitFac_phaeo     (nl-1,node_size))

  TqlimitFac_phyto      (:,:) = 0.d0
  TqlimitFac_diatoms    (:,:) = 0.d0
  TqlimitFac_cocco      (:,:) = 0.d0
  TqlimitFac_phaeo      (:,:) = 0.d0


  allocate(TCphotLigLim_diatoms    (nl-1,node_size))
  allocate(TCphotLigLim_phyto      (nl-1,node_size))
  allocate(TCphotLigLim_cocco      (nl-1,node_size))
  allocate(TCphotLigLim_phaeo      (nl-1,node_size))


  TCphotLigLim_diatoms  (:,:) = 0.d0
  TCphotLigLim_phyto    (:,:) = 0.d0
  TCphotLigLim_cocco    (:,:) = 0.d0
  TCphotLigLim_phaeo    (:,:) = 0.d0

  allocate(TCphot_diatoms       (nl-1,node_size))
  allocate(TCphot_phyto         (nl-1,node_size))
  allocate(TCphot_cocco         (nl-1,node_size))
  allocate(TCphot_phaeo         (nl-1,node_size))

  TCphot_diatoms        (:,:) = 0.d0
  TCphot_phyto          (:,:) = 0.d0
  TCphot_cocco          (:,:) = 0.d0
  TCphot_phaeo          (:,:) = 0.d0

  allocate(TSi_assimDia         (nl-1,node_size))

  TSi_assimDia          (:,:) = 0.d0

    end if

!! *** Allocate 3D mocsy ***
    allocate(CO23D        ( nl-1, node_size ))
    allocate(pH3D         ( nl-1, node_size ))
    allocate(pCO23D       ( nl-1, node_size ))
    allocate(HCO33D       ( nl-1, node_size ))
    allocate(CO33D        ( nl-1, node_size ))
    allocate(OmegaC3D     ( nl-1, node_size ))
    allocate(kspc3D       ( nl-1, node_size ))
    allocate(rhoSW3D      ( nl-1, node_size ))
  
    CO23D(:,:)          = 0.d0
    pH3D(:,:)           = 0.d0
    pCO23D(:,:)         = 0.d0
    HCO33D(:,:)         = 0.d0
    CO33D(:,:)          = 0.d0
    OmegaC3D(:,:)       = 0.d0
    kspc3D(:,:)         = 0.d0
    rhoSW3D(:,:)        = 0.d0

!! *** Allocate ballasting ***
    allocate(rho_particle1       ( nl-1, node_size ))
    allocate(rho_particle2       ( nl-1, node_size ))
    allocate(scaling_density1_3D ( nl,   node_size ))
    allocate(scaling_density2_3D ( nl,   node_size ))
    allocate(scaling_visc_3D     ( nl,   node_size ))
    allocate(seawater_visc_3D    ( nl-1, node_size ))
    rho_particle1       = 0.d0
    rho_particle2       = 0.d0
    scaling_density1_3D = 0.d0
    scaling_density2_3D = 0.d0
    scaling_visc_3D     = 0.d0
    seawater_visc_3D    = 0.d0

    allocate(Sinkingvel1(nl,node_size), Sinkingvel2(nl,node_size))
    Sinkingvel1(:,:)      = 0.d0
    Sinkingvel2(:,:)      = 0.d0

    allocate(Sinkvel1_tr(nl,node_size,num_tracers), Sinkvel2_tr(nl,node_size,num_tracers))  ! OG 16.03.23
    Sinkvel1_tr(:,:,:)    = 0.0d0
    Sinkvel2_tr(:,:,:)    = 0.0d0

    if (use_MEDUSA) then
        allocate(GloSed(node_size,sedflx_num))
        allocate(SinkFlx(node_size,bottflx_num))
        allocate(SinkFlx_tr(node_size,bottflx_num,num_tracers)) ! kh 25.03.22 buffer sums per tracer index

        SinkFlx(:,:)      = 0.d0
        SinkFlx_tr(:,:,:) = 0.0d0 ! kh 25.03.22
        GloSed(:,:)       = 0.d0
        allocate(lb_flux(node_size,9))
        lb_flux(:,:)      = 0.d0
    end if

    ! After reading parecomsetup namelist
    call initialize_tracer_indices

    ! Validation check here
    call validate_recom_tracers(num_tracers, mype)

    ! ... populate tracer_id_array from namelist ...
    tracer_id_array = tracers%data(1:tracers%num_tracers)%ID
    call validate_tracer_id_sequence(tracer_id_array, num_tracers, mype)

    !===============================================================================
    ! Model Configuration Summary
    !===============================================================================
    ! Configuration 1: Base model (enable_3zoo2det=F, enable_coccos=F)
    !   - 2 Phytoplankton: General Phy, Diatoms
    !   - 1 Zooplankton: Heterotrophs
    !   - 1 Detritus pool
    !
    ! Configuration 2: 3Zoo2Det (enable_3zoo2det=T, enable_coccos=F)
    !   - 2 Phytoplankton: General Phy, Diatoms
    !   - 3 Zooplankton: Het, Zoo2, Zoo3
    !   - 2 Detritus pools
    !
    ! Configuration 3: Coccos (enable_3zoo2det=F, enable_coccos=T)
    !   - 4 Phytoplankton: General Phy, Diatoms, Coccos, Phaeo
    !   - 1 Zooplankton: Heterotrophs
    !   - 1 Detritus pool
    !
    ! Configuration 4: Full model (enable_3zoo2det=T, enable_coccos=T)
    !   - 4 Phytoplankton: General Phy, Diatoms, Coccos, Phaeo
    !   - 3 Zooplankton: Het, Zoo2, Zoo3
    !   - 2 Detritus pools
    !===============================================================================

    DO i=num_tracers-bgc_num+1, num_tracers
        id=tracers%data(i)%ID

        SELECT CASE (id)

        !---------------------------------------------------------------------------
        ! Base Model: 2 Phytoplankton + 1 Zooplankton + 1 Detritus
        !---------------------------------------------------------------------------
        ! Skip: DIN, DIC, Alk, DSi and O2 are read from files
        ! Fe [mol/L] => [umol/m3] Check the units again!

        ! --- Small Phytoplankton
        CASE (1004)  ! PhyN - Phytoplankton Nitrogen
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max

        CASE (1005)  ! PhyC - Phytoplankton Carbon
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max/NCmax

        CASE (1006)  ! PhyChl - Phytoplankton Chlorophyll
            tracers%data(i)%values(:,:) = tiny_chl

        ! --- Detritus (Non-living organic matter) ---
        CASE (1007)  ! DetN - Detrital Nitrogen
            tracers%data(i)%values(:,:) = tiny

        CASE (1008)  ! DetC - Detrital Carbon
            tracers%data(i)%values(:,:) = tiny

        ! --- Mesozooplankton (Heterotrophs) ---
        CASE (1009)  ! HetN - Heterotroph Nitrogen
            tracers%data(i)%values(:,:) = tiny

        CASE (1010)  ! HetC - Heterotroph Carbon (using Redfield ratio)
            tracers%data(i)%values(:,:) = tiny * Redfield

        ! --- Dissolved Organic Matter ---
        CASE (1011)  ! DON - Dissolved Organic Nitrogen
            tracers%data(i)%values(:,:) = tiny

        CASE (1012)  ! DOC - Dissolved Organic Carbon
            tracers%data(i)%values(:,:) = tiny

        ! --- Diatoms ---
        CASE (1013)  ! DiaN - Diatom Nitrogen
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max

        CASE (1014)  ! DiaC - Diatom Carbon
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max/NCmax

        CASE (1015)  ! DiaChl - Diatom Chlorophyll
            tracers%data(i)%values(:,:) = tiny_chl

        CASE (1016)  ! DiaSi - Diatom Silica
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max_d/NCmax_d/SiCmax

        CASE (1017)  ! DetSi - Detrital Silica
            tracers%data(i)%values(:,:) = tiny

        ! --- Iron (micronutrient) ---
        CASE (1019)  ! Fe - Iron (unit conversion: mol/L => umol/m3)
            tracers%data(i)%values(:,:) = tracers%data(i)%values(:,:)* 1.e9

        ! --- Calcium Carbonate (Calcite) ---
        CASE (1020)  ! PhyCalc - Phytoplankton Calcite
            tracers%data(i)%values(:,:) = tiny * Redfield

        CASE (1021)  ! DetCalc - Detrital Calcite
            tracers%data(i)%values(:,:) = tiny

        !---------------------------------------------------------------------------
        ! Extended Model: Additional Zooplankton and Detritus (enable_3zoo2det)
        !---------------------------------------------------------------------------

        CASE (1023)
            IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                ! Zoo2N - Macrozooplankton Nitrogen
                tracers%data(i)%values(:,:) = tiny
            ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                ! CoccoN - Coccolithophore Nitrogen
                tracers%data(i)%values(:,:) = tiny_chl / chl2N_max
            END IF

        CASE (1024)
            IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                ! Zoo2C - Macrozooplankton Carbon
                tracers%data(i)%values(:,:) = tiny * Redfield
            ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                ! CoccoC - Coccolithophore Carbon
                tracers%data(i)%values(:,:) = tiny_chl / chl2N_max / NCmax
            END IF

        CASE (1025)
            IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                ! DetZ2N - Macrozooplankton Detrital Nitrogen
                tracers%data(i)%values(:,:) = tiny
            ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                ! CoccoChl - Coccolithophore Chlorophyll
                tracers%data(i)%values(:,:) = tiny_chl
            END IF

        CASE (1026)
            IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                ! DetZ2C - Macrozooplankton Detrital Carbon
                tracers%data(i)%values(:,:) = tiny
            ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                ! PhaeoN - Phaeocystis Nitrogen
                tracers%data(i)%values(:,:) = tiny_chl / chl2N_max
            END IF

        CASE (1027)
            IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                ! DetZ2Si - Zooplankton 2 Detrital Silica
                tracers%data(i)%values(:,:) = tiny
            ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                ! PhaeoC - Phaeocystis Carbon
                tracers%data(i)%values(:,:) = tiny_chl / chl2N_max / NCmax
            END IF

        CASE (1028)
            IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                ! DetZ2Calc - Macrozooplankton Detrital Calcite
                tracers%data(i)%values(:,:) = tiny
            ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                ! PhaeoChl - Phaeocystis Chlorophyll
                tracers%data(i)%values(:,:) = tiny_chl
            END IF

        !---------------------------------------------------------------------------
        ! Extended Model: Coccolithophores with 3Zoo2Det
        !---------------------------------------------------------------------------

        CASE (1029)
            IF (enable_coccos .AND. enable_3zoo2det) THEN
                ! CoccoN - Coccolithophore Nitrogen
                tracers%data(i)%values(:,:) = tiny_chl / chl2N_max
            ELSE IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                ! Zoo3N - Microzooplankton Nitrogen
                tracers%data(i)%values(:,:) = tiny
            END IF

        CASE (1030)
            IF (enable_coccos .AND. enable_3zoo2det) THEN
                ! CoccoC - Coccolithophore Carbon
                tracers%data(i)%values(:,:) = tiny_chl / chl2N_max / NCmax
            ELSE IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                ! Zoo3C - Microzooplankton Carbon
                tracers%data(i)%values(:,:) = tiny * Redfield
            END IF

        CASE (1031)
            IF (enable_coccos .AND. enable_3zoo2det) THEN
                ! CoccoChl - Coccolithophore Chlorophyll
                tracers%data(i)%values(:,:) = tiny_chl
            END IF

        CASE (1032)
            IF (enable_coccos .AND. enable_3zoo2det) THEN
                ! PhaeoN - Phaeocystis Nitrogen
                tracers%data(i)%values(:,:) = tiny_chl / chl2N_max
            END IF

        CASE (1033)
            IF (enable_coccos .AND. enable_3zoo2det) THEN
                ! PhaeoC - Phaeocystis Carbon
                tracers%data(i)%values(:,:) = tiny_chl / chl2N_max / NCmax
            END IF

        CASE (1034)
            IF (enable_coccos .AND. enable_3zoo2det) THEN
                ! PhaeoChl - Phaeocystis Chlorophyll
                tracers%data(i)%values(:,:) = tiny_chl
            END IF

        CASE (1035)
            IF (enable_coccos .AND. enable_3zoo2det) THEN
                ! Zoo3N - Zooplankton 3 Nitrogen
                tracers%data(i)%values(:,:) = tiny
            END IF

        CASE (1036)
            IF (enable_coccos .AND. enable_3zoo2det) THEN
                ! Zoo3C - Zooplankton 3 Carbon
                tracers%data(i)%values(:,:) = tiny * Redfield
            END IF

        END SELECT
    END DO
!------------------------------------------

    !< Mask hydrothermal vent in Eastern Equatorial Pacific GO
    do row=1, myDim_nod2D+eDim_nod2D
        !if (ulevels_nod2D(row)>1) cycle 
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)-1
        do k=nzmin, nzmax
            ! do not take regions shallower than 2000 m into account
            if (((geo_coord_nod2D(2,row) > -12.5*rad) .and. (geo_coord_nod2D(2,row) < 9.5*rad))&
                .and.((geo_coord_nod2D(1,row)> -106.0*rad) .and. (geo_coord_nod2D(1,row) < -72.0*rad))) then
                if (abs(Z_3d_n(k,row))<2000.0_WP) cycle
                tracers%data(21)%values(k,row) = min(0.3, tracers%data(21)%values(k,row)) ! OG todo: try 0.6 
            end if
        end do
    end do

    !< Mask negative values
    tracers%data(21)%values(:,:) = max(tiny, tracers%data(21)%values(:,:))
!------------------------------------------

    if(mype==0) write(*,*),'Tracers have been initialized as spinup from WOA/glodap netcdf files'
        locDINmax = -66666
        locDINmin = 66666
        locDICmax = locDINmax
        locDICmin = locDINmin
        locAlkmax = locDINmax
        locAlkmin = locDINmin
        locDSimax = locDINmax
        locDSimin = locDINmin
        locDFemax = locDINmax
        locDFemin = locDINmin
        locO2max  = locDINmax
        locO2min  = locDINmin
        locDICremax = locDICremax ! init DIC remin (added by Sina)
        locDICremin = locDICremin ! init DIC remin (added by Sina) 
        do n=1, myDim_nod2d
            locDINmax = max(locDINmax,maxval(tracers%data(3)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDINmin = min(locDINmin,minval(tracers%data(3)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDICmax = max(locDICmax,maxval(tracers%data(4)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDICmin = min(locDICmin,minval(tracers%data(4)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locAlkmax = max(locAlkmax,maxval(tracers%data(5)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locAlkmin = min(locAlkmin,minval(tracers%data(5)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDSimax = max(locDSimax,maxval(tracers%data(20)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDSimin = min(locDSimin,minval(tracers%data(20)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDFemax = max(locDFemax,maxval(tracers%data(21)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDFemin = min(locDFemin,minval(tracers%data(21)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locO2max  = max(locO2max,maxval(tracers%data(24)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locO2min  = min(locO2min,minval(tracers%data(24)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDICremax  = min(locDICremax,minval(tracers%data(37)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) ) ! init DIC remin (added by Sina)
            locDICremin  = min(locDICremin,minval(tracers%data(37)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) ) ! init DIC remin (added by Sina)
        end do

        if (mype==0) write(*,*) "Sanity check for REcoM variables after recom_init call"
        call MPI_AllREDUCE(locDINmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DIN. =', glo
        call MPI_AllREDUCE(locDINmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DIN. =', glo

        call MPI_AllREDUCE(locDICmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DIC. =', glo
        call MPI_AllREDUCE(locDICmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DIC. =', glo
        call MPI_AllREDUCE(locAlkmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. Alk. =', glo
        call MPI_AllREDUCE(locAlkmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. Alk. =', glo
        call MPI_AllREDUCE(locDSimax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DSi. =', glo
        call MPI_AllREDUCE(locDSimin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DSi. =', glo
        call MPI_AllREDUCE(locDFemax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DFe. =', glo
        call MPI_AllREDUCE(locDFemin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  `-> gobal min init. DFe. =', glo
        call MPI_AllREDUCE(locO2max , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. O2. =', glo
        call MPI_AllREDUCE(locO2min , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  `-> gobal min init. O2. =', glo
        call MPI_AllREDUCE(locDICremax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DICremin =', glo
        call MPI_AllREDUCE(locDICremin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DICremin =', glo

        if (enable_3zoo2det) then
            is_3zoo2det=1.0_WP
        else
            is_3zoo2det=0.0_WP
        endif

        if (enable_coccos) then
            is_coccos=1.0_WP
        else
            is_coccos=0.0_WP
        endif
    end subroutine recom_init

