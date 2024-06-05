! CONTENT:
! ------------
!    subroutine recom_init
!
! written by V. Schourup-Kristensen, 
! adapted to fesom2.0 by ozgur gurses, 22.05.2020
!    
!===============================================================================
! allocate & initialise arrays for REcoM
subroutine recom_init(mesh)
    USE MOD_MESH
    use g_clock
    use REcoM_declarations
    use REcoM_GloVar
    use REcoM_locVar
    use recom_config
    use REcoM_ciso
    ! fesom modules
    use o_ARRAYS
    use g_PARSUP
    use o_MESH
    implicit none
#include "netcdf.inc"
 
    integer                            :: status, ncid, dimid_rec, nrec, varid
    integer                            :: tra_varid(bgc_num), benthos_varid(benthos_num), hplus_varid
    integer                            :: din_varid, dsi_varid, dfe_varid, alk_varid, dic_varid
    integer                            :: istart(2), icount(2), size2D, size3D     ! Number of 2D and 3D nodes on current processor
    integer                            :: i, j, k, row, nzmin, nzmax
    character(100)                     :: filename
    character(2)                       :: trind
    character(4)                       :: tr_name
    real(kind=8), allocatable          :: ncdata(:)
    character(100)                     :: CO2filename
    character(20)                      :: CO2vari
    integer                            :: CO2start, CO2count
    integer                            :: astat !MB allocation error flag

    integer     :: elem_size, node_size
    integer     :: n
    type(t_mesh), intent(in) , target :: mesh
#include "../associate_mesh.h"
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D+eDim_nod2D

    if (.not. use_REcoM) return
    
    ! *** allocate and initialize ***

    ! GloFeDust and AtmFeInput: atm dep of iron and arrays for diagnostics.
    ! When felimit is not used, arrays are set to zero.

    allocate(GloFeDust(node_size))
    GloFeDust = 0.d0

    allocate(AtmFeInput(node_size))
    AtmFeInput = 0.d0

    ! GloNDust AtmNInput: atm dep of nitrogen sources and sinks. 
    ! They are allocated even nitrogen sources and sinks are not used

    allocate(GloNDust(node_size))
    GloNDust = 0.d0

    allocate(AtmNInput(node_size))
    AtmNInput = 0.d0
    
    ! cosAI: cos of angle of incidence
    allocate(cosAI(node_size))
    cosAI = 0.d0

    allocate(GloPCO2surf(node_size))
    allocate(GloCO2flux(node_size))
    allocate(GloCO2flux_seaicemask(node_size))
    allocate(GloO2flux_seaicemask(node_size))
    allocate(GlodPCO2surf(node_size))
    allocate(GlodecayBenthos(node_size,benthos_num))
    allocate(PAR3D(nl-1,node_size))
    allocate(DenitBen(node_size)) 
    allocate(Benthos(node_size,benthos_num))
    allocate(Benthos_tr(node_size,benthos_num,num_tracers)) ! kh 25.03.22 buffer per tracer index
    allocate(LocBenthos(benthos_num))
    allocate(decayBenthos(benthos_num)) ! [1/day] Decay rate of detritus in the benthic layer

    allocate(wFluxPhy(benthos_num))     ! [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of phytoplankton
    allocate(wFluxDia(benthos_num))     ! [mmol/(m2 * day)] Flux of N,C, Si and chl through sinking of diatoms 	


#if defined (__3Zoo2Det)
    allocate(GlowFluxDet(node_size,benthos_num+4))
    allocate(wFluxDet(benthos_num+4))     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
#else

    allocate(GlowFluxDet(node_size,benthos_num))
    allocate(wFluxDet(benthos_num))     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
#endif

    allocate(GlowFluxPhy(node_size,benthos_num))
    allocate(GlowFluxDia(node_size,benthos_num))

    GlowFluxDet=0.0d0
    GlowFluxPhy=0.0d0
    GlowFluxDia=0.0d0    

    allocate(GloHplus(node_size))

!  allocate(Benthos(benthos_num,node_size), stat=astat)
!  if (astat /= 0) then
!    print *, "Allocation error: Benthos!"
!    stop
!  end if
!  allocate(LocBenthos(benthos_num), stat=astat)
!  if (astat /= 0) then
!    print *, "Allocation error: LocBenthos!"
!    stop
!  end if
!  allocate(decayBenthos(benthos_num), stat=astat)
!  if (astat /= 0) then
!    print *, "Allocation error: decayBenthos!"
!    stop
!  end if
!  allocate(wFluxDet(benthos_num), stat=astat)
!  if (astat /= 0) then
!    print *, "Allocation error: wFluxDet"
!    stop
!  end if
!  allocate(wFluxPhy(benthos_num), stat=astat)
!  if (astat /= 0) then
!    print *, "Allocation error: wFluxPhy"
!    stop
!  end if
!  allocate(wFluxDia(benthos_num), stat=astat)
!  if (astat /= 0) then
!    print *, "Allocation error: wFluxDia"
!    stop
!  end if

    allocate(RiverDIN2D(node_size))
    allocate(RiverDON2D(node_size))
    allocate(RiverDOC2D(node_size))
    allocate(RiverDSi2D(node_size))
    allocate(RiverDIC2D(node_size))
    allocate(RiverAlk2D(node_size))

    allocate(ErosionTON2D(node_size))
    allocate(ErosionTOC2D(node_size))
    allocate(ErosionTSi2D(node_size))


    ! initialise 2d field of CO2 related diagnostics
    GloPCO2surf = 0.d0
    GloCO2flux = 0.d0
    GloCO2flux_seaicemask = 0.0d0 
    GloO2flux_seaicemask = 0.0d0 
    GlodPCO2surf = 0.d0
    GlodecayBenthos = 0.0d0

    ! initialise 2d field of CO2 related variables 
    ! atmospheric CO2 is read from a file if constant_CO2 is false 
    AtmCO2 = 0.d0 
    Hplus = 0.d0
    pco2surf = 0.d0
    dflux = 0.d0
    co2flux_seaicemask = 0.d0
    o2flux_seaicemask = 0.d0
    dpco2surf= 0.d0
    co2      = 0.d0                                    ! NEW

    if (Diags) then
        allocate(diags2D(node_size,12))   ! NEW  (8 -> 12 added 3rd zoo and ballasting)
        diags2D(:,:)      = 0.d0

!--- Allocate 2D diagnostics
  allocate(NPPn(node_size))
  NPPn = 0.d0
  allocate(NPPd(node_size))
  NPPd = 0.d0
  allocate(GPPn(node_size))
  GPPn = 0.d0
  allocate(GPPd(node_size))
  GPPd = 0.d0
  allocate(NNAn(node_size))
  NNAn = 0.d0
  allocate(NNAd(node_size))
  NNAd = 0.d0
  allocate(Chldegn(node_size))
  Chldegn = 0.d0
  allocate(Chldegd(node_size))
  Chldegd = 0.d0
  allocate(NPPc(node_size))
  NPPc = 0.d0
  allocate(GPPc(node_size))
  GPPc = 0.d0
  allocate(NNAc(node_size))
  NNAc = 0.d0
  allocate(Chldegc(node_size))
  Chldegc = 0.d0
  allocate(grazmeso_tot(node_size))
  grazmeso_tot = 0.d0
  allocate(grazmeso_n(node_size))
  grazmeso_n   = 0.d0
  allocate(grazmeso_d(node_size))
  grazmeso_d   = 0.d0
  allocate(grazmeso_c(node_size))
  grazmeso_c   = 0.d0
  allocate(grazmeso_det(node_size))
  grazmeso_det = 0.d0
  allocate(grazmeso_mic(node_size))
  grazmeso_mic = 0.d0
  allocate(grazmeso_det2(node_size))
  grazmeso_det2= 0.d0
  allocate(grazmacro_tot(node_size))
  grazmacro_tot = 0.d0
  allocate(grazmacro_n(node_size))
  grazmacro_n = 0.d0
  allocate(grazmacro_d(node_size))
  grazmacro_d = 0.d0
  allocate(grazmacro_c(node_size))
  grazmacro_c = 0.d0
  allocate(grazmacro_mes(node_size))
  grazmacro_mes = 0.d0
  allocate(grazmacro_det(node_size))
  grazmacro_det = 0.d0
  allocate(grazmacro_mic(node_size))
  grazmacro_mic = 0.d0
  allocate(grazmacro_det2(node_size))
  grazmacro_det2= 0.d0
  allocate(grazmicro_tot(node_size))
  grazmicro_tot = 0.d0
  allocate(grazmicro_n(node_size))
  grazmicro_n = 0.d0
  allocate(grazmicro_d(node_size))
  grazmicro_d = 0.d0
  allocate(grazmicro_c(node_size))
  grazmicro_c = 0.d0
!--- Allocate 3D diagnostics
!  allocate(grazmeso_tot(nl-1,node_size))  ! Comment Miriam (02/2024): changed grazing output from 3D to 2D diagnostics
!  grazmeso_tot(:,:) = 0.d0
!  allocate(grazmeso_n(nl-1,node_size))
!  grazmeso_n(:,:) = 0.d0
!  allocate(grazmeso_d(nl-1,node_size))
!  grazmeso_d(:,:) = 0.d0
!  allocate(grazmeso_c(nl-1,node_size))
!  grazmeso_c(:,:) = 0.d0
  allocate(respmeso(nl-1,node_size))
  respmeso(:,:) = 0.d0
  allocate(respmacro(nl-1,node_size))
  respmacro(:,:) = 0.d0
  allocate(respmicro(nl-1,node_size))
  respmicro(:,:) = 0.d0
  allocate(calcdiss(nl-1,node_size))
  calcdiss(:,:) = 0.d0
  allocate(calcif(nl-1,node_size))
  calcif(:,:) = 0.d0
  allocate(aggn(nl-1,node_size))
  aggn(:,:) = 0.d0
  allocate(aggd(nl-1,node_size))
  aggd(:,:) = 0.d0
  allocate(aggc(nl-1,node_size))
  aggc(:,:) = 0.d0
  allocate(docexn(nl-1,node_size))
  docexn(:,:) = 0.d0
  allocate(docexd(nl-1,node_size))
  docexd(:,:) = 0.d0
  allocate(docexc(nl-1,node_size))
  docexc(:,:) = 0.d0
  allocate(respn(nl-1,node_size))
  respn(:,:) = 0.d0
  allocate(respd(nl-1,node_size))
  respd(:,:) = 0.d0
  allocate(respc(nl-1,node_size))
  respc(:,:) = 0.d0
  allocate(NPPn3D(nl-1,node_size))
  NPPn3D(:,:) = 0.d0
  allocate(NPPd3D(nl-1,node_size))
  NPPd3D(:,:) = 0.d0
  allocate(NPPc3D(nl-1,node_size))
  NPPc3D(:,:) = 0.d0

    end if  

    allocate(CO23D(nl-1,node_size))                   !NEW MOCSY
    CO23D(:,:)          = 0.d0

    allocate(pH3D(nl-1,node_size))                    !NEW MOCSY 
    pH3D(:,:)           = 0.d0

    allocate(pCO23D(nl-1,node_size))                  !NEW MOCSY 
    pCO23D(:,:)         = 0.d0

    allocate(HCO33D(nl-1,node_size))                  !NEW MOCSY 
    HCO33D(:,:)         = 0.d0

    allocate(CO33D(nl-1,node_size))                   !NEW DISS
    CO33D(:,:)          = 0.d0

    allocate(OmegaC3D(nl-1,node_size))                !NEW DISS
    OmegaC3D(:,:)       = 0.d0

    allocate(kspc3D(nl-1,node_size))                  !NEW DISS
    kspc3D(:,:)         = 0.d0

    allocate(rhoSW3D(nl-1,node_size))                 !NEW DISS
    rhoSW3D(:,:)        = 0.d0

    allocate(rho_particle1(nl-1,node_size))           !NEW BALL
    rho_particle1(:,:)  = 0.d0

    allocate(rho_particle2(nl-1,node_size))           !NEW BALL
    rho_particle2(:,:)  = 0.d0

    allocate(scaling_density1_3D(nl,node_size))     !NEW BALL nl-1 -> nl
    scaling_density1_3D(:,:) = 0.d0
 
    allocate(scaling_density2_3D(nl,node_size))     !NEW BALL nl-1 -> nl
    scaling_density2_3D(:,:) = 0.d0

    allocate(scaling_visc_3D(nl,node_size))         !NEW BALL nl-1 -> nl
    scaling_visc_3D(:,:)     = 0.d0

    allocate(seawater_visc_3D(nl-1,node_size))        !NEW BALL
    seawater_visc_3D(:,:)     = 0.d0

    PAR3D(:,:) = 0.d0
    DenitBen(:) = 0.d0

  if (ciso) then
    allocate(Cphot_z(nl-1))
    allocate(Cphot_dia_z(nl-1))
    allocate(GloPCO2surf_13(node_size))
    allocate(GloPCO2surf_14(node_size))
    allocate(GloCO2flux_13(node_size))
    allocate(GloCO2flux_14(node_size))
    allocate(GloCO2flux_seaicemask_13(node_size))
    allocate(GloCO2flux_seaicemask_14(node_size))

    Cphot_z          = 0.d0
    Cphot_dia_z      = 0.d0
    GloPCO2surf_13   = 0.d0
    GloPCO2surf_14   = 0.d0
    GloCO2flux_13    = 0.d0
    GloCO2flux_14    = 0.d0
    GloCO2flux_seaicemask_13 = 0.0d0
    GloCO2flux_seaicemask_14 = 0.0d0
  end if


    RiverDIN2D = 0.d0
    RiverDON2D = 0.d0
    RiverDOC2D = 0.d0
    RiverDSi2D = 0.d0
    RiverDIC2D = 0.d0
    RiverAlk2D = 0.d0

    ErosionTON2D = 0.d0
    ErosionTON2D = 0.d0
    ErosionTSi2D = 0.d0

        allocate(Sinkingvel1(nl,node_size), Sinkingvel2(nl,node_size))  ! OG 16.03.23
        allocate(Sinkvel1_tr(nl,node_size,num_tracers), Sinkvel2_tr(nl,node_size,num_tracers))  ! OG 16.03.23
    if (use_MEDUSA) then
        allocate(GloSed(node_size,sedflx_num))
        allocate(SinkFlx(node_size,bottflx_num))
        allocate(SinkFlx_tr(node_size,bottflx_num,num_tracers)) ! kh 25.03.22 buffer sums per tracer index

Sinkingvel1(:,:)      = 0.d0  ! OG 16.03.23
Sinkingvel2(:,:)      = 0.d0  ! OG 16.03.23
Sinkvel1_tr(:,:,:)    = 0.0d0 ! OG 16.03.23
Sinkvel2_tr(:,:,:)    = 0.0d0 ! OG 16.03.23
!        allocate(LocSinkFlx(bottflx_num))
!        allocate(LocSed(sedflx_num))
        SinkFlx(:,:)      = 0.d0
        SinkFlx_tr(:,:,:) = 0.0d0 ! kh 25.03.22
        GloSed(:,:)       = 0.d0
!        LocSed(:)         = 0.d0
        allocate(lb_flux(node_size,9))
        lb_flux(:,:)      = 0.d0
    end if


    if (useRivFe) then
       allocate(RiverFe(node_size))
       RiverFe(:)   = 0.d0
    end if

! Atmospheric box model
    if (use_atbox) then
      if (mype==0 .and. my_fesom_group == 0) print *, "Initializing the atmospheric isoCO2 box model ..." !OG
      allocate(x_co2atm(node_size))
      x_co2atm    = CO2_for_spinup
      if (ciso) then
        allocate(x_co2atm_13(node_size))
        r_atm_spinup_13 = 1. + 0.001 * delta_co2_13
        x_co2atm_13 = CO2_for_spinup * r_atm_spinup_13
        if (ciso_14) then
          allocate(x_co2atm_14(node_size))
          allocate(cosmic_14(node_size))
          if (ciso_organic_14) then
            delta_co2_14 = (big_delta_co2_14(1) + 2. * delta_co2_13 + 50.) / (0.95 - 0.002 * delta_co2_13)
          else
            delta_co2_14 = big_delta_co2_14(1)
          end if
          r_atm_spinup_14 = 1. + 0.001 * delta_co2_14
          x_co2atm_14    = CO2_for_spinup * r_atm_spinup_14
!         Conversion of initial cosmogenic 14C production rates (mol / s) to fluxes (atoms / s / cm**2) 
!         Since 14C values are scaled to 12C, we need to include the standard 14C / 12C ratio here:
!         1.176e-12 (Karlen et al., 1964) * 6.0221e23 (Avogadro constant) * 1.e-4 (cm**2 / m**2) 
!         = 7.0820e7 cm**2 / m**2
          production_rate_to_flux_14 = 7.0820e7 / ocean_area
          cosmic_14 = cosmic_14_init / production_rate_to_flux_14
        end if
      end if
    end if  ! use_atbox


    if (ciso) then
! Define ciso variables assigning additional ciso tracer indices
      idic_13    = bgc_base_num + 1
      iphyc_13   = bgc_base_num + 2
      idetc_13   = bgc_base_num + 3
      ihetc_13   = bgc_base_num + 4
      idoc_13    = bgc_base_num + 5
      idiac_13   = bgc_base_num + 6
      iphycal_13 = bgc_base_num + 7
      idetcal_13 = bgc_base_num + 8
      idic_14    = bgc_base_num + 9
      iphyc_14   = bgc_base_num + 10
      idetc_14   = bgc_base_num + 11
      ihetc_14   = bgc_base_num + 12
      idoc_14    = bgc_base_num + 13
      idiac_14   = bgc_base_num + 14
      iphycal_14 = bgc_base_num + 15
      idetcal_14 = bgc_base_num + 16


! Allocate 13CO2 surface fields
      allocate(GloPCO2surf_13(node_size))
      allocate(GloCO2flux_13(node_size))
      allocate(GloCO2flux_seaicemask_13(node_size))

      GloPCO2surf_13   = 0.d0
      GloCO2flux_13    = 0.d0
      GloCO2flux_seaicemask_13 = 0.0d0
      
! Allocate auxiliary inital delta13C_DIC field
      allocate(delta_dic_13_init(nl-1,nod2D))

      if (ciso_14) then
! Allocate 14CO2 surface fields
        allocate(GloPCO2surf_14(node_size))
        allocate(GloCO2flux_14(node_size))
        allocate(GloCO2flux_seaicemask_14(node_size))    

        GloPCO2surf_14   = 0.d0
        GloCO2flux_14    = 0.d0
        GloCO2flux_seaicemask_14 = 0.0d0

! Allocate auxiliary inital d|Delta14C_DIC fields
        allocate(delta_dic_14_init(nl-1,nod2D))
        allocate(big_delta_dic_14_init(nl-1,nod2D))
      end if    ! ciso_14

    end if      ! ciso

    !___________________________________________________________________________
    ! Initialization of benthos
    ! Benthic layer consists of Benthos(1) = N, Benthos(2)=C, Benthos(3)=Si, Benthos(4)=calc
  
    Benthos(:,:) = 0.d0 !tiny
    Benthos_tr(:,:,:) = 0.0d0 ! kh 25.03.22

    GloHplus = exp(-8.d0 * log(10.d0)) ! = 10**(-8)
    !___________________________________________________________________________
    if(mype==0) write(*,*) 'Benthic layers are set'

! change to fesom tracer counter, not hard coded counters:

    !tr_arr(:,:,3)                         ! tracer 3  = DIN
    !tr_arr(:,:,4)                         ! tracer 4  = DIC
    !tr_arr(:,:,5)                         ! tracer 5  = Alk

    tr_arr(:,:,6)  = tiny_chl/chl2N_max       !tiny                  ! tracer 6  = PhyN   -> Intracellular conc of Nitrogen in small phytoplankton
    tr_arr(:,:,7)  = tiny_chl/chl2N_max/NCmax !tiny * Redfield       ! tracer 7  = PhyC   -> Intracellular conc of Carbon in small phytoplankton
    tr_arr(:,:,8)  = tiny_chl                 !tiny * 1.56d0         ! tracer 8  = PhyChl -> Current intracellular ChlA conc

    tr_arr(:,:,9)  = tiny                  ! tracer 9  = DetN
    tr_arr(:,:,10) = tiny                  ! tracer 10 = DetC

    tr_arr(:,:,11) = tiny                  ! tracer 11 = HetN
    tr_arr(:,:,12) = tiny * Redfield       ! tracer 12 = HetC

    tr_arr(:,:,13) = tiny                  ! tracer 13 = DON
    tr_arr(:,:,14) = tiny                  ! tracer 14 = DOC

    tr_arr(:,:,15) = tiny_chl/chl2N_max !tiny                  ! tracer 15 = DiaN
    tr_arr(:,:,16) = tiny_chl/chl2N_max/NCmax !tiny * Redfield       ! tracer 16 = DiaC
    tr_arr(:,:,17) = tiny_chl !tiny * 1.56d0         ! tracer 17 = DiaChl

    tr_arr(:,:,18) = tiny_chl/chl2N_max_d/NCmax_d/SiCmax !tiny                   ! tracer 18 = DiaSi

    tr_arr(:,:,19) = tiny                  ! tracer 19 = DetSi 
    !tr_arr(:,:,20)                        ! tracer 20 = DSi     
    tr_arr(:,:,21) = tr_arr(:,:,21) * 1.e9 ! Fe [mol/L] => [umol/m3] Check the units again!

!    tr_arr(:,:,22) = tiny !cPhyN * 0.25d0        ! tracer 22 = PhyCalc 
    tr_arr(:,:,22) = tiny * Redfield       ! tracer 22 = PhyCalc ! NEW now dependent on Cocco biomass
    tr_arr(:,:,23) = tiny                  ! tracer 23 = DetCalc
    !tr_arr(:,:,24)                        ! tracer 24 = Oxy     ! read from the file

#if defined (__3Zoo2Det)
    tr_arr(:,:,25) = tiny                   ! tracer 25 = Zoo2N
    tr_arr(:,:,26) = tiny * Redfield        ! tracer 26 = Zoo2C
    tr_arr(:,:,27) = tiny                   ! tracer 26 = DetZ2N                              
    tr_arr(:,:,28) = tiny                   ! tracer 27 = DetZ2C                                    
    tr_arr(:,:,29) = tiny                   ! tracer 28 = DetZ2Si                            
    tr_arr(:,:,30) = tiny                   ! tracer 29 = DetZ2Calc 
#endif

#if defined (__coccos) & defined (__3Zoo2Det)
    tr_arr(:,:,31) = tiny_chl/chl2N_max        ! tiny             ! tracer 31 = CoccoN
    tr_arr(:,:,32) = tiny_chl/chl2N_max/NCmax  ! tiny * Redfield  ! tracer 32 = CoccoC 
    tr_arr(:,:,33) = tiny_chl                  ! tiny * 1.56d0    ! tracer 33 = CoccoChl
#elif defined (__coccos) & !defined (__3Zoo2Det)
   if (mype==0 .and. my_fesom_group == 0)  print *, "case_3p1z1d"
    tr_arr(:,:,25) = tiny_chl/chl2N_max        ! tracer 25 = CoccoN
    tr_arr(:,:,26) = tiny_chl/chl2N_max/NCmax  ! tracer 26 = CoccoC 
    tr_arr(:,:,27) = tiny_chl  
#endif

#if defined (__coccos) & defined (__3Zoo2Det)
    tr_arr(:,:,34) = tiny                      ! tracer 34 = Zoo3N
    tr_arr(:,:,35) = tiny * Redfield           ! tracer 35 = Zoo3C
#elif !defined (__coccos) & defined (__3Zoo2Det)
    tr_arr(:,:,31) = tiny                      ! tracer 31 = Zoo3N
    tr_arr(:,:,32) = tiny * Redfield           ! tracer 32 = Zoo3C

#endif

  if (ciso) then

!   DIC_13
    if (ciso_init) then
!     delta13C_DIC according to GLODAP for depths > 500 m
      delta_dic_13_init = (2.3 - 0.06 * tr_arr(:,:,3))
    else
      delta_dic_13_init = 0.
    end if
    tr_arr(:,:,idic_13 + 2) = (1. + 0.001 * delta_dic_13_init) * tr_arr(:,:,4) 

!   POC_13
    tr_arr(:,:,iphyc_13 + 2)   = tr_arr(:,:,7)
    tr_arr(:,:,idetc_13 + 2)   = tr_arr(:,:,10)
    tr_arr(:,:,ihetc_13 + 2)   = tr_arr(:,:,12)
    tr_arr(:,:,idoc_13  + 2)   = tr_arr(:,:,14)
    tr_arr(:,:,idiac_13 + 2)   = tr_arr(:,:,16)
    tr_arr(:,:,iphycal_13 + 2) = tr_arr(:,:,22)
    tr_arr(:,:,idetcal_13 + 2) = tr_arr(:,:,23)

    if (ciso_14) then
!     DIC_14
      if (ciso_init) then
!       Delta14C_DIC according to Broecker et al. (1995):
        big_delta_dic_14_init = -70. - tr_arr(:,:,20)
      else
!       global-mean value (GLODAP-1 BkgC14 ~ -146) permil
        big_delta_dic_14_init = -150.
      end if
      if (ciso_organic_14) then
!       Stuiver & Pollach (1977, eq. (2)):
        delta_dic_14_init = (big_delta_dic_14_init + 2. * (delta_dic_13_init + 25.)) / &
                            (0.95 - 0.002 * delta_dic_13_init)
      else
!       simplified "inorganic" radiocarbon
        delta_dic_14_init = big_delta_dic_14_init
      end if
      tr_arr(1:16,:,idic_14 + 2)    = 0.95 * tr_arr(1:16,:,4)
      tr_arr(17:nl-1,:,idic_14 + 2) = (1. + 0.001 * delta_dic_14_init) * tr_arr(17:nl-1,:,4)

!     POC_14
      if (ciso_organic_14) then
        tr_arr(:,:,iphyc_14 + 2)   = tr_arr(:,:,7)
        tr_arr(:,:,idetc_14 + 2)   = tr_arr(:,:,10)
        tr_arr(:,:,ihetc_14 + 2)   = tr_arr(:,:,12)
        tr_arr(:,:,idoc_14 + 2)    = tr_arr(:,:,14)
        tr_arr(:,:,idiac_14 +2)    = tr_arr(:,:,16)
        tr_arr(:,:,iphycal_14 + 2) = tr_arr(:,:,22)
        tr_arr(:,:,idetcal_14 + 2) = tr_arr(:,:,23)
      end if
    end if
!   ciso_14
  end if
! ciso

! Mask hydrothermal vent in Eastern Equatorial Pacific GO
if (1) then 
         do row=1, myDim_nod2D+eDim_nod2D
         if (ulevels_nod2D(row)>1) cycle 
         nzmin = ulevels_nod2D(row)
         nzmax = nlevels_nod2D(row)-1
            !do k=1, nlevels_nod2D(row)
            do k=nzmin, nzmax
     ! do not take regions deeper than 2000 m into account

               if     (((geo_coord_nod2D(2,row)>-12.5*rad) .and. (geo_coord_nod2D(2,row)<9.5*rad))&
                  .and.((geo_coord_nod2D(1,row)>-106.0  *rad) .and. (geo_coord_nod2D(1,row)<-72.0*rad))) then
!                  if (Z_3d_n(k,row)<-2000.0_WP) cycle
                  tr_arr(k,row,21) = min(0.3,tr_arr(k,row,21))  ! Fe [mol/L] => [umol/m3] Check the units again!
               end if
            end do
         end do
end if
!---------------------------------------------------------------------------------------------------------
if(mype==0) write(*,*),'Tracers have been initialized as spinup from WOA/glodap netcdf files'

end subroutine recom_init

