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
    ! fesom modules
!    use o_ARRAYS
!    use o_MESH
    implicit none
#include "netcdf.inc"
    !___________________________________________________________________________
    ! pointer on necessary derived types 
    integer                                 :: n, k, row, nzmin, nzmax
    integer                                 :: elem_size, node_size, num_tracers

    real(kind=WP)                           :: locDINmax, locDINmin, locDICmax, locDICmin, locAlkmax, glo
    real(kind=WP)                           :: locAlkmin, locDSimax, locDSimin, locDFemax, locDFemin

    type(t_tracer), intent(inout), target   :: tracers
    type(t_partit), intent(inout), target   :: partit
    type(t_mesh),   intent(in) ,   target   :: mesh
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D+eDim_nod2D
    num_tracers=tracers%num_tracers

    if (.not. use_REcoM) return
    
    ! *** allocate and initialize ***

    allocate(GloFeDust             ( node_size ))
    allocate(AtmFeInput            ( node_size ))
    allocate(GloNDust              ( node_size ))
    allocate(AtmNInput             ( node_size ))

    allocate(cosAI                 ( node_size ))
    allocate(GloPCO2surf           ( node_size ))
    allocate(GloCO2flux            ( node_size ))
    allocate(GloCO2flux_seaicemask ( node_size ))
    allocate(GloO2flux_seaicemask  ( node_size ))
    allocate(GlodPCO2surf          ( node_size ))
    allocate(GlodecayBenthos       ( node_size, benthos_num ))
    allocate(Benthos               ( node_size, benthos_num ))
!    allocate(Benthos_tr            ( node_size, benthos_num, num_tracers )) ! kh 25.03.22 buffer per tracer index
    allocate(GloHplus              ( node_size ))
    allocate(DenitBen              ( node_size )) 
    allocate(RiverDIN2D            ( node_size ))
    allocate(RiverDON2D            ( node_size ))
    allocate(RiverDOC2D            ( node_size ))
    allocate(RiverDSi2D            ( node_size ))
    allocate(RiverDIC2D            ( node_size ))
    allocate(RiverAlk2D            ( node_size ))
    allocate(ErosionTON2D          ( node_size ))
    allocate(ErosionTOC2D          ( node_size ))
    allocate(ErosionTSi2D          ( node_size ))
    allocate(LocBenthos            ( benthos_num ))
    allocate(decayBenthos          ( benthos_num )) ! [1/day] Decay rate of detritus in the benthic layer
    allocate(wFluxPhy              ( benthos_num ))     ! [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of phytoplankton
    allocate(wFluxDia              ( benthos_num ))     ! [mmol/(m2 * day)] Flux of N,C, Si and chl through sinking of diatoms 	
    allocate(PAR3D                 (nl-1, node_size ))

    GloFeDust             = 0.d0
    AtmFeInput            = 0.d0
    GloNDust              = 0.d0
    AtmNInput             = 0.d0
  
    cosAI                 = 0.d0
    GloPCO2surf           = 0.d0
    GloCO2flux            = 0.d0
    GloCO2flux_seaicemask = 0.d0 
    GloO2flux_seaicemask  = 0.d0 
    GlodPCO2surf          = 0.d0
    GlodecayBenthos       = 0.d0
    Benthos               = 0.d0
!    Benthos_tr(:,:,:) = 0.0d0 ! kh 25.03.22
    GloHplus              = exp(-8.d0 * log(10.d0)) ! = 10**(-8)
    DenitBen              = 0.d0

    RiverDIN2D            = 0.d0
    RiverDON2D            = 0.d0
    RiverDOC2D            = 0.d0
    RiverDSi2D            = 0.d0
    RiverDIC2D            = 0.d0
    RiverAlk2D            = 0.d0
    ErosionTON2D          = 0.d0
    ErosionTON2D          = 0.d0
    ErosionTSi2D          = 0.d0

    LocBenthos            = 0.d0
    decayBenthos          = 0.d0
    wFluxPhy              = 0.d0
    wFluxDia              = 0.d0
    PAR3D                 = 0.d0

#if defined (__3Zoo2Det)
    allocate(GlowFluxDet ( node_size, benthos_num+4 ))
    GlowFluxDet=0.0d0

    allocate(wFluxDet    ( benthos_num+4 ))     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
    wFluxDet=0.0d0

#else
    allocate(GlowFluxDet(node_size,benthos_num))
    GlowFluxDet=0.0d0
    allocate(wFluxDet(benthos_num))     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
    wFluxDet=0.0d0

#endif

    allocate(GlowFluxPhy(node_size,benthos_num))
    GlowFluxPhy=0.0d0
    allocate(GlowFluxDia(node_size,benthos_num))
    GlowFluxDia=0.0d0    

    AtmCO2             = 0.d0 
    Hplus              = 0.d0
    pco2surf           = 0.d0
    dflux              = 0.d0
    co2flux_seaicemask = 0.d0
    o2flux_seaicemask  = 0.d0
    dpco2surf          = 0.d0
    co2                = 0.d0

    if (Diags) then
!        allocate(diags2D(node_size,12))   ! NEW  (8 -> 12 added 3rd zoo and ballasting)
!        diags2D(:,:)      = 0.d0

!--- Allocate 2D diagnostics
    allocate(NPPn    ( node_size ))
    allocate(NPPd    ( node_size ))
    allocate(NPPc    ( node_size ))
    allocate(GPPn    ( node_size ))
    allocate(GPPd    ( node_size ))
    allocate(GPPc    ( node_size ))
    allocate(NNAn    ( node_size ))
    allocate(NNAd    ( node_size ))
    allocate(NNAc    ( node_size ))
    allocate(Chldegn ( node_size ))
    allocate(Chldegd ( node_size ))
    allocate(Chldegc ( node_size ))

    NPPn    = 0.d0
    NPPd    = 0.d0
    NPPc    = 0.d0
    GPPn    = 0.d0
    GPPd    = 0.d0
    GPPc    = 0.d0
    NNAn    = 0.d0
    NNAd    = 0.d0
    NNAc    = 0.d0
    Chldegn = 0.d0
    Chldegd = 0.d0
    Chldegc = 0.d0

!--- Allocate 3D diagnostics
    allocate(grazmeso_tot ( nl-1, node_size ))
    allocate(grazmeso_n   ( nl-1, node_size ))
    allocate(grazmeso_d   ( nl-1, node_size ))
    allocate(grazmeso_c   ( nl-1, node_size ))
    allocate(respmeso     ( nl-1, node_size ))
    allocate(respmacro    ( nl-1, node_size ))
    allocate(respmicro    ( nl-1, node_size ))
    allocate(calcdiss     ( nl-1, node_size ))
    allocate(calcif       ( nl-1, node_size ))
    allocate(aggn         ( nl-1, node_size ))
    allocate(aggd         ( nl-1, node_size ))
    allocate(aggc         ( nl-1, node_size ))
    allocate(docexn       ( nl-1, node_size ))
    allocate(docexd       ( nl-1, node_size ))
    allocate(docexc       ( nl-1, node_size ))
    allocate(respn        ( nl-1, node_size ))
    allocate(respd        ( nl-1, node_size ))
    allocate(respc        ( nl-1, node_size ))
    allocate(NPPn3D       ( nl-1, node_size ))
    allocate(NPPd3D       ( nl-1, node_size ))
    allocate(NPPc3D       ( nl-1, node_size ))

    grazmeso_tot = 0.d0
    grazmeso_n   = 0.d0
    grazmeso_d   = 0.d0
    grazmeso_c   = 0.d0
    respmeso     = 0.d0
    respmacro    = 0.d0
    respmicro    = 0.d0
    calcdiss     = 0.d0
    calcif       = 0.d0
    aggn         = 0.d0
    aggd         = 0.d0
    aggc         = 0.d0
    docexn       = 0.d0
    docexd       = 0.d0
    docexc       = 0.d0
    respn        = 0.d0
    respd        = 0.d0
    respc        = 0.d0
    NPPn3D       = 0.d0
    NPPd3D       = 0.d0
    NPPc3D       = 0.d0

    end if  

    allocate(CO23D               ( nl-1, node_size ))
    allocate(pH3D                ( nl-1, node_size ))
    allocate(pCO23D              ( nl-1, node_size ))
    allocate(HCO33D              ( nl-1, node_size ))
    allocate(CO33D               ( nl-1, node_size ))
    allocate(OmegaC3D            ( nl-1, node_size ))
    allocate(kspc3D              ( nl-1, node_size ))
    allocate(rhoSW3D             ( nl-1, node_size ))
    allocate(rho_particle1       ( nl-1, node_size ))
    allocate(rho_particle2       ( nl-1, node_size ))
    allocate(scaling_density1_3D ( nl,   node_size ))
    allocate(scaling_density2_3D ( nl,   node_size ))
    allocate(scaling_visc_3D     ( nl,   node_size ))
    allocate(seawater_visc_3D    ( nl-1, node_size ))

    CO23D               = 0.d0
    pH3D                = 0.d0
    pCO23D              = 0.d0
    HCO33D              = 0.d0
    CO33D               = 0.d0
    OmegaC3D            = 0.d0
    kspc3D              = 0.d0
    rhoSW3D             = 0.d0
    rho_particle1       = 0.d0
    rho_particle2       = 0.d0
    scaling_density1_3D = 0.d0
    scaling_density2_3D = 0.d0
    scaling_visc_3D     = 0.d0
    seawater_visc_3D    = 0.d0

    if (ciso) then
    allocate(Cphot_z                  ( nl-1 ))
    allocate(Cphot_dia_z              ( nl-1 ))
    allocate(GloPCO2surf_13           ( node_size ))
    allocate(GloPCO2surf_14           ( node_size ))
    allocate(GloCO2flux_13            ( node_size ))
    allocate(GloCO2flux_14            ( node_size ))
    allocate(GloCO2flux_seaicemask_13 ( node_size ))
    allocate(GloCO2flux_seaicemask_14 ( node_size ))

    Cphot_z                  = 0.d0
    Cphot_dia_z              = 0.d0
    GloPCO2surf_13           = 0.d0
    GloPCO2surf_14           = 0.d0
    GloCO2flux_13            = 0.d0
    GloCO2flux_14            = 0.d0
    GloCO2flux_seaicemask_13 = 0.0d0
    GloCO2flux_seaicemask_14 = 0.0d0
   end if

!        allocate(Sinkingvel1(nl,node_size), Sinkingvel2(nl,node_size))  ! OG 16.03.23
!        allocate(Sinkvel1_tr(nl,node_size,num_tracers), Sinkvel2_tr(nl,node_size,num_tracers))  ! OG 16.03.23
!    if (use_MEDUSA) then
!        allocate(GloSed(node_size,sedflx_num))
!        allocate(SinkFlx(node_size,bottflx_num))
!        allocate(SinkFlx_tr(node_size,bottflx_num,num_tracers)) ! kh 25.03.22 buffer sums per tracer index

!Sinkingvel1(:,:)      = 0.d0  ! OG 16.03.23
!Sinkingvel2(:,:)      = 0.d0  ! OG 16.03.23
!Sinkvel1_tr(:,:,:)    = 0.0d0 ! OG 16.03.23
!Sinkvel2_tr(:,:,:)    = 0.0d0 ! OG 16.03.23
!        allocate(LocSinkFlx(bottflx_num))
!        allocate(LocSed(sedflx_num))
!        SinkFlx(:,:)      = 0.d0
!        SinkFlx_tr(:,:,:) = 0.0d0 ! kh 25.03.22
!        GloSed(:,:)       = 0.d0
!!        LocSed(:)         = 0.d0
!        allocate(lb_flux(node_size,9))
!        lb_flux(:,:)      = 0.d0
!    end if


!    if (useRivFe) then
!       allocate(RiverFe(node_size))
!       RiverFe(:)   = 0.d0
!    end if

! Atmospheric box model
    if (use_atbox) then
!      if (mype==0 .and. my_fesom_group == 0) print *, "Initializing the atmospheric isoCO2 box model ..." !OG
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
      allocate(GloPCO2surf_13           ( node_size ))
      allocate(GloCO2flux_13            ( node_size ))
      allocate(GloCO2flux_seaicemask_13 ( node_size ))

      GloPCO2surf_13           = 0.d0
      GloCO2flux_13            = 0.d0
      GloCO2flux_seaicemask_13 = 0.0d0
      
! Allocate auxiliary inital delta13C_DIC field
      allocate(delta_dic_13_init (nl-1, nod2D ))

      if (ciso_14) then
! Allocate 14CO2 surface fields
      allocate(GloPCO2surf_14           ( node_size ))
      allocate(GloCO2flux_14            ( node_size ))
      allocate(GloCO2flux_seaicemask_14 ( node_size ))    

      GloPCO2surf_14           = 0.d0
      GloCO2flux_14            = 0.d0
      GloCO2flux_seaicemask_14 = 0.0d0

! Allocate auxiliary inital d|Delta14C_DIC fields
      allocate(delta_dic_14_init     ( nl-1, nod2D ))
      allocate(big_delta_dic_14_init ( nl-1, nod2D ))
      end if    ! ciso_14

    end if      ! ciso

    !___________________________________________________________________________
    if(mype==0) write(*,*) 'Benthic layers are set'

    !tracers%data(3)%values ! DIN
    !tracers%data(4)%values ! DIC
    !tracers%data(5)%values ! Alk

    tracers%data(6)%values  = tiny_chl/chl2N_max       ! PhyN
    tracers%data(7)%values  = tiny_chl/chl2N_max/NCmax ! PhyC
    tracers%data(8)%values  = tiny_chl                 ! PhyChl


    tracers%data(9)%values  = tiny            ! DetN
    tracers%data(10)%values = tiny            ! DetC

    tracers%data(11)%values = tiny            ! HetN
    tracers%data(12)%values = tiny * Redfield ! HetC

    tracers%data(13)%values = tiny            ! DON
    tracers%data(14)%values = tiny            ! DOC

    tracers%data(15)%values = tiny_chl/chl2N_max       ! DiaN
    tracers%data(16)%values = tiny_chl/chl2N_max/NCmax ! DiaC
    tracers%data(17)%values = tiny_chl                 ! DiaChl

    tracers%data(18)%values = tiny_chl/chl2N_max_d/NCmax_d/SiCmax ! DiaSi

    tracers%data(19)%values = tiny ! DetSi 
    !tracers%data(20)%values ! DSi     
    tracers%data(21)%values = tracers%data(21)%values * 1.e9 ! Fe [mol/L] => [umol/m3] Check the units again!

    tracers%data(22)%values = tiny * Redfield ! PhyCalc
    tracers%data(23)%values = tiny            ! DetCalc
    !tracers%data(24)%values ! Oxy

#if defined (__3Zoo2Det)
    tracers%data(25)%values = tiny            ! Zoo2N
    tracers%data(26)%values = tiny * Redfield ! Zoo2C
    tracers%data(27)%values = tiny            ! DetZ2N                              
    tracers%data(28)%values = tiny            ! DetZ2C                                    
    tracers%data(29)%values = tiny            ! DetZ2Si                            
    tracers%data(30)%values = tiny            ! DetZ2Calc 
#endif

#if defined (__coccos)
    tracers%data(31)%values = tiny_chl/chl2N_max       ! CoccoN
    tracers%data(32)%values = tiny_chl/chl2N_max/NCmax ! CoccoC 
    tracers%data(33)%values = tiny_chl                 ! CoccoChl
#endif

#if defined (__3Zoo2Det)
    tracers%data(34)%values = tiny            ! Zoo3N
    tracers%data(35)%values = tiny * Redfield ! Zoo3C
#endif

!  if (ciso) then
!
!!   DIC_13
!    if (ciso_init) then
!!     delta13C_DIC according to GLODAP for depths > 500 m
!      delta_dic_13_init = (2.3 - 0.06 * tr_arr(:,:,3))
!    else
!      delta_dic_13_init = 0.
!    end if
!    tr_arr(:,:,idic_13 + 2) = (1. + 0.001 * delta_dic_13_init) * tr_arr(:,:,4) 
!
!!   POC_13
!    tr_arr(:,:,iphyc_13 + 2)   = tr_arr(:,:,7)
!    tr_arr(:,:,idetc_13 + 2)   = tr_arr(:,:,10)
!    tr_arr(:,:,ihetc_13 + 2)   = tr_arr(:,:,12)
!    tr_arr(:,:,idoc_13  + 2)   = tr_arr(:,:,14)
!    tr_arr(:,:,idiac_13 + 2)   = tr_arr(:,:,16)
!    tr_arr(:,:,iphycal_13 + 2) = tr_arr(:,:,22)
!    tr_arr(:,:,idetcal_13 + 2) = tr_arr(:,:,23)
!
!    if (ciso_14) then
!!     DIC_14
!      if (ciso_init) then
!!       Delta14C_DIC according to Broecker et al. (1995):
!        big_delta_dic_14_init = -70. - tr_arr(:,:,20)
!      else
!!       global-mean value (GLODAP-1 BkgC14 ~ -146) permil
!        big_delta_dic_14_init = -150.
!      end if
!      if (ciso_organic_14) then
!!       Stuiver & Pollach (1977, eq. (2)):
!        delta_dic_14_init = (big_delta_dic_14_init + 2. * (delta_dic_13_init + 25.)) / &
!                            (0.95 - 0.002 * delta_dic_13_init)
!      else
!!       simplified "inorganic" radiocarbon
!        delta_dic_14_init = big_delta_dic_14_init
!      end if
!      tr_arr(1:16,:,idic_14 + 2)    = 0.95 * tr_arr(1:16,:,4)
!      tr_arr(17:nl-1,:,idic_14 + 2) = (1. + 0.001 * delta_dic_14_init) * tr_arr(17:nl-1,:,4)
!
!!     POC_14
!      if (ciso_organic_14) then
!        tr_arr(:,:,iphyc_14 + 2)   = tr_arr(:,:,7)
!        tr_arr(:,:,idetc_14 + 2)   = tr_arr(:,:,10)
!        tr_arr(:,:,ihetc_14 + 2)   = tr_arr(:,:,12)
!        tr_arr(:,:,idoc_14 + 2)    = tr_arr(:,:,14)
!        tr_arr(:,:,idiac_14 +2)    = tr_arr(:,:,16)
!        tr_arr(:,:,iphycal_14 + 2) = tr_arr(:,:,22)
!        tr_arr(:,:,idetcal_14 + 2) = tr_arr(:,:,23)
!      end if
!    end if
!!   ciso_14
!  end if
! ciso

! Mask hydrothermal vent in Eastern Equatorial Pacific GO
     do row=1, myDim_nod2D+eDim_nod2D
     if (ulevels_nod2D(row)>1) cycle 
         nzmin = ulevels_nod2D(row)
         nzmax = nlevels_nod2D(row)-1
         do k=nzmin, nzmax
     ! do not take regions deeper than 2000 m into account ! NOT for now OG 23.03.2022
             if (((geo_coord_nod2D(2,row)>-12.5*rad) .and. (geo_coord_nod2D(2,row)<9.5*rad))&
                  .and.((geo_coord_nod2D(1,row)>-106.0  *rad) .and. (geo_coord_nod2D(1,row)<-72.0*rad))) then
!                  if (Z_3d_n(k,row)<-2000.0_WP) cycle
                  tracers%data(21)%values(k,row) = min(0.3,tracers%data(21)%values(k,row))  ! Fe
             end if
         end do
     end do



! Forcing arrays
    allocate(Alk_surf    ( node_size ))
    allocate(relax_alk   ( node_size ))
    allocate(virtual_alk ( node_size ))
    Alk_surf=0.0_WP
    relax_alk=0.0_WP
    virtual_alk=0.0_WP

    if (restore_alkalinity) then
      if (mype==0)  print *, achar(27)//'[46;1m'//'--> Set surface field for alkalinity restoring'//achar(27)//'[0m' 
      Alk_surf=tracers%data(2+ialk)%values(1,:)
      if(mype==0) write(*,*),'Alkalinity restoring = true. Field read in.' 
    endif

!---------------------------------------------------------------------------------------------------------
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

end subroutine recom_init

