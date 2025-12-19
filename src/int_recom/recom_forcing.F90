!===============================================================================
! REcoM_Forcing
!===============================================================================
subroutine REcoM_Forcing(zNodes, n, Nn, state, SurfSW, Loc_slp, Temp, Sali, Sali_depth &
            , CO2_watercolumn                                          &
            , pH_watercolumn                                           &
            , pCO2_watercolumn                                         &
            , HCO3_watercolumn                                         &
            , CO3_watercolumn                                          &
            , OmegaC_watercolumn                                       &
            , kspc_watercolumn                                         &
            , rhoSW_watercolumn                                        &
            , PAR, ice, dynamics, tracers, partit, mesh)

    use recom_declarations
    use recom_locvar
    use recom_config
    use recom_glovar
    use gasx
    use recom_ciso
    use g_clock
    use o_PARAM
    use g_rotate_grid
    use g_config
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use mod_tracer
    use MOD_DYN
    use MOD_ICE

    use o_param
    use o_arrays
    use g_forcing_arrays
    use g_comm_auto
    use g_support
#ifdef RECOM_WAVEBANDS    
    use REcoM_spectral
#endif /* RECOM_WAVEBANDS  */    
    implicit none

    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_ice)   , intent(inout), target :: ice
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh

    real(kind=8)                              :: Latr          
    integer                                   :: n, Nn         ! Nn -Total number of nodes
    real(kind=8),dimension(mesh%nl-1)	      :: zNodes	       ! Depth of nodes   zr(1:nzmax) = Z_3d_n(1:nzmax,n)
    real(kind=8),dimension(mesh%nl-1,bgc_num) :: state         	
    real(kind=8)                              :: SurfSW        ! [W/m2] ShortWave radiation at surface
    real(kind=8)                              :: Loc_slp       ! [Pa] sea-level pressure
    real(kind=8),dimension(mesh%nl-1)         :: Temp          ! [degrees C] Ocean temperature
    real(kind=8),dimension(mesh%nl-1)         :: Sali_depth    ! Salinity for the whole water column

    !!---- Watercolumn carbonate chemistry
    real(kind=8),dimension(mesh%nl-1)         :: CO2_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: pH_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: pCO2_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: HCO3_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: CO3_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: OmegaC_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: kspc_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: rhoSW_watercolumn

    real(kind=8),dimension(mesh%nl-1)         :: PAR
!SL following code lines are added for spectral light
!SL could be still rearranged later
#ifdef RECOM_WAVEBANDS
    !!---- Spectral Light related
    integer :: ilam, Nr
    integer :: iday,iyr,imon,isec,lp,wd,mydate(4)
    INTEGER :: idiscEs,jdiscEs,kdiscEs,ldiscEs
    INTEGER :: idiscEu,jdiscEu,kdiscEu,ldiscEu
    real(kind=8)                                  :: solz
    Real(kind=8),dimension(tlam)                  :: PARwup
    Real(kind=8),dimension(tlam)                  :: PARwdn
    Real(kind=8),dimension(tlam,mesh%nl-1)        :: PARw_k ! or (mesh%nl-1,tlam)
    Real(kind=8),dimension(mesh%nl-1)             :: PARl
    Real(kind=8),dimension(tlam)                  :: PARwup_diag
    Real(kind=8)                                  :: PARwup_total
    Real(kind=8),dimension(tlam)                  :: Edwsf
    Real(kind=8),dimension(tlam)                  :: Eswsf
#ifdef RECOM_CALC_REFLEC
    INTEGER :: index
    Real(kind=8),dimension(mesh%nl-1)             :: PARw_kwb
#endif /* RECOM_CALC_REFLEC */
    Real(kind=8),dimension(tlam)                  :: C_phot_nl
    Real(kind=8),dimension(tlam)                  :: C_phot_nl_dia
    Real(kind=8),dimension(tlam)                  :: Ek_nl
    Real(kind=8),dimension(tlam)                  :: Ek_nl_dia
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: a_k
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: acdom_k
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: aphy_chl_k
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: aphy_chl_dia_k
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: apart_k
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: actot
!    Real(kind=8),dimension(2,tlam)                :: aclocal
    Real(kind=8),dimension(4,tlam)                :: aclocal
    Real(kind=8)                                  :: atten
    Real(kind=8)                                  :: discEs
    Real(kind=8)                                  :: discEu
    Real(kind=8),dimension(mesh%nl-1)             :: dz_k
! Biomass    
    Real(kind=8),dimension(mesh%nl-1)             :: part_k
!    Real(kind=8),dimension(2,mesh%nl-1)           :: Phy_k
!    Real(kind=8),dimension(2,mesh%nl-1)           :: phychl_k
    Real(kind=8),dimension(4,mesh%nl-1)           :: Phy_k
    Real(kind=8),dimension(4,mesh%nl-1)           :: phychl_k
    if (RECOM_CDOM) then
    Real(kind=8),dimension(mesh%nl-1)             :: cdom_k
    endif  !/* RECOM_CDOM */
    if (RECOM_CALC_APHYT=.true. .and. RECOM_MARSHALL=.true.) then
    Real(kind=8),dimension(mesh%nl-1)             :: phyD1_k
    Real(kind=8),dimension(mesh%nl-1)             :: diaD1_k
    endif  !/* defined(RECOM_CALC_APHYT) && defined(RECOM_MARSHALL)*/
! for diagnostics
    Real(kind=8),dimension(mesh%nl-1)             :: a_kave
    Real(kind=8),dimension(mesh%nl-1)             :: acdom_kave
    Real(kind=8),dimension(mesh%nl-1)             :: apart_kave
    Real(kind=8),dimension(mesh%nl-1)             :: actot_ave
! iops surface
    Real(kind=8),dimension(tlam)                  :: a_ksur
    Real(kind=8),dimension(tlam)                  :: acdom_ksur
    Real(kind=8),dimension(tlam)                  :: apart_ksur
    Real(kind=8),dimension(tlam)                  :: actot_sur
    if (RECOM_CALC_REFLEC) then
    Real(kind=8),dimension(mesh%nl-1)             :: a_kwb
    Real(kind=8),dimension(mesh%nl-1)             :: acdom_kwb
    Real(kind=8),dimension(mesh%nl-1)             :: apart_kwb
    Real(kind=8),dimension(mesh%nl-1)             :: actot_wb
    endif   !/* RECOM_CALC_REFLEC */
    if (RECOM_RADTRANS) then
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: bt_k
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: bb_k
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: bpart_k
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: bbpart_k
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: bctot
    Real(kind=8),dimension(mesh%nl-1,tlam)        :: bbctot
    Real(kind=8),dimension(tlam,mesh%nl-1)     :: Edz      ! or (tlam, mesh%nl-1)
    Real(kind=8),dimension(tlam,mesh%nl-1)     :: Esz      ! or (tlam, mesh%nl-1)
    Real(kind=8),dimension(tlam,mesh%nl-1)     :: Euz      ! or (tlam, mesh%nl-1)
    Real(kind=8),dimension(tlam,mesh%nl-1)     :: Estop    ! or (tlam, mesh%nl-1)
    Real(kind=8),dimension(tlam,mesh%nl-1)     :: Eutop    ! or (tlam, mesh%nl-1)
    Real(kind=8),dimension(mesh%nl-1)          :: tirrq
    Real(kind=8),dimension(tlam,mesh%nl-1)     :: tirrwq   ! or (tlam, mesh%nl-1)
    Real(kind=8),dimension(tlam,mesh%nl-1)     :: amp1     ! or (tlam, mesh%nl-1)
    Real(kind=8),dimension(tlam,mesh%nl-1)     :: amp2     ! or (tlam, mesh%nl-1)
    Real(kind=8),dimension(mesh%nl-1)             :: rmud
    Real(kind=8)                                  :: rn=1.341d0     !refractive index of seawater
    Real(kind=8)                                  :: sinszaw
    Real(kind=8)                                  :: szaw
    Real(kind=8)                                  :: rmudl
    Real(kind=8)                                  :: rn
!  only for diagnostics
    Real(kind=8),dimension(mesh%nl-1)          :: bt_kave
    Real(kind=8),dimension(mesh%nl-1)          :: bpart_kave
    Real(kind=8),dimension(mesh%nl-1)          :: bctot_ave
    Real(kind=8),dimension(mesh%nl-1)          :: bb_kave
    Real(kind=8),dimension(mesh%nl-1)          :: bbpart_kave
    Real(kind=8),dimension(mesh%nl-1)          :: bbctot_ave
    Real(kind=8),dimension(tlam)               :: Eupwel
    Real(kind=8),dimension(tlam)               :: Reflec
!  iops surface     
    Real(kind=8),dimension(tlam)        :: bt_ksur
    Real(kind=8),dimension(tlam)        :: bpart_ksur
    Real(kind=8),dimension(tlam)        :: bctot_sur
    Real(kind=8),dimension(tlam)        :: bb_ksur
    Real(kind=8),dimension(tlam)        :: bbpart_ksur
    Real(kind=8),dimension(tlam)        :: bbctot_sur
    if (RECOM_CALC_REFLEC) then
    Real(kind=8),dimension(mesh%nl-1)             :: bt_kwb
    Real(kind=8),dimension(mesh%nl-1)             :: bpart_kwb
    Real(kind=8),dimension(mesh%nl-1)             :: bctot_wb
    Real(kind=8),dimension(mesh%nl-1)             :: bb_kwb
    Real(kind=8),dimension(mesh%nl-1)             :: bbpart_kwb
    Real(kind=8),dimension(mesh%nl-1)             :: bbctot_wb
    Real(kind=8),dimension(mesh%nl-1)             :: Edz_wb
    Real(kind=8),dimension(mesh%nl-1)             :: Esz_wb
    Real(kind=8),dimension(mesh%nl-1)             :: Euz_wb
    Real(kind=8),dimension(mesh%nl-1)             :: Estop_wb
    Real(kind=8),dimension(mesh%nl-1)             :: Eutop_wb
    Real(kind=8),dimension(mesh%nl-1)             :: amp1_wb
    Real(kind=8),dimension(mesh%nl-1)             :: amp2_wb
    endif !/* RECOM_CALC_REFLEC */
    endif !/* RADTRANS */
!if (enable_coccos) then
!! remember to introduce and declare new _cocco and _phaeo related
!endif
#endif /* RECOM_WAVEBANDS */    


    !!---- Subroutine Depth

    real(kind=8),dimension(mesh%nl)           :: zF                   ! [m] Depth of fluxes
    real(kind=8),dimension(mesh%nl,6)         :: SinkVel              ! [m/day]
    real(kind=8),dimension(mesh%nl-1)         :: thick                ! [m] Vertical distance between two nodes = Thickness 
    real(kind=8),dimension(mesh%nl-1)         :: recipthick           ! [1/m] reciprocal of thick

    !!---- Subroutine CO2Flux /mocsy
    real(kind=8)                              :: REcoM_DIC(1)         ! [mol/m3] Conc of DIC in the surface water, used to calculate CO2 flux
    real(kind=8)                              :: REcoM_Alk(1)         ! [mol/m3] Conc of Alk in the surface water, used to calculate CO2 flux
    real(kind=8)                              :: REcoM_Si(1)          ! [mol/m3] Conc of Si in the surface water, used to calculate CO2 flux
    real(kind=8)                              :: REcoM_Phos(1)        ! [mol/m3] Conc of Phos in the surface water, used to calculate the CO2 flux
    real(kind=8)                              :: Sali(1)              ! Salinity of current surface layer
    real(kind=8)                              :: Latd(1)              ! latitude in degree
    real(kind=8)                              :: Lond(1)              ! longitude in degree
    real(kind=8)                              :: REcoM_T(1)           ! temperature again, for mocsy minimum defined as -2
    real(kind=8)                              :: REcoM_S(1)           ! temperature again, for mocsy minimum defined as 21
! atm pressure, now read in as forcing!!
    !!---- atm pressure
    real(kind=8)                              :: Patm(1)              ! atmospheric pressure [atm]

    !!---- Subroutine o2flux /mocsy 
    real(kind=8)                              :: ppo(1)               ! atmospheric pressure, divided by 1 atm 
    real(kind=8)                              :: REcoM_O2(1)          ! [mmol/m3] Conc of O2 in the surface water, used to calculate O2 flux

    !!---- Subroutine REcoM_sms
    real(kind=8),dimension(mesh%nl-1,bgc_num) :: sms                  ! matrix that entail changes in tracer concentrations

    !!---- Diagnostics
    integer                                   :: idiags,k
    integer                                   :: tr_num

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"


    tiny_N   = tiny_chl/chl2N_max   ! 0.00001/ 3.15d0   Chl2N_max [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
    tiny_N_d = tiny_chl/chl2N_max_d ! 0.00001/ 4.2d0

    tiny_C   = tiny_N  /NCmax       ! NCmax   = 0.2d0   [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
    tiny_C_d = tiny_N_d/NCmax_d     ! NCmax_d = 0.2d0 

    tiny_Si  = tiny_C_d/SiCmax      ! SiCmax = 0.8d0

if (enable_coccos) then
    tiny_N_c = tiny_chl/chl2N_max_c ! 0.00001/ 3.5d0
    tiny_C_c = tiny_N_c/NCmax_c     ! NCmax_c = 0.15d0

    tiny_N_p = tiny_chl/chl2N_max_p ! 0.00001/ 3.5d0 
    tiny_C_p = tiny_N_p/NCmax_p     ! NCmax_c = 0.15d0
endif

    call Cobeta(partit, mesh)  

    call Depth_calculations(n, Nn,SinkVel,zF,thick,recipthick, partit, mesh)

!------------SPECTRAL LIGHT--------------------------------------------
!======================================================================
#ifdef RECOM_WAVEBANDS
       Nr = mesh%nl-1
       idiscEs = 0
       jdiscEs = 0
       kdiscEs = 0
       ldiscEs = 0
       idiscEu = 0
       jdiscEu = 0
       kdiscEu = 0
       ldiscEu = 0
       discEs = 0.
       discEu = 0.
! ------ GET constant acdom_k -------
     DO k=1,Nr
     if (RECOM_CALC_ACDOM =.false.) then
       do ilam = 1,tlam
        acdom_k(k,ilam) = acdom(ilam)
       enddo
     endif      !/* no RECOM_CALC_ACDOM */
! ------ GET constant aphy_chl_k & aphy_chl_dia_k -------
     if (RECOM_CALC_APHYT =.false.) then
       do ilam = 1,tlam
        aphy_chl_k(k,ilam) = aphy_chl(ilam)
        aphy_chl_dia_k(k,ilam) = aphy_chl_dia(ilam)
       enddo
      ENDDO
      endif     !/* no RECOM_CALC_APHYT */
     if (RECOM_CALC_ACDOM) then
!------------ COMPUTE ACDOM_k  ----------------------------------
        if (RECOM_CDOM) then
        call MONOD_ACDOM(Nr,cdom_k(1:Nr),                       &
                         acdom_k(1:Nr,1:tlam),                  &
                         myThid)
        else
        aclocal(1,1:tlam)=aphy_chl
        aclocal(2,1:tlam)=aphy_chl_dia
        aclocal(3,1:tlam)=aphy_chl_cocco
        aclocal(4,1:tlam)=aphy_chl_phaeo
         call MONOD_ACDOM(Nr,phychl_k(4,1:Nr), aclocal, aw,     &
                          acdom_k(1:Nr,1:tlam),                 &
                          myThid)
        endif   !/* RECOM_CDOM */
     endif      !/* RECOM_CALC_ACDOM */
! ------------ COMPUTE aphy_chl_k & aphy_chl_dia_k ----------------
     if (RECOM_CALC_APHYT =.true. .and. RECOM_MARSHALL=.true.) then
        call RECOM_APHYTO(Nr,phyD1_k(1:Nr),QYmax,Drel,          & 
                          aphy_chl_ps(1:tlam),aphy_chl(1:tlam), &
                          aphy_chl_k(1:Nr, 1:tlam),             & 
                          myThid)
        call RECOM_APHYTO(Nr,diaD1_k(1:Nr),QYmax_d,Drel,                &
                          aphy_chl_ps_dia(1:tlam),aphy_chl_dia(1:tlam), &
                          aphy_chl_dia_k(1:Nr, 1:tlam),                 &
                          myThid)
if (enable_coccos) then
        call RECOM_APHYTO(Nr, coccoD1_k(1:Nr),QYmax,Drel,               &
                   aphy_chl_ps_cocco(1:tlam),aphy_chl_cocco(1:tlam),    &
                   aphy_chl_k_cocco(1:Nr, 1:tlam),                      &
                          myThid)
        call RECOM_APHYTO(Nr, phaeoD1_k(1:Nr),QYmax_d,Drel,             &
                   aphy_chl_ps_phaeo(1:tlam),aphy_chl_phaeo(1:tlam),    &
                   aphy_chl_phaeo_k(1:Mr, 1:tlam),                      &
                          myThid)
endif
     endif      !/* RECOM_CALC_APHYT */
! ------------ GET PART_k FOR WAVEBANDS_3D and RADTRANS ----------------
!     In Darwin-MONOD particulate matter is calculated in P units.
!     We use either detC and a C per particle factor (Stramski 2001) or
!     detC with a biomass-specific spectrum (Gallegos 2011).

      DO ilam = 1,tlam
       DO k=1,Nr

      if (RECOM_CALC_APART) then
      apart_k(k,ilam) = part_k(k) * aparcoeff * exapar(ilam)
         if (RECOM_RADTRANS) then
         bpart_k(k,ilam) = part_k(k) * bparcoeff * exbpar(ilam)
         bbpart_k(k,ilam) = part_k(k) * bparcoeff                         &
                       * exbpar(ilam) * bb_to_b
         endif
      else
      apart_k(k,ilam) = part_k(k)*apart_P(ilam)
         if (RECOM_RADTRANS)
         bpart_k(k,ilam) = part_k(k)*bpart_P(ilam)
         bbpart_k(k,ilam) = part_k(k)*bbpart_P(ilam)
         endif
      endif     !/* RECOM_CALC_APART*/
        ENDDO   !k
       ENDDO    !ilam
!C ------------- GET SPECTRAL LIGHT -----------------
       do ilam = 1,tlam
          if (OASIM)
! add direct and diffuse, convert to uEin/m2/s/nm
          Edwsf(ilam) = oasim_ed(ilam) ! * 1000
          Eswsf(ilam) = oasim_es(ilam) ! * 1000
          PARwup(ilam) = WtouEins(ilam) * (Edwsf(ilam)                   &
                               + Eswsf(ilam))
          else
! sf is per nm; convert to per waveband
          Edwsf(ilam) = wb_width(ilam) * sf(ilam)                        &
                      * PARadiation
          Eswsf(ilam) = tiny
          PARwup(ilam) = Edwsf(ilam) * WtouEins(ilam)
          endif
          PARwup_diag(ilam) = PARwup(ilam)
       enddo   ! ilam
          PARwup_total = 0.
             do ilam = 1,tlam
             PARwup_total = PARwup_total + PARwup(ilam)
             enddo ! ilam
if (RECOM_RADTRANS) then
!     Compute 1/cos(zenith) for direct light below surface given solar zenith
!     angle (in radians) at surface (solz=zenith_deg from RECOM_INSOLATION.F)
         sinszaw = sin(solz)/rn
         szaw = asin(sinszaw)
         rmudl = 1.0/cos(szaw)    !avg cosine direct (1 over)
         rmud = min(rmudl,1.5)
         rmud = max(rmud,0.0)
end if   
if (RECOM_RADTRANS =.false.) then
! ------------ WAVEBANDS W/O RADTRANS ----------------------------------
!SL drF is thick
!SL dz_k = drF*hFacC for now lets assume hFacC = 1;  dz_k = thick
!SL Confirm with Sergey or Dima about the grid design
!SL mind kSurface
         dz_k = thick
         hFacC = 1.d0 
         do k=1,Nr
           do ilam = 1,tlam
               if (hFacC(k).gt.0. _d 0) then
! get total attenuation (absorption) by phyto at each wavelength
              actot(k,ilam) = 0.
              actot(k,ilam) = actot(k,ilam)                         &
                        + (phychl_k(1,k)*aphy_chl_k(k,ilam))        &
                        + (phychl_k(2,k)*aphy_chl_dia_k(k,ilam))
if (enable_coccos) then
              actot(k,ilam) = actot(k,ilam)                         &
                        + (phychl_k(3,k)*aphy_chl_cocco_k(k,ilam))  &
                        + (phychl_k(4,k)*aphy_chl_phaeo_k(k,ilam))
endif
              a_k(k,ilam) = aw(ilam)                                &
                           + actot(k,ilam)                          &
                           + acdom_k(k,ilam)                        &
                           + apart_k(k,ilam)
              atten = a_k(k,ilam) * thick(k)    !find drF(k) it is now thick why not dz_k ?
              PARwdn(ilam) = PARwup(ilam)*exp(-atten)
              endif
           enddo !ilam
! find for the midpoint of the gridcell (gridcell mean)
! what could it be instead of hFacC dch w.r.t. drF(k)=:thick
           do ilam = 1,tlam
                   if (hFacC(k).gt.0. _d 0) then
                   PARw_k(ilam,k)=sqrt(PARwup(ilam)                 &
                          * PARwdn(ilam))
                   else
                   PARw_k(ilam,k) = 0. _d 0
                   endif
           enddo
! cycle
           do ilam=1,tlam
                if (hFacC(k).gt.0. _d 0) then
                   PARwup(ilam) = PARwdn(ilam)
                endif
           enddo   !ilam
! sum wavebands for total PAR at the mid point of the gridcell (PARl)
           PARl(k) = 0.
            do ilam = 1,tlam
              PARl(k) = PARl(k) + PARw_k(ilam,k)
            enddo !ilam
!CCEA compute averages and wb=index for exporting
             a_kave(k) = 0.d0
             acdom_kave(k) = 0.d0
             apart_kave(k) = 0.d0
             actot_ave(k) = 0.d0
              do ilam = 1,tlam
                 a_kave(k)     = a_kave(k)                       &
                               + wb_width(ilam) * a_k(k,ilam)
                 acdom_kave(k) = acdom_kave(k)                   &
                               + wb_width(ilam) * acdom_k(k,ilam)
                 apart_kave(k) = apart_kave(k)                   &
                               + wb_width(ilam) * apart_k(k,ilam)
                 actot_ave(k)  = actot_ave(k)                    &
                               + wb_width(ilam) * actot(k,ilam)
              enddo   !ilam
              a_kave(k) = a_kave(k) / wb_totalWidth
             acdom_kave(k) = acdom_kave(k) / wb_totalWidth
             apart_kave(k) = apart_kave(k) / wb_totalWidth
             actot_ave(k) = actot_ave(k) / wb_totalWidth
if (RECOM_CALC_REFLEC) then
             index = darwin_diag_acdom_ilam
             PARw_kwb(k) = PARw_k(index,k)
             a_kwb(k) = a_k(k,index)
             acdom_kwb(k) = acdom_k(k,index)
             apart_kwb(k) = apart_k(k,index)
             actot_wb(k) = actot(k,index)
endif
      enddo       !k
! iops surface for diagnostics
        DO ilam = 1,tlam
              IF( dz_k(kSurface) .GT. 0.0 )THEN
                IF(Eswsf(ilam).GE.tiny .OR.            &
                   Edwsf(ilam).GE.tiny ) THEN
             a_ksur(ilam) = a_k(kSurface,ilam)
             acdom_ksur(ilam) = acdom_k(kSurface,ilam)
             apart_ksur(ilam) = apart_k(kSurface,ilam)
             actot_sur(ilam) = actot(kSurface,ilam)
                ENDIF !light
              ENDIF !depth
        ENDDO  !ilam

else        !/* RECOM_RADTRANS */
!c ------------ FULL RADIATIVE TRANSFER CODE ----------------------------
!CEA Direct and difusse irradiance: both with OASIM, only Ed without
!
! Compute total absorption/scattering coefficients
      DO k=1,Nr
       DO ilam = 1,tlam
!   absorption by phyto: if RADTRANS actot is computed twice
            actot(k,ilam) = 0.0
            bctot(k,ilam) = 0.0
            bbctot(k,ilam) = 0.0
!            DO np = 1,npmax
            actot(k,ilam)  = actot(k,ilam)                       &
                           + phychl_k(1,k) * aphy_chl_k(k,ilam)  &
                      + phychl_k(2,k) * aphy_chl_dia_k(k,ilam)
if (enable_coccos) then
            actot(k,ilam)  = actot(k,ilam)                       &
                      + phychl_k(3,k) * aphy_chl_cocco_k(k,ilam) &
                      + phychl_k(4,k) * aphy_chl_phaeo_k(k,ilam)
endif
if (RECOM_BMASS) then
            bctot(k,ilam)  = bctot(k,ilam)                       &
                           + Phy_k(1,k) * bphy_chl(ilam)         &
                           + Phy_k(2,k) * bphy_chl_dia(ilam)
if (enable_coccos) then
            bctot(k,ilam)  = bctot(k,ilam)                       &
                           + Phy_k(3,k) * bphy_chl_cocco(ilam)   &
                           + Phy_k(4,k) * bphy_chl_phaeo(ilam)
endif
            bbctot(k,ilam) = bbctot(k,ilam)                      &
                           + Phy_k(1,k) * bbphy_chl(ilam)        &
                           + Phy_k(2,k) * bbphy_chl_dia(ilam)
if (enable_coccos) then
            bbctot(k,ilam) = bbctot(k,ilam)                      &
                           + Phy_k(3,k) * bbphy_chl_cocco(ilam)  &
                           + Phy_k(4,k) * bbphy_chl_phaeo(ilam)
endif
else
            bctot(k,ilam)  = bctot(k,ilam)                       &
                           + phychl_k(1,k) * bphy_chl(ilam)      &
                         + phychl_k(2,k) * bphy_chl_dia(ilam)
if (enable_coccos) then
            bctot(k,ilam)  = bctot(k,ilam)                       &
                         + phychl_k(3,k) * bphy_chl_cocco(ilam)  &
                         + phychl_k(4,k) * bphy_chl_phaeo(ilam)
endif
            bbctot(k,ilam) = bbctot(k,ilam)                      &
                         + phychl_k(1,k) * bbphy_chl(ilam)       &
                         + phychl_k(2,k) * bbphy_chl_dia(ilam)
if (enable_coccos) then
            bbctot(k,ilam) = bbctot(k,ilam)                      &
                         + phychl_k(3,k) * bbphy_chl_cocco(ilam) &
                         + phychl_k(4,k) * bbphy_chl_phaeo(ilam)
endif
endif
!            ENDDO
!   total: water, CDOM, phyto, particles
            a_k(k,ilam) = aw(ilam) + acdom_k(k,ilam)             &
                        + actot(k,ilam) + apart_k(k,ilam)
            bt_k(k,ilam) = bw(ilam)                              &
                         + bctot(k,ilam) + bpart_k(k,ilam)
            bb_k(k,ilam) = darwin_bbw * bw(ilam)                 &
                         + bbctot(k,ilam) + bbpart_k(k,ilam)
            bb_k(k,ilam) = MAX(darwin_bbmin, bb_k(k,ilam))
!   initialize output variables
            Edz(ilam,k) = 0.0
            Esz(ilam,k) = 0.0
            Euz(ilam,k) = 0.0
            Estop(ilam,k) = 0.0
            Eutop(ilam,k) = 0.0
            amp1(ilam,k) = 0.0
            amp2(ilam,k) = 0.0
        ENDDO    !ilam
       ENDDO     !k

! ------ Propagate three-beam light in the water column -------
!CEA Some of the routines use drF and others dz_k, why?
         IF (darwin_radtrans_niter.GE.0) THEN
           call MONOD_RADTRANS_ITER(                             &
                    Nr,                                          &
                    dz_k(1:Nr),rmud,                             &
                    Edwsf(1:tlam),                               &
                    Eswsf(1:tlam),                               &
                    a_k(1:Nr,1:tlam),                            &
                    bt_k(1:Nr,1:tlam),                           &
                    bb_k(1:Nr,1:tlam),                           &
                    darwin_radtrans_kmax,darwin_radtrans_niter,  &
                    Edz(1:tlam,1:Nr),                            &
                    Esz(1:tlam,1:Nr),                            &
                    Euz(1:tlam,1:Nr),                            &
                    Eutop(1:tlam,1:Nr),                          &
                    tirrq(1:Nr),                                 &
                    tirrwq(1:tlam,1:Nr),                         &
                    amp1(1:tlam,1:Nr),amp2(1:tlam,1:Nr),         &
                    myThid)
         ELSEIF (darwin_radtrans_niter.EQ.-1) THEN
           call MONOD_RADTRANS(                                  &
                    Nr,                                          &
                    thick(1:Nr),rmud,                            &
                    Edwsf(1:tlam),Eswsf(1:tlam),                 &
                    a_k(1:Nr,1:tlam),                            &
                    bt_k(1:Nr,1:tlam),                           &
                    bb_k(1:Nr,1:tlam),                           &
                    Edz(1:tlam,1:Nr),Esz(1:tlam,1:Nr),           &
                    Euz(1:tlam,1:Nr),Eutop(1:tlam,1:Nr),         &
                    tirrq(1:Nr),                                 &
                    tirrwq(1:tlam,1:Nr),                         &
                    myThid)
         ELSE
            call MONOD_RADTRANS_DIRECT(                          &
                    Nr,                                          &
                    dz_k(1:Nr),rmud,                             &
                    Edwsf(1:tlam),Eswsf(1:tlam),                 &
                    a_k(1:Nr,1:tlam),                            &
                    bt_k(1:Nr,1:tlam),                           &
                    bb_k(1:Nr,1:tlam),                           &
                    darwin_radtrans_kmax,                        &
                    Edz(1:tlam,1:Nr),Esz(1:tlam,1:Nr),           &
                    Euz(1:tlam,1:Nr),                            &
                    Estop(1:tlam,1:Nr),Eutop(1:tlam,1:Nr),       &
                    tirrq(1:Nr),                                 &
                    tirrwq(1:tlam,1:Nr),                         &
                    amp1(1:tlam,1:Nr),amp2(1:tlam,1:Nr),         &
                    myThid)

         ENDIF
!     Uses chl from prev timestep (as wavebands does) keep like this in case
!     need to consider upwelling irradiance as affecting the grid box above
!     Pass to sms: PARw_k only, but will be for this timestep for RADTRANST
!     and previous timestep for WAVEBANDS.
!
!     Now copy
       DO k=1,Nr
                       PARl(k) = tirrq(k)
           DO ilam = 1,tlam
                       PARw_k(ilam,k) = tirrwq(ilam,k)
           ENDDO  !ilam
!CCEA Compute averages and wb=index for exporting
             a_kave(k) = 0.d0
             acdom_kave(k) = 0.d0
             apart_kave(k) = 0.d0
             actot_ave(k) = 0.d0

             bt_kave(k) = 0.d0
             bpart_kave(k) = 0.d0
             bctot_ave(k) = 0.d0
             bb_kave(k) = 0.d0
             bbpart_kave(k) = 0.d0
             bbctot_ave(k) = 0.d0
              do ilam = 1,tlam
                 a_kave(k) = a_kave(k)                         &
                           + wb_width(ilam) * a_k(k,ilam)
                 acdom_kave(k) = acdom_kave(k)                 &
                               + wb_width(ilam) * acdom_k(k,ilam)
                 apart_kave(k) = apart_kave(k)                 &
                               + wb_width(ilam) * apart_k(k,ilam)
                 actot_ave(k) = actot_ave(k)                   &
                              + wb_width(ilam) * actot(k,ilam)

                 bt_kave(k) = bt_kave(k)                       &
                            + wb_width(ilam) * bt_k(k,ilam)
                 bpart_kave(k) = bpart_kave(k)                 &
                               + wb_width(ilam) * bpart_k(k,ilam)
                 bctot_ave(k) = bctot_ave(k)                   &
                              + wb_width(ilam) * bctot(k,ilam)

                 bb_kave(k) = bb_kave(k)                       &
                            + wb_width(ilam) * bb_k(k,ilam)
                 bbpart_kave(k) = bbpart_kave(k)               &
                                + wb_width(ilam) * bbpart_k(k,ilam)
                 bbctot_ave(k) = bbctot_ave(k)                 &
                               + wb_width(ilam) * bbctot(k,ilam)

              enddo   !ilam
             a_kave(k) = a_kave(k) / wb_totalWidth
             acdom_kave(k) = acdom_kave(k) / wb_totalWidth
             apart_kave(k) = apart_kave(k) / wb_totalWidth
             actot_ave(k) = actot_ave(k) / wb_totalWidth

             bt_kave(k) = bt_kave(k) / wb_totalWidth
             bpart_kave(k) = bpart_kave(k) / wb_totalWidth
             bctot_ave(k) = bctot_ave(k) / wb_totalWidth

             bb_kave(k) = bb_kave(k) / wb_totalWidth
             bbpart_kave(k) = bbpart_kave(k) / wb_totalWidth
             bbctot_ave(k) = bbctot_ave(k) / wb_totalWidth

if (RECOM_CALC_REFLEC)
             index = darwin_diag_acdom_ilam
             PARw_kwb(k) = PARw_k(index,k)
             a_kwb(k) = a_k(k,index)
             acdom_kwb(k) = acdom_k(k,index)
             apart_kwb(k) = apart_k(k,index)
             actot_wb(k) = actot(k,index)

             bt_kwb(k) = bt_k(k,index)
             bpart_kwb(k) = bpart_k(k,index)
             bctot_wb(k) = bctot(k,index)

             bb_kwb(k) = bb_k(k,index)
             bbpart_kwb(k) = bbpart_k(k,index)
             bbctot_wb(k) = bbctot(k,index)

             Edz_wb(k) = Edz(index,k)
             Esz_wb(k) = Esz(index,k)
             Euz_wb(k) = Euz(index,k)
             Estop_wb(k) = Estop(index,k)
             Eutop_wb(k) = Eutop(index,k)
             amp1_wb(k) = amp1(index,k)
             amp2_wb(k) = amp2(index,k)
             amp2_wb(k) = amp2(index,k)
endif            
      ENDDO        !k
!     PARw and PARwup from non-spectral RECOM are from previous timestep
!     (attenuation done in recom_sms) but PARw and PARwup from WAVEBANDS
!     and RADTRANS are for the current timestep.


! ----------------- Reflectance -------------------------         
!douple check w.r.t. existance name dz_k(kSurface) 
        DO ilam = 1,tlam
              IF( dz_k(kSurface) .GT. 0.0 )THEN
                IF(Eswsf(ilam).GE.darwin_radmodThresh .OR.    &
                   Edwsf(ilam).GE.darwin_radmodThresh ) THEN
                Eupwel(ilam) = Eutop(ilam,kSurface)
                Reflec(ilam) = Eupwel(ilam) /                 &
                        (Eswsf(ilam) + Edwsf(ilam))
!     iops surface
             a_ksur(ilam) = a_k(kSurface,ilam)
             acdom_ksur(ilam) = acdom_k(kSurface,ilam)
             apart_ksur(ilam) = apart_k(kSurface,ilam)
             actot_sur(ilam) = actot(kSurface,ilam)
             bt_ksur(ilam) = bt_k(kSurface,ilam)
             bpart_ksur(ilam) = bpart_k(kSurface,ilam)
             bctot_sur(ilam) = bctot(kSurface,ilam)
             bb_ksur(ilam) = bb_k(kSurface,ilam)
             bbpart_ksur(ilam) = bbpart_k(kSurface,ilam)
             bbctot_sur(ilam) = bbctot(kSurface,ilam)
                ENDIF !light
              ENDIF !depth         
          ENDDO   !ilam   
endif        !/* RECOM_RADTRANS */       
#endif /* RECOM_WAVEBANDS */
!======================================================================
    

    !! *** Mocsy ***

    !!---- convert from mmol/m3 to mol/m3
    REcoM_DIC  = max(tiny*1e-3, state(one,idic)*1e-3)
    REcoM_Alk  = max(tiny*1e-3, state(one,ialk)*1e-3)
    REcoM_Si   = max(tiny*1e-3, state(one,isi) *1e-3)

    !!---- convert N to P with Redfield ratio
    REcoM_Phos = max(tiny*1e-3, state(one,idin)*1e-3) /16.

    !!---- minimum set to 2 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
    REcoM_T    = max(2.d0, Temp(1))
    !!---- maximum set to 40 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
    REcoM_T    = min(REcoM_T, 40.d0) 

    !!---- minimum set to 21: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble in regions with S between 19 and 21 and ice conc above 97%
    REcoM_S    = max(21.d0, Sali(1)) 
    !!---- maximum set to 43: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble   REcoM_S    = min(REcoM_S, 43.d0)  !!!!!!!!

    !!---- convert from Pa to atm.
    Patm = Loc_slp/Pa2atm  

    !!---- lon
    Lond=geo_coord_nod2D(1,n)/rad !! convert from rad to degree
    !!---- lat
    Latr=geo_coord_nod2D(2,n)
    Latd=geo_coord_nod2D(2,n)/rad !! convert from rad to degree

    !!---- calculate piston velocity kw660, which is an input to the flxco2 calculation
    !!---- pistonvel already scaled for ice-free area
    !!---- compute piston velolicty kw660 (at 25 C) from wind speed
    !!---- BUT without Schmidt number temperature correction (Sc differs each gas)
    !! ULoc: wind speed at 10-m height
    !! Loc_ice_conc: modeled sea-ice cover: fraction of grid cell, varying between 0.0 (no ice) and 1.0 (full cover)
    !! kw660: piston velocity at 25°C [m/s], uncorrected by the Schmidt number for different temperatures

    call pistonvel(ULoc, Loc_ice_conc, Nmocsy, kw660)

    !! *** check ***

    if((REcoM_DIC(1) > 10000.d0)) then               ! NEW: added this entire print statement (if to endif)
        print*, 'NEW ERROR: DIC !'  
        print*, 'pco2surf: ',pco2surf
        print*, 'co2: ',co2
        print*, 'rhoSW: ', rhoSW
        print*, 'temp: ',REcoM_T
        print*, 'tempis: ',tempis
        print*, 'REcoM_S: ', REcoM_S
        print*, 'REcoM_Alk: ', REcom_Alk
        print*, 'REcoM_DIC: ', REcoM_DIC
        print*, 'REcoM_Si: ', REcoM_Si
        print*, 'REcoM_Phos: ', REcoM_Phos
        print*, 'kw660: ',kw660
        print*, 'LocAtmCO2: ', LocAtmCO2
        print*, 'Patm: ', Patm
        print*, 'thick(One): ',thick(One)
        print*, 'Nmocsy: ', Nmocsy
        print*, 'Lond: ', Lond
        print*, 'Latd: ', Latd
        print*, 'ULoc: ', ULoc
        print*, 'Loc_ice_conc: ', Loc_ice_conc
        stop
    endif

    call flxco2(co2flux, co2ex, dpco2surf,                                                       &
                ph, pco2surf, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis, K0, &
                REcoM_T, REcoM_S, REcoM_Alk, REcoM_DIC, REcoM_Si, REcoM_Phos, kw660, LocAtmCO2, Patm, thick(One), Nmocsy, Lond,Latd, &
                optCON='mol/m3',optT='Tpot   ',optP='m ',optB='u74',optK1K2='l  ',optKf='dg',optGAS='Pinsitu',optS='Sprc')

! changed optK1K2='l  ' to 'm10'
  if((co2flux(1)>1.e10) .or. (co2flux(1)<-1.e10)) then
!     co2flux(1)=0.0  
      print*, 'ERROR: co2 flux !'
      print*, 'pco2surf: ',pco2surf
      print*, 'co2: ',co2
      print*, 'rhoSW: ', rhoSW
      print*, 'temp: ',REcoM_T
      print*, 'tempis: ',tempis
      print*, 'REcoM_S: ', REcoM_S
      print*, 'REcoM_Alk: ', REcom_Alk
      print*, 'REcoM_DIC: ', REcoM_DIC
      print*, 'REcoM_Si: ', REcoM_Si
      print*, 'REcoM_Phos: ', REcoM_Phos
      print*, 'kw660: ',kw660
      print*, 'LocAtmCO2: ', LocAtmCO2
      print*, 'Patm: ', Patm
      print*, 'thick(One): ',thick(One) 
      print*, 'Nmocsy: ', Nmocsy
      print*, 'Lond: ', Lond
      print*, 'Latd: ', Latd   
      print*, 'ULoc: ', ULoc
      print*, 'Loc_ice_conc: ', Loc_ice_conc
      stop
    endif

! use ice-free area and also convert from mol/m2/s to mmol/m2/d
!   if(mype==0) write(*,*), 'co2flux (mol/m2/s) =',co2flux

! ice-fraction is already considered in piston-velocity, so don't apply it here
   dflux     = co2flux * 1.e3 *SecondsPerDay  !* (1.d0 - Loc_ice_conc)
!   if(mype==0) write(*,*), 'dflux (mmol/m2/d) =',dflux

   co2flux_seaicemask = co2flux * 1.e3 !  [mmol/m2/s]  * (1.d0 - Loc_ice_conc)
!   if(mype==0) write(*,*), 'co2flux_seaicemask (mmol/m2/s) =',co2flux_seaicemask

! then oxygen
   ppo = Loc_slp/Pa2atm !1 !slp divided by 1 atm
   REcoM_O2 = max(tiny*1e-3,state(one,ioxy)*1e-3) ! convert from mmol/m3 to mol/m3 for mocsy

   call  o2flux(REcoM_T, REcoM_S, kw660, ppo, REcoM_O2, Nmocsy, o2ex)
   oflux     = o2ex * 1.e3 *SecondsPerDay  !* (1.d0 - Loc_ice_conc) [mmol/m2/d]
   o2flux_seaicemask = o2ex * 1.e3 ! back to mmol here [mmol/m2/s] 

! Source-Minus-Sinks

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> REcoM_sms'//achar(27)//'[0m'

!  call REcoM_sms(n, Nn, state, thick, recipthick, SurfSW, sms, Temp ,zF, PAR, mesh)

  call REcoM_sms(n, Nn, state, thick, recipthick, SurfSW, sms, Temp, Sali_depth &
        , CO2_watercolumn                                              & ! MOCSY [mol/m3]
        , pH_watercolumn                                               & ! MOCSY on total scale
        , pCO2_watercolumn                                             & ! MOCSY [uatm]
        , HCO3_watercolumn                                             & ! MOCSY [mol/m3]
        , CO3_watercolumn                                              & ! DISS [mol/m3]
        , OmegaC_watercolumn                                           & ! DISS calcite saturation state
        , kspc_watercolumn                                             & ! DISS stoichiometric solubility product [mol^2/kg^2]
        , rhoSW_watercolumn                                            & ! DISS in-situ density of seawater [kg/m3]
        , Loc_slp                                                      &
        , zF, PAR                                                      &
#ifdef RECOM_WAVEBANDS   
        ,PARlocal                                                            &
        ,PARwlocal                                                           &
        ,C_phot_nl                                                           &
        ,C_phot_nl_dia                                                       &
if (enable_coccos) then        
        ,C_phot_nl_cocco                                                     &
        ,C_phot_nl_phaeo                                                     &
endif                                                        
        ,Ek_nl                                                               &
        ,Ek_nl_dia
if (enable_coccos) then                                                      &
        ,Ek_nl_cocco                                                         &
        ,Ek_nl_phaeo
endif
#endif        
        , Lond, Latd, ice, dynamics, tracers, partit, mesh)

  state(1:nn,:)      = max(tiny,state(1:nn,:) + sms(1:nn,:))

  state(1:nn,ipchl)  = max(tiny_chl,state(1:nn,ipchl))
  state(1:nn,iphyn)  = max(tiny_N,  state(1:nn,iphyn))
  state(1:nn,iphyc)  = max(tiny_C,  state(1:nn,iphyc))
  state(1:nn,idchl)  = max(tiny_chl,state(1:nn,idchl))
  state(1:nn,idian)  = max(tiny_N_d,state(1:nn,idian))
  state(1:nn,idiac)  = max(tiny_C_d,state(1:nn,idiac))
  state(1:nn,idiasi) = max(tiny_Si, state(1:nn,idiasi))

if (enable_coccos) then
  state(1:nn,icchl)  = max(tiny_chl,state(1:nn,icchl))
  state(1:nn,icocn)  = max(tiny_N_c,state(1:nn,icocn))
  state(1:nn,icocc)  = max(tiny_C_c,state(1:nn,icocc))

  state(1:nn,iphachl)  = max(tiny_chl,state(1:nn,iphachl))
  state(1:nn,iphan)  = max(tiny_N_p,state(1:nn,iphan))
  state(1:nn,iphac)  = max(tiny_C_p,state(1:nn,iphac))
endif

if (enable_3zoo2det) then
  state(1:nn,imiczoon)  = max(tiny,state(1:nn,imiczoon))
  state(1:nn,imiczooc)  = max(tiny,state(1:nn,imiczooc))
endif

#ifdef RECOM_WAVEBANDS
!SL------CDOM and Marshal related --------------------------
if (RECOM_CDOM) then
  state(1:nn,idcdom)  = max(tiny,state(1:nn,idcdom))
endif
if (RECOM_CALC_APHYT=.true. .and. RECOM_MARSHALL =.true.)
  state(1:nn,id1)  = max(tiny,state(1:nn,id1)) !rel
  state(1:nn,id1d)  = max(tiny,state(1:nn,id1d)) !rel
endif
#endif /* RECOM_WAVEBANDS */


if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> ciso after REcoM_Forcing'//achar(27)//'[0m'

  if (ciso) then
!   Calculate carbon-isotopic fractionation, radioactive decay is calculated in oce_ale_tracer.F90

!   Fractionation due to air-sea exchange and chemical speciation of CO2
    call recom_ciso_airsea(recom_t(1), co3(1), recom_dic(1)) ! -> alpha_aq, alpha_dic. CO3 is taken from mocsy

!   Isotopic ratios of dissolved CO2, also needed to calculate biogenic fractionation
    r_dic_13     = max(tiny*1e-3,state(1,idic_13)*1e-3) / recom_dic(1)
    r_co2s_13    = alpha_aq_13 / alpha_dic_13 * r_dic_13
!   Calculate air-sea fluxes of 13|14CO2 in mmol / m**2 / s
    kwco2  = kw660(1) * (660/scco2(REcoM_T(1)))**0.5  ! Piston velocity (via mocsy)
    co2sat = co2flux(1) / (kwco2 + tiny) + co2(1)     ! Saturation concentration of CO2 (via mocsy)
!   co2flux_13   = kwco2 * alpha_k_13 * (alpha_aq_13 * r_atm_13 * co2sat - r_co2s_13 * co2(1))
!   co2flux_13   = alpha_k_13 * alpha_aq_13 * kwco2 * (r_atm_13 * co2sat - r_dic_13 * co2(1) / alpha_dic_13)
!   Fractionation factors were determined for freshwater, include a correction for enhanced fractionation in seawater
    co2flux_13   = (alpha_k_13 * alpha_aq_13 - 0.0002) * kwco2 * (r_atm_13 * co2sat - r_dic_13 * co2(1) / alpha_dic_13)
    co2flux_seaicemask_13 = co2flux_13 * 1.e3

!   Biogenic fractionation due to photosynthesis of plankton
!   phyc_13|14 and diac_13|14 are only used in REcoM_sms to calculate DIC_13|14, DOC_13|14 and DetC_13|14

    call recom_ciso_photo(co2(1)) ! -> alpha_p
    r_phyc_13 = r_co2s_13 / alpha_p_13
    r_diac_13 = r_co2s_13 / alpha_p_dia_13
    state(1:nn,iphyc_13)   = max((tiny_C   * r_phyc_13), (state(1:nn,iphyc) * r_phyc_13))
    state(1:nn,idiac_13)   = max((tiny_C_d * r_diac_13), (state(1:nn,idiac) * r_diac_13))

!   The same for radiocarbon, fractionation factors have been already derived above
    if (ciso_14) then
!   Air-sea exchange
      r_dic_14   = max(tiny*1e-3,state(1,idic_14)*1e-3) / recom_dic(1)
      r_co2s_14  = alpha_aq_14 / alpha_dic_14 * r_dic_14
!     co2flux_14 = kwco2 * alpha_k_14 * (alpha_aq_14 * r_atm_14 * co2sat - r_co2s_14 * co2(1))
!     Fractionation factors were determined for freshwater, include a correction for enhanced fractionation seawater
      co2flux_14 = (alpha_k_14 * alpha_aq_14 - 0.0004) * kwco2 * (r_atm_14 * co2sat - r_dic_14 * co2(1) / alpha_dic_14)
      co2flux_seaicemask_14 = co2flux_14 * 1.e3
!   Biogenic fractionation
      if (ciso_organic_14) then
        r_phyc_14 = r_co2s_14 / alpha_p_14
        r_diac_14 = r_co2s_14 / alpha_p_dia_14
        state(1:nn,iphyc_14) = max((tiny_C   * r_phyc_14), (state(1:nn,iphyc) * r_phyc_14))
        state(1:nn,idiac_14) = max((tiny_C_d * r_diac_14), (state(1:nn,idiac) * r_diac_14))
      end if
    end if
!   Radiocarbon
  end if
! ciso

!-------------------------------------------------------------------------------
! Diagnostics
  if (Diags) then

!    logical, optional                 :: lNPPn

!    if (present(lNPPn))then
!        locNPPn = sum(diags3Dloc(1:nn,idiags) * thick(1:nn))
!    endif
     locNPPn = sum(vertNPPn(1:nn) * thick(1:nn))
     locGPPn = sum(vertGPPn(1:nn) * thick(1:nn))
     locNNAn = sum(vertNNAn(1:nn) * thick(1:nn))
     locChldegn = sum(vertChldegn(1:nn) * thick(1:nn))

     locNPPd = sum(vertNPPd(1:nn) * thick(1:nn))
     locGPPd = sum(vertGPPd(1:nn) * thick(1:nn))
     locNNAd = sum(vertNNAd(1:nn) * thick(1:nn))
     locChldegd = sum(vertChldegd(1:nn) * thick(1:nn))

if (enable_coccos) then
     locNPPc = sum(vertNPPc(1:nn) * thick(1:nn))
     locGPPc = sum(vertGPPc(1:nn) * thick(1:nn))
     locNNAc = sum(vertNNAc(1:nn) * thick(1:nn))
     locChldegc = sum(vertChldegc(1:nn) * thick(1:nn))

     locNPPp = sum(vertNPPp(1:nn) * thick(1:nn))
     locGPPp = sum(vertGPPp(1:nn) * thick(1:nn))
     locNNAp = sum(vertNNAp(1:nn) * thick(1:nn))
     locChldegp = sum(vertChldegp(1:nn) * thick(1:nn))
endif

  end if
end subroutine REcoM_Forcing
