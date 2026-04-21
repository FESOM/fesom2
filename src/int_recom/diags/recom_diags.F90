! ==============================================================================
! MODULE: recom_diags_management
! Purpose: Allocation, initialization, update, and deallocation of all REcoM
!          diagnostic arrays (per-node water column and global 2D/3D fields).
! ==============================================================================
module recom_diags_management
    use recom_config
    implicit none

contains

! ==============================================================================
! SUBROUTINE: allocate_and_init_diags
! Purpose: Allocate and initialize all local (per-column) diagnostic arrays
! ==============================================================================
subroutine allocate_and_init_diags(nl)
    use recom_locvar
    use REcoM_declarations
    implicit none
    
    integer, intent(in) :: nl  ! Number of vertical levels

    ! --------------------------------------------------------------------------
    ! Small Phytoplankton & Diatoms (always active)
    ! --------------------------------------------------------------------------
    call alloc_init_phyto_diags(nl)

    ! --------------------------------------------------------------------------
    ! Coccolithophores and Phaeocystis (optional)
    ! --------------------------------------------------------------------------
    if (enable_coccos) call alloc_init_cocco_diags(nl)

    ! Calcification arrays (always needed)
    allocate(vertcalcdiss(nl-1), vertcalcif(nl-1))
    vertcalcdiss = 0.d0
    vertcalcif   = 0.d0

    ! --------------------------------------------------------------------------
    ! Zooplankton grazing (optional)
    ! --------------------------------------------------------------------------
    if (Grazing_detritus) call alloc_init_zoo_diags(nl)
 
    ! --------------------------------------------------------------------------
    ! Dissolution, remineralization, and carbonate chemistry
    ! --------------------------------------------------------------------------
    call alloc_init_bgc_diags(nl)
 
end subroutine allocate_and_init_diags

! ==============================================================================
! SUBROUTINE: alloc_init_phyto_diags
! Purpose: Allocate and initialize small phytoplankton and diatom diagnostics
! ==============================================================================
subroutine alloc_init_phyto_diags(nl)
    use recom_locvar
    use recom_declarations
    implicit none
 
    integer, intent(in) :: nl

    ! --------------------------------------------------------------------------
    ! Small Phytoplankton
    ! --------------------------------------------------------------------------
    allocate(vertNPPn(nl-1), vertGPPn(nl-1), vertNNAn(nl-1), vertChldegn(nl-1))
    allocate(vertrespn(nl-1), vertdocexn(nl-1), vertaggn(nl-1))

    vertNPPn    = 0.d0
    vertGPPn    = 0.d0
    vertNNAn    = 0.d0
    vertChldegn = 0.d0
    vertrespn   = 0.d0
    vertdocexn  = 0.d0
    vertaggn    = 0.d0

    ! --------------------------------------------------------------------------
    ! Diatoms
    ! --------------------------------------------------------------------------
    allocate(vertNPPd(nl-1), vertGPPd(nl-1), vertNNAd(nl-1), vertChldegd(nl-1))
    allocate(vertrespd(nl-1), vertdocexd(nl-1), vertaggd(nl-1))

    vertNPPd    = 0.d0
    vertGPPd    = 0.d0
    vertNNAd    = 0.d0
    vertChldegd = 0.d0
    vertrespd   = 0.d0
    vertdocexd  = 0.d0
    vertaggd    = 0.d0

    ! --------------------------------------------------------------------------
    ! Photosynthesis and Nitrogen Assimilation
    ! --------------------------------------------------------------------------
    allocate(vertphotn(nl-1), vertphotd(nl-1))
    allocate(vertNassimn(nl-1), vertNassimd(nl-1))
    vertphotn = 0.d0
    vertphotd = 0.d0
    vertNassimn = 0.d0
    vertNassimd = 0.d0

    ! --------------------------------------------------------------------------
    ! DOM Remineralisation
    ! --------------------------------------------------------------------------
    allocate(vertDONremin(nl-1), vertDOCremin(nl-1))
    vertDONremin = 0.d0
    vertDOCremin = 0.d0

    ! --------------------------------------------------------------------------
    ! Temperature and Photosynthesis Tracking Variables
    ! --------------------------------------------------------------------------

    allocate(VTPhyCO2(nl-1),  VTDiaCO2(nl-1))
    allocate(VTCphotLigLim_phyto(nl-1),   VTCphotLigLim_diatoms(nl-1))
    allocate(VTCphot_phyto(nl-1),         VTCphot_diatoms(nl-1))
    allocate(VTTemp_phyto(nl-1),          VTTemp_diatoms(nl-1))
    allocate(VTqlimitFac_phyto(nl-1),     VTqlimitFac_diatoms(nl-1))
    allocate(VTSi_assimDia(nl-1))

    VTPhyCO2  = 0.d0
    VTDiaCO2  = 0.d0
    VTCphotLigLim_phyto   = 0.d0
    VTCphotLigLim_diatoms = 0.d0
    VTCphot_phyto   = 0.d0
    VTCphot_diatoms = 0.d0
    VTTemp_phyto   = 0.d0
    VTTemp_diatoms = 0.d0
    VTqlimitFac_phyto   = 0.d0
    VTqlimitFac_diatoms = 0.d0
    VTSi_assimDia = 0.d0
end subroutine alloc_init_phyto_diags

! ==============================================================================
! SUBROUTINE: alloc_init_cocco_diags
! Purpose: Allocate and initialize coccolithophore and Phaeocystis diagnostics
! ==============================================================================
subroutine alloc_init_cocco_diags(nl)
    use recom_locvar
    use recom_declarations
    implicit none
 
    integer, intent(in) :: nl

    ! ----------------------------------------------------------------------
    ! Coccolithophores
    ! ----------------------------------------------------------------------
    allocate(vertNPPc(nl-1), vertGPPc(nl-1), vertNNAc(nl-1), vertChldegc(nl-1))
    allocate(vertrespc(nl-1), vertdocexc(nl-1), vertaggc(nl-1))
    vertNPPc     = 0.d0
    vertGPPc     = 0.d0
    vertNNAc     = 0.d0
    vertChldegc  = 0.d0
    vertrespc    = 0.d0
    vertdocexc   = 0.d0
    vertaggc     = 0.d0
 
    ! ----------------------------------------------------------------------
    ! Phaeocystis
    ! ----------------------------------------------------------------------
    allocate(vertNPPp(nl-1), vertGPPp(nl-1), vertNNAp(nl-1), vertChldegp(nl-1))
    allocate(vertrespp(nl-1), vertdocexp(nl-1), vertaggp(nl-1))
    vertNPPp    = 0.d0
    vertGPPp    = 0.d0
    vertNNAp    = 0.d0
    vertChldegp = 0.d0
    vertrespp   = 0.d0
    vertdocexp  = 0.d0
    vertaggp    = 0.d0

    ! --------------------------------------------------------------------------
    ! Photosynthesis and Nitrogen Assimilation
    ! --------------------------------------------------------------------------
    allocate(vertphotc(nl-1), vertphotp(nl-1))
    allocate(vertNassimc(nl-1), vertNassimp(nl-1))
    vertphotc = 0.d0
    vertphotp = 0.d0
    vertNassimc = 0.d0
    vertNassimp = 0.d0

    ! --------------------------------------------------------------------------
    ! Temperature / photosynthesis tracking
    ! --------------------------------------------------------------------------
    allocate(VTTemp_cocco(nl-1),       VTTemp_phaeo(nl-1))
    allocate(VTCoccoCO2(nl-1),         VTPhaeoCO2(nl-1))
    allocate(VTqlimitFac_cocco(nl-1),  VTqlimitFac_phaeo(nl-1))
    allocate(VTCphotLigLim_cocco(nl-1),VTCphotLigLim_phaeo(nl-1))
    allocate(VTCphot_cocco(nl-1),      VTCphot_phaeo(nl-1))
    VTTemp_cocco = 0.d0
    VTTemp_phaeo = 0.d0
    VTCoccoCO2 = 0.d0
    VTPhaeoCO2 = 0.d0
    VTqlimitFac_cocco = 0.d0
    VTqlimitFac_phaeo = 0.d0
    VTCphotLigLim_cocco = 0.d0
    VTCphotLigLim_phaeo = 0.d0
    VTCphot_cocco = 0.d0
    VTCphot_phaeo = 0.d0
end subroutine alloc_init_cocco_diags

! ==============================================================================
! SUBROUTINE: alloc_init_zoo_diags
! Purpose: Allocate and initialize zooplankton grazing diagnostics
! ==============================================================================
subroutine alloc_init_zoo_diags(nl)
    use recom_locvar
    use recom_declarations
    implicit none
 
    integer, intent(in) :: nl

    ! --------------------------------------------------------------------------
    ! Mesozooplankton (always present when Grazing_detritus is on)
    ! --------------------------------------------------------------------------
    allocate(vertgrazmeso_tot(nl-1), vertgrazmeso_n(nl-1), vertgrazmeso_d(nl-1))
    allocate(vertgrazmeso_det(nl-1), vertrespmeso(nl-1))
    vertgrazmeso_tot  = 0.d0
    vertgrazmeso_n    = 0.d0
    vertgrazmeso_d    = 0.d0
    vertgrazmeso_det  = 0.d0
    vertrespmeso      = 0.d0

    ! --------------------------------------------------------------------------
    ! Mesozooplankton Respiration
    ! --------------------------------------------------------------------------
    allocate(vertmesocdis(nl-1))
    vertmesocdis = 0.d0

    if (.not. enable_3zoo2det) return

    ! --------------------------------------------------------------------------
    ! Microzooplankton
    ! --------------------------------------------------------------------------
    allocate(vertgrazmicro_tot(nl-1), vertgrazmicro_n(nl-1), vertgrazmicro_d(nl-1))
    allocate(vertrespmicro(nl-1))             
    vertgrazmicro_tot = 0.d0
    vertgrazmicro_n   = 0.d0
    vertgrazmicro_d   = 0.d0
    vertrespmicro     = 0.d0

    ! --------------------------------------------------------------------------
    ! Microzooplankton Dissolution
    ! --------------------------------------------------------------------------
    allocate( vertmicrocdis(nl-1))
    vertmicrocdis   = 0.d0

    ! --------------------------------------------------------------------------
    ! Additional mesozooplankton arrays (3-zoo config)
    ! --------------------------------------------------------------------------
    allocate(vertgrazmeso_mic(nl-1), vertgrazmeso_det2(nl-1))
    vertgrazmeso_mic  = 0.d0
    vertgrazmeso_det2 = 0.d0

    ! --------------------------------------------------------------------------
    ! Macrozooplankton
    ! --------------------------------------------------------------------------
    allocate(vertgrazmacro_tot(nl-1), vertgrazmacro_n(nl-1),   vertgrazmacro_d(nl-1))
    allocate(vertgrazmacro_mes(nl-1), vertgrazmacro_det(nl-1))
    allocate(vertgrazmacro_mic(nl-1), vertgrazmacro_det2(nl-1))
    allocate(vertrespmacro(nl-1))
    vertgrazmacro_tot  = 0.d0
    vertgrazmacro_n    = 0.d0
    vertgrazmacro_d    = 0.d0
    vertgrazmacro_mes  = 0.d0
    vertgrazmacro_det  = 0.d0
    vertgrazmacro_mic  = 0.d0
    vertgrazmacro_det2 = 0.d0
    vertrespmacro      = 0.d0

    ! --------------------------------------------------------------------------
    ! Macrozooplankton Dissolution
    ! --------------------------------------------------------------------------
    allocate(vertmacrocdis(nl-1))
    allocate(vertfastcdis(nl-1))
    vertmacrocdis = 0.d0
    vertfastcdis = 0.d0

    if (.not. enable_coccos) return
 
    allocate(vertgrazmicro_c(nl-1), vertgrazmicro_p(nl-1))
    allocate(vertgrazmeso_c(nl-1),  vertgrazmeso_p(nl-1))
    allocate(vertgrazmacro_c(nl-1), vertgrazmacro_p(nl-1))
    vertgrazmicro_c = 0.d0
    vertgrazmicro_p = 0.d0
    vertgrazmeso_c  = 0.d0
    vertgrazmeso_p  = 0.d0
    vertgrazmacro_c = 0.d0
    vertgrazmacro_p = 0.d0

end subroutine alloc_init_zoo_diags

! ==============================================================================
! SUBROUTINE: alloc_init_bgc_diags
! Purpose: Allocate and initialize dissolution/remineralization diagnostics
! ==============================================================================
subroutine alloc_init_bgc_diags(nl)
    use REcoM_declarations
    implicit none
 
    integer, intent(in) :: nl
 
    allocate(vertDISSOC(nl-1), vertDISSON(nl-1), vertDISSOSi(nl-1))
    allocate(vertREMOC(nl-1),  vertREMOCt(nl-1), vertREMON(nl-1))
    vertDISSOC  = 0.d0
    vertDISSON  = 0.d0
    vertDISSOSi = 0.d0
    vertREMOC   = 0.d0
    vertREMOCt  = 0.d0
    vertREMON   = 0.d0
end subroutine alloc_init_bgc_diags

! ==============================================================================
! SUBROUTINE: update_2d_diags
! Purpose: Transfer local diagnostic values to 2D global arrays
! ==============================================================================
subroutine update_2d_diags(n)
    use recom_locvar
    use recom_glovar
    use REcoM_declarations
    implicit none
    
    integer, intent(in) :: n  ! Node index
    
    ! --------------------------------------------------------------------------
    ! Small Phytoplankton
    ! --------------------------------------------------------------------------
    NPPn(n)    = locNPPn
    GPPn(n)    = locGPPn
    NNAn(n)    = locNNAn
    Chldegn(n) = locChldegn
    
    ! --------------------------------------------------------------------------
    ! Diatoms
    ! --------------------------------------------------------------------------
    NPPd(n)    = locNPPd
    GPPd(n)    = locGPPd
    NNAd(n)    = locNNAd
    Chldegd(n) = locChldegd
    
    ! --------------------------------------------------------------------------
    ! Coccolithophores and Phaeocystis (if enabled)
    ! --------------------------------------------------------------------------
    if (enable_coccos) then
        NPPc(n)    = locNPPc
        GPPc(n)    = locGPPc
        NNAc(n)    = locNNAc
        Chldegc(n) = locChldegc
        
        NPPp(n)    = locNPPp
        GPPp(n)    = locGPPp
        NNAp(n)    = locNNAp
        Chldegp(n) = locChldegp
    endif
    
    ! --------------------------------------------------------------------------
    ! Zooplankton Grazing (if enabled)
    ! --------------------------------------------------------------------------
    if (Grazing_detritus) then
        ! Mesozooplankton
        grazmeso_tot(n) = locgrazmeso_tot
        grazmeso_n(n)   = locgrazmeso_n
        grazmeso_d(n)   = locgrazmeso_d
        grazmeso_det(n) = locgrazmeso_det
        
        if (enable_coccos) then
            grazmeso_c(n) = locgrazmeso_c
            grazmeso_p(n) = locgrazmeso_p
        endif
        
        if (enable_3zoo2det) then
            grazmeso_mic(n)  = locgrazmeso_mic
            grazmeso_det2(n) = locgrazmeso_det2
            
            ! Macrozooplankton
            grazmacro_tot(n)  = locgrazmacro_tot
            grazmacro_n(n)    = locgrazmacro_n
            grazmacro_d(n)    = locgrazmacro_d
            grazmacro_mes(n)  = locgrazmacro_mes
            grazmacro_det(n)  = locgrazmacro_det
            grazmacro_mic(n)  = locgrazmacro_mic
            grazmacro_det2(n) = locgrazmacro_det2
            
            if (enable_coccos) then
                grazmacro_c(n) = locgrazmacro_c
                grazmacro_p(n) = locgrazmacro_p
            endif
            
            ! Microzooplankton
            grazmicro_tot(n) = locgrazmicro_tot
            grazmicro_n(n)   = locgrazmicro_n
            grazmicro_d(n)   = locgrazmicro_d
            
            if (enable_coccos) then
                grazmicro_c(n) = locgrazmicro_c
                grazmicro_p(n) = locgrazmicro_p
            endif
        endif
    endif

    ! Dissolution and remineralization ! R2OMIP
     DISSOC(n)  = locDISSOC
     DISSON(n)  = locDISSON
     DISSOSi(n) = locDISSOSi
     REMOC(n)   = locREMOC
     REMOCt(n)  = locREMOCt
     REMON(n)   = locREMON
    
end subroutine update_2d_diags

! ==============================================================================
! SUBROUTINE: update_3d_diags
! Purpose: Transfer vertical profile diagnostic values to 3D global arrays
! ==============================================================================
subroutine update_3d_diags(n, nzmax)
    use recom_locvar
    use recom_glovar
    use REcoM_declarations
    implicit none
    
    integer, intent(in) :: n      ! Node index
    integer, intent(in) :: nzmax  ! Maximum vertical level for this node
    
    ! --------------------------------------------------------------------------
    ! Small Phytoplankton
    ! --------------------------------------------------------------------------
    aggn(1:nzmax,n)   = vertaggn(1:nzmax)
    docexn(1:nzmax,n) = vertdocexn(1:nzmax)
    respn(1:nzmax,n)  = vertrespn(1:nzmax)
    NPPn3D(1:nzmax,n) = vertNPPn(1:nzmax)
    
    ! --------------------------------------------------------------------------
    ! Diatoms
    ! --------------------------------------------------------------------------
    aggd(1:nzmax,n)   = vertaggd(1:nzmax)
    docexd(1:nzmax,n) = vertdocexd(1:nzmax)
    respd(1:nzmax,n)  = vertrespd(1:nzmax)
    NPPd3D(1:nzmax,n) = vertNPPd(1:nzmax)
    
    ! --------------------------------------------------------------------------
    ! Coccolithophores and Phaeocystis (if enabled)
    ! --------------------------------------------------------------------------
    if (enable_coccos) then
        aggc(1:nzmax,n)   = vertaggc(1:nzmax)
        docexc(1:nzmax,n) = vertdocexc(1:nzmax)
        respc(1:nzmax,n)  = vertrespc(1:nzmax)
        NPPc3D(1:nzmax,n) = vertNPPc(1:nzmax)
        
        aggp(1:nzmax,n)   = vertaggp(1:nzmax)
        docexp(1:nzmax,n) = vertdocexp(1:nzmax)
        respp(1:nzmax,n)  = vertrespp(1:nzmax)
        NPPp3D(1:nzmax,n) = vertNPPp(1:nzmax)
    endif
    
    ! --------------------------------------------------------------------------
    ! Calcification
    ! --------------------------------------------------------------------------
    calcdiss(1:nzmax,n) = vertcalcdiss(1:nzmax)
    calcif(1:nzmax,n)   = vertcalcif(1:nzmax)
    
    ! --------------------------------------------------------------------------
    ! Zooplankton Respiration
    ! --------------------------------------------------------------------------
    respmeso(1:nzmax,n) = vertrespmeso(1:nzmax)
    
    if (enable_3zoo2det) then
        respmacro(1:nzmax,n) = vertrespmacro(1:nzmax)
        respmicro(1:nzmax,n) = vertrespmicro(1:nzmax)
    endif

    ! --------------------------------------------------------------------------
    ! Zooplankton Mortality / Dissolution
    ! --------------------------------------------------------------------------
    mesocdis(1:nzmax,n) = vertmesocdis(1:nzmax)

    if (enable_3zoo2det) then
        microcdis(1:nzmax,n) = vertmicrocdis(1:nzmax)
        macrocdis(1:nzmax,n) = vertmacrocdis(1:nzmax)
        fastcdis(1:nzmax,n)  = vertfastcdis(1:nzmax)
    endif

    ! --------------------------------------------------------------------------
    ! Photosynthesis and Nitrogen Assimilation
    ! --------------------------------------------------------------------------
    photn(1:nzmax,n)    = vertphotn(1:nzmax)
    photd(1:nzmax,n)    = vertphotd(1:nzmax)
    Nassimn(1:nzmax,n)  = vertNassimn(1:nzmax)
    Nassimd(1:nzmax,n)  = vertNassimd(1:nzmax)

    if (enable_coccos) then
        photc(1:nzmax,n)   = vertphotc(1:nzmax)
        photp(1:nzmax,n)   = vertphotp(1:nzmax)
        Nassimc(1:nzmax,n) = vertNassimc(1:nzmax)
        Nassimp(1:nzmax,n) = vertNassimp(1:nzmax)
    endif

    ! --------------------------------------------------------------------------
    ! DOM Remineralisation
    ! --------------------------------------------------------------------------
    DONremin(1:nzmax,n) = vertDONremin(1:nzmax)
    DOCremin(1:nzmax,n) = vertDOCremin(1:nzmax)

    ! --------------------------------------------------------------------------
    ! Temperature and Photosynthesis Tracking - Phytoplankton
    ! Always active: global arrays are allocated unconditionally in recom_init.
    ! --------------------------------------------------------------------------

    TPhyCO2(1:nzmax,n)             = VTPhyCO2(1:nzmax)
    TDiaCO2(1:nzmax,n)             = VTDiaCO2(1:nzmax)
    TCphotLigLim_phyto(1:nzmax,n)  = VTCphotLigLim_phyto(1:nzmax)
    TCphot_phyto(1:nzmax,n)        = VTCphot_phyto(1:nzmax)
    TCphotLigLim_diatoms(1:nzmax,n)= VTCphotLigLim_diatoms(1:nzmax)
    TCphot_diatoms(1:nzmax,n)      = VTCphot_diatoms(1:nzmax)
    TTemp_phyto(1:nzmax,n)         = VTTemp_phyto(1:nzmax)
    TqlimitFac_phyto(1:nzmax,n)    = VTqlimitFac_phyto(1:nzmax)
    TTemp_diatoms(1:nzmax,n)       = VTTemp_diatoms(1:nzmax)
    TqlimitFac_diatoms(1:nzmax,n)  = VTqlimitFac_diatoms(1:nzmax)
    TSi_assimDia(1:nzmax,n)        = VTSi_assimDia(1:nzmax)

    ! --------------------------------------------------------------------------
    ! Temperature and Photosynthesis Tracking - Coccos/Phaeo (if enabled)
    ! --------------------------------------------------------------------------
    if (enable_coccos) then
        TTemp_cocco(1:nzmax,n)         = VTTemp_cocco(1:nzmax)
        TCoccoCO2(1:nzmax,n)           = VTCoccoCO2(1:nzmax)
        TqlimitFac_cocco(1:nzmax,n)    = VTqlimitFac_cocco(1:nzmax)
        TCphotLigLim_cocco(1:nzmax,n)  = VTCphotLigLim_cocco(1:nzmax)
        TCphot_cocco(1:nzmax,n)        = VTCphot_cocco(1:nzmax)
        TTemp_phaeo(1:nzmax,n)         = VTTemp_phaeo(1:nzmax)
        TPhaeoCO2(1:nzmax,n)           = VTPhaeoCO2(1:nzmax)
        TqlimitFac_phaeo(1:nzmax,n)    = VTqlimitFac_phaeo(1:nzmax)
        TCphotLigLim_phaeo(1:nzmax,n)  = VTCphotLigLim_phaeo(1:nzmax)
        TCphot_phaeo(1:nzmax,n)        = VTCphot_phaeo(1:nzmax)
    endif
    
end subroutine update_3d_diags

! ==============================================================================
! SUBROUTINE: deallocate_diags
! Purpose: Deallocate all diagnostic arrays
! ==============================================================================
subroutine deallocate_diags()
    use recom_locvar
    use REcoM_declarations
    implicit none

        ! --------------------------------------------------------------------------
        ! Small Phytoplankton
        ! --------------------------------------------------------------------------
        deallocate(vertNPPn, vertGPPn, vertNNAn, vertChldegn)
        deallocate(vertaggn, vertdocexn, vertrespn)
        deallocate(VTPhyCO2, VTCphotLigLim_phyto, VTCphot_phyto)
    
        ! --------------------------------------------------------------------------
        ! Diatoms
        ! --------------------------------------------------------------------------
        deallocate(vertNPPd, vertGPPd, vertNNAd, vertChldegd)
        deallocate(vertaggd, vertdocexd, vertrespd)
        deallocate(VTDiaCO2, VTCphotLigLim_diatoms, VTCphot_diatoms)

        deallocate(VTTemp_phyto, VTqlimitFac_phyto)
        deallocate(VTTemp_diatoms, VTqlimitFac_diatoms)
        deallocate(VTSi_assimDia)
        deallocate(vertcalcdiss, vertcalcif)

        ! --------------------------------------------------------------------------
        ! Zooplankton Mortality / Dissolution
        ! --------------------------------------------------------------------------
        deallocate(vertmesocdis)

        if (enable_3zoo2det) then
            deallocate(vertfastcdis)
            deallocate(vertmicrocdis,vertmacrocdis)
        endif

        ! --------------------------------------------------------------------------
        ! Photosynthesis and Nitrogen Assimilation
        ! --------------------------------------------------------------------------
        deallocate(vertphotn, vertNassimn)
        deallocate(vertphotd, vertNassimd)


        ! --------------------------------------------------------------------------
        ! DOM Remineralisation
        ! --------------------------------------------------------------------------
        deallocate(vertDONremin, vertDOCremin)

    if (enable_coccos) then
        ! --------------------------------------------------------------------------
        ! Coccolithophores and Phaeocystis (if enabled)
        ! --------------------------------------------------------------------------
        deallocate(vertNPPc, vertGPPc, vertNNAc, vertChldegc)
        deallocate(vertaggc, vertdocexc, vertrespc)
        deallocate(VTTemp_cocco, VTCoccoCO2, VTqlimitFac_cocco)
        deallocate(VTCphotLigLim_cocco, VTCphot_cocco)

        deallocate(vertNPPp, vertGPPp, vertNNAp, vertChldegp)
        deallocate(vertaggp, vertdocexp, vertrespp)
        deallocate(VTTemp_phaeo, VTPhaeoCO2, VTqlimitFac_phaeo)
        deallocate(VTCphotLigLim_phaeo, VTCphot_phaeo)

        ! --------------------------------------------------------------------------
        ! Photosynthesis and Nitrogen Assimilation
        ! --------------------------------------------------------------------------
        deallocate(vertphotc, vertNassimc)
        deallocate(vertphotp, vertNassimp)
    endif
    
    ! --------------------------------------------------------------------------
    ! Zooplankton Grazing (if enabled)
    ! --------------------------------------------------------------------------
    if (Grazing_detritus) then
        deallocate(vertgrazmeso_tot, vertgrazmeso_n, vertgrazmeso_d)
        deallocate(vertgrazmeso_det, vertrespmeso)
        
        if (enable_coccos) then
            deallocate(vertgrazmeso_c, vertgrazmeso_p)
        endif
        
        if (enable_3zoo2det) then
            deallocate(vertgrazmeso_mic, vertgrazmeso_det2)
            
            deallocate(vertgrazmacro_tot, vertgrazmacro_n, vertgrazmacro_d)
            deallocate(vertgrazmacro_mes, vertgrazmacro_det)
            deallocate(vertgrazmacro_mic, vertgrazmacro_det2)
            deallocate(vertrespmacro)
            
            if (enable_coccos) then
                deallocate(vertgrazmacro_c, vertgrazmacro_p)
            endif
            
            deallocate(vertgrazmicro_tot, vertgrazmicro_n, vertgrazmicro_d)
            deallocate(vertrespmicro)
            
            if (enable_coccos) then
                deallocate(vertgrazmicro_c, vertgrazmicro_p)
            endif
        endif
    endif

    ! Dissolution and remineralization ! R2OMIP
    deallocate(vertDISSOC, vertDISSON, vertDISSOSi, vertREMOC, vertREMOCt, vertREMON)
    
end subroutine deallocate_diags

! ==============================================================================
! SUBROUTINE: exchange_diags
! Purpose: MPI halo exchange for all active 2D/3D diagnostic fields
! ==============================================================================
subroutine exchange_diags(partit)
    use recom_glovar
    use recom_declarations
    use MOD_PARTIT
    use g_comm_auto
    implicit none
 
    type(t_partit), intent(inout), target :: partit
 
    ! Net Primary / Gross Primary Production and N assimilation
    call exchange_nod(NPPn,    partit); call exchange_nod(GPPn,    partit)
    call exchange_nod(NNAn,    partit); call exchange_nod(Chldegn, partit)
    call exchange_nod(NPPd,    partit); call exchange_nod(GPPd,    partit)
    call exchange_nod(NNAd,    partit); call exchange_nod(Chldegd, partit)
 
    if (enable_coccos) then
        call exchange_nod(NPPc,    partit); call exchange_nod(GPPc,    partit)
        call exchange_nod(NNAc,    partit); call exchange_nod(Chldegc, partit)
        call exchange_nod(NPPp,    partit); call exchange_nod(GPPp,    partit)
        call exchange_nod(NNAp,    partit); call exchange_nod(Chldegp, partit)
    endif
 
    ! Mesozooplankton
    call exchange_nod(grazmeso_tot, partit); call exchange_nod(grazmeso_n,   partit)
    call exchange_nod(grazmeso_d,   partit); call exchange_nod(grazmeso_det, partit)
 
    if (enable_coccos) then
        call exchange_nod(grazmeso_c, partit); call exchange_nod(grazmeso_p, partit)
    endif
 
    if (enable_3zoo2det) then
        call exchange_nod(grazmeso_mic,   partit); call exchange_nod(grazmeso_det2,  partit)
        call exchange_nod(grazmacro_tot,  partit); call exchange_nod(grazmacro_n,    partit)
        call exchange_nod(grazmacro_d,    partit); call exchange_nod(grazmacro_mes,  partit)
        call exchange_nod(grazmacro_det,  partit); call exchange_nod(grazmacro_mic,  partit)
        call exchange_nod(grazmacro_det2, partit)
        call exchange_nod(grazmicro_tot,  partit); call exchange_nod(grazmicro_n,    partit)
        call exchange_nod(grazmicro_d,    partit)
 
        if (enable_coccos) then
            call exchange_nod(grazmacro_c, partit); call exchange_nod(grazmacro_p, partit)
            call exchange_nod(grazmicro_c, partit); call exchange_nod(grazmicro_p, partit)
        endif
    endif
 
end subroutine exchange_diags

end module recom_diags_management
