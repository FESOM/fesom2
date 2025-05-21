!=======================================================================
!
! Contains Icepack component driver routines common to all drivers.
!
! Authors: Lorenzo Zampieri ( lorenzo.zampieri@awi.de )
! Adapted for Icepack 1.4.1 by F. Kauker (frank.kauker@awi.de)
!
!=======================================================================

submodule (icedrv_main) icedrv_step

  use icedrv_kinds
  use icedrv_system, only: icedrv_system_abort
  use icepack_intfc, only: &
       icepack_warnings_flush, icepack_warnings_aborted, icepack_query_tracer_flags,      &
       icepack_query_tracer_indices, icepack_query_tracer_sizes, icepack_query_parameters
      
  implicit none

  !=======================================================================

contains

  !=======================================================================
  !
  ! Scales radiation fields computed on the previous time step.
  !

  subroutine prep_radiation ()

    ! column package includes
    use icepack_intfc, only: icepack_prep_radiation

    implicit none
    
    ! local variables    
    integer (kind=int_kind) :: i               ! horizontal indices
    
    character(len=*), parameter :: subname='(prep_radiation)'
    
    !-----------------------------------------------------------------
    ! Compute netsw scaling factor (new netsw / old netsw)
    !-----------------------------------------------------------------
    
    do i = 1, nx
 
       alvdr_init(i) = alvdr_ai(i)
       alvdf_init(i) = alvdf_ai(i)
       alidr_init(i) = alidr_ai(i)
       alidf_init(i) = alidf_ai(i)
 
       call icepack_prep_radiation( &
            ncat=ncat,                      nilyr=nilyr,                    &
            nslyr=nslyr,                                                    &
            aice=aice(i),                   aicen=aicen(i,:),               &
            swvdr=swvdr(i),                 swvdf=swvdf(i),                 &
            swidr=swidr(i),                 swidf=swidf(i),                 &
            alvdr_ai=alvdr_ai(i),           alvdf_ai=alvdf_ai(i),           &
            alidr_ai=alidr_ai(i),           alidf_ai=alidf_ai(i),           & 
            scale_factor=scale_factor(i),                                   &
            fswsfcn=fswsfcn(i,:),           fswintn=fswintn(i,:),           &
            fswthrun=fswthrun(i,:),         fswthrun_vdr=fswthrun_vdr(i,:), &
            fswthrun_vdf=fswthrun_vdf(i,:), fswthrun_idr=fswthrun_idr(i,:), &
            fswthrun_idf=fswthrun_idf(i,:), fswpenln=fswpenln(i,:,:),       &
            Sswabsn=Sswabsn(i,:,:),         Iswabsn=Iswabsn(i,:,:)) 
    enddo               ! i

    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)

  end subroutine prep_radiation

  !=======================================================================
  !
  ! Driver for updating ice and snow internal temperatures and
  ! computing thermodynamic growth rates and coupler fluxes.
  !

  subroutine step_therm1 (dt)

    ! column packge includes
    use icepack_intfc, only: icepack_step_therm1
    
    implicit none
    
    logical (kind=log_kind) :: prescribed_ice   ! if .true., use prescribed ice instead of computed
    
    real (kind=dbl_kind), intent(in) :: dt      ! time step
    
    ! local variables
    
    integer (kind=int_kind) :: &
         i,                                     & ! horizontal indices
         n,                                     & ! thickness category index
         k, kk                                    ! indices for aerosols
    
    integer (kind=int_kind) :: &
         ntrcr, nt_apnd, nt_hpnd, nt_ipnd, nt_alvl, nt_vlvl, nt_Tsfc,   &
         nt_iage, nt_FY, nt_qice, nt_sice, nt_qsno, nt_aero, nt_isosno, &
         nt_isoice, nt_rsnw, nt_smice, nt_smliq
    
    logical (kind=log_kind) :: &
         tr_iage, tr_FY, tr_aero, tr_iso, calc_Tsfc, snwgrain, tr_pond, tr_pond_lvl, tr_pond_topo
    
    real (kind=dbl_kind), dimension(n_aero,2,ncat) :: aerosno,  aeroice ! kg/m^2
        
    real (kind=dbl_kind), dimension(n_iso,ncat) :: &
         isosno, isoice    ! kg/m^2

    real (kind=dbl_kind), dimension(nslyr,ncat) :: &
         rsnwn, smicen, smliqn
    
    real (kind=dbl_kind) :: puny
    
    character(len=*), parameter :: subname='(step_therm1)'
    
    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------
    
    call icepack_query_parameters( &
         puny_out=puny, calc_Tsfc_out=calc_Tsfc,  snwgrain_out=snwgrain)

    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
    
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_aero_out=tr_aero, tr_pond_out=tr_pond,    &
         tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)

    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    call icepack_query_tracer_indices( &
         nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, nt_alvl_out=nt_alvl, &
         nt_vlvl_out=nt_vlvl, nt_Tsfc_out=nt_Tsfc, nt_iage_out=nt_iage, nt_FY_out=nt_FY,     &
         nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_aero_out=nt_aero, nt_qsno_out=nt_qsno, &
         nt_rsnw_out=nt_rsnw, nt_smice_out=nt_smice, nt_smliq_out=nt_smliq,                  &
         nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice )

    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    !-----------------------------------------------------------------
    
    prescribed_ice = .false.
    aerosno(:,:,:) = c0
    aeroice(:,:,:) = c0
    isosno (:,:)   = c0
    isoice (:,:)   = c0
    rsnwn  (:,:)   = c0
    smicen (:,:)   = c0
    smliqn (:,:)   = c0
    
    do i = 1, nx
    
       !-----------------------------------------------------------------
       ! Save the ice area passed to the coupler (so that history fields
       !  can be made consistent with coupler fields).
       ! Save the initial ice area and volume in each category.
       !-----------------------------------------------------------------
    
       aice_init (i) = aice (i)
    
       do n = 1, ncat
          aicen_init(i,n) = aicen(i,n)
          vicen_init(i,n) = vicen(i,n)
          vsnon_init(i,n) = vsnon(i,n)
       enddo
       
    enddo ! i

    do i = 1, nx
       if(mype==7.and.i==5) then
          l_print_point=.true.
       else
          l_print_point=.false.
       endif
       if (tr_aero) then
          ! trcrn(nt_aero) has units kg/m^3
          do n=1,ncat
             do k=1,n_aero
                aerosno (k,:,n) = &
                     trcrn(i,nt_aero+(k-1)*4  :nt_aero+(k-1)*4+1,n) &
                     * vsnon_init(i,n)
                aeroice (k,:,n) = &
                     trcrn(i,nt_aero+(k-1)*4+2:nt_aero+(k-1)*4+3,n) &
                     * vicen_init(i,n)
             enddo
          enddo
       endif ! tr_aero

       if (tr_iso) then
          ! trcrn(nt_isosno/ice) has units kg/m^3
          do n = 1, ncat
             do k = 1, n_iso
                isosno(k,n) = trcrn(i,nt_isosno+k-1,n) * vsnon_init(i,n)
                isoice(k,n) = trcrn(i,nt_isoice+k-1,n) * vicen_init(i,n)
             enddo
          enddo
       endif ! tr_iso

       if (snwgrain) then
          do n = 1, ncat
             do k = 1, nslyr
                rsnwn (k,n) = trcrn(i,nt_rsnw +k-1,n)
                smicen(k,n) = trcrn(i,nt_smice+k-1,n)
                smliqn(k,n) = trcrn(i,nt_smliq+k-1,n)
             enddo
          enddo
       endif ! snwgrain

       call icepack_step_therm1( &            
            l_print_point=l_print_point,                          &
            dt=dt,            ncat=ncat,                          &
            nilyr=nilyr,      nslyr=nslyr,                        &
            aicen_init = aicen_init(i,:),                         &
            vicen_init = vicen_init(i,:),                         &
            vsnon_init = vsnon_init(i,:),                         &
            aice = aice(i),   aicen = aicen(i,:),                 &
            vice = vice(i),   vicen = vicen(i,:),                 &
            vsno = vsno(i),   vsnon = vsnon(i,:),                 &
            uvel = uvel(i),   vvel  = vvel(i),                    &
            Tsfc = trcrn(i,nt_Tsfc,:),                            &
            zqsn = trcrn(i,nt_qsno:nt_qsno+nslyr-1,:),            &
            zqin = trcrn(i,nt_qice:nt_qice+nilyr-1,:),            &
            zSin = trcrn(i,nt_sice:nt_sice+nilyr-1,:),            &
	    alvl = trcrn(i,nt_alvl,:),                            &
            vlvl = trcrn(i,nt_vlvl,:),                            &
            apnd = trcrn(i,nt_apnd,:),                            &
            hpnd = trcrn(i,nt_hpnd,:),                            &
            ipnd = trcrn(i,nt_ipnd,:),                            &
            iage = trcrn(i,nt_iage,:),                            &
            FY   = trcrn(i,nt_FY,:),                              &
            rsnwn  = rsnwn (:,:),                                 &
            smicen = smicen(:,:),                                 &
            smliqn = smliqn(:,:),                                 &
            aerosno = aerosno(:,:,:),                             &
            aeroice = aeroice(:,:,:),                             &
            isosno  = isosno(:,:),                                &
            isoice  = isoice(:,:),                                &
            uatm = uatm(i), vatm = vatm(i),                       &
            wind = wind(i), zlvl = zlvl_v,                        &
            Qa   = Qa(i),   rhoa = rhoa(i),                       &
            Qa_iso = Qa_iso(i,:),                                 &
            Tair = T_air(i), Tref = Tref(i),                      &
            Qref = Qref(i), Uref = Uref(i),                       &
            Qref_iso = Qref_iso(i,:),                             &
            Cdn_atm_ratio = Cdn_atm_ratio(i),                     &
            Cdn_ocn       = Cdn_ocn(i),                           &
            Cdn_ocn_skin  = Cdn_ocn_skin(i),                      &
            Cdn_ocn_floe  = Cdn_ocn_floe(i),                      &
            Cdn_ocn_keel  = Cdn_ocn_keel(i),                      &
            Cdn_atm       = Cdn_atm(i),                           &
            Cdn_atm_skin  = Cdn_atm_skin(i),                      &
            Cdn_atm_floe  = Cdn_atm_floe(i),                      &
            Cdn_atm_pond  = Cdn_atm_pond(i),                      &
            Cdn_atm_rdg   = Cdn_atm_rdg(i),                       &
            hfreebd  = hfreebd(i),    hkeel     = hkeel(i),       &
            hdraft   = hdraft(i),     hridge    = hridge(i),      &
            distrdg  = distrdg(i),    dkeel     = dkeel(i),       &
            lfloe    = lfloe(i),      dfloe     = dfloe(i),       &
            strax    = strax(i),      stray     = stray(i),       &
            strairxT = strairxT(i),   strairyT  = strairyT(i),    &
            potT     = potT(i),       sst       = sst(i),         &
            sss      = sss(i),        Tf        = Tf(i),          &
            strocnxT = strocnxT(i),   strocnyT  = strocnyT(i),    &
            fbot     = fbot(i),       frzmlt    = frzmlt(i),      &
            Tbot     = Tbot(i),       Tsnice    = Tsnice(i),      &
            rside    = rside(i),      fside     = fside(i),       &
            wlat     = wlat(i),                                   &
            fsnow    = fsnow(i),      frain     = frain(i),       &
            fpond    = fpond(i),      fsloss    = fsloss(i),      &
            fsurf    = fsurf(i),      fsurfn    = fsurfn(i,:),    &
            fcondtop = fcondtop(i),   fcondtopn = fcondtopn(i,:), &
            fcondbot = fcondbot(i),   fcondbotn = fcondbotn(i,:), &
            fswsfcn  = fswsfcn(i,:),  fswintn   = fswintn(i,:),   &
            fswthrun = fswthrun(i,:),                             &
            fswthrun_vdr = fswthrun_vdr(i,:),                     &
            fswthrun_vdf = fswthrun_vdf(i,:),                     &
            fswthrun_idr = fswthrun_idr(i,:),                     &
            fswthrun_idf = fswthrun_idf(i,:),                     &
            fswabs    = fswabs(i),                                &
            flwout   = flwout(i),     flw       = flw(i),         &
            fsens    = fsens(i),      fsensn    = fsensn(i,:),    &
            flat     = flat(i),       flatn     = flatn(i,:),     &
            fresh    = fresh(i),      fsalt     = fsalt(i),       &
            fhocn    = fhocn(i),                                  &
            fswthru   = fswthru(i),                               &
            fswthru_vdr= fswthru_vdr(i),                          &
            fswthru_vdf= fswthru_vdf(i),                          &
            fswthru_idr= fswthru_idr(i),                          &
            fswthru_idf= fswthru_idf(i),                          &
            flatn_f  = flatn_f(i,:),  fsensn_f  = fsensn_f(i,:),  &
            fsurfn_f = fsurfn_f(i,:),                             &
            fcondtopn_f = fcondtopn_f(i,:),                       &
            faero_atm   = faero_atm(i,1:n_aero),                  &
            faero_ocn   = faero_ocn(i,1:n_aero),                  &
            fiso_atm    = fiso_atm   (i,:),                       &
            fiso_ocn    = fiso_ocn   (i,:),                       &
            fiso_evap   = fiso_evap  (i,:),                       &
            HDO_ocn     = HDO_ocn (i),                            &
            H2_16O_ocn  = H2_16O_ocn (i),                         &
            H2_18O_ocn  = H2_18O_ocn (i),                         &
            Sswabsn  = Sswabsn(i,:,:),Iswabsn   = Iswabsn(i,:,:), &
            evap = evap(i), evaps = evaps(i), evapi = evapi(i),   &
            dhsn     = dhsn(i,:),     ffracn    = ffracn(i,:),    &
            meltt    = meltt(i),      melttn    = melttn(i,:),    &
            meltb    = meltb(i),      meltbn    = meltbn(i,:),    &
            melts    = melts(i),      meltsn    = meltsn(i,:),    &
            congel   = congel(i),     congeln   = congeln(i,:),   &
            snoice   = snoice(i),     snoicen   = snoicen(i,:),   &
            dsnow    = dsnow(i),      dsnown    = dsnown(i,:),    &
            meltsliqn= meltsliqn(i,:),                            &
            lmask_n  = lmask_n(i),    lmask_s   = lmask_s(i),     &
            mlt_onset=mlt_onset(i),   frz_onset = frz_onset(i),   &
            yday = yday,  prescribed_ice = prescribed_ice)
    
       if (tr_aero) then
          do n = 1, ncat
             if (vicen(i,n) > puny) &
                  aeroice(:,:,n) = aeroice(:,:,n)/vicen(i,n)
             if (vsnon(i,n) > puny) &
                  aerosno(:,:,n) = aerosno(:,:,n)/vsnon(i,n)
             do k = 1, n_aero
                do kk = 1, 2
                   trcrn(i,nt_aero+(k-1)*4+kk-1,n)=aerosno(k,kk,n)
                   trcrn(i,nt_aero+(k-1)*4+kk+1,n)=aeroice(k,kk,n)
                enddo
             enddo
          enddo
       endif ! tr_aero
            
       if (tr_iso) then
          do n = 1, ncat
             if (vicen(i,n) > puny) isoice(:,n) = isoice(:,n)/vicen(i,n)
             if (vsnon(i,n) > puny) isosno(:,n) = isosno(:,n)/vsnon(i,n)
             do k = 1, n_iso
                trcrn(i,nt_isosno+k-1,n) = isosno(k,n)
                trcrn(i,nt_isoice+k-1,n) = isoice(k,n)
             enddo
          enddo
       endif ! tr_iso
           
       if (snwgrain) then
          do n = 1, ncat
             do k = 1, nslyr
                trcrn(i,nt_rsnw +k-1,n) = rsnwn (k,n)
                trcrn(i,nt_smice+k-1,n) = smicen(k,n)
                trcrn(i,nt_smliq+k-1,n) = smliqn(k,n)
             enddo
          enddo
       endif ! snwgrain
    enddo ! i

    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)
      
  end subroutine step_therm1

  !=======================================================================
  ! Driver for thermodynamic changes not needed for coupling:
  ! transport in thickness space, lateral growth and melting.
  !
  ! authors: William H. Lipscomb, LANL
  !          Elizabeth C. Hunke, LANL

  subroutine step_therm2 (dt)

    ! column package_includes
    use icepack_intfc, only: icepack_step_therm2

    implicit none
    
    real (kind=dbl_kind), intent(in) :: dt      ! time step
    
    ! local variables
    
    integer (kind=int_kind) :: &
         i,                                   & ! horizontal index
         ntrcr, nbtrcr
    
    logical (kind=log_kind) :: &
         tr_fsd,                              & ! floe size distribution tracers
         update_ocn_f
          
    character(len=*), parameter :: subname='(step_therm2)'
    
    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------
    
    call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
    call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
    call icepack_query_parameters(update_ocn_f_out=update_ocn_f)

    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    !-----------------------------------------------------------------
    
    do i = 1, nx
    
       ! wave_sig_ht - compute here to pass to add new ice
       if (tr_fsd) &
            wave_sig_ht(i) = c4*SQRT(SUM(wave_spectrum(i,:)*dwavefreq(:)))
    
       call icepack_step_therm2( &
            dt=dt, ncat=ncat,                            &
            nltrcr=nltrcr, nilyr=nilyr, nslyr=nslyr,     &
            hin_max=hin_max(:), nblyr=nblyr,             &   
            aicen=aicen(i,:),                            &
            vicen=vicen(i,:),                            &
            vsnon=vsnon(i,:),                            &
            aicen_init=aicen_init(i,:),                  &
            vicen_init=vicen_init(i,:),                  &
            trcrn=trcrn(i,1:ntrcr,:),                    &
            aice0=aice0(i),                              &
            aice =aice(i),                               &
            trcr_depend=trcr_depend(1:ntrcr),            &
            trcr_base=trcr_base(1:ntrcr,:),              &
            n_trcr_strata=n_trcr_strata(1:ntrcr),        &
            nt_strata=nt_strata(1:ntrcr,:),              &
            Tf=Tf(i), sss=sss(i),                        &
            salinz=salinz(i,:), fside=fside(i),          &
            rside=rside(i),   meltl=meltl(i),            &
            frzmlt=frzmlt(i), frazil=frazil(i),          &
            frain=frain(i),   fpond=fpond(i),            &
            fresh=fresh(i),   fsalt=fsalt(i),            &
            fhocn=fhocn(i),   update_ocn_f=update_ocn_f, &
            bgrid=bgrid,      cgrid=cgrid,               &
            igrid=igrid,      faero_ocn=faero_ocn(i,:),  &
            first_ice=first_ice(i,:),                    &
            fzsal=fzsal(i),                              &
            flux_bio=flux_bio(i,1:nbtrcr),               &
            ocean_bio=ocean_bio(i,1:nbtrcr),             &
            frazil_diag=frazil_diag(i),                  &
            frz_onset=frz_onset(i),                      &
            yday=yday,                                   &
            nfsd=nfsd,   wave_sig_ht=wave_sig_ht(i),     &
            wave_spectrum=wave_spectrum(i,:),            &
            wavefreq=wavefreq(:),                        &
            dwavefreq=dwavefreq(:),                      &
            d_afsd_latg=d_afsd_latg(i,:),                &
            d_afsd_newi=d_afsd_newi(i,:),                &
            d_afsd_latm=d_afsd_latm(i,:),                &
            d_afsd_weld=d_afsd_weld(i,:),                &
            floe_rad_c=floe_rad_c(:),                    &
            floe_binwidth=floe_binwidth(:))

    enddo ! i

    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)
         
  end subroutine step_therm2

  !=======================================================================
  !
  ! finalize thermo updates
  !

  subroutine update_state (dt, daidt, dvidt, dagedt, offset)

    ! column package includes
    use icepack_intfc, only: icepack_aggregate
    
    implicit none
    
    real (kind=dbl_kind), intent(in) :: &
         dt,                           & ! time step
         offset                          ! d(age)/dt time offset = dt for thermo, 0 for dyn
    
    real (kind=dbl_kind), dimension(:), intent(inout) :: &
         daidt,                        & ! change in ice area per time step
         dvidt,                        & ! change in ice volume per time step
         dagedt                          ! change in ice age per time step
    
    integer (kind=int_kind) :: & 
         i,                            & ! horizontal indices
         ntrcr, nt_iage
    
    logical (kind=log_kind) :: tr_iage   ! ice age tracer
    
    character(len=*), parameter :: subname='(update_state)'
    
    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------
    
    call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
    call icepack_warnings_flush(ice_stderr)

    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    call icepack_query_tracer_indices(nt_iage_out=nt_iage)
    call icepack_warnings_flush(ice_stderr)
    
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    call icepack_query_tracer_flags(tr_iage_out=tr_iage)

    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    do i = 1, nx
    
       !-----------------------------------------------------------------
       ! Aggregate the updated state variables (includes ghost cells). 
       !----------------------------------------------------------------- 
     
      call icepack_aggregate( &
            ncat=ncat,                    &
            trcrn=trcrn(i,1:ntrcr,:),     &
            aicen=aicen(i,:),             &
            vicen=vicen(i,:),             &
            vsnon=vsnon(i,:),             &
            trcr=trcr(i,1:ntrcr),         &
            aice=aice(i),                 &
            vice=vice(i),                 &
            vsno =vsno(i),                &
            aice0=aice0(i),               &
            ntrcr=ntrcr,                  &
            trcr_depend=trcr_depend(1:ntrcr),     &
            trcr_base=trcr_base(1:ntrcr,:),       &
            n_trcr_strata=n_trcr_strata(1:ntrcr), &
            nt_strata=nt_strata(1:ntrcr,:),       &
            Tf=Tf(i)          )

       !-----------------------------------------------------------------
       ! Compute thermodynamic area and volume tendencies.
       !-----------------------------------------------------------------
    
       daidt(i) = (aice(i) - daidt(i)) / dt
       dvidt(i) = (vice(i) - dvidt(i)) / dt
       if (tr_iage) then
          if (offset > c0) then                 ! thermo
             if (trcr(i,nt_iage) > c0) &
                  dagedt(i) = (trcr(i,nt_iage) &
                  - dagedt(i) - offset) / dt
          else                                  ! dynamics
             dagedt(i) = (trcr(i,nt_iage) &
                  - dagedt(i)) / dt
          endif
       endif
       
    enddo ! i
    
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)
    
  end subroutine update_state

!=======================================================================
!
! Run one time step of wave-fracturing the floe size distribution
!

  subroutine step_dyn_wave (dt)
 
    ! column package includes
    use icepack_intfc, only: icepack_step_wavefracture
    
    implicit none
    
    real (kind=dbl_kind), intent(in) :: dt      ! time step
    
    ! local variables
    
    integer (kind=int_kind) :: i, j, ntrcr, nbtrcr
    
    character (len=char_len) :: wave_spec_type
    
    character(len=*), parameter :: subname = '(step_dyn_wave)'
    
    call icepack_query_parameters(wave_spec_type_out=wave_spec_type)
    
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    do i = 1, nx
       d_afsd_wave(i,:) = c0
       call icepack_step_wavefracture (wave_spec_type=wave_spec_type, &
            dt=dt, ncat=ncat, nfsd=nfsd, nfreq=nfreq, &
            aice          = aice         (i),      &
            vice          = vice         (i),      &
            aicen         = aicen        (i,:),    &
            floe_rad_l    = floe_rad_l     (:),    &
            floe_rad_c    = floe_rad_c     (:),    &
            wave_spectrum = wave_spectrum(i,:),    &
            wavefreq      = wavefreq       (:),    &
            dwavefreq     = dwavefreq      (:),    &
            trcrn         = trcrn        (i,:,:),  &
            d_afsd_wave   = d_afsd_wave  (i,:))
    end do ! i
    
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
  end subroutine step_dyn_wave

  !=======================================================================
  !
  ! Run one time step of ridging.
  !

  subroutine step_dyn_ridge (dt, ndtd)

    ! column package includes
    use icepack_intfc, only: icepack_step_ridge
    
    implicit none

    real (kind=dbl_kind), intent(in) :: dt         ! time step
    
    integer (kind=int_kind), intent(in) :: ndtd    ! number of dynamics subcycles
    
    ! local variables
    
    integer (kind=int_kind) :: i, ntrcr, nbtrcr
    
    character(len=*), parameter :: subname='(step_dyn_ridge)'

    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------
 
    call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
 
    !-----------------------------------------------------------------
    ! Ridging
    !-----------------------------------------------------------------

    do i = 1, nx
       call icepack_step_ridge ( &
            dt=dt,                    ndtd=ndtd,                &
            nilyr=nilyr,              nslyr=nslyr,              &
            nblyr=nblyr,                                        &
            ncat=ncat,                hin_max=hin_max(:),       &
            rdg_conv=rdg_conv(i),     rdg_shear=rdg_shear(i),   &
            aicen=aicen(i,:),                                   &
            trcrn=trcrn(i,1:ntrcr,:),                           &
            vicen=vicen(i,:),         vsnon=vsnon(i,:),         &
            aice0=aice0(i),                                     &
            trcr_depend=trcr_depend(1:ntrcr),                   &
            trcr_base=trcr_base(1:ntrcr,:),                     &
            n_trcr_strata=n_trcr_strata(1:ntrcr),               &
            nt_strata=nt_strata(1:ntrcr,:),                     &
            dardg1dt=dardg1dt(i),     dardg2dt=dardg2dt(i),     &
            dvirdgdt=dvirdgdt(i),     opening=opening(i),       &
            fpond=fpond(i),                                     &
            fresh=fresh(i),           fhocn=fhocn(i),           &
            n_aero=n_aero,                                      &
            faero_ocn=faero_ocn(i,:), fiso_ocn=fiso_ocn(i,:),   &
            aparticn=aparticn(i,:),   krdgn=krdgn(i,:),         &
            aredistn=aredistn(i,:),   vredistn=vredistn(i,:),   &
            dardg1ndt=dardg1ndt(i,:), dardg2ndt=dardg2ndt(i,:), &
            dvirdgndt=dvirdgndt(i,:),                           &
            araftn=araftn(i,:),       vraftn=vraftn(i,:),       &
            aice=aice(i),             fsalt=fsalt(i),           &
            first_ice=first_ice(i,:),                           &
            flux_bio=flux_bio(i,1:nbtrcr), Tf = Tf(i))
    enddo
    
    call cut_off_icepack
  
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)
          
  end subroutine step_dyn_ridge

  !=======================================================================
  !
  ! Updates snow tracers
  !
  ! authors: Elizabeth C. Hunke, LANL
  !          Nicole Jeffery, LANL
  
  subroutine step_snow (dt)
    
    use icepack_intfc, only: icepack_step_snow
    
    real (kind=dbl_kind), intent(in) :: dt        ! time step
    
    ! local variables
    
    integer (kind=int_kind) :: &
         nt_smice, nt_smliq, nt_rsnw, nt_Tsfc, nt_qice, nt_sice, nt_qsno, &
         nt_alvl, nt_vlvl, nt_rhos
    
    integer (kind=int_kind) :: i, n

    character(len=*), parameter :: subname='(step_snow)'

    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------
    
    call icepack_query_tracer_indices( &
         nt_smice_out=nt_smice, nt_smliq_out=nt_smliq, &
         nt_rsnw_out=nt_rsnw, nt_Tsfc_out=nt_Tsfc, &
         nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_qsno_out=nt_qsno, &
         nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_rhos_out=nt_rhos)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    !-----------------------------------------------------------------
    ! Snow redistribution and metamorphosis
    !-----------------------------------------------------------------

    do i = 1, nx
       call icepack_step_snow ( &
            dt,                 nilyr,              &
            nslyr,              ncat,               &
            wind(i),            aice(i),            &
            aicen(i,:),         vicen(i,:),         &
            vsnon(i,:),         trcrn(i,nt_Tsfc,:), &
            trcrn(i,nt_qice,:),    & ! top layer only
            trcrn(i,nt_sice,:),    & ! top layer only
            trcrn(i,nt_qsno:nt_qsno+nslyr-1,:),     &
            trcrn(i,nt_alvl,:), trcrn(i,nt_vlvl,:), &
            trcrn(i,nt_smice:nt_smice+nslyr-1,:),   &
            trcrn(i,nt_smliq:nt_smliq+nslyr-1,:),   &
            trcrn(i,nt_rsnw:nt_rsnw+nslyr-1,:),     &
            trcrn(i,nt_rhos:nt_rhos+nslyr-1,:),     &
            fresh(i),           fhocn(i),           &
            fsloss(i),          fsnow(i))
    enddo

  end subroutine step_snow
    
  !=======================================================================
  !
  ! Computes radiation fields
  !
  
  subroutine step_radiation (dt)

    ! column package includes
    use icepack_intfc, only: icepack_step_radiation

    implicit none
    
    real (kind=dbl_kind), intent(in) :: dt                 ! time step
    
    ! local variables
    
    integer (kind=int_kind) :: i, n, k
    
    integer (kind=int_kind) :: &
         max_aero, max_algae, nt_Tsfc, nt_alvl, &
         nt_apnd, nt_hpnd, nt_ipnd, nt_aero, nlt_chl_sw, &
         ntrcr, nbtrcr_sw, nt_fbri, nt_rsnw
    
    integer (kind=int_kind), dimension(:), allocatable :: &
         nlt_zaero_sw, nt_zaero, nt_bgc_N
    
    logical (kind=log_kind) :: &
         tr_bgc_N, tr_zaero, tr_brine, snwgrain
    
    real (kind=dbl_kind), dimension(ncat) :: fbri          ! brine height to ice thickness
    
    real(kind= dbl_kind), dimension(:,:), allocatable :: &
             rsnow,                                      & ! snow grain radius
             ztrcr_sw                                      ! BGC tracers affecting radiation
    
    logical (kind=log_kind) :: &
         l_print_point,                                  & ! flag for printing debugging information
         dEdd_algae,                                     & ! from icepack
         modal_aero                                        ! from icepack
    
    character(len=*), parameter :: subname='(step_radiation)'
    
    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------
    
    call icepack_query_tracer_sizes( &
         max_aero_out=max_aero, max_algae_out=max_algae)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)

    allocate(nlt_zaero_sw(max_aero))
    allocate(nt_zaero(max_aero))
    allocate(nt_bgc_N(max_algae))
    call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_sw_out=nbtrcr_sw)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    call icepack_query_tracer_flags( &
         tr_brine_out=tr_brine, tr_bgc_N_out=tr_bgc_N, tr_zaero_out=tr_zaero)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__,line= __LINE__)
    
    call icepack_query_tracer_indices( &
         nt_Tsfc_out=nt_Tsfc, nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
         nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero, &
         nlt_chl_sw_out=nlt_chl_sw, nlt_zaero_sw_out=nlt_zaero_sw, &
         nt_rsnw_out=nt_rsnw, &
         nt_fbri_out=nt_fbri, nt_zaero_out=nt_zaero, nt_bgc_N_out=nt_bgc_N)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    call icepack_query_parameters( &
         dEdd_algae_out=dEdd_algae, modal_aero_out=modal_aero, snwgrain_out=snwgrain )
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    !-----------------------------------------------------------------

    allocate(rsnow(nslyr,ncat))
    allocate(ztrcr_sw(nbtrcr_sw,ncat))
    
    l_print_point = .false.
    
    do i = 1, nx
       fbri(:) = c0
       rsnow(:,:) = c0
       ztrcr_sw(:,:) = c0
       do n = 1, ncat
          if (tr_brine) fbri(n) = trcrn(i,nt_fbri,n)
          if (snwgrain) rsnow(:,n) = trcrn(i,nt_rsnw:nt_rsnw+nslyr-1,n)
       enddo
    
       call icepack_step_radiation( &
            dt=dt,                     ncat=ncat,                &
            nblyr=nblyr,               nilyr=nilyr,              &
            nslyr=nslyr,               dEdd_algae=dEdd_algae,    &
            modal_aero=modal_aero,                               &
            sw_grid=swgrid(:),         i_grid=igrid(:),          &
            fbri=fbri(:),              aicen=aicen(i,:),         &
            vicen=vicen(i,:),          vsnon=vsnon(i,:),         &
            Tsfcn=trcrn(i,nt_Tsfc,:),                            &
            alvln=trcrn(i,nt_alvl,:),                            &
            apndn=trcrn(i,nt_apnd,:),                            &
            hpndn=trcrn(i,nt_hpnd,:),                            &
            ipndn=trcrn(i,nt_ipnd,:),                            &
            aeron=trcrn(i,nt_aero:nt_aero+4*n_aero-1,:),         &
            bgcNn=trcrn(i,nt_bgc_N(1):nt_bgc_N(1)+n_algae*(nblyr+3)-1,:), &
            zaeron=trcrn(i,nt_zaero(1):nt_zaero(1)+n_zaero*(nblyr+3)-1,:), &
            trcrn_bgcsw=ztrcr_sw,                                &
            TLAT=lat_val(i),           TLON=lon_val(i),          &
            calendar_type=calendar_type,                         &
            days_per_year=days_per_year, sec=sec,                &
            nextsw_cday=nextsw_cday,   yday=yday,                &
            swvdr=swvdr(i),            swvdf=swvdf(i),           &
            swidr=swidr(i),            swidf=swidf(i),           &
            coszen=cos_zen(i),         fsnow=fsnow(i),           &
            alvdrn=alvdrn(i,:),        alvdfn=alvdfn(i,:),       &
            alidrn=alidrn(i,:),        alidfn=alidfn(i,:),       &
            fswsfcn=fswsfcn(i,:),      fswintn=fswintn(i,:),     &
            fswthrun=fswthrun(i,:),                              &
            fswthrun_vdr=fswthrun_vdr(i,:),                      &
            fswthrun_vdf=fswthrun_vdf(i,:),                      &
            fswthrun_idr=fswthrun_idr(i,:),                      &
            fswthrun_idf=fswthrun_idf(i,:),                      &
            fswpenln=fswpenln(i,:,:),                            &
            Sswabsn=Sswabsn(i,:,:),    Iswabsn=Iswabsn(i,:,:),   &
            albicen=albicen(i,:),      albsnon=albsnon(i,:),     &
            albpndn=albpndn(i,:),      apeffn=apeffn(i,:),       &
            snowfracn=snowfracn(i,:),                            &
            dhsn=dhsn(i,:),            ffracn=ffracn(i,:),       &
            rsnow=rsnow(:,:),                                    &
            l_print_point=l_print_point)
    
       if (dEdd_algae .and. (tr_zaero .or. tr_bgc_N)) then
          do n = 1, ncat
             do k = 1, nbtrcr_sw
                trcrn_sw(i,k,n) = ztrcr_sw(k,n)
             enddo
          enddo
       endif
       
    enddo ! i
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)
    
        
    deallocate(rsnow)
    deallocate(ztrcr_sw)
    deallocate(nlt_zaero_sw)
    deallocate(nt_zaero)
    deallocate(nt_bgc_N)

  end subroutine step_radiation

  !=======================================================================
  !
  ! Ocean mixed layer calculation (internal to sea ice model).
  ! Allows heat storage in ocean for uncoupled runs.
  !

  subroutine ocean_mixed_layer (ice, dt)
          
    use icepack_intfc, only: icepack_atm_boundary ! icepack_ocn_mixed_layer removed frank.kauker@awi.de
    use MOD_ICE
    
    implicit none

    type(t_ice), target, intent(inout) :: ice
          
    real (kind=dbl_kind), intent(in) :: dt      ! time step
    
    ! local variables
    
    integer (kind=int_kind) :: i                ! horizontal indices
    
    real (kind=dbl_kind) :: albocn
    
    real (kind=dbl_kind), dimension(nx) :: &
         delt,                                & ! potential temperature difference   (K)
         delq,                                & ! specific humidity difference   (kg/kg)
         shcoef,                              & ! transfer coefficient for sensible heat
         lhcoef                                 ! transfer coefficient for latent heat
    
    character(len=*), parameter :: subname='(ocean_mixed_layer)'
    character (len=100) :: value
    character (len=*), parameter :: fmt = "(f6.2)" 
    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------
    
    call icepack_query_parameters(albocn_out=albocn)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)
    
    !----------------------------------------------------------------- 
    ! Compute boundary layer quantities
    !-----------------------------------------------------------------

    do i = 1, nx
       
       call icepack_atm_boundary( &
            sfctype = 'ocn',          &
            Tsf     = sst(i),         &
            potT    = potT(i),        &
            uatm    = uatm(i),        &   
            vatm    = vatm(i),        &   
            wind    = wind(i),        &   
            zlvl    = zlvl_v,         &   
            zlvs    = zlvl_s,         &   
            Qa      = Qa(i),          &     
            rhoa    = rhoa(i),        &
            strx    = strairx_ocn(i), & 
            stry    = strairy_ocn(i), & 
            Tref    = Tref_ocn(i),    & 
            Qref    = Qref_ocn(i),    & 
            delt    = delt(i),        &    
            delq    = delq(i),        &
            lhcoef  = lhcoef(i),      &
            shcoef  = shcoef(i),      &
            Cdn_atm = Cdn_atm(i),     & 
            Cdn_atm_ratio_n = Cdn_atm_ratio(i))    
    enddo ! i
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)
    
    !-----------------------------------------------------------------
    ! Ocean albedo
    ! For now, assume albedo = albocn in each spectral band.
    !-----------------------------------------------------------------
    
    alvdr_ocn(:) = albocn
    alidr_ocn(:) = albocn
    alvdf_ocn(:) = albocn
    alidf_ocn(:) = albocn
    
    !-----------------------------------------------------------------
    ! Compute ocean fluxes and update SST
    !-----------------------------------------------------------------
    
    do i = 1, nx
       call ocn_mixed_layer_icepack( &
            ice,                                             &
            alvdr_ocn=alvdr_ocn(i),  swvdr=swvdr(i),         & 
            alidr_ocn=alidr_ocn(i),  swidr=swidr(i),         &
            alvdf_ocn=alvdf_ocn(i),  swvdf=swvdf(i),         &
            alidf_ocn=alidf_ocn(i),  swidf=swidf(i),         &
            flwout_ocn=flwout_ocn(i),sst=sst(i),             &
            fsens_ocn=fsens_ocn(i),  shcoef=shcoef(i),       &
            flat_ocn=flat_ocn(i),    lhcoef=lhcoef(i),       &
            evap_ocn=evap_ocn(i),    flw=flw(i),             &
            delt=delt(i),            delq=delq(i),           & 
            aice=aice(i),            fhocn=fhocn(i),         &
            fswthru=fswthru(i),      hmix=hmix(i),           &
            Tf=Tf(i),                fresh=fresh(i),         &
            frain=frain(i),          fsnow=fsnow(i),         &
            fhocn_tot=fhocn_tot(i),  fresh_tot=fresh_tot(i), &
            frzmlt=frzmlt(i),        fsalt=fsalt(i),         &
            sss=sss(i)          ) 
    enddo                    ! i
    
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)

  end subroutine ocean_mixed_layer

  !=======================================================================

  subroutine ocn_mixed_layer_icepack( &
       ice,                   &
       alvdr_ocn, swvdr,      &
       alidr_ocn, swidr,      &
       alvdf_ocn, swvdf,      &
       alidf_ocn, swidf,      &
       sst,       flwout_ocn, &
       fsens_ocn, shcoef,     &
       flat_ocn,  lhcoef,     &
       evap_ocn,  flw,        &
       delt,      delq,       &
       aice,      fhocn,      &
       fswthru,   hmix,       &
       Tf,        fresh,      &
       frain,     fsnow,      &
       fhocn_tot, fresh_tot,  &
       frzmlt,    fsalt,      &
       sss)

    use g_forcing_param,  only: use_virt_salt
    use MOD_ICE
          
    implicit none
          
    real (kind=dbl_kind), intent(in) :: &
         alvdr_ocn, alidr_ocn, alvdf_ocn, alidf_ocn, & ! visible/near-IR direct/diffusive (fraction)
         swvdr, swvdf, swidr, swidf,                 & ! visible/near-IR direct/diffusive sw down (W/m^2)
         flw,                                        & ! incoming longwave radiation (W/m^2)
         Tf,                                         & ! freezing temperature (C)
         hmix,                                       & ! mixed layer depth (m)
         delt,                                       & ! potential temperature difference   (K)
         delq,                                       & ! specific humidity difference   (kg/kg)
         shcoef,                                     & ! transfer coefficient for sensible heat
         lhcoef,                                     & ! transfer coefficient for latent heat
         fswthru,                                    & ! shortwave penetrating to ocean (W/m^2)
         aice,                                       & ! ice area fraction
         sst,                                        & ! sea surface temperature (C)
         sss,                                        & ! sea surface salinity
         frain,                                      & ! rainfall rate (kg/m^2/s)
         fsnow,                                      & ! snowfall rate (kg/m^2/s)
         fsalt                                         ! salt flux from ice to the ocean (kg/m^2/s) 

    real (kind=dbl_kind), intent(inout) :: &
         flwout_ocn,                                 & ! outgoing longwave radiation (W/m^2)
         fsens_ocn,                                  & ! sensible heat flux (W/m^2)
         flat_ocn,                                   & ! latent heat flux   (W/m^2)
         evap_ocn,                                   & ! evaporative water flux (kg/m^2/s)
         fhocn,                                      & ! net heat flux to ocean (W/m^2)
         fresh,                                      & ! fresh water flux to ocean (kg/m^2/s)
         frzmlt,                                     & ! freezing/melting potential (W/m^2)
         fhocn_tot,                                  & ! net total heat flux to ocean (W/m^2)
         fresh_tot                                     ! fresh total water flux to ocean (kg/m^2/s)

    real (kind=dbl_kind), parameter :: frzmlt_max = c1000 ! max magnitude of frzmlt (W/m^2)

    real (kind=dbl_kind) :: &
         TsfK,                                       & ! surface temperature (K)
         swabs,                                      & ! surface absorbed shortwave heat flux (W/m^2)
         Tffresh,                                    & ! 0 C in K
         Lfresh,                                     &   
         Lvap,                                       & 
         lfs_corr,                                   & ! fresh water correction for linear free surface      
         stefan_boltzmann,                           &
         ice_ref_salinity

    character(len=*),parameter :: subname='(icepack_ocn_mixed_layer)'
    
    type(t_ice), target, intent(inout) :: ice  
    real(kind=WP), pointer  :: emiss_wat
    emiss_wat => ice%thermo%emiss_wat

    call icepack_query_parameters( &
         Tffresh_out=Tffresh, Lfresh_out=Lfresh, &
         stefan_boltzmann_out=stefan_boltzmann,  &
         ice_ref_salinity_out=ice_ref_salinity,  &
         Lvap_out=Lvap                           )

    ! shortwave radiative flux ! Visible is absorbed by clorophil
    ! afterwards
    swabs = (c1-alidr_ocn) * swidr + (c1-alidf_ocn) * swidf + (c1-alvdr_ocn) * swvdr + (c1-alvdf_ocn) * swvdf

    ! ocean surface temperature in Kelvin
    TsfK = sst + Tffresh

    ! longwave radiative flux
    ! Water emissivity added to be consistent
    ! with the standard FESOM2 version
    flwout_ocn = - emiss_wat * stefan_boltzmann * TsfK**4
    
    ! downward latent and sensible heat fluxes
    fsens_ocn =  shcoef * delt
    flat_ocn  =  lhcoef * delq
    evap_ocn  = -flat_ocn / Lvap

    ! Compute heat change due to exchange between ocean and atmosphere

    fhocn_tot = fhocn + fswthru                                  &  ! these are *aice already
         + (fsens_ocn + flat_ocn + flwout_ocn + flw + swabs &
         + Lfresh*fsnow) * (c1-aice) + max(c0,frzmlt)*aice
         
    if (use_virt_salt) then
       lfs_corr = fsalt/ice_ref_salinity/p001
       fresh = fresh - lfs_corr * ice_ref_salinity / sss
    endif

    !PS fresh_tot = fresh + (-evap_ocn + frain + fsnow)*(c1-aice)
    fresh_tot = fresh + frain + (-evap_ocn + fsnow)*(c1-aice)
    !            |         |         |        |--> depending on ice concentration eiter snow adds to the                  
    !            |         |         |             freshwater in the ocean or accumulates on the ice as snoe layer
    !            |         |         |
    !            |         |         |--> evaporation ocean-->atmosphere
    !            |         |
    !            |         |--> add here the total rain, at the end all the rain
    !            |              drains through the ice, therefor comment the line
    !            |              in ice_pack_therm_itd.F90, subroutine icepack_step_therm2()
    !            |              Line: 1999
    !            |              !!! If i dont do this here im not able to balance
    !            |              the ocean volume under zstar to maschine precision !!!
    !            | 
    !            |--> at that point fresh contains the freshwater flux contributions
    !                 from the thermodynamic growth rates of ice and snow but also 
    !                 the contributions from the sublimation of ice-->atmos (see.
    !                 icepack_therm_vertical.F90-->subroutine thermo_vertical(...)
    !                 Line: 453 --> freshn = freshn + evapn - (rhoi*dhi + rhos*dhs) / dt
    !                 evapn == evaporative water flux (kg/m^2/s)  from sublimation

  end subroutine ocn_mixed_layer_icepack

  !=======================================================================

  subroutine coupling_prep(ice, dt)
            
    use MOD_ICE  

    implicit none
          
    type(t_ice), target, intent(inout) :: ice
    
    real (kind=dbl_kind), intent(in) :: dt         ! time step
    
    ! local variables
    integer (kind=int_kind) :: &
         n,                                      & ! thickness category index
         i,                                      & ! horizontal index
         k,                                      & ! tracer index
         nbtrcr
    
    real (kind=dbl_kind) :: &
         netsw,                                  & ! flag for shortwave radiation presence
         rhofresh,                               &
         puny
    
    character(len=*), parameter :: subname='(coupling_prep)'
    
    !-----------------------------------------------------------------
    ! Save current value of frzmlt for diagnostics.
    ! Update mixed layer with heat and radiation from ice.
    !-----------------------------------------------------------------
    
    call icepack_query_parameters(puny_out=puny, rhofresh_out=rhofresh)
    call icepack_query_tracer_sizes(nbtrcr_out=nbtrcr)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    do i = 1, nx
       frzmlt_init(i) = frzmlt(i)
    enddo
    
    call ocean_mixed_layer (ice, dt) ! ocean surface fluxes and sst
    
    !-----------------------------------------------------------------
    ! Aggregate albedos
    !-----------------------------------------------------------------
    
    do i = 1, nx
       alvdf(i) = c0
       alidf(i) = c0
       alvdr(i) = c0
       alidr(i) = c0
       
       albice(i) = c0
       albsno(i) = c0
       albpnd(i) = c0
       apeff_ai(i) = c0
       snowfrac(i) = c0
    enddo
    do n = 1, ncat
       do i = 1, nx
          if (aicen(i,n) > puny) then
             
             alvdf(i) = alvdf(i) + alvdfn(i,n)*aicen(i,n)
             alidf(i) = alidf(i) + alidfn(i,n)*aicen(i,n)
             alvdr(i) = alvdr(i) + alvdrn(i,n)*aicen(i,n)
             alidr(i) = alidr(i) + alidrn(i,n)*aicen(i,n)
             
             netsw = swvdr(i) + swidr(i) + swvdf(i) + swidf(i)
             if (netsw > puny) then ! sun above horizon
                albice(i) = albice(i) + albicen(i,n)*aicen(i,n)
                albsno(i) = albsno(i) + albsnon(i,n)*aicen(i,n)
                albpnd(i) = albpnd(i) + albpndn(i,n)*aicen(i,n)
             endif
             
             apeff_ai(i) = apeff_ai(i) + apeffn(i,n)*aicen(i,n) ! for history
             snowfrac(i) = snowfrac(i) + snowfracn(i,n)*aicen(i,n) ! for history
             
          endif ! aicen > puny
       enddo
    enddo
    
    do i = 1, nx
       
       !-----------------------------------------------------------------
       ! reduce fresh by fpond for coupling
       !-----------------------------------------------------------------
       
       if (l_mpond_fresh) then
          fpond(i) = fpond(i) * rhofresh/dt
          fresh(i) = fresh(i) - fpond(i)
       endif
       
       !----------------------------------------------------------------
       ! Store grid box mean albedos and fluxes before scaling by aice
       !----------------------------------------------------------------
       
       alvdf_ai  (i) = alvdf  (i)
       alidf_ai  (i) = alidf  (i)
       alvdr_ai  (i) = alvdr  (i)
       alidr_ai  (i) = alidr  (i)
       fresh_ai  (i) = fresh  (i)
       fsalt_ai  (i) = fsalt  (i)
       fhocn_ai  (i) = fhocn  (i)
       fswthru_ai(i) = fswthru(i)
       fzsal_ai  (i) = fzsal  (i)
       fzsal_g_ai(i) = fzsal_g(i)
       
       if (nbtrcr > 0) then
          do k = 1, nbtrcr
             flux_bio_ai  (i,k) = flux_bio  (i,k)
          enddo
       endif
       
       !-----------------------------------------------------------------
       ! Save net shortwave for scaling factor in scale_factor
       !-----------------------------------------------------------------
       scale_factor(i) = &
            swvdr(i)*(c1 - alvdr_ai(i)) &
            + swvdf(i)*(c1 - alvdf_ai(i)) &
            + swidr(i)*(c1 - alidr_ai(i)) &
            + swidf(i)*(c1 - alidf_ai(i))
    enddo
    
  end subroutine coupling_prep

  !=======================================================================

  module subroutine step_icepack(flag_debug, ice, mesh, time_evp, time_advec, time_therm)

    use icepack_intfc, only: icepack_ice_strength
    use g_config, only: dt
    use MOD_MESH    
    use MOD_ICE    
    use ice_EVPdynamics_interface
    use ice_maEVPdynamics_interface

    implicit none
    
    logical (kind=log_kind), intent(in) :: flag_debug
    logical (kind=log_kind) :: &
         calc_Tsfc, skl_bgc, solve_zsal, z_tracers, tr_brine, &  ! from icepack
         tr_fsd, wave_spec
    integer (kind=int_kind) :: k, i                              ! for loop indexes   
    real (kind=dbl_kind) :: offset, t1, t2, t3, t4               ! d(age)/dt time offset
    real (kind=dbl_kind), intent(out) :: time_therm, time_advec, time_evp

    type(t_ice), target, intent(inout) :: ice  
    type(t_mesh), target, intent(in) :: mesh    

    character(len=*), parameter :: subname='(ice_step)'
    
    t1 = c0
    t2 = c0
    t3 = c0
    t4 = c0
    
    t1 = MPI_Wtime()    
    
    !-----------------------------------------------------------------
    ! query Icepack values
    !-----------------------------------------------------------------

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call icepack_query_...'//achar(27)//'[0m'
    call icepack_query_parameters(skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)
    call icepack_query_parameters(solve_zsal_out=solve_zsal, calc_Tsfc_out=calc_Tsfc, &
         wave_spec_out=wave_spec)
    call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_fsd_out=tr_fsd)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)
    
    ! TODO: Add appropriate timing
    
    !-----------------------------------------------------------------
    ! copy variables from fesom2 (also ice velocities)
    !-----------------------------------------------------------------

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call fesom_to_icepack'//achar(27)//'[0m'
    call fesom_to_icepack(ice, mesh)
    
    !-----------------------------------------------------------------
    ! tendencies needed by fesom
    !-----------------------------------------------------------------
    
    dhi_t_dt(:) = vice(:)
    dhs_t_dt(:) = vsno(:)
    
    !-----------------------------------------------------------------
    ! initialize diagnostics
    !-----------------------------------------------------------------
    
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call init_history'//achar(27)//'[0m'
    call init_history_therm
    call init_history_bgc
    
    !-----------------------------------------------------------------
    ! Scale radiation fields
    !-----------------------------------------------------------------

    if (calc_Tsfc) then
       if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call prep_radiation'//achar(27)//'[0m'
       call prep_radiation ()
    endif
    
    !-----------------------------------------------------------------
    ! thermodynamics and biogeochemistry
    !-----------------------------------------------------------------

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call step_therm1'//achar(27)//'[0m'
    call step_therm1     (dt) ! vertical thermodynamics
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call step_therm2'//achar(27)//'[0m'
    call step_therm2     (dt) ! ice thickness distribution thermo
    
    !-----------------------------------------------------------------
    ! clean up, update tendency diagnostics
    !-----------------------------------------------------------------
    offset = dt
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call update_state'//achar(27)//'[0m'
    call update_state (dt, daidtt, dvidtt, dagedtt, offset)
    
    !-----------------------------------------------------------------
    ! tendencies needed by fesom
    !-----------------------------------------------------------------
    
    dhi_t_dt(:) = ( vice(:) - dhi_t_dt(:) ) / dt
    dhs_t_dt(:) = ( vsno(:) - dhs_t_dt(:) ) / dt
    
    !-----------------------------------------------------------------
    ! dynamics, transport, ridging
    !-----------------------------------------------------------------

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call init_history_dyn'//achar(27)//'[0m'
    call init_history_dyn
    
    ! Compute sea-ice internal stress (immediately before EVP)
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call icepack_ice_strength'//achar(27)//'[0m'
    do i = 1, nx
       call icepack_ice_strength(ncat,                     &
            aice(i),     vice(i),     &
            aice0(i),    aicen(i,:),  &
            vicen(i,:),  strength(i))
    end do
    ! wave fracture of the floe size distribution
    ! note this is called outside of the dynamics subcycling loop
    if (tr_fsd .and. wave_spec) then
       if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call step_dyn_wave'//achar(27)//'[0m'
       call step_dyn_wave(dt)
    endif
    
    do k = 1, ndtd
       
       !-----------------------------------------------------------------
       ! EVP 
       !-----------------------------------------------------------------
       
       t2 = MPI_Wtime()
       
       select case (ice%whichEVP)
       case (0)
          if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics'//achar(27)//'[0m'
          call EVPdynamics  (ice, p_partit, mesh)
       case (1)
          if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics_m'//achar(27)//'[0m'
          call EVPdynamics_m(ice, p_partit, mesh)
       case (2)
          if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics_a'//achar(27)//'[0m'
          call EVPdynamics_a(ice, p_partit, mesh)
       case default
          if (mype==0) write(*,*) 'A non existing EVP scheme specified!'
          call par_ex(p_partit%MPI_COMM_FESOM, p_partit%mype)
          stop
       end select
       
       t3 = MPI_Wtime()
       time_evp = t3 - t2
       
       !-----------------------------------------------------------------
       ! update ice velocities
       !-----------------------------------------------------------------
       
       if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call fesom_to_icepack'//achar(27)//'[0m'
       call fesom_to_icepack(ice, mesh)
       
       !-----------------------------------------------------------------
       ! advect tracers
       !-----------------------------------------------------------------
       
       t2 = MPI_Wtime()
       if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call tracer_advection_icepack'//achar(27)//'[0m'
       call tracer_advection_icepack(ice, mesh)

       t3 = MPI_Wtime()
       time_advec = t3 - t2
       
       !-----------------------------------------------------------------
       ! initialize tendencies needed by fesom
       !-----------------------------------------------------------------
       dhi_r_dt(:) = vice(:)
       dhs_r_dt(:) = vsno(:)
       
       !-----------------------------------------------------------------
       ! ridging
       !-----------------------------------------------------------------

       if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call step_dyn_ridge'//achar(27)//'[0m'
       call step_dyn_ridge (dt_dyn, ndtd)
       
       ! clean up, update tendency diagnostics
       offset = c0
       if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call update_stat'//achar(27)//'[0m'
       call update_state (dt_dyn, daidtd, dvidtd, dagedtd, offset)
       
       !-----------------------------------------------------------------
       ! tendencies needed by fesom
       !-----------------------------------------------------------------
       ! --> ridging adds fresh water need to know for compensation of thermodynamic
       !     growth rate of ice and snow in fesom
       dhi_r_dt(:) = ( vice(:) - dhi_r_dt(:) ) / dt
       dhs_r_dt(:) = ( vsno(:) - dhs_r_dt(:) ) / dt
       
    enddo
    
    !-----------------------------------------------------------------
    ! total tendencies of thermodynamic and ridging needed by fesom
    !-----------------------------------------------------------------
    ! --> needed for total compenstion of fresh of thgr and thgrsn
    dhi_dt(:) = dhi_r_dt(:) + dhi_t_dt(:) 
    dhs_dt(:) = dhs_r_dt(:) + dhs_t_dt(:)
    
    !-----------------------------------------------------------------
    ! albedo, shortwave radiation
    !-----------------------------------------------------------------
    
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call step_radiation'//achar(27)//'[0m'
    call step_radiation (dt)
    
    !-----------------------------------------------------------------
    ! get ready for coupling and the next time step
    !-----------------------------------------------------------------
    
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call coupling_prep'//achar(27)//'[0m'
    call coupling_prep (ice, dt)
    
    !-----------------------------------------------------------------
    ! icepack timing
    !-----------------------------------------------------------------  
    
    t4 = MPI_Wtime()
    time_therm = t4 - t1 - time_advec - time_evp
    
    !time_advec = c0
    !time_therm = c0
    !time_evp   = c0
    
  end subroutine step_icepack
  
  !=======================================================================
  
end submodule icedrv_step

!=======================================================================
