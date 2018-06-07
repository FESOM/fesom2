SUBROUTINE nemogcmcoup_lim2_update( mype, npes, icomm, &
   &                                npoints,  &
   &                                taux_oce, tauy_oce, taux_ice, tauy_ice, &
   &                                qs___oce, qs___ice, qns__oce, qns__ice, dqdt_ice, &
   &                                evap_tot, evap_ice, prcp_liq, prcp_sol, &
   &                                runoff, ocerunoff, tcc, lcc, tice_atm, &
   &                                kt, ldebug, loceicemix, lqnsicefilt )

   ! Update fluxes in nemogcmcoup_data by parallel
   ! interpolation of the input gaussian grid data
   
   USE par_kind

   IMPLICIT NONE

   ! Arguments

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Fluxes on the Gaussian grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wp), DIMENSION(npoints), INTENT(IN) :: &
      & taux_oce, tauy_oce, taux_ice, tauy_ice, &
      & qs___oce, qs___ice, qns__oce, qns__ice, &
      & dqdt_ice, evap_tot, evap_ice, prcp_liq, prcp_sol, &
      & runoff, ocerunoff, tcc, lcc, tice_atm

   ! Current time step
   INTEGER, INTENT(in) :: kt
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug
   ! QS/QNS mixed switch
   LOGICAL, INTENT(IN) :: loceicemix
   ! QNS ice filter switch (requires tice_atm to be sent)
   LOGICAL, INTENT(IN) :: lqnsicefilt

   ! Local variables

#ifdef FESOM_TODO

   ! Packed receive buffer
   REAL(wp), DIMENSION((nlei-nldi+1)*(nlej-nldj+1)) :: zrecv
   ! Unpacked fields on ORCA grids
   REAL(wp), DIMENSION(jpi,jpj) :: zqs___oce, zqs___ice, zqns__oce, zqns__ice
   REAL(wp), DIMENSION(jpi,jpj) :: zdqdt_ice, zevap_tot, zevap_ice, zprcp_liq, zprcp_sol
   REAL(wp), DIMENSION(jpi,jpj) :: zrunoff, zocerunoff
   REAL(wp), DIMENSION(jpi,jpj) :: ztmp, zicefr
   ! Arrays for rotation
   REAL(wp), DIMENSION(jpi,jpj) :: zuu,zvu,zuv,zvv,zutau,zvtau 
   ! Lead fraction for both LIM2/LIM3
   REAL(wp), DIMENSION(jpi,jpj) :: zfrld
   ! Mask for masking for I grid
   REAL(wp) :: zmsksum
   ! For summing up LIM3 contributions to ice temperature
   REAL(wp) :: zval,zweig

   ! Loop variables
   INTEGER :: ji,jj,jk,jl
   ! netCDF debugging output variables
   CHARACTER(len=128) :: cdoutfile
   INTEGER :: inum
   REAL(wp) :: zhook_handle ! Dr Hook handle

   IF(lhook) CALL dr_hook('nemogcmcoup_lim2_update',0,zhook_handle)
   IF(nn_timing == 1) CALL timing_start('nemogcmcoup_lim2_update')
   
   ! Allocate the storage data

   IF (.NOT.lallociceflx) THEN
      ALLOCATE( &
         & zsqns_tot(jpi,jpj),   &
         & zsqns_ice(jpi,jpj),   &
         & zsqsr_tot(jpi,jpj),   &
         & zsqsr_ice(jpi,jpj),   &
         & zsemp_tot(jpi,jpj),   &
         & zsemp_ice(jpi,jpj),   &
	 & zsevap_ice(jpi,jpj),  &
         & zsdqdns_ice(jpi,jpj), &
         & zssprecip(jpi,jpj),   &
	 & zstprecip(jpi,jpj),   &
         & zstcc(jpi,jpj),       &
         & zslcc(jpi,jpj),       &
         & zsatmist(jpi,jpj),    &
         & zsqns_ice_add(jpi,jpj)&
         & )
      lallociceflx = .TRUE.
   ENDIF
   IF (.NOT.lallocstress) THEN
      ALLOCATE( &
         & zsutau(jpi,jpj),     &
         & zsvtau(jpi,jpj),     &
         & zsutau_ice(jpi,jpj), &
         & zsvtau_ice(jpi,jpj)  &
         & )
      lallocstress = .TRUE.
   ENDIF

   ! Sort out incoming arrays from the IFS and put them on the ocean grid
   
   !1. Interpolate ocean solar radiation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, qs___oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean solar radiation

   zqs___oce(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zqs___oce(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !2. Interpolate ice solar radiation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, qs___ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ice solar radiation

   zqs___ice(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zqs___ice(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !3. Interpolate ocean non-solar radiation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, qns__oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean non-solar radiation

   zqns__oce(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zqns__oce(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !4. Interpolate ice non-solar radiation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, qns__ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ice non-solar radiation

   zqns__ice(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zqns__ice(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !5. Interpolate  D(q)/dT to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, dqdt_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack D(q)/D(T) 

   zdqdt_ice(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zdqdt_ice(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !6. Interpolate total evaporation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, evap_tot,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack total evaporation

   zevap_tot(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zevap_tot(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !7. Interpolate evaporation over ice to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, evap_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack evaporation over ice

   zevap_ice(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zevap_ice(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
 
   !8. Interpolate liquid precipitation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, prcp_liq,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack liquid precipitation

   zprcp_liq(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zprcp_liq(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !9. Interpolate solid precipitation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, prcp_sol,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack precipitation over ice

   zprcp_sol(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zprcp_sol(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !10. Interpolate runoff to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, runoff,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack runoff

   zrunoff(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zrunoff(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !11. Interpolate ocean runoff to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, ocerunoff,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean runoff

   zocerunoff(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zocerunoff(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !12. Interpolate total cloud fractions to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, tcc,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean runoff

   zstcc(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zstcc(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !13. Interpolate low cloud fractions to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, lcc,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean runoff

   zslcc(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zslcc(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! get sea ice fraction and lead fraction

#if defined key_lim2
   zfrld(:,:) = frld(:,:)
   zicefr(:,:) = 1 - zfrld(:,:)
#else
   zicefr(:,:) = 0.0_wp
   DO jl = 1, jpl
      zicefr(:,:) = zicefr(:,:) + a_i(:,:,jl)
   ENDDO
   zfrld(:,:) = 1 - zicefr(:,:)
#endif
    
   zsemp_tot(:,:) = zevap_tot(:,:) - zprcp_liq(:,:) - zprcp_sol(:,:)
   zstprecip(:,:) = zprcp_liq(:,:) + zprcp_sol(:,:)
   ! More consistent with NEMO, but does changes the results, so
   ! we don't do it for now.
   ! zsemp_tot(:,:) = zevap_tot(:,:) - zstprecip(:,:)
   zsemp_ice(:,:) = zevap_ice(:,:) - zprcp_sol(:,:)
   zssprecip(:,:) = - zsemp_ice(:,:)
   zsemp_tot(:,:) = zsemp_tot(:,:) - zrunoff(:,:)
   zsemp_tot(:,:) = zsemp_tot(:,:) - zocerunoff(:,:)
   zsevap_ice(:,:) = zevap_ice(:,:)

   !   non solar heat fluxes   !  (qns)
   IF (loceicemix) THEN
      zsqns_tot(:,:) =  zqns__oce(:,:)
   ELSE
      zsqns_tot(:,:) =  zfrld(:,:) * zqns__oce(:,:) + zicefr(:,:) * zqns__ice(:,:)
   ENDIF
   zsqns_ice(:,:) =  zqns__ice(:,:)
   ztmp(:,:) = zfrld(:,:) * zprcp_sol(:,:) * lfus  ! add the latent heat of solid precip. melting

   zsqns_tot(:,:) = zsqns_tot(:,:) - ztmp(:,:)    ! over free ocean 
   !      solar heat fluxes    !   (qsr)
  
   IF (loceicemix) THEN
      zsqsr_tot(:,:) =  zqs___oce(:,:)
   ELSE
      zsqsr_tot(:,:) =  zfrld(:,:) * zqs___oce(:,:) + zicefr(:,:) * zqs___ice(:,:)
   ENDIF
   zsqsr_ice(:,:) =  zqs___ice(:,:)
   
   IF( ln_dm2dc ) THEN   ! modify qsr to include the diurnal cycle
      zsqsr_tot(:,:) = sbc_dcy( zsqsr_tot(:,:) )
      zsqsr_ice(:,:) = sbc_dcy( zsqsr_ice(:,:) )
   ENDIF
  
   zsdqdns_ice(:,:) = zdqdt_ice(:,:)

   ! Apply lateral boundary condition
   
   CALL lbc_lnk(zsqns_tot, 'T', 1.0)
   CALL lbc_lnk(zsqns_ice, 'T', 1.0)
   CALL lbc_lnk(zsqsr_tot, 'T', 1.0)
   CALL lbc_lnk(zsqsr_ice, 'T', 1.0)
   CALL lbc_lnk(zsemp_tot, 'T', 1.0)
   CALL lbc_lnk(zsemp_ice, 'T', 1.0)
   CALL lbc_lnk(zsdqdns_ice, 'T', 1.0)
   CALL lbc_lnk(zssprecip, 'T', 1.0)
   CALL lbc_lnk(zstprecip, 'T', 1.0)
   CALL lbc_lnk(zstcc, 'T', 1.0)
   CALL lbc_lnk(zslcc, 'T', 1.0)

   ! Interpolate  atmospheric ice temperature to T grid
      
   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, tice_atm,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack atmospheric ice temperature
      
   zsatmist(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zsatmist(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   CALL lbc_lnk(zsatmist, 'T', 1.0)
   
   zsqns_ice_add(:,:) = 0.0_wp

   ! Use the dqns_ice filter

   IF (lqnsicefilt) THEN

      ! Add filtr to qns_ice
      
#if defined key_lim2 
      ztmp(:,:) = tn_ice(:,:,1)
#else
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zval=0.0
            zweig=0.0
            DO jl = 1, jpl
               zval = zval + tn_ice(ji,jj,jl) * a_i(ji,jj,jl)
               zweig = zweig + a_i(ji,jj,jl)
            ENDDO
            IF ( zweig > 0.0 ) THEN
               ztmp(ji,jj) = zval /zweig
            ELSE
               ztmp(ji,jj) = rt0
            ENDIF
         ENDDO
      ENDDO
      CALL lbc_lnk(ztmp, 'T', 1.0)
#endif

      WHERE ( zicefr(:,:) > .001_wp )
         zsqns_ice_add(:,:) = zsdqdns_ice(:,:) * ( ztmp(:,:) - zsatmist(:,:) )
      END WHERE

      zsqns_ice(:,:) = zsqns_ice(:,:) + zsqns_ice_add(:,:)
      
   ENDIF

   ! Interpolate u-stress to U grid

   CALL parinter_fld( mype, npes, icomm, gausstoU, npoints,taux_oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack u stress on U grid

   zuu(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zuu(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate v-stress to U grid

   CALL parinter_fld( mype, npes, icomm, gausstoU, npoints, tauy_oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack v stress on U grid

   zvu(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zvu(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate u-stress to V grid

   CALL parinter_fld( mype, npes, icomm, gausstoV, npoints,taux_oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack u stress on V grid

   zuv(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zuv(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate v-stress to V grid

   CALL parinter_fld( mype, npes, icomm, gausstoV, npoints, tauy_oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack v stress on V grid

   zvv(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zvv(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   ! Rotate stresses from en to ij and put u,v stresses on U,V grids
   
   CALL repcmo( zuu, zvu, zuv, zvv, zsutau, zsvtau )

   ! Apply lateral boundary condition on u,v stresses on the U,V grids

   CALL lbc_lnk( zsutau, 'U', -1.0 )
   CALL lbc_lnk( zsvtau, 'V', -1.0 )

   ! Interpolate ice u-stress to U grid

   CALL parinter_fld( mype, npes, icomm, gausstoU, npoints,taux_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack ice u stress on U grid

   zuu(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zuu(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate ice v-stress to U grid

   CALL parinter_fld( mype, npes, icomm, gausstoU, npoints, tauy_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ice v stress on U grid

   zvu(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zvu(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate ice u-stress to V grid

   CALL parinter_fld( mype, npes, icomm, gausstoV, npoints,taux_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack ice u stress on V grid

   zuv(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zuv(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate ice v-stress to V grid

   CALL parinter_fld( mype, npes, icomm, gausstoV, npoints, tauy_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack ice v stress on V grid

   zvv(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zvv(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   ! Rotate stresses from en to ij and put u,v stresses on U,V grids

   CALL repcmo( zuu, zvu, zuv, zvv, zutau, zvtau )

   ! Apply lateral boundary condition on u,v stresses on the U,V grids

   CALL lbc_lnk( zutau, 'U', -1.0 )
   CALL lbc_lnk( zvtau, 'V', -1.0 )

#if defined key_lim2_vp

   ! Convert to I grid for LIM2 for key_lim_vp
   DO jj = 2, jpjm1                                   ! (U,V) ==> I
      DO ji = 2, jpim1   ! NO vector opt.
         zmsksum = umask(ji-1,jj,1) + umask(ji-1,jj-1,1) 
         zsutau_ice(ji,jj) = ( umask(ji-1,jj,1) * zutau(ji-1,jj) + &
            &                  umask(ji-1,jj-1,1) * zutau(ji-1,jj-1) )
         IF ( zmsksum > 0.0 ) THEN
            zsutau_ice(ji,jj) = zsutau_ice(ji,jj) / zmsksum
         ENDIF
         zmsksum = vmask(ji,jj-1,1) + vmask(ji-1,jj-1,1) 
         zsvtau_ice(ji,jj) = ( vmask(ji,jj-1,1) * zvtau(ji,jj-1) + &
            &                  vmask(ji-1,jj-1,1) * zvtau(ji-1,jj-1) )
         IF ( zmsksum > 0.0 ) THEN
            zsvtau_ice(ji,jj) = zsvtau_ice(ji,jj) / zmsksum
         ENDIF
      END DO
   END DO

#else
   
   zsutau_ice(:,:) = zutau(:,:)
   zsvtau_ice(:,:) = zvtau(:,:)

#endif

   CALL lbc_lnk( zsutau_ice, 'I', -1.0 )
   CALL lbc_lnk( zsvtau_ice, 'I', -1.0 )
   
   ! Optionally write files write the data on the ORCA grid via IOM.

   IF (ldebug) THEN      
      WRITE(cdoutfile,'(A,I8.8)') 'zsutau_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsutau' , zsutau )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsvtau_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsvtau' , zsvtau )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsutau_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsutau_ice' , zsutau_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsvtau_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsvtau_ice' , zsvtau_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqns_tot_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqns_tot' , zsqns_tot )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqns_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqns_ice' , zsqns_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqsr_tot_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqsr_tot' , zsqsr_tot )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqsr_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqsr_ice' , zsqsr_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsemp_tot_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsemp_tot' , zsemp_tot )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsemp_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsemp_ice' , zsemp_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsdqdns_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsdqdns_ice' , zsdqdns_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zssprecip_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zssprecip' , zssprecip )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zstprecip_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zstprecip' , zstprecip )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsevap_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsevap_ice' , zsevap_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zstcc_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zstcc' , zstcc )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zslcc_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zslcc' , zslcc )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsatmist_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsatmist' , zsatmist )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqns_ice_add_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqns_ice_add' , zsqns_ice_add )
      CALL iom_close( inum )
   ENDIF

   IF(nn_timing == 1) CALL timing_stop('nemogcmcoup_lim2_update')
   IF(lhook) CALL dr_hook('nemogcmcoup_lim2_update',1,zhook_handle)

#else

   WRITE(0,*)'nemogcmcoup_lim2_update not done for FESOM yet'
   CALL abort

#endif

END SUBROUTINE nemogcmcoup_lim2_update

   
