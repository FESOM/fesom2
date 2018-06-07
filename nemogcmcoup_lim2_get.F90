SUBROUTINE nemogcmcoup_lim2_get( mype, npes, icomm, &
   &                             nopoints, pgsst, pgist, pgalb, &
   &                             pgifr, pghic, pghsn, pgucur, pgvcur, &
   &                             pgistl, licelvls )

   ! Interpolate sst, ice: surf T; albedo; concentration; thickness,
   ! snow thickness and currents from the ORCA grid to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind

   IMPLICIT NONE
   
   ! Arguments
   REAL(wp), DIMENSION(nopoints) :: pgsst, pgist, pgalb, pgifr, pghic, pghsn, pgucur, pgvcur
   REAL(wp), DIMENSION(nopoints,3) :: pgistl
   LOGICAL :: licelvls

   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number Gaussian grid points
   INTEGER, INTENT(IN) :: nopoints

   ! Local variables

#ifdef FESOM_TODO

   ! Temporary array for packing of input data without halos.
   REAL(wp), DIMENSION((nlei-nldi+1)*(nlej-nldj+1)) :: zsend
   ! Arrays for rotation of current vectors from ij to ne.
   REAL(wp), DIMENSION(jpi,jpj) :: zotx1, zoty1, ztmpx, ztmpy
   ! Array for fraction of leads (i.e. ocean)
   REAL(wp), DIMENSION(jpi,jpj) :: zfr_l 
   REAL(wp) :: zt
   ! Loop variables
   INTEGER :: ji, jj, jk, jl
   REAL(wp) :: zhook_handle ! Dr Hook handle
   
   IF(lhook) CALL dr_hook('nemogcmcoup_lim2_get',0,zhook_handle)
   IF(nn_timing == 1) CALL timing_start('nemogcmcoup_lim2_get')

   zfr_l(:,:) = 1.- fr_i(:,:)
   
   IF (.NOT.ALLOCATED(zscplsst)) THEN
      ALLOCATE(zscplsst(jpi,jpj))
   ENDIF

   ! Pack SST data and convert to K.

   IF ( nsstlvl(1) == nsstlvl(2) ) THEN
      jk = 0 
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            jk = jk + 1
            zsend(jk) = tsn(ji,jj,nsstlvl(1),jp_tem) + rt0
            zscplsst(ji,jj) = zsend(jk) - rt0
         ENDDO
      ENDDO
   ELSE
      jk = 0 
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            jk = jk + 1
            zsend(jk) = SUM(&
               & tsn(ji,jj,nsstlvl(1):nsstlvl(2),jp_tem) * &
               & tmask(ji,jj,nsstlvl(1):nsstlvl(2)) * &
               & fse3t(ji,jj,nsstlvl(1):nsstlvl(2)) ) / &
               & MAX( SUM( &
               & tmask(ji,jj,nsstlvl(1):nsstlvl(2)) * &
               & fse3t(ji,jj,nsstlvl(1):nsstlvl(2))) , 1.0 ) + rt0
            zscplsst(ji,jj) = zsend(jk) - rt0
         ENDDO
      ENDDO
   ENDIF
   CALL lbc_lnk( zscplsst, 'T', 1. )

   ! Interpolate SST

   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
      &               nopoints, pgsst )

   ! Pack ice temperature data (already in K)

#if defined key_lim2
   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = tn_ice(ji,jj,1)
      ENDDO
   ENDDO
#else
   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = 0
         zt=0.0
         DO jl = 1, jpl
            zsend(jk) = zsend(jk) + tn_ice(ji,jj,jl) * a_i(ji,jj,jl)
            zt = zt + a_i(ji,jj,jl)
         ENDDO
         IF ( zt > 0.0 ) THEN
            zsend(jk) = zsend(jk) / zt
         ELSE
            zsend(jk) = rt0
         ENDIF
      ENDDO
   ENDDO
#endif
   
   ! Interpolate ice temperature 

   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
      &               nopoints, pgist )

   ! Ice level temperatures

   IF (licelvls) THEN

#if defined key_lim2

      DO jl = 1, 3
         
         ! Pack ice temperatures data at level jl(already in K)
         
         jk = 0 
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               jk = jk + 1
               zsend(jk) = tbif (ji,jj,jl)
            ENDDO
         ENDDO
         
         ! Interpolate ice temperature  at level jl
         
         CALL parinter_fld( mype, npes, icomm, Ttogauss, &
            &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
            &               nopoints, pgistl(:,jl) )
         
      ENDDO

#else
      WRITE(0,*)'licelvls needs to be sorted for LIM3'
      CALL abort
#endif     

   ENDIF

   ! Pack ice albedo data 

#if defined key_lim2
   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = alb_ice(ji,jj,1)
      ENDDO
   ENDDO
#else
   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = 0
         zt=0.0
         DO jl = 1, jpl
            zsend(jk) = zsend(jk) + alb_ice(ji,jj,jl) * a_i(ji,jj,jl)
            zt = zt + a_i(ji,jj,jl)
         ENDDO
         IF ( zt > 0.0_wp ) THEN
            zsend(jk) = zsend(jk) / zt
         ELSE
            zsend(jk) = albedo_oce_mix(ji,jj)
         ENDIF
      ENDDO
   ENDDO
#endif
   
   ! Interpolate ice albedo

   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
      &               nopoints, pgalb )

   ! Pack ice fraction data

   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = fr_i(ji,jj)
      ENDDO
   ENDDO

   ! Interpolation of ice fraction.

   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
      &               nopoints, pgifr )

   ! Pack ice thickness data

#if defined key_lim2
   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = hicif(ji,jj)
      ENDDO
   ENDDO
#else
   ! LIM3
   ! Average over categories (to be revised).
   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = 0
         DO jl = 1, jpl
            zsend(jk) = zsend(jk) + ht_i(ji,jj,jl) * a_i(ji,jj,jl)
         ENDDO
      ENDDO
   ENDDO
#endif

   ! Interpolation of ice thickness

   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
      &               nopoints, pghic )  

   ! Pack snow thickness data

#if defined key_lim2
   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = hsnif(ji,jj)
      ENDDO
   ENDDO
#else
   ! LIM3
   ! Average over categories (to be revised).
   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = 0
         DO jl = 1, jpl
            zsend(jk) = zsend(jk) + ht_s(ji,jj,jl) * a_i(ji,jj,jl)
         ENDDO
      ENDDO
   ENDDO
#endif

   ! Interpolation of snow thickness

   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
      &               nopoints, pghsn )

   ! Currents needs to be rotated from ij to ne first

   DO jj = 2, jpjm1
      DO ji = 2, jpim1
         zotx1(ji,jj) = 0.5 * ( un(ji,jj,1) + un(ji-1,jj  ,1) )
         zoty1(ji,jj) = 0.5 * ( vn(ji,jj,1) + vn(ji  ,jj-1,1) ) 
      END DO
   END DO
   CALL lbc_lnk( zotx1, 'T', -1. )
   CALL lbc_lnk( zoty1, 'T', -1. )
   CALL rot_rep( zotx1, zoty1, 'T', 'ij->e', ztmpx )
   CALL rot_rep( zotx1, zoty1, 'T', 'ij->n', ztmpy )

   ! Pack U current

   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = ztmpx(ji,jj)
      ENDDO
   ENDDO

   ! Interpolate U current

   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
      &               nopoints, pgucur )

   ! Pack V current

   jk = 0 
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = jk + 1
         zsend(jk) = ztmpy(ji,jj)
      ENDDO
   ENDDO

   ! Interpolate V current

   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
      &               nopoints, pgvcur )

   IF(nn_timing == 1) CALL timing_stop('nemogcmcoup_lim2_get')
   IF(lhook) CALL dr_hook('nemogcmcoup_lim2_get',1,zhook_handle)

#else

   WRITE(0,*)'nemogcmcoup_lim2_get not done for FESOM yet'
   CALL abort

#endif

END SUBROUTINE nemogcmcoup_lim2_get
   
