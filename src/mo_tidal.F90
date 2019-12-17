!> Tidal Module
!> Real time computation of luni-solar gravitational potential
!> @author Maik Thomas, Ann-Margrit Hellmich
!> @date Last modified 9.9.2009 Malte Mueller
!! Included subroutines are:
!! alloc_mem_tidal   -- Allocation of memory
!! foreph_ini        -- Initialization of tidal module
!!    eph            -- determination of the phase with respect to 1.1.2000
!! foreph            -- main routine for calculation of realtime lunisolar
!!                      gravitational potential
!!    ephvsop87      -- ephemeredes of sun and moon
!!        sidt2      -- computes sidereal time
!!        obliq      -- calculation of obliquity of ecliptic
!!        Sun_n      -- calculation of position of sun
!!          anomaly  -- calculation of true anomaly AT and eccentric anomaly AE
!!        moon      -- calculation of position of moon
!!        aufb2
!!    eqecl          -- conversion of ecliptic into equatorial coordinates
!! tipouv            -- horizontal gradients of gravitational potential


MODULE mo_tidal
      USE o_PARAM
      USE o_ARRAYS, only : ssh_gp
      USE o_MESH,   only : geo_coord_nod2D
      USE g_config, only : dt
      USE g_PARSUP
      USE g_clock
      IMPLICIT NONE
      !Earth Tides ( maik thomas, emr  pers. comm )
      REAL(wp), PARAMETER :: sal=0.69_wp

      INTEGER :: mmccdt      ! mmccdt is the number of timesteps since 01.01.2000 00:00:00

      CONTAINS


      SUBROUTINE foreph_ini(lyear,lmonth)
      ! Initialization of tidal module
      ! Determination of Julian Day of first time step
      ! Projection of mpiom grid on tidal module internal coordinates

      IMPLICIT NONE

      INTEGER,INTENT(IN)::lyear,lmonth
      INTEGER :: i, j, jcc, moph

      mmccdt = 0; jcc = 0; moph = 0
      CALL eph(lyear,lmonth,jcc,moph)
      ! FIXME: why not trot from mo_planetary_constants?
      mmccdt = (jcc + moph - 0.5_wp) * NINT(86400._wp / dt)
      ! mmccdt is the number of timesteps since 01.01.2000 00:00:00
      ! bugfix by M.Mueller : substracting 0.5_wp shifts the reference time by -12h to 01.01.2000 12:00:00
      ! This is needed to get the correct sideral time (sbr sidt2)
      ! and calculate the correct hour angle for the tidal potential (->sbr aufb2)
      ! FIXME : mmccdt = ldtmon-1 + int(jcc + moph - 0.5_wp) * NINT(86400._wp / dt)
      ! FIXME : replace eph by a some to code that directly calculates julian days and
      ! centuries as needed by siderial time and ephemerides

      WRITE(*,*)'tidal: phase relative to 2000 :'    &
      ,'year= ',lyear, 'month= ',lmonth, 'yearoff= ',jcc,' monoff= ',moph ,'mmccdt= ',mmccdt

      END SUBROUTINE foreph_ini

      SUBROUTINE eph(jul,mon,jahrph,moph)
      !    Berechnet Phase [in Tagen]  bezueglich des 1.1.2000
      !       jahrph=phase zu beginn des Jahres
      !       gesamtphase=jahrph+moph

      IMPLICIT NONE

      INTEGER,INTENT(IN)::jul,mon
      INTEGER,INTENT(OUT)::jahrph,moph
      INTEGER::jpl,jj,m

      ! Before year 2000
      IF(jul.LT.2000) THEN
        jahrph=0
        DO jj=jul,1999
          jpl=365

          IF ((MOD(jj,4).EQ.0  .AND. MOD(jj,100).NE.0)  .OR. MOD(jj,400).EQ.0  ) jpl = 366

          jahrph=jahrph-jpl
        ENDDO
      ENDIF
      ! year 2000
      IF(jul.EQ.2000) jahrph=0
      ! After year 2000
      IF(jul.GT.2000) THEN
        jahrph=0
        DO jj=2000,jul-1
          jpl=365

          IF ((MOD(jj,4).EQ.0  .AND. MOD(jj,100).NE.0)  .OR. MOD(jj,400).EQ.0  ) jpl = 366

          jahrph=jahrph+jpl
        ENDDO
      ENDIF

      moph=0
      IF (mon.GT.1)THEN
        DO m=1,mon-1
!         moph=moph+monlen(m,jul)
          moph=moph+num_day_in_month(fleapyear,m)
        ENDDO
      ENDIF

      END SUBROUTINE eph

      SUBROUTINE foreph
!     calculates the realtime gravitational potential of sun & moon
!     input : lon_rad,lat_rad
!     output: tipoto

      IMPLICIT NONE

      REAL(wp) :: dres(3,2),crim3,rkomp,erdrad,rekts,dekls
      REAL(wp) :: cris3,rektm,deklm,deklm2,dekls2,sidm,sidmq
      REAL(wp) :: rkosp,codm,codmq,sids,sidsq,cods,codsq,sidm2
      REAL(wp) :: sids2,hamp,hasp
      INTEGER :: i,j

      mmccdt = mmccdt + 1

      rkomp = -4.113e-07_wp ! factor of the tidal potential due to the moon
                            ! attention the factor is defined negative (contrary to the standard).
      rkosp = 0.46051_wp * rkomp
      erdrad = 6371000._wp

      CALL ephvsop87(dres)

      rekts=dres(1,1)
      dekls=dres(2,1)
      cris3=dres(3,1)

      rektm=dres(1,2)
      deklm=dres(2,2)
      crim3=dres(3,2)


      deklm2 = deklm * 2._wp
      dekls2 = dekls * 2._wp
      sidm   = Sin(deklm)
      sidmq  = sidm*sidm
      codm   = Cos(deklm)
      codmq  = codm*codm
      sids   = Sin(dekls)
      sidsq  = sids*sids
      cods   = Cos(dekls)
      codsq  = cods*cods
      sidm2=Sin(deklm2)
      sids2=Sin(dekls2)


      DO i=1, myDim_nod2D+eDim_nod2D

         hamp = rektm + geo_coord_nod2D(1,i)
         hasp = rekts + geo_coord_nod2D(1,i)

         ! attention : mpiom uses a negative tidal potential due to negative factor rkomp
         ssh_gp(i) = erdrad * rkomp * crim3 &
              * (3._wp * (SIN(geo_coord_nod2D(2,i))**2 - 1._wp/3._wp) * (sidmq - 1._wp/3._wp)&
              &  + SIN(2._wp * geo_coord_nod2D(2,i)) * sidm2 * COS(hamp) &
              &  + COS(geo_coord_nod2D(2,i))**2 * codmq * COS(2._wp * hamp))        &
              &  + erdrad * rkosp * cris3 &
              &    * (3._wp * (SIN(geo_coord_nod2D(2,i))**2 - 1._wp/3._wp) &
              &       * (sidsq - 1._wp/3._wp) &
              &       + SIN(2._wp * geo_coord_nod2D(2,i)) * sids2 * COS(hasp) &
              &       + COS(geo_coord_nod2D(2,i))**2 * codsq * COS(2._wp * hasp))

      END DO

!     no exchange needed
!     CALL bounds_exch(tipoto)

      END SUBROUTINE foreph

      Subroutine ephvsop87(res2)
!     calcuates the ephemeredes of sun/moon
!     input: mmccdt : time step information
!     output: res(3,2) right ascension, declination and geocentric distance
!             of Sun and Moon
!
      IMPLICIT NONE
!
!  optional preparation for calculation of potential
!  according to ephaufb.f by Maik Thomas
!
      INTEGER fnut
      REAL(wp) :: pi2,pic,T,sidt,ecl,nutob,nutl,res(3,2),res2(3,2)

      pi2 = dacos(-1.d0) * 2._wp
      pic = dacos(-1.d0)/DBLE(180.d0)

!     Transformation of Time into fractional julian centuries T
      t = (REAL(mmccdt-1, wp) * dt / 86400._wp)/36525._wp

!     corresponding sidereal time Greenwich
      Call sidt2(pic,pi2,T,sidt)

!  set fnut (perform nutation -> 1; don't -> 0)
      fnut=0

! C obliquity of the ecliptic
      CALL obliq(fnut,pic,T,ecl,nutob,nutl)

! C Sun
      CALL sun_n(fnut,pic,pi2,T,ecl,nutl,res)

! C Moon
      CALL moon(fnut,pic,pi2,T,ecl,nutl,res)

! C modifications in preparaion of calculation of potentials
      CALL aufb2(sidt,res,res2)
!
      END SUBROUTINE ephvsop87

      SUBROUTINE sidt2(pic,pi2,T,sidt)
!  convertion of Julian Centuries since J2000 T to siderical time sidt
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      REAL(WP) :: pic,pi2,T,sidt,JD,T2,T3
!
      T2=T*T
      T3=T2*T
!
!  Julian days
      jd = t * 36525.0_wp + 2451545.0_wp
!
!  mean siderial time of Greenwich in rad
      sidt = (280.46061837_WP + 360.98564736629_WP * (jd - 2451545.0_WP) &
           + 0.000387933_wp * T2 - T3/REAL(38710000, WP))*pic
!  between 0 and 2pi :
      IF (sidt .LT. 0._wp) THEN
        Call negangle2(pi2,sidt)
      ENDIF
      IF(sidt.Ge.pi2) THEN
        CALL langle2(pi2,sidt)
      ENDIF
!
      END SUBROUTINE sidt2

      Subroutine obliq(fnut,pic,T,ecl,nutob,nutl)
!  calculation of obliquity of ecliptic
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      INTEGER :: fnut
      REAL(wp) :: pic,T,ecl,nutob,nutl,A,B,C,T1,T2,T3
      REAL(wp) :: L1,L2,D1,D2,M1,M2,N1,N2
!
!  see page 57
! correction terms if requested
      If(fnut.Eq.1) Then
        t1 = t + 1._wp     ! adding one century to "shift" reference time to J1900
        T2=T1*T1
!
        a = 100.0021358_wp * t1
        b = 360._wp * (a - AINT(a, wp))
        l1 = 279.6967_wp + 0.000303_wp * t2 + b
        l2 = 2.0_wp * l1 * pic
!
        a = 1336.855231_wp * t1
        b = 360._wp * (a - AINT(a, wp))
        d1 = 270.4342_wp - 0.001133_wp * t2 + b
        d2 = 2.0_wp * d1 * pic
!
        a = 99.99736056_wp * t1
        b = 360._wp * (a - AINT(a, wp))
        m1 = (358.4758_wp - 0.00015_wp * t2 + b) * pic
!
        a = 1325.552359_wp * t1
        b = 360._wp * (a - AINT(a, wp))
        m2 = (296.1046_wp + 0.009192_wp * t2 + b) * pic
!
        a = 5.372616667_wp * t1
        b = 360._wp * (a - AINT(a, wp))
        n1 = (259.1833_wp + 0.002078_wp * t2 - b) * pic
        n2 = 2._wp * n1
! correction term for nutation in longitude
        nutl = ((-17.2327_wp - 0.01737_wp * t1) * SIN(n1)               &
             + (-1.2729_wp - 0.00013_wp * t1) * SIN(l2) + 0.2088_wp * SIN(n2) &
             - 0.2037_wp * SIN(d2) + (0.1261_wp - 0.00031_wp * t1) * SIN(m1)  &
             + 0.0675_wp * SIN(m2) - (0.0497_wp - 0.00012_wp * t1) * SIN(l2 + m1) &
             - 0.0342_wp * SIN(d2 - n1) - 0.0261_wp * SIN(d2 + m2)      &
             + 0.0214_wp * SIN(l2 - m1) - 0.0149_wp * SIN(l2 - d2 + m2) &
             + 0.0124_wp * SIN(l2 - n1) + 0.0114_wp * SIN(d2 - m2)) &
             / 3600.0_WP * pic
! correction term for nutation in obliquity of the ecliptic
        nutob = ((9.21_wp + 0.00091_wp * t1) * COS(n1)                  &
             + (0.5522_wp - 0.00029_wp * t1) * COS(l2) - 0.0904_wp * COS(n2) &
             + 0.0884_wp * COS(d2) + 0.0216_wp * COS(l2 + m1)           &
             + 0.0183_wp * COS(d2 - n1) + 0.0113_wp * COS(d2 + m2)      &
             - 0.0093_wp * COS(l2 - m1) - 0.0066_wp * COS(l2 - n1)) &
             / 3600.0_WP * pic
!
      Else
        nutob = 0.0_wp
        nutl = 0.0_wp
      Endif
! obliquity of the ecliptic
      t1 = t + 1._wp     ! adding one century to "shift" reference time to J1900
      T2=T1*T1
      T3=T2*T1
      c = 46.815_wp * t1 + 0.0006_wp * t2 - 0.00181_wp * t3
      ecl = (23.43929167_wp - c/3600.0_WP) * pic + nutob
!
      END SUBROUTINE obliq

      SUBROUTINE Sun_n(fnut,pic,pi2,T,ecl,nutl,res)
!  calculation of position of the Sun
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      Integer fnut
      REAL(wp) :: pic,pi2,T,T1,T2,T3,ecl,A,B,nutl
      REAL(wp) :: L,M1,EC,AT,AE
      REAL(wp) :: A1,B1,C1,D1,E1,H1,D2,D3,S1,S2,S3,SW,X1,X2,res(3,2)
!
!  see page 116
      t1 = t + 1._wp ! adding one century to "shift" reference time to J1900
      T2=T1*T1
      T3=T2*T1
!
      a = 100.0021359_wp * t1
      b = 360._wp * (a - AINT(a, wp))
      l = (279.69668_wp + 0.0003025_wp * t2 + b) * pic
!
      a = 99.99736042_wp * t1
      b = 360._wp * (a - AINT(a, wp))
      m1 = (358.47583_wp - 0.00015_wp * t2 + 0.0000033_wp * t3 + b) * pic
      ec = 0.01675104_wp - 0.0000418_wp * t1 - 0.000000126_wp * t2
!
!  true and eccentric anomaly in rad
      Call anomaly(pi2,M1,EC,AT,AE)
!
!  various arguments in rad
      a = 62.55209472_wp * t1
      b = 360._wp * (a - AINT(a, wp))
      a1 = (153.23_wp + b) * pic
!
      a = 125.1041894_wp * t1
      b = 360._wp * (a - AINT(a, wp))
      b1 = (216.57_wp + b) * pic
!
      a = 91.56766028_wp * t1
      b = 360._wp * (a - AINT(a, wp))
      c1 = (312.69_wp + b) * pic
!
      a = 1236.853095_wp * t1
      b = 360._wp * (a - AINT(a, wp))
      d1 = (350.74_wp - 0.00144_wp * t2 + b) * pic
      e1 = (231.19_wp + 20.2_wp * t1) * pic
!
      a = 183.1353208_wp * t1
      b = 360._wp * (a - AINT(a, wp))
      h1 = (353.4_wp + b) * pic
!
      d2 = (0.00134_wp * COS(a1) + 0.00154_wp * COS(b1) + 0.002_wp * COS(c1) &
           + 0.00179_wp * SIN(d1) + 0.00178_wp * SIN(e1)) * pic
      d3 = 0.00000543_wp * SIN(a1) + 0.00001575_wp * SIN(b1) &
           + 0.00001627_wp * SIN(c1) + 0.00003076_wp * COS(d1) &
           + 0.00000927_wp * SIN(h1)
!
!  geocentric ecliptic coordinates of the Sun
      S1=AT+L-M1+D2
      If(fnut.Eq.1) Then
        S1=S1+nutl
      Endif
      IF (s1 .LT. 0._wp) CALL negangle2(pi2,S1)
      If(S1.Ge.pi2) Call langle2(pi2,S1)
      s2 = 0.0_wp
      s3 = 1.0000002_wp * (1.0_wp - ec * COS(ae)) + d3
!
!  geocentric equatorial coordinates of the Sun
      SW  =  -1._wp
      Call eqecl(pi2,S1,S2,X1,X2,ecl,SW)
      res(1,1)=X1
      res(2,1)=X2
      res(3,1)=S3
!
      End Subroutine sun_n

      Subroutine moon(fnut,pic,pi2,T,ecl,nutl,res)
!  calculation of position of the Moon
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      Integer fnut
      REAL(wp) :: pic,pi2,T,T1,T2,T3,ecl,nutl,A,B,C,SW
      REAL(wp) :: Q,M1,M2,M3,M4,M5,M6,ML,MS,MD,ME,MF,NA,S1,S2,S3,S4,E,E2
      REAL(wp) :: L,G,W1,W2,PM,MO1,MO2,MO3,X1,X2,res(3,2)
!
!  see page 157
      t1 = t + 1._wp ! adding one century to "shift" reference time to J1900
      T2=T1*T1
      T3=T2*T1
!
      q = t1 * 36525._wp
      m1 = q/27.32158213_wp
      m1 = 360._wp * (m1 - AINT(m1, wp))
      m2 = q / 365.2596407_wp
      m2 = 360._wp * (m2 - AINT(m2, wp))
      m3 = q / 27.55455094_wp
      m3 = 360._wp * (m3 - AINT(m3, wp))
      m4 = q / 29.53058868_wp
      m4 = 360._wp * (m4 - AINT(m4, wp))
      m5 = q / 27.21222039_wp
      m5 = 360._wp * (m5 - AINT(m5, wp))
      m6 = q / 6798.363307_wp
      m6 = 360._wp * (m6 - AINT(m6, wp))
      ml = 270.434164_wp + m1 - 0.001133_wp * t2 + 0.0000019_wp * t3
      ms = 358.475833_wp + m2 - 0.00015_wp * t2 + 0.0000033_wp * t3
      md = 296.104608_wp + m3 + 0.009192_wp * t2 + 0.0000144_wp * t3
      me = 350.737486_wp + m4 - 0.001436_wp * t2 + 0.0000019_wp * t3
      mf = 11.250889_wp + m5 - 0.003211_wp * t2 - 0.0000003_wp * t3
      na = (259.183275_wp - m6 + 0.002078_wp * t2 + 0.0000022_wp * t3) * pic
      s2 = SIN(na)
      a = (51.2_wp + 20.2_wp * t1) * pic
      s1 = SIN(a)
      b = (346.56_wp + 132.87_wp * t1 - 0.0091731_wp * t2) * pic
      s3 = 0.003964_wp * SIN(b)
      c = na + (275.05_wp - 2.3_wp * t1) * pic
      s4 = SIN(c)
      ml = (ml + 0.000233_wp * s1 + s3 + 0.001964_wp * s2) * pic
      ms = (ms - 0.001778_wp * s1) * pic
      md = (md + 0.000817_wp * s1 + s3 + 0.002541_wp * s2) * pic
      mf = (mf + s3 - 0.024691_wp * s2 - 0.004328_wp * s4) * pic
      me = (me + 0.002011_wp * s1 + s3 + 0.001964_wp * s2) * pic
      e = 1._wp - 0.002495_wp * t1 + 0.00000752_wp * t2
      e2 = e * e

!  ecliptic longitude MO1
      l = 6.28875_wp * SIN(md) + 1.274018_wp * SIN(2._wp * me - md)        &
           + 0.658309_wp * SIN(2._wp * me) + 0.213616_wp * SIN(2._wp * md) &
           - e * 0.185596_wp * SIN(ms) - 0.114336_wp * SIN(2._wp * mf)     &
           + 0.058793_wp * SIN(2._wp * (me - md))                          &
           + 0.057212_wp * e * SIN(2._wp * me - ms - md) + 0.05332_wp      &
           & * SIN(2._wp * me + md)                                        &
           + 0.045874_wp * e * SIN(2._wp * me - ms)                        &
           + 0.041024_wp * e * SIN(md - ms)                                &
           - 0.034718_wp * SIN(me) - e * 0.030465_wp * SIN(md + ms)        &
           + 0.015326_wp * SIN(2._wp * (me - mf))                          &
           - 0.012528_wp * SIN(2._wp * mf + md)                            &
           - 0.01098_wp * SIN(2._wp * mf - md)                             &
           + 0.010674_wp * SIN(4._wp * me - md)                            &
           + 0.010034_wp * SIN(3._wp * md)                                 &
           + 0.008548_wp * SIN(4._wp * me - 2._wp * md)                    &
           - e * 0.00791_wp * SIN(ms - md + 2._wp * me)                    &
           - e * 0.006783_wp * SIN(2._wp * me + ms)                        &
           + 0.005162_wp * SIN(md - me) + e * 0.005_wp * SIN(me + ms)      &
           + 0.003862_wp * SIN(4._wp * me)                                 &
           + e * 0.004049_wp * SIN(md - ms + 2._wp * me)                   &
           + 0.003996_wp * SIN(2._wp * (md + me))                          &
           + 0.003665_wp * SIN(2._wp * me - 3._wp * md)                    &
           + e * 0.002695_wp * SIN(2._wp * md - ms)                        &
           + 0.002602_wp * SIN(md - 2._wp * (mf + me))                     &
           + e * 0.002396_wp * SIN(2._wp * (me - md) - ms)                 &
           - 0.002349_wp * SIN(me + md)                                    &
           + e2 * 0.002249_wp * SIN(2._wp * (me - ms))                     &
           - e * 0.002125_wp * SIN(ms + 2._wp * md)                        &
           - e2 * 0.002079_wp * SIN(2._wp * ms)                            &
           + e2 * 0.002059_wp * SIN(2._wp * (me - ms) - md)                &
           - 0.001773_wp * SIN(2._wp * (me - mf) + md)                     &
           - 0.001595_wp * SIN(2._wp * (me + mf))                          &
           + e * 0.00122_wp * SIN(4._wp * me - ms - md)                    &
           - 0.00111_wp * SIN(2._wp * (md + mf))                           &
           + 0.000892_wp * SIN(md - 3._wp * me)                            &
           - e * 0.000811_wp * SIN(ms + md + 2._wp * me)                   &
           + e * 0.000761_wp * SIN(4._wp * me - ms - 2._wp * md)           &
           + e2 * 0.000704_wp * SIN(md - 2._wp * (ms + me))                &
           + e * 0.000693_wp * SIN(ms - 2._wp * (md - me))                 &
           + e * 0.000598_wp * SIN(2._wp * (me - mf) - ms)                 &
           + 0.00055_wp * SIN(md + 4._wp * me)                             &
           + 0.000538_wp * SIN(4._wp * md)                                 &
           + e * 0.000521_wp * SIN(4._wp * me - ms)                        &
           + 0.000486_wp * SIN(2._wp * md - me)                            &
           + e2 * 0.000717_wp * SIN(md - 2._wp * ms)
      MO1=ML+L*pic
      If(fnut.Eq.1) Then
        MO1=MO1+nutl
      Endif
      IF (MO1 .LT. 0._wp) CALL negangle2(pi2,MO1)
      IF(MO1 .GE. pi2) CALL langle2(pi2,MO1)

!  ecliptic latitude MO2
      g = 5.128189_wp * SIN(mf) + 0.280606_wp * SIN(md + mf)                 &
           + 0.277693_wp * SIN(md - mf) + 0.173238_wp * SIN(2._wp * me - mf) &
           + 0.055413_wp * SIN(2._wp * me + mf - md)                         &
           + 0.046272_wp * SIN(2._wp * me - mf - md)                         &
           + 0.032573_wp * SIN(2._wp * me + mf)                              &
           + 0.017198_wp * SIN(2._wp * md + mf)                              &
           + 0.009267_wp * SIN(2._wp * me - mf + md)                         &
           + 0.008823_wp * SIN(2._wp * md - mf)                              &
           + e * 0.008247_wp * SIN(2._wp * me - ms - mf)                     &
           + 0.004323_wp * SIN(2._wp * (me + md) - mf)                       &
           + 0.0042_wp * SIN(2._wp * me + md + mf)                           &
           + e * 0.003372_wp * SIN(mf - ms - 2._wp * me)                     &
           + e * 0.002472_wp * SIN(2._wp * me - md + mf - ms)                &
           + e * 0.002222_wp * SIN(2._wp * me + mf - ms)                     &
           + e * 0.002072_wp * SIN(2._wp * me - md - mf - ms)                &
           + e * 0.001877_wp * SIN(mf - ms + md)                             &
           + 0.001828_wp * SIN(4._wp * me - md - mf)                         &
           - e * 0.001803_wp * SIN(ms + mf) - 0.00175_wp * SIN(3._wp * mf)   &
           + e * 0.00157_wp * SIN(md - mf - ms) - 0.001487_wp * SIN(me + mf) &
           - e * 0.001481_wp * SIN(mf + ms + md)                             &
           + e * 0.001417_wp * SIN(mf - ms - md)                             &
           + e * 0.00135_wp * SIN(mf - ms) + 0.00133_wp * SIN(mf - me)       &
           + 0.001106_wp * SIN(mf + 3._wp * md)                              &
           + 0.00102_wp * SIN(4._wp * me - mf)                               &
           + 0.000833_wp * SIN(mf + 4._wp * me - md)                         &
           + 0.000781_wp * SIN(md - 3._wp * mf)                              &
           + 0.00067_wp * SIN(mf + 3._wp * me - 2._wp * md)                  &
           + 0.000606_wp * SIN(2._wp * me - 3._wp * mf)                      &
           + 0.000597_wp * SIN(2._wp * (me + md) - mf)                       &
           + e * 0.000492_wp * SIN(2._wp * me + md - ms - mf)                &
           + 0.00045_wp * SIN(2._wp * (md - me) - mf)                        &
           + 0.000439_wp * SIN(3._wp * me - mf)                              &
           + 0.000423_wp * SIN(mf + 2._wp * (me + md))                       &
           + 0.000422_wp * SIN(2._wp * me - 3._wp * md - mf)                 &
           - e * 0.000367_wp * SIN(mf + ms + 2._wp * me - md)                &
           - e * 0.000353_wp * SIN(mf + ms + 2._wp * me)                     &
           + 0.000331_wp * SIN(mf + 4._wp * me)                              &
           + e * 0.000317_wp * SIN(2._wp * me + md - ms + mf)                &
           + e2 * 0.000306_wp * SIN(2._wp * (me - ms) - mf)                  &
           - 0.000283_wp *SIN(md + 3._wp * mf)
      w1 = 0.0004664_wp * COS(na)
      w2 = 0.0000754_wp * COS(c)
      mo2 = g * pic * (1.0_wp - w1 - w2)

!  horizontal parallax PM
      pm = 0.950724_wp + 0.051818_wp * COS(md)                             &
           + 0.009531_wp * COS(2._wp * me - md)                            &
           + 0.007843_wp * COS(2._wp * me) + 0.002824_wp * COS(2._wp * md) &
           + 0.000857_wp * COS(2._wp * me + md)                            &
           + e * 0.000533_wp * COS(2._wp * me - ms)                        &
           + e * 0.000401_wp * COS(2._wp * me - md - ms)                   &
           + e * 0.00032_wp * COS(md - ms) - 0.000271_wp * COS(me)         &
           - e * 0.000264_wp * COS(md + ms)                                &
           - 0.000198_wp * COS(2._wp * mf - md)                            &
           + 0.000173_wp * COS(3._wp * md)                                 &
           + 0.000167_wp * COS(4._wp * me - md)- e * 0.000111_wp * COS(ms) &
           + 0.000103_wp * COS(4._wp * me - 2._wp * md)                    &
           - 0.000084_wp * COS(2._wp * md - 2._wp * me)                    &
           - e * 0.000083_wp * COS(2._wp * me + ms)                        &
           + 0.000079_wp * COS(2._wp * me + 2._wp * md)                    &
           + 0.000072_wp * COS(4._wp * me)                                 &
           + e * 0.000064_wp * COS(2._wp * me - ms + md)                   &
           - e * 0.000063_wp * COS(2._wp * me + ms - md)                   &
           + e * 0.000041_wp * COS(ms + me)                                &
           + e * 0.000035_wp * COS(2._wp * md - ms)                        &
           - 0.000033_wp * COS(3._wp * md - 2._wp * me)                    &
           - 0.00003_wp * COS(md + me)                                     &
           - 0.000029_wp * COS(2._wp * (mf - me))                          &
           - e * 0.000029_wp * COS(2._wp * md + ms)                        &
           + e2 * 0.000026_wp * COS(2._wp * (me - ms))                     &
           - 0.000023_wp * COS(2._wp * (mf - me) + md)                     &
           + e * 0.000019_wp * COS(4._wp * me - md - ms)
      PM=PM*pic

!  geocentric distance MO3 in km
      mo3 = 6378.14_wp/SIN(pm)

!  geocentric equatorial coordinates of the Moon
      sw = -1._wp
      Call eqecl(pi2,MO1,MO2,X1,X2,ecl,SW)
      res(1,2)=X1
      res(2,2)=X2
      res(3,2)=MO3
!
      END SUBROUTINE moon

      SUBROUTINE anomaly(pi2,AM,EC,AT,AE)
!  calculation of true anomaly AT and eccentric anomaly AE
!  given mean anomaly AM and eccentricity EC
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      REAL(wp) :: pi2,AM,EC,AT,AE,M,D,A
!
!  see page 113
      m = am - pi2 * AINT((am / pi2), wp)
      AE=M
  1   D=AE-(EC*Sin(AE))-M
      IF (ABS(d) .GE. 0.000006_wp) THEN
        d = d/(1.0_wp - ec * COS(ae))
        AE=AE-D
        Goto 1
      Else
        a = SQRT((1.0_wp + ec)/(1.0_wp - ec))*TAN(ae/2.0_wp)
        at = 2.0_wp * ATAN(a)
      Endif
!
      END SUBROUTINE anomaly

      SUBROUTINE eqecl(pi2,X,Y,P,Q,ecl,SW)
!  conversion of ecliptic into equatorial coordinates
!  according to Duffett, 1990
!
! if SW=+1: equatorial (X,Y..alpha,delta) to ecliptic (P,Q..lambda,beta)
! if SW=-1: equatorial (X,Y..lambda,beta) to ecliptic (P,Q..alpha,delta)
!
      IMPLICIT NONE
!
      REAL(wp) :: pi2,ecl,P,Q,X,Y,SW
!
!  see page 62
      p = ATAN2((SIN(x) * COS(ecl) + TAN(y) * SIN(ecl) * sw), COS(X))
      IF (p .LT. 0._wp) CALL negangle2(pi2,p)
      IF (P .GE. pi2) CALL langle2(pi2,P)
      Q=Asin(Sin(Y)*Cos(ecl)-Cos(Y)*Sin(ecl)*Sin(X)*SW)
!
      END SUBROUTINE eqecl

      SUBROUTINE aufb2(sidt,res,res2)
!
! modifications according to "ephaufb.f" by Maik Thomas
! (rekt(rad)->sid.time.green.-r.asc.; dekl(rad); cri3->(a/r)^3
! for Sun and Moon)
!
      IMPLICIT NONE
!
      REAL(wp) :: sidt,h(3)
      REAL(wp) :: res(3,2),res2(3,2)
!
      res2(1,1)=sidt-res(1,1)
      res2(1,2)=sidt-res(1,2)
      res2(2,1)=res(2,1)
      res2(2,2)=res(2,2)
      h(1) = 1._wp / res(3,1)
      h(2) = 384400._wp / res(3,2)
      res2(3,1)=h(1)*h(1)*h(1)
      res2(3,2)=h(2)*h(2)*h(2)
!
      Return
      END SUBROUTINE aufb2

      SUBROUTINE negangle2(pi2,x)
!  transformation of negative angles
!  to angles in interval [0\B0;360\B0)
!
      IMPLICIT NONE
!
      Logical endvar
      REAL(wp) :: pi2,x
!
      endvar=.False.
!
      Do while (.Not.endvar)
        x=x+pi2
        endvar = (x .GE. 0._wp)
      Enddo
!
      Return
      END  SUBROUTINE  negangle2

      SUBROUTINE langle2(pi2,x)
!  transformation of large angles
!  to angles in interval [0\B0;360\B0)
!
      IMPLICIT NONE
!
      Logical endvar
      REAL(wp) :: pi2,x
!
      endvar=.False.
!
      Do while (.Not.endvar)
        x=x-pi2
        endvar=(x.Lt.pi2)
      Enddo
!
      Return
      END SUBROUTINE    langle2

   END MODULE mo_tidal


