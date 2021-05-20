!> Tidal Module
!> Real time computation of luni-solar gravitational potential
!> @author Maik Thomas, Ann-Margrit Hellmich
!> @date Last modified 9.9.2009 Malte Mueller
!! Included subroutines are:
!! alloc_mem_tidal   -- Allocation of memory
!! foreph_ini        -- Initialization of tidal module
!!    eph            -- determination of the phase with respect to 1.1.2000
!! foreph            -- main routine for calculation of realtime lunisolar
!!                      gravitational potential, with effective Earth Elasticity
!!                      Factor (EEF) of Body Earth Tide considered (Kantha-1995-JGR)
!!             FIXME::  Self-Attraction and Loading effect need to be fixed ?
!!    ephvsop87      -- ephemeredes of sun and moon
!!        sidt2      -- computes sidereal time
!!        obliq      -- calculation of obliquity of ecliptic
!!        Sun_n      -- calculation of position of sun
!!          anomaly  -- calculation of true anomaly AT and eccentric anomaly AE
!!        moon      -- calculation of position of moon
!!        aufb2
!!    eqecl          -- conversion of ecliptic into equatorial coordinates


MODULE mo_tidal
      USE o_PARAM
      USE o_ARRAYS, only : ssh_gp
      USE mod_mesh
      USE g_config, only : dt
      USE g_PARSUP
      USE g_clock
      IMPLICIT NONE
      !Earth Tides ( maik thomas, emr  pers. comm )
      REAL(WP), PARAMETER :: eef=0.69_WP

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
      mmccdt = (jcc + moph - 0.5_WP) * NINT(86400._WP / dt)
      ! mmccdt is the number of timesteps since 01.01.2000 00:00:00
      ! bugfix by M.Mueller : substracting 0.5_WP shifts the reference time by -12h to 01.01.2000 12:00:00
      ! This is needed to get the correct sideral time (sbr sidt2)
      ! and calculate the correct hour angle for the tidal potential (->sbr aufb2)
      ! FIXME : mmccdt = ldtmon-1 + int(jcc + moph - 0.5_WP) * NINT(86400._WP / dt)
      ! FIXME : replace eph by a some to code that directly calculates julian days and
      ! centuries as needed by siderial time and ephemerides

      if (mype==0) WRITE(*,*)'tidal: phase relative to 2000 :'    &
      ,'year= ',lyear, 'month= ',lmonth, 'yearoff= ',jcc,' monoff= ',moph ,'mmccdt= ',mmccdt

      END SUBROUTINE foreph_ini

      SUBROUTINE eph(jul,mon,jahrph,moph)
      !    Berechnet Phase [in Tagen]  bezueglich des 1.1.2000
      !       jahrph=phase zu beginn des Jahres
      !       gesamtphase=jahrph+moph

      IMPLICIT NONE

      INTEGER,INTENT(IN)::jul,mon
      INTEGER,INTENT(OUT)::jahrph,moph
      INTEGER::jpl,jj,m,iflag

      ! Before year 2000
      IF(jul.LT.2000) THEN
        jahrph=0
        DO jj=jul,1999
          call check_fleapyr(jj,iflag)
          jpl=365+iflag

          jahrph=jahrph-jpl
        ENDDO
      ENDIF

      ! year 2000
      IF(jul.EQ.2000) jahrph=0

      ! After year 2000
      IF(jul.GT.2000) THEN
        jahrph=0
        DO jj=2000,jul-1
          call check_fleapyr(jj,iflag)
          jpl=365+iflag

          jahrph=jahrph+jpl
        ENDDO
      ENDIF

      moph=0
      IF (mon.GT.1)THEN
        DO m=1,mon-1
          moph=moph+num_day_in_month(fleapyear,m)
        ENDDO
      ENDIF

      END SUBROUTINE eph

      SUBROUTINE foreph(mesh)
!     calculates the realtime gravitational potential of sun & moon
!     output: ssh_gp (with Body Earth Tide effect)

      IMPLICIT NONE
      type(t_mesh), intent(in) , target :: mesh
      REAL(WP) :: dres(3,2),crim3,rkomp,erdrad,rekts,dekls
      REAL(WP) :: cris3,rektm,deklm,deklm2,dekls2,sidm,sidmq
      REAL(WP) :: rkosp,codm,codmq,sids,sidsq,cods,codsq,sidm2
      REAL(WP) :: sids2,hamp,hasp
      INTEGER :: i,j
#include  "associate_mesh.h"

      mmccdt = mmccdt + 1

      rkomp = -4.113e-07_WP ! factor of the tidal potential due to the moon
                            ! attention the factor is defined negative (contrary to the standard).
      rkosp = 0.46051_WP * rkomp
      erdrad = 6371000._WP

      CALL ephvsop87(dres)

      rekts=dres(1,1)
      dekls=dres(2,1)
      cris3=dres(3,1)

      rektm=dres(1,2)
      deklm=dres(2,2)
      crim3=dres(3,2)


      deklm2 = deklm * 2._WP
      dekls2 = dekls * 2._WP
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
         ssh_gp(i) = eef * erdrad * rkomp * crim3 &
              * (3._WP * (SIN(geo_coord_nod2D(2,i))**2 - 1._WP/3._WP) * (sidmq - 1._WP/3._WP)&
              &  + SIN(2._WP * geo_coord_nod2D(2,i)) * sidm2 * COS(hamp) &
              &  + COS(geo_coord_nod2D(2,i))**2 * codmq * COS(2._WP * hamp))        &
              &  + erdrad * rkosp * cris3 &
              &    * (3._WP * (SIN(geo_coord_nod2D(2,i))**2 - 1._WP/3._WP) &
              &       * (sidsq - 1._WP/3._WP) &
              &       + SIN(2._WP * geo_coord_nod2D(2,i)) * sids2 * COS(hasp) &
              &       + COS(geo_coord_nod2D(2,i))**2 * codsq * COS(2._WP * hasp))

      END DO

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
      REAL(WP) :: T,sidt,ecl,nutob,nutl,res(3,2),res2(3,2)

!     Transformation of Time into fractional julian centuries T
      t = (REAL(mmccdt-1, WP) * dt / 86400._WP)/36525._WP

!     corresponding sidereal time Greenwich
      Call sidt2(T,sidt)

!  set fnut (perform nutation -> 1; don't -> 0)
      fnut=0

! C obliquity of the ecliptic
      CALL obliq(fnut,T,ecl,nutob,nutl)

! C Sun
      CALL sun_n(fnut,T,ecl,nutl,res)

! C Moon
      CALL moon(fnut,T,ecl,nutl,res)

! C modifications in preparaion of calculation of potentials
      CALL aufb2(sidt,res,res2)
!
      END SUBROUTINE ephvsop87

      SUBROUTINE sidt2(T,sidt)
!  convertion of Julian Centuries since J2000 T to siderical time sidt
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      REAL(WP) :: T,sidt,JD,T2,T3
!
      T2=T*T
      T3=T2*T
!
!  Julian days
      jd = t * 36525.0_WP + 2451545.0_WP
!
!  mean siderial time of Greenwich in rad
      sidt = (280.46061837_WP + 360.98564736629_WP * (jd - 2451545.0_WP) &
           + 0.000387933_WP * T2 - T3/REAL(38710000, WP))*rad
!  between 0 and 2pi :
      IF (sidt .LT.  0._WP) THEN
        Call negangle2(sidt)
      ENDIF
      IF (sidt .Ge. (pi*2)) THEN
        CALL langle2(sidt)
      ENDIF
!
      END SUBROUTINE sidt2

      Subroutine obliq(fnut,T,ecl,nutob,nutl)
!  calculation of obliquity of ecliptic
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      INTEGER :: fnut
      REAL(WP) :: T,ecl,nutob,nutl,A,B,C,T1,T2,T3
      REAL(WP) :: L1,L2,D1,D2,M1,M2,N1,N2
!
!  see page 57
! correction terms if requested
      If(fnut.Eq.1) Then
        t1 = t + 1._WP     ! adding one century to "shift" reference time to J1900
        T2=T1*T1
!
        a = 100.0021358_WP * t1
        b = 360._WP * (a - AINT(a, WP))
        l1 = 279.6967_WP + 0.000303_WP * t2 + b
        l2 = 2.0_WP * l1 * rad
!
        a = 1336.855231_WP * t1
        b = 360._WP * (a - AINT(a, WP))
        d1 = 270.4342_WP - 0.001133_WP * t2 + b
        d2 = 2.0_WP * d1 * rad
!
        a = 99.99736056_WP * t1
        b = 360._WP * (a - AINT(a, WP))
        m1 = (358.4758_WP - 0.00015_WP * t2 + b) * rad
!
        a = 1325.552359_WP * t1
        b = 360._WP * (a - AINT(a, WP))
        m2 = (296.1046_WP + 0.009192_WP * t2 + b) * rad
!
        a = 5.372616667_WP * t1
        b = 360._WP * (a - AINT(a, WP))
        n1 = (259.1833_WP + 0.002078_WP * t2 - b) * rad
        n2 = 2._WP * n1
! correction term for nutation in longitude
        nutl = ((-17.2327_WP - 0.01737_WP * t1) * SIN(n1)               &
             + (-1.2729_WP - 0.00013_WP * t1) * SIN(l2) + 0.2088_WP * SIN(n2) &
             - 0.2037_WP * SIN(d2) + (0.1261_WP - 0.00031_WP * t1) * SIN(m1)  &
             + 0.0675_WP * SIN(m2) - (0.0497_WP - 0.00012_WP * t1) * SIN(l2 + m1) &
             - 0.0342_WP * SIN(d2 - n1) - 0.0261_WP * SIN(d2 + m2)      &
             + 0.0214_WP * SIN(l2 - m1) - 0.0149_WP * SIN(l2 - d2 + m2) &
             + 0.0124_WP * SIN(l2 - n1) + 0.0114_WP * SIN(d2 - m2)) &
             / 3600.0_WP * rad
! correction term for nutation in obliquity of the ecliptic
        nutob = ((9.21_WP + 0.00091_WP * t1) * COS(n1)                  &
             + (0.5522_WP - 0.00029_WP * t1) * COS(l2) - 0.0904_WP * COS(n2) &
             + 0.0884_WP * COS(d2) + 0.0216_WP * COS(l2 + m1)           &
             + 0.0183_WP * COS(d2 - n1) + 0.0113_WP * COS(d2 + m2)      &
             - 0.0093_WP * COS(l2 - m1) - 0.0066_WP * COS(l2 - n1)) &
             / 3600.0_WP * rad
!
      Else
        nutob = 0.0_WP
        nutl = 0.0_WP
      Endif
! obliquity of the ecliptic
      t1 = t + 1._WP     ! adding one century to "shift" reference time to J1900
      T2=T1*T1
      T3=T2*T1
      c = 46.815_WP * t1 + 0.0006_WP * t2 - 0.00181_WP * t3
      ecl = (23.43929167_WP - c/3600.0_WP) * rad + nutob
!
      END SUBROUTINE obliq

      SUBROUTINE Sun_n(fnut,T,ecl,nutl,res)
!  calculation of position of the Sun
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      Integer fnut
      REAL(WP) :: T,T1,T2,T3,ecl,A,B,nutl
      REAL(WP) :: L,M1,EC,AT,AE
      REAL(WP) :: A1,B1,C1,D1,E1,H1,D2,D3,S1,S2,S3,SW,X1,X2,res(3,2)
!
!  see page 116
      t1 = t + 1._WP ! adding one century to "shift" reference time to J1900
      T2=T1*T1
      T3=T2*T1
!
      a = 100.0021359_WP * t1
      b = 360._WP * (a - AINT(a, WP))
      l = (279.69668_WP + 0.0003025_WP * t2 + b) * rad
!
      a = 99.99736042_WP * t1
      b = 360._WP * (a - AINT(a, WP))
      m1 = (358.47583_WP - 0.00015_WP * t2 + 0.0000033_WP * t3 + b) * rad
      ec = 0.01675104_WP - 0.0000418_WP * t1 - 0.000000126_WP * t2
!
!  true and eccentric anomaly in rad
      Call anomaly(M1,EC,AT,AE)
!
!  various arguments in rad
      a = 62.55209472_WP * t1
      b = 360._WP * (a - AINT(a, WP))
      a1 = (153.23_WP + b) * rad
!
      a = 125.1041894_WP * t1
      b = 360._WP * (a - AINT(a, WP))
      b1 = (216.57_WP + b) * rad
!
      a = 91.56766028_WP * t1
      b = 360._WP * (a - AINT(a, WP))
      c1 = (312.69_WP + b) * rad
!
      a = 1236.853095_WP * t1
      b = 360._WP * (a - AINT(a, WP))
      d1 = (350.74_WP - 0.00144_WP * t2 + b) * rad
      e1 = (231.19_WP + 20.2_WP * t1) * rad
!
      a = 183.1353208_WP * t1
      b = 360._WP * (a - AINT(a, WP))
      h1 = (353.4_WP + b) * rad
!
      d2 = (0.00134_WP * COS(a1) + 0.00154_WP * COS(b1) + 0.002_WP * COS(c1) &
           + 0.00179_WP * SIN(d1) + 0.00178_WP * SIN(e1)) * rad
      d3 = 0.00000543_WP * SIN(a1) + 0.00001575_WP * SIN(b1) &
           + 0.00001627_WP * SIN(c1) + 0.00003076_WP * COS(d1) &
           + 0.00000927_WP * SIN(h1)
!
!  geocentric ecliptic coordinates of the Sun
      S1=AT+L-M1+D2
      If(fnut.Eq.1) Then
        S1=S1+nutl
      Endif
      IF (s1 .LT.  0._WP) CALL negangle2(S1)
      If (S1 .Ge. (pi*2)) Call langle2(S1)
      s2 = 0.0_WP
      s3 = 1.0000002_WP * (1.0_WP - ec * COS(ae)) + d3
!
!  geocentric equatorial coordinates of the Sun
      SW  =  -1._WP
      Call eqecl(S1,S2,X1,X2,ecl,SW)
      res(1,1)=X1
      res(2,1)=X2
      res(3,1)=S3
!
      End Subroutine sun_n

      Subroutine moon(fnut,T,ecl,nutl,res)
!  calculation of position of the Moon
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      Integer fnut
      REAL(WP) :: T,T1,T2,T3,ecl,nutl,A,B,C,SW
      REAL(WP) :: Q,M1,M2,M3,M4,M5,M6,ML,MS,MD,ME,MF,NA,S1,S2,S3,S4,E,E2
      REAL(WP) :: L,G,W1,W2,PM,MO1,MO2,MO3,X1,X2,res(3,2)
!
!  see page 157
      t1 = t + 1._WP ! adding one century to "shift" reference time to J1900
      T2=T1*T1
      T3=T2*T1
!
      q = t1 * 36525._WP
      m1 = q/27.32158213_WP
      m1 = 360._WP * (m1 - AINT(m1, WP))
      m2 = q / 365.2596407_WP
      m2 = 360._WP * (m2 - AINT(m2, WP))
      m3 = q / 27.55455094_WP
      m3 = 360._WP * (m3 - AINT(m3, WP))
      m4 = q / 29.53058868_WP
      m4 = 360._WP * (m4 - AINT(m4, WP))
      m5 = q / 27.21222039_WP
      m5 = 360._WP * (m5 - AINT(m5, WP))
      m6 = q / 6798.363307_WP
      m6 = 360._WP * (m6 - AINT(m6, WP))
      ml = 270.434164_WP + m1 - 0.001133_WP * t2 + 0.0000019_WP * t3
      ms = 358.475833_WP + m2 - 0.00015_WP * t2 + 0.0000033_WP * t3
      md = 296.104608_WP + m3 + 0.009192_WP * t2 + 0.0000144_WP * t3
      me = 350.737486_WP + m4 - 0.001436_WP * t2 + 0.0000019_WP * t3
      mf = 11.250889_WP + m5 - 0.003211_WP * t2 - 0.0000003_WP * t3
      na = (259.183275_WP - m6 + 0.002078_WP * t2 + 0.0000022_WP * t3) * rad
      s2 = SIN(na)
      a = (51.2_WP + 20.2_WP * t1) * rad
      s1 = SIN(a)
      b = (346.56_WP + 132.87_WP * t1 - 0.0091731_WP * t2) * rad
      s3 = 0.003964_WP * SIN(b)
      c = na + (275.05_WP - 2.3_WP * t1) * rad
      s4 = SIN(c)
      ml = (ml + 0.000233_WP * s1 + s3 + 0.001964_WP * s2) * rad
      ms = (ms - 0.001778_WP * s1) * rad
      md = (md + 0.000817_WP * s1 + s3 + 0.002541_WP * s2) * rad
      mf = (mf + s3 - 0.024691_WP * s2 - 0.004328_WP * s4) * rad
      me = (me + 0.002011_WP * s1 + s3 + 0.001964_WP * s2) * rad
      e = 1._WP - 0.002495_WP * t1 + 0.00000752_WP * t2
      e2 = e * e

!  ecliptic longitude MO1
      l = 6.28875_WP * SIN(md) + 1.274018_WP * SIN(2._WP * me - md)        &
           + 0.658309_WP * SIN(2._WP * me) + 0.213616_WP * SIN(2._WP * md) &
           - e * 0.185596_WP * SIN(ms) - 0.114336_WP * SIN(2._WP * mf)     &
           + 0.058793_WP * SIN(2._WP * (me - md))                          &
           + 0.057212_WP * e * SIN(2._WP * me - ms - md) + 0.05332_WP      &
           & * SIN(2._WP * me + md)                                        &
           + 0.045874_WP * e * SIN(2._WP * me - ms)                        &
           + 0.041024_WP * e * SIN(md - ms)                                &
           - 0.034718_WP * SIN(me) - e * 0.030465_WP * SIN(md + ms)        &
           + 0.015326_WP * SIN(2._WP * (me - mf))                          &
           - 0.012528_WP * SIN(2._WP * mf + md)                            &
           - 0.01098_WP * SIN(2._WP * mf - md)                             &
           + 0.010674_WP * SIN(4._WP * me - md)                            &
           + 0.010034_WP * SIN(3._WP * md)                                 &
           + 0.008548_WP * SIN(4._WP * me - 2._WP * md)                    &
           - e * 0.00791_WP * SIN(ms - md + 2._WP * me)                    &
           - e * 0.006783_WP * SIN(2._WP * me + ms)                        &
           + 0.005162_WP * SIN(md - me) + e * 0.005_WP * SIN(me + ms)      &
           + 0.003862_WP * SIN(4._WP * me)                                 &
           + e * 0.004049_WP * SIN(md - ms + 2._WP * me)                   &
           + 0.003996_WP * SIN(2._WP * (md + me))                          &
           + 0.003665_WP * SIN(2._WP * me - 3._WP * md)                    &
           + e * 0.002695_WP * SIN(2._WP * md - ms)                        &
           + 0.002602_WP * SIN(md - 2._WP * (mf + me))                     &
           + e * 0.002396_WP * SIN(2._WP * (me - md) - ms)                 &
           - 0.002349_WP * SIN(me + md)                                    &
           + e2 * 0.002249_WP * SIN(2._WP * (me - ms))                     &
           - e * 0.002125_WP * SIN(ms + 2._WP * md)                        &
           - e2 * 0.002079_WP * SIN(2._WP * ms)                            &
           + e2 * 0.002059_WP * SIN(2._WP * (me - ms) - md)                &
           - 0.001773_WP * SIN(2._WP * (me - mf) + md)                     &
           - 0.001595_WP * SIN(2._WP * (me + mf))                          &
           + e * 0.00122_WP * SIN(4._WP * me - ms - md)                    &
           - 0.00111_WP * SIN(2._WP * (md + mf))                           &
           + 0.000892_WP * SIN(md - 3._WP * me)                            &
           - e * 0.000811_WP * SIN(ms + md + 2._WP * me)                   &
           + e * 0.000761_WP * SIN(4._WP * me - ms - 2._WP * md)           &
           + e2 * 0.000704_WP * SIN(md - 2._WP * (ms + me))                &
           + e * 0.000693_WP * SIN(ms - 2._WP * (md - me))                 &
           + e * 0.000598_WP * SIN(2._WP * (me - mf) - ms)                 &
           + 0.00055_WP * SIN(md + 4._WP * me)                             &
           + 0.000538_WP * SIN(4._WP * md)                                 &
           + e * 0.000521_WP * SIN(4._WP * me - ms)                        &
           + 0.000486_WP * SIN(2._WP * md - me)                            &
           + e2 * 0.000717_WP * SIN(md - 2._WP * ms)
      MO1=ML+L*rad
      If(fnut.Eq.1) Then
        MO1=MO1+nutl
      Endif
      IF (MO1 .LT.  0._WP) CALL negangle2(MO1)
      IF (MO1 .GE. (pi*2)) CALL langle2(MO1)

!  ecliptic latitude MO2
      g = 5.128189_WP * SIN(mf) + 0.280606_WP * SIN(md + mf)                 &
           + 0.277693_WP * SIN(md - mf) + 0.173238_WP * SIN(2._WP * me - mf) &
           + 0.055413_WP * SIN(2._WP * me + mf - md)                         &
           + 0.046272_WP * SIN(2._WP * me - mf - md)                         &
           + 0.032573_WP * SIN(2._WP * me + mf)                              &
           + 0.017198_WP * SIN(2._WP * md + mf)                              &
           + 0.009267_WP * SIN(2._WP * me - mf + md)                         &
           + 0.008823_WP * SIN(2._WP * md - mf)                              &
           + e * 0.008247_WP * SIN(2._WP * me - ms - mf)                     &
           + 0.004323_WP * SIN(2._WP * (me + md) - mf)                       &
           + 0.0042_WP * SIN(2._WP * me + md + mf)                           &
           + e * 0.003372_WP * SIN(mf - ms - 2._WP * me)                     &
           + e * 0.002472_WP * SIN(2._WP * me - md + mf - ms)                &
           + e * 0.002222_WP * SIN(2._WP * me + mf - ms)                     &
           + e * 0.002072_WP * SIN(2._WP * me - md - mf - ms)                &
           + e * 0.001877_WP * SIN(mf - ms + md)                             &
           + 0.001828_WP * SIN(4._WP * me - md - mf)                         &
           - e * 0.001803_WP * SIN(ms + mf) - 0.00175_WP * SIN(3._WP * mf)   &
           + e * 0.00157_WP * SIN(md - mf - ms) - 0.001487_WP * SIN(me + mf) &
           - e * 0.001481_WP * SIN(mf + ms + md)                             &
           + e * 0.001417_WP * SIN(mf - ms - md)                             &
           + e * 0.00135_WP * SIN(mf - ms) + 0.00133_WP * SIN(mf - me)       &
           + 0.001106_WP * SIN(mf + 3._WP * md)                              &
           + 0.00102_WP * SIN(4._WP * me - mf)                               &
           + 0.000833_WP * SIN(mf + 4._WP * me - md)                         &
           + 0.000781_WP * SIN(md - 3._WP * mf)                              &
           + 0.00067_WP * SIN(mf + 3._WP * me - 2._WP * md)                  &
           + 0.000606_WP * SIN(2._WP * me - 3._WP * mf)                      &
           + 0.000597_WP * SIN(2._WP * (me + md) - mf)                       &
           + e * 0.000492_WP * SIN(2._WP * me + md - ms - mf)                &
           + 0.00045_WP * SIN(2._WP * (md - me) - mf)                        &
           + 0.000439_WP * SIN(3._WP * me - mf)                              &
           + 0.000423_WP * SIN(mf + 2._WP * (me + md))                       &
           + 0.000422_WP * SIN(2._WP * me - 3._WP * md - mf)                 &
           - e * 0.000367_WP * SIN(mf + ms + 2._WP * me - md)                &
           - e * 0.000353_WP * SIN(mf + ms + 2._WP * me)                     &
           + 0.000331_WP * SIN(mf + 4._WP * me)                              &
           + e * 0.000317_WP * SIN(2._WP * me + md - ms + mf)                &
           + e2 * 0.000306_WP * SIN(2._WP * (me - ms) - mf)                  &
           - 0.000283_WP *SIN(md + 3._WP * mf)
      w1 = 0.0004664_WP * COS(na)
      w2 = 0.0000754_WP * COS(c)
      mo2 = g * rad * (1.0_WP - w1 - w2)

!  horizontal parallax PM
      pm = 0.950724_WP + 0.051818_WP * COS(md)                             &
           + 0.009531_WP * COS(2._WP * me - md)                            &
           + 0.007843_WP * COS(2._WP * me) + 0.002824_WP * COS(2._WP * md) &
           + 0.000857_WP * COS(2._WP * me + md)                            &
           + e * 0.000533_WP * COS(2._WP * me - ms)                        &
           + e * 0.000401_WP * COS(2._WP * me - md - ms)                   &
           + e * 0.00032_WP * COS(md - ms) - 0.000271_WP * COS(me)         &
           - e * 0.000264_WP * COS(md + ms)                                &
           - 0.000198_WP * COS(2._WP * mf - md)                            &
           + 0.000173_WP * COS(3._WP * md)                                 &
           + 0.000167_WP * COS(4._WP * me - md)- e * 0.000111_WP * COS(ms) &
           + 0.000103_WP * COS(4._WP * me - 2._WP * md)                    &
           - 0.000084_WP * COS(2._WP * md - 2._WP * me)                    &
           - e * 0.000083_WP * COS(2._WP * me + ms)                        &
           + 0.000079_WP * COS(2._WP * me + 2._WP * md)                    &
           + 0.000072_WP * COS(4._WP * me)                                 &
           + e * 0.000064_WP * COS(2._WP * me - ms + md)                   &
           - e * 0.000063_WP * COS(2._WP * me + ms - md)                   &
           + e * 0.000041_WP * COS(ms + me)                                &
           + e * 0.000035_WP * COS(2._WP * md - ms)                        &
           - 0.000033_WP * COS(3._WP * md - 2._WP * me)                    &
           - 0.00003_WP * COS(md + me)                                     &
           - 0.000029_WP * COS(2._WP * (mf - me))                          &
           - e * 0.000029_WP * COS(2._WP * md + ms)                        &
           + e2 * 0.000026_WP * COS(2._WP * (me - ms))                     &
           - 0.000023_WP * COS(2._WP * (mf - me) + md)                     &
           + e * 0.000019_WP * COS(4._WP * me - md - ms)
      PM=PM*rad

!  geocentric distance MO3 in km
      mo3 = 6378.14_WP/SIN(pm)

!  geocentric equatorial coordinates of the Moon
      sw = -1._WP
      Call eqecl(MO1,MO2,X1,X2,ecl,SW)
      res(1,2)=X1
      res(2,2)=X2
      res(3,2)=MO3
!
      END SUBROUTINE moon

      SUBROUTINE anomaly(AM,EC,AT,AE)
!  calculation of true anomaly AT and eccentric anomaly AE
!  given mean anomaly AM and eccentricity EC
!  according to Duffett, 1990
!
      IMPLICIT NONE
!
      REAL(WP) :: AM,EC,AT,AE,M,D,A
!
!  see page 113
      m = am - pi*2 * AINT((am / (pi*2)), WP)
      AE=M
  1   D=AE-(EC*Sin(AE))-M
      IF (ABS(d) .GE. 0.000006_WP) THEN
        d = d/(1.0_WP - ec * COS(ae))
        AE=AE-D
        Goto 1
      Else
        a = SQRT((1.0_WP + ec)/(1.0_WP - ec))*TAN(ae/2.0_WP)
        at = 2.0_WP * ATAN(a)
      Endif
!
      END SUBROUTINE anomaly

      SUBROUTINE eqecl(X,Y,P,Q,ecl,SW)
!  conversion of ecliptic into equatorial coordinates
!  according to Duffett, 1990
!
! if SW=+1: equatorial (X,Y..alpha,delta) to ecliptic (P,Q..lambda,beta)
! if SW=-1: equatorial (X,Y..lambda,beta) to ecliptic (P,Q..alpha,delta)
!
      IMPLICIT NONE
!
      REAL(WP) :: ecl,P,Q,X,Y,SW
!
!  see page 62
      p = ATAN2((SIN(x) * COS(ecl) + TAN(y) * SIN(ecl) * sw), COS(X))
      IF (p .LT.  0._WP) CALL negangle2(p)
      IF (P .GE. (pi*2)) CALL langle2(P)
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
      REAL(WP) :: sidt,h(3)
      REAL(WP) :: res(3,2),res2(3,2)
!
      res2(1,1)=sidt-res(1,1)
      res2(1,2)=sidt-res(1,2)
      res2(2,1)=res(2,1)
      res2(2,2)=res(2,2)
      h(1) = 1._WP / res(3,1)
      h(2) = 384400._WP / res(3,2)
      res2(3,1)=h(1)*h(1)*h(1)
      res2(3,2)=h(2)*h(2)*h(2)
!
      Return
      END SUBROUTINE aufb2

      SUBROUTINE negangle2(x)
!  transformation of negative angles
!  to angles in interval [0;2pi)
!
      IMPLICIT NONE
!
      Logical endvar
      REAL(WP) :: x
!
      endvar=.False.
!
      Do while (.Not.endvar)
        x=x+pi*2
        endvar = (x .GE. 0._WP)
      Enddo
!
      Return
      END  SUBROUTINE  negangle2

      SUBROUTINE langle2(x)
!  transformation of large angles
!  to angles in interval [0;2pi)
!
      IMPLICIT NONE
!
      Logical endvar
      REAL(WP) :: x
!
      endvar=.False.
!
      Do while (.Not.endvar)
        x=x-pi*2
        endvar=(x.Lt.(pi*2))
      Enddo
!
      Return
      END SUBROUTINE    langle2

   END MODULE mo_tidal


