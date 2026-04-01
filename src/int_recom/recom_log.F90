!=============================================================================
! REcoM LOGGING MODULE  —  centralised, colourful log statements
!
! All print / write statements that were scattered inside
!   #if defined (__recom) ... #endif
! blocks of gen_surface_forcing.F90 are collected here.
!
! ANSI colour legend used throughout:
!   CYAN    [36m  — debug / verbose trace (recom_debug only)
!   GREEN   [32m  — success / file opened / data loaded OK
!   YELLOW  [33m  — informational update (monthly data, year change)
!   MAGENTA [35m  — CO2 / carbon cycle values
!   RED     [31m  — errors (abort follows immediately)
!   RESET   [0m   — always appended after every coloured line
!=============================================================================

#if defined (__recom)

   !==========================================================================
   !  Helper: emit one coloured line from rank-0 only
   !  Usage:  call recom_log(mype, COLOR_GREEN, 'message text')
   !
   !  Colour constants:
   !    COLOR_CYAN    = achar(27)//'[36m'
   !    COLOR_GREEN   = achar(27)//'[32m'
   !    COLOR_YELLOW  = achar(27)//'[33m'
   !    COLOR_MAGENTA = achar(27)//'[35m'
   !    COLOR_RED     = achar(27)//'[31m'
   !    COLOR_RESET   = achar(27)//'[0m'
   !==========================================================================
   SUBROUTINE recom_log(mype, color_code, msg)
      implicit none
      integer,          intent(in) :: mype
      character(len=*), intent(in) :: color_code
      character(len=*), intent(in) :: msg
      if (mype == 0) then
         write(*, '(3A)') trim(color_code), trim(msg), achar(27)//'[0m'
      end if
   END SUBROUTINE recom_log

   !==========================================================================
   !  SECTION 1 — sbc_ini  (called once at initialisation)
   !  Original location: end of sbc_ini, inside  #if defined(__recom)
   !==========================================================================
   SUBROUTINE log_recom_namelist_open(mype, iostat_val)
      ! Replaces the open/error block for namelist.recom
      implicit none
      integer, intent(in) :: mype, iostat_val
      if (iostat_val == 0) then
         call recom_log(mype, achar(27)//'[32m', &
              '     file   : namelist.recom for sbc  open ok')
      else
         call recom_log(mype, achar(27)//'[31m', &
              'ERROR: --> bad opening file   : namelist.recom for sbc')
         ! caller must invoke par_ex / stop after this
      end if
   END SUBROUTINE log_recom_namelist_open

   !==========================================================================
   !  SECTION 2 — sbc_do : debug entry trace
   !  Original: if (recom_debug .and. mype==0) print *, '... --> Atm_input ...'
   !==========================================================================
   SUBROUTINE log_recom_atm_input_entry(mype, recom_debug)
      implicit none
      integer, intent(in) :: mype
      logical, intent(in) :: recom_debug
      if (recom_debug .and. mype == 0) then
         call recom_log(mype, achar(27)//'[36m', &
              '  --> Atm_input  (entering atmospheric forcing update)')
      end if
   END SUBROUTINE log_recom_atm_input_entry

   !==========================================================================
   !  SECTION 3 — Atmospheric CO2
   !  Covers: constant_CO2, transient file, atbox, control output
   !==========================================================================

   ! 3a  Constant CO2 spinup value
   SUBROUTINE log_co2_constant(mype, CO2_for_spinup)
      implicit none
      integer,  intent(in) :: mype
      real(kind=8), intent(in) :: CO2_for_spinup
      if (mype == 0) then
         call recom_log(mype, achar(27)//'[35m', &
              '  [CO2] Mode: constant spinup value')
         write(*, '(A,F10.4)') achar(27)//'[35m'// &
              '         Constant_CO2 = ', CO2_for_spinup
         write(*, '(A)') achar(27)//'[0m'
      end if
   END SUBROUTINE log_co2_constant

   ! 3b  Transient CO2 file read
   SUBROUTINE log_co2_transient_update(mype, month_idx, currentCO2year, filename)
      implicit none
      integer,          intent(in) :: mype, month_idx, currentCO2year
      character(len=*), intent(in) :: filename
      if (mype == 0) then
         write(*, '(2A,I0,A,I0,2A)') achar(27)//'[35m', &
              '  [CO2] Updating climatology for month ', month_idx, &
              '  |  carbon year = ', currentCO2year, &
              '  from ', trim(filename)
         write(*, '(A)') achar(27)//'[0m'
      end if
   END SUBROUTINE log_co2_transient_update

   ! 3c  CO2 error: cannot open file
   SUBROUTINE log_co2_file_error(mype, filename)
      implicit none
      integer,          intent(in) :: mype
      character(len=*), intent(in) :: filename
      if (mype == 0) then
         call recom_log(mype, achar(27)//'[31m', &
              'ERROR: CANNOT READ CO2 FILE CORRECTLY !!!!!')
         call recom_log(mype, achar(27)//'[31m', &
              'Error in opening netcdf file: '//trim(filename))
      end if
   END SUBROUTINE log_co2_file_error

   ! 3d  Control output of atmospheric CO2 values (AtmCO2, AtmCO2_13, AtmCO2_14)
   SUBROUTINE log_co2_control_output(mype, AtmCO2_val, &
                                     ciso, AtmCO2_13_val, &
                                     ciso_14, AtmCO2_14_val, &
                                     use_atbox)
      implicit none
      integer,  intent(in) :: mype
      real(kind=8), intent(in) :: AtmCO2_val, AtmCO2_13_val, AtmCO2_14_val
      logical,  intent(in) :: ciso, ciso_14, use_atbox
      if (mype == 0) then
         write(*, '(2A,ES12.5)') achar(27)//'[35m', &
              '  [CO2] AtmCO2     = ', AtmCO2_val
         if (ciso) then
            write(*, '(A,ES12.5)') &
                 '        AtmCO2_13 = ', AtmCO2_13_val
            if (ciso_14) &
               write(*, '(A,ES12.5)') &
                    '        AtmCO2_14 = ', AtmCO2_14_val
         end if
         if (use_atbox) &
            write(*, '(A)') '        use_atbox = .true.'
         write(*, '(A)') achar(27)//'[0m'
      end if
   END SUBROUTINE log_co2_control_output

   !==========================================================================
   !  SECTION 4 — Iron (Fe) deposition
   !==========================================================================

   ! 4a  Monthly Albani Fe update
   SUBROUTINE log_fe_update(mype, month_idx, filename)
      implicit none
      integer,          intent(in) :: mype, month_idx
      character(len=*), intent(in) :: filename
      if (mype == 0) then
         write(*, '(2A,I3,2A)') achar(27)//'[33m', &
              '  [Fe]  Updating iron climatology for month ', month_idx, &
              '  from ', trim(filename)
         write(*, '(A)') achar(27)//'[0m'
      end if
   END SUBROUTINE log_fe_update

   ! 4b  Albani switched off
   SUBROUTINE log_fe_disabled(mype)
      implicit none
      integer, intent(in) :: mype
      call recom_log(mype, achar(27)//'[33m', &
           '  [Fe]  Albani dust is switched off --> Check namelist.recom')
   END SUBROUTINE log_fe_disabled

   !==========================================================================
   !  SECTION 5 — Nitrogen deposition
   !==========================================================================

   ! 5a  Monthly N update
   SUBROUTINE log_nitrogen_update(mype, month_idx, filename, Nvari)
      implicit none
      integer,          intent(in) :: mype, month_idx
      character(len=*), intent(in) :: filename, Nvari
      if (mype == 0) then
         write(*, '(2A,I3,4A)') achar(27)//'[33m', &
              '  [N]   Updating nitrogen climatology for month ', month_idx, &
              '  variable=', trim(Nvari), &
              '  from ', trim(filename)
         write(*, '(A)') achar(27)//'[0m'
      end if
   END SUBROUTINE log_nitrogen_update

   ! 5b  N deposition disabled
   SUBROUTINE log_nitrogen_disabled(mype)
      implicit none
      integer, intent(in) :: mype
      call recom_log(mype, achar(27)//'[33m', &
           '  [N]   useAeolianN = .false.')
   END SUBROUTINE log_nitrogen_disabled

   !==========================================================================
   !  SECTION 6 — River inputs
   !==========================================================================

   ! 6a  R2OMIP: constant pre-industrial
   SUBROUTINE log_river_r2omip_pi(mype, filename)
      implicit none
      integer,          intent(in) :: mype
      character(len=*), intent(in) :: filename
      if (mype == 0) then
         call recom_log(mype, achar(27)//'[32m', &
              '  [RIV] Mode: constant pre-industrial riverine inputs (R2OMIP PI)')
         call recom_log(mype, achar(27)//'[32m', &
              '        Opening: '//trim(filename))
      end if
   END SUBROUTINE log_river_r2omip_pi

   ! 6b  R2OMIP: transient year-specific
   SUBROUTINE log_river_r2omip_transient(mype, cyearnew, filename)
      implicit none
      integer,          intent(in) :: mype
      character(len=*), intent(in) :: cyearnew, filename
      if (mype == 0) then
         write(*, '(4A)') achar(27)//'[32m', &
              '  [RIV] Mode: transient R2OMIP  year=', trim(cyearnew), &
              '  from '//trim(filename)
         write(*, '(A)') achar(27)//'[0m'
      end if
   END SUBROUTINE log_river_r2omip_transient

   ! 6c  R2OMIP file open error
   SUBROUTINE log_river_r2omip_error(mype, filename)
      implicit none
      integer,          intent(in) :: mype
      character(len=*), intent(in) :: filename
      if (mype == 0) then
         call recom_log(mype, achar(27)//'[31m', &
              'ERROR: Failed to open river input file: '//trim(filename))
      end if
   END SUBROUTINE log_river_r2omip_error

   ! 6d  R2OMIP sanity check start
   SUBROUTINE log_river_sanity_start(mype)
      implicit none
      integer, intent(in) :: mype
      call recom_log(mype, achar(27)//'[32m', &
           '  [RIV] Sanity-checking R2OMIP river variables...')
   END SUBROUTINE log_river_sanity_start

   ! 6e  Standard monthly river update
   SUBROUTINE log_river_monthly_update(mype, month_idx, filename)
      implicit none
      integer,          intent(in) :: mype, month_idx
      character(len=*), intent(in) :: filename
      if (mype == 0) then
         write(*, '(2A,I3,2A)') achar(27)//'[33m', &
              '  [RIV] Updating riverine data for month ', month_idx, &
              '  from ', trim(filename)
         write(*, '(A)') achar(27)//'[0m'
      end if
   END SUBROUTINE log_river_monthly_update

   ! 6f  Rivers disabled
   SUBROUTINE log_river_disabled(mype)
      implicit none
      integer, intent(in) :: mype
      call recom_log(mype, achar(27)//'[33m', &
           '  [RIV] Riverine input disabled  (useRivers = .false.)')
   END SUBROUTINE log_river_disabled

   !==========================================================================
   !  SECTION 7 — Erosion inputs
   !==========================================================================

   ! 7a  Debug entry trace
   SUBROUTINE log_erosion_entry(mype, recom_debug)
      implicit none
      integer, intent(in) :: mype
      logical, intent(in) :: recom_debug
      if (recom_debug .and. mype == 0) then
         call recom_log(mype, achar(27)//'[36m', &
              '  --> Erosion_input  (entering erosion update)')
      end if
   END SUBROUTINE log_erosion_entry

   ! 7b  Monthly erosion update
   SUBROUTINE log_erosion_monthly_update(mype, month_idx, filename)
      implicit none
      integer,          intent(in) :: mype, month_idx
      character(len=*), intent(in) :: filename
      if (mype == 0) then
         write(*, '(2A,I3,2A)') achar(27)//'[33m', &
              '  [ERO] Updating erosion data for month ', month_idx, &
              '  from ', trim(filename)
         write(*, '(A)') achar(27)//'[0m'
      end if
   END SUBROUTINE log_erosion_monthly_update

   ! 7c  Erosion disabled
   SUBROUTINE log_erosion_disabled(mype)
      implicit none
      integer, intent(in) :: mype
      call recom_log(mype, achar(27)//'[33m', &
           '  [ERO] Erosion input disabled  (useErosion = .false.)')
   END SUBROUTINE log_erosion_disabled

#endif
!=============================================================================
!
!  HOW TO REPLACE CALLS IN gen_surface_forcing.F90
!  ─────────────────────────────────────────────────
!
!  sbc_ini
!  ───────
!  BEFORE:
!      if (iost == 0) then
!          if (mype==0) WRITE(*,*) '     file   : ', 'namelist.recom for sbc',' open ok'
!      else
!          if (mype==0) WRITE(*,*) 'ERROR: --> bad opening file   : namelist.recom ...'
!
!  AFTER:
!      call log_recom_namelist_open(mype, iost)
!      if (iost /= 0) then; call par_ex(...); stop; end if
!
!  ─────────────────────────────────────────────────
!  sbc_do  —  debug entry
!  ───────
!  BEFORE:
!      if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> Atm_input'//...
!
!  AFTER:
!      call log_recom_atm_input_entry(mype, recom_debug)
!
!  ─────────────────────────────────────────────────
!  sbc_do  —  CO2
!  ──────────────
!  BEFORE:
!      if (mype==0) write(*,*) 'Updating CO2 climatology for month  ', i,' from ', trim(filename)
!      ...
!      if (mype==0) write(*,*),'Current carbon year=',currentCO2year
!      if (mype==0) write(*,*),'Atm CO2=', AtmCO2
!      print *, "In atm_input: AtmCO2    = ", AtmCO2(1)
!      ...
!
!  AFTER:
!      call log_co2_transient_update(mype, i, currentCO2year, filename)
!      ...
!      call log_co2_control_output(mype, AtmCO2(1), ciso, AtmCO2_13(1), &
!                                  ciso_14, AtmCO2_14(1,1), use_atbox)
!
!  ─────────────────────────────────────────────────
!  sbc_do  —  Fe
!  ─────────────
!  BEFORE:
!      if (mype==0) write(*,*) 'Updating iron climatology for month  ', i,' from ', trim(filename)
!
!  AFTER:
!      call log_fe_update(mype, i, filename)
!
!  ─────────────────────────────────────────────────
!  sbc_do  —  N
!  ────────────
!  BEFORE:
!      if (mype==0) write(*,*) 'Updating nitrogen climatology for month  ', i,' from ', trim(filename)
!      ...
!      GloNDust = 0.0_WP
!      if (mstep==1 .and. mype==0) write(*,*) 'useAeolianN is switched off'
!
!  AFTER:
!      call log_nitrogen_update(mype, i, filename, Nvari)
!      ...
!      call log_nitrogen_disabled(mype)
!
!  ─────────────────────────────────────────────────
!  sbc_do  —  Rivers
!  ─────────────────
!  BEFORE (R2OMIP PI):
!      write(*,'(A)') 'INFO: Using constant pre-industrial riverine inputs'
!
!  AFTER:
!      call log_river_r2omip_pi(mype, filename)
!
!  BEFORE (R2OMIP transient):
!      write(*,'(4A)') 'INFO: Using transient riverine inputs ...', trim(cyearnew), ...
!
!  AFTER:
!      call log_river_r2omip_transient(mype, cyearnew, filename)
!
!  BEFORE (sanity check):
!      if (mype == 0) write(*,*) 'Sanity-checking R2OMIP river variables...'
!
!  AFTER:
!      call log_river_sanity_start(mype)
!
!  BEFORE (standard monthly):
!      write(*,'(A,I2,2A)') 'Updating riverine data for month ', i, ' from ', trim(filename)
!
!  AFTER:
!      call log_river_monthly_update(mype, i, filename)
!
!  BEFORE (disabled):
!      if (mype == 0 .and. mstep == 1) write(*,*) 'INFO: Riverine input disabled ...'
!
!  AFTER:
!      if (mstep == 1) call log_river_disabled(mype)
!
!  ─────────────────────────────────────────────────
!  sbc_do  —  Erosion
!  ──────────────────
!  BEFORE:
!      if (recom_debug .and. mype == 0) print *, achar(27)//'[36m'//' --> Erosion_input'//...
!
!  AFTER:
!      call log_erosion_entry(mype, recom_debug)
!
!  BEFORE:
!      write(*,'(A,I2,2A)') 'Updating erosion data for month ', i, ' from ', trim(filename)
!
!  AFTER:
!      call log_erosion_monthly_update(mype, i, filename)
!
!  BEFORE:
!      if (mype == 0 .and. mstep == 1) write(*,*) 'INFO: Erosion input disabled ...'
!
!  AFTER:
!      if (mstep == 1) call log_erosion_disabled(mype)
!
!=============================================================================

