module recom_diag

!>
!! @brief Module computing the total carbon and silicate budget
!!
!! $ID: n/a$
!!
!
  use g_config
  use mod_mesh
  use g_parsup
  use g_clock
  use g_comm_auto
  use o_ARRAYS
  use g_forcing_arrays
  use i_ARRAYS
  use o_mixing_KPP_mod
  use g_rotate_grid
  use g_support
  use REcoM_GloVar
  use io_mesh_info
  use recom_config

  implicit none
#include "netcdf.inc"
  private

  public :: ldiag_carbon, ldiag_silicate, recom_diag_freq, recom_diag_freq_unit, recom_logfile_outfreq, total_carbon, total_silicate, &
            compute_carbon_diag, compute_silicate_diag, write_recom_diag, compute_recom_diagnostics, precom_diag_list

  real(kind=WP),  save,  target                 :: total_carbon
  real(kind=WP),  save,  target                 :: valDIC
  real(kind=WP),  save,  target                 :: valDOC
  real(kind=WP),  save,  target                 :: valPhyC
  real(kind=WP),  save,  target                 :: valDetC
  real(kind=WP),  save,  target                 :: valHetC
  real(kind=WP),  save,  target                 :: valDiaC
  real(kind=WP),  save,  target                 :: valPhyCalc
  real(kind=WP),  save,  target                 :: valDetCalc
  real(kind=WP),  save,  target                 :: valDSi
  real(kind=WP),  save,  target                 :: valDiaSi
  real(kind=WP),  save,  target                 :: valDetSi
  real(kind=WP),  save,  target                 :: valDetz2Si
  real(kind=WP),  save,  target                 :: valBenSi

  real(kind=WP),  save,  target                 :: total_silicate

  logical                                       :: ldiag_carbon        =.true.
  logical                                       :: ldiag_silicate      =.true.
  integer                                       :: recom_diag_freq       = 1         !only required for d,h,s cases,  y, m take 1
  character                                     :: recom_diag_freq_unit  = 'm'        !output period: y,  d, h, s 
  integer                                       :: recom_logfile_outfreq = 120         !in logfile info. output frequency, # steps

  real(kind=WP)                                 :: ctime !current time in seconds from the beginning of the year
  integer                                       :: row 
  
  namelist /precom_diag_list/ ldiag_carbon, ldiag_silicate, recom_diag_freq, recom_diag_freq_unit, recom_logfile_outfreq

  contains

  !---------------------------------------------------------------------------
  !>
  !! @brief Computes the interactive total carbon and total silicate
  !! 
  !! @remarks This routine reads mesh 
  !


! ==============================================================
subroutine compute_carbon_diag(mode,mesh)

  implicit none
  integer, intent(in)        :: mode
  type(t_mesh), intent(in)  , target :: mesh
  logical, save              :: firstcall=.true.

  if (firstcall) then  !allocate the stuff at the first call
    total_carbon=0.0
    firstcall=.false.
    if (mode==0) return
  end if

        ! DIC
        call integrate_nod(tr_arr(:,:,4), valDIC, mesh)
        total_carbon=total_carbon+valDIC
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DIC at timestep :', mstep, valDIC
        end if

        ! DOC
        call integrate_nod(tr_arr(:,:,14), valDOC, mesh)
        total_carbon=total_carbon+valDOC
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DOC at timestep :', mstep, valDOC
        end if

        !PhyC
        call integrate_nod(tr_arr(:,:,7), valPhyC, mesh)
        total_carbon=total_carbon+valPhyC
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of PhyC at timestep :', mstep, valPhyC
        end if

        !DetC
        call integrate_nod(tr_arr(:,:,10), valDetC, mesh)
        total_carbon=total_carbon+valDetC
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DetC at timestep :', mstep, valDetC
        end if

        !HetC
        call integrate_nod(tr_arr(:,:,12), valHetC, mesh)
        total_carbon=total_carbon+valHetC
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of HetC at timestep :', mstep, valHetC
        end if

        !DiaC
        call integrate_nod(tr_arr(:,:,16), valDiaC, mesh)
        total_carbon=total_carbon+valDiaC
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DiaC at timestep :', mstep, valDiaC
        end if

       !PhyCalc
        call integrate_nod(tr_arr(:,:,22), valPhyCalc, mesh)
        total_carbon=total_carbon+valPhyCalc
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of PhyCalc at timestep :', mstep, valPhyCalc
        end if

        !DetCalc
        call integrate_nod(tr_arr(:,:,23), valDetCalc, mesh)
        total_carbon=total_carbon+valDetCalc
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DetCalc at timestep :', mstep, valDetCalc
        end if

        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of carbon at timestep :', mstep, total_carbon
        end if

end subroutine compute_carbon_diag


subroutine compute_silicate_diag(mode,mesh)

  implicit none
  integer, intent(in)        :: mode
  type(t_mesh), intent(in)  , target :: mesh
  logical, save              :: firstcall=.true.

  if (firstcall) then  !allocate the stuff at the first call
    firstcall=.false.
    if (mode==0) return
  end if
    total_silicate=0.0 ! snapshot

        !DSi
        call integrate_nod(tr_arr(:,:,20), valDSi, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral of DSi at timestep :', mstep, valDSi
        total_silicate=total_silicate+valDSi

        !DiaSi
        call integrate_nod(tr_arr(:,:,18), valDiaSi, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral of DiaSi at timestep :', mstep, valDiaSi
        total_silicate=total_silicate+valDiaSi

        !DetSi
        call integrate_nod(tr_arr(:,:,19), valDetSi, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral of DetSi at timestep :', mstep, valDetSi
        total_silicate=total_silicate+valDetSi

!if (REcoM_Second_Zoo) then
        !Detz2Si
        call integrate_nod(tr_arr(:,:,29), valDetz2Si, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral of Detz2Si at timestep :', mstep, valDetSi
        total_silicate=total_silicate+valDetz2Si
!end if 
        !BenSi
!        call integrate_nod(Benthos(:,3), valBenSi, mesh)
         call integrate_bottom(valBenSi,mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral of BenSi at timestep :', mstep, valBenSi
        total_silicate=total_silicate+valBenSi

        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral of silicate at timestep :', mstep, total_silicate

end subroutine compute_silicate_diag

! ==============================================================
subroutine write_recom_diag(mode, mesh)

  implicit none
  integer                            :: status, ncid, j, k
  character(2000)                    :: filename
  integer                            :: recID, tID, tcID, tsID, tsdelID
  integer                            :: valDICID, valDOCID, valPhyCID, valDetCID, valHetCID, valDiaCID, valPhyCalcID, valDetCalcID
  integer                            :: valDSiID, valDiaSiID, valDetSiID, valDetz2SiID, valBenSiID
  integer                            :: rec_count=0
  character(2000)                    :: att_text
  real(real64)                       :: rtime !timestamp of the record
  logical                            :: do_output
  integer, intent(in)                :: mode
  logical, save                      :: firstcall=.true.
  type(t_mesh), intent(in)  , target :: mesh


! only master (rank=0) writes the output out
  if (mype/=0) return

  ctime=timeold+(dayold-1.)*86400

! control for output frequency 
  do_output=.false.

  if (recom_diag_freq_unit.eq.'y') then
     call annual_event(do_output)

  else if (recom_diag_freq_unit == 'm') then 
     call monthly_event(do_output) 

  else if (recom_diag_freq_unit == 'd') then
     call daily_event(do_output, recom_diag_freq)

  else if (recom_diag_freq_unit == 'h') then
     call hourly_event(do_output, recom_diag_freq)

  else if (recom_diag_freq_unit == 's') then
     call step_event(do_output, mstep, recom_diag_freq)

  else
     write(*,*) 'You did not specify a supported outputflag.'
     write(*,*) 'The program will stop to give you opportunity to do it.'
     call par_ex(1)
     stop
  endif

! create output file

  if (firstcall) then  !create the stuff at the first call

     filename=trim(ResultPath)//trim(runid)//'.'//cyearnew//'.recom.diag.nc'
     status = nf_create(filename, IOR(NF_CLOBBER,IOR(NF_NETCDF4,NF_CLASSIC_MODEL)), ncid)
!! define dimensions (time unlimited in this case)
     status = nf_def_dim(ncid, 'time', NF_UNLIMITED, recID)
!! define variables
     status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, recID, tID)
     status = nf_def_var(ncid, 'total_carbon', NF_DOUBLE, 1, recID, tcID)
     status = nf_def_var(ncid, 'total_DIC', NF_DOUBLE, 1, recID, valDICID)
     status = nf_def_var(ncid, 'total_DOC', NF_DOUBLE, 1, recID, valDOCID)
     status = nf_def_var(ncid, 'total_PhyC', NF_DOUBLE, 1, recID, valPhyCID)
     status = nf_def_var(ncid, 'total_DetC', NF_DOUBLE, 1, recID, valDetCID)
     status = nf_def_var(ncid, 'total_HetC', NF_DOUBLE, 1, recID, valHetCID)
     status = nf_def_var(ncid, 'total_DiaC', NF_DOUBLE, 1, recID, valDiaCID)
     status = nf_def_var(ncid, 'total_PhyCalc', NF_DOUBLE, 1, recID, valPhyCalcID)
     status = nf_def_var(ncid, 'total_DetCalc', NF_DOUBLE, 1, recID, valDetCalcID)

     status = nf_def_var(ncid, 'total_silicate', NF_DOUBLE, 1, recID, tsID)

     status = nf_def_var(ncid, 'total_DSi', NF_DOUBLE, 1, recID, valDSiID)
     status = nf_def_var(ncid, 'total_DiaSi', NF_DOUBLE, 1, recID, valDiaSiID)
     status = nf_def_var(ncid, 'total_DetSi', NF_DOUBLE, 1, recID, valDetSiID)
     status = nf_def_var(ncid, 'total_Detz2Si', NF_DOUBLE, 1, recID, valDetz2SiID)
     status = nf_def_var(ncid, 'total_BenSi', NF_DOUBLE, 1, recID, valBenSiID)


!! add attributes
     att_text='time'
     status = nf_put_att_text(ncid, tID, 'long_name', len_trim(att_text), trim(att_text))
     write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
     status = nf_put_att_text(ncid, tID, 'units', len_trim(att_text), trim(att_text))
     status = nf_close(ncid)

     firstcall=.false.
     if (mode==0) return
  end if
  

if (do_output) then

! fill the output file
  filename=trim(ResultPath)//trim(runid)//'.'//cyearnew//'.recom.diag.nc'
  status = nf_open(filename, nf_write, ncid)

  status = nf_inq_dimid (ncid, 'time', recID)
  status = nf_inq_dimlen(ncid, recID, rec_count)

  status = nf_inq_varid(ncid, 'time', tID)
  status = nf_inq_varid(ncid, 'total_carbon', tcID)

  status = nf_inq_varid(ncid, 'total_DIC', valDICID)
  status = nf_inq_varid(ncid, 'total_DOC', valDOCID)
  status = nf_inq_varid(ncid, 'total_PhyC', valPhyCID)
  status = nf_inq_varid(ncid, 'total_DetC', valDetCID)
  status = nf_inq_varid(ncid, 'total_HetC', valHetCID)
  status = nf_inq_varid(ncid, 'total_DiaC', valDiaCID)
  status = nf_inq_varid(ncid, 'total_PhyCalc', valPhyCalcID)
  status = nf_inq_varid(ncid, 'total_DetCalc', valDetCalcID)

  status = nf_inq_varid(ncid, 'total_silicate', tsID)

  status = nf_inq_varid(ncid, 'total_DSi', valDSiID)
  status = nf_inq_varid(ncid, 'total_DiaSi', valDiaSiID)
  status = nf_inq_varid(ncid, 'total_DetSi', valDetSiID)
  status = nf_inq_varid(ncid, 'total_Detz2Si', valDetz2SiID)
  status = nf_inq_varid(ncid, 'total_BenSi', valBenSiID)

  do k=rec_count, 1, -1
     status=nf_get_vara_double(ncid, tID, k, 1, rtime, 1);
     if (ctime > rtime) then
        rec_count=k+1
!       write(*,*) 'I/O '//trim(entry%name)//' : current record = ', entry%rec_count, '; ', entry%rec_count, ' records in the file;'
        exit ! a proper rec_count detected, exit the loop
     end if
     if (k==1) then
        write(*,*) 'I/O '//'time'//' WARNING: the existing output file will be overwritten'//'; ', rec_count, ' records in the file;'
        rec_count=1
        exit ! no appropriate rec_count detected
     end if
  end do

  rec_count=max(rec_count, 1)

  status = nf_put_vara_double(ncid, tID, rec_count, 1, ctime, 1)
  status = nf_put_vara_double(ncid, tcID, rec_count, 1, total_carbon, 1)
  status = nf_put_vara_double(ncid, valDICID, rec_count, 1, valDIC, 1)
  status = nf_put_vara_double(ncid, valDOCID, rec_count, 1, valDOC, 1)
  status = nf_put_vara_double(ncid, valPhyCID, rec_count, 1, valPhyC, 1)
  status = nf_put_vara_double(ncid, valDetCID, rec_count, 1, valDetC, 1)
  status = nf_put_vara_double(ncid, valHetCID, rec_count, 1, valHetC, 1)
  status = nf_put_vara_double(ncid, valDiaCID, rec_count, 1, valDiaC, 1)
  status = nf_put_vara_double(ncid, valPhyCalcID, rec_count, 1, valPhyCalc, 1)
  status = nf_put_vara_double(ncid, valDetCalcID, rec_count, 1, valDetCalc, 1)

  status = nf_put_vara_double(ncid, tsID, rec_count, 1, total_silicate, 1)

  status = nf_put_vara_double(ncid, valDSiID, rec_count, 1, valDSi, 1)
  status = nf_put_vara_double(ncid, valDiaSiID, rec_count, 1, valDiaSi, 1)
  status = nf_put_vara_double(ncid, valDetSiID, rec_count, 1, valDetSi, 1)
  status = nf_put_vara_double(ncid, valDetz2SiID, rec_count, 1, valDetz2Si, 1)
  status = nf_put_vara_double(ncid, valBenSiID, rec_count, 1, valBenSi, 1)

  status=nf_close(ncid)

end if

end subroutine write_recom_diag

! ==============================================================
subroutine compute_recom_diagnostics(mode, mesh)

  implicit none
  integer, intent(in)                :: mode !constructor mode (0=only allocation; any other=do diagnostic)
  type(t_mesh), intent(in)  , target :: mesh

  !1. carbon diagnostic
  if (ldiag_carbon)      call compute_carbon_diag(mode,mesh)
  !2. silicate diagnostic
  if (ldiag_silicate)    call compute_silicate_diag(mode,mesh)
  !3. write total carbon and silicate out into recom.diag.nc
  if (ldiag_carbon .or. ldiag_silicate) call write_recom_diag(mode, mesh)

end subroutine compute_recom_diagnostics


end module recom_diag
