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

  public :: ldiag_carbon, ldiag_silicate, ldiag_nitrate, recom_diag_freq, recom_diag_freq_unit, recom_logfile_outfreq, total_carbon_66N, total_silicate_66N, total_nitrate_66N, total_carbon_80N, total_silicate_80N, total_nitrate_80N, compute_carbon_diag, compute_silicate_diag, compute_nitrate_diag, write_recom_diag, compute_recom_diagnostics, precom_diag_list

  real(kind=WP),  save,  target                 :: total_carbon_66N
  real(kind=WP),  save,  target                 :: total_silicate_66N
  real(kind=WP),  save,  target                 :: total_nitrate_66N
  real(kind=WP),  save,  target                 :: valDIC_66N
  real(kind=WP),  save,  target                 :: valDOC_66N
  real(kind=WP),  save,  target                 :: valPhyC_66N
  real(kind=WP),  save,  target                 :: valDetC_66N
  real(kind=WP),  save,  target                 :: valHetC_66N
  real(kind=WP),  save,  target                 :: valDiaC_66N
  real(kind=WP),  save,  target                 :: valBenC_66N
  real(kind=WP),  save,  target                 :: valPhyCalc_66N
  real(kind=WP),  save,  target                 :: valDetCalc_66N
  real(kind=WP),  save,  target                 :: valBenCalc_66N
  real(kind=WP),  save,  target                 :: valDSi_66N
  real(kind=WP),  save,  target                 :: valDiaSi_66N
  real(kind=WP),  save,  target                 :: valDetSi_66N
  real(kind=WP),  save,  target                 :: valBenSi_66N
  real(kind=WP),  save,  target                 :: valDIN_66N
  real(kind=WP),  save,  target                 :: valDON_66N
  real(kind=WP),  save,  target                 :: valDiaN_66N
  real(kind=WP),  save,  target                 :: valPhyN_66N
  real(kind=WP),  save,  target                 :: valHetN_66N
  real(kind=WP),  save,  target                 :: valDetN_66N
  real(kind=WP),  save,  target                 :: valBenN_66N
  
  real(kind=WP),  save,  target                 :: total_silicate_80N
  real(kind=WP),  save,  target                 :: total_nitrate_80N
  real(kind=WP),  save,  target                 :: total_carbon_80N
  real(kind=WP),  save,  target                 :: valDIC_80N
  real(kind=WP),  save,  target                 :: valDOC_80N
  real(kind=WP),  save,  target                 :: valPhyC_80N
  real(kind=WP),  save,  target                 :: valDetC_80N
  real(kind=WP),  save,  target                 :: valHetC_80N
  real(kind=WP),  save,  target                 :: valDiaC_80N
  real(kind=WP),  save,  target                 :: valBenC_80N
  real(kind=WP),  save,  target                 :: valPhyCalc_80N
  real(kind=WP),  save,  target                 :: valDetCalc_80N
  real(kind=WP),  save,  target                 :: valBenCalc_80N
  real(kind=WP),  save,  target                 :: valDSi_80N
  real(kind=WP),  save,  target                 :: valDiaSi_80N
  real(kind=WP),  save,  target                 :: valDetSi_80N
  real(kind=WP),  save,  target                 :: valBenSi_80N
  real(kind=WP),  save,  target                 :: valDIN_80N
  real(kind=WP),  save,  target                 :: valDON_80N
  real(kind=WP),  save,  target                 :: valDiaN_80N
  real(kind=WP),  save,  target                 :: valPhyN_80N
  real(kind=WP),  save,  target                 :: valHetN_80N
  real(kind=WP),  save,  target                 :: valDetN_80N
  real(kind=WP),  save,  target                 :: valBenN_80N

  logical                                       :: ldiag_carbon        =.true.
  logical                                       :: ldiag_silicate      =.true.
  logical                                       :: ldiag_nitrate       =.true.
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
! ==============================================================
subroutine compute_carbon_diag(mode,mesh)

  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use g_clock
  use o_PARAM
  USE o_ARRAYS
  use g_PARSUP
  use mod_MESH
  use g_comm_auto
  implicit none
  integer, intent(in)        :: mode
  type(t_mesh), intent(in)  , target :: mesh
  logical, save              :: firstcall=.true.
  integer                    :: k
#include  "../associate_mesh.h"
  
  ind_arctic_66_3D = 1.0_WP
  ind_arctic_66_2D = 1.0_WP
  do k=1, myDim_nod2D
      if (geo_coord_nod2D(2,k) < 66*rad) then
          ind_arctic_66_3D(:,k) = 0.0_WP
          ind_arctic_66_2D(k) = 0.0_WP
      end if
  end do
  
  ind_arctic_80_3D = 1.0_WP
  ind_arctic_80_2D = 1.0_WP
  do k=1, myDim_nod2D
      if (geo_coord_nod2D(2,k) < 80*rad) then
          ind_arctic_80_3D(:,k) = 0.0_WP
          ind_arctic_80_2D(k) = 0.0_WP
      end if
  end do

  if (firstcall) then  !allocate the stuff at the first call
    total_carbon_80N=0.0
    total_carbon_66N=0.0
    firstcall=.false.
    if (mode==0) return
  end if

        ! DIC
        call integrate_nod(tr_arr(:,:,4)*ind_arctic_80_3D, valDIC_80N, mesh)
        total_carbon_80N=total_carbon_80N+valDIC_80N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DIC (>80N) at timestep :', mstep, valDIC_80N
        end if

        ! DOC
        call integrate_nod(tr_arr(:,:,14)*ind_arctic_80_3D, valDOC_80N, mesh)
        total_carbon_80N=total_carbon_80N+valDOC_80N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DOC (>80N) at timestep :', mstep, valDOC_80N
        end if

        !PhyC
        call integrate_nod(tr_arr(:,:,7)*ind_arctic_80_3D, valPhyC_80N, mesh)
        total_carbon_80N=total_carbon_80N+valPhyC_80N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of PhyC (>80N) at timestep :', mstep, valPhyC_80N
        end if

        !DetC
        call integrate_nod(tr_arr(:,:,10)*ind_arctic_80_3D+tr_arr(:,:,28)*ind_arctic_80_3D, valDetC_80N, mesh)
        total_carbon_80N=total_carbon_80N+valDetC_80N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DetC+DetZ2C (>80N) at timestep :', mstep, valDetC_80N
        end if

        !HetC
        call integrate_nod(tr_arr(:,:,12)*ind_arctic_80_3D+tr_arr(:,:,26)*ind_arctic_80_3D, valHetC_80N, mesh)
        total_carbon_80N=total_carbon_80N+valHetC_80N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of HetC+Zoo2C (>80N) at timestep :', mstep, valHetC_80N
        end if

        !DiaC
        call integrate_nod(tr_arr(:,:,16)*ind_arctic_80_3D, valDiaC_80N, mesh)
        total_carbon_80N=total_carbon_80N+valDiaC_80N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DiaC (>80N) at timestep :', mstep, valDiaC_80N
        end if

       !PhyCalc
        call integrate_nod(tr_arr(:,:,22)*ind_arctic_80_3D, valPhyCalc_80N, mesh)
        total_carbon_80N=total_carbon_80N+valPhyCalc_80N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of PhyCalc (>80N) at timestep :', mstep, valPhyCalc_80N
        end if

        !DetCalc
        call integrate_nod(tr_arr(:,:,23)*ind_arctic_80_3D+tr_arr(:,:,30)*ind_arctic_80_3D, valDetCalc_80N, mesh)
        total_carbon_80N=total_carbon_80N+valDetCalc_80N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DetCalc+DetZ2Calc (>80N) at timestep :', mstep, valDetCalc_80N
        end if
        
        !BenC
!        call integrate_nod(Benthos(:,3), valBenSi, mesh)
         call integrate_bottom(Benthos(:,2)*ind_arctic_80_2D,valBenC_80N,mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of BenC at timestep :', mstep, valBenC_80N
        total_carbon_80N=total_carbon_80N+valBenC_80N
        
        !BenCalc
!        call integrate_nod(Benthos(:,3), valBenSi, mesh)
         call integrate_bottom(Benthos(:,4)*ind_arctic_80_2D,valBenCalc_80N,mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of BenCalc at timestep :', mstep, valBenCalc_80N
        total_carbon_80N=total_carbon_80N+valBenCalc_80N

        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of carbon (>80N) at timestep :', mstep, total_carbon_80N
        end if
        
        
        
        
        
        
        ! DIC
        call integrate_nod(tr_arr(:,:,4)*ind_arctic_66_3D, valDIC_66N, mesh)
        total_carbon_66N=total_carbon_66N+valDIC_66N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DIC (>66N) at timestep :', mstep, valDIC_66N
        end if

        ! DOC
        call integrate_nod(tr_arr(:,:,14)*ind_arctic_66_3D, valDOC_66N, mesh)
        total_carbon_66N=total_carbon_66N+valDOC_66N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DOC (>66N) at timestep :', mstep, valDOC_66N
        end if

        !PhyC
        call integrate_nod(tr_arr(:,:,7)*ind_arctic_66_3D, valPhyC_66N, mesh)
        total_carbon_66N=total_carbon_66N+valPhyC_66N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of PhyC (>66N) at timestep :', mstep, valPhyC_66N
        end if

        !DetC
        call integrate_nod(tr_arr(:,:,10)*ind_arctic_66_3D+tr_arr(:,:,28)*ind_arctic_66_3D, valDetC_66N, mesh)
        total_carbon_66N=total_carbon_66N+valDetC_66N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DetC+DetZ2C (>66N) at timestep :', mstep, valDetC_66N
        end if

        !HetC
        call integrate_nod(tr_arr(:,:,12)*ind_arctic_66_3D+tr_arr(:,:,26)*ind_arctic_66_3D, valHetC_66N, mesh)
        total_carbon_66N=total_carbon_66N+valHetC_66N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of HetC+Zoo2C (>66N) at timestep :', mstep, valHetC_66N
        end if

        !DiaC
        call integrate_nod(tr_arr(:,:,16)*ind_arctic_66_3D, valDiaC_66N, mesh)
        total_carbon_66N=total_carbon_66N+valDiaC_66N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DiaC (>66N) at timestep :', mstep, valDiaC_66N
        end if

       !PhyCalc
        call integrate_nod(tr_arr(:,:,22)*ind_arctic_66_3D, valPhyCalc_66N, mesh)
        total_carbon_66N=total_carbon_66N+valPhyCalc_66N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of PhyCalc (>80N) at timestep :', mstep, valPhyCalc_66N
        end if

        !DetCalc
        call integrate_nod(tr_arr(:,:,23)*ind_arctic_66_3D+tr_arr(:,:,30)*ind_arctic_66_3D, valDetCalc_66N, mesh)
        total_carbon_66N=total_carbon_66N+valDetCalc_66N
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of DetCalc+DetZ2Calc (>66N) at timestep :', mstep, valDetCalc_66N
        end if
        
        !BenC
!        call integrate_nod(Benthos(:,3), valBenSi, mesh)
         call integrate_bottom(Benthos(:,2)*ind_arctic_66_2D,valBenC_66N,mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of BenC at timestep :', mstep, valBenC_66N
        total_carbon_66N=total_carbon_66N+valBenC_66N
        
        !BenCalc
!        call integrate_nod(Benthos(:,3), valBenSi, mesh)
         call integrate_bottom(Benthos(:,4)*ind_arctic_66_2D,valBenCalc_66N,mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of BenCalc at timestep :', mstep, valBenCalc_66N
        total_carbon_66N=total_carbon_66N+valBenCalc_66N

        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) then
           write(*,*) 'total integral of carbon (>66N) at timestep :', mstep, total_carbon_66N
        end if

end subroutine compute_carbon_diag


subroutine compute_silicate_diag(mode,mesh)
  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use g_clock
  use o_PARAM
  USE o_ARRAYS
  use g_PARSUP
  use mod_MESH
  use g_comm_auto
  implicit none
  integer, intent(in)        :: mode
  type(t_mesh), intent(in)  , target :: mesh
  logical, save              :: firstcall=.true.
  integer                    :: k
#include  "../associate_mesh.h"
  
  ind_arctic_66_3D = 1.0_WP
  ind_arctic_66_2D = 1.0_WP
  do k=1, myDim_nod2D
      if (geo_coord_nod2D(2,k) < 66*rad) then
          ind_arctic_66_3D(:,k) = 0.0_WP
          ind_arctic_66_2D(k) = 0.0_WP
      end if
  end do
  
  ind_arctic_80_3D = 1.0_WP
  ind_arctic_80_2D = 1.0_WP
  do k=1, myDim_nod2D
      if (geo_coord_nod2D(2,k) < 80*rad) then
          ind_arctic_80_3D(:,k) = 0.0_WP
          ind_arctic_80_2D(k) = 0.0_WP
      end if
  end do

  if (firstcall) then  !allocate the stuff at the first call
    total_silicate_80N=0.0 ! accumulates
    total_silicate_66N=0.0 ! accumulates
    firstcall=.false.
    if (mode==0) return
  end if

        !DSi
        call integrate_nod(tr_arr(:,:,20)*ind_arctic_80_3D, valDSi_80N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of DSi at timestep :', mstep, valDSi_80N
        total_silicate_80N=total_silicate_80N+valDSi_80N

        !DiaSi
        call integrate_nod(tr_arr(:,:,18)*ind_arctic_80_3D, valDiaSi_80N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of DiaSi at timestep :', mstep, valDiaSi_80N
        total_silicate_80N=total_silicate_80N+valDiaSi_80N

        !DetSi
        call integrate_nod(tr_arr(:,:,19)*ind_arctic_80_3D+tr_arr(:,:,29)*ind_arctic_80_3D, valDetSi_80N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of DetSi+Detz2Si at timestep :', mstep, valDetSi_80N
        total_silicate_80N=total_silicate_80N+valDetSi_80N

        !BenSi
!        call integrate_nod(Benthos(:,3), valBenSi, mesh)
         call integrate_bottom(Benthos(:,3)*ind_arctic_80_2D,valBenSi_80N,mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of BenSi at timestep :', mstep, valBenSi_80N
        total_silicate_80N=total_silicate_80N+valBenSi_80N

        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of silicate at timestep :', mstep, total_silicate_80N
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !DSi
        call integrate_nod(tr_arr(:,:,20)*ind_arctic_66_3D, valDSi_66N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of DSi at timestep :', mstep, valDSi_80N
        total_silicate_80N=total_silicate_80N+valDSi_80N

        !DiaSi
        call integrate_nod(tr_arr(:,:,18)*ind_arctic_66_3D, valDiaSi_66N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of DiaSi at timestep :', mstep, valDiaSi_66N
        total_silicate_66N=total_silicate_66N+valDiaSi_66N

        !DetSi
        call integrate_nod(tr_arr(:,:,19)*ind_arctic_66_3D+tr_arr(:,:,29)*ind_arctic_66_3D, valDetSi_66N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of DetSi+Detz2Si at timestep :', mstep, valDetSi_66N
        total_silicate_66N=total_silicate_66N+valDetSi_66N

        !BenSi
!        call integrate_nod(Benthos(:,3), valBenSi, mesh)
         call integrate_bottom(Benthos(:,3)*ind_arctic_66_2D,valBenSi_66N,mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of BenSi at timestep :', mstep, valBenSi_66N
        total_silicate_66N=total_silicate_66N+valBenSi_66N

        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of silicate at timestep :', mstep, total_silicate_66N

end subroutine compute_silicate_diag


subroutine compute_nitrate_diag(mode,mesh)
  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use g_clock
  use o_PARAM
  USE o_ARRAYS
  use g_PARSUP
  use mod_MESH
  use g_comm_auto
  implicit none
  integer, intent(in)        :: mode
  type(t_mesh), intent(in)  , target :: mesh
  logical, save              :: firstcall=.true.
  integer                    :: k
#include  "../associate_mesh.h"
  
  ind_arctic_66_3D = 1.0_WP
  ind_arctic_66_2D = 1.0_WP
  do k=1, myDim_nod2D
      if (geo_coord_nod2D(2,k) < 66*rad) then
          ind_arctic_66_3D(:,k) = 0.0_WP
          ind_arctic_66_2D(k) = 0.0_WP
      end if
  end do
  
  ind_arctic_80_3D = 1.0_WP
  ind_arctic_80_2D = 1.0_WP
  do k=1, myDim_nod2D
      if (geo_coord_nod2D(2,k) < 80*rad) then
          ind_arctic_80_3D(:,k) = 0.0_WP
          ind_arctic_80_2D(k) = 0.0_WP
      end if
  end do

  if (firstcall) then  !allocate the stuff at the first call
    total_nitrate_80N=0.0 ! accumulates
    total_nitrate_66N=0.0 ! accumulates
    firstcall=.false.
    if (mode==0) return
  end if

        !DIN
        call integrate_nod(tr_arr(:,:,3)*ind_arctic_80_3D, valDIN_80N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of DIN at timestep :', mstep, valDIN_80N
        total_nitrate_80N=total_nitrate_80N+valDIN_80N
        
        ! DON
        call integrate_nod(tr_arr(:,:,13)*ind_arctic_80_3D, valDON_80N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of DON at timestep :', mstep, valDON_80N
        total_nitrate_80N=total_nitrate_80N+valDON_80N

        !DiaN
        call integrate_nod(tr_arr(:,:,15)*ind_arctic_80_3D, valDiaN_80N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of DiaN at timestep :', mstep, valDiaN_80N
        total_nitrate_80N=total_nitrate_80N+valDiaN_80N
        
        !PhyN
        call integrate_nod(tr_arr(:,:,6)*ind_arctic_80_3D, valPhyN_80N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of PhyN at timestep :', mstep, valPhyN_80N
        total_nitrate_80N=total_nitrate_80N+valPhyN_80N
        
        !HetN
        call integrate_nod(tr_arr(:,:,11)*ind_arctic_80_3D+tr_arr(:,:,25)*ind_arctic_80_3D, valHetN_80N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of Zoo2N+HetN at timestep :', mstep, valHetN_80N
        total_nitrate_80N=total_nitrate_80N+valHetN_80N

        !DetN
        call integrate_nod(tr_arr(:,:,9)*ind_arctic_80_3D+tr_arr(:,:,27)*ind_arctic_80_3D, valDetN_80N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of DetZ2N+DetN at timestep :', mstep, valDetN_80N
        total_nitrate_80N=total_nitrate_80N+valDetN_80N

        !BenN
!        call integrate_nod(Benthos(:,3), valBenSi, mesh)
         call integrate_bottom(Benthos(:,1)*ind_arctic_80_2D,valBenN_80N,mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of BenN at timestep :', mstep, valBenN_80N
        total_nitrate_80N=total_nitrate_80N+valBenN_80N

        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of nitrate at timestep :', mstep, total_nitrate_80N
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !DIN
        call integrate_nod(tr_arr(:,:,3)*ind_arctic_66_3D, valDIN_66N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of DIN at timestep :', mstep, valDIN_66N
        total_nitrate_66N=total_nitrate_66N+valDIN_66N
        
        ! DON
        call integrate_nod(tr_arr(:,:,13)*ind_arctic_66_3D, valDON_66N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of DON at timestep :', mstep, valDON_66N
        total_nitrate_66N=total_nitrate_66N+valDON_66N

        !DiaN
        call integrate_nod(tr_arr(:,:,15)*ind_arctic_66_3D, valDiaN_66N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of DiaN at timestep :', mstep, valDiaN_66N
        total_nitrate_66N=total_nitrate_66N+valDiaN_66N

        !PhyN
        call integrate_nod(tr_arr(:,:,6)*ind_arctic_66_3D, valPhyN_66N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of PhyN at timestep :', mstep, valPhyN_66N
        total_nitrate_66N=total_nitrate_66N+valPhyN_66N
        
        !HetN
        call integrate_nod(tr_arr(:,:,11)*ind_arctic_66_3D+tr_arr(:,:,25)*ind_arctic_66_3D, valHetN_66N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of Zoo2N+HetN at timestep :', mstep, valHetN_66N
        total_nitrate_66N=total_nitrate_66N+valHetN_66N

        !DetN
        call integrate_nod(tr_arr(:,:,9)*ind_arctic_66_3D+tr_arr(:,:,27)*ind_arctic_66_3D, valDetN_66N, mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>80N) of DetZ2N+DetN at timestep :', mstep, valDetN_80N
        total_nitrate_66N=total_nitrate_66N+valDetN_66N
!end if 
        !BenN
!        call integrate_nod(Benthos(:,3), valBenSi, mesh)
         call integrate_bottom(Benthos(:,3)*ind_arctic_66_2D,valBenN_66N,mesh)
        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of BenN at timestep :', mstep, valBenN_66N
        total_nitrate_66N=total_nitrate_66N+valBenN_66N

        if (mype==0 .and. mod(mstep,recom_logfile_outfreq)==0) write(*,*) 'total integral (>66N) of nitrate at timestep :', mstep, total_nitrate_66N

end subroutine compute_nitrate_diag

! ==============================================================
subroutine write_recom_diag(mode, mesh)

  implicit none
  integer                            :: status, ncid, j, k
  character(2000)                    :: filename
  integer                            :: recID, tID, tcID_66N, tsID_66N, tsdelID, tnID_66N,  tcID_80N, tsID_80N, tnID_80N
  integer                            :: valDICID_66N, valDOCID_66N, valPhyCID_66N, valDetCID_66N, valHetCID_66N, valDiaCID_66N, valPhyCalcID_66N, valDetCalcID_66N, valBenCID_66N, valBenCalcID_66N
  integer                            :: valDICID_80N, valDOCID_80N, valPhyCID_80N, valDetCID_80N, valHetCID_80N, valDiaCID_80N, valPhyCalcID_80N, valDetCalcID_80N, valBenCID_80N, valBenCalcID_80N
  integer                            :: valDSiID_66N, valDiaSiID_66N, valDetSiID_66N, valBenSiID_66N
  integer                            :: valDSiID_80N, valDiaSiID_80N, valDetSiID_80N, valBenSiID_80N
  integer                            :: valDINID_66N, valDONID_66N, valDiaNID_66N, valDetNID_66N, valDetz2NID_66N, valBenNID_66N, valPhyNID_66N, valHetNID_66N
  integer                            :: valDINID_80N, valDONID_80N, valDiaNID_80N, valDetNID_80N, valDetz2NID_80N, valBenNID_80N, valPhyNID_80N, valHetNID_80N
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
     status = nf_def_var(ncid, 'total_carbon_66N', NF_DOUBLE, 1, recID, tcID_66N)
     status = nf_def_var(ncid, 'total_DIC_66N', NF_DOUBLE, 1, recID, valDICID_66N)
     status = nf_def_var(ncid, 'total_DOC_66N', NF_DOUBLE, 1, recID, valDOCID_66N)
     status = nf_def_var(ncid, 'total_PhyC_66N', NF_DOUBLE, 1, recID, valPhyCID_66N)
     status = nf_def_var(ncid, 'total_DetC_66N', NF_DOUBLE, 1, recID, valDetCID_66N)
     status = nf_def_var(ncid, 'total_HetC_66N', NF_DOUBLE, 1, recID, valHetCID_66N)
     status = nf_def_var(ncid, 'total_DiaC_66N', NF_DOUBLE, 1, recID, valDiaCID_66N)
     status = nf_def_var(ncid, 'total_PhyCalc_66N', NF_DOUBLE, 1, recID, valPhyCalcID_66N)
     status = nf_def_var(ncid, 'total_DetCalc_66N', NF_DOUBLE, 1, recID, valDetCalcID_66N)
     status = nf_def_var(ncid, 'total_BenCalc_66N', NF_DOUBLE, 1, recID, valBenCalcID_66N)
     status = nf_def_var(ncid, 'total_BenC_66N', NF_DOUBLE, 1, recID, valBenCID_66N)
     
     status = nf_def_var(ncid, 'total_carbon_80N', NF_DOUBLE, 1, recID, tcID_80N)
     status = nf_def_var(ncid, 'total_DIC_80N', NF_DOUBLE, 1, recID, valDICID_80N)
     status = nf_def_var(ncid, 'total_DOC_80N', NF_DOUBLE, 1, recID, valDOCID_80N)
     status = nf_def_var(ncid, 'total_PhyC_80N', NF_DOUBLE, 1, recID, valPhyCID_80N)
     status = nf_def_var(ncid, 'total_DetC_80N', NF_DOUBLE, 1, recID, valDetCID_80N)
     status = nf_def_var(ncid, 'total_HetC_80N', NF_DOUBLE, 1, recID, valHetCID_80N)
     status = nf_def_var(ncid, 'total_DiaC_80N', NF_DOUBLE, 1, recID, valDiaCID_80N)
     status = nf_def_var(ncid, 'total_PhyCalc_80N', NF_DOUBLE, 1, recID, valPhyCalcID_80N)
     status = nf_def_var(ncid, 'total_DetCalc_80N', NF_DOUBLE, 1, recID, valDetCalcID_80N)
     status = nf_def_var(ncid, 'total_BenCalc_80N', NF_DOUBLE, 1, recID, valBenCalcID_80N)
     status = nf_def_var(ncid, 'total_BenC_80N', NF_DOUBLE, 1, recID, valBenCID_80N)
     
     status = nf_def_var(ncid, 'total_silicate_66N', NF_DOUBLE, 1, recID, tsID_66N)
     status = nf_def_var(ncid, 'total_DSi_66N', NF_DOUBLE, 1, recID, valDSiID_66N)
     status = nf_def_var(ncid, 'total_DiaSi_66N', NF_DOUBLE, 1, recID, valDiaSiID_66N)
     status = nf_def_var(ncid, 'total_DetSi_66N', NF_DOUBLE, 1, recID, valDetSiID_66N)
     status = nf_def_var(ncid, 'total_BenSi_66N', NF_DOUBLE, 1, recID, valBenSiID_66N)
     status = nf_def_var(ncid, 'total_silicate_80N', NF_DOUBLE, 1, recID, tsID_80N)
     status = nf_def_var(ncid, 'total_DSi_80N', NF_DOUBLE, 1, recID, valDSiID_80N)
     status = nf_def_var(ncid, 'total_DiaSi_80N', NF_DOUBLE, 1, recID, valDiaSiID_80N)
     status = nf_def_var(ncid, 'total_DetSi_80N', NF_DOUBLE, 1, recID, valDetSiID_80N)
     status = nf_def_var(ncid, 'total_BenSi_80N', NF_DOUBLE, 1, recID, valBenSiID_80N)
     
     status = nf_def_var(ncid, 'total_nitrate_66N', NF_DOUBLE, 1, recID, tnID_66N)
     status = nf_def_var(ncid, 'total_DIN_66N', NF_DOUBLE, 1, recID, valDINID_66N)
     status = nf_def_var(ncid, 'total_DON_66N', NF_DOUBLE, 1, recID, valDONID_66N)
     status = nf_def_var(ncid, 'total_DiaN_66N', NF_DOUBLE, 1, recID, valDiaNID_66N)
     status = nf_def_var(ncid, 'total_DetN_66N', NF_DOUBLE, 1, recID, valDetNID_66N)
     status = nf_def_var(ncid, 'total_BenN_66N', NF_DOUBLE, 1, recID, valBenNID_66N)
     status = nf_def_var(ncid, 'total_PhyN_66N', NF_DOUBLE, 1, recID, valPhyNID_66N)
     status = nf_def_var(ncid, 'total_HetN_66N', NF_DOUBLE, 1, recID, valHetNID_66N)
     status = nf_def_var(ncid, 'total_nitrate_80N', NF_DOUBLE, 1, recID, tnID_80N)
     status = nf_def_var(ncid, 'total_DIN_80N', NF_DOUBLE, 1, recID, valDINID_80N)
     status = nf_def_var(ncid, 'total_DON_80N', NF_DOUBLE, 1, recID, valDONID_80N)
     status = nf_def_var(ncid, 'total_DiaN_80N', NF_DOUBLE, 1, recID, valDiaNID_80N)
     status = nf_def_var(ncid, 'total_DetN_80N', NF_DOUBLE, 1, recID, valDetNID_80N)
     status = nf_def_var(ncid, 'total_BenN_80N', NF_DOUBLE, 1, recID, valBenNID_80N)
     status = nf_def_var(ncid, 'total_PhyN_80N', NF_DOUBLE, 1, recID, valPhyNID_80N)
     status = nf_def_var(ncid, 'total_HetN_80N', NF_DOUBLE, 1, recID, valHetNID_80N)

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
  
  status = nf_inq_varid(ncid, 'total_carbon_80N', tcID_80N)
  status = nf_inq_varid(ncid, 'total_DIC_80N', valDICID_80N)
  status = nf_inq_varid(ncid, 'total_DOC_80N', valDOCID_80N)
  status = nf_inq_varid(ncid, 'total_PhyC_80N', valPhyCID_80N)
  status = nf_inq_varid(ncid, 'total_DetC_80N', valDetCID_80N)
  status = nf_inq_varid(ncid, 'total_HetC_80N', valHetCID_80N)
  status = nf_inq_varid(ncid, 'total_DiaC_80N', valDiaCID_80N)
  status = nf_inq_varid(ncid, 'total_BenC_80N', valBenCID_80N)
  status = nf_inq_varid(ncid, 'total_BenCalc_80N', valBenCalcID_80N)
  status = nf_inq_varid(ncid, 'total_PhyCalc_80N', valPhyCalcID_80N)
  status = nf_inq_varid(ncid, 'total_DetCalc_80N', valDetCalcID_80N)
  
  status = nf_inq_varid(ncid, 'total_carbon_66N', tcID_66N)
  status = nf_inq_varid(ncid, 'total_DIC_66N', valDICID_66N)
  status = nf_inq_varid(ncid, 'total_DOC_66N', valDOCID_66N)
  status = nf_inq_varid(ncid, 'total_PhyC_66N', valPhyCID_66N)
  status = nf_inq_varid(ncid, 'total_DetC_66N', valDetCID_66N)
  status = nf_inq_varid(ncid, 'total_HetC_66N', valHetCID_66N)
  status = nf_inq_varid(ncid, 'total_DiaC_66N', valDiaCID_66N)
  status = nf_inq_varid(ncid, 'total_BenC_66N', valBenCID_66N)
  status = nf_inq_varid(ncid, 'total_BenCalc_66N', valBenCalcID_66N)
  status = nf_inq_varid(ncid, 'total_PhyCalc_66N', valPhyCalcID_66N)
  status = nf_inq_varid(ncid, 'total_DetCalc_66N', valDetCalcID_66N)

  status = nf_inq_varid(ncid, 'total_silicate_80N', tsID_80N)
  status = nf_inq_varid(ncid, 'total_DSi_80N', valDSiID_80N)
  status = nf_inq_varid(ncid, 'total_DiaSi_80N', valDiaSiID_80N)
  status = nf_inq_varid(ncid, 'total_DetSi_80N', valDetSiID_80N)
  status = nf_inq_varid(ncid, 'total_BenSi_80N', valBenSiID_80N)
  
  status = nf_inq_varid(ncid, 'total_silicate_66N', tsID_66N)
  status = nf_inq_varid(ncid, 'total_DSi_66N', valDSiID_66N)
  status = nf_inq_varid(ncid, 'total_DiaSi_66N', valDiaSiID_66N)
  status = nf_inq_varid(ncid, 'total_DetSi_66N', valDetSiID_66N)
  status = nf_inq_varid(ncid, 'total_BenSi_66N', valBenSiID_66N)
  
  status = nf_inq_varid(ncid, 'total_nitrate_66N', tnID_66N)
  status = nf_inq_varid(ncid, 'total_DIN_66N', valDINID_66N)
  status = nf_inq_varid(ncid, 'total_DON_66N', valDONID_66N)
  status = nf_inq_varid(ncid, 'total_DiaN_66N', valDiaNID_66N)
  status = nf_inq_varid(ncid, 'total_DetN_66N', valDetNID_66N)
  status = nf_inq_varid(ncid, 'total_BenN_66N', valBenNID_66N)
  status = nf_inq_varid(ncid, 'total_PhyN_66N', valPhyNID_66N)
  status = nf_inq_varid(ncid, 'total_HetN_66N', valHetNID_66N)
  
  status = nf_inq_varid(ncid, 'total_nitrate_80N', tnID_80N)
  status = nf_inq_varid(ncid, 'total_DIN_80N', valDINID_80N)
  status = nf_inq_varid(ncid, 'total_DON_80N', valDONID_80N)
  status = nf_inq_varid(ncid, 'total_DiaN_80N', valDiaNID_80N)
  status = nf_inq_varid(ncid, 'total_DetN_80N', valDetNID_80N)
  status = nf_inq_varid(ncid, 'total_BenN_80N', valBenNID_80N)
  status = nf_inq_varid(ncid, 'total_PhyN_80N', valPhyNID_80N)
  status = nf_inq_varid(ncid, 'total_HetN_80N', valHetNID_80N)

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
  
  status = nf_put_vara_double(ncid, tcID_66N, rec_count, 1, total_carbon_66N, 1)
  status = nf_put_vara_double(ncid, valDICID_66N, rec_count, 1, valDIC_66N, 1)
  status = nf_put_vara_double(ncid, valDOCID_66N, rec_count, 1, valDOC_66N, 1)
  status = nf_put_vara_double(ncid, valPhyCID_66N, rec_count, 1, valPhyC_66N, 1)
  status = nf_put_vara_double(ncid, valDetCID_66N, rec_count, 1, valDetC_66N, 1)
  status = nf_put_vara_double(ncid, valHetCID_66N, rec_count, 1, valHetC_66N, 1)
  status = nf_put_vara_double(ncid, valDiaCID_66N, rec_count, 1, valDiaC_66N, 1)
  status = nf_put_vara_double(ncid, valPhyCalcID_66N, rec_count, 1, valPhyCalc_66N, 1)
  status = nf_put_vara_double(ncid, valDetCalcID_66N, rec_count, 1, valDetCalc_66N, 1)
  status = nf_put_vara_double(ncid, valBenCalcID_66N, rec_count, 1, valBenCalc_66N, 1)
  status = nf_put_vara_double(ncid, valBenCID_66N, rec_count, 1, valBenC_66N, 1)
  
  status = nf_put_vara_double(ncid, tcID_80N, rec_count, 1, total_carbon_80N, 1)
  status = nf_put_vara_double(ncid, valDICID_80N, rec_count, 1, valDIC_80N, 1)
  status = nf_put_vara_double(ncid, valDOCID_80N, rec_count, 1, valDOC_80N, 1)
  status = nf_put_vara_double(ncid, valPhyCID_80N, rec_count, 1, valPhyC_80N, 1)
  status = nf_put_vara_double(ncid, valDetCID_80N, rec_count, 1, valDetC_80N, 1)
  status = nf_put_vara_double(ncid, valHetCID_80N, rec_count, 1, valHetC_80N, 1)
  status = nf_put_vara_double(ncid, valDiaCID_80N, rec_count, 1, valDiaC_80N, 1)
  status = nf_put_vara_double(ncid, valPhyCalcID_80N, rec_count, 1, valPhyCalc_80N, 1)
  status = nf_put_vara_double(ncid, valDetCalcID_80N, rec_count, 1, valDetCalc_80N, 1)
  status = nf_put_vara_double(ncid, valBenCalcID_80N, rec_count, 1, valBenCalc_80N, 1)
  status = nf_put_vara_double(ncid, valBenCID_80N, rec_count, 1, valBenC_80N, 1)

  status = nf_put_vara_double(ncid, tsID_66N, rec_count, 1, total_silicate_66N, 1)
  status = nf_put_vara_double(ncid, valDSiID_66N, rec_count, 1, valDSi_66N, 1)
  status = nf_put_vara_double(ncid, valDiaSiID_66N, rec_count, 1, valDiaSi_66N, 1)
  status = nf_put_vara_double(ncid, valDetSiID_66N, rec_count, 1, valDetSi_66N, 1)
  status = nf_put_vara_double(ncid, valBenSiID_66N, rec_count, 1, valBenSi_66N, 1)
  
  status = nf_put_vara_double(ncid, tsID_80N, rec_count, 1, total_silicate_80N, 1)
  status = nf_put_vara_double(ncid, valDSiID_80N, rec_count, 1, valDSi_80N, 1)
  status = nf_put_vara_double(ncid, valDiaSiID_80N, rec_count, 1, valDiaSi_80N, 1)
  status = nf_put_vara_double(ncid, valDetSiID_80N, rec_count, 1, valDetSi_80N, 1)
  status = nf_put_vara_double(ncid, valBenSiID_80N, rec_count, 1, valBenSi_80N, 1)
  
  status = nf_put_vara_double(ncid, tnID_66N, rec_count, 1, total_nitrate_66N, 1)
  status = nf_put_vara_double(ncid, valDINID_66N, rec_count, 1, valDIN_66N, 1)
  status = nf_put_vara_double(ncid, valDONID_66N, rec_count, 1, valDON_66N, 1)
  status = nf_put_vara_double(ncid, valDiaNID_66N, rec_count, 1, valDiaN_66N, 1)
  status = nf_put_vara_double(ncid, valDetNID_66N, rec_count, 1, valDetN_66N, 1)
  status = nf_put_vara_double(ncid, valBenNID_66N, rec_count, 1, valBenN_66N, 1)
  status = nf_put_vara_double(ncid, valPhyNID_66N, rec_count, 1, valPhyN_66N, 1)
  status = nf_put_vara_double(ncid, valHetNID_66N, rec_count, 1, valHetN_66N, 1)
  
  status = nf_put_vara_double(ncid, tnID_80N, rec_count, 1, total_nitrate_80N, 1)
  status = nf_put_vara_double(ncid, valDINID_80N, rec_count, 1, valDIN_80N, 1)
  status = nf_put_vara_double(ncid, valDONID_80N, rec_count, 1, valDON_80N, 1)
  status = nf_put_vara_double(ncid, valDiaNID_80N, rec_count, 1, valDiaN_80N, 1)
  status = nf_put_vara_double(ncid, valDetNID_80N, rec_count, 1, valDetN_80N, 1)
  status = nf_put_vara_double(ncid, valBenNID_80N, rec_count, 1, valBenN_80N, 1)
  status = nf_put_vara_double(ncid, valPhyNID_80N, rec_count, 1, valPhyN_80N, 1)
  status = nf_put_vara_double(ncid, valHetNID_80N, rec_count, 1, valHetN_80N, 1)

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
  !3. nitrate diagnostic
  if (ldiag_nitrate)    call compute_nitrate_diag(mode,mesh)
  !4. write total carbon and silicate out into recom.diag.nc
  if (ldiag_carbon .or. ldiag_silicate .or. ldiag_nitrate) call write_recom_diag(mode, mesh)

end subroutine compute_recom_diagnostics


end module recom_diag