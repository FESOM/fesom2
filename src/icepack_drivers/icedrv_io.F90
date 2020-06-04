!=======================================================================
!
! This submodule initializes the IO subroutines
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!
!=======================================================================

      submodule (icedrv_main) icedrv_io

      use icepack_intfc,    only: icepack_query_parameters
      use icepack_intfc,    only: icepack_query_tracer_flags
      use icepack_intfc,    only: icepack_query_tracer_sizes
      use icepack_intfc,    only: icepack_query_tracer_indices
      use icedrv_system,    only: icedrv_system_abort
      use io_meandata,      only: def_stream3D 

      contains

      module subroutine init_io_icepack(mesh)

          use mod_mesh
          use g_parsup

          implicit none

          type(t_mesh), target, intent(in) :: mesh
        
          integer           :: i, j, k,                            &
                               nt_Tsfc, nt_sice, nt_qice, nt_qsno, &
                               nt_apnd, nt_hpnd, nt_ipnd, nt_alvl, &
                               nt_vlvl, nt_iage, nt_FY,   nt_aero, &
                               ktherm,  nt_fbri
        
          integer, save     :: nm_io_unit  = 102       ! unit to open namelist file
          integer, save     :: nm_icepack_unit = 103
          integer           :: iost
          character(len=10) :: id_string
          character(500)    :: longname, trname, units
        
          logical (kind=log_kind)   ::                        &
               solve_zsal, skl_bgc, z_tracers,                &
               tr_iage, tr_FY, tr_lvl, tr_aero, tr_pond_cesm, &
               tr_pond_topo, tr_pond_lvl, tr_brine,           &
               tr_bgc_N, tr_bgc_C, tr_bgc_Nit,                &
               tr_bgc_Sil,  tr_bgc_DMS,                       &
               tr_bgc_chl,  tr_bgc_Am,                        &
               tr_bgc_PON,  tr_bgc_DON,                       &
               tr_zaero,    tr_bgc_Fe,                        &
               tr_bgc_hum

          integer, save                  :: io_listsize=0
          type io_entry
                character(len=10)        :: id        ='unknown   '
                integer                  :: freq      =0
                character                :: unit      =''
                integer                  :: precision =0
          end type
        
          type(io_entry), save, allocatable, target   :: io_list_icepack(:)

#include "../associate_mesh.h"

          ! Get the tracers information from icepack
          call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
               nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
          call icepack_query_tracer_indices(                                          &
               nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd,         &
               nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_Tsfc_out=nt_Tsfc,         &
               nt_iage_out=nt_iage, nt_FY_out=nt_FY,                                  &
               nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,         &
               nt_aero_out=nt_aero, nt_qsno_out=nt_qsno)
          call icepack_query_parameters(solve_zsal_out=solve_zsal,                    &
               skl_bgc_out=skl_bgc, z_tracers_out=z_tracers, ktherm_out=ktherm)
          call icepack_query_tracer_flags(                                            &
               tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_lvl_out=tr_lvl,               &
               tr_aero_out=tr_aero, tr_pond_cesm_out=tr_pond_cesm,                    &
               tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl,            &
               tr_brine_out=tr_brine, tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C,   &
               tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_Sil_out=tr_bgc_Sil,                  &
               tr_bgc_DMS_out=tr_bgc_DMS,                                             &
               tr_bgc_chl_out=tr_bgc_chl, tr_bgc_Am_out=tr_bgc_Am,                    &
               tr_bgc_PON_out=tr_bgc_PON, tr_bgc_DON_out=tr_bgc_DON,                  &
               tr_zaero_out=tr_zaero,     tr_bgc_Fe_out=tr_bgc_Fe,                    &
               tr_bgc_hum_out=tr_bgc_hum)
        
          namelist /nml_listsize      / io_listsize
          namelist /nml_list_icepack  / io_list_icepack
        
          ! OPEN and read namelist for icepack I/O
          open( unit=nm_io_unit, file='namelist.io', form='formatted', access='sequential', status='old', iostat=iost )
          if (iost == 0) then
          if (mype==0) write(*,*) '     file   : ', 'namelist.io',' open ok'
             else
          if (mype==0) write(*,*) 'ERROR: --> bad opening file   : ','namelist.io',' ;    iostat=',iost
             call par_ex
             stop
          end if
          open( unit=nm_icepack_unit, file='namelist.icepack', form='formatted', access='sequential', status='old', iostat=iost )
          if (iost == 0) then
          if (mype==0) write(*,*) '     file   : ', 'namelist.icepack',' open ok' 
             else
          if (mype==0) write(*,*) 'ERROR: --> bad opening file   : ','namelist.icepack',' ;    iostat=',iost
             call par_ex
             stop
          end if

          read(nm_io_unit, nml=nml_listsize, iostat=iost )
          allocate(io_list_icepack(io_listsize))
          read(nm_icepack_unit, nml=nml_list_icepack, iostat=iost )
          close(nm_icepack_unit)
        
        
          do i=1, io_listsize
             if (trim(io_list_icepack(i)%id)=='unknown   ') then
                if (mype==0) write(*,*) 'io_listsize will be changed from ', io_listsize, ' to ', i-1, '!'
                io_listsize=i-1
                exit
             end if
          end do
        
          do i=1, io_listsize
             select case (trim(io_list_icepack(i)%id))
             case ('aicen     ')
                 call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'aicen', 'sea ice concentration',     'none', aicen(:,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh) 
             case ('vicen     ')
                 call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'vicen', 'volume per unit area of ice',  'm', vicen(:,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             case ('vsnon     ')
                 call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'vsnon', 'volume per unit area of snow', 'm', vsnon(:,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             case ('Tsfc      ')
                 call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'Tsfc',  'sea ice surf. temperature', 'degC', trcrn(:,nt_Tsfc,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             ! If the following tracers are not defined they will not be outputed
             case ('iage      ')
                if (tr_iage) then
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'iage', 'sea ice age', 's', trcrn(:,nt_iage,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('FY        ')
                if (tr_FY) then 
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'FY', 'first year ice', 'none', trcrn(:,nt_FY,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('lvl       ')
                if (tr_lvl) then
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'alvl', 'ridged sea ice area',   'none', trcrn(:,nt_alvl,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'vlvl', 'ridged sea ice volume', 'm',    trcrn(:,nt_vlvl,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('pond_cesm ')
                if (tr_pond_cesm) then
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'apnd', 'melt ponds area',   'none', trcrn(:,nt_apnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh) 
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'vpnd', 'melt ponds volume', 'm',    trcrn(:,nt_hpnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh) 
                end if
             case ('pond_topo ')
                if (tr_pond_topo) then
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'apnd', 'melt ponds area',                   'none', trcrn(:,nt_apnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'vpnd', 'melt ponds volume',                 'm',    trcrn(:,nt_hpnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'ipnd', 'melt ponds refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('pond_lvl  ')
                if (tr_pond_lvl) then
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'apnd', 'melt ponds area',                   'none', trcrn(:,nt_apnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'vpnd', 'melt ponds volume',                 'm',    trcrn(:,nt_hpnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'ipnd', 'melt ponds refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('brine     ')
                if (tr_brine) then
                  call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'fbri', 'volume fraction of ice with dynamic salt', 'none',    trcrn(:,nt_fbri,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('qicen     ')
                 do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'qicen_', k
                    write(longname,'(A22,i1)') 'sea ice enthalpy lyr: ', k 
                    units='J/m3'
                    call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/), trim(trname), trim(longname), trim(units), trcrn(:,nt_qice+k-1,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                 end do
             case ('sicen     ')
                 do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'sicen_', k
                    write(longname,'(A22,i1)') 'sea ice salinity lyr: ', k
                    units='psu'
                    call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/), trim(trname), trim(longname), trim(units), trcrn(:,nt_sice+k-1,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                 end do
             case ('qsnon     ')
                 do k = 1,nslyr  ! Separate variable for each snow layer
                    write(trname,'(A6,i1)') 'qsnon_', k
                    write(longname,'(A19,i1)') 'snow enthalpy lyr: ', k
                    units='J/m3'
                    call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/), trim(trname), trim(longname), trim(units), trcrn(:,nt_qsno+k-1,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                 end do
             case default
                 if (mype==0) write(*,*) 'stream ', io_list_icepack(i)%id, ' is not defined !'
             end select
          end do

      end subroutine init_io_icepack

      end submodule icedrv_io
