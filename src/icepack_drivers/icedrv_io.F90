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
      use icepack_intfc,    only: icepack_warnings_flush
      use icepack_intfc,    only: icepack_warnings_aborted                  
      use icedrv_system,    only: icedrv_system_abort
 
      contains

    !
    !
    !_________________________________________________________________________
    ! define mean IO output of icepack 
    module subroutine ini_mean_icepack_io(mesh)

        use mod_mesh
        use io_meandata,      only: def_stream

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

        namelist /nml_general       / io_listsize
        namelist /nml_list_icepack  / io_list_icepack        

#include "associate_mesh.h"

        ! Get the tracers information from icepack
        call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
            nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
        cAll icepack_query_tracer_indices(                                          &
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
        
        !_______________________________________________________________________
        ! OPEN and read namelist.io --> need to extract variable io_listsize
        open( unit=nm_io_unit, file='namelist.io', form='formatted', access='sequential', status='old', iostat=iost )
        if (iost == 0) then
            if (mype==0) write(*,*) '     file   : ', 'namelist.io',' open ok'
        else
            if (mype==0) write(*,*) 'ERROR: --> bad opening file   : ','namelist.io',' ;    iostat=',iost
            call par_ex(p_partit%MPI_COMM_FESOM, p_partit%mype)
            stop
        end if
        
        ! read list_size from namelist.io for allocation
        read(nm_io_unit, nml=nml_general, iostat=iost )
        close(nm_io_unit)
        allocate(io_list_icepack(io_listsize))
        
        !_______________________________________________________________________
        ! OPEN and read namelist.icepack --> need to extract io_list_icepack
        open( unit=nm_icepack_unit, file='namelist.icepack', form='formatted', access='sequential', status='old', iostat=iost )
        if (iost == 0) then
            if (mype==0) write(*,*) '     file   : ', 'namelist.icepack',' open ok' 
        else
            if (mype==0) write(*,*) 'ERROR: --> bad opening file   : ','namelist.icepack',' ;    iostat=',iost
            call par_ex(p_partit%MPI_COMM_FESOM, p_partit%mype)
            stop
        end if
        
        ! read io_list_icepack from namelist to fill up what has been previously 
        ! allocated --> allocate(io_list_icepack(io_listsize))
        read(nm_icepack_unit, nml=nml_list_icepack, iostat=iost )
        close(nm_icepack_unit)
        
        !_______________________________________________________________________
        ! reduce running index to the number that is actually filt up 
        do i=1, io_listsize
            if (trim(io_list_icepack(i)%id)=='unknown   ') then
                if (mype==0) write(*,*) 'io_listsize will be changed from ', io_listsize, ' to ', i-1, '!'
                io_listsize=i-1
                exit
            end if
        end do
        
        !_______________________________________________________________________
        ! define output streams
        do i=1, io_listsize
            select case (trim(io_list_icepack(i)%id))
            case ('aice0     ')
                call def_stream(nod2D          , nx_nh          , 'aice0'   , 'open water fraction'         , 'none', aice0(:)          , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh) 
            case ('aicen     ')
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'aicen'   , 'sea ice concentration'       , 'none', aicen(:,:)        , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.) 
            case ('vicen     ')
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'vicen'   , 'volume per unit area of ice' , 'm'   , vicen(:,:)        , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
            case ('vsnon     ')
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'vsnon'   , 'volume per unit area of snow', 'm'   , vsnon(:,:)        , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
            case ('aice      ')
                call def_stream(nod2D          , nx_nh          , 'aice'    , 'sea ice concentration'       , 'none', aice(:)           , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh) 
            case ('vice      ')
                call def_stream(nod2D          , nx_nh          , 'vice'    , 'volume per unit area of ice' , 'm'   , vice(:)           , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
            case ('vsno      ')
                call def_stream(nod2D          , nx_nh          , 'vsno'    , 'volume per unit area of snow', 'm'   , vsno(:)           , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
            ! Sea ice velocity components
            case ('uvel      ')
                call def_stream(nod2D          , nx_nh          , 'uvel'    , 'x-component of ice velocity' , 'm/s' , uvel(:)           , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
            case ('vvel      ')
                call def_stream(nod2D          , nx_nh          , 'vvel'    , 'y-component of ice velocity' , 'm/s' , vvel(:)           , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
            ! Sea ice or snow surface temperature
            case ('Tsfc      ')
                call def_stream(nod2D          , nx_nh          , 'Tsfc'    , 'sea ice surf. temperature'   , 'degC', trcr(:,nt_Tsfc)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
            case ('Tsfcn     ')
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'Tsfcn'   , 'sea ice surf. temperature'   , 'degC', trcrn(:,nt_Tsfc,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
            case ('strength  ') 
                call def_stream(nod2D          , nx_nh          , 'strength', 'sea ice strength'            , 'N'   , strength(:)       , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh) 
            ! If the following tracers are not defined they will not be outputed
            case ('iagen     ')
                if (tr_iage) then
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'iage'    , 'sea ice age'                 , 's'   , trcrn(:,nt_iage,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                end if
            case ('FYn       ')
                if (tr_FY) then 
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'FY'      , 'first year ice'              , 'none', trcrn(:,nt_FY,:)  , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                end if
            case ('lvln      ')
                if (tr_lvl) then
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'alvl'    , 'ridged sea ice area'         , 'none', trcrn(:,nt_alvl,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'vlvl'    , 'ridged sea ice volume'       , 'm'   , trcrn(:,nt_vlvl,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                end if
            case ('pond_cesmn')
                if (tr_pond_cesm) then
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'apnd'    , 'melt pond area fraction'     , 'none', trcrn(:,nt_apnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.) 
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'hpnd'    , 'melt pond depth'             , 'm'   , trcrn(:,nt_hpnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.) 
                end if
            case ('pond_topon')
                if (tr_pond_topo) then
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'apnd'    , 'melt pond area fraction'     , 'none', trcrn(:,nt_apnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'hpnd'    , 'melt pond depth'             , 'm'   , trcrn(:,nt_hpnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'ipnd'    , 'melt pond refrozen lid thickness', 'm',trcrn(:,nt_ipnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                end if
            case ('pond_lvln ')
                if (tr_pond_lvl) then
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'apnd'    , 'melt pond area fraction'     , 'none', trcrn(:,nt_apnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'hpnd'    , 'melt pond depth'             , 'm'   , trcrn(:,nt_hpnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                call def_stream((/ncat, nod2D/), (/ncat, nx_nh/), 'ipnd'    , 'melt pond refrozen lid thickness','m', trcrn(:,nt_ipnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                end if
            case ('brinen    ')
                if (tr_brine) then
                call def_stream((/ncat, nod2D/),  (/ncat, nx_nh/), 'fbri'   , 'volume fraction of ice with dynamic salt', 'none', trcrn(:,nt_fbri,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                end if
            case ('qicen     ')
                do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'qicen_', k
                    write(longname,'(A22,i1)') 'sea ice enthalpy lyr: ', k 
                    units='J/m3'
                    call def_stream((/ncat, nod2D/),  (/ncat, nx_nh/), trim(trname), trim(longname), trim(units), trcrn(:,nt_qice+k-1,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                end do
            case ('sicen     ')
                do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'sicen_', k
                    write(longname,'(A22,i1)') 'sea ice salinity lyr: ', k
                    units='psu'
                    call def_stream((/ncat, nod2D/),  (/ncat, nx_nh/), trim(trname), trim(longname), trim(units), trcrn(:,nt_sice+k-1,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                end do
            case ('qsnon     ')
                do k = 1,nslyr  ! Separate variable for each snow layer
                    write(trname,'(A6,i1)') 'qsnon_', k
                    write(longname,'(A19,i1)') 'snow enthalpy lyr: ', k
                    units='J/m3'
                    call def_stream((/ncat, nod2D/),  (/ncat, nx_nh/), trim(trname), trim(longname), trim(units), trcrn(:,nt_qsno+k-1,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh, .true.)
                end do
            ! Average over categories
            case ('iage      ')
                if (tr_iage) then
                    call def_stream(nod2D      , nx_nh          , 'iage'    , 'sea ice age'                 , 's'   , trcr(:,nt_iage)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                end if
            case ('FY        ')
                if (tr_FY) then 
                    call def_stream(nod2D      , nx_nh          , 'FY'      , 'first year ice'              , 'none', trcr(:,nt_FY)     , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                end if
            case ('lvl       ')
                if (tr_lvl) then
                    call def_stream(nod2D      , nx_nh          , 'alvl'    , 'ridged sea ice area'         , 'none', trcr(:,nt_alvl)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                    call def_stream(nod2D      , nx_nh          , 'vlvl'    , 'ridged sea ice volume'       , 'm'   , trcr(:,nt_vlvl)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                end if
            case ('pond_cesm ')
                if (tr_pond_cesm) then
                    call def_stream(nod2D      , nx_nh          , 'apnd'    , 'melt pond area fraction'     , 'none', trcr(:,nt_apnd)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh) 
                    call def_stream(nod2D      , nx_nh          , 'hpnd'    , 'melt pond depth'             , 'm'   , trcr(:,nt_hpnd)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh) 
                end if
            case ('pond_topo ')
                if (tr_pond_topo) then
                    call def_stream(nod2D      , nx_nh          , 'apnd'    , 'melt pond area fraction'     , 'none', trcr(:,nt_apnd)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                    call def_stream(nod2D      , nx_nh          , 'hpnd'    , 'melt pond depth'             , 'm'   , trcr(:,nt_hpnd)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                    call def_stream(nod2D      , nx_nh          , 'ipnd'    , 'melt pond refrozen lid thickness', 'm', trcr(:,nt_ipnd)  , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                end if
            case ('pond_lvl  ')
                if (tr_pond_lvl) then
                    call def_stream(nod2D      , nx_nh          , 'apnd'    , 'melt pond area fraction'     , 'none', trcr(:,nt_apnd)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                    call def_stream(nod2D      , nx_nh          , 'hpnd'    , 'melt pond depth'             , 'm'   , trcr(:,nt_hpnd)   , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                    !call def_stream(nod2D,  nx_nh,  'ipnd', 'melt pond refrozen lid thickness', 'm',    trcr(:,nt_ipnd), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                end if
            case ('brine     ')
                if (tr_brine) then
                    call def_stream(nod2D      , nx_nh          , 'fbri'    , 'volume fraction of ice with dynamic salt', 'none', trcr(:,nt_fbri), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                end if
            case ('qice      ')
                do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'qicen_', k
                    write(longname,'(A22,i1)') 'sea ice enthalpy lyr: ', k 
                    units='J/m3'
                    call def_stream(nod2D      , nx_nh          , trim(trname), trim(longname)              , trim(units), trcr(:,nt_qice+k-1), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                end do
            case ('sice      ')
                do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'sicen_', k
                    write(longname,'(A22,i1)') 'sea ice salinity lyr: ', k
                    units='psu'
                    call def_stream(nod2D      , nx_nh          , trim(trname), trim(longname)              , trim(units), trcr(:,nt_sice+k-1), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                end do
            case ('qsno      ')
                do k = 1,nslyr  ! Separate variable for each snow layer
                    write(trname,'(A6,i1)') 'qsnon_', k
                    write(longname,'(A19,i1)') 'snow enthalpy lyr: ', k
                    units='J/m3'
                    call def_stream(nod2D      , nx_nh          , trim(trname), trim(longname)              , trim(units), trcr(:,nt_qsno+k-1), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
                end do
            case ('rdg_conv  ')
                call def_stream(nod2D          , nx_nh          , 'rdg_conv' , 'Convergence term for ridging', '1/s', rdg_conv(:)       , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
            case ('rdg_shear ')
                call def_stream(nod2D          , nx_nh          , 'rdg_shear', 'Shear term for ridging'     , '1/s' , rdg_shear(:)      , io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, p_partit, mesh)
            case default
                if (mype==0) write(*,*) 'stream ', io_list_icepack(i)%id, ' is not defined !'
            end select
        end do ! --> do i=1, io_listsize

    end subroutine ini_mean_icepack_io

    !
    !
    !___________________________________________________________________________
    ! define mean IO output of icepack 
    module subroutine ini_icepack_io(year, partit, mesh)
        
        use mod_mesh
        use mod_partit
        use mod_parsup
        use g_config,     only: runid, ResultPath, RestartOutPath
        use io_restart,   only: icepack_files, icepack_path, nc_restart_path
    
        implicit none
    
        type(t_mesh)  , intent(in)   , target :: mesh
        type(t_partit), intent(inout), target :: partit
        integer       , intent(in)            :: year
        logical       , save                  :: has_been_called = .false.
        
        integer (kind=int_kind)   :: i, j, k, iblk, &     ! counting indices
                                     nt_Tsfc, nt_sice, nt_qice, nt_qsno,    &
                                     nt_apnd, nt_hpnd, nt_ipnd, nt_alvl,    &
                                     nt_vlvl, nt_iage, nt_FY,   nt_aero,    &
                                     ktherm,  nt_fbri
        character(500)            :: longname
        character(500)            :: trname, units
        character(4)              :: cyear
    
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
    
#include "associate_mesh.h"
    
        ! Get the tracers information from icepack
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
        call icepack_warnings_flush(nu_diag)
        ! The following error message needs to be fixed
        !if (icepack_warnings_aborted()) call abort_ice(error_message=subname,       &
        !    file=__FILE__, line=__LINE__)
      
        write(cyear,'(i4)') year
        ! Create an icepack restart file
        ! Only serial output implemented so far
        icepack_path = nc_restart_path('icepack', year, RestartOutPath)
        
        if(has_been_called) return
        has_been_called = .true.
        
        ! Define the netCDF variables for surface
        ! and vertically constant fields
      
        !-----------------------------------------------------------------
        ! 3D restart fields (ncat)
        !-----------------------------------------------------------------
        call icepack_files%def_node_var('aice'     , 'sea ice concentration'                        , 'none', aice(:)           , mesh, partit)
        call icepack_files%def_node_var('vice'     , 'volum per unit area of ice'                   , 'm'   , vice(:)           , mesh, partit)
        call icepack_files%def_node_var('vsno'     , 'volum per unit area of snow'                  , 'm'   , vsno(:)           , mesh, partit)
        
        call icepack_files%def_node_var('aicen'     , 'sea ice concentration per class'             , 'none', aicen(:,:)        , mesh, partit, ncat)
        call icepack_files%def_node_var('vicen'     , 'volum per unit area of ice per class'        , 'm'   , vicen(:,:)        , mesh, partit, ncat)
        call icepack_files%def_node_var('vsnon'     , 'volum per unit area of snow per class'       , 'm'   , vsnon(:,:)        , mesh, partit, ncat)
        call icepack_files%def_node_var('Tsfc'      , 'sea ice surf. temperature'                   , 'degC', trcrn(:,nt_Tsfc,:), mesh, partit, ncat)
        call icepack_files%def_node_var('uvel'      , 'zonal component of ice velocity'             , 'm/s' , uvel(:)           , mesh, partit)
        call icepack_files%def_node_var('vvel'      , 'meridional component of ice velocity'        , 'm/s' , vvel(:)           , mesh, partit)
      
        if (tr_iage) then
            call icepack_files%def_node_var('iage'  , 'sea ice age'                                 , 's'   , trcrn(:,nt_iage,:), mesh, partit, ncat)
        end if
      
        if (tr_FY) then
            call icepack_files%def_node_var('FY'    , 'first year ice'                              , 'none', trcrn(:,nt_FY,:)  , mesh, partit, ncat)
        end if
      
        if (tr_lvl) then
            call icepack_files%def_node_var('alvl'  , 'ridged sea ice area'                         , 'none', trcrn(:,nt_alvl,:), mesh, partit, ncat)
            call icepack_files%def_node_var('vlvl'  , 'ridged sea ice volume'                       , 'm'   , trcrn(:,nt_vlvl,:), mesh, partit, ncat)
        end if
      
        if (tr_pond_cesm) then
            call icepack_files%def_node_var('apnd'  , 'melt pond area fraction'                     , 'none', trcrn(:,nt_apnd,:), mesh, partit, ncat)
            call icepack_files%def_node_var('hpnd'  , 'melt pond depth'                             , 'm'   , trcrn(:,nt_hpnd,:), mesh, partit, ncat)
        end if
      
        if (tr_pond_topo) then
            call icepack_files%def_node_var('apnd'  , 'melt pond area fraction'                     , 'none', trcrn(:,nt_apnd,:), mesh, partit, ncat)
            call icepack_files%def_node_var('hpnd'  , 'melt pond depth'                             , 'm'   , trcrn(:,nt_hpnd,:), mesh, partit, ncat)
            call icepack_files%def_node_var('ipnd'  , 'melt pond refrozen lid thickness'            , 'm'   , trcrn(:,nt_ipnd,:), mesh, partit, ncat)
        end if
      
        if (tr_pond_lvl) then
            call icepack_files%def_node_var('apnd'  , 'melt pond area fraction'                     , 'none', trcrn(:,nt_apnd,:), mesh, partit, ncat)
            call icepack_files%def_node_var('hpnd'  , 'melt pond depth'                             , 'm'   , trcrn(:,nt_hpnd,:), mesh, partit, ncat)
            call icepack_files%def_node_var('ipnd'  , 'melt pond refrozen lid thickness'            , 'm'   , trcrn(:,nt_ipnd,:), mesh, partit, ncat)
            call icepack_files%def_node_var('ffracn', 'fraction of fsurfn over pond used to melt ipond'     , 'none', ffracn, mesh, partit);
            call icepack_files%def_node_var('dhsn'  ,  'depth difference for snow on sea ice and pond ice'  ,  'm'  , dhsn  , mesh, partit);
        end if
      
        if (tr_brine) then
            call icepack_files%def_node_var('fbri'  ,     'volume fraction of ice with dynamic salt', 'none', trcrn(:,nt_fbri,:), mesh, partit, ncat)
            call icepack_files%def_node_var('first_ice', 'distinguishes ice that disappears'        , 'logical', first_ice_real(:,:), mesh, partit, ncat)
        end if
      
        !-----------------------------------------------------------------
        ! 4D restart fields, written as layers of 3D
        !-----------------------------------------------------------------
      
        ! Ice
      
        do k = 1,nilyr
           write(trname,'(A6,i1)') 'sicen_', k
           write(longname,'(A21,i1)') 'sea ice salinity lyr:', k
           units='psu'
           call icepack_files%def_node_var(trim(trname), trim(longname), trim(units), trcrn(:,nt_sice+k-1,:), mesh, partit, ncat)
           write(trname,'(A6,i1)') 'qicen_', k
           write(longname,'(A21,i1)') 'sea ice enthalpy lyr:', k
           units='J/m3'
           call icepack_files%def_node_var(trim(trname), trim(longname), trim(units), trcrn(:,nt_qice+k-1,:), mesh, partit, ncat)
        end do
      
        ! Snow
      
        do k = 1,nslyr
           write(trname,'(A6,i1)') 'qsnon_', k
           write(longname,'(A18,i1)') 'snow enthalpy lyr:', k
           units='J/m3'
           call icepack_files%def_node_var(trim(trname), trim(longname), trim(units), trcrn(:,nt_qsno+k-1,:), mesh, partit, ncat)
        end do
      
        !
        ! All the other 4D tracers (linked to aerosols and biogeochemistry) are at the
        ! moment not supported for restart. This might change if someone is interested
        ! in using the biogeochemistry modules. At this stage, I do not know the model
        ! enough to use these options. Lorenzo Zampieri - 16/10/2019.
        !
      
    end subroutine ini_icepack_io

end submodule icedrv_io
