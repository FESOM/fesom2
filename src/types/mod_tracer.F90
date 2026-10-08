module mod_tracer
    ! Tracer type (decision D5: EVOLVING). Structure from tracer_dwarf MOD_TRACER.
    ! Per-tracer scheme config rides inside t_tracer_data (D4/D6). DVD-diagnostic
    ! fields dropped (not v1). Tracer DATA = WP (precision under test); tracer WORK
    ! = MP (D3); MP == WP at the DP anchor.
    use mod_precision, only: WP, MP
    use, intrinsic :: iso_fortran_env, only: int32
    use mod_binary_arrays, only: write_bin_array, read_bin_array
    implicit none
    save

    type t_tracer_data
        real(kind=WP), allocatable, dimension(:,:)   :: values    ! (nl-1,nod2D) instant
        real(kind=WP), allocatable, dimension(:,:)   :: valuesAB  ! AB interpolation
        real(kind=WP), allocatable, dimension(:,:,:) :: valuesold ! previous steps
        logical       :: smooth_bh_tra = .false.
        real(kind=WP) :: gamma0_tra = 0.0_WP, gamma1_tra = 0.0_WP, gamma2_tra = 0.0_WP
        logical       :: i_vert_diff = .false.
        character(20) :: tra_adv_hor = 'NONE', tra_adv_ver = 'NONE', tra_adv_lim = 'NONE'
        real(kind=WP) :: tra_adv_ph = 1.0_WP   ! horiz. adv parameter (MUSCL 4th-order fraction)
        real(kind=WP) :: tra_adv_pv = 1.0_WP   ! vert. adv parameter (QR4C 4th-order fraction)
        integer       :: AB_order = 2
        integer       :: ID = 0
    contains
        procedure :: write_td => write_t_tracer_data
        procedure :: read_td  => read_t_tracer_data
    end type t_tracer_data

    type t_tracer_work
        real(kind=MP), allocatable, dimension(:,:) :: del_ttf
        real(kind=MP), allocatable, dimension(:,:) :: del_ttf_advhoriz, del_ttf_advvert
        ! FCT (Zalesak)
        real(kind=MP), allocatable, dimension(:,:) :: fct_LO            ! low-order solution
        real(kind=MP), allocatable, dimension(:,:) :: adv_flux_hor      ! antidiffusive horiz.
        real(kind=MP), allocatable, dimension(:,:) :: adv_flux_ver      ! antidiffusive vert.
        real(kind=MP), allocatable, dimension(:,:) :: fct_ttf_max, fct_ttf_min
        real(kind=MP), allocatable, dimension(:,:) :: fct_plus, fct_minus
        ! MUSCL reconstruction
        integer,       allocatable, dimension(:,:) :: edge_up_dn_tri
        ! Elemental gradient of values over the FULL element halo, rebuilt per tracer by
        ! init_tracers_AB: the MUSCL kernels form the up/downwind reconstruction from it on
        ! the fly (no stored per-edge edge_up_dn_grad; docs/plans/completed/
        ! 2026-10-06-muscl-onthefly.md), the horizontal diffusion / Redi of the same tracer
        ! reuses it. Scratch: NOT serialized.
        real(kind=WP), allocatable, dimension(:,:,:) :: tr_xy             ! (2, nl-1, nElemF)
        ! M4d Redi: vertical tracer gradient (tracer_gradient_z), recomputed per tracer; NOT
        ! serialized (allocated only when Redi). Feeds diff_part_hor_redi's K13/K23 slope terms.
        real(kind=WP), allocatable, dimension(:,:) :: tr_z              ! (nl, nod2D)
    contains
        procedure :: write_tw => write_t_tracer_work
        procedure :: read_tw  => read_t_tracer_work
    end type t_tracer_work

    ! Auxiliary type for reading namelist.tra tracer list.
    type nml_tracer_list_type
        integer       :: ID      = -1
        character(20) :: adv_hor = 'NONE'
        character(20) :: adv_ver = 'NONE'
        character(20) :: adv_lim = 'NONE'
        real(kind=WP) :: adv_ph  = 1.0_WP
        real(kind=WP) :: adv_pv  = 1.0_WP
    end type nml_tracer_list_type

    type t_tracer
        integer :: num_tracers = 2
        type(t_tracer_data), allocatable :: data(:)
        type(t_tracer_work)              :: work
    contains
        procedure :: write_unformatted => write_t_tracer
        procedure :: read_unformatted  => read_t_tracer
        generic   :: write(unformatted) => write_unformatted
        generic   :: read(unformatted)  => read_unformatted
    end type t_tracer

contains

    subroutine write_t_tracer_data(td, unit)
        class(t_tracer_data), intent(in) :: td
        integer,              intent(in) :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        call write_bin_array(td%values,    unit, iostat, iomsg)
        call write_bin_array(td%valuesold, unit, iostat, iomsg)
        call write_bin_array(td%valuesAB,  unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) td%smooth_bh_tra
        write(unit, iostat=iostat, iomsg=iomsg) td%gamma0_tra, td%gamma1_tra, td%gamma2_tra
        write(unit, iostat=iostat, iomsg=iomsg) td%i_vert_diff
        write(unit, iostat=iostat, iomsg=iomsg) td%tra_adv_hor, td%tra_adv_ver, td%tra_adv_lim
        write(unit, iostat=iostat, iomsg=iomsg) td%tra_adv_ph, td%tra_adv_pv
        write(unit, iostat=iostat, iomsg=iomsg) td%AB_order, td%ID
    end subroutine write_t_tracer_data

    subroutine read_t_tracer_data(td, unit)
        class(t_tracer_data), intent(inout) :: td
        integer,              intent(in)    :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        call read_bin_array(td%values,    unit, iostat, iomsg)
        call read_bin_array(td%valuesold, unit, iostat, iomsg)
        call read_bin_array(td%valuesAB,  unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) td%smooth_bh_tra
        read(unit, iostat=iostat, iomsg=iomsg) td%gamma0_tra, td%gamma1_tra, td%gamma2_tra
        read(unit, iostat=iostat, iomsg=iomsg) td%i_vert_diff
        read(unit, iostat=iostat, iomsg=iomsg) td%tra_adv_hor, td%tra_adv_ver, td%tra_adv_lim
        read(unit, iostat=iostat, iomsg=iomsg) td%tra_adv_ph, td%tra_adv_pv
        read(unit, iostat=iostat, iomsg=iomsg) td%AB_order, td%ID
    end subroutine read_t_tracer_data

    subroutine write_t_tracer_work(tw, unit)
        class(t_tracer_work), intent(in) :: tw
        integer,              intent(in) :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        call write_bin_array(tw%del_ttf,          unit, iostat, iomsg)
        call write_bin_array(tw%del_ttf_advhoriz, unit, iostat, iomsg)
        call write_bin_array(tw%del_ttf_advvert,  unit, iostat, iomsg)
        call write_bin_array(tw%fct_LO,           unit, iostat, iomsg)
        call write_bin_array(tw%adv_flux_hor,     unit, iostat, iomsg)
        call write_bin_array(tw%adv_flux_ver,     unit, iostat, iomsg)
        call write_bin_array(tw%fct_ttf_max,      unit, iostat, iomsg)
        call write_bin_array(tw%fct_ttf_min,      unit, iostat, iomsg)
        call write_bin_array(tw%fct_plus,         unit, iostat, iomsg)
        call write_bin_array(tw%fct_minus,        unit, iostat, iomsg)
        call write_bin_array(tw%edge_up_dn_tri,   unit, iostat, iomsg)
    end subroutine write_t_tracer_work

    subroutine read_t_tracer_work(tw, unit)
        class(t_tracer_work), intent(inout) :: tw
        integer,              intent(in)    :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        call read_bin_array(tw%del_ttf,          unit, iostat, iomsg)
        call read_bin_array(tw%del_ttf_advhoriz, unit, iostat, iomsg)
        call read_bin_array(tw%del_ttf_advvert,  unit, iostat, iomsg)
        call read_bin_array(tw%fct_LO,           unit, iostat, iomsg)
        call read_bin_array(tw%adv_flux_hor,     unit, iostat, iomsg)
        call read_bin_array(tw%adv_flux_ver,     unit, iostat, iomsg)
        call read_bin_array(tw%fct_ttf_max,      unit, iostat, iomsg)
        call read_bin_array(tw%fct_ttf_min,      unit, iostat, iomsg)
        call read_bin_array(tw%fct_plus,         unit, iostat, iomsg)
        call read_bin_array(tw%fct_minus,        unit, iostat, iomsg)
        call read_bin_array(tw%edge_up_dn_tri,   unit, iostat, iomsg)
    end subroutine read_t_tracer_work

    subroutine write_t_tracer(tracer, unit, iostat, iomsg)
        class(t_tracer), intent(in)    :: tracer
        integer,         intent(in)    :: unit
        integer,         intent(out)   :: iostat
        character(*),    intent(inout) :: iomsg
        integer :: i
        write(unit, iostat=iostat, iomsg=iomsg) tracer%num_tracers
        do i = 1, tracer%num_tracers
            call tracer%data(i)%write_td(unit)
        end do
        call tracer%work%write_tw(unit)
    end subroutine write_t_tracer

    subroutine read_t_tracer(tracer, unit, iostat, iomsg)
        class(t_tracer), intent(inout) :: tracer
        integer,         intent(in)    :: unit
        integer,         intent(out)   :: iostat
        character(*),    intent(inout) :: iomsg
        integer :: i
        read(unit, iostat=iostat, iomsg=iomsg) tracer%num_tracers
        if (.not. allocated(tracer%data)) allocate(tracer%data(tracer%num_tracers))
        do i = 1, tracer%num_tracers
            call tracer%data(i)%read_td(unit)
        end do
        call tracer%work%read_tw(unit)
    end subroutine read_t_tracer

end module mod_tracer
