! Inputs and outputs of the sea-ice step. ice_import copies the atmospheric forcing and the
! surface flux fields the ice step uses from g_forcing_arrays / o_arrays into ice%forc,
! ice_export copies the fields the ice step changes back. Between the two the ice step
! (ice_timestep) works only on the ice derived type, its mesh and its partition, which is what
! running it on another decomposition needs. The copies are exact: results do not change.
! Coupled builds (OASIS, IFS interface, YAC) use their own thermodynamics on the global
! arrays, and Icepack its own driver; there only the mass budget diagnostics are copied, or
! nothing.
module ice_import_export_module
    use o_param, only: WP
    use MOD_ICE
    use MOD_PARTIT
    use g_forcing_arrays
    use o_arrays, only: fw_ice, fw_snw
    implicit none
    private
    public :: ice_import, ice_export

contains

    subroutine copy_in(dst, src)
        real(kind=WP), allocatable, intent(inout) :: dst(:)
        real(kind=WP),              intent(in)    :: src(:)
        if (.not. allocated(dst)) allocate(dst(size(src)))
        dst = src
    end subroutine copy_in

    subroutine ice_import(ice, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_partit), intent(in),    target :: partit
#if !defined (__icepack)
        associate(f => ice%forc)
#if !defined (__oasis) && !defined (__ifsinterface) && !defined (__yac)
        call copy_in(f%shortwave,      shortwave)
        call copy_in(f%longwave,       longwave)
        call copy_in(f%tair,           tair)
        call copy_in(f%shum,           shum)
        call copy_in(f%u_wind,         u_wind)
        call copy_in(f%v_wind,         v_wind)
        call copy_in(f%runoff,         runoff)
        call copy_in(f%ch_atm_oce_arr, ch_atm_oce_arr)
        call copy_in(f%ce_atm_oce_arr, ce_atm_oce_arr)
        call copy_in(f%prec_rain,      prec_rain)
        call copy_in(f%prec_snow,      prec_snow)
        call copy_in(f%evaporation,    evaporation)
        call copy_in(f%ice_sublimation, ice_sublimation)
        call copy_in(f%flice,          flice)
        call copy_in(f%fw_ice,         fw_ice)
        call copy_in(f%fw_snw,         fw_snw)
        call copy_in(f%real_salt_flux, real_salt_flux)
        call copy_in(f%hf_qlat,        hf_qlat)
        call copy_in(f%hf_qsen,        hf_qsen)
        call copy_in(f%hf_qradtot,     hf_qradtot)
        call copy_in(f%hf_qswr,        hf_qswr)
        call copy_in(f%hf_qlwr,        hf_qlwr)
        call copy_in(f%hf_qlwrout,     hf_qlwrout)
#endif
        call copy_in(f%sidm_dyn,  sidm_dyn)
        call copy_in(f%sidm_cut,  sidm_cut)
        call copy_in(f%sidm_mpre, sidm_mpre)
        call copy_in(f%sisn_dyn,  sisn_dyn)
        call copy_in(f%sisn_cut,  sisn_cut)
        call copy_in(f%sisn_spre, sisn_spre)
        call copy_in(f%sidc_dyn,  sidc_dyn)
        call copy_in(f%sidc_cut,  sidc_cut)
        call copy_in(f%sidc_apre, sidc_apre)
        end associate
#endif
    end subroutine ice_import

    subroutine ice_export(ice, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_partit), intent(in),    target :: partit
#if !defined (__icepack)
        associate(f => ice%forc)
#if !defined (__oasis) && !defined (__ifsinterface) && !defined (__yac)
        prec_rain       = f%prec_rain
        prec_snow       = f%prec_snow
        evaporation     = f%evaporation
        ice_sublimation = f%ice_sublimation
        flice           = f%flice
        fw_ice          = f%fw_ice
        fw_snw          = f%fw_snw
        real_salt_flux  = f%real_salt_flux
        hf_qlat         = f%hf_qlat
        hf_qsen         = f%hf_qsen
        hf_qradtot      = f%hf_qradtot
        hf_qswr         = f%hf_qswr
        hf_qlwr         = f%hf_qlwr
        hf_qlwrout      = f%hf_qlwrout
#endif
        sidm_dyn  = f%sidm_dyn
        sidm_cut  = f%sidm_cut
        sidm_mpre = f%sidm_mpre
        sisn_dyn  = f%sisn_dyn
        sisn_cut  = f%sisn_cut
        sisn_spre = f%sisn_spre
        sidc_dyn  = f%sidc_dyn
        sidc_cut  = f%sidc_cut
        sidc_apre = f%sidc_apre
        end associate
#endif
    end subroutine ice_export
end module ice_import_export_module
