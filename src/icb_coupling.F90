!===========================================
subroutine icb2ice
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param
    use o_arrays
    use o_mesh
    use o_mesh
    use g_parsup
    use g_config
    use i_arrays
    use iceberg_params

    do ib=1, ib_num
        ! ... iterate over all icebergs and get freshwater and heat fluxes
        do m=1, 3
            node2_id = elem2D_nodes(m, ib)
        
            fresh_wa_flux(node2_id) = fresh_wa_flux(node2_id) - fw_flux_ib(ib)/3
            net_heat_flux(node2_id) = net_heat_flux(node2d_id) - heat_flux_ib(ib)/3
        enddo
    enddo

end subroutine icb2ice
