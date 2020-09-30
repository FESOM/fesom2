subroutine icb2fesom(ib,i_have_element,localelement,depth_ib)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param
    use o_arrays
    use o_mesh
    use g_config
    use g_parsup
    use i_arrays
    use iceberg_params

    integer   			:: iceberg_node  
    real                        :: scaling_factor

    if(i_have_element) then 
        iceberg_node=elem2D_nodes(1,localelement)
        scaling_factor = 1/area(1,iceberg_node)

        if (iceberg_node<=mydim_nod2d) then
            ibfwbv(iceberg_node) = ibfwbv(iceberg_node) - fwbv_flux_ib(ib) * scaling_factor
            ibfwb(iceberg_node) = ibfwb(iceberg_node) - fwb_flux_ib(ib) * scaling_factor
            ibfwl(iceberg_node) = ibfwl(iceberg_node) - fwl_flux_ib(ib) * scaling_factor
            ibfwe(iceberg_node) = ibfwe(iceberg_node) - fwe_flux_ib(ib) * scaling_factor
            ibhf(iceberg_node) = ibhf(iceberg_node) - heat_flux_ib(ib) * scaling_factor
            water_flux(iceberg_node) = water_flux(iceberg_node)-ibfwbv(iceberg_node)-ibfwb(iceberg_node)-ibfwl(iceberg_node)-ibfwe(iceberg_node)
        else
            write(*,*) 'iceberg_node only communication node'
        end if
    end if
end subroutine icb2fesom
