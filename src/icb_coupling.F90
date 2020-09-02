!===========================================
subroutine icb2ice
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param
    use o_arrays
    use o_mesh
    use g_parsup
    use g_config
    use i_arrays
    use iceberg_params

    do ib=1, ib_num
        ! ... iterate over all icebergs and get freshwater and heat fluxes
        do m=1, 3
            node2_id = elem2D_nodes(m, ib)
        
            !fresh_wa_flux(node2_id) = fresh_wa_flux(node2_id) - fw_flux_ib(ib)/3
            !net_heat_flux(node2_id) = net_heat_flux(node2d_id) - heat_flux_ib(ib)/3
        enddo
    enddo

end subroutine icb2ice










subroutine icb2fesom_grid(ib,i_have_element,localelement)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param
    use o_arrays
    use o_mesh
    use g_config
    use i_arrays
    use iceberg_params

    integer, dimension(:), save, allocatable :: local_idx_of
    integer   			:: iceberg_node  
  
      if(i_have_element) then 
       !write(*,*) '*** LA debugging ***'
       !write(*,*) 'IB ',ib,' in element ',iceberg_elem
       !write(*,*) 'IB ',ib,' rot. coords:', lon_deg, lat_deg !,lon_rad, lat_rad
       do i=1,3
        iceberg_node=elem2D_nodes(i,localelement)
        !write(*,*) 'node ', i, 's coords lon:', geo_coord_nod2d(1,iceberg_node)/rad,' lat: ', geo_coord_nod2d(2,iceberg_node)/rad
        !write(*,*) 'fwb_flux_ib(ib)=',fwb_flux_ib(ib)
        !write(*,*) 'fwl_flux_ib(ib)=',fwl_flux_ib(ib)
        !write(*,*) 'fwe_flux_ib(ib)=',fwe_flux_ib(ib)
        
        ibfwb(iceberg_node) = fwb_flux_ib(ib)/3
        ibfwl(iceberg_node) = fwl_flux_ib(ib)/3
        ibfwe(iceberg_node) = fwe_flux_ib(ib)/3
        ibhf(iceberg_node) = heat_flux_ib(ib)/3
        water_flux(iceberg_node) = water_flux(iceberg_node)-(fwb_flux_ib(ib)+fwl_flux_ib(ib)+fwe_flux_ib(ib))/3
        !heat_flux(iceberg_node) = heat_flux(iceberg_node)-heat_flux_ib(iceberg_node)/3
       end do
      end if
end subroutine icb2fesom_grid
