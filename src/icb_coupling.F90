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
       !write(*,*) '*** LA debugging ***'
       !write(*,*) 'IB ',ib,' in element ',iceberg_elem
       !write(*,*) 'IB ',ib,' rot. coords:', lon_deg, lat_deg !,lon_rad, lat_rad
        !do i=1,3
        iceberg_node=elem2D_nodes(1,localelement)
        scaling_factor = dt/REAL(steps_per_FESOM_step)/area(1,iceberg_node)
        !scaling_factor = 1/area(1,iceberg_node)
        !write(*,*) 'LA DEBUG: localelment=',localelement
        !write(*,*) 'LA DEBUG: iceberg_node=',iceberg_node
        !write(*,*) 'LA DEBUG: mydim_nod2d=',mydim_nod2d
        !write(*,*) 'node ', i, 's coords lon:', geo_coord_nod2d(1,iceberg_node)/rad,' lat: ', geo_coord_nod2d(2,iceberg_node)/rad
        !write(*,*) 'fwb_flux_ib(ib)=',fwb_flux_ib(ib)
        !write(*,*) 'fwl_flux_ib(ib)=',fwl_flux_ib(ib)
        !write(*,*) 'fwe_flux_ib(ib)=',fwe_flux_ib(ib)
        !
        !write(*,*) 'LA DEBUG: dt=',dt,', steps_per_FESOM_step=',steps_per_FESOM_step
        !write(*,*) 'LA DEBUG: length(ib)=',length_ib(ib)
        !write(*,*) 'LA DEBUG: height(ib)=',height_ib(ib)
        !write(*,*) 'LA DEBUG: width(ib)=',width_ib(ib)
        !write(*,*) 'LA DEBUG: depth(ib)=',depth_ib
        !write(*,*) 'LA DEBUG: area(1,iceberg_node)=',area(1,iceberg_node)
       
        !dh_b = fwb_flux_ib(ib)*dt/REAL(steps_per_FESOM_step)    !change of height..
        !dh_v = fwl_flux_ib(ib)*dt/REAL(steps_per_FESOM_step)    !..and length due to melting..
        !dh_e = fwe_flux_ib(ib)*dt/REAL(steps_per_FESOM_step)    !..and due to wave erosion [m].
        !dh_bv = fwbv_flux_ib(ib)*dt/REAL(steps_per_FESOM_step)  !change of length due to 'basal meltrate'
        !CALCULATION OF WORKING SURFACES AS IN BIGG (1997) & SILVA (2010)
        !basal volume loss
        bvl = fwb_flux_ib(ib)*length_ib(ib)**2
        !lateral volume loss due to ...
        !erosion
        lvl_e = fwe_flux_ib(ib)*length_ib(ib)*height_ib(ib) + fwe_flux_ib(ib)*width_ib(ib)*height_ib(ib)
        !basal melting
        lvl_b = fwbv_flux_ib(ib)*2*length_ib(ib)*abs(depth_ib) + fwbv_flux_ib(ib)*2*width_ib(ib)*abs(depth_ib)
        !lateral melting
        lvl_v = fwl_flux_ib(ib)*2*length_ib(ib)*abs(depth_ib) + fwl_flux_ib(ib)*2*width_ib(ib)*abs(depth_ib)
        !total volume loss
        tvl = bvl + lvl_b + lvl_v + lvl_e       ![m^3] per timestep, for freshwater flux convert somehow to [m/s]
        !write(*,*) 'LA DEBUG 1: tvl=',tvl
   
        !tvl = tvl*dt/REAL(steps_per_FESOM_step)
        !write(*,*) 'LA DEBUG 2: tvl=',tvl

        if (iceberg_node<=mydim_nod2d) then
            ibfwb(iceberg_node) = ibfwb(iceberg_node) - fwb_flux_ib(ib) * scaling_factor
            ibfwl(iceberg_node) = ibfwl(iceberg_node) - fwl_flux_ib(ib) * scaling_factor
            ibfwe(iceberg_node) = ibfwe(iceberg_node) - fwe_flux_ib(ib) * scaling_factor
            ibhf(iceberg_node) = ibhf(iceberg_node) - heat_flux_ib(ib) * scaling_factor
            water_flux(iceberg_node) = water_flux(iceberg_node) - tvl * scaling_factor
        else
            write(*,*) 'iceberg_node only communication node'
        end if
        !end do
    end if
end subroutine icb2fesom
