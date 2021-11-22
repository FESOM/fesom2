subroutine reset_ib_fluxes()
    use o_param

! kh 18.03.21 not really used here
!   use o_arrays

    use o_mesh
    use g_config
    use g_parsup
    use i_arrays
    use iceberg_params

    !alles wieder auf null setzen
    ibfwbv = 0
    ibfwb = 0
    ibfwl = 0
    ibfwe = 0
    ibhf = 0
end subroutine


subroutine prepare_icb2fesom(mesh, ib,i_have_element,localelement,depth_ib)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param

! kh 18.03.21 not really used here
!   use o_arrays

    use o_mesh
    use MOD_MESH
    use g_config
    use g_parsup
    use i_arrays
    use iceberg_params

    integer   			:: iceberg_node  
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"

    if(i_have_element) then 
        iceberg_node=elem2D_nodes(1,localelement)

        if (iceberg_node<=mydim_nod2d) then
            ibfwbv(iceberg_node) = ibfwbv(iceberg_node) - fwbv_flux_ib(ib) / area(1,iceberg_node)
            ibfwb(iceberg_node) = ibfwb(iceberg_node) - fwb_flux_ib(ib) / area(1,iceberg_node)
            ibfwl(iceberg_node) = ibfwl(iceberg_node) - fwl_flux_ib(ib) / area(1,iceberg_node)
            ibfwe(iceberg_node) = ibfwe(iceberg_node) - fwe_flux_ib(ib) / area(1,iceberg_node)
            ibhf(iceberg_node) = ibhf(iceberg_node) - heat_flux_ib(ib) / area(1,iceberg_node)
        else
            write(*,*) 'iceberg_node only communication node'
        end if
    end if
end subroutine prepare_icb2fesom


subroutine icb2fesom(mesh)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param

! kh 18.03.21 specification of structure used
    use o_arrays, only: water_flux, heat_flux

    use o_mesh
    use MOD_MESH
    use g_config
    use g_parsup
    use i_arrays
    use iceberg_params
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"

    do n=1, myDim_nod2d
        water_flux(n) = water_flux(n)-ibfwbv(n)-ibfwb(n)-ibfwl(n)-ibfwe(n)
        heat_flux(n) = heat_flux(n) - ibhf(n)
    end do
end subroutine icb2fesom
