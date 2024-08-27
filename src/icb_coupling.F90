subroutine reset_ib_fluxes()
    use o_param

! kh 18.03.21 not really used here
    use g_config
    use iceberg_params

    ibfwbv = 0
    ibfwb = 0
    ibfwl = 0
    ibfwe = 0
    ibhf = 0
end subroutine


subroutine prepare_icb2fesom(mesh, partit, ib,i_have_element,localelement,depth_ib)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param
    use MOD_MESH
    use g_config
    use MOD_PARTIT
    use iceberg_params

    logical                 :: i_have_element
    real                    :: depth_ib
    integer                 :: localelement
    integer                 :: iceberg_node  
    integer, dimension(3)   :: ib_nods_in_ib_elem
    integer                 :: num_ib_nods_in_ib_elem
    real                    :: tot_area_nods_in_ib_elem
    integer                 :: i, ib
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    if(i_have_element) then 
        num_ib_nods_in_ib_elem=0
        tot_area_nods_in_ib_elem=0

        do i=1,3
            iceberg_node=elem2d_nodes(i,localelement)

            if (iceberg_node<=mydim_nod2d) then
                ib_nods_in_ib_elem(i)   = iceberg_node
                num_ib_nods_in_ib_elem  = num_ib_nods_in_ib_elem + 1
                tot_area_nods_in_ib_elem= tot_area_nods_in_ib_elem + mesh%area(1,iceberg_node)
            else
                ib_nods_in_ib_elem(i)   = 0
            end if
        end do
        
        do i=1, 3
            iceberg_node=ib_nods_in_ib_elem(i)

            if (iceberg_node>0) then
                ibfwbv(iceberg_node) = ibfwbv(iceberg_node) - fwbv_flux_ib(ib) / tot_area_nods_in_ib_elem
                ibfwb(iceberg_node) = ibfwb(iceberg_node) - fwb_flux_ib(ib) / tot_area_nods_in_ib_elem
                ibfwl(iceberg_node) = ibfwl(iceberg_node) - fwl_flux_ib(ib) / tot_area_nods_in_ib_elem
                ibfwe(iceberg_node) = ibfwe(iceberg_node) - fwe_flux_ib(ib) / tot_area_nods_in_ib_elem
                ibhf(iceberg_node) = ibhf(iceberg_node) - heat_flux_ib(ib) / tot_area_nods_in_ib_elem
            end if
        end do
    end if
end subroutine prepare_icb2fesom

subroutine icb2fesom(mesh, partit, ice)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param

! kh 18.03.21 specification of structure used
    use o_arrays, only: water_flux, heat_flux, wiso_flux_oce
    use MOD_MESH
    use g_config
    use MOD_PARTIT
    use MOD_ICE
    use iceberg_params
    integer                                 :: n
    real(kind=WP), dimension(:)  , pointer :: fresh_wa_flux, net_heat_flux
type(t_ice)   , intent(inout), target :: ice
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
    fresh_wa_flux => ice%flx_fw(:)
    net_heat_flux => ice%flx_h(:)

    do n=1, myDim_nod2d+eDim_nod2D
        if (.not.turn_off_hf) then
                net_heat_flux(n)   = net_heat_flux(n) + ibhf(n) !* steps_per_ib_step
        end if
        if (.not.turn_off_fw) then
                fresh_wa_flux(n)   = fresh_wa_flux(n) + (ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n)) !* steps_per_ib_step
        end if
    end do
!---wiso-code-begin
    if(lwiso) then
      do n=1, myDim_nod2D+eDim_nod2D
      wiso_flux_oce(n,1)=wiso_flux_oce(n,1)+(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))*1000.0*wiso_smow(1)*(1-30.0/1000.0)
      wiso_flux_oce(n,2)=wiso_flux_oce(n,2)+(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))*1000.0*wiso_smow(2)*(1-240.0/1000.0)
      wiso_flux_oce(n,3)=wiso_flux_oce(n,3)+(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))*1000.0
      end do
    end if
!---wiso-code-end
end subroutine icb2fesom
