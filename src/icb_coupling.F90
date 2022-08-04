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


! kh 17.12.21 iceberg indices are at first just buffered in prepare_icb2fesom_part_1 and processed deferred in prepare_icb2fesom_part_2, 
! so that bit identical results are achieved even when the loop over all icebergs is processed in parallel in a "random order"
subroutine prepare_icb2fesom_part_1(mesh, ib,i_have_element,localelement,depth_ib)
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

    logical                 :: i_have_element
    real                    :: depth_ib
    integer                 :: localelement
    integer                 :: iceberg_node  
    integer, dimension(3)   :: ib_nods_in_ib_elem
    integer                 :: num_ib_nods_in_ib_elem
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"

    if(i_have_element) then 

        ! ##############################################################
        ! LA: spread fluxes over all nodes of element, 25.04.2022
        call get_iceberg_nodes_for_element(mesh, localelement, ib_nods_in_ib_elem, num_ib_nods_in_ib_elem)
        !write(*,*) "LA DEBUG: ib_nods_in_ib_elem=",ib_nods_in_ib_elem,", num_ib_nods_in_ib_elem=",num_ib_nods_in_ib_elem

! kh 17.12.21 just buffer the number of nodes
        flux_iceberg_num_nodes_ib(ib) = num_ib_nods_in_ib_elem
        do i=1, 3
            iceberg_node=ib_nods_in_ib_elem(i)

            if (iceberg_node>0) then

! kh 17.12.21 just buffer the index
                flux_iceberg_node_ib(i, ib) = iceberg_node

!               ibfwbv(iceberg_node) = ibfwbv(iceberg_node) - fwbv_flux_ib(ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
!               ibfwb (iceberg_node) = ibfwb (iceberg_node) - fwb_flux_ib (ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
!               ibfwl (iceberg_node) = ibfwl (iceberg_node) - fwl_flux_ib (ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
!               ibfwe (iceberg_node) = ibfwe (iceberg_node) - fwe_flux_ib (ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
!               ibhf  (iceberg_node) = ibhf  (iceberg_node) - heat_flux_ib(ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
            !else
            !   write(*,*) 'iceberg_node only communication node. iceberg_node=',iceberg_node,', mydim_nod2d=',mydim_nod2d
            end if
        end do
    end if
end subroutine prepare_icb2fesom_part_1

subroutine get_iceberg_nodes_for_element(mesh, localelement, ib_nods_in_ib_elem, num_ib_nods_in_ib_elem)
    use o_param
    use o_mesh
    use MOD_MESH
    use g_config
    use g_parsup
    use i_arrays
    use iceberg_params

    integer, intent(in)                     :: localelement
    integer, dimension(3), intent(inout)    :: ib_nods_in_ib_elem
    integer, intent(inout)                  :: num_ib_nods_in_ib_elem
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"

    num_ib_nods_in_ib_elem=0

    do i=1,3
        iceberg_node=elem2d_nodes(i,localelement)

        if (iceberg_node<=mydim_nod2d) then
            ib_nods_in_ib_elem(i)   = iceberg_node
            num_ib_nods_in_ib_elem  = num_ib_nods_in_ib_elem + 1
        else
            ib_nods_in_ib_elem(i)   = 0
        end if
    end do
end subroutine get_iceberg_nodes_for_element

! kh 17.12.21
subroutine prepare_icb2fesom_part_2(mesh)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param
    use o_mesh
    use MOD_MESH
    use g_config
    use g_parsup
    use i_arrays
    use iceberg_params

    integer :: iceberg_node
    integer :: num_ib_nods_in_ib_elem
    integer :: i

type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
    
! kh 17.12.21
    do ib = 1, ib_num
        num_ib_nods_in_ib_elem = flux_iceberg_num_nodes_ib(ib)
        do i = 1, 3
            iceberg_node = flux_iceberg_node_ib(i, ib)
            if(iceberg_node > 0) then
                ibfwbv(iceberg_node) = ibfwbv(iceberg_node) - fwbv_flux_ib(ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
                ibfwb (iceberg_node) = ibfwb (iceberg_node) - fwb_flux_ib (ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
                ibfwl (iceberg_node) = ibfwl (iceberg_node) - fwl_flux_ib (ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
                ibfwe (iceberg_node) = ibfwe (iceberg_node) - fwe_flux_ib (ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
                ibhf  (iceberg_node) = ibhf  (iceberg_node) - heat_flux_ib(ib) / area(1,iceberg_node) / num_ib_nods_in_ib_elem
            end if
        end do
    end do
end subroutine prepare_icb2fesom_part_2


subroutine icb2fesom(mesh)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param

! kh 18.03.21 specification of structure used
    use o_arrays, only: water_flux, heat_flux, wiso_flux_oce

    use o_mesh
    use MOD_MESH
    use g_config
    use g_parsup
    use i_arrays
    use iceberg_params
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"

    do n=1, myDim_nod2d+eDim_nod2D
        water_flux(n)   = water_flux(n) - (ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n)) * steps_per_ib_step
        heat_flux(n)    = heat_flux(n)  - ibhf(n) * steps_per_ib_step
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
