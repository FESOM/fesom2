subroutine reset_ib_fluxes()
    use o_param

! kh 18.03.21 not really used here
!   use o_arrays

!    use o_mesh
    use g_config
!    use MOD_PARTIT
!    use i_arrays
    use iceberg_params

    !alles wieder auf null setzen
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

! kh 18.03.21 not really used here
!   use o_arrays

!    use o_mesh
    use MOD_MESH
    use g_config
    use MOD_PARTIT
!    use i_arrays
    use iceberg_params

    logical                 :: i_have_element
    real                    :: depth_ib
    integer                 :: localelement
    integer                 :: iceberg_node  
    integer, dimension(3)   :: ib_nods_in_ib_elem
    integer                 :: num_ib_nods_in_ib_elem
    integer                 :: i, ib
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    if(i_have_element) then 
        ! ##############################################################
        ! LA: spread fluxes over all nodes of element, 25.04.2022
        ! LA: Das geht nicht auf, wenn die areas der verschiedenen nodes unterscheidlich groß sind:
        !                       1/3 * flux1 * area1 + 1/3 * flux2 * area2 + 1/3 * flux3 * area3 != flux_tot
        ! richtig waere:        area1 / area_tot * flux1 + area2 / area_tot * flux2 + area3 / area_tot * flux3 = flux_tot
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
        ! ##############################################################
        
        !write(*,*) "LA DEBUG: ib_nods_in_ib_elem=",ib_nods_in_ib_elem,", num_ib_nods_in_ib_elem=",num_ib_nods_in_ib_elem 
        do i=1, 3
            iceberg_node=ib_nods_in_ib_elem(i)

            if (iceberg_node>0) then
                ibfwbv(iceberg_node) = ibfwbv(iceberg_node) - fwbv_flux_ib(ib) / mesh%area(1,iceberg_node) / num_ib_nods_in_ib_elem
                ibfwb(iceberg_node) = ibfwb(iceberg_node) - fwb_flux_ib(ib) / mesh%area(1,iceberg_node) / num_ib_nods_in_ib_elem
                ibfwl(iceberg_node) = ibfwl(iceberg_node) - fwl_flux_ib(ib) / mesh%area(1,iceberg_node) / num_ib_nods_in_ib_elem
                ibfwe(iceberg_node) = ibfwe(iceberg_node) - fwe_flux_ib(ib) / mesh%area(1,iceberg_node) / num_ib_nods_in_ib_elem
                ibhf(iceberg_node) = ibhf(iceberg_node) - heat_flux_ib(ib) / mesh%area(1,iceberg_node) / num_ib_nods_in_ib_elem
                write(*,*) "LA DEBUG: ibhf(iceberg_node)=",ibhf(iceberg_node) 
            !else                                                             
            !    write(*,*) 'iceberg_node only communication node. iceberg_node=',iceberg_node,', mydim_nod2d=',mydim_nod2d
            end if
        end do
    end if
end subroutine prepare_icb2fesom

!subroutine get_iceberg_nodes_for_element(elem2d_nodes, localelement, ib_nods_in_ib_elem, num_ib_nods_in_ib_elem)
!    use o_param
!    use o_mesh
!    use MOD_MESH
!    use g_config
!    use MOD_PARTIT
!    use i_arrays
!    use iceberg_params
!
!    integer, intent(in)                     :: localelement
!    integer, dimension(3), intent(inout)    :: ib_nods_in_ib_elem
!    integer, intent(inout)                  :: num_ib_nods_in_ib_elem
!    integer                                 :: i, iceberg_node
!    integer, dimension(:,:), intent(in)     :: elem2d_nodes
!!    integer, intent(in)                     :: mydim_nod2d
!!type(t_mesh), intent(in)                :: mesh
!!#include "associate_mesh.h"
!
!    num_ib_nods_in_ib_elem=0
!
!    do i=1,3
!        iceberg_node=elem2d_nodes(i,localelement)
!
!        if (iceberg_node<=mydim_nod2d) then
!            ib_nods_in_ib_elem(i)   = iceberg_node
!            num_ib_nods_in_ib_elem  = num_ib_nods_in_ib_elem + 1
!        else
!            ib_nods_in_ib_elem(i)   = 0
!        end if
!    end do
!end subroutine get_iceberg_nodes_for_element


subroutine icb2fesom(mesh, partit, ice)
    !transmits the relevant fields from the iceberg to the ocean model
    !Lars Ackermann, 17.03.2020

    use o_param

! kh 18.03.21 specification of structure used
    use o_arrays, only: water_flux, heat_flux !, wiso_flux_oce

!    use o_mesh
    use MOD_MESH
    use g_config
    use MOD_PARTIT
    use MOD_ICE
    !use i_arrays
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
        !if (ibhf(n).ne.0) then
        !        write(*,*) "LA DEBUG: *** n = ",n
        !        write(*,*) "LA DEBUG: fresh_wa_flux(n) = ", fresh_wa_flux(n), "(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n)) = ", (ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))
        !        write(*,*) "LA DEBUG: net_heat_flux(n) = ", net_heat_flux(n), "ibhf(n) = ", (ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))
        !end if
        fresh_wa_flux(n)   = fresh_wa_flux(n) + (ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n)) * steps_per_ib_step !sign needs to be negative
        net_heat_flux(n)   = net_heat_flux(n) + ibhf(n) * steps_per_ib_step !LA: reversed sign!!! needs to be negative
    end do
!!---wiso-code-begin
!   if(lwiso) then
!     do n=1, myDim_nod2D+eDim_nod2D
!     wiso_flux_oce(n,1)=wiso_flux_oce(n,1)+(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))*1000.0*wiso_smow(1)*(1-30.0/1000.0)
!     wiso_flux_oce(n,2)=wiso_flux_oce(n,2)+(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))*1000.0*wiso_smow(2)*(1-240.0/1000.0)
!     wiso_flux_oce(n,3)=wiso_flux_oce(n,3)+(ibfwb(n)+ibfwl(n)+ibfwe(n)+ibfwbv(n))*1000.0
!     end do
!   end if
!!---wiso-code-end

end subroutine icb2fesom
