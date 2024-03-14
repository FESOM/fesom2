! ------------
! 23.03.2023
! OG   
!===============================================================================
! allocate & initialise arrays for REcoM
module recom_init_interface
    interface
        subroutine recom_init(tracers, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in)  ,  target :: mesh
        end subroutine
    end interface
end module
!
!
!_______________________________________________________________________________
subroutine recom_init(tracers, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    use g_clock
    use REcoM_declarations
    use REcoM_GloVar
    use REcoM_locVar
    use recom_config
    use REcoM_ciso
    ! fesom modules
!    use o_ARRAYS
!    use o_MESH
    implicit none
#include "netcdf.inc"
    !___________________________________________________________________________
    ! pointer on necessary derived types 
    integer                                 :: n, k, row, nzmin, nzmax, i, id
    integer                                 :: elem_size, node_size, num_tracers

    real(kind=WP)                           :: locDINmax, locDINmin, locDICmax, locDICmin, locAlkmax, glo
    real(kind=WP)                           :: locAlkmin, locDSimax, locDSimin, locDFemax, locDFemin

    type(t_tracer), intent(inout), target   :: tracers
    type(t_partit), intent(inout), target   :: partit
    type(t_mesh),   intent(in) ,   target   :: mesh

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    elem_size   = myDim_elem2D+eDim_elem2D
    node_size   = myDim_nod2D+eDim_nod2D
    num_tracers = tracers%num_tracers


!! *** Allocate and initialize ***

    allocate(GloFeDust             ( node_size ))
    allocate(AtmFeInput            ( node_size ))
    allocate(GloNDust              ( node_size ))
    allocate(AtmNInput             ( node_size ))

    allocate(cosAI                 ( node_size ))
    allocate(GloPCO2surf           ( node_size ))
    allocate(GloCO2flux            ( node_size ))
    allocate(GloCO2flux_seaicemask ( node_size ))
    allocate(GloO2flux_seaicemask  ( node_size ))
    allocate(GlodPCO2surf          ( node_size ))
    allocate(GlodecayBenthos       ( node_size, benthos_num ))
    allocate(Benthos               ( node_size, benthos_num ))
    allocate(Benthos_tr            ( node_size, benthos_num, num_tracers )) ! kh 25.03.22 buffer per tracer index
    allocate(GloHplus              ( node_size ))
    allocate(DenitBen              ( node_size )) 
    allocate(RiverDIN2D            ( node_size ))
    allocate(RiverDON2D            ( node_size ))
    allocate(RiverDOC2D            ( node_size ))
    allocate(RiverDSi2D            ( node_size ))
    allocate(RiverDIC2D            ( node_size ))
    allocate(RiverAlk2D            ( node_size ))
    allocate(ErosionTON2D          ( node_size ))
    allocate(ErosionTOC2D          ( node_size ))
    allocate(ErosionTSi2D          ( node_size ))
    allocate(relax_alk             ( node_size ))
    allocate(virtual_alk           ( node_size ))
    allocate(LocBenthos            ( benthos_num ))
    allocate(decayBenthos          ( benthos_num ))     ! [1/day] Decay rate of detritus in the benthic layer
    allocate(wFluxPhy              ( benthos_num ))     ! [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of phytoplankton
    allocate(wFluxDia              ( benthos_num ))     ! [mmol/(m2 * day)] Flux of N,C, Si and chl through sinking of diatoms 	
    allocate(PAR3D                 ( nl-1, node_size ))

    GloFeDust             = 0.d0
    AtmFeInput            = 0.d0
    GloNDust              = 0.d0
    AtmNInput             = 0.d0
  
    cosAI                 = 0.d0
    GloPCO2surf           = 0.d0
    GloCO2flux            = 0.d0
    GloCO2flux_seaicemask = 0.d0 
    GloO2flux_seaicemask  = 0.d0 
    GlodPCO2surf          = 0.d0
    GlodecayBenthos       = 0.d0
    Benthos               = 0.d0
    Benthos_tr(:,:,:) = 0.0d0 ! kh 25.03.22
    GloHplus              = exp(-8.d0 * log(10.d0)) ! = 10**(-8)
    DenitBen              = 0.d0

    RiverDIN2D            = 0.d0
    RiverDON2D            = 0.d0
    RiverDOC2D            = 0.d0
    RiverDSi2D            = 0.d0
    RiverDIC2D            = 0.d0
    RiverAlk2D            = 0.d0
    ErosionTON2D          = 0.d0
    ErosionTON2D          = 0.d0
    ErosionTSi2D          = 0.d0
    relax_alk             = 0.d0
    virtual_alk           = 0.d0

    LocBenthos            = 0.d0
    decayBenthos          = 0.d0
    wFluxPhy              = 0.d0
    wFluxDia              = 0.d0
    PAR3D                 = 0.d0
    end subroutine recom_init

