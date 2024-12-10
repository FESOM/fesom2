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
    implicit none
#include "netcdf.inc"
    !___________________________________________________________________________
    ! pointer on necessary derived types 
    integer                                 :: n, k, row, nzmin, nzmax, i, id
    integer                                 :: elem_size, node_size, num_tracers

    real(kind=WP)                           :: locDINmax, locDINmin, locDICmax, locDICmin, locAlkmax, glo
    real(kind=WP)                           :: locAlkmin, locDSimax, locDSimin, locDFemax, locDFemin
    real(kind=WP)                           :: locO2max, locO2min

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

    !! * Fe and N deposition as surface boundary condition *
    allocate(GloFeDust             ( node_size ))
    allocate(AtmFeInput            ( node_size ))
    allocate(GloNDust              ( node_size ))
    allocate(AtmNInput             ( node_size ))

    !! * River nutrients as surface boundary condition *
    allocate(RiverDIN2D            ( node_size ))
    allocate(RiverDON2D            ( node_size ))
    allocate(RiverDOC2D            ( node_size ))
    allocate(RiverDSi2D            ( node_size ))
    allocate(RiverDIC2D            ( node_size ))
    allocate(RiverAlk2D            ( node_size ))

    !! * Erosion nutrients as surface boundary condition *
    allocate(ErosionTON2D          ( node_size ))
    allocate(ErosionTOC2D          ( node_size ))
    allocate(ErosionTSi2D          ( node_size ))

    !! * Alkalinity restoring to climatology *
    allocate(relax_alk             ( node_size ))
    allocate(virtual_alk           ( node_size ))

    allocate(cosAI                 ( node_size ))
    allocate(GloPCO2surf           ( node_size ))
    allocate(GloCO2flux            ( node_size ))
    allocate(GloO2flux             ( node_size ))
    allocate(GloCO2flux_seaicemask ( node_size ))
    allocate(GloO2flux_seaicemask  ( node_size ))
    allocate(GlodPCO2surf          ( node_size ))
    allocate(GlodecayBenthos       ( node_size, benthos_num ))
    allocate(Benthos               ( node_size, benthos_num ))
    allocate(Benthos_tr            ( node_size, benthos_num, num_tracers )) ! kh 25.03.22 buffer per tracer index
    allocate(GloHplus              ( node_size ))
    allocate(DenitBen              ( node_size ))

    allocate(LocBenthos            ( benthos_num ))
    allocate(decayBenthos          ( benthos_num ))     ! [1/day] Decay rate of detritus in the benthic layer
    allocate(PAR3D                 ( nl-1, node_size ))

    GloFeDust             = 0.d0
    AtmFeInput            = 0.d0
    GloNDust              = 0.d0
    AtmNInput             = 0.d0

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

    cosAI                 = 0.d0
    GloPCO2surf           = 0.d0
    GloCO2flux            = 0.d0
    GloCO2flux_seaicemask = 0.d0
    GloO2flux_seaicemask  = 0.d0
    GlodPCO2surf          = 0.d0
    GlodecayBenthos       = 0.d0
    Benthos               = 0.d0
    Benthos_tr(:,:,:)     = 0.0d0 ! kh 25.03.22
    GloHplus              = exp(-8.d0 * log(10.d0)) ! = 10**(-8)
    DenitBen              = 0.d0

    LocBenthos            = 0.d0
    decayBenthos          = 0.d0
    wFluxPhy              = 0.d0
    wFluxDia              = 0.d0
    PAR3D                 = 0.d0

!    pco2surf           = 0.d0
!    dflux              = 0.d0
!    oflux              = 0.d0
!    co2flux_seaicemask = 0.d0
!    o2flux_seaicemask  = 0.d0
!    dpco2surf          = 0.d0
!    co2                = 0.d0

    if (Diags) then

!! *** Allocate 2D diagnostics ***
    allocate(NPPn    ( node_size ))
    allocate(NPPd    ( node_size ))
    allocate(NPPc    ( node_size ))
    allocate(GPPn    ( node_size ))
    allocate(GPPd    ( node_size ))
    allocate(GPPc    ( node_size ))
    allocate(NNAn    ( node_size ))
    allocate(NNAd    ( node_size ))
    allocate(NNAc    ( node_size ))
    allocate(Chldegn ( node_size ))
    allocate(Chldegd ( node_size ))
    allocate(Chldegc ( node_size ))

    NPPn    = 0.d0
    NPPd    = 0.d0
    NPPc    = 0.d0
    GPPn    = 0.d0
    GPPd    = 0.d0
    GPPc    = 0.d0
    NNAn    = 0.d0
    NNAd    = 0.d0
    NNAc    = 0.d0
    Chldegn = 0.d0
    Chldegd = 0.d0
    Chldegc = 0.d0

!! *** Allocate 3D diagnostics ***
    allocate(respmeso     ( nl-1, node_size ))
    allocate(respmacro    ( nl-1, node_size ))
    allocate(respmicro    ( nl-1, node_size ))
    allocate(calcdiss     ( nl-1, node_size ))
    allocate(calcif       ( nl-1, node_size ))
    allocate(aggn         ( nl-1, node_size ))
    allocate(aggd         ( nl-1, node_size ))
    allocate(aggc         ( nl-1, node_size ))
    allocate(docexn       ( nl-1, node_size ))
    allocate(docexd       ( nl-1, node_size ))
    allocate(docexc       ( nl-1, node_size ))
    allocate(respn        ( nl-1, node_size ))
    allocate(respd        ( nl-1, node_size ))
    allocate(respc        ( nl-1, node_size ))
    allocate(NPPn3D       ( nl-1, node_size ))
    allocate(NPPd3D       ( nl-1, node_size ))
    allocate(NPPc3D       ( nl-1, node_size ))

    respmeso     = 0.d0
    respmacro    = 0.d0
    respmicro    = 0.d0
    calcdiss     = 0.d0
    calcif       = 0.d0
    aggn         = 0.d0
    aggd         = 0.d0
    aggc         = 0.d0
    docexn       = 0.d0
    docexd       = 0.d0
    docexc       = 0.d0
    respn        = 0.d0
    respd        = 0.d0
    respc        = 0.d0
    NPPn3D       = 0.d0
    NPPd3D       = 0.d0
    NPPc3D       = 0.d0
    end if

!! *** Allocate 3D mocsy ***
    allocate(CO23D        ( nl-1, node_size ))
    allocate(pH3D         ( nl-1, node_size ))
    allocate(pCO23D       ( nl-1, node_size ))
    allocate(HCO33D       ( nl-1, node_size ))
    allocate(CO33D        ( nl-1, node_size ))
    allocate(OmegaC3D     ( nl-1, node_size ))
    allocate(kspc3D       ( nl-1, node_size ))
    allocate(rhoSW3D      ( nl-1, node_size ))
  
    CO23D(:,:)          = 0.d0
    pH3D(:,:)           = 0.d0
    pCO23D(:,:)         = 0.d0
    HCO33D(:,:)         = 0.d0
    CO33D(:,:)          = 0.d0
    OmegaC3D(:,:)       = 0.d0
    kspc3D(:,:)         = 0.d0
    rhoSW3D(:,:)        = 0.d0

!! *** Allocate ballasting ***
    allocate(rho_particle1       ( nl-1, node_size ))
    allocate(rho_particle2       ( nl-1, node_size ))
    allocate(scaling_density1_3D ( nl,   node_size ))
    allocate(scaling_density2_3D ( nl,   node_size ))
    allocate(scaling_visc_3D     ( nl,   node_size ))
    allocate(seawater_visc_3D    ( nl-1, node_size ))
    rho_particle1       = 0.d0
    rho_particle2       = 0.d0
    scaling_density1_3D = 0.d0
    scaling_density2_3D = 0.d0
    scaling_visc_3D     = 0.d0
    seawater_visc_3D    = 0.d0

    allocate(Sinkingvel1(nl,node_size), Sinkingvel2(nl,node_size))
    Sinkingvel1(:,:)      = 0.d0
    Sinkingvel2(:,:)      = 0.d0

    if (use_MEDUSA) then
        allocate(GloSed(node_size,sedflx_num))
        allocate(SinkFlx(node_size,bottflx_num))
        allocate(SinkFlx_tr(node_size,bottflx_num,num_tracers)) ! kh 25.03.22 buffer sums per tracer index

        SinkFlx(:,:)      = 0.d0
        SinkFlx_tr(:,:,:) = 0.0d0 ! kh 25.03.22
        GloSed(:,:)       = 0.d0
        allocate(lb_flux(node_size,9))
        lb_flux(:,:)      = 0.d0
    end if

    DO i=num_tracers-bgc_num+1, num_tracers
        id=tracers%data(i)%ID

        SELECT CASE (id)

! *******************
! CASE 2phy 1zoo 1det
! *******************
! Skip: DIN, DIC, Alk, DSi and O2 are read from files
! Fe [mol/L] => [umol/m3] Check the units again!

        CASE (1004)
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max       ! PhyN

        CASE (1005)
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max/NCmax ! PhyC

        CASE (1006)
            tracers%data(i)%values(:,:) = tiny_chl                 ! PhyChl

        CASE (1007)
            tracers%data(i)%values(:,:) = tiny                     ! DetN

        CASE (1008)
            tracers%data(i)%values(:,:) = tiny                     ! DetC

        CASE (1009)
            tracers%data(i)%values(:,:) = tiny                     ! HetN

        CASE (1010)
            tracers%data(i)%values(:,:) = tiny * Redfield          ! HetC

        CASE (1011)
            tracers%data(i)%values(:,:) = tiny                     ! DON

        CASE (1012)
            tracers%data(i)%values(:,:) = tiny                     ! DOC

        CASE (1013)
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max       ! DiaN

        CASE (1014)
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max/NCmax ! DiaC

        CASE (1015)
            tracers%data(i)%values(:,:) = tiny_chl                 ! DiaChl

        CASE (1016)
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max_d/NCmax_d/SiCmax ! DiaSi

        CASE (1017)
            tracers%data(i)%values(:,:) = tiny                     ! DetSi

        CASE (1019)
            tracers%data(i)%values(:,:) = tracers%data(i)%values(:,:)* 1.e9 ! Fe [mol/L] => [umol/m3] Check the units again!

        CASE (1020)
            tracers%data(i)%values(:,:) = tiny * Redfield          ! PhyCalc

        CASE (1021)
            tracers%data(i)%values(:,:) = tiny                     ! DetCalc

! *******************
! CASE 2phy 2zoo 2det
! *******************
#if defined (__3Zoo2Det)
        CASE (1023)
            tracers%data(i)%values(:,:) = tiny                     ! Zoo2N
        CASE (1024)
            tracers%data(i)%values(:,:) = tiny * Redfield          ! Zoo2C
        CASE (1025)
            tracers%data(i)%values(:,:) = tiny                     ! DetZ2N 
        CASE (1026)
            tracers%data(i)%values(:,:) = tiny                     ! DetZ2C
        CASE (1027)
            tracers%data(i)%values(:,:) = tiny                     ! DetZ2Si
        CASE (1028)
            tracers%data(i)%values(:,:) = tiny                     ! DetZ2Calc
#endif

! *******************
! CASE 3phy 2zoo 2det
! *******************
#if defined (__coccos) & defined (__3Zoo2Det)
        CASE (1029)
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max       ! CoccoN

        CASE (1030)
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max/NCmax ! CoccoC

        CASE (1031)
            tracers%data(i)%values(:,:) = tiny_chl                 ! CoccoChl

! *******************
! CASE 3phy 1zoo 1det
! *******************
#elif defined (__coccos) & !defined (__3Zoo2Det)
        CASE (1023)
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max       ! CoccoN

        CASE (1024)
            tracers%data(i)%values(:,:) = tiny_chl/chl2N_max/NCmax ! CoccoC

        CASE (1025)
            tracers%data(i)%values(:,:) = tiny_chl                 ! CoccoChl

#endif

! *******************
! CASE 3phy 3zoo 2det
! *******************
#if defined (__coccos) & defined (__3Zoo2Det)
        CASE (1032)
            tracers%data(i)%values(:,:) = tiny                     ! Zoo3N

        CASE (1033)
            tracers%data(i)%values(:,:) = tiny * Redfield          ! Zoo3C

#elif !defined (__coccos) & defined (__3Zoo2Det)
! *******************
! CASE 2phy 3zoo 2det
! *******************
        CASE (1029)
            tracers%data(i)%values(:,:) = tiny                     ! Zoo3N

        CASE (1030)
            tracers%data(i)%values(:,:) = tiny * Redfield          ! Zoo3C

#endif

        END SELECT
    END DO
!------------------------------------------

    !< Mask hydrothermal vent in Eastern Equatorial Pacific GO
    do row=1, myDim_nod2D+eDim_nod2D
        !if (ulevels_nod2D(row)>1) cycle 
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)-1
        do k=nzmin, nzmax
            ! do not take regions shallower than 2000 m into account
            if (((geo_coord_nod2D(2,row) > -12.5*rad) .and. (geo_coord_nod2D(2,row) < 9.5*rad))&
                .and.((geo_coord_nod2D(1,row)> -106.0*rad) .and. (geo_coord_nod2D(1,row) < -72.0*rad))) then
                if (abs(Z_3d_n(k,row))<2000.0_WP) cycle
                tracers%data(21)%values(k,row) = min(0.3, tracers%data(21)%values(k,row)) ! OG todo: try 0.6 
            end if
        end do
    end do

    !< Mask negative values
    tracers%data(21)%values(:,:) = max(tiny, tracers%data(21)%values(:,:))
!------------------------------------------

    if(mype==0) write(*,*),'Tracers have been initialized as spinup from WOA/glodap netcdf files'
        locDINmax = -66666
        locDINmin = 66666
        locDICmax = locDINmax
        locDICmin = locDINmin
        locAlkmax = locDINmax
        locAlkmin = locDINmin
        locDSimax = locDINmax
        locDSimin = locDINmin
        locDFemax = locDINmax
        locDFemin = locDINmin
        locO2max  = locDINmax
        locO2min  = locDINmin

        do n=1, myDim_nod2d
            locDINmax = max(locDINmax,maxval(tracers%data(3)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDINmin = min(locDINmin,minval(tracers%data(3)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDICmax = max(locDICmax,maxval(tracers%data(4)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDICmin = min(locDICmin,minval(tracers%data(4)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locAlkmax = max(locAlkmax,maxval(tracers%data(5)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locAlkmin = min(locAlkmin,minval(tracers%data(5)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDSimax = max(locDSimax,maxval(tracers%data(20)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDSimin = min(locDSimin,minval(tracers%data(20)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDFemax = max(locDFemax,maxval(tracers%data(21)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDFemin = min(locDFemin,minval(tracers%data(21)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locO2max  = max(locO2max,maxval(tracers%data(24)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locO2min  = min(locO2min,minval(tracers%data(24)%values(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
        end do

        if (mype==0) write(*,*) "Sanity check for REcoM variables after recom_init call"
        call MPI_AllREDUCE(locDINmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DIN. =', glo
        call MPI_AllREDUCE(locDINmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DIN. =', glo

        call MPI_AllREDUCE(locDICmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DIC. =', glo
        call MPI_AllREDUCE(locDICmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DIC. =', glo
        call MPI_AllREDUCE(locAlkmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. Alk. =', glo
        call MPI_AllREDUCE(locAlkmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. Alk. =', glo
        call MPI_AllREDUCE(locDSimax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DSi. =', glo
        call MPI_AllREDUCE(locDSimin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DSi. =', glo
        call MPI_AllREDUCE(locDFemax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DFe. =', glo
        call MPI_AllREDUCE(locDFemin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  `-> gobal min init. DFe. =', glo
        call MPI_AllREDUCE(locO2max , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. O2. =', glo
        call MPI_AllREDUCE(locO2min , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  `-> gobal min init. O2. =', glo

    end subroutine recom_init

