module mod_geom_dump
    ! Write FESOM3's computed mesh geometry in the SAME self-describing binary
    ! format as the FESOM2 oracle's src/fesom_geom_dump.F90, so tools/geom_diff.py
    ! compares them field-by-field for the M1 geometry byte-gate (closes deferred
    ! M0.7). Single-rank only (local == global), matching FESOM2's 1-rank geom run.
    !
    ! Hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Tasks 7/14): the byte gate is
    ! retired (docs/HANDOFF.md "Gate migration"); the dump is REMAPPED to the FESOM2 layout so
    ! geom_diff.py --signal=<tol> stays usable as a roundoff diagnostic on TRIANGLE meshes:
    ! elem2D_nodes(1:3), gradient_sca x = slots 1:3 / y = slots 5:7 (6 rows), edge_tri with
    ! FESOM2's 0 for a missing 2nd element (the model stores the 1st). The FESOM2 layout has no
    ! room for a 4th vertex, so a mesh with quads stops with a message.
    !
    ! Format (stream, little-endian):
    !   char8  "FGEOMDMP" | int32 nod2D,elem2D,edge2D,nl
    !   records: char24 name | int32 dtype(0=r64,1=i32) | int32 d1,d2 | data
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision, only: WP, MP
    use mod_mesh,      only: t_mesh, MAX_NV
    use mod_partit,    only: t_partit
    implicit none
    private
    public :: geom_dump_write

contains

    subroutine geom_dump_write(mesh, partit, path)
        ! M2.12a: dump per-rank LOCAL OWNED arrays (1..myDim_*), mirroring the FESOM2
        ! oracle shim (src/fesom_geom_dump.F90). At npes==1 the owned counts equal the
        ! global counts and the path is unchanged, so the proven 1-rank geom gate is
        ! byte-for-byte unaffected; at npes>1 each rank writes <path>.<mype5>.
        type(t_mesh),     intent(in) :: mesh
        type(t_partit),   intent(in) :: partit
        character(len=*), intent(in) :: path
        character(len=6)  :: rsuf
        character(len=:), allocatable :: fpath
        integer :: u, ios, ne, nn, n2, nl
        real(kind=MP), allocatable :: g6(:,:)
        integer, allocatable :: tri2(:,:)
        ne = partit%myDim_elem2D; nn = partit%myDim_nod2D; n2 = partit%myDim_edge2D; nl = mesh%nl
        if (any(mesh%elem2D_nnodes(1:ne) == MAX_NV)) then
            write(*,'(a)') 'geom_dump_write: the FESOM2 geometry layout is triangle-only; this mesh has quads'
            error stop 1
        end if
        fpath = trim(path)
        if (partit%npes /= 1) then
            write(rsuf,'(i5.5)') partit%mype
            fpath = trim(path)//'.'//rsuf
        end if
        open(newunit=u, file=fpath, status='replace', form='unformatted', &
             access='stream', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') 'geom_dump_write: cannot open '//fpath; error stop 1
        end if
        write(u) 'FGEOMDMP'
        write(u) int(nn,int32), int(ne,int32), int(n2,int32), int(nl,int32)
        call wr_r2(u, 'coord_nod2D',       mesh%coord_nod2D(1:2, 1:nn))
        call wr_r2(u, 'geo_coord_nod2D',   mesh%geo_coord_nod2D(1:2, 1:nn))
        call wr_i2(u, 'elem2D_nodes',      mesh%elem2D_nodes(1:3, 1:ne))
        call wr_r1(u, 'elem_area',         mesh%elem_area(1:ne))
        call wr_r1(u, 'elem_cos',          mesh%elem_cos(1:ne))
        call wr_r1(u, 'metric_factor',     mesh%metric_factor(1:ne))
        ! FESOM2 layout (x 1:3, y 4:6), remapped from the 4-slot one (x 1:4, y 5:8)
        allocate(g6(6, ne))
        g6(1:3, :) = mesh%gradient_sca(1:3, 1:ne); g6(4:6, :) = mesh%gradient_sca(5:7, 1:ne)
        call wr_r2(u, 'gradient_sca',      g6)
        deallocate(g6)
        call wr_r2(u, 'edge_dxdy',         mesh%edge_dxdy(1:2, 1:n2))    ! R7: METRES
        call wr_r1(u, 'edge_len',          mesh%edge_len(1:n2))          ! R7: METRES
        call wr_r2(u, 'edge_cross_dxdy',   mesh%edge_cross_dxdy(1:4, 1:n2))
        call wr_r1(u, 'area',              mesh%area(1:nn))
        call wr_r1(u, 'areasvol',          mesh%areasvol(1:nn))
        call wr_r1(u, 'area_inv',          mesh%area_inv(1:nn))
        call wr_r1(u, 'areasvol_inv',      mesh%areasvol_inv(1:nn))
        call wr_r1(u, 'mesh_resolution',   mesh%mesh_resolution(1:nn))    ! M4 GM: scalar cell resolution
        call wr_i2(u, 'edges',             mesh%edges(1:2, 1:n2))
        ! FESOM2 encoding: a boundary edge has edge_tri(2) = 0 (the model repeats edge_tri(1))
        tri2 = mesh%edge_tri(1:2, 1:n2)
        where (tri2(2,:) == tri2(1,:)) tri2(2,:) = 0
        call wr_i2(u, 'edge_tri',          tri2)
        deallocate(tri2)
        call wr_i1(u, 'nlevels',           mesh%nlevels(1:ne))
        call wr_i1(u, 'ulevels',           mesh%ulevels(1:ne))
        call wr_i1(u, 'nlevels_nod2D',     mesh%nlevels_nod2D(1:nn))
        call wr_i1(u, 'nlevels_nod2D_min', mesh%nlevels_nod2D_min(1:nn))
        close(u)
        write(*,'(a,i0,a)') 'geom_dump_write: rank ', partit%mype, ' wrote '//fpath
    end subroutine geom_dump_write

    subroutine wr_r1(u, name, a)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        real(kind=MP), intent(in) :: a(:)
        character(len=24) :: nm
        nm = name
        write(u) nm, int(0,int32), int(size(a),int32), int(1,int32)
        write(u) real(a, real64)
    end subroutine wr_r1

    subroutine wr_r2(u, name, a)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        real(kind=MP), intent(in) :: a(:,:)
        character(len=24) :: nm
        nm = name
        write(u) nm, int(0,int32), int(size(a,1),int32), int(size(a,2),int32)
        write(u) real(a, real64)
    end subroutine wr_r2

    subroutine wr_i1(u, name, a)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        integer, intent(in) :: a(:)
        character(len=24) :: nm
        nm = name
        write(u) nm, int(1,int32), int(size(a),int32), int(1,int32)
        write(u) int(a, int32)
    end subroutine wr_i1

    subroutine wr_i2(u, name, a)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        integer, intent(in) :: a(:,:)
        character(len=24) :: nm
        nm = name
        write(u) nm, int(1,int32), int(size(a,1),int32), int(size(a,2),int32)
        write(u) int(a, int32)
    end subroutine wr_i2

end module mod_geom_dump
