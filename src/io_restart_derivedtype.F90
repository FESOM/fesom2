module restart_derivedtype_module
    interface
        subroutine write_all_bin_restarts(ctarr, path_in, pathi_in, partit, mesh, ice, dynamics, tracers)
            use MOD_ICE
            use MOD_DYN
            use MOD_TRACER
            use MOD_PARTIT
            use MOD_MESH
            integer, dimension(3) , intent(in)                   :: ctarr
            character(len=*),       intent(in)                   :: path_in
            character(len=*),       intent(in)                   :: pathi_in
            type(t_partit),         intent(in), target           :: partit
            type(t_mesh)  ,         intent(in), target           :: mesh
            type(t_ice)   ,         intent(in), target, optional :: ice
            type(t_dyn)   ,         intent(in), target, optional :: dynamics
            type(t_tracer),         intent(in), target, optional :: tracers
        end subroutine write_all_bin_restarts

        subroutine read_all_bin_restarts(path_in, partit, mesh, ice, dynamics, tracers)
            use MOD_ICE
            use MOD_DYN
            use MOD_TRACER
            use MOD_PARTIT
            use MOD_MESH
            character(len=*), intent(in)                    :: path_in
            type(t_partit), intent(inout), target           :: partit
            type(t_mesh)  , intent(inout), target           :: mesh
            type(t_ice)   , intent(inout), target, optional :: ice
            type(t_dyn)   , intent(inout), target, optional :: dynamics
            type(t_tracer), intent(inout), target, optional :: tracers
        end subroutine read_all_bin_restarts
    end interface
end module restart_derivedtype_module
!
!
!_______________________________________________________________________________
subroutine write_all_bin_restarts(ctarr, path_in, pathi_in, partit, mesh, ice, dynamics, tracers)
    use MOD_ICE
    use MOD_DYN
    use MOD_TRACER
    use MOD_PARTIT
    use MOD_MESH
    use fortran_utils
    implicit none

    integer, dimension(3) , intent(in)           :: ctarr ! //cstep,ctime,cyear//
    character(len=*)      , intent(in)           :: path_in
    character(len=*)      , intent(in)           :: pathi_in
    type(t_partit), target, intent(in)           :: partit
    type(t_mesh)  , target, intent(in)           :: mesh
    type(t_ice)   , target, intent(in), optional :: ice
    type(t_dyn)   , target, intent(in), optional :: dynamics
    type(t_tracer), target, intent(in), optional :: tracers

    ! EO parameters
    integer fileunit, fileunit_i

    !___________________________________________________________________________
    ! write info file
    if(partit%mype == 0) then
        print *, achar(27)//'[1;33m'//' --> writing derived type binary restarts to '//trim(path_in)//achar(27)//'[0m'
        ! store metadata about the raw restart
        fileunit_i = 299
        open(newunit = fileunit_i, file = trim(pathi_in))
        write(fileunit_i, '(g0)') ctarr(1)
        write(fileunit_i, '(g0)') ctarr(2)
        write(fileunit_i, '(2(g0))') "! year: ",ctarr(3)
    end if

    !___________________________________________________________________________
    ! mesh derived type
    fileunit = partit%mype+300
    open(newunit = fileunit, &
        file     = trim(path_in)//'/'//'t_mesh.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
        status   = 'replace', &
        form     = 'unformatted')
    write(fileunit) mesh
    close(fileunit)
    if(partit%mype == 0) then
        write(fileunit_i, '(1(g0))') "!   t_mesh"
        print *, achar(27)//'[33m'//'     > write derived type t_mesh'//achar(27)//'[0m'
    end if

    !___________________________________________________________________________
    ! partit derived type
    fileunit = partit%mype+300
    open(newunit = fileunit, &
        file     = trim(path_in)//'/'//'t_partit.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
        status   = 'replace', &
        form     = 'unformatted')
    write(fileunit) partit
    close(fileunit)
    if(partit%mype == 0) then
        write(fileunit_i, '(1(g0))') "!   t_partit"
        print *, achar(27)//'[33m'//'     > write derived type t_partit'//achar(27)//'[0m'
    end if

    !___________________________________________________________________________
    ! tracer derived type
    if (present(tracers)) then
        fileunit = partit%mype+300
        open(newunit = fileunit, &
            file     = trim(path_in)//'/'//'t_tracer.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
            status   = 'replace', &
            form     = 'unformatted')
        write(fileunit) tracers
        close(fileunit)
        if(partit%mype == 0) then
            write(fileunit_i, '(1(g0))') "!   t_tracer"
            print *, achar(27)//'[33m'//'     > write derived type t_tracer'//achar(27)//'[0m'
        end if
    end if

    !___________________________________________________________________________
    ! dynamics derived type
    if (present(dynamics)) then
        fileunit = partit%mype+300
        open(newunit = fileunit, &
            file     = trim(path_in)//'/'//'t_dynamics.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
            status   = 'replace', &
            form     = 'unformatted')
        write(fileunit) dynamics
        close(fileunit)
        if(partit%mype == 0) then
            write(fileunit_i, '(1(g0))') "!   t_dynamics"
            print *, achar(27)//'[33m'//'     > write derived type t_dynamics'//achar(27)//'[0m'
        end if
    end if

    !___________________________________________________________________________
    ! ice derived type
    if (present(ice)) then
        fileunit = partit%mype+300
        open(newunit = fileunit, &
            file    = trim(path_in)//'/'//'t_ice.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
            status  = 'replace', &
            form    = 'unformatted')
        write(fileunit) ice
        close(fileunit)
        if(partit%mype == 0) then
            write(fileunit_i, '(1(g0))') "!   t_ice"
            print *, achar(27)//'[33m'//'     > write derived type t_ice'//achar(27)//'[0m'
        end if
    end if
    !___________________________________________________________________________
    if(partit%mype == 0) close(fileunit_i)
end subroutine write_all_bin_restarts
!
!
!_______________________________________________________________________________
! read derived type binary restart files, depending on input (see optional) not
! all derived type binaries are read --> functionalitiy for dwarfs !
subroutine read_all_bin_restarts(path_in, partit, mesh, ice, dynamics, tracers)
    use MOD_ICE
    use MOD_DYN
    use MOD_TRACER
    use MOD_PARTIT
    use MOD_MESH
    use fortran_utils
    implicit none

    ! do optional here for the usage with dwarfs, since there only specific derived
    ! types will be needed
    character(len=*), intent(in)                    :: path_in
    type(t_partit), intent(inout), target           :: partit
    type(t_mesh)  , intent(inout), target           :: mesh
    type(t_ice)   , intent(inout), target, optional :: ice
    type(t_dyn)   , intent(inout), target, optional :: dynamics
    type(t_tracer), intent(inout), target, optional :: tracers
    integer fileunit

    !___________________________________________________________________________
    if (partit%mype==0) print *, achar(27)//'[1;33m'//' --> read restarts from derived type binary'//achar(27)//'[0m'
    !___________________________________________________________________________
    ! mesh derived type
    fileunit = partit%mype+300
    open( newunit = fileunit, &
          file     = trim(path_in)//'/'//'t_mesh.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
          status   = 'old', &
          form     = 'unformatted')
    read(fileunit) mesh
    close(fileunit)
    if (partit%mype==0) print *, achar(27)//'[33m'//'     > read derived type t_mesh'//achar(27)//'[0m'

    !___________________________________________________________________________
    ! partit derived type
    fileunit = partit%mype+300
    open(newunit = fileunit, &
         file     = trim(path_in)//'/'//'t_partit.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
         status   = 'old', &
         form     = 'unformatted')
    read(fileunit) partit
    close(fileunit)
    if (partit%mype==0) print *, achar(27)//'[33m'//'     > read derived type t_partit'//achar(27)//'[0m'

    !___________________________________________________________________________
    ! tracer derived type
    if (present(tracers)) then
        fileunit = partit%mype+300
        open(newunit = fileunit, &
            file     = trim(path_in)//'/'//'t_tracer.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
            status   = 'old', &
            form     = 'unformatted')
        read(fileunit) tracers
        close(fileunit)
        if (partit%mype==0) print *, achar(27)//'[33m'//'     > read derived type t_tracer'//achar(27)//'[0m'
    end if

    !___________________________________________________________________________
    ! dynamics derived type
    if (present(dynamics)) then
        fileunit = partit%mype+300
        open(newunit = fileunit, &
            file     = trim(path_in)//'/'//'t_dynamics.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
            status   = 'old', &
            form     = 'unformatted')
        read(fileunit) dynamics
        close(fileunit)
        if (partit%mype==0) print *, achar(27)//'[33m'//'     > read derived type t_dynamics'//achar(27)//'[0m'
    end if

    !___________________________________________________________________________
    ! ice derived type
    if (present(ice)) then
        fileunit = partit%mype+300
        open(newunit = fileunit, &
            file    = trim(path_in)//'/'//'t_ice.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
            status  = 'old', &
            form    = 'unformatted')
        read(fileunit) ice
        close(fileunit)
        if (partit%mype==0) print *, achar(27)//'[33m'//'     > read derived type t_ice'//achar(27)//'[0m'
    end if
end subroutine read_all_bin_restarts
