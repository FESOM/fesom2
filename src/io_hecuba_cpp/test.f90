program test
    use hecuba_exporter_lib
    implicit none
    type(hecuba_exporter) :: hb

    ! Create an object of type foo
    hb = hecuba_exporter(3, 4)

    ! Call bound procedures (member functions)
    !write(*,*) hb%load_datamodel(10d0), " should be ", 14.0d0
    write(*,*) hb%retrieve_data("tablesst")
    ! The destructor should be called automatically here, but this is not yet
    ! implemented in gfortran. So let's do it manually.
#ifdef __GNUC__
    call hb%delete
#endif
end program
