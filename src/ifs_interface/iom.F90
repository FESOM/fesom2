!=====================================================
! Input/Output manager :  Library to write output files
!
! -Original code for NEMOv40 by ECMWF.
! -Adapted to FESOM2 by Razvan Aguridan, ECMWF, 2023.
!-----------------------------------------------------

MODULE iom
#if defined(__MULTIO)
    USE multio_api
    USE, INTRINSIC :: iso_fortran_env, only: real64

    IMPLICIT NONE
    PRIVATE

    TYPE(multio_handle) :: mio_handle
    INTEGER(8), PRIVATE :: mio_parent_comm

    PUBLIC iom_enable_multio
    PUBLIC iom_initialize, iom_init_server, iom_finalize
    PUBLIC iom_send_fesom_domains
    PUBLIC iom_field_request, iom_send_fesom_data
    PUBLIC iom_flush

    LOGICAL :: lnomultio = .TRUE.

    PRIVATE ctl_stop
    !!----------------------------------------------------------------------
    !! NEMO/OCE 4.0 , NEMO Consortium (2018)
    !! $Id: iom.F90 13297 2020-07-13 08:01:58Z andmirek $
    !! Software governed by the CeCILL license (see ./LICENSE)
    !!----------------------------------------------------------------------

    TYPE iom_field_request
        CHARACTER(100)                          :: name     = REPEAT(" ", 100)
        CHARACTER(100)                          :: category = REPEAT(" ", 100)
        CHARACTER(6)                            :: gridType = REPEAT(" ", 6)
        REAL(real64), DIMENSION(:), POINTER     :: values => NULL()
        INTEGER                                 :: globalSize = 0
        INTEGER                                 :: sampleInterval=0
        INTEGER                                 :: level = 0
        INTEGER                                 :: step = 0
        INTEGER                                 :: currentDate,  currentTime
        INTEGER                                 :: previousDate, previousTime
        INTEGER                                 :: startDate,    startTime
        INTEGER                                 :: lastcounter
    END TYPE

CONTAINS

    SUBROUTINE iom_enable_multio()
        IMPLICIT NONE
        lnomultio = .FALSE.
    END SUBROUTINE

    SUBROUTINE multio_custom_error_handler(context, err, info)
        USE mpi
        USE, intrinsic :: iso_fortran_env, ONLY: int64
        USE :: multio_api_constants_mod, ONLY: multio_failure_info

        IMPLICIT NONE

        INTEGER(int64),                 INTENT(INOUT)   :: context  ! Use mpi communicator as context
        INTEGER,                        INTENT(IN)      :: err
        TYPE(multio_failure_info),      INTENT(in)      :: info
        INTEGER                                         :: mpierr

        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop( 'MULTIO ERROR: ', multio_error_string(err, info))
            IF (context /= MPI_UNDEFINED) THEN
                CALL mpi_abort(int(context), MPI_ERR_OTHER, mpierr)
                context = MPI_UNDEFINED
            ENDIF
        ENDIF
    END SUBROUTINE

    SUBROUTINE iom_initialize(client_id, local_comm, return_comm, global_comm )
        USE mpi
        USE :: multio_api, ONLY: failure_handler_t

        IMPLICIT NONE

        CHARACTER(LEN=*),               INTENT(IN)                  :: client_id
        INTEGER,                        INTENT(IN),     OPTIONAL    :: local_comm
        INTEGER,                        INTENT(OUT),    OPTIONAL    :: return_comm
        INTEGER,                        INTENT(IN),     OPTIONAL    :: global_comm

        TYPE(multio_configuration)                                  :: conf_ctx
        INTEGER                                                     :: err
        CHARACTER(len=16)                                           :: err_str
        PROCEDURE(failure_handler_t),   POINTER                     :: pf

        IF (lnomultio) RETURN

        mio_parent_comm = mpi_comm_world

        err = multio_initialise()
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('Initializing multio failed: ', multio_error_string(err))
        END IF

        IF (PRESENT(global_comm)) THEN
            mio_parent_comm = global_comm
        ENDIF

        ! Prepare context and check errors explicitly until everything is set up - then failure handler is used
        BLOCK
            CHARACTER(:), allocatable :: config_file
            INTEGER :: config_file_length

            CALL get_environment_variable('MULTIO_FESOM_CONFIG_FILE', length=config_file_length)
            IF (config_file_length == 0) THEN
                call ctl_stop('The fesom plan file is not correctly set!')
                err = conf_ctx%new()
            ELSE
                ALLOCATE(character(len=config_file_length + 1) :: config_file)

                CALL get_environment_variable('MULTIO_FESOM_CONFIG_FILE', config_file)
                err = conf_ctx%new(config_file)

                DEALLOCATE(config_file)
            ENDIF
        END BLOCK

        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('Creating multio configuration context failed: ', multio_error_string(err))
        END IF

        ! Setting a failure handler that reacts on interface problems or exceptions that are not handled within the interface
        pf => multio_custom_error_handler
        err = conf_ctx%set_failure_handler(pf, mio_parent_comm)
        if (err /= MULTIO_SUCCESS) then
            CALL ctl_stop( 'setting multio failure handler failed: ', multio_error_string(err))
        end if

        err = conf_ctx%mpi_allow_world_default_comm(.FALSE.)
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('conf_ctx%mpi_allow_world_default_comm(.FALSE.) failed: ', multio_error_string(err))
        END IF

        err = conf_ctx%mpi_return_client_comm(return_comm)
        IF (err /= MULTIO_SUCCESS) THEN
            WRITE (err_str, "(I10)") return_comm
            CALL ctl_stop('conf_ctx%mpi_return_client_comm(', err_str,') failed: ', multio_error_string(err))
        END IF

        err = conf_ctx%mpi_parent_comm(int(mio_parent_comm))
        IF (err /= MULTIO_SUCCESS) THEN
            WRITE (err_str, "(I10)") mio_parent_comm
            CALL ctl_stop('conf_ctx%mpi_parent_comm(', err_str,') failed: ', multio_error_string(err))
        END IF

        err = mio_handle%new(conf_ctx)
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('mio_handle%new(conf_ctx) failed: ', multio_error_string(err))
        END IF

        err = conf_ctx%delete()
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('conf_ctx%delete() failed: ', multio_error_string(err))
        END IF

        err = mio_handle%open_connections();
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('mio_handle%open_connections failed: ', multio_error_string(err))
        END IF
    END SUBROUTINE iom_initialize

    SUBROUTINE iom_finalize()
        IMPLICIT NONE
        INTEGER :: err

        IF (lnomultio) RETURN

        err = mio_handle%close_connections();
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('mio_handle%close_connections failed: ', multio_error_string(err))
        END IF

        err = mio_handle%delete();
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('mio_handle%delete failed: ', multio_error_string(err))
        END IF
    END SUBROUTINE iom_finalize

    SUBROUTINE iom_init_server(server_comm)
        USE :: multio_api, ONLY: failure_handler_t

        IMPLICIT NONE

        INTEGER,                         INTENT(IN)  :: server_comm

        TYPE(multio_configuration)                   :: conf_ctx
        INTEGER                                      :: err
        CHARACTER(len=16)                            :: err_str
        PROCEDURE(failure_handler_t),    POINTER     :: pf

        IF (lnomultio) RETURN

        mio_parent_comm = server_comm

        err = multio_initialise()
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('Initializing multio failed: ', multio_error_string(err))
        END IF

        ! Prepare context and check errors explicitly until everything is set up - then failure handler is used

        BLOCK
            CHARACTER(:), allocatable :: config_file
            INTEGER :: config_file_length

            CALL get_environment_variable('MULTIO_FESOM_CONFIG_FILE', length=config_file_length)
            IF (config_file_length == 0) THEN
                err = conf_ctx%new()
            ELSE
                ALLOCATE(character(len=config_file_length + 1) :: config_file)

                CALL get_environment_variable('MULTIO_FESOM_CONFIG_FILE', config_file)
                err = conf_ctx%new(config_file)

                DEALLOCATE(config_file)
            ENDIF
        END BLOCK

        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('Creating multio server configuration context failed: ', multio_error_string(err))
        END IF

        ! Setting a failure handler that reacts on interface problems or exceptions that are not handled within the interface
        ! Set handler before invoking blocking start server call
        pf => multio_custom_error_handler
        err = conf_ctx%set_failure_handler(pf, mio_parent_comm)
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('setting multio failure handler failed: ', multio_error_string(err))
        END IF

        err = conf_ctx%mpi_allow_world_default_comm(.FALSE.)
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('conf_ctx%mpi_allow_world_default_comm(.FALSE.) failed: ', multio_error_string(err))
        END IF

        err = conf_ctx%mpi_parent_comm(int(mio_parent_comm))
        IF (err /= MULTIO_SUCCESS) THEN
            WRITE (err_str, "(I10)") mio_parent_comm
            CALL ctl_stop('conf_ctx%mpi_parent_comm(', err_str,') failed: ', multio_error_string(err))
        END IF
        ! Blocking call
        err = multio_start_server(conf_ctx)
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('multio_start_server(conf_ctx) failed: ', multio_error_string(err))
        END IF

        err = conf_ctx%delete()
        IF (err /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('conf_ctx%delete() failed: ', multio_error_string(err))
        END IF
    END SUBROUTINE iom_init_server

    SUBROUTINE iom_send_fesom_domains(partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT

        IMPLICIT NONE

        TYPE(multio_metadata)              :: md
        INTEGER                            :: cerr
        INTEGER                            :: elem, elnodes(3), aux
        TYPE(t_partit), INTENT(IN), TARGET :: partit
        TYPE(t_mesh),   intent(in), TARGET :: mesh
        INTEGER, DIMENSION(:), POINTER     :: temp

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

        IF (lnomultio) RETURN

        cerr = md%new(mio_handle)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: ngrid, md%new() failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("name", "N grid")
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: ngrid, md%set_string(name) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("category", "fesom-domain-nodemap")
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: ngrid, md%set_string(category) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("representation", "unstructured")
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: ngrid, md%set_string(representation) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("globalSize", mesh%nod2D)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: ngrid, md%set_int(globalSize) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("toAllServers", .TRUE._1)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: ngrid, md%set_bool(toAllServers) failed: ', multio_error_string(cerr))
        END IF

        temp => partit%myList_nod2D(1:partit%myDim_nod2D)
        cerr = mio_handle%write_domain(md, temp - 1)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: ngrid, mio_handle%write_domain() failed: ', multio_error_string(cerr))
        END IF

        cerr = md%delete()
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: ngrid, md%delete() failed: ', multio_error_string(cerr))
        END IF

        !declare grid at elements
        cerr = md%new(mio_handle)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: egrid, md%new() failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("name", "C grid")
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: egrid, md%set_string(name) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("category", "fesom-domain-elemmap")
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: egrid, md%set_string(category) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("representation", "unstructured")
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: egrid, md%set_string(representation) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("globalSize", mesh%elem2D)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: egrid, md%set_int(globalSize) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("toAllServers", .TRUE._1)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: egrid, md%set_bool(toAllServers) failed: ', multio_error_string(cerr))
        END IF

        cerr = mio_handle%write_domain(md, partit%myList_elem2D(partit%myInd_elem2D_shrinked) - 1)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: egrid, mio_handle%write_domain() failed: ', multio_error_string(cerr))
        END IF

        cerr = md%delete()
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_domains: egrid, md%delete() failed: ', multio_error_string(cerr))
        END IF
    END SUBROUTINE iom_send_fesom_domains

    SUBROUTINE iom_send_fesom_data(data)
        USE g_clock
        USE g_config, only: MeshId
        IMPLICIT NONE

        TYPE(iom_field_request), INTENT(INOUT)  :: data
        INTEGER                                 :: cerr
        TYPE(multio_metadata)                   :: md

        IF (lnomultio) RETURN

        cerr = md%new(mio_handle)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%new() failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("category", data%category)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_string(category) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("globalSize", data%globalSize)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_int(globalSize) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("level", data%level)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_int(level) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("toAllServers", .FALSE._1)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_bool(toAllServers) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("name", trim(data%name))
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_string(name) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("gridType", "unstructured_grid")
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_string(gridType) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("unstructuredGridType", MeshId)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_string(unstructuredGridType) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("unstructuredGridSubtype", data%gridType(1:1))
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_string(unstructuredGridSubtype) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("operation", "average")
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_string(operation) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("domain", data%gridType)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%set_string(domain) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set("currentDate",    data%currentDate)
        cerr = md%set("currentTime",    data%currentTime)
        cerr = md%set("previousDate",   data%previousDate)
        cerr = md%set("previousTime",   data%previousTime)
        cerr = md%set("startDate",      data%startDate)
        cerr = md%set("startTime",      data%startTime)
        cerr = md%set("sampleInterval", data%sampleInterval)
!       cerr = md%set_int("sampleIntervalInSeconds", data%sampleInterval)
        cerr = md%set("sampleIntervalUnit", 'S')
        cerr = md%set("sampleIntervalInSeconds", data%sampleInterval)
        cerr = md%set("timeStep",                data%sampleInterval) !we do not distinguish between the timestep & sampling interval legacy code for MULTIO
        cerr = md%set("step-frequency",          data%lastcounter)
        cerr = md%set("step",                    data%step)
        IF (cerr /= MULTIO_SUCCESS) THEN
           CALL ctl_stop('send_fesom_data: md%set_int(date) failed: ', multio_error_string(cerr))
        END IF

        cerr = mio_handle%write_field(md, data%values)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: mio_handle%write_field failed: ', multio_error_string(cerr))
        END IF

        cerr = md%delete()
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('send_fesom_data: md%delete failed: ', multio_error_string(cerr))
        END IF
    END SUBROUTINE

    SUBROUTINE iom_flush(domain, step)
        IMPLICIT NONE

        CHARACTER(6), INTENT(IN)                :: domain
        INTEGER, INTENT(IN)                     :: step

        INTEGER                                 :: cerr
        TYPE(multio_metadata)                   :: md

        IF (lnomultio) RETURN

        cerr = md%new(mio_handle)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('iom_flush: md%new() failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set_bool("toAllServers", .TRUE._1)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('iom_flush: md%set_bool(toAllServers) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set_string("domain", domain)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('iom_flush: md%set_string(domain) failed: ', multio_error_string(cerr))
        END IF

        cerr = md%set_int("step", step)
        IF (cerr /= MULTIO_SUCCESS) THEN
           CALL ctl_stop('iom_flush: md%set_int(step) failed: ', multio_error_string(cerr))
        END IF

        cerr = mio_handle%flush(md)
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('iom_flush: mio_handle%multio_flush failed: ', multio_error_string(cerr))
        END IF

        cerr = md%delete()
        IF (cerr /= MULTIO_SUCCESS) THEN
            CALL ctl_stop('iom_flush: md%delete failed: ', multio_error_string(cerr))
        END IF
    END SUBROUTINE

    SUBROUTINE ctl_stop(m1, m2, m3, m4)
        USE mpi

        IMPLICIT NONE
        CHARACTER(len=*), INTENT(in), OPTIONAL :: m1, m2, m3, m4
        INTEGER :: dummy

        IF ( PRESENT(m1) ) WRITE(*,*) m1
        IF ( PRESENT(m2) ) WRITE(*,*) m2
        IF ( PRESENT(m3) ) WRITE(*,*) m3
        IF ( PRESENT(m4) ) WRITE(*,*) m4

        CALL mpi_abort(mpi_comm_world, 1, dummy)
    END SUBROUTINE ctl_stop

    !!======================================================================
#endif
END MODULE iom
