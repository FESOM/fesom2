! This is just a test to see whether I can add extra source files and
! FESOM2 remains compilable!

!==========================================================
MODULE MOD_NEURALNET
  USE o_PARAM
  USE MOD_READ_BINARY_ARRAYS !, ONLY: read1d_real, read2d_real, read1d_char, read1d_int
  USE MOD_MESH
  USE, INTRINSIC :: ieee_arithmetic
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  PRIVATE
  PUBLIC :: assign_input_values, read_nn_architecture, read_nn_weights, read_nn_biases, read_nn_activation 
  PUBLIC :: forward_pass_single, forward_pass_full
  PUBLIC :: neuralnet_init, get_neural_net
  PUBLIC :: t_neural_net

  TYPE t_neural_net
      INTEGER :: nlayers
      INTEGER, ALLOCATABLE :: layer_sizes(:)
      REAL(kind=WP), ALLOCATABLE :: weights(:,:,:), biases(:,:)
      CHARACTER(LEN=16), ALLOCATABLE :: activations(:)
  END TYPE t_neural_net
  
  ! Module-level variable - initialized once, shared across all code
  TYPE(t_neural_net), TARGET :: nn_gm_module
  LOGICAL :: nn_initialized = .FALSE.

  CONTAINS

    SUBROUTINE assign_input_values(input_var_names, num_input_var_names, z, n, ordered_nb_list, input_data, n_input_data)
      ! This is a placeholder for the logic to extract the input data from the FESOM data structure based on the variable name. 
      ! This will likely involve a series of IF statements or a SELECT CASE statement to match the variable names to the corresponding data extraction logic.
      CHARACTER(LEN=*), INTENT(IN) :: input_var_names(:)
      REAL(kind=WP), INTENT(OUT) :: input_data(:), n_input_data(:) 
      INTEGER, INTENT(IN) :: z, n ! These are the vertical level and node index for which we want to extract the data. They are needed to know which value to extract from the FESOM data structure.
      INTEGER :: i ! Current index in loop over input variables. var_name and var_name_nb... are successive in the input_var_names.txt, so we can just keep incrementing j to fill the input_data array in the correct order.
      INTEGER :: nb, u
      INTEGER :: max_num_nb = 6 ! Maximum number of neighbours in dbgyre mesh.
      INTEGER :: num_input_var_names
      INTEGER, INTENT(IN) :: ordered_nb_list(:) ! This is the list of neighbour node indices, ordered ascending.
      CHARACTER(LEN=32) :: name, fname
      REAL(kind=WP), ALLOCATABLE :: means(:), stds(:) ! Arrays to hold the means and stds for all input variables, to be read from ./normalization_params. This is needed for normalization of the input data.

      ALLOCATE(means(num_input_var_names), stds(num_input_var_names))

      OPEN(NEWUNIT=u, FILE=TRIM(fname), STATUS='old', ACTION='read', FORM='unformatted', ACCESS='stream')
      DO i = 1, num_input_var_names
            ! Loop over input variable names. For each variable name, extract the corresponding data from the FESOM data structure and fill the input_data array. 
            ! Also apply normalization to the input data using the mean and std from the training data, which are stored in ./normalization_params.
          name = input_var_names(i)
          WRITE(fname, '(A,A,A)') './normalization_params/mean/', input_var_names(i), '.bin'
          READ(u) means(i)
          WRITE(fname, '(A,A,A)') './normalization_params/std/', input_var_names(i), '.bin'
          READ(u) stds(i)
          SELECT CASE (name)
              CASE ('temp')
                  ! input_data(i) = tracers%data(1)%values(z,n)
                  DO nb=1, max_num_nb
                      IF ( ordered_nb_list(nb) == -1 ) THEN
                          input_data(i+nb) = means(i) ! or, equivalently, the mean value from the training data after normalization
                      ELSE
                          !  input_data(i+nb) = tracers%data(1)%values(z, ordered_nb_list(nb))
                      END IF
                  END DO
              CASE ('unod')
                  ! input_data(i) = dynamics%uvnode(1,z,n)
                  DO nb=1, max_num_nb
                      IF ( ordered_nb_list(nb) == -1 ) THEN
                          input_data(i+nb) = means(i) ! or, equivalently, the mean value from the training data after normalization
                      ELSE
                          !  input_data(i+nb) = dynamics%uvnode(1,z, ordered_nb_list(nb))
                      END IF
                  END DO
              CASE ('vnod')
                  ! input_data(i) = dynamics%uvnode(2,z,n)
                  DO nb=1, max_num_nb
                      IF ( ordered_nb_list(nb) == -1 ) THEN
                          input_data(i+nb) = means(i) ! or, equivalently, the mean value from the training data after normalization
                      ELSE
                          !  input_data(i+nb) = dynamics%uvnode(2,z, ordered_nb_list(nb))
                      END IF
                  END DO
              CASE ('curl_u_nn')
                  ! Call customized subroutine to compute curl_u_nn as it is not computed by default
                  ! curl_u_nn = ...
                  DO nb=1, max_num_nb
                      IF ( ordered_nb_list(nb) == -1 ) THEN
                          input_data(i+nb) = means(i) ! or, equivalently, the mean value from the training data after normalization
                      ELSE
                          !  input_data(i+nb) = curl_u_nn(z, ordered_nb_list(nb))
                      END IF
                  END DO
              CASE ('slope_x')
                  ! input_data(i) = neutral_slope(1,z,n)
                  DO nb=1, max_num_nb
                      IF ( ordered_nb_list(nb) == -1 ) THEN
                          input_data(i+nb) = means(i) ! or, equivalently, the mean value from the training data after normalization
                      ELSE
                          !  input_data(i+nb) = neutral_slope(1,z, ordered_nb_list(nb))
                      END IF
                  END DO
              CASE ('slope_y')
                  ! input_data(i) = neutral_slope(2,z,n)
                  DO nb=1, max_num_nb
                      IF ( ordered_nb_list(nb) == -1 ) THEN
                          input_data(i+nb) = means(i) ! or, equivalently, the mean value from the training data after normalization
                      ELSE
                          !  input_data(i+nb) = neutral_slope(2,z, ordered_nb_list(nb))
                      END IF
                  END DO
              CASE ('N2')
                  ! Watch out: N2 is defined at vertical interfaces, so it needs to be interpolated to the vertical levels of the nodes.
                  ! input_data(i) = bvfreq(z,n)
                  DO nb=1, max_num_nb
                      IF ( ordered_nb_list(nb) == -1 ) THEN
                          input_data(i+nb) = means(i) ! or, equivalently, the mean value from the training data after normalization
                      ELSE
                          !  input_data(i+nb) = bvfreq(z, ordered_nb_list(nb))
                      END IF
                  END DO
              CASE ('ld_baroc1')
                  ! input_data(i) = rosb(n)
              CASE DEFAULT
                  ! WRITE(*,(A,A,A)) 'This variable (' // var_name // ') probably ends with _nb and its value is already filled. Skipping...'
          END SELECT

          WRITE(*,*) 'Extracted input variable: ', name, ' with original value: ', input_data(i), 'and normalized value: '

      END DO
      CLOSE(u)
      ! Now normalize by doing array computation
      n_input_data = (input_data - means) / stds
      DEALLOCATE(means, stds)
    END SUBROUTINE assign_input_values

    SUBROUTINE neuralnet_init(mype, mpi_comm)
        ! Initialize the global neural network (only once)
        ! Subsequent calls are no-ops
        INTEGER, INTENT(IN), OPTIONAL :: mype, mpi_comm
        INTEGER :: i, nlayers, my_rank, comm, ierr
        CHARACTER(LEN=256) :: path, weights_path, biases_path, act_path
        CHARACTER(LEN=4) :: i_str
        LOGICAL :: use_mpi
        
        ! Check if already initialized
        IF (nn_initialized) RETURN
        
        ! Determine if we should use MPI
        use_mpi = PRESENT(mype) .AND. PRESENT(mpi_comm)
        
        IF (use_mpi) THEN
            my_rank = mype
            comm = mpi_comm
        ELSE
            my_rank = 0
            comm = -1
        END IF
        
        ! 1. Read network architecture (only on rank 0)
        IF (my_rank == 0) THEN
            CALL read_nn_architecture(nn_gm_module%nlayers, nn_gm_module%layer_sizes)
        END IF
        
        ! Broadcast nlayers to all ranks
        IF (use_mpi) THEN
            CALL MPI_BCAST(nn_gm_module%nlayers, 1, MPI_INTEGER, 0, comm, ierr)
        END IF
        
        nlayers = nn_gm_module%nlayers
        
        ! Allocate layer_sizes on all ranks (needed before broadcast)
        IF (my_rank /= 0) THEN
            ALLOCATE(nn_gm_module%layer_sizes(nlayers+1))
        END IF
        
        ! Broadcast layer_sizes
        IF (use_mpi) THEN
            CALL MPI_BCAST(nn_gm_module%layer_sizes, nlayers+1, MPI_INTEGER, 0, comm, ierr)
        END IF
        
        ! 2. Allocate arrays for weights, biases, activations on all ranks
        ALLOCATE(nn_gm_module%weights(nlayers, maxval(nn_gm_module%layer_sizes), maxval(nn_gm_module%layer_sizes)))
        ALLOCATE(nn_gm_module%biases(nlayers, maxval(nn_gm_module%layer_sizes)))
        ALLOCATE(nn_gm_module%activations(nlayers))
        
        ! 3. Initialize with NaN on rank 0 (others will receive via broadcast)
        IF (my_rank == 0) THEN
            nn_gm_module%weights = ieee_value(nn_gm_module%weights, ieee_quiet_nan)
            nn_gm_module%biases = ieee_value(nn_gm_module%biases, ieee_quiet_nan)
        END IF
        
        ! 4. Load each layer's weights, biases, and activation functions (only on rank 0)
        path = '/albedo/home/rjuhrban/fesom2/src/neuralnet_params'
        
        IF (my_rank == 0) THEN
            DO i = 1, nlayers
                WRITE(i_str, '(I0)') i - 1
                
                ! Build file paths
                weights_path = TRIM(path) // '/layer_' // TRIM(i_str) // '_weights.bin'
                biases_path = TRIM(path) // '/layer_' // TRIM(i_str) // '_biases.bin'
                act_path = TRIM(path) // '/layer_' // TRIM(i_str) // '_act.txt'
                
                ! Read layer parameters
                CALL read_nn_weights(weights_path, nn_gm_module%weights(i,:,:), &
                                nn_gm_module%layer_sizes(i), nn_gm_module%layer_sizes(i+1))
                CALL read_nn_biases(biases_path, nn_gm_module%biases(i,:), &
                                nn_gm_module%layer_sizes(i+1))
                CALL read_nn_activation(act_path, nn_gm_module%activations(i))
            END DO
        END IF
        
        ! 5. Broadcast weights, biases, and activations to all ranks
        IF (use_mpi) THEN
            CALL MPI_BCAST(nn_gm_module%weights, &
                          nlayers * maxval(nn_gm_module%layer_sizes) * maxval(nn_gm_module%layer_sizes), &
                          MPI_DOUBLE_PRECISION, 0, comm, ierr)
            CALL MPI_BCAST(nn_gm_module%biases, &
                          nlayers * maxval(nn_gm_module%layer_sizes), &
                          MPI_DOUBLE_PRECISION, 0, comm, ierr)
            ! For character array, broadcast as bytes
            CALL MPI_BCAST(nn_gm_module%activations, &
                          nlayers * 16, &
                          MPI_CHARACTER, 0, comm, ierr)
        END IF
        
        nn_initialized = .TRUE.

        ! Load means and stds for later normalization
        
    END SUBROUTINE neuralnet_init
    
    FUNCTION get_neural_net() RESULT(nn)
        ! Return pointer to the globally-initialized neural network
        TYPE(t_neural_net), POINTER :: nn
        nn => nn_gm_module
    END FUNCTION get_neural_net

    SUBROUTINE read_nn_architecture(nl, layer_sizes)
        INTEGER, INTENT(OUT) :: nl
        INTEGER, ALLOCATABLE, INTENT(OUT) :: layer_sizes(:)
        CHARACTER(LEN=256) :: nlname, nnname
        INTEGER :: iunit, iostat
        CHARACTER(256) :: iomsg
        
        ! Use module parameter or environment variable for path
        nlname = '/albedo/home/rjuhrban/fesom2/src/neuralnet_params/nlayers.bin'
        
        OPEN(NEWUNIT=iunit, FILE=TRIM(nlname), STATUS='old', ACTION='read', &
            FORM='unformatted', ACCESS='stream', IOSTAT=iostat, IOMSG=iomsg)
        
        IF (iostat /= 0) THEN
            WRITE(*,'(A)') 'ERROR opening nlayers.bin: ', TRIM(iomsg)
            STOP 1
        END IF
        
        READ(iunit, IOSTAT=iostat) nl
        CLOSE(iunit)
        
        IF (iostat /= 0) THEN
            WRITE(*,'(A)') 'ERROR reading nlayers.bin'
            STOP 1
        END IF
        
        ALLOCATE(layer_sizes(nl+1))
        
        nnname = '/albedo/home/rjuhrban/fesom2/src/neuralnet_params/nneurons.bin'
        
        OPEN(NEWUNIT=iunit, FILE=TRIM(nnname), STATUS='old', ACTION='read', &
            FORM='unformatted', ACCESS='stream', IOSTAT=iostat, IOMSG=iomsg)
        READ(iunit, IOSTAT=iostat) layer_sizes
        CLOSE(iunit)
        
        IF (iostat /= 0) THEN
            WRITE(*,'(A)') 'ERROR reading nneurons.bin'
            STOP 1
        END IF
    END SUBROUTINE read_nn_architecture

    SUBROUTINE read_nn_weights(filename_weights, weights, input_neurons, output_neurons)
      CHARACTER(LEN=*), INTENT(IN) :: filename_weights
      REAL(kind=WP), INTENT(INOUT) :: weights(:,:)
      INTEGER :: iunit, iostat
      CHARACTER(512) :: iomsg
      INTEGER, INTENT(IN) :: input_neurons, output_neurons
      REAL(kind=WP), ALLOCATABLE :: buffer(:)
      INTEGER :: total_elements, bytes_per_real
      
      total_elements = input_neurons * output_neurons
      ALLOCATE(buffer(total_elements))
      
      ! Calculate record length in bytes based on actual REAL kind
      bytes_per_real = STORAGE_SIZE(buffer(1)) / 8  ! bits to bytes
      
      ! Use direct access with record length = total bytes
      OPEN(NEWUNIT=iunit, FILE=filename_weights, STATUS='old', &
           FORM='unformatted', ACCESS='direct', RECL=total_elements*bytes_per_real, &
           IOSTAT=iostat, IOMSG=iomsg)
      
      IF (iostat /= 0) THEN
          WRITE(*,'(A)') 'ERROR opening weights file: ', TRIM(iomsg)
          DEALLOCATE(buffer)
          STOP 1
      END IF
      
      READ(iunit, REC=1, IOSTAT=iostat, IOMSG=iomsg) buffer
      
      IF (iostat /= 0) THEN
          WRITE(*,'(A)') 'ERROR reading weights: ', TRIM(iomsg)
          DEALLOCATE(buffer)
          CLOSE(iunit)
          STOP 1
      END IF
      
      ! Reshape the 1D buffer into 2D array in column-major order
      weights(1:input_neurons, 1:output_neurons) = &
          RESHAPE(buffer, [input_neurons, output_neurons])
      
      CLOSE(iunit)
      DEALLOCATE(buffer)
      
    END SUBROUTINE read_nn_weights

    SUBROUTINE read_nn_biases(filename_biases, biases, output_neurons)
      CHARACTER(LEN=*), INTENT(IN) :: filename_biases
      REAL(kind=WP), INTENT(INOUT) :: biases(:)
      INTEGER, INTENT(IN) :: output_neurons
      INTEGER :: iunit, iostat
      CHARACTER(512) :: iomsg
      REAL(kind=WP), ALLOCATABLE :: buffer(:)
      INTEGER :: bytes_per_real
      
      ALLOCATE(buffer(output_neurons))
      
      ! Calculate record length in bytes based on actual REAL kind
      bytes_per_real = STORAGE_SIZE(buffer(1)) / 8  ! bits to bytes
      
      ! Use direct access with record length = total bytes
      OPEN(NEWUNIT=iunit, FILE=filename_biases, STATUS='old', &
           FORM='unformatted', ACCESS='direct', RECL=output_neurons*bytes_per_real, &
           IOSTAT=iostat, IOMSG=iomsg)
      
      IF (iostat /= 0) THEN
          WRITE(*,'(A)') 'ERROR opening biases file: ', TRIM(iomsg)
          DEALLOCATE(buffer)
          STOP 1
      END IF
      
      READ(iunit, REC=1, IOSTAT=iostat, IOMSG=iomsg) buffer
      
      IF (iostat /= 0) THEN
          WRITE(*,'(A)') 'ERROR reading biases: ', TRIM(iomsg)
          DEALLOCATE(buffer)
          CLOSE(iunit)
          STOP 1
      END IF
      
      biases(1:output_neurons) = buffer
      
      CLOSE(iunit)
      DEALLOCATE(buffer)
      
    END SUBROUTINE read_nn_biases

    SUBROUTINE read_nn_activation(filename, activation)
      ! Reads activations. Routine is called layerwise, so this is only for one layer. (single string)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      CHARACTER(LEN=*), INTENT(OUT) :: activation
      INTEGER :: iunit, iostat, i
      CHARACTER(512) :: iomsg
      CHARACTER(LEN=16) :: temp_activation, valid_act
      CHARACTER(LEN=*), PARAMETER :: valid_activations(2) = ['relu   ', 'id     ']
      LOGICAL :: is_valid
      
      OPEN(NEWUNIT=iunit, FILE=filename, STATUS='old', ACTION='read', &
           IOSTAT=iostat, IOMSG=iomsg)
      
      IF (iostat /= 0) THEN
          WRITE(*,'(A)') 'ERROR opening activation file: ', TRIM(iomsg)
          STOP 1
      END IF
      
      READ(iunit, '(A)', IOSTAT=iostat) temp_activation
      CLOSE(iunit)
      
      IF (iostat /= 0) THEN
          WRITE(*,'(A)') 'ERROR reading activation file: ', TRIM(filename)
          STOP 1
      END IF
      
      ! Trim and validate
      activation = TRIM(ADJUSTL(temp_activation))
      
      ! Check if valid activation function
      is_valid = .false.
      DO i = 1, SIZE(valid_activations)
          IF (TRIM(activation) == TRIM(valid_activations(i))) THEN
              is_valid = .true.
              EXIT
          END IF
      END DO
      
      IF (.NOT. is_valid) THEN
          WRITE(*,'(A)') 'ERROR: Unknown activation function: ', TRIM(activation)
          WRITE(*,'(A)') 'Valid options: relu, sigmoid, id'
          STOP 1
      END IF

    END SUBROUTINE read_nn_activation

    SUBROUTINE forward_pass_single(inputs, layer_idx, nn, outputs)
        ! Forward pass through a single layer with activation function
        ! inputs:   input vector to the layer
        ! layer_idx: which layer (1 to nlayers)
        ! nn:       neural network data structure
        ! outputs:  result after matrix mult and activation
        REAL(kind=WP), INTENT(IN) :: inputs(:)
        INTEGER, INTENT(IN) :: layer_idx
        TYPE(t_neural_net), POINTER, INTENT(IN) :: nn
        REAL(kind=WP), ALLOCATABLE, INTENT(OUT) :: outputs(:)
        
        INTEGER :: in_size, out_size
        
        IF (.NOT. ASSOCIATED(nn)) THEN
            WRITE(*,*) 'ERROR: NN not initialized in forward_pass_single'
            STOP 1
        END IF
        
        IF (layer_idx < 1 .OR. layer_idx > nn%nlayers) THEN
            WRITE(*,*) 'ERROR: Invalid layer index ' , layer_idx, ' for NN with ', nn%nlayers, ' layers'
            STOP 1
        END IF
        
        in_size = nn%layer_sizes(layer_idx)
        out_size = nn%layer_sizes(layer_idx + 1)
        
        IF (SIZE(inputs) /= in_size) THEN
            WRITE(*,*) 'ERROR: Input size mismatch. Expected ', in_size, ' but got ', SIZE(inputs)
            STOP 1
        END IF
        
        ALLOCATE(outputs(out_size))
        
        ! Matrix multiplication: output = W^T . input + bias
        ! Weights are stored as (in_size, out_size)
        outputs = MATMUL(TRANSPOSE(nn%weights(layer_idx, 1:in_size, 1:out_size)), inputs) &
                + nn%biases(layer_idx, 1:out_size)
        
        ! Apply activation function
        SELECT CASE (TRIM(nn%activations(layer_idx)))
            CASE ('relu')
                outputs = MAX(outputs, 0.0_WP)
            CASE ('id')
                ! Identity activation - no operation
            CASE DEFAULT
                WRITE(*,'(A)') 'ERROR: Unknown activation function: ' // TRIM(nn%activations(layer_idx))
                STOP 1
        END SELECT
        
    END SUBROUTINE forward_pass_single

    SUBROUTINE forward_pass_full(input, nn, output)
        ! Full forward pass through entire neural network
        ! Uses the global neural network data structure
        REAL(kind=WP), INTENT(IN) :: input(:)
        TYPE(t_neural_net), POINTER, INTENT(IN) :: nn
        REAL(kind=WP), ALLOCATABLE, INTENT(OUT) :: output(:)
        
        INTEGER :: i, nlayers
        REAL(kind=WP), ALLOCATABLE :: layer_input(:), layer_output(:)
        
        IF (.NOT. ASSOCIATED(nn)) THEN
            WRITE(*,*) 'ERROR: NN not initialized in forward_pass_full'
            STOP 1
        END IF
        
        nlayers = nn%nlayers
        ALLOCATE(layer_input(SIZE(input)))
        layer_input = input
        
        ! Loop through each layer
        DO i = 1, nlayers
            CALL forward_pass_single(layer_input, i, nn, layer_output)
            DEALLOCATE(layer_input)
            ALLOCATE(layer_input(SIZE(layer_output)))
            layer_input = layer_output
            DEALLOCATE(layer_output)
        END DO
        
        ! Final output
        ALLOCATE(output(SIZE(layer_input)))
        output = layer_input
        DEALLOCATE(layer_input)
        
    END SUBROUTINE forward_pass_full

    ! SUBROUTINE nn_inference_for_node_layer(node, layer, dynamics, tracers, mesh, partit, &
    !                                         curl_u_nn, ld_baroc, N2_nn, slope_x_nn, slope_y_nn, &
    !                                         nn, nn_output)
    !     ! Extract inputs and run NN inference for a single (node, vertical layer) pair
    !     ! Returns NN output (e.g., flux components)
    !     INTEGER, INTENT(IN) :: node, layer  
    !     TYPE(t_dyn), INTENT(IN), TARGET :: dynamics
    !     TYPE(t_tracer), INTENT(IN), TARGET :: tracers
    !     TYPE(t_mesh), INTENT(IN), TARGET :: mesh
    !     TYPE(t_partit), INTENT(IN), TARGET :: partit
    !     TYPE(t_neural_net), POINTER, INTENT(IN) :: nn
    !     REAL(kind=WP), ALLOCATABLE, INTENT(OUT) :: nn_output(:)
        
    !     ! Optional pre-computed fields (if not provided, will be accessed from structures)
    !     REAL(kind=WP), INTENT(IN), OPTIONAL :: curl_u_nn(:,:), ld_baroc(:), N2_nn(:,:), &
    !                                             slope_x_nn(:,:), slope_y_nn(:,:)
        
    !     REAL(kind=WP), ALLOCATABLE :: nn_input(:)
    !     INTEGER :: n_features
        
    !     IF (.NOT. ASSOCIATED(nn)) RETURN
        
    !     n_features = nn%layer_sizes(1)  ! Input layer size
    !     ALLOCATE(nn_input(n_features))
        
    !     ! Extract input features from FESOM data at (node, layer)
    !     CALL extract_nn_features(node, layer, dynamics, tracers, mesh, partit, &
    !                              curl_u_nn, ld_baroc, N2_nn, slope_x_nn, slope_y_nn, &
    !                              nn_input)
        
    !     ! Run full NN forward pass
    !     CALL forward_pass_full(nn_input, nn, nn_output)
        
    !     DEALLOCATE(nn_input)
        
    ! END SUBROUTINE nn_inference_for_node_layer

    ! ============================================================================================
    ! EXTRACT INPUT FEATURE VALUES DURING RUN
    ! ============================================================================================
    SUBROUTINE extract_nn_features(nod2, nz1, dynamics, tracers, mesh, partit, &
        curl_u_nn, ld_baroc, N2_nn, slope_x_nn, slope_y_nn, nn_input)
        ! Extract and normalize NN input features for one (nod2, nz1) pair
        ! Feature order (MUST match training data order):
        !   1-7:   curl_u_nn (node + 6 neighbors)
        !   8:     ld_baroc1 (rosb, 2D only)
        !   9-15:  N2 (node + 6 neighbors)
        !   16-22: slope_x (node + 6 neighbors)
        !   23-29: slope_y (node + 6 neighbors)
        !   30-36: temperature (node + 6 neighbors)
        !   37-43: unod (node + 6 neighbors)
        !   44-50: vnod (node + 6 neighbors)
        ! Total: 50 features

        ! WATCH OUT: For neighbor value access, use LOCAL node indices and only if != 0
        
        INTEGER, INTENT(IN) :: nod2, nz1
        TYPE(t_dyn), INTENT(IN), TARGET :: dynamics
        TYPE(t_tracer), INTENT(IN), TARGET :: tracers
        TYPE(t_mesh), INTENT(IN), TARGET :: mesh
        TYPE(t_partit), INTENT(IN), TARGET :: partit
        REAL(KIND=WP), INTENT(IN) :: curl_u_nn(:,:), ld_baroc(:), N2_nn(:,:), &
                                                slope_x_nn(:,:), slope_y_nn(:,:)
        REAL(kind=WP), INTENT(OUT) :: nn_input(:)
        
    !     ! Cached normalization parameters (load once per simulation)
    !     REAL(kind=WP), ALLOCATABLE, SAVE :: means(:), stds(:)
    !     LOGICAL, SAVE :: params_loaded = .FALSE.
        
    !     ! Local variables
    !     INTEGER :: idx, nb_idx, nb_node, i, j
    !     REAL(kind=WP) :: raw_value
    !     INTEGER :: neighbors(6)
        
    !     ! Load normalization parameters once per simulation
    !     IF (.NOT. params_loaded) THEN
    !         CALL load_nn_normalization_params(means, stds)
    !         params_loaded = .TRUE.
    !     END IF
        
    !     IF (SIZE(nn_input) /= SIZE(means)) THEN
    !         WRITE(*,*) 'ERROR: NN input size mismatch. Expected ', SIZE(means), ' got ', SIZE(nn_input)
    !         STOP 1
    !     END IF
        
    !     ! Get neighbors of current node
    !     CALL get_node_neighbors(node, mesh, neighbors)
        
    !     idx = 1
        
    !     ! ===== 1-7: curl_u_nn (node + 6 neighbors) =====
    !     IF (PRESENT(curl_u_nn)) THEN
    !         raw_value = curl_u_nn(layer, node)
    !     ELSE
    !         raw_value = 0.0_WP  ! TODO: Compute or fetch from available structures
    !     END IF
    !     nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !     idx = idx + 1
        
    !     DO nb_idx = 1, 6
    !         IF (neighbors(nb_idx) > 0) THEN
    !             IF (PRESENT(curl_u_nn)) THEN
    !                 raw_value = curl_u_nn(layer, neighbors(nb_idx))
    !             ELSE
    !                 raw_value = 0.0_WP
    !             END IF
    !         ELSE
    !             raw_value = means(idx)  ! Pad with zero (normalize to 0)
    !         END IF
    !         nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !         idx = idx + 1
    !     END DO
        
    !     ! ===== 8: ld_baroc1 (rosb, 2D only, no neighbors) =====
    !     IF (PRESENT(ld_baroc)) THEN
    !         raw_value = ld_baroc(node)
    !     ELSE
    !         raw_value = 0.0_WP  ! TODO: Compute or fetch from available structures
    !     END IF
    !     nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !     idx = idx + 1
        
    !     ! ===== 9-15: N2 (node + 6 neighbors) =====
    !     IF (PRESENT(N2_nn)) THEN
    !         raw_value = N2_nn(layer, node)
    !     ELSE
    !         raw_value = 0.0_WP  ! TODO: Use bvfreq from pressure/stability module
    !     END IF
    !     nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !     idx = idx + 1
        
    !     DO nb_idx = 1, 6
    !         IF (neighbors(nb_idx) > 0) THEN
    !             IF (PRESENT(N2_nn)) THEN
    !                 raw_value = N2_nn(layer, neighbors(nb_idx))
    !             ELSE
    !                 raw_value = 0.0_WP
    !             END IF
    !         ELSE
    !             raw_value = means(idx)
    !         END IF
    !         nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !         idx = idx + 1
    !     END DO
        
    !     ! ===== 16-22: slope_x (node + 6 neighbors) =====
    !     IF (PRESENT(slope_x_nn)) THEN
    !         raw_value = slope_x_nn(layer, node)
    !     ELSE
    !         raw_value = 0.0_WP  ! TODO: Use neutral_slope(1,:,:)
    !     END IF
    !     nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !     idx = idx + 1
        
    !     DO nb_idx = 1, 6
    !         IF (neighbors(nb_idx) > 0) THEN
    !             IF (PRESENT(slope_x_nn)) THEN
    !                 raw_value = slope_x_nn(layer, neighbors(nb_idx))
    !             ELSE
    !                 raw_value = 0.0_WP
    !             END IF
    !         ELSE
    !             raw_value = means(idx)
    !         END IF
    !         nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !         idx = idx + 1
    !     END DO
        
    !     ! ===== 23-29: slope_y (node + 6 neighbors) =====
    !     IF (PRESENT(slope_y_nn)) THEN
    !         raw_value = slope_y_nn(layer, node)
    !     ELSE
    !         raw_value = 0.0_WP  ! TODO: Use neutral_slope(2,:,:)
    !     END IF
    !     nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !     idx = idx + 1
        
    !     DO nb_idx = 1, 6
    !         IF (neighbors(nb_idx) > 0) THEN
    !             IF (PRESENT(slope_y_nn)) THEN
    !                 raw_value = slope_y_nn(layer, neighbors(nb_idx))
    !             ELSE
    !                 raw_value = 0.0_WP
    !             END IF
    !         ELSE
    !             raw_value = means(idx)
    !         END IF
    !         nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !         idx = idx + 1
    !     END DO
        
    !     ! ===== 30-36: temperature (node + 6 neighbors) =====
    !     raw_value = tracers%data(1)%values(layer, node)
    !     nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !     idx = idx + 1
        
    !     DO nb_idx = 1, 6
    !         IF (neighbors(nb_idx) > 0) THEN
    !             raw_value = tracers%data(1)%values(layer, neighbors(nb_idx))
    !         ELSE
    !             raw_value = means(idx)
    !         END IF
    !         nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !         idx = idx + 1
    !     END DO
        
    !     ! ===== 37-43: unod (node + 6 neighbors) =====
    !     raw_value = dynamics%uvnode(1, layer, node)
    !     nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !     idx = idx + 1
        
    !     DO nb_idx = 1, 6
    !         IF (neighbors(nb_idx) > 0) THEN
    !             raw_value = dynamics%uvnode(1, layer, neighbors(nb_idx))
    !         ELSE
    !             raw_value = means(idx)
    !         END IF
    !         nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !         idx = idx + 1
    !     END DO
        
    !     ! ===== 44-50: vnod (node + 6 neighbors) =====
    !     raw_value = dynamics%uvnode(2, layer, node)
    !     nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !     idx = idx + 1
        
    !     DO nb_idx = 1, 6
    !         IF (neighbors(nb_idx) > 0) THEN
    !             raw_value = dynamics%uvnode(2, layer, neighbors(nb_idx))
    !         ELSE
    !             raw_value = means(idx)
    !         END IF
    !         nn_input(idx) = (raw_value - means(idx)) / stds(idx)
    !         idx = idx + 1
    !     END DO
        
    END SUBROUTINE extract_nn_features

    SUBROUTINE load_nn_normalization_params(means, stds)
        ! Load normalization parameters (means and stds) from disk
        ! Called once per simulation, cached via SAVE in extract_nn_features
        ! Total features: 50 (curl_u_nn + ld_baroc1 + N2 + slope_x + slope_y + temp + unod + vnod, 
        !                      all with 6 neighbors each except ld_baroc1)
        REAL(kind=WP), ALLOCATABLE, INTENT(OUT) :: means(:), stds(:)
        
        CHARACTER(LEN=32), PARAMETER :: base_path = './neuralnet_params'
        INTEGER :: iunit, iostat
        CHARACTER(256) :: iomsg
        INTEGER :: n_features
        LOGICAL :: normalization_params_loaded
        
        ! Total number of features: 50
        n_features = 50
        
        ALLOCATE(means(n_features))
        ALLOCATE(stds(n_features))
        
        ! Load normalization constants from individual files
        ! Expected files structure (confirm with your training setup):
        ! /albedo/home/rjuhrban/fesom2/src/neuralnet_params/
        !   mean/curl_u_nn.bin, mean/ld_baroc1.bin, mean/N2.bin, etc.
        !   std/curl_u_nn.bin, std/ld_baroc1.bin, std/N2.bin, etc.
        !
        ! If neighbors are stored separately, adjust accordingly
        
        ! Initialize with placeholder values (will fail obviously if not replaced)
        means = 0.0_WP
        stds = 1.0_WP

        normalization_params_loaded = .FALSE.

        ! Do normalization with float64 precision. But WP in FESOM is already 8, so nothing to do here
        
        ! TODO: Implement actual file reading
        ! CALL read_nn_feature_means_stds(base_path, means, stds)
        
        ! Placeholder warning
        WRITE(*,'(A)') 'WARNING: Normalization parameters not loaded (using defaults: means=0, stds=1)'
        
    END SUBROUTINE load_nn_normalization_params

    ! SUBROUTINE apply_nn_corrections_tracer(dynamics, tracers, partit, mesh, tr_num, &
    !                                       curl_u_nn, ld_baroc, N2_nn, slope_x_nn, slope_y_nn)
    !     ! Apply NN-based corrections to tracer tendencies
    !     ! Call this after advection but with access to del_ttf for adding flux divergences
    !     ! This is where unresolved flux contributions from the NN are added
    !     !
    !     ! Optional arguments (pre-computed fields for efficiency):
    !     !   - curl_u_nn: Curl of velocity from diag_curl_vel3, shape(nlev, nod2D)
    !     !   - ld_baroc: Rossby radius from oce_fer_gm, shape(nod2D)
    !     !   - N2_nn: Buoyancy frequency BEFORE smoothing from oce_ale_pressure_bv, shape(nlev, nod2D)
    !     !   - slope_x_nn: Neutral slope x-component from oce_ale_pressure_bv, shape(nlev, nod2D)
    !     !   - slope_y_nn: Neutral slope y-component from oce_ale_pressure_bv, shape(nlev, nod2D)
        
    !     TYPE(t_dyn), INTENT(INOUT), TARGET :: dynamics
    !     TYPE(t_tracer), INTENT(INOUT), TARGET :: tracers
    !     TYPE(t_partit), INTENT(INOUT), TARGET :: partit
    !     TYPE(t_mesh), INTENT(IN), TARGET :: mesh
    !     INTEGER, INTENT(IN) :: tr_num
    !     REAL(kind=WP), INTENT(IN), OPTIONAL :: curl_u_nn(:,:), ld_baroc(:), N2_nn(:,:), &
    !                                             slope_x_nn(:,:), slope_y_nn(:,:)
        
    !     TYPE(t_neural_net), POINTER :: nn
    !     REAL(kind=WP), ALLOCATABLE :: nn_output(:)
    !     INTEGER :: n, nz, nzmax, nzmin
        
    !     nn => get_neural_net()
    !     IF (.NOT. ASSOCIATED(nn)) RETURN
        
    !     ! Parallelize over nodes and vertical levels
    !     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmax, nzmin, nn_output)
    !     DO n = 1, myDim_nod2D
    !         nzmin = ulevels_nod2D(n)
    !         nzmax = kmt(n)
            
    !         DO nz = nzmin, nzmax
    !             ! Run NN inference for this (node, level)
    !             ! CALL nn_inference_for_node_layer(n, nz, dynamics, tracers, mesh, partit, &
    !             !                                  curl_u_nn, ld_baroc, N2_nn, slope_x_nn, slope_y_nn, &
    !             !                                  nn, nn_output)
                
    !             ! TODO: Add NN output to flux divergence or tracer tendency
    !             ! Currently: nn_output contains NN predictions (likely flux divergences or flux components)
    !             ! To be integrated: 
    !             !   del_ttf(nz, n) = del_ttf(nz, n) + nn_output_contribution
                
    !             IF (ALLOCATED(nn_output)) DEALLOCATE(nn_output)
    !         END DO
    !     END DO
    !     !$OMP END PARALLEL DO
        
    ! END SUBROUTINE apply_nn_corrections_tracer

END MODULE MOD_NEURALNET