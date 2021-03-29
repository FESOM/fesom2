module openacc_params

    integer, parameter :: z_vector_length = 64 ! vector length of vertical loops

    integer, parameter :: stream_hor_adv_tra = 3
    integer, parameter :: stream_ver_adv_tra = 4
    integer, parameter :: stream_hnode_update = 5
    integer, parameter :: stream_redi = 6

end module
