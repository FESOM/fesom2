module openacc_params

    integer, parameter :: z_vector_length = 64 ! vector length of vertical loops

    integer, parameter :: stream_hor_adv_tra = 3
    integer, parameter :: stream_ver_adv_tra = 4
    integer, parameter :: stream_hnode_update = 5
    integer, parameter :: stream_hor_diff_tra = 6
    integer, parameter :: stream_ver_diff_tra = 7

end module
