MODULE o_LU_ARRAYS

! LU arrays ...

 USE o_PARAM
 IMPLICIT NONE
 save
  real(kind=WP), allocatable    :: sdbt_x_lu(:,:,:)
  real(kind=WP), allocatable    :: sdbt_y_lu(:,:,:)
  real(kind=WP), allocatable    :: sdbt_z_lu(:,:,:)

  real(kind=WP), allocatable    :: us_lu(:,:,:)
  real(kind=WP), allocatable    :: vs_lu(:,:,:)
  real(kind=WP), allocatable    :: ws_lu(:,:,:)

END MODULE o_LU_ARRAYS
