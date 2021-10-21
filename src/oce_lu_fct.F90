subroutine create_noise(mesh)
 use MOD_MESH
 use o_PARAM
 use o_ARRAYS

 implicit none
 integer         :: nobs = 20
 integer         :: patchsize = 9
 integer         :: nz, elem, i
! real            :: mat_W(3*myDim_elem2D,nobs)
 type(t_mesh), intent(in)              , target :: mesh


 

   ! fill matrix W with u and v random observation for each layer
              
   !compute the associated w velocity

   ! compute SVD of W giving matirx U and vector S

   ! rescale S

   ! compute variance tensor a 

   ! rescale a

   ! compute sdbt 

   ! project on spherical coordinates ???

 

end subroutine create_noise
