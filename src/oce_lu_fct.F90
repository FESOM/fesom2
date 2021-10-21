subroutine create_noise(mesh)
 use MOD_MESH
 use o_PARAM

 implicit none
 integer         :: nobs = 20
 integer         :: patchsize = 9
 integer         :: nz, elem, i
 real            :: mat_W(3*myDim_elem2D,nobs)
 type(t_mesh), intent(in)              , target :: mesh


DO nz=nzmin, nzmax-1
    
        mean_u=sum(UV(1,nz,:))/real(myDim_elem2D)
        mean_v=sum(UV(2,nz,:))/real(myDim_elem2D)
        
        DO elem=1, myDim_elem2D


         ! step 1:  fill matrix W with u and v random observation for each layer
                                        
            call random_number(rnd_ary)
            rnd_obs_id = 1 + FLOOR((patchsize+1-1)*rnd_ary)  ! random integer in [1,...,patchsize]

            DO i=1, nobs
                mat_W(elem,i) = UV(1,nz,uv_neighbour_set(nz,elem,rnd_obs_id(i)))-mean_u
                mat_W(elem+myDim_elem2D,i) = UV(2,nz,uv_neighbour_set(nz,elem,rnd_obs_id(i)))-mean_v
            END
        
        END DO

   ! compute the associated w velocity

   ! compute SVD of W giving matirx U and vector S

   ! rescale S

   ! compute variance tensor a 

   ! rescale a

   ! compute sdbt 

   ! project on spherical coordinates ???

END DO

end subroutine create_noise
