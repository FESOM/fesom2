subroutine create_noise(mesh)
 use MOD_MESH
 use o_PARAM
 use o_ARRAYS
 use g_PARSUP

 implicit none
 integer,parameter :: nobs = 20
 integer,parameter :: patchsize = 9
 real              :: rnd_ary(nobs)
 integer           :: rnd_obs_id(nobs)
 integer           :: nz, elem, i, elem_size
 real              :: mean_u, mean_v, scaling, tau, dx, dy, dz, dt
 type(t_mesh), intent(in)              , target :: mesh

 
 dx=1.
 dy=1.
 dz=1.
 dt=1.
 
 scaling = SQRT(real(patchsize)**(-2./3.))/SQRT(real(nObs) - 1.)
 
 elem_size=myDim_elem2D+eDim_elem2D

   Do nz=1, 10!nl-1 ! TODOLU where is nl defined ...?
    
    mean_u=sum(UV(1,nz,:))/real(elem_size)
    mean_v=sum(UV(2,nz,:))/real(elem_size)
    
    ! set a_.. to 0
    
    Do elem=1, myDim_elem2D
      
       call random_number(rnd_ary)
       rnd_obs_id = 1 + FLOOR((patchsize+1-1)*rnd_ary)  ! random integer in [1,...,patchsize]

       Do i=1, nobs  ! fill matrix W with u and v random observation for each layer
           
           lu_mat_W(elem,i)              = UV(1,nz,uv_neighbour_set(nz,elem,rnd_obs_id(i)))-mean_u
           lu_mat_W(elem + myDim_elem2D,i)  = UV(2,nz,uv_neighbour_set(nz,elem,rnd_obs_id(i)))-mean_v

!compute the associated w velocity
           
           lu_mat_W(elem + 2*myDim_elem2D,i)  = lu_mat_W(elem,i)/2 + lu_mat_W(elem + myDim_elem2D,i)/2 ! TODOLU just a dummy here    
           
       end do
       
    end do


   ! compute SVD of W giving matirx U and vector S
    
      call SVD(lu_mat_W,lu_mat_U,lu_mat_S,myDim_elem2D,20)
    
   ! rescale S
       
       lu_mat_S = lu_mat_S*scaling
       
   ! compute variance tensor a 
       Do i=1, nobs
           lu_a_xx(:,nz) =lu_a_xx(:,nz) + (lu_mat_S(i)**2)*(lu_mat_U(1:myDim_elem2D,i)**2)
           lu_a_yy(:,nz) =lu_a_yy(:,nz) + (lu_mat_S(i)**2)*(lu_mat_U(myDim_elem2D+1:myDim_elem2D,i)**2)
           lu_a_zz(:,nz) =lu_a_zz(:,nz) + (lu_mat_S(i)**2)*(lu_mat_U(2*myDim_elem2D+1:3*myDim_elem2D,i)**2)
           lu_a_xy(:,nz) =lu_a_xy(:,nz) + (lu_mat_S(i)**2)*(lu_mat_U(1:myDim_elem2D,i)*lu_mat_U(myDim_elem2D+1:2*myDim_elem2D,i))
           lu_a_xz(:,nz) =lu_a_xz(:,nz) + (lu_mat_S(i)**2)*(lu_mat_U(1:myDim_elem2D,i)*lu_mat_U(2*myDim_elem2D+1:3*myDim_elem2D,i))
           lu_a_yz(:,nz) =lu_a_yz(:,nz) + (lu_mat_S(i)**2)*(lu_mat_U(myDim_elem2D+1:myDim_elem2D,i)*lu_mat_U(2*myDim_elem2D+1:3*myDim_elem2D,i))
       end do
   ! rescale a
  
       tau=min(  ( dx/sqrt(maxval(lu_a_xx(:,nz))) )**(2./3.)*dt**(1./3.), ( dy/sqrt(maxval(lu_a_yy(:,nz))) )**(2./3.)*dt**(1./3.), ( dz/sqrt(maxval(lu_a_zz(:,nz))) )**(2./3.)*dt**(1./3.)  )
    
       lu_a_xx(:,nz) =tau*lu_a_xx(:,nz)
       lu_a_yy(:,nz) =tau*lu_a_yy(:,nz)
       lu_a_zz(:,nz) =tau*lu_a_zz(:,nz)
       lu_a_xy(:,nz) =tau*lu_a_xy(:,nz)
       lu_a_xz(:,nz) =tau*lu_a_xz(:,nz)
       lu_a_yz(:,nz) =tau*lu_a_yz(:,nz)
       
   ! compute sdbt 
       call random_number(rnd_ary)
       lu_mat_S=lu_mat_S*(2.*rnd_ary-1.)
       
       Do elem=1, myDim_elem2D
           lu_sdbt_UV(1,nz,elem) = sum(lu_mat_U(elem,:)*lu_mat_S)
           lu_sdbt_UV(2,nz,elem) = sum(lu_mat_U(elem + myDim_elem2D,:)*lu_mat_S)
           !lu_sdbt_w(nz,node)   = sum(lu_mat_U(elem + 2* myDim_elem2D,:)*lu_mat_S) ! TODOLU needs to be mapped to the nodes
       end do    
       
   ! project on spherical coordinates ???
    
    
   end do 

end subroutine create_noise


subroutine SVD(A,U,S,M,N)

! Program computes the matrix singular value decomposition. 
! Using Lapack library.

 DOUBLE PRECISION A(M,N),U(M,M),S(N)
 DOUBLE PRECISION,ALLOCATABLE :: WORK(:)
 INTEGER LDA,M,N,LWORK,INFO

 LDA=M
 LDU=M

 LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))

 ALLOCATE(work(lwork))

 CALL DGESVD('A','N', M, N, A, LDA, S, U, LDU, 0,1, WORK, LWORK, INFO )
        
end subroutine SVD
