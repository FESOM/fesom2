! A set of routines to implement quadratic reconstruction. Additional 
! global arrays that are necessary for implementation are announced 
! in module mesh.
!
! They are allocated and filled in quadratic_reconstruction_init
! 
! solve_tracer_second_order_rec is the advection routine. The rest are 
! auxiliary routines called either by quadratic_reconstruction_init or 
! the advection routine. 
! 
! The procedure (in its horizontal part) corresponds to  

! Ollivier-Gooh, C., Van Altena, M., 2002. A high-order-accurate
! unstructured mesh finite-volume scheme for the advection/diffusion
! equation. J. Comput. Phys. 181, 729--752.

! Ouvrard, H., Kozubskaya, T., Abalakin, I., Koobus, B., Dervieux, A., 2009. 
! Advective vertex-centered reconstruction scheme on unstructured meshes. INRIA,
! Rapport de recherche 7033.
! 
! I use linear interpolation where the number of neighbors is insufficient.
! sergey.danilov@awi.de, 2012
!Contains:
!	quadratic_reconstruction_init
!	quadratic_reconstruction
!	linear_reconstruction
!	compute_moments
!	quad_int
!	lin_int
!	FINDInv
!	exchange_nod_quad_int
!	adv_tracer_second_order_rec
!===========================================================================
SUBROUTINE quadratic_reconstruction_init
! allocate and fill in arrays needed for quadratic interpolation
! They include
! quad_int_mat(mn,6,myDim_nod2D)    Matrices (for each of nodes) 
! quad_int_coef(6, nl-1, nod2D)     Place to store interpolation coefficients
! nlevels_nod2D_min(myDim_nod2D)    Full levels
! nn_num(myDim_nod2D),nn_pos(mn,myDim_nod2D) Neighboring nodes. 
!                                            We have encountered them in ssh_stiff,
!                                            but as temporary arrays.        
! 

USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
integer    :: n, k, mn, n1, n2
real(kind=WP), allocatable  :: aux(:,:)
real(kind=WP)               :: mcos

! =======
! a. Max number of neighbors 
! =======
  
  mn=0
  DO n=1, myDim_nod2D
  k=SSH_stiff%rowptr(n+1)-SSH_stiff%rowptr(n)
  if(k>mn) mn=k
  END DO
 
  allocate(quad_int_mat(mn,6,myDim_nod2D))
  allocate(quad_int_coef(6, nl-1, nod2D))
  allocate(nlevels_nod2D_min(myDim_nod2D))
  allocate(nn_num(myDim_nod2D),nn_pos(mn,myDim_nod2D))
  quad_int_coef=0.0_8

  DO n=1,myDim_nod2d
     nn_num(n)=1
     nn_pos(1,n)=n
  end do   
  Do n=1, myDim_edge2D
     n1=edges(1,n)
     n2=edges(2,n)
     if(n1<=myDim_nod2D) then
     nn_pos(nn_num(n1)+1,n1)=n2
     nn_num(n1)=nn_num(n1)+1
     end if
     if(n2<=myDim_nod2D) then
     nn_pos(nn_num(n2)+1,n2)=n1
     nn_num(n2)=nn_num(n2)+1
     end if
  END DO 

! =======
! b. Compute moments
! =======
  allocate(aux(5,myDim_nod2D+eDim_nod2D))

  call compute_moments(aux)

  
! =======
! c. Fill matrices
! =======
  
  
  DO n=1, myDim_nod2D
                               !! n=myList_nod2D(m)
  k=SSH_stiff%rowptr(n+1)-SSH_stiff%rowptr(n)
  if(k>=7) then
  call quad_int(n, aux)
  end if
  END DO
  
  DO n=1, myDim_nod2D
                               !! n=myList_nod2D(m)
  k=SSH_stiff%rowptr(n+1)-SSH_stiff%rowptr(n)
  if(k<7) then
  call lin_int(n, aux)
  end if
  END DO
  deallocate(aux)  
! =======
! Reduction to physical coordinates:
! (It would be better to keep distances in radian measure,
! but it is not done yet)
! =======
  DO n=1, myDim_nod2D
                               !! n=myList_nod2D(m)
  mcos=cos(coord_nod2D(2,n))
  if(cartesian) mcos=1.0_8
  quad_int_mat(:,2,n)=quad_int_mat(:,2,n)/mcos/r_earth
  quad_int_mat(:,3,n)=quad_int_mat(:,3,n)/r_earth
  quad_int_mat(:,4,n)=quad_int_mat(:,4,n)/(mcos*r_earth)**2
  quad_int_mat(:,5,n)=quad_int_mat(:,5,n)/mcos/r_earth/r_earth
  quad_int_mat(:,6,n)=quad_int_mat(:,6,n)/r_earth/r_earth	
  END DO
! ======
! d. levels with full control volumes
! ======
 
  
  DO n=1,myDim_nod2D
    nlevels_nod2D_min(n)=minval(nlevels_nod2D(nn_pos(1:nn_num(n),n))) 
  END DO
 
END SUBROUTINE quadratic_reconstruction_init
! ==========================================================================
SUBROUTINE quadratic_reconstruction(ttf)
! Find coefficients of quadratic reconstruction
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
IMPLICIT NONE
real(kind=WP)   :: ttf(nl-1,myDim_nod2D+eDIm_nod2D)
integer        :: n, nz, k, n_num, n_ind(20)
  ! =========
  ! compute coefficients of quadratic reconstruction
  ! =========
 DO n=1,myDim_nod2D
                      !! n=myList_nod2D(m)
    n_num=nn_num(n)   !! SSH_stiff%rowptr(n+1)-SSH_stiff%rowptr(n)  
    n_ind(1:n_num)=nn_pos(1:n_num, n)   !!SSH_stiff%colind(SSH_stiff%rowptr(n):SSH_stiff%rowptr(n+1)-1)
    DO nz=1,nlevels_nod2d_min(n)-1
       DO k=1,6
       quad_int_coef(k,nz,n)=sum(quad_int_mat(1:n_num,k,n)*ttf(nz,n_ind(1:n_num)))
       END DO
    END DO 
    DO nz=nlevels_nod2d_min(n), nlevels_nod2D(n)-1 
       quad_int_coef(1,nz,n)=ttf(nz,n)   ! will be used in linear rec. (biased)
    END DO     
 END DO
  ! ======== 
  ! Linear reconstruction is used for incomplete 
  ! control volumes. It has to be called before this 
  ! subroutine
  ! ======== 
  
  ! ========
  ! Exchange halo nodes
  ! ========
  call exchange_nod_quad_int(quad_int_coef)
  
END SUBROUTINE quadratic_reconstruction
!=======================================================================
SUBROUTINE linear_reconstruction (tt_xy)
! Biased linear reconstruction
! ttx, tty  elemental gradient of tracer 
! ttxnodes, ttynodes nodal gradients obtained by averaging
USE o_PARAM
USE o_MESH
!USE o_ARRAYS
USE g_PARSUP
IMPLICIT NONE
real(kind=WP)      :: tt_xy(2,nl-1,myDim_elem2D+eDim_elem2D)
integer            :: n, nz, elem, k
real(kind=WP)      :: tvol, tx, ty
  DO n=1, myDim_nod2D        
     DO nz=nlevels_nod2D_min(n), nlevels_nod2D(n)-1
        tvol=0.0
	tx=0.0
	ty=0.0
	DO k=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(k,n)
           if(nlevels(elem)-1<nz) cycle
           tvol=tvol+elem_area(elem)
           tx=tx+tt_xy(1,nz,elem)*elem_area(elem)
	   ty=ty+tt_xy(2,nz,elem)*elem_area(elem)
        END DO
	quad_int_coef(2,nz,n)=tx/tvol
	quad_int_coef(3,nz,n)=ty/tvol
     END DO
   END DO 
   ! quad_int_coef(1,nz, n) is filled in quadratic reconstruction routine
   
END SUBROUTINE linear_reconstruction
 !===========================================================================
SUBROUTINE compute_moments(aux)
! Compute  integrals from y^nx^m over control volume 
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
  use g_comm_auto
IMPLICIT NONE
real(kind=WP)      :: aux(5,myDim_nod2D+eDim_nod2D)
real(kind=WP)      :: aa, b(2), c(2)
integer            :: edge, el(2), ednodes(2), n
real(kind=WP), allocatable :: temp_array(:)
aux=0.0

 DO edge=1, myDim_edge2D    !!!+eDim_edge2D
    el=edge_tri(:,edge)
    ! ====
    ! left element (always exist)
    ! ====
    ednodes=edges(:,edge)
    c=coord_nod2D(:,ednodes(2))-coord_nod2D(:,ednodes(1))    
    call elem_center(el(1),b(1),b(2))
    b=b-coord_nod2D(:,ednodes(1))
	  if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
          if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
	  if(c(1)>cyclic_length/2.) c(1)=c(1)-cyclic_length
          if(c(1)<-cyclic_length/2.) c(1)=c(1)+cyclic_length
    c=c/2.0_8    ! we need distance to midedge
    
    aa=elem_area(el(1))/6.0_8	  
    aux(1,ednodes(1))=aux(1,ednodes(1))+aa*(b(1)+c(1))/3.0_8
    aux(2,ednodes(1))=aux(2,ednodes(1))+aa*(b(2)+c(2))/3.0_8
    aux(3,ednodes(1))=aux(3,ednodes(1))+aa*(b(1)*b(1)+c(1)*c(1)+b(1)*c(1))/6.0_8
    aux(4,ednodes(1))=aux(4,ednodes(1))+aa*(2*b(1)*b(2)+2*c(1)*c(2)+b(1)*c(2)+b(2)*c(1))/12.0_8
    aux(5,ednodes(1))=aux(5,ednodes(1))+aa*(b(2)*b(2)+c(2)*c(2)+b(2)*c(2))/6.0_8

    c=-c
    call elem_center(el(1),b(1),b(2))
    b=b-coord_nod2D(:,ednodes(2))
	  if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
          if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
    aa=elem_area(el(1))/6.0_8	  
    aux(1,ednodes(2))=aux(1,ednodes(2))+aa*(b(1)+c(1))/3.0_8
    aux(2,ednodes(2))=aux(2,ednodes(2))+aa*(b(2)+c(2))/3.0_8
    aux(3,ednodes(2))=aux(3,ednodes(2))+aa*(b(1)*b(1)+c(1)*c(1)+b(1)*c(1))/6.0_8
    aux(4,ednodes(2))=aux(4,ednodes(2))+aa*(2*b(1)*b(2)+2*c(1)*c(2)+b(1)*c(2)+b(2)*c(1))/12.0_8
    aux(5,ednodes(2))=aux(5,ednodes(2))+aa*(b(2)*b(2)+c(2)*c(2)+b(2)*c(2))/6.0_8
	
    ! ====
    ! right element 
    ! ====
    if(el(2)<1) cycle

    c=-c
    call elem_center(el(2),b(1),b(2))
    b=b-coord_nod2D(:,ednodes(1))
    if (b(1)> cyclic_length/2.)  b(1)=b(1)-cyclic_length
    if (b(1)<-cyclic_length/2.)  b(1)=b(1)+cyclic_length
    aa=elem_area(el(2))/6.0_8	  
    aux(1,ednodes(1))=aux(1,ednodes(1))+aa*(b(1)+c(1))/3.0_8
    aux(2,ednodes(1))=aux(2,ednodes(1))+aa*(b(2)+c(2))/3.0_8
    aux(3,ednodes(1))=aux(3,ednodes(1))+aa*(b(1)*b(1)+c(1)*c(1)+b(1)*c(1))/6.0_8
    aux(4,ednodes(1))=aux(4,ednodes(1))+aa*(2*b(1)*b(2)+2*c(1)*c(2)+b(1)*c(2)+b(2)*c(1))/12.0_8
    aux(5,ednodes(1))=aux(5,ednodes(1))+aa*(b(2)*b(2)+c(2)*c(2)+b(2)*c(2))/6.0_8

    c=-c
    call elem_center(el(2),b(1),b(2))
    b=b-coord_nod2D(:,ednodes(2))
	  if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
          if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
    aa=elem_area(el(2))/6.0_8	  
    aux(1,ednodes(2))=aux(1,ednodes(2))+aa*(b(1)+c(1))/3.0_8
    aux(2,ednodes(2))=aux(2,ednodes(2))+aa*(b(2)+c(2))/3.0_8
    aux(3,ednodes(2))=aux(3,ednodes(2))+aa*(b(1)*b(1)+c(1)*c(1)+b(1)*c(1))/6.0_8
    aux(4,ednodes(2))=aux(4,ednodes(2))+aa*(2*b(1)*b(2)+2*c(1)*c(2)+b(1)*c(2)+b(2)*c(1))/12.0_8
    aux(5,ednodes(2))=aux(5,ednodes(2))+aa*(b(2)*b(2)+c(2)*c(2)+b(2)*c(2))/6.0_8
  END DO 	

  DO n=1,myDim_nod2D
    aux(:,n)=aux(:,n)/area(1,n)
  END DO   
  ! ===
  ! Exchange halo information
  ! ===
  allocate(temp_array(myDim_nod2D+eDim_nod2D))
  DO n=1,5
  temp_array=aux(n,:)
  call exchange_nod(temp_array)
  aux(n,:)=temp_array
  END DO
  deallocate(temp_array)
end SUBROUTINE compute_moments    
! ================================================================================    
SUBROUTINE quad_int(n,aux)
! Fill in matrix quad_int_mat(:,:,n)
! n is the node number
! aux contains moments in radian measure
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
integer      :: k, q, n, n_num, n_ind(20)
real(kind=8) :: aux(5,myDim_nod2D+eDim_nod2D), b(2), weights(20)
real(kind=8), allocatable :: aux_mat(:,:), aux_r(:,:), moments(:,:),S_mat(:,:) 
real(kind=8), allocatable :: R_mat(:,:), R_aux(:,:)  
  
  n_num=nn_num(n)
  n_ind(1:n_num)=nn_pos(1:n_num,n)
  DO k=2,n_num
  b=coord_nod2D(:,n_ind(k))-coord_nod2D(:,n)
    if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
    if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
  
  weights(1:n_num)=1.0_8/sum(b*b)
  END DO
  weights(1)=1.0_8

allocate(aux_mat(7,7))
allocate(aux_r(7,n_num))
allocate(moments(6,n_num))

 ! =========
 ! fill moments
 ! The first column will be that of node n, the rest 
 ! neighbors of n
 ! ========= 
 DO k=1,n_num
    b=coord_nod2D(:, n_ind(k))-coord_nod2D(:, n)
    if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
    if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
    moments(1,k)=1.0_8
    moments(2,k)=b(1)+aux(1, n_ind(k))   
    moments(3,k)=b(2)+aux(2, n_ind(k))   
    moments(4,k)=b(1)*b(1)+2*b(1)*aux(1, n_ind(k))+aux(3, n_ind(k))
    moments(5,k)=b(1)*b(2)+b(1)*aux(2, n_ind(k))+b(2)*aux(1, n_ind(k))+aux(4, n_ind(k))      
    moments(6,k)=b(2)*b(2)+2*b(2)*aux(2, n_ind(k))+aux(5, n_ind(k))
 END DO
 ! ==========
 ! fill aux_mat
 ! T=u0+alpha x+ beta y +gamma x^2 +delta xy =epsilon y^2
 ! Minimization problem is solved
 ! ==========
   DO k=2,7
      DO q=2,7
      aux_mat(k,q)=sum(weights(2:n_num)*moments(k-1,2:n_num)*moments(q-1,2:n_num))
      END DO
   END DO
   aux_mat(1,2:7)=moments(:,1)
   aux_mat(1,1)=0.0
   aux_mat(2:7,1)=moments(:,1) 
 ! Test : aux_mat should be symmetric
 
 ! ==========
 ! fill aux_r
 ! (matrix of rhs operator) 
 ! ==========  
   aux_r=0.0
   DO k=2,n_num
      DO q=2,7
         aux_r(q,k)=weights(k)*moments(q-1,k)
      END DO
   END DO
   aux_r(1,1)=1.0_8
 ! ==========
 ! Remove the row with Lagrangian multiplier
 ! ==========   
  
   DO k=3,7  
      aux_mat(k,:)=aux_mat(k,:)-aux_mat(2,:)*moments(k-1,1)
      aux_r(k,:)=aux_r(k,:)-aux_r(2,:)*moments(k-1,1) 
   END DO
   ! aux_mat([1,3:7],2:7]) is now our S-matrix
   ! aux_r([1,3:7],:) is our rhs matrix
 ! ========== 
 ! Fill proper matrices 
 ! ========== 
   allocate(S_mat(6,6),R_aux(6,n_num),R_mat(6,n_num))  
   S_mat(1,:)=aux_mat(1,2:7)
   S_mat(2:6,:)=aux_mat(3:7,2:7)
   R_aux(1,:)=aux_r(1,:)
   R_aux(2:6,:)=aux_r(3:7,:)
   
   call FindInv(S_mat,R_aux,R_mat,6,n_num,k)
   DO k=1,6
      DO q=1,n_num 
         quad_int_mat(q,k,n)=R_mat(k,q)
      END DO
   END DO 
   Deallocate(R_mat,R_aux,S_mat,moments, aux_r,aux_mat)
END SUBROUTINE quad_int
!===========================================================================
SUBROUTINE lin_int(n,aux)
! Fill in matrix quad_int_mat(:,:,n), but only with linear 
! interpolation coefficients
! n is the node number
! aux contains moments in radian measure
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
integer      :: k, q, n, n_num, n_ind(20)
real(kind=8) :: aux(5,myDim_nod2D+eDim_nod2D), b(2), weights(20)
real(kind=8), allocatable :: aux_mat(:,:), aux_r(:,:), moments(:,:),S_mat(:,:) 
real(kind=8), allocatable :: R_mat(:,:), R_aux(:,:)  
  
  n_num=nn_num(n)
  n_ind(1:n_num)=nn_pos(1:n_num,n)
  DO k=2,n_num
  b=coord_nod2D(:,n_ind(k))-coord_nod2D(:,n)
    if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
    if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
  
  weights(1:n_num)=1.0_8/sum(b*b)
  END DO
  weights(1)=1.0_8
allocate(aux_mat(4,4))
allocate(aux_r(4,n_num))
allocate(moments(3,n_num))

 ! =========
 ! fill moments
 ! The first column will be that of node n, the rest 
 ! neighbors of n
 ! ========= 
 DO k=1,n_num
    b=coord_nod2D(:, n_ind(k))-coord_nod2D(:, n)
    if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
    if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
    moments(1,k)=1.0_8
    moments(2,k)=b(1)+aux(1, n_ind(k))   
    moments(3,k)=b(2)+aux(2, n_ind(k))   
 END DO
 ! ==========
 ! fill aux_mat
 ! T=u0+alpha x+ beta y 
 ! Minimization problem is solved
 ! ==========
   DO k=2,4
      DO q=2,4
      aux_mat(k,q)=sum(weights(2:n_num)*moments(k-1,2:n_num)*moments(q-1,2:n_num))
      END DO
   END DO
   aux_mat(1,2:4)=moments(:,1)
   aux_mat(1,1)=0.0
   aux_mat(2:4,1)=moments(:,1) 
 ! Test : aux_mat should be symmetric
 ! ==========
 ! fill aux_r
 ! (matrix of rhs operator) 
 ! ==========  
   aux_r=0.0
   DO k=2,n_num
      DO q=2,4
         aux_r(q,k)=weights(k)*moments(q-1,k)
      END DO
   END DO
   aux_r(1,1)=1.0_8
 ! ==========
 ! Remove the row with Lagrangian multiplier
 ! ==========   
  
   DO k=3,4  
      aux_mat(k,:)=aux_mat(k,:)-aux_mat(2,:)*moments(k-1,1)
      aux_r(k,:)=aux_r(k,:)-aux_r(2,:)*moments(k-1,1) 
   END DO
   ! aux_mat([1,3:4],2:4]) is now our S-matrix
   ! aux_r([1,3:4],:) is our rhs matrix
 ! ========== 
 ! Fill proper matrices 
 ! ========== 
   allocate(S_mat(3,3),R_aux(3,n_num),R_mat(3,n_num))  
   S_mat(1,:)=aux_mat(1,2:4)
   S_mat(2:3,:)=aux_mat(3:4,2:4)
   R_aux(1,:)=aux_r(1,:)
   R_aux(2:3,:)=aux_r(3:4,:)
   call FindInv(S_mat,R_aux,R_mat,3,n_num,k)
   DO k=1,3
      DO q=1,n_num 
         quad_int_mat(q,k,n)=R_mat(k,q)
      END DO
   END DO 
   DO k=4,6
      DO q=1,n_num 
         quad_int_mat(q,k,n)=0.
      END DO
   END DO 
   Deallocate(R_mat,R_aux,S_mat,moments, aux_r,aux_mat)
END SUBROUTINE lin_int
!============================================================== 
SUBROUTINE FINDInv(matrix, Rm, Cm, n, n_num, errorflag)
! Matrix inverse
USE o_MESH
IMPLICIT NONE

INTEGER, INTENT(IN) :: n, n_num
INTEGER, INTENT(OUT) :: errorflag  ! -1 for error, 0 otherwise
REAL(kind=WP), INTENT(IN), DIMENSION(n,n) :: matrix  
REAL(kind=WP) , DIMENSION(n,n_num) :: Rm, Cm 
	
LOGICAL :: FLAG = .TRUE.
INTEGER :: i, j, k, l
REAL(kind=WP) :: m
REAL(kind=WP), DIMENSION(n,2*n) :: augmatrix !augmented matrix
Real(kind=WP) :: inverse(n,n)
	
	
	!Augment input matrix with an identity matrix
	DO i = 1, n
           DO j = 1, 2*n
	      IF(j <= n ) THEN
		 augmatrix(i,j) = matrix(i,j)
	      ELSE IF ((i+n) == j) THEN
                 augmatrix(i,j) = 1
	      Else
		 augmatrix(i,j) = 0
	      ENDIF
	   END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k=1, n-1
	  IF(augmatrix(k,k) == 0) THEN
	     FLAG = .FALSE.
	     DO i = k+1, n
		IF (augmatrix(i,k) /= 0) THEN
		   DO j = 1,2*n
		      augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
		   END DO
		   FLAG = .TRUE.
		   EXIT
		ENDIF
		IF (FLAG .EQV. .FALSE.) THEN
		   PRINT*, "Matrix is non - invertible"
		   inverse = 0
		   errorflag = -1
		   return
		ENDIF
	      END DO
	    ENDIF
	    DO j = k+1, n			
	       m = augmatrix(j,k)/augmatrix(k,k)
	       DO i = k, 2*n
		  augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
	       END DO
	    END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
	  IF (augmatrix(i,i) == 0) THEN
	   PRINT*, "Matrix is non - invertible"
           inverse = 0
	   errorflag = -1
	   return
	  ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
	     m = augmatrix(i,i)
	   DO j = i , (2 * n)				
	     augmatrix(i,j) = (augmatrix(i,j) / m)
	   END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
	   DO i =1, k
	      m = augmatrix(i,k+1)
	      DO j = k, (2*n)
		augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
	      END DO
	  END DO
	END DO				
	
	!store answer
	DO i =1, n
	   DO j = 1, n
	      inverse(i,j) = augmatrix(i,j+n)
	  END DO
	END DO
	errorflag = 0

         DO j=1,n
	    DO i=1,n_num
	    Cm(j,i)=sum(inverse(j,:)*Rm(:,i))
	    END DO   
	 END DO
   

END SUBROUTINE FINDinv
!=========================================================================
subroutine exchange_nod_quad_int(int_array)
! A temporary exchange solution (just through cycle). 
! Larger exchange buffers will be introduced later.
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
 
 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 integer  :: nz, nh, nc, q
 real(kind=WP) :: int_array(6,nl-1,myDim_nod2D+eDim_nod2D) 
  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum
  ! Put data to be communicated into send buffer 

DO q=1,6  ! Coefficients are exchanged in serial way. To be replaced later.
  Do n=1,com_nod2D%sPEnum
   nini=com_nod2D%sptr(n)
   nend=com_nod2D%sptr(n+1)-1
   nc=0
   DO nh=nini, nend 
      DO nz=1,nl-1
      nc=nc+1
      s_buff_nod3D(n)%array(nc)=int_array(q,nz,com_nod2D%slist(nh))
      END DO
   END DO   
  end do  
 
  DO n=1, sn
     dest=com_nod2D%sPE(n)
     nini=com_nod2D%sptr(n)
     offset=(com_nod2D%sptr(n+1) - nini)*(nl-1)
     
     
     call MPI_ISEND(s_buff_nod3D(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1,rn   
     source=com_nod2D%rPE(n)
     nini=com_nod2D%rptr(n)
     offset=(com_nod2D%rptr(n+1) - nini)*(nl-1)
     
      
     call MPI_IRECV(r_buff_nod3D(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
  Do n=1,com_nod2D%rPEnum
   nini=com_nod2D%rptr(n)
   nend=com_nod2D%rptr(n+1)-1
   nc=0
   DO nh=nini, nend
      DO nz=1, nl-1
      nc=nc+1
      int_array(q,nz,com_nod2D%rlist(nh))=r_buff_nod3D(n)%array(nc)
      END DO
   END DO   
  end do 
END DO   
END SUBROUTINE exchange_nod_quad_int
!===========================================================================
SUBROUTINE adv_tracer_second_order_rec(ttf,dttf,tr_num,tt_xy,tt_xynodes)

! Direct time-space quadratic reconstruction advection scheme
! SImilar to the Miura one in technical imlementation

USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
integer      :: el(2), enodes(2), n, nz, edge, m
integer      :: nl1, nl2
real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2, flux=0.0 
real(kind=WP) :: tvert(nl), a, b, c, d, da, db, dg, wm, dbm, dgm, lx, ly
real(kind=WP) :: Tx, Ty, Tmean, rdata=0.0_WP
real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D), dttf(nl-1, myDim_nod2D+eDim_nod2D)
integer :: tr_num 
real*8      :: tt_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
real*8      :: tt_xynodes(2,nl-1,myDim_nod2D+eDim_nod2D)
ttrhs=0d0

! =================
! Horizontal advection
! =================
DO edge=1, myDim_edge2D
   enodes=edges(:,edge)   
   el=edge_tri(:,edge)
   c1=0.0_WP
   c2=0.0_WP
   nl1=nlevels(el(1))-1
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   nl2=0
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   nl2=nlevels(el(2))-1
   end if     
   ! ============
   ! the first column
   ! ============
   DO nz=1, nl1 
   ! ============
   ! Upwind implementation (quadratic reconstruction)
   ! ============
   
   
   if(UV(2,nz,el(1))*deltaX1- UV(1,nz,el(1))*deltaY1>0) then   
   
   lx=0.5_WP*(-edge_dxdy(1,edge)*r_earth*elem_cos(el(1))+ deltaX1-UV(1,nz,el(1))*dt)
   ly=0.5_WP*(-edge_dxdy(2,edge)*r_earth + deltaY1-UV(2,nz,el(1))*dt)
   
   Tmean=quad_int_coef(1,nz,enodes(2))+ &
         quad_int_coef(2,nz,enodes(2))*lx+quad_int_coef(3,nz,enodes(2))*ly + & 
         quad_int_coef(4,nz,enodes(2))* &
	 (lx*lx+(deltaX1*deltaX1+UV(1,nz,el(1))*UV(1,nz,el(1))*dt*dt)/12.0_WP)+ &
	 quad_int_coef(5,nz,enodes(2))* & 
         (lx*ly+(deltaX1*deltaY1+UV(1,nz,el(1))*UV(2,nz,el(1))*dt*dt)/12.0_WP)+ &
	 quad_int_coef(6,nz,enodes(2))* & 
         (ly*ly+(deltaY1*deltaY1+UV(2,nz,el(1))*UV(2,nz,el(1))*dt*dt)/12.0_WP)
   else
   lx=0.5_WP*(edge_dxdy(1,edge)*r_earth*elem_cos(el(1))+ deltaX1-UV(1,nz,el(1))*dt)
   ly=0.5_WP*(edge_dxdy(2,edge)*r_earth+ deltaY1-UV(2,nz,el(1))*dt)
   
   Tmean=quad_int_coef(1,nz,enodes(1))+ &
         quad_int_coef(2,nz,enodes(1))*lx+quad_int_coef(3,nz,enodes(1))*ly + & 
         quad_int_coef(4,nz,enodes(1))* &
         (lx*lx+(deltaX1*deltaX1+UV(1,nz,el(1))*UV(1,nz,el(1))*dt*dt)/12.0_WP)+ &
	 quad_int_coef(5,nz,enodes(1))* & 
         (lx*ly+(deltaX1*deltaY1+UV(1,nz,el(1))*UV(2,nz,el(1))*dt*dt)/12.0_WP)+ &
	 quad_int_coef(6,nz,enodes(1))* & 
         (ly*ly+(deltaY1*deltaY1+UV(2,nz,el(1))*UV(2,nz,el(1))*dt*dt)/12.0_WP)
   end if
   c1=UV(2,nz,el(1))*Tmean*deltaX1- UV(1,nz,el(1))*Tmean*deltaY1
   ttrhs(nz,enodes(1))=ttrhs(nz,enodes(1))+c1
   ttrhs(nz,enodes(2))=ttrhs(nz,enodes(2))-c1
   END DO
   
   ! ============
   ! the second column
   ! ============
   if(el(2)>0)then
   DO nz=1, nl2
   ! ============
   ! Upwind implementation (quadratic reconstruction)
   ! ============
   if(UV(2,nz,el(2))*deltaX2- UV(1,nz,el(2))*deltaY2<0) then   
   
   lx=0.5_WP*(-edge_dxdy(1,edge)*r_earth*elem_cos(el(2))+ deltaX2-UV(1,nz,el(2))*dt)
   ly=0.5_WP*(-edge_dxdy(2,edge)*r_earth+ deltaY2-UV(2,nz,el(2))*dt)
   
   
   Tmean=quad_int_coef(1,nz,enodes(2))+ &
         quad_int_coef(2,nz,enodes(2))*lx+quad_int_coef(3,nz,enodes(2))*ly + & 
         quad_int_coef(4,nz,enodes(2))* &
	  (lx*lx+(deltaX2*deltaX2+UV(1,nz,el(2))*UV(1,nz,el(2))*dt*dt)/12.0_WP)+ &
	 quad_int_coef(5,nz,enodes(2))* & 
          (lx*ly+(deltaX2*deltaY2+UV(1,nz,el(2))*UV(2,nz,el(2))*dt*dt)/12.0_WP)+ &
	 quad_int_coef(6,nz,enodes(2))* & 
          (ly*ly+(deltaY2*deltaY2+UV(2,nz,el(2))*UV(2,nz,el(2))*dt*dt)/12.0_WP)
   else
   lx=0.5_WP*(edge_dxdy(1,edge)*r_earth*elem_cos(el(2))+ deltaX2-UV(1,nz,el(2))*dt)
   ly=0.5_WP*(edge_dxdy(2,edge)*r_earth+ deltaY2-UV(2,nz,el(2))*dt)
   
   
   Tmean=quad_int_coef(1,nz,enodes(1))+ &
         quad_int_coef(2,nz,enodes(1))*lx+quad_int_coef(3,nz,enodes(1))*ly + & 
         quad_int_coef(4,nz,enodes(1))* &
	  (lx*lx+(deltaX2*deltaX2+UV(1,nz,el(2))*UV(1,nz,el(2))*dt*dt)/12.0_WP)+ &
	 quad_int_coef(5,nz,enodes(1))* & 
          (lx*ly+(deltaX2*deltaY2+UV(1,nz,el(2))*UV(2,nz,el(2))*dt*dt)/12.0_WP)+ &
	 quad_int_coef(6,nz,enodes(1))* & 
          (ly*ly+(deltaY2*deltaY2+UV(2,nz,el(2))*UV(2,nz,el(2))*dt*dt)/12.0_WP)
   
   end if
   c2=-UV(2,nz,el(2))*Tmean*deltaX2+ UV(1,nz,el(2))*Tmean*deltaY2
   ttrhs(nz,enodes(1))=ttrhs(nz,enodes(1))+c2
   ttrhs(nz,enodes(2))=ttrhs(nz,enodes(2))-c2
   END DO
   end if
 END DO
! ===================
! Vertical advection and diffusion
! ===================
 DO n=1, myDim_nod2D                  !! P (d) n=1, nod2D
                                      !! n=myList_nod2D(m)	    
  tvert(1) = -Wvel(1,n)*ttf(1,n)*area(1,n)		    
  ! Bottom conditions	  
  tvert(nlevels_nod2D(n))=0.	 
   DO nz=2, nlevels_nod2D(n)-1
      ! ============
      ! QUICK upwind 
      ! ============

      if(Wvel(nz,n)>0) then
        if(nz==nlevels_nod2D(n)-1) then
	  Tmean=0.5_WP*(ttf(nz-1,n)+ttf(nz,n))  ! or replace this with first
	                                       ! order upwind  ttf(nz,n)
	else
!	a=Z(nz-1)-zbar(nz) ; b=zbar(nz)-Z(nz) ; c=zbar(nz)-Z(nz+1)
	wm=Wvel(nz,n)*dt*0.5_WP
	dgm=( 2.*zbar(nz)-Z(nz-1)-Z(nz)  -wm)*wm !(b-a-wm)*wm
	dbm=(-2.*zbar(nz)+Z(nz-1)+Z(nz+1)+wm)*wm !(a-c+wm)*wm
        dg = ( (Z(nz-1)-zbar(nz))*(zbar(nz)-Z(nz  )) + dgm)  &
             / ((Z(nz-1)-Z(nz+1))*(Z(nz+1)-Z(nz))) 
        db = (-(Z(nz-1)-zbar(nz))*(zbar(nz)-Z(nz+1)) + dbm)  &
             / ((Z(nz-1)-Z(nz  ))*(Z(nz+1)-Z(nz))) 
	!dg=(a*b+dgm)/(c+a)/(b-c) ; db=(-a*c+dbm)/(b+a)/(b-c)
	da=1.0_WP-dg-db
	
	Tmean=ttf(nz-1,n)*da+ttf(nz,n)*db+ttf(nz+1,n)*dg
	end if
      end if

      if(Wvel(nz,n)<0) then
        if(nz==2) then
	  Tmean=0.5_WP*(ttf(nz-1,n)+ttf(nz,n))  
	else  
	!a=zbar(nz)-Z(nz) ; b=Z(nz-1)-zbar(nz) ;c=Z(nz-2)-zbar(nz)
	wm=-Wvel(nz,n)*dt*0.5_WP
	dgm=(-2.*zbar(nz)+Z(nz-1)+Z(nz)-wm)*wm !(b-a-wm)*wm
	dbm=( 2.*zbar(nz)-Z(nz)-Z(nz-2)+wm)*wm !(a-c+wm)*wm
	dg=( (zbar(nz)-Z(nz))*(Z(nz-1)-zbar(nz))+dgm)  &
             /((Z(nz-2)-Z(nz))*(Z(nz-1)-Z(nz-2)))
	db=(-(zbar(nz)-Z(nz))*(Z(nz-2)-zbar(nz))+dbm)  &
             /((Z(nz-1)-Z(nz))*(Z(nz-1)-Z(nz-2)))
	!dg=(a*b+dgm)/(c+a)/(b-c) ; db=(-a*c+dbm)/(b+a)/(b-c)
	da=1.0_WP-dg-db
	Tmean=ttf(nz,n)*da+ttf(nz-1,n)*db+ttf(nz-2,n)*dg
	end if
      end if
      tvert(nz)= -Tmean*Wvel(nz,n)*area(nz,n)
   END DO
 
   DO nz=1,nlevels_nod2D(n)-1
      ttrhs(nz,n)=(ttrhs(nz,n)+ &
                    (tvert(nz)-tvert(nz+1))/(zbar(nz)-zbar(nz+1)))
   END DO
END DO

!Update
DO n=1, myDim_nod2D
     DO nz=1,nlevels_nod2D(n)-1
        dttf(nz,n)=dttf(nz,n)+ttrhs(nz,n)*dt/area(nz,n)
     END DO
END DO
end subroutine adv_tracer_second_order_rec
!===========================================================================

