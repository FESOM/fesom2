! ===================================================================
! Contains routines needed for computations of dynamics.
! ===================================================================

! ===================================================================
! Stiffness matrix for the elevation
! 
! We first use local numbering and assemble the matrix
! Having completed this step we substitute global contiguous numbers.
!
! Our colind cannot be used to access local node neighbors
! This is a reminder from FESOM
!        do q=1, nghbr_nod2D(row)%nmb
!           col_pos(nghbr_nod2D(row)%addresses(q)) = q
!        enddo
!       ipos=sshstiff%rowptr(row)+col_pos(col)-sshstiff%rowptr(1)
! 
! To achive it we should use global arrays n_num and n_pos.
! Reserved for future. 


subroutine stiff_mat
  use o_PARAM
  use o_MESH
  use g_PARSUP
  use g_CONFIG
 ! use o_arrays
  !
  implicit none
  integer                             :: n, n1, n2, i, j,  row, ed
  integer                             :: enodes(2), elnodes(3), el(2)
  integer                             :: elem, npos(3), offset, nini, nend
  real(kind=WP)                       :: dmean, ff, factor, a 
  real(kind=WP)                       :: fx(3), fy(3), ax, ay
  integer, allocatable                :: n_num(:), n_pos(:,:), pnza(:), rpnza(:)
  integer, allocatable                :: mapping(:)
  logical :: flag
  character*10   :: mype_string,npes_string
  character*80   :: dist_mesh_dir,file_name

  !
  ! a)
  ssh_stiff%dim=nod2D
  allocate(ssh_stiff%rowptr(myDim_nod2D+1))
  ssh_stiff%rowptr(1)=1               ! This has to be updated as
                                      ! contiguous numbering is required
  
  allocate(n_num(myDim_nod2D+eDim_nod2D),n_pos(12,myDim_nod2D))
  n_pos=0
  
  ! b) Neighbourhood information
  DO n=1,myDim_nod2d
     n_num(n)=1
     n_pos(1,n)=n
  end do   
  Do n=1, myDim_edge2D
#ifdef ELEM_BOUNDARY
     if (boundary_nod(n)/=0) cycle
#endif
     n1=edges(1,n)
     n2=edges(2,n)

#ifdef ELEM_BOUNDARY
     if ((n1<=myDim_nod2D).and.(boundary_nod(n)==0)) then
#else
     if(n1<=myDim_nod2D) then
#endif
     n_pos(n_num(n1)+1,n1)=n2
     n_num(n1)=n_num(n1)+1
     end if
#ifdef ELEM_BOUNDARY
     if ((n2<=myDim_nod2D).and.(boundary_nod(n)==0)) then
#else
     if(n2<=myDim_nod2D) then
#endif

     n_pos(n_num(n2)+1,n2)=n1
     n_num(n2)=n_num(n2)+1
     end if
  END DO 
  
         ! n_num contains the number of neighbors
	 ! n_pos -- their indices 
  do n=1,myDim_nod2D
     ssh_stiff%rowptr(n+1)=ssh_stiff%rowptr(n)+n_num(n)
  end do
  ssh_stiff%nza=ssh_stiff%rowptr(myDim_nod2D+1)-1								 
  ! c) 
  allocate(ssh_stiff%colind(ssh_stiff%nza))
  allocate(ssh_stiff%values(ssh_stiff%nza))
  ssh_stiff%values=0.0_WP
  
  ! d) 
  do n=1,myDim_nod2D
     nini=ssh_stiff%rowptr(n)-ssh_stiff%rowptr(1)+1
     nend=ssh_stiff%rowptr(n+1)-ssh_stiff%rowptr(1)
     ssh_stiff%colind(nini:nend)= &
          n_pos(1:n_num(n),n)
  end do
  ! Thus far everything is in local numbering.
  ! We will update it later when the values are 
  ! filled in 
  !
  ! e) fill in 
  !M/dt-alpha*theta*g*dt*\nabla(H\nabla\eta))
  !  
  n_num=0
  factor=g*dt*alpha*theta
  DO ed=1,myDim_edge2D   !! Attention
   enodes=edges(:,ed)
   el=edge_tri(:,ed)
    DO j=1,2 
      row=enodes(j)
      if(row>myDim_nod2D) cycle    !! Attention
#ifdef ELEM_BOUNDARY
      if (boundary_nod(n)/=0) cycle
#endif
      ! ==========  
      offset=SSH_stiff%rowptr(row)-ssh_stiff%rowptr(1)
      DO n=1,SSH_stiff%rowptr(row+1)-SSH_stiff%rowptr(row)
	 n2=SSH_stiff%colind(offset+n)
	 n_num(n2)=offset+n
      END DO
      ! ==========
      DO i=1,2  ! Two elements related to the edge
                ! It should be just grad on elements 
      elem=el(i)
      if(elem<1) cycle
      elnodes=elem2D_nodes(:,elem)
      call elem_center(elem, ax, ay)
      ay=cos(ay)
      if(cartesian) ay=1.0_WP
      dmean=-zbar(nlevels(elem))
      dmean=0.5_WP*dmean*r_earth
      a=dmean*(coord_nod2D(2,elnodes(3))-coord_nod2D(2,elnodes(1)))
      fx(2)=a
      fx(1)=-a
      a=-dmean*(coord_nod2D(2,elnodes(2))-coord_nod2D(2,elnodes(1)))
      fx(3)=a
      fx(1)=fx(1)-a
      a=(coord_nod2D(1,elnodes(3))-coord_nod2D(1,elnodes(1)))
        if(a>cyclic_length/2.) a=a-cyclic_length
        if(a<-cyclic_length/2.) a=a+cyclic_length
      a=-ay*dmean*a
      fy(2)=a
      fy(1)=-a
      a=(coord_nod2D(1,elnodes(2))-coord_nod2D(1,elnodes(1)))
        if(a>cyclic_length/2.) a=a-cyclic_length
        if(a<-cyclic_length/2.) a=a+cyclic_length
      a=ay*dmean*a
      fy(3)=a
      fy(1)=fy(1)-a
      
      fx=fx/elem_area(elem)
      fy=fy/elem_area(elem)
      
      ! This should work too !!
      !fx=zbar(nlevels(elem))*gradient_sca(1:3,elem)
      !fy=zbar(nlevels(elem))*gradient_sca(4:6,elem)
      !
      ! fx, fy are contribution to -velocity from elem   (-h\nabla\eta)
      !
      fy=fy*(edge_cross_dxdy(2*i-1,ed))- &
         fx*(edge_cross_dxdy(2*i,ed))
   
      if(i==2) fy=-fy
      if(j==2) fy=-fy
      ! 
      ! In the computation above, I've used rules from ssh_rhs (where it is 
      ! on the rhs. So the sign is changed in the expression below.
      !
      npos=n_num(elnodes)
      SSH_stiff%values(npos)=SSH_stiff%values(npos)- fy*factor
      END DO
    END DO
  END DO 
  ! Mass matrix part
  DO row=1, myDim_nod2D
#ifdef ELEM_BOUNDARY
      if (boundary_nod(n)/=0) then
        SSH_stiff%values(row)=1d0
        cycle
      endif
#endif
     offset=ssh_stiff%rowptr(row)
     SSH_stiff%values(offset)=SSH_stiff%values(offset)+ &
                              area(1,row)/dt
  END DO
  deallocate(n_pos,n_num)

  ! =================================
  ! Global contiguous numbers:
  ! =================================
  
  ! =================================
  ! Now we need to exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! ================================= 
  allocate(pnza(npes), rpnza(npes))    
  pnza(1:npes)=0
  rpnza=0
  
  pnza(mype+1)=ssh_stiff%nza
  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
       
  if (mype==0) then
  offset=0
  else
  offset=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  ssh_stiff%rowptr=ssh_stiff%rowptr+offset   ! pointers are global
  
  ! =================================
  ! replace local nza with a global one
  ! =================================
  ssh_stiff%nza=sum(rpnza(1:npes))
  deallocate(rpnza, pnza)

  ! ==================================
  ! colindices are now local to PE. We need to make them local
  ! contiguous
  ! ==================================
  allocate(mapping(nod2d))
  write(npes_string,"(I10)") npes
  dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'
  file_name=trim(dist_mesh_dir)//'rpart.out'
  open(10,file=trim(dist_mesh_dir)//'rpart.out')
  read(10,*) n
  read(10,*) mapping(1:npes)
  read(10,*) mapping
  close(10) 
      ! (i) global natural: 
    do n=1,ssh_stiff%rowptr(myDim_nod2D+1)-ssh_stiff%rowptr(1)  
       ssh_stiff%colind(n)=myList_nod2D(ssh_stiff%colind(n))    
    end do
    ! (ii) global PE contiguous: 
    do n=1,ssh_stiff%rowptr(myDim_nod2D+1)-ssh_stiff%rowptr(1)  
       ssh_stiff%colind(n)=mapping(ssh_stiff%colind(n))    
    end do
  ! ==================================
  deallocate(mapping)

end subroutine stiff_mat
! =================================================================
! =================================================================
subroutine solve_ssh
  use o_PARAM
  use o_MESH
  use o_ARRAYS
  use g_PARSUP
  use g_comm_auto
  !
#ifdef PETSC
  implicit none
#include "petscf.h"
  
  
  integer                              :: myrows
  integer                              :: Pmode
  real(kind=8)                         :: rinfo(20,20)
  integer                              :: maxiter=2000
  integer                              :: restarts=15
  integer                              :: fillin=3
  integer                              :: lutype=2
  integer                              :: nrhs=1
  real(kind=8)                         :: droptol=1.e-7
  real(kind=8)                         :: soltol =1e-10  !1.e-10
  !integer, allocatable                 :: part(:)    
  logical, save                        :: lfirst=.true.
  real(kind=WP), allocatable           :: arr_nod2D(:),arr_nod2D2(:,:),arr_nod2D3(:)
  real*8 :: cssh1,cssh2,crhs
  integer i
     Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB +PET_REPORT + PET_QUIET+ PET_RCM+PET_PCBJ
     !Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB +PET_REPORT + PET_QUIET
     if (lfirst) then   
        Pmode = Pmode+PET_STRUCT+PET_PMVALS + PET_PCASM+PET_OVL_2 !+PET_PCBJ+PET_ILU
        lfirst=.false.
     end if
     call PETSC_S(Pmode, 1, ssh_stiff%dim, ssh_stiff%nza, myrows, &
          maxiter, & 
          restarts, &
          fillin,  &
          droptol, &  
          soltol,  &
          part, ssh_stiff%rowptr, ssh_stiff%colind, ssh_stiff%values, &
          ssh_rhs, d_eta, &
          rinfo, MPI_COMM_WORLD)
#elif defined(PARMS)
  !use o_solver

  implicit none

#include "fparms.h"

  logical, save                        :: lfirst=.true.
  integer    :: ident
  integer               :: n3, comm, reuse
  integer               :: maxiter, restart, lutype, fillin
  real(kind=8)          :: droptol, soltol
  COMM=MPI_COMM_WORLD
  ident=1
  maxiter=2000
  restart=15
  fillin=3
  lutype=2
  droptol=1.e-8
  soltol=1.e-10
  reuse=0
  call MPI_Barrier(MPI_COMM_WORLD, MPIERR)

    if(lfirst) then
!      call solver_init(ident, ssh_stiff, SOLBICGS, PCBJ, PCILUK, &
!                       lutype, fillin, droptol, maxiter, restart, soltol, 0)
!subroutine solver_init(ident, stiff, stype, pctype, pcilutype, &
!            lutype, fillin, droptol, maxiter, restart, soltol, reuse)

       call psolver_init(ident, SOLBICGS, PCBJ, PCILUK, lutype, &
              fillin, droptol, maxiter, restart, soltol, &
              part-1, ssh_stiff%rowptr(:)-ssh_stiff%rowptr(1), &
              ssh_stiff%colind-1, ssh_stiff%values, reuse, MPI_COMM_WORLD)
	lfirst=.false.
    end if
    call psolve(ident, ssh_rhs, 0., d_eta, 0, MPI_COMM_WORLD)

#endif
 call exchange_nod(d_eta)
end subroutine solve_ssh
!===========================================================================
SUBROUTINE compute_ssh_rhs
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
  use g_comm_auto
IMPLICIT NONE

! In the semiimplicit method: 
! ssh_rhs=-\nabla\int(U_n+alpha*U_rhs) dz

integer      :: ed, el(2), enodes(2),  nz
real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2 
ssh_rhs=0.0_WP
DO ed=1, myDim_edge2D                     !!  m=1, myDim_edge2D
                                          !!ed=myList_edge2D(m)
   enodes=edges(:,ed)
   el=edge_tri(:,ed)
   
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,ed)
   deltaY1=edge_cross_dxdy(2,ed)
   if(el(2)>0)  then
    deltaX2=edge_cross_dxdy(3,ed)
    deltaY2=edge_cross_dxdy(4,ed)
   end if
   DO nz=1, nlevels(el(1))-1
    c1=c1+((UV(2,nz,el(1))+alpha*UV_rhs(2,nz,el(1)))*deltaX1- &
          (UV(1,nz,el(1))+alpha*UV_rhs(1,nz,el(1)))*deltaY1)*(zbar(nz)-zbar(nz+1))
   END DO
 
   if(el(2)>0) then
    DO nz=1, nlevels(el(2))-1
     c2=c2-((UV(2,nz,el(2))+alpha*UV_rhs(2,nz,el(2)))*deltaX2- &
         (UV(1,nz,el(2))+alpha*UV_rhs(1,nz,el(2)))*deltaY2)*(zbar(nz)-zbar(nz+1))
    END DO
   end if
   ssh_rhs(enodes(1))=ssh_rhs(enodes(1))+(c1+c2)
   ssh_rhs(enodes(2))=ssh_rhs(enodes(2))-(c1+c2)
   
END DO

! The rhs is the total transport into the column
! associated with the scalar cell
  
  !ssh_rhs(1:myDim_nod2D)=eta_n(1:myDim_nod2D)*area(1,:)/dt+ssh_rhs(1:myDim_nod2D)
!! P-E should be just added as it is to the rhs.  
!call exchange_nod(ssh_rhs)
!call test_ssh_rhs ! delete me

END subroutine compute_ssh_rhs

! ===================================================================================

subroutine test_ssh_rhs
USE g_PARSUP
use o_MESH
use o_ARRAYS
real*8 c1,c2
  c1=sum(ssh_rhs(1:mydim_nod2d))
  c2=0d0
  call MPI_Allreduce (c1, c2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
  if (mype==0) write(*,*) 'DBG:sum(ssh_rhs)=',c2
end subroutine test_ssh_rhs
! ===================================================================================
subroutine compute_vel_nodes
USE o_MESH
USE o_PARAM
USE o_ARRAYS
USE g_PARSUP
use g_comm_auto
IMPLICIT NONE
integer            :: n, nz, k, elem
real(kind=WP)      :: tx, ty, tvol
 DO n=1, myDim_nod2D 
    DO nz=1, nlevels_nod2D(n)-1
       tvol=0.0_WP
       tx=0.0_WP
       ty=0.0_WP
       DO k=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(k,n)
           if(nlevels(elem)-1<nz) cycle
           tvol=tvol+elem_area(elem)
           tx=tx+UV(1,nz,elem)*elem_area(elem)
           ty=ty+UV(2,nz,elem)*elem_area(elem)
        END DO
        Unode(1,nz,n)=tx/tvol
        Unode(2,nz,n)=ty/tvol
    END DO
 END DO
   call exchange_nod(Unode(1,:,:))
   call exchange_nod(Unode(2,:,:))
end subroutine compute_vel_nodes
!===========================================================================
subroutine impl_vert_visc
USE o_MESH
USE o_PARAM
USE o_ARRAYS
USE g_PARSUP
USE g_CONFIG,only: dt
IMPLICIT NONE

real(kind=WP)              ::  a(nl-1), b(nl-1), c(nl-1), ur(nl-1), vr(nl-1)
real(kind=WP)              ::  cp(nl-1), up(nl-1), vp(nl-1)
integer                    ::  nz, elem, nzmax, elnodes(3)
real(kind=WP)              ::  zinv, m, friction, wu, wd
 DO elem=1,myDim_elem2D
    elnodes=elem2D_nodes(:,elem)
    nzmax=nlevels(elem)
    ! ==========================
    ! Operator
    ! ==========================
    ! Regular part of coefficients:
    DO nz=2, nzmax-2
       zinv=1.0_WP*dt/(zbar(nz)-zbar(nz+1))
       a(nz)=-Av(nz,elem)/(Z(nz-1)-Z(nz))*zinv
       c(nz)=-Av(nz+1,elem)/(Z(nz)-Z(nz+1))*zinv
       b(nz)=-a(nz)-c(nz)+1.0_WP
    !Update from the vertical advection
       wu=sum(Wvel_i(nz,   elnodes))/3.
       wd=sum(Wvel_i(nz+1, elnodes))/3.
       a(nz)=a(nz)+min(0._WP, wu)*zinv
       b(nz)=b(nz)+max(0._WP, wu)*zinv

       b(nz)=b(nz)-min(0._WP, wd)*zinv
       c(nz)=c(nz)-max(0._WP, wd)*zinv
    END DO
   ! The last row
       zinv=1.0_WP*dt/(zbar(nzmax-1)-zbar(nzmax))
       a(nzmax-1)=-Av(nzmax-1,elem)/(Z(nzmax-2)-Z(nzmax-1))*zinv
       b(nzmax-1)=-a(nzmax-1)+1.0_WP
       c(nzmax-1)=0.0_WP

    !Update from the vertical advection
       wu=sum(Wvel_i(nzmax-1, elnodes))/3.
       a(nzmax-1)=a(nzmax-1)+min(0._WP, wu)*zinv
       b(nzmax-1)=b(nzmax-1)+max(0._WP, wu)*zinv

   ! The first row
       zinv=1.0_WP*dt/(zbar(1)-zbar(2))
       c(1)=-Av(2,elem)/(Z(1)-Z(2))*zinv
       a(1)=0.0_WP
       b(1)=-c(1)+1.0_WP
    !Update from the vertical advection
       wu=sum(Wvel_i(1, elnodes))/3.
       wd=sum(Wvel_i(2, elnodes))/3.

       b(1)=b(1)+wu*zinv
       b(1)=b(1)-min(0._WP, wd)*zinv
       c(1)=c(1)-max(0._WP, wd)*zinv
   ! ===========================
   ! The rhs:
   ! ===========================
    ur(1:nzmax-1)=UV_rhs(1,1:nzmax-1,elem)
    vr(1:nzmax-1)=UV_rhs(2,1:nzmax-1,elem)
   ! The first row contains surface forcing
    ur(1)= ur(1)+zinv*stress_surf(1,elem)/density_0
    vr(1)= vr(1)+zinv*stress_surf(2,elem)/density_0
   ! The last row contains bottom friction
    zinv=1.0_WP*dt/(zbar(nzmax-1)-zbar(nzmax))
    friction=-C_d*sqrt(UV(1,nlevels(elem)-1,elem)**2+ &
                      UV(2,nlevels(elem)-1,elem)**2)
    ur(nzmax-1)=ur(nzmax-1)+zinv*friction*UV(1,nzmax-1,elem)
    vr(nzmax-1)=vr(nzmax-1)+zinv*friction*UV(2,nzmax-1,elem)
    ! Model solves for the difference to the timestep N and therefore we need to 
    ! update the RHS for advective and diffusive contributions at the previous time step
    DO nz=2, nzmax-2
	ur(nz)=ur(nz)-a(nz)*UV(1,nz-1,elem)-(b(nz)-1.0_WP)*UV(1,nz,elem)-c(nz)*UV(1,nz+1,elem)
	vr(nz)=vr(nz)-a(nz)*UV(2,nz-1,elem)-(b(nz)-1.0_WP)*UV(2,nz,elem)-c(nz)*UV(2,nz+1,elem)
    END DO
    ur(1)=ur(1)-(b(1)-1.0_WP)*UV(1,1,elem)-c(1)*UV(1,2,elem)
    vr(1)=vr(1)-(b(1)-1.0_WP)*UV(2,1,elem)-c(1)*UV(2,2,elem)

    ur(nzmax-1)=ur(nzmax-1)-a(nzmax-1)*UV(1,nzmax-2,elem)-(b(nzmax-1)-1.0_WP)*UV(1,nzmax-1,elem)
    vr(nzmax-1)=vr(nzmax-1)-a(nzmax-1)*UV(2,nzmax-2,elem)-(b(nzmax-1)-1.0_WP)*UV(2,nzmax-1,elem)
   ! ===========================
   ! The sweep algorithm
   ! ===========================
   ! initialize c-prime and s,t-prime
    cp(1) = c(1)/b(1)
    up(1) = ur(1)/b(1)
    vp(1) = vr(1)/b(1)
   ! solve for vectors c-prime and t, s-prime
    do nz = 2,nzmax-1
       m = b(nz)-cp(nz-1)*a(nz)
       cp(nz) = c(nz)/m
       up(nz) = (ur(nz)-up(nz-1)*a(nz))/m
       vp(nz) = (vr(nz)-vp(nz-1)*a(nz))/m
    enddo
   ! initialize x
    ur(nzmax-1) = up(nzmax-1)
    vr(nzmax-1) = vp(nzmax-1)
   ! solve for x from the vectors c-prime and d-prime
    do nz = nzmax-2, 1, -1
       ur(nz) = up(nz)-cp(nz)*ur(nz+1)
       vr(nz) = vp(nz)-cp(nz)*vr(nz+1)
    end do
   ! ===========================
   ! RHS update
   ! ===========================
    DO nz=1,nzmax-1
       UV_rhs(1,nz,elem)=ur(nz)
       UV_rhs(2,nz,elem)=vr(nz)
    END DO
 END DO   !!! cycle over elements
END subroutine impl_vert_visc        
!===========================================================================
SUBROUTINE vel_gradients
! Compute derivatives of velocity by least square interpolation.
! The interpolation coefficients are already saved
! For the no-slip case, it is assumed that velocity at 
! the boundary edge == 0. For the free-slip case, there are only 2
! neighbours
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
  use g_comm_auto
IMPLICIT NONE
real(kind=WP)    :: x(3), y(3), u, v, r1, r2
real(kind=WP)    :: zc1, zc2, zc_inv, un
integer         :: elem, el, j, nz, m
DO elem=1, myDim_elem2D     !! m=1, myDim_elem2D
   !! elem=myList_elem2D(m)
   vel_grad(1:nlevels(elem)-1,:,elem) = 0.0_WP

   DO j=1,3
      el=elem_neighbors(j,elem)

      if (el>0) then
         ! ======================
         ! fill in virtual values in the velocity array
         ! ======================
         ! Filling in velocity gradients:
         ! ======================
         if(nlevels(el)<nlevels(elem)) then
            DO nz=1, nlevels(el)-1
               u=UV(1,nz,el)-UV(1,nz,elem)
               v=UV(2,nz,el)-UV(2,nz,elem)
               vel_grad(nz,1,elem)=vel_grad(nz,1,elem)+gradient_vec(j,elem)*u
               vel_grad(nz,2,elem)=vel_grad(nz,2,elem)+gradient_vec(j+3,elem)*u
               vel_grad(nz,3,elem)=vel_grad(nz,3,elem)+gradient_vec(j,elem)*v
               vel_grad(nz,4,elem)=vel_grad(nz,4,elem)+gradient_vec(j+3,elem)*v
            END DO
            if(free_slip) then
               zc1=edge_dxdy(1,elem_edges(j,elem))*elem_cos(elem)
               zc2=edge_dxdy(2,elem_edges(j,elem))
               zc_inv = 1./(zc1*zc1+zc2*zc2)
               DO nz=nlevels(el),nlevels(elem)-1
                  un = -2.*(UV(1,nz,elem)*zc2 - UV(2,nz,elem)*zc1)*zc_inv
                  u = un*zc2
                  v =-un*zc1
                  vel_grad(nz,1,elem)=vel_grad(nz,1,elem)+gradient_vec(j,elem)*u
                  vel_grad(nz,2,elem)=vel_grad(nz,2,elem)+gradient_vec(j+3,elem)*u
                  vel_grad(nz,3,elem)=vel_grad(nz,3,elem)+gradient_vec(j,elem)*v
                  vel_grad(nz,4,elem)=vel_grad(nz,4,elem)+gradient_vec(j+3,elem)*v
               END DO
            else     ! noslip
               DO nz=nlevels(el),nlevels(elem)-1
                  u= -2.*UV(1,nz,elem)
                  v= -2.*UV(2,nz,elem)
                  vel_grad(nz,1,elem)=vel_grad(nz,1,elem)+gradient_vec(j,elem)*u
                  vel_grad(nz,2,elem)=vel_grad(nz,2,elem)+gradient_vec(j+3,elem)*u
                  vel_grad(nz,3,elem)=vel_grad(nz,3,elem)+gradient_vec(j,elem)*v
                  vel_grad(nz,4,elem)=vel_grad(nz,4,elem)+gradient_vec(j+3,elem)*v    
               END DO
            end if
         
         else
            DO nz=1, nlevels(elem)-1
               u=UV(1,nz,el)-UV(1,nz,elem)
               v=UV(2,nz,el)-UV(2,nz,elem)
               vel_grad(nz,1,elem)=vel_grad(nz,1,elem)+gradient_vec(j,elem)*u
               vel_grad(nz,2,elem)=vel_grad(nz,2,elem)+gradient_vec(j+3,elem)*u
               vel_grad(nz,3,elem)=vel_grad(nz,3,elem)+gradient_vec(j,elem)*v
               vel_grad(nz,4,elem)=vel_grad(nz,4,elem)+gradient_vec(j+3,elem)*v

            END DO
         end if
      else

         ! ===============
         ! Boundary element
         ! ===============
         !    (Here we do not have place for virtual velocities
         !     in the velocity array so we use auxiliary array)
         ! ======================
         if(free_slip) then
            zc1=edge_dxdy(1,elem_edges(j,elem))*elem_cos(elem)
            zc2=edge_dxdy(2,elem_edges(j,elem))
            zc_inv = 1./(zc1*zc1+zc2*zc2)
            DO nz=1,nlevels(elem)-1
               un=-2.*(UV(1,nz,elem)*zc2-UV(2,nz,elem)*zc1)*zc_inv
               u = un*zc2
               v = -un*zc1
               vel_grad(nz,1,elem)=vel_grad(nz,1,elem)+gradient_vec(j,elem)*u
               vel_grad(nz,2,elem)=vel_grad(nz,2,elem)+gradient_vec(j+3,elem)*u
               vel_grad(nz,3,elem)=vel_grad(nz,3,elem)+gradient_vec(j,elem)*v
               vel_grad(nz,4,elem)=vel_grad(nz,4,elem)+gradient_vec(j+3,elem)*v
            END DO
         else     ! noslip
            DO nz=1,nlevels(elem)-1
               u=-2.*UV(1,nz,elem)
               v=-2.*UV(2,nz,elem)
               vel_grad(nz,1,elem)=vel_grad(nz,1,elem)+gradient_vec(j,elem)*u
               vel_grad(nz,2,elem)=vel_grad(nz,2,elem)+gradient_vec(j+3,elem)*u
               vel_grad(nz,3,elem)=vel_grad(nz,3,elem)+gradient_vec(j,elem)*v
               vel_grad(nz,4,elem)=vel_grad(nz,4,elem)+gradient_vec(j+3,elem)*v    
            END DO
         end if
      end if
   end do   ! cycle over neighbor elements
END DO

do j=1,4
 call exchange_elem(vel_grad(:,j,:))
enddo
END SUBROUTINE vel_gradients
!===========================================================================
SUBROUTINE compute_vel_rhs
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
  use g_comm_auto
IMPLICIT NONE
integer          :: elem, elnodes(3), nz 
real(kind=WP)    :: eta(3), ff, gg, mm 
real(kind=WP)    :: Fx, Fy, pre(3)
logical, save    :: lfirst=.true.
real(kind=WP)    :: t1, t2, t3, t4, t5, t6

t1=MPI_Wtime()
  
 ! =================
 ! Take care of the AB part
 ! =================
 Do elem=1, myDim_elem2D    !!m=1, myDim_elem2D           (a)
                            !!elem=myList_elem2D(m)
    DO nz=1,nl-1 
       UV_rhs(1,nz,elem)=-(0.5_WP+epsilon)*UV_rhsAB(1,nz,elem)   
       UV_rhs(2,nz,elem)=-(0.5_WP+epsilon)*UV_rhsAB(2,nz,elem)
    END DO
 END DO 
 ! ====================
 ! Sea level and pressure contribution   -\nabla(\eta +hpressure/rho_0)
 ! and the Coriolis force + metric terms
 ! ====================
 DO elem=1,  myDim_elem2D           !! m=1,  myDim_elem2D          
                                    !! elem=myList_elem2D(m)
    elnodes=elem2D_nodes(:,elem)
    !eta=g*eta_n(elnodes)*(1-theta)        !! this place needs update (1-theta)!!!
    eta=g*eta_n(elnodes)
    gg=elem_area(elem)
    ff=coriolis(elem)*gg
    !mm=metric_factor(elem)*gg
    if(mom_adv==4) mm=0.0_WP
    DO nz=1,nlevels(elem)-1
      pre=-(eta+hpressure(nz,elnodes)/density_0)
      Fx=sum(gradient_sca(1:3,elem)*pre)
      Fy=sum(gradient_sca(4:6,elem)*pre)
      UV_rhs(1,nz,elem)=UV_rhs(1,nz,elem)+Fx*gg 
      UV_rhs(2,nz,elem)=UV_rhs(2,nz,elem)+Fy*gg 
      UV_rhsAB(1,nz,elem)= UV(2,nz,elem)*ff! + mm*UV(1,nz,elem)*UV(2,nz,elem)
      UV_rhsAB(2,nz,elem)=-UV(1,nz,elem)*ff! - mm*UV(1,nz,elem)*UV(2,nz,elem)
    END DO
 END DO
t2=MPI_Wtime() 
 ! ====================
 ! Compute velocity gradients
 ! (to be used in viscosity operator
 ! and in flux estimates)
 ! ====================
 call vel_gradients
t3=MPI_Wtime() 
 ! ====================
 ! Horizontal advection
 ! ====================
 if(mom_adv==1) call momentum_adv_upwind   
 If(mom_adv==2) call momentum_adv_p1
 if(mom_adv==3) call momentum_adv_scalar

t4=MPI_Wtime() 
 !======
 ! Horizontal viscosity part
 !======
  if(laplacian.or.biharmonic) call h_viscosity2

  if(laplacian.and.biharmonic) then    ! Both harmonic and biharmonic 
  call biPlusharmonic_viscosity
  elseif(laplacian) then               ! Only harmonic
  call laplacian_viscosity 
  elseif (biharmonic) then
  call biharmonic_viscosity            ! Only biharmonic
  end if
 ! ======================
 ! subroutine biharmonic redefines gradients!
 ! Do not move it from here
 ! ======================
 
t5=MPI_Wtime() 
 !======
 ! Vertical advection and viscosity
 !======
 if(mom_adv==3) then
    !call momentum_vert_visc
 else
    call momentum_vert_adv_visc
 end if 

 ! =======================
 ! Update the rhs   
 ! =======================
 ff=(1.5_WP+epsilon)
  if(lfirst.and.(.not.r_restart)) then
    ff=1.0_WP
    lfirst=.false.
  end if
  
 DO elem=1, myDim_elem2D       !! m=1, myDim_elem2D 
                               !! elem=myList_elem2D(m)
    DO nz=1,nlevels(elem)-1   
       UV_rhs(1,nz,elem)=dt*(UV_rhs(1,nz,elem)+UV_rhsAB(1,nz,elem)*ff)/elem_area(elem)
       UV_rhs(2,nz,elem)=dt*(UV_rhs(2,nz,elem)+UV_rhsAB(2,nz,elem)*ff)/elem_area(elem)
    END DO 

 END DO

 ! =======================  
 ! U_rhs contains all contributions to velocity from old time steps   
 ! =======================
 t6=MPI_Wtime() 
call MPI_Barrier(MPI_COMM_WORLD, MPIERR)

 if(mype==0) then 
   write(*,*) 'Momentum:   ', t6-t1
   write(*,*) 'pres., Cor: ', t2-t1
   write(*,*) 'grads:      ', t3-t2
   write(*,*) 'h adv       ', t4-t3
   write(*,*) 'h visc      ', t5-t4
   write(*,*) 'vert. part  ', t6-t5
 end if     
END SUBROUTINE compute_vel_rhs
!==========================================================================
SUBROUTINE update_vel
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
  use g_comm_auto
#ifdef BTR_SPLIT
use o_split
#endif
IMPLICIT NONE
integer   :: elem, elnodes(3), nz, m
real(kind=WP) :: eta(3) 
real(kind=WP) :: Fx, Fy

 DO elem=1, myDim_elem2D                 !! m=1, myDim_elem2D  
                                         !! elem=myList_elem2D(m)
    elnodes=elem2D_nodes(:,elem)
    !eta=-g*theta*dt*eta_n(elnodes)
    eta=-g*theta*dt*d_eta(elnodes)
    Fx=sum(gradient_sca(1:3,elem)*eta)
    Fy=sum(gradient_sca(4:6,elem)*eta)
#ifdef BTR_SPLIT
    DO nz=1, nlevels(elem)-1
    UV(1,nz,elem)= UV(1,nz,elem) + UV_rhs(1,nz,elem) + Fx
    UV(2,nz,elem)= UV(2,nz,elem) + UV_rhs(2,nz,elem) + Fy
    END DO
    rhs_btr=0d0
!    depth0(elem)=depth0(elem)+sum(elev_mean(elem2D_nodes(:,elem))-elev(elem2D_nodes(:,elem)))/3d0!TODO Should we do it?
    DO nz=1, nlevels(elem)-1
      rhs_btr(elem,1)=rhs_btr(elem,1)+UV(1,nz,elem)*(zbar(nz)-zbar(nz+1))
      rhs_btr(elem,2)=rhs_btr(elem,2)+UV(2,nz,elem)*(zbar(nz)-zbar(nz+1)) 
    ENDDO
    rhs_btr(elem,:)=rhs_btr(elem,:)/depth0(elem)
    DO nz=1, nlevels(elem)-1
    UV(1,nz,elem)=  UV(1,nz,elem) - ( rhs_btr(elem,1) - vel_mean(elem,1))
    UV(2,nz,elem)=  UV(2,nz,elem) - ( rhs_btr(elem,2) - vel_mean(elem,2))
    ENDDO
#else
    DO nz=1, nlevels(elem)-1
     UV(1,nz,elem)= UV(1,nz,elem) + UV_rhs(1,nz,elem) + Fx
     UV(2,nz,elem)= UV(2,nz,elem) + UV_rhs(2,nz,elem) + Fy
    END DO
#endif
 END DO
 eta_n=eta_n+d_eta
 call exchange_elem(UV(1,:,:))
 call exchange_elem(UV(2,:,:))

end subroutine update_vel
!==========================================================================
SUBROUTINE vert_vel
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
  use g_comm_auto
IMPLICIT NONE
integer       :: el(2), enodes(2), n, nz, edge
real(kind=WP) :: ff, c1, c2, deltaX1, deltaY1, deltaX2, deltaY2 

! =================
! Set w to zero    ???? This is seemingly not needed
! =================
   Do n=1, myDim_nod2D        !! m=1, myDim_nod2D       
                              !! n=myList_nod2D(m)
      DO nz=1, nl
         wvel(nz,n)=0.0_8
         if (Fer_GM) then
            fer_wvel(nz,n)=0.0_8
         end if
      END DO
   END DO
! =================
! Contributions from levels in divergence
! =================
 DO edge=1, myDim_edge2D          !! m=1, myDim_edge2D 
                                  !! edge=myList_edge2D(m)
   enodes=edges(:,edge)   
   el=edge_tri(:,edge)
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   end if     
   DO nz=nlevels(el(1))-1,1,-1
   c1=(UV(2,nz,el(1))*deltaX1- &
      UV(1,nz,el(1))*deltaY1)*(zbar(nz)-zbar(nz+1))
      
   Wvel(nz,enodes(1))=Wvel(nz,enodes(1))+c1
   Wvel(nz,enodes(2))=Wvel(nz,enodes(2))-c1
   if (Fer_GM) then
      c1=(fer_UV(2,nz,el(1))*deltaX1- &
         fer_UV(1,nz,el(1))*deltaY1)*(zbar(nz)-zbar(nz+1))
      
      fer_Wvel(nz,enodes(1))=fer_Wvel(nz,enodes(1))+c1
      fer_Wvel(nz,enodes(2))=fer_Wvel(nz,enodes(2))-c1
   end if   
   END DO

   if(el(2)>0)then
   DO nz=nlevels(el(2))-1,1,-1
   c2=-(UV(2,nz,el(2))*deltaX2- &
      UV(1,nz,el(2))*deltaY2)*(zbar(nz)-zbar(nz+1))
   Wvel(nz,enodes(1))=Wvel(nz,enodes(1))+c2
   Wvel(nz,enodes(2))=Wvel(nz,enodes(2))-c2

   if (Fer_GM) then
      c2=-(fer_UV(2,nz,el(2))*deltaX2- &
         fer_UV(1,nz,el(2))*deltaY2)*(zbar(nz)-zbar(nz+1))
      fer_Wvel(nz,enodes(1))=fer_Wvel(nz,enodes(1))+c2
      fer_Wvel(nz,enodes(2))=fer_Wvel(nz,enodes(2))-c2
   end if
   
   END DO
   end if

END DO

! ===================
! Sum up to get W
! ===================
 Do n=1, myDim_nod2D                   !! m=1, myDim_nod2D 
                                       !! n=myList_nod2D(m)
   DO nz=nl-1,1,-1
   Wvel(nz,n)=Wvel(nz,n)+Wvel(nz+1,n)
   if (Fer_GM) then 
      fer_Wvel(nz,n)=fer_Wvel(nz,n)+fer_Wvel(nz+1,n)
   end if
   END DO
END DO

 Do n=1, myDim_nod2D                   !! m=1, myDim_nod2D
                                       !! n=myList_nod2D(m)
   DO nz=1,nlevels_nod2D(n)-1
      Wvel(nz,n)=Wvel(nz,n)/area(nz,n)
      if (Fer_GM) then 
         fer_Wvel(nz,n)=fer_Wvel(nz,n)/area(nz,n)          
      end if
   END DO
END DO

call exchange_nod(Wvel)
call exchange_nod(fer_Wvel)

! Split implicit vertical velocity onto implicit and explicit components
if (w_split) then
      Do n=1, myDim_nod2D+eDim_nod2D
         DO nz=1,nlevels_nod2D(n)-1
            Wvel_e(nz,n)=min(max(Wvel(nz,n), -w_exp_max), w_exp_max)
         END DO
      END DO
else
      Wvel_e=Wvel
end if
Wvel_i=Wvel-Wvel_e
end subroutine vert_vel
! =================================================================
SUBROUTINE momentum_adv_upwind
! 
! Linear reconstruction upwind horizontal momentum advection.
! It is typically too damping or too noisy, so other options are
! recommended instead.
! 
USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP

IMPLICIT NONE
integer          :: ed, nodes(2), el(2), nz, nl1, nl2
real(kind=WP)    :: x1, y1, x2, y2, u1, u2, v1, v2
real(kind=WP)    :: un, xe, ye

! ======================
! Horizontal momentum advection   -\int div(uu)dS=-\sum u(un)l
! (assembly over edges)        
! ======================


 DO ed=1, myDim_edge2D+eDim_edge2D  !! m=1, myDim_edge2D+eDim_edge2D   
                                     !! ed=myList_edge2D(m)
   if(myList_edge2D(ed)>edge2D_in) cycle
   
   nodes=edges(:,ed)   
   el=edge_tri(:,ed)
   nl1=nlevels(el(1))-1
   nl2=nlevels(el(2))-1
   
   x1=-edge_cross_dxdy(1,ed)
   y1=-edge_cross_dxdy(2,ed)
   x2=-edge_cross_dxdy(3,ed)
   y2=-edge_cross_dxdy(4,ed)
   xe=edge_dxdy(1,ed)
   ye=edge_dxdy(2,ed)
   
   DO nz=1, min(nl1,nl2)   ! Only faces that do not belong to
                           ! vertical walls can contribute to
			   ! the momentum advection 
   !====== 
   ! The piece below gives second order spatial accuracy for
   ! the momentum fluxes. 
   !======
   
   u1=UV(1,nz,el(1))+ vel_grad(nz,1,el(1))*x1+vel_grad(nz,2,el(1))*y1
   v1=UV(2,nz,el(1))+ vel_grad(nz,3,el(1))*x1+vel_grad(nz,4,el(1))*y1
   u2=UV(1,nz,el(2))+ vel_grad(nz,1,el(2))*x2+vel_grad(nz,2,el(2))*y2
   v2=UV(2,nz,el(2))+ vel_grad(nz,3,el(2))*x2+vel_grad(nz,4,el(2))*y2
   
   !======
   ! Normal velocity at edge ed directed to el(2)
   ! (outer to el(1)) multiplied with the length of the edge
   ! and mean depth dmean
   !======
   un=0.5_WP*r_earth*((u1+u2)*ye-(v1*elem_cos(el(1))+v2*elem_cos(el(2)))*xe)
   !======
   ! If it is positive, take velocity in the left element (el(1)),
   ! and use the velocity at el(2) otherwise.
   !======  
   
  if(un>=0) then    
  u2=u1
  v2=v1
  end if   
   
   UV_rhsAB(1,nz,el(1))=UV_rhsAB(1,nz,el(1))-u2*un
   UV_rhsAB(2,nz,el(1))=UV_rhsAB(2,nz,el(1))-v2*un
   UV_rhsAB(1,nz,el(2))=UV_rhsAB(1,nz,el(2))+u2*un
   UV_rhsAB(2,nz,el(2))=UV_rhsAB(2,nz,el(2))+v2*un
 
   END DO
 END DO  

END SUBROUTINE momentum_adv_upwind
! ============================================================================
SUBROUTINE momentum_vert_adv_visc
! 
! Vertical momentum advection and viscosity
! For advection, quadratic upwind reconstruction is used.
! 
 
USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE

integer          :: elem, elnodes(3), nz 
real(kind=WP)    :: friction,temp
real(kind=WP)    :: w, uvertAB(2, nl)
real(kind=WP)    :: uvert(2,nl), umean, vmean, a, b, c, d, dg, da, db

! =======================
! Vertical momentum advection 
! and vertical viscosity
! =======================
 uvert=0.0_WP
 uvertAB=0.0_WP

 DO elem=1, myDim_elem2D         !! m=1, myDim_elem2D        !! P (d) elem=1,elem2D
                                 !! elem=myList_elem2D(m)
   elnodes=elem2D_nodes(:,elem)
   DO nz=2, nlevels(elem)-1
          w=sum(Wvel(nz,elnodes))*elem_area(elem)/3.0_WP
	  umean=0.5_WP*(UV(1,nz-1,elem)+UV(1,nz,elem))
  	  vmean=0.5_WP*(UV(2,nz-1,elem)+UV(2,nz,elem))
      if(w>0) then
        if(nz==nlevels(elem)-1) then
	    umean=0.5_WP*(UV(1,nz-1,elem)+UV(1,nz,elem))
  	    vmean=0.5_WP*(UV(2,nz-1,elem)+UV(2,nz,elem))
	                                       ! or replace this with first
	                                       ! order upwind  
	else
	a=Z(nz-1)-zbar(nz)
	b=zbar(nz)-Z(nz)
	c=zbar(nz)-Z(nz+1)
	d=(c+a)*(b**2-a**2)-(c**2-a**2)*(b+a)
        dg=a*b*(a+b)/d
	db=-a*c*(a+c)/d
	da=1.0_WP-dg-db
	umean=UV(1,nz-1,elem)*da+UV(1,nz,elem)*db+UV(1,nz+1,elem)*dg
	vmean=UV(2,nz-1,elem)*da+UV(2,nz,elem)*db+UV(2,nz+1,elem)*dg
	end if
      end if

      if(w<0) then
        if(nz==2) then
	    umean=0.5_WP*(UV(1,nz-1,elem)+UV(1,nz,elem))
  	    vmean=0.5_WP*(UV(2,nz-1,elem)+UV(2,nz,elem))
	else  
	a=zbar(nz)-Z(nz)
	b=Z(nz-1)-zbar(nz)
	c=Z(nz-2)-zbar(nz)
	d=(c+a)*(b**2-a**2)-(c**2-a**2)*(b+a)
        dg=a*b*(a+b)/d
	db=-a*c*(a+c)/d
	da=1.0_WP-dg-db
	umean=UV(1,nz,elem)*da+UV(1,nz-1,elem)*db+UV(1,nz-2,elem)*dg
	vmean=UV(2,nz,elem)*da+UV(2,nz-1,elem)*db+UV(2,nz-2,elem)*dg
	end if
      end if
          temp=Av(nz,elem)/(Z(nz-1)-Z(nz))*elem_area(elem)
          uvertAB(1,nz)= -umean*w 
	  uvert(1,nz)= temp*(UV(1,nz-1,elem)-UV(1,nz,elem))
	  uvertAB(2,nz)= -vmean*w 
	  uvert(2,nz)= temp*(UV(2,nz-1,elem)-UV(2,nz,elem))
   END DO
   ! Wind stress and bottom drag
   if (i_vert_visc) then
          uvert(1,1)=0d0
          uvert(2,1)=0d0
          uvert(1,nlevels(elem))=0d0
          uvert(2,nlevels(elem))=0d0
   else
	  uvert(1,1)= stress_surf(1,elem)*elem_area(elem)/density_0
	  uvert(2,1)= stress_surf(2,elem)*elem_area(elem)/density_0
	  friction=C_d*sqrt(UV(1,nlevels(elem)-1,elem)**2+ &
	               UV(2,nlevels(elem)-1,elem)**2)*elem_area(elem)
          
	  uvert(1,nlevels(elem))=friction*UV(1,nlevels(elem)-1,elem)
	  uvert(2,nlevels(elem))=friction*UV(2,nlevels(elem)-1,elem)
    endif
          w=sum(Wvel(1,elnodes))*elem_area(elem)/3.0_WP 
	  uvertAB(1,1)= -w*UV(1,1,elem)
	  uvertAB(2,1)= -w*UV(2,1,elem) 
	  ! + sign here because it is subtracted!
	  uvertAB(1,nlevels(elem))=0.0_WP 
	  uvertAB(2,nlevels(elem))=0.0_WP 
   DO nz=1,nlevels(elem)-1
      UV_rhs(1,nz,elem)=UV_rhs(1,nz,elem)+(uvert(1,nz)-uvert(1,nz+1))/(zbar(nz)-zbar(nz+1))
      UV_rhs(2,nz,elem)=UV_rhs(2,nz,elem)+(uvert(2,nz)-uvert(2,nz+1))/(zbar(nz)-zbar(nz+1))
      UV_rhsAB(1,nz,elem)=UV_rhsAB(1,nz,elem)+(uvertAB(1,nz)-uvertAB(1,nz+1))/(zbar(nz)-zbar(nz+1))
      UV_rhsAB(2,nz,elem)=UV_rhsAB(2,nz,elem)+(uvertAB(2,nz)-uvertAB(2,nz+1))/(zbar(nz)-zbar(nz+1))
   END DO
END DO
END SUBROUTINE momentum_vert_adv_visc
! ===================================================================================
SUBROUTINE momentum_vert_visc
! 
! Vertical viscosity, explicit 
! 
USE o_PARAM
USE o_ARRAYS
USE o_MESH
USE g_PARSUP
use g_comm_auto
USE g_CONFIG
IMPLICIT NONE

integer          :: elem, elnodes(3), nz 
real(kind=WP)    :: friction 
real(kind=WP)    :: uvert(2,nl)
real(kind=WP)    :: s_coeff

! =======================
! Vertical viscosity
! =======================
 DO elem=1, myDim_elem2D         !! m=1, myDim_elem2D        !! P (d) elem=1,elem2D

    if (mom_adv/=4) then
       s_coeff=elem_area(elem)
    else
       s_coeff=dt
    end if
                                 !! elem=myList_elem2D(m)
   elnodes=elem2D_nodes(:,elem)
   DO nz=2, nlevels(elem)-1
          uvert(1,nz)= Av(nz,elem)*(UV(1,nz-1,elem)-UV(1,nz,elem))/(Z(nz-1)-Z(nz))
	  uvert(2,nz)= Av(nz,elem)*(UV(2,nz-1,elem)-UV(2,nz,elem))/(Z(nz-1)-Z(nz))
   END DO
   if (i_vert_visc) then
          uvert(1,1)=0d0
          uvert(2,1)=0d0
          uvert(1,nlevels(elem))=0d0
          uvert(2,nlevels(elem))=0d0
   else
   ! Wind stress and bottom drag
          uvert(1,1)= stress_surf(1,elem)*elem_area(elem)/density_0
	  uvert(2,1)= stress_surf(2,elem)*elem_area(elem)/density_0
	  friction=C_d*sqrt(UV(1,nlevels(elem)-1,elem)**2+ &
	               UV(2,nlevels(elem)-1,elem)**2)*elem_area(elem)
          
	  uvert(1,nlevels(elem))=friction*UV(1,nlevels(elem)-1,elem)
	  uvert(2,nlevels(elem))=friction*UV(2,nlevels(elem)-1,elem)
	  ! + sign here because it is subtracted!
	 
   endif
   DO nz=1, nlevels(elem)-1
      UV_rhs(1,nz,elem)=UV_rhs(1,nz,elem)+(uvert(1,nz)-uvert(1,nz+1))/(zbar(nz)-zbar(nz+1))*s_coeff
      UV_rhs(2,nz,elem)=UV_rhs(2,nz,elem)+(uvert(2,nz)-uvert(2,nz+1))/(zbar(nz)-zbar(nz+1))*s_coeff
   END DO
 END DO

END SUBROUTINE momentum_vert_visc
! ===================================================================================
