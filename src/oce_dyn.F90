! ===================================================================
! Contains routines needed for computations of dynamics.
! includes: solve_ssh, vert_vel, update_vel, compute_vel_nodes, compute_ssh_rhs, stiff_mat, 
!           test_ssh_rhs, impl_vert_visc, viscosity_filt2x
! ===================================================================
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
integer                         :: myrows
integer                         :: Pmode
real(kind=8)                    :: rinfo(20,20)
integer                         :: maxiter=2000
integer                         :: restarts=15
integer                         :: fillin=3
integer                         :: lutype=2
integer                         :: nrhs=1
real(kind=8)                    :: droptol=1.e-7
real(kind=8)                    :: soltol =1e-10  !1.e-10
logical, save                   :: lfirst=.true.
real(kind=WP), allocatable      :: arr_nod2D(:),arr_nod2D2(:,:),arr_nod2D3(:)
real(kind=8)                    :: cssh1,cssh2,crhs
integer                         :: i
Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB +PET_REPORT + PET_QUIET+ PET_RCM+PET_PCBJ
if (lfirst) then   
   Pmode = Pmode+PET_STRUCT+PET_PMVALS + PET_PCASM+PET_OVL_2 !+PET_PCBJ+PET_ILU
   lfirst=.false.
end if
call PETSC_S(Pmode, 1, ssh_stiff%dim, ssh_stiff%nza, myrows, &
     maxiter,  & 
     restarts, &
     fillin,   &
     droptol,  &  
     soltol,   &
     part, ssh_stiff%rowptr, ssh_stiff%colind, ssh_stiff%values, &
     ssh_rhs, d_eta, &
     rinfo, MPI_COMM_WORLD)

#elif defined(PARMS)

  use iso_c_binding, only: C_INT, C_DOUBLE
  implicit none
#include "fparms.h"
logical, save        :: lfirst=.true.
integer(kind=C_INT)  :: ident
integer(kind=C_INT)  :: n3, reuse
integer(kind=C_INT)  :: maxiter, restart, lutype, fillin
real(kind=C_DOUBLE)  :: droptol, soltol
integer :: n

interface
   subroutine psolver_init(ident, SOL, PCGLOB, PCLOC, lutype, &
        fillin, droptol, maxiter, restart, soltol, &
        part, rowptr, colind, values, reuse, MPI_COMM) bind(C)
     use iso_c_binding, only: C_INT, C_DOUBLE
     integer(kind=C_INT) :: ident, SOL, PCGLOB, PCLOC, lutype, &
                            fillin,  maxiter, restart, &
                            part(*), rowptr(*), colind(*), reuse, MPI_COMM
     real(kind=C_DOUBLE) :: droptol,  soltol, values(*)
   end subroutine psolver_init
end interface
interface
   subroutine psolve(ident, ssh_rhs, zero_r, d_eta, zero_i) bind(C)

     use iso_c_binding, only: C_INT, C_DOUBLE
     integer(kind=C_INT) :: ident, zero_i
     real(kind=C_DOUBLE) :: zero_r, ssh_rhs(*), d_eta(*)

   end subroutine psolve
end interface

ident=1
maxiter=2000
restart=15
fillin=2 !3
lutype=1 !2
droptol=1.e-8
soltol=1.e-10
reuse=0

if (lfirst) then
   ! Set SOLCG for CG solver (symmetric, positiv definit matrices only!!)
   !     SOLBICGS for BiCGstab solver (arbitrary matrices)
   call psolver_init(ident, SOLCG, PCRAS, PCILUK, lutype, &
        fillin, droptol, maxiter, restart, soltol, &
        part-1, ssh_stiff%rowptr(:)-ssh_stiff%rowptr(1), &
        ssh_stiff%colind-1, ssh_stiff%values, reuse, MPI_COMM_WORLD)
   lfirst=.false.
end if

   call psolve(ident, ssh_rhs, real(0., C_DOUBLE), d_eta, 0)
#endif

call exchange_nod(d_eta)

end subroutine solve_ssh
! ===================================================================================
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
Do n=1, myDim_nod2D
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
DO edge=1, myDim_edge2D
   enodes=edges(:,edge)   
   el=edge_tri(:,edge)
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   if (el(2)>0) then
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

   if (el(2)>0) then
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
Do n=1, myDim_nod2D
   DO nz=nl-1,1,-1
      Wvel(nz,n)=Wvel(nz,n)+Wvel(nz+1,n)
      if (Fer_GM) then 
         fer_Wvel(nz,n)=fer_Wvel(nz,n)+fer_Wvel(nz+1,n)
      end if
   END DO
END DO

Do n=1, myDim_nod2D
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
!==========================================================================
SUBROUTINE update_vel
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
integer       :: elem, elnodes(3), nz, m
real(kind=WP) :: eta(3) 
real(kind=WP) :: Fx, Fy

 DO elem=1, myDim_elem2D
    elnodes=elem2D_nodes(:,elem)
    eta=-g*theta*dt*d_eta(elnodes)
    Fx=sum(gradient_sca(1:3,elem)*eta)
    Fy=sum(gradient_sca(4:6,elem)*eta)
    DO nz=1, nlevels(elem)-1
     UV(1,nz,elem)= UV(1,nz,elem) + UV_rhs(1,nz,elem) + Fx
     UV(2,nz,elem)= UV(2,nz,elem) + UV_rhs(2,nz,elem) + Fy
    END DO
 END DO
eta_n=eta_n+d_eta
call exchange_elem(UV)
end subroutine update_vel
!==========================================================================
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
         if (nlevels(elem)-1<nz) cycle
            tvol=tvol+elem_area(elem)
            tx=tx+UV(1,nz,elem)*elem_area(elem)
            ty=ty+UV(2,nz,elem)*elem_area(elem)
      END DO
      Unode(1,nz,n)=tx/tvol
      Unode(2,nz,n)=ty/tvol
   END DO
END DO
call exchange_nod(Unode)
end subroutine compute_vel_nodes
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

integer       :: ed, el(2), enodes(2),  nz
real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2 

ssh_rhs=0.0_WP
DO ed=1, myDim_edge2D
   enodes=edges(:,ed)
   el=edge_tri(:,ed)   
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,ed)
   deltaY1=edge_cross_dxdy(2,ed)
   if (el(2)>0)  then
      deltaX2=edge_cross_dxdy(3,ed)
      deltaY2=edge_cross_dxdy(4,ed)
   end if
   DO nz=1, nlevels(el(1))-1
      c1=c1+((UV(2,nz,el(1))+alpha*UV_rhs(2,nz,el(1)))*deltaX1- &
         (UV(1,nz,el(1))+alpha*UV_rhs(1,nz,el(1)))*deltaY1)*(zbar(nz)-zbar(nz+1))
   END DO
 
   if (el(2)>0) then
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
! ssh_rhs(1:myDim_nod2D)=eta_n(1:myDim_nod2D)*area(1,:)/dt+ssh_rhs(1:myDim_nod2D)
!! P-E should be just added as it is to the rhs.  
! call exchange_nod(ssh_rhs)
! call test_ssh_rhs ! delete me
END subroutine compute_ssh_rhs
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
implicit none

integer                             :: n, n1, n2, i, j,  row, ed
integer                             :: enodes(2), elnodes(3), el(2)
integer                             :: elem, npos(3), offset, nini, nend
real(kind=WP)                       :: dmean, ff, factor, a 
real(kind=WP)                       :: fx(3), fy(3), ax, ay
integer, allocatable                :: n_num(:), n_pos(:,:), pnza(:), rpnza(:)
integer, allocatable                :: mapping(:)
logical                             :: flag
character*10                        :: mype_string,npes_string
character*80                        :: dist_mesh_dir,file_name
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
   n1=edges(1,n)
   n2=edges(2,n)
   if (n1<=myDim_nod2D) then
      n_pos(n_num(n1)+1,n1)=n2
      n_num(n1)=n_num(n1)+1
   end if
   if (n2<=myDim_nod2D) then
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
factor = g*dt*alpha*theta
DO ed=1,myDim_edge2D   !! Attention
   enodes=edges(:,ed)
   el=edge_tri(:,ed)
   DO j=1,2 
      row=enodes(j)
      if (row>myDim_nod2D) cycle    !! Attention
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
         if (elem<1) cycle
         elnodes=elem2D_nodes(:,elem)
         call elem_center(elem, ax, ay)
         ay=cos(ay)
         if (cartesian) ay=1.0_WP
         dmean=-zbar(nlevels(elem))
         dmean=0.5_WP*dmean*r_earth
         a=dmean*(coord_nod2D(2,elnodes(3))-coord_nod2D(2,elnodes(1)))
         fx(2)=a
         fx(1)=-a
         a=-dmean*(coord_nod2D(2,elnodes(2))-coord_nod2D(2,elnodes(1)))
         fx(3)=a
         fx(1)=fx(1)-a
         a=(coord_nod2D(1,elnodes(3))-coord_nod2D(1,elnodes(1)))
         if (a>cyclic_length/2.) a=a-cyclic_length
         if (a<-cyclic_length/2.) a=a+cyclic_length
         a=-ay*dmean*a
         fy(2)=a
         fy(1)=-a
         a=(coord_nod2D(1,elnodes(2))-coord_nod2D(1,elnodes(1)))
         if (a>cyclic_length/2.) a=a-cyclic_length
         if (a<-cyclic_length/2.) a=a+cyclic_length
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
         fy = fy*(edge_cross_dxdy(2*i-1,ed)) &
            - fx*(edge_cross_dxdy(2*i  ,ed))
   
         if(i==2) fy=-fy
         if(j==2) fy=-fy
         ! 
         ! In the computation above, I have used rules from ssh_rhs (where it is 
         ! on the rhs. So the sign is changed in the expression below.
         !
         npos=n_num(elnodes)
         SSH_stiff%values(npos)=SSH_stiff%values(npos)- fy*factor
      END DO
   END DO
END DO 
! Mass matrix part
DO row=1, myDim_nod2D
   offset=ssh_stiff%rowptr(row)
   SSH_stiff%values(offset)=SSH_stiff%values(offset)+ &
                            area(1,row)/dt
END DO
deallocate(n_pos,n_num)

! =================================
! Global contiguous numbers:
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
subroutine test_ssh_rhs
USE g_PARSUP
use o_MESH
use o_ARRAYS
IMPLICIT NONE

real(kind=WP) :: c1,c2
c1=sum(ssh_rhs(1:mydim_nod2d))
c2=0d0
call MPI_Allreduce (c1, c2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
if (mype==0) write(*,*) 'DBG:sum(ssh_rhs)=',c2
end subroutine test_ssh_rhs
! ===================================================================================
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
   ! Update from the vertical advection
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

! Update from the vertical advection
   wu=sum(Wvel_i(nzmax-1, elnodes))/3.
   a(nzmax-1)=a(nzmax-1)+min(0._WP, wu)*zinv
   b(nzmax-1)=b(nzmax-1)+max(0._WP, wu)*zinv

! The first row
   zinv=1.0_WP*dt/(zbar(1)-zbar(2))
   c(1)=-Av(2,elem)/(Z(1)-Z(2))*zinv
   a(1)=0.0_WP
   b(1)=-c(1)+1.0_WP
! Update from the vertical advection
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
SUBROUTINE viscosity_filt2x
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_config
USE g_comm_auto
IMPLICIT NONE

real(kind=8)               :: u1, v1, tau_inv, s
integer                    :: ed, el(2), nz
real(kind=8), allocatable  :: UV_c(:,:,:), UV_f(:,:,:)
 
! Filter is applied twice. It should be approximately 
! equivalent to biharmonic operator with the coefficient
! (tau_c/day)a^3/9. Scaling inside is found to help 
! with smoothness in places of mesh transition. *(it makes a^3 from a^4) 
ed=myDim_elem2D+eDim_elem2D
allocate(UV_c(2,nl-1,ed)) ! to store the filtered velocity
allocate(UV_f(2,nl-1,ed)) ! to store the contributions into the RHS


UV_c=0.0_8
UV_f=0.0_8
tau_inv=dt*tau_c/3600.0/24.0     ! SET IT experimentally 
  
DO ed=1, myDim_edge2D+eDim_edge2D
   if(myList_edge2D(ed)>edge2D_in) cycle
   el=edge_tri(:,ed)
   DO nz=1,minval(nlevels(el))-1
      u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
      v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))

      UV_c(1,nz,el(1))=UV_c(1,nz,el(1))-u1
      UV_c(1,nz,el(2))=UV_c(1,nz,el(2))+u1
      UV_c(2,nz,el(1))=UV_c(2,nz,el(1))-v1
      UV_c(2,nz,el(2))=UV_c(2,nz,el(2))+v1
   END DO 
END DO
 
! ============ 
! Contribution from boundary edges (Dirichlet boundary conditions)
! ============
! DO ed=1, myDim_edge2D+eDim_edge2D
!    if(myList_edge2D(ed)<=edge2D_in) cycle
!    el=edge_tri(:, ed)
!    DO  nz=1, nlevels(el(1))-1
!        UV_c(1,nz,el(1))=UV_c(1,nz,el(1))-2.0_WP*UV(1,nz,el(1))
!        UV_c(2,nz,el(1))=UV_c(2,nz,el(1))-2.0_WP*UV(2,nz,el(1))
!    END DO
! END DO

Do ed=1,myDim_elem2D
   Do nz=1,nlevels(ed)-1
       UV_c(1,nz,ed)=-UV_c(1,nz,ed)*tau_inv*sqrt(scale_area/elem_area(ed))
       UV_c(2,nz,ed)=-UV_c(2,nz,ed)*tau_inv*sqrt(scale_area/elem_area(ed))
   END DO
end do

call exchange_elem(UV_c)
! call exchange_elem(UV_c(1,:,:))
! call exchange_elem(UV_c(2,:,:))

DO ed=1, myDim_edge2D+eDim_edge2D
   if(myList_edge2D(ed)>edge2D_in) cycle
   el=edge_tri(:,ed)
   DO nz=1,minval(nlevels(el))-1 
      u1=(UV_c(1,nz,el(1))-UV_c(1,nz,el(2)))
      v1=(UV_c(2,nz,el(1))-UV_c(2,nz,el(2)))

      UV_f(1,nz,el(1))=UV_f(1,nz,el(1))-u1
      UV_f(1,nz,el(2))=UV_f(1,nz,el(2))+u1
      UV_f(2,nz,el(1))=UV_f(2,nz,el(1))-v1
      UV_f(2,nz,el(2))=UV_f(2,nz,el(2))+v1
   END DO 
END DO
 
DO ed=1, myDim_elem2D
   DO nz=1, nlevels(ed)-1 
      u1=sqrt(UV_f(1,nz,ed)**2+UV_f(2,nz,ed)**2)+1.e-5
      v1=sqrt(UV(1,nz,ed)**2+UV(2,nz,ed)**2)
   ! we limit the maximum contribution from the filter such, that the update is less than the N (N=2 currently) times velocity
   ! this is done to force the CFL, which is otherwise exceeded in some points
   ! some other criteria is welcome (i.e. like computing the eigenvalues from filtering)
      UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+UV_f(1,nz,ed)*min(1.0_WP, 2.0_WP*v1/u1)
      UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+UV_f(2,nz,ed)*min(1.0_WP, 2.0_WP*v1/u1)
   END DO 
END DO
deallocate(UV_f, UV_c)    
end subroutine viscosity_filt2x
