module o_split
implicit none
save
real*8,allocatable,dimension(:,:) :: vel_btr,rhs_btr,vel_mean
real*8,allocatable,dimension(:) :: depth0,elev,elev_mean
integer :: num_btr_steps=100
real*8 :: dt_btr

!-----------------------------------------------------
contains
subroutine init_split
use g_parsup
use o_mesh
use g_config
implicit none
integer :: nod_size,elem_size
integer i
nod_size=myDim_nod2D+eDim_nod2D
elem_size=myDim_elem2D+eDim_elem2D
allocate(elev(nod_size),elev_mean(nod_size),depth0(elem_size))
allocate(vel_btr(elem_size,2),rhs_btr(elem_size,2),vel_mean(elem_size,2))
dt_btr=2d0*dt/dble(num_btr_steps)
end subroutine init_split
!------------------------------------------------------
subroutine bar_split
use g_parsup
use g_config
use o_mesh
use o_arrays
use g_comm_auto
use o_param
!Here baroclinic split is treated as a computational trick.
!eta_new=A*eta_old, to make it stable we need to filter out components that we cannot possibly resolve (i.e. barotropic motion)
!In case of implicit time step this done automatically
!Here fast mode is filtered via averaging
implicit none
Real*8  :: ff, pre(3), Fx, Fy,c1,c2
real*8  :: deltaX1, deltaY1, deltaX2, deltaY2
integer :: n,nz,elem,elnodes(3),el(2),enodes(2),ed
real*8 :: t1,t2,t3,t4,t5
!1. Find baroclinic velocity and its RHS
t1=MPI_Wtime()
vel_btr=0d0
rhs_btr=0d0
do elem=1,myDim_elem2D
   elnodes=elem2D_nodes(:, elem)
!   depth0(elem)=sum(eta_n(elnodes))/3d0 !  TODO  I'm not sure that's entirely correct // KKV
   depth0(elem)=0d0 !  TODO We don't really use free surface, so it seems like a more logical choice // KKV
   ff=sum(coriolis_node(elnodes))/3d0
   do nz=1,nlevels(elem)-1
      vel_btr(elem,1)=vel_btr(elem,1)-(zbar(nz+1)-zbar(nz))*UV(1,nz,elem)
      vel_btr(elem,2)=vel_btr(elem,2)-(zbar(nz+1)-zbar(nz))*UV(2,nz,elem)
      depth0(elem)=depth0(elem) - zbar(nz+1) + zbar(nz)
   enddo
   vel_btr(elem,:)=vel_btr(elem,:)/depth0(elem)
      !RHS is w/o Coriolis force, so it can be split up. grad(SSH) is also removed(for debugging period).
   do nz=1,nlevels(elem)-1
      pre=-g*eta_n(elnodes)
      Fx=sum(gradient_sca(1:3,elem)*pre)
      Fy=sum(gradient_sca(4:6,elem)*pre)
      UV_rhs(1,nz,elem)=UV_rhs(1,nz,elem)-dt*vel_btr(elem,2)*ff-Fx*dt
      UV_rhs(2,nz,elem)=UV_rhs(2,nz,elem)+dt*vel_btr(elem,1)*ff-Fy*dt
      UV(1,nz,elem)=UV(1,nz,elem)-vel_btr(elem,1)
      UV(2,nz,elem)=UV(2,nz,elem)-vel_btr(elem,2)
      rhs_btr(elem,1)=rhs_btr(elem,1)-2d0*(zbar(nz+1)-zbar(nz))*(UV_rhs(1,nz,elem)+UV(1,nz,elem))/dble(num_btr_steps)
      rhs_btr(elem,2)=rhs_btr(elem,2)-2d0*(zbar(nz+1)-zbar(nz))*(UV_rhs(2,nz,elem)+UV(2,nz,elem))/dble(num_btr_steps)
   enddo
   rhs_btr(elem,:)=rhs_btr(elem,:)/depth0(elem)
   
   do nz=1,nlevels(elem)-1
      UV_rhs(1,nz,elem)=UV_rhs(1,nz,elem)-rhs_btr(elem,1)*real(num_btr_steps)*5e-1
      UV_rhs(2,nz,elem)=UV_rhs(2,nz,elem)-rhs_btr(elem,2)*real(num_btr_steps)*5e-1
   enddo
enddo
vel_mean=vel_btr
elev_mean=eta_n
elev=eta_n
ssh_rhs=0.0_WP
t2=MPI_Wtime()

!2. Advance b. trop. velocity and elevation
!Forward Euler scheme is used
do n=1, num_btr_steps
   do elem=1,myDim_elem2D
      elnodes=elem2D_nodes(:, elem)
!      depth0(elem)=depth0(elem)+sum(ssh_rhs(elem2D_nodes(:,elem)))/3d0!TODO Should we do it?
      pre=-g*elev(elnodes)
      Fx=sum(gradient_sca(1:3,elem)*pre)
      Fy=sum(gradient_sca(4:6,elem)*pre)
      ff=sum(coriolis_node(elnodes))/3d0
      c1=vel_btr(elem,1)
      c2=vel_btr(elem,2)

      vel_btr(elem,1)=c1+rhs_btr(elem,1)+( ff*c2+Fx )*dt_btr
      vel_btr(elem,2)=c2+rhs_btr(elem,2)+(-ff*c1+Fy )*dt_btr
   enddo
   call exchange_elem(vel_btr(:,1))
   call exchange_elem(vel_btr(:,2))
   ssh_rhs=0.0_WP
   DO ed=1, myDim_edge2D
     enodes=edges(:,ed)
     el=edge_tri(:,ed)
     deltaX1=edge_cross_dxdy(1,ed) ; deltaY1=edge_cross_dxdy(2,ed)
     c1=0.0_WP ; c2=0.0_WP
     if(el(2)>0)  then
       deltaX2=edge_cross_dxdy(3,ed) ; deltaY2=edge_cross_dxdy(4,ed)
     end if
                 c1=( vel_btr(el(1),2)*deltaX1 - vel_btr(el(1),1)*deltaY1)*depth0(el(1))
     if(el(2)>0) c2=(-vel_btr(el(2),2)*deltaX2 + vel_btr(el(2),1)*deltaY2)*depth0(el(2))
     ssh_rhs(enodes(1))=ssh_rhs(enodes(1))+(c1+c2)*dt_btr/area(1,enodes(1))
     ssh_rhs(enodes(2))=ssh_rhs(enodes(2))-(c1+c2)*dt_btr/area(1,enodes(2))
   END DO
   call exchange_nod(ssh_rhs(:))
   elev=elev+ssh_rhs !TODO Check the sign
   
   vel_mean=vel_mean+vel_btr
   elev_mean=elev_mean+elev

enddo
t3=MPI_Wtime()
!3.Average barotropic vars, advance U_n,V_n
elev_mean=elev_mean/dble(num_btr_steps+1)
vel_mean=vel_mean/dble(num_btr_steps+1)
d_eta=elev_mean-eta_n
do elem=1,myDim_elem2D
 do nz=1,nlevels(elem)-1
     UV_rhs(1,nz,elem)=UV_rhs(1,nz,elem)+vel_mean(elem,1)
     UV_rhs(2,nz,elem)=UV_rhs(2,nz,elem)+vel_mean(elem,2)
 enddo
enddo


t4=MPI_Wtime()
!if (mype==0) write(*,*) 'SSH DBG', t2-t1,t3-t2,t4-t3
end subroutine bar_split
end module o_split
