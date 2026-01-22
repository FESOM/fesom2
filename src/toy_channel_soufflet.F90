MODULE Toy_Channel_Soufflet
  use mod_mesh
  USE o_ARRAYS
  USE o_PARAM
  USE g_PARSUP
  USE g_config

  implicit none
  SAVE 
  private
  public :: relax_zonal_vel, relax_zonal_temp, compute_zonal_mean_ini, compute_zonal_mean, initial_state_soufflet, energy_out_soufflet, soufflet_forc_update

  real(kind=WP), allocatable  :: zvel(:,:), Uclim(:,:)
  real(kind=WP), allocatable  :: ztem(:,:)
  integer,       allocatable  :: bpos(:)
  real(kind=WP), allocatable  :: znum(:,:)
  integer           ::   soufflet_forc_update      =10
  real(kind=WP)     ::   tau_inv     =1.0/50.0/24.0/3600.0 ! piston velocity for zonal restoring [1/s] (1/50 days in Stouflet et al. 2016)
  real(kind=WP)     ::   lat0        =0.0                  ! the Southern boundary of the mesh (in radians)
  real(kind=WP)     ::   ysize       =2000000.0            ! the meridional lenght of the channel [m]
  real(kind=WP)     ::   xsize       =90018410.49779853    ! the zonal lenght of the channel [m] !4.5*pi*r_earth=90018410.49779853
  integer           ::   nybins      =100                  ! number of meridional bins for zonal restoring [m]
  real(kind=WP)     ::   dy                                ! ysize/nybins/r_earth (size of bins in radians)
  real(kind=WP)     ::   Ljet        =1600000.0   ! m
  real(kind=WP)     ::   rhomax      =27.75       ! background density 
  real(kind=WP)     ::   Sb          =9.8e-6      ! background density vertical profile
  real(kind=WP)     ::   zsize       =4000.0      ! m  The depth
  real(kind=WP)     ::   drho_No     =1.41
  real(kind=WP)     ::   drho_So     =1.4
  real(kind=WP)     ::   z_No        =-400.0
  real(kind=WP)     ::   z_So        =-1000.0
  real(kind=WP)     ::   dz_No       =300.0
  real(kind=WP)     ::   dz_So       =700.0
  real(kind=WP)     ::   drhosurf_No =0.0
  real(kind=WP)     ::   drhosurf_So =1.5
  real(kind=WP)     ::   zsurf       =-300.0
!
!--------------------------------------------------------------------------------------------
!
  contains
!
!--------------------------------------------------------------------------------------------
!
subroutine relax_zonal_vel(mesh)
  implicit none
  integer        :: elem,  nz, nn, nn1
  real(kind=WP)  :: a, yy, uzon 
  type(t_mesh), intent(in)     , target :: mesh
#include  "associate_mesh.h"

  DO elem=1, myDim_elem2D
     ! ========
     ! Interpolation
     ! ========
     yy=sum(coord_nod2D(2,elem2D_nodes(:,elem)))/3.0-lat0
     a=0  !interpolation coefficient
     if (yy<dy/2) then  ! southward from the center of the first bin
        nn=1
        nn1=1
     else 
        nn=floor(yy/dy-0.5)+1
        nn1=nn+1
        if (nn1>100) nn1=nn  ! northward of the center of the last bin 
                             ! Linear interpolation (nearest if close to boundary)
        a=yy/dy+0.5-real(nn)
     end if
     ! ========
     DO nz=1, nlevels(elem)-1
        Uzon              = (1.0-a)*zvel(nz,nn)+a*zvel(nz,nn1)
        UV_rhs(1,nz,elem) = UV_rhs(1,nz,elem)+dt*tau_inv*(Uclim(nz,elem)-Uzon)
     END DO
  END DO
end subroutine relax_zonal_vel
!==========================================================================
subroutine relax_zonal_temp(mesh)
  implicit none
  integer                           :: n, nz, nn, nn1
  real(kind=WP)                     :: yy, a, Tzon
  type(t_mesh), intent(in) , target :: mesh
#include  "associate_mesh.h"

  do n=1, myDim_nod2D
     yy=coord_nod2D(2,n)-lat0
     a=0 
    if (yy<dy/2) then  ! southward from the center of the first bin
       nn=1
       nn1=1
    else 
       nn=floor(yy/dy-0.5)+1
       nn1=nn+1
    if (nn1>100) nn1=nn  ! northward of the center of the last bin 
                         ! Linear interpolation (nearest if close to boundary)
       a=yy/dy+0.5-nn
    end if
    do nz=1, nlevels_nod2D(n)-1
       Tzon=(1.0-a)*ztem(nz,nn)+a*ztem(nz,nn1)
       tr_arr(nz,n,1)=  tr_arr(nz,n,1)+dt*tau_inv*(Tclim(nz,n)-Tzon)
    end do
  end do
end subroutine relax_zonal_temp
!==========================================================================
subroutine compute_zonal_mean_ini(mesh)
  implicit none 
  real(kind=8)                      :: ymean, Ly
  integer                           :: elem, nz, m, elnodes(3)
  real(kind=8), allocatable         :: zvel1D(:), znum1D(:)
  type(t_mesh), intent(in) , target :: mesh
#include  "associate_mesh.h"

Ly=ysize/r_earth       ! The meridional lenght in radians
dy=Ly/real(nybins)

allocate(bpos(myDim_elem2D))
allocate(zvel(nl-1, nybins), ztem(nl-1,nybins), znum(nl-1,nybins))
! =======
! Only channel geometry with flat bottom is supported
! =======
!
 
znum=0.0
zvel=0.0
ztem=0.0
! find element's  position in bins:
 DO elem=1,myDim_elem2D
    elnodes=elem2D_nodes(:,elem)
    ymean=sum(coord_nod2D(2,elnodes))/3.0
    bpos(elem)=floor((ymean-lat0)/dy)+1
 END DO
 if(maxval(bpos)>100) write(*,*) mype, 'BPOS1'
 if(minval(bpos)<1) write(*,*) mype, 'BPOS2' 
   
 DO elem=1,myDim_elem2D
    if(elem2D_nodes(1,elem)>myDim_nod2D) cycle
    znum(:,bpos(elem))=znum(:,bpos(elem))+1.0
 END DO
! =============
! Communication:
! =============
 allocate(zvel1D(nybins*(nl-1)), znum1D(nybins*(nl-1)))
 DO m=1,nybins
    zvel1D((m-1)*(nl-1)+1:m*(nl-1))=znum(:,m)
 END DO
    znum1D=0.0
 call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
 call MPI_AllREDUCE( zvel1D, znum1D, nybins*(nl-1), MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr) 
 ! fill them back in
  DO m=1,nybins
    znum(:,m)=znum1D((m-1)*(nl-1)+1:m*(nl-1))
  END DO
  deallocate(znum1d, zvel1d)      
 ! ==========
 ! Depending on partitioning, znum may contain 0 at some places. Take care that 
 ! no division by 0 is occurring 
end subroutine compute_zonal_mean_ini
!==========================================================================
subroutine compute_zonal_mean(mesh)
  implicit none 
  integer       :: elem, nz, m, elnodes(3)
  real(kind=8), allocatable  :: zvel1D(:), znum1D(:)
  type(t_mesh), intent(in) , target :: mesh
#include  "associate_mesh.h"


 ztem=0.
 zvel=0.
 DO elem=1,myDim_elem2D
    if(elem2D_nodes(1,elem)>myDim_nod2D) cycle
    Do nz=1,nlevels(elem)-1
    ztem(nz,bpos(elem))=ztem(nz,bpos(elem))+sum(tr_arr(nz,elem2D_nodes(:,elem),1))/3.0_8
    zvel(nz,bpos(elem))=zvel(nz,bpos(elem))+UV(1,nz,elem)
    END DO
 END DO
! =============
! Communication:
! =============
! Velocity:
! put vel result into 1D array
 allocate(zvel1D(nybins*(nl-1)), znum1D(nybins*(nl-1)))
 DO m=1,nybins
    zvel1D((m-1)*(nl-1)+1:m*(nl-1))=zvel(:,m)
 END DO
    znum1D=0.0
 call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
 call MPI_AllREDUCE( zvel1D, znum1D, nybins*(nl-1), MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr) 
 ! fill in back
  DO m=1,nybins
    zvel(:,m)=znum1D((m-1)*(nl-1)+1:m*(nl-1))
  END DO
! temperature
 DO m=1,nybins
    zvel1D((m-1)*(nl-1)+1:m*(nl-1))=ztem(:,m)
 END DO
    znum1D=0.0
 call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
 call MPI_AllREDUCE( zvel1D, znum1D, nybins*(nl-1), MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr) 
 ! fill in back
 DO m=1,nybins
    ztem(:,m)=znum1D((m-1)*(nl-1)+1:m*(nl-1))
 END DO
 ! zvel = mean zonal velocity 
  zvel=zvel/(znum+0.001)
 ! ztem = mean zonal temperature
  ztem=ztem/(znum+0.001)
  deallocate(znum1d, zvel1d)      
 ! ==========
 !if(mype==0) then
 !open(20, file='zonal.out')
 !   write(20,'(e13.5)') zvel, ztem
 !close(20)
 !endif

end subroutine compute_zonal_mean
! ====================================================================================
subroutine initial_state_soufflet(mesh)
 ! Profiles Soufflet 2016 (OM)
  implicit none
  type(t_mesh), intent(in) , target  :: mesh
  integer                            :: n, nz, elnodes(3)
  real(kind=8)                       :: dst, yn, Fy, Lx
! real(kind=8)                       :: Ljet,rhomax,Sb, drho_No, drho_So
! real(kind=8)                       :: z_No, z_So,dz_No,dz_So, drhosurf_No, drhosurf_So, zsurf 
  real(kind=8)                       :: d_No(mesh%nl-1), d_So(mesh%nl-1), rho_No(mesh%nl-1), rho_So(mesh%nl-1)
#include  "associate_mesh.h"

  dy=ysize/nybins/r_earth

  ! Default values
  stress_surf   = 0.0
  heat_flux     = 0.0_8
  tr_arr(:,:,2) = 35.0_8
  Ssurf         = tr_arr(1,:,1)
  water_flux    = 0.0_8
  relax2clim    = 0.0_8

! Have to set density_0=1028._WP in oce_modules.F90
! ========
! The north and south density profiles
! ========
DO nz=1,nl-1
    d_No(nz)=z_No+(Z(nz)-z_No)*sqrt(1+0.5*(((Z(nz)-z_No)+abs(Z(nz)-z_No))/1.3/dz_No)**2)
    d_So(nz)=z_So+(Z(nz)-z_So)*sqrt(1+0.5*(((Z(nz)-z_So)+abs(Z(nz)-z_So))/1.3/dz_So)**2)
end DO

DO nz=1,nl-1
    rho_No(nz)=rhomax-Sb*(Z(nz)+zsize)-0.5*drho_No*(1+tanh((d_No(nz)-z_No)/dz_No))- & 
        1.0/(2*tanh(1.0))*drhosurf_No*(1+tanh((zsurf-Z(nz))/zsurf))
    
    rho_So(nz)=rhomax-Sb*(Z(nz)+zsize)-0.5*drho_So*(1+tanh((d_So(nz)-z_So)/dz_So))- & 
        1.0/(2*tanh(1.0))*drhosurf_So*(1+tanh((zsurf-Z(nz))/zsurf))
end DO
! We need to convert density profiles to temperature profiles
! T=T_0-(rho-rho_0)/alpha/rho_0; alpha*rho_0=0.00025*1030 kg/K/m^3
! T_0=10K

rho_No=10.0-(rho_No-rhomax)/(0.00025_WP*density_0)        ! T_0 it the bottom temperature
rho_So=10.0-(rho_So-rhomax)/(0.00025_WP*density_0)

! ========
!  2D profile:
! ========

do n=1, myDim_nod2D+eDim_nod2D
     ! The north and south profiles are combined as 
     dst=(coord_nod2D(2, n)-lat0)*r_earth
     yn=pi*(ysize/Ljet)*(dst/ysize-0.5)+pi/2.0
     if(yn<0) then
       Fy=1.0;
     else
        if(yn>pi) then
            Fy=0.0
        else
            Fy=1.0-(yn-sin(yn)*cos(yn))/pi;
        end if
     end if
     do nz=1, nlevels_nod2D(n)-1
        tr_arr(nz, n,1)=rho_So(nz)+(rho_No(nz)-rho_So(nz))*(1.0-Fy)
     end do
  end do
 
  ! ========
  ! Make consistent
  ! ========
  Tsurf=tr_arr(1,:,1)
  Tclim=tr_arr(:,:,1)
  ! ========
  ! add small perturbation:
  do n=1, myDim_nod2D+eDim_nod2D
     dst=(coord_nod2D(2, n)-lat0)*r_earth
     do nz=1, nlevels_nod2D(n)-1
        tr_arr(nz,n,1)=tr_arr(nz,n,1)-0.1*sin(2*pi*dst/ysize)*exp(2*Z(nz)/zsize) &
	       *(sin(8*pi*coord_nod2D(1,n)*r_earth/xsize)+ &
	       0.5*sin(3*pi*coord_nod2D(1,n)*r_earth/xsize))
     end do
  end do
  ! =======
  ! Compute geostrophically balanced flow
  ! =======
  write(*,*) mype, 'T', maxval(tr_arr(:,:,1)), minval(tr_arr(:,:,1))
  ! Redefine Coriolis (to agree with the Soufflet paper) 
  DO n=1,myDim_elem2D
  elnodes=elem2D_nodes(:,n)
  dst=(sum(coord_nod2D(2, elnodes))/3.0-lat0)*r_earth-ysize/2
  coriolis(n)=1.0e-4+dst*1.6e-11	 
  END DO
  write(*,*) mype, 'COR', maxval(coriolis*10000.0), minval(coriolis*10000.0) 
  DO n=1,myDim_elem2D 
  elnodes=elem2D_nodes(:,n)
  ! Thermal wind \partial_z UV(1,:,:)=(g/rho_0/f)\partial_y rho 
  DO nz=1,nlevels(n)-1
     d_No(nz)=(-(0.00025_WP*density_0)*g/density_0/coriolis(n))*sum(gradient_sca(4:6,n)*Tclim(nz, elnodes))
  !  d_N is used here as a placeholder
  ! -(sw_alpha*density_0) here is from the equation of state d\rho=-(sw_alpha*density_0) dT 
  END DO
  ! Vertical integration
  nz=nlevels(n)-1
  UV(1,nz,n)=d_No(nz)*(Z(nz)-zbar(nz+1)) 
  DO nz=nlevels(n)-2,1,-1
  UV(1,nz,n)=UV(1,nz+1,n)+d_No(nz+1)*(zbar(nz+1)-Z(nz+1))+d_No(nz)*(Z(nz)-zbar(nz+1)) 
  END DO 
  END DO
  allocate(Uclim(nl-1,myDim_elem2D))
  Uclim=UV(1,:,:)
  write(*,*) mype, 'Vel', maxval(UV(1,:,:)), minval(UV(1,:,:))
 END subroutine initial_state_soufflet
! ===============================================================================
subroutine energy_out_soufflet(mesh)
  implicit none 
  real(kind=8)                      :: tke(2), aux(2), ww, wwaux 
  integer                           :: elem, nz, m, elnodes(3), nybins
  real(kind=8), allocatable         :: zvel1D(:), znum1D(:)
  type(t_mesh), intent(in) , target :: mesh
#include  "associate_mesh.h"

 nybins=100
 zvel=0.
 DO elem=1,myDim_elem2D
    if(elem2D_nodes(1,elem)>myDim_nod2D) cycle
    Do nz=1,nlevels(elem)
    zvel(nz,bpos(elem))=zvel(nz,bpos(elem))+UV(1,nz,elem)
    END DO
 END DO
! =============
! Communication:
! =============
! Velocity:
! put vel result into 1D array
 allocate(zvel1D(nybins*(nl-1)), znum1D(nybins*(nl-1)))
 DO m=1,nybins
    zvel1D((m-1)*(nl-1)+1:m*(nl-1))=zvel(:,m)
 END DO
    znum1D=0.0
 call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
 call MPI_AllREDUCE( zvel1D, znum1D, nybins*(nl-1), MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr) 
 ! fill in back
  DO m=1,nybins
    zvel(:,m)=znum1D((m-1)*(nl-1)+1:m*(nl-1))
  END DO
 ! zvel = mean zonal velocity 
  zvel=zvel/(znum+0.001)
  deallocate(znum1d, zvel1d)      
 ! ==========
 ! total kinetic energy
 ! ========== 
tke=0.
aux=0.
Do elem=1,myDim_elem2D
   if(elem2D_nodes(1,elem)>myDim_nod2D) cycle
   DO nz=1, nlevels(elem)-1
      tke(1)=tke(1)+(UV(1,nz,elem)**2+UV(2,nz,elem)**2)*elem_area(elem)*(zbar(nz)-zbar(nz+1))
      tke(2)=tke(2)+elem_area(elem)*(zbar(nz)-zbar(nz+1))
   END DO
END DO
 ! =========
 ! Vertical velocity variance
 ! =========
 ww=0.0_8
Do m=1, myDim_nod2D
 ww=ww+Wvel(1,m)**2*(zbar(1)-Z(1))*area(1,m)
 Do nz=2, nlevels_nod2D(m)-1
 ww=ww+Wvel(nz,m)**2*area(nz,m)*(Z(nz-1)-Z(nz))
 end do
End do
 call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
 call MPI_AllREDUCE( tke, aux, 2, MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
 call MPI_Allreduce(ww,wwaux,1, MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
 if(mype==0) then
 open(5, file='tke.out', position='append')
 write(5,'(e13.5)') aux(1)/aux(2) 
 close(5)
 open(15, file='wvar.out', position='append')
 write(15,'(e13.5)') wwaux/aux(2)
 close(15)
 end if
 ! ===========
 ! Eddy kinetic energy (zonal mean is subtracted)
 ! ===========
tke=0.
aux=0.
Do elem=1,myDim_elem2D
   if(elem2D_nodes(1,elem)>myDim_nod2D) cycle
   DO nz=1, nlevels(elem)-1
      tke(1)=tke(1)+((UV(1,nz,elem)-zvel(nz,bpos(elem)))**2+UV(2,nz,elem)**2)*elem_area(elem)*(zbar(nz)-zbar(nz+1))
      tke(2)=tke(2)+elem_area(elem)*(zbar(nz)-zbar(nz+1))
   END DO
END DO

 call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
 call MPI_AllREDUCE( tke, aux, 2, MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
 if(mype==0) then
 open(5, file='eke.out', position='append')
 write(5,'(e13.5)') aux(1)/aux(2)
 close(5)
 end if
 end subroutine energy_out_soufflet
END MODULE Toy_Channel_Soufflet

