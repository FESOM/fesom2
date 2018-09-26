module diagnostics

  use g_config
  use o_mesh
  use g_parsup
  use g_clock
  use g_comm_auto
  use o_ARRAYS
  use g_forcing_arrays
  use i_ARRAYS
  use o_mixing_KPP_mod
  use g_rotate_grid
  use g_support
  implicit none

  private
  public :: compute_diagnostics, ldiag_solver, rhs_diag, lcurt_stress_surf, curl_stress_surf, ldiag_curl_vel3, curl_vel3, ldiag_energy, wzmid, wzmidrho, rho, &
            u_x_u, u_x_v, v_x_v, v_x_w, u_x_w, dudx, dudy, dvdx, dvdy, dudz, dvdz, utau_surf, utau_bott, av_dudz_sq, &
            taux_nod, tauy_nod, tbotx_nod, tboty_nod, av_nod, u_bott, v_bott

  ! Arrays used for diagnostics, some shall be accessible to the I/O
  ! 1. solver diagnostics: A*x=rhs? 
  ! A=ssh_stiff, x=d_eta, rhs=ssh_rhs; rhs_diag=A*x;
  real(kind=WP),  save, allocatable, target      :: rhs_diag(:)
  real(kind=WP),  save, allocatable, target      :: curl_stress_surf(:)
  real(kind=WP),  save, allocatable, target      :: curl_vel3(:,:)
  real(kind=WP),  save, allocatable, target      :: wzmid(:,:), wzmidrho(:,:), rho(:,:)
  real(kind=WP),  save, allocatable, target      :: u_x_u(:,:), u_x_v(:,:), v_x_v(:,:), v_x_w(:,:), u_x_w(:,:)
  real(kind=WP),  save, allocatable, target      :: dudx(:,:), dudy(:,:), dvdx(:,:), dvdy(:,:), dudz(:,:), dvdz(:,:)
  real(kind=WP),  save, allocatable, target      :: utau_surf(:), utau_bott(:), av_dudz_sq(:)
  real(kind=WP),  save, allocatable, target      :: taux_nod(:), tauy_nod(:), tbotx_nod(:), tboty_nod(:), av_nod(:,:), u_bott(:), v_bott(:)

  logical                                       :: ldiag_solver     =.false.
  logical                                       :: lcurt_stress_surf=.false.
  logical                                       :: ldiag_curl_vel3  =.false.
  logical                                       :: ldiag_energy     =.false.
  logical                                       :: ldiag_salt3D     =.true.
  contains

! ==============================================================
!rhs_diag=ssh_rhs?
subroutine diag_solver(mode)
  implicit none
  integer, intent(in)           :: mode
  integer                       :: n, is, ie
  logical, save                 :: firstcall=.true.

!=====================

  if (firstcall) then !allocate the stuff at the first call
     allocate(rhs_diag(mydim_nod2D))
     firstcall=.false.
     if (mode==0) return
  end if

  do n=1, myDim_nod2D
     is=ssh_stiff%rowptr_loc(n)
     ie=ssh_stiff%rowptr_loc(n+1)-1
     rhs_diag(n)=sum(ssh_stiff%values(is:ie)*d_eta(ssh_stiff%colind_loc(is:ie)))
  end do
end subroutine diag_solver
! ==============================================================
!curt(stress_surf)
subroutine diag_curl_stress_surf(mode)
  implicit none
  integer, intent(in)           :: mode
  logical, save                 :: firstcall=.true.
  integer         :: enodes(2), el(2), ed, n
  real(kind=WP)   :: deltaX1, deltaY1, deltaX2, deltaY2, c1

!=====================

  if (firstcall) then  !allocate the stuff at the first call
     allocate(curl_stress_surf(myDim_nod2D+eDim_nod2D))
     firstcall=.false.
     if (mode==0) return
  end if

  curl_stress_surf=0.

  DO ed=1, myDim_edge2D
     enodes=edges(:,ed)
     el=edge_tri(:,ed)
     deltaX1=edge_cross_dxdy(1,ed)
     deltaY1=edge_cross_dxdy(2,ed)
     if (el(2)>0) then
        deltaX2=edge_cross_dxdy(3,ed)
        deltaY2=edge_cross_dxdy(4,ed)
     end if
     c1=deltaX1*stress_surf(1,el(1))+deltaY1*stress_surf(2,el(1))
     curl_stress_surf(enodes(1))=curl_stress_surf(enodes(1))+c1
     curl_stress_surf(enodes(2))=curl_stress_surf(enodes(2))-c1
     if (el(2)>0) then
        c1= -deltaX2*stress_surf(1,el(2))-deltaY2*stress_surf(2,el(2))
        curl_stress_surf(enodes(1))=curl_stress_surf(enodes(1))+c1
        curl_stress_surf(enodes(2))=curl_stress_surf(enodes(2))-c1
     end if
  END DO
  DO n=1, myDim_nod2D+eDim_nod2D
     curl_stress_surf(n)=curl_stress_surf(n)/area(1,n)
  END DO
end subroutine diag_curl_stress_surf
! ==============================================================
!3D curl(velocity)
subroutine diag_curl_vel3(mode)
  implicit none
  integer, intent(in)           :: mode
  logical, save                 :: firstcall=.true.
  integer         :: enodes(2), el(2), ed, n, nz, nl1, nl2
  real(kind=WP)   :: deltaX1, deltaY1, deltaX2, deltaY2, c1

!=====================
  if (firstcall) then  !allocate the stuff at the first call
     allocate(curl_vel3(nl-1, myDim_nod2D+eDim_nod2D))
     firstcall=.false.
     if (mode==0) return
  end if

  curl_vel3=0.

  DO ed=1,myDim_edge2D
     enodes=edges(:,ed)
     el=edge_tri(:,ed)
     nl1=nlevels(el(1))-1
     deltaX1=edge_cross_dxdy(1,ed)
     deltaY1=edge_cross_dxdy(2,ed)
     nl2=0
     if (el(2)>0) then
        deltaX2=edge_cross_dxdy(3,ed)
        deltaY2=edge_cross_dxdy(4,ed)
        nl2=nlevels(el(2))-1
      end if     
      DO nz=1,min(nl1,nl2)
         c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))- &
         deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
         curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
         curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
      END DO
      DO nz=min(nl1,nl2)+1,nl1
         c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
         curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
         curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
      END DO
      DO nz=min(nl1,nl2)+1,nl2
         c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
         curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
         curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
      END DO
   END DO

   DO n=1, myDim_nod2D
      DO nz=1, nlevels_nod2D(n)-1
         curl_vel3(nz,n)=curl_vel3(nz,n)/area(nz,n)
      END DO
   END DO    
end subroutine diag_curl_vel3
! ==============================================================
!energy budget
subroutine diag_energy(mode)
  implicit none
  integer, intent(in)           :: mode
  logical, save                 :: firstcall=.true.
  integer            :: n, nz, k, i, elem, nzmax, elnodes(3)
  real(kind=WP)      :: ux, vx, uy, vy, tvol, rval(2)
  real(kind=WP)      :: stress_surf_n(2), stress_bott_n(2)
  real(kind=WP)      :: geo_grad_x(3), geo_grad_y(3), geo_u(3), geo_v(3)

!=====================
  if (firstcall) then  !allocate the stuff at the first call
     allocate(wzmid(nl-1, myDim_nod2D), wzmidrho(nl-1, myDim_nod2D), rho(nl-1, myDim_nod2D))
     allocate(u_x_u(nl-1, myDim_nod2D), u_x_v(nl-1, myDim_nod2D), v_x_v(nl-1, myDim_nod2D), v_x_w(nl-1, myDim_nod2D), u_x_w(nl-1, myDim_nod2D))
     allocate(dudx(nl-1, myDim_nod2D), dudy(nl-1, myDim_nod2D), dvdx(nl-1, myDim_nod2D), dvdy(nl-1, myDim_nod2D), dudz(nl-1, myDim_nod2D), dvdz(nl-1, myDim_nod2D))
     allocate(utau_surf(myDim_nod2D), utau_bott(myDim_nod2D), av_dudz_sq(myDim_nod2D))
     allocate(Av_nod(nl-1, myDim_nod2D), u_bott(myDim_nod2D), v_bott(myDim_nod2D))
     allocate(taux_nod(myDim_nod2D), tauy_nod(myDim_nod2D), tbotx_nod(myDim_nod2D), tboty_nod(myDim_nod2D))
     wzmid=0.
     rho  =0.
     wzmidrho=0.
     u_x_u=0.
     u_x_v=0.
     v_x_v=0.
     u_x_w=0.
     v_x_w=0.
     dudx=0.
     dudy=0.
     dvdx=0.
     dvdy=0.
     dudz=0.
     dvdz=0.
     av_dudz_sq=0.

     taux_nod=0.
     tauy_nod=0.
     tbotx_nod=0.
     tboty_nod=0.

     av_nod=0.
     u_bott=0.
     v_bott=0.


     firstcall=.false.
     if (mode==0) return
  end if
  
  !vertical velocity at the mid depths
  wzmid=0.5*(Wvel(1:nl-1, 1:myDim_nod2D)+Wvel(2:nl, 1:myDim_nod2D))
  rho  =density_m_rho0(1:nl-1, 1:myDim_nod2D)
  where(abs(rho)>0.)
     rho=rho+density_0
  end where

  !vertical velocity times density
  wzmidrho=rho(1:nl-1, 1:myDim_nod2D)*wzmid(1:nl-1, 1:myDim_nod2D)
  
  ! etc.
  u_x_u=Unode(1,1:nl-1,1:myDim_nod2D)*Unode(1,1:nl-1,1:myDim_nod2D)
  u_x_v=Unode(1,1:nl-1,1:myDim_nod2D)*Unode(2,1:nl-1,1:myDim_nod2D)
  v_x_v=Unode(2,1:nl-1,1:myDim_nod2D)*Unode(2,1:nl-1,1:myDim_nod2D)
  u_x_w=Unode(1,1:nl-1,1:myDim_nod2D)*wzmid(1:nl-1, 1:myDim_nod2D)
  v_x_w=Unode(2,1:nl-1,1:myDim_nod2D)*wzmid(1:nl-1, 1:myDim_nod2D)

  utau_surf=0.
  utau_bott=0.
  taux_nod =0.
  tauy_nod =0.
  tbotx_nod =0.
  tboty_nod =0.
  u_bott   =0.
  v_bott   =0.
  Av_nod   =0.
  ! this loop might be very expensive
  DO n=1, myDim_nod2D
     nzmax=nlevels_nod2D(n)
     stress_surf_n=0.
     stress_bott_n=0.
     DO nz=1, nzmax-1
        tvol=0.0_WP
        ux  =0.0_WP
        uy  =0.0_WP
        vx  =0.0_WP
        vy  =0.0_WP
        DO k=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(k,n)
	   elnodes=elem2D_nodes(:, elem)
           if (nz==1) then
              rval=stress_surf(:,elem)/density_0*elem_area(elem)
              stress_surf_n = stress_surf_n+rval/Av(2, elem)      !accumulation of: du/dz at the surface (du/dz=tau/Av)
              utau_surf(n)=utau_surf(n)+sum(rval(:)*UV(:,1,elem)) !a scalar product tau_surface times u 
                                                                  !shall rotation be considered? (scalar product)
              taux_nod(n)=taux_nod(n)+rval(1)
              tauy_nod(n)=tauy_nod(n)+rval(2)
           end if
           if (nlevels(elem)-1 == nz) then
              rval=C_d*sqrt(UV(1, nz,elem)**2+ UV(2, nz, elem)**2)*UV(:,nz,elem)*elem_area(elem)
              stress_bott_n = stress_bott_n-rval/Av(nz, elem)  !accumulation of: du/dz at the bottom (du/dz=tau/Av)
              utau_bott(n)=utau_bott(n)-sum(rval*UV(:,nz,elem))!a scalar product tau_bottom times u 
                                                               !shall rotation be considered? (scalar product)
              u_bott(n)=u_bott(n)+UV(1,nz,elem)*elem_area(elem)
              v_bott(n)=v_bott(n)+UV(2,nz,elem)*elem_area(elem)
               
              tbotx_nod(n)=tbotx_nod(n)-rval(1)
              tboty_nod(n)=tboty_nod(n)-rval(2)
           end if
           if (nlevels(elem)-1 < nz) cycle
           Av_nod(nz, n)=Av_nod(nz, n)+Av(nz+1, elem)*elem_area(elem)!we will need Av at nodes later
           geo_grad_x=gradient_sca(1:3,elem)                       !reassign vector arrays to make the rotation possible (if required)
           geo_grad_y=gradient_sca(4:6,elem)
           geo_u=Unode(1,nz,elnodes)
           geo_v=Unode(2,nz,elnodes)
           if (rotated_grid) then
              DO i=1, 3
                 call vector_r2g(geo_grad_x(i), geo_grad_y(i), coord_nod2D(1, n), coord_nod2D(2, n), 0)
                 call vector_r2g(geo_u(i), geo_v(i), coord_nod2D(1, n), coord_nod2D(2, n), 0)
              END DO
           end if
           tvol=tvol+elem_area(elem)
           ux=ux+sum(geo_grad_x*geo_u)*elem_area(elem)         !accumulate tensor of velocity derivatives
           vx=vx+sum(geo_grad_x*geo_v)*elem_area(elem)
           uy=uy+sum(geo_grad_y*geo_u)*elem_area(elem)
           vy=vy+sum(geo_grad_y*geo_v)*elem_area(elem)
        END DO
        Av_nod(nz, n)=Av_nod(nz, n)/tvol !convert the accumulation over layers to actual nodal values
        dudx(nz,n)=ux/tvol!/area(nz, n)/3.
        dvdx(nz,n)=vx/tvol
        dudy(nz,n)=uy/tvol
        dvdy(nz,n)=vy/tvol
     END DO

!    Av_nod(1, n)   =Av_nod(2, n)              !convert the accumulation over surface to actual surface values at nodes
     stress_surf_n=stress_surf_n/area(1, n)/3. !*3. because of full areas were used in the loop
     stress_bott_n=stress_bott_n/area(1, n)/3.
     utau_surf(n)=utau_surf(n)/area(1, n)/3.
     utau_bott(n)=utau_bott(n)/area(1, n)/3.
     
     taux_nod(n)=taux_nod(n)/area(1, n)/3.
     tauy_nod(n)=tauy_nod(n)/area(1, n)/3.

     tbotx_nod(n)=tbotx_nod(n)/area(1, n)/3.
     tboty_nod(n)=tboty_nod(n)/area(1, n)/3.

     u_bott(n)=u_bott(n)/area(1, n)/3.
     v_bott(n)=v_bott(n)/area(1, n)/3.

     if (rotated_grid) then
        call vector_r2g(stress_surf_n(1), stress_surf_n(2), coord_nod2D(1, n), coord_nod2D(2, n), 0)
        call vector_r2g(stress_bott_n(1), stress_bott_n(2), coord_nod2D(1, n), coord_nod2D(2, n), 0)
        call vector_r2g(taux_nod(n),  tauy_nod(n),  coord_nod2D(1, n), coord_nod2D(2, n), 0)
        call vector_r2g(tbotx_nod(n), tboty_nod(n), coord_nod2D(1, n), coord_nod2D(2, n), 0)
        call vector_r2g(u_bott(n), v_bott(n), coord_nod2D(1, n), coord_nod2D(2, n), 0)
     end if
         
     ! compute Z_n which is the mid depth of prisms (ALE supports changing layer thicknesses)
     zbar_n=0.0_WP
     Z_n=0.0_WP
     zbar_n(nzmax)=zbar_n_bot(n)
     Z_n(nzmax-1)=zbar_n(nzmax) + hnode_new(nzmax-1,n)/2.0_WP
     do nz=nzmax-1,2,-1
        zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
        Z_n(nz-1) = zbar_n(nz) + hnode_new(nz-1,n)/2.0_WP
     end do
     zbar_n(1) = zbar_n(2) + hnode_new(1,n)

     !compute du/dz & dv/dz
     dudz(2:nzmax-2, n)=(Unode(1, 1:nzmax-3, n)-Unode(1, 3:nzmax-1, n))/(Z_n(1:nzmax-3)-Z_n(3:nzmax-1)) !central differences in vertical
     dvdz(2:nzmax-2, n)=(Unode(2, 1:nzmax-3, n)-Unode(2, 3:nzmax-1, n))/(Z_n(1:nzmax-3)-Z_n(3:nzmax-1))

     dudz(1, n)=stress_surf_n(1)
     dvdz(1, n)=stress_surf_n(2)
     dudz(nl-1, n)    =stress_bott_n(1)
     dvdz(nl-1, n)    =stress_bott_n(2)

     !compute int(Av * (du/dz)^2)
     av_dudz_sq(n)=sum((dudz(1:nzmax-1, n)**2+dvdz(1:nzmax-1, n)**2)*Av_nod(1:nzmax-1, n)* hnode_new(1:nzmax-1,n))
  END DO  
end subroutine diag_energy
! ==============================================================
subroutine compute_diagnostics(mode)
  implicit none
  integer, intent(in)           :: mode !constructor mode (0=only allocation; any other=do diagnostic)
  real(kind=WP)                 :: val  

  !1. solver diagnostic
  if (ldiag_solver)      call diag_solver(mode)
  !2. compute curl(stress_surf)
  if (lcurt_stress_surf) call diag_curl_stress_surf(mode)
  !3. compute curl(velocity)
  if (ldiag_curl_vel3)   call diag_curl_vel3(mode)
  !4. compute energy budget
  if (ldiag_energy)      call diag_energy(mode)
  !5. print integrated temperature 
  if (ldiag_salt3d) then
     if (mod(mstep,logfile_outfreq)==0) then
        call integrate_nod(tr_arr(:,:,2), val)
        if (mype==0) then
           write(*,*) 'total integral of salinity at timestep :', mstep, val
        end if
     end if
  end if
end subroutine compute_diagnostics
end module diagnostics
