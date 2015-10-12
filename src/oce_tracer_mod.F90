!============================================================================================
module o_tracers
implicit none
contains
!=======================================================================
subroutine pp_diff(nz, node, diffcoeff)
! Computes the Richardson number dependent diffusivity as in Gent & Cane 
! parameterization of mixing (see Large and Gent JPO 1999, 449-464).
! (slightly modified Pacanowski-Philander)
!
USE o_PARAM
USE o_MESH
USE o_ARRAYS
IMPLICIT NONE
integer, intent(in)         :: nz, node
real*8, intent(out)   :: diffcoeff

real*8                :: dz_inv, bv, shear, a, rho_up, rho_dn, t, s
integer                     :: nn, elem, ne 
 
  dz_inv=1./(Z(nz-1)-Z(nz))
  t=tr_arr(nz-1, node,1)
  s=tr_arr(nz-1, node,2)
  call densityJM(t, s, zbar(nz), rho_up)
  t=tr_arr(nz, node,1)
  s=tr_arr(nz, node,2)
  call densityJM(t, s, zbar(nz), rho_dn)  
  ! both densities are referenced to the plane 
  ! where flux is computed (zbar(nz)
  bv  = -g*dz_inv*(rho_up-rho_dn)/density_0
  
    if (bv<0) then
        a=1_WP   
    else
        shear=0.0_WP
	nn=0
	DO ne=1, nod_in_elem2D_num(node)
	   elem=nod_in_elem2D(ne,node)
	   if(nz>nlevels(elem)-1) cycle
	   shear=shear+(UV(1,nz-1,elem)-UV(1,nz,elem))**2 +&
	               (UV(2,nz-1,elem)-UV(2,nz,elem))**2 
           nn=nn+1
	END DO   
        shear=shear*dz_inv*dz_inv/nn
	a=shear/(shear+10*bv)
    end if
    
    diffcoeff=K_ver+0.01_WP*a*a*a	    	
    
    
end subroutine pp_diff 
!=======================================================================
SUBROUTINE tracer_gradient_nodes (tt_xy, tt_xynodes)

! tt_xy      elemental gradient of tracer 
! tt_xynodes nodal gradients obtained by averaging
USE o_PARAM
USE o_MESH
!USE o_ARRAYS
USE g_PARSUP
use g_comm_auto
IMPLICIT NONE
real*8      :: tt_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
real*8      :: tt_xynodes(2,nl-1,myDim_nod2D+eDim_nod2D)
integer         :: n, nz, elem, k
real*8   :: tvol, tx, ty
  DO n=1, myDim_nod2D !! m=1, myDim_nod2D
                      !! n=myList_nod2D(m)
     DO nz=1, nlevels_nod2D(n)-1
        tvol=0.0_WP
	tx=0.0_WP
	ty=0.0_WP
	DO k=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(k,n)
           if(nlevels(elem)-1<nz) cycle
           tvol=tvol+elem_area(elem)
           tx=tx+tt_xy(1,nz,elem)*elem_area(elem)
	   ty=ty+tt_xy(2,nz,elem)*elem_area(elem)
        END DO
	tt_xynodes(1,nz,n)=tx/tvol
	tt_xynodes(2,nz,n)=ty/tvol
     END DO
   END DO 
   call exchange_nod(tt_xynodes(1,:,:))
   call exchange_nod(tt_xynodes(2,:,:))
END SUBROUTINE tracer_gradient_nodes
!=======================================================================
SUBROUTINE tracer_gradient_elements (ttf, tt_xy)
! ttf - tracer field
! ttx, tty  elemental gradient of tracer 

USE o_PARAM
USE o_MESH
!USE o_ARRAYS
USE g_PARSUP
IMPLICIT NONE
real*8      :: ttf(nl-1,myDim_nod2D+eDim_nod2D)
real*8      :: tt_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
integer           :: elem,  elnodes(3)
integer           :: n, nz

 DO elem=1, myDim_elem2D       !! m=1, myDim_elem2D
                               !! elem=myList_elem2D(m)
   elnodes=elem2D_nodes(:,elem)
   DO nz=1, nlevels(elem)-1   
      tt_xy(1,nz, elem)=sum(gradient_sca(1:3,elem)*ttf(nz,elnodes))
      tt_xy(2,nz, elem)=sum(gradient_sca(4:6,elem)*ttf(nz,elnodes))
   END DO
 END DO
     
END SUBROUTINE tracer_gradient_elements
!========================================================================================
SUBROUTINE init_tracers(tr_num)
use g_parsup
use o_PARAM, only: tracer_adv,Redi_GM
use o_arrays
use o_mesh
IMPLICIT NONE
integer :: tr_num,n,nz
!Allocate stuff if needed
!Fill work arrays
! elemental gradients are in ttx, tty
! nodal gradients are in ttxnodes, ttynodes.
!TS grads are stored separately from the rest of tracers
!they are loaded in work arrays ttx,tty,ttxnodes,ttynodes,tsv
del_ttf=0d0

if (tr_num<=2) then
  if (tr_num==1) then
  do n=1,2
   call tracer_gradient_elements(tr_arr(:,:,n), tt_xy_stored(:,:,:,n))
   call tracer_gradient_nodes(tt_xy_stored(:,:,:,n), tt_xynodes_stored(:,:,:,n))
  enddo

  if (Redi_GM) then
      call compute_neutral_slope(tr_arr(:,:,1),tr_arr(:,:,2))
      call redi_gm_coef  
  endif
   
  endif
else
  call tracer_gradient_elements(tr_arr(:,:,tr_num), tt_xy_stored(:,:,:,tr_num))
  call tracer_gradient_nodes(tt_xy_stored(:,:,:,tr_num), tt_xynodes_stored(:,:,:,tr_num))
endif
select case (tracer_adv)
case(1) 
	continue
case(2) !Quadratic reconstr.
! Get coefs
  	call linear_reconstruction(tt_xy_stored(:,:,:,tr_num))
        call quadratic_reconstruction(tr_arr(:,:,tr_num)) 
CASE DEFAULT !unknown
	if (mype==0) write(*,*) 'Unknown advection type. Check your namelists.'
        call par_ex(1)
END SELECT
end subroutine init_tracers
!=================================================================================
SUBROUTINE init_tracers_AB(tr_num)
use g_parsup
use o_PARAM, only: tracer_adv,Redi_GM
use o_arrays
use o_mesh
use g_comm
IMPLICIT NONE
integer :: tr_num,n,nz 
!Allocate stuff if needed
!Fill work arrays
! elemental gradients are in ttx, tty
! nodal gradients are in ttxnodes, ttynodes.
!TS grads are stored separately from the rest of tracers
!they are loaded in work arrays ttx,tty,ttxnodes,ttynodes,tsv

del_ttf=0d0


if (tr_num<=2) then
  if (tr_num==1) then
  tr_arr_old(:,:,1:2)=-(0.5+epsilon)*tr_arr_old(:,:,1:2)+(1.5+epsilon)*tr_arr(:,:,1:2)

 do n=1,2
  call tracer_gradient_elements(tr_arr_old(:,:,n), tt_xy_stored(:,:,:,n))
  call tracer_gradient_nodes(tt_xy_stored(:,:,:,n), tt_xynodes_stored(:,:,:,n))
 enddo
  if (Redi_GM) then
      call compute_neutral_slope(tr_arr_old(:,:,1),tr_arr_old(:,:,2))
      call redi_gm_coef 
  endif
  endif
else
! =================
! AB interpolation
! =================
  tr_arr_old(:,:,tr_num)=-(0.5+epsilon)*tr_arr_old(:,:,tr_num)+(1.5+epsilon)*tr_arr(:,:,tr_num)
  call tracer_gradient_elements(tr_arr_old(:,:,tr_num), tt_xy_stored(:,:,:,tr_num))
  call tracer_gradient_nodes(tt_xy_stored(:,:,:,tr_num), tt_xynodes_stored(:,:,:,tr_num))
endif

select case (tracer_adv)
case(3,4,5)

        call exchange_elem3D_full(tt_xy_stored(1,:,:,tr_num))
	call exchange_elem3D_full(tt_xy_stored(2,:,:,tr_num))
	call fill_up_dn_grad(tt_xy_stored(:,:,:,tr_num))
CASE DEFAULT !Unknown
	if (mype==0) write(*,*) 'Unknown advection type. Check your namelists.'
        call par_ex(1)
END SELECT
end subroutine init_tracers_AB

!==============================================================================
SUBROUTINE adv_tracers(tr_num)
use g_parsup
use o_mesh,only: nl,nod2D
use o_PARAM, only: tracer_adv
use o_arrays
IMPLICIT NONE
integer :: tr_num 
select case (tracer_adv)
case(1) !Miura
	call adv_tracer_miura(tr_arr(:,:,tr_num), del_ttf, tr_num,tt_xy_stored(:,:,:,tr_num),tt_xynodes_stored(:,:,:,tr_num) )
case(2) !Quadratic reconstr.
	call adv_tracer_second_order_rec(tr_arr(:,:,tr_num),del_ttf,tr_num,tt_xy_stored(:,:,:,tr_num),tt_xynodes_stored(:,:,:,tr_num) )
case(3) !MUSCL
        call adv_tracer_muscl(tr_arr(:,:,tr_num), del_ttf, tr_arr_old(:,:,tr_num), tr_num,tt_xy_stored(:,:,:,tr_num),tt_xynodes_stored(:,:,:,tr_num) )
case(4) !MUSCL+FCT(3D)
        call adv_tracer_fct(tr_arr(:,:,tr_num),del_ttf,tr_arr_old(:,:,tr_num), 0.0_WP)
case(5) !MUSCL+FCT(2D+1D)
        call adv_tracer_fct2p1_34(tr_arr(:,:,tr_num),del_ttf,tr_arr_old(:,:,tr_num),tr_num,tt_xy_stored(:,:,:,tr_num),tt_xynodes_stored(:,:,:,tr_num) )
CASE DEFAULT !unknown
	if (mype==0) write(*,*) 'Unknown advection type. Check your namelists.'
        call par_ex(1)
END SELECT
end subroutine adv_tracers
!==================================
SUBROUTINE diff_tracers(tr_num)
use o_mesh
USE g_PARSUP
use o_arrays
IMPLICIT NONE
integer, intent(in) :: tr_num
integer             :: n, nl1

call diff_part_hor(tr_arr(:,:,tr_num),del_ttf,tr_num,tt_xy_stored(:,:,:,tr_num),tt_xynodes_stored(:,:,:,tr_num))
if (Redi_GM)          call ver_redi_gm(tr_arr(:,:,tr_num),del_ttf,tr_num, tt_xynodes_stored(:,:,:,tr_num))
if (not(i_vert_diff)) call diff_ver_part_expl(tr_arr(:,:,tr_num),del_ttf,tr_num)

!Update tracer berofe implicit operator is applied
DO n=1, myDim_nod2D !! m=1, myDim_nod2D
      nl1=nlevels_nod2D(n)-1
      tr_arr(1:nl1,n,tr_num)=tr_arr(1:nl1,n,tr_num)+del_ttf(1:nl1,n)
END DO
!We do not set del_ttf to zero because it will not be used in this timestep anymore
!init_tracers will set it to zero for the next timestep

if (i_vert_diff)      call diff_ver_part_impl(tr_arr(:,:,tr_num), tr_num)
end subroutine diff_tracers
!===========================================================================
subroutine diff_part_hor(ttf, dttf, tr_num,tt_xy,tt_xynodes)
use o_ARRAYS
use g_PARSUP
use o_MESH
USE o_param
use g_config
implicit none
real*8::deltaX1,deltaY1,deltaX2,deltaY2
integer::edge,tr_num
integer::nl1,nz,el(2),elnodes(3),n,nl2,enodes(2)
real*8::temp(3),d,c1,Tx,Ty,Kh
real*8 :: ttf(nl-1, myDim_nod2D+eDim_nod2D), dttf(nl-1, myDim_nod2D+eDim_nod2D)
real*8 :: tsv(nl,myDim_nod2D+eDIm_nod2D),tvol
real*8      :: tt_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
real*8      :: tt_xynodes(2,nl-1,myDim_nod2D+eDim_nod2D)
ttrhs=0d0
DO n=1,myDim_nod2D+eDim_nod2D
     DO nz=2,nlevels_nod2D(n)-1
        tvol=Z(nz-1)-Z(nz)
        tsv(nz,n)=(ttf(nz-1,n)-ttf(nz,n))/tvol
     END DO
        tsv(1,n)=0.0_WP
	nz=nlevels_nod2D(n)
	tsv(nz,n)=0.0_WP
END DO
DO edge=1, myDim_edge2D
 deltaX1=edge_cross_dxdy(1,edge)
 deltaY1=edge_cross_dxdy(2,edge)
 el=edge_tri(:,edge)
 enodes=edges(:,edge)
 nl1=nlevels(el(1))-1
 elnodes=elem2d_nodes(:,el(1))
 Kh=elem_area(el(1))
 if (el(2)>0) then 
   Kh=0.5*(Kh+elem_area(el(2)))
   nl2=nlevels(el(2))-1
 endif
 Kh=K_hor*Kh/scale_area

 do nz=1,nl1
  Tx=0d0
  Ty=0d0

  if (Redi_GM) then
    d=sum(Kd(1,nz,elnodes))/3.0_WP
    temp=(Kd(2,nz,elnodes)-Kd(3,nz,elnodes))*neutral_slope(1,nz,elnodes)
    c1=sum(temp*(tsv(nz,elnodes)+tsv(nz+1,elnodes)))/6.0_WP
    Tx=d*tt_xy(1,nz,el(1))+c1
    temp=(Kd(2,nz,elnodes)-Kd(3,nz,elnodes))*neutral_slope(2,nz,elnodes)
    c1=sum(temp*(tsv(nz,elnodes)+tsv(nz+1,elnodes)))/6.0_WP
    Ty=d*tt_xy(2,nz,el(1))+c1
  endif !NO Redi_GM

  if ((el(2)>0).and.(nz<=nl2)) then
    Tx=Tx+0.5*Kh*(tt_xy(1,nz,el(1))+tt_xy(1,nz,el(2)))
    Ty=Ty+0.5*Kh*(tt_xy(2,nz,el(1))+tt_xy(2,nz,el(2)))
  else
    Tx=Tx+Kh*tt_xy(1,nz,el(1))
    Ty=Ty+Kh*tt_xy(2,nz,el(1))
  end if
  
  ttrhs(nz,enodes(1)) = ttrhs(nz,enodes(1)) - Ty*deltaX1 + Tx*deltaY1
  ttrhs(nz,enodes(2)) = ttrhs(nz,enodes(2)) + Ty*deltaX1 - Tx*deltaY1
 enddo

 if (el(2)>0) then
  deltaX2=edge_cross_dxdy(3,edge)
  deltaY2=edge_cross_dxdy(4,edge)
  do nz=1,nl2
  Tx=0d0
  Ty=0d0
  if (Redi_GM) then
   elnodes=elem2D_nodes(:,el(2))
   d=sum(Kd(1,nz,elnodes))/3.0_WP                                                                     
   temp=(Kd(2,nz,elnodes)-Kd(3,nz,elnodes))*neutral_slope(1,nz,elnodes)                               
   c1=sum(temp*(tsv(nz,elnodes)+tsv(nz+1,elnodes)))/6.0_WP
   Tx=d*tt_xy(1,nz,el(2))+c1
   temp=(Kd(2,nz,elnodes)-Kd(3,nz,elnodes))*neutral_slope(2,nz,elnodes)                               
   c1=sum(temp*(tsv(nz,elnodes)+tsv(nz+1,elnodes)))/6.0_WP                                    
   Ty=d*tt_xy(2,nz,el(2))+c1
  endif

  if (nz<=nl1) then
      Tx=Tx+0.5*Kh*(tt_xy(1,nz,el(1))+tt_xy(1,nz,el(2)))
      Ty=Ty+0.5*Kh*(tt_xy(2,nz,el(1))+tt_xy(2,nz,el(2)))
  else
      Tx=Tx+Kh*tt_xy(1,nz,el(2))
      Ty=Ty+Kh*tt_xy(2,nz,el(2))
  end if

   ttrhs(nz,enodes(1)) = ttrhs(nz,enodes(1)) + Ty*deltaX2 - Tx*deltaY2
   ttrhs(nz,enodes(2)) = ttrhs(nz,enodes(2)) - Ty*deltaX2 + Tx*deltaY2
  enddo
 endif

ENDDO

DO n=1, myDim_nod2D
     DO nz=1,nlevels_nod2D(n)-1
        dttf(nz,n)=dttf(nz,n)+ttrhs(nz,n)*dt/area(nz,n)
     END DO
END DO
end subroutine diff_part_hor
!========================================================
subroutine ver_redi_gm(ttf, dttf, tr_num,tt_xynodes)
! Vertical diffusive flux from redi_GM:                                                                            
use o_ARRAYS,only: neutral_slope,Kd,Kv
use o_MESH
use g_PARSUP
use g_config,only: dt
implicit none
real*8::vd_flux(nl-1)
real*8 :: Tx
integer:: nz,nl1,n,tr_num
real*8 :: ttf(nl-1, myDim_nod2D+eDim_nod2D), dttf(nl-1, myDim_nod2D+eDim_nod2D)
real*8      :: tt_xynodes(2,nl-1,myDim_nod2D+eDim_nod2D)
DO n=1, myDim_nod2D
nl1=nlevels_nod2D(n)-1
vd_flux=0d0

do nz=2,nl1
  vd_flux(nz)=  ((Kd(2,nz-1,n)+Kd(3,nz-1,n))*(Z(nz-1)-zbar(nz)) \
      *( neutral_slope(1,nz-1,n)*tt_xynodes(1,nz-1,n)+neutral_slope(2,nz-1,n)*tt_xynodes(2,nz-1,n) ) \
      + \
       (Kd(2,nz,n)+Kd(3,nz,n))*(zbar(nz)-Z(nz)) \
      *( neutral_slope(1,nz,n)*tt_xynodes(1,nz,n)+neutral_slope(2,nz,n)*tt_xynodes(2,nz,n) )) \
      / (Z(nz-1)-Z(nz))*area(nz,n)
enddo
do nz=1,nl1-1
 dttf(nz,n) = dttf(nz,n)+(vd_flux(nz) - vd_flux(nz+1))/(zbar(nz)-zbar(nz+1))*dt/area(nz,n)
enddo
dttf(nl1,n) = dttf(nl1,n)+vd_flux(nl1)/(zbar(nl1)-zbar(nl1+1))*dt/area(nl1,n)

ENDDO
end subroutine ver_redi_gm
!===================================================================================
subroutine diff_ver_part_expl(ttf, dttf, tr_num)
! Vertical diffusive flux(explicit scheme):                                                                            
use o_ARRAYS
use o_MESH
use g_PARSUP
use g_config,only: dt
implicit none
real*8::vd_flux(nl-1)
real*8 :: rdata,flux,rlx
integer:: nz,nl1,tr_num,n
real*8 :: ttf(nl-1, myDim_nod2D+eDim_nod2D), dttf(nl-1, myDim_nod2D+eDim_nod2D)
real*8 :: zinv1,Ty
DO n=1, myDim_nod2D
nl1=nlevels_nod2D(n)-1
vd_flux=0d0

if (tr_num==1) then
  flux  = -heat_flux(n)/density_0/4200.0
  rdata =  Tsurf(n)
  rlx   =  surf_relax_T
elseif (tr_num==2) then
  flux  =  water_flux(n)*tr_arr(1,n,2)
  rdata =  Ssurf(n)
  rlx   =  surf_relax_S
else
  flux  = 0d0
  rdata = 0d0
  rlx=0d0
endif
!Surface forcing

vd_flux(1)= flux + &
                    rlx*(rdata-ttf(1,n))

do nz=2,nl1
 zinv1=1.0_WP/(Z(nz-1)-Z(nz))

 Ty= Kd(4,nz-1,n)*(Z(nz-1)-zbar(nz))*zinv1 *neutral_slope(3,nz-1,n)**2 + \
              Kd(4,nz,n)*(zbar(nz)-Z(nz))*zinv1 *neutral_slope(3,nz,n)**2

 vd_flux(nz) = (Kv(nz,n)+Ty)*(ttf(nz-1,n)-ttf(nz,n))*zinv1*area(nz,n)
enddo
do nz=1,nl1-1
 dttf(nz,n) =dttf(nz,n) + (vd_flux(nz) - vd_flux(nz+1))/(zbar(nz)-zbar(nz+1))*dt/area(nz,n)
enddo
dttf(nl1,n) = dttf(nl1,n) + (vd_flux(nl1)/(zbar(nl1)-zbar(nl1+1)))*dt/area(nl1,n)
ENDDO

end subroutine diff_ver_part_expl
!===========================================================================
subroutine diff_ver_part_impl(ttf, tr_num)
USE o_MESH
USE o_PARAM
USE o_ARRAYS
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
real*8, intent(inout) :: ttf(nl-1, myDim_nod2D+eDim_nod2D)
integer, intent(in)   :: tr_num
real*8                :: a(nl), b(nl), c(nl), tr(nl)
real*8                ::  cp(nl), tp(nl)
integer               ::  nz, n, nzmax
real*8                ::  m, zinv, dt_inv
real*8                ::  rsss, Ty,Ty1,c1,zinv1,zinv2, v_adv
   dt_inv=1.0_WP/dt
   DO n=1,myDim_nod2D      !! m=1,myDim_nod2D 
                           !! n=myList_nod2D(m)
          nzmax=nlevels_nod2D(n)
          ! The first row
          nz=1
          zinv2=1.0_WP/(Z(nz)-Z(nz+1))
          zinv=1.0_WP*dt/(zbar(1)-zbar(2))
          Ty1= Kd(4,nz,n)*(Z(nz)-zbar(nz+1))*zinv2 *neutral_slope(3,nz,n)**2 + \
              Kd(4,nz+1,n)*(zbar(nz+1)-Z(nz+1))*zinv2 *neutral_slope(3,nz+1,n)**2
          c(1)=-(Kv(2,n)+Ty1)*zinv2*zinv*area(2,n)/area(1,n)
          a(1)=0.0_WP
          b(1)=-c(1)+1.0_WP
	  ! update from the vertical advection
          v_adv=zinv*area(2,n)/area(1,n)
          b(1)=b(1)+Wvel_i(1, n)*zinv-min(0._WP, Wvel_i(2, n))*v_adv
          c(1)=c(1)-max(0._WP, Wvel_i(2, n))*v_adv

          zinv1=zinv2
          ! Regular part of coefficients:
          DO nz=2, nzmax-2
             zinv2=1.0_WP/(Z(nz)-Z(nz+1))
             Ty= Kd(4,nz-1,n)*(Z(nz-1)-zbar(nz))*zinv1 *neutral_slope(3,nz-1,n)**2 + \
              Kd(4,nz,n)*(zbar(nz)-Z(nz))*zinv1 *neutral_slope(3,nz,n)**2
             Ty1= Kd(4,nz,n)*(Z(nz)-zbar(nz+1))*zinv2 *neutral_slope(3,nz,n)**2 + \
              Kd(4,nz+1,n)*(zbar(nz+1)-Z(nz+1))*zinv2 *neutral_slope(3,nz+1,n)**2
        
             zinv=1.0_WP*dt/(zbar(nz)-zbar(nz+1))
             a(nz)=-(Kv(nz,n)+Ty)*zinv1*zinv
             c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv*area(nz+1,n)/area(nz,n)
             b(nz)=-a(nz)-c(nz)+1.0_WP
             zinv1=zinv2

  	     ! update from the vertical advection
             v_adv=zinv
	     a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv
             b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv

  	     v_adv=v_adv*area(nz+1,n)/area(nz,n)
             b(nz)=b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
             c(nz)=c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
          END DO
          ! The last row
          nz=nzmax-1
          zinv=1.0_WP*dt/(zbar(nzmax-1)-zbar(nzmax))
          Ty= Kd(4,nz-1,n)*(Z(nz-1)-zbar(nz))*zinv1 *neutral_slope(3,nz-1,n)**2 + \
              Kd(4,nz,n)*(zbar(nz)-Z(nz))*zinv1 *neutral_slope(3,nz,n)**2
          a(nzmax-1)=-(Kv(nzmax-1,n)+Ty)*zinv1*zinv
          b(nzmax-1)=-a(nzmax-1)+1.0_WP
          c(nzmax-1)=0.0_WP

	  ! update from the vertical advection
          v_adv=zinv
          a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv       
          b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv

          ! ===========================================
          ! The rhs:
          DO nz=2,nzmax-2
                tr(nz)=-a(nz)*ttf(nz-1,n)-(b(nz)-1.0_WP)*ttf(nz,n)-c(nz)*ttf(nz+1,n)
          END DO
          nz=nzmax-1
          tr(nz)=-a(nz)*ttf(nz-1,n)-(b(nz)-1.0_WP)*ttf(nz,n)
          !  The first row contains also surface forcing
          tr(1)=-(b(1)-1.0_WP)*ttf(1,n)-c(1)*ttf(2,n)

          zinv=1.0_WP*dt/(zbar(1)-zbar(2))
          if (tr_num==1) then
          tr(1)=tr(1)  -  &
                 zinv*(1d-3*heat_flux(n)/(4.2_WP)/density_0 - surf_relax_T*(Tsurf(n)-tr_arr(1,n,1)))
          elseif (tr_num==2) then
          rsss=ref_sss
          if(ref_sss_local) rsss = tr_arr(1,n,2)
          tr(1)=tr(1)  +  &
                 zinv*(rsss*water_flux(n) + surf_relax_S*(Ssurf(n)-tr_arr(1,n,2)))
          endif
          ! =============================================
          ! The sweep algorithm
          ! initialize c-prime and s,t-prime
          cp(1) = c(1)/b(1)
          tp(1) = tr(1)/b(1)
! solve for vectors c-prime and t, s-prime
          DO nz = 2,nzmax-1
           m = b(nz)-cp(nz-1)*a(nz)
           cp(nz) = c(nz)/m
           tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
          END DO
! initialize x
          tr(nzmax-1) = tp(nzmax-1)
         ! solve for x from the vectors c-prime and d-prime
          do nz = nzmax-2, 1, -1
           tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
          end do
        
          DO nz=1,nzmax-1
            ttf(nz,n)=ttf(nz,n)+tr(nz)
          END DO
        
   END DO   !!! cycle over nodes
   END subroutine diff_ver_part_impl
!====================================================================
Subroutine relax_to_clim(tr_num)
use o_mesh
use g_config,only: dt
USE g_PARSUP
use o_arrays
use o_PARAM, only: tracer_adv
implicit none
integer :: tr_num,n,nz

if ((tracer_adv==3).or.(tracer_adv==4).or.(tracer_adv==5).or.(Redi_GM)) then
     tr_arr_old(:,:,tr_num)=tr_arr(:,:,tr_num)
endif

if((clim_relax>1.0e-8).and.(tr_num==1)) then
  DO n=1, myDim_nod2D !! m=1, myDim_nod2D
                      !! n=myList_nod2D(m)
tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)=tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)+relax2clim(n)*dt*(Tclim(1:nlevels_nod2D(n)-1,n)-tr_arr(1:nlevels_nod2D(n)-1,n,tr_num))
  END DO
end if
if((clim_relax>1.0e-8).and.(tr_num==2)) then
  DO n=1, myDim_nod2D
tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)=tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)+relax2clim(n)*dt*(Sclim(1:nlevels_nod2D(n)-1,n)-tr_arr(1:nlevels_nod2D(n)-1,n,tr_num))
  END DO
end if 
end subroutine relax_to_clim
end module o_tracers
