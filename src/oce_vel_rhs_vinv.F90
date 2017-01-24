! Vector invariant momentum advection:
! (curl u+f)\times u+grad(u^2/2)+w du/dz
!
! ===================================================================
subroutine relative_vorticity
USE o_ARRAYS
USE o_MESH
USE g_PARSUP
  use g_comm_auto
IMPLICIT NONE
integer        :: n, nz, el(2), enodes(2), nl1, nl2, edge
real(kind=WP)  :: deltaX1, deltaY1, deltaX2, deltaY2, c1

DO n=1,myDim_nod2D
                                 !! n=myList_nod2D(m)
   DO nz=1, nlevels_nod2D(n)-1
      vorticity(nz,n)=0.0_WP
   END DO
END DO      

DO edge=1,myDim_edge2D
                                 !! edge=myList_edge2D(m)
   enodes=edges(:,edge)
   el=edge_tri(:,edge)
   nl1=nlevels(el(1))-1
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   nl2=0
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   nl2=nlevels(el(2))-1
   end if     
   DO nz=1,min(nl1,nl2)
      c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))- &
      deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
      vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
      vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
   END DO
   DO nz=min(nl1,nl2)+1,nl1
      c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
      vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
      vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
   END DO
   DO nz=min(nl1,nl2)+1,nl2
      c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
      vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
      vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
   END DO
 END DO
      	 
! vorticity = vorticity*area at this stage
! It is correct only on myDim nodes
DO n=1,myDim_nod2D
                              !! n=myList_nod2D(m)
   DO nz=1,nlevels_nod2D(n)-1
      vorticity(nz,n)=vorticity(nz,n)/area(nz,n)
   END DO
END DO      
 call exchange_nod(vorticity)
! Now it the relative vorticity known on neighbors too

end subroutine relative_vorticity
! ==========================================================================
subroutine compute_vel_rhs_vinv !vector invariant
USE o_PARAM
USE o_ARRAYS
USE o_MESH
USE g_PARSUP
USE g_CONFIG
  use g_comm_auto
IMPLICIT NONE
integer           :: n, n1, nz, elem, elnodes(3), nl1, j
real(kind=WP)     :: a, b, c, da, db, dc, dg, ff(3), gg, eta(3), pre(3), Fx, Fy,w
real(kind=WP)     :: uvert(nl,2), umean, vmean, friction
logical, save     :: lfirst=.true.
real(kind=WP)     :: KE_node(nl-1,myDim_nod2D+eDim_nod2D)
real(kind=WP)     :: dZ_inv(2:nl-1), dzbar_inv(nl-1), elem_area_inv
real(kind=WP)     :: density0_inv = 1./density_0

uvert=0.0_WP 

! ======================
! Kinetic energy at nodes:
! ======================
 

KE_node(:,:)=0d0 

DO elem=1, myDim_elem2D
                            !! elem=myList_elem2D(m)
   elnodes=elem2D_nodes(:,elem)
   DO j=1,3                !NR interchange loops => nz-loop vectorizes
      DO nz=1,nlevels(elem)-1 
         KE_node(nz,elnodes(j)) = KE_node(nz,elnodes(j))+(UV(1,nz,elem)*UV(1,nz,elem) &
              +UV(2,nz,elem)*UV(2,nz,elem))*elem_area(elem)  !NR/6.0_8 below
      END DO
   END DO
END DO

DO n=1,myDim_nod2D
                            !! n=myList_nod2D(m)
   DO nz=1, nlevels_nod2D(n)-1
   !DO nz=1, nl-1
      KE_node(nz,n)=KE_node(nz,n)/(6._WP*area(nz,n))  !NR divide by 6 here
   END DO
END DO   
! Set the kinetic energy to zero at lateral walls:
DO n=1,myDim_edge2D
                            !! n=myList_edge2D(m)
   if(myList_edge2D(n) > edge2D_in) then
      elnodes(1:2)=edges(:,n)
      KE_node(:,elnodes(1:2))=0.0
   endif
end DO   
 
 call exchange_nod(KE_node)
! Now gradients of KE will be correct on myDim_elem2D

! ==================
! AB contribution from the old time step
! ==================
 Do elem=1, myDim_elem2D          !! P (a)
                                  !! elem=myList_elem2D(m)
    DO nz=1,nl-1 
       UV_rhs(1,nz,elem)=-(0.5+epsilon)*UV_rhsAB(1,nz,elem)   
       UV_rhs(2,nz,elem)=-(0.5+epsilon)*UV_rhsAB(2,nz,elem)
    END DO
 END DO 

 call relative_vorticity 
! ====================
! Sea level and pressure contribution   -\nabla(g\eta +hpressure/rho_0+V^2/2)
! and the Coriolis force (elemental part)
! ====================

!DS KE_node=0.		!DS
!DS vorticity=0. 	!DS
DO elem=1,  myDim_elem2D          !! P (b)  elem=1,elem2D 
                                  !! elem=myList_elem2D(m)
   elnodes = elem2D_nodes(:,elem)
   eta = g*eta_n(elnodes)
   gg = elem_area(elem)
   ff = coriolis_node(elnodes)

   DO nz=1,nlevels(elem)-1
      pre = -(eta + hpressure(nz,elnodes)*density0_inv)
      Fx = sum(gradient_sca(1:3,elem)*pre)
      Fy = sum(gradient_sca(4:6,elem)*pre)
      UV_rhs(1,nz,elem) = UV_rhs(1,nz,elem)+Fx*gg 
      UV_rhs(2,nz,elem) = UV_rhs(2,nz,elem)+Fy*gg 

      pre = -KE_node(nz,elnodes)
      Fx = sum(gradient_sca(1:3,elem)*pre)
      Fy = sum(gradient_sca(4:6,elem)*pre)
   
      da = UV(2,nz,elem)*sum(ff+vorticity(nz,elnodes))/3.0_8
      db =-UV(1,nz,elem)*sum(ff+vorticity(nz,elnodes))/3.0_8

      UV_rhsAB(1,nz,elem)=(da+Fx)*gg
      UV_rhsAB(2,nz,elem)=(db+Fy)*gg

   END DO   
END DO
! =======================
! Compute w du/dz at elements: wdu/dz=d(wu)/dz-udw/dz
! The central estimate of u in the flux term will correspond  to energy
! conservation
! =======================

!NR precompute
DO nz=2,nl-1
   dZ_inv(nz) = 1./(Z(nz-1)-Z(nz))
ENDDO
DO nz=1,nl-1
   dzbar_inv(nz) = 1./(zbar(nz)-zbar(nz+1))
END DO

!DO elem=1, myDim_elem2D                  
!                                  !! elem=myList_elem2D(m)
!  elnodes=elem2D_nodes(:,elem)
!  nl1=nlevels(elem)-1
!
!  uvert(1,1:2)=0d0
!  uvert(nl1+1,1:2)=0d0
!
!  DO nz=2, nl1
!     w=sum(Wvel(nz,elnodes))/3.0_8
!     umean=0.5_8*(UV(1,nz-1,elem)+UV(1,nz,elem))
!     vmean=0.5_8*(UV(2,nz-1,elem)+UV(2,nz,elem))
!     uvert(nz,1)=-umean*w
!     uvert(nz,2)=-vmean*w
!  END DO
!  DO nz=1,nl1
!     da=sum(Wvel(nz,elnodes)-Wvel(nz+1,elnodes))/3.0_8
!     UV_rhsAB(1,nz,elem) = UV_rhsAB(1,nz,elem) + (uvert(nz,1)-uvert(nz+1,1)+&
!          da*UV(1,nz,elem))*elem_area(elem)*dzbar_inv(nz) !/(zbar(nz)-zbar(nz+1))
!     UV_rhsAB(2,nz,elem)=UV_rhsAB(2,nz,elem)+(uvert(nz,2)-uvert(nz+1,2)+&
!          da*UV(2,nz,elem))*elem_area(elem)*dzbar_inv(nz) !/(zbar(nz)-zbar(nz+1))
!     
!  END DO
!END DO


DO elem=1, myDim_elem2D                  
                                  !! elem=myList_elem2D(m)
  elnodes=elem2D_nodes(:,elem)
  nl1=nlevels(elem)-1

!  w=sum(Wvel(2, elnodes))/3.0_8
!  w=min(abs(w), 0.0001)*sign(1.0_8, w)
  uvert(1,1)=w*(UV(1,1,elem)-UV(1,2,elem))*dZ_inv(2)*0.5_8
  uvert(1,2)=w*(UV(2,1,elem)-UV(2,2,elem))*dZ_inv(2)*0.5_8

!  w=sum(Wvel(nl1, elnodes))/3.0_8
!  w=min(abs(w), 0.0001)*sign(1.0_8, w)
  uvert(nl1,1)=w*(UV(1,nl1-1,elem)-UV(1,nl1,elem))*dZ_inv(nl1)*0.5_8
  uvert(nl1,2)=w*(UV(2,nl1-1,elem)-UV(2,nl1,elem))*dZ_inv(nl1)*0.5_8


  DO nz=2, nl1-1
!     w=sum(Wvel(nz,elnodes)+Wvel(nz+1,elnodes))/6.0_8
!     w=min(abs(w), 0.0001)*sign(1.0_8, w)
     if (w >= 0.) then
        uvert(nz,1)=w*(UV(1,nz,elem)-UV(1,nz+1,elem))*dZ_inv(nz+1)
        uvert(nz,2)=w*(UV(2,nz,elem)-UV(2,nz+1,elem))*dZ_inv(nz+1)
     else
        uvert(nz,1)=w*(UV(1,nz-1,elem)-UV(1,nz,elem))*dZ_inv(nz)
        uvert(nz,2)=w*(UV(2,nz-1,elem)-UV(2,nz,elem))*dZ_inv(nz)
     end if
  END DO
  UV_rhsAB(1,1:nl1,elem) = UV_rhsAB(1,1:nl1,elem) - uvert(1:nl1,1)*elem_area(elem)
  UV_rhsAB(2,1:nl1,elem) = UV_rhsAB(2,1:nl1,elem) - uvert(1:nl1,2)*elem_area(elem)

END DO

! =======================
! Update the rhs   
! =======================
gg=(1.5_8+epsilon)
if(lfirst.and.(.not.r_restart)) then
   gg=1.0
   lfirst=.false.
end if

DO elem=1, myDim_elem2D                    !! P(e) elem=1, elem2D    
   !! elem=myList_elem2D(m)
   elem_area_inv = dt/elem_area(elem)
   DO nz=1,nlevels(elem)-1   
      UV_rhs(1,nz,elem)= (UV_rhs(1,nz,elem)+UV_rhsAB(1,nz,elem)*gg) *elem_area_inv
      UV_rhs(2,nz,elem)= (UV_rhs(2,nz,elem)+UV_rhsAB(2,nz,elem)*gg) *elem_area_inv
   END DO
END DO
! U_rhs contains all contributions to velocity from old time steps   

end subroutine compute_vel_rhs_vinv
