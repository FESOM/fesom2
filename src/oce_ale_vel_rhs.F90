module momentum_adv_scalar_interface
  interface
    subroutine momentum_adv_scalar(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module

!
!
!_______________________________________________________________________________
subroutine compute_vel_rhs(mesh)
    use MOD_MESH
    use o_ARRAYS
    use i_ARRAYS
    use i_therm_param
    use o_PARAM
    use g_PARSUP
    use g_CONFIG
    use g_forcing_param, only: use_virt_salt
    use g_comm_auto
    use momentum_adv_scalar_interface
    
    implicit none 
    type(t_mesh), intent(in) , target :: mesh   
    integer                  :: elem, elnodes(3), nz 
    real(kind=WP)            :: eta(3), ff, mm 
    real(kind=WP)            :: Fx, Fy, pre(3)
    logical, save            :: lfirst=.true.
    real(kind=WP)            :: t1, t2, t3, t4
    real(kind=WP)            :: p_ice(3)
    integer                  :: use_pice

#include "associate_mesh.h"

    t1=MPI_Wtime()
    use_pice=0
    if (use_floatice .and.  .not. trim(which_ale)=='linfs') use_pice=1
    
    do elem=1, myDim_elem2D
        !___________________________________________________________________________
        ! Take care of the AB part
        do nz=1,nl-1 
            UV_rhs(1,nz,elem)=-(0.5_WP+epsilon)*UV_rhsAB(1,nz,elem)   
            UV_rhs(2,nz,elem)=-(0.5_WP+epsilon)*UV_rhsAB(2,nz,elem)
        end do
        
        !___________________________________________________________________________
        ! Sea level and pressure contribution   -\nabla(\eta +hpressure/rho_0)
        ! and the Coriolis force + metric terms
        elnodes=elem2D_nodes(:,elem)
        
        !  eta=g*eta_n(elnodes)*(1-theta)        !! this place needs update (1-theta)!!!
        eta = g*eta_n(elnodes)   
        
        ff  = coriolis(elem)*elem_area(elem)
        !mm=metric_factor(elem)*gg
        
        !___________________________________________________________________________
        ! in case of ALE zlevel and zstar add pressure from ice to atmospheric pressure
        ! to account for floating ice
        if (use_pice > 0) then
            p_ice = 0.0_WP
            p_ice = (m_ice(elnodes)*rhoice+m_snow(elnodes)*rhosno)*inv_rhowat
            ! limit maximum ice loading like in FESOM1.4
            p_ice = g*min(p_ice,max_ice_loading)
        else
            p_ice = 0.0_WP
        endif
        
        !____________________________________________________________________________
        ! apply pressure gradient force, as well as contributions from gradient of 
        ! the sea surface height as well as ice pressure in case of floating sea ice
        ! to velocity rhs
        pre = -(eta+p_ice)!+atmospheric pressure etc.
        Fx  = sum(gradient_sca(1:3,elem)*pre)
        Fy  = sum(gradient_sca(4:6,elem)*pre)
        do nz=1,nlevels(elem)-1
            ! add pressure gradient terms
            UV_rhs(1,nz,elem)   = UV_rhs(1,nz,elem) + (Fx-pgf_x(nz,elem))*elem_area(elem) 
            UV_rhs(2,nz,elem)   = UV_rhs(2,nz,elem) + (Fy-pgf_y(nz,elem))*elem_area(elem)
            
            ! add coriolis force
            UV_rhsAB(1,nz,elem) = UV(2,nz,elem)*ff! + mm*UV(1,nz,elem)*UV(2,nz,elem)
            UV_rhsAB(2,nz,elem) =-UV(1,nz,elem)*ff! - mm*UV(1,nz,elem)*UV(2,nz,elem)
        end do
    end do

    t2=MPI_Wtime() 
    !___________________________________________________________________________
    ! advection
    if (mom_adv==1) then
       if (mype==0) write(*,*) 'in moment not adapted mom_adv advection typ for ALE, check your namelist'
       call par_ex(1)
    elseif (mom_adv==2) then
       call momentum_adv_scalar(mesh)
    end if
    t3=MPI_Wtime() 

    !___________________________________________________________________________
    ! Update the rhs   
    ff=(1.5_WP+epsilon)
    if (lfirst.and.(.not.r_restart)) then
        ff=1.0_WP
        lfirst=.false.
    end if

    do elem=1, myDim_elem2D
        do nz=1,nlevels(elem)-1   
            UV_rhs(1,nz,elem)=dt*(UV_rhs(1,nz,elem)+UV_rhsAB(1,nz,elem)*ff)/elem_area(elem)
            UV_rhs(2,nz,elem)=dt*(UV_rhs(2,nz,elem)+UV_rhsAB(2,nz,elem)*ff)/elem_area(elem)
        end do
    end do
    ! =======================  
    ! U_rhs contains all contributions to velocity from old time steps   
    ! =======================
    t4=MPI_Wtime() 
    ! if (mod(mstep,logfile_outfreq)==0 .and. mype==0) then
    !    write(*,*) 'Momentum:   ', t4-t1
    !    write(*,*) 'pres., Cor: ', t2-t1
    !    write(*,*) 'h adv       ', t3-t2
    !    write(*,*) 'vert. part  ', t4-t3
    ! end if     
END SUBROUTINE compute_vel_rhs
! ===================================================================
!
! Momentum advection on scalar control volumes with ALE adaption--> exchange zinv(nz)
! against hnode(nz,node)
!_______________________________________________________________________________
subroutine momentum_adv_scalar(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
use g_comm_auto
IMPLICIT NONE

type(t_mesh), intent(in) , target :: mesh
integer                  :: n, nz, el1, el2
integer                  :: nl1, nl2, nod(2), el, ed, k, nle
real(kind=WP)            :: un1(1:mesh%nl-1), un2(1:mesh%nl-1)
real(kind=WP)            :: wu(1:mesh%nl), wv(1:mesh%nl)
real(kind=WP)            :: Unode_rhs(2,mesh%nl-1,myDim_nod2d+eDim_nod2D)

#include "associate_mesh.h"

!_______________________________________________________________________________
do n=1,myDim_nod2d
   nl1 = nlevels_nod2D(n)-1
   wu(1:nl1+1) = 0._WP
   wv(1:nl1+1) = 0._WP

   do k=1,nod_in_elem2D_num(n)
      el = nod_in_elem2D(k,n)
      nle = nlevels(el)-1
   !___________________________________________________________________________
   ! The vertical part for each element is collected
      wu(1) = wu(1) + UV(1,1,el)*elem_area(el)
      wv(1) = wv(1) + UV(2,1,el)*elem_area(el)
             
      wu(2:nle) = wu(2:nle) + 0.5_WP*(UV(1,2:nle,el)+UV(1,1:nle-1,el))*elem_area(el)
      wv(2:nle) = wv(2:nle) + 0.5_WP*(UV(2,2:nle,el)+UV(2,1:nle-1,el))*elem_area(el)
   enddo

   wu(1:nl1) = wu(1:nl1)*Wvel_e(1:nl1,n)
   wv(1:nl1) = wv(1:nl1)*Wvel_e(1:nl1,n)

   do nz=1,nl1
      ! Here 1/3 because 1/3 of the area is related to the node
      Unode_rhs(1,nz,n) = - (wu(nz) - wu(nz+1) ) / (3._WP*hnode(nz,n)) 
      Unode_rhs(2,nz,n) = - (wv(nz) - wv(nz+1) ) / (3._WP*hnode(nz,n)) 
      
   enddo

   ! To get a clean checksum, set the remaining values to zero
   Unode_rhs(1:2,nl1+1:nl-1,n) = 0._WP

end do


!_______________________________________________________________________________
DO ed=1, myDim_edge2D
   nod = edges(:,ed)   
   el1  = edge_tri(1,ed)   
   el2  = edge_tri(2,ed)
   nl1 = nlevels(el1)-1

   !___________________________________________________________________________
   ! The horizontal part
   un1(1:nl1) = UV(2,1:nl1,el1)*edge_cross_dxdy(1,ed)   &
              - UV(1,1:nl1,el1)*edge_cross_dxdy(2,ed)  
   
   !___________________________________________________________________________
   if (el2>0) then
      nl2 = nlevels(el2)-1
      
      un2(1:nl2) = - UV(2,1:nl2,el2)*edge_cross_dxdy(3,ed) &
                   + UV(1,1:nl2,el2)*edge_cross_dxdy(4,ed)

      ! fill with zeros to combine the loops
      ! Usually, no or only a very few levels have to be filled. In this case, 
      ! computing "zeros" is cheaper than the loop overhead.
      
      un1(nl1+1:max(nl1,nl2)) = 0._WP
      un2(nl2+1:max(nl1,nl2)) = 0._WP

      ! first edge node
      ! Do not calculate on Halo nodes, as the result will not be used. 
      ! The "if" is cheaper than the avoided computiations.
      if (nod(1) <= myDim_nod2d) then
         do nz=1, max(nl1,nl2)
            
            Unode_rhs(1,nz,nod(1)) = Unode_rhs(1,nz,nod(1)) + un1(nz)*UV(1,nz,el1) + un2(nz)*UV(1,nz,el2) 
            Unode_rhs(2,nz,nod(1)) = Unode_rhs(2,nz,nod(1)) + un1(nz)*UV(2,nz,el1) + un2(nz)*UV(2,nz,el2)
         end do
      endif
      
      if (nod(2) <= myDim_nod2d) then
         do nz=1, max(nl1,nl2)
            Unode_rhs(1,nz,nod(2)) = Unode_rhs(1,nz,nod(2)) - un1(nz)*UV(1,nz,el1) - un2(nz)*UV(1,nz,el2)
            Unode_rhs(2,nz,nod(2)) = Unode_rhs(2,nz,nod(2)) - un1(nz)*UV(2,nz,el1) - un2(nz)*UV(2,nz,el2)
         end do
      endif


   else  ! ed is a boundary edge, there is only the contribution from el1
      if (nod(1) <= myDim_nod2d) then
         do nz=1, nl1
            
            Unode_rhs(1,nz,nod(1)) = Unode_rhs(1,nz,nod(1)) + un1(nz)*UV(1,nz,el1)
            Unode_rhs(2,nz,nod(1)) = Unode_rhs(2,nz,nod(1)) + un1(nz)*UV(2,nz,el1)
         end do
      endif
      ! second edge node
      if  (nod(2) <= myDim_nod2d) then
         do nz=1, nl1
            Unode_rhs(1,nz,nod(2)) = Unode_rhs(1,nz,nod(2)) - un1(nz)*UV(1,nz,el1)
            Unode_rhs(2,nz,nod(2)) = Unode_rhs(2,nz,nod(2)) - un1(nz)*UV(2,nz,el1)
         end do
      endif
   endif

end do

!_______________________________________________________________________________
do n=1,myDim_nod2d
   nl1 = nlevels_nod2D(n)-1
   
   Unode_rhs(1,1:nl1,n) = Unode_rhs(1,1:nl1,n) *area_inv(1:nl1,n)
   Unode_rhs(2,1:nl1,n) = Unode_rhs(2,1:nl1,n) *area_inv(1:nl1,n)
   
end do

!_______________________________________________________________________________
call exchange_nod(Unode_rhs)

!_______________________________________________________________________________
do el=1, myDim_elem2D
   nl1 = nlevels(el)-1
   UV_rhsAB(1:2,1:nl1,el) = UV_rhsAB(1:2,1:nl1,el) &
        + elem_area(el)*(Unode_rhs(1:2,1:nl1,elem2D_nodes(1,el)) &
        + Unode_rhs(1:2,1:nl1,elem2D_nodes(2,el)) & 
        + Unode_rhs(1:2,1:nl1,elem2D_nodes(3,el))) / 3.0_WP
   
end do
end subroutine momentum_adv_scalar


! ===================================================================

