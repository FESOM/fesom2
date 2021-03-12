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
    use g_forcing_arrays, only: press_air
    use g_comm_auto
    use g_sbf, only: l_mslp
    use momentum_adv_scalar_interface
    
    implicit none 
    type(t_mesh), intent(in) , target :: mesh   
    integer                  :: elem, elnodes(3), nz, nzmax, nzmin 
    real(kind=WP)            :: ff, mm 
    real(kind=WP)            :: Fx, Fy, pre(3)
    logical, save            :: lfirst=.true.
    real(kind=WP)            :: t1, t2, t3, t4
    real(kind=WP)            :: p_ice(3), p_air(3), p_eta(3)
    integer                  :: use_pice

#include "associate_mesh.h"

    t1=MPI_Wtime()
    use_pice=0
    if (use_floatice .and.  .not. trim(which_ale)=='linfs') use_pice=1
    
    do elem=1, myDim_elem2D
        nzmax = nlevels(elem)
        nzmin = ulevels(elem)
        !___________________________________________________________________________
        ! Take care of the AB part
        !!PS do nz=1,nl-1
        do nz=nzmin,nzmax-1
            UV_rhs(1,nz,elem)=-(0.5_WP+epsilon)*UV_rhsAB(1,nz,elem)   
            UV_rhs(2,nz,elem)=-(0.5_WP+epsilon)*UV_rhsAB(2,nz,elem)
        end do
        
        !___________________________________________________________________________
        ! Sea level and pressure contribution   -\nabla(\eta +hpressure/rho_0)
        ! and the Coriolis force + metric terms
        elnodes=elem2D_nodes(:,elem)
        
        !  p_eta=g*eta_n(elnodes)*(1-theta)        !! this place needs update (1-theta)!!!
        p_eta = g*eta_n(elnodes)   
        
        ff  = coriolis(elem)*elem_area(elem)
        !mm=metric_factor(elem)*gg
        
        !___________________________________________________________________________
        ! contribution from sea level pressure
        if (l_mslp) then
            p_air = press_air(elnodes)/1000
        else                   !|-> convert press_air from: Pa--> bar)
            p_air = 0.0_WP
        end if
        
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

        pre = -(p_eta+p_ice+p_air)

        if (use_global_tides) then
           pre=pre-ssh_gp(elnodes)
        end if

        Fx  = sum(gradient_sca(1:3,elem)*pre)
        Fy  = sum(gradient_sca(4:6,elem)*pre)
        !!PS do nz=1,nlevels(elem)-1
        do nz=nzmin,nzmax-1
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
        nzmax = nlevels(elem)
        nzmin = ulevels(elem)
        !!PS do nz=1,nlevels(elem)-1
        do nz=nzmin,nzmax-1
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
integer                  :: nl1, nl2, ul1, ul2, nod(2), el, ed, k, nle, ule
real(kind=WP)            :: un1(1:mesh%nl-1), un2(1:mesh%nl-1)
real(kind=WP)            :: wu(1:mesh%nl), wv(1:mesh%nl)
real(kind=WP)            :: Unode_rhs(2,mesh%nl-1,myDim_nod2d+eDim_nod2D)

#include "associate_mesh.h"

    !___________________________________________________________________________
    ! 1st. compute vertical momentum advection component: w * du/dz, w*dv/dz
    do n=1,myDim_nod2d
        nl1 = nlevels_nod2D(n)-1
        ul1 = ulevels_nod2D(n)
        wu(1:nl1+1) = 0._WP
        wv(1:nl1+1) = 0._WP
        
        !_______________________________________________________________________
        ! loop over adjacent elements of vertice n 
        do k=1,nod_in_elem2D_num(n)
            el = nod_in_elem2D(k,n)
            !___________________________________________________________________
            nle = nlevels(el)-1
            ule = ulevels(el)
            
            !___________________________________________________________________
            ! accumulate horizontal velocities at full depth levels (top and 
            ! bottom faces of prism) 
            ! account here also for boundary condition below cavity --> 
            ! horizontal velocity at cavity-ocean interce ule (if ule>1) must be  
            ! zero ???
            if (ule==1) then
                wu(ule) = wu(ule) + UV(1,ule,el)*elem_area(el)
                wv(ule) = wv(ule) + UV(2,ule,el)*elem_area(el)
            end if 
            
            ! interpolate horizontal velocity from mid-depth levels to full
            ! depth levels of upper and lower prism faces and average over adjacent
            ! elements of vertice n
            wu(ule+1:nle) = wu(ule+1:nle) + 0.5_WP*(UV(1,ule+1:nle,el)+UV(1,ule:nle-1,el))*elem_area(el)
            wv(ule+1:nle) = wv(ule+1:nle) + 0.5_WP*(UV(2,ule+1:nle,el)+UV(2,ule:nle-1,el))*elem_area(el)
        enddo
        
        !_______________________________________________________________________
        ! multiply w*du and w*dv
        wu(ul1:nl1) = wu(ul1:nl1)*Wvel_e(ul1:nl1,n)
        wv(ul1:nl1) = wv(ul1:nl1)*Wvel_e(ul1:nl1,n)
        
        !_______________________________________________________________________
        ! compute w*du/dz, w*dv/dz
        do nz=ul1,nl1
!!PS             if (ul1>1) write(*,*) mype, wu(ul1:nl1)
            ! Here 1/3 because 1/3 of the area is related to the node --> comes from
            ! averaging the elemental velocities
            Unode_rhs(1,nz,n) = - (wu(nz) - wu(nz+1) ) / (3._WP*hnode(nz,n)) 
            Unode_rhs(2,nz,n) = - (wv(nz) - wv(nz+1) ) / (3._WP*hnode(nz,n)) 
            
        enddo
        
        !_______________________________________________________________________
        ! To get a clean checksum, set the remaining values to zero
        Unode_rhs(1:2,nl1+1:nl-1,n) = 0._WP
        Unode_rhs(1:2,1:ul1-1   ,n) = 0._WP
    end do


    !___________________________________________________________________________
    ! 2nd. compute horizontal advection component: u*du/dx, u*dv/dx & v*du/dy, v*dv/dy
    ! loop over triangle edges
    do ed=1, myDim_edge2D
        nod = edges(:,ed)   
        el1 = edge_tri(1,ed)   
        el2 = edge_tri(2,ed)
        nl1 = nlevels(el1)-1
        ul1 = ulevels(el1)
        
        !_______________________________________________________________________
        ! compute horizontal normal velocity with respect to the edge from triangle 
        ! centroid towards triangel edge mid-pointe for element el1
        !                     .o.    
        !                   ./   \.                
        !                 ./  el1  \.   
        !               ./     x     \. 
        !             ./       |-------\.-----------------edge_cross_dxdy(1:2,ed) --> (dx,dy)
        !            /         |->n_vec  \                  
        !    nod(1) o----------O----------o nod(2)   
        !            \.        |->n_vec ./
        !              \.      |------./------------------edge_cross_dxdy(3:4,ed) --> (dx,dy)
        !                \.    x    ./
        !                  \. el2 ./
        !                    \. ./  
        !                      °
        un1(ul1:nl1) =   UV(2,ul1:nl1,el1)*edge_cross_dxdy(1,ed)   &
                       - UV(1,ul1:nl1,el1)*edge_cross_dxdy(2,ed)  
                       
        !_______________________________________________________________________
        ! compute horizontal normal velocity with respect to the edge from triangle 
        ! centroid towards triangel edge mid-pointe for element el2 when it is valid
        ! --> if its a boundary triangle el2 will be not valid
        if (el2>0) then ! --> el2 is valid element
            nl2 = nlevels(el2)-1
            ul2 = ulevels(el2)
            
            un2(ul2:nl2) = - UV(2,ul2:nl2,el2)*edge_cross_dxdy(3,ed) &
                            + UV(1,ul2:nl2,el2)*edge_cross_dxdy(4,ed)
            
            ! fill with zeros to combine the loops
            ! Usually, no or only a very few levels have to be filled. In this case, 
            ! computing "zeros" is cheaper than the loop overhead.
            un1(nl1+1:max(nl1,nl2)) = 0._WP
            un2(nl2+1:max(nl1,nl2)) = 0._WP
            un1(1:ul1-1)            = 0._WP
            un2(1:ul2-1)            = 0._WP
            
            ! first edge node
            ! Do not calculate on Halo nodes, as the result will not be used. 
            ! The "if" is cheaper than the avoided computiations.
            if (nod(1) <= myDim_nod2d) then
                do nz=min(ul1,ul2), max(nl1,nl2)
                    ! add w*du/dz+(u*du/dx+v*du/dy) & w*dv/dz+(u*dv/dx+v*dv/dy)
                    Unode_rhs(1,nz,nod(1)) = Unode_rhs(1,nz,nod(1)) + un1(nz)*UV(1,nz,el1) + un2(nz)*UV(1,nz,el2) 
                    Unode_rhs(2,nz,nod(1)) = Unode_rhs(2,nz,nod(1)) + un1(nz)*UV(2,nz,el1) + un2(nz)*UV(2,nz,el2)
                end do
            endif
            
            ! second edge node
            if (nod(2) <= myDim_nod2d) then
                do nz=min(ul1,ul2), max(nl1,nl2)
                    ! add w*du/dz+(u*du/dx+v*du/dy) & w*dv/dz+(u*dv/dx+v*dv/dy)
                    Unode_rhs(1,nz,nod(2)) = Unode_rhs(1,nz,nod(2)) - un1(nz)*UV(1,nz,el1) - un2(nz)*UV(1,nz,el2)
                    Unode_rhs(2,nz,nod(2)) = Unode_rhs(2,nz,nod(2)) - un1(nz)*UV(2,nz,el1) - un2(nz)*UV(2,nz,el2)
                end do
            endif
            
        else  ! el2 is not a valid element --> ed is a boundary edge, there is only the contribution from el1
            ! first edge node
            if (nod(1) <= myDim_nod2d) then
                do nz=ul1, nl1
                    ! add w*du/dz+(u*du/dx+v*du/dy) & w*dv/dz+(u*dv/dx+v*dv/dy)
                    Unode_rhs(1,nz,nod(1)) = Unode_rhs(1,nz,nod(1)) + un1(nz)*UV(1,nz,el1)
                    Unode_rhs(2,nz,nod(1)) = Unode_rhs(2,nz,nod(1)) + un1(nz)*UV(2,nz,el1)
                end do ! --> do nz=ul1, nl1
            endif 
            
            ! second edge node
            if  (nod(2) <= myDim_nod2d) then
                !!PS do nz=1, nl1
                do nz=ul1, nl1
                    ! add w*du/dz+(u*du/dx+v*du/dy) & w*dv/dz+(u*dv/dx+v*dv/dy)
                    Unode_rhs(1,nz,nod(2)) = Unode_rhs(1,nz,nod(2)) - un1(nz)*UV(1,nz,el1)
                    Unode_rhs(2,nz,nod(2)) = Unode_rhs(2,nz,nod(2)) - un1(nz)*UV(2,nz,el1)
                end do ! --> do nz=ul1, nl1
            endif
        endif ! --> if (el2>0) then
    end do ! --> do ed=1, myDim_edge2D

    !___________________________________________________________________________
    ! divide total nodal advection by scalar area
    do n=1,myDim_nod2d
        nl1 = nlevels_nod2D(n)-1
        ul1 = ulevels_nod2D(n)
!!PS         Unode_rhs(1,ul1:nl1,n) = Unode_rhs(1,ul1:nl1,n) *area_inv(ul1:nl1,n) ! --> TEST_cavity
!!PS         Unode_rhs(2,ul1:nl1,n) = Unode_rhs(2,ul1:nl1,n) *area_inv(ul1:nl1,n) ! --> TEST_cavity
        Unode_rhs(1,ul1:nl1,n) = Unode_rhs(1,ul1:nl1,n) *areasvol_inv(ul1:nl1,n)
        Unode_rhs(2,ul1:nl1,n) = Unode_rhs(2,ul1:nl1,n) *areasvol_inv(ul1:nl1,n)
    end do !-->do n=1,myDim_nod2d

    !___________________________________________________________________________
    call exchange_nod(Unode_rhs)

    !___________________________________________________________________________
    ! convert total nodal advection from vertice --> elements
    do el=1, myDim_elem2D
        nl1 = nlevels(el)-1
        ul1 = ulevels(el)
        UV_rhsAB(1:2,ul1:nl1,el) = UV_rhsAB(1:2,ul1:nl1,el) &
                + elem_area(el)*(Unode_rhs(1:2,ul1:nl1,elem2D_nodes(1,el)) &
                + Unode_rhs(1:2,ul1:nl1,elem2D_nodes(2,el)) & 
                + Unode_rhs(1:2,ul1:nl1,elem2D_nodes(3,el))) / 3.0_WP     
    
    end do ! --> do el=1, myDim_elem2D
end subroutine momentum_adv_scalar


! ===================================================================

