module diff_part_hor_redi_interface
    interface
        subroutine diff_part_hor_redi(tracer, partit, mesh) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module diff_ver_part_expl_ale_interface
    interface
        subroutine diff_ver_part_expl_ale(tr_num, tracer, partit, mesh) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        integer       , intent(in)   , target :: tr_num
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module diff_ver_part_redi_expl_interface
    interface
        subroutine diff_ver_part_redi_expl(tracer, partit, mesh) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module diff_ver_part_impl_ale_interface
    interface
        subroutine diff_ver_part_impl_ale(tr_num, dynamics,  tracer, partit, mesh) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        use MOD_DYN
        integer       , intent(in)   , target :: tr_num
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module diff_tracers_ale_interface
    interface
        subroutine diff_tracers_ale(tr_num, dynamics, tracer, partit, mesh) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        use MOD_DYN
        integer       , intent(in),    target :: tr_num
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module bc_surface_interface
    interface
        function bc_surface(n, id, sval, nzmin, partit) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        integer , intent(in)                  :: n, id, nzmin
        type(t_partit), intent(inout), target :: partit
        real(kind=WP)                         :: bc_surface
        real(kind=WP), intent(in)             :: sval
        end function
    end interface
end module

module diff_part_bh_interface
    interface
        subroutine diff_part_bh(tr_num, dynamics, tracer, partit, mesh) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        use MOD_DYN
        integer       , intent(in)   , target :: tr_num
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module solve_tracers_ale_interface
    interface
        subroutine solve_tracers_ale(ice, dynamics, tracers, partit, mesh) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer     
        use MOD_DYN
        USE MOD_ICE
        type(t_ice)   , intent(in),    target :: ice
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module
!
!
!===============================================================================
! Driving routine    Here with ALE changes!!!
subroutine solve_tracers_ale(ice, dynamics, tracers, partit, mesh) 
    use g_config
    use o_PARAM, only: SPP, Fer_GM
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE MOD_ICE
    use mod_tracer
    use g_comm_auto
    use o_tracers
    use Toy_Channel_Soufflet
    use diff_tracers_ale_interface
    use oce_adv_tra_driver_interfaces    
    implicit none
    type(t_ice)   , intent(in)   , target    :: ice
    type(t_dyn)   , intent(inout), target    :: dynamics
    type(t_tracer), intent(inout), target    :: tracers
    type(t_partit), intent(inout), target    :: partit
    type(t_mesh)  , intent(in)   , target    :: mesh
    !___________________________________________________________________________
    integer                                  :: tr_num, node, elem, nzmax, nzmin
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, fer_UV
    real(kind=WP), dimension(:,:)  , pointer :: Wvel, Wvel_e, Wvel_i, fer_Wvel
    real(kind=WP), dimension(:,:)  , pointer :: del_ttf
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    UV     => dynamics%uv(:,:,:)
    Wvel   => dynamics%w(:,:)
    Wvel_e => dynamics%w_e(:,:)
    Wvel_i => dynamics%w_i(:,:)
    if (Fer_GM) then
        fer_UV     => dynamics%fer_uv(:,:,:)
        fer_Wvel   => dynamics%fer_w(:,:)
    end if 
    del_ttf => tracers%work%del_ttf

    !___________________________________________________________________________
    if (SPP) call cal_rejected_salt(ice, partit, mesh)
    if (SPP) call app_rejected_salt(tracers%data(2)%values, partit, mesh)
    
    !___________________________________________________________________________
    ! update 3D velocities with the bolus velocities:
    ! 1. bolus velocities are computed according to GM implementation after R. Ferrari et al., 2010
    ! 2. bolus velocities are used only for advecting tracers and shall be subtracted back afterwards
    if (Fer_GM) then
!$OMP PARALLEL DO
        do elem=1, myDim_elem2D+eDim_elem2D
           UV(:, :, elem)    =UV(:, :, elem) + fer_UV(:, :, elem)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do node=1, myDim_nod2D+eDim_nod2D
           Wvel_e(:, node)=Wvel_e(:, node)+fer_Wvel(:, node)
           Wvel  (:, node)=Wvel  (:, node)+fer_Wvel(:, node)
        end do
!$OMP END PARALLEL DO
    end if
    
    !___________________________________________________________________________
    ! loop over all tracers 
    do tr_num=1, tracers%num_tracers
        ! do tracer AB (Adams-Bashfort) interpolation only for advectiv part 
        ! needed
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call init_tracers_AB'//achar(27)//'[0m'
        call init_tracers_AB(tr_num, tracers, partit, mesh) 
        
        ! advect tracers
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call adv_tracers_ale'//achar(27)//'[0m'
        ! it will update del_ttf with contributions from horizontal and vertical advection parts (del_ttf_advhoriz and del_ttf_advvert)
        call do_oce_adv_tra(dt, UV, Wvel, Wvel_i, Wvel_e, tr_num, dynamics, tracers, partit, mesh)
        
        !___________________________________________________________________________
        ! update array for total tracer flux del_ttf with the fluxes from horizontal
        ! and vertical advection
!$OMP PARALLEL DO
        do node=1, myDim_nod2d
           tracers%work%del_ttf(:, node)=tracers%work%del_ttf(:, node)+tracers%work%del_ttf_advhoriz(:, node)+tracers%work%del_ttf_advvert(:, node)
        end do
!$OMP END PARALLEL DO
           !___________________________________________________________________________
           ! AB is not needed after the advection step. Initialize it with the current tracer before it is modified.
           ! call init_tracers_AB at the beginning of this loop will compute AB for the next time step then.
!$OMP PARALLEL DO
        do node=1, myDim_nod2d+eDim_nod2D
           tracers%data(tr_num)%valuesAB(:, node)=tracers%data(tr_num)%values(:, node) !DS: check that this is the right place!
        end do
!$OMP END PARALLEL DO
        ! diffuse tracers 
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call diff_tracers_ale'//achar(27)//'[0m'
        call diff_tracers_ale(tr_num, dynamics, tracers, partit, mesh) 

        ! relax to salt and temp climatology
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call relax_to_clim'//achar(27)//'[0m'
!       if ((toy_ocean) .AND. ((tr_num==1) .AND. (TRIM(which_toy)=="soufflet"))) then
        if ((toy_ocean) .AND. ((TRIM(which_toy)=="soufflet"))) then
            call relax_zonal_temp(tracers%data(1), partit, mesh) 
        else
            call relax_to_clim(tr_num, tracers, partit, mesh) 
        end if
        call exchange_nod(tracers%data(tr_num)%values(:,:), partit)
!$OMP BARRIER
    end do
    
    !___________________________________________________________________________
    ! 3D restoring for "passive" tracers
    !!!$OMPTODO: add OpenMP later, not needed right now!
    do tr_num=1, ptracers_restore_total
       tracers%data(ptracers_restore(tr_num)%locid)%values(:, ptracers_restore(tr_num)%ind2)=1.0_WP
    end do
    
    !___________________________________________________________________________
    ! subtract the the bolus velocities back from 3D velocities:
    if (Fer_GM) then
!$OMP PARALLEL DO
        do elem=1, myDim_elem2D+eDim_elem2D
           UV(:, :, elem)    =UV(:, :, elem) - fer_UV(:, :, elem)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do node=1, myDim_nod2D+eDim_nod2D
           Wvel_e(:, node)=Wvel_e(:, node)-fer_Wvel(:, node)
           Wvel  (:, node)=Wvel  (:, node)-fer_Wvel(:, node)
        end do
!$OMP END PARALLEL DO
    end if
    
    !___________________________________________________________________________
    ! to avoid crash with high salinities when coupled to atmosphere
    ! --> if we do only where (tr_arr(:,:,2) < 3._WP ) we also fill up the bottom 
    !     topogrpahy with values which are then writte into the output --> thats why
    !     do node=1,.... and tr_arr(node,1:nzmax,2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nzmin, nzmax)
    do node=1,myDim_nod2D+eDim_nod2D
        nzmax=nlevels_nod2D(node)-1
        nzmin=ulevels_nod2D(node)
        where (tracers%data(2)%values(nzmin:nzmax,node) > 45._WP)
               tracers%data(2)%values(nzmin:nzmax,node)=45._WP
        end where

        where (tracers%data(2)%values(nzmin:nzmax,node) < 3._WP )
               tracers%data(2)%values(nzmin:nzmax,node) = 3._WP
        end where        
    end do
!$OMP END PARALLEL DO
end subroutine solve_tracers_ale
!
!
!===============================================================================
subroutine diff_tracers_ale(tr_num, dynamics, tracers, partit, mesh) 
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use mod_tracer
    use MOD_DYN
    use o_arrays
    use o_tracers
    use diff_part_hor_redi_interface
    use diff_ver_part_expl_ale_interface
    use diff_ver_part_redi_expl_interface
    use diff_ver_part_impl_ale_interface
    use diff_part_bh_interface
    implicit none
    integer       , intent(in)   , target :: tr_num
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: n, nzmax, nzmin
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), pointer                :: del_ttf(:,:)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    del_ttf => tracers%work%del_ttf
    
    !___________________________________________________________________________
    ! do horizontal diffusiion
    ! write there also horizontal diffusion rhs to del_ttf which is equal the R_T^n 
    ! in danilovs srcipt
    ! includes Redi diffusivity if Redi=.true.
    call diff_part_hor_redi(tracers, partit, mesh)  ! seems to be ~9% faster than diff_part_hor
    
    !___________________________________________________________________________
    ! do vertical diffusion: explicit
    if (.not. tracers%data(tr_num)%i_vert_diff) call diff_ver_part_expl_ale(tr_num, tracers, partit, mesh) 
    ! A projection of horizontal Redi diffussivity onto vertical. This par contains horizontal
    ! derivatives and has to be computed explicitly!
    if (Redi) call diff_ver_part_redi_expl(tracers, partit, mesh)     

    !___________________________________________________________________________
    ! Update tracers --> calculate T* see Danilov et al. (2017)
    ! T* =  (dt*R_T^n + h^(n-0.5)*T^(n-0.5))/h^(n+0.5)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nzmin, nzmax)
    do n=1, myDim_nod2D 
        nzmax=nlevels_nod2D(n)-1
        nzmin=ulevels_nod2D(n)
        del_ttf(nzmin:nzmax,n)=del_ttf(nzmin:nzmax,n)+tracers%data(tr_num)%values(nzmin:nzmax,n)* &
                                    (hnode(nzmin:nzmax,n)-hnode_new(nzmin:nzmax,n))
        tracers%data(tr_num)%values(nzmin:nzmax,n)=tracers%data(tr_num)%values(nzmin:nzmax,n)+ &
                                    del_ttf(nzmin:nzmax,n)/hnode_new(nzmin:nzmax,n)
        ! WHY NOT ??? --> whats advantage of above --> tested it --> the upper 
        ! equation has a 30% smaller nummerical drift
        ! tr_arr(1:nzmax,n,tr_num)=(hnode(1:nzmax,n)*tr_arr(1:nzmax,n,tr_num)+ &
        !                           del_ttf(1:nzmax,n))/hnode_new(1:nzmax,n)
    end do
!$OMP END PARALLEL DO
    !___________________________________________________________________________
    if (tracers%data(tr_num)%i_vert_diff) then
        ! do vertical diffusion: implicite 
        call diff_ver_part_impl_ale(tr_num, dynamics, tracers, partit, mesh) 
    end if
    !We DO not set del_ttf to zero because it will not be used in this timestep anymore
    !init_tracers_AB will set it to zero for the next timestep
    if (tracers%data(tr_num)%smooth_bh_tra) then
       call diff_part_bh(tr_num, dynamics, tracers, partit, mesh)  ! alpply biharmonic diffusion (implemented as filter)                                                
    end if
end subroutine diff_tracers_ale
!
!
!===============================================================================
!Vertical diffusive flux(explicit scheme):                                                                            
subroutine diff_ver_part_expl_ale(tr_num, tracers, partit, mesh) 
    use o_ARRAYS
    use g_forcing_arrays
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_TRACER
    use g_config,only: dt
    implicit none 
    integer       , intent(in)   , target :: tr_num
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: n, nz, nl1, ul1
    real(kind=WP)            :: vd_flux(mesh%nl-1)
    real(kind=WP)            :: rdata, flux, rlx
    real(kind=WP)            :: zinv1
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), pointer   :: del_ttf(:,:)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    del_ttf => tracers%work%del_ttf
    
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nl1, ul1, vd_flux, rdata, flux, rlx, zinv1)
    !___________________________________________________________________________    
    do n=1, myDim_nod2D
        nl1=nlevels_nod2D(n)-1
        ul1=ulevels_nod2D(n)
        
        vd_flux=0._WP
        if (tracers%data(tr_num)%ID==1) then
            flux  = -heat_flux(n)/vcpw
            rdata =  Tsurf(n)
            rlx   =  surf_relax_T
        elseif (tracers%data(tr_num)%ID==2) then
            flux  =  virtual_salt(n)+relax_salt(n)- real_salt_flux(n)*is_nonlinfs
        else
            flux  = 0._WP
            rdata = 0._WP
            rlx=0._WP
        endif        
        !_______________________________________________________________________
        !Surface forcing
        vd_flux(ul1)= flux        
        do nz=ul1+1,nl1
            !___________________________________________________________________
            zinv1=1.0_WP/(Z_3d_n(nz-1,n)-Z_3d_n(nz,n))                       
            vd_flux(nz) = Kv(nz,n)*(tracers%data(tr_num)%values(nz-1,n)-tracers%data(tr_num)%values(nz,n))*zinv1*area(nz,n)            
        end do
        
        do nz=ul1,nl1-1
            del_ttf(nz,n) = del_ttf(nz,n) + (vd_flux(nz) - vd_flux(nz+1))/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))*dt/areasvol(nz,n)
        end do
        del_ttf(nl1,n) = del_ttf(nl1,n) + (vd_flux(nl1)/(zbar_3d_n(nl1,n)-zbar_3d_n(nl1+1,n)))*dt/areasvol(nl1,n)
        
    end do ! --> do n=1, myDim_nod2D
!$OMP END PARALLEL DO
end subroutine diff_ver_part_expl_ale
!
!
!===============================================================================
! vertical diffusivity augmented with Redi contribution [vertical flux of K(3,3)*d_zT]
subroutine diff_ver_part_impl_ale(tr_num, dynamics, tracers, partit, mesh) 
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_TRACER
    use MOD_DYN
    use o_PARAM
    use o_ARRAYS, only: Ki, Kv, heat_flux, water_flux, slope_tapered
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_CONFIG
    use g_forcing_arrays
    use o_mixing_KPP_mod !for ghats _GO_   
    use g_cvmix_kpp, only: kpp_nonlcltranspT, kpp_nonlcltranspS, kpp_oblmixc
    use bc_surface_interface
    implicit none
    integer       , intent(in)   , target :: tr_num
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=WP)            :: a(mesh%nl), b(mesh%nl), c(mesh%nl), tr(mesh%nl)
    real(kind=WP)            :: cp(mesh%nl), tp(mesh%nl)
    integer                  :: nz, n, nzmax, nzmin
    real(kind=WP)            :: m, zinv, dz
    real(kind=WP)            :: rsss, Ty, Ty1, c1, zinv1, zinv2, v_adv
    real(kind=WP)            :: isredi=0._WP
    logical                  :: do_wimpl=.true.
    real(kind=WP)            :: zbar_n(mesh%nl), z_n(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:), pointer :: trarr
    real(kind=WP), dimension(:,:), pointer :: Wvel_i
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    trarr  => tracers%data(tr_num)%values(:,:)
    Wvel_i => dynamics%w_i(:,:)
    
    !___________________________________________________________________________
    if ((trim(tracers%data(tr_num)%tra_adv_lim)=='FCT') .OR. (.not. dynamics%use_wsplit)) do_wimpl=.false.    
    if (Redi) isredi=1._WP
    Ty    =0.0_WP
    Ty1   =0.0_WP
    
    ! solve equation diffusion equation implicite part: 
    ! -->   h^(n+0.5)* (T^(n+0.5)-Tstar) = dt*( K_33*d/dz*(T^(n+0.5)-Tstar) + K_33*d/dz*Tstar )
    ! -->   Tnew = T^(n+0.5)-Tstar
    ! -->   h^(n+0.5)* (Tnew) = dt*(K_33*d/dz*Tnew) + K_33*dt*d/dz*Tstar 
    ! -->   h^(n+0.5)* (Tnew) = dt*(K_33*d/dz*Tnew) + RHS        
    ! -->   solve for T_new                                      
    ! -->   V_1 (Skalar Volume), A_1 (Area of edge),              .                 
    !       no Cavity A1==V1, yes Cavity A1 !=V1                 /I\ nvec_up (+1)
    !                                                             I
    !    ----------- zbar_1, A_1                             *----I----*
    ! Z_1 o T_1, V1                                          |\   I  ./|
    !    ----------- zbar_2, A_2                             | \   ./  |   Gaus Theorem:
    ! Z_2 o T_2, V2                                          |  \ /    |    --> Flux form
    !    ----------- zbar_3, A_3                             |   |     |    --> normal vec outwards facing
    ! Z_3 o T_3, V3                                          *---|-----*
    !    ----------- zbar_4                                   \  | I ./
    !        :                                                 \ | I/ 
    !                                                           \|/I   
    !                                                            * I
    !                                                             \I/
    !                                                              *  nvec_dwn (-1)
    ! --> 1st. solve homogenouse part:
    ! f(Tnew) = h^(n+0.5)* (Tnew) - dt*(K_33*dTnew/dz) = 0
    !
    ! --> 2nd. Difference Quotient at Tnew_i in flux form (Gaus Theorem, dont forget normal vectors!!!):  
    ! V_i*Tnew_i *h_i = -dt * [ K_33 * (Tnew_i-1 - Tnew_i)/(Z_i-1 - Z_i) * A_i * nvec_up 
    !                          +K_33 * (Tnew_i - Tnew_i+1)/(Z_i - Z_i+1) * A_i+1 * nvec_dwn ]
    !     Tnew_i *h_i = -dt * [ K_33 * (Tnew_i-1 - Tnew_i)/(Z_i-1 - Z_i) * A_i  /V_i * nvec_up 
    !                          +K_33 * (Tnew_i - Tnew_i+1)/(Z_i - Z_i+1) * A_i+1/V_i * nvec_dwn ]
    !
    ! --> 3rd. solve for coefficents a, b, c:
    ! f(Tnew) = [ a*dTnew_i-1 + b*dTnew_i + c*dTnew_i+1 ] 
    !
    !     df(Tnew)/dTnew_i-1 = a = -dt*K_33/(Z_i-1 - Z_i) * A_i/V_i * (nvec_up =1)
    ! 
    !     df(Tnew)/dTnew_i+1 = c =  dt * K_33 * 1/(Z_i - Z_i+1) * A_i+1/V_i * (nvec_dwn=-1)
    !                            = -dt * K_33 * 1/(Z_i - Z_i+1) * A_i+1/V_i
    !
    !     df(Tnew)/dTnew_i   = b = h_i + dt*K_33/(Z_i-1 - Z_i) * A_i/V_i   * (nvec_up=+1)
    !                                  - dt*K_33/(Z_i - Z_i+1) * A_i+1/V_i * (nvec_dwn=-1)
    !                            = h_i + dt*K_33/(Z_i-1 - Z_i) * A_i/V_i
    !                                  + dt*K_33/(Z_i - Z_i+1) * A_i+1/V_i
    !                            = h_i -(a+c)
    !
    ! --> 4th. solve inhomogenous part:
    ! [ a*dTnew_i-1 + b*dTnew_i + c*dTnew_i+1 ] = RHS/V_i
    !
    ! RHS     = K_33*dt*d/dz*Tstar
    !
    ! --> write as Difference Quotient in flux form
    ! RHS/V_i =  K_33 * dt * (Tstar_i-1 - Tstar_i)/(Z_i-1 - Z_i) * A_i/V_i   * (nvec_up=1)
    !          + K_33 * dt * (Tstar_i - Tstar_i+1)/(Z_i - Z_i+1) * A_i+1/V_i * (nvec_dwn=-1)
    !
    !         =  K_33*dt/(Z_i-1 - Z_i) * A_i/V_i   * Tstar_i-1
    !          - K_33*dt/(Z_i-1 - Z_i) * A_i/V_i   * Tstar_i
    !          - K_33*dt/(Z_i - Z_i+1) * A_i+1/V_i * Tstar_i
    !          + K_33*dt/(Z_i - Z_i+1) * A_i+1/V_i * Tstar_i+1
    !         
    !         = -a*Tstar_i-1 + (a+c)*Tstar_i - c * Tstar_i+1
    !                            |-> b = h_i - (a+c), a+c = h_i-b
    !
    !         = -a*Tstar_i-1 - (b-h_i)*Tstar_i - c * Tstar_i+1
    !
    ! --> 5th. solve for Tnew_i --> forward sweep algorithm --> see lower
    !  | b_1 c_1 ...            |   |dTnew_1|
    !  | a_2 b_2 c_2 ...        |   |dTnew_2|
    !  |     a_3 b_3 c_3 ...    | * |dTnew_3| = RHS/V_i
    !  |         a_4 b_4 c_4 ...|   |dTnew_4| 
    !  |              :         |   |   :   |   
    !
    ! --> a = -dt*K_33 / (Z_i-1 - Z_i) * A_i/V_i
    ! 
    ! --> c = -dt*K_33 / (Z_i - Z_i+1) * A_i+1/V_i
    !
    ! --> b = h^(n+0.5) -[ dt*K_33/(Z_i-1 - Z_i)*A_i/V_i + dt*K_33/(Z_i - Z_i+1) * A_i+1/V_i ] = -(a+c) + h^(n+0.5)
    
    !___________________________________________________________________________
    ! loop over local nodes

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, nzmax, nzmin, a, b, c, tr, cp, tp, m, zinv, dz, &
!$OMP                                         rsss, Ty, Ty1, c1, zinv1, zinv2, v_adv, zbar_n, z_n)
!$OMP DO
    do n=1,myDim_nod2D 
        ! initialise
        a  = 0.0_WP
        b  = 0.0_WP
        c  = 0.0_WP
        tr = 0.0_WP
        tp = 0.0_WP
        cp = 0.0_WP
        ! max. number of levels at node n
        nzmax=nlevels_nod2D(n)
        nzmin=ulevels_nod2D(n)
        !___________________________________________________________________________
        ! Here can not exchange zbar_n & Z_n with zbar_3d_n & Z_3d_n because  
        ! they be calculate from the actualized mesh with hnode_new
        ! calculate new zbar (depth of layers) and Z (mid depths of layers) 
        ! depending on layer thinkness over depth at node n
        ! Be carefull here vertical operation have to be done on NEW vertical mesh !!!
        zbar_n=0.0_WP
        Z_n=0.0_WP
        zbar_n(nzmax)=zbar_n_bot(n)
        Z_n(nzmax-1)=zbar_n(nzmax) + hnode_new(nzmax-1,n)/2.0_WP
        do nz=nzmax-1,nzmin+1,-1
            zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
            Z_n(nz-1)  = zbar_n(nz) + hnode_new(nz-1,n)/2.0_WP
        end do
        zbar_n(nzmin) = zbar_n(nzmin+1) + hnode_new(nzmin,n)
        !_______________________________________________________________________
        ! Regular part of coefficients: --> surface layer 
        nz=nzmin
        
        ! 1/dz(nz)
        zinv2=1.0_WP/(Z_n(nz)-Z_n(nz+1))
        zinv=1.0_WP*dt    ! no .../(zbar(1)-zbar(2)) because of  ALE
        
        ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
        Ty1= (Z_n(nz)     -zbar_n(nz+1))*zinv2 *slope_tapered(3,nz  ,n)**2*Ki(nz  ,n) + &
             (zbar_n(nz+1)-Z_n(   nz+1))*zinv2 *slope_tapered(3,nz+1,n)**2*Ki(nz+1,n)
        Ty1=Ty1*isredi
        
        ! layer dependent coefficients for for solving dT(1)/dt+d/dz*K_33*d/dz*T(1) = ...
        a(nz)=0.0_WP
        !!PS c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv * (area(nz+1,n)/areasvol(nz,n))
        c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv * area(nz+1,n)/areasvol(nz,n) 
        b(nz)=-c(nz)+hnode_new(nz,n)      ! ale
        ! update from the vertical advection --> comes from splitting of vert 
        ! velocity into explicite and implicite contribution
        if (do_wimpl) then
            !___________________________________________________________________
            ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for 
            ! numerical reasons, to gurante that area/areasvol in case of no 
            ! cavity is ==1.0_WP
            v_adv =zinv * ( area(nz  ,n)/areasvol(nz,n) )
            b(nz) =b(nz)+Wvel_i(nz, n)*v_adv
            
            !!PS v_adv =zinv * ( area(nz+1,n)/areasvol(nz,n) )
            v_adv =zinv * area(nz+1,n)/areasvol(nz,n)
            b(nz) =b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
            c(nz) =c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
        end if        
        ! backup zinv2 for next depth level
        zinv1=zinv2
        !_______________________________________________________________________
        ! Regular part of coefficients: --> 2nd...nl-2 layer
        do nz=nzmin+1, nzmax-2
        
            ! 1/dz(nz)
            zinv2=1.0_WP/(Z_n(nz)-Z_n(nz+1))
            ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
            Ty = (Z_n(nz-1   )-zbar_n(nz  ))*zinv1 *slope_tapered(3,nz-1,n)**2*Ki(nz-1,n)+ &
                 (zbar_n(nz  )-Z_n(nz     ))*zinv1 *slope_tapered(3,nz  ,n)**2*Ki(nz  ,n)
            Ty1= (Z_n(nz     )-zbar_n(nz+1))*zinv2 *slope_tapered(3,nz  ,n)**2*Ki(nz  ,n)+ &
                 (zbar_n(nz+1)-Z_n(nz+1   ))*zinv2 *slope_tapered(3,nz+1,n)**2*Ki(nz+1,n)
            Ty =Ty *isredi
            Ty1=Ty1*isredi
            
            ! layer dependent coefficients for for solving dT(nz)/dt+d/dz*K_33*d/dz*T(nz) = ...
            !___________________________________________________________________
            ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for 
            ! numerical reasons, to gurante that area/areasvol in case of no 
            ! cavity is ==1.0_WP   
            a(nz)=-(Kv(nz,n)  +Ty )*zinv1*zinv * ( area(nz  ,n)/areasvol(nz,n) ) 
            !!PS c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv * ( area(nz+1,n)/areasvol(nz,n) )
            c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv * area(nz+1,n)/areasvol(nz,n)
            b(nz)=-a(nz)-c(nz)+hnode_new(nz,n)
            
            ! backup zinv2 for next depth level
            zinv1=zinv2
            ! update from the vertical advection
            if (do_wimpl) then
                !_______________________________________________________________
                ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for 
                ! numerical reasons, to gurante that area/areasvol in case of no 
                ! cavity is ==1.0_WP   
                v_adv=zinv * ( area(nz  ,n)/areasvol(nz,n) )
                a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv
                b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv
                !!PS v_adv=v_adv*areasvol(nz+1,n)/areasvol(nz,n)
                !!PS v_adv=zinv * ( area(nz+1,n)/areasvol(nz,n) )
                v_adv=zinv * area(nz+1,n)/areasvol(nz,n) 
                b(nz)=b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
                c(nz)=c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
            end if
        end do ! --> do nz=nzmin+1, nzmax-2
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> nl-1 layer
        nz=nzmax-1 
        
        zinv=1.0_WP*dt   ! no ... /(zbar(nzmax-1)-zbar(nzmax)) because of ale
        
        ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
        Ty= (Z_n(nz-1) -zbar_n(nz)) * zinv1 * slope_tapered(3,nz-1,n)**2 * Ki(nz-1,n) + &
            (zbar_n(nz)-Z_n(nz)   ) * zinv1 * slope_tapered(3,nz  ,n)**2 * Ki(nz,n)
        Ty =Ty *isredi
        ! layer dependent coefficients for for solving dT(nz)/dt+d/dz*K_33*d/dz*T(nz) = ...
        
        !___________________________________________________________________
        ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for 
        ! numerical reasons, to gurante that area/areasvol in case of no 
        ! cavity is ==1.0_WP
        a(nz)=-(Kv(nz,n)+Ty)*zinv1*zinv* ( area(nz  ,n)/areasvol(nz,n) )
        c(nz)=0.0_WP
        b(nz)=-a(nz)+hnode_new(nz,n)
        
        ! update from the vertical advection
        if (do_wimpl) then
            !___________________________________________________________________
            ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for 
            ! numerical reasons, to gurante that area/areasvol in case of no 
            ! cavity is ==1.0_WP
            v_adv=zinv* ( area(nz  ,n)/areasvol(nz,n) )
            a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv       
            b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv
        end if
        
        !_______________________________________________________________________
        ! the rhs (inhomogene part): --> rhs = K_33*dt*d/dz*Tstar --> Tstar...trarr
        ! solve difference quotient for rhs --> tr
        !  RHS at Volume_2:
        !  
        !  RHS*V_2 = K_33*dt*(T_1-T_2)/(Z_1-Z_2)*V_2 - K_33*dt*(T_2-T_3)/(Z_2-Z_3)*V_3
        !          = -a*T_1 + (a+c)*T_2 - c*T_3
        !
        ! -+--> tr(1) =(a(1)+c(1))*trarr(1,n)-c(1)*trarr(2,n)
        !  |--> a(1)=0
        nz=nzmin
        dz=hnode_new(nz,n) ! It would be (zbar(nz)-zbar(nz+1)) if not ALE
        tr(nz)=-(b(nz)-dz)*trarr(nz,n)-c(nz)*trarr(nz+1,n)
        
        do nz=nzmin+1,nzmax-2
            dz=hnode_new(nz,n)
            tr(nz)= -a(nz)     * trarr(nz-1,n) &
                    -(b(nz)-dz)* trarr(nz,n) &
                    -c(nz)     * trarr(nz+1,n)
        end do
        
        nz=nzmax-1
        dz=hnode_new(nz,n)
        tr(nz)=-a(nz)*trarr(nz-1,n)-(b(nz)-dz)*trarr(nz,n)       
        !_______________________________________________________________________
        ! Add KPP nonlocal fluxes to the rhs (only T and S currently)
        ! use here blmc or kpp_oblmixc instead of Kv, since Kv already contains 
        ! at this point the mixing enhancments from momix, instable
        ! mixing or windmixing which are to much for nonlocal 
        ! transports and lead to instability of the model
        if (use_kpp_nonlclflx) then
            if (tracers%data(tr_num)%ID==2) then 
                rsss=ref_sss
                if (ref_sss_local) rsss=tracers%data(tr_num)%values(1,n)
            end if
            
            !___________________________________________________________________
            ! use fesom1.4 KPP
            if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then
                if     (tracers%data(tr_num)%ID==1) then ! temperature
                    ! --> no fluxes to the top out of the surface, no fluxes 
                    !     downwards out of the bottom
                    !___surface_________________________________________________
                    nz = nzmin
                    tr(nz)=tr(nz) &
                               +(-MIN(ghats(nz+1,n)*blmc(nz+1,n,2), 1.0_WP)*(area(nz+1,n)/areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                    !___bulk____________________________________________________
                    do nz=nzmin+1, nzmax-2
                        tr(nz)=tr(nz) &
                               +( MIN(ghats(nz  ,n)*blmc(nz  ,n,2), 1.0_WP)*(area(nz  ,n)/areasvol(nz,n)) &
                                 -MIN(ghats(nz+1,n)*blmc(nz+1,n,2), 1.0_WP)*(area(nz+1,n)/areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                    end do
                    !___bottom__________________________________________________
                    nz = nzmax-1
                    tr(nz)=tr(nz) &
                               +( MIN(ghats(nz  ,n)*blmc(nz  ,n,2), 1.0_WP)*(area(nz  ,n)/areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                                
                elseif (tracers%data(tr_num)%ID==2) then ! salinity
                    ! --> no fluxes to the top out of the surface, no fluxes 
                    !     downwards out of the bottom
                    !___surface_________________________________________________
                    nz = nzmin
                    tr(nz)=tr(nz) &
                               -(-MIN(ghats(nz+1,n)*blmc(nz+1,n,3), 1.0_WP)*(area(nz+1,n)/areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                    !___bulk____________________________________________________            
                    do nz=nzmin+1, nzmax-2
                        tr(nz)=tr(nz) &
                               -( MIN(ghats(nz  ,n)*blmc(nz  ,n,3), 1.0_WP)*(area(nz  ,n)/areasvol(nz,n)) &
                                 -MIN(ghats(nz+1,n)*blmc(nz+1,n,3), 1.0_WP)*(area(nz+1,n)/areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                    end do
                    !___bottom__________________________________________________
                    nz = nzmax-1
                    tr(nz)=tr(nz) &
                               -( MIN(ghats(nz  ,n)*blmc(nz  ,n,3), 1.0_WP)*(area(nz  ,n)/areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                end if 
            !___________________________________________________________________
            ! use cvmix KPP
            elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then
                if     (tracers%data(tr_num)%ID==1) then ! temperature
                    !___surface_________________________________________________
                    nz = nzmin
                    tr(nz)=tr(nz) &
                               +(-MIN(kpp_nonlcltranspT(nz+1,n)*kpp_oblmixc(nz+1,n,2), 1.0_WP)*(area(nz+1,n)/areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                    !___bulk____________________________________________________
                    do nz=nzmin+1, nzmax-2
                        tr(nz)=tr(nz) &
                               +( MIN(kpp_nonlcltranspT(nz  ,n)*kpp_oblmixc(nz  ,n,2), 1.0_WP)*(area(nz  ,n)/areasvol(nz,n)) &
                                 -MIN(kpp_nonlcltranspT(nz+1,n)*kpp_oblmixc(nz+1,n,2), 1.0_WP)*(area(nz+1,n)/areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                    end do
                    !___bottom__________________________________________________
                    nz = nzmax-1
                    tr(nz)=tr(nz) &
                               +( MIN(kpp_nonlcltranspT(nz  ,n)*kpp_oblmixc(nz  ,n,2), 1.0_WP)*(area(nz  ,n)/areasvol(nz,n)) &
                                ) * heat_flux(n) / vcpw * dt
                                
                elseif (tracers%data(tr_num)%ID==2) then ! salinity
                    !___surface_________________________________________________
                    nz = nzmin
                    tr(nz)=tr(nz) &
                               -(-MIN(kpp_nonlcltranspS(nz+1,n)*kpp_oblmixc(nz+1,n,3), 1.0_WP)*(area(nz+1,n)/areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                    !___bulk____________________________________________________            
                    do nz=nzmin+1, nzmax-2
                        tr(nz)=tr(nz) &
                               -( MIN(kpp_nonlcltranspS(nz  ,n)*kpp_oblmixc(nz  ,n,3), 1.0_WP)*(area(nz  ,n)/areasvol(nz,n)) &
                                 -MIN(kpp_nonlcltranspS(nz+1,n)*kpp_oblmixc(nz+1,n,3), 1.0_WP)*(area(nz+1,n)/areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                    end do
                    !___bottom__________________________________________________
                    nz = nzmax-1
                    tr(nz)=tr(nz) &
                               -( MIN(kpp_nonlcltranspS(nz  ,n)*kpp_oblmixc(nz  ,n,3), 1.0_WP)*(area(nz  ,n)/areasvol(nz,n)) &
                                ) * rsss * water_flux(n) * dt
                end if
            end if
        end if ! --> if (use_kpp_nonlclflx) then
        
        !_______________________________________________________________________
        ! case of activated shortwave penetration into the ocean, ad 3d contribution
        if (use_sw_pene .and. tracers%data(tr_num)%ID==1) then
            do nz=nzmin, nzmax-1
                zinv=1.0_WP*dt  !/(zbar(nz)-zbar(nz+1)) ale!
                !!PS tr(nz)=tr(nz)+(sw_3d(nz, n)-sw_3d(nz+1, n) * ( area(nz+1,n)/areasvol(nz,n)) ) * zinv
                tr(nz)=tr(nz)+(sw_3d(nz, n)-sw_3d(nz+1, n) * area(nz+1,n)/areasvol(nz,n)) * zinv
            end do
        end if
        
        !_______________________________________________________________________
        !  The first row contains also the boundary condition from heatflux, 
        !  freshwaterflux and relaxation terms
        !  --> trarr(1,n)*water_flux(n) : latent heatflux contribution due to 
        !      cell volume. If Volume decreases --> temp has to raise, if volume 
        !      expended --> temp has to decrease
        !                           (-)   ^                        (-)   ^ 
        !                            |    |                         |    | 
        !   IN MOMENT: heat_flux ~~~~|~~~~|~~~~   ,  water_flux ~~~~|~~~~|~~~~
        !  (BUT CHECK!)              |    |                         |    |
        !                            v   (+)                        v   (+) 
        !                            
        tr(nzmin)= tr(nzmin)+bc_surface(n, tracers%data(tr_num)%ID, trarr(nzmin,n), nzmin, partit)  
        
        !_______________________________________________________________________
        ! The forward sweep algorithm to solve the three-diagonal matrix 
        ! problem
        ! 
        !  | b_1 c_1 ...            |   |dTnew_1|
        !  | a_2 b_2 c_2 ...        |   |dTnew_2|
        !  |     a_3 b_3 c_3 ...    | * |dTnew_3| = RHS
        !  |         a_4 b_4 c_4 ...|   |dTnew_3| 
        !  |              :         |   |   :   |
        ! 
        ! 1st: define new coefficents:
        !      --> c'_i = c_i/b_i                               ; i=1
        !          c'_i = c_i/(b_i-a_i*c'_i-1)                  ; i = 2,3,...,n-1
        !      --> rhs'_i = rhs_i/b_i                           ; i=1
        !          rhs'_i = (rhs_i-a_i*d'_i-1)/(b_i-a_i*c'_i-1) ; i = 2,3,...,n-1
        !
        ! 2nd: solution is optained by back substitution
        !      --> dTnew_n = rhs'_n
        !      --> dTnew_i = rhs'_i-c'_i*dTnew_i+1 ; i = n-1,n-2,...,1
        !
        ! initialize c-prime and s,t-prime
        cp(nzmin) = c(nzmin)/b(nzmin)
        tp(nzmin) = tr(nzmin)/b(nzmin)
        
        ! solve for vectors c-prime and t, s-prime
        do nz = nzmin+1,nzmax-1
            m = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
        end do
        
        ! start with back substitution 
        tr(nzmax-1) = tp(nzmax-1)
        
        ! solve for x from the vectors c-prime and d-prime
        do nz = nzmax-2, nzmin, -1
            tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
        end do
        
        !_______________________________________________________________________
        ! update tracer
        ! tr ... dTnew = T^(n+0.5) - T*
        do nz=nzmin,nzmax-1
            ! trarr - before ... T*
            trarr(nz,n)=trarr(nz,n)+tr(nz)
        end do
        
    end do ! --> do n=1,myDim_nod2D
!$OMP END DO
!$OMP END PARALLEL
end subroutine diff_ver_part_impl_ale
!
!
!===============================================================================
subroutine diff_ver_part_redi_expl(tracers, partit, mesh) 
    use o_ARRAYS
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_TRACER
    USE o_param
    use g_config
    use g_comm_auto
    IMPLICIT NONE
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    !___________________________________________________________________________
    integer                  :: n, k, elem, nz
    integer                  :: n2, nl1, ul1, nl2
    real(kind=WP)            :: Tx, Ty, vd_flux(mesh%nl)
    real(kind=WP)            :: tr_xynodes(2,mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP)            :: zbar_n(mesh%nl), z_n(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), pointer   :: del_ttf(:,:)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    del_ttf => tracers%work%del_ttf

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, k, elem, nz, n2, nl1, ul1, nl2, Tx, Ty, vd_flux, zbar_n, z_n)
!$OMP DO
    !___________________________________________________________________________
    do n=1, myDim_nod2D
        nl1=nlevels_nod2D(n)-1
        ul1=ulevels_nod2D(n)
        !!PS do nz=1, nl1
        do nz=ul1, nl1
            Tx=0.0_WP
            Ty=0.0_WP
            do k=1, nod_in_elem2D_num(n)
                elem=nod_in_elem2D(k,n)
                !!PS if(nz.LE.(nlevels(elem)-1)) then
                if( nz.LE.(nlevels(elem)-1) .and. nz.GE.(ulevels(elem))) then
                    Tx=Tx+tr_xy(1,nz,elem)*elem_area(elem)
                    Ty=Ty+tr_xy(2,nz,elem)*elem_area(elem)
                endif
            end do
            tr_xynodes(1,nz,n)=tx/3.0_WP/areasvol(nz,n)
            tr_xynodes(2,nz,n)=ty/3.0_WP/areasvol(nz,n)
        end do
    end do
    
!$OMP END DO
    ! no halo exchange of tr_xynodes is needed !
!$OMP DO
    do n=1, myDim_nod2D
        nl1=nlevels_nod2D(n)-1
        ul1=ulevels_nod2D(n)
        vd_flux=0._WP
        
        !_______________________________________________________________________
        zbar_n(1:mesh%nl  )=0.0_WP
        z_n   (1:mesh%nl-1)=0.0_WP
        zbar_n(nl1+1)=zbar_n_bot(n)
        z_n(nl1)=zbar_n(nl1+1) + hnode_new(nl1,n)/2.0_WP
        do nz=nl1, ul1+1, -1
            zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
            z_n(nz-1)  = zbar_n(nz)   + hnode_new(nz-1,n)/2.0_WP
        end do
        zbar_n(ul1) = zbar_n(ul1+1)   + hnode_new(ul1,n)
        
        !_______________________________________________________________________
        do nz=ul1+1,nl1
            vd_flux(nz)=(z_n(nz-1)-zbar_n(nz))*(slope_tapered(1,nz-1,n)*tr_xynodes(1,nz-1,n)+slope_tapered(2,nz-1,n)*tr_xynodes(2,nz-1,n))*Ki(nz-1,n)
            vd_flux(nz)=vd_flux(nz)+&
                        (zbar_n(nz)-z_n(nz))  *(slope_tapered(1,nz,n)  *tr_xynodes(1,nz,n)  +slope_tapered(2,nz,n)  *tr_xynodes(2,nz,n))  *Ki(nz,n)
            vd_flux(nz)=vd_flux(nz)/(z_n(nz-1)-z_n(nz))*area(nz,n)
        enddo
        do nz=ul1,nl1
            del_ttf(nz,n) = del_ttf(nz,n)+(vd_flux(nz) - vd_flux(nz+1))*dt/areasvol(nz,n)
        enddo
    end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine diff_ver_part_redi_expl
!
!
!===============================================================================
subroutine diff_part_hor_redi(tracers, partit, mesh) 
    use o_ARRAYS
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_TRACER
    use o_param
    use g_config
    IMPLICIT NONE
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: edge
    real(kind=WP)            :: deltaX1, deltaY1, deltaX2, deltaY2
    integer                  :: nl1, ul1, nl2, ul2, nl12, ul12, nz, el(2), elnodes(3), enodes(2)
    real(kind=WP)            :: c, Fx, Fy, Tx, Ty, Tx_z, Ty_z, SxTz, SyTz, Tz(2)
    real(kind=WP)            :: rhs1(mesh%nl-1), rhs2(mesh%nl-1), Kh, dz
    real(kind=WP)            :: isredi=0._WP
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), pointer   :: del_ttf(:,:)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    del_ttf => tracers%work%del_ttf

    !___________________________________________________________________________
    if (Redi) isredi=1._WP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, deltaX1, deltaY1, deltaX2, deltaY2, &
!$OMP                   nl1, ul1, nl2, ul2, nl12, ul12, nz, el, elnodes, enodes, &
!$OMP                             c, Fx, Fy, Tx, Ty, Tx_z, Ty_z, SxTz, SyTz, Tz, &
!$OMP                                                          rhs1, rhs2, Kh, dz)
!$OMP DO
    do edge=1, myDim_edge2D
        rhs1=0.0_WP
        rhs2=0.0_WP
        !_______________________________________________________________________
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        el=edge_tri(:,edge)
        enodes=edges(:,edge)
        nl1=nlevels(el(1))-1
        ul1=ulevels(el(1))
        elnodes=elem2d_nodes(:,el(1))
        !_______________________________________________________________________
        nl2=0
        ul2=0
        if (el(2)>0) then 
            nl2=nlevels(el(2))-1
            ul2=ulevels(el(2))
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
        endif
        !_______________________________________________________________________
        nl12=min(nl1,nl2)
        ul12=max(ul1,ul2)
        !_______________________________________________________________________
        ! (A)
        do nz=ul1, ul12-1
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=helem(nz, el(1))
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=tr_xy(1,nz,el(1))
            Ty=tr_xy(2,nz,el(1))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=(-deltaX1*Fy+deltaY1*Fx)*dz
            rhs1(nz) = rhs1(nz) + c
            rhs2(nz) = rhs2(nz) - c
        end do
        !_______________________________________________________________________
        ! (B)
        if (ul2>0) then
            do nz=ul2,ul12-1
                Kh=sum(Ki(nz, enodes))/2.0_WP
                dz=helem(nz, el(2))
                Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
                SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
                SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
                Tx=tr_xy(1,nz,el(2))
                Ty=tr_xy(2,nz,el(2))
                Fx=Kh*(Tx+SxTz*isredi)
                Fy=Kh*(Ty+SyTz*isredi)
                c=(deltaX2*Fy-deltaY2*Fx)*dz
                rhs1(nz) = rhs1(nz) + c
                rhs2(nz) = rhs2(nz) - c
            end do
        end if
        !_______________________________________________________________________
        ! (C)
        do nz=ul12, nl12
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=sum(helem(nz, el))/2.0_WP
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=0.5_WP*(tr_xy(1,nz,el(1))+tr_xy(1,nz,el(2)))
            Ty=0.5_WP*(tr_xy(2,nz,el(1))+tr_xy(2,nz,el(2)))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=((deltaX2-deltaX1)*Fy-(deltaY2-deltaY1)*Fx)*dz
            rhs1(nz) = rhs1(nz) + c
            rhs2(nz) = rhs2(nz) - c
        enddo
        !_______________________________________________________________________
        ! (D)
        do nz=nl12+1, nl1
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=helem(nz, el(1))
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=tr_xy(1,nz,el(1))
            Ty=tr_xy(2,nz,el(1))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=(-deltaX1*Fy+deltaY1*Fx)*dz
            rhs1(nz) = rhs1(nz) + c
            rhs2(nz) = rhs2(nz) - c
        end do
        !_______________________________________________________________________
        ! (E)
        do nz=nl12+1, nl2
            Kh=sum(Ki(nz, enodes))/2.0_WP
            dz=helem(nz, el(2))
            Tz=0.5_WP*(tr_z(nz,enodes)+tr_z(nz+1,enodes))
            SxTz=sum(Tz*slope_tapered(1,nz,enodes))/2.0_WP
            SyTz=sum(Tz*slope_tapered(2,nz,enodes))/2.0_WP
            Tx=tr_xy(1,nz,el(2))
            Ty=tr_xy(2,nz,el(2))
            Fx=Kh*(Tx+SxTz*isredi)
            Fy=Kh*(Ty+SyTz*isredi)
            c=(deltaX2*Fy-deltaY2*Fx)*dz
            rhs1(nz) = rhs1(nz) + c
            rhs2(nz) = rhs2(nz) - c
        end do
        !_______________________________________________________________________
        nl12=max(nl1,nl2)
        ul12 = ul1
        if (ul2>0) ul12=min(ul1,ul2)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_set_lock(partit%plock(enodes(1)))
#else
!$OMP ORDERED
#endif
        del_ttf(ul12:nl12,enodes(1))=del_ttf(ul12:nl12,enodes(1))+rhs1(ul12:nl12)*dt/areasvol(ul12:nl12,enodes(1))
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(enodes(1)))
        call omp_set_lock  (partit%plock(enodes(2)))
#endif
        del_ttf(ul12:nl12,enodes(2))=del_ttf(ul12:nl12,enodes(2))+rhs2(ul12:nl12)*dt/areasvol(ul12:nl12,enodes(2))
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(enodes(2)))
#else
!$OMP END ORDERED
#endif
    end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine diff_part_hor_redi
!
!
!===============================================================================
SUBROUTINE diff_part_bh(tr_num, dynamics, tracers, partit, mesh) 
    use o_ARRAYS, only: 
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_TRACER
    use MOD_DYN
    use o_param
    use g_config
    use g_comm_auto
    IMPLICIT NONE
    integer,        intent(in),    target    :: tr_num
    type(t_dyn)   , intent(inout), target    :: dynamics
    type(t_tracer), intent(inout), target    :: tracers
    type(t_mesh)  , intent(in)   , target    :: mesh
    type(t_partit), intent(inout), target    :: partit
    integer                                  :: n, nz, ed, el(2), en(2), k, elem, nzmin, nzmax
    integer                                  :: elnodes1(3), elnodes2(3)
    real(kind=WP)                            :: u1, v1, len, vi, ww, tt(mesh%nl-1)
    real(kind=WP), pointer                   :: temporary_ttf(:,:)
    real(kind=WP), pointer                   :: UV(:,:,:)
    real(kind=WP), pointer                   :: ttf(:,:)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 


    UV            => dynamics%uv(:,:,:)
    ttf           => tracers%data(tr_num)%values
    temporary_ttf => tracers%work%del_ttf !use already allocated working array. could be fct_LO instead etc.

!$OMP PARALLEL DO
        do n=1, myDim_nod2D+eDim_nod2D
           temporary_ttf(:, n)=0.0_8
        end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, ed, el, en, k, elem, nzmin, nzmax, u1, v1, len, vi, tt, ww, &
!$OMP elnodes1, elnodes2)
!$OMP DO
    DO ed=1, myDim_edge2D!+eDim_edge2D
       if (myList_edge2D(ed) > edge2D_in) cycle
       el=edge_tri(:,ed)
       en=edges(:,ed)
       len=sqrt(sum(elem_area(el)))
       nzmax = minval(nlevels(el))
       nzmin = maxval(ulevels(el))
       elnodes1=elem2d_nodes(:,el(1))
       elnodes2=elem2d_nodes(:,el(2))
       DO  nz=nzmin, nzmax-1
           u1=maxval(ttf(nz, elnodes1))-minval(ttf(nz, elnodes2))
           v1=minval(ttf(nz, elnodes1))-maxval(ttf(nz, elnodes2))
           vi=u1*u1+v1*v1
           tt(nz)=ttf(nz,en(1))-ttf(nz,en(2))
           vi=sqrt(max(tracers%data(tr_num)%gamma0_tra,            &
                   max(tracers%data(tr_num)%gamma1_tra*sqrt(vi),   & 
                       tracers%data(tr_num)%gamma2_tra*     vi)    &
                                                         )*len)       
           tt(nz)=tt(nz)*vi
       END DO
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_set_lock  (partit%plock(en(1)))
#else
!$OMP ORDERED
#endif
       temporary_ttf(nzmin:nzmax-1,en(1))=temporary_ttf(nzmin:nzmax-1,en(1))-tt(nzmin:nzmax-1)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_unset_lock(partit%plock(en(1)))
       call omp_set_lock  (partit%plock(en(2)))
#endif
       temporary_ttf(nzmin:nzmax-1,en(2))=temporary_ttf(nzmin:nzmax-1,en(2))+tt(nzmin:nzmax-1)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_unset_lock(partit%plock(en(2)))
#else
!$OMP END ORDERED
#endif
    END DO
!$OMP END DO
!$OMP MASTER
    call exchange_nod(temporary_ttf, partit)
!$OMP END MASTER
!$OMP BARRIER
    ! ===========
    ! Second round: 
    ! ===========
!$OMP DO
    DO ed=1, myDim_edge2D!+eDim_edge2D
       if (myList_edge2D(ed)>edge2D_in) cycle
          el=edge_tri(:,ed)
          en=edges(:,ed)
          len=sqrt(sum(elem_area(el)))
          nzmax = minval(nlevels(el))
          nzmin = maxval(ulevels(el))
          elnodes1=elem2d_nodes(:,el(1))
          elnodes2=elem2d_nodes(:,el(2))
          DO  nz=nzmin, nzmax-1
              u1=maxval(ttf(nz, elnodes1))-minval(ttf(nz, elnodes2))
              v1=minval(ttf(nz, elnodes1))-maxval(ttf(nz, elnodes2))
              vi=u1*u1+v1*v1
              tt(nz)=temporary_ttf(nz,en(1))-temporary_ttf(nz,en(2))
              vi=sqrt(max(tracers%data(tr_num)%gamma0_tra,            &
                      max(tracers%data(tr_num)%gamma1_tra*sqrt(vi),   & 
                          tracers%data(tr_num)%gamma2_tra*     vi)    &
                                                            )*len)                    
              tt(nz)=-tt(nz)*vi*dt
          END DO 
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
          call omp_set_lock  (partit%plock(en(1)))
#else
!$OMP ORDERED
#endif
          ttf(nzmin:nzmax-1,en(1))=ttf(nzmin:nzmax-1,en(1))-tt(nzmin:nzmax-1)/area(nzmin:nzmax-1,en(1))
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
          call omp_unset_lock(partit%plock(en(1)))
          call omp_set_lock  (partit%plock(en(2)))
#endif
          ttf(nzmin:nzmax-1,en(2))=ttf(nzmin:nzmax-1,en(2))+tt(nzmin:nzmax-1)/area(nzmin:nzmax-1,en(2))
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
          call omp_unset_lock(partit%plock(en(2)))
#else
!$OMP END ORDERED
#endif
    END DO
!$OMP END DO
!$OMP END PARALLEL
call exchange_nod(ttf, partit)
!$OMP BARRIER
end subroutine diff_part_bh
!
!
!===============================================================================
! this function returns a boundary conditions for a specified tracer ID and surface node
! ID = 0 and 1 are reserved for temperature and salinity
FUNCTION bc_surface(n, id, sval, nzmin, partit) 
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE o_ARRAYS
  USE g_forcing_arrays
  USE g_config
  implicit none
  
  integer,       intent(in)            :: n, id, nzmin 
  real(kind=WP), intent(in)            :: sval
  type(t_partit),intent(inout), target :: partit
  REAL(kind=WP)                        :: bc_surface
  character(len=10)                    :: id_string

  !  --> is_nonlinfs=1.0 for zelvel,zstar ....                            
  !  --> is_nonlinfs=0.0 for linfs
  SELECT CASE (id)
    CASE (1)
        bc_surface=-dt*(heat_flux(n)/vcpw + sval*water_flux(n)*is_nonlinfs)    
    CASE (2)
        ! --> real_salt_flux(:): salt flux due to containment/releasing of salt
        !     by forming/melting of sea ice
        bc_surface= dt*(virtual_salt(n) & !--> is zeros for zlevel/zstar
                    + relax_salt(n) - real_salt_flux(n)*is_nonlinfs)
    CASE (101) ! apply boundary conditions to tracer ID=101
        bc_surface= dt*(prec_rain(n))! - real_salt_flux(n)*is_nonlinfs)
    CASE (301)
        bc_surface=0.0_WP
    CASE (302)
        bc_surface=0.0_WP
    CASE (303)
        bc_surface=0.0_WP
    CASE (304)
        bc_surface= dt*(hosing_flux(n))
    CASE (501) ! ice-shelf water due to basal melting
        if (nzmin==1) then
           bc_surface = 0.0_WP
        else
           bc_surface = dt*water_flux(n)*(sval-1.0)
        end if
    CASE DEFAULT
      if (partit%mype==0) then
         write (id_string, "(I3)") id
         if (partit%mype==0) write(*,*) 'invalid ID '//trim(id_string)//' specified in boundary conditions'
         if (partit%mype==0) write(*,*) 'the model will stop!'
      end if
      call par_ex(partit%MPI_COMM_FESOM, partit%mype)
      stop
  END SELECT
  RETURN
END FUNCTION
