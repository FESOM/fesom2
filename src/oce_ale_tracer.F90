module diff_part_hor_redi_interface
  interface
    subroutine diff_part_hor_redi(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module adv_tracers_muscle_ale_interface
  interface
    subroutine adv_tracers_muscle_ale(ttfAB, num_ord, do_Xmoment, mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh), intent(in) , target :: mesh
      integer                  :: do_Xmoment
      real(kind=WP)            :: ttfAB(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP)            :: num_ord
    end subroutine
  end interface
end module
module adv_tracers_vert_ppm_ale_interface
  interface
    subroutine adv_tracers_vert_ppm_ale(ttf, do_Xmoment, mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh), intent(in) , target :: mesh
      integer                  :: do_Xmoment
      real(kind=WP)            :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    end subroutine
  end interface
end module
module adv_tracers_ale_interface
  interface
    subroutine adv_tracers_ale(tr_num, mesh)
      use mod_mesh
      integer :: tr_num
      type(t_mesh), intent(in) , target :: mesh    
    end subroutine
  end interface
end module
module diff_ver_part_expl_ale_interface
  interface
    subroutine diff_ver_part_expl_ale(tr_num, mesh)
      use MOD_MESH
      type(t_mesh), intent(in) , target :: mesh    
      integer                  :: tr_num
    end subroutine
  end interface
end module
module diff_ver_part_redi_expl_interface
  interface
    subroutine diff_ver_part_redi_expl(mesh)
      use MOD_MESH
      type(t_mesh), intent(in) , target :: mesh    
    end subroutine
  end interface
end module
module diff_ver_part_impl_ale_interface
  interface
    subroutine diff_ver_part_impl_ale(tr_num, mesh)
      use MOD_MESH
      type(t_mesh), intent(in) , target :: mesh
      integer                  :: tr_num
    end subroutine
  end interface
end module
module diff_tracers_ale_interface
  interface
    subroutine diff_tracers_ale(tr_num, mesh)
      use mod_mesh
      integer, intent(in)      :: tr_num
      type(t_mesh), intent(in) , target :: mesh
    end subroutine
  end interface
end module
module bc_surface_interface
  interface
    function bc_surface(n, id, mesh)
      use mod_mesh
      integer , intent(in)              :: n, id
      type(t_mesh), intent(in) , target :: mesh
      real(kind=WP)                     :: bc_surface
    end function
  end interface
end module
module diff_part_bh_interface
  interface
    subroutine diff_part_bh(ttf, mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh) , intent(in),    target :: mesh
      real(kind=WP), intent(inout), target :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    end subroutine
  end interface
end module
module diff_ver_recom_expl_interface
  interface
    subroutine diff_ver_recom_expl(tr_num, mesh)
      use mod_mesh
      use g_PARSUP
      integer, intent(in)      :: tr_num
      type(t_mesh), intent(in) , target :: mesh
    end subroutine
  end interface
end module
module ver_sinking_recom_benthos_interface
  interface
    subroutine ver_sinking_recom_benthos(tr_num, mesh)
      use mod_mesh
      use g_PARSUP
      type(t_mesh), intent(in) , target :: mesh
      integer, intent(in)      :: tr_num
!      real(kind=WP), intent(inout), target :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    end subroutine
  end interface
end module
module integrate_bottom_interface
  interface
    subroutine integrate_bottom(tflux, mesh)
      use mod_mesh
      use g_PARSUP
      type(t_mesh), intent(in) , target :: mesh
      integer, intent(inout)      :: tflux
    end subroutine
  end interface
end module
!
!
!===============================================================================
! Driving routine    Here with ALE changes!!!
subroutine solve_tracers_ale(mesh)
    use g_config
    use g_parsup
    use o_PARAM, only: num_tracers, SPP, Fer_GM
    use o_arrays
    use mod_mesh
    use g_comm_auto
    use o_tracers
    use Toy_Channel_Soufflet
    use adv_tracers_ale_interface
    use diff_tracers_ale_interface
    use REcom_config, only: ciso     ! to calculation radioactive decay of 14C
    use REcoM_ciso, only: lambda_14  ! decay constant of 14C
    
    implicit none
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: tr_num, node, nzmax, nzmin
    real(kind=WP)            :: aux_tr(mesh%nl-1,myDim_nod2D+eDim_nod2D)

#include "associate_mesh.h"
    !___________________________________________________________________________
    if (SPP) call cal_rejected_salt(mesh)
    if (SPP) call app_rejected_salt(mesh)
    
    !___________________________________________________________________________
    ! update 3D velocities with the bolus velocities:
    ! 1. bolus velocities are computed according to GM implementation after R. Ferrari et al., 2010
    ! 2. bolus velocities are used only for advecting tracers and shall be subtracted back afterwards
    if (Fer_GM) then
        UV    =UV    +fer_UV
        Wvel_e=Wvel_e+fer_Wvel
        Wvel  =Wvel  +fer_Wvel
    end if
    !___________________________________________________________________________
    ! loop over all tracers 
    do tr_num=1,num_tracers
        ! do tracer AB (Adams-Bashfort) interpolation only for advectiv part 
        ! needed
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call init_tracers_AB'//achar(27)//'[0m'
        call init_tracers_AB(tr_num, mesh)
        
        ! advect tracers
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call adv_tracers_ale'//achar(27)//'[0m'
        call adv_tracers_ale(tr_num, mesh)
        
        ! diffuse tracers 
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call diff_tracers_ale'//achar(27)//'[0m'
        call diff_tracers_ale(tr_num, mesh)
        
        ! relax to salt and temp climatology
        if (flag_debug .and. mype==0)  print *, achar(27)//'[37m'//'         --> call relax_to_clim'//achar(27)//'[0m'
        call relax_to_clim(tr_num, mesh)
        if ((toy_ocean) .AND. (TRIM(which_toy)=="soufflet")) call relax_zonal_temp(mesh)
        call exchange_nod(tr_arr(:,:,tr_num))
    end do
    
    !___________________________________________________________________________
    do tr_num=1, ptracers_restore_total           
        tr_arr(:,ptracers_restore(tr_num)%ind2,ptracers_restore(tr_num)%locid)=1.0_WP
    end do
    
    !___________________________________________________________________________
    ! subtract the the bolus velocities back from 3D velocities:
    if (Fer_GM) then
        UV    =UV    -fer_UV
        Wvel_e=Wvel_e-fer_Wvel
        Wvel  =Wvel  -fer_Wvel
    end if
    
    !___________________________________________________________________________
    ! to avoid crash with high salinities when coupled to atmosphere
    ! --> if we do only where (tr_arr(:,:,2) < 3._WP ) we also fill up the bottom 
    !     topogrpahy with values which are then writte into the output --> thats why
    !     do node=1,.... and tr_arr(node,1:nzmax,2)
    do node=1,myDim_nod2D+eDim_nod2D
        nzmax=nlevels_nod2D(node)-1
        nzmin=ulevels_nod2D(node)
        !!PS where (tr_arr(1:nzmax,node,2) > 45._WP)
        !!PS     tr_arr(1:nzmax,node,2)=45._WP
        !!PS end where
        where (tr_arr(nzmin:nzmax,node,2) > 45._WP)
            tr_arr(nzmin:nzmax,node,2)=45._WP
        end where

        !!PS where (tr_arr(1:nzmax,node,2) < 3._WP )
        !!PS     tr_arr(1:nzmax,node,2)=3._WP
        !!PS end where
        where (tr_arr(nzmin:nzmax,node,2) < 3._WP )
            tr_arr(nzmin:nzmax,node,2)=3._WP
        end where
        
!!PS         if (nzmin>15 .and. mype==0) then 
!!PS             write(*,*) ' tr_arr(:,node,1) = ',tr_arr(:,node,1)
!!PS             write(*,*)
!!PS             write(*,*) ' tr_arr(:,node,2) = ',tr_arr(:,node,2)
!!PS         end if 
    end do
end subroutine solve_tracers_ale
!
!
!===============================================================================
subroutine adv_tracers_ale(tr_num, mesh)
    use g_config, only: flag_debug
    use g_parsup
    use mod_mesh
    use o_arrays
    use diagnostics, only: ldiag_DVD, compute_diag_dvd_2ndmoment_klingbeil_etal_2014, & 
                           compute_diag_dvd_2ndmoment_burchard_etal_2008, compute_diag_dvd
    use adv_tracers_muscle_ale_interface
    use adv_tracers_vert_ppm_ale_interface
    use oce_adv_tra_driver_interfaces
    implicit none
    integer :: tr_num, node, nz
    type(t_mesh), intent(in) , target :: mesh    
    ! del_ttf ... initialised and setted to zero in call init_tracers_AB(tr_num)
    ! --> del_ttf ... equivalent to R_T^n in Danilov etal FESOM2: "from finite element
    !     to finite volume". At the end R_T^n should contain all advection therms and 
    !     the terms due to diffusion.
    ! del_ttf=0d0
    
    !___________________________________________________________________________
    ! if ldiag_DVD=.true. --> compute tracer second moments for the calcualtion 
    ! of discret variance decay
    if (ldiag_DVD .and. tr_num <= 2) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call compute_diag_dvd_2ndmoment'//achar(27)//'[0m'
        call compute_diag_dvd_2ndmoment_klingbeil_etal_2014(tr_num, mesh)
        !!PS call compute_diag_dvd_2ndmoment_burchard_etal_2008(tr_num)
    end if    
    
    !___________________________________________________________________________
    ! horizontal ale tracer advection 
    ! here --> add horizontal advection part to del_ttf(nz,n) = del_ttf(nz,n) + ...
    del_ttf_advhoriz = 0.0_WP
    del_ttf_advvert  = 0.0_WP
    call do_oce_adv_tra(tr_arr(:,:,tr_num), tr_arr_old(:,:,tr_num), UV, wvel, wvel_i, wvel_e, 1, del_ttf_advhoriz, del_ttf_advvert, tra_adv_ph, tra_adv_pv, mesh)   
    !___________________________________________________________________________
    ! update array for total tracer flux del_ttf with the fluxes from horizontal
    ! and vertical advection
    del_ttf=del_ttf+del_ttf_advhoriz+del_ttf_advvert
    
    !___________________________________________________________________________
    ! compute discrete variance decay after Burchard and Rennau 2008
    if (ldiag_DVD .and. tr_num <= 2) then
        if (flag_debug .and. mype==0)  print *, achar(27)//'[38m'//'             --> call compute_diag_dvd'//achar(27)//'[0m'
        call compute_diag_dvd(tr_num, mesh)
    end if     
    
end subroutine adv_tracers_ale
!
!
!===============================================================================
subroutine diff_tracers_ale(tr_num, mesh)
    use mod_mesh
    use g_PARSUP
    use o_arrays
    use o_tracers
    use diff_part_hor_redi_interface
    use diff_ver_part_expl_ale_interface
    use diff_ver_part_redi_expl_interface
    use diff_ver_part_impl_ale_interface
    use diff_part_bh_interface
    use diff_ver_recom_expl_interface
    use ver_sinking_recom_benthos_interface
#if defined(__recom)
    USE REcoM_GloVar
    use recom_config !, recom_debug
use g_support
#endif
    implicit none
    
    integer, intent(in)      :: tr_num
    integer                  :: n, nzmax, nzmin
    type(t_mesh), intent(in) , target :: mesh
real :: net

#include "associate_mesh.h"

#if defined(__recom)
    dtr_bf         = 0.0_WP
    str_bf         = 0.0_WP
    vert_sink      = 0.0_WP
#endif

    !___________________________________________________________________________
    ! convert tr_arr_old(:,:,tr_num)=ttr_n-0.5   --> prepare to calc ttr_n+0.5
    ! eliminate AB (adams bashfort) interpolates tracer, which is only needed for 
    ! tracer advection. For diffusion only need tracer from previouse time step
    tr_arr_old(:,:,tr_num)=tr_arr(:,:,tr_num) !DS: check that this is the right place!
    
    !___________________________________________________________________________
    ! do horizontal diffusiion
    ! write there also horizontal diffusion rhs to del_ttf which is equal the R_T^n 
    ! in danilovs srcipt
    ! includes Redi diffusivity if Redi=.true.
    call diff_part_hor_redi(mesh) ! seems to be ~9% faster than diff_part_hor
    
    !___________________________________________________________________________
    ! do vertical diffusion: explicite 
    if (.not. i_vert_diff) call diff_ver_part_expl_ale(tr_num, mesh)
    ! A projection of horizontal Redi diffussivity onto vertical. This par contains horizontal
    ! derivatives and has to be computed explicitly!
    if (Redi) call diff_ver_part_redi_expl(mesh)

! OG recom bottom boundary layer contribution
#if defined(__recom)

! 1) Remineralization from the benthos
!    Nutrient fluxes come from the bottom boundary
!    Unit [mmol/m2/s]

if (1) then

    if (tracer_id(tr_num) == 1001 .or.    &   ! DIN
        tracer_id(tr_num) == 1002 .or.    &   ! DIC
        tracer_id(tr_num) == 1003 .or.    &   ! Alk
        tracer_id(tr_num) == 1018 .or.    &   ! Si
        tracer_id(tr_num) == 1019 .or.    &   ! Fe

#if defined(__ciso)
        tracer_id(tr_num) == 1033 .or.    &   ! DIC_13
        tracer_id(tr_num) == 1034 .or.    &   ! DIC_14
#endif
        tracer_id(tr_num) == 1022     ) then  ! Oxy

        call diff_ver_recom_expl(tr_num,mesh)

! update tracer fields
        do n=1, myDim_nod2D 
            nzmax=nlevels_nod2D(n)-1
            nzmin=ulevels_nod2D(n)
            tr_arr(nzmin:nzmax,n,tr_num)=tr_arr(nzmin:nzmax,n,tr_num)+ &
                                        dtr_bf(nzmin:nzmax,n)
        end do
    end if

! 2) Sinking in water column 

    if (tracer_id(tr_num) == 1007 .or.    &   ! idetn
        tracer_id(tr_num) == 1008 .or.    &   ! idetc
        tracer_id(tr_num) == 1017 .or.    &   ! idetsi
        tracer_id(tr_num) == 1021 .or.    &   ! idetcal

        tracer_id(tr_num) == 1004 .or.    &   !iphyn
        tracer_id(tr_num) == 1005 .or.    &   !iphyc
        tracer_id(tr_num) == 1020 .or.    &   !iphycal
        tracer_id(tr_num) == 1006 .or.    &   !ipchl

        tracer_id(tr_num) == 1013 .or.    &   !idian
        tracer_id(tr_num) == 1014 .or.    &   !idiac
        tracer_id(tr_num) == 1016 .or.    &   !idiasi
        tracer_id(tr_num) == 1015 .or.    &   !idchl

! if (REcoM_Second_Zoo) then
        tracer_id(tr_num) == 1025 .or.    &   !idetz2n
        tracer_id(tr_num) == 1026 .or.    &   !idetz2c
        tracer_id(tr_num) == 1027 .or.    &   !idetz2si
        tracer_id(tr_num) == 1028 ) then      !idetz2calc
! endif

! sinking
       call recom_sinking_new(tr_num,mesh) !vert_sink

! sinking into the benthos
        call ver_sinking_recom_benthos(tr_num,mesh) !str_bf
       
! update tracer fields
        do n=1, myDim_nod2D 
            nzmax=nlevels_nod2D(n)-1
            nzmin=ulevels_nod2D(n)
tr_arr(nzmin:nzmax,n,tr_num)=tr_arr(nzmin:nzmax,n,tr_num)+ &
                                        vert_sink(nzmin:nzmax,n)
tr_arr(nzmin:nzmax,n,tr_num)=tr_arr(nzmin:nzmax,n,tr_num)+ &
                                        str_bf(nzmin:nzmax,n)
        end do                             
    end if
end if

#endif
    
    !___________________________________________________________________________
    ! Update tracers --> calculate T* see Danilov etal "FESOM2 from finite elements
    ! to finite volume" 
    ! T* =  (dt*R_T^n + h^(n-0.5)*T^(n-0.5))/h^(n+0.5)
    do n=1, myDim_nod2D 
        nzmax=nlevels_nod2D(n)-1
        nzmin=ulevels_nod2D(n)
        !!PS del_ttf(1:nzmax,n)=del_ttf(1:nzmax,n)+tr_arr(1:nzmax,n,tr_num)* &
        !!PS                             (hnode(1:nzmax,n)-hnode_new(1:nzmax,n))
        !!PS tr_arr(1:nzmax,n,tr_num)=tr_arr(1:nzmax,n,tr_num)+ &
        !!PS                             del_ttf(1:nzmax,n)/hnode_new(1:nzmax,n)
        
        del_ttf(nzmin:nzmax,n)=del_ttf(nzmin:nzmax,n)+tr_arr(nzmin:nzmax,n,tr_num)* &
                                    (hnode(nzmin:nzmax,n)-hnode_new(nzmin:nzmax,n))
        tr_arr(nzmin:nzmax,n,tr_num)=tr_arr(nzmin:nzmax,n,tr_num)+ &
                                    del_ttf(nzmin:nzmax,n)/hnode_new(nzmin:nzmax,n)
        ! WHY NOT ??? --> whats advantage of above --> tested it --> the upper 
        ! equation has a 30% smaller nummerical drift
        !tr_arr(1:nzmax,n,tr_num)=(hnode(1:nzmax,n)*tr_arr(1:nzmax,n,tr_num)+ &
        !                        del_ttf(1:nzmax,n))/hnode_new(1:nzmax,n)
    end do
    
    !___________________________________________________________________________
    if (i_vert_diff) then
        ! do vertical diffusion: implicite 
        call diff_ver_part_impl_ale(tr_num, mesh)
        
    end if
    
    !We DO not set del_ttf to zero because it will not be used in this timestep anymore
    !init_tracers will set it to zero for the next timestep
    !init_tracers will set it to zero for the next timestep
    if (smooth_bh_tra) then
       call diff_part_bh(tr_arr(:,:,tr_num), mesh) ! alpply biharmonic diffusion (implemented as filter)                                                
    end if
end subroutine diff_tracers_ale
!
!
!===============================================================================
!Vertical diffusive flux(explicit scheme):                                                                            
subroutine diff_ver_part_expl_ale(tr_num, mesh)
    use o_ARRAYS
    use g_forcing_arrays
    use MOD_MESH
    use g_PARSUP
    use g_config,only: dt
    
    implicit none 
    type(t_mesh), intent(in) , target :: mesh    
    real(kind=WP)            :: vd_flux(mesh%nl-1)
    real(kind=WP)            :: rdata,flux,rlx
    integer                  :: nz,nl1,ul1, tr_num,n
    real(kind=WP)            :: zinv1,Ty

#include "associate_mesh.h"
    !___________________________________________________________________________    
    do n=1, myDim_nod2D
        nl1=nlevels_nod2D(n)-1
        ul1=ulevels_nod2D(n)
        
        vd_flux=0._WP
        if (tr_num==1) then
            flux  = -heat_flux(n)/vcpw
            rdata =  Tsurf(n)
            rlx   =  surf_relax_T
        elseif (tr_num==2) then
            flux  =  virtual_salt(n)+relax_salt(n)- real_salt_flux(n)*is_nonlinfs
        else
            flux  = 0._WP
            rdata = 0._WP
            rlx=0._WP
        endif
        
        !_______________________________________________________________________
        !Surface forcing
        !!PS vd_flux(1)= flux
        vd_flux(ul1)= flux
        
        !_______________________________________________________________________
        !!PS do nz=2,nl1
        do nz=ul1+1,nl1
            !___________________________________________________________________
            zinv1=1.0_WP/(Z_3d_n(nz-1,n)-Z_3d_n(nz,n))
            
            !___________________________________________________________________
!            Ty= Kd(4,nz-1,n)*(Z_3d_n(nz-1,n)-zbar_3d_n(nz,n))*zinv1 *neutral_slope(3,nz-1,n)**2 + &
!                Kd(4,nz,n)*(zbar_3d_n(nz,n)-Z_3d_n(nz,n))*zinv1 *neutral_slope(3,nz,n)**2
            
            vd_flux(nz) = (Kv(nz,n)+Ty)*(tr_arr(nz-1,n,tr_num)-tr_arr(nz,n,tr_num))*zinv1*area(nz,n)
            
        end do
        
        !_______________________________________________________________________
        !!PS do nz=1,nl1-1
        do nz=ul1,nl1-1
            del_ttf(nz,n) = del_ttf(nz,n) + (vd_flux(nz) - vd_flux(nz+1))/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))*dt/areasvol(nz,n)
        end do
        del_ttf(nl1,n) = del_ttf(nl1,n) + (vd_flux(nl1)/(zbar_3d_n(nl1,n)-zbar_3d_n(nl1+1,n)))*dt/areasvol(nl1,n)
        
    end do ! --> do n=1, myDim_nod2D
end subroutine diff_ver_part_expl_ale
!
!
!===============================================================================
! vertical diffusivity augmented with Redi contribution [vertical flux of K(3,3)*d_zT]
subroutine diff_ver_part_impl_ale(tr_num, mesh)
    use MOD_MESH
    use o_PARAM
    use o_ARRAYS
    use i_ARRAYS
    use g_PARSUP
    use g_CONFIG
    use g_forcing_arrays
    use o_mixing_KPP_mod !for ghats _GO_    
    use bc_surface_interface
        
    implicit none
    type(t_mesh), intent(in) , target :: mesh
!!PS     real(kind=WP)            :: bc_surface    
    real(kind=WP)            :: a(mesh%nl), b(mesh%nl), c(mesh%nl), tr(mesh%nl)
    real(kind=WP)            :: cp(mesh%nl), tp(mesh%nl)
    integer                  :: nz, n, nzmax,nzmin, tr_num
    real(kind=WP)            :: m, zinv, dt_inv, dz
    real(kind=WP)            :: rsss, Ty,Ty1, c1,zinv1,zinv2,v_adv
    real(kind=WP), external  :: TFrez  ! Sea water freeze temperature.
    real(kind=WP)            :: isredi=0._WP
    logical                  :: do_wimpl=.true.
#include "associate_mesh.h"

    !___________________________________________________________________________
    if ((trim(tra_adv_lim)=='FCT') .OR. (.not. w_split)) do_wimpl=.false.
    
    if (Redi) isredi=1._WP
    dt_inv=1.0_WP/dt
    Ty    =0.0_WP
    Ty1   =0.0_WP
    
    ! solve equation diffusion equation implicite part: 
    ! --> h^(n+0.5)* (T^(n+0.5)-Tstar) = dt*( K_33*d/dz*(T^(n+0.5)-Tstar) + K_33*d/dz*Tstar )
    ! -->   dTnew = T^(n+0.5)-Tstar
    ! -->   h^(n+0.5)* (dTnew) = dt*(K_33*d/dz*dTnew) + K_33*dt*d/dz*Tstar 
    ! -->   h^(n+0.5)* (dTnew) = dt*(K_33*d/dz*dTnew) + RHS 
    ! -->   solve for dT_new
    !    
    !    ----------- zbar_1, V_1 (Skalar Volume), A_1 (Area of edge),  no Cavity A1==V1, yes Cavity A1 !=V1
    ! Z_1 o T_1
    !    ----------- zbar_2, V_2
    ! Z_2 o T_2
    !    ----------- zbar_3, V_3
    ! Z_3 o T_3
    !    ----------- zbar_4
    !        :
    ! --> Difference Quotient at Volume _2:  ddTnew_2/dt + d/dz*K_33 d/dz*dTnew_2 = 0 --> homogene solution 
    ! V2*dTnew_2 *h^(n+0.5) = -dt * [ (dTnew_1-dTnew_2)/(Z_1-Z_2)*A_2 + (dTnew_2-dTnew_3)/(Z_2-Z_3)*A_3 ] + RHS
    !    dTnew_2 *h^(n+0.5) = -dt * [ (dTnew_1-dTnew_2)/(Z_1-Z_2)*A_2/V_2 + (dTnew_2-dTnew_3)/(Z_2-Z_3)*A_3/V_2 ] + RHS
    !                                                  |                                 |
    !                                                  v                                 v
    !                                         diffusive flux towards             diffusive flux towards
    !                                         T_2 trough boundary V2             T_2 trough boundary V3 
    !    
    ! --> solve coefficents for homogene part   
    !    dTnew_2 *h^(n+0.5) = -dt * [ a*dTnew_1 + b*dTnew_2 + c*dTnew_3 ] 
    !
    ! --> a = -dt*K_33/(Z_1-Z_2)*A_2/V_2
    ! 
    ! --> c = -dt*K_33/(Z_2-Z_3)*A_3/V_2
    !
    ! --> b = h^(n+0.5) -[ dt*K_33/(Z_1-Z_2)*A_2/V_2 + dt*K_33/(Z_2-Z_3)*A_3/V_2 ] = -(a+c) + h^(n+0.5)
    
    !___________________________________________________________________________
    ! loop over local nodes
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
!         zbar_n(nzmax)=zbar(nzmax)
        zbar_n(nzmax)=zbar_n_bot(n)
        Z_n(nzmax-1)=zbar_n(nzmax) + hnode_new(nzmax-1,n)/2.0_WP
        !!PS do nz=nzmax-1,2,-1
        do nz=nzmax-1,nzmin+1,-1
            zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
            Z_n(nz-1)  = zbar_n(nz) + hnode_new(nz-1,n)/2.0_WP
        end do
        !!PS zbar_n(1) = zbar_n(2) + hnode_new(1,n)
        zbar_n(nzmin) = zbar_n(nzmin+1) + hnode_new(nzmin,n)
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> surface layer 
        !!PS nz=1
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
        !!PS c(nz)=-(Kv(2,n)+Ty1)*zinv2*zinv*area(nz+1,n)/area(nz,n)
        c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv*area(nz+1,n)/areasvol(nz,n)
        b(nz)=-c(nz)+hnode_new(nz,n)      ! ale
        
        ! update from the vertical advection --> comes from splitting of vert 
        ! velocity into explicite and implicite contribution
        if (do_wimpl) then
            !!PS v_adv =zinv*area(nz+1,n)/areasvol(nz,n)
            !!PS b(nz) =b(nz)+Wvel_i(nz, n)*zinv-min(0._WP, Wvel_i(nz+1, n))*v_adv
            !!PS c(nz) =c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
            
            !___________________________________________________________________
            ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for 
            ! numerical reasons, to gurante that area/areasvol in case of no 
            ! cavity is ==1.0_WP
            v_adv =zinv* ( area(nz  ,n)/areasvol(nz,n) )
            b(nz) =b(nz)+Wvel_i(nz, n)*v_adv
            
            v_adv =zinv*area(nz+1,n)/areasvol(nz,n)
            b(nz) =b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
            c(nz) =c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
        end if        
        ! backup zinv2 for next depth level
        zinv1=zinv2
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> 2nd...nl-2 layer
        !!PS do nz=2, nzmax-2
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
            a(nz)=-(Kv(nz,n)  +Ty )*zinv1*zinv* ( area(nz  ,n)/areasvol(nz,n) ) 
            c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv*area(nz+1,n)/areasvol(nz,n)
            b(nz)=-a(nz)-c(nz)+hnode_new(nz,n)
            
            ! backup zinv2 for next depth level
            zinv1=zinv2
            
            ! update from the vertical advection
            if (do_wimpl) then
                !_______________________________________________________________
                ! use brackets when computing ( area(nz  ,n)/areasvol(nz,n) ) for 
                ! numerical reasons, to gurante that area/areasvol in case of no 
                ! cavity is ==1.0_WP   
                v_adv=zinv* ( area(nz  ,n)/areasvol(nz,n) )
                a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv
                b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv
                !!PS v_adv=v_adv*areasvol(nz+1,n)/areasvol(nz,n)
                v_adv=zinv*area(nz+1,n)/areasvol(nz,n)
                b(nz)=b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
                c(nz)=c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
            end if
        end do ! --> do nz=2, nzmax-2
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> nl-1 layer
        nz=nzmax-1 
        
        zinv=1.0_WP*dt   ! no ... /(zbar(nzmax-1)-zbar(nzmax)) because of ale
        
        ! calculate isoneutral diffusivity : Kd*s^2 --> K_33 = Kv + Kd*s^2
        Ty= (Z_n(nz-1)-zbar_n(nz))   *zinv1 *slope_tapered(3,nz-1,n)**2*Ki(nz-1,n) + &
            (zbar_n(nz)-Z_n(nz)) *zinv1 *slope_tapered(3,nz,n)**2  *Ki(nz,n)
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
        ! the rhs (inhomogene part): --> rhs = K_33*dt*d/dz*Tstar --> Tstar...tr_arr
        ! solve difference quotient for rhs --> tr
        !  RHS at Volume_2:
        !  
        !  RHS*V_2 = K_33*dt*(T_1-T_2)/(Z_1-Z_2)*V_2 - K_33*dt*(T_2-T_3)/(Z_2-Z_3)*V_3
        !          = -a*T_1 + (a+c)*T_2 - c*T_3
        !
        ! -+--> tr(1) =(a(1)+c(1))*tr_arr(1,n,tr_num)-c(1)*tr_arr(2,n,tr_num)
        !  |--> a(1)=0
        !!PS nz=1
        nz=nzmin
        dz=hnode_new(nz,n) ! It would be (zbar(nz)-zbar(nz+1)) if not ALE
        tr(nz)=-(b(nz)-dz)*tr_arr(nz,n,tr_num)-c(nz)*tr_arr(nz+1,n,tr_num)
        !tr(nz)=c(nz)*(tr_arr(nz,n,tr_num) - tr_arr(nz+1,n,tr_num))
        
        
        ! *******************************************************************
        ! nonlocal transport to the rhs (only T and S currently) _GO_
        ! *******************************************************************
        ! rsss will be used later to compute:
        ! 1. the virtual salinity flux 
        ! 2. the contribution from the nonlocal term in KPP for salinity
        if (tr_num==2) then 
            rsss=ref_sss
                if (ref_sss_local) rsss=tr_arr(1,n,2)
        end if
        
        !!PS do nz=2,nzmax-2
        do nz=nzmin+1,nzmax-2
            dz=hnode_new(nz,n)
            tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num)-(b(nz)-dz)*tr_arr(nz,n,tr_num)-c(nz)*tr_arr(nz+1,n,tr_num)
            !tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num) &
            !       -c(nz)*tr_arr(nz+1,n,tr_num) &
            !       +(a(nz)+c(nz))*tr_arr(nz,n,tr_num)
            
            ! *******************************************************************
            ! nonlocal transport to the rhs (only T and S currently) _GO_
            ! *******************************************************************
!leads to non conservation in 8th digit. needs to be checked!
!            if (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then
!                if (tr_num==1) then ! T
!                    tr(nz)=tr(nz)+(MIN(ghats(nz,n)*Kv(nz,n), 1.0_WP)-MIN(ghats(nz+1,n)*Kv(nz+1,n), 1.0_WP)*area(nz+1,n)/area(nz,n))*heat_flux(n)/vcpw
!                elseif (tr_num==2) then ! S
!                    tr(nz)=tr(nz)-(MIN(ghats(nz,n)*Kv(nz,n), 1.0_WP)-MIN(ghats(nz+1,n)*Kv(nz+1,n), 1.0_WP)*area(nz+1,n)/area(nz,n))*rsss*water_flux(n)
!                end if
!            end if 
        end do
        nz=nzmax-1
        dz=hnode_new(nz,n)
        tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num)-(b(nz)-dz)*tr_arr(nz,n,tr_num)
        !tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num)+a(nz)*tr_arr(nz,n,tr_num)
        
        !_______________________________________________________________________
        ! case of activated shortwave penetration into the ocean, ad 3d contribution
        if (use_sw_pene .and. tr_num==1) then
            !!PS do nz=1, nzmax-1
            do nz=nzmin, nzmax-1
                zinv=1.0_WP*dt  !/(zbar(nz)-zbar(nz+1)) ale!
                tr(nz)=tr(nz)+(sw_3d(nz, n)-sw_3d(nz+1, n)*area(nz+1,n)/areasvol(nz,n))*zinv
            end do
        end if
        
        !_______________________________________________________________________
        !  The first row contains also the boundary condition from heatflux, 
        !  freshwaterflux and relaxation terms
        !  --> tr_arr(1,n,1)*water_flux(n) : latent heatflux contribution due to 
        !      cell volume. If Volume decreases --> temp has to raise, if volume 
        !      expended --> temp has to decrease
        !                           (-)   ^                        (-)   ^ 
        !                            |    |                         |    | 
        !   IN MOMENT: heat_flux ~~~~|~~~~|~~~~   ,  water_flux ~~~~|~~~~|~~~~
        !  (BUT CHECK!)              |    |                         |    |
        !                            v   (+)                        v   (+) 
        !                            
        !!PS tr(1)= tr(1)+bc_surface(n, tracer_id(tr_num))        
        tr(nzmin)= tr(nzmin)+bc_surface(n, tracer_id(tr_num),mesh) 
        
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
        !!PS cp(1) = c(1)/b(1)
        !!PS tp(1) = tr(1)/b(1)
        cp(nzmin) = c(nzmin)/b(nzmin)
        tp(nzmin) = tr(nzmin)/b(nzmin)
        
        ! solve for vectors c-prime and t, s-prime
        !!PS do nz = 2,nzmax-1
        do nz = nzmin+1,nzmax-1
            m = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
        end do
        
        ! start with back substitution 
        tr(nzmax-1) = tp(nzmax-1)
        
        ! solve for x from the vectors c-prime and d-prime
        !!PS do nz = nzmax-2, 1, -1
        do nz = nzmax-2, nzmin, -1
            tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
        end do
        
        !_______________________________________________________________________
        ! update tracer
        ! tr ... dTnew = T^(n+0.5) - T*
        !!PS do nz=1,nzmax-1
        do nz=nzmin,nzmax-1
            ! tr_arr - before ... T*
            tr_arr(nz,n,tr_num)=tr_arr(nz,n,tr_num)+tr(nz)
            ! tr_arr - after ... T^(n+0.5) = dTnew + T* = T^(n+0.5) - T* + T* 
            
        end do
        
    end do ! --> do n=1,myDim_nod2D   
end subroutine diff_ver_part_impl_ale

!
!
!===============================================================================
subroutine ver_sinking_recom_benthos(tr_num,mesh)
    use o_ARRAYS
    use g_PARSUP
    use MOD_MESH
    USE o_param
    use g_config
    use g_comm_auto
    USE O_MESH
    use g_forcing_arrays
    use g_support
use ver_sinking_recom_benthos_interface
#if defined(__recom)
    USE REcoM_GloVar
    use recom_config !, recom_debug
    use g_support
#endif
    IMPLICIT NONE
    type(t_mesh), intent(in) , target  :: mesh
    integer                   :: elem,k, tr_num
    integer                   :: nl1,ul1,nz,n,nzmin, nzmax, net
    real(kind=WP)             :: Vben(mesh%nl),  aux(mesh%nl-1), flux(mesh%nl), add_benthos_2d(myDim_nod2D)
    integer                   :: nlevels_nod2D_minimum
    real(kind=WP)             :: tv
#include "associate_mesh.h"

   do n=1, myDim_nod2D ! needs exchange_nod in the end
        nl1=nlevels_nod2D(n)-1
        ul1=ulevels_nod2D(n)
        
        aux=0._WP
        Vben=0._WP
        add_benthos_2d=0._WP

! 1) Calculate sinking velociy for vertical sinking case
! ******************************************************
        if (tracer_id(tr_num)==1007 .or. &  !idetn
            tracer_id(tr_num)==1008 .or. &  !idetc
            tracer_id(tr_num)==1017 .or. &  !idetsi
            tracer_id(tr_num)==1021 ) then  !idetcal
            
            Vben = VDet
	    if (allow_var_sinking) Vben = Vdet_a * abs(zbar_3d_n(:,n)) + VDet

        elseif(tracer_id(tr_num)==1004 .or. &  !iphyn
               tracer_id(tr_num)==1005 .or. &  !iphyc
               tracer_id(tr_num)==1020 .or. &  !iphycal
               tracer_id(tr_num)==1006 ) then  !ipchl

            Vben = VPhy
	    if (allow_var_sinking) Vben = Vdet_a * abs(zbar_3d_n(:,n)) + VPhy

        elseif(tracer_id(tr_num)==1013 .or. &  !idian
               tracer_id(tr_num)==1014 .or. &  !idiac
               tracer_id(tr_num)==1016 .or. &  !idiasi
               tracer_id(tr_num)==1015 ) then  !idchl

            Vben = VDia
	    if (allow_var_sinking) Vben = Vdet_a * abs(zbar_3d_n(:,n)) + VDia
      
! Constant vertical sinking for the second detritus class
! *******************************************************

!   if (REcoM_Second_Zoo) then ! No variable sinking
        elseif(tracer_id(tr_num)==1025 .or. &  !idetz2n
               tracer_id(tr_num)==1026 .or. &  !idetz2c
               tracer_id(tr_num)==1027 .or. &  !idetz2si
               tracer_id(tr_num)==1028 ) then  !idetz2calc

               Vben = VDet_zoo2
        endif

        Vben= Vben/SecondsPerDay ! conversion [m/d] --> [m/s] (vertical velocity, note that it is positive here)

        k=nod_in_elem2D_num(n)
        ! Screening minimum depth in neigbouring nodes around node n
        nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

        do nz=nlevels_nod2D_minimum, nl1
           tv = tr_arr(nz,n,tr_num)*Vben(nz)
           aux(nz)= - tv*(area(nz,n)-area(nz+1,n))
        end do
        nz=nl1
        tv = tr_arr(nz,n,tr_num)*Vben(nz)
        aux(nz)= - tv*(area(nz+1,n))
        
        do nz=ul1,nl1
           str_bf(nz,n) = str_bf(nz,n) + (aux(nz))*dt/area(nz,n)/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
           add_benthos_2d(n) = add_benthos_2d(n) - (aux(nz+1))*dt
        end do

            ! N
            if( tracer_id(tr_num)==1004 .or. &  !iphyn
                tracer_id(tr_num)==1007 .or. &  !idetn
                tracer_id(tr_num)==1013 .or. &  !idian
                tracer_id(tr_num)==1025 ) then  !idetz2n
                Benthos(n,1)= Benthos(n,1) +  add_benthos_2d(n) ![mmol]
            endif
         
            ! C
            if( tracer_id(tr_num)==1005 .or. &  !iphyc
                tracer_id(tr_num)==1008 .or. &  !idetc
                tracer_id(tr_num)==1014 .or. &  !idiac
                tracer_id(tr_num)==1026 ) then  !idetz2c
                Benthos(n,2)= Benthos(n,2) + add_benthos_2d(n)
            endif

            ! Si
            if( tracer_id(tr_num)==1016 .or. &  !idiasi
                tracer_id(tr_num)==1017 .or. &  !idetsi
                tracer_id(tr_num)==1027 ) then  !idetz2si
                Benthos(n,3)= Benthos(n,3) + add_benthos_2d(n)
            endif

            ! Cal
            if( tracer_id(tr_num)==1020 .or. &  !iphycal
                tracer_id(tr_num)==1021 .or. &  !idetcal
                tracer_id(tr_num)==1028 ) then  !idetz2cal
                Benthos(n,4)= Benthos(n,4) + add_benthos_2d(n) 
            endif
end do

    do n=1, benthos_num
      call exchange_nod(Benthos(:,n))
    end do
end subroutine ver_sinking_recom_benthos

subroutine integrate_bottom(tflux,mesh)
    use o_ARRAYS
    use g_PARSUP
    use MOD_MESH
    USE o_param
    use g_config
    use g_comm_auto
    USE O_MESH
    use g_forcing_arrays
    use integrate_bottom_interface
#if defined(__recom)
    USE REcoM_GloVar
    use recom_config !, recom_debug
#endif
    IMPLICIT NONE
    type(t_mesh), intent(in) , target  :: mesh
    integer                            :: elem,k, tr_num
    integer                            :: nl1,ul1,nz,n
    real(kind=WP)                      :: tf, aux(mesh%nl-1)
    real(kind=WP), intent(inout)       :: tflux
    integer                            :: nlevels_nod2D_minimum
#include "associate_mesh.h"

   tf =0._WP
   do n=1, myDim_nod2D 
         tf=tf+Benthos(n,3)
   end do
   tflux=0.0_WP
   call MPI_AllREDUCE(tf, tflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
        MPI_COMM_FESOM, MPIerr)
end subroutine integrate_bottom
!
!
!===============================================================================
subroutine diff_ver_recom_expl(tr_num,mesh)
    use o_ARRAYS
    use g_PARSUP
    use MOD_MESH
    USE o_param
    use g_config
    use g_comm_auto
    USE O_MESH
    use g_forcing_arrays
use diff_ver_recom_expl_interface
#if defined(__recom)
    USE REcoM_GloVar
    use recom_config !, recom_debug
#endif
    IMPLICIT NONE
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: elem,k,tr_num
    integer                  :: n2,nl1,nl2,nz,n,id,ul1
    real(kind=WP)            :: vd_flux(mesh%nl)
    integer                  :: nlevels_nod2D_minimum
    real(kind=WP)            :: bottom_flux(myDim_nod2D+eDim_nod2D)

#include "associate_mesh.h"

bottom_flux = 0._WP
id = tracer_id(tr_num)

  SELECT CASE (id)
    CASE (1001)
      bottom_flux = GlodecayBenthos(:,1) ! DIN [mmolN/m^2/s]
    CASE (1002)
      bottom_flux = GlodecayBenthos(:,2) + GlodecayBenthos(:,4) ! DIC + calcification
    CASE (1003)
      bottom_flux = GlodecayBenthos(:,4) * 2.0_WP ! Alk
    CASE (1018)
      bottom_flux = GlodecayBenthos(:,3) ! Si
    CASE (1019)
      if(use_Fe2N) then 
        bottom_flux = GlodecayBenthos(:,1) * Fe2N_benthos ! Fe
      else
        bottom_flux = GlodecayBenthos(:,2) * Fe2C_benthos
      end if
    CASE (1022)
      bottom_flux = -GlodecayBenthos(:,2) * redO2C ! Oxy
    CASE (1033)
      if (ciso) then
        bottom_flux = GlodecayBenthos(:,5) + GlodecayBenthos(:,7) ! DIC_13 and Calc: DIC_13
      end if
    CASE (1034)
      if (ciso) then
        bottom_flux = GlodecayBenthos(:,6) + GlodecayBenthos(:,8) ! DIC_14 and Calc: DIC_14
      end if
    CASE DEFAULT
      if (mype==0) then
         if (mype==0) write(*,*) 'check specified in boundary conditions'
         if (mype==0) write(*,*) 'the model will stop!'
      end if
      call par_ex
      stop
  END SELECT

   do n=1, myDim_nod2D

        nl1=nlevels_nod2D(n)-1
        ul1=ulevels_nod2D(n)

        vd_flux=0._WP

        k=nod_in_elem2D_num(n)
        ! Screening minimum depth in neigbouring nodes around node n
        nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

        !_______________________________________________________________________
        ! Bottom flux
        do nz=nlevels_nod2D_minimum, nl1
            vd_flux(nz)=(area(nz,n)-area(nz+1,n))* bottom_flux(n)/(area(1,n))           
        end do
        nz=nl1
        vd_flux(nz+1)= (area(nz+1,n))* bottom_flux(n)/(area(1,n))
        !_______________________________________________________________________
        ! writing flux into rhs
        do nz=ul1,nl1
            ! flux contribute only the cell through its bottom !!!
!            dtr_bf(nz,n) = dtr_bf(nz,n) + vd_flux(nz+1)*dt/area(nz,n)/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
            dtr_bf(nz,n) = dtr_bf(nz,n) + vd_flux(nz+1)*dt/areasvol(nz,n)/hnode_new(nz,n)
        end do
    end do
end subroutine diff_ver_recom_expl
!
!
!===============================================================================
subroutine diff_ver_part_redi_expl(mesh)
    use o_ARRAYS
    use g_PARSUP
    use MOD_MESH
    USE o_param
    use g_config
    use g_comm_auto
    IMPLICIT NONE
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: elem,k
    integer                  :: n2,nl1,ul1,nl2,nz,n
    real(kind=WP)            :: Tx, Ty
    real(kind=WP)            :: tr_xynodes(2,mesh%nl-1,myDim_nod2D+eDim_nod2D), vd_flux(mesh%nl)

#include "associate_mesh.h"

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
    
    ! call exchange_nod_begin(tr_xynodes)  !NR the halo is not needed
    
    do n=1, myDim_nod2D
        nl1=nlevels_nod2D(n)-1
        ul1=ulevels_nod2D(n)
        vd_flux=0._WP
        
        !_______________________________________________________________________
        zbar_n=0.0_WP
        Z_n   =0.0_WP
!         zbar_n(nl1+1)=zbar(nl1+1)
        zbar_n(nl1+1)=zbar_n_bot(n)
        Z_n(nl1)=zbar_n(nl1+1) + hnode_new(nl1,n)/2.0_WP
        !!PS do nz=nl1, 2, -1
        do nz=nl1, ul1+1, -1
            zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
            Z_n(nz-1)  = zbar_n(nz) + hnode_new(nz-1,n)/2.0_WP
        end do
        !!PS zbar_n(1) = zbar_n(2) + hnode_new(1,n)
        zbar_n(ul1) = zbar_n(ul1+1) + hnode_new(ul1,n)
        
        !_______________________________________________________________________
        !!PS do nz=2,nl1
        do nz=ul1+1,nl1
            vd_flux(nz)=(Z_n(nz-1)-zbar_n(nz))*(slope_tapered(1,nz-1,n)*tr_xynodes(1,nz-1,n)+slope_tapered(2,nz-1,n)*tr_xynodes(2,nz-1,n))*Ki(nz-1,n)
            vd_flux(nz)=vd_flux(nz)+&
                        (zbar_n(nz)-Z_n(nz))  *(slope_tapered(1,nz,n)  *tr_xynodes(1,nz,n)  +slope_tapered(2,nz,n)  *tr_xynodes(2,nz,n))  *Ki(nz,n)
            vd_flux(nz)=vd_flux(nz)/(Z_n(nz-1)-Z_n(nz))*area(nz,n)
        enddo
        !!PS do nz=1,nl1
        do nz=ul1,nl1
            del_ttf(nz,n) = del_ttf(nz,n)+(vd_flux(nz) - vd_flux(nz+1))*dt/areasvol(nz,n)
        enddo
    end do
end subroutine diff_ver_part_redi_expl!
!
!
!===============================================================================
subroutine diff_part_hor_redi(mesh)
    use o_ARRAYS
    use g_PARSUP
    use MOD_MESH
    use o_param
    use g_config
    IMPLICIT NONE
    type(t_mesh), intent(in) , target :: mesh
    real(kind=WP)            :: deltaX1,deltaY1,deltaX2,deltaY2
    integer                  :: edge
    integer                  :: n2,nl1,ul1,nl2,ul2,nl12,ul12,nz,el(2),elnodes(3),n,enodes(2)
    real(kind=WP)            :: c, Fx, Fy,Tx, Ty, Tx_z, Ty_z, SxTz, SyTz, Tz(2)
    real(kind=WP)            :: rhs1(mesh%nl-1), rhs2(mesh%nl-1), Kh, dz
    real(kind=WP)            :: isredi=0._WP

#include "associate_mesh.h"

    if (Redi) isredi=1._WP
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
        !Kh=elem_area(el(1))
        !_______________________________________________________________________
        nl2=0
        ul2=0
        if (el(2)>0) then 
            !Kh=0.5_WP*(Kh+elem_area(el(2)))
            nl2=nlevels(el(2))-1
            ul2=ulevels(el(2))
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
        endif
        !Kh=K_hor*Kh/scale_area
        !_______________________________________________________________________
        nl12=min(nl1,nl2)
        ul12=max(ul1,ul2)
        
        !_______________________________________________________________________
        ! (A)
        do nz=ul1,ul12-1
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
        !!PS do nz=1,nl12
        do nz=ul12,nl12
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
        do nz=nl12+1,nl1
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
        do nz=nl12+1,nl2
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
        !!PS del_ttf(1:nl12,enodes(1))=del_ttf(1:nl12,enodes(1))+rhs1(1:nl12)*dt/area(1:nl12,enodes(1))
        !!PS del_ttf(1:nl12,enodes(2))=del_ttf(1:nl12,enodes(2))+rhs2(1:nl12)*dt/area(1:nl12,enodes(2))
        del_ttf(ul12:nl12,enodes(1))=del_ttf(ul12:nl12,enodes(1))+rhs1(ul12:nl12)*dt/areasvol(ul12:nl12,enodes(1))
        del_ttf(ul12:nl12,enodes(2))=del_ttf(ul12:nl12,enodes(2))+rhs2(ul12:nl12)*dt/areasvol(ul12:nl12,enodes(2))
        
    end do
end subroutine diff_part_hor_redi
!
!
!===============================================================================
SUBROUTINE diff_part_bh(ttf, mesh)
    use o_ARRAYS
    use g_PARSUP
    use MOD_MESH
    use O_MESH
    use o_param
    use g_config
    use g_comm_auto

    IMPLICIT NONE
    type(t_mesh),  intent(in),    target :: mesh
    real(kind=WP), intent(inout), target :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP)                        :: u1, v1, len, vi, tt, ww 
    integer                              :: nz, ed, el(2), en(2), k, elem, nl1
    real(kind=WP), allocatable           :: temporary_ttf(:,:)

#include "associate_mesh.h"

    ed=myDim_nod2D+eDim_nod2D
    allocate(temporary_ttf(nl-1, ed))

    temporary_ttf=0.0_8
    DO ed=1, myDim_edge2D+eDim_edge2D
       if (myList_edge2D(ed)>edge2D_in) cycle
       el=edge_tri(:,ed)
       en=edges(:,ed)
       len=sqrt(sum(elem_area(el)))
       nl1=maxval(nlevels_nod2D_min(en))-1
       DO  nz=1,nl1
           u1=UV(1, nz,el(1))-UV(1, nz,el(2))
           v1=UV(2, nz,el(1))-UV(2, nz,el(2))
           vi=u1*u1+v1*v1
           tt=ttf(nz,en(1))-ttf(nz,en(2))
           vi=sqrt(max(gamma0, max(gamma1*sqrt(vi), gamma2*vi))*len)
           !vi=sqrt(max(sqrt(u1*u1+v1*v1),0.04)*le)  ! 10m^2/s for 10 km (0.04 h/50)
           !vi=sqrt(10.*le)
           tt=tt*vi
           temporary_ttf(nz,en(1))=temporary_ttf(nz,en(1))-tt
           temporary_ttf(nz,en(2))=temporary_ttf(nz,en(2))+tt
       END DO 
    END DO
    call exchange_nod(temporary_ttf)
    ! ===========
    ! Second round: 
    ! ===========
    DO ed=1, myDim_edge2D+eDim_edge2D
       if (myList_edge2D(ed)>edge2D_in) cycle
          el=edge_tri(:,ed)
          en=edges(:,ed)
          len=sqrt(sum(elem_area(el)))
          nl1=maxval(nlevels_nod2D_min(en))-1
          DO  nz=1,nl1
              u1=UV(1, nz,el(1))-UV(1, nz,el(2))
              v1=UV(2, nz,el(1))-UV(2, nz,el(2))
              vi=u1*u1+v1*v1
              tt=temporary_ttf(nz,en(1))-temporary_ttf(nz,en(2))
              vi=sqrt(max(gamma0, max(gamma1*sqrt(vi), gamma2*vi))*len)
              !vi=sqrt(max(sqrt(u1*u1+v1*v1),0.04)*le)  ! 10m^2/s for 10 km (0.04 h/50)
              !vi=sqrt(10.*le) 
              tt=-tt*vi*dt
              ttf(nz,en(1))=ttf(nz,en(1))-tt/area(nz,en(1))
              ttf(nz,en(2))=ttf(nz,en(2))+tt/area(nz,en(2))
          END DO 
    END DO  
    deallocate(temporary_ttf)
end subroutine diff_part_bh
!
!
!===============================================================================
! this function returns a boundary conditions for a specified thacer ID and surface node
! ID = 0 and 1 are reserved for temperature and salinity
FUNCTION bc_surface(n, id, mesh)
  use MOD_MESH
  USE o_ARRAYS
  USE g_forcing_arrays
  USE g_PARSUP, only: mype, par_ex
  USE g_config
#if defined(__recom)
USE REcoM_GloVar
use recom_config, only: ciso, recom_debug
use REcoM_declarations
use REcoM_ciso
#endif

  implicit none
  
  type(t_mesh), intent(in) , target :: mesh  
  REAL(kind=WP)       :: bc_surface
  integer, intent(in) :: n, id
  character(len=10)   :: id_string

  !  --> is_nonlinfs=1.0 for zelvel,zstar ....                            
  !  --> is_nonlinfs=0.0 for linfs
  SELECT CASE (id)
    CASE (0)
        bc_surface=-dt*(heat_flux(n)/vcpw + tr_arr(mesh%ulevels_nod2D(n),n,1)*water_flux(n)*is_nonlinfs)
    CASE (1)
        ! --> real_salt_flux(:): salt flux due to containment/releasing of salt
        !     by forming/melting of sea ice
        bc_surface= dt*(virtual_salt(n) & !--> is zeros for zlevel/zstar
                    + relax_salt(n) - real_salt_flux(n)*is_nonlinfs)
#if defined(__recom)
    CASE (1001) ! DIN
        bc_surface= dt*(AtmNInput(n)  + RiverDIN2D(n) * is_riverinput + ErosionTON2D(n) * is_erosioninput)
!  if (mype==0) then
!     write(*,*) '____________________________________________________________'
!     write(*,*) ' --> DIN_surface,  = ', bc_surface
!  endif
    CASE (1002) ! DIC
        bc_surface= dt*(GloCO2flux_seaicemask(n) + RiverDIC2D(n) * is_riverinput + ErosionTOC2D(n) * is_erosioninput)
!  if (mype==0) then
!     write(*,*) '____________________________________________________________'
!     write(*,*) ' --> DIC_surface,  = ', bc_surface
!  endif
    CASE (1003) ! Alk
        ! --> Here we need the alkalinity flux
        bc_surface= dt*(virtual_alk(n) &  
                    + relax_alk(n) + RiverAlk2D(n) * is_riverinput)
!  if (mype==0) then
!     write(*,*) '____________________________________________________________'
!     write(*,*) ' --> Alk_surface,  = ', bc_surface
!  endif
    CASE (1004:1010)
        bc_surface=0.0_WP
    CASE (1011) ! DON
        bc_surface= dt*RiverDON2D(n) * is_riverinput
    CASE (1012) ! DOC
        bc_surface= dt*RiverDOC2D(n) * is_riverinput
    CASE (1013:1017)
        bc_surface=0.0_WP
    CASE (1018) ! DSi
        bc_surface=dt*(RiverDSi2D(n) * is_riverinput + ErosionTSi2D(n) * is_erosioninput)
    CASE (1019) ! Fe
        bc_surface= dt*AtmFeInput(n)
!  if (mype==0) then
!     write(*,*) '____________________________________________________________'
!     write(*,*) ' --> Fe_surface,  = ', bc_surface
!  endif
    CASE (1020:1021) ! Si
        bc_surface=0.0_WP  ! OG added bc for recom fields 
    CASE (1022) ! OXY
        bc_surface= dt*GloO2flux_seaicemask(n)
!  if (mype==0) then
!     write(*,*) '____________________________________________________________'
!     write(*,*) ' --> DIC_surface,  = ', bc_surface
!  endif
    CASE (1023:1032)
        bc_surface=0.0_WP  ! OG added bc for recom fields 
!ciso adapted by MB
    CASE (1033) ! DIC_13
         if (ciso) then
           bc_surface= dt*GloCO2flux_seaicemask_13(n)
         else
           bc_surface=0.0_WP
         end if
         if (recom_debug .and. mype==0) then
             write(*,*) '____________________________________________________________'
             write(*,*) ' --> DIC_13_surface,  = ', bc_surface
         endif
    CASE (1034) ! DIC_14
         if (ciso) then
           bc_surface= dt*GloCO2flux_seaicemask_14(n)
         else
           bc_surface=0.0_WP
         end if
         if (recom_debug .and. mype==0) then
             write(*,*) '____________________________________________________________'
             write(*,*) ' --> DIC_14_surface,  = ', bc_surface
         endif
    CASE (1035:1099)
        bc_surface=0.0_WP  ! OG added bc for recom fields - adapted to ciso by MB 
    CASE (1102:1299)
        bc_surface=0.0_WP  ! added by MB for ciso
!ciso adapted by MB
#endif 
    CASE (101) ! apply boundary conditions to tracer ID=101
        bc_surface= dt*(prec_rain(n))! - real_salt_flux(n)*is_nonlinfs)
    CASE (301)
        bc_surface=0.0_WP
    CASE (302)
        bc_surface=0.0_WP
    CASE (303)
        bc_surface=0.0_WP
    CASE DEFAULT
      if (mype==0) then
         write (id_string, "(I3)") id
         if (mype==0) write(*,*) 'invalid ID '//trim(id_string)//' specified in boundary conditions'
         if (mype==0) write(*,*) 'the model will stop!'
      end if
      call par_ex
      stop
  END SELECT
  RETURN
END FUNCTION
