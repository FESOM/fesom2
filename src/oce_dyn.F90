
! A set of routines for computing the horizonlal viscosity
! the control parameters (their default values) are:
! gamma0 (0.01 [m/s]), gamma1 (0.1 [no dim.]), gamma2 (10.[s/m]), Div_c [1.], Leith_c[1.?]
! 1. gamma0 has the dimension of velocity. It should be as small as possible, but in any case smaller than 0.01 m/s. 
!    All major ocean circulation models are stable with harmonic viscosity 0.01*len.
! 2. gamma1 is nondimensional. In commonly used Leith or Smagorinsky parameterizations it is C/pi^2=0.1 (C is about 1). 
!    We therefore try to follow this, allowing some adjustments (because our mesh is triangular, our resolution is different, etc.). 
!    We however, try to keep gamma1<0.1
! 3. gamma2 is dimensional (1/velocity). If it is 10, then the respective term dominates starting from |u|=0.1 m/s an so on. It is only used in: 
!    (5) visc_filt_bcksct, (6) visc_filt_bilapl, (7) visc_filt_bidiff
! 4. Div_c  =1.    should be default
! 5. Leith_c=?    (need to be adjusted)
module visc_filt_bcksct_interface
  interface
    subroutine visc_filt_bcksct(dynamics, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      USE MOD_DYN
      type(t_dyn)   , intent(inout), target :: dynamics
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
      
    end subroutine
  end interface
end module
module visc_filt_bilapl_interface
  interface
    subroutine visc_filt_bilapl(dynamics, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      USE MOD_DYN
      type(t_dyn)   , intent(inout), target :: dynamics
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
      
    end subroutine
  end interface
end module
module visc_filt_bidiff_interface
  interface
    subroutine visc_filt_bidiff(dynamics, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      USE MOD_DYN
      type(t_dyn)   , intent(inout), target :: dynamics
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
      
    end subroutine
  end interface
end module
!!PS module visc_filt_dbcksc_interface
!!PS   interface
!!PS     subroutine visc_filt_dbcksc(dynamics, partit, mesh)
!!PS       use mod_mesh
!!PS       USE MOD_PARTIT
!!PS       USE MOD_PARSUP
!!PS       USE MOD_DYN
!!PS       type(t_dyn)   , intent(inout), target :: dynamics
!!PS       type(t_partit), intent(inout), target :: partit
!!PS       type(t_mesh)  , intent(in)   , target :: mesh
!!PS       
!!PS     end subroutine
!!PS   end interface
!!PS end module
!!PS module backscatter_coef_interface
!!PS   interface
!!PS     subroutine backscatter_coef(dynamics, partit, mesh)
!!PS       use mod_mesh
!!PS       USE MOD_PARTIT
!!PS       USE MOD_PARSUP
!!PS       USE MOD_DYN
!!PS       type(t_dyn)   , intent(inout), target :: dynamics
!!PS       type(t_partit), intent(inout), target :: partit
!!PS       type(t_mesh)  , intent(in)   , target :: mesh
!!PS       
!!PS     end subroutine
!!PS   end interface
!!PS end module
!!PS module uke_update_interface
!!PS   interface
!!PS     subroutine uke_update(dynamics, partit, mesh)
!!PS       use mod_mesh
!!PS       USE MOD_PARTIT
!!PS       USE MOD_PARSUP
!!PS       USE MOD_DYN
!!PS       type(t_dyn)   , intent(inout), target :: dynamics
!!PS       type(t_partit), intent(inout), target :: partit
!!PS       type(t_mesh)  , intent(in)   , target :: mesh
!!PS       
!!PS     end subroutine
!!PS   end interface
!!PS end module

module relative_vorticity_interface
  interface
    subroutine relative_vorticity(dynamics, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      use MOD_DYN
      type(t_dyn)   , intent(inout), target :: dynamics
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
      
    end subroutine
  end interface
end module

! 
! Contains routines needed for computations of dynamics.
! includes: update_vel, compute_vel_nodes
!_______________________________________________________________________________
SUBROUTINE update_vel(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    use g_comm_auto
    IMPLICIT NONE
    integer       :: elem, elnodes(3), nz, m, nzmax, nzmin
    real(kind=WP) :: eta(3) 
    real(kind=WP) :: Fx, Fy
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:), pointer :: eta_n, d_eta

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV=>dynamics%uv(:,:,:)
    UV_rhs=>dynamics%uv_rhs(:,:,:)
    eta_n=>dynamics%eta_n(:)
    d_eta=>dynamics%d_eta(:)
        
    DO elem=1, myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        eta=-g*theta*dt*d_eta(elnodes)
        Fx=sum(gradient_sca(1:3,elem)*eta)
        Fy=sum(gradient_sca(4:6,elem)*eta)
        nzmin = ulevels(elem)
        nzmax = nlevels(elem)
        !!PS DO nz=1, nlevels(elem)-1
        DO nz=nzmin, nzmax-1
            UV(1,nz,elem)= UV(1,nz,elem) + UV_rhs(1,nz,elem) + Fx
            UV(2,nz,elem)= UV(2,nz,elem) + UV_rhs(2,nz,elem) + Fy
        END DO
    END DO
    eta_n=eta_n+d_eta
    call exchange_elem(UV, partit)
end subroutine update_vel
!
!
!_______________________________________________________________________________
subroutine compute_vel_nodes(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE o_PARAM
    use g_comm_auto
    IMPLICIT NONE
    integer            :: n, nz, k, elem, nln, uln, nle, ule
    real(kind=WP)      :: tx, ty, tvol
    
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    real(kind=WP), dimension(:,:,:), pointer :: UV, UVnode
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV=>dynamics%uv(:,:,:)
    UVnode=>dynamics%uvnode(:,:,:)

    DO n=1, myDim_nod2D 
        uln = ulevels_nod2D(n)
        nln = nlevels_nod2D(n) 
        !!PS DO nz=1, nlevels_nod2D(n)-1
        DO nz=uln, nln-1
            tvol=0.0_WP
            tx  =0.0_WP
            ty  =0.0_WP
            DO k=1, nod_in_elem2D_num(n)
                elem=nod_in_elem2D(k,n)
                ule = ulevels(elem)
                nle = nlevels(elem)
                !!PS if (nlevels(elem)-1<nz) cycle
                if (nle-1<nz .or. nz<ule) cycle
                tvol=tvol+elem_area(elem)
                tx=tx+UV(1,nz,elem)*elem_area(elem)
                ty=ty+UV(2,nz,elem)*elem_area(elem)
            END DO
            UVnode(1,nz,n)=tx/tvol
            UVnode(2,nz,n)=ty/tvol
        END DO
    END DO
    call exchange_nod(UVnode, partit)
end subroutine compute_vel_nodes
!
!
!_______________________________________________________________________________
subroutine viscosity_filter(option, dynamics, partit, mesh)
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    use visc_filt_bcksct_interface
    use visc_filt_bilapl_interface
    use visc_filt_bidiff_interface
!!PS     use visc_filt_dbcksc_interface
!!PS     use backscatter_coef_interface
    use g_backscatter
    IMPLICIT NONE 
    integer                               :: option
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit

    ! Driving routine 
    ! Background viscosity is selected in terms of Vl, where V is 
    ! background velocity scale and l is the resolution. V is 0.005 
    ! or 0.01, perhaps it would be better to pass it as a parameter.

    ! h_viscosity_leiht needs vorticity, so vorticity array should be 
    ! allocated. At present, there are two rounds of smoothing in 
    ! h_viscosity. 
    SELECT CASE (option)
        CASE (5)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call visc_filt_bcksct'//achar(27)//'[0m'
            call visc_filt_bcksct(dynamics, partit, mesh)
        CASE (6)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call visc_filt_bilapl'//achar(27)//'[0m'
            call visc_filt_bilapl(dynamics, partit, mesh)
        CASE (7)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call visc_filt_bidiff'//achar(27)//'[0m'
            call visc_filt_bidiff(dynamics, partit, mesh)
        CASE (8)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call backscatter_coef'//achar(27)//'[0m'
            call backscatter_coef(partit, mesh)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call visc_filt_dbcksc'//achar(27)//'[0m'
            call visc_filt_dbcksc(dynamics, partit, mesh) 
        CASE DEFAULT
            if (partit%mype==0) write(*,*) 'mixing scheme with option ' , option, 'has not yet been implemented'
            call par_ex(partit%MPI_COMM_FESOM, partit%mype)
            stop
    END SELECT
end subroutine viscosity_filter  
!
!
!_______________________________________________________________________________
SUBROUTINE visc_filt_bcksct(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    USE g_comm_auto
    IMPLICIT NONE

    real(kind=8)  :: u1, v1, len, vi 
    integer       :: nz, ed, el(2), nelem(3),k, elem, nzmin, nzmax
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c, U_b, V_b
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV     => dynamics%uv(    :,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)
    U_b    => dynamics%work%u_b(:,:)
    V_b    => dynamics%work%v_b(:,:)

    ! An analog of harmonic viscosity operator.
    ! Same as visc_filt_h, but with the backscatter. 
    ! Here the contribution from squared velocities is added to the viscosity.    
    ! The contribution from boundary edges is neglected (free slip). 

    ed=myDim_elem2D+eDim_elem2D
    allocate(U_b(nl-1,ed), V_b(nl-1, ed))
    ed=myDim_nod2D+eDim_nod2D
    allocate(U_c(nl-1,ed), V_c(nl-1,ed))
    U_b=0.0_WP
    V_b=0.0_WP
    U_c=0.0_WP
    V_c=0.0_WP
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el=edge_tri(:,ed)
        len=sqrt(sum(elem_area(el)))
        nzmax = minval(nlevels(el))
        nzmin = maxval(ulevels(el))
        !!PS DO  nz=1,minval(nlevels(el))-1
        DO  nz=nzmin,nzmax-1
            u1=UV(1,nz,el(1))-UV(1,nz,el(2))
            v1=UV(2,nz,el(1))-UV(2,nz,el(2))
            vi=dt*max(gamma0, max(gamma1*sqrt(u1*u1+v1*v1), gamma2*(u1*u1+v1*v1)))*len
!            vi=dt*max(gamma0, gamma1*max(sqrt(u1*u1+v1*v1), gamma2*(u1*u1+v1*v1)))*len 
            !here gamma2 is dimensional (1/velocity). If it is 10, then the respective term dominates starting from |u|=0.1 m/s an so on.
            u1=u1*vi
            v1=v1*vi
            U_b(nz,el(1))=U_b(nz,el(1))-u1/elem_area(el(1))
            U_b(nz,el(2))=U_b(nz,el(2))+u1/elem_area(el(2))
            V_b(nz,el(1))=V_b(nz,el(1))-v1/elem_area(el(1))
            V_b(nz,el(2))=V_b(nz,el(2))+v1/elem_area(el(2))
        END DO 
    END DO
    call exchange_elem(U_b, partit)
    call exchange_elem(V_b, partit)
    ! ===========
    ! Compute smoothed viscous term: 
    ! ===========
    DO ed=1, myDim_nod2D 
        nzmin = ulevels_nod2D(ed)
        nzmax = nlevels_nod2D(ed)
        !!PS DO nz=1, nlevels_nod2D(ed)-1
        DO nz=nzmin, nzmax-1
            vi=0.0_WP
            u1=0.0_WP
            v1=0.0_WP
            DO k=1, nod_in_elem2D_num(ed)
                elem=nod_in_elem2D(k,ed)
                vi=vi+elem_area(elem)
                u1=u1+U_b(nz,elem)*elem_area(elem)
                v1=v1+V_b(nz,elem)*elem_area(elem)
            END DO
            U_c(nz,ed)=u1/vi
            V_c(nz,ed)=v1/vi
        END DO
    END DO
    call exchange_nod(U_c, partit)
    call exchange_nod(V_c, partit)
    do ed=1, myDim_elem2D
        nelem=elem2D_nodes(:,ed)
        nzmin = ulevels(ed)
        nzmax = nlevels(ed)
        !!PS Do nz=1, nlevels(ed)-1
        Do nz=nzmin, nzmax-1
            UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+U_b(nz,ed) -easy_bs_return*sum(U_c(nz,nelem))/3.0_WP
            UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+V_b(nz,ed) -easy_bs_return*sum(V_c(nz,nelem))/3.0_WP
        END DO
    end do
end subroutine visc_filt_bcksct
!
!
!_______________________________________________________________________________
! Strictly energy dissipative and momentum conserving version
! Viscosity depends on velocity Laplacian, i.e., on an analog of
! the Leith viscosity (Lapl==second derivatives)
! \nu=|3u_c-u_n1-u_n2-u_n3|*sqrt(S_c)/100. There is an additional term
! in viscosity that is proportional to the velocity amplitude squared.
! The coefficient has to be selected experimentally.
SUBROUTINE visc_filt_bilapl(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    USE g_comm_auto
    IMPLICIT NONE
    real(kind=8)  :: u1, v1, vi, len
    integer       :: ed, el(2), nz, nzmin, nzmax
    
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV => dynamics%uv(:,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)

    ed=myDim_elem2D+eDim_elem2D
    allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
    U_c=0.0_WP
    V_c=0.0_WP
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el=edge_tri(:,ed)
        nzmin = maxval(ulevels(el))
        nzmax = minval(nlevels(el))
        !!PS DO  nz=1,minval(nlevels(el))-1
        DO  nz=nzmin,nzmax-1
            u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
            v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
            U_c(nz,el(1))=U_c(nz,el(1))-u1
            U_c(nz,el(2))=U_c(nz,el(2))+u1
            V_c(nz,el(1))=V_c(nz,el(1))-v1
            V_c(nz,el(2))=V_c(nz,el(2))+v1
        END DO 
    END DO
 
    Do ed=1,myDim_elem2D
        len=sqrt(elem_area(ed))
        nzmin = ulevels(ed)
        nzmax = nlevels(ed)
        !!PS Do nz=1,nlevels(ed)-1
        Do nz=nzmin,nzmax-1
            ! vi has the sense of harmonic viscosity coef. because of 
            ! division by area in the end 
            u1=U_c(nz,ed)**2+V_c(nz,ed)**2
            vi=max(gamma0, max(gamma1*sqrt(u1), gamma2*u1))*len*dt
            U_c(nz,ed)=-U_c(nz,ed)*vi                             
            V_c(nz,ed)=-V_c(nz,ed)*vi
        END DO
    end do
    
    call exchange_elem(U_c, partit)
    call exchange_elem(V_c, partit)
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el=edge_tri(:,ed)
        nzmin = maxval(ulevels(el))
        nzmax = minval(nlevels(el))
        !!PS DO  nz=1,minval(nlevels(el))-1
        DO  nz=nzmin,nzmax-1
            u1=(U_c(nz,el(1))-U_c(nz,el(2)))
            v1=(V_c(nz,el(1))-V_c(nz,el(2)))
            UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
            UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
            UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
            UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
        END DO 
    END DO  
    deallocate(V_c,U_c)

end subroutine visc_filt_bilapl
!
!
!_______________________________________________________________________________
! Strictly energy dissipative and momentum conserving version
! Viscosity depends on velocity differences, and is introduced symmetrically 
! into both stages of biharmonic operator
! On each edge, \nu=sqrt(|u_c1-u_c2|*sqrt(S_c1+S_c2)/100)
! The effect is \nu^2
! Quadratic in velocity term can be introduced if needed.
SUBROUTINE visc_filt_bidiff(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    USE g_comm_auto
    IMPLICIT NONE
    real(kind=8)  :: u1, v1, vi, len
    integer       :: ed, el(2), nz, nzmin, nzmax
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV => dynamics%uv(:,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)
    !
    ed=myDim_elem2D+eDim_elem2D
    allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
    U_c=0.0_WP
    V_c=0.0_WP
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el=edge_tri(:,ed)
        len=sqrt(sum(elem_area(el)))
        nzmin = maxval(ulevels(el))
        nzmax = minval(nlevels(el))
        !!PS DO  nz=1,minval(nlevels(el))-1
        DO  nz=nzmin,nzmax-1
            u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
            v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
            vi=u1*u1+v1*v1
            vi=sqrt(max(gamma0, max(gamma1*sqrt(vi), gamma2*vi))*len)
            ! vi=sqrt(max(gamma0, gamma1*max(sqrt(vi), gamma2*vi))*len)
            u1=u1*vi
            v1=v1*vi
            U_c(nz,el(1))=U_c(nz,el(1))-u1
            U_c(nz,el(2))=U_c(nz,el(2))+u1
            V_c(nz,el(1))=V_c(nz,el(1))-v1
            V_c(nz,el(2))=V_c(nz,el(2))+v1
        END DO 
    END DO
    
    call exchange_elem(U_c, partit)
    call exchange_elem(V_c, partit)
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el=edge_tri(:,ed)
        len=sqrt(sum(elem_area(el)))
        nzmin = maxval(ulevels(el))
        nzmax = minval(nlevels(el))
        !!PS DO  nz=1,minval(nlevels(el))-1
        DO  nz=nzmin,nzmax-1
            u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
            v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
            vi=u1*u1+v1*v1
            vi=-dt*sqrt(max(gamma0, max(gamma1*sqrt(vi), gamma2*vi))*len)
            ! vi=-dt*sqrt(max(gamma0, gamma1*max(sqrt(vi), gamma2*vi))*len)
            u1=vi*(U_c(nz,el(1))-U_c(nz,el(2)))
            v1=vi*(V_c(nz,el(1))-V_c(nz,el(2)))
            UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
            UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
            UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
            UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
        END DO 
    END DO
    deallocate(V_c, U_c)

end subroutine visc_filt_bidiff
!!PS !
!!PS !
!!PS !_______________________________________________________________________________
!!PS SUBROUTINE visc_filt_dbcksc(dynamics, partit, mesh)
!!PS USE MOD_MESH
!!PS USE MOD_PARTIT
!!PS USE MOD_PARSUP
!!PS use MOD_DYN
!!PS USE o_ARRAYS, only: v_back, UV_dis_tend, UV_total_tend, UV_back_tend, &
!!PS                     uke, uke_dif
!!PS USE o_PARAM
!!PS USE g_CONFIG
!!PS USE g_comm_auto
!!PS USE g_support
!!PS USE uke_update_interface
!!PS IMPLICIT NONE
!!PS 
!!PS real(kind=8)  :: u1, v1, le(2), len, crosslen, vi, uke1 
!!PS integer       :: nz, ed, el(2)
!!PS !!PS real(kind=8), allocatable  :: U_c(:,:), V_c(:,:)
!!PS real(kind=8)  , allocatable  :: UV_back(:,:,:), UV_dis(:,:,:), uke_d(:,:) 
!!PS real(kind=8)  , allocatable  :: uuu(:)
!!PS type(t_dyn)   , intent(inout), target :: dynamics
!!PS type(t_partit), intent(inout), target :: partit
!!PS type(t_mesh)  , intent(in)   , target :: mesh
!!PS real(kind=WP) , dimension(:,:,:), pointer :: UV, UV_rhs
!!PS real(kind=WP) , dimension(:,:)  , pointer :: U_c, V_c
!!PS #include "associate_part_def.h"
!!PS #include "associate_mesh_def.h"
!!PS #include "associate_part_ass.h"
!!PS #include "associate_mesh_ass.h"
!!PS UV     => dynamics%uv(:,:,:)
!!PS UV_rhs => dynamics%uv_rhs(:,:,:)
!!PS U_c    => dynamics%work%u_c(:,:)
!!PS V_c    => dynamics%work%v_c(:,:)
!!PS 
!!PS  ! An analog of harmonic viscosity operator.  
!!PS  ! It adds to the rhs(0) Visc*(u1+u2+u3-3*u0)/area
!!PS  ! on triangles, which is Visc*Laplacian/4 on equilateral triangles. 
!!PS  ! The contribution from boundary edges is neglected (free slip). 
!!PS  ! Filter is applied twice. 
!!PS 
!!PS ed=myDim_elem2D+eDim_elem2D
!!PS allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
!!PS allocate(UV_back(2,nl-1,ed), UV_dis(2,nl-1, ed)) 
!!PS allocate(uke_d(nl-1,ed)) 
!!PS allocate(uuu(ed))
!!PS  
!!PS  U_c=0.0_8
!!PS  V_c=0.0_8
!!PS  UV_back=0.0_8
!!PS  UV_dis=0.0_8
!!PS  uke_d=0.0_8
!!PS 
!!PS  DO ed=1, myDim_edge2D+eDim_edge2D
!!PS     if(myList_edge2D(ed)>edge2D_in) cycle
!!PS     el=edge_tri(:,ed)
!!PS     DO  nz=1,minval(nlevels(el))-1
!!PS      u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
!!PS      v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
!!PS      
!!PS      U_c(nz,el(1))=U_c(nz,el(1))-u1
!!PS      U_c(nz,el(2))=U_c(nz,el(2))+u1
!!PS      V_c(nz,el(1))=V_c(nz,el(1))-v1
!!PS      V_c(nz,el(2))=V_c(nz,el(2))+v1
!!PS     END DO 
!!PS  END DO
!!PS  
!!PS 
!!PS  Do ed=1,myDim_elem2D
!!PS     len=sqrt(elem_area(ed))                     
!!PS     len=dt*len/30.0_8
!!PS     Do nz=1,nlevels(ed)-1
!!PS      ! vi has the sense of harmonic viscosity coefficient because of 
!!PS      ! the division by area in the end 
!!PS      ! ====
!!PS      ! Case 1 -- an analog to the third-order upwind (vi=|u|l/12)
!!PS      ! ====
!!PS      vi=max(0.2_8,sqrt(UV(1,nz,ed)**2+UV(2,nz,ed)**2))*len 
!!PS      U_c(nz,ed)=-U_c(nz,ed)*vi                             
!!PS      V_c(nz,ed)=-V_c(nz,ed)*vi
!!PS     END DO
!!PS  end do
!!PS 
!!PS 
!!PS  call exchange_elem(U_c, partit)
!!PS  call exchange_elem(V_c, partit) 
!!PS 
!!PS  DO ed=1, myDim_edge2D+eDim_edge2D
!!PS     if(myList_edge2D(ed)>edge2D_in) cycle
!!PS     el=edge_tri(:,ed)
!!PS     le=edge_dxdy(:,ed)
!!PS     le(1)=le(1)*sum(elem_cos(el))*0.25_8
!!PS     len=sqrt(le(1)**2+le(2)**2)*r_earth
!!PS     le(1)=edge_cross_dxdy(1,ed)-edge_cross_dxdy(3,ed)
!!PS     le(2)=edge_cross_dxdy(2,ed)-edge_cross_dxdy(4,ed)
!!PS     crosslen=sqrt(le(1)**2+le(2)**2) 
!!PS     DO  nz=1,minval(nlevels(el))-1
!!PS      vi=dt*len*(v_back(nz,el(1))+v_back(nz,el(2)))/crosslen
!!PS      !if(mype==0) write(*,*) 'vi ', vi , ' and ed' , ed
!!PS      !if(mype==0) write(*,*) 'dt*len/crosslen ', dt*len/crosslen, ' and ed' , ed
!!PS      !vi=max(vi,0.005*len*dt) ! This helps to reduce noise in places where 
!!PS                               ! Visc is small and decoupling might happen 
!!PS      !Backscatter contribution
!!PS      u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))*vi
!!PS      v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))*vi
!!PS      
!!PS      !UKE diffusion
!!PS      vi=dt*len*(K_back*sqrt(elem_area(el(1))/scale_area)+K_back*sqrt(elem_area(el(2))/scale_area))/crosslen
!!PS 
!!PS      uke1=(uke(nz,el(1))-uke(nz,el(2)))*vi
!!PS 
!!PS      
!!PS      UV_back(1,nz,el(1))=UV_back(1,nz,el(1))-u1/elem_area(el(1))
!!PS      UV_back(1,nz,el(2))=UV_back(1,nz,el(2))+u1/elem_area(el(2))
!!PS      UV_back(2,nz,el(1))=UV_back(2,nz,el(1))-v1/elem_area(el(1))
!!PS      UV_back(2,nz,el(2))=UV_back(2,nz,el(2))+v1/elem_area(el(2))  
!!PS      
!!PS      !Correct scaling for the diffusion?
!!PS      uke_d(nz,el(1))=uke_d(nz,el(1))-uke1/elem_area(el(1))
!!PS      uke_d(nz,el(2))=uke_d(nz,el(2))+uke1/elem_area(el(2))
!!PS      
!!PS      
!!PS      
!!PS      !Biharmonic contribution
!!PS      u1=(U_c(nz,el(1))-U_c(nz,el(2)))
!!PS      v1=(V_c(nz,el(1))-V_c(nz,el(2)))
!!PS 
!!PS      UV_dis(1,nz,el(1))=UV_dis(1,nz,el(1))-u1/elem_area(el(1))
!!PS      UV_dis(1,nz,el(2))=UV_dis(1,nz,el(2))+u1/elem_area(el(2))
!!PS      UV_dis(2,nz,el(1))=UV_dis(2,nz,el(1))-v1/elem_area(el(1))
!!PS      UV_dis(2,nz,el(2))=UV_dis(2,nz,el(2))+v1/elem_area(el(2))
!!PS      
!!PS     END DO 
!!PS  END DO
!!PS  
!!PS call exchange_elem(UV_back, partit)
!!PS 
!!PS DO  nz=1, nl-1
!!PS     uuu=0.0_8
!!PS     uuu=UV_back(1,nz,:)
!!PS     call smooth_elem(uuu,smooth_back_tend, partit, mesh)
!!PS     UV_back(1,nz,:)=uuu
!!PS     uuu=0.0_8
!!PS     uuu=UV_back(2,nz,:)
!!PS     call smooth_elem(uuu,smooth_back_tend, partit, mesh)
!!PS     UV_back(2,nz,:)=uuu 
!!PS END DO
!!PS 
!!PS  DO ed=1, myDim_elem2D
!!PS      DO  nz=1,nlevels(ed)-1
!!PS      UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+UV_dis(1,nz,ed)+UV_back(1,nz,ed)
!!PS      UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+UV_dis(2,nz,ed)+UV_back(2,nz,ed)               
!!PS      END DO 
!!PS  END DO
!!PS 
!!PS  UV_dis_tend=UV_dis!+UV_back
!!PS  UV_total_tend=UV_dis+UV_back
!!PS  UV_back_tend=UV_back
!!PS  uke_dif=uke_d
!!PS  
!!PS  call uke_update(dynamics, partit, mesh)
!!PS  deallocate(V_c,U_c)
!!PS  deallocate(UV_dis,UV_back) 
!!PS  deallocate(uke_d)
!!PS  deallocate(uuu)
!!PS 
!!PS end subroutine visc_filt_dbcksc
!!PS !
!!PS !
!!PS !_______________________________________________________________________________
!!PS SUBROUTINE backscatter_coef(partit, mesh)
!!PS USE MOD_MESH
!!PS USE MOD_PARTIT
!!PS USE MOD_PARSUP
!!PS USE o_ARRAYS
!!PS USE o_PARAM
!!PS USE g_CONFIG
!!PS use g_comm_auto
!!PS IMPLICIT NONE
!!PS type(t_mesh),   intent(in),    target :: mesh
!!PS type(t_partit), intent(inout), target :: partit
!!PS integer                               :: elem, nz
!!PS #include "associate_part_def.h"
!!PS #include "associate_mesh_def.h"
!!PS #include "associate_part_ass.h"
!!PS #include "associate_mesh_ass.h"
!!PS 
!!PS !Potentially add the Rossby number scaling to the script...
!!PS !check if sign is right! Different in the Jansen paper
!!PS !Also check with the normalization by area; as before we use element length sqrt(2*elem_area(ed))
!!PS 
!!PS v_back=0.0_8
!!PS DO  elem=1, myDim_elem2D 
!!PS      DO  nz=1,nlevels(elem)-1
!!PS !v_back(1,ed)=c_back*sqrt(2.0_WP*elem_area(ed))*sqrt(max(2.0_WP*uke(1,ed),0.0_WP))*(3600.0_WP*24.0_WP/tau_c)*4.0_WP/sqrt(2.0_WP*elem_area(ed))**2 !*sqrt(max(2.0_WP*uke(1,ed),0.0_WP))
!!PS !v_back(nz,elem)=-c_back*sqrt(4._8/sqrt(3.0_8)*elem_area(elem))*sqrt(max(2.0_8*uke(nz,elem),0.0_8)) !Is the scaling correct
!!PS v_back(nz,elem)=min(-c_back*sqrt(elem_area(elem))*sqrt(max(2.0_8*uke(nz,elem),0.0_8)),0.2*elem_area(elem)/dt) !Is the scaling correct
!!PS !Scaling by sqrt(2*elem_area) or sqrt(elem_area)?
!!PS      END DO
!!PS END DO
!!PS 
!!PS call exchange_elem(v_back, partit)
!!PS 
!!PS end subroutine backscatter_coef
!!PS !
!!PS !
!!PS !_______________________________________________________________________________
!!PS SUBROUTINE uke_update(dynamics, partit, mesh)
!!PS USE MOD_MESH
!!PS USE MOD_PARTIT
!!PS USE MOD_PARSUP
!!PS use MOD_DYN
!!PS USE o_ARRAYS, only: uke_rhs, uke_dif, uke_back, uke_dis, uke, UV_dis_tend, uv_back_tend, uke_rhs_old, &
!!PS                     bvfreq, coriolis_node
!!PS USE o_PARAM
!!PS USE g_CONFIG
!!PS use g_comm_auto
!!PS USE g_support
!!PS USE g_rotate_grid
!!PS IMPLICIT NONE
!!PS 
!!PS !I had to change uke(:) to uke(:,:) to make output and restart work!!
!!PS 
!!PS !Why is it necessary to implement the length of the array? It doesn't work without!
!!PS !integer, intent(in)        :: t_levels
!!PS type(t_dyn)   , intent(inout), target :: dynamics
!!PS type(t_partit), intent(inout), target :: partit
!!PS type(t_mesh)  , intent(in)   , target :: mesh
!!PS 
!!PS real(kind=8)		   :: hall, h1_eta, hnz, vol
!!PS integer			   :: elnodes(3), nz, ed, edi, node, j, elem, q
!!PS real(kind=8), allocatable  :: uuu(:), work_array(:), U_work(:,:), V_work(:,:), rosb_array(:,:), work_uv(:)
!!PS integer			   :: kk, nzmax, el
!!PS real(kind=8)       :: c1, rosb, vel_u, vel_v, vel_uv, scaling, reso
!!PS real*8             :: c_min=0.5, f_min=1.e-6, r_max=200000., ex, ey, a1, a2, len_reg, dist_reg(2) ! Are those values still correct?
!!PS real(kind=WP), dimension(:,:,:), pointer :: UV
!!PS #include "associate_part_def.h"
!!PS #include "associate_mesh_def.h"
!!PS #include "associate_part_ass.h"
!!PS #include "associate_mesh_ass.h" 
!!PS UV => dynamics%uv(:,:,:)
!!PS 
!!PS !rosb_dis=1._8 !Should be variable to control how much of the dissipated energy is backscattered
!!PS !rossby_num=2
!!PS 
!!PS ed=myDim_elem2D+eDim_elem2D
!!PS allocate(uuu(ed)) 
!!PS 
!!PS uke_back=0.0_8
!!PS uke_dis=0.0_8
!!PS DO ed=1, myDim_elem2D
!!PS DO nz=1, nlevels(ed)-1  
!!PS    uke_dis(nz,ed)=(UV(1,nz,ed)*UV_dis_tend(1,nz,ed)+UV(2,nz,ed)*UV_dis_tend(2,nz,ed))   
!!PS    uke_back(nz,ed)=(UV(1,nz,ed)*UV_back_tend(1,nz,ed)+UV(2,nz,ed)*UV_back_tend(2,nz,ed))
!!PS END DO
!!PS END DO
!!PS 
!!PS DO  nz=1,nl-1
!!PS     uuu=0.0_8
!!PS     uuu=uke_back(nz,:)
!!PS     call smooth_elem(uuu,smooth_back, partit, mesh) !3) ?
!!PS     uke_back(nz,:)=uuu
!!PS END DO
!!PS 
!!PS 
!!PS 
!!PS !Timestepping use simple backward timestepping; all components should have dt in it, unless they need it twice
!!PS !Amplitudes should be right given the correction of the viscosities; check for all, also for biharmonic
!!PS !uke(1,ed)=uke(1,ed)-uke_dis(1,ed)-uke_back(1,ed)+uke_dif(1,ed)
!!PS ed=myDim_elem2D+eDim_elem2D
!!PS allocate(U_work(nl-1,myDim_nod2D+eDim_nod2D),V_work(nl-1,myDim_nod2D+eDim_nod2D))
!!PS allocate(work_uv(myDim_nod2D+eDim_nod2D))
!!PS allocate(rosb_array(nl-1,ed))
!!PS call exchange_elem(UV, partit)
!!PS rosb_array=0._8
!!PS DO nz=1, nl-1
!!PS     work_uv=0._WP
!!PS     DO node=1, myDim_nod2D
!!PS        vol=0._WP
!!PS        U_work(nz,node)=0._WP 
!!PS        V_work(nz,node)=0._WP 
!!PS        DO j=1, nod_in_elem2D_num(node)
!!PS            elem=nod_in_elem2D(j, node)
!!PS            U_work(nz,node)=U_work(nz,node)+UV(1,nz,elem)*elem_area(elem)
!!PS            V_work(nz,node)=V_work(nz,node)+UV(2,nz,elem)*elem_area(elem)
!!PS            vol=vol+elem_area(elem)
!!PS        END DO
!!PS        U_work(nz,node)=U_work(nz,node)/vol
!!PS        V_work(nz,node)=U_work(nz,node)/vol
!!PS     END DO
!!PS     work_uv=U_work(nz,:)
!!PS     call exchange_nod(work_uv, partit)
!!PS     U_work(nz,:)=work_uv
!!PS     work_uv=V_work(nz,:)
!!PS     call exchange_nod(work_uv, partit)
!!PS     V_work(nz,:)=work_uv    
!!PS END DO
!!PS 
!!PS     DO el=1,myDim_elem2D
!!PS      DO nz=1, nlevels(el)-1     
!!PS         rosb_array(nz,el)=sqrt((sum(gradient_sca(1:3,el)*U_work(nz,elem2D_nodes(1:3,el)))-&
!!PS               sum(gradient_sca(4:6, el)*V_work(nz,elem2D_nodes(1:3,el))))**2+&
!!PS               (sum(gradient_sca(4:6, el)*U_work(nz,elem2D_nodes(1:3,el)))+&
!!PS               sum(gradient_sca(1:3, el)*V_work(nz,elem2D_nodes(1:3,el))))**2)
!!PS !        hall=hall+hnz
!!PS      END DO
!!PS !     rosb_array(el)=rosb_array(el)/hall
!!PS     END DO
!!PS DO ed=1, myDim_elem2D
!!PS     scaling=1._WP
!!PS     IF(uke_scaling) then
!!PS       reso=sqrt(elem_area(ed)*4._wp/sqrt(3._wp))
!!PS       rosb=0._wp
!!PS       elnodes=elem2D_nodes(:, ed)   
!!PS       DO kk=1,3
!!PS         c1=0._wp
!!PS         nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(elnodes(kk)), elnodes(kk))), 1)
!!PS         !Vertical average; same scaling in the vertical
!!PS         DO nz=1, nzmax-1
!!PS          c1=c1+hnode_new(nz,elnodes(kk))*(sqrt(max(bvfreq(nz,elnodes(kk)), 0._WP))+sqrt(max(bvfreq(nz+1,elnodes(kk)), 0._WP)))/2.
!!PS         END DO
!!PS         c1=max(c_min, c1/pi) !ca. first baroclinic gravity wave speed limited from below by c_min
!!PS         !Cutoff K_GM depending on (Resolution/Rossby radius) ratio
!!PS         rosb=rosb+min(c1/max(abs(coriolis_node(elnodes(kk))), f_min), r_max)
!!PS       END DO
!!PS       rosb=rosb/3._8
!!PS       scaling=1._WP/(1._WP+(uke_scaling_factor*reso/rosb))!(4._wp*reso/rosb))
!!PS     END IF
!!PS     
!!PS     DO nz=1, nlevels(ed)-1  
!!PS     elnodes=elem2D_nodes(:,ed)
!!PS     
!!PS     !Taking out that one place where it is always weird (Pacific Southern Ocean)
!!PS     !Should not really be used later on, once we fix the issue with the 1/4 degree grid
!!PS     if(.not. (TRIM(which_toy)=="soufflet")) then
!!PS       call elem_center(ed, ex, ey)
!!PS       !a1=-104.*rad
!!PS       !a2=-49.*rad
!!PS       call g2r(-104.*rad, -49.*rad, a1, a2)
!!PS       dist_reg(1)=ex-a1
!!PS       dist_reg(2)=ey-a2
!!PS       call trim_cyclic(dist_reg(1))
!!PS       dist_reg(1)=dist_reg(1)*elem_cos(ed)
!!PS       dist_reg=dist_reg*r_earth
!!PS       len_reg=sqrt(dist_reg(1)**2+dist_reg(2)**2)
!!PS     
!!PS     
!!PS       !if(mype==0) write(*,*) 'len_reg ', len_reg , ' and dist_reg' , dist_reg, ' and ex, ey', ex, ey, ' and a ', a1, a2
!!PS       rosb_array(nz,ed)=rosb_array(nz,ed)/max(abs(sum(coriolis_node(elnodes(:)))), f_min)
!!PS       !uke_dif(nz, ed)=scaling*(1-exp(-len_reg/300000))*1._8/(1._8+rosb_array(nz,ed)/rosb_dis)!UV_dif(1,ed)
!!PS       uke_dis(nz,ed)=scaling*(1-exp(-len_reg/300000))*1._8/(1._8+rosb_array(nz,ed)/rosb_dis)*uke_dis(nz,ed)
!!PS     else
!!PS       rosb_array(nz,ed)=rosb_array(nz,ed)/max(abs(sum(coriolis_node(elnodes(:)))), f_min)
!!PS       !uke_dif(nz, ed)=scaling*1._8/(1._8+rosb_array(nz,ed)/rosb_dis)!UV_dif(1,ed)
!!PS       uke_dis(nz,ed)=scaling*1._8/(1._8+rosb_array(nz,ed)/rosb_dis)*uke_dis(nz,ed)
!!PS     end if
!!PS     
!!PS     END DO
!!PS END DO
!!PS deallocate(U_work, V_work)
!!PS deallocate(rosb_array)
!!PS deallocate(work_uv)
!!PS call exchange_elem(uke_dis, partit)
!!PS DO nz=1, nl-1  
!!PS     uuu=uke_dis(nz,:)
!!PS     call smooth_elem(uuu,smooth_dis, partit, mesh)
!!PS     uke_dis(nz,:)=uuu
!!PS END DO
!!PS DO ed=1, myDim_elem2D
!!PS     DO  nz=1,nlevels(ed)-1
!!PS     uke_rhs_old(nz,ed)=uke_rhs(nz,ed)
!!PS     uke_rhs(nz,ed)=-uke_dis(nz,ed)-uke_back(nz,ed)+uke_dif(nz,ed)
!!PS     uke(nz,ed)=uke(nz,ed)+1.5_8*uke_rhs(nz,ed)-0.5_8*uke_rhs_old(nz,ed)
!!PS     END DO
!!PS END DO
!!PS call exchange_elem(uke, partit)
!!PS 
!!PS deallocate(uuu)
!!PS end subroutine uke_update
!
!
!_______________________________________________________________________________
subroutine relative_vorticity(dynamics, partit, mesh)
    USE o_ARRAYS, only: vorticity
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use g_comm_auto
    IMPLICIT NONE
    integer        :: n, nz, el(2), enodes(2), nl1, nl2, edge, ul1, ul2, nl12, ul12
    real(kind=WP)  :: deltaX1, deltaY1, deltaX2, deltaY2, c1
    
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    real(kind=WP), dimension(:,:,:), pointer :: UV
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    UV => dynamics%uv(:,:,:)

    !!PS DO n=1,myDim_nod2D
    !!PS    nl1 = nlevels_nod2D(n)-1
    !!PS    ul1 = ulevels_nod2D(n)
    !!PS    vorticity(ul1:nl1,n)=0.0_WP
    !!PS    !!PS DO nz=1, nlevels_nod2D(n)-1
    !!PS    !!PS    vorticity(nz,n)=0.0_WP
    !!PS    !!PS END DO
    !!PS END DO      
    vorticity(:,1:myDim_nod2D) = 0.0_WP
    DO edge=1,myDim_edge2D
                                    !! edge=myList_edge2D(m)
        enodes=edges(:,edge)
        el=edge_tri(:,edge)
        nl1=nlevels(el(1))-1
        ul1=ulevels(el(1))
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        nl2=0
        ul2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
            nl2=nlevels(el(2))-1
            ul2=ulevels(el(2))
        end if  
        nl12 = min(nl1,nl2)
        ul12 = max(ul1,ul2)
        
        DO nz=ul1,ul12-1
            c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
            vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
            vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
        END DO
        if (ul2>0) then
            DO nz=ul2,ul12-1
                c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
                vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
                vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
            END DO
        endif 
        !!PS DO nz=1,min(nl1,nl2)
        DO nz=ul12,nl12
            c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))- &
            deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
            vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
            vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
        END DO
        !!PS DO nz=min(nl1,nl2)+1,nl1
        DO nz=nl12+1,nl1
            c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
            vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
            vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
        END DO
        !!PS DO nz=min(nl1,nl2)+1,nl2
        DO nz=nl12+1,nl2
            c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
            vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
            vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
        END DO
    END DO
    
    ! vorticity = vorticity*area at this stage
    ! It is correct only on myDim nodes
    DO n=1,myDim_nod2D
                                !! n=myList_nod2D(m)
        ul1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        !!PS DO nz=1,nlevels_nod2D(n)-1
        DO nz=ul1,nl1-1
            vorticity(nz,n)=vorticity(nz,n)/areasvol(nz,n)
        END DO
    END DO      
    
    call exchange_nod(vorticity, partit)
    
! Now it the relative vorticity known on neighbors too
end subroutine relative_vorticity

