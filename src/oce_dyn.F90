
! A set of routines for computing the horizonlal viscosity
! the control parameters (their default values) are:
! dynamics%visc_gamma0 (0.01 [m/s]), dynamics%visc_gamma1 (0.1 [no dim.]), dynamics%visc_gamma2 (10.[s/m]), Div_c [1.], Leith_c[1.?]
! 1. dynamics%visc_gamma0 has the dimension of velocity. It should be as small as possible, but in any case smaller than 0.01 m/s. 
!    All major ocean circulation models are stable with harmonic viscosity 0.01*len.
! 2. dynamics%visc_gamma1 is nondimensional. In commonly used Leith or Smagorinsky parameterizations it is C/pi^2=0.1 (C is about 1). 
!    We therefore try to follow this, allowing some adjustments (because our mesh is triangular, our resolution is different, etc.). 
!    We however, try to keep dynamics%visc_gamma1<0.1
! 3. dynamics%visc_gamma2 is dimensional (1/velocity). If it is 10, then the respective term dominates starting from |u|=0.1 m/s an so on. It is only used in: 
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
            vi=dt*max(dynamics%visc_gamma0,                         &
                      max(dynamics%visc_gamma1*sqrt(u1*u1+v1*v1),   &
                      dynamics%visc_gamma2*(u1*u1+v1*v1))           &
                    )*len
!            vi=dt*max(dynamics%visc_gamma0, dynamics%visc_gamma1*max(sqrt(u1*u1+v1*v1), dynamics%visc_gamma2*(u1*u1+v1*v1)))*len 
            !here dynamics%visc_gamma2 is dimensional (1/velocity). If it is 10, then the respective term dominates starting from |u|=0.1 m/s an so on.
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
            UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+U_b(nz,ed) -dynamics%visc_easybsreturn*sum(U_c(nz,nelem))/3.0_WP
            UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+V_b(nz,ed) -dynamics%visc_easybsreturn*sum(V_c(nz,nelem))/3.0_WP
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
            vi=max(dynamics%visc_gamma0,                &
                   max(dynamics%visc_gamma1*sqrt(u1),   &
                   dynamics%visc_gamma2*u1)             &
                  )*len*dt
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
            vi=sqrt(max(dynamics%visc_gamma0,           &
                    max(dynamics%visc_gamma1*sqrt(vi),  &
                    dynamics%visc_gamma2*vi)            &
                   )*len)
            ! vi=sqrt(max(dynamics%visc_gamma0, dynamics%visc_gamma1*max(sqrt(vi), dynamics%visc_gamma2*vi))*len)
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
            vi=-dt*sqrt(max(dynamics%visc_gamma0,           &
                        max(dynamics%visc_gamma1*sqrt(vi),  &
                        dynamics%visc_gamma2*vi)            &
                       )*len)
            ! vi=-dt*sqrt(max(dynamics%visc_gamma0, dynamics%visc_gamma1*max(sqrt(vi), dynamics%visc_gamma2*vi))*len)
            u1=vi*(U_c(nz,el(1))-U_c(nz,el(2)))
            v1=vi*(V_c(nz,el(1))-V_c(nz,el(2)))
            UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
            UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
            UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
            UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
        END DO 
    END DO
end subroutine visc_filt_bidiff

