module oce_adv_tra_ver_interfaces
  interface
! implicit 1st order upwind vertical advection with to solve for fct_LO
! updates the input tracer ttf
    subroutine adv_tra_vert_impl(ttf, w, mesh)
      use mod_mesh
      use g_PARSUP
      type(t_mesh),  intent(in), target  :: mesh
      real(kind=WP), intent(inout)       :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)          :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
    end subroutine
!===============================================================================
! 1st order upwind (explicit)
! returns flux given at vertical interfaces of scalar volumes
! IF init_zero=.TRUE.  : flux will be set to zero before computation
! IF init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_ver_upw1(ttf, w, do_Xmoment, mesh, flux, init_zero)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh
      integer,       intent(in)  :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
      real(kind=WP), intent(in)  :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)  :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout) :: flux(mesh%nl,  myDim_nod2D)
      logical, optional          :: init_zero
    end subroutine
!===============================================================================
! QR (4th order centerd)
! returns flux given at vertical interfaces of scalar volumes
! IF init_zero=.TRUE.  : flux will be set to zero before computation
! IF init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_ver_qr4c(ttf, w, do_Xmoment, mesh, num_ord, flux, init_zero)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh
      integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
      real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: flux(mesh%nl,  myDim_nod2D)
      logical, optional                 :: init_zero
    end subroutine
!===============================================================================
! Vertical advection with PPM reconstruction (5th order)
! returns flux given at vertical interfaces of scalar volumes
! IF init_zero=.TRUE.  : flux will be set to zero before computation
! IF init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
   subroutine adv_tra_vert_ppm(ttf, w, do_Xmoment, mesh, flux, init_zero)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh
      integer                           :: n, nz, nl1
      integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
      real(kind=WP)                     :: tvert(mesh%nl), tv
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: flux(mesh%nl,  myDim_nod2D)
      logical, optional                 :: init_zero
    end subroutine
! central difference reconstruction (2nd order, use only with FCT)
! returns flux given at vertical interfaces of scalar volumes
! IF init_zero=.TRUE.  : flux will be set to zero before computation
! IF init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_ver_cdiff(ttf, w, do_Xmoment, mesh, flux, init_zero)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh
      integer                           :: n, nz, nl1
      integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
      real(kind=WP)                     :: tvert(mesh%nl), tv
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: flux(mesh%nl,  myDim_nod2D)
      logical, optional                 :: init_zero
    end subroutine
  end interface
end module
!===============================================================================
subroutine adv_tra_vert_impl(ttf, w, mesh)
    use MOD_MESH
    use O_MESH
    use o_PARAM
    use o_ARRAYS
    use i_ARRAYS
    use g_PARSUP
    use g_CONFIG
    use g_forcing_arrays
    use o_mixing_KPP_mod !for ghats _GO_        
    
    implicit none
    type(t_mesh),  intent(in) , target :: mesh
    real(kind=WP), intent(inout)       :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)          :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
    real(kind=WP)                      :: a(mesh%nl), b(mesh%nl), c(mesh%nl), tr(mesh%nl)
    real(kind=WP)                      :: cp(mesh%nl), tp(mesh%nl)
    integer                            :: nz, n, nzmax, nzmin, tr_num
    real(kind=WP)                      :: m, zinv, dt_inv, dz
    real(kind=WP)                      :: c1, v_adv

#include "associate_mesh.h"

    dt_inv=1.0_WP/dt
    
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
        
        ! upper surface index, in case of cavity !=1
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
        Z_n(nzmax-1) =zbar_n(nzmax) + hnode_new(nzmax-1,n)/2.0_WP
        do nz=nzmax-1,nzmin+1,-1
            zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
            Z_n(nz-1)  = zbar_n(nz)   + hnode_new(nz-1,n)/2.0_WP
        end do
        zbar_n(nzmin) = zbar_n(nzmin+1) + hnode_new(nzmin,n)
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> surface layer 
        nz=nzmin
        
        ! 1/dz(nz)
        zinv=1.0_WP*dt    ! no .../(zbar(1)-zbar(2)) because of  ALE
        
        !!PS a(nz)=0.0_WP
        !!PS v_adv=zinv*areasvol(nz+1,n)/areasvol(nz,n)
        !!PS b(nz)= hnode_new(nz,n)+W(nz, n)*zinv-min(0._WP, W(nz+1, n))*v_adv
        !!PS c(nz)=-max(0._WP, W(nz+1, n))*v_adv
        
        a(nz)=0.0_WP
        v_adv=zinv*area(nz  ,n)/areasvol(nz,n)
        b(nz)= hnode_new(nz,n)+W(nz, n)*v_adv
        
        v_adv=zinv*area(nz+1,n)/areasvol(nz,n)
        b(nz)= b(nz)-min(0._WP, W(nz+1, n))*v_adv
        c(nz)=-max(0._WP, W(nz+1, n))*v_adv
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> 2nd...nl-2 layer
        do nz=nzmin+1, nzmax-2
            ! update from the vertical advection
            v_adv=zinv*area(nz  ,n)/areasvol(nz,n)
            a(nz)=min(0._WP, W(nz, n))*v_adv
            b(nz)=hnode_new(nz,n)+max(0._WP, W(nz, n))*v_adv
            
            v_adv=zinv*area(nz+1,n)/areasvol(nz,n)
            b(nz)=b(nz)-min(0._WP, W(nz+1, n))*v_adv
            c(nz)=     -max(0._WP, W(nz+1, n))*v_adv
        end do ! --> do nz=2, nzmax-2
        
        !_______________________________________________________________________
        ! Regular part of coefficients: --> nl-1 layer
        nz=nzmax-1
        ! update from the vertical advection
        !!PS a(nz)=                min(0._WP, W(nz, n))*zinv
        !!PS b(nz)=hnode_new(nz,n)+max(0._WP, W(nz, n))*zinv
        !!PS c(nz)=0.0_WP
        v_adv=zinv*area(nz  ,n)/areasvol(nz,n)
        a(nz)=                min(0._WP, W(nz, n))*v_adv
        b(nz)=hnode_new(nz,n)+max(0._WP, W(nz, n))*v_adv
        c(nz)=0.0_WP
        
        !_______________________________________________________________________
        nz=nzmin
        dz=hnode_new(nz,n) ! It would be (zbar(nz)-zbar(nz+1)) if not ALE
        tr(nz)=-(b(nz)-dz)*ttf(nz,n)-c(nz)*ttf(nz+1,n)
        
        do nz=nzmin+1,nzmax-2
            dz=hnode_new(nz,n)
            tr(nz)=-a(nz)*ttf(nz-1,n)-(b(nz)-dz)*ttf(nz,n)-c(nz)*ttf(nz+1,n)
        end do
        nz=nzmax-1
        dz=hnode_new(nz,n)
        tr(nz)=-a(nz)*ttf(nz-1,n)-(b(nz)-dz)*ttf(nz,n)
        
        !_______________________________________________________________________
        nz = nzmin
        cp(nz) = c(nz)/b(nz)
        tp(nz) = tr(nz)/b(nz)
        
        ! solve for vectors c-prime and t, s-prime
        do nz = nzmin+1,nzmax-1
            m = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
        end do
        
        !_______________________________________________________________________
        ! start with back substitution 
        tr(nzmax-1) = tp(nzmax-1)
        
        ! solve for x from the vectors c-prime and d-prime
        do nz = nzmax-2, nzmin, -1
            tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
        end do
        
        !_______________________________________________________________________
        ! update tracer
        do nz=nzmin,nzmax-1
            ttf(nz,n)=ttf(nz,n)+tr(nz)
        end do
    end do ! --> do n=1,myDim_nod2D
end subroutine adv_tra_vert_impl
!
!
!===============================================================================
subroutine adv_tra_ver_upw1(ttf, w, do_Xmoment, mesh, flux, init_zero)
    use g_config
    use MOD_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_forcing_arrays
    implicit none
    type(t_mesh),  intent(in), target :: mesh
    integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP)                     :: tvert(mesh%nl)
    integer                           :: n, nz, nzmax, nzmin
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl,  myDim_nod2D)
    logical, optional                 :: init_zero
#include "associate_mesh.h"

    if (present(init_zero))then
       if (init_zero) flux=0.0_WP
    else
       flux=0.0_WP
    end if

    do n=1, myDim_nod2D
       !_______________________________________________________________________
       nzmax=nlevels_nod2D(n)
       nzmin=ulevels_nod2D(n)
       
       !_______________________________________________________________________
       ! vert. flux at surface layer
       nz=nzmin
       flux(nz,n)=-W(nz,n)*ttf(nz,n)*area(nz,n)-flux(nz,n)
       
       !_______________________________________________________________________
       ! vert. flux at bottom layer --> zero bottom flux
       nz=nzmax
       flux(nz,n)= 0.0_WP-flux(nz,n)
       
       !_______________________________________________________________________
       ! Be carefull have to do vertical tracer advection here on old vertical grid
       ! also horizontal advection is done on old mesh (see helem contains old 
       ! mesh information)
       !_______________________________________________________________________
       ! vert. flux at remaining levels    
       do nz=nzmin+1,nzmax-1
          flux(nz,n)=-0.5*(                                                        &
                      (ttf(nz  ,n)**do_Xmoment)*(W(nz,n)+abs(W(nz,n)))+ &
                      (ttf(nz-1,n)**do_Xmoment)*(W(nz,n)-abs(W(nz,n))))*area(nz,n)-flux(nz,n)
       end do
    end do
end subroutine adv_tra_ver_upw1
!
!
!===============================================================================
subroutine adv_tra_ver_qr4c(ttf, w, do_Xmoment, mesh, num_ord, flux, init_zero)
    use g_config
    use MOD_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_forcing_arrays
    implicit none
    type(t_mesh),  intent(in), target :: mesh
    integer,       intent(in)    :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP), intent(in)    :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
    real(kind=WP), intent(in)    :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)    :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout) :: flux(mesh%nl,  myDim_nod2D)
    logical, optional            :: init_zero
    real(kind=WP)                :: tvert(mesh%nl)
    integer                      :: n, nz, nzmax, nzmin
    real(kind=WP)                :: Tmean, Tmean1, Tmean2
    real(kind=WP)                :: qc, qu, qd

#include "associate_mesh.h"

    if (present(init_zero))then
       if (init_zero) flux=0.0_WP
    else
       flux=0.0_WP
    end if

    do n=1, myDim_nod2D
       !_______________________________________________________________________
       nzmax=nlevels_nod2D(n)
       nzmin=ulevels_nod2D(n)
       !_______________________________________________________________________
       ! vert. flux at surface layer
       nz=nzmin
       flux(nz,n)=-ttf(nz,n)*W(nz,n)*area(nz,n)-flux(nz,n)
       
       !_______________________________________________________________________
       ! vert. flux 2nd layer --> centered differences
       nz=nzmin+1
       flux(nz,n)=-0.5_WP*(ttf(nz-1,n)+ttf(nz,n))*W(nz,n)*area(nz,n)-flux(nz,n)
       
       !_______________________________________________________________________
       ! vert. flux at bottom - 1 layer --> centered differences
       nz=nzmax-1
       flux(nz,n)=-0.5_WP*(ttf(nz-1,n)+ttf(nz,n))*W(nz,n)*area(nz,n)-flux(nz,n)
       
       !_______________________________________________________________________
       ! vert. flux at bottom layer --> zero bottom flux
       nz=nzmax
       flux(nz,n)= 0.0_WP-flux(nz,n)
       
       !_______________________________________________________________________
       ! Be carefull have to do vertical tracer advection here on old vertical grid
       ! also horizontal advection is done on old mesh (see helem contains old 
       ! mesh information)
       !_______________________________________________________________________
       ! vert. flux at remaining levels    
       do nz=nzmin+2,nzmax-2
            !centered (4th order)
            qc=(ttf(nz-1,n)-ttf(nz  ,n))/(Z_3d_n(nz-1,n)-Z_3d_n(nz  ,n))
            qu=(ttf(nz  ,n)-ttf(nz+1,n))/(Z_3d_n(nz  ,n)-Z_3d_n(nz+1,n))    
            qd=(ttf(nz-2,n)-ttf(nz-1,n))/(Z_3d_n(nz-2,n)-Z_3d_n(nz-1,n))
                    
            Tmean1=ttf(nz  ,n)+(2*qc+qu)*(zbar_3d_n(nz,n)-Z_3d_n(nz  ,n))/3.0_WP
            Tmean2=ttf(nz-1,n)+(2*qc+qd)*(zbar_3d_n(nz,n)-Z_3d_n(nz-1,n))/3.0_WP
            Tmean =(W(nz,n)+abs(W(nz,n)))*(Tmean1**do_Xmoment)+(W(nz,n)-abs(W(nz,n)))*(Tmean2**do_Xmoment)
    !         flux(nz,n)=-0.5_WP*(num_ord*(Tmean1+Tmean2)*W(nz,n)+(1.0_WP-num_ord)*Tmean)*area(nz,n)-flux(nz,n)
            flux(nz,n)=(-0.5_WP*(1.0_WP-num_ord)*Tmean - num_ord*((0.5_WP*(Tmean1+Tmean2))**do_Xmoment)*W(nz,n))*area(nz,n)-flux(nz,n)
       end do
    end do
end subroutine adv_tra_ver_qr4c
!
!
!===============================================================================
subroutine adv_tra_vert_ppm(ttf, w, do_Xmoment, mesh, flux, init_zero)
    use g_config
    use MOD_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_forcing_arrays
    implicit none
    type(t_mesh),  intent(in) , target :: mesh
    integer,       intent(in)          :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP), intent(in)          :: ttf (mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)          :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)       :: flux(mesh%nl,   myDim_nod2D)
    logical, optional                  :: init_zero
    real(kind=WP)                      :: tvert(mesh%nl), tv(mesh%nl), aL, aR, aj, x
    real(kind=WP)                      :: dzjm1, dzj, dzjp1, dzjp2, deltaj, deltajp1
    integer                            :: n, nz, nzmax, nzmin
    integer                            :: overshoot_counter, counter

#include "associate_mesh.h"

    if (present(init_zero))then
       if (init_zero) flux=0.0_WP
    else
       flux=0.0_WP
    end if

    ! --------------------------------------------------------------------------
    ! Vertical advection
    ! --------------------------------------------------------------------------
    ! A piecewise parabolic scheme for uniformly-spaced layers.
    ! See Colella and Woodward, JCP, 1984, 174-201. It can be coded so as to to take 
    ! non-uniformity into account, but this is more cumbersome. This is the version for AB
    ! time stepping
    ! --------------------------------------------------------------------------
    overshoot_counter=0
    counter          =0
    do n=1, myDim_nod2D
        !_______________________________________________________________________
        !Interpolate to zbar...depth levels --> all quantities (tracer ...) are 
        ! calculated on mid depth levels 
        ! nzmax ... number of depth levels at node n
        nzmax=nlevels_nod2D(n)
        nzmin=ulevels_nod2D(n)
        
        ! tracer at surface level
        tv(nzmin)=ttf(nzmin,n)        
        ! tracer at surface+1 level
!       tv(2)=-ttf(1,n)*min(sign(1.0, W(2,n)), 0._WP)+ttf(2,n)*max(sign(1.0, W(2,n)), 0._WP)
!       tv(3)=-ttf(2,n)*min(sign(1.0, W(3,n)), 0._WP)+ttf(3,n)*max(sign(1.0, W(3,n)), 0._WP)
        tv(nzmin+1)=0.5*(ttf(nzmin,  n)+ttf(nzmin+1,n))
        ! tacer at bottom-1 level
        tv(nzmax-1)=-ttf(nzmax-2,n)*min(sign(1.0, W(nzmax-1,n)), 0._WP)+ttf(nzmax-1,n)*max(sign(1.0, W(nzmax-1,n)), 0._WP)
!       tv(nzmax-1)=0.5_WP*(ttf(nzmax-2,n)+ttf(nzmax-1,n))
        ! tracer at bottom level
        tv(nzmax)=ttf(nzmax-1,n)
        
        !_______________________________________________________________________
        ! calc tracer for surface+2 until depth-2 layer
        ! see Colella and Woodward, JCP, 1984, 174-201 --> equation (1.9)
        ! loop over layers (segments)
        !!PS do nz=3, nzmax-3
        do nz=nzmin+1, nzmax-3
            !___________________________________________________________________
            ! for uniform spaced vertical grids --> piecewise parabolic method (ppm)
            ! equation (1.9)
            ! tv(nz)=(7.0_WP*(ttf(nz-1,n)+ttf(nz,n))-(ttf(nz-2,n)+ttf(nz+1,n)))/12.0_WP
            
            !___________________________________________________________________
            ! for non-uniformity spaced vertical grids --> piecewise parabolic 
            ! method (ppm) see see Colella and Woodward, JCP, 1984, 174-201 
            ! --> full equation (1.6), (1.7) and (1.8)
            dzjm1    = hnode_new(nz-1,n)
            dzj      = hnode_new(nz  ,n)
            dzjp1    = hnode_new(nz+1,n)
            dzjp2    = hnode_new(nz+2,n)
            ! Be carefull here vertical operation have to be done on NEW vertical mesh !!!
            
            !___________________________________________________________________
            ! equation (1.7)
            ! --> Here deltaj is the average slope in the jth zone of the parabola 
            !     with zone averages a_(j-1) and a_j, a_(j+1)
            ! --> a_j^n
            deltaj   = dzj/(dzjm1+dzj+dzjp1)* &
                      ( &
                       (2._WP*dzjm1+dzj    )/(dzjp1+dzj)*(ttf(nz+1,n)-ttf(nz  ,n)) +  &
                       (dzj    +2._WP*dzjp1)/(dzjm1+dzj)*(ttf(nz  ,n)-ttf(nz-1,n)) &
                      )
            ! --> a_(j+1)^n          
            deltajp1 = dzjp1/(dzj+dzjp1+dzjp2)* &
                      ( &
                       (2._WP*dzj+dzjp1  )/(dzjp2+dzjp1)*(ttf(nz+2,n)-ttf(nz+1,n)) +  &
                       (dzjp1+2._WP*dzjp2)/(dzj  +dzjp1)*(ttf(nz+1,n)-ttf(nz  ,n)) &
                      )
            !___________________________________________________________________
            ! condition (1.8)
            ! --> This modification leads to a somewhat steeper representation of 
            !     discontinuities in the solution. It also guarantees that a_(j+0.5)
            !     lies in the range of values defined by a_j; and a_(j+1);
            if ( (ttf(nz+1,n)-ttf(nz  ,n))*(ttf(nz  ,n)-ttf(nz-1,n)) > 0._WP ) then
                deltaj = min(  abs(deltaj), &
                             2._WP*abs(ttf(nz+1,n)-ttf(nz  ,n)),&
                             2._WP*abs(ttf(nz  ,n)-ttf(nz-1,n)) &
                             )*sign(1.0_WP,deltaj)
            else
                deltaj = 0.0_WP
            endif
            if ( (ttf(nz+2,n)-ttf(nz+1,n))*(ttf(nz+1,n)-ttf(nz  ,n)) > 0._WP ) then
                deltajp1 = min(  abs(deltajp1),&
                               2._WP*abs(ttf(nz+2,n)-ttf(nz+1,n)),&
                               2._WP*abs(ttf(nz+1,n)-ttf(nz,n)) &
                               )*sign(1.0_WP,deltajp1)
            else
                deltajp1 = 0.0_WP
            endif
            !___________________________________________________________________
            ! equation (1.6)
            ! --> calcualte a_(j+0.5)
            ! nz+1 is the interface betweel layers (segments) nz and nz+1
            tv(nz+1)=    ttf(nz,n) &
                        + dzj/(dzj+dzjp1)*(ttf(nz+1,n)-ttf(nz,n)) &
                        + 1._WP/(dzjm1+dzj+dzjp1+dzjp2) * &
                        ( &
                            (2._WP*dzjp1*dzj)/(dzj+dzjp1)* &
                                ((dzjm1+dzj)/(2._WP*dzj+dzjp1) - (dzjp2+dzjp1)/(2._WP*dzjp1+dzj))*(ttf(nz+1,n)-ttf(nz,n)) &
                        - dzj*(dzjm1+dzj)/(2._WP*dzj+dzjp1)*deltajp1 &
                        + dzjp1*(dzjp1+dzjp2)/(dzj+2._WP*dzjp1)*deltaj &
                        )
                       !tv(nz+1)=max(min(ttf(nz, n), ttf(nz+1, n)), min(max(ttf(nz, n), ttf(nz+1, n)), tv(nz+1)))
        end do ! --> do nz=2,nzmax-3
        
        tvert(1:nzmax)=0._WP
        ! loop over layers (segments)
        do nz=nzmin, nzmax-1
            if ((W(nz,n)<=0._WP) .AND. (W(nz+1,n)>=0._WP)) CYCLE
            counter=counter+1
            aL=tv(nz)
            aR=tv(nz+1)
            if ((aR-ttf(nz, n))*(ttf(nz, n)-aL)<=0._WP) then
                !   write(*,*) aL, ttf(nz, n), aR
                overshoot_counter=overshoot_counter+1
                aL =ttf(nz, n)
                aR =ttf(nz, n)
            end if
            if ((aR-aL)*(ttf(nz, n)-0.5_WP*(aL+aR))> (aR-aL)**2/6._WP) then
                aL =3._WP*ttf(nz, n)-2._WP*aR
            end if
            if ((aR-aL)*(ttf(nz, n)-0.5_WP*(aR+aL))<-(aR-aL)**2/6._WP) then
                aR =3._WP*ttf(nz, n)-2._WP*aL
            end if
            
            dzj   = hnode(nz,n)
            aj=6.0_WP*(ttf(nz, n)-0.5_WP*(aL+aR))
            
            if (W(nz,n)>0._WP) then
                x=min(W(nz,n)*dt/dzj, 1._WP)
                tvert(nz  )=(-aL-0.5_WP*x*(aR-aL+(1._WP-2._WP/3._WP*x)*aj))
                tvert(nz  )=( tvert(nz)**do_Xmoment ) ! compute 2nd moment for DVD
                tvert(nz  )=tvert(nz)*area(nz,n)*W(nz,n)
            end if
            
            if (W(nz+1,n)<0._WP) then
                x=min(-W(nz+1,n)*dt/dzj, 1._WP)
                tvert(nz+1)=(-aR+0.5_WP*x*(aR-aL-(1._WP-2._WP/3._WP*x)*aj))
                tvert(nz+1)=( tvert(nz+1)**do_Xmoment ) ! compute 2nd moment for DVD
                tvert(nz+1)=tvert(nz+1)*area(nz+1,n)*W(nz+1,n)
            end if
        end do
        
        !_______________________________________________________________________
        ! Surface flux
        tvert(nzmin)= -( tv(nzmin)**do_Xmoment )*W(nzmin,n)*area(nzmin,n)
        ! Zero bottom flux
        tvert(nzmax)=0.0_WP        
        flux(nzmin:nzmax, n)=tvert(nzmin:nzmax)-flux(nzmin:nzmax, n)
    end do ! --> do n=1, myDim_nod2D
!       if (mype==0) write(*,*) 'PPM overshoot statistics:', real(overshoot_counter)/real(counter)
end subroutine adv_tra_vert_ppm
!
!
!===============================================================================
subroutine adv_tra_ver_cdiff(ttf, w, do_Xmoment, mesh, flux, init_zero)
    use g_config
    use MOD_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_forcing_arrays
    implicit none
    type(t_mesh),  intent(in), target :: mesh
    integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: W  (mesh%nl,   myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: flux(mesh%nl,  myDim_nod2D)
    logical, optional                 :: init_zero
    integer                           :: n, nz, nzmax, nzmin
    real(kind=WP)                     :: tvert(mesh%nl), tv
#include "associate_mesh.h"

    if (present(init_zero))then
       if (init_zero) flux=0.0_WP
    else
       flux=0.0_WP
    end if

    do n=1, myDim_nod2D
        !_______________________________________________________________________
        nzmax=nlevels_nod2D(n)-1
        nzmin=ulevels_nod2D(n)
        
        !_______________________________________________________________________
        ! Surface flux
        tvert(nzmin)= -W(nzmin,n)*(ttf(nzmin,n)**do_Xmoment)*area(nzmin,n)
        
        !_______________________________________________________________________
        ! Zero bottom flux
        tvert(nzmax+1)=0.0_WP        
        
        !_______________________________________________________________________
        ! Other levels
        do nz=nzmin+1, nzmax
            tv=0.5_WP*(ttf(nz-1,n)+ttf(nz,n))
            tv=tv**do_Xmoment
            tvert(nz)= -tv*W(nz,n)*area(nz,n)
        end do
        
        !_______________________________________________________________________
        flux(nzmin:nzmax, n)=tvert(nzmin:nzmax)-flux(nzmin:nzmax, n)
    end do ! --> do n=1, myDim_nod2D
end subroutine adv_tra_ver_cdiff
