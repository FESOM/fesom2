module oce_adv_tra_driver_interfaces
  interface
   subroutine do_oce_adv_tra(ttf, ttfAB, vel, w, wi, we, do_Xmoment, dttf_h, dttf_v, opth, optv, mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, myDim_elem2D+eDim_elem2D)
      real(kind=WP), intent(in), target :: W(mesh%nl,   myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in), target :: WI(mesh%nl,   myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in), target :: WE(mesh%nl,   myDim_nod2D+eDim_nod2D)
      integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
      real(kind=WP), intent(in)         :: ttf   (mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: ttfAB (mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: dttf_h(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: dttf_v(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: opth, optv
    end subroutine
  end interface
end module

module oce_tra_adv_flux2dtracer_interface
  interface
    subroutine oce_tra_adv_flux2dtracer(dttf_h, dttf_v, flux_h, flux_v, mesh, use_lo, ttf, lo)
    !update the solution for vertical and horizontal flux contributions
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP), intent(inout)      :: dttf_h(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: dttf_v(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: flux_h(mesh%nl-1, myDim_edge2D)
      real(kind=WP), intent(inout)      :: flux_v(mesh%nl,  myDim_nod2D)
      logical,       optional           :: use_lo
      real(kind=WP), optional           :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), optional           :: lo (mesh%nl-1, myDim_nod2D+eDim_nod2D)
    end subroutine
  end interface
end module
!
!
!===============================================================================
subroutine do_oce_adv_tra(ttf, ttfAB, vel, w, wi, we, do_Xmoment, dttf_h, dttf_v, opth, optv, mesh)
    use MOD_MESH
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_CONFIG
    use g_comm_auto
    use oce_adv_tra_hor_interfaces
    use oce_adv_tra_ver_interfaces
    use oce_adv_tra_fct_interfaces
    use oce_tra_adv_flux2dtracer_interface
    implicit none
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, myDim_elem2D+eDim_elem2D)
    real(kind=WP), intent(in), target :: W(mesh%nl,   myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in), target :: WI(mesh%nl,   myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in), target :: WE(mesh%nl,   myDim_nod2D+eDim_nod2D)
    integer,       intent(in)         :: do_Xmoment !--> = [1,2] compute 1st & 2nd moment of tracer transport
    real(kind=WP), intent(in)         :: ttf  (mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: ttfAB(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: dttf_h(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: dttf_v(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: opth, optv
    real(kind=WP), pointer, dimension (:,:) :: pwvel
 
    integer       :: el(2), enodes(2), nz, n, e
    integer       :: nl12, nu12, nl1, nl2, nu1, nu2, tr_num
    real(kind=WP) :: cLO, cHO, deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=WP) :: qc, qu, qd
    real(kind=WP) :: tvert(mesh%nl), tvert_e(mesh%nl), a, b, c, d, da, db, dg, vflux, Tupw1
    real(kind=WP) :: Tmean, Tmean1, Tmean2, num_ord
    logical       :: do_zero_flux

#include "associate_mesh.h"
    !___________________________________________________________________________
    ! compute FCT horzontal and vertical low order solution as well as lw order 
    ! part of antidiffusive flux
    if (trim(tra_adv_lim)=='FCT') then 
        ! compute the low order upwind horizontal flux
        ! init_zero=.true.  : zero the horizontal flux before computation
        ! init_zero=.false. : input flux will be substracted
        call adv_tra_hor_upw1(ttf, vel, do_Xmoment, mesh, adv_flux_hor, init_zero=.true.)
        
        ! update the LO solution for horizontal contribution
        fct_LO=0.0_WP
        do e=1, myDim_edge2D
            enodes=edges(:,e)
            el=edge_tri(:,e)        
            nl1=nlevels(el(1))-1
            nu1=ulevels(el(1))
            nl2=0
            nu2=0
            if(el(2)>0) then
                nl2=nlevels(el(2))-1
                nu2=ulevels(el(2))
            end if     
            
            nl12 = max(nl1,nl2)
            nu12 = nu1
            if (nu2>0) nu12 = min(nu1,nu2)
            
            !!PS do  nz=1, max(nl1, nl2)
            do nz=nu12, nl12
                fct_LO(nz, enodes(1))=fct_LO(nz, enodes(1))+adv_flux_hor(nz, e)
                fct_LO(nz, enodes(2))=fct_LO(nz, enodes(2))-adv_flux_hor(nz, e)
            end do
        end do
        
        ! compute the low order upwind vertical flux (explicit part only)
        ! zero the input/output flux before computation
        call adv_tra_ver_upw1(ttf, we, do_Xmoment, mesh, adv_flux_ver, init_zero=.true.)
        
        ! update the LO solution for vertical contribution
        do n=1, myDim_nod2D
            nu1 = ulevels_nod2D(n)
            nl1 = nlevels_nod2D(n)
            !!PS do  nz=1, nlevels_nod2D(n)-1
            do  nz= nu1, nl1-1
                fct_LO(nz,n)=(ttf(nz,n)*hnode(nz,n)+(fct_LO(nz,n)+(adv_flux_ver(nz, n)-adv_flux_ver(nz+1, n)))*dt/areasvol(nz,n))/hnode_new(nz,n)
            end do
        end do
        
        if (w_split) then !wvel/=wvel_e
            ! update for implicit contribution (w_split option)
            call adv_tra_vert_impl(fct_LO, wi, mesh)
            ! compute the low order upwind vertical flux (full vertical velocity)
            ! zero the input/output flux before computation
            ! --> compute here low order part of vertical anti diffusive fluxes, 
            !     has to be done on the full vertical velocity w
            call adv_tra_ver_upw1(ttf, w, do_Xmoment, mesh, adv_flux_ver, init_zero=.true.)
        end if
        
        call exchange_nod(fct_LO)
    end if

    do_zero_flux=.true.
    if (trim(tra_adv_lim)=='FCT') do_zero_flux=.false.
   
    !___________________________________________________________________________
    ! do horizontal tracer advection, in case of FCT high order solution 
    SELECT CASE(trim(tra_adv_hor))
        CASE('MUSCL')
            ! compute the untidiffusive horizontal flux (init_zero=.false.: input is the LO horizontal flux computed above)
            call adv_tra_hor_muscl(ttfAB, uv,   do_Xmoment, mesh, opth,  adv_flux_hor, init_zero=do_zero_flux)
        CASE('MFCT')
             call adv_tra_hor_mfct(ttfAB, uv,   do_Xmoment, mesh, opth,  adv_flux_hor, init_zero=do_zero_flux)
        CASE('UPW1')
             call adv_tra_hor_upw1(ttfAB, uv,   do_Xmoment, mesh,        adv_flux_hor, init_zero=do_zero_flux)
        CASE DEFAULT !unknown
            if (mype==0) write(*,*) 'Unknown horizontal advection type ',  trim(tra_adv_hor), '! Check your namelists!'
            call par_ex(1)
    END SELECT
   
    if (trim(tra_adv_lim)=='FCT') then
       pwvel=>w
    else
       pwvel=>we
    end if
 
    !___________________________________________________________________________
    ! do vertical tracer advection, in case of FCT high order solution 
    SELECT CASE(trim(tra_adv_ver))
        CASE('QR4C')
            ! compute the untidiffusive vertical flux   (init_zero=.false.:input is the LO vertical flux computed above)
            call adv_tra_ver_qr4c (ttfAB, pwvel,   do_Xmoment, mesh, optv, adv_flux_ver, init_zero=do_zero_flux)
        CASE('CDIFF')
            call adv_tra_ver_cdiff(ttfAB, pwvel,   do_Xmoment, mesh,       adv_flux_ver, init_zero=do_zero_flux)
        CASE('PPM')
            call adv_tra_vert_ppm (ttfAB, pwvel,   do_Xmoment, mesh,       adv_flux_ver, init_zero=do_zero_flux)
        CASE('UPW1')
        call adv_tra_ver_upw1 (ttfAB, pwvel,   do_Xmoment, mesh,       adv_flux_ver, init_zero=do_zero_flux)
        CASE DEFAULT !unknown
            if (mype==0) write(*,*) 'Unknown vertical advection type ',  trim(tra_adv_ver), '! Check your namelists!'
            call par_ex(1)
        ! --> be aware the vertical implicite part in case without FCT is done in 
        !     oce_ale_tracer.F90 --> subroutine diff_ver_part_impl_ale(tr_num, mesh)
        !     for do_wimpl=.true.
    END SELECT
    
    !___________________________________________________________________________
    !
!if (mype==0) then
!   write(*,*) 'check new:'
!   write(*,*) '1:', minval(fct_LO),       maxval(fct_LO),       sum(fct_LO)
!   write(*,*) '2:', minval(adv_flux_hor), maxval(adv_flux_hor), sum(adv_flux_hor)
!   write(*,*) '3:', minval(adv_flux_ver), maxval(adv_flux_ver), sum(adv_flux_ver)
!end if
    if (trim(tra_adv_lim)=='FCT') then
!if (mype==0) write(*,*) 'before:', sum(abs(adv_flux_ver)), sum(abs(adv_flux_hor))
       call oce_tra_adv_fct(dttf_h, dttf_v, ttf, fct_LO, adv_flux_hor, adv_flux_ver, mesh)
!if (mype==0) write(*,*) 'after:', sum(abs(adv_flux_ver)), sum(abs(adv_flux_hor))
       call oce_tra_adv_flux2dtracer(dttf_h, dttf_v, adv_flux_hor, adv_flux_ver, mesh, use_lo=.TRUE., ttf=ttf, lo=fct_LO)
    else
       call oce_tra_adv_flux2dtracer(dttf_h, dttf_v, adv_flux_hor, adv_flux_ver, mesh)
    end if
end subroutine do_oce_adv_tra
!
!
!===============================================================================
subroutine oce_tra_adv_flux2dtracer(dttf_h, dttf_v, flux_h, flux_v, mesh, use_lo, ttf, lo)
    use MOD_MESH
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_CONFIG
    use g_comm_auto
    implicit none
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(inout)      :: dttf_h(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: dttf_v(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: flux_h(mesh%nl-1, myDim_edge2D)
    real(kind=WP), intent(inout)      :: flux_v(mesh%nl,  myDim_nod2D)
    logical,       optional           :: use_lo
    real(kind=WP), optional           :: lo (mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), optional           :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    integer                           :: n, nz, k, elem, enodes(3), num, el(2), nu12, nl12, nu1, nu2, nl1, nl2, edge
#include "associate_mesh.h"
    !___________________________________________________________________________
    ! c. Update the solution
    ! Vertical
    if (present(use_lo)) then
       if (use_lo) then
          do n=1, myDim_nod2d
             nu1 = ulevels_nod2D(n)
             nl1 = nlevels_nod2D(n)
             !!PS do nz=1,nlevels_nod2D(n)-1
             do nz=nu1, nl1-1  
                dttf_v(nz,n)=dttf_v(nz,n)-ttf(nz,n)*hnode(nz,n)+LO(nz,n)*hnode_new(nz,n)
             end do
           end do
       end if
    end if

    do n=1, myDim_nod2d
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        do nz=nu1,nl1-1  
            dttf_v(nz,n)=dttf_v(nz,n) + (flux_v(nz,n)-flux_v(nz+1,n))*dt/areasvol(nz,n)
        end do
    end do

    
    ! Horizontal
    do edge=1, myDim_edge2D
        enodes(1:2)=edges(:,edge)
        el=edge_tri(:,edge)
        nl1=nlevels(el(1))-1
        nu1=ulevels(el(1))
        
        nl2=0
        nu2=0
        if(el(2)>0) then
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
        end if 
        
        nl12 = max(nl1,nl2)
        nu12 = nu1
        if (nu2>0) nu12 = min(nu1,nu2)
            
        !!PS do  nz=1, max(nl1, nl2)
        do nz=nu12, nl12
            dttf_h(nz,enodes(1))=dttf_h(nz,enodes(1))+flux_h(nz,edge)*dt/areasvol(nz,enodes(1))
            dttf_h(nz,enodes(2))=dttf_h(nz,enodes(2))-flux_h(nz,edge)*dt/areasvol(nz,enodes(2))
        end do
    end do
end subroutine oce_tra_adv_flux2dtracer
