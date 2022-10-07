!---------------------------------------------------------------------------
!Implementation of Gent & McWiliams parameterization after R. Ferrari et al., 2010
!Contains:
!  fer_solve_Gamma
!  fer_gamma2vel
!  fer_compute_C_K ! this subroutine shall be a subject of future tuning (with respect to fer_k)
!===========================================================================
subroutine fer_solve_Gamma(mesh)
    USE MOD_MESH
    USE o_PARAM
    USE o_ARRAYS, ONLY: sigma_xy, fer_gamma, bvfreq, fer_c, fer_K, zbar_n, Z_n, hnode_new, zbar_n_bot
    USE g_PARSUP
    USE g_CONFIG
    use g_comm_auto
    IMPLICIT NONE
    type(t_mesh), intent(in)               , target :: mesh	
    integer                                :: nz, n, nzmax, nzmin
    real(kind=WP)                          :: zinv1,zinv2, zinv, m, r
    real(kind=WP)                          :: a(mesh%nl), b(mesh%nl), c(mesh%nl)
    real(kind=WP)                          :: cp(mesh%nl), tp(2,mesh%nl)
    real(kind=WP), dimension(:,:), pointer :: tr

#include "associate_mesh.h"

    DO n=1,myDim_nod2D
        tr=>fer_gamma(:,:,n)
!       !_____________________________________________________________________
!       ! minimum number of levels below elements containing node n
!       nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
        nzmax=nlevels_nod2D(n)
        nzmin=ulevels_nod2D(n)
        
        !_______________________________________________________________________
        ! compute Z_n(:) and zbar_n(:) for the current ALE step
        zbar_n=0.0_WP
        Z_n   =0.0_WP
        zbar_n(nzmax)=zbar_n_bot(n)                                ! depth of the deepest level with respect to partial cell
        Z_n(nzmax-1) =zbar_n(nzmax) + hnode_new(nzmax-1,n)/2.0_WP  ! depth of the deepest layer
        !!PS do nz=nzmax-1,2,-1
        do nz=nzmax-1,nzmin+1,-1
            zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)            ! depth of the nz level
            Z_n(nz-1)  = zbar_n(nz)   + hnode_new(nz-1,n)/2.0_WP   ! depth of the nz layer
        end do
        !!PS zbar_n(1) = zbar_n(2) + hnode_new(1,n)                ! surface level height (/depth)
        zbar_n(nzmin) = zbar_n(nzmin+1) + hnode_new(nzmin,n)       ! surface level height (/depth)
        
        !_____________________________________________________________________
        ! minimum number of levels below elements containing node n
        !!PS nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
        !!PS nzmin=maxval(ulevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
        nzmax=nlevels_nod2D_min(n)
        nzmin=ulevels_nod2D_max(n)
        
        ! The first row
        !!PS c(1)=0.0_WP
        !!PS a(1)=0.0_WP
        !!PS b(1)=1.0_WP
        c(nzmin)=0.0_WP
        a(nzmin)=0.0_WP
        b(nzmin)=1.0_WP
        
        !!PS zinv2=1.0_WP/(zbar_n(1)-zbar_n(2))
        zinv2=1.0_WP/(zbar_n(nzmin)-zbar_n(nzmin+1))
        
        !!PS DO nz=2, nzmax-1
        DO nz=nzmin+1, nzmax-1
            zinv1=zinv2
            zinv2=1.0_WP/(zbar_n(nz)-zbar_n(nz+1))
            zinv =1.0_WP/(Z_n(nz-1)-Z_n(nz))
            a(nz)= fer_c(n)*zinv1*zinv
            c(nz)= fer_c(n)*zinv2*zinv
            b(nz)=-a(nz)-c(nz)-max(bvfreq(nz,n), 1.e-8)
        END DO
        
        ! The last row
        nz=nzmax
        c(nz)=0.0_WP
        a(nz)=0.0_WP
        b(nz)=1.0_WP
        
        ! ===========================================
        ! The rhs:
        !!PS tr(:, 1)=0.
        tr(:, nzmin)=0.0_WP
        tr(:, nzmax)=0.0_WP
        !!PS DO nz=2, nzmax-1
        DO nz=nzmin+1, nzmax-1
            r=g/density_0
            tr(1, nz)=r*0.5_WP*sum(sigma_xy(1,nz-1:nz,n))*fer_K(nz, n)
            tr(2, nz)=r*0.5_WP*sum(sigma_xy(2,nz-1:nz,n))*fer_K(nz, n)
        END DO
        
        ! =============================================
        ! The sweep algorithm
        ! initialize c-prime and s,t-prime
        !!PS cp(1) = c(1)/b(1)
        !!PS tp(:,1) = tr(:,1)/b(1)
        cp(nzmin) = c(nzmin)/b(nzmin)
        tp(:,nzmin) = tr(:,nzmin)/b(nzmin)
        
        ! solve for vectors c-prime and t, s-prime
        !!PS DO nz = 2, nzmax
        DO nz = nzmin+1, nzmax
            m = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            tp(:,nz) = (tr(:,nz)-tp(:,nz-1)*a(nz))/m
        END DO
        
        ! initialize x
        tr(:,nzmax) = tp(:,nzmax)
        
        ! solve for x from the vectors c-prime and d-prime
        !!PS do nz = nzmax-1, 1, -1
        do nz = nzmax-1, nzmin, -1
            tr(:,nz) = tp(:,nz)-cp(nz)*tr(:,nz+1)
        end do
    END DO   !!! cycle over nodes
    
    call exchange_nod(fer_gamma)
END subroutine fer_solve_Gamma
!
!
!
!====================================================================
subroutine fer_gamma2vel(mesh)
    USE MOD_MESH
    USE o_PARAM
    USE o_ARRAYS, ONLY: fer_gamma, fer_uv, helem
    USE g_PARSUP
    USE g_CONFIG
    use g_comm_auto
    IMPLICIT NONE

    integer                                :: nz, nzmax, el, elnod(3), nzmin
    real(kind=WP)                          :: zinv
    real(kind=WP)                          :: onethird=1._WP/3._WP
    type(t_mesh), intent(in)               , target :: mesh

#include  "associate_mesh.h"

    DO el=1, myDim_elem2D
        elnod=elem2D_nodes(:,el)
        ! max. number of levels at element el
        nzmax=nlevels(el)
        nzmin=ulevels(el)
        !!PS DO nz=1, nzmax-1
        DO nz=nzmin, nzmax-1
            zinv=onethird/helem(nz,el)
            fer_uv(1,nz,el)=sum(fer_gamma(1,nz,elnod)-fer_gamma(1,nz+1,elnod))*zinv
            fer_uv(2,nz,el)=sum(fer_gamma(2,nz,elnod)-fer_gamma(2,nz+1,elnod))*zinv
        END DO
    END DO
    call exchange_elem(fer_uv)
end subroutine fer_gamma2vel
!
!
!
!===============================================================================
subroutine init_Redi_GM(mesh) !fer_compute_C_K_Redi
    USE MOD_MESH
    USE o_PARAM
    USE o_ARRAYS, ONLY: fer_c, fer_k, fer_scal, Ki, bvfreq, MLD1_ind, neutral_slope, coriolis_node, hnode_new, Z_3d_n
    USE g_PARSUP
    USE g_CONFIG
    use g_comm_auto
    IMPLICIT NONE
    type(t_mesh), intent(in) , target :: mesh
    integer                  :: n, nz, nzmax, nzmin
    real(kind=WP)            :: reso, c1, rosb, scaling, rr_ratio, aux_zz(mesh%nl)
    real(kind=WP)            :: x0=1.5_WP, sigma=.15_WP ! Fermi function parameters to cut off GM where Rossby radius is resolved
    real(kind=WP)            :: c_min=0.5_WP, f_min=1.e-6_WP, r_max=200000._WP
    real(kind=WP)            :: zscaling(mesh%nl)
    real(kind=WP)            :: bvref

#include "associate_mesh.h"

    ! fill arrays for 3D Redi and GM coefficients: F1(xy)*F2(z)
    !******************************* F1(x,y) ***********************************
    do n=1, myDim_nod2D
        nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
        nzmin=maxval(ulevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
        reso=mesh_resolution(n)
        if (Fer_GM) then
            c1=0._wp
            
            !!PS do nz=1, nzmax-1
            do nz=nzmin, nzmax-1
                !!PS c1=c1+hnode_new(nz,n)*(sqrt(max(bvfreq(nz,n), 0._WP))+sqrt(max(bvfreq(nz+1,n), 0._WP)))/2._WP
                c1=c1+hnode_new(nz,n)*(sqrt(abs(max(bvfreq(nz,n), 0._WP)))+sqrt(abs(max(bvfreq(nz+1,n), 0._WP))))/2._WP ! add abs() for -0 case, cray
            end do
            c1=max(c_min, c1/pi) !ca. first baroclinic gravity wave speed limited from below by c_min
            scaling=1._WP
            
            !___________________________________________________________________
            ! Cutoff K_GM depending on (Resolution/Rossby radius) ratio
            if (scaling_Rossby) then
                rosb=min(c1/max(abs(coriolis_node(n)), f_min), r_max)
                rr_ratio=min(reso/rosb, 5._WP)
                scaling=1._WP/(1._WP+exp(-(rr_ratio-x0)/sigma))
            end if
            
            !___________________________________________________________________
            ! Scale K_GM with resolution (referenced to 100,000m)
            if (scaling_resolution) then
                scaling=scaling*(reso/100000._WP)**K_GM_resscalorder !put to repo
            end if
            
            !___________________________________________________________________
            ! resolution ramp function for the switch off of GM 
            ! default: 
            !     ^GM_scaling
            !   1-|           .-----------
            !     |         ./ |
            !     |       ./   |
            !   0-|------/-----|-------------->Resolution
            !            |     |
            !           30km   40km (FESOM1.4/2.0 default, MPAS: 20km...30km )
!!PS             if (reso < 40000.0_WP) then
!!PS                 scaling=scaling*max((reso/10000.0_WP-3.0_WP), 0._WP) !no GM below 30km resolution
!!PS             end if
            if (reso/1000.0_WP < K_GM_rampmax) then
                scaling=scaling*max((reso/1000.0_WP-K_GM_rampmin)/(K_GM_rampmax-K_GM_rampmin), 0._WP) !no GM below 30km resolution
            end if
            
            !___________________________________________________________________
            ! apply KGM scaling paramter (K_GM_scal)
            fer_scal(n) = min(scaling,1.0_WP)
             ! set maximum amplitude to K_GM_max
            !!PS fer_k(1,n)  = fer_scal(n)*K_GM_max
            !!PS ! limit lower values to K_GM_min
            !!PS fer_k(1,n)  = max(fer_k(1,n),K_GM_min)
            fer_k(nzmin,n)  = fer_scal(n)*K_GM_max
            ! limit lower values to K_GM_min
            fer_k(nzmin,n)  = max(fer_k(nzmin,n),K_GM_min)
            fer_c(n)    = c1*c1                          !put to repo
        end if
        
        !_______________________________________________________________________
        ! note, Redi can be used without GM and vise versa!
        ! if both are used it will be reset below
        if (Redi) then
            !!PS Ki(1,n)=K_hor*(reso/100000.0_WP)**2
            Ki(nzmin,n)=K_hor*(reso/100000.0_WP)**2
        end if
    end do
 
    !Like in FESOM 1.4 we make Redi equal GM
    if (Redi .and. Fer_GM) then
        !!PS Ki(1,:)=fer_k(1,:)
        Ki(nzmin,:)=fer_k(nzmin,:)
    end if
   
!******************************* F2(z) (e.g. Ferreira et al., 2005) *********************************
!Ferreira, D., Marshall, J. and Heimbach, P.: Estimating Eddy Stresses by Fitting Dynamics to Observations Using a
!Residual-Mean Ocean Circulation Model and Its Adjoint, Journal of Physical Oceanography, 35(10), 1891–
!1910, doi:10.1175/jpo2785.1, 2005.

    do n=1,myDim_nod2D
        nzmax=nlevels_nod2D(n)
        nzmin=ulevels_nod2D(n)
        !_______________________________________________________________________
        ! Allpy vertical scaling after Ferreira et al.(2005)
        if (scaling_Ferreira) then
            !___________________________________________________________________
            ! choose reference buoyancy
            if (K_GM_bvref==0) then
                ! ferreira bvref value surface (original Ferreira)
                !!PS bvref=max(bvfreq(1, n), 1.e-6_WP)
                bvref=max(bvfreq(nzmin, n), 1.e-6_WP)
            elseif (K_GM_bvref==1) then
                ! ferreira bvref value bottom mixed layer (Dima)
                bvref=max(bvfreq(MLD1_ind(n)+1, n), 1.e-6_WP)
            elseif (K_GM_bvref==2) then
                ! ferreira bvref value mean over mixed layer
                !!PS bvref=max(sum(bvfreq(1:MLD1_ind(n), n))/(MLD1_ind(n)), 1.e-6_WP)
                bvref=max(sum(bvfreq(nzmin:MLD1_ind(n), n))/(MLD1_ind(n)), 1.e-6_WP)
            elseif (K_GM_bvref==3) then
                ! ferreira bvref value depth weighted mean over mixed layer
                aux_zz=0.0_WP
                !!PS aux_zz(1:MLD1_ind(n)-1) = Z_3d_n(2:MLD1_ind(n),n)-Z_3d_n(1:MLD1_ind(n)-1,n)
                !!PS bvref=max(sum(bvfreq(2:MLD1_ind(n),n)*aux_zz(1:MLD1_ind(n)-1))/sum(aux_zz(1:MLD1_ind(n)-1)), 1.e-6_WP)
                aux_zz(nzmin:MLD1_ind(n)-1) = Z_3d_n(nzmin+1:MLD1_ind(n),n)-Z_3d_n(nzmin:MLD1_ind(n)-1,n)
                bvref=max(sum(bvfreq(nzmin+1:MLD1_ind(n),n)*aux_zz(nzmin:MLD1_ind(n)-1))/sum(aux_zz(nzmin:MLD1_ind(n)-1)), 1.e-6_WP)    
            end if 
            
            !___________________________________________________________________
            ! compute scaling with respect to reference buoyancy
            !!PS do nz=1, nzmax
            do nz=nzmin, nzmax
                zscaling(nz)=max(bvfreq(nz, n)/bvref, 0.2_WP)
                zscaling(nz)=min(zscaling(nz), 1.0_WP)
            end do
        else
            zscaling=1.0_WP
        end if
        
        !_______________________________________________________________________
        ! Switch off GM and Redi within a BL in NH (a strategy following FESOM 1.4)
        if (scaling_FESOM14) then
            !zscaling(1:MLD1_ind(n)+1)=0.0_WP
            !!PS do nz=1, nzmax
            do nz=nzmin, nzmax
                if (neutral_slope(3, min(nz, nl-1), n) > 5.e-3_WP) zscaling(nz)=0.0_WP
            end do
        end if
      
        !_______________________________________________________________________
        ! do vertical Ferreira scaling and limiting for GM diffusivity
        if (Fer_GM) then
            ! start with index 2 to not alternate fer_k(1,n) which contains here 
            ! the surface template for the scaling 
            !!PS do nz=2, nzmax
            do nz=nzmin+1, nzmax
                fer_k(nz,n)=fer_k(1,n)*zscaling(nz)
            end do 
            ! after vertical Ferreira scaling is done also scale surface template
            !!PS fer_k(1,n)=fer_k(1,n)*zscaling(1)
            fer_k(nzmin,n)=fer_k(nzmin,n)*zscaling(nzmin)
        end if
        
        !_______________________________________________________________________
        ! do vertical Ferreira scaling and limiting for Redi diffusivity
        if (Redi) then
            ! start with index 2 to not alternate fer_k(1,n) which contains here 
            ! the surface template for the scaling 
            !!PS do nz=2, nzmax-1
            do nz=nzmin+1, nzmax-1
                !!PS Ki(nz,n)= Ki(1,n)*0.5_WP*(zscaling(nz)+zscaling(nz+1))
                Ki(nz,n)= Ki(nzmin,n)*0.5_WP*(zscaling(nz)+zscaling(nz+1))
            end do
            ! after vertical Ferreira scaling is done also scale surface template
            !!PS Ki(1,n)=Ki(1,n)*0.5_WP*(zscaling(1)+zscaling(2))
            Ki(nzmin,n)=Ki(nzmin,n)*0.5_WP*(zscaling(nzmin)+zscaling(nzmin+1))
        end if
   end do

   if (Fer_GM) call exchange_nod(fer_c)
   if (Fer_GM) call exchange_nod(fer_k)
   if (Redi)   call exchange_nod(Ki)
end subroutine init_Redi_GM
!====================================================================
