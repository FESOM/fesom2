module densityJM_components_interface
  interface
    subroutine densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)
      USE MOD_MESH
      real(kind=WP), intent(IN)  :: t,s
      real(kind=WP), intent(OUT) :: bulk_0, bulk_pz, bulk_pz2, rhopot
      type(t_mesh), intent(in) , target :: mesh
    end subroutine
  end interface
end module

module density_linear_interface
  interface
    subroutine density_linear(t, s, bulk_0, bulk_pz, bulk_pz2, rho_out, mesh)
      USE MOD_MESH
      real(kind=WP), intent(IN)  :: t,s
      real(kind=WP), intent(OUT) :: bulk_0, bulk_pz, bulk_pz2, rho_out
      type(t_mesh), intent(in) , target :: mesh
    end subroutine
  end interface
end module

module pressure_force_4_linfs_fullcell_interface
  interface
    subroutine pressure_force_4_linfs_fullcell(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module pressure_force_4_linfs_nemo_interface
  interface
    subroutine pressure_force_4_linfs_nemo(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module pressure_force_4_linfs_shchepetkin_interface
  interface
    subroutine pressure_force_4_linfs_shchepetkin(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module pressure_force_4_linfs_easypgf_interface
  interface
    subroutine pressure_force_4_linfs_easypgf(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module pressure_force_4_linfs_cubicspline_interface
  interface
    subroutine pressure_force_4_linfs_cubicspline(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module pressure_force_4_linfs_cavity_interface
  interface
    subroutine pressure_force_4_linfs_cavity(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module pressure_force_4_zxxxx_shchepetkin_interface
  interface
    subroutine pressure_force_4_zxxxx_shchepetkin(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module pressure_force_4_zxxxx_easypgf_interface
  interface
    subroutine pressure_force_4_zxxxx_easypgf(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module pressure_force_4_zxxxx_cubicspline_interface
  interface
    subroutine pressure_force_4_zxxxx_cubicspline(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module init_ref_density_interface
  interface
    subroutine init_ref_density(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
!
!
!===============================================================================
subroutine pressure_bv(mesh)
! fill in the hydrostatic pressure and the Brunt-Vaisala frequency 
! in a single pass the using split form of the equation of state
! as proposed by NR
    use g_config
    USE o_PARAM
    USE MOD_MESH
    USE o_ARRAYS
    USE g_PARSUP
    use i_arrays
    USE o_mixing_KPP_mod, only: dbsfc
    USE diagnostics,      only: ldiag_dMOC
    use densityJM_components_interface
    use density_linear_interface
    IMPLICIT NONE
    type(t_mesh), intent(in) , target :: mesh    
    real(kind=WP)            :: dz_inv, bv,  a, rho_up, rho_dn, t, s
    integer                  :: node, nz, nl1, nzmax, nzmin
    real(kind=WP)            :: rhopot(mesh%nl), bulk_0(mesh%nl), bulk_pz(mesh%nl), bulk_pz2(mesh%nl), rho(mesh%nl), dbsfc1(mesh%nl), db_max
    real(kind=WP)            :: bulk_up, bulk_dn, smallvalue, buoyancy_crit, rho_surf, aux_rho, aux_rho1
    real(kind=WP)            :: sigma_theta_crit= 0.03 !Boyer-Montegut (2004) !0.125_WP   !kg/m3, Levitus threshold for computing MLD2
    logical                  :: flag1, flag2, mixing_kpp
#include "associate_mesh.h"
    smallvalue=1.0e-20
    buoyancy_crit=0.0003_WP
    mixing_kpp = (mix_scheme_nmb==1 .or. mix_scheme_nmb==17)  ! NR Evaluate string comparison outside the loop. It is expensive.
!!PS     mixing_kpp = (trim(mix_scheme)=='KPP' .or. trim(mix_scheme)=='cvmix_KPP')  ! NR Evaluate string comparison outside the loop. It is expensive.
    !___________________________________________________________________________
    ! Screen salinity
    a=0.0_WP
    do node=1, myDim_nod2D+eDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)
        !!PS do nz=1,nlevels_nod2d(node)-1
        do nz=nzmin,nzmax-1
            a=min(a,tr_arr(nz,node,2))
        enddo
    enddo
    
    !___________________________________________________________________________
    if(a<0.0_WP) then
        write (*,*)' --> pressure_bv: s<0 happens!', a
        pe_status=1
        do node=1, myDim_nod2D+eDim_nod2D
            nzmin = ulevels_nod2D(node)
            nzmax = nlevels_nod2D(node)
            !!PS do nz=1, nlevels_nod2d(node)-1
            do nz=nzmin, nzmax-1
                if (tr_arr(nz, node, 2) < 0) write (*,*) 'the model blows up at n=', mylist_nod2D(node), ' ; ', 'nz=', nz
            end do
        end do
    endif
    
    !___________________________________________________________________________
    do node=1, myDim_nod2D+eDim_nod2D
        nzmin = ulevels_nod2D(node)
        nzmax = nlevels_nod2D(node)
            
        !!PS nl1= nlevels_nod2d(node)-1
        rho      = 0.0_WP
        bulk_0   = 0.0_WP 
        bulk_pz  = 0.0_WP
        bulk_pz2 = 0.0_WP
        rhopot   = 0.0_WP
        ! also compute the maximum buoyancy gradient between the surface and any depth
        ! it will be used for computing MLD according to FESOM 1.4 implementation (after Large et al. 1997)
        db_max   = 0.0_WP
        dbsfc1   = 0.0_WP
        
        !_______________________________________________________________________
        ! apply equation of state
        do nz=nzmin, nzmax-1
            t=tr_arr(nz, node,1)
            s=tr_arr(nz, node,2)
            select case(state_equation)
                case(0)
                    call density_linear(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz), mesh)
                case(1)
                    call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz), mesh)
                case default !unknown
                    if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                    call par_ex(1)
            end select
        end do
            
        !NR split the loop here. The Intel compiler could not resolve that there is no dependency 
        !NR and did not vectorize the full loop. 
        !_______________________________________________________________________
        ! calculate density for MOC
        if (ldiag_dMOC) then
            !!PS do nz=1, nl1
            do nz=nzmin, nzmax-1
                rho(nz)              = bulk_0(nz) - 2000._WP*(bulk_pz(nz)   -2000._WP*bulk_pz2(nz))
                density_dmoc(nz,node)= rho(nz)*rhopot(nz)/(rho(nz)-200._WP)
                        !           density_dmoc(nz,node)   = rhopot(nz)
            end do
        end if 
            
        !_______________________________________________________________________    
        ! compute density for PGF 
        !!PS do nz=1, nl1
        do nz=nzmin, nzmax-1
            !___________________________________________________________________
            rho(nz) = bulk_0(nz) + Z_3d_n(nz,node)*(bulk_pz(nz)   + Z_3d_n(nz,node)*bulk_pz2(nz))
            !!PS rho(nz)=rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*Z_3d_n(nz,node))-density_0
            rho(nz) = rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*Z_3d_n(nz,node)*real(state_equation))-density_ref(nz,node)
            density_m_rho0(nz,node) = rho(nz)
            
            !___________________________________________________________________
            ! buoyancy difference between the surface and the grid points blow (adopted from FESOM 1.4)
            ! --> bring density of surface point adiabatically to the same 
            !     depth level as the deep point --> than calculate bouyancy 
            !     difference
            rho_surf=bulk_0(nzmin)   + Z_3d_n(nz,node)*(bulk_pz(nzmin)   + Z_3d_n(nz,node)*bulk_pz2(nzmin))
            !!PS rho_surf=rho_surf*rhopot(nzmin)/(rho_surf+0.1_WP*Z_3d_n(nz,node))-density_0
            rho_surf=rho_surf*rhopot(nzmin)/(rho_surf+0.1_WP*Z_3d_n(nz,node)*real(state_equation)) !-density_ref(nzmin,node)
            
            !!PS dbsfc1(nz) = -g * ( rho_surf - rho(nz) ) / (rho(nz)+density_0)      ! this is also required when KPP is ON
            !!PS dbsfc1(nz) = -g * density_0_r * ( rho_surf - rho(nz) )
            dbsfc1(nz) = -g * ( rho_surf - (rho(nz)+density_ref(nz,node)) ) / (rho(nz)+density_ref(nz,node))      ! this is also required when KPP is ON
            
            db_max = max(dbsfc1(nz)/abs(Z_3d_n(nzmin,node)-Z_3d_n(max(nz, nzmin+1),node)), db_max)
        end do
        
        dbsfc1(nzmax)=dbsfc1(nzmax-1)
        if (mixing_kpp) then ! in case KPP is ON store the buoyancy difference with respect to the surface (m/s2)
            dbsfc(nzmin:nzmax, node )=dbsfc1(nzmin:nzmax)
        end if
        
        !_______________________________________________________________________
        ! fill density levels that are occupied by the cavity with the density 
        ! that corresponds to a water mass that has the same temperature and salinity 
        ! like at the cavity-ocean interface --> compute water mass density that
        ! is replaced by the cavity
        if (nzmin>1) then 
            t=tr_arr(nzmin, node,1)
            s=tr_arr(nzmin, node,2)
            do nz=1, nzmin-1
                select case(state_equation)
                    case(0)
                        call density_linear(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz), mesh)
                    case(1)
                        call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz), mesh)
                    case default !unknown
                        if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                        call par_ex(1)
                end select
                !_______________________________________________________________
                rho(nz)= bulk_0(nz)   + Z_3d_n(nz,node)*(bulk_pz(nz)   + Z_3d_n(nz,node)*bulk_pz2(nz))
                rho(nz)=rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*Z_3d_n(nz,node)*real(state_equation))-density_ref(nz,node)
                density_m_rho0(nz,node) = rho(nz)
            end do
        end if
            
        !_______________________________________________________________________
        ! calculate pressure 
        if (trim(which_ale)=='linfs' .or. use_cavity .eqv. .true.) then
            !!PS hpressure(1, node)=-Z_3d_n(1,node)*rho(1)*g
            !!PS hpressure(1, node)=0.5_WP*hnode(1,node)*rho(1)*g
            !___________________________________________________________________
            ! compute pressure boundary condition at cavity-ocean interface
            ! hpressure(nzmin, node)
            if (nzmin>1) then ! cavity case 
                
                ! --> this as a upper pressure boundary condition incase of a 
                ! homogenous ocean, with no fluxes and boundary condition creates
                ! no pressure gradient errors 
                hpressure(nzmin, node)=0.5_WP*(zbar_3d_n(1,node)-zbar_3d_n(2,node))*rho(1)*g
                do nz=2,nzmin
                    a=0.5_WP*g*(rho(nz-1)*(zbar_3d_n(nz-1,node)-zbar_3d_n(nz,node))+rho(nz)*(zbar_3d_n(nz,node)-zbar_3d_n(nz+1,node)))
                    hpressure(nzmin, node)=hpressure(nzmin, node)+a                
                end do
                
                !  --> this as a upper pressure boundary condition incase of a 
                ! homogenous ocean, with no fluxes and boundary condition creates
                ! pgf errors 
                !!PS hpressure(nzmin, node)=-Z_3d_n(nzmin,node)*rho(nzmin)*g
            else
                hpressure(nzmin, node)=-Z_3d_n(nzmin,node)*rho(nzmin)*g
            end if 
            
            !___________________________________________________________________
            ! compute pressure below surface boundary
            do nz=nzmin+1,nzmax-1
                ! why 0.5 ... integrate g*rho*dz vertically, integrate half layer 
                ! thickness of previouse layer and half layer thickness of actual 
                ! layer to integrate pressure on mid depth level of actual layer
                a=0.5_WP*g*(rho(nz-1)*hnode(nz-1,node)+rho(nz)*hnode(nz,node))
                hpressure(nz, node)=hpressure(nz-1, node)+a
            end do
        end if    
            
        !___________________________________________________________________
        ! calculate mixed layer depth after Monterey and Levitus, (1997) who 
        ! compute MLD as the depth at which the density over depth differs 
        ! by 0.125 sigma units from the surface density (Griffies et al., 2009). 
        ! This MLD definition was also supported in FESOM1.4 (-->MLD2) and after 
        ! the definition of Large et al. (1997), who suggest to compute MLD 
        ! as the shallowest depth where the vertical derivative of buoyancy 
        ! is equal to a local critical buoyancy gradient (Griffies et al., 
        ! 2009) (-->MLD1). 
        ! BV frequency:  bvfreq(nl,:), squared value is stored   
        MLD1(node)=Z_3d_n(nzmin+1,node)
        MLD2(node)=Z_3d_n(nzmin+1,node)
        MLD1_ind(node)=nzmin+1
        MLD2_ind(node)=nzmin+1
        
        flag1=.true.
        flag2=.true.
        do nz=nzmin+1,nzmax-1
            bulk_up = bulk_0(nz-1) + zbar_3d_n(nz,node)*(bulk_pz(nz-1) + zbar_3d_n(nz,node)*bulk_pz2(nz-1)) 
            bulk_dn = bulk_0(nz)   + zbar_3d_n(nz,node)*(bulk_pz(nz)   + zbar_3d_n(nz,node)*bulk_pz2(nz))
            rho_up = bulk_up*rhopot(nz-1) / (bulk_up + 0.1_WP*zbar_3d_n(nz,node)*real(state_equation))  
            rho_dn = bulk_dn*rhopot(nz)   / (bulk_dn + 0.1_WP*zbar_3d_n(nz,node)*real(state_equation))  
            dz_inv=1.0_WP/(Z_3d_n(nz-1,node)-Z_3d_n(nz,node))  
            
            !_______________________________________________________________
            ! squared brunt väisälä frequence N^2 --> N^2>0 stratification is 
            ! stable, vertical elongated parcel is accelaratedtowards 
            ! initial point --> does oscillation with frequency N. 
            ! N^2<0 stratification is unstable vertical elongated parcel is 
            ! accelerated away from initial point 
            bvfreq(nz,node)  = -g*dz_inv*(rho_up-rho_dn)/density_0
!!PS             bvfreq(nz,node)  = -g*dz_inv*(rho_up-rho_dn)/density_ref(nz,node)
            
            !!PS !--> Why not like this ?
            !!PS bvfreq(nz,node)  = -g*dz_inv*(rho_up-rho_dn)/(rho_dn)
            
            !_______________________________________________________________
            ! define MLD following Large et al. 1997
            ! MLD is the shallowest depth where the local buoyancy gradient matches the maximum buoyancy gradient 
            ! between the surface and any discrete depth within the water column.
            if (bvfreq(nz, node) > db_max .and. flag1) then
                MLD1(node)    =Z_3d_n(nz, node)
                MLD1_ind(node)=nz
                flag1=.false.
            end if
            
            ! another definition of MLD after Levitus
            if ((rhopot(nz)-rhopot(nzmin) > sigma_theta_crit) .and. flag2) then
                MLD2(node)=MLD2(node)+(Z_3d_n(nz,node)-MLD2(node))/(rhopot(nz)-rhopot(nz-1)+1.e-20)*(rhopot(1)+sigma_theta_crit-rhopot(nz-1))
                MLD2_ind(node)=nz
                flag2=.false.
            elseif (flag2) then
                MLD2(node)=Z_3d_n(nz,node)
            end if
        end do
        
        if (flag2) MLD2_ind(node)=nzmax-1
                
        bvfreq(nzmin,node)=bvfreq(nzmin+1,node)
        bvfreq(nzmax,node)=bvfreq(nzmax-1,node) 
        !___________________________________________________________________
        ! The mixed layer depth 
        ! mixlay_depth    
        ! bv_ref
    end do
    !_______________________________________________________________________
    ! BV is defined on full levels except for the first and the last ones.
end subroutine pressure_bv
!
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) for linear free surface case
subroutine pressure_force_4_linfs(mesh)
    use g_config
    use g_PARSUP
    use mod_mesh
    use pressure_force_4_linfs_fullcell_interface
    use pressure_force_4_linfs_nemo_interface
    use pressure_force_4_linfs_shchepetkin_interface
    use pressure_force_4_linfs_cubicspline_interface
    use pressure_force_4_linfs_cavity_interface
    use pressure_force_4_linfs_easypgf_interface
    implicit none
    type(t_mesh), intent(in) , target :: mesh    
    !___________________________________________________________________________
    ! calculate pressure gradient force (PGF) for linfs with full cells
    if ( .not. use_partial_cell .and. .not. use_cavity_partial_cell) then
        call pressure_force_4_linfs_fullcell(mesh)
        
    elseif (use_cavity .and. use_cavity_partial_cell ) then
        if     (trim(which_pgf)=='sergey') then
            call pressure_force_4_linfs_cavity(mesh)  
        elseif (trim(which_pgf)=='shchepetkin') then    
            call pressure_force_4_linfs_shchepetkin(mesh)
        elseif (trim(which_pgf)=='easypgf') then
            call pressure_force_4_linfs_easypgf(mesh)    
        else
            write(*,*) '________________________________________________________'
            write(*,*) ' --> ERROR: the choosen form of pressure gradient       '
            write(*,*) '            calculation (PGF) in case of CAVITY is not supported for linfs !!!'
            write(*,*) '            see in namelist.oce --> which_pgf = sergey,   '
            write(*,*) '            shchepetkin, easypgf '
            write(*,*) '________________________________________________________'
            call par_ex(1)
        end if     
        
    !___________________________________________________________________________
    ! calculate pressure gradient force (PGF) for linfs with partiall cells
    else ! --> (trim(which_ale)=='linfs' .and. use_partial_cell )
        if     (trim(which_pgf)=='nemo') then
            call pressure_force_4_linfs_nemo(mesh)
        elseif (trim(which_pgf)=='shchepetkin') then
            call pressure_force_4_linfs_shchepetkin(mesh)
        elseif (trim(which_pgf)=='cubicspline') then
            call pressure_force_4_linfs_cubicspline(mesh)
        elseif (trim(which_pgf)=='easypgf') then
            call pressure_force_4_linfs_easypgf(mesh)    
        else
            write(*,*) '________________________________________________________'
            write(*,*) ' --> ERROR: the choosen form of pressure gradient       '
            write(*,*) '            calculation (PGF) is not supported for linfs !!!'
            write(*,*) '            see in namelist.oce --> which_pgf = nemo,   '
            write(*,*) '            shchepetkin, cubicspline '
            write(*,*) '________________________________________________________'
            call par_ex(1)
        end if 
    end if 
end subroutine pressure_force_4_linfs
!
!
!
!===============================================================================
! calculate pressure gradient force for linfs in case full cells
subroutine pressure_force_4_linfs_fullcell(mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    use g_PARSUP
    use g_config
    implicit none
    
    integer                  :: elem, elnodes(3), nle, ule,  nlz
    type(t_mesh), intent(in) , target :: mesh

#include  "associate_mesh.h"
    
    !___________________________________________________________________________
    ! loop over triangular elemments
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        ! number of levels at elem
        nle=nlevels(elem)-1
        ule=ulevels(elem)
            
        !_______________________________________________________________________
        ! node indices of elem 
        elnodes = elem2D_nodes(:,elem)
            
        !_______________________________________________________________________
        ! loop over mid-depth levels to calculate the pressure gradient 
        ! force (pgf) --> from top to bottom
        !!PS do nlz=1,nle
        do nlz=ule,nle
            pgf_x(nlz,elem) = sum(gradient_sca(1:3,elem)*hpressure(nlz,elnodes)/density_0)
            pgf_y(nlz,elem) = sum(gradient_sca(4:6,elem)*hpressure(nlz,elnodes)/density_0)
        end do 
    end do !-->do elem=1, myDim_elem2D
end subroutine pressure_force_4_linfs_fullcell   
!
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) like in NEMO based on NEMO ocean engine 
! Gurvan Madec, and the NEMO team gurvan.madec@locean-ipsl.umpc.fr, nemo 
! st@locean-ipsl.umpc.fr calculate vertical center index for linear 
! interpolation, densities are interpolated to shallowest nodal mid depth level 
! that contribute to element --> advantage no extrapolation neccessary
! Calculate pressure gradient force (PGF) like in NEMO based on NEMO ocean engine
! Gurvan Madec, and the NEMO team gurvan.madec@locean-ipsl.umpc.fr, nemo st@locean-ipsl.umpc.fr
! November 2015, – version 3.6 stable –
subroutine pressure_force_4_linfs_nemo(mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    use g_PARSUP
    use g_config
    use densityJM_components_interface
    use density_linear_interface
    implicit none
    
    logical             :: do_interpTS=.true.
    integer             :: elem, elnodes(3), nle, ule, nlz, nln(3), uln(3), ni, nlc, nlce
    real(kind=WP)       :: hpress_n_bottom(3)
    real(kind=WP)       :: interp_n_dens(3), interp_n_temp, interp_n_salt, &
                           dZn, dZn_i, dh, dval, mean_e_rho,dZn_rho_grad(2)
    real(kind=WP)       :: rhopot, bulk_0, bulk_pz, bulk_pz2
    type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
    !___________________________________________________________________________
    ! loop over triangular elemments
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        ! nle...number of mid-depth levels at elem
        nle          = nlevels(elem)-1
        ule          = ulevels(elem)
        
        ! node indices of elem 
        elnodes = elem2D_nodes(:,elem)
        
        !_______________________________________________________________________
        ! 1) loop over mid-depth levels to calculate the pressure gradient 
        !    force (pgf) --> goes until one layer above the bottom 
        !        
        !     :             :      --> nle-2
        !     :             :
        ! ----------   ----------
        ! 
        !     o            o       --> nle-1
        ! T,S,rho(nle-1)
        ! ----------   ----------
        ! 
        !     x <--- nle --o linear interpolate densities to shallowest mid depth level
        !     o
        !              ----------
        ! T,S,rho(nle) //////////
        ! ----------   //////////
        ! //////////   //////////
        ! //////////   //////////
        ! //////////   //////////
        !!PS do nlz=1,nle-1
        do nlz=ule,nle-1 
            pgf_x(nlz,elem) = sum(gradient_sca(1:3,elem)*hpressure(nlz,elnodes)/density_0)
            pgf_y(nlz,elem) = sum(gradient_sca(4:6,elem)*hpressure(nlz,elnodes)/density_0)
        end do
                
        !_______________________________________________________________________
        ! 2) interpolate/extrapolate nodal density linearly to elemental 
        !    bottom depth
        !_______________________________________________________________________
        ! nln...number of mid-depth levels at node
        nln     = nlevels_nod2d(elnodes)-1
        uln     = ulevels_nod2d(elnodes)
        !!PS uln     = ulevels_nod2d(elnodes)
        
        !_______________________________________________________________________
        ! calculate mid depth element level --> Z_e
        zbar_n       = 0.0_WP
        Z_n          = 0.0_WP
        zbar_n(nle+1)= zbar_e_bot(elem)
        Z_n(nle)     = zbar_n(nle+1) + helem(nle,elem)/2.0_WP
        !!PS do nlz=nle,2,-1
        do nlz=nle,ule+1,-1
            zbar_n(nlz) = zbar_n(nlz+1) + helem(nlz,elem)
            Z_n(nlz-1)  = zbar_n(nlz)   + helem(nlz-1,elem)/2.0_WP
        end do
        !!PS zbar_n(1)    = zbar_n(2) + helem(1,elem)
        zbar_n(ule)    = zbar_n(ule+1) + helem(ule,elem)
        ! --> zbar_n    ... depth of level at elements
        ! --> Z_n       ... depth of mid-depth level at elements
        ! --> zbar_n_3d ... depth of level at nodes
        ! --> Z_n_3d    ... depth of mid-depth level at nodes
        
        !_______________________________________________________________________
        ! loop over nodal indices of element
        do ni=1,3
            !___________________________________________________________________
            ! calculate vertical center index for linear interpolation, densities 
            ! are interpolated to shallowest nodal mid depth level that contribute
            ! to element --> advantage no extrapolation neccessary
            !!PS nlc   = minloc( Z_3d_n(1:nln(ni),elnodes(ni))-maxval(Z_3d_n(nle,elnodes)),1,& 
            !!PS                 Z_3d_n(1:nln(ni),elnodes(ni))-maxval(Z_3d_n(nle,elnodes))>0.0_WP &  ! mask for index selection
            !!PS                 )+1
            nlc   = minloc( Z_3d_n(uln(ni):nln(ni),elnodes(ni))-maxval(Z_3d_n(nle,elnodes)),1,& 
                            Z_3d_n(uln(ni):nln(ni),elnodes(ni))-maxval(Z_3d_n(nle,elnodes))>0.0_WP &  ! mask for index selection
                            )+1
            nlc   = min(nlc,nln(ni)) 
            dZn   = Z_3d_n(nlc,elnodes(ni))    -Z_3d_n(nlc-1,elnodes(ni))
            dZn_i = maxval(Z_3d_n(nle,elnodes))-Z_3d_n(nlc-1,elnodes(ni))
            dh    = minval(hnode(nle  ,elnodes)) 
            
            !___________________________________________________________________
            !! Option (A): interpolate density directly ...
            !if (.not. do_interpTS ) then 
            !    dval              = density_m_rho0(nlc  ,elnodes(ni))-density_m_rho0(nlc-1,elnodes(ni))
            !    interp_n_dens(ni) = density_m_rho0(nlc-1,elnodes(ni)) + (dval/dZn*dZn_i)
            !    
            !! Option (B): NEMO ocean engine Gurvan Madec, and the NEMO team suggest 
            !! not to linearly interpolate the density towards the bottom rather to 
            !! interpolate temperature and salinity and calculate from them the 
            !! bottom density to account for the non linearities in the equation of 
            !! state ...
            !else
                ! ... interpolate temperature and saltinity ...
                dval          = tr_arr(nlc,  elnodes(ni),1) - tr_arr(nlc-1,elnodes(ni),1)
                interp_n_temp = tr_arr(nlc-1,elnodes(ni),1) + (dval/dZn*dZn_i)
                dval          = tr_arr(nlc  ,elnodes(ni),2) - tr_arr(nlc-1,elnodes(ni),2)
                interp_n_salt = tr_arr(nlc-1,elnodes(ni),2) + (dval/dZn*dZn_i)
                
                ! calculate density at element mid-depth bottom depth via 
                ! equation of state from linear interpolated temperature and 
                ! salinity
                select case(state_equation)
                    case(0)
                        call density_linear(interp_n_temp, interp_n_salt, bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)
                    case(1)
                        call densityJM_components(interp_n_temp, interp_n_salt, bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)
                    case default !unknown
                        if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                        call par_ex(1)
                end select
                interp_n_dens(ni) = bulk_0 + Z_n(nle)*(bulk_pz + Z_n(nle)*bulk_pz2)
                !!PS interp_n_dens(ni) = interp_n_dens(ni)*rhopot/(interp_n_dens(ni)+0.1_WP*Z_n(nle))*real(state_equation))-density_0
                interp_n_dens(ni) = interp_n_dens(ni)*rhopot/(interp_n_dens(ni)+0.1_WP*Z_n(nle)*real(state_equation))-density_ref(nle,elnodes(ni))
            !end if 
                
            ! calculate proper starting level for extrapolation of 
            ! bottom pressure
            nlce = min(nlc,nle)
                
            !___________________________________________________________________
            ! integrate nodal density until mid-depth level of bottom 
            ! element using linearly interpolated nodal bottom density 
            ! (interp_n_dens) to obtain bottom pressure value
            hpress_n_bottom(ni) = hpressure(nlce-1, elnodes(ni))        &
                                    + 0.5_WP*g*                         & 
                                    (density_m_rho0(nlce-1,elnodes(ni))*&
                                              hnode(nlce-1,elnodes(ni)) & 
                                    +                                   &  
                                    interp_n_dens(ni)*dh                &
                                    )                             
        end do ! --> do ni=1,3
        !_______________________________________________________________________
        pgf_x(nle,elem) = sum(gradient_sca(1:3,elem)*hpress_n_bottom)/density_0
        pgf_y(nle,elem) = sum(gradient_sca(4:6,elem)*hpress_n_bottom)/density_0
        
    end do ! --> do elem=1, myDim_elem2D
end subroutine pressure_force_4_linfs_nemo
!
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) based on "Shchepetkin 
! and McWilliams," A method for computing horizontal pressure-gradient
! force in an oceanic model with a nonaligned vertical coordinate",
! Journal of geophysical research, vol 108, no C3, 3090
! --> based on density jacobian method ...
! calculate PGF for linfs with partiell cell on/off
! First coded by P. Scholz for FESOM2.0, 08.02.2019
subroutine pressure_force_4_linfs_shchepetkin(mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    use g_PARSUP
    use g_config
    implicit none
    
    integer             :: elem, elnodes(3), nle, ule, nlz, idx(3),ni
    real(kind=WP)       :: int_dp_dx(2), drho_dx, dz_dx, aux_sum
    real(kind=WP)       :: dx10(3), dx20(3), dx21(3), df10(3), df21(3), drho_dz(3)
    type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
    !___________________________________________________________________________
    ! loop over triangular elemments
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        ! nle...number of mid-depth levels at elem
        nle     = nlevels(elem)-1
        ule     = ulevels(elem)
        
        ! node indices of elem 
        elnodes = elem2D_nodes(:,elem)
        
        !_______________________________________________________________________
        ! calculate mid depth element level --> Z_e
        ! nle...number of mid-depth levels at elem
        zbar_n       = 0.0_WP
        Z_n          = 0.0_WP
        zbar_n(nle+1)= zbar_e_bot(elem)
        Z_n(nle)     = zbar_n(nle+1) + helem(nle,elem)*0.5_WP
        do nlz=nle,ule+1,-1
            zbar_n(nlz) = zbar_n(nlz+1) + helem(nlz,elem)
            Z_n(nlz-1)  = zbar_n(nlz)   + helem(nlz-1,elem)*0.5_WP
        end do
        zbar_n(ule) = zbar_n(ule+1) + helem(ule,elem)
        
        !_______________________________________________________________________
        ! Calculate pressure gradient force (PGF) based on "Shchepetkin 
        ! and McWilliams," A method for computing horizontal pressure-gradient
        ! force in an oceanic model with a nonaligned vertical coordinate",
        ! Journal of geophysical research, vol 108, no C3, 3090
        ! --> based on density jacobian method ...
        ! --> equation 1.7
        ! 
        ! -1/rho_0* dP/dx|_z = -1/rho_0*dP/dx|_z=eta - g/rho_0*int_z^eta( drho/dx|_z)*dz'
        !
        !                    = g*rho(eta)/rho_0*deta/dx
        !                      - g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
        !
        ! --> eta...free surface elevation for linfs start at zero
        
        
        !_______________________________________________________________________
        ! calculate pressure gradient for surface layer
        nlz=ule
        if (ule==1) then 
            !___________________________________________________________________
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(1)    = aux_sum
                
            !___________________________________________________________________
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(2)    = aux_sum
        else
            !___________________________________________________________________
            ! vertical gradient --> with average density and average 
            ! mid-depth level on element
            !       x0        x1               x2
            !  -----o---------o----------------o-------
            !       f0        f1               f2
            ! --> use Newtonsche interpolation polynom of second order to calculate 
            ! central difference quotient for not equidestant distributed values
            ! f(x) = a0 + a1*(x-x0) + a2*(x-x0)*(x-x1)
            !  | |
            !  | +--> derivative : f'(x) = a1 + a2*[(x-x1)+(x-x0)]
            !  |
            !  +--> determine Coefficients a0, a1, a2 for f(x)=f1,f2,f3
            !  |
            !  +--> a0 = f0 , 
            !       a1 = (f1-f0)/(x1-x0) , 
            !       a2 = ((x1-x0)*(f2-f1)-(x2-x1)*(f1-f0))/((x2-x0)*(x2-x1)*(x1-x0))
            !
            ! f'(x) = (f1-f0)/(x1-x0) 
            !         - [(x1-x0)*(f2-f1)-(x2-x1)*(f1-f0)]/[(x2-x1)*(x1-x0)*(x2-x0)]
            !           *[(x-x1)+(x-x0)]
            !
            idx = (/1, 1, 1/)*nlz
            idx = idx - ulevels_nod2D(elnodes)
            do ni=1,3
                if (idx(ni)==0) then
                    ! elemental surface index nlz ist also surface index of vertice ni
                    !_______________________________________________________________
                    ! interpolation coefficients 
                    dx10(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    dx21(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz+1,elnodes(ni))
                    dx20(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    
                    df10(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz+2,elnodes(ni))-density_m_rho0(nlz+1,elnodes(ni))
                    !_______________________________________________________________
                    ! compute vertical density gradient
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                *((Z_n(nlz)-Z_3d_n(nlz+1,elnodes(ni))) + &
                                (Z_n(nlz)-Z_3d_n(nlz  ,elnodes(ni))))
                else
                    ! elemental surface index nlz ist deeper than surface index of vertices ni
                    !_______________________________________________________________
                    ! interpolation coefficients 
                    dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                    dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                    
                    df10(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                    !_______________________________________________________________
                    ! compute vertical density gradient
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                *((Z_n(nlz)-Z_3d_n(nlz  ,elnodes(ni))) + &
                                (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))))
                end if 
            end do
            
            !___________________________________________________________________
            ! -->...-g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            dz_dx           = sum(gradient_sca(1:3,elem)*Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0               
            pgf_x(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(1)    = aux_sum
            
            !___________________________________________________________________
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            dz_dx           = sum(gradient_sca(4:6,elem)*Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0             
            pgf_y(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(2)    = aux_sum
        end if 
        
        !_______________________________________________________________________
        ! calculate pressure gradient for surface layer +1 until one layer above bottom
        do nlz=ule+1,nle-1
            !___________________________________________________________________
            ! - g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
            ! --> in case linfs: dz_dx == 0.0
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
            int_dp_dx(1)    = int_dp_dx(1) + aux_sum
            
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
            int_dp_dx(2)    = int_dp_dx(2) + aux_sum
            
        end do ! --> do nlz=1,nle-1
        
        !_______________________________________________________________________
        ! calculate pressure gradient for bottom layer
        ! vertical gradient --> with average density on element
        !       x0        x1               x2
        !  -----o---------o----------------o-------
        !       f0        f1               f2
        ! --> use Newtonsche interpolation polynom of second order to calculate 
        ! central difference quotient for not equidestant distributed values
        ! f(x) = a0 + a1*(x-x0) + a2*(x-x0)*(x-x1)
        !  | |
        !  | +--> derivative : f'(x) = a1 + a2*[(x-x1)+(x-x0)]
        !  |
        !  +--> determine Coefficients a0, a1, a2 for f(x)=f1,f2,f3
        !  |
        !  +--> a0 = f0 , 
        !       a1 = (f1-f0)/(x1-x0) , 
        !       a2 = ((x1-x0)*(f2-f1)-(x2-x1)*(f1-f0))/((x2-x0)*(x2-x1)*(x1-x0))
        !
        ! f'(x) = (f1-f0)/(x1-x0) 
        !         - [(x1-x0)*(f2-f1)-(x2-x1)*(f1-f0)]/[(x2-x1)*(x1-x0)*(x2-x0)]
        !           *[(x-x1)+(x-x0)]
        !
        nlz = nle 
        idx = (/1, 1, 1/)*nlz
        idx = nlevels_nod2D(elnodes)-1 - idx
        do ni=1,3
            if (idx(ni)==0) then
                ! elemental bottom index nlz ist also bottom index of vertice ni
                !_______________________________________________________________
                ! interpolation coefficients 
                dx10(ni)   = Z_3d_n(nlz-1,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                
                df10(ni)   = density_m_rho0(nlz-1,elnodes(ni))-density_m_rho0(nlz-2,elnodes(ni))
                df21(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                !_______________________________________________________________
                ! compute vertical density gradient
                drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                             *((Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) + &
                               (Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni))))
            else
                ! elemental bottom index nlz ist shallower than bottom index of vertice ni
                !_______________________________________________________________
                ! interpolation coefficients 
                dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                
                df10(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                df21(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                !_______________________________________________________________
                ! compute vertical density gradient
                drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                             *((Z_n(nlz)-Z_3d_n(nlz  ,elnodes(ni))) + &
                               (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))))
            end if 
        end do
        
        !_______________________________________________________________________
        ! - g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
        ! zonal gradients
        drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
        dz_dx           = sum(gradient_sca(1:3,elem)*Z_3d_n(nlz,elnodes))
        !!PS aux_sum         = (drho_dx-drho_dz)*dz_dx)*helem(nlz,elem)*g/density_0
        aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0               
        pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
        
        !_______________________________________________________________________
        ! meridional gradients
        drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
        dz_dx           = sum(gradient_sca(4:6,elem)*Z_3d_n(nlz,elnodes))
        !!PS aux_sum         = (drho_dx-drho_dz*dz_dx)*helem(nlz,elem)*g/density_0
        aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0             
        pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
        
    end do ! --> do elem=1, myDim_elem2D
end subroutine pressure_force_4_linfs_shchepetkin
!
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) 
! First coded by P. Scholz for FESOM2.0, 08.02.2019
subroutine pressure_force_4_linfs_easypgf(mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    use g_PARSUP
    use g_config
    use densityJM_components_interface
    use density_linear_interface
    implicit none
    
    integer             :: elem, elnodes(3), nle, ule, nlz, idx(3),ni
    real(kind=WP)       :: int_dp_dx(2), drho_dx, aux_sum
    real(kind=WP)       :: dx10(3), dx20(3), dx21(3)
    real(kind=WP)       :: t0(3), dt10(3), dt21(3)
    real(kind=WP)       :: s0(3), ds10(3), ds21(3)
    real(kind=WP)       :: rho_at_Zn(3), temp_at_Zn(3), salt_at_Zn(3), drho_dz(3), aux_dref
    real(kind=WP)       :: rhopot(3), bulk_0(3), bulk_pz(3), bulk_pz2(3)
    real(kind=WP)       :: dref_rhopot, dref_bulk_0, dref_bulk_pz, dref_bulk_pz2
    
    type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
    !___________________________________________________________________________
    ! loop over triangular elemments
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        ! nle...number of mid-depth levels at elem
        nle     = nlevels(elem)-1
        ule     = ulevels(elem)
        
        ! node indices of elem 
        elnodes = elem2D_nodes(:,elem)
        
        !_______________________________________________________________________
        ! calculate mid depth element level --> Z_e
        ! nle...number of mid-depth levels at elem
        zbar_n       = 0.0_WP
        Z_n          = 0.0_WP
        zbar_n(nle+1)= zbar_e_bot(elem)
        Z_n(nle)     = zbar_n(nle+1) + helem(nle,elem)*0.5_WP
        do nlz=nle,ule+1,-1
            zbar_n(nlz) = zbar_n(nlz+1) + helem(nlz,elem)
            Z_n(nlz-1)  = zbar_n(nlz)   + helem(nlz-1,elem)*0.5_WP
        end do
        zbar_n(ule) = zbar_n(ule+1) + helem(ule,elem)
        
        !_______________________________________________________________________
        aux_dref     = density_0
        if (use_cavity .and. .not. use_density_ref) then
            select case(state_equation)
                case(0)
                    call density_linear(density_ref_T, density_ref_S, dref_bulk_0, dref_bulk_pz, dref_bulk_pz2, dref_rhopot, mesh)
                case(1)
                    call densityJM_components(density_ref_T, density_ref_S, dref_bulk_0, dref_bulk_pz, dref_bulk_pz2, dref_rhopot, mesh)
                case default !unknown
                    if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                    call par_ex(1)
            end select
        end if
        
        !_______________________________________________________________________
        ! calculate pressure gradient for surface layer
        nlz=ule
        if (ule==1) then 
            !___________________________________________________________________
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(1)    = aux_sum
                
            !___________________________________________________________________
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(2)    = aux_sum
        else
            !___________________________________________________________________
            ! vertical gradient --> with average density and average 
            ! mid-depth level on element
            !       x0        x1               x2
            !  -----o---------o----------------o-------
            !       f0        f1               f2
            ! --> use Newtonsche interpolation polynom of second order to calculate 
            ! central difference quotient for not equidestant distributed values
            ! f(x) = a0 + a1*(x-x0) + a2*(x-x0)*(x-x1)
            !  | |
            !  | +--> derivative : f'(x) = a1 + a2*[(x-x1)+(x-x0)]
            !  |
            !  +--> determine Coefficients a0, a1, a2 for f(x)=f1,f2,f3
            !  |
            !  +--> a0 = f0 , 
            !       a1 = (f1-f0)/(x1-x0) , 
            !       a2 = ((x1-x0)*(f2-f1)-(x2-x1)*(f1-f0))/((x2-x0)*(x2-x1)*(x1-x0))
            !
            ! f(x) = f0 
            !        + (f1-f0)/(x1-x0)*(x-x0) 
            !        + [(x1-x0)*(f2-f1)-(x2-x1)*(f1-f0)]/[(x2-x0)*(x2-x1)*(x1-x0)]
            !          *[(x-x1)+(x-x0)]
            !
            !
            
            ! account for reference density when using cavities
            if (use_cavity .and. .not. use_density_ref) then
                aux_dref = dref_bulk_0   + Z_n(nlz)*dref_bulk_pz   + Z_n(nlz)*dref_bulk_pz2
                aux_dref = aux_dref*dref_rhopot/(aux_dref+0.1_WP*Z_n(nlz))
            end if 
            
            idx = (/1, 1, 1/)*nlz
            idx = idx - ulevels_nod2D(elnodes)
            do ni=1,3
                if (idx(ni)==0) then
                    ! elemental surface index nlz ist also surface index of vertice ni
                    !___________________________________________________________________
                    ! compute interpolation coefficients
                    dx10(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    dx21(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz+1,elnodes(ni))
                    dx20(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    
                    t0(ni)     = tr_arr(nlz  ,elnodes(ni),1)
                    dt10(ni)   = tr_arr(nlz+1,elnodes(ni),1)-tr_arr(nlz  ,elnodes(ni),1)
                    dt21(ni)   = tr_arr(nlz+2,elnodes(ni),1)-tr_arr(nlz+1,elnodes(ni),1)
                    
                    s0(ni)     = tr_arr(nlz  ,elnodes(ni),2)
                    ds10(ni)   = tr_arr(nlz+1,elnodes(ni),2)-tr_arr(nlz  ,elnodes(ni),2)
                    ds21(ni)   = tr_arr(nlz+2,elnodes(ni),2)-tr_arr(nlz+1,elnodes(ni),2)
                    !___________________________________________________________________
                    ! interpoalte vertice temp and salinity to elemental level Z_n
                    temp_at_Zn(ni) = t0(ni) &
                                    + dt10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni))) & 
                                    + (dx10(ni)*dt21(ni)-dx21(ni)*dt10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                    (Z_n(nlz)-Z_3d_n(nlz+1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))
                    salt_at_Zn(ni) = s0(ni) &
                                    + ds10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni))) & 
                                    + (dx10(ni)*ds21(ni)-dx21(ni)*ds10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                    (Z_n(nlz)-Z_3d_n(nlz+1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))
                    !___________________________________________________________________
                    ! compute density from state equation 
                    select case(state_equation)
                        case(0)
                            call density_linear(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)   
                        case(1)
                            call densityJM_components(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)   
                        case default !unknown
                            if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                            call par_ex(1)
                    end select
                    rho_at_Zn(ni) = bulk_0(ni) + Z_n(nlz)*(bulk_pz(ni) + Z_n(nlz)*bulk_pz2(ni))
                    rho_at_Zn(ni) = rho_at_Zn(ni)*rhopot(ni)/(rho_at_Zn(ni)+0.1_WP*Z_n(nlz)*real(state_equation))-aux_dref
                    
                else
                    ! elemental surface index nlz ist deeper than surface index of vertices ni
                    !___________________________________________________________________
                    ! compute interpolation coefficients
                    dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                    dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                    
                    t0(ni)     = tr_arr(nlz-1,elnodes(ni),1)
                    dt10(ni)   = tr_arr(nlz  ,elnodes(ni),1)-tr_arr(nlz-1,elnodes(ni),1)
                    dt21(ni)   = tr_arr(nlz+1,elnodes(ni),1)-tr_arr(nlz  ,elnodes(ni),1)
                    
                    s0(ni)     = tr_arr(nlz-1,elnodes(ni),2)
                    ds10(ni)   = tr_arr(nlz  ,elnodes(ni),2)-tr_arr(nlz-1,elnodes(ni),2)
                    ds21(ni)   = tr_arr(nlz+1,elnodes(ni),2)-tr_arr(nlz  ,elnodes(ni),2)
                    !___________________________________________________________________               
                    ! interpoalte vertice temp and salinity to elemental level Z_n
                    temp_at_Zn(ni) = t0(ni) &
                                    + dt10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                                    + (dx10(ni)*dt21(ni)-dx21(ni)*dt10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                    (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
                    salt_at_Zn(ni) = s0(ni) &
                                    + ds10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                                    + (dx10(ni)*ds21(ni)-dx21(ni)*ds10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                    (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
                    !___________________________________________________________________
                    ! compute density from state equation 
                    select case(state_equation)
                        case(0)
                            call density_linear(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)    
                        case(1)
                            call densityJM_components(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)    
                        case default !unknown
                            if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                            call par_ex(1)
                    end select
                    rho_at_Zn(ni) = bulk_0(ni) + Z_n(nlz)*(bulk_pz(ni) + Z_n(nlz)*bulk_pz2(ni))
                    rho_at_Zn(ni) = rho_at_Zn(ni)*rhopot(ni)/(rho_at_Zn(ni)+0.1_WP*Z_n(nlz)*real(state_equation))-aux_dref
                end if 
            end do
            
            !___________________________________________________________________
            ! -->...-g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*rho_at_Zn)
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(1)    = aux_sum
            
            !___________________________________________________________________
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*rho_at_Zn)
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0            
            pgf_y(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(2)    = aux_sum
        end if 
        
        !_______________________________________________________________________
        ! calculate pressure gradient for surface layer +1 until one layer above bottom
        do nlz=ule+1,nle-1
            !___________________________________________________________________
            ! - g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
            ! --> in case linfs: dz_dx == 0.0
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
            int_dp_dx(1)    = int_dp_dx(1) + aux_sum
            
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
            int_dp_dx(2)    = int_dp_dx(2) + aux_sum
            
        end do ! --> do nlz=1,nle-1
        
        !_______________________________________________________________________
        ! calculate pressure gradient for bottom layer
        ! vertical gradient --> with average density on element
        !       x0        x1               x2
        !  -----o---------o----------------o-------
        !       f0        f1               f2
        ! --> use Newtonsche interpolation polynom of second order to calculate 
        ! central difference quotient for not equidestant distributed values
        ! f(x) = a0 + a1*(x-x0) + a2*(x-x0)*(x-x1)
        !  | |
        !  | +--> derivative : f'(x) = a1 + a2*[(x-x1)+(x-x0)]
        !  |
        !  +--> determine Coefficients a0, a1, a2 for f(x)=f1,f2,f3
        !  |
        !  +--> a0 = f0 , 
        !       a1 = (f1-f0)/(x1-x0) , 
        !       a2 = ((x1-x0)*(f2-f1)-(x2-x1)*(f1-f0))/((x2-x0)*(x2-x1)*(x1-x0))
        !
        ! f'(x) = (f1-f0)/(x1-x0) 
        !         - [(x1-x0)*(f2-f1)-(x2-x1)*(f1-f0)]/[(x2-x1)*(x1-x0)*(x2-x0)]
        !           *[(x-x1)+(x-x0)]
        !
        nlz = nle 
        idx = (/1, 1, 1/)*nlz
        idx = nlevels_nod2D(elnodes)-1 - idx
        do ni=1,3
            if (idx(ni)==0) then
                ! elemental bottom index nlz ist also bottom index of vertice ni
                !___________________________________________________________________
                ! compute interpolation coefficients
                dx10(ni)   = Z_3d_n(nlz-1,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                
                t0(ni)     = tr_arr(nlz-2,elnodes(ni),1)
                dt10(ni)   = tr_arr(nlz-1,elnodes(ni),1)-tr_arr(nlz-2,elnodes(ni),1)
                dt21(ni)   = tr_arr(nlz  ,elnodes(ni),1)-tr_arr(nlz-1,elnodes(ni),1)
                
                s0(ni)     = tr_arr(nlz-2,elnodes(ni),2)
                ds10(ni)   = tr_arr(nlz-1,elnodes(ni),2)-tr_arr(nlz-2,elnodes(ni),2)
                ds21(ni)   = tr_arr(nlz  ,elnodes(ni),2)-tr_arr(nlz-1,elnodes(ni),2)
                !___________________________________________________________________
                ! interpoalte vertice temp and salinity to elemental level Z_n
                temp_at_Zn(ni) = t0(ni) &
                                + dt10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni))) & 
                                + (dx10(ni)*dt21(ni)-dx21(ni)*dt10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni)))
                salt_at_Zn(ni) = s0(ni) &
                                + ds10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni))) & 
                                + (dx10(ni)*ds21(ni)-dx21(ni)*ds10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni)))
                !___________________________________________________________________
                ! compute density from state equation 
                select case(state_equation)
                    case(0)
                        call density_linear(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case(1)
                        call densityJM_components(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case default !unknown
                        if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                        call par_ex(1)
                end select
                rho_at_Zn(ni) = bulk_0(ni) + Z_n(nlz)*(bulk_pz(ni) + Z_n(nlz)*bulk_pz2(ni))
                rho_at_Zn(ni) = rho_at_Zn(ni)*rhopot(ni)/(rho_at_Zn(ni)+0.1_WP*Z_n(nlz)*real(state_equation))-aux_dref
                
            else
                ! elemental bottom index nlz ist shallower than bottom index of vertice ni
                !___________________________________________________________________
                ! compute interpolation coefficients
                dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                
                t0(ni)     = tr_arr(nlz-1,elnodes(ni),1)
                dt10(ni)   = tr_arr(nlz  ,elnodes(ni),1)-tr_arr(nlz-1,elnodes(ni),1)
                dt21(ni)   = tr_arr(nlz+1,elnodes(ni),1)-tr_arr(nlz  ,elnodes(ni),1)
                
                s0(ni)     = tr_arr(nlz-1,elnodes(ni),2)
                ds10(ni)   = tr_arr(nlz  ,elnodes(ni),2)-tr_arr(nlz-1,elnodes(ni),2)
                ds21(ni)   = tr_arr(nlz+1,elnodes(ni),2)-tr_arr(nlz  ,elnodes(ni),2)
                !___________________________________________________________________               
                ! interpoalte vertice temp and salinity to elemental level Z_n
                temp_at_Zn(ni) = t0(ni) &
                                + dt10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                                + (dx10(ni)*dt21(ni)-dx21(ni)*dt10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
                salt_at_Zn(ni) = s0(ni) &
                                + ds10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                                + (dx10(ni)*ds21(ni)-dx21(ni)*ds10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
                !___________________________________________________________________
                ! compute density from state equation 
                select case(state_equation)
                    case(0)
                        call density_linear(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case(1)
                        call densityJM_components(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case default !unknown
                        if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                        call par_ex(1)
                end select
                rho_at_Zn(ni) = bulk_0(ni) + Z_n(nlz)*(bulk_pz(ni) + Z_n(nlz)*bulk_pz2(ni))
                rho_at_Zn(ni) = rho_at_Zn(ni)*rhopot(ni)/(rho_at_Zn(ni)+0.1_WP*Z_n(nlz)*real(state_equation))-aux_dref
            end if 
        end do
        
        !_______________________________________________________________________
        ! - g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
        ! zonal gradients
        drho_dx         = sum(gradient_sca(1:3,elem)*rho_at_Zn)
        aux_sum         = drho_dx*helem(nlz,elem)*g/density_0            
        pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
        
        !_______________________________________________________________________
        ! meridional gradients
        drho_dx         = sum(gradient_sca(4:6,elem)*rho_at_Zn)
        aux_sum         = drho_dx*helem(nlz,elem)*g/density_0            
        pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
        
    end do ! --> do elem=1, myDim_elem2D
end subroutine pressure_force_4_linfs_easypgf
!
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) via cubicspline used in FEOSM1.4
! First coded by Q. Wang for FESOM1.4, adapted by P. Scholz for FESOM2.0, 08.02.2019
subroutine pressure_force_4_linfs_cubicspline(mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    use g_PARSUP
    use g_config
    implicit none
    
    integer             :: elem, elnodes(3), nle, ule, nlz, nlc, ni, node, nln(3), uln(3), dd
    real(kind=WP)       :: int_dp_dx(2), drho_dx, dz_dx, drho_dz, auxp
    real(kind=WP)       :: dx10, dx20, dx21, df10, df21
    real(kind=WP)       :: interp_n_dens(3) 
    integer             :: s_ind(4)
    real(kind=WP)       :: s_z(4), s_dens(4), s_H, aux1, aux2, s_dup, s_dlo
    real(kind=WP)       :: a, b, c, d, dz 
    type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
    !___________________________________________________________________________
    ! loop over triangular elemments
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        ! nle...number of mid-depth levels at elem
        nle     = nlevels(elem)-1
        ule     = ulevels(elem)
        
        ! node indices of elem 
        elnodes = elem2D_nodes(:,elem)
        
        ! calculate mid depth element level --> Z_e
        ! nle...number of mid-depth levels at elem
        zbar_n       = 0.0_WP
        Z_n          = 0.0_WP
        zbar_n(nle+1)= zbar_e_bot(elem)
        Z_n(nle)     = zbar_n(nle+1) + helem(nle,elem)/2.0_WP
        do nlz=nle,2,-1
            zbar_n(nlz) = zbar_n(nlz+1) + helem(nlz,elem)
            Z_n(nlz-1)  = zbar_n(nlz)   + helem(nlz-1,elem)/2.0_WP
        end do
        zbar_n(1) = zbar_n(2) + helem(1,elem)
        
        ! nln...number of mid depth levels at each node that corresponds to 
        ! element elem
        nln=nlevels_nod2D(elnodes)-1
        uln=ulevels_nod2D(elnodes)
        
        !_______________________________________________________________________
        ! calculate pressure gradient for surface layer
        int_dp_dx     = 0.0_WP
        nlz=ule
        
        if (ule==1) then ! in case no cavity
            !___________________________________________________________________
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            auxp            = drho_dx*helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = int_dp_dx(1) + auxp*0.5_WP
            int_dp_dx(1)    = int_dp_dx(1) + auxp
            
            !___________________________________________________________________
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            auxp            = drho_dx*helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = int_dp_dx(2) + auxp*0.5_WP
            int_dp_dx(2)    = int_dp_dx(2) + auxp
            
        else! in case of cavity (includes pressure bnd condition under cavity)
            !___________________________________________________________________
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            
            ! cavity-ocean interface pressure boundary condition 
            auxp            = drho_dx*(-zbar_e_srf(elem))*g/density_0
            int_dp_dx(1)    = int_dp_dx(1) + auxp
            
            auxp            = drho_dx*helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = int_dp_dx(1) + auxp*0.5_WP
            int_dp_dx(1)    = int_dp_dx(1) + auxp
                
            !___________________________________________________________________
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            
            ! cavity-ocean interface pressure boundary condition 
            auxp            = drho_dx*(-zbar_e_srf(elem))*g/density_0
            int_dp_dx(2)    = int_dp_dx(2) + auxp
            
            auxp            = drho_dx*helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = int_dp_dx(2) + auxp*0.5_WP
            int_dp_dx(2)    = int_dp_dx(2) + auxp
        endif 
        
        !_______________________________________________________________________
        ! calculate pressure gradient for surface layer+1 until one layer above bottom
        !!PS do nlz=1,nle-1
        do nlz=ule+1,nle-1
            !___________________________________________________________________
            ! - g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
            ! --> in case linfs: dz_dx == 0.0
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            auxp            = drho_dx*helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = int_dp_dx(1) + auxp*0.5_WP
            int_dp_dx(1)    = int_dp_dx(1) + auxp
            
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            auxp            = drho_dx*helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = int_dp_dx(2) + auxp*0.5_WP
            int_dp_dx(2)    = int_dp_dx(2) + auxp
            
        end do ! --> do nlz=1,nle-1
        
        !_______________________________________________________________________
        ! calculate pressure gradient for bottom layer
        nlz = nle 
        
        ! loop over the three node points that span up element
        do ni =1, 3
            ! node...
            node = elnodes(ni)
            
            !___________________________________________________________________
            ! calculate vertical center index
            nlc=nln(ni)-1
            !!PS do dd=1,nln(ni)
            do dd=uln(ni),nln(ni)
                if (Z_3d_n(dd,elnodes(ni))<=Z_n(nlz)) then
                    nlc = dd-1
                    if (dd==1) nlc=1
                    exit
                end if 
            end do
            ! --> this is 10% slower than the above one
            ! nlc   = minloc( Z_3d_n(1:nln(ni),elnodes(ni))-Z_n(nlz),1,& 
            !                 Z_3d_n(1:nln(ni),elnodes(ni))-Z_n(nlz)>0.0_WP &  ! mask for index selection
            !                )
            ! nlc   = min(nlc,nln(ni)-1) 
            
            !___________________________________________________________________
            ! "Cubic spline interpolation computes a third order polynomial only 
            ! from two data points with the additional constraint that the first 
            ! and second derivative at the interpolation points are continuous. 
            ! So if you have 4 points, then you compute 3 different polynomials 
            !(between points 1-2, 2-3, and 3-4), and these polynomials are 
            ! smoothly connected in the sense that their first and second 
            ! derivatives are equal at the given data points." --> as a comment
            ! prepare cubic spline interpolation 
            s_ind=(/nlc-1,nlc,nlc+1,nlc+2/)
            
            ! we are at the bottom 
            s_ind(4)=nlc+1
            
            s_z    = Z_3d_n(s_ind,node)
            s_dens = density_m_rho0(s_ind,node)
            s_H    = s_z(3)-s_z(2)
            aux1   = (s_dens(3)-s_dens(2))/s_H
            
            !___________________________________________________________________
            ! calculate derivatives in a way to get monotonic profile 
            ! --> bottom case (see FESOM1.4)
            aux2=(s_dens(2)-s_dens(1))/(s_z(2)-s_z(1))
            s_dup=0.0_WP
            if(aux1*aux2>0._WP)  s_dup=2.0_WP*aux1*aux2/(aux1+aux2)
            s_dlo=1.5_WP*aux1-0.5_WP*s_dup
            
            !___________________________________________________________________
            ! cubic polynomial coefficients
            a=s_dens(2)
            b=s_dup
            c=-(2.0_WP*s_dup+s_dlo)/s_H + 3.0_WP*(s_dens(3)-s_dens(2))/s_H**2
            d=(s_dup+s_dlo)/s_H**2 - 2.0_WP*(s_dens(3)-s_dens(2))/s_H**3
            
            !___________________________________________________________________
            ! interpolate
            dz=Z_n(nlz)-s_z(2)
            interp_n_dens(ni)=a+b*dz+c*dz**2+d*dz**3
            
        end do ! --> do ni=1,3
        
        ! zonal gradients
        drho_dx         = sum(gradient_sca(1:3,elem)*interp_n_dens)
        auxp            = drho_dx*helem(nlz,elem)*g/density_0
        pgf_x(nlz,elem) = int_dp_dx(1) + auxp*0.5_WP
        int_dp_dx(1)    = int_dp_dx(1) + auxp
            
        ! meridional gradients
        drho_dx         = sum(gradient_sca(4:6,elem)*interp_n_dens)
        auxp            = drho_dx*helem(nlz,elem)*g/density_0
        pgf_y(nlz,elem) = int_dp_dx(2) + auxp*0.5_WP
        int_dp_dx(2)    = int_dp_dx(2) + auxp
        
    end do ! --> do elem=1, myDim_elem2D
end subroutine pressure_force_4_linfs_cubicspline
!
!
!
!===============================================================================
! calculate pressure gradient force for linfs in case cavities are used with 
! surface partial cells or bottom partial cells
subroutine pressure_force_4_linfs_cavity(mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    use g_PARSUP
    use g_config
    implicit none
    
    integer             :: elem, elnodes(3), nle, ule,  nlz, idx(3), ni
    real(kind=WP)       :: int_dp_dx(2), drho_dx, dz_dx, aux_sum
    real(kind=WP)       :: dx10(3), dx20(3), dx21(3), df10(3), df21(3), drho_dz(3)
    type(t_mesh), intent(in) , target :: mesh

#include  "associate_mesh.h"
    
    !___________________________________________________________________________
    ! loop over triangular elemments
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        ! number of levels at elem
        nle=nlevels(elem)-1
        ule=ulevels(elem)
            
        !_______________________________________________________________________
        ! node indices of elem 
        elnodes = elem2D_nodes(:,elem)
        
        !_______________________________________________________________________
        ! calculate mid depth element level --> Z_e
        ! nle...number of mid-depth levels at elem
        zbar_n       = 0.0_WP
        Z_n          = 0.0_WP
        zbar_n(nle+1)= zbar_e_bot(elem)
        Z_n(nle)     = zbar_n(nle+1) + helem(nle,elem)*0.5_WP
        do nlz=nle,ule+1,-1
            zbar_n(nlz) = zbar_n(nlz+1) + helem(nlz,elem)
            Z_n(nlz-1)  = zbar_n(nlz)   + helem(nlz-1,elem)*0.5_WP
        end do
        zbar_n(ule) = zbar_n(ule+1) + helem(ule,elem)
        
        !_______________________________________________________________________
        ! calculate pressure gradient for surface layer 
        nlz=ule
        if (ule==1) then 
            pgf_x(nlz,elem) = sum(gradient_sca(1:3,elem)*hpressure(nlz,elnodes)/density_0)
            pgf_y(nlz,elem) = sum(gradient_sca(4:6,elem)*hpressure(nlz,elnodes)/density_0)
        else
            !___________________________________________________________________
            ! compute vertical density gradient 
            !       x0        x1               x2
            !  -----o---------o----------------o-------
            !       f0        f1               f2
            ! --> use Newtonsche interpolation polynom of second order to calculate 
            ! central difference quotient for not equidestant distributed values
            ! f(x) = a0 + a1*(x-x0) + a2*(x-x0)*(x-x1)
            !  | |
            !  | +--> derivative : f'(x) = a1 + a2*[(x-x1)+(x-x0)]
            !  |
            !  +--> determine Coefficients a0, a1, a2 for f(x)=f1,f2,f3
            !  |
            !  +--> a0 = f0 , 
            !       a1 = (f1-f0)/(x1-x0) , 
            !       a2 = ((x1-x0)*(f2-f1)-(x2-x1)*(f1-f0))/((x2-x0)*(x2-x1)*(x1-x0))
            !
            ! f'(x) = (f1-f0)/(x1-x0) 
            !         - [(x1-x0)*(f2-f1)-(x2-x1)*(f1-f0)]/[(x2-x1)*(x1-x0)*(x2-x0)]
            !           *[(x-x1)+(x-x0)]
            !
            nlz = ule
            idx = (/1, 1, 1/)*nlz
            idx = idx - ulevels_nod2D(elnodes)
            do ni=1,3
                if (idx(ni)==0) then
                    ! elemental surface index nlz ist also surface index of vertice ni
                    !_______________________________________________________________
                    ! interpolation coefficients 
                    dx10(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    dx21(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz+1,elnodes(ni))
                    dx20(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    
                    df10(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz+2,elnodes(ni))-density_m_rho0(nlz+1,elnodes(ni))
                    !_______________________________________________________________
                    ! compute vertical density gradient
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                *((Z_n(nlz)-Z_3d_n(nlz+1,elnodes(ni))) + &
                                (Z_n(nlz)-Z_3d_n(nlz  ,elnodes(ni))))
                else
                    ! elemental surface index nlz ist deeper than surface index of vertices ni
                    !_______________________________________________________________
                    ! interpolation coefficients 
                    dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                    dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                    
                    df10(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                    !_______________________________________________________________
                    ! compute vertical density gradient
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                *((Z_n(nlz)-Z_3d_n(nlz  ,elnodes(ni))) + &
                                (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))))
                end if 
            end do
            
            !_______________________________________________________________________
            ! -->...-g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            dz_dx           = sum(gradient_sca(1:3,elem)*Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0               
            pgf_x(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(1)    = aux_sum
             
            !_______________________________________________________________________
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            dz_dx           = sum(gradient_sca(4:6,elem)*Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0             
            pgf_y(nlz,elem) = aux_sum*0.5_WP
            int_dp_dx(2)    = aux_sum
        end if 
        
        !_______________________________________________________________________
        ! calculate pressure gradient for the ocean bulk based on pressure
        do nlz=ule+1,nle-1
            pgf_x(nlz,elem) = sum(gradient_sca(1:3,elem)*hpressure(nlz,elnodes)/density_0)
            pgf_y(nlz,elem) = sum(gradient_sca(4:6,elem)*hpressure(nlz,elnodes)/density_0)
        end do 
        
        !_______________________________________________________________________
        ! calculate pressure gradient for the ocean bottom
        nlz=nle
        if (.not. use_partial_cell) then 
            pgf_x(nlz,elem) = sum(gradient_sca(1:3,elem)*hpressure(nlz,elnodes)/density_0)
            pgf_y(nlz,elem) = sum(gradient_sca(4:6,elem)*hpressure(nlz,elnodes)/density_0)
        else
            !___________________________________________________________________
            ! (1) compute pgf_xy onto zbar_n(nle-1) --> needed for vertical 
            ! integration of bottom density gradient from zbar_n(nle-1) until Z_n(nle)
            int_dp_dx(1) = sum(gradient_sca(1:3,elem)*(hpressure(nlz-1,elnodes) & 
                           + 0.5_WP*g*(density_m_rho0(nlz-1,elnodes)*hnode(nlz-1,elnodes)))/density_0)
            int_dp_dx(2) = sum(gradient_sca(4:6,elem)*(hpressure(nlz-1,elnodes) & 
                           + 0.5_WP*g*(density_m_rho0(nlz-1,elnodes)*hnode(nlz-1,elnodes)))/density_0)
                           
            !___________________________________________________________________
            ! vertical gradient --> with average density and average mid-depth level 
            ! on element
            ! f'(x) = (f1-f0)/(x1-x0) 
            !         - [(x1-x0)*(f2-f1)-(x2-x1)*(f1-f0)]/[(x2-x1)*(x1-x0)*(x2-x0)]
            !           *[(x-x1)+(x-x0)]
            idx = (/1, 1, 1/)*nlz
            idx = nlevels_nod2D(elnodes)-1 - idx
            do ni=1,3
                if (idx(ni)==0) then
                    ! elemental bottom index nlz ist also bottom index of vertice ni
                    !_______________________________________________________________
                    ! interpolation coefficients 
                    dx10(ni)   = Z_3d_n(nlz-1,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                    dx21(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                    dx20(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                    
                    df10(ni)   = density_m_rho0(nlz-1,elnodes(ni))-density_m_rho0(nlz-2,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                    !_______________________________________________________________
                    ! compute vertical density gradient
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                *((Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) + &
                                (Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni))))
                else
                    ! elemental bottom index nlz ist shallower than bottom index of vertice ni
                    !_______________________________________________________________
                    ! interpolation coefficients 
                    dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                    dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                    dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                    
                    df10(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                    df21(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                    !_______________________________________________________________
                    ! compute vertical density gradient
                    drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                                *((Z_n(nlz)-Z_3d_n(nlz  ,elnodes(ni))) + &
                                (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))))
                end if 
            end do
            
            !___________________________________________________________________
            ! -->...-g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            dz_dx           = sum(gradient_sca(1:3,elem)*Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0               
            pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
            
            !_______________________________________________________________________
            ! meridional gradients
            drho_dx         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            dz_dx           = sum(gradient_sca(4:6,elem)*Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0             
            pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
        end if 
        
    end do !-->do elem=1, myDim_elem2D
end subroutine pressure_force_4_linfs_cavity
!
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) for full free surface case zlevel and zstar
subroutine pressure_force_4_zxxxx(mesh)
    use g_PARSUP
    use g_config
    use mod_mesh
    use pressure_force_4_zxxxx_shchepetkin_interface
    use pressure_force_4_zxxxx_cubicspline_interface
    use pressure_force_4_zxxxx_easypgf_interface
    implicit none
    type(t_mesh), intent(in) , target :: mesh    
    !___________________________________________________________________________
    if     (trim(which_pgf)=='shchepetkin') then
        call pressure_force_4_zxxxx_shchepetkin(mesh)
    elseif (trim(which_pgf)=='cubicspline') then
        call pressure_force_4_zxxxx_cubicspline(mesh)
    elseif (trim(which_pgf)=='easypgf'    ) then
        call pressure_force_4_zxxxx_easypgf(mesh)    
    else
        write(*,*) '________________________________________________________'
        write(*,*) ' --> ERROR: the choosen form of pressure gradient       '
        write(*,*) '            calculation (PGF) is not supported for      '
        write(*,*) '            zlevel/zstar !!!                            '
        write(*,*) '            see in namelist.oce --> which_pgf =         '
        write(*,*) '            shchepetkin, cubicspline, easypgf                   '
        write(*,*) '________________________________________________________'
        call par_ex(1)
    end if 
end subroutine pressure_force_4_zxxxx
!
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) on elements by interpolating the nodal 
! density on nodal depthlayers on the depth of the elments using a cubic-spline 
! interpolation.
! First coded by Q. Wang for FESOM1.4, adapted by P. Scholz for FESOM2.0
! 26.04.2018
subroutine pressure_force_4_zxxxx_cubicspline(mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    use g_PARSUP
    use g_config
    implicit none
    
    integer             :: elem, elnodes(3), nle, ule, nln(3), uln(3), nlz, nlc,dd
    integer             :: ni, node, dens_ind,kk
    real(kind=WP)       :: ze
    integer             :: s_ind(4)
    real(kind=WP)       :: s_z(4), s_dens(4), s_H, aux1, aux2, aux(2), s_dup, s_dlo
    real(kind=WP)       :: a, b, c, d, dz, rho_n(3), rhograd_e(2), p_grad(2)
    type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
    !___________________________________________________________________________
    ! loop over triangular elemments
    do elem=1, myDim_elem2D
        ule          = ulevels(elem)
        nle          = nlevels(elem)-1
        !_______________________________________________________________________
        ! calculate mid depth element level --> Z_e
        ! nle...number of mid-depth levels at elem
        zbar_n       = 0.0_WP
        Z_n          = 0.0_WP
        zbar_n(nle+1)= zbar_e_bot(elem)
        Z_n(nle)     = zbar_n(nle+1) + helem(nle,elem)/2.0_WP
        !!PS do nlz=nle,2,-1
        do nlz=nle,ule+1,-1
            zbar_n(nlz) = zbar_n(nlz+1) + helem(nlz,elem)
            Z_n(nlz-1)  = zbar_n(nlz)   + helem(nlz-1,elem)/2.0_WP
        end do
        !!PS zbar_n(1) = zbar_n(2) + helem(1,elem)
        zbar_n(ule) = zbar_n(ule+1) + helem(ule,elem)
        
        !_______________________________________________________________________
        ! node indices of elem 
        elnodes = elem2D_nodes(:,elem)
        
        !_______________________________________________________________________
        ! nln...number of mid depth levels at each node that corresponds to 
        ! element elem
        nln=nlevels_nod2D(elnodes)-1
        uln=ulevels_nod2D(elnodes)
        
        !_______________________________________________________________________
        ! loop over mid-depth levels at element elem
        p_grad=0.0_WP
        !!PS do nlz=1,nle
        do nlz=ule,nle
            
            !___________________________________________________________________
            rho_n = 0.0_WP
            ! loop over the three node points that span up element elem
            do ni=1,3
                ! node...
                node = elnodes(ni)
                
                !_______________________________________________________________
                ! calculate vertical center index 
                nlc=nln(ni)-1
                !!PS do dd=1,nln(ni)
                do dd=uln(ni),nln(ni)
                    if (Z_3d_n(dd,elnodes(ni))<=Z_n(nlz)) then
                        nlc = dd-1
                        if (dd==1) nlc=1
                        exit
                    end if 
                end do
                ! --> this is 10% slower than the above one
                ! nlc   = minloc( Z_3d_n(1:nln(ni),elnodes(ni))-Z_n(nlz),1,& 
                !                 Z_3d_n(1:nln(ni),elnodes(ni))-Z_n(nlz)>0.0_WP)
                ! nlc   = max(nlc,1)                
                ! nlc   = min(nlc,nln(ni)-1) 
                
                !_______________________________________________________________
                ! prepare cubic spline interpolation 
                s_ind=(/nlc-1, nlc, nlc+1, nlc+2/)
                
                !_______________________________________________________________
                ! calculate derivatives in a way to get monotonic profile
                ! distinguish between surface, bottom and bulg
                !!PS if      (nlc==1) then ! surface case
                if      (nlc==uln(ni)) then ! surface case
                    !___________________________________________________________
!!PS                     s_ind(1)=1
                    s_ind(1)=uln(ni)
                    
                    !___________________________________________________________
                    s_z    = Z_3d_n(s_ind,node)
                    s_dens = density_m_rho0(s_ind,node)
                    s_H    = s_z(3)-s_z(2)
                    aux1   = (s_dens(3)-s_dens(2))/s_H
                    
                    !___________________________________________________________
                    aux2=(s_dens(4)-s_dens(3))/(s_z(4)-s_z(3))
                    s_dlo=0.0_WP
                    if(aux1*aux2>0._WP) s_dlo=2.0_WP*aux1*aux2/(aux1+aux2)
                    s_dup=1.5_WP*aux1-0.5_WP*s_dlo
                    
                elseif (nlc == nln(ni)-1) then! bottom case
                    !___________________________________________________________
                    s_ind(4)=nlc+1
                    
                    !___________________________________________________________
                    s_z    = Z_3d_n(s_ind,node)
                    s_dens = density_m_rho0(s_ind,node)
                    s_H    = s_z(3)-s_z(2)
                    aux1   = (s_dens(3)-s_dens(2))/s_H
                    
                    !___________________________________________________________
                    aux2=(s_dens(2)-s_dens(1))/(s_z(2)-s_z(1))
                    s_dup=0.0_WP
                    if(aux1*aux2>0._WP)  s_dup=2.0_WP*aux1*aux2/(aux1+aux2)
                    s_dlo=1.5_WP*aux1-0.5_WP*s_dup
                    
                else ! bulk, subsurface/above bottom case
                    !___________________________________________________________
                    s_z    = Z_3d_n(s_ind,node)
                    s_dens = density_m_rho0(s_ind,node)
                    s_H    = s_z(3)-s_z(2)
                    aux1   = (s_dens(3)-s_dens(2))/s_H
                    
                    !___________________________________________________________
                    aux2=(s_dens(2)-s_dens(1))/(s_z(2)-s_z(1))
                    s_dup=0.0_WP
                    if(aux1*aux2>0._WP)  s_dup=2.0_WP*aux1*aux2/(aux1+aux2)
                    aux2=(s_dens(4)-s_dens(3))/(s_z(4)-s_z(3))
                    s_dlo=0.0_WP
                    if(aux1*aux2>0._WP) s_dlo=2.0_WP*aux1*aux2/(aux1+aux2)
                    
                end if
                
                !_______________________________________________________________
                ! cubic polynomial coefficients
                a=s_dens(2)
                b=s_dup
                c=-(2.0_WP*s_dup+s_dlo)/s_H + 3.0_WP*(s_dens(3)-s_dens(2))/s_H**2
                d=(s_dup+s_dlo)/s_H**2 - 2.0_WP*(s_dens(3)-s_dens(2))/s_H**3
                
                !_______________________________________________________________
                ! interpolate
                dz=Z_n(nlz)-s_z(2)
                rho_n(ni)=a+b*dz+c*dz**2+d*dz**3
                
            end do ! --> do ni=1,3
            
            !___________________________________________________________________
            ! calculate element wise density gradient
            rhograd_e(1) = sum(gradient_sca(1:3,elem)*rho_n)
            rhograd_e(2) = sum(gradient_sca(4:6,elem)*rho_n)
            
            !___________________________________________________________________
            ! calculate element wise pressure gradient force 
            ! helem ... here because of vertical integral 
            aux             = g*helem(nlz,elem)*rhograd_e/density_0
            
            ! *0.5 because pgf_xy is calculated at mid depth levels but at 
            ! this point p_grad is integrated pressure gradient force until
            ! full depth  layers of previouse depth layer
            pgf_x(nlz,elem) = p_grad(1) + aux(1)*0.5_WP
            pgf_y(nlz,elem) = p_grad(2) + aux(2)*0.5_WP
            
            ! integration to full depth levels
            p_grad          = p_grad    + aux
            
        end do ! --> do nlz=1,nle
    end do ! --> do elem=1, myDim_elem2D
end subroutine pressure_force_4_zxxxx_cubicspline
!
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) based on "Shchepetkin 
! and McWilliams," A method for computing horizontal pressure-gradient
! force in an oceanic model with a nonaligned vertical coordinate",
! Journal of geophysical research, vol 108, no C3, 3090
! --> based on density jacobian method ...
! calculate PGF for linfs with partiell cell on/off
! First coded by P. Scholz for FESOM2.0, 08.02.2019
subroutine pressure_force_4_zxxxx_shchepetkin(mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    use g_PARSUP
    use g_config
    use densityJM_components_interface
    use density_linear_interface
    implicit none
    
    integer             :: elem, elnodes(3), nle,ule, nlz, nln(3), ni, nlc, nlce, idx(3)
    real(kind=WP)       :: int_dp_dx(2), drho_dx, drho_dy, drho_dz(3), dz_dx, dz_dy, aux_sum 
    real(kind=WP)       :: dx10(3), dx20(3), dx21(3), df10(3), df21(3)
    real(kind=WP)       :: rhopot(3), bulk_0(3), bulk_pz(3), bulk_pz2(3)
    type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
    
    !___________________________________________________________________________
    ! loop over triangular elemments
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        ! nle...number of mid-depth levels at elem
        nle          = nlevels(elem)-1
        ule          = ulevels(elem)
        
        ! node indices of elem 
        elnodes      = elem2D_nodes(:,elem) 
        
        !_______________________________________________________________________
        ! calculate mid depth element level --> Z_e
        ! nle...number of mid-depth levels at elem
        zbar_n       = 0.0_WP
        Z_n          = 0.0_WP
        zbar_n(nle+1)= zbar_e_bot(elem)
        Z_n(nle)     = zbar_n(nle+1) + helem(nle,elem)*0.5_WP
        do nlz=nle,ule+1,-1
            zbar_n(nlz) = zbar_n(nlz+1) + helem(nlz,elem)
            Z_n(nlz-1)  = zbar_n(nlz)   + helem(nlz-1,elem)*0.5_WP
        end do
        zbar_n(ule) = zbar_n(ule+1) + helem(ule,elem)
        
        !_______________________________________________________________________
        ! Calculate pressure gradient force (PGF) based on "Shchepetkin 
        ! and McWilliams," A method for computing horizontal pressure-gradient
        ! force in an oceanic model with a nonaligned vertical coordinate",
        ! Journal of geophysical research, vol 108, no C3, 3090
        ! --> based on density jacobian method ...
        ! --> equation 1.7
        ! 
        ! -1/rho_0* dP/dx|_z = -1/rho_0*dP/dx|_z=eta - g/rho_0*int_z^eta( drho/dx|_z)*dz'
        !
        !                    = g*rho(eta)/rho_0*deta/dx
        !                      - g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
        !
        ! --> g*rho(eta)/rho_0*deta/dx
        
        !_______________________________________________________________________
        ! vertical gradient --> with average density and average 
        ! mid-depth level on element
        !       x0        x1               x2
        !  -----o---------o----------------o-------
        !       f0        f1               f2
        ! --> use Newtonsche interpolation polynom of second order to calculate 
        ! central difference quotient for not equidestant distributed values
        ! f(x) = a0 + a1*(x-x0) + a2*(x-x0)*(x-x1)
        !  | |
        !  | +--> derivative : f'(x) = a1 + a2*[(x-x1)+(x-x0)]
        !  |
        !  +--> determine Coefficients a0, a1, a2 for f(x)=f1,f2,f3
        !  |
        !  +--> a0 = f0 , 
        !       a1 = (f1-f0)/(x1-x0) , 
        !       a2 = ((x1-x0)*(f2-f1)-(x2-x1)*(f1-f0))/((x2-x0)*(x2-x1)*(x1-x0))
        !
        ! f'(x) = (f1-f0)/(x1-x0) 
        !         - [(x1-x0)*(f2-f1)-(x2-x1)*(f1-f0)]/[(x2-x1)*(x1-x0)*(x2-x0)]
        !           *[(x-x1)+(x-x0)]
        !_______________________________________________________________________
        ! calculate pressure gradient for surface layer 
        nlz = ule
        idx = (/1, 1, 1/)*nlz
        idx = idx - ulevels_nod2D(elnodes)
        do ni=1,3
            if (idx(ni)==0) then
                ! elemental surface index nlz ist also surface index of vertice ni
                !_______________________________________________________________
                ! interpolation coefficients 
                dx10(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz+1,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                
                df10(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                df21(ni)   = density_m_rho0(nlz+2,elnodes(ni))-density_m_rho0(nlz+1,elnodes(ni))
                !_______________________________________________________________
                ! compute vertical density gradient
                drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                             *((Z_n(nlz)-Z_3d_n(nlz+1,elnodes(ni))) + &
                               (Z_n(nlz)-Z_3d_n(nlz  ,elnodes(ni))))
            else
                ! elemental surface index nlz ist deeper than surface index of vertices ni
                !_______________________________________________________________
                ! interpolation coefficients 
                dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                
                df10(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                df21(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                !_______________________________________________________________
                ! compute vertical density gradient
                drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                             *((Z_n(nlz)-Z_3d_n(nlz  ,elnodes(ni))) + &
                               (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))))
            end if 
        end do
        
        !_______________________________________________________________________
        ! -->...-g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
        ! zonal surface pressure gradients
        drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
        dz_dx           = sum(gradient_sca(1:3,elem)*Z_3d_n(nlz,elnodes))
        aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0               
        pgf_x(nlz,elem) = aux_sum*0.5_WP
        int_dp_dx(1)    = aux_sum
        
        !_______________________________________________________________________
        ! meridional surface pressure gradients
        drho_dy         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
        dz_dy           = sum(gradient_sca(4:6,elem)*Z_3d_n(nlz,elnodes))
        aux_sum         = (drho_dy-sum(drho_dz)/3.0_WP*dz_dy)*helem(nlz,elem)*g/density_0             
        pgf_y(nlz,elem) = aux_sum*0.5_WP
        int_dp_dx(2)    = aux_sum
        
        !_______________________________________________________________________
        ! calculate pressure gradient for subsurface layer until one layer above 
        ! the bottom 
        do nlz=ule+1,nle-1
            !___________________________________________________________________
            ! compute coefficients
            dx10            = Z_3d_n(nlz  ,elnodes)-Z_3d_n(nlz-1,elnodes)
            dx21            = Z_3d_n(nlz+1,elnodes)-Z_3d_n(nlz  ,elnodes)
            dx20            = Z_3d_n(nlz+1,elnodes)-Z_3d_n(nlz-1,elnodes)
            
            df10            = density_m_rho0(nlz  ,elnodes)-density_m_rho0(nlz-1,elnodes)
            df21            = density_m_rho0(nlz+1,elnodes)-density_m_rho0(nlz  ,elnodes)
            !_______________________________________________________________
            ! compute vertical density gradient
            drho_dz         = df10/dx10 + (dx10*df21-dx21*df10)/(dx20*dx21*dx10)&
                              *(((/1.0,1.0,1.0/)*Z_n(nlz)-Z_3d_n(nlz,elnodes)) + &
                                ((/1.0,1.0,1.0/)*Z_n(nlz)-Z_3d_n(nlz-1,elnodes)))
            
            !___________________________________________________________________
            ! -->...-g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
            ! zonal bulk pressure gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
            dz_dx           = sum(gradient_sca(1:3,elem)*Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0               
            pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
            int_dp_dx(1)    = int_dp_dx(1) + aux_sum
            
            !___________________________________________________________________
            ! meridional bulk pressure gradients
            drho_dy         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
            dz_dy           = sum(gradient_sca(4:6,elem)*Z_3d_n(nlz,elnodes))
            aux_sum         = (drho_dy-sum(drho_dz)/3.0_WP*dz_dy)*helem(nlz,elem)*g/density_0             
            pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
            int_dp_dx(2)    = int_dp_dx(2) + aux_sum
            
        end do ! --> do nlz=2,nle-1
        
        !_______________________________________________________________________
        ! calculate pressure gradient for bottom layer
        nlz = nle 
        idx = (/1, 1, 1/)*nlz
        idx = nlevels_nod2D(elnodes)-1 - idx
        do ni=1,3
            if (idx(ni)==0) then
                ! elemental bottom index nlz ist also bottom index of vertice ni
                !_______________________________________________________________
                ! interpolation coefficients 
                dx10(ni)   = Z_3d_n(nlz-1,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                
                df10(ni)   = density_m_rho0(nlz-1,elnodes(ni))-density_m_rho0(nlz-2,elnodes(ni))
                df21(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                !_______________________________________________________________
                ! compute vertical density gradient
                drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                             *((Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) + &
                               (Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni))))
            else
                ! elemental bottom index nlz ist shallower than bottom index of vertice ni
                !_______________________________________________________________
                ! interpolation coefficients 
                dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                
                df10(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                df21(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                !_______________________________________________________________
                ! compute vertical density gradient
                drho_dz(ni)= df10(ni)/dx10(ni) + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))&
                             *((Z_n(nlz)-Z_3d_n(nlz  ,elnodes(ni))) + &
                               (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))))
            end if 
        end do
        
        !_______________________________________________________________________
        ! -->...-g/rho*int_z^eta( drho/dx|_s - drho/dz'*dz'/dx|_s )*dz'
        ! zonal bottom pressure gradients
        drho_dx         = sum(gradient_sca(1:3,elem)*density_m_rho0(nlz,elnodes))
        dz_dx           = sum(gradient_sca(1:3,elem)*(Z_3d_n(nlz,elnodes)))
        aux_sum         = (drho_dx-sum(drho_dz)/3.0_WP*dz_dx)*helem(nlz,elem)*g/density_0
        pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
        
        !_______________________________________________________________________
        ! meridional bottom pressure gradients
        drho_dy         = sum(gradient_sca(4:6,elem)*density_m_rho0(nlz,elnodes))
        dz_dy           = sum(gradient_sca(4:6,elem)*(Z_3d_n(nlz,elnodes)))
        aux_sum         = (drho_dy-sum(drho_dz)/3.0_WP*dz_dy)*helem(nlz,elem)*g/density_0
        pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
        
    end do ! --> do elem=1, myDim_elem2D
    
end subroutine pressure_force_4_zxxxx_shchepetkin
!
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) based on "Shchepetkin 
! and McWilliams," A method for computing horizontal pressure-gradient
! force in an oceanic model with a nonaligned vertical coordinate",
! Journal of geophysical research, vol 108, no C3, 3090
! --> based on density jacobian method ...
! calculate PGF for linfs with partiell cell on/off
! First coded by P. Scholz for FESOM2.0, 08.02.2019
subroutine pressure_force_4_zxxxx_easypgf(mesh)
    use o_PARAM
    use MOD_MESH
    use o_ARRAYS
    use g_PARSUP
    use g_config
    use densityJM_components_interface
    use density_linear_interface
    implicit none
    
    integer             :: elem, elnodes(3), nle,ule, nlz, nln(3), ni, nlc, nlce, idx(3)
    real(kind=WP)       :: int_dp_dx(2), drho_dx, dz_dx, drho_dy, dz_dy,aux_sum 
    real(kind=WP)       :: dx10(3), dx20(3), dx21(3)
    real(kind=WP)       :: f0(3), df10(3), df21(3)
    real(kind=WP)       :: t0(3), dt10(3), dt21(3)
    real(kind=WP)       :: s0(3), ds10(3), ds21(3)
    real(kind=WP)       :: rho_at_Zn(3), temp_at_Zn(3), salt_at_Zn(3), drho_dz(3), aux_dref
    real(kind=WP)       :: rhopot(3), bulk_0(3), bulk_pz(3), bulk_pz2(3)
    real(kind=WP)       :: dref_rhopot, dref_bulk_0, dref_bulk_pz, dref_bulk_pz2
    type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
    !___________________________________________________________________________
    ! loop over triangular elemments
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        ! nle...number of mid-depth levels at elem
        nle          = nlevels(elem)-1
        ule          = ulevels(elem)
        
        ! node indices of elem 
        elnodes      = elem2D_nodes(:,elem) 
        
        !_______________________________________________________________________
        ! calculate mid depth element level --> Z_e
        ! nle...number of mid-depth levels at elem
        zbar_n       = 0.0_WP
        Z_n          = 0.0_WP
        zbar_n(nle+1)= zbar_e_bot(elem)
        Z_n(nle)     = zbar_n(nle+1) + helem(nle,elem)*0.5_WP
        do nlz=nle,ule+1,-1
            zbar_n(nlz) = zbar_n(nlz+1) + helem(nlz,elem)
            Z_n(nlz-1)  = zbar_n(nlz)   + helem(nlz-1,elem)*0.5_WP
        end do
        zbar_n(ule)  = zbar_n(ule+1) + helem(ule,elem)
        
        !_______________________________________________________________________
        aux_dref     = density_0
        if (use_cavity .and. .not. use_density_ref) then
            select case(state_equation)
                case(0)
                    call density_linear(density_ref_T, density_ref_S, dref_bulk_0, dref_bulk_pz, dref_bulk_pz2, dref_rhopot, mesh)
                case(1)
                    call densityJM_components(density_ref_T, density_ref_S, dref_bulk_0, dref_bulk_pz, dref_bulk_pz2, dref_rhopot, mesh)
                case default !unknown
                    if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                    call par_ex(1)
            end select
        end if     
        
        !_______________________________________________________________________
        ! vertical gradient --> with average density and average 
        ! mid-depth level on element
        !       x0        x1               x2
        !  -----o---------o----------------o-------
        !       f0        f1               f2
        ! --> use Newtonsche interpolation polynom of second order to calculate 
        ! central difference quotient for not equidestant distributed values
        ! f(x) = a0 + a1*(x-x0) + a2*(x-x0)*(x-x1)
        !  | |
        !  | +--> derivative : f'(x) = a1 + a2*[(x-x1)+(x-x0)]
        !  |
        !  +--> determine Coefficients a0, a1, a2 for f(x)=f1,f2,f3
        !  |
        !  +--> a0 = f0 , 
        !       a1 = (f1-f0)/(x1-x0) , 
        !       a2 = ((x1-x0)*(f2-f1)-(x2-x1)*(f1-f0))/((x2-x0)*(x2-x1)*(x1-x0))
        !
        ! f(x) = f0 
        !        + (f1-f0)/(x1-x0)*(x-x0) 
        !        + [(x1-x0)*(f2-f1)-(x2-x1)*(f1-f0)]/[(x2-x0)*(x2-x1)*(x1-x0)]
        !          *[(x-x1)+(x-x0)]
        !
        !_______________________________________________________________________
        ! calculate pressure gradient for surface layer 
        
        ! account for reference density when using cavities
        if (use_cavity .and. .not. use_density_ref) then
            aux_dref = dref_bulk_0   + Z_n(nlz)*dref_bulk_pz   + Z_n(nlz)*dref_bulk_pz2
            aux_dref = aux_dref*dref_rhopot/(aux_dref+0.1_WP*Z_n(nlz))
        end if 
        
        nlz = ule
        idx = (/1, 1, 1/)*nlz
        idx = idx - ulevels_nod2D(elnodes)
        do ni=1,3
            if (idx(ni)==0) then
                ! elemental surface index nlz ist also surface index of vertice ni
                ! --> do second order newtonian surface boundary interpolation
                !___________________________________________________________________
                ! compute interpolation coefficients
                dx10(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz+1,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz+2,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                
                !!PS f0(ni)     = density_m_rho0(nlz  ,elnodes(ni))
                !!PS df10(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                !!PS df21(ni)   = density_m_rho0(nlz+2,elnodes(ni))-density_m_rho0(nlz+1,elnodes(ni))
                
                t0(ni)     = tr_arr(nlz  ,elnodes(ni),1)
                dt10(ni)   = tr_arr(nlz+1,elnodes(ni),1)-tr_arr(nlz  ,elnodes(ni),1)
                dt21(ni)   = tr_arr(nlz+2,elnodes(ni),1)-tr_arr(nlz+1,elnodes(ni),1)
                
                s0(ni)     = tr_arr(nlz  ,elnodes(ni),2)
                ds10(ni)   = tr_arr(nlz+1,elnodes(ni),2)-tr_arr(nlz  ,elnodes(ni),2)
                ds21(ni)   = tr_arr(nlz+2,elnodes(ni),2)-tr_arr(nlz+1,elnodes(ni),2)
                !___________________________________________________________________
                ! interpoalte vertice temp and salinity to elemental level Z_n
                temp_at_Zn(ni) = t0(ni) &
                                + dt10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni))) & 
                                + (dx10(ni)*dt21(ni)-dx21(ni)*dt10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz+1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))
                salt_at_Zn(ni) = s0(ni) &
                                + ds10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni))) & 
                                + (dx10(ni)*ds21(ni)-dx21(ni)*ds10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz+1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))
                !___________________________________________________________________
                ! compute density from state equation 
                select case(state_equation)
                    case(0)
                        call density_linear(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case(1)
                        call densityJM_components(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case default !unknown
                        if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                        call par_ex(1)
                end select
                rho_at_Zn(ni) = bulk_0(ni) + Z_n(nlz)*(bulk_pz(ni) + Z_n(nlz)*bulk_pz2(ni))
                rho_at_Zn(ni) = rho_at_Zn(ni)*rhopot(ni)/(rho_at_Zn(ni)+0.1_WP*Z_n(nlz)*real(state_equation))-aux_dref
                
                !!PS rho_at_Zn(ni) = f0(ni) &
                !!PS                 + df10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni))) & 
                !!PS                 + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                !!PS                   (Z_n(nlz)-Z_3d_n(nlz+1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))
                
            elseif (idx(ni)>0) then
                ! elemental surface index nlz ist deeper than surface index of vertices ni
                ! --> do normal interpolation
                !___________________________________________________________________
                ! compute interpolation coefficients
                dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                
                !!PS f0(ni)     = density_m_rho0(nlz-1,elnodes(ni))
                !!PS df10(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                !!PS df21(ni)   = density_m_rho0(nlz+1,elnodes(ni))-density_m_rho0(nlz  ,elnodes(ni))
                
                t0(ni)     = tr_arr(nlz-1,elnodes(ni),1)
                dt10(ni)   = tr_arr(nlz  ,elnodes(ni),1)-tr_arr(nlz-1,elnodes(ni),1)
                dt21(ni)   = tr_arr(nlz+1,elnodes(ni),1)-tr_arr(nlz  ,elnodes(ni),1)
                
                s0(ni)     = tr_arr(nlz-1,elnodes(ni),2)
                ds10(ni)   = tr_arr(nlz  ,elnodes(ni),2)-tr_arr(nlz-1,elnodes(ni),2)
                ds21(ni)   = tr_arr(nlz+1,elnodes(ni),2)-tr_arr(nlz  ,elnodes(ni),2)
                !___________________________________________________________________               
                ! interpoalte vertice temp and salinity to elemental level Z_n
                temp_at_Zn(ni) = t0(ni) &
                                + dt10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                                + (dx10(ni)*dt21(ni)-dx21(ni)*dt10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
                salt_at_Zn(ni) = s0(ni) &
                                + ds10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                                + (dx10(ni)*ds21(ni)-dx21(ni)*ds10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
                !___________________________________________________________________
                ! compute density from state equation 
                select case(state_equation)
                    case(0)
                        call density_linear(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case(1)
                        call densityJM_components(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case default !unknown
                        if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                        call par_ex(1)
                end select
                rho_at_Zn(ni) = bulk_0(ni) + Z_n(nlz)*(bulk_pz(ni) + Z_n(nlz)*bulk_pz2(ni))
                rho_at_Zn(ni) = rho_at_Zn(ni)*rhopot(ni)/(rho_at_Zn(ni)+0.1_WP*Z_n(nlz)*real(state_equation))-aux_dref
                
                !!PS rho_at_Zn(ni) = f0(ni) &
                !!PS                 + df10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                !!PS                 + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                !!PS                   (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
            else
                ! would do second order Newtonian surface boundary extrapolation 
                ! --> this is not wanted !!!
                write(*,*) ' --> would do second order surface boundary density extrapolation'
                write(*,*) '     This is not wanted, model stops here'
                call par_ex(0)    
            end if 
        end do
        !_______________________________________________________________________
        ! zonal surface pressure gradients
        drho_dx         = sum(gradient_sca(1:3,elem)*rho_at_Zn)
        aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
        pgf_x(nlz,elem) = aux_sum*0.5_WP
        int_dp_dx(1)    = aux_sum
        !_______________________________________________________________________
        ! meridional surface pressure gradients
        drho_dy         = sum(gradient_sca(4:6,elem)*rho_at_Zn)
        aux_sum         = drho_dy*helem(nlz,elem)*g/density_0
        pgf_y(nlz,elem) = aux_sum*0.5_WP
        int_dp_dx(2)    = aux_sum
        
        !_______________________________________________________________________
        ! calculate pressure gradient for subsurface layer until one layer above 
        ! the bottom 
        do nlz=ule+1,nle-1
        
            !___________________________________________________________________
            ! elemental bottom index nlz ist shallower than bottom index of vertice ni
            if (use_cavity .and. .not. use_density_ref) then
                aux_dref = dref_bulk_0   + Z_n(nlz)*dref_bulk_pz   + Z_n(nlz)*dref_bulk_pz2
                aux_dref = aux_dref*dref_rhopot/(aux_dref+0.1_WP*Z_n(nlz))
            end if 
            
            !___________________________________________________________________
            ! compute interpolation coefficients
            dx10   = Z_3d_n(nlz  ,elnodes)-Z_3d_n(nlz-1,elnodes)
            dx21   = Z_3d_n(nlz+1,elnodes)-Z_3d_n(nlz  ,elnodes)
            dx20   = Z_3d_n(nlz+1,elnodes)-Z_3d_n(nlz-1,elnodes)
            
            !!PS f0     = density_m_rho0(nlz-1,elnodes)
            !!PS df10   = density_m_rho0(nlz  ,elnodes)-density_m_rho0(nlz-1,elnodes)
            !!PS df21   = density_m_rho0(nlz+1,elnodes)-density_m_rho0(nlz  ,elnodes)
            
            t0     = tr_arr(nlz-1,elnodes,1)
            dt10   = tr_arr(nlz  ,elnodes,1)-tr_arr(nlz-1,elnodes,1)
            dt21   = tr_arr(nlz+1,elnodes,1)-tr_arr(nlz  ,elnodes,1)
            
            s0     = tr_arr(nlz-1,elnodes,2)
            ds10   = tr_arr(nlz  ,elnodes,2)-tr_arr(nlz-1,elnodes,2)
            ds21   = tr_arr(nlz+1,elnodes,2)-tr_arr(nlz  ,elnodes,2)
            !___________________________________________________________________               
            ! interpoalte vertice temp and salinity to elemental level Z_n
            temp_at_Zn = t0 &
                            + dt10/dx10*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes)) & 
                            + (dx10*dt21-dx21*dt10)/(dx20*dx21*dx10)* &
                                ((/1.0_WP,1.0_WP,1.0_WP/)*Z_n(nlz)-Z_3d_n(nlz,elnodes))*&
                                ((/1.0_WP,1.0_WP,1.0_WP/)*Z_n(nlz)-Z_3d_n(nlz-1,elnodes))
            salt_at_Zn = s0 &
                            + ds10/dx10*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes)) & 
                            + (dx10*ds21-dx21*ds10)/(dx20*dx21*dx10)* &
                                ((/1.0_WP,1.0_WP,1.0_WP/)*Z_n(nlz)-Z_3d_n(nlz,elnodes))*&
                                ((/1.0_WP,1.0_WP,1.0_WP/)*Z_n(nlz)-Z_3d_n(nlz-1,elnodes))
            !___________________________________________________________________
            ! compute density from state equation 
            select case(state_equation)
                case(0)
                    call density_linear(temp_at_Zn(1), salt_at_Zn(1), bulk_0(1), bulk_pz(1), bulk_pz2(1), rhopot(1), mesh)             
                    call density_linear(temp_at_Zn(2), salt_at_Zn(2), bulk_0(2), bulk_pz(2), bulk_pz2(2), rhopot(2), mesh)             
                    call density_linear(temp_at_Zn(3), salt_at_Zn(3), bulk_0(3), bulk_pz(3), bulk_pz2(3), rhopot(3), mesh)             
                case(1)
                    call densityJM_components(temp_at_Zn(1), salt_at_Zn(1), bulk_0(1), bulk_pz(1), bulk_pz2(1), rhopot(1), mesh)             
                    call densityJM_components(temp_at_Zn(2), salt_at_Zn(2), bulk_0(2), bulk_pz(2), bulk_pz2(2), rhopot(2), mesh)             
                    call densityJM_components(temp_at_Zn(3), salt_at_Zn(3), bulk_0(3), bulk_pz(3), bulk_pz2(3), rhopot(3), mesh)             
                case default !unknown
                    if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                    call par_ex(1)
            end select
            rho_at_Zn = bulk_0 + Z_n(nlz)*(bulk_pz + Z_n(nlz)*bulk_pz2)
            rho_at_Zn = rho_at_Zn*rhopot/(rho_at_Zn+0.1_WP*Z_n(nlz)*real(state_equation))-aux_dref
            
            !!PS rho_at_Zn = f0 &
            !!PS                 + df10/dx10*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes)) & 
            !!PS                 + (dx10*df21-dx21*df10)/(dx20*dx21*dx10)* &
            !!PS                     (Z_n(nlz)-Z_3d_n(nlz,elnodes))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes))
            
            !___________________________________________________________________
            ! zonal gradients
            drho_dx         = sum(gradient_sca(1:3,elem)*rho_at_Zn)
            aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
            pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
            int_dp_dx(1)    = int_dp_dx(1) + aux_sum
            
            !___________________________________________________________________
            ! meridional gradients
            drho_dy         = sum(gradient_sca(4:6,elem)*rho_at_Zn)
            aux_sum         = drho_dy*helem(nlz,elem)*g/density_0
            pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
            int_dp_dx(2)    = int_dp_dx(2) + aux_sum
        end do ! --> do nlz=2,nle-1
        
        !_______________________________________________________________________
        ! calculate pressure gradient for bottom layer
        ! vertical gradient --> with average density and average mid-depth level 
        ! on element
        ! f(x) = a0 + a1(x-x0) + a2*(x-x1)(x-x0)
        ! f(x) = f0 + df10/dx10*(x-x0) + (dx10*df21-dx21*df10)/(dx20*dx21*dx10)*(x-x1)(x-x0)
        
        ! account for reference density when using cavities
        if (use_cavity .and. .not. use_density_ref) then
            aux_dref = dref_bulk_0 + Z_n(nlz)*dref_bulk_pz + Z_n(nlz)*dref_bulk_pz2
            aux_dref = aux_dref*dref_rhopot/(aux_dref+0.1_WP*Z_n(nlz))
        end if 
        
        nlz = nle 
        idx = (/1, 1, 1/)*nlz
        idx = nlevels_nod2D(elnodes)-1 - idx
        do ni=1,3
            if (idx(ni)==0) then
                ! elemental bottom index nlz ist also bottom index of vertice ni
                ! --> do second order newtonian bottom boundary interpolation 
                !___________________________________________________________________
                ! compute interpolation coefficients
                dx10(ni)   = Z_3d_n(nlz-1,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-2,elnodes(ni))
                
                !!PS f0(ni)     = density_m_rho0(nlz-2,elnodes(ni))
                !!PS df10(ni)   = density_m_rho0(nlz-1,elnodes(ni))-density_m_rho0(nlz-2,elnodes(ni))
                !!PS df21(ni)   = density_m_rho0(nlz  ,elnodes(ni))-density_m_rho0(nlz-1,elnodes(ni))
                
                t0(ni)     = tr_arr(nlz-2,elnodes(ni),1)
                dt10(ni)   = tr_arr(nlz-1,elnodes(ni),1)-tr_arr(nlz-2,elnodes(ni),1)
                dt21(ni)   = tr_arr(nlz  ,elnodes(ni),1)-tr_arr(nlz-1,elnodes(ni),1)
                
                s0(ni)     = tr_arr(nlz-2,elnodes(ni),2)
                ds10(ni)   = tr_arr(nlz-1,elnodes(ni),2)-tr_arr(nlz-2,elnodes(ni),2)
                ds21(ni)   = tr_arr(nlz  ,elnodes(ni),2)-tr_arr(nlz-1,elnodes(ni),2)
                !___________________________________________________________________
                ! interpoalte vertice temp and salinity to elemental level Z_n
                temp_at_Zn(ni) = t0(ni) &
                                + dt10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni))) & 
                                + (dx10(ni)*dt21(ni)-dx21(ni)*dt10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni)))
                salt_at_Zn(ni) = s0(ni) &
                                + ds10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni))) & 
                                + (dx10(ni)*ds21(ni)-dx21(ni)*ds10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni)))
                !___________________________________________________________________
                ! compute density from state equation 
                select case(state_equation)
                    case(0)
                        call density_linear(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case(1)
                        call densityJM_components(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case default !unknown
                        if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                        call par_ex(1)
                end select
                rho_at_Zn(ni) = bulk_0(ni) + Z_n(nlz)*(bulk_pz(ni) + Z_n(nlz)*bulk_pz2(ni))
                rho_at_Zn(ni) = rho_at_Zn(ni)*rhopot(ni)/(rho_at_Zn(ni)+0.1_WP*Z_n(nlz)*real(state_equation))-aux_dref
                
                !!PS rho_at_Zn(ni) = f0(ni) &
                !!PS                 + df10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni))) & 
                !!PS                 + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                !!PS                  (Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-2,elnodes(ni)))
                
            elseif (idx(ni)>0) then
                ! elemental bottom index nlz ist shallower than bottom index of vertice ni
                ! --> do normal interpolation
                !___________________________________________________________________
                ! compute interpolation coefficients
                dx10(ni)   = Z_3d_n(nlz  ,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                dx21(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz  ,elnodes(ni))
                dx20(ni)   = Z_3d_n(nlz+1,elnodes(ni))-Z_3d_n(nlz-1,elnodes(ni))
                
                !!PS f0     = density_m_rho0(nlz-1,elnodes)
                !!PS df10   = density_m_rho0(nlz  ,elnodes)-density_m_rho0(nlz-1,elnodes)
                !!PS df21   = density_m_rho0(nlz+1,elnodes)-density_m_rho0(nlz  ,elnodes)
            
                t0(ni)     = tr_arr(nlz-1,elnodes(ni),1)
                dt10(ni)   = tr_arr(nlz  ,elnodes(ni),1)-tr_arr(nlz-1,elnodes(ni),1)
                dt21(ni)   = tr_arr(nlz+1,elnodes(ni),1)-tr_arr(nlz  ,elnodes(ni),1)
                
                s0(ni)     = tr_arr(nlz-1,elnodes(ni),2)
                ds10(ni)   = tr_arr(nlz  ,elnodes(ni),2)-tr_arr(nlz-1,elnodes(ni),2)
                ds21(ni)   = tr_arr(nlz+1,elnodes(ni),2)-tr_arr(nlz  ,elnodes(ni),2)
                !___________________________________________________________________               
                ! interpoalte vertice temp and salinity to elemental level Z_n
                temp_at_Zn(ni) = t0(ni) &
                                + dt10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                                + (dx10(ni)*dt21(ni)-dx21(ni)*dt10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
                salt_at_Zn(ni) = s0(ni) &
                                + ds10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                                + (dx10(ni)*ds21(ni)-dx21(ni)*ds10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                                  (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
                !___________________________________________________________________
                ! compute density from state equation 
                select case(state_equation)
                    case(0)
                        call density_linear(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case(1)
                        call densityJM_components(temp_at_Zn(ni), salt_at_Zn(ni), bulk_0(ni), bulk_pz(ni), bulk_pz2(ni), rhopot(ni), mesh)             
                    case default !unknown
                        if (mype==0) write(*,*) 'Wrong type of the equation of state. Check your namelists.'
                        call par_ex(1)
                end select
                rho_at_Zn(ni) = bulk_0(ni) + Z_n(nlz)*(bulk_pz(ni) + Z_n(nlz)*bulk_pz2(ni))
                rho_at_Zn(ni) = rho_at_Zn(ni)*rhopot(ni)/(rho_at_Zn(ni)+0.1_WP*Z_n(nlz)*real(state_equation))-aux_dref
                
                !!PS rho_at_Zn(ni) = f0(ni) &
                !!PS                 + df10(ni)/dx10(ni)*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni))) & 
                !!PS                 + (dx10(ni)*df21(ni)-dx21(ni)*df10(ni))/(dx20(ni)*dx21(ni)*dx10(ni))* &
                !!PS                   (Z_n(nlz)-Z_3d_n(nlz,elnodes(ni)))*(Z_n(nlz)-Z_3d_n(nlz-1,elnodes(ni)))
            else
                ! would do second order Newtonian bottom boundary extrapolation 
                ! --> this is not wanted !!!
                write(*,*) ' --> would do second order bottom boundary density extrapolation'
                write(*,*) '     This is not wanted, model stops here'
                write(*,*) ' idx = ', idx
                write(*,*) ' nle = ', nle
                write(*,*) ' nln = ', nlevels_nod2D(elnodes)-1
                call par_ex(0)
            end if 
        end do
        !_______________________________________________________________________
        ! zonal bottom pressure gradients
        drho_dx         = sum(gradient_sca(1:3,elem)*rho_at_Zn)
        aux_sum         = drho_dx*helem(nlz,elem)*g/density_0
        pgf_x(nlz,elem) = int_dp_dx(1) + aux_sum*0.5_WP
        
        !_______________________________________________________________________
        ! meridional bottom pressure gradients
        drho_dy         = sum(gradient_sca(4:6,elem)*rho_at_Zn)
        aux_sum         = drho_dy*helem(nlz,elem)*g/density_0
        pgf_y(nlz,elem) = int_dp_dx(2) + aux_sum*0.5_WP
       
    end do ! --> do elem=1, myDim_elem2D
end subroutine pressure_force_4_zxxxx_easypgf
!
!
!
!===============================================================================
SUBROUTINE densityJM_local(t, s, pz, rho_out, mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP !, only: par_ex,pe_status
use densityJM_components_interface
IMPLICIT NONE

  !
  ! - calculates in-situ density as a function of potential temperature
  !   (relative to the surface)
  !   using the Jackett and McDougall equation of state
  !   (Copyright (c) 1992, CSIRO, Australia)
  ! - has been derived from the SPEM subroutine rhocal
  !
  !---------------------------------------------------------------------------

  real(kind=WP), intent(IN)  :: t,s,pz
  real(kind=WP), intent(OUT) :: rho_out                 
  real(kind=WP)              :: rhopot, bulk
  real(kind=WP)              :: bulk_0, bulk_pz, bulk_pz2
  type(t_mesh), intent(in)   , target :: mesh
#include "associate_mesh.h"
  !compute secant bulk modulus

  call densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)

  bulk = bulk_0 + pz*(bulk_pz + pz*bulk_pz2) 

  rho_out = bulk*rhopot / (bulk + 0.1_WP*pz) - density_0
end subroutine densityJM_local
!
!
!
!===============================================================================
SUBROUTINE densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP !, only: par_ex,pe_status
IMPLICIT NONE

  !
  ! - calculates in-situ density as a function of potential temperature
  !   (relative to the surface)
  !   using the Jackett and McDougall equation of state
  !   (Copyright (c) 1992, CSIRO, Australia)
  ! - has been derived from the SPEM subroutine rhocal
  !
  ! Ralph Timmermann, August 2005
  !---------------------------------------------------------------------------
  ! N. Rakowski 2014 the split form
  !---------------------------------------------------------------------------
  real(kind=WP), intent(IN)  :: t,s
  real(kind=WP), intent(OUT) :: bulk_0, bulk_pz, bulk_pz2, rhopot
  real(kind=WP)              :: s_sqrt

  real(kind=WP), parameter   :: a0    = 19092.56,     at   = 209.8925
  real(kind=WP), parameter   :: at2   = -3.041638,    at3  = -1.852732e-3
  real(kind=WP), parameter   :: at4   = -1.361629e-5
  real(kind=WP), parameter   :: as    = 104.4077,     ast  = -6.500517
  real(kind=WP), parameter   :: ast2  = .1553190,     ast3 = 2.326469e-4 
  real(kind=WP), parameter   :: ass   = -5.587545,    asst = 0.7390729 
  real(kind=WP), parameter   :: asst2 = -1.909078e-2
  real(kind=WP), parameter   :: ap    = -4.721788e-1, apt  = -1.028859e-2
  real(kind=WP), parameter   :: apt2  = 2.512549e-4,  apt3 = 5.939910e-7 
  real(kind=WP), parameter   :: aps   = 1.571896e-2,  apst = 2.598241e-4 
  real(kind=WP), parameter   :: apst2 = -7.267926e-6, apss = -2.042967e-3
  real(kind=WP), parameter   :: ap2   = 1.045941e-5,  ap2t = -5.782165e-10 
  real(kind=WP), parameter   :: ap2t2 = 1.296821e-7
  real(kind=WP), parameter   :: ap2s  = -2.595994e-7,ap2st = -1.248266e-9 
  real(kind=WP), parameter   :: ap2st2= -3.508914e-9
  
  real(kind=WP), parameter   :: b0 = 999.842594,    bt  = 6.793952e-2
  real(kind=WP), parameter   :: bt2 = -9.095290e-3, bt3 = 1.001685e-4
  real(kind=WP), parameter   :: bt4 = -1.120083e-6, bt5 = 6.536332e-9
  real(kind=WP), parameter   :: bs = 0.824493,      bst = -4.08990e-3
  real(kind=WP), parameter   :: bst2 = 7.64380e-5,  bst3 = -8.24670e-7		
  real(kind=WP), parameter   :: bst4 = 5.38750e-9	
  real(kind=WP), parameter   :: bss = -5.72466e-3,  bsst = 1.02270e-4
  real(kind=WP), parameter   :: bsst2 = -1.65460e-6,bss2 = 4.8314e-4
  type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

  !compute secant bulk modulus

  s_sqrt = sqrt(s)

  bulk_0 =  a0      + t*(at   + t*(at2  + t*(at3 + t*at4)))      &
          + s* (as  + t*(ast  + t*(ast2 + t*ast3))               &
               + s_sqrt*(ass  + t*(asst + t*asst2)))                

  bulk_pz =  ap  + t*(apt  + t*(apt2 + t*apt3))                  &
                  + s*(aps + t*(apst + t*apst2) + s_sqrt*apss)

  bulk_pz2 = ap2 + t*(ap2t + t*ap2t2)		                 &
                + s *(ap2s + t*(ap2st + t*ap2st2))

  rhopot =  b0 + t*(bt + t*(bt2 + t*(bt3  + t*(bt4  + t*bt5))))	 &
               + s*(bs + t*(bst + t*(bst2 + t*(bst3 + t*bst4)))  &
                  + s_sqrt*(bss + t*(bsst + t*bsst2))            &
                       + s* bss2)
end subroutine densityJM_components
!
!
!
!===============================================================================
function ptheta(s,t,p,pr)
  ! Compute local potential temperature at pr
  ! using bryden 1973 polynomial for adiabatic lapse rate
  ! and runge-kutta 4-th order integration algorithm.
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! fofonoff,n.,1977,deep-sea res.,24,489-491
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       reference prs   pr       decibars
  !       potential tmp.  theta    deg celsius
  ! checkvalue: theta= 36.89073 c,s=40 (ipss-78),t=40 deg c,
  ! p=10000 decibars,pr=0 decibars
  !
  ! Coded by ??
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_param, only: WP
  implicit none
  real(kind=WP) 			:: ptheta, s, t, p, pr
  real(kind=WP) 			:: h, xk, q
  real(kind=WP), external	        :: atg

  h = pr - p
  xk = h*atg(s,t,p)
  t = t + 0.5_WP*xk
  q = xk
  p = p + 0.5_WP*h
  xk = h*atg(s,t,p)
  t = t + 0.29289322_WP*(xk-q)
  q = 0.58578644_WP*xk + 0.121320344_WP*q
  xk = h*atg(s,t,p)
  t = t + 1.707106781_WP*(xk-q)
  q = 3.414213562_WP*xk - 4.121320344_WP*q
  p = p + 0.5_WP*h
  xk = h*atg(s,t,p)
  ptheta = t + (xk-2.0_WP*q)/6.0_WP
  return
end function ptheta
!
!
!
!===============================================================================
function atg(s,t,p)
  ! adiabatic temperature gradient deg c per decibar
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       adiabatic       atg      deg. c/decibar
  ! checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
  ! t=40 deg c,p0=10000 decibars
  !
  ! Coded by ??
  ! Reviewed by ??
  !--------------------------------------------------------

  use o_param, only: WP
  implicit none
  real(kind=WP)  atg, s, t, p, ds

  ds = s - 35.0_WP
  atg = (((-2.1687e-16_WP*t+1.8676e-14_WP)*t-4.6206e-13_WP)*p   &
       +((2.7759e-12_WP*t-1.1351e-10_WP)*ds+((-5.4481e-14_WP*t        &
       +8.733e-12_WP)*t-6.7795e-10_WP)*t+1.8741e-8_WP))*p             &
       +(-4.2393e-8_WP*t+1.8932e-6_WP)*ds                          &
       +((6.6228e-10_WP*t-6.836e-8_WP)*t+8.5258e-6_WP)*t+3.5803e-5_WP

  return
end function atg
!
!
!
!===============================================================================
subroutine sw_alpha_beta(TF1,SF1, mesh)
  ! DESCRIPTION:
  !   A function to calculate the thermal expansion coefficient
  !   and saline contraction coefficient. (elementwise)
  !
  ! INPUT:
  !   tracer(:,2) = salinity              [psu      (PSS-78)]
  !   tracer(:,1) = potential temperature [degree C (ITS-90)]
  !   z           = pressure (or -depth)  [db]
  !
  ! OUTPUT:
  !   sw_alpha = Thermal expansion coeff (alpha) [degree_C.^-1]
  !   sw_beta  = Saline contraction coeff (beta) [psu.^-1]
  !
  ! Qiang Wang, 25,11,2004
  !
  ! REFERENCE:
  !    McDougall, T.J. 1987.  Neutral Surfaces
  !    Journal of Physical Oceanography, vol 17, 1950-1964,
  !-----------------------------------------------------------------
  ! CHECK VALUE:
  !    sw_beta=0.72088e-3 psu^-1 @ S=40.0psu, ptmp=10.0C (ITS-90), p=4000db
  !    a_over_b=0.34765 psu*C^-1 @ S=40.0psu, ptmp=10.0C, p=4000db
  !-----------------------------------------------------------------
  use mod_mesh
  use o_arrays
  use g_parsup
  use o_param
  use g_comm_auto
  implicit none
  !
  type(t_mesh), intent(in) , target :: mesh
  integer        :: n, nz, nzmin, nzmax
  real(kind=WP)  :: t1,t1_2,t1_3,t1_4,p1,p1_2,p1_3,s1,s35,s35_2 
  real(kind=WP)  :: a_over_b    
  real(kind=WP)  :: TF1(mesh%nl-1, myDim_nod2D+eDim_nod2D),SF1(mesh%nl-1, myDim_nod2D+eDim_nod2D)

#include "associate_mesh.h"

  do n = 1,myDim_nod2d
     nzmin = ulevels_nod2d(n)   
     nzmax = nlevels_nod2d(n)   
     !!PS do nz=1, nlevels_nod2d(n) -1
     do nz=nzmin, nzmax-1
     
     t1 = TF1(nz,n)*1.00024_WP
     s1 = SF1(nz,n)
    !!PS      p1 = abs(Z(nz)) 
     p1 = abs(Z_3d_n(nz,n)) 
     
     t1_2 = t1*t1
     t1_3 = t1_2*t1
     t1_4 = t1_3*t1
     p1_2 = p1*p1
     p1_3 = p1_2*p1
     s35 = s1-35.0_WP
     s35_2 = s35*s35

     ! calculate beta
     sw_beta(nz,n) = 0.785567e-3_WP - 0.301985e-5_WP*t1 &
          + 0.555579e-7_WP*t1_2 - 0.415613e-9_WP*t1_3 &
          + s35*(-0.356603e-6_WP + 0.788212e-8_WP*t1 &
          + 0.408195e-10_WP*p1 - 0.602281e-15_WP*p1_2) &
          + s35_2*(0.515032e-8_WP) & 
          + p1*(-0.121555e-7_WP + 0.192867e-9_WP*t1 - 0.213127e-11_WP*t1_2) &
          + p1_2*(0.176621e-12_WP - 0.175379e-14_WP*t1) &
          + p1_3*(0.121551e-17_WP)

     ! calculate the thermal expansion / saline contraction ratio
     a_over_b = 0.665157e-1_WP + 0.170907e-1_WP*t1 &
          - 0.203814e-3_WP*t1_2 + 0.298357e-5_WP*t1_3 &
          - 0.255019e-7_WP*t1_4 &
          + s35*(0.378110e-2_WP - 0.846960e-4_WP*t1 &
          - 0.164759e-6_WP*p1 - 0.251520e-11_WP*p1_2) &
          + s35_2*(-0.678662e-5_WP) &
          + p1*(0.380374e-4_WP - 0.933746e-6_WP*t1 + 0.791325e-8_WP*t1_2) &
          + p1_2*t1_2*(0.512857e-12_WP) &
          - p1_3*(0.302285e-13_WP)

     ! calculate alpha
     sw_alpha(nz,n) = a_over_b*sw_beta(nz,n)
   end do
 end do
call exchange_nod(sw_alpha)
call exchange_nod(sw_beta)
end subroutine sw_alpha_beta
!
!
!
!===============================================================================
subroutine compute_sigma_xy(TF1,SF1, mesh)
  !--------------------------------------------------------------------
  ! DESCRIPTION:
  !   computes density gradient
  !
  ! INPUT:
  !   SF          = salinity              [psu      (PSS-78)]
  !   TF          = potential temperature [degree C (ITS-90)]
  !
  ! OUTPUT:
  ! based on thermal expansion and saline contraction coefficients
  ! computes density gradient sigma_xy
  !-------------------------------------------------------------------
  use mod_mesh
  use o_param
  use o_arrays
  use g_parsup
  use g_comm_auto
  implicit none
  !
  type(t_mesh),  intent(in)   , target :: mesh
  real(kind=WP), intent(IN)   :: TF1(mesh%nl-1, myDim_nod2D+eDim_nod2D), SF1(mesh%nl-1, myDim_nod2D+eDim_nod2D)
  real(kind=WP)               :: tx(mesh%nl-1), ty(mesh%nl-1), sx(mesh%nl-1), sy(mesh%nl-1), vol(mesh%nl-1), testino(2)
  integer                     :: n, nz, elnodes(3),el, k, nln, uln, nle, ule

#include "associate_mesh.h"
  !
  DO n=1, myDim_nod2D
        nln = nlevels_nod2D(n)-1
        uln = ulevels_nod2D(n)
        !!PS vol(1:nl1) = 0.0_WP
        !!PS tx(1:nl1)  = 0.0_WP
        !!PS ty(1:nl1)  = 0.0_WP
        !!PS sx(1:nl1)  = 0.0_WP
        !!PS sy(1:nl1)  = 0.0_WP
        vol(uln:nln) = 0.0_WP
        tx(uln:nln)  = 0.0_WP
        ty(uln:nln)  = 0.0_WP
        sx(uln:nln)  = 0.0_WP
        sy(uln:nln)  = 0.0_WP
        DO k=1, nod_in_elem2D_num(n)
           el=nod_in_elem2D(k, n)
           nle = nlevels(el)-1 
           ule = ulevels(el)
           !!PS DO nz=1, nlevels(el)-1 
           DO nz=ule, nle           
              vol(nz) = vol(nz)+elem_area(el)

              !NR  writing the sum over elem2D_nodes explicitly helps the compiler to vectorize the nz-loop

              tx(nz) = tx(nz)+(gradient_sca(1,el)*TF1(nz,elem2D_nodes(1,el)) &
                             + gradient_sca(2,el)*TF1(nz,elem2D_nodes(2,el)) &
                             + gradient_sca(3,el)*TF1(nz,elem2D_nodes(3,el)))*elem_area(el)

              ty(nz) = ty(nz)+(gradient_sca(4,el)*TF1(nz,elem2D_nodes(1,el)) &
                             + gradient_sca(5,el)*TF1(nz,elem2D_nodes(2,el)) &
                             + gradient_sca(6,el)*TF1(nz,elem2D_nodes(3,el)))*elem_area(el)

              sx(nz) = sx(nz)+(gradient_sca(1,el)*SF1(nz,elem2D_nodes(1,el)) &
                             + gradient_sca(2,el)*SF1(nz,elem2D_nodes(2,el)) &
                             + gradient_sca(3,el)*SF1(nz,elem2D_nodes(3,el)))*elem_area(el)

              sy(nz) = sy(nz)+(gradient_sca(4,el)*SF1(nz,elem2D_nodes(1,el)) &
                             + gradient_sca(5,el)*SF1(nz,elem2D_nodes(2,el)) &
                             + gradient_sca(6,el)*SF1(nz,elem2D_nodes(3,el)))*elem_area(el)
           END DO
        enddo
        !!PS sigma_xy(1,1:nl1,n) = (-sw_alpha(1:nl1,n)*tx(1:nl1)+sw_beta(1:nl1,n)*sx(1:nl1))/vol(1:nl1)*density_0
        !!PS sigma_xy(2,1:nl1,n) = (-sw_alpha(1:nl1,n)*ty(1:nl1)+sw_beta(1:nl1,n)*sy(1:nl1))/vol(1:nl1)*density_0
        sigma_xy(1,uln:nln,n) = (-sw_alpha(uln:nln,n)*tx(uln:nln)+sw_beta(uln:nln,n)*sx(uln:nln))/vol(uln:nln)*density_0
        sigma_xy(2,uln:nln,n) = (-sw_alpha(uln:nln,n)*ty(uln:nln)+sw_beta(uln:nln,n)*sy(uln:nln))/vol(uln:nln)*density_0
  END DO 

  call exchange_nod(sigma_xy)
end subroutine compute_sigma_xy
!
!
!
!===============================================================================
subroutine compute_neutral_slope(mesh)
    use o_ARRAYS
    use g_PARSUP
    use MOD_MESH
    use o_param
    use g_config
    use g_comm_auto
    IMPLICIT NONE
    real(kind=WP)   :: deltaX1,deltaY1,deltaX2,deltaY2
    integer         :: edge
    integer         :: n,nz,nl1,ul1,el(2),elnodes(3),enodes(2)
    real(kind=WP)   :: c, ro_z_inv,eps,S_cr,S_d
    type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"
    !if sigma_xy is not computed
    eps=5.0e-6_WP
    S_cr=1.0e-2_WP
    S_d=1.0e-3_WP
    slope_tapered=0._WP
    do n=1, myDim_nod2D
        nl1=nlevels_nod2d(n)-1
        ul1=ulevels_nod2d(n)
        !!PS do nz = 2,nl1
        do nz = ul1+1,nl1
            ro_z_inv=2._WP*g/density_0/max(bvfreq(nz,n)+bvfreq(nz+1,n), eps**2) !without minus, because neutral slope S=-(nabla\rho)/(d\rho/dz)
            neutral_slope(1,nz,n)=sigma_xy(1,nz,n)*ro_z_inv
            neutral_slope(2,nz,n)=sigma_xy(2,nz,n)*ro_z_inv
            neutral_slope(3,nz,n)=sqrt(neutral_slope(1,nz,n)**2+neutral_slope(2,nz,n)**2)
            !tapering
            c=1.0_WP
            c=0.5_WP*(1.0_WP + tanh((S_cr - neutral_slope(3,nz,n))/S_d))
            if ((bvfreq(nz,n) <= 0.0_WP) .or. (bvfreq(nz+1,n) <= 0.0_WP)) c=0.0_WP
            slope_tapered(:,nz,n)=neutral_slope(:,nz,n)*c
!                       slope_tapered(:,nl1-1:nl1,n)=0.
!                       slope_tapered(:,1:2,n)      =0.
        enddo
    enddo

    call exchange_nod(neutral_slope)
    call exchange_nod(slope_tapered)
end subroutine compute_neutral_slope
!
!
!
!===============================================================================
!converts insitu temperature to a potential one
!               tr_arr(:,:,1) will be modified!
subroutine insitu2pot(mesh)
  use mod_mesh
  use o_param
  use o_arrays
  use g_config
  use g_PARSUP
  implicit none
  real(kind=WP), external     :: ptheta
  real(kind=WP)               :: pp, pr, tt, ss
  integer                     :: n, nz, nzmin,nzmax
  type(t_mesh), intent(in) , target :: mesh

#include  "associate_mesh.h"
 
  ! Convert in situ temperature into potential temperature
  pr=0.0_WP
  do n=1,myDim_nod2d+eDim_nod2D
     nzmin = ulevels_nod2D(n)
     nzmax = nlevels_nod2D(n)
     !!PS do nz=1, nlevels_nod2D(n)-1
     do nz=nzmin, nzmax-1    
        tt=tr_arr(nz,n,1)
        ss=tr_arr(nz,n,2)
        
        !!PS ___________________________________________________________________
        !!PS using here Z_3d_n at the beginning makes the model very instable after 
        !!PS the initialisation when using partial cell, requires smaller time step
        !!PS when using partial cell --> i would stick here with Z --> since we 
        !!PS anyway do a spinup and it its only used at initialisation time
        !!PS pp=abs(Z_3d_n(nz,n))
        pp=abs(Z(nz))
        tr_arr(nz,n,1)=ptheta(ss, tt, pp, pr)
     end do
  end do
end subroutine insitu2pot
!
!
!
!===============================================================================
SUBROUTINE density_linear(t, s, bulk_0, bulk_pz, bulk_pz2, rho_out, mesh)
!coded by Margarita Smolentseva, 21.05.2020
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP !, only: par_ex,pe_status
use g_config !, only: which_toy, toy_ocean
IMPLICIT NONE

  real(kind=WP), intent(IN)  :: t,s
  real(kind=WP), intent(OUT) :: rho_out                 
  real(kind=WP)              :: rhopot, bulk
  real(kind=WP), intent(OUT) :: bulk_0, bulk_pz, bulk_pz2
  type(t_mesh), intent(in)   , target :: mesh
#include "associate_mesh.h"
  !compute secant bulk modulus
  
  bulk_0   = 1
  bulk_pz  = 0
  bulk_pz2 = 0

  IF((toy_ocean) .AND. (TRIM(which_toy)=="soufflet")) THEN
      rho_out  = density_0 - 0.00025_WP*(t - 10.0_WP)*density_0
  ELSE
      rho_out  = density_0 + 0.8_WP*(s - 34.0_WP) - 0.2*(t - 20.0_WP)
  END IF
  
end subroutine density_linear
!
!
!
!===============================================================================
subroutine init_ref_density(mesh)
    ! compute reference density
    ! Coded by Qiang Wang
    ! Reviewed by ??
    !___________________________________________________________________________
    USE MOD_MESH
    use o_PARAM
    use o_ARRAYS
    use g_PARSUP
    use densityJM_components_interface
    implicit none
    
    !___________________________________________________________________________
    type(t_mesh), intent(in) , target :: mesh    
    integer         :: node, nz, nzmin, nzmax
    real(kind=WP)   :: rhopot, bulk_0, bulk_pz, bulk_pz2, rho 
    real(kind=8)    :: T, S, auxz

#include "associate_mesh.h"

    !___________________________________________________________________________
!!PS     S=34.  
!!PS     T=2.0 
    
    !___________________________________________________________________________
    density_ref = 0.0_WP
    do node=1,myDim_nod2d+eDim_nod2d
        !!PS nzmin = ulevels_nod2d(node)
        nzmin = 1
        nzmax = nlevels_nod2d(node)-1
        auxz=min(0.0,Z_3d_n(nzmin,node))
        
        !_______________________________________________________________________
        call densityJM_components(density_ref_T, density_ref_S, bulk_0, bulk_pz, bulk_pz2, rhopot, mesh)
        rho = bulk_0   + auxz*bulk_pz   + auxz*bulk_pz2
        density_ref(nzmin,node) = rho*rhopot/(rho+0.1_WP*auxz)
        
        !_______________________________________________________________________
        do nz=nzmin+1,nzmax
            auxz=Z_3d_n(nz,node)
            rho = bulk_0   + auxz*bulk_pz   + auxz*bulk_pz2
            density_ref(nz,node) = rho*rhopot/(rho+0.1_WP*auxz)
        end do
    end do
    if(mype==0) write(*,*) ' --> compute reference density'
end subroutine init_ref_density

