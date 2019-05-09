!
!
!===============================================================================
! module interface to FESOM2.0 for the CVMIX library for the calculation 
! of vertical mixing: PP scheme 
!
! @see  Pacanowski R.C. and Philander S.G.H. --> PP81
!       Parameterizations of Vertical Mixing in Numerical Models of Tropical Oceans.
!       J. Phys. Oceanogr., XX, XX–XX, doi:, 1981.
!
!
! @see  Large W.G., McWilliams J.C., Doney S.C. -> KPD94
!       Oceanic Vertical Mixing: A Review and a Model with a Nonlocal
!       Boundary Layer Parameterizations.
!       Rev. of Geophys., XX,XX–XX. doi: , 1994.
!
module g_cvmix_pp
    !___________________________________________________________________________
    ! module calls from cvmix library
    use cvmix_shear,    only: cvmix_init_shear, cvmix_coeffs_shear  ! calculates avo,dvo (PP81)
    
    !___________________________________________________________________________
    ! module calls from FESOM
    use g_config
    use o_param           
    use o_mesh
    use g_parsup
    use o_arrays
    use g_comm_auto 
    use i_arrays
    implicit none
   
    ! namelist parameters associated with the PP scheme
    ! replace later with updates from namelist values
    real(kind=WP)      :: pp_Av0             = 0.01
    real(kind=WP)      :: pp_alpha           = 5.0
    real(kind=WP)      :: pp_exp             = 2.0 
    real(kind=WP)      :: pp_Avbckg          = 1.0e-4
    real(kind=WP)      :: pp_Kvbckg          = 1.0e-5
    logical            :: pp_use_monob       = .true.
    real(kind=WP)      :: pp_monob_Kv        = 0.01
    logical            :: pp_use_nonconstKvb = .true.
    logical            :: pp_use_windmix     = .true.
    logical            :: pp_use_instabmix   = .true.
    
    namelist /param_pp/ pp_Av0, pp_alpha, pp_exp, pp_Avbckg, pp_Kvbckg, pp_use_monob, &
                        pp_monob_Kv, pp_use_nonconstKvb, pp_use_windmix, pp_use_instabmix
    
    real(kind=WP), allocatable, dimension(:,:) :: pp_Av, pp_Kv
    real(kind=WP), allocatable, dimension(:,:) :: pp_richardnmb ! store Ridchardsen number
    real(kind=WP), allocatable, dimension(:)   :: pp_monob_mixl  ! store Monin-Obukov mixing length
    
    contains
    !
    !
    !
    !===========================================================================
    ! allocate and initialize CVMIX PP variables --> call initialisation 
    ! routine from cvmix library
    subroutine init_cvmix_pp
        character(len=100) :: nmlfile
        logical            :: nmlfile_exist=.False.
        integer            :: node_size
        
        !_______________________________________________________________________
        ! allocate + initialse all tke arrays
        node_size=myDim_nod2D+eDim_nod2D
        
        ! initialize input fields 
        allocate(pp_Av(nl,node_size),pp_Kv(nl,node_size))
        pp_Av         = 0.0_WP
        pp_Kv         = 0.0_WP
        allocate(pp_richardnmb(nl-1,node_size))
        pp_richardnmb = 0.0_WP
        allocate(pp_monob_mixl(node_size))
        pp_monob_mixl = 0.0_WP
        
        !_______________________________________________________________________
        ! read cvmix namelist file 
        nmlfile ='namelist.cvmix'    ! name of ocean namelist file
        ! check if cvmix namelist file exists if not use default values 
        inquire(file=trim(nmlfile),exist=nmlfile_exist) 
        if (nmlfile_exist) then
            open(20,file=trim(nmlfile))
                read(20,nml=param_pp)
            close(20)
        end if
        
        ! if non-constant background diffusivity is supposed to be used --> set 
        ! normaly used cvmix background diffusivity to zero independet of what 
        ! is written in namelist.cvmix
        if (pp_use_nonconstKvb==.False.) then
            pp_Kvbckg = 0.0_WP
        endif 
        
        !_______________________________________________________________________
        ! write info to log file 
        if (mype==0) then
            write(*,*) "     pp_Av0             = ", pp_Av0
            write(*,*) "     pp_alpha           = ", pp_alpha
            write(*,*) "     pp_exp             = ", pp_exp
            write(*,*) "     pp_Avbck           = ", pp_Avbckg
            write(*,*) "     pp_Kvbck           = ", pp_Kvbckg
            write(*,*) "     pp_use_monob       = ", pp_use_monob
            write(*,*) "     pp_monob_Kv        = ", pp_monob_Kv
            write(*,*) "     pp_use_nonconstKvb = ", pp_use_nonconstKvb
            write(*,*) "     pp_use_windmix     = ", pp_use_windmix
            write(*,*) "     pp_use_instabmix   = ", pp_use_instabmix
            write(*,*)
        end if
        
        !_______________________________________________________________________
        ! Pacanowski and Philander 1981
        ! eq.1 ... visc  = PP_nu_zero/(1+PP_alpha*Ri)^PP_exp + PP_nu_b
        ! eq.2 ... kappa = visc/(1+PP_alpha*Ri) + PP_kappa_b  
        call cvmix_init_shear(mix_scheme  = 'PP',           &
                              PP_nu_zero  = pp_Av0,         &
                              PP_alpha    = pp_alpha,       &
                              PP_exp      = pp_exp,         &
                              PP_nu_b     = pp_Avbckg,      &
                              PP_kappa_b  = pp_Kvbckg)      
        
    end subroutine init_cvmix_pp
    !
    !
    !
    !===========================================================================
    ! calculate PP vertrical mixing coefficients from CVMIX library
    subroutine calc_cvmix_pp
        
        integer       :: node, elem, nz, nln, elnodes(3)
        real(kind=WP) :: vshear2, dz2, Kvb
        real(kind=WP) :: wndmix=1.e-3, wndnl=2, kv_conv=0.1_WP, av_conv=0.1_WP
        
        !_______________________________________________________________________
        pp_richardnmb = 0.0_WP
        
        do node = 1,myDim_nod2D
            !___________________________________________________________________
            ! number of above bottom levels at node
            nln = nlevels_nod2D(node)-1
            
            !___________________________________________________________________
            ! calculate Richardsen number
            do nz=2,nln
                dz2     = (Z_3d_n(nz-1,node)-Z_3d_n(nz,node))**2
                vshear2 = (Unode(1,nz-1,node)-Unode(1,nz,node))**2 +&
                          (Unode(2,nz-1,node)-Unode(2,nz,node))**2 
                vshear2 = vshear2/dz2
                pp_richardnmb(nz,node) = bvfreq(nz,node)/vshear2
            end do
            
            !___________________________________________________________________
            ! use cvmix library function 
            call cvmix_coeffs_shear(Mdiff_out = pp_Av(:,node),         &
                                    Tdiff_out = pp_Kv(:,node),         &
                                    RICH      = pp_richardnmb(:,node), &
                                    nlev      = nln,                   &
                                    max_nlev  = nl-1)
            !___________________________________________________________________
            pp_Av(1,node)=0.0_WP
            pp_Kv(1,node)=0.0_WP
            
        end do !--> do node = 1,myDim_nod2D
        
        !_______________________________________________________________________
        ! add vertical mixing scheme of Timmermann and Beckmann, 2004,"Parameterization 
        ! of vertical mixing in the Weddell Sea!
        ! Computes the mixing length derived from the Monin-Obukhov length
        ! --> in FESOM1.4 refered as TB04 mixing scheme
        if (pp_use_monob) then
            do node = 1,myDim_nod2D
                !_______________________________________________________________
                ! calcualte monin obukov length
                call mo_length(water_flux(node),heat_flux(node), &         
                               stress_atmoce_x(node),stress_atmoce_y(node), &    
                               u_ice(node),v_ice(node),a_ice(node), &                             
                               dt, pp_monob_mixl(node))
                !_______________________________________________________________
                ! increase vertical diffusion within monin obukov length to namelist
                ! parameter pp_monob_Kv. pp_monob_Kv in moment set to 0.01 --> that means 
                ! very strong vertical mixing within mixlength
                do nz = 2,nlevels_nod2D(node)-1
                    if(abs(zbar_3d_n(nz,node)) <= pp_monob_mixl(node)) then
                        pp_Kv(nz,node) = pp_Kv(nz,node) + pp_monob_Kv
                    else 
                        exit    
                    end if 
                end do 
            end do    
        end if 
           
        !_______________________________________________________________________
        ! calculate and add latitudinal and depth dependend background 
        ! diffusivity of Q. Wang from FESOM1.4
        if (pp_use_nonconstKvb) then   
            do node = 1,myDim_nod2D
                do nz = 2,nlevels_nod2D(node)-1
                    call Kv0_background_qiang(Kvb,geo_coord_nod2D(2,node)/rad,abs(zbar_3d_n(nz,node)))
                    pp_Kv(nz,node) = pp_Kv(nz,node) + Kvb
                end do
            end do
        end if 
        
        !_______________________________________________________________________
        ! enhance mixing in case of instable stratification  --> (N^2>0)
        if (pp_use_instabmix) then
            do node = 1,myDim_nod2D
                do nz = 2,nlevels_nod2D(node)-1
                    if (bvfreq(nz,node) < 0.0_WP) then 
                        pp_Kv(nz,node) = kv_conv
                        pp_Av(nz,node) = av_conv
                    end if 
                end do
            end do
        end if
        
        !_______________________________________________________________________
        ! add additional wind mixing for upper layers 
        if (pp_use_windmix) then
            do node = 1,myDim_nod2D
                do nz = 2,nlevels_nod2D(node)-1
                    if (nz <= wndnl+1) then
                        pp_Kv(nz,node)=max(pp_Kv(nz,node), wndmix)
                        pp_Av(nz,node)=max(pp_Av(nz,node), wndmix)
                    else
                        exit
                    end if 
                end do
            end do
        end if
        
        !_______________________________________________________________________
        ! write out diffusivity
        call exchange_nod(pp_Kv)
        Kv = pp_Kv
           
        !_______________________________________________________________________
        ! write out viscosity -->interpolate therefor from nodes to elements
        call exchange_nod(pp_Av) !Warning: don't forget to communicate before averaging on elements!!!
        Av = 0.0_WP
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            do nz=2,nlevels(elem)-1
                Av(nz,elem) = sum(pp_Av(nz,elnodes))/3.0_WP    ! (elementwise)                
            end do
        end do
        call exchange_elem(Av)
        
    end subroutine calc_cvmix_pp
end module g_cvmix_pp
