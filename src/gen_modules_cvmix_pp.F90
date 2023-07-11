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
! written by Patrick Scholz, 10.05.2019
module g_cvmix_pp
    !___________________________________________________________________________
    ! module calls from cvmix library
    use cvmix_shear,    only: cvmix_init_shear, cvmix_coeffs_shear
    
    !___________________________________________________________________________
    ! module calls from FESOM
    use g_config
    use o_param           
    use MOD_MESH
    use g_parsup
    use o_arrays
    use g_comm_auto 
    use i_arrays
    implicit none
   
    !___________________________________________________________________________
    ! setting default namelist parameters associated with the PP scheme
    ! replace later with updates from namelist values
    real(kind=WP)      :: pp_Av0             = 0.01
    real(kind=WP)      :: pp_alpha           = 5.0
    real(kind=WP)      :: pp_exp             = 2.0 
    real(kind=WP)      :: pp_Avbckg          = 1.0e-4
    real(kind=WP)      :: pp_Kvbckg          = 1.0e-5
    
    ! pp_use_fesompp ... general flag to use fesom flavour of PP mixing if false 
    ! original pp mixing of Pacanowski and Philnader 1981 is used
    logical            :: pp_use_fesompp     = .true. 
    
    ! single elements of fesom flavour PP mixing can be switched of
    logical            :: pp_use_AvbinKv     = .true. ! use term Kv=...+AvB/(1+5Ri)+...
    logical            :: pp_use_nonconstKvb = .true. ! use nonkonst Kvbckg of qiang 
    
    
    namelist /param_pp/ pp_Av0, pp_alpha, pp_exp, pp_Avbckg, pp_Kvbckg, &
                        pp_use_nonconstKvb, pp_use_fesompp, pp_use_AvbinKv
    
    !___________________________________________________________________________
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
    subroutine init_cvmix_pp(mesh)
        use MOD_MESH
        implicit none
        type(t_mesh), intent(in), target :: mesh        
        character(len=MAX_PATH)       :: nmlfile
        logical                  :: nmlfile_exist=.False.
        integer                  :: node_size
#include "associate_mesh.h"
        !_______________________________________________________________________
        if(mype==0) then
            write(*,*) '____________________________________________________________'
            write(*,*) ' --> initialise CVMIX_PP'
            write(*,*)
        end if
        
        !_______________________________________________________________________
        ! allocate + initialse all pp arrays
        node_size=myDim_nod2D+eDim_nod2D
        allocate(pp_Av(nl,node_size),pp_Kv(nl,node_size))
        pp_Av         = 0.0_WP
        pp_Kv         = 0.0_WP
        
        allocate(pp_richardnmb(nl,node_size))
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
        else
            write(*,*) '     could not find namelist.cvmix, will use default values !'
        end if
        
        ! if non-constant background diffusivity is supposed to be used --> set 
        ! normaly used cvmix background diffusivity to zero independet of what 
        ! is written in namelist.cvmix
        if (pp_use_fesompp .and. pp_use_nonconstKvb) then
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
            write(*,*) "     pp_use_fesompp     = ", pp_use_fesompp
            write(*,*) "     pp_use_nonconstKvb = ", pp_use_nonconstKvb
            write(*,*)
        end if
        
        !_______________________________________________________________________
        ! Initialise CVMIX Pacanowski and Philander 1981 vertical mixing 
        ! parameterisation
        ! eq.1 ... visc  = PP_nu_zero/(1+PP_alpha*Ri)^PP_exp + PP_nu_b
        ! eq.2 ... kappa = visc/(1+PP_alpha*Ri) + PP_kappa_b  
        !                = PP_nu_zero/(1+PP_alpha*Ri)^(PP_exp+1) + 
        !                  PP_nu_b/(1+PP_alpha*Ri) +
        !                  PP_kappa_b
        if (pp_use_fesompp .and. pp_use_AvbinKv .eqv. .false.) then
            ! ommit the term PP_nu_b/(1+PP_alpha*Ri) in kappa this can make an
            ! already diffusive model even more diffusive --> it was done in 
            ! FESOM1.4 like this. In this case set pp_Avbckg and pp_Kvbckg by
            ! hand
            call cvmix_init_shear(mix_scheme  = 'PP',         &
                                PP_nu_zero  = pp_Av0,         &
                                PP_alpha    = pp_alpha,       &
                                PP_exp      = pp_exp,         &
                                PP_nu_b     = 0.0_WP,         &
                                PP_kappa_b  = 0.0_WP)      
        else
            call cvmix_init_shear(mix_scheme  = 'PP',         &
                                PP_nu_zero  = pp_Av0,         &
                                PP_alpha    = pp_alpha,       &
                                PP_exp      = pp_exp,         &
                                PP_nu_b     = pp_Avbckg,      &
                                PP_kappa_b  = pp_Kvbckg)      
        end if
    end subroutine init_cvmix_pp
    !
    !
    !
    !===========================================================================
    ! calculate PP vertrical mixing coefficients from CVMIX library
    subroutine calc_cvmix_pp(mesh)
        use MOD_MESH
        implicit none
        type(t_mesh), intent(in), target :: mesh                
        integer       :: node, elem, nz, nln, nun, elnodes(3), windnl=2, node_size
        real(kind=WP) :: vshear2, dz2, Kvb
#include "associate_mesh.h"
        node_size = myDim_nod2D
        !_______________________________________________________________________
        do node = 1,node_size
            !___________________________________________________________________
            ! number of above bottom levels at node
            nln = nlevels_nod2D(node)-1
            nun = ulevels_nod2D(node)
            
            !___________________________________________________________________
            ! calculate Richardson number
            !!PS do nz=2,nln
            do nz=nun+1,nln
                dz2     = (Z_3d_n( nz-1,node)-Z_3d_n( nz,node))**2
                vshear2 = (Unode(1,nz-1,node)-Unode(1,nz,node))**2 +&
                          (Unode(2,nz-1,node)-Unode(2,nz,node))**2 
                vshear2 = vshear2/dz2
                ! WIKIPEDIA: The Richardson number is always 
                ! considered positive. A negative value of N² (i.e. complex N) 
                ! indicates unstable density gradients with active convective 
                ! overturning. Under such circumstances the magnitude of negative 
                ! Ri is not generally of interest. It can be shown that Ri < 1/4 
                ! is a necessary condition for velocity shear to overcome the 
                ! tendency of a stratified fluid to remain stratified, and 
                ! some mixing (turbulence) will generally occur. When Ri is 
                ! large, turbulent mixing across the stratification is 
                ! generally suppressed
                !!PS pp_richardnmb(nz,node) = bvfreq(nz,node)/vshear2
                pp_richardnmb(nz,node) = max(bvfreq(nz,node),0.0_WP)/vshear2
                !                                            ______
                !          takes care that Ri stays positiv <--|
            end do
            
            
            !___________________________________________________________________
            ! use cvmix library function 
            call cvmix_coeffs_shear(Mdiff_out = pp_Av(:,node),         &
                                    Tdiff_out = pp_Kv(:,node),         &
                                    RICH      = pp_richardnmb(:,node), &
                                    nlev      = nln,                   &
                                    max_nlev  = nl-1)
            
            !___________________________________________________________________
            ! In the fesom flavour of PP mixing is the term from the background
            ! viscosity in the diffusivity Avb/(1+5Ri) ommited --> therefor set 
            ! cvmix PP_nu_b and PP_kappa_b to zero and add them from hand
            ! Av = Av0/(1+5Ri)^2 + Avb
            ! Kv = Av/(1+5Ri) + Kvb
            !    = Av0/(1+5Ri)^3 + Avb/(1+5Ri) + Kvb
            !                      -----------
            ! This can make an already diffusive model even more diffusive --> it 
            ! needs to be tested if this is an advantage for FESOM2.0 or not 
            if (pp_use_fesompp .and. pp_use_AvbinKv .eqv. .false.) then
                !!PS pp_Av(2:nln,node) = pp_Av(2:nln,node) + pp_Avbckg
                pp_Av(nun+1:nln,node) = pp_Av(nun+1:nln,node) + pp_Avbckg
            end if  
            
            !_______________________________________________________________________
            ! calculate and add latitudinal and depth dependend background 
            ! diffusivity of Q. Wang from FESOM1.4
            if (pp_use_fesompp .and. pp_use_nonconstKvb) then   
                !!PS do nz = 2,nlevels_nod2D(node)-1
                do nz = nun+1,nlevels_nod2D(node)-1
                    call Kv0_background_qiang(Kvb,                       &
                                                geo_coord_nod2D(2,node)/rad, &
                                                abs(zbar_3d_n(nz,node))      &
                                                )
                    pp_Kv(nz,node) = pp_Kv(nz,node) + Kvb
                end do
            else
                !!PS pp_Kv(2:nln,node) = pp_Kv(2:nln,node) + pp_Kvbckg
                pp_Kv(nun+1:nln,node) = pp_Kv(nun+1:nln,node) + pp_Kvbckg
            end if 
            
            !___________________________________________________________________
            !!PS pp_Av(1,node)=0.0_WP
            !!PS pp_Kv(1,node)=0.0_WP
            pp_Av(nun,node)=0.0_WP
            pp_Kv(nun,node)=0.0_WP
            
        end do !--> do node = 1,myDim_nod2D

        !_______________________________________________________________________
        ! write out diffusivities to FESOM2.0 --> diffusivities remain on nodes
        call exchange_nod(pp_Kv)
        Kv = pp_Kv
           
        !_______________________________________________________________________
        ! write out viscosities to FESOM2.0 --> viscosities for FESOM2.0 are 
        ! defined on elements --> interpolate therefor from nodes to elements
        call exchange_nod(pp_Av)
        Av = 0.0_WP
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            !!PS do nz=2,nlevels(elem)-1
            do nz=ulevels(elem)+1,nlevels(elem)-1
                Av(nz,elem) = sum(pp_Av(nz,elnodes))/3.0_WP    ! (elementwise)                
            end do
        end do
    end subroutine calc_cvmix_pp
end module g_cvmix_pp
