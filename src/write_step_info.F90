!===============================================================================
module write_step_info_module
  use g_config, only: dt, use_ice, use_icebergs, ib_num, logfile_outfreq, which_ALE, toy_ocean
  use MOD_MESH
  use MOD_PARTIT
  use MOD_PARSUP
  use MOD_TRACER
  use MOD_DYN
  use MOD_ICE
  use o_PARAM
  use o_ARRAYS, only: water_flux, heat_flux, pgf_x, pgf_y, Av, Kv, stress_surf
  use g_comm_auto
  use g_support
  use iceberg_params
  use io_BLOWUP
  use g_forcing_arrays
  use diagnostics
  implicit none
  private
  public :: write_step_info, check_blowup, write_enegry_info
contains
!===============================================================================
subroutine write_step_info(istep, outfreq, ice, dynamics, tracers, partit, mesh)
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  use MOD_TRACER
  use MOD_DYN
  use MOD_ICE
  integer                               :: istep,outfreq
  type(t_mesh),   intent(in)   , target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_tracer), intent(in)   , target :: tracers
  type(t_dyn)   , intent(in)   , target :: dynamics
  type(t_ice)   , intent(in)   , target :: ice
  real(kind=WP), dimension(:,:,:), pointer :: UV, UVnode
  real(kind=WP), dimension(:,:)  , pointer :: Wvel, CFL_z
  real(kind=WP), dimension(:)    , pointer :: eta_n, d_eta, m_ice
  
  ! Variable declarations
  real(kind=WP) :: int_eta, int_hbar, int_deta, int_dhbar, int_wflux, int_hflux, int_temp, int_salt
  real(kind=WP) :: loc_eta, loc_hbar, loc_deta, loc_dhbar, loc_wflux, loc
  real(kind=WP) :: min_eta, min_hbar, min_deta, min_wflux, min_hflux, min_temp, min_salt
  real(kind=WP) :: min_wvel, min_wvel2, min_uvel, min_uvel2, min_vvel, min_vvel2, min_hnode, min_hnode2
  real(kind=WP) :: max_eta, max_hbar, max_deta, max_wflux, max_hflux, max_temp, max_salt
  real(kind=WP) :: max_wvel, max_wvel2, max_uvel, max_uvel2, max_vvel, max_vvel2, max_hnode, max_hnode2
  real(kind=WP) :: max_cfl_z, max_pgfx, max_pgfy, max_m_ice, max_av, max_kv
  integer :: n
  
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
  UV     => dynamics%uv(:,:,:)
  UVnode => dynamics%uvnode(:,:,:)
  Wvel   => dynamics%w(:,:)
  CFL_z  => dynamics%cfl_z(:,:)
  eta_n  => dynamics%eta_n(:)
  if ( .not. dynamics%use_ssh_se_subcycl) d_eta => dynamics%d_eta(:)
  m_ice  => ice%data(2)%values(:)
  
  if (mod(istep,outfreq)==0) then
  
  !_______________________________________________________________________
  int_eta   =0.
  int_hbar  =0.
  int_deta  =0.
  int_dhbar =0.
  int_wflux =0.
  int_hflux =0.
  int_temp  =0.
  int_salt  =0.
  loc_eta   =0.
  loc_hbar  =0.
  loc_deta  =0.
  loc_dhbar =0.
  loc_wflux =0.
  loc       =0.
  
  !_______________________________________________________________________
#if !defined(__openmp_reproducible)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n) REDUCTION(+:loc_eta, loc_hbar, loc_deta, loc_dhbar, loc_wflux)
#endif
  do n=1, myDim_nod2D
     loc_eta   = loc_eta   + areasvol(ulevels_nod2D(n), n)*eta_n(n)
     loc_hbar  = loc_hbar  + areasvol(ulevels_nod2D(n), n)*hbar(n)
     loc_dhbar = loc_dhbar + areasvol(ulevels_nod2D(n), n)*(hbar(n)-hbar_old(n))
     if ( .not. dynamics%use_ssh_se_subcycl) then 
          loc_deta  = loc_deta  + areasvol(ulevels_nod2D(n), n)*d_eta(n)
     end if 
     loc_wflux = loc_wflux + areasvol(ulevels_nod2D(n), n)*water_flux(n)
  end do
#if !defined(__openmp_reproducible)
!$OMP END PARALLEL DO     
#endif
  if (dynamics%use_ssh_se_subcycl) then
      loc_deta=loc_dhbar
  end if
  
  !_______________________________________________________________________
  call MPI_AllREDUCE(loc_eta  , int_eta  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  call MPI_AllREDUCE(loc_hbar , int_hbar , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
!PS     call MPI_AllREDUCE(loc_deta , int_deta , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  call MPI_AllREDUCE(loc_dhbar, int_dhbar, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  call MPI_AllREDUCE(loc_wflux, int_wflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)     

  int_eta  = int_eta  /ocean_areawithcav
  int_hbar = int_hbar /ocean_areawithcav
!PS     int_deta = int_deta /ocean_areawithcav
  int_dhbar= int_dhbar/ocean_areawithcav
  int_wflux= int_wflux/ocean_areawithcav      
  !_______________________________________________________________________
  loc=omp_min_max_sum1(eta_n,       1, myDim_nod2D, 'min', partit) 
  call MPI_AllREDUCE(loc , min_eta  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(hbar,        1, myDim_nod2D, 'min', partit) 
  call MPI_AllREDUCE(loc , min_hbar , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(water_flux,  1, myDim_nod2D, 'min', partit) 
  call MPI_AllREDUCE(loc , min_wflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(heat_flux,   1, myDim_nod2D, 'min', partit) 
  call MPI_AllREDUCE(loc , min_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum2(tracers%data(1)%values, 1, nl-1, 1, myDim_nod2D, 'min', partit, 0.0_WP) 
  call MPI_AllREDUCE(loc , min_temp , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum2(tracers%data(2)%values, 1, nl-1, 1, myDim_nod2D, 'min', partit, 0.0_WP) 
  call MPI_AllREDUCE(loc , min_salt , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(Wvel(1,:), 1, myDim_nod2D, 'min', partit)
  call MPI_AllREDUCE(loc , min_wvel , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(Wvel(2,:), 1, myDim_nod2D, 'min', partit)
  call MPI_AllREDUCE(loc , min_wvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(UVnode(1,1,:), 1, myDim_nod2D, 'min', partit)
  call MPI_AllREDUCE(loc , min_uvel , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(UVnode(1,2,:), 1, myDim_nod2D, 'min', partit)
  call MPI_AllREDUCE(loc , min_uvel2, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(UVnode(2,1,:), 1, myDim_nod2D, 'min', partit)
  call MPI_AllREDUCE(loc , min_vvel , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(UVnode(2,2,:), 1, myDim_nod2D, 'min', partit)
  call MPI_AllREDUCE(loc , min_vvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  if ( .not. dynamics%use_ssh_se_subcycl) then
      loc=omp_min_max_sum1(d_eta, 1, myDim_nod2D, 'min', partit)
  else
      loc=omp_min_max_sum1(hbar-hbar_old, 1, myDim_nod2D, 'min', partit)
  end if 
  call MPI_AllREDUCE(loc , min_deta  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(hnode(1,:), 1, myDim_nod2D, 'min', partit)
  call MPI_AllREDUCE(loc , min_hnode , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(hnode(2,:), 1, myDim_nod2D, 'min', partit)
  call MPI_AllREDUCE(loc , min_hnode2 , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  
  !_______________________________________________________________________
  loc=omp_min_max_sum1(eta_n, 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_eta  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(hbar, 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_hbar , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(water_flux, 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_wflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(heat_flux, 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum2(tracers%data(1)%values, 1, nl-1, 1, myDim_nod2D, 'max', partit, 0.0_WP) 
  call MPI_AllREDUCE(loc , max_temp , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum2(tracers%data(2)%values, 1, nl-1, 1, myDim_nod2D, 'max', partit, 0.0_WP)
  call MPI_AllREDUCE(loc , max_salt , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(Wvel(1,:), 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_wvel , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(Wvel(2,:), 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_wvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(UVnode(1,1,:), 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_uvel , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(UVnode(1,2,:), 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_uvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(UVnode(2,1,:), 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_vvel , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(UVnode(2,2,:), 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_vvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  if ( .not. dynamics%use_ssh_se_subcycl) then
      loc=omp_min_max_sum1(d_eta, 1, myDim_nod2D, 'max', partit)
  else
      loc=omp_min_max_sum1(hbar-hbar_old, 1, myDim_nod2D, 'max', partit)
  end if 
  call MPI_AllREDUCE(loc , max_deta  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(hnode(1, :), 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_hnode , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum1(hnode(2, :), 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_hnode2 , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum2(CFL_z, 1, nl-1, 1, myDim_nod2D, 'max', partit) 
  call MPI_AllREDUCE(loc , max_cfl_z , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum2(pgf_x, 1, nl-1, 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_pgfx , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum2(pgf_y, 1, nl-1, 1, myDim_nod2D, 'max', partit)
  call MPI_AllREDUCE(loc , max_pgfy , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  if (use_ice) then
     loc=omp_min_max_sum1(m_ice, 1, myDim_nod2D, 'max', partit)
     call MPI_AllREDUCE(loc , max_m_ice , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  end if
  loc=omp_min_max_sum2(Av, 1, nl, 1, myDim_elem2D, 'max', partit) 
  call MPI_AllREDUCE(loc , max_av , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  loc=omp_min_max_sum2(Kv, 1, nl, 1, myDim_nod2D, 'max', partit) 
  call MPI_AllREDUCE(loc , max_kv , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  !_______________________________________________________________________
  if (mype==0) then
     write(*,*) '___CHECK GLOBAL OCEAN VARIABLES --> mstep=',mstep
     write(*,*) '  ___global estimat of eta & hbar____________________'
     write(*,*) '   int(eta), int(hbar)      =', int_eta, int_hbar
     write(*,*) '    --> error(eta-hbar)     =', int_eta-int_hbar
     write(*,*) '   min(eta) , max(eta)      =', min_eta, max_eta
     write(*,*) '   max(hbar), max(hbar)     =', min_hbar, max_hbar
     write(*,*)
     write(*,*) '   int(deta), int(dhbar)    =', int_deta, int_dhbar
     write(*,*) '    --> error(deta-dhbar)   =', int_deta-int_dhbar
     write(*,*) '    --> error(deta-wflux)   =', int_deta-int_wflux
     write(*,*) '    --> error(dhbar-wflux)  =', int_dhbar-int_wflux
     write(*,*)
     write(*,*) '   -int(wflux)*dt          =', int_wflux*dt*(-1.0)
     write(*,*) '   int(deta )-int(wflux)*dt =', int_deta-int_wflux*dt*(-1.0)
     write(*,*) '   int(dhbar)-int(wflux)*dt =', int_dhbar-int_wflux*dt*(-1.0)
     write(*,*)
     write(*,*) '  ___global min/max/mean  --> mstep=',mstep,'____________'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '       eta= ', min_eta  ,' | ',max_eta  ,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '      deta= ', min_deta ,' | ',max_deta ,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '      hbar= ', min_hbar ,' | ',max_hbar ,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") '     wflux= ', min_wflux,' | ',max_wflux,' | ',int_wflux
     write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") '     hflux= ', min_hflux,' | ',max_hflux,' | ',int_hflux
     write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") '      temp= ', min_temp ,' | ',max_temp ,' | ',int_temp
     write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") '      salt= ', min_salt ,' | ',max_salt ,' | ',int_salt
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '    wvel(1,:)= ', min_wvel ,' | ',max_wvel ,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '    wvel(2,:)= ', min_wvel2,' | ',max_wvel2,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '    uvel(1,:)= ', min_uvel ,' | ',max_uvel ,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '    uvel(2,:)= ', min_uvel2,' | ',max_uvel2,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '    vvel(1,:)= ', min_vvel ,' | ',max_vvel ,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '    vvel(2,:)= ', min_vvel2,' | ',max_vvel2,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '   hnode(1,:)= ', min_hnode,' | ',max_hnode,' | ','N.A.'
     write(*,"(A, ES10.3, A, ES10.3, A, A     )") '   hnode(2,:)= ', min_hnode2,' | ',max_hnode2,' | ','N.A.'
     write(*,"(A, A     , A, ES10.3, A, A     )") '     cfl_z= ',' N.A.     ',' | ',max_cfl_z  ,' | ','N.A.'
     write(*,"(A, A     , A, ES10.3, A, A     )") '     pgf_x= ',' N.A.     ',' | ',max_pgfx  ,' | ','N.A.'
     write(*,"(A, A     , A, ES10.3, A, A     )") '     pgf_y= ',' N.A.     ',' | ',max_pgfy  ,' | ','N.A.'
     write(*,"(A, A     , A, ES10.3, A, A     )") '          Av= ',' N.A.     ',' | ',max_av    ,' | ','N.A.'
     write(*,"(A, A     , A, ES10.3, A, A     )") '          Kv= ',' N.A.     ',' | ',max_kv    ,' | ','N.A.'
     if (use_ice)  then
     write(*,"(A, A     , A, ES10.3, A, A)")      '     m_ice= ',' N.A.     ',' | ',max_m_ice  ,' | ','N.A.'
     end if
   end if
   endif ! --> if (mod(istep,logfile_outfreq)==0) then
end subroutine write_step_info
!
!
!===============================================================================
subroutine check_blowup(istep, ice, dynamics, tracers, partit, mesh)
    USE MOD_ICE
    USE MOD_DYN
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use g_config, only: logfile_outfreq, which_ALE, toy_ocean, use_ice, use_icebergs, ib_num
    use o_PARAM
    use o_ARRAYS, only: water_flux, stress_surf, &
                    heat_flux, Kv, Av
    use g_comm_auto
    use io_BLOWUP
    use g_forcing_arrays
    use diagnostics

    use iceberg_params
    implicit none
  
    type(t_ice)   , intent(in)   , target :: ice
    type(t_dyn)   , intent(in)   , target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(in)   , target :: tracers
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                       :: n, nz, istep, found_blowup_loc=0, found_blowup=0, ib
    integer                       :: el, elidx
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV
    real(kind=WP), dimension(:,:)  , pointer :: Wvel, CFL_z
    real(kind=WP), dimension(:)    , pointer :: ssh_rhs, ssh_rhs_old
    real(kind=WP), dimension(:)    , pointer :: eta_n, d_eta
    real(kind=WP), dimension(:)    , pointer :: u_ice, v_ice
    real(kind=WP), dimension(:)    , pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:)    , pointer :: a_ice_old, m_ice_old, m_snow_old
    real(kind=WP), dimension(:), allocatable, target :: dhbar
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    UV          => dynamics%uv(:,:,:)
    Wvel        => dynamics%w(:,:)
    CFL_z       => dynamics%cfl_z(:,:)
    
    eta_n       => dynamics%eta_n(:)
    u_ice       => ice%uice(:)
    v_ice       => ice%vice(:)
    a_ice       => ice%data(1)%values(:)
    m_ice       => ice%data(2)%values(:)
    m_snow      => ice%data(3)%values(:)
    a_ice_old   => ice%data(1)%values_old(:)
    m_ice_old   => ice%data(2)%values_old(:)
    m_snow_old  => ice%data(3)%values_old(:)
    if ( .not. dynamics%use_ssh_se_subcycl) then 
        d_eta       => dynamics%d_eta(:)
        ssh_rhs     => dynamics%ssh_rhs(:)
        ssh_rhs_old => dynamics%ssh_rhs_old(:)
    else
        allocate(dhbar(myDim_nod2D+eDim_nod2D))
        dhbar = hbar-hbar_old
        d_eta => dhbar
    end if 
    
    !___________________________________________________________________________
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz)
    do n=1, myDim_nod2d       
       !___________________________________________________________________
       ! check ssh
       if ( ((eta_n(n) /= eta_n(n)) .or. eta_n(n)<-10.0 .or. eta_n(n)>10.0 .or. (d_eta(n) /= d_eta(n)) ) ) then
!$OMP CRITICAL
          found_blowup_loc=1
          write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
          write(*,*) ' --STOP--> found eta_n become NaN or <-10.0, >10.0'
          write(*,*) 'mype     = ',mype
          write(*,*) 'mstep    = ',istep
          write(*,*) 'node     = ',n
          write(*,*) 'uln, nln    = ',ulevels_nod2D(n), nlevels_nod2D(n)
          write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
          write(*,*)
          write(*,*) 'eta_n(n)    = ',eta_n(n)
          write(*,*) 'd_eta(n)    = ',d_eta(n)
          write(*,*)
          write(*,*) 'zbar_3d_n   = ',zbar_3d_n(:,n)
          write(*,*) 'Z_3d_n      = ',Z_3d_n(:,n)
          write(*,*)
          if ( .not. dynamics%use_ssh_se_subcycl) then 
            write(*,*) 'ssh_rhs = ',ssh_rhs(n),', ssh_rhs_old = ',ssh_rhs_old(n)
          end if
          write(*,*)
          write(*,*) 'hbar = ',hbar(n),', hbar_old = ',hbar_old(n)
          write(*,*)
          write(*,*) 'wflux = ',water_flux(n)
          write(*,*)
          if (.not. toy_ocean) then
            write(*,*) 'u_wind = ',u_wind(n),', v_wind = ',v_wind(n)
            write(*,*)
            do nz=1,nod_in_elem2D_num(n)
                    write(*,*) 'stress_surf(1:2,',nz,') = ',stress_surf(:,nod_in_elem2D(nz,n))
            end do
          end if
          if (use_ice) then
          write(*,*)
          write(*,*) 'm_ice = ',m_ice(n),', m_ice_old = ',m_ice_old(n)
          write(*,*) 'a_ice = ',a_ice(n),', a_ice_old = ',a_ice_old(n)
          end if 
          write(*,*)
          write(*,*) 'Wvel(:, n)  = ',Wvel(ulevels_nod2D(n):nlevels_nod2D(n),n)
          write(*,*)
          write(*,*) 'CFL_z(:,n)  = ',CFL_z(ulevels_nod2D(n):nlevels_nod2D(n),n)
          write(*,*)
          write(*,*) 'hnode(:, n)  = ',hnode(ulevels_nod2D(n):nlevels_nod2D(n),n)
          write(*,*)
          if (use_icebergs) then
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                if (mesh%elem2d_nodes(1, iceberg_elem(ib)) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
!$OMP END CRITICAL
       endif
       
       !___________________________________________________________________
       ! check surface vertical velocity --> in case of zlevel and zstar 
       ! vertical coordinate its indicator if Volume is conserved  for 
       ! Wvel(1,n)~maschine preccision
       if ( .not. trim(which_ALE)=='linfs' .and. ( Wvel(1, n) /= Wvel(1, n)  )) then
!$OMP CRITICAL
          found_blowup_loc=1
          write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
          write(*,*) ' --STOP--> found surface layer vertical velocity becomes NaN or >1e-12'
          write(*,*) 'mype     = ',mype
          write(*,*) 'mstep    = ',istep
          write(*,*) 'node     = ',n
          write(*,*) 'uln, nln    = ',ulevels_nod2D(n), nlevels_nod2D(n)
          write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
          write(*,*)
          write(*,*) 'Wvel(1, n)  = ',Wvel(1,n)
          write(*,*) 'Wvel(:, n)  = ',Wvel(:,n)
          write(*,*)
          write(*,*) 'hnode(1, n) = ',hnode(1,n)
          write(*,*) 'hnode(:, n) = ',hnode(:,n)
          write(*,*) 'hflux    = ',heat_flux(n)
          write(*,*) 'wflux    = ',water_flux(n)
          write(*,*)
          write(*,*) 'eta_n    = ',eta_n(n)
          write(*,*) 'd_eta(n)    = ',d_eta(n)
          write(*,*) 'hbar     = ',hbar(n)
          write(*,*) 'hbar_old    = ',hbar_old(n)
          if ( .not. dynamics%use_ssh_se_subcycl) then 
            write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
            write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
          end if
          write(*,*)
          write(*,*) 'CFL_z(:,n)  = ',CFL_z(:,n)
          write(*,*)
          if (use_icebergs) then
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                if (mesh%elem2d_nodes(1, iceberg_elem(ib)) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
!$OMP END CRITICAL
       end if ! --> if ( .not. trim(which_ALE)=='linfs' .and. ...
          
       !___________________________________________________________________
       ! check surface layer thinknesss
       if ( .not. trim(which_ALE)=='linfs' .and. ( hnode(1, n) /= hnode(1, n)  .or. hnode(1,n)< 0 )) then
!$OMP CRITICAL
          found_blowup_loc=1
          write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
          write(*,*) ' --STOP--> found surface layer thickness becomes NaN or <0'
          write(*,*) 'mype     = ',mype
          write(*,*) 'mstep    = ',istep
          write(*,*) 'node     = ',n
          write(*,*)
          write(*,*) 'hnode(1, n)  = ',hnode(1,n)
          write(*,*) 'hnode(:, n)  = ',hnode(:,n)
          write(*,*)
          write(*,*) 'eta_n    = ',eta_n(n)
          write(*,*) 'd_eta(n)    = ',d_eta(n)
          write(*,*) 'hbar     = ',hbar(n)
          write(*,*) 'hbar_old    = ',hbar_old(n)
          if ( .not. dynamics%use_ssh_se_subcycl) then 
            write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
            write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
          end if
          write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
          write(*,*)
          if (use_icebergs) then
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                if (mesh%elem2d_nodes(1, iceberg_elem(ib)) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
!$OMP END CRITICAL
       end if ! --> if ( .not. trim(which_ALE)=='linfs' .and. ...
          
       
       do nz=ulevels_nod2D(n),nlevels_nod2D(n)-1
          !_______________________________________________________________
          ! check temp
          if ( (tracers%data(1)%values(nz, n) /= tracers%data(1)%values(nz, n)) .or. &
             tracers%data(1)%values(nz, n) < -5.0 .or. tracers%data(1)%values(nz, n)>60) then
!$OMP CRITICAL
             found_blowup_loc=1
             write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
             write(*,*) ' --STOP--> found temperture becomes NaN or <-5.0, >60'
             write(*,*) 'mype     = ',mype
             write(*,*) 'mstep    = ',istep
             write(*,*) 'node     = ',n
             write(*,*) 'lon,lat     = ',geo_coord_nod2D(:,n)/rad
             write(*,*) 'nz       = ',nz
             write(*,*) 'nzmin, nzmax= ',ulevels_nod2D(n),nlevels_nod2D(n)
             write(*,*) 'x=', geo_coord_nod2D(1,n)/rad, ' ; ', 'y=', geo_coord_nod2D(2,n)/rad
             write(*,*) 'temp(nz, n) = ',tracers%data(1)%values(nz, n)
             write(*,*) 'temp(: , n) = ',tracers%data(1)%values(:, n)
             write(*,*) 'temp_old(nz,n)= ',tracers%data(1)%valuesAB(nz, n)
             write(*,*) 'temp_old(: ,n)= ',tracers%data(1)%valuesAB(:, n)
             write(*,*)
             write(*,*) 'hflux    = ',heat_flux(n)
             write(*,*) 'wflux    = ',water_flux(n)
             write(*,*)
             write(*,*) 'eta_n    = ',eta_n(n)
             write(*,*) 'd_eta(n)    = ',d_eta(n)
             write(*,*) 'hbar     = ',hbar(n)
             write(*,*) 'hbar_old    = ',hbar_old(n)
             if ( .not. dynamics%use_ssh_se_subcycl) then 
                write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
                write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
             end if
             write(*,*)
             if (use_ice) then
                write(*,*) 'm_ice    = ',m_ice(n)
                write(*,*) 'm_ice_old   = ',m_ice_old(n)
                write(*,*) 'm_snow      = ',m_snow(n)
                write(*,*) 'm_snow_old  = ',m_snow_old(n)
                write(*,*)
             end if 
             write(*,*) 'hnode    = ',hnode(:,n)
             write(*,*) 'hnode_new   = ',hnode_new(:,n)
             write(*,*)
             write(*,*) 'Kv       = ',Kv(:,n)
             write(*,*)
             write(*,*) 'W          = ',Wvel(:,n)
             write(*,*)
             write(*,*) 'CFL_z(:,n)  = ',CFL_z(:,n)
             write(*,*)
          if (use_icebergs) then
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                if (mesh%elem2d_nodes(1, iceberg_elem(ib)) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
             write(*,*)
!$OMP END CRITICAL
          endif ! --> if ( (tracers%data(1)%values(nz, n) /= tracers%data(1)%values(nz, n)) .or. & ...
          
          !_______________________________________________________________
          ! check salt
          if ( (tracers%data(2)%values(nz, n) /= tracers%data(2)%values(nz, n)) .or.  &
             tracers%data(2)%values(nz, n) <3.0_WP .or. tracers%data(2)%values(nz, n) >45.0_WP ) then
!$OMP CRITICAL
             found_blowup_loc=1
             write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
             write(*,*) ' --STOP--> found salinity becomes NaN or <=3.0, >=45.0'
             write(*,*) 'mype     = ',mype
             write(*,*) 'mstep    = ',istep
             write(*,*) 'node     = ',n
             write(*,*) 'nz       = ',nz
             write(*,*) 'nzmin, nzmax= ',ulevels_nod2D(n),nlevels_nod2D(n)
             write(*,*) 'x=', geo_coord_nod2D(1,n)/rad, ' ; ', 'y=', geo_coord_nod2D(2,n)/rad
             write(*,*) 'salt(nz, n) = ',tracers%data(2)%values(nz, n)
             write(*,*) 'salt(: , n) = ',tracers%data(2)%values(:, n)
             write(*,*)
             write(*,*) 'temp(nz, n) = ',tracers%data(1)%values(nz, n)
             write(*,*) 'temp(: , n) = ',tracers%data(1)%values(:, n)
             write(*,*)
             write(*,*) 'hflux    = ',heat_flux(n)
             write(*,*)
             write(*,*) 'wflux    = ',water_flux(n)
             write(*,*) 'eta_n    = ',eta_n(n)
             write(*,*) 'd_eta(n)    = ',d_eta(n)
             write(*,*) 'hbar     = ',hbar(n)
             write(*,*) 'hbar_old    = ',hbar_old(n)
             if ( .not. dynamics%use_ssh_se_subcycl) then 
                write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
                write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
             end if 
             write(*,*)
             write(*,*) 'hnode    = ',hnode(:,n)
             write(*,*) 'hnode_new   = ',hnode_new(:,n)
             write(*,*)
             write(*,*) 'zbar_3d_n   = ',zbar_3d_n(:,n)
             write(*,*) 'Z_3d_n      = ',Z_3d_n(:,n)
             write(*,*)
             write(*,*) 'Kv       = ',Kv(:,n)
             write(*,*)
             do el=1,nod_in_elem2d_num(n)
                elidx = nod_in_elem2D(el,n)
                 write(*,*) ' elem#=',el,', elemidx=',elidx
                 write(*,*) '    Av =',Av(:,elidx)
             enddo
             write(*,*) 'Wvel     = ',Wvel(:,n)
             write(*,*)
             write(*,*) 'CFL_z(:,n)  = ',CFL_z(:,n)
             write(*,*)
             write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
             write(*,*)
          if (use_icebergs) then
            write(*,*) 'ibhf_n(:, n) = ',ibhf_n(ulevels_nod2D(n):nlevels_nod2D(n),n)
            write(*,*) 'ibfwb(n) = ',ibfwb(n)
            write(*,*) 'ibfwl(n) = ',ibfwl(n)
            write(*,*) 'ibfwe(n) = ',ibfwe(n)
            write(*,*) 'ibfwbv(n) = ',ibfwbv(n)
            do ib=1, ib_num
                if (mesh%elem2d_nodes(1, iceberg_elem(ib)) == n) then
                    write(*,*) 'ib = ',ib, ', length = ',length_ib(ib), ', height = ', height_ib(ib), ', scaling = ', scaling(ib) 
                    write(*,*) 'hfb_flux_ib(ib) = ',hfb_flux_ib(ib)
                    write(*,*) 'hfl_flux_ib(ib,n) = ',hfl_flux_ib(ib,n)
                    write(*,*) 'hfe_flux_ib(ib) = ',hfe_flux_ib(ib)
                    write(*,*) 'hfbv_flux_ib(ib,n) = ',hfbv_flux_ib(ib,n)
                end if
            end do
            write(*,*)
          end if
!$OMP END CRITICAL
          endif ! --> if ( (tracers%data(2)%values(nz, n) /= tracers%data(2)%values(nz, n)) .or.  & ...
       end do ! --> do nz=1,nlevels_nod2D(n)-1
    end do ! --> do n=1, myDim_nod2d
!$OMP END PARALLEL DO
    !_______________________________________________________________________
    ! check globally if one of the cpus hat a blowup situation. if its the
    ! case CPU mype==0 needs to write out the stuff. Write out occurs in 
    ! moment only over CPU mype==0
    call MPI_AllREDUCE(found_blowup_loc  , found_blowup  , 1, MPI_INTEGER, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    if (found_blowup==1) then
        call write_step_info(istep, 1, ice, dynamics, tracers, partit, mesh)
        if (mype==0) then
            call sleep(1)
            write(*,*)
            write(*,*) '                      ,-*                 ,-*             '
            write(*,*) '                     (_)  MODEL BLOW UP  (_)              '
            write(*,*) '                              ____                        '
            write(*,*) '                       __,-~~/~   `---.                   '
            write(*,*) '                     _/_,---(      ,   )                  '
            write(*,*) '                 __ /        <   /   )   \___             '
            write(*,*) ' - -- ----===;;;`====------------------===;;;===---- -- - '
            write(*,*) '                    \/  ~"~"~"~"~"~\~"~)~"/               '
            write(*,*) '                    (_ (   \  (     >    \)               '
            write(*,*) '                     \_( _ <         >_>`                 '
            write(*,*) '                        ~ `-i` ::>|--"                    '
            write(*,*) '                            I;|.|.|                       '
            write(*,*) '                           <|i::|i|`                      '
            write(*,*) '                          (` ^`"`- ")                     '
            write(*,*) ' _______________________.,-#%&$@%#&#~,.__________________ '
            write(*,*) '                                                          '
            write(*,*) '            (`- ́)  _ (`- ́).->          <-. (`- ́)          ' 
            write(*,*) '   <-.      ( OO).-/ ( OO)_      .->      \(OO )_         '
            write(*,*) '(`- ́)-----.(,------.(_)--\_)(`- ́)----. ,--./  ,-.) .----. '
            write(*,*) '(OO|(_\--- ́ |  .--- ́/    _ /( OO).-.  `|   `. ́   |\_,-.  |'
            write(*,*) ' / |  `--. (|  `--. \_..`--.( _) | |  ||  |`. ́|  |   . ́ . ́'
            write(*,*) ' \_)  .-- ́  |  .-- ́ .-._)   \\|  |)|  ||  |   |  | . ́  /_ '
            write(*,*) '  `|  |_)   |  `---.\       / `  `- ́   ́|  |   |  ||      |'
            write(*,*) '   `-- ́     `------ ́ `----- ́   `----- ́ `-- ́   `-- ́`------ ́'
            write(*,*)
        end if
        call blowup(istep, ice, dynamics, tracers, partit, mesh)
        if (mype==0) write(*,*) ' --> finished writing blow up file'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    endif 
end subroutine check_blowup
!===============================================================================
subroutine write_enegry_info(dynamics, partit, mesh)
   use MOD_MESH
   USE MOD_PARTIT
   USE MOD_PARSUP
   use MOD_DYN
   use g_support
   IMPLICIT NONE
   type(t_mesh),   intent(in)   , target :: mesh
   type(t_partit), intent(inout), target :: partit
   type(t_dyn)   , intent(in)   , target :: dynamics
   real(kind=WP)                         :: budget(2)
   integer, pointer                      :: mype

   mype            => partit%mype

   if (mype==0) write(*,*) '*******KE budget analysis...*******'
   if (mype==0) write(*,*) '     U     |     V     |     TOTAL'
   call integrate_elem(dynamics%ke_du2(1,:,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_du2(2,:,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. du2=', budget(1),' | ',  budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_pre_xVEL(1,:,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_pre_xVEL(2,:,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. pre=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_nod(dynamics%ke_wrho,  budget(1), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7)") 'w * rho=', budget(1)

   call integrate_elem(dynamics%ke_adv_xVEL(1,:,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_adv_xVEL(2,:,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. adv=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_hvis_xVEL(1,:,:), budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_hvis_xVEL(2,:,:), budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke.  ah=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_vvis_xVEL(1,:,:), budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_vvis_xVEL(2,:,:), budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke.  av=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_cor_xVEL(1,:,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_cor_xVEL(2,:,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. cor=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_wind_xVEL(1,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_wind_xVEL(2,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. wind=', budget(1), ' | ', budget(2), ' | ', sum(budget)

   call integrate_elem(dynamics%ke_drag_xVEL(1,:),  budget(1), partit, mesh)
   call integrate_elem(dynamics%ke_drag_xVEL(2,:),  budget(2), partit, mesh)
   if (mype==0) write(*,"(A, ES14.7, A, ES14.7, A, ES14.7)") 'ke. drag=', budget(1), ' | ', budget(2), ' | ', sum(budget)
   if (mype==0) write(*,*) '***********************************'   
end subroutine write_enegry_info

end module write_step_info_module
