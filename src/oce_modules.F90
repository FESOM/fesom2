! Modules of cell-vertex ocean model
! S. Danilov, 2012 (sergey.danilov@awi.de)
! SI units are used
  
MODULE o_PARAM
integer, parameter            :: WP=8        ! Working precision
real(kind=WP), parameter      :: pi=3.14159265358979
real(kind=WP), parameter      :: rad=pi/180.0_WP
real(kind=WP), parameter      :: density_0=1030.0_WP
real(kind=WP), parameter      :: g=9.81_WP
real(kind=WP), parameter      :: r_earth=6367500.0_WP
real(kind=WP), parameter      :: omega=2*pi/(3600.0_WP*24.0_WP)
real(kind=WP), parameter      :: vcpw=4.2e6   ![J/m^3/K] water heat cap
real(kind=WP), parameter      :: small=1.0e-8 !small value

real(kind=WP)                 :: C_d= 0.0025_WP ! Bottom drag coefficient
real(kind=WP)	              :: kappa=0.4      !von Karman's constant
real(kind=WP)                 :: mix_coeff_PP=0.01_WP   ! mixing coef for PP scheme
real(kind=WP)                 :: A_hor=100.0_WP		! Horizontal harm. visc.    
real(kind=WP)                 :: A_hor_max=1500.0_WP	! Maximum viscosity allowed (to limit Smag and Leith
							! contributions when they are too large
real(kind=WP)                 :: Leith_c=0.		! Leith viscosity. It needs vorticity, which is only computed for 
							! the vector invariant form of momentum advection (mom_adv=4)
real(kind=WP)                 :: tau_c=0.4		! Controls the strength of filters. Should be about 0.4
real(kind=WP)                 :: A_ver=0.001_WP ! Vertical harm. visc.
real(kind=WP)                 :: Div_c=0.5_WP
real(kind=WP)                 :: Smag_c=0.0_WP  ! 0.2   ! (C/pi)^2
real(kind=WP)                 :: Abh0=8.0e12    ! 
logical                       :: laplacian=.false.
logical                       :: biharmonic=.true.
real(kind=WP)                 :: K_hor=10._WP
real(kind=WP)                 :: K_ver=0.00001_WP
real(kind=WP)                 :: scale_area=2.0e8
real(kind=WP)                 :: surf_relax_T= 0.0_WP
real(kind=WP)                 :: surf_relax_S= 10.0_WP/(60*3600.0_WP*24)
real(kind=WP)                 :: clim_relax= 1.0_WP/(10*3600.0_WP*24)
real(kind=WP)                 :: clim_decay, clim_growth
                                 ! set to 0.0 if no relaxation
logical                       :: ref_sss_local=.false.
real(kind=WP)                 :: ref_sss=34.7
logical                       :: Fer_GM =.true.   !flag for Ferrari et al. (2010) GM scheme
! Time stepping                               
#ifdef BTR_SPLIT
real(kind=WP)                 :: alpha=0.0_WP, theta=0.0_WP ! implicitness for
#else 
real(kind=WP)                 :: alpha=1.0_WP, theta=1.0_WP ! implicitness for
#endif
                                                 ! elevation and divergence
real(kind=WP)                 :: epsilon=0.1_WP ! AB2 offset 
! Tracers
logical                       :: i_vert_diff= .true.
logical                       :: i_vert_visc= .true.
integer                       :: tracer_adv=1

logical                       :: w_split  =.false.
real(kind=WP)                 :: w_exp_max=1.e-5_WP

                                ! 1 Muira
				! 2 quadratic reconstruction Miura
				! 3 MUSCL
				! 4 MUSCL-FCT
integer	                       :: num_tracers=2
! Momentum
logical                       :: free_slip=.false.
                                ! false=no slip 
integer                       :: mom_adv=4 
                                ! 1 vector control volumes
				! 2 vector c. v., p1 velocities
				! 3 scalar control volumes  
				! 4 vector invariant 

logical                       :: open_b=.false.   ! Reserved    

logical                       :: mo_on=.false. !Monin-Obukhov
real*8 :: modiff=0.01                   !for PP, mixing coefficient within MO length

  ! *** active tracer cutoff
logical          :: limit_salinity=.true.         !set an allowed range for salinity
real(kind=WP)    :: salinity_min=5.0              !minimal salinity 
real(kind=WP)    :: coeff_limit_salinity=0.0023   !m/s, coefficient to restore s to s_min

  namelist /tracer_cutoff/ limit_salinity, salinity_min, coeff_limit_salinity

! *** others ***
 integer                       :: num_tracer


 NAMELIST /oce_dyn/ C_d, A_ver, laplacian, A_hor, A_hor_max, Leith_c, tau_c, Div_c, Smag_c, &
                    biharmonic, Abh0, scale_area, mom_adv, free_slip, i_vert_visc, w_split, w_exp_max
 NAMELIST /oce_tra/ K_ver, K_hor, surf_relax_T, surf_relax_S, clim_relax, &
		    ref_sss_local, ref_sss, i_vert_diff, &
		    tracer_adv
END MODULE o_PARAM  
!==========================================================
MODULE o_MESH
USE o_PARAM
! All variables used to keep the mesh structure +
! auxiliary variables involved in implementation 
! of open boundaries and advection schemes
! 

integer                                    ::   nod2D      ! the number of 2D nodes
real(kind=WP)                              ::   ocean_area
real(kind=WP), allocatable, dimension(:,:) ::   coord_nod2D, geo_coord_nod2D
integer                                    ::   edge2D     ! the number of 2D edges
integer                                    ::   edge2D_in  
                                              ! the number of internal 2D edges
integer                                    ::   elem2D     ! the number of 2D elements
integer, allocatable, dimension(:,:)       ::   elem2D_nodes
                                              ! elem2D_nodes(:,n) lists
				              ! 3 nodes of element n   
integer, allocatable, dimension(:,:)       ::   edges
                                              ! edge(:,n) lists 2 nodes
				              ! edge n
integer, allocatable, dimension(:,:)       ::   edge_tri
                                              ! edge_tri(:,n) lists 2 
				              ! elements containing edge n
				              ! The first one is to left 
				              ! of the line directed
				              ! to the second node
integer, allocatable, dimension(:,:)       ::   elem_edges
                                              ! elem_edges(:,n) are edges of 
                                              ! element n.  
real(kind=WP), allocatable, dimension(:)   ::   elem_area
real(kind=WP), allocatable, dimension(:,:) ::   edge_dxdy, edge_cross_dxdy
real(kind=WP), allocatable, dimension(:)   ::   elem_cos, metric_factor
integer,allocatable,dimension(:,:)         ::   elem_neighbors
integer,allocatable,dimension(:,:)         ::   nod_in_elem2D
integer,allocatable,dimension(:)           ::   nod_in_elem2D_num
real(kind=WP),allocatable,dimension(:)     ::   depth
                                              ! depth(n) is the depths at 
				              ! node n 
real(kind=WP),allocatable,dimension(:,:)    ::   gradient_vec 
                                              ! Coefficients of linear reconstruction
					      ! of velocities on elements
real(kind=WP),allocatable,dimension(:,:)    ::   gradient_sca
                                              ! Coefficients to compute
					      ! gradient of scalars on elements
! Vertical structure             
integer                                    :: nl
real(kind=8), allocatable, dimension(:)    :: zbar, Z,elem_depth
integer, allocatable, dimension(:)         :: nlevels, nlevels_nod2D
real(kind=8), allocatable, dimension(:,:)  :: area
real(kind=WP), allocatable, dimension(:)   :: mesh_resolution


  type sparse_matrix 
     integer :: nza
     integer :: dim
     real(kind=8), allocatable, dimension(:)      :: values
     integer(KIND=4), allocatable,   dimension(:) :: colind
     integer(KIND=4), allocatable,   dimension(:) :: rowptr
  end type sparse_matrix
! Elevation stiffness matrix
type(sparse_matrix)                           :: ssh_stiff

! Auxiliary arrays. They are not related to mesh structure, but are 
! kept here because they are just used for temporary storage in computations

! Open boundary:
integer                                       :: ob_num  ! number of OB fragments

TYPE ob_type
    integer      :: len
    integer, allocatable, dimension(:)       :: list
END TYPE ob_type

TYPE ob_rhs_type
    integer      :: len
    real(kind=WP), allocatable, dimension(:) :: list
END TYPE ob_rhs_type

type(ob_type), allocatable                    ::  ob_info(:)
type(ob_rhs_type), allocatable                ::  ob_2rhs(:)
!
! The fct part
real(kind=WP),allocatable,dimension(:,:)      :: fct_aec  ! Antidif. elem. contrib.
real(kind=WP),allocatable,dimension(:,:)      :: fct_LO, fct_HO  ! Low-order solution
real(kind=WP),allocatable,dimension(:,:)      :: fct_aec_ver  ! Antidif. vert. fluxes
real(kind=WP),allocatable,dimension(:,:)      :: fct_ttf_max,fct_ttf_min
real(kind=WP),allocatable,dimension(:,:)      :: fct_plus,fct_minus
! Quadratic reconstruction part
integer,allocatable,dimension(:)              :: nlevels_nod2D_min, nn_num
real(kind=WP),allocatable,dimension(:,:,:)    :: quad_int_mat, quad_int_coef
integer,allocatable,dimension(:,:)            :: nn_pos
! MUSCL type reconstruction
integer,allocatable,dimension(:,:)            :: edge_up_dn_tri
real(kind=WP),allocatable,dimension(:,:,:)    :: edge_up_dn_grad


end module o_MESH
!==========================================================
!==========================================================
module g_PARSUP
USE o_PARAM
! Variables to organize parallel work  
implicit none
save

#ifdef PETSC
#include "finclude/petsc.h"
#else
  include 'mpif.h'
#endif

  type com_struct
     integer    :: rPEnum                    ! the number of PE I receive info from 
     integer, dimension(:), allocatable :: rPE   ! their list
     integer, dimension(:), allocatable :: rptr  ! allocatables to the list of nodes
     integer, dimension(:), allocatable :: rlist ! the list of nodes
     integer    :: sPEnum                    ! send part 
     integer, dimension(:), allocatable :: sPE
     integer, dimension(:), allocatable :: sptr
     integer, dimension(:), allocatable :: slist
  end type com_struct
  type com_array_i
     integer, dimension(:), allocatable :: array
  end  type com_array_i

  type(com_struct)   :: com_nod2D
  type(com_struct)   :: com_edge2D
  type(com_struct)   :: com_elem2D
  type(com_struct)   :: com_elem2D_full
  ! Buffer arrays to store information to be communicated
  type com_array
     real(kind=WP), dimension(:), allocatable :: array
  end  type com_array
  type(com_array), allocatable             :: s_buff_nod2D(:),  r_buff_nod2D(:)
  type(com_array), allocatable             :: s_buff_nod3D(:),  r_buff_nod3D(:)
  type(com_array), allocatable             :: s_buff_elem2D(:), r_buff_elem2D(:)
  type(com_array), allocatable             :: s_buff_elem3D(:), r_buff_elem3D(:)
  type(com_array), allocatable             :: s_buff_elem3D_full(:),s_buff_elem2D_full(:)
  type(com_array), allocatable             :: r_buff_elem3D_full(:),r_buff_elem2D_full(:)
  type(com_array), allocatable             :: s_buff_edge2D(:), r_buff_edge2D(:)
  type(com_array), allocatable             :: s_buff_edge3D(:), r_buff_edge3D(:)
  type(com_array_i), allocatable :: s_buff_elem2D_full_i(:), r_buff_elem2D_full_i(:) 
  type(com_array_i), allocatable :: s_buff_nod2D_i(:),  r_buff_nod2D_i(:)
  ! general MPI part
  integer            :: MPIERR
  integer            :: npes
  integer            :: mype
  integer            :: maxPEnum=100
  integer, allocatable, dimension(:)  :: part

  ! Mesh partition
  integer                             :: myDim_nod2D, eDim_nod2D
  integer, allocatable, dimension(:)  :: myList_nod2D
  integer                             :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
  integer, allocatable, dimension(:)  :: myList_elem2D
  integer                             :: myDim_edge2D, eDim_edge2D
  integer, allocatable, dimension(:)  :: myList_edge2D

  integer :: pe_status = 0 ! if /=0 then something is wrong 

   integer, allocatable ::  remPtr_nod2D(:),  remList_nod2D(:)
   integer, allocatable ::  remPtr_elem2D(:), remList_elem2D(:)

contains
subroutine par_init    ! initializes MPI

  implicit none


  integer :: i

  call MPI_INIT(i);
  call MPI_Comm_Size(MPI_COMM_WORLD,npes,i)
  call MPI_Comm_Rank(MPI_COMM_WORLD,mype,i)
#ifdef PETSC
  call PETSCInitialize(PETSC_NULL_CHARACTER,i);
#endif
  if(mype==0) then
  write(*,*) 'MPI has been initialized'
  write(*, *) 'Running on ', npes, ' PEs'
  end if
end subroutine par_init
!=================================================================
subroutine par_ex(abort)       ! finalizes MPI
  implicit none
  integer,optional :: abort
  if (present(abort)) then
   write(*,*) 'Run finished unexpectedly!'
   call MPI_ABORT( MPI_COMM_WORLD, 1 )
  else
  call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call  MPI_Finalize(MPIerr)
  endif
end subroutine par_ex
!=================================================================
subroutine set_par_support_ini
  use o_MESH
  implicit none

  integer   n, j, k, nini, nend

  ! Construct partitioning vector
  allocate(part(nod2D))
  part=0  
  !do n=0, npes-1            ! this works only on geometrically simple meshes
  !   nini=(nod2D/npes)*n+1
  !   nend=(nod2D/npes)*(n+1)
  !   if (n==npes-1) nend=nod2D
  !   part(nini:nend)=n
  !end do
  if(npes>1) then
#ifdef PART_WEIGHTED
  call partit2(ssh_stiff%dim,ssh_stiff%rowptr,ssh_stiff%colind,nlevels_nod2D,npes, part)
#else
  call partit (ssh_stiff%dim,ssh_stiff%rowptr,ssh_stiff%colind, npes, part)
#endif
  if(mype==0) write(*,*) 'partitioning is done'
  end if
  call communication_edgen
  call communication_nodn
  call communication_elemn
  call communication_elem_fulln
  call mymesh
  if(mype==0) write(*,*) 'Communication arrays have been set up'   
end subroutine set_par_support_ini

!=======================================================================

subroutine set_par_support
  use o_MESH
  implicit none

  integer   n, offset
  !
  ! In the distributed memory version, most of the job is already done 
  ! at the initialization phase and is taken into account in read_mesh
  ! routine. Here only communication buffers are set. 
  ! 
  !=============
  ! Allocate communication buffers: 
  !=============
  if (npes>1) then
      !=========
      ! Send and receive buffers for nodal arrays 
      !=========
      ! 2D
     allocate(s_buff_nod2D(com_nod2D%sPEnum),r_buff_nod2D(com_nod2D%rPEnum))
     do n=1, com_nod2D%sPEnum
        offset=com_nod2D%sptr(n+1) - com_nod2D%sptr(n)
        allocate(s_buff_nod2D(n)%array(offset))
     end do
     do n=1, com_nod2D%rPEnum
        offset=com_nod2D%rptr(n+1) - com_nod2D%rptr(n)
        allocate(r_buff_nod2D(n)%array(offset))
     end do
      allocate(s_buff_nod2D_i(com_nod2D%sPEnum),r_buff_nod2D_i(com_nod2D%rPEnum))                     
      do n=1, com_nod2D%sPEnum                                                                        
        offset=com_nod2D%sptr(n+1) - com_nod2D%sptr(n)                                                
        allocate(s_buff_nod2D_i(n)%array(offset))                                                     
      end do                                                                                          
      do n=1, com_nod2D%rPEnum                                                                        
        offset=com_nod2D%rptr(n+1) - com_nod2D%rptr(n)                                                
        allocate(r_buff_nod2D_i(n)%array(offset))                                                     
      end do                                                                                          
              
      ! 3D
     allocate(s_buff_nod3D(com_nod2D%sPEnum),r_buff_nod3D(com_nod2D%rPEnum))
     do n=1, com_nod2D%sPEnum
        offset=com_nod2D%sptr(n+1) - com_nod2D%sptr(n)
	!allocate(s_buff_nod3D(n)%array(offset*(nl-1)))
	allocate(s_buff_nod3D(n)%array(offset*nl))
     end do
     do n=1, com_nod2D%rPEnum
        offset=com_nod2D%rptr(n+1) - com_nod2D%rptr(n)
	!allocate(r_buff_nod3D(n)%array(offset*(nl-1)))
	allocate(r_buff_nod3D(n)%array(offset*nl))
     end do
      !=========
      ! Send and receive buffers for elemental arrays
      !========= 
        ! 2D
     allocate(s_buff_elem2D(com_elem2D%sPEnum),r_buff_elem2D(com_elem2D%rPEnum))
     do n=1, com_elem2D%sPEnum
        offset=com_elem2D%sptr(n+1) - com_elem2D%sptr(n)
	allocate(s_buff_elem2D(n)%array(offset))
     end do
     do n=1, com_elem2D%rPEnum
        offset=com_elem2D%rptr(n+1) - com_elem2D%rptr(n)
	allocate(r_buff_elem2D(n)%array(offset))
     end do
      allocate(s_buff_elem2D_full(com_elem2D_full%sPEnum))
      allocate(r_buff_elem2D_full(com_elem2D_full%rPEnum))
      do n=1, com_elem2D_full%sPEnum
        offset=com_elem2D_full%sptr(n+1) - com_elem2D_full%sptr(n)
        allocate(s_buff_elem2D_full(n)%array(offset))
      end do
      do n=1, com_elem2D_full%rPEnum
        offset=com_elem2D_full%rptr(n+1) - com_elem2D_full%rptr(n)
        allocate(r_buff_elem2D_full(n)%array(offset))
      end do
       ! 3D
     allocate(s_buff_elem3D(com_elem2D%sPEnum),r_buff_elem3D(com_elem2D%rPEnum))
     do n=1, com_elem2D%sPEnum
        offset=com_elem2D%sptr(n+1) - com_elem2D%sptr(n)
	allocate(s_buff_elem3D(n)%array(offset*(nl-1)))
	!allocate(s_buff_elem3D(n)%array(offset*nl))
     end do
     do n=1, com_elem2D%rPEnum
        offset=com_elem2D%rptr(n+1) - com_elem2D%rptr(n)
	allocate(r_buff_elem3D(n)%array(offset*(nl-1)))
	!allocate(r_buff_elem3D(n)%array(offset*nl))
     end do
      allocate(s_buff_elem2D_full_i(com_elem2D_full%sPEnum))
      allocate(r_buff_elem2D_full_i(com_elem2D_full%rPEnum))
      do n=1, com_elem2D_full%sPEnum
        offset=com_elem2D_full%sptr(n+1) - com_elem2D_full%sptr(n)
        allocate(s_buff_elem2D_full_i(n)%array(offset))
      end do
      do n=1, com_elem2D_full%rPEnum
        offset=com_elem2D_full%rptr(n+1) - com_elem2D_full%rptr(n)
        allocate(r_buff_elem2D_full_i(n)%array(offset))
      end do

       ! 3D full
    
     allocate(s_buff_elem3D_full(com_elem2D_full%sPEnum), &
              r_buff_elem3D_full(com_elem2D_full%rPEnum))
     do n=1, com_elem2D_full%sPEnum
        offset=com_elem2D_full%sptr(n+1) - com_elem2D_full%sptr(n)
	allocate(s_buff_elem3D_full(n)%array(offset*(nl-1)))
     end do
     do n=1, com_elem2D_full%rPEnum
        offset=com_elem2D_full%rptr(n+1) - com_elem2D_full%rptr(n)
	allocate(r_buff_elem3D_full(n)%array(offset*(nl-1)))
     end do
     
      !=========
      ! Send and receive buffers for edge arrays
      !========= 
      ! 2D
     
     allocate(s_buff_edge2D(com_edge2D%sPEnum),r_buff_edge2D(com_edge2D%rPEnum))
     do n=1, com_edge2D%sPEnum
        offset=com_edge2D%sptr(n+1) - com_edge2D%sptr(n)
	allocate(s_buff_edge2D(n)%array(offset))
     end do
     do n=1, com_edge2D%rPEnum
        offset=com_edge2D%rptr(n+1) - com_edge2D%rptr(n)
        allocate(r_buff_edge2D(n)%array(offset))
     end do
     
   end if

   call init_gatherLists

  if(mype==0) write(*,*) 'Communication arrays are set' 

end subroutine set_par_support


!===================================================================
subroutine init_gatherLists

  use o_MESH
  implicit none
  
  integer :: n2D, e2D, sum_loc_elem2D
  integer :: n, estart, nstart

  if (mype==0) then

     if (npes > 1) then

        allocate(remPtr_nod2D(npes))
        allocate(remPtr_elem2D(npes))

        remPtr_nod2D(1) = 1
        remPtr_elem2D(1) = 1
        
        do n=1, npes-1
           call MPI_RECV(n2D, 1, MPI_INTEGER, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPIerr )
           call MPI_RECV(e2D, 1, MPI_INTEGER, n, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPIerr )

           remPtr_nod2D(n+1)  = remPtr_nod2D(n)  + n2D
           remPtr_elem2D(n+1) = remPtr_elem2D(n) + e2D 
        enddo



        allocate(remList_nod2D(remPtr_nod2D(npes)))   ! this should be nod2D - myDim_nod2D
        allocate(remList_elem2D(remPtr_elem2D(npes))) ! this is > elem2D, because the elements overlap.
                                                      ! Consider optimization: avoid multiple communication
                                                      ! of the same elem from different PEs.

        do n=1, npes-1
           nstart = remPtr_nod2D(n)
           n2D    = remPtr_nod2D(n+1) - remPtr_nod2D(n) 
           call MPI_RECV(remList_nod2D(nstart), n2D, MPI_INTEGER, n, 2, MPI_COMM_WORLD, &
                                                           MPI_STATUS_IGNORE, MPIerr ) 
           estart = remPtr_elem2D(n)
           e2D    = remPtr_elem2D(n+1) - remPtr_elem2D(n)
           call MPI_RECV(remList_elem2D(estart),e2D, MPI_INTEGER, n, 3, MPI_COMM_WORLD, &
                                                           MPI_STATUS_IGNORE, MPIerr ) 

        enddo
     end if
  else

     call MPI_SEND(myDim_nod2D,   1,            MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
     call MPI_SEND(myDim_elem2D,  1,            MPI_INTEGER, 0, 1, MPI_COMM_WORLD, MPIerr )
     call MPI_SEND(myList_nod2D,  myDim_nod2D,  MPI_INTEGER, 0, 2, MPI_COMM_WORLD, MPIerr )
     call MPI_SEND(myList_elem2D, myDim_elem2D, MPI_INTEGER, 0, 3, MPI_COMM_WORLD, MPIerr )

  endif

end subroutine init_gatherLists


end module g_PARSUP

!==========================================================
MODULE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
! Arrays are described in subroutine array_setup  
real(kind=WP), allocatable    :: UV(:,:,:)
real(kind=WP), allocatable    :: UV_rhs(:,:,:), UV_rhsAB(:,:,:)
real(kind=WP), allocatable    :: eta_n(:), d_eta(:)
real(kind=WP), allocatable    :: ssh_rhs(:), Wvel(:,:), hpressure(:,:)
real(kind=WP), allocatable    :: Wvel_e(:,:), Wvel_i(:,:)
real(kind=WP), allocatable    :: stress_surf(:,:)
real(kind=WP), allocatable    :: T_rhs(:,:) 
real(kind=WP), allocatable    :: heat_flux(:), Tsurf(:) 
real(kind=WP), allocatable    :: S_rhs(:,:)
real(kind=WP), allocatable    :: tr_arr(:,:,:),tr_arr_old(:,:,:)
real(kind=WP), allocatable    :: del_ttf(:,:)

real(kind=WP), allocatable    :: water_flux(:), Ssurf(:)
real(kind=WP), allocatable    :: Tclim(:,:), Sclim(:,:)
real(kind=WP), allocatable    :: Visc(:,:)
real(kind=WP), allocatable    :: Tsurf_t(:,:), Ssurf_t(:,:)
real(kind=WP), allocatable    :: tau_x_t(:,:), tau_y_t(:,:) 
real(kind=WP), allocatable    :: heat_flux_t(:,:), heat_rel_t(:,:), heat_rel(:) 
real(kind=WP), allocatable    :: coriolis(:), coriolis_node(:)
real(kind=WP), allocatable    :: relax2clim(:)

! Passive and age tracers
 character(4), allocatable    :: prog_tracer_name(:)
real(kind=WP), allocatable    :: tracer(:,:,:), tracer_rhs(:,:,:)   
!Tracer gradients&RHS      
real(kind=8), allocatable :: ttrhs(:,:)
real(kind=8), allocatable :: tr_xy(:,:,:)
! Auxiliary arrays to store velocity gradients
real(kind=WP),allocatable,dimension(:,:,:)    ::   vel_grad

! Auxiliary arrays for vector-invariant form of momentum advection
real(kind=WP), allocatable,dimension(:,:)   :: vorticity

!Viscosity and diff coefs
real(kind=WP), allocatable,dimension(:,:)   :: Av,Kv
!Velocities interpolated to nodes
real(kind=WP), allocatable,dimension(:,:,:)   :: Unode

! Auxiliary arrays to store Redi-GM fields
real(kind=WP), allocatable,dimension(:,:,:) :: neutral_slope
real(kind=WP), allocatable,dimension(:,:,:) :: sigma_xy
real(kind=WP), allocatable,dimension(:,:)   :: sw_beta, sw_alpha
!real(kind=WP), allocatable,dimension(:,:,:) :: tsh, tsv, tsh_nodes
!real(kind=WP), allocatable,dimension(:,:)   :: hd_flux,vd_flux
!Array for Redi-GM coefs
real(kind=WP), allocatable :: Kd(:,:,:)


!Mean arrays
real(kind=WP), allocatable    :: UV_mean(:,:,:)
real(kind=WP), allocatable    :: eta_n_mean(:),Wvel_mean(:,:)
real(kind=WP), allocatable    :: tr_arr_mean(:,:,:)
real(kind=WP), allocatable    :: fer_UV_mean(:,:,:), fer_wvel_mean(:,:)

!Monin-Obukhov correction
real*8,allocatable :: mo(:,:),mixlength(:)
!GM_stuff
real*8,allocatable :: bvfreq(:,:),mixlay_dep(:),bv_ref(:)

real(kind=WP), allocatable    :: fer_UV(:,:,:), fer_wvel(:,:)
real(kind=WP), target, allocatable    :: fer_c(:), fer_K(:), fer_gamma(:,:,:)
END MODULE o_ARRAYS
!==========================================================
