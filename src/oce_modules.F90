! Modules of cell-vertex ocean model
! S. Danilov, 2012 (sergey.danilov@awi.de)
! SI units are used
  
MODULE o_PARAM
integer, parameter            :: WP=8        ! Working precision
integer		                  :: mstep
real(kind=WP), parameter      :: pi=3.14159265358979
real(kind=WP), parameter      :: rad=pi/180.0_WP
real(kind=WP), parameter      :: density_0=1030.0_WP
real(kind=WP), parameter      :: density_0_r=1.0_WP/density_0 ! [m^3/kg]         
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
logical                       :: Fer_GM =.false.   !flag for Ferrari et al. (2010) GM scheme
real(kind=WP)                 :: visc_sh_limit=5.0e-3      !for KPP, max visc due to shear instability
real(kind=WP)                 :: diff_sh_limit=5.0e-3      !for KPP, max diff due to shear instability
logical                       :: Kv0_const=.true.		    !use Kv0 varying with depth and latitude 
logical                       :: double_diffusion=.false.  !for KPP,dd switch
                                 ! KPP parametrization
 character(5)                 :: mix_scheme='KPP'	   !'KPP','PP'
real(KIND=WP)                 :: Ricr   = 0.3_WP  ! critical bulk Richardson Number
real(KIND=WP)                 :: concv  = 1.6_WP  ! constant for pure convection (eqn. 23) (Large 1.5-1.6; MOM default 1.8)

logical                       :: AvKv =.false.   ! write Av, Kv
logical                       :: hbl_diag =.false.   ! write boundary layer depth

! Time stepping                               
! real(kind=WP)                 :: alpha=1.0_WP, theta=1.0_WP ! implicitness for
real(kind=WP)                 :: alpha=1.0_WP, theta=1.0_WP ! implicitness for
                                                 ! elevation and divergence
real(kind=WP)                 :: epsilon=0.01_WP ! AB2 offset 
! Tracers
logical                       :: i_vert_diff= .true.
logical                       :: i_vert_visc= .true.
integer                       :: tracer_adv=1

logical                       :: w_split  =.false.
real(kind=WP)                 :: w_exp_max=1.e-5_WP

				! 1 MUSCL
				! 2 MUSCL-FCT
integer	                       :: num_tracers=2
! Momentum
logical                       :: free_slip=.false.
                                ! false=no slip 
integer                       :: mom_adv=2
                                ! 1 vector control volumes, p1 velocities
				! 2 scalar control volumes  
				! 3 vector invariant 

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
 real*8                        :: time_sum=0.0 ! for runtime estimate


 NAMELIST /oce_dyn/ C_d, A_ver, laplacian, A_hor, A_hor_max, Leith_c, tau_c, Div_c, Smag_c, &
                    biharmonic, Abh0, scale_area, mom_adv, free_slip, i_vert_visc, w_split, w_exp_max, &
                    Fer_GM, visc_sh_limit, mix_scheme, Ricr, concv, AvKv,hbl_diag
 NAMELIST /oce_tra/ diff_sh_limit, Kv0_const, double_diffusion, K_ver, K_hor, surf_relax_T, surf_relax_S, clim_relax, &
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
real(kind=8), allocatable, dimension(:,:)  :: area, area_inv
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
     integer, dimension(:), allocatable :: req  ! request for MPI_Wait
  end type com_struct

  type(com_struct)   :: com_nod2D
  type(com_struct)   :: com_edge2D
  type(com_struct), target :: com_elem2D
  type(com_struct), target :: com_elem2D_full
 
  ! MPI Datatypes for interface exchange

  ! Edge fields (2D)
  integer, allocatable       :: s_mpitype_edge2D(:),         r_mpitype_edge2D(:)   

  ! Element fields (2D; 2D integer; 3D with nl-1 or nl levels, 1 - 4 values)
  !                 small halo and / or full halo
  integer, allocatable, target :: s_mpitype_elem2D(:,:),       r_mpitype_elem2D(:,:) 
  integer, allocatable         :: s_mpitype_elem2D_full_i(:),  r_mpitype_elem2D_full_i(:) 
  integer, allocatable, target :: s_mpitype_elem2D_full(:,:),  r_mpitype_elem2D_full(:,:) 
  integer, allocatable, target :: s_mpitype_elem3D(:,:,:),     r_mpitype_elem3D(:,:,:) 
  integer, allocatable, target :: s_mpitype_elem3D_full(:,:,:),r_mpitype_elem3D_full(:,:,:) 

  ! Nodal fields (2D; 2D integer; 3D with nl-1 or nl levels, one, two, or three values)
  integer, allocatable       :: s_mpitype_nod2D(:),     r_mpitype_nod2D(:) 
  integer, allocatable       :: s_mpitype_nod2D_i(:),   r_mpitype_nod2D_i(:)
  integer, allocatable       :: s_mpitype_nod3D(:,:,:), r_mpitype_nod3D(:,:,:) 

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

   logical :: elem_full_flag  
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

  use iso_c_binding, only: idx_t=>C_INT32_T
  use o_MESH
  implicit none

  integer   n, j, k, nini, nend, ierr

  interface 
     subroutine partit(n,ptr,adj,wgt,np,part) bind(C)
       use iso_c_binding, only: idx_t=>C_INT32_T
       integer(idx_t), intent(in)  :: n, ptr(*), adj(*), wgt(*),np
       integer(idx_t), intent(out) :: part(*)
     end subroutine partit
  end interface

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

     if (mype==0) then
        write(*,*) 'Calling partit'
        call partit(ssh_stiff%dim, ssh_stiff%rowptr, ssh_stiff%colind, &
             nlevels_nod2D, npes, part)
  	
        call check_partitioning

        write(*,*) 'partitioning is done.'

! The stiffness matrix is no longer needed. 
        deallocate(ssh_stiff%rowptr)
        deallocate(ssh_stiff%colind)
        
     end if
     !NR No longer needed - last use was as weight for partitioning
     deallocate(nlevels_nod2D)

     call MPI_BCAST(part,nod2D,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  endif
  
  call communication_edgen
  call communication_nodn
  call communication_elemn
  deallocate(elem_neighbors,elem_edges)
  call communication_elem_fulln
  call mymesh
  if(mype==0) write(*,*) 'Communication arrays have been set up'   
end subroutine set_par_support_ini

!=======================================================================

subroutine check_partitioning

! In general, METIS 5 has several advantages compared to METIS 4, e.g.,
!   * neighbouring tasks get neighbouring partitions (important for multicore computers!)
!   * lower maximum of weights per partition (better load balancing)
!   * lower memory demand
!
! BUT: there might be outliers, single nodes connected to their partition by
!      only one edge or even completely isolated. This spoils everything :-(
!
! This routine checks for isolated nodes and moves them to an adjacent partition,
! trying not to spoil the load balance.

  use o_MESH
  integer :: i, j, k, n, node, n_iso, n_iter, is, ie, kmax, np
  integer :: nod_per_partition(2,0:npes-1)
  integer :: max_nod_per_part(2), min_nod_per_part(2)
  integer :: average_nod_per_part(2)
  logical :: already_counted, found_part
  
  integer :: max_adjacent
  integer, allocatable :: ne_part(:), ne_part_num(:), ne_part_load(:,:)

! call partit(ssh_stiff%dim, ssh_stiff%rowptr, ssh_stiff%colind, &
!             nlevels_nod2D, npes, part)

! Check load balancing
  nod_per_partition(:,:)=0
  do n=1,nod2D
     nod_per_partition(1,part(n)) = nod_per_partition(1,part(n)) +1
     nod_per_partition(2,part(n)) = nod_per_partition(2,part(n)) +nlevels_nod2D(n)
  enddo

  min_nod_per_part(1) = minval( nod_per_partition(1,:))
  min_nod_per_part(2) = minval( nod_per_partition(2,:))

  max_nod_per_part(1) = maxval( nod_per_partition(1,:))
  max_nod_per_part(2) = maxval( nod_per_partition(2,:))

  average_nod_per_part(1) = nod2D / npes
  average_nod_per_part(2) = sum(nlevels_nod2D(:)) / npes
  
! Now check for isolated nodes (connect by one or even no edge to other
! nodes of its partition) and repair, if possible

  max_adjacent = maxval(ssh_stiff%rowptr(2:nod2D+1) - ssh_stiff%rowptr(1:nod2D))
  allocate(ne_part(max_adjacent), ne_part_num(max_adjacent), ne_part_load(2,max_adjacent))

  isolated_nodes_check: do n_iter = 1, 10
     print *,' '
     print *,'Check for isolated nodes, new iteration ========'
     n_iso = 0
     do n=1,nod2D
        is = ssh_stiff%rowptr(n)
        ie = ssh_stiff%rowptr(n+1) -1

        if (count(part(ssh_stiff%colind(is:ie)) == part(n)) <= 1) then

           n_iso = n_iso+1
           print *,'Isolated node',n, 'in partition', part(n)
           print *,'Neighbouring nodes are in partitions', part(ssh_stiff%colind(is:ie))

        ! count the adjacent nodes of the other PEs
        
           np=1
           ne_part(1) = part(ssh_stiff%colind(is))
           ne_part_num(1) = 1
           ne_part_load(1,1) = nod_per_partition(1,ne_part(1)) + 1
           ne_part_load(2,1) = nod_per_partition(2,ne_part(1)) + nlevels_nod2D(n)
           
           do i=is+1,ie
              node = ssh_stiff%colind(i)
              if (part(node)==part(n)) cycle
              already_counted = .false.
              do k=1,np
                 if (part(node) == ne_part(k)) then
                    ne_part_num(k) = ne_part_num(k) + 1
                    already_counted = .true.
                    exit
                 endif
              enddo
              if (.not. already_counted) then
                 np = np+1
                 ne_part(np) = part(node)
                 ne_part_num(np) = 1
                 ne_part_load(1,np) = nod_per_partition(1,ne_part(np)) + 1
                 ne_part_load(2,np) = nod_per_partition(2,ne_part(np)) + nlevels_nod2D(n)
              endif
           enddo
        
        ! Now, check for two things: The load balance, and if 
        ! there is more than one node of that partition.
        ! Otherwise, it would become isolated again.

        ! Best choice would be the partition with most adjacent nodes (edgecut!)
        ! Choose, if it does not decrease the load balance. 
        !        (There might be two partitions with the same number of adjacent
        !         nodes. Don't care about this here)

           kmax = maxloc(ne_part_num(1:np),1)

           if (ne_part_num(kmax) <= 1) then
              print *,'Sorry, no chance to solve an isolated node problem'
              exit isolated_nodes_check
           endif

           if  (ne_part_load(1,kmax) <= max_nod_per_part(1) .and. &
                ne_part_load(2,kmax) <= max_nod_per_part(2) ) then
              k = kmax
           else
           ! Don't make it too compicated. Reject partitions that have only one
           ! adjacent node. Take the next not violating the load balance.
              found_part = .false.
              do k=1,np
                 if (ne_part_num(k)==1 .or. k==kmax) cycle
                 
                 if  (ne_part_load(1,k) <= max_nod_per_part(1) .and. &
                      ne_part_load(2,k) <= max_nod_per_part(2) ) then
                    
                    found_part = .true.
                    exit
                 endif
              enddo

              if (.not. found_part) then
              ! Ok, don't think to much. Simply go for minimized edge cut.
                 k = kmax
              endif
           endif
        
! Adjust the load balancing
        
           nod_per_partition(1,ne_part(k)) = nod_per_partition(1,ne_part(k)) + 1 
           nod_per_partition(2,ne_part(k)) = nod_per_partition(2,ne_part(k)) + nlevels_nod2D(n)
           nod_per_partition(1,part(n))    = nod_per_partition(1,part(n)) - 1
           nod_per_partition(2,part(n))    = nod_per_partition(2,part(n)) - nlevels_nod2D(n)

! And, finally, move nod n to other partition        
           part(n) = ne_part(k)
           print *,'Node',n,'is moved to part',part(n)
        endif
  enddo  
  
  if (n_iso==0) then
     print *,'No isolated nodes found'
     exit isolated_nodes_check 
  endif
  ! Check for isolated nodes again
end do isolated_nodes_check

  deallocate(ne_part, ne_part_num, ne_part_load)

  if (n_iso > 0) then
     print *,'++++++++++++++++++++++++++++++++++++++++++++'
     print *,'+  WARNING: PARTITIONING LOOKS REALLY BAD  +'
     print *,'+  It was not possible to remove all       +'
     print *,'+  isolated nodes. Consider to rerun with  +'
     print *,'+  new metis random seed (see Makefile.in) +'
     print *,'++++++++++++++++++++++++++++++++++++++++++++'
  endif

  print *,'=== LOAD BALANCING ==='
  print *,'2D nodes: min, aver, max per part',min_nod_per_part(1), &
          average_nod_per_part(1),max_nod_per_part(1)

  write(*,"('2D nodes: percent min, aver, max ',f8.3,'%, 100%, ',f8.3,'%')") &
          100.*real(min_nod_per_part(1)) / real(average_nod_per_part(1)), &
          100.*real(max_nod_per_part(1)) / real(average_nod_per_part(1))

  print *,'3D nodes: Min, aver, max per part',min_nod_per_part(2), &
          average_nod_per_part(2),max_nod_per_part(2)
  write(*,"('3D nodes: percent min, aver, max ',f8.3,'%, 100%, ',f8.3,'%')") &
          100.*real(min_nod_per_part(2)) / real(average_nod_per_part(2)), &
          100.*real(max_nod_per_part(2)) / real(average_nod_per_part(2))




end subroutine check_partitioning
!=======================================================================

subroutine set_par_support
  use o_MESH
  implicit none

  integer   n, offset
  integer :: i, max_nb, nb, nini, nend, nl1, n_val
  integer, allocatable :: blocklen(:),     displace(:)

  !
  ! In the distributed memory version, most of the job is already done 
  ! at the initialization phase and is taken into account in read_mesh
  ! routine. Here, MPI datatypes are built and buffers for MPI wait requests
  ! are allocated. 

   if (npes > 1) then

!================================================
! MPI REQUEST BUFFERS
!================================================
      allocate(com_edge2D%req(          com_edge2D%rPEnum +      com_edge2D%sPEnum))
      allocate(com_nod2D%req(            com_nod2D%rPEnum +       com_nod2D%sPEnum))
      allocate(com_elem2D%req(          com_elem2D%rPEnum +      com_elem2D%sPEnum))
      allocate(com_elem2D_full%req(com_elem2D_full%rPEnum + com_elem2D_full%sPEnum))

!================================================
! MPI DATATYPES
!================================================
      ! Build MPI Data types for halo exchange: Edges
      allocate(r_mpitype_edge2D(com_edge2D%rPEnum))  ! 2D
      allocate(s_mpitype_edge2D(com_edge2D%sPEnum))  

      ! Upper limit for the length of the local interface between the neighbor PEs 
      max_nb = max(maxval(com_edge2D%rptr(2:com_edge2D%rPEnum+1) - com_edge2D%rptr(1:com_edge2D%rPEnum)), &
                   maxval(com_edge2D%sptr(2:com_edge2D%sPEnum+1) - com_edge2D%sptr(1:com_edge2D%sPEnum)))

      allocate(displace(max_nb), blocklen(max_nb))

      do n=1,com_edge2D%rPEnum
         nb = 1
         nini = com_edge2D%rptr(n)
         nend = com_edge2D%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_edge2D%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_edge2D%rlist(i) /= com_edge2D%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_edge2D%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_DOUBLE_PRECISION, r_mpitype_edge2D(n), MPIerr)

         call MPI_TYPE_COMMIT(r_mpitype_edge2D(n),   MPIerr) 
      enddo

      do n=1,com_edge2D%sPEnum
         nb = 1
         nini = com_edge2D%sptr(n)
         nend = com_edge2D%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_edge2D%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_edge2D%slist(i) /= com_edge2D%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_edge2D%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_DOUBLE_PRECISION, s_mpitype_edge2D(n), MPIerr)

         call MPI_TYPE_COMMIT(s_mpitype_edge2D(n),   MPIerr) 

      enddo

      deallocate(displace, blocklen)


      ! Build MPI Data types for halo exchange: Elements
      allocate(r_mpitype_elem2D(com_elem2D%rPEnum,4))     ! 2D, small halo
      allocate(s_mpitype_elem2D(com_elem2D%sPEnum,4))
      allocate(r_mpitype_elem2D_full_i(com_elem2D_full%rPEnum))   ! 2D, wide halo, integer
      allocate(s_mpitype_elem2D_full_i(com_elem2D_full%sPEnum))

      allocate(r_mpitype_elem2D_full(com_elem2D_full%rPEnum,4))     ! 2D, wide halo 
      allocate(s_mpitype_elem2D_full(com_elem2D_full%sPEnum,4))

      allocate(r_mpitype_elem3D(com_elem2D%rPEnum, nl-1:nl,4))     ! 3D, small halo 
      allocate(s_mpitype_elem3D(com_elem2D%sPEnum, nl-1:nl,4))

      allocate(r_mpitype_elem3D_full(com_elem2D_full%rPEnum, nl-1:nl,4))     ! 3D, wide halo
      allocate(s_mpitype_elem3D_full(com_elem2D_full%sPEnum, nl-1:nl,4))
      
      
      ! Upper limit for the length of the local interface between the neighbor PEs 
      max_nb = max(  &
           maxval(com_elem2D%rptr(2:com_elem2D%rPEnum+1) - com_elem2D%rptr(1:com_elem2D%rPEnum)), &
           maxval(com_elem2D%sptr(2:com_elem2D%sPEnum+1) - com_elem2D%sptr(1:com_elem2D%sPEnum)), &
           maxval(com_elem2D_full%rptr(2:com_elem2D_full%rPEnum+1) - com_elem2D_full%rptr(1:com_elem2D_full%rPEnum)), &
           maxval(com_elem2D_full%sptr(2:com_elem2D_full%sPEnum+1) - com_elem2D_full%sptr(1:com_elem2D_full%sPEnum)))
      
      allocate(displace(max_nb), blocklen(max_nb))

      
      do n=1,com_elem2D%rPEnum
         nb = 1
         nini = com_elem2D%rptr(n)
         nend = com_elem2D%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_elem2D%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_elem2D%rlist(i) /= com_elem2D%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_elem2D%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         DO n_val=1,4
            call MPI_TYPE_INDEXED(nb, blocklen*n_val, displace*n_val, MPI_DOUBLE_PRECISION, &
                 r_mpitype_elem2D(n,n_val), MPIerr)

            call MPI_TYPE_COMMIT(r_mpitype_elem2D(n,n_val), MPIerr) 

            DO nl1=nl-1, nl
               call MPI_TYPE_INDEXED(nb, blocklen*nl1*n_val, displace*nl1*n_val, MPI_DOUBLE_PRECISION, & 
                    r_mpitype_elem3D(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(r_mpitype_elem3D(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO  
      enddo

      do n=1,com_elem2D%sPEnum
         nb = 1
         nini = com_elem2D%sptr(n)
         nend = com_elem2D%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_elem2D%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_elem2D%slist(i) /= com_elem2D%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_elem2D%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
                  
         DO n_val=1,4
            call MPI_TYPE_INDEXED(nb, blocklen*n_val, displace*n_val, MPI_DOUBLE_PRECISION, &
                 s_mpitype_elem2D(n, n_val), MPIerr)

            call MPI_TYPE_COMMIT(s_mpitype_elem2D(n, n_val),   MPIerr) 
 
            DO nl1=nl-1, nl
               call MPI_TYPE_INDEXED(nb, blocklen*nl1*n_val, displace*nl1*n_val, MPI_DOUBLE_PRECISION, & 
                    s_mpitype_elem3D(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(s_mpitype_elem3D(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO  
      enddo
      
      do n=1,com_elem2D_full%rPEnum
         nb = 1
         nini = com_elem2D_full%rptr(n)
         nend = com_elem2D_full%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_elem2D_full%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_elem2D_full%rlist(i) /= com_elem2D_full%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_elem2D_full%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         call MPI_TYPE_INDEXED(nb, blocklen,displace,MPI_INTEGER, r_mpitype_elem2D_full_i(n),MPIerr)

         call MPI_TYPE_COMMIT(r_mpitype_elem2D_full_i(n), MPIerr)

         DO n_val=1,4

            call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_DOUBLE_PRECISION, &
                 r_mpitype_elem2D_full(n,n_val), MPIerr)
            call MPI_TYPE_COMMIT(r_mpitype_elem2D_full(n, n_val),   MPIerr)

            DO nl1=nl-1, nl
               call MPI_TYPE_INDEXED(nb, blocklen*nl1*n_val, displace*nl1*n_val, MPI_DOUBLE_PRECISION, & 
                    r_mpitype_elem3D_full(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(r_mpitype_elem3D_full(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO   
      enddo

      do n=1,com_elem2D_full%sPEnum
         nb = 1
         nini = com_elem2D_full%sptr(n)
         nend = com_elem2D_full%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_elem2D_full%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_elem2D_full%slist(i) /= com_elem2D_full%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_elem2D_full%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo
         
         call MPI_TYPE_INDEXED(nb, blocklen,displace,MPI_INTEGER, s_mpitype_elem2D_full_i(n), MPIerr)

         call MPI_TYPE_COMMIT(s_mpitype_elem2D_full_i(n), MPIerr)  
 
         DO n_val=1,4
            call MPI_TYPE_INDEXED(nb, blocklen, displace, MPI_DOUBLE_PRECISION, &
                 s_mpitype_elem2D_full(n,n_val),MPIerr)
            call MPI_TYPE_COMMIT(s_mpitype_elem2D_full(n,n_val),   MPIerr)
  
            DO nl1=nl-1, nl
               call MPI_TYPE_INDEXED(nb, blocklen*nl1*n_val, displace*nl1*n_val, MPI_DOUBLE_PRECISION, & 
                    s_mpitype_elem3D_full(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(s_mpitype_elem3D_full(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO
      enddo

      deallocate(displace,blocklen)


   ! Build MPI Data types for halo exchange: Nodes

      allocate(r_mpitype_nod2D(com_nod2D%rPEnum))     ! 2D
      allocate(s_mpitype_nod2D(com_nod2D%sPEnum))
      allocate(r_mpitype_nod2D_i(com_nod2D%rPEnum))   ! 2D integer
      allocate(s_mpitype_nod2D_i(com_nod2D%sPEnum))   

      allocate(r_mpitype_nod3D(com_nod2D%rPEnum,nl-1:nl,3))  ! 3D with nl-1 or nl layers, 1-3 values 
      allocate(s_mpitype_nod3D(com_nod2D%sPEnum,nl-1:nl,3))
  
      ! Upper limit for the length of the local interface between the neighbor PEs 
      max_nb = max(maxval(com_nod2D%rptr(2:com_nod2D%rPEnum+1) - com_nod2D%rptr(1:com_nod2D%rPEnum)), &
                   maxval(com_nod2D%sptr(2:com_nod2D%sPEnum+1) - com_nod2D%sptr(1:com_nod2D%sPEnum)))

      allocate(displace(max_nb), blocklen(max_nb))

      do n=1,com_nod2D%rPEnum
         nb = 1
         nini = com_nod2D%rptr(n)
         nend = com_nod2D%rptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_nod2D%rlist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_nod2D%rlist(i) /= com_nod2D%rlist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_nod2D%rlist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_DOUBLE_PRECISION, & 
              r_mpitype_nod2D(n),     MPIerr)

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_INTEGER, & 
              r_mpitype_nod2D_i(n),   MPIerr)

         call MPI_TYPE_COMMIT(r_mpitype_nod2D(n),     MPIerr)    
         call MPI_TYPE_COMMIT(r_mpitype_nod2D_i(n),   MPIerr)     

         DO nl1=nl-1, nl
            DO n_val=1,3
               call MPI_TYPE_INDEXED(nb, blocklen*nl1*n_val, displace*nl1*n_val, MPI_DOUBLE_PRECISION, & 
                    r_mpitype_nod3D(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(r_mpitype_nod3D(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO
      enddo

      do n=1,com_nod2D%sPEnum
         nb = 1
         nini = com_nod2D%sptr(n)
         nend = com_nod2D%sptr(n+1) - 1
         displace(:) = 0
         displace(1) = com_nod2D%slist(nini) -1  ! C counting, start at 0
         blocklen(:) = 1
         do i=nini+1, nend
            if (com_nod2D%slist(i) /= com_nod2D%slist(i-1) + 1) then  
               ! New block
               nb = nb+1
               displace(nb) = com_nod2D%slist(i) -1
            else
               blocklen(nb) = blocklen(nb)+1
            endif
         enddo

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_DOUBLE_PRECISION, & 
              s_mpitype_nod2D(n),     MPIerr)

         call MPI_TYPE_INDEXED(nb, blocklen,      displace,      MPI_INTEGER, & 
              s_mpitype_nod2D_i(n),   MPIerr)

         call MPI_TYPE_COMMIT(s_mpitype_nod2D(n),     MPIerr)    
         call MPI_TYPE_COMMIT(s_mpitype_nod2D_i(n),   MPIerr)     

         DO nl1=nl-1, nl
            DO n_val=1,3
               call MPI_TYPE_INDEXED(nb, blocklen*nl1*n_val, displace*nl1*n_val, MPI_DOUBLE_PRECISION, & 
                    s_mpitype_nod3D(n,nl1,n_val),  MPIerr)

               call MPI_TYPE_COMMIT(s_mpitype_nod3D(n,nl1,n_val),  MPIerr)  
            ENDDO
         ENDDO
      enddo

      deallocate(blocklen,     displace)

   endif

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
real(kind=WP), allocatable    :: heat_flux_old(:), Tsurf_old(:)  !PS
real(kind=WP), allocatable    :: S_rhs(:,:)
real(kind=WP), allocatable    :: tr_arr(:,:,:),tr_arr_old(:,:,:)
real(kind=WP), allocatable    :: del_ttf(:,:)

real(kind=WP), allocatable    :: water_flux(:), Ssurf(:)
real(kind=WP), allocatable    :: water_flux_old(:), Ssurf_old(:) !PS
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

! Auxiliary arrays for vector-invariant form of momentum advection
real(kind=WP), allocatable,dimension(:,:)   :: vorticity

!Viscosity and diff coefs
real(kind=WP), allocatable,dimension(:,:)   :: Av,Kv
real(kind=WP), allocatable,dimension(:,:,:)   :: Kv2 
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
real(kind=WP), allocatable    :: Av_mean(:,:),Kv_mean(:,:,:) 
real(kind=WP), allocatable    :: hbl_mean(:)
!_______________________________________________________________________________
! Arrays added for ALE implementation:
! --> layer thinkness at node and depthlayer for t=n and t=n+1
real(kind=WP), allocatable,dimension(:,:)   :: hnode, hnode_new, zbar_3d_n, Z_3d_n

! --> layer thinkness at elements, interpolated from hnode
real(kind=WP), allocatable,dimension(:,:)   :: helem

! --> thinkness of bottom elem (important for partial cells)
real(kind=WP), allocatable,dimension(:)     :: bottom_elem_thickness 

! --> The increment of total fluid depth on elements. It is used to update the matrix
real(kind=WP), allocatable,dimension(:)     :: dhe

! --> hbar, hbar_old: correspond to the elevation, but on semi-integer time steps.
real(kind=WP), allocatable,dimension(:)     :: hbar, hbar_old

! --> auxiliary array to store an intermediate part of the rhs computations.
real(kind=WP), allocatable,dimension(:)     :: ssh_rhs_old !, ssh_rhs_old2 !PS

! --> auxiliary array to store depth of layers and depth of mid level due to changing 
!     layer thinkness at every node
real(kind=WP), allocatable,dimension(:)     :: zbar_n, Z_n

! --> multiplication factor for surface boundary condition in 
!     diff_ver_part_impl_ale(tr_num) between linfs -->=0.0 and noninfs 
!     (zlevel,zstar...) --> = 1.0
real(kind=WP)                               :: is_nonlinfs

!_______________________________________________________________________________
!Monin-Obukhov correction
real*8,allocatable :: mo(:,:),mixlength(:)
!GM_stuff
real*8,allocatable :: bvfreq(:,:),mixlay_dep(:),bv_ref(:)

real(kind=WP), allocatable    :: fer_UV(:,:,:), fer_wvel(:,:)
real(kind=WP), target, allocatable    :: fer_c(:), fer_K(:), fer_gamma(:,:,:)
END MODULE o_ARRAYS
!==========================================================
