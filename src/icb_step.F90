module iceberg_step
 USE MOD_MESH
 use MOD_PARTIT
 use MOD_ICE
 USE MOD_DYN
 use iceberg_params
 use iceberg_dynamics
 use iceberg_element

implicit none

  public    ::  iceberg_calculation
  public    ::  iceberg_step1
  public    ::  iceberg_step2
  public    ::  initialize_velo
  public    ::  trajectory
  public    ::  depth_bathy
  public    ::  parallel2coast
  public    ::  projection
  public    ::  iceberg_restart
  public    ::  iceberg_restart_with_icesheet
  public    ::  iceberg_out
  public    ::  init_buoys
  public    ::  init_icebergs
  public    ::  init_icebergs_with_icesheet
  public    ::  determine_save_count
  public    ::  init_buoy_output
  public    ::  write_buoy_props_netcdf

contains 

subroutine iceberg_calculation(ice, mesh, partit, dynamics, istep) 
 !======================================================================!
 !									!
 !                     ICEBERG MODULE FOR FESOM				!
 !									!
 !======================================================================!
 ! last update: 07.10.2015, T. Rackow (AWI)				
 !======================================================================!  

 !==================== MODULES & DECLARATIONS ==========================!=
 									!=
 use o_param 		!for ?						!=
 use g_config		!for istep, step_per_day, logfile_outfreq	!=
!=
 implicit none								!=
									!=
 integer	:: ib, times, istep
 integer	:: istep_end_synced
 integer:: req, status(MPI_STATUS_SIZE)
 logical:: completed
 real(kind=8) 	:: t0, t1, t2, t3, t4, t0_restart, t1_restart   	!=
 logical	:: firstcall=.true. 					!=
 logical	:: lastsubstep  					!=

 real		:: arr_from_block(15)					!=
 integer	:: elem_from_block					!=  
 real		:: vl_from_block(4)					!=	
 real,dimension(15*ib_num):: arr_block_red				!=
 integer,dimension(ib_num):: elem_block_red				!=
 integer,dimension(ib_num):: pe_block_red				!=
 integer    :: n
 real, dimension(4*ib_num):: vl_block_red				!=

type(t_ice), intent(inout), target :: ice
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
type(t_dyn)   , intent(inout), target :: dynamics
!==================== MODULES & DECLARATIONS ==========================!= 
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

! kh 16.03.21 (asynchronous) iceberg computation starts with the content in common arrays at istep and will merge its results at istep_end_synced
 istep_end_synced = istep + steps_per_ib_step - 1
 if(firstcall) then
  !overwrite icb_modules if restart, initialize netcdf output if no restart:
  
  t0_restart=MPI_Wtime()
  if (use_icesheet_coupling) then
    call iceberg_restart_with_icesheet(partit)
  else
    call iceberg_restart(partit)
  end if
  t1_restart=MPI_Wtime()
  firstcall = .false.
  !call init_global_tides
  !call tides_distr
 end if  

 t0=MPI_Wtime()
 
 !call update_global_tides !for each timestep istep once

 !write ib values in 2 larger arrays for
 !faster communication via ALLREDUCE
 arr_block = 0.0
 elem_block = 0
 pe_block = 0
 !for communication of averaged volume losses
 vl_block = 0.0

 !the original routine iceberg_step has been splitted 
 !in the part before (step1) and after (step2) communication
 !
 !this results in two do-loops BUT the expensive routine
 !com_values can be replaced by MPI_ALLREDUCE
 
 
 !============================= STEP 1 =================================!
 
 
 do ib=1, ib_num
  lastsubstep = .false.
  if( real(istep) > real(step_per_day)*calving_day(ib) ) then !iceberg calved

    lastsubstep = .true. !do output every timestep

    if( .not.melted(ib) ) then
        call iceberg_step1(ice, mesh, partit, dynamics, ib, height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),&
                            Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
                            conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
                            u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib), lastsubstep,&
                            f_u_ib_old(ib), f_v_ib_old(ib), l_semiimplicit,   &
                            semiimplicit_coeff, AB_coeff, istep)			
    end if 
  end if
 end do

 t1=MPI_Wtime()
 
 
 !========================== COMMUNICATION =============================!

 !all PEs need the array arr(15) and the iceberg element
 !in step2
 !ALLREDUCE: arr_block, elem_block
 
 arr_block_red = 0.0
 elem_block_red= 0
 pe_block_red= 0
 vl_block_red = 0.0

!$omp critical 
 call MPI_IAllREDUCE(arr_block, arr_block_red, 15*ib_num, MPI_DOUBLE_PRECISION, MPI_SUM, partit%MPI_COMM_FESOM_IB, req, partit%MPIERR_IB)
!$omp end critical

 completed = .false.
 do while (.not. completed)
!$omp critical
CALL MPI_TEST(req, completed, status, partit%MPIERR_IB)
!$omp end critical
 end do

!$omp critical 
 call MPI_IAllREDUCE(elem_block, elem_block_red, ib_num, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM_IB, req, partit%MPIERR_IB)  
!$omp end critical

completed = .false.
 do while (.not. completed)
!$omp critical
  CALL MPI_TEST(req, completed, status, partit%MPIERR_IB)
!$omp end critical
 end do

!$omp critical
 call MPI_IAllREDUCE(pe_block, pe_block_red, ib_num, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM_IB, req, partit%MPIERR_IB)
!$omp end critical

completed = .false.
 do while (.not. completed)
!$omp critical
  CALL MPI_TEST(req, completed, status, partit%MPIERR_IB)
!$omp end critical
 end do

!!$omp critical
! call MPI_IAllREDUCE(elem_area_block, elem_area_block_red, ib_num, MPI_REAL, MPI_SUM, partit%MPI_COMM_FESOM_IB, req, partit%MPIERR_IB)
!!$omp end critical
!
!completed = .false.
! do while (.not. completed)
!!$omp critical
!  CALL MPI_TEST(req, completed, status, partit%MPIERR_IB)
!!$omp end critical
! end do


!$omp critical 
 call MPI_IAllREDUCE(vl_block, vl_block_red, 4*ib_num, MPI_DOUBLE_PRECISION, MPI_SUM, partit%MPI_COMM_FESOM_IB, req, partit%MPIERR_IB)
!$omp end critical

 completed = .false.
 do while (.not. completed)
!$omp critical
  CALL MPI_TEST(req, completed, status, partit%MPIERR_IB)
!$omp end critical
 end do

 buoy_props=0.

 t2=MPI_Wtime()   
  
 
 !============================= STEP 2 =================================!

 do ib=1, ib_num

  !get the smaller array arr(15) and 
  !the iceberg element for iceberg ib
  !as before
  arr_from_block = arr_block_red( (ib-1)*15+1 : ib*15)
  elem_from_block= elem_block_red(ib)
  !averaged volume losses
  vl_from_block = vl_block_red( (ib-1)*4+1 : ib*4)
  bvl_mean(ib)=vl_from_block(1)
  lvlv_mean(ib)=vl_from_block(2)
  lvle_mean(ib)=vl_from_block(3)
  lvlb_mean(ib)=vl_from_block(4)

  lastsubstep = .false.
  if( real(istep) > real(step_per_day)*calving_day(ib) ) then !iceberg calved
  
    !substeps don't work anymore with new communication
    lastsubstep = .true. !do output every timestep
   
    if( .not.melted(ib) ) then
        call iceberg_step2(mesh, partit, arr_from_block, elem_from_block, ib, height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),&
                            Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
                            conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
                            u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib), lastsubstep,&
                            f_u_ib_old(ib), f_v_ib_old(ib), l_semiimplicit,   &
                            semiimplicit_coeff, AB_coeff, istep, elem_block_red, &
                            pe_block_red)
    end if
  end if
end do

 t3=MPI_Wtime() 


 !========================== VECTOR OUTPUT =============================!
 
 !call iceberg_vector_ncout !look in routines what is really written out!
 !istep, lon_deg_geo, lat_deg_geo, u_ib_geo, v_ib_geo, volume losses (set to zero)
 !introduce force_last_output(ib)?

 if (mod(istep_end_synced,icb_outfreq)==0 .AND. .not.ascii_out) then

   if (mype==0 .AND. (real(istep) > real(step_per_day)*calving_day(1) ) ) call write_buoy_props_netcdf(partit)
       
   ! all PEs: set back to zero for next round
   bvl_mean=0.0
   lvlv_mean=0.0
   lvle_mean=0.0      
   lvlb_mean=0.0

 end if

 t4=MPI_Wtime() 
 
 if (mod(istep,logfile_outfreq)==0 .and. mype==0) then 
   write(*,*) 'icebergs took', t4-t0 
   write(*,*) 'iceberg step1 took', t1-t0
   write(*,*) 'NEW comvalues took', t2-t1
   write(*,*) 'iceberg step2 took', t3-t2
   write(*,*) 'vector output took', t4-t3
   write(*,*) 'reading restart took', t1_restart-t0_restart
   write(*,*) '*************************************************************'
 end if

end subroutine iceberg_calculation


!****************************************************************************************************************************
!****************************************************************************************************************************


subroutine iceberg_step1(ice, mesh, partit, dynamics, ib, height_ib_single,length_ib_single,width_ib_single, lon_deg,lat_deg, &
			Co,Ca,Ci, Cdo_skin,Cda_skin, rho_icb, 		   &
			conc_sill,P_sill, rho_h2o,rho_air,rho_ice,	   & 
			u_ib,v_ib, iceberg_elem, find_iceberg_elem, 	   &
			lastsubstep, f_u_ib_old,	   &
			f_v_ib_old, l_semiimplicit, semiimplicit_coeff,    &
			AB_coeff, istep)
			
 !============================= MODULES & DECLARATIONS =========================================!=
 												!=
 use o_param 		!for rad								!=
 use g_rotate_grid	!for subroutine g2r, logfile_outfreq					!=
 use g_config, only: steps_per_ib_step
 !=
use iceberg_params, only: length_ib, width_ib, scaling, elem_block, elem_area_glob !smallestvol_icb, arr_block, elem_block, l_geo_out, icb_outfreq, l_allowgrounding, draft_scale, reject_elem, melted, grounded, scaling !, length_ib, width_ib, scaling
!#else
! use iceberg_params, only: smallestvol_icb, arr_block, elem_block, l_geo_out, icb_outfreq, l_allowgrounding, draft_scale, melted, grounded, scaling !, length_ib, width_ib, scaling
!#endif
 												!=
 implicit none											!=
 
 logical                :: reject_tmp 
 integer, intent(in)	:: ib, istep
 real,    intent(inout)	:: height_ib_single,length_ib_single,width_ib_single
 real,    intent(inout)	:: lon_deg,lat_deg
 real, 	  intent(in)	:: Co,Ca,Ci, Cdo_skin,Cda_skin
 real,	  intent(in)	:: rho_icb, conc_sill,P_sill, rho_h2o,rho_air,rho_ice
 real,	  intent(inout)	:: u_ib,v_ib
 integer, intent(inout)	:: iceberg_elem !global
 logical, intent(inout)	:: find_iceberg_elem
 logical, intent(in)	:: lastsubstep 
 real,    intent(inout)	:: f_u_ib_old, f_v_ib_old
 logical, intent(in)	:: l_semiimplicit
 real,    intent(in)	:: semiimplicit_coeff
 real,    intent(in)	:: AB_coeff
 real, dimension(:), pointer     :: coriolis
 integer :: istep_end_synced
 
 integer, dimension(:), save, allocatable :: local_idx_of
 real      			:: depth_ib, volume_ib, mass_ib
 real				:: lon_rad, lat_rad, new_u_ib, new_v_ib
 real	   			:: old_lon,old_lat, frozen_in, P_ib, conci_ib
 real				:: lon_rad_out, lat_rad_out  !for unrotated output
 real				:: lon_deg_out, lat_deg_out  !for unrotated output
 integer   			:: i, iceberg_node  
 real 				:: dudt, dvdt
!! LA: add threshold for number of icebergs in one elemt
 integer                        :: num_ib_in_elem, idx
 real                           :: area_ib_tot

 !iceberg output 
 character 			:: ib_char*10
 character 			:: file_track*80
 character 			:: file_forces_u*80
 character 			:: file_forces_v*80
 character 			:: file_meltrates*80	
 logical			:: l_output
 							
 !MPI											
 real	   			:: arr(15)						!=
 logical   			:: i_have_element					!=
 real	   			:: left_mype						!=
 integer   			:: old_element						!=
 real(kind=8) 			:: t0, t1, t2, t3, t4, t5, t6, t7, t8                   !=
 											!=
 !for restart										!=
 logical, save   		:: firstcall=.true.					!=
 !for grounding										!=
 real, dimension(3)		:: Zdepth3						!=
 real				:: Zdepth						!=
 											!=
 real, dimension(2)             :: coords_tmp
! integer, pointer  :: mype
type(t_ice),  intent(inout), target :: ice
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
type(t_dyn)   , intent(inout), target :: dynamics
!========================= MODULES & DECLARATIONS =====================================!=
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"


 mype          =>partit%mype

 istep_end_synced = istep + steps_per_ib_step - 1

 depth_ib = -height_ib_single * rho_icb/rho_h2o
 volume_ib= length_ib_single * width_ib_single * height_ib_single
 mass_ib = volume_ib * rho_icb	 !less mass 
 lon_rad = lon_deg*rad
 lat_rad = lat_deg*rad
 
 if(volume_ib .le. smallestvol_icb) then
  melted(ib) = .true.

  if (mod(istep_end_synced,logfile_outfreq)==0 .and. mype==0 .and. lastsubstep) then
   write(*,*) 'iceberg ', ib,' melted'
  end if

  return
 end if 
 
 if (firstcall) then
  if(mype==0) write(*,*) 'Preparing local_idx_of array...'
  allocate(local_idx_of(elem2D))
  !creates mapping
  call global2local(mesh, partit, local_idx_of, elem2D)
  firstcall=.false.
  if(mype==0) write(*,*) 'Preparing local_idx_of done.' 
 end if 
 
 if (find_iceberg_elem) then
  lon_rad = lon_deg*rad
  lat_rad = lat_deg*rad
  call g2r(lon_rad, lat_rad, lon_rad, lat_rad)
  lat_deg=lat_rad/rad !rotated lat in degree
  lon_deg=lon_rad/rad !rotated lon in degree   
  
  !find LOCAL element where the iceberg starts:
  coords_tmp = [lon_deg, lat_deg]
  call point_in_triangle(mesh, partit, iceberg_elem, coords_tmp)
  !call point_in_triangle(mesh, iceberg_elem, (/lon_deg, lat_deg/))
  i_have_element= (iceberg_elem .ne. 0) !up to 3 PEs possible
  
  if(i_have_element) then
   i_have_element= mesh%elem2D_nodes(1,iceberg_elem) <= partit%myDim_nod2D !1 PE still .true.
   
   
   
   if (use_cavity) then
      reject_tmp = all( (mesh%cavity_depth(mesh%elem2D_nodes(:,iceberg_elem))/=0.0) .OR. (mesh%bc_index_nod2D(mesh%elem2D_nodes(:,iceberg_elem))==0.0) )
      if(reject_tmp) then
!       write(*,*) " * set IB elem ",iceberg_elem,"to zero for IB=",ib
!       write(*,*) " cavity: ",all((mesh%cavity_depth(mesh%elem2D_nodes(:,iceberg_elem))/=0.0))
!       write(*,*) " boundary: ", all(mesh%bc_index_nod2D(mesh%elem2D_nodes(:,iceberg_elem))==1)
       iceberg_elem=0 !reject element
       i_have_element=.false.
      else 
       iceberg_elem=partit%myList_elem2D(iceberg_elem) !global now
      end if
   else
      iceberg_elem=partit%myList_elem2D(iceberg_elem) !global now
   endif   
  end if
  call com_integer(partit, i_have_element,iceberg_elem)
 
  if(iceberg_elem .EQ. 0) then
        write(*,*) 'IB ',ib,' rot. coords:', lon_deg, lat_deg !,lon_rad, lat_rad
   	call par_ex
   	stop 'ICEBERG OUTSIDE MODEL DOMAIN OR IN ICE SHELF REGION'
  end if
  
  ! initialize the iceberg velocity
  if(local_idx_of(iceberg_elem) <= partit%myDim_elem2D ) then
    call initialize_velo(mesh, partit, dynamics, i_have_element, ib, u_ib, v_ib, lon_rad, lat_rad, depth_ib, local_idx_of(iceberg_elem))
  else
    write(*,*) " * skip initialize_velo"
  end if
  !iceberg elem of ib is found
  find_iceberg_elem = .false.
 
 coriolis => mesh%coriolis(:)
  !for AB method

! kh 06.08.21 observed via -check bounds: forrtl: severe (408): fort: (3): Subscript #1 of the array CORIOLIS has value 0 which is less than the lower bound of 1
  if(local_idx_of(iceberg_elem) > 0) then
    if(local_idx_of(iceberg_elem) <= partit%myDim_elem2D ) then

      f_u_ib_old = coriolis(local_idx_of(iceberg_elem))*u_ib
      f_v_ib_old = coriolis(local_idx_of(iceberg_elem))*v_ib
      
    endif
  endif
 end if
 
 
 ! ================== START ICEBERG CALCULATION ====================
 
 arr=0.
 frozen_in = 0.
 i_have_element=.false.
 !if the first node belongs to this processor.. (just one processor enters here!)
 !if( local_idx_of(iceberg_elem) > 0 .and. elem2D_nodes(1,local_idx_of(iceberg_elem)) <= myDim_nod2D ) then
if( local_idx_of(iceberg_elem) > 0 ) then 

  if( elem2D_nodes(1,local_idx_of(iceberg_elem)) <= partit%myDim_nod2D ) then

  i_have_element=.true. 

! kh 16.03.21 (asynchronous) iceberg calculation starts with the content in common arrays at istep and will merge its results at istep_end_synced
! l_output =  lastsubstep .and. mod(istep,icb_outfreq)==0
  l_output =  lastsubstep .and. mod(istep_end_synced,icb_outfreq)==0
  
  !===========================DYNAMICS===============================
  

  call iceberg_dyn(mesh, partit, ice, dynamics, ib, new_u_ib, new_v_ib, u_ib, v_ib, lon_rad,lat_rad, depth_ib, &
                   height_ib_single, length_ib_single, width_ib_single, local_idx_of(iceberg_elem), &
  		   mass_ib, Ci, Ca, Co, Cda_skin, Cdo_skin, &
  		   rho_ice, rho_air, rho_h2o, P_sill,conc_sill, frozen_in, &
  		   file_forces_u, file_forces_v, P_ib, conci_ib, &
		   dt*REAL(steps_per_ib_step), l_output, f_u_ib_old, &
		   f_v_ib_old, l_semiimplicit, semiimplicit_coeff, &
		   AB_coeff, file_meltrates, rho_icb)

  dudt = (new_u_ib-u_ib)/REAL(steps_per_ib_step) / dt
  dvdt = (new_v_ib-v_ib)/REAL(steps_per_ib_step) / dt
		   
  !=======================END OF DYNAMICS============================
 
  call depth_bathy(mesh,partit, Zdepth3, local_idx_of(iceberg_elem))
  !interpolate depth to location of iceberg (2 times because FEM_3eval expects a 2 component vector...)
  call FEM_3eval(mesh,partit, Zdepth,Zdepth,lon_rad,lat_rad,Zdepth3,Zdepth3,local_idx_of(iceberg_elem))
  !write(*,*) 'nodal depth in iceberg ', ib,'s element:', Zdepth3
  !write(*,*) 'depth at iceberg ', ib, 's location:', Zdepth
  
  !=================CHECK IF ICEBERG IS GROUNDED...===================
 old_element = iceberg_elem !save if iceberg left model domain
 if((draft_scale(ib)*abs(depth_ib) .gt. Zdepth) .and. l_allowgrounding ) then 
 !if((draft_scale(ib)*abs(depth_ib) .gt. minval(Zdepth3)) .and. l_allowgrounding ) then 
   !icebergs remains stationary (iceberg can melt above in iceberg_dyn!)
    left_mype = 0.0 
    u_ib = 0.0
    v_ib = 0.0
    old_lon = lon_rad
    old_lat = lat_rad
 
! kh 16.03.21 (asynchronous) iceberg calculation starts with the content in common arrays at istep and will merge its results at istep_end_synced
    grounded(ib) = .true.
    !if (mod(istep_end_synced,logfile_outfreq)==0) then
    if (lverbose_icb) write(*,*) 'iceberg ib ', ib, 'is grounded'
    !end if
 	
 else 
  !===================...ELSE CALCULATE TRAJECTORY====================

 
 t0=MPI_Wtime()
  call trajectory( lon_rad,lat_rad, u_ib,v_ib, new_u_ib,new_v_ib, &
		   lon_deg,lat_deg,old_lon,old_lat, dt*REAL(steps_per_ib_step))
  	   
 t1=MPI_Wtime()
  iceberg_elem=local_idx_of(iceberg_elem)  	!local
  
 t2=MPI_Wtime()
  call find_new_iceberg_elem(mesh, partit, iceberg_elem, [lon_deg, lat_deg], left_mype)
 t3=MPI_Wtime()
  iceberg_elem=partit%myList_elem2D(iceberg_elem)  	!global
  
  if(left_mype > 0.) then
   lon_rad = old_lon
   lat_rad = old_lat
 t4=MPI_Wtime()
   call parallel2coast(mesh,partit, new_u_ib, new_v_ib, lon_rad,lat_rad, local_idx_of(iceberg_elem))
 t5=MPI_Wtime()
   call trajectory( lon_rad,lat_rad, new_u_ib,new_v_ib, new_u_ib,new_v_ib, &
		   lon_deg,lat_deg,old_lon,old_lat, dt*REAL(steps_per_ib_step))
 t6=MPI_Wtime()
   u_ib = new_u_ib
   v_ib = new_v_ib
		   
   iceberg_elem=local_idx_of(iceberg_elem)  	!local
 t7=MPI_Wtime()
   call find_new_iceberg_elem(mesh,partit, iceberg_elem, (/lon_deg, lat_deg/), left_mype)

 t8=MPI_Wtime()
   iceberg_elem=partit%myList_elem2D(iceberg_elem)  	!global
  end if		   
  !================END OF TRAJECTORY CALCULATION=====================
 end if ! iceberg stationary?

  !-----------------------------
  ! LA 2022-11-30
  if(cell_saturation > 0) then
    if( lverbose_icb) then
     write(*,*) " * checking for cell saturation - ib: ", ib, ", old_elem: ", old_element, ", new_elem: ", iceberg_elem
    end if
    select case(cell_saturation) !num of coastal points
     case(1) 
      area_ib_tot = 0.0
     case(2) 
      area_ib_tot = length_ib_single*width_ib_single*scaling(ib)
    end select
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(idx, area_ib_tot)
!$OMP DO
    do idx = 1, size(elem_block)
        if (elem_block(idx) == iceberg_elem) then
            area_ib_tot = area_ib_tot + length_ib(idx) * width_ib(idx) * scaling(idx)
        end if
    end do
!$OMP END DO
!$OMP END PARALLEL
  !-----------------------------

    if ((area_ib_tot > elem_area_glob(iceberg_elem)) .and. (old_element.ne.0) .and. (iceberg_elem.ne.old_element)) then ! .and. (left_mype == 0)) then 
        if( lverbose_icb) then
            write(*,*) " *******************************************************"
            write(*,*) " * set iceberg ", ib, " back to elem ", old_element, " from elem ", iceberg_elem
            write(*,*) " * area_ib_tot = ", area_ib_tot, "; elem_area = ", elem_area(local_idx_of(iceberg_elem))
        
        end if
        i_have_element = .true.
        left_mype = 0.0
        lon_rad = old_lon
        lat_rad = old_lat
        lon_deg = lon_rad/rad
        lat_deg = lat_rad/rad
        iceberg_elem = old_element
        u_ib    = 0.
        v_ib    = 0.
    end if
  end if
  !###########################################
 
  !values for communication
  arr= (/ height_ib_single,length_ib_single,width_ib_single, u_ib,v_ib, lon_rad,lat_rad, &
          left_mype, old_lon,old_lat, frozen_in, dudt, dvdt, P_ib, conci_ib/) 

  !save in larger array	  
  arr_block((ib-1)*15+1 : ib*15)=arr
  elem_block(ib)=iceberg_elem
  pe_block(ib)=mype

  volume_ib=height_ib_single*length_ib_single*width_ib_single
  call prepare_icb2fesom(mesh,partit,ib,i_have_element,local_idx_of(iceberg_elem),depth_ib,height_ib_single)
 end if !processor has element?
end if !... and first node belongs to processor?

 !t1=MPI_Wtime()
 !if (mod(istep,logfile_outfreq)==0 .and. i_have_element .and. lastsubstep) write(*,*) 'dynamics  took', t1-t0
 !if (mod(istep,logfile_outfreq)==0 .and. i_have_element .and. lastsubstep) then
 !   write(*,*) 'trajectory 1  took', t1-t0
 !   write(*,*) 'find_new_iceberg_elem 1  took', t3-t2
 !   write(*,*) 'parallel2coast  took', t5-t4
 !   write(*,*) 'trajectory 2  took', t6-t5
 !   write(*,*) 'find_new_iceberg_elem 2  took', t8-t7
 !end if
 ! =================== END OF ICEBERG CALCULATION ==================
 
 end subroutine iceberg_step1


 
subroutine iceberg_step2(mesh, partit,arr, elem_from_block, ib, height_ib_single, length_ib_single, width_ib_single, lon_deg,lat_deg, &
			Co,Ca,Ci, Cdo_skin,Cda_skin, rho_icb, 		   &
			conc_sill,P_sill, rho_h2o,rho_air,rho_ice,	   & 
			u_ib,v_ib, iceberg_elem, find_iceberg_elem, 	   &
			lastsubstep, f_u_ib_old,	   &
			f_v_ib_old, l_semiimplicit, semiimplicit_coeff,    &
			AB_coeff, istep, elem_block_red, pe_block_red)
			
 !============================= MODULES & DECLARATIONS =========================================!=
 												!=
 use o_param 		!for rad								!=
 use g_rotate_grid	!for subroutine g2r, logfile_outfreq					!=
 use g_config, only: steps_per_ib_step
 use g_comm_auto
!=
use g_comm
use iceberg_params, only: length_ib, width_ib, scaling !smallestvol_icb, arr_block, elem_block, l_geo_out, icb_outfreq, l_allowgrounding, draft_scale, reject_elem, melted, grounded, scaling !, length_ib, width_ib, scaling
!#else
! use iceberg_params, only: smallestvol_icb, buoy_props, bvl_mean, lvlv_mean, lvle_mean, lvlb_mean, ascii_out, l_geo_out, icb_outfreq, l_allowgrounding, draft_scale, elem_block
!#endif
 												!=
 implicit none											!=
 
 real, 	  intent(in)	:: arr(15)
 integer, intent(in)	:: elem_from_block
 integer, intent(in)	:: ib
 real,    intent(inout)	:: height_ib_single, length_ib_single, width_ib_single
 real,    intent(inout)	:: lon_deg,lat_deg
 real, 	  intent(in)	:: Co,Ca,Ci, Cdo_skin,Cda_skin
 real,	  intent(in)	:: rho_icb, conc_sill,P_sill, rho_h2o,rho_air,rho_ice
 real,	  intent(inout)	:: u_ib,v_ib
 integer, intent(inout)	:: iceberg_elem !global
 logical, intent(inout)	:: find_iceberg_elem
 logical, intent(in)	:: lastsubstep 
 real,    intent(inout)	:: f_u_ib_old, f_v_ib_old
 logical, intent(in)	:: l_semiimplicit
 real,    intent(in)	:: semiimplicit_coeff
 real,    intent(in)	:: AB_coeff						
 integer, intent(in), dimension(ib_num):: elem_block_red				!=
 integer, intent(in), dimension(ib_num):: pe_block_red				!=
 
 integer, dimension(:), save, allocatable :: local_idx_of
 real      			:: depth_ib, volume_ib, mass_ib
 real				:: lon_rad, lat_rad, new_u_ib, new_v_ib
 real	   			:: old_lon,old_lat, frozen_in, P_ib, conci_ib
 real				:: lon_rad_out, lat_rad_out  !for unrotated output
 real				:: lon_deg_out, lat_deg_out  !for unrotated output
 real				:: u_ib_out, v_ib_out  	     !for unrotated output
 real				:: dudt_out, dvdt_out  	     !for unrotated output
 integer   			:: i, iceberg_node, istep 
 real 				:: dudt, dvdt

! kh 16.03.21
 integer :: istep_end_synced
 
! LA: add threshold for number of icebergs in one elemt
 integer status(MPI_STATUS_SIZE)
 integer                        :: num_ib_in_elem, idx
 real                           :: area_ib_tot
 !real(real64), dimension(:), allocatable    :: rbuffer, local_elem_area
 real(real64)                   :: elem_area_tmp

 !iceberg output 
 character 			:: ib_char*10
 character 			:: file_track*80
 character 			:: file_forces_u*80
 character 			:: file_forces_v*80
 character 			:: file_meltrates*80	
 logical			:: l_output
 							
 !MPI											
 logical   			:: i_have_element					!=
 real	   			:: left_mype						!=
 integer   			:: old_element						!=
 real(kind=8) 			:: t0, t1, t2, t3, t4					!=
 											!=
 !for restart										!=
 logical, save   		:: firstcall=.true.					!=
 !for grounding										!=
 real, dimension(3)		:: Zdepth3						!=
 real				:: Zdepth						!=

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
 !========================= MODULES & DECLARATIONS =====================================!=
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
 
 !all PEs enter here with identical array arr
  
 !**** check if iceberg melted in step 1 ****! 
 !mass_ib = arr(1) * arr(2) * arr(3) * rho_icb
 !if(mass_ib .le. 1.0e-6) then
 ! return
 !end if 
 
 !**** check if iceberg melted in step 1 ****!
 !call com_values(i_have_element, arr, iceberg_elem) 

! kh 16.03.21 (asynchronous) iceberg calculation starts with the content in common arrays at istep and will merge its results at istep_end_synced
 istep_end_synced = istep + steps_per_ib_step - 1

 iceberg_elem= elem_from_block !update element as before in com_values
 old_element = elem_from_block !save if iceberg left model domain
 height_ib_single= arr(1)
 length_ib_single= arr(2)
 width_ib_single = arr(3)
 u_ib     = arr(4)
 v_ib     = arr(5)
 lon_rad  = arr(6)
 lat_rad  = arr(7) 
 lon_deg  = lon_rad/rad
 lat_deg  = lat_rad/rad
 left_mype= arr(8)
 old_lon  = arr(9)
 old_lat  = arr(10)
 frozen_in= arr(11)
 dudt	  = arr(12)
 dvdt	  = arr(13)
 P_ib	  = arr(14)
 conci_ib = arr(15) 

 !**** check if iceberg melted in step 1 ****! 
 volume_ib = height_ib_single * length_ib_single * width_ib_single ! * rho_icb
 if(volume_ib .le. smallestvol_icb) then
  buoy_props(ib, :) = 0. ! for output: NaN or MissVal could be written here
  return
 end if 
 !**** check if iceberg melted in step 1 ****!

 t2=MPI_Wtime()

  !!**** LA: check if iceberg changed element and new element is too full 
  !!if(local_idx_of(iceberg_elem) > 0 .and. iceberg_elem .ne. old_element) then !IB left model domain
  !if(iceberg_elem .ne. old_element) then
    if (firstcall) then
      allocate(local_idx_of(elem2D))
      !creates mapping
      call global2local(mesh, partit, local_idx_of, elem2D)
      firstcall=.false.
    end if 

 if(left_mype > 0.) then
   call iceberg_elem4all(mesh, partit, iceberg_elem, lon_deg, lat_deg) !Just PE changed?
   if(iceberg_elem == 0 ) then
           left_mype = 0.0
           lon_rad = old_lon
           lat_rad = old_lat 
           lon_deg = lon_rad/rad
           lat_deg = lat_rad/rad
           iceberg_elem = old_element
           u_ib    = 0.
           v_ib    = 0.
   else if(cell_saturation > 0) then
     if (mype==0) write(*,*) 'iceberg ',ib, ' changed PE or was very fast'
     elem_area_tmp = elem_area_glob(iceberg_elem)
     select case(cell_saturation) !num of coastal points
      case(1) 
       area_ib_tot = 0.0
      case(2) 
       area_ib_tot = length_ib_single*width_ib_single*scaling(ib)
     end select
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(idx, area_ib_tot)
!$OMP DO
     do idx = 1, size(elem_block_red)
         if (elem_block_red(idx) == iceberg_elem) then
     !        write(*,*) " * Found element ",elem_block_red(idx), " for ib ",idx, ", elem_area=",elem_area_tmp
             area_ib_tot = area_ib_tot + length_ib(idx) * width_ib(idx) * scaling(idx)
         end if
     end do
!$OMP END DO
!$OMP END PARALLEL
     if((area_ib_tot > elem_area_tmp) .and. (elem_area_tmp > 0.0) .and. (iceberg_elem.ne.old_element)) then
         if(mype==pe_block_red(ib) .and. lverbose_icb) then
            write(*,*) " *******************************************************"
            write(*,*) " * iceberg changed PE and saturation"
            write(*,*) " * set iceberg ", ib, " back to elem ", old_element, " from elem ", iceberg_elem
            write(*,*) " * area_ib = ", length_ib_single * width_ib_single, "; area_ib_tot = ", area_ib_tot, "; elem_area = ", elem_area_tmp
         end if
         left_mype = 0.0
         lon_rad = old_lon
         lat_rad = old_lat
         lon_deg = lon_rad/rad
         lat_deg = lat_rad/rad
         iceberg_elem = old_element
         u_ib    = 0.
         v_ib    = 0.
     end if
   else 
     if (mype==0) write(*,*) 'iceberg ',ib, ' changed PE or was very fast'
   end if
 end if
 
 t3=MPI_Wtime()
 

 !if(mype==0 .and. lastsubstep .and. mod(istep,icb_outfreq)==0 .and. ascii_out) then
 !if(mype==mod(ib,npes-1) .and. lastsubstep .and. mod(istep,icb_outfreq)==0 .and. ascii_out) then

! kh 16.03.21 (asynchronous) iceberg calculation starts with the content in common arrays at istep and will merge its results at istep_end_synced
! if(mype==0 .and. lastsubstep .and. mod(istep,icb_outfreq)==0) then
 if(mype==0 .and. lastsubstep .and. mod(istep_end_synced,icb_outfreq)==0) then

   !output in 1. unrotated or 2. rotated coordinates      
   u_ib_out = u_ib
   v_ib_out = v_ib
   dudt_out = dudt
   dvdt_out = dvdt
   if(l_geo_out) then
      call r2g(lon_rad_out, lat_rad_out, lon_rad, lat_rad)
      lon_deg_out = lon_rad_out/rad
      lat_deg_out = lat_rad_out/rad
      call vector_r2g(u_ib_out, v_ib_out, lon_rad, lat_rad, 0)
      call vector_r2g(dudt_out, dvdt_out, lon_rad, lat_rad, 0)
   else
      lon_rad_out = lon_rad
      lat_rad_out = lat_rad
      lon_deg_out = lon_deg
      lat_deg_out = lat_deg
   end if

! if(ascii_out) then !use old ASCII output
!  
! file_track='/work/ollie/lackerma/iceberg/iceberg_ICBref_'
! !convert ib integer to string
! write(ib_char,'(I10)') ib
! !left-adjust the string..
! ib_char = adjustl(ib_char)
! !.. and trim while concatenating:
! file_track	=  trim(file_track) // trim(ib_char) // '.dat'
!   
!   open(unit=42,file=file_track,position='append')
!
!! kh 16.03.21 (asynchronous) iceberg calculation starts with the content in common arrays at istep and will merge its results at istep_end_synced
!!  write(42,'(I,12e15.7)') 	istep, lon_rad_out, lat_rad_out, lon_deg_out, lat_deg_out, &
!!  				u_ib_out, v_ib_out, frozen_in, P_sill, P_ib, conci_ib, dudt_out, dvdt_out
!   write(42,'(I,12e15.7)') 	istep_end_synced, lon_rad_out, lat_rad_out, lon_deg_out, lat_deg_out, &
!   				u_ib_out, v_ib_out, frozen_in, P_sill, P_ib, conci_ib, dudt_out, dvdt_out
!   close(42)
!
! else !write in array for faster netcdf output

  buoy_props(ib, 1) = lon_rad_out
  buoy_props(ib, 2) = lat_rad_out
  buoy_props(ib, 3) = lon_deg_out
  buoy_props(ib, 4) = lat_deg_out
  buoy_props(ib, 5) = frozen_in
  buoy_props(ib, 6) = dudt_out
  buoy_props(ib, 7) = dvdt_out
  buoy_props(ib, 8) = u_ib_out
  buoy_props(ib, 9) = v_ib_out
  buoy_props(ib,10) = height_ib_single
  buoy_props(ib,11) = length_ib_single
  buoy_props(ib,12) = width_ib_single
  buoy_props(ib,13) = iceberg_elem

! end if

 end if 
 
 t4=MPI_Wtime()
end subroutine iceberg_step2


!****************************************************************************************************************************
!****************************************************************************************************************************

subroutine initialize_velo(mesh,partit,dynamics, i_have_element, ib, u_ib, v_ib, lon_rad, lat_rad, depth_ib, localelem)			
 
 use g_rotate_grid,  only: vector_g2r
! use iceberg_params, only: l_initial, l_iniuser, ini_u, ini_v
implicit none
 
 logical, intent(in)	:: i_have_element
 integer, intent(in)	:: ib
 real, intent(inout) 	:: u_ib, v_ib
 real, intent(in)	:: lon_rad, lat_rad, depth_ib
 integer, intent(in)	:: localelem
 
 real, dimension(3)	:: startu, startv
 real 			:: ini_u_rot, ini_v_rot	

type(t_mesh), intent(in) , target :: mesh
type(t_dyn), intent(in) , target :: dynamics
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"



 !initialize ZERO for all PEs
 u_ib=0.
 v_ib=0. 
  
 !use initial velocities? 
 if(l_initial .AND. i_have_element) then
 
 	if(l_iniuser) then
   		ini_u_rot = ini_u(ib)	!still in geo. coord.
   		ini_v_rot = ini_v(ib)	!still in geo. coord.
		call vector_g2r(ini_u_rot, ini_v_rot, lon_rad, lat_rad, 0)
		u_ib = ini_u_rot
		v_ib = ini_v_rot	
	else
   		!OCEAN VELOCITY uo_ib, voib is start velocity
   		call iceberg_avvelo(mesh, partit, dynamics, startu,startv,depth_ib,localelem)
        call FEM_3eval(mesh, partit,u_ib,v_ib,lon_rad,lat_rad,startu,startv,localelem)
	end if
 end if
 
end subroutine initialize_velo


!****************************************************************************************************************************
!****************************************************************************************************************************

subroutine trajectory( lon_rad,lat_rad, old_u,old_v, new_u,new_v, &
		       lon_deg,lat_deg,old_lon,old_lat, dt_ib )			
 use o_param 	!for dt, r_earth
 implicit none
 
 real, intent(inout) 	:: lon_rad,lat_rad
 real, intent(inout)    :: old_u,old_v
 real, intent(in)	:: new_u,new_v
 real, intent(out)	:: lon_deg,lat_deg,old_lon,old_lat
 real, intent(in)	:: dt_ib
 
 real :: deltax1, deltay1, deltax2, deltay2	
 
 !save old position in case the iceberg leaves the domain
 old_lon = lon_rad
 old_lat = lat_rad
  
 !displacement vectors
 deltax1 = old_u * dt_ib
 deltay1 = old_v * dt_ib
 deltax2 = new_u * dt_ib
 deltay2 = new_v * dt_ib
   
 !heun method
 lon_rad = lon_rad + (0.5*(deltax1 + deltax2) / (r_earth*cos(lat_rad)) )
 lat_rad = lat_rad + (0.5*(deltay1 + deltay2) /  r_earth )
 lon_deg=lon_rad/rad
 lat_deg=lat_rad/rad
   
 !update velocity here (old value was needed for heun method)
 old_u=new_u
 old_v=new_v

end subroutine trajectory


!****************************************************************************************************************************
!****************************************************************************************************************************

subroutine depth_bathy(mesh, partit,Zdepth3, elem)  
  use o_param
  use g_clock
  use g_forcing_arrays
  use g_rotate_grid
  
  implicit none
  
  real, dimension(3), intent(OUT) 	:: Zdepth3 !depth in column below element
  integer, intent(IN)			:: elem	  !local element
  integer				:: m, n2, k, n_low
  
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
  Zdepth3=0.0
  
  ! loop over all nodes of the iceberg element
  do m=1, 3
   !for each 2D node of the iceberg element..
   n2=mesh%elem2D_nodes(m,elem)

   !k=num_layers_below_nod2d(n2)+1
   !n_low= nod3d_below_nod2d(k,   n2)	!deepest node below n2
   k=mesh%nlevels_nod2D(n2)

   !..compute depth below this node: 
   !Zdepth3(m) = abs(coord_nod3D(3, n_low))
   !Zdepth3(m) = abs(mesh%Z_3d_n_ib(k, n2))
   Zdepth3(m) = abs(mesh%zbar_n_bot(n2))
   !Zdepth3(m) = abs(mesh%zbar(k))
   !if (!Zdepth3(m)<0.0) then
   !    Zdepth3(m) = -Zdepth3(m)
   !end if
 end do
  
end subroutine depth_bathy


!****************************************************************************************************************************
!****************************************************************************************************************************

subroutine parallel2coast(mesh, partit,u, v, lon,lat, elem)
 use iceberg_params, only: coastal_nodes
 implicit none
 
 real, intent(inout) 	:: u, v 	!velocity
 real, intent(in)	:: lon, lat 	!radiant
 integer, intent(in)	:: elem
 
 integer :: fld_tmp
 integer, dimension(3) :: n
 integer :: node, m, i
 real, dimension(2) :: velocity, velocity1, velocity2
 real :: d1, d2
 
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
  if( use_cavity ) then
    !fld_tmp = coastal_nodes(mesh, partit, elem)
    fld_tmp =  count( (mesh%cavity_depth(mesh%elem2D_nodes(:,elem))/=0.0) .OR. (mesh%bc_index_nod2D(mesh%elem2D_nodes(:,elem))==0.0) )
  else
    fld_tmp =  count( (mesh%bc_index_nod2D(mesh%elem2D_nodes(:,elem))==0.0) )
  end if
  
  SELECT CASE ( fld_tmp ) !num of coastal points
   CASE (0) !...coastal points: do nothing
    return
    
   CASE (1) !...coastal point
   n = 0
   i = 1
   velocity = [ u, v ]
    do m = 1, 3
      node = mesh%elem2D_nodes(m,elem)
      !write(*,*) 'index ', m, ':', index_nod2D(node)
      if( use_cavity ) then
        !if( mesh%bc_index_nod2D(node)==1 .OR. cavity_flag_nod2d(node)==1 ) then
        if( mesh%bc_index_nod2D(node)==0.0 .OR.  (mesh%cavity_depth(node)/=0.0) ) then
         n(i) = node
         exit
        end if 
      else
        if( mesh%bc_index_nod2D(node)==1 ) then
         n(i) = node
         exit
        end if 
      end if
    end do 
    
   !write(*,*) 'one coastal node ', n(1)
  
  !LA comment for testing
   !i = 2 
   !if ( n(1) <= myDim_nod2D ) then	!all neighbours known

   ! do m = 1, nghbr_nod2D(n(1))%nmb
   !   node = nghbr_nod2D(n(1))%addresses(m) 
!#ifdef use_cavity
   !   if ( (node /= n(1)) .and. ( (bc_index_nod2D(node)==1) .OR. (cavity_flag_nod2d(node)==1) ) ) then   
!#else
   !   if ( (node /= n(1)) .and. (bc_index_nod2D(node)==1)) then
!#endif
   !    n(i) = node
   !    i = i+1
   !    if(i==4) exit
   !   end if
   ! end do
    
    !write(*,*) 'nodes n(i) ', n
    
    d1 = sqrt( (lon - coord_nod2D(1, n(2)))**2 + (lat - coord_nod2D(2, n(2)))**2 )
    d2 = sqrt( (lon - coord_nod2D(1, n(3)))**2 + (lat - coord_nod2D(2, n(3)))**2 )
    !write(*,*) 'distances :' , d1, d2
    !write(*,*) 'velocity vor :' , velocity
    if (d1 < d2) then
      call projection(mesh,partit,  velocity, n(2), n(1))
    else
      call projection(mesh,partit,  velocity, n(3), n(1))
    end if
    !write(*,*) 'velocity nach:', velocity
    !call projection(velocity, n(3), n(2))
      
   !else
   ! !if coastal point is not first node of element, the coastal point could be in eDim_nod2D,
   ! !so not all neighbours of this node are known to PE. WHAT SHOULD BE DONE?
   !end if    
    
    
   CASE (2) !...coastal points
    n = 0
    i = 1
    velocity = [ u, v ]
    do m = 1, 3
      node = mesh%elem2D_nodes(m,elem) 
      if( use_cavity ) then
        !if( (mesh%bc_index_nod2D(node)==1) .OR. (cavity_flag_nod2d(node)==1)) then
        if( mesh%bc_index_nod2D(node)==0.0 .OR.  (mesh%cavity_depth(node)/=0.0) ) then
         n(i) = node
         i = i+1
        end if
      else
        if( mesh%bc_index_nod2D(node)==1 ) then
         n(i) = node
         i = i+1
        end if
      end if
    end do   
    call projection(mesh,partit,  velocity, n(1), n(2))
    
   
   CASE DEFAULT 
    return  	!mesh element MUST NOT have 3 coastal points!

 END SELECT
 
 u = velocity(1)
 v = velocity(2)

end subroutine parallel2coast


!****************************************************************************************************************************
!****************************************************************************************************************************


subroutine projection(mesh, partit, velocity, n1, n2)
implicit none
 
 real, dimension(2), intent(inout) :: velocity
 integer, intent(in) :: n1, n2  
 
 real, dimension(2) :: direction
 real :: length, sp  

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

 ! direction: node1 - node2 (pointing from 2 to 1)
 direction(1) = coord_nod2D(1, n1) - coord_nod2D(1, n2)
 direction(2) = coord_nod2D(2, n1) - coord_nod2D(2, n2)
 length = sqrt( direction(1)**2 + direction(2)**2 )
 direction(1) = direction(1)/length
 direction(2) = direction(2)/length
    
 sp = sum( velocity(:) * direction(:) )	!ab = b'|a|, b' scalar projection of vec b on vec a
 velocity = sp * direction

end subroutine projection


!****************************************************************************************************************************
!****************************************************************************************************************************


subroutine iceberg_restart(partit)
! use iceberg_params 
 use g_config, only : ib_num

 implicit none
 integer :: icbID, ib
 LOGICAL :: file_exists
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_part_ass.h"
 INQUIRE(FILE=IcebergRestartPath, EXIST=file_exists) 
 icbID = mype+10
 
 !call allocate_icb
 
 if(file_exists) then
  open(unit=icbID,file=IcebergRestartPath,status='old', form='formatted')
  
  do ib=1, ib_num 
  
   !read all parameters that icb_step needs:			
   read(icbID,'(18e15.7,I8,L,3e15.7,L,I5,L)')						&
   	height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),	&
	Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
	conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
	u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib),		&
	f_u_ib_old(ib), f_v_ib_old(ib), calving_day(ib), grounded(ib), scaling(ib), melted(ib)

  end do
  close(icbID)

  if(mype==0) then
  write(*,*) 'read iceberg restart file'

  !if(.NOT.ascii_out) call determine_save_count ! computed from existing records in netcdf file
  if(.NOT.ascii_out) call init_buoy_output(partit)
  !call init_icebergs_with_icesheet ! all PEs read LON,LAT,LENGTH from files

  !write(*,*) '*************************************************************'
  end if
 else

  if(mype==0) then
  write(*,*) 'no iceberg restart'

  if(.NOT.ascii_out) call init_buoy_output(partit)

  end if

  !call init_buoys ! all PEs read LON,LAT from files  
  !write(*,*) 'initialized positions from file'
  call init_icebergs ! all PEs read LON,LAT,LENGTH from files
  !write(*,*) 'initialized positions and length/width from file'
  !write(*,*) '*************************************************************'
 end if

end subroutine iceberg_restart


!****************************************************************************************************************************
!****************************************************************************************************************************


subroutine iceberg_restart_with_icesheet(partit)
! use iceberg_params 
 use g_config, only : ib_num

 implicit none
 integer :: icbID_ISM, icbID_non_melted_icb, ib, st
 LOGICAL :: file_exists, file_exists_non_melted
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_part_ass.h"
 
 INQUIRE(FILE=num_non_melted_icb_file, EXIST=file_exists_non_melted) 
 INQUIRE(FILE=IcebergRestartPath_ISM, EXIST=file_exists) 
 icbID_ISM = mype+10
 icbID_non_melted_icb = mype+11

 !call allocate_icb
 
 open(unit=icbID_non_melted_icb,file=num_non_melted_icb_file,status='old', form='formatted')
    read(icbID_non_melted_icb,*) num_non_melted_icb
 close(icbID_non_melted_icb)

 if(file_exists) then
  open(unit=icbID_ISM,file=IcebergRestartPath_ISM,status='old', form='formatted')
  do ib=1, num_non_melted_icb 
   !read all parameters that icb_step needs:			
   read(icbID_ISM,'(18e15.7,I8,L,3e15.7,L,I5,L)',iostat=st)						&
   	height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),	&
	Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
	conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
	u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib),		&
	f_u_ib_old(ib), f_v_ib_old(ib), calving_day(ib), grounded(ib), scaling(ib), melted(ib)
  end do
  close(icbID_ISM)

  if(mype==0) then
  write(*,*) 'read iceberg restart file'

  !if(.NOT.ascii_out) call determine_save_count ! computed from existing records in netcdf file
  if(.NOT.ascii_out) call init_buoy_output(partit)
  end if
  call init_icebergs_with_icesheet
  !write(*,*) 'initialized positions and length/width from file'
  !write(*,*) '*************************************************************'
 else

  if(mype==0) then
  write(*,*) 'no iceberg restart'

  if(.NOT.ascii_out) call init_buoy_output(partit)

  end if

  !call init_buoys ! all PEs read LON,LAT from files  
  !write(*,*) 'initialized positions from file'
  call init_icebergs !_with_icesheet ! all PEs read LON,LAT,LENGTH from files
  !write(*,*) 'initialized positions and length/width from file'
  !write(*,*) '*************************************************************'
 end if

end subroutine iceberg_restart_with_icesheet


!****************************************************************************************************************************
!****************************************************************************************************************************


subroutine iceberg_out(partit)
! use iceberg_params
 use g_clock		!for dayold
 implicit none
 integer :: icbID, icbID_ISM, ib, istep
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_part_ass.h"
 
 icbID = 42
 icbID_ISM = 43
 
 !calving_day has to be adjusted for restarts because calving_day gives the amount
 !of days (since the model FIRST has been started) after which icebergs are released
 !Criterion for calving is:
 !if( real(istep) > real(step_per_day)*calving_day(ib) -1 ) then !iceberg calved

! kh 10.02.21 istep is not initialized
!calving_day = calving_day - REAL(istep/step_per_day)
 where(calving_day <= 0.0)
 calving_day = 0.0	!to avoid negative calving_days
 end where

 if(mype==0) then
  open(unit=icbID,file=IcebergRestartPath,position='append', status='replace', form='formatted')
  open(unit=icbID_ISM,file=IcebergRestartPath_ISM,position='append', status='replace', form='formatted')
  
  do ib=1, ib_num 
  
   !write all parameters that icb_step needs:
   write(icbID,'(18e15.7,I8,L,3e15.7,L,I5,L)')						&
   	height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),	&
	Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
	conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
	u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib),		&
	f_u_ib_old(ib), f_v_ib_old(ib), calving_day(ib), grounded(ib), scaling(ib), melted(ib)
   
   !***************************************************************
   !write new restart file with only non melted icebergs
   if(.not.melted(ib)) then
       !write all parameters that icb_step needs:
       write(icbID_ISM,'(18e15.7,I8,L,3e15.7,L,I5,L)')						&
            height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),	&
            Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
            conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
            u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib),		&
            f_u_ib_old(ib), f_v_ib_old(ib), calving_day(ib), grounded(ib), scaling(ib), melted(ib)
   end if

  end do
  close(icbID_ISM)
  close(icbID)
 end if
end subroutine iceberg_out

!========================================================================
! reads lon and lat values for buoys start position
! from files
! written by Madlen Kimmritz, 25.07.2015
!========================================================================
subroutine init_buoys
! use iceberg_params
 use g_config

 implicit none
 integer :: i
 integer :: io_error

!buoys_xlon_file > lon_deg
 open(unit=97, file=buoys_xlon_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file buoys_xlon_file'
 do i = 1, ib_num
    read(97,*) lon_deg(i)
 end do
 close (97)
!buoys_ylat_file > lat_deg
 open(unit=98, file=buoys_ylat_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file buoys_ylat_file'
 do i = 1, ib_num
    read(98,*) lat_deg(i)
 end do
 close(98)
end subroutine init_buoys
!
!--------------------------------------------------------------------------------------------
!

!========================================================================
! reads lon and lat values for iceberg starting positions
! from files as well as their length and width
! written by Madlen Kimmritz, 25.07.2015
! added length for iceberg case, 07.10.2015
!========================================================================
subroutine init_icebergs
! use iceberg_params
 use g_config
! use MOD_PARTIT

 implicit none
 integer :: i, myunit
 integer :: io_error, io_read
!type(t_partit), intent(inout), target :: partit
!#include "associate_part_def.h"
!#include "associate_part_ass.h"

!buoys_xlon_file > lon_deg
 open(newunit=myunit, file=buoys_xlon_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file buoys_xlon_file'
 do i = 1, ib_num
    read(myunit,*,iostat=io_read) lon_deg(i)
 end do
 close (myunit)
!buoys_ylat_file > lat_deg
 open(unit=98, file=buoys_ylat_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file buoys_ylat_file'
 do i = 1, ib_num
    read(98,*) lat_deg(i)
 end do
 close(98)
!length_icb_file > length_ib
 open(unit=97, file=length_icb_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file length_icb_file'
 do i = 1, ib_num
    read(97,*) length_ib(i)
 end do
 close (97)
!width_icb_file > lat_deg
 open(unit=98, file=width_icb_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file width_icb_file'
 do i = 1, ib_num
    read(98,*) width_ib(i)
 end do
 close(98)
!height_icb_file > height_ib
 open(unit=97, file=height_icb_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file width_icb_file'
 do i = 1, ib_num
    read(97,*) height_ib(i)
 end do
 close(97)
!scaling_file
 open(unit=98, file=scaling_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file scaling_file'
 do i = 1, ib_num
    read(98,*) scaling(i)
 end do
 close(98)

end subroutine init_icebergs
!
!--------------------------------------------------------------------------------------------
!

subroutine init_icebergs_with_icesheet
! use iceberg_params
 use g_config
! use MOD_PARTIT

 implicit none
 integer :: i
 integer :: io_error
!type(t_partit), intent(inout), target :: partit
!#include "associate_part_def.h"
!#include "associate_part_ass.h"

!buoys_xlon_file > lon_deg
 open(unit=97, file=buoys_xlon_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file buoys_xlon_file'
 do i = 1+num_non_melted_icb, ib_num
    read(97,*) lon_deg(i)
 end do
 close (97)
!buoys_ylat_file > lat_deg
 open(unit=98, file=buoys_ylat_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file buoys_ylat_file'
 do i = 1+num_non_melted_icb, ib_num
    read(98,*) lat_deg(i)
 end do
 close(98)
!length_icb_file > length_ib
 open(unit=97, file=length_icb_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file length_icb_file'
 do i = 1+num_non_melted_icb, ib_num
    read(97,*) length_ib(i)
 end do
 close (97)
!width_icb_file > lat_deg
 open(unit=98, file=width_icb_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file width_icb_file'
 do i = 1+num_non_melted_icb, ib_num
    read(98,*) width_ib(i)
 end do
 close(98)
!height_icb_file > height_ib
 open(unit=97, file=height_icb_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file width_icb_file'
 do i = 1+num_non_melted_icb, ib_num
    read(97,*) height_ib(i)
 end do
 close(97)
!scaling_file
 open(unit=98, file=scaling_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file scaling_file'
 do i = 1+num_non_melted_icb, ib_num
    read(98,*) scaling(i)
 end do
 close(98)

end subroutine init_icebergs_with_icesheet
!
!--------------------------------------------------------------------------------------------
!

subroutine determine_save_count(partit)
  ! computes save_count_buoys and prev_sec_in_year from records in existing netcdf file
  !-----------------------------------------------------------  
  use g_clock
!  use iceberg_params, only : file_icb_netcdf, save_count_buoys, prev_sec_in_year
  !use iceberg_params, only : save_count_buoys, prev_sec_in_year
  implicit none

#include "netcdf.inc" 
  integer		    :: buoy_nrec
  integer                   :: status, ncid
  integer                   :: dimid_rec
  integer                   :: time_varid
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_part_ass.h"

 if (mype==0) then

  ! open file
  status = nf_open(file_icb_netcdf, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! inquire time dimension ID and its length
  status = nf_inq_dimid(ncid, 'time', dimid_rec)
  if(status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_dimlen(ncid, dimid_rec, buoy_nrec)
  if(status .ne. nf_noerr) call handle_err(status)

  ! the next buoy/iceberg record to be saved
  save_count_buoys=buoy_nrec+1

  write(*,*) 'next record is #',save_count_buoys

  ! load sec_in_year up to now in 'prev_sec_in_year', time axis will be continued
  status=nf_inq_varid(ncid, 'time', time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_get_vara_double(ncid, time_varid, save_count_buoys-1, 1, prev_sec_in_year) 
  if (status .ne. nf_noerr) call handle_err(status)

  write(*,*) 'seconds passed up to now: ',prev_sec_in_year

  !close file
  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

end if

end subroutine determine_save_count
!=============================================================================================

!
!--------------------------------------------------------------------------------------------
!
subroutine init_buoy_output(partit)
  ! Initialize output file for buoys/icebergs
  ! written by Madlen Kimmritz, 24.07.2015
  ! reviewed by T. Rackow, 14.08.2015
  !-----------------------------------------------------------  
  use g_clock
  use g_config, only : ib_num
!  use iceberg_params, only : file_icb_netcdf, save_count_buoys !ggf in namelist
  implicit none

#include "netcdf.inc" 
  integer                   :: status,ncid,year_start,month_start,day_start
  integer                   :: dimid_ib, dimid_rec, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: lonrad_id, latrad_id, londeg_id, latdeg_id
  integer                   :: frozen_id, dudt_id, dvdt_id
  integer                   :: uib_id, vib_id
  integer                   :: height_id, length_id, width_id
  integer                   :: bvl_id, lvlv_id, lvle_id, lvlb_id, felem_id
  character(100)            :: longname, att_text, description
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_part_ass.h"

 if (mype==0) then
  write(*,*) 'initialize new buoy/iceberg output file'
 

 ! create a file
  status = nf_create(file_icb_netcdf, nf_clobber, ncid)
  if (status.ne.nf_noerr)  call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'number_tracer', ib_num, dimid_ib)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'time', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)


  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for the tracers.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_ib
  dimids(2) = dimid_rec


  status = nf_def_var(ncid, 'pos_lon_rad', NF_DOUBLE, 2, dimids, lonrad_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'pos_lat_rad', NF_DOUBLE, 2, dimids, latrad_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'pos_lon_deg', NF_DOUBLE, 2, dimids, londeg_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'pos_lat_deg', NF_DOUBLE, 2, dimids, latdeg_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'frozen_in', NF_DOUBLE, 2, dimids, frozen_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'du_dt', NF_DOUBLE, 2, dimids, dudt_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'dv_dt', NF_DOUBLE, 2, dimids, dvdt_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'icb_vel_u', NF_DOUBLE, 2, dimids, uib_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'icb_vel_v', NF_DOUBLE, 2, dimids, vib_id)
  if (status .ne. nf_noerr) call handle_err(status)

  ! 3 dimensions of the iceberg, comment for buoy case

  status = nf_def_var(ncid, 'height', NF_DOUBLE, 2, dimids, height_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'length', NF_DOUBLE, 2, dimids, length_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'width', NF_DOUBLE, 2, dimids, width_id)
  if (status .ne. nf_noerr) call handle_err(status)

  ! 4 additional iceberg variables (meltrates), comment for buoy case

  status = nf_def_var(ncid, 'bvl', NF_DOUBLE, 2, dimids, bvl_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'lvlv', NF_DOUBLE, 2, dimids, lvlv_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'lvle', NF_DOUBLE, 2, dimids, lvle_id)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'lvlb', NF_DOUBLE, 2, dimids, lvlb_id)
  if (status .ne. nf_noerr) call handle_err(status)

  ! LA: add felem
  status = nf_def_var(ncid, 'felem', NF_DOUBLE, 2, dimids, felem_id)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Assign long_name and units attributes to variables.
  longname='time' ! use NetCDF Climate and Forecast (CF) Metadata Convention
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)'), 'seconds since ', year_start, '-', month_start, '-', day_start, ' 00:00:00'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'units', len_trim(att_text), trim(att_text))
  if (status .ne. nf_noerr) call handle_err(status)
  if (include_fleapyear) then
      att_text='standard'
  else
      att_text='noleap'
  end if
  status = nf_put_att_text(ncid, time_varid, 'calendar', len_trim(att_text), trim(att_text))
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='longitude of buoy/iceberg position in radiant'
  status = nf_PUT_ATT_TEXT(ncid, lonrad_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, lonrad_id, 'units', 7, 'radiant')
  if (status .ne. nf_noerr) call handle_err(status)
  !rotated or not rotated due to setting in iceberg module
  description='(un)rotated according to setting of l_geo_out in iceberg module'
  status = nf_put_att_text(ncid, lonrad_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='latitude of buoy/iceberg position in radiant'
  status = nf_PUT_ATT_TEXT(ncid, latrad_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, latrad_id, 'units', 7, 'radiant')
  if (status .ne. nf_noerr) call handle_err(status)
  !rotated or not rotated due to setting in iceberg module
  description='(un)rotated according to setting of l_geo_out in iceberg module'
  status = nf_put_att_text(ncid, latrad_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='longitude of buoy/iceberg position in degree'
  status = nf_PUT_ATT_TEXT(ncid, londeg_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, londeg_id, 'units', 12, 'degrees_east')
  if (status .ne. nf_noerr) call handle_err(status)
  !rotated or not rotated due to setting in iceberg module
  description='(un)rotated according to setting of l_geo_out in iceberg module'
  status = nf_put_att_text(ncid, londeg_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='latitude of buoy/iceberg position in degree'
  status = nf_PUT_ATT_TEXT(ncid, latdeg_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, latdeg_id, 'units', 13, 'degrees_north')
  if (status .ne. nf_noerr) call handle_err(status)
  !rotated or not rotated due to setting in iceberg module
  description='(un)rotated according to setting of l_geo_out in iceberg module'
  status = nf_put_att_text(ncid, latdeg_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='status of buoy/iceberg (frozen in/not frozen in)'
  status = nf_PUT_ATT_TEXT(ncid, frozen_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  description='1 = frozen, 0 = not frozen, else partially frozen'
  status = nf_put_att_text(ncid, frozen_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, frozen_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='du/dt of buoy/iceberg in last time step'
  status = nf_PUT_ATT_TEXT(ncid, dudt_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, dudt_id, 'units', 8, 'm s^(-2)')
  if (status .ne. nf_noerr) call handle_err(status)  
  !rotated or not rotated due to setting in iceberg module
  description='(un)rotated according to setting of l_geo_out in iceberg module'
  status = nf_put_att_text(ncid, dudt_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, dudt_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='dv/dt of buoy/iceberg in last time step'
  status = nf_PUT_ATT_TEXT(ncid, dvdt_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, dvdt_id, 'units', 8, 'm s^(-2)')
  if (status .ne. nf_noerr) call handle_err(status)  
  !rotated or not rotated due to setting in iceberg module
  description='(un)rotated according to setting of l_geo_out in iceberg module'
  status = nf_put_att_text(ncid, dvdt_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, dvdt_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='velocity of buoy/iceberg, u component'
  status = nf_PUT_ATT_TEXT(ncid, uib_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, uib_id, 'units', 8, 'm s^(-1)')
  if (status .ne. nf_noerr) call handle_err(status)  
  !rotated or not rotated due to setting in iceberg module
  description='(un)rotated according to setting of l_geo_out in iceberg module'
  status = nf_put_att_text(ncid, uib_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, uib_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='velocity of buoy/iceberg, v component'
  status = nf_PUT_ATT_TEXT(ncid, vib_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, vib_id, 'units', 8, 'm s^(-1)')
  if (status .ne. nf_noerr) call handle_err(status)  
  !rotated or not rotated due to setting in iceberg module
  description='(un)rotated according to setting of l_geo_out in iceberg module'
  status = nf_put_att_text(ncid, vib_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, vib_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  ! 3 dimensions of the iceberg, comment for buoy case

  longname='height of the iceberg'
  status = nf_PUT_ATT_TEXT(ncid, height_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, height_id, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)  
  description='freeboard + draft'
  status = nf_put_att_text(ncid, height_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, height_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='length of the iceberg'
  status = nf_PUT_ATT_TEXT(ncid, length_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, length_id, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)  
  description='open'
  status = nf_put_att_text(ncid, length_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, length_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='width of the iceberg'
  status = nf_PUT_ATT_TEXT(ncid, width_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, width_id, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)  
  description='open'
  status = nf_put_att_text(ncid, width_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, width_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)


  ! 4 additional iceberg variables (meltrates), comment for buoy case

  longname='basal volume loss'
  status = nf_PUT_ATT_TEXT(ncid, bvl_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, bvl_id, 'units', 18, 'm^3 (ice) day^(-1)')
  if (status .ne. nf_noerr) call handle_err(status)
  description='losses are averaged over the preceding output interval'
  status = nf_put_att_text(ncid, bvl_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, bvl_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='lateral volume loss due to 1) bouyant convection'
  status = nf_PUT_ATT_TEXT(ncid, lvlv_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, lvlv_id, 'units', 18, 'm^3 (ice) day^(-1)')
  if (status .ne. nf_noerr) call handle_err(status)
  description='losses are averaged over the preceding output interval'
  status = nf_put_att_text(ncid, lvlv_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, lvlv_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='lateral volume loss due to 2) wave erosion'
  status = nf_PUT_ATT_TEXT(ncid, lvle_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, lvle_id, 'units', 18, 'm^3 (ice) day^(-1)')
  if (status .ne. nf_noerr) call handle_err(status)
  description='losses are averaged over the preceding output interval'
  status = nf_put_att_text(ncid, lvle_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, lvle_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='lateral volume loss due to 3) "basal" formulation'
  status = nf_PUT_ATT_TEXT(ncid, lvlb_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, lvlb_id, 'units', 18, 'm^3 (ice) day^(-1)')
  if (status .ne. nf_noerr) call handle_err(status)
  description='losses are averaged over the preceding output interval'
  status = nf_put_att_text(ncid, lvlb_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, lvlb_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  ! LA: add felem
  longname='fesom element'
  status = nf_PUT_ATT_TEXT(ncid, felem_id, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, felem_id, 'units', 18, '')
  if (status .ne. nf_noerr) call handle_err(status)
  description=''
  status = nf_put_att_text(ncid, felem_id, 'description', len_trim(description), trim(description)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_PUT_ATT_TEXT(ncid, felem_id, 'Coordinates', 23, 'pos_lon_deg pos_lat_deg') ! arcGIS 
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status) 

 ! initialize the counter for saving results
  save_count_buoys=1
  write(*,*) 'initialize new buoy/iceberg output file: done'

 end if
end subroutine init_buoy_output
!=============================================================================================

!=============================================================================================

subroutine write_buoy_props_netcdf(partit)
  ! write output file for buoys/icebergs for current time step
  ! written by Madlen Kimmritz, 25.07.2015
  ! reviewed by T. Rackow, 17.08.2015
  !-----------------------------------------------------------  

 use g_config
 use g_clock
! use iceberg_params, only : buoy_props, file_icb_netcdf, save_count_buoys, prev_sec_in_year, bvl_mean, lvlv_mean, lvle_mean, lvlb_mean
 use g_forcing_param

 implicit none

#include "netcdf.inc" 

  integer                   :: status,ncid, istep
  integer                   :: dimid_ib, dimid_rec, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: lonrad_id, latrad_id, londeg_id, latdeg_id
  integer                   :: frozen_id, dudt_id, dvdt_id
  integer                   :: uib_id, vib_id  
  integer                   :: height_id, length_id, width_id
  integer                   :: bvl_id, lvlv_id, lvle_id, lvlb_id, felem_id
  integer                   :: start(2), count(2)
  real(kind=8)              :: sec_in_year
type(t_partit), intent(inout), target :: partit
!type(t_ice),    intent(inout), target :: ice
#include "associate_part_def.h"
#include "associate_part_ass.h"

!   /gfs1/work/hbkkim15/output/slabt/buoys_track.nc
  
  if (mype==0) then 

! kh 16.03.21 ?! istep is not initialized, intitialize to 0 here
     istep = 0
     sec_in_year=dt*istep

    ! open files
     status = nf_open(trim(file_icb_netcdf), nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id

     status = nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'pos_lon_rad', lonrad_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'pos_lat_rad', latrad_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'pos_lon_deg', londeg_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'pos_lat_deg', latdeg_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'frozen_in',  frozen_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'du_dt',  dudt_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'dv_dt',  dvdt_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'icb_vel_u',  uib_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'icb_vel_v',  vib_id)
     if (status .ne. nf_noerr) call handle_err(status)  

     ! inquire 3 additional IDs for iceberg case, comment for buoy case

     status = nf_inq_varid(ncid, 'height',  height_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'length',  length_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'width',  width_id)
     if (status .ne. nf_noerr) call handle_err(status)  

     ! inquire 4 additional IDs for iceberg case, comment for buoy case

     status = nf_inq_varid(ncid, 'bvl',  bvl_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'lvlv',  lvlv_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'lvle',  lvle_id)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_inq_varid(ncid, 'lvlb',  lvlb_id)
     if (status .ne. nf_noerr) call handle_err(status) 

     ! * LA: include fesom elemt in output
     status = nf_inq_varid(ncid, 'felem', felem_id)
     if (status .ne. nf_noerr) call handle_err(status) 

     !buoy_props(ib, 1) = lon_rad_out
     !buoy_props(ib, 2) = lat_rad_out
     !buoy_props(ib, 3) = lon_deg_out
     !buoy_props(ib, 4) = lat_deg_out
     !buoy_props(ib, 5) = frozen_in
     !buoy_props(ib, 6) = dudt_out
     !buoy_props(ib, 7) = dvdt_out
     !buoy_props(ib, 8) = u_ib_out
     !buoy_props(ib, 9) = v_ib_out
     !buoy_props(ib,10) = height_ib
     !buoy_props(ib,11) = length_ib 
     !buoy_props(ib,12) = width_ib

     !bvl_mean(ib)*step_per_day
     !lvlv_mean(ib)*step_per_day
     !lvle_mean(ib)*step_per_day
     !lvlb_mean(ib)*step_per_day

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count_buoys, 1, prev_sec_in_year+sec_in_year) 
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count_buoys, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)

     !variables
     start=(/1,save_count_buoys/)
     count=(/ib_num, 1/)
     status=nf_put_vara_double(ncid, lonrad_id, start, count, buoy_props(:, 1)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, latrad_id, start, count, buoy_props(:, 2)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, londeg_id, start, count, buoy_props(:, 3)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, latdeg_id, start, count, buoy_props(:, 4)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, frozen_id, start, count, buoy_props(:, 5)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, dudt_id, start, count, buoy_props(:, 6)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, dvdt_id, start, count, buoy_props(:, 7)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, uib_id, start, count, buoy_props(:, 8)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, vib_id, start, count, buoy_props(:, 9)) 
     if (status .ne. nf_noerr) call handle_err(status)

     ! write 3 additional variables for iceberg case, comment for buoy case

     status=nf_put_vara_double(ncid, height_id, start, count, buoy_props(:,10)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, length_id, start, count, buoy_props(:,11)) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, width_id, start, count, buoy_props(:,12)) 
     if (status .ne. nf_noerr) call handle_err(status)

     ! write 4 additional variables for iceberg case, comment for buoy case

     status=nf_put_vara_double(ncid, bvl_id, start, count, bvl_mean(:)*step_per_day) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, lvlv_id, start, count, lvlv_mean(:)*step_per_day) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, lvle_id, start, count, lvle_mean(:)*step_per_day) 
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_put_vara_double(ncid, lvlb_id, start, count, lvlb_mean(:)*step_per_day) 
     if (status .ne. nf_noerr) call handle_err(status)

     ! LA: add felem
     status=nf_put_vara_double(ncid, felem_id, start, count, buoy_props(:,13)) 
     if (status .ne. nf_noerr) call handle_err(status)

     !close file
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     save_count_buoys=save_count_buoys+1
!==========================================================

  end if!mype==0


end subroutine write_buoy_props_netcdf
end module iceberg_step
