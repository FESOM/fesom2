subroutine iceberg_calculation  
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
 use g_parsup		!for mype					!=
 use iceberg_params	!for ..all others				!=
 									!=
 implicit none								!=
									!=
 integer	:: ib, times						!=	
 real(kind=8) 	:: t0, t1, t2, t3, t4					!=
 logical	:: firstcall=.true. 					!=
 logical	:: lastsubstep  					!=
 									!=
 real		:: arr_from_block(15)					!=
 integer	:: elem_from_block					!=  
 real		:: vl_from_block(4)					!=	
 real,dimension(15*ib_num):: arr_block_red				!=
 integer,dimension(ib_num):: elem_block_red				!=
 real, dimension(4*ib_num):: vl_block_red				!=
									!=
 !==================== MODULES & DECLARATIONS ==========================!= 

 if(firstcall) then
  !overwrite icb_modules if restart, initialize netcdf output if no restart:
  call iceberg_restart
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
  
    !substeps don't work anymore with new communication
    !do times=1, steps_per_FESOM_step
    !if(times == steps_per_FESOM_step) lastsubstep = .true. !do output at last substep
    lastsubstep = .true. !do output every timestep

    call iceberg_step1(	ib, height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),&
			Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
			conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
			u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib), lastsubstep,&
			steps_per_FESOM_step, f_u_ib_old(ib), f_v_ib_old(ib), l_semiimplicit,   &
			semiimplicit_coeff, AB_coeff)			
    !call MPI_Barrier(MPI_COMM_WORLD, MPIERR) !necessary?
    !end do
   
  end if
 end do
 
 t1=MPI_Wtime()
 
 
 !========================== COMMUNICATION =============================!

 !all PEs need the array arr(15) and the iceberg element
 !in step2
 !ALLREDUCE: arr_block, elem_block

 call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
 arr_block_red = 0.0
 elem_block_red= 0
 vl_block_red = 0.0
 call MPI_AllREDUCE(arr_block, arr_block_red, 15*ib_num, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
 call MPI_AllREDUCE(elem_block, elem_block_red, ib_num, MPI_INTEGER, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)  
 !ALLREDUCE: vl_block, containing the volume losses (IBs may switch PE during the output interval)
 call MPI_AllREDUCE(vl_block, vl_block_red, 4*ib_num, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
 
! call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
! arr_block_red = 0.0
! elem_block_red= 0
! vl_block_red = 0.0
! call MPI_AllREDUCE(arr_block, arr_block_red, 15*ib_num, MPI_DOUBLE_PRECISION, MPI_SUM, &
!       MPI_COMM_WORLD, MPIerr)
! call MPI_AllREDUCE(elem_block, elem_block_red, ib_num, MPI_INTEGER, MPI_SUM, &
!       MPI_COMM_WORLD, MPIerr)  
! !ALLREDUCE: vl_block, containing the volume losses (IBs may switch PE during the output interval)
! call MPI_AllREDUCE(vl_block, vl_block_red, 4*ib_num, MPI_DOUBLE_PRECISION, MPI_SUM, &
!       MPI_COMM_WORLD, MPIerr)
 
 !is filled in STEP 2
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
    !do times=1, steps_per_FESOM_step
    !if(times == steps_per_FESOM_step) lastsubstep = .true. !do output at last substep
    lastsubstep = .true. !do output every timestep
    
    call iceberg_step2(	arr_from_block, elem_from_block, ib, height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),&
			Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
			conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
			u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib), lastsubstep,&
			steps_per_FESOM_step, f_u_ib_old(ib), f_v_ib_old(ib), l_semiimplicit,   &
			semiimplicit_coeff, AB_coeff)	
    !call MPI_Barrier(MPI_COMM_WORLD, MPIERR) !necessary?
    !end do
  end if
end do
 
 !call MPI_Barrier(MPI_COMM_WORLD, MPIERR) !MK: necessary for output !commented 1 JAN 2000

 t3=MPI_Wtime() 


 !========================== VECTOR OUTPUT =============================!
 
 !call iceberg_vector_ncout !look in routines what is really written out!
 !istep, lon_deg_geo, lat_deg_geo, u_ib_geo, v_ib_geo, volume losses (set to zero)
 !introduce force_last_output(ib)?

 if (mod(istep,icb_outfreq)==0 .AND. ascii_out==.false.) then

   if (mype==0 .AND. (real(istep) > real(step_per_day)*calving_day(1) ) ) call write_buoy_props_netcdf
       
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
 write(*,*) '*************************************************************'
 end if
 
end subroutine iceberg_calculation


!****************************************************************************************************************************
!****************************************************************************************************************************


subroutine iceberg_step1(ib, height_ib,length_ib,width_ib, lon_deg,lat_deg, &
			Co,Ca,Ci, Cdo_skin,Cda_skin, rho_icb, 		   &
			conc_sill,P_sill, rho_h2o,rho_air,rho_ice,	   & 
			u_ib,v_ib, iceberg_elem, find_iceberg_elem, 	   &
			lastsubstep, steps_per_FESOM_step, f_u_ib_old,	   &
			f_v_ib_old, l_semiimplicit, semiimplicit_coeff,    &
			AB_coeff)
			
 !============================= MODULES & DECLARATIONS =========================================!=
 												!=
 use o_param 		!for rad								!=
 use o_arrays		!for coriolis_param_elem2D						!=
 use o_mesh		!for nod2D, (cavities: for cavity_flag_nod2d)				!=
 use g_parsup		!for myDim_elem2D, myList_nod2D						!=
 use g_rotate_grid	!for subroutine g2r, logfile_outfreq					!=
 												!=
#ifdef use_cavity
 use iceberg_params, only: smallestvol_icb, arr_block, elem_block, l_geo_out, icb_outfreq, l_allowgrounding, draft_scale, reject_elem	!=
#else
 use iceberg_params, only: smallestvol_icb, arr_block, elem_block, l_geo_out, icb_outfreq, l_allowgrounding, draft_scale			!=
#endif
 												!=
 implicit none											!=
 
 
 integer, intent(in)	:: ib
 real,    intent(inout)	:: height_ib,length_ib,width_ib
 real,    intent(inout)	:: lon_deg,lat_deg
 real, 	  intent(in)	:: Co,Ca,Ci, Cdo_skin,Cda_skin
 real,	  intent(in)	:: rho_icb, conc_sill,P_sill, rho_h2o,rho_air,rho_ice
 real,	  intent(inout)	:: u_ib,v_ib
 integer, intent(inout)	:: iceberg_elem !global
 logical, intent(inout)	:: find_iceberg_elem
 logical, intent(in)	:: lastsubstep 
 real, intent(in)	:: steps_per_FESOM_step
 real,    intent(inout)	:: f_u_ib_old, f_v_ib_old
 logical, intent(in)	:: l_semiimplicit
 real,    intent(in)	:: semiimplicit_coeff
 real,    intent(in)	:: AB_coeff
 
 
 integer, dimension(:), save, allocatable :: local_idx_of
 real      			:: depth_ib, volume_ib, mass_ib
 real				:: lon_rad, lat_rad, new_u_ib, new_v_ib
 real	   			:: old_lon,old_lat, frozen_in, P_ib, conci_ib
 real				:: lon_rad_out, lat_rad_out  !for unrotated output
 real				:: lon_deg_out, lat_deg_out  !for unrotated output
 integer   			:: i, iceberg_node  
 real 				:: dudt, dvdt
 
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
 real(kind=8) 			:: t0, t1, t2, t3, t4					!=
 											!=
 !for restart										!=
 logical, save   		:: firstcall=.true.					!=
 !for grounding										!=
 real, dimension(3)		:: Zdepth3						!=
 real				:: Zdepth						!=
 											!=
 !========================= MODULES & DECLARATIONS =====================================!=
 
 
 t0=MPI_Wtime()

 depth_ib = -height_ib * rho_icb/rho_h2o
 volume_ib= length_ib * width_ib * height_ib  
 mass_ib = volume_ib * rho_icb	 !less mass 
 lon_rad = lon_deg*rad
 lat_rad = lat_deg*rad
 
 if(volume_ib .le. smallestvol_icb) then
 
  if (mod(istep,logfile_outfreq)==0 .and. mype==0 .and. lastsubstep) then
   write(*,*) 'iceberg ', ib,' melted'
  end if
   
  return
 end if 
 
 if (firstcall) then
  if(mype==0) write(*,*) 'Preparing local_idx_of array...'
  allocate(local_idx_of(elem2D))
  !creates mapping
  call global2local(local_idx_of)
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
  call point_in_triangle(iceberg_elem, (/lon_deg, lat_deg/))
  i_have_element= (iceberg_elem .ne. 0) !up to 3 PEs possible
  
  if(i_have_element) then
   i_have_element= elem2D_nodes(1,iceberg_elem) <= myDim_nod2D !1 PE still .true.
#ifdef use_cavity
   if(reject_elem(iceberg_elem)) then
    iceberg_elem=0 !reject element
    i_have_element=.false.
   else 
    iceberg_elem=myList_elem2D(iceberg_elem) !global now
   end if
#else
   iceberg_elem=myList_elem2D(iceberg_elem) !global now
#endif 
  end if
  call com_integer(i_have_element,iceberg_elem)
  
  if(iceberg_elem .EQ. 0) then
   	call par_ex
   	stop 'ICEBERG OUTSIDE MODEL DOMAIN OR IN ICE SHELF REGION'
  end if
  
  if(i_have_element) then 
   write(*,*) 'IB ',ib,' in element ',iceberg_elem
   write(*,*) 'IB ',ib,' rot. coords:', lon_deg, lat_deg !,lon_rad, lat_rad
   !do i=1,3
   ! iceberg_node=elem2D_nodes(i,local_idx_of(iceberg_elem))
   ! write(*,*) 'node ', i, 's coords lon:', coord_nod2d(1,iceberg_node)/rad,' lat: ', coord_nod2d(2,iceberg_node)/rad
   !end do
  end if
  
  ! initialize the iceberg velocity
  call initialize_velo(i_have_element, ib, u_ib, v_ib, lon_rad, lat_rad, depth_ib, local_idx_of(iceberg_elem))

  !iceberg elem of ib is found
  find_iceberg_elem = .false.
  
  !for AB method
  f_u_ib_old = coriolis_param_elem2D(local_idx_of(iceberg_elem))*u_ib
  f_v_ib_old = coriolis_param_elem2D(local_idx_of(iceberg_elem))*v_ib
 end if
 
 file_track='/work/ollie/lackerma/iceberg/iceberg_ICBref_'
 file_forces_u='/work/ollie/lackerma/iceberg/iceberg_ICBref_forces_u_'
 file_forces_v='/work/ollie/lackerma/iceberg/iceberg_ICBref_forces_v_'
 file_meltrates='/work/ollie/lackerma/iceberg/iceberg_ICBref_melt_'

 !convert ib integer to string
 write(ib_char,'(I10)') ib
 
 !left-adjust the string..
 ib_char = adjustl(ib_char)
 
 !.. and trim while concatenating:
 file_track	=  trim(file_track) // trim(ib_char) // '.dat'
 file_forces_u	=  trim(file_forces_u) // trim(ib_char) // '.dat'
 file_forces_v	=  trim(file_forces_v) // trim(ib_char) // '.dat'
 file_meltrates	=  trim(file_meltrates) // trim(ib_char) // '.dat' 
 
 
 ! ================== START ICEBERG CALCULATION ====================
 
 arr=0.
 frozen_in = 0.
 i_have_element=.false.
 !if the first node belongs to this processor.. (just one processor enters here!)
 !if( local_idx_of(iceberg_elem) > 0 .and. elem2D_nodes(1,local_idx_of(iceberg_elem)) <= myDim_nod2D ) then
if( local_idx_of(iceberg_elem) > 0 ) then 
  if( elem2D_nodes(1,local_idx_of(iceberg_elem)) <= myDim_nod2D ) then

  i_have_element=.true. 
  l_output =  lastsubstep .and. mod(istep,icb_outfreq)==0
  
  !===========================DYNAMICS===============================
  
  call iceberg_dyn(ib, new_u_ib, new_v_ib, u_ib, v_ib, lon_rad,lat_rad, depth_ib, &
                   height_ib, length_ib, width_ib, local_idx_of(iceberg_elem), &
  		   mass_ib, Ci, Ca, Co, Cda_skin, Cdo_skin, &
  		   rho_ice, rho_air, rho_h2o, P_sill,conc_sill, frozen_in, &
  		   file_forces_u, file_forces_v, P_ib, conci_ib, &
		   dt/REAL(steps_per_FESOM_step), l_output, f_u_ib_old, &
		   f_v_ib_old, l_semiimplicit, semiimplicit_coeff, &
		   AB_coeff, file_meltrates, rho_icb)

  !new_u_ib = 2.0
  !new_v_ib = 0.0

  dudt = (new_u_ib-u_ib)*REAL(steps_per_FESOM_step) / dt
  dvdt = (new_v_ib-v_ib)*REAL(steps_per_FESOM_step) / dt		   
		   
  !new_u_ib = 0.
  !new_v_ib = -0.20
  !P_ib = 20000.
  !conci_ib = 0.3
		   
  !=======================END OF DYNAMICS============================
  
  call depth_bathy(Zdepth3, local_idx_of(iceberg_elem))
  !interpolate depth to location of iceberg (2 times because FEM_3eval expects a 2 component vector...)
  call FEM_3eval(Zdepth,Zdepth,lon_rad,lat_rad,Zdepth3,Zdepth3,local_idx_of(iceberg_elem))
  
  !write(*,*) 'nodal depth in iceberg ', ib,'s element:', Zdepth3
  !write(*,*) 'depth at iceberg ', ib, 's location:', Zdepth
  
  !=================CHECK IF ICEBERG IS GROUNDED...===================
 if((draft_scale(ib)*abs(depth_ib) .gt. Zdepth) .and. l_allowgrounding ) then 
   !icebergs remains stationary (iceberg can melt above in iceberg_dyn!)
    left_mype = 0.0 
    u_ib = 0.0
    v_ib = 0.0
    old_lon = lon_rad
    old_lat = lat_rad
    if (mod(istep,logfile_outfreq)==0) then 
    write(*,*) 'iceberg ib ', ib, 'is grounded'
    end if
 	
 else 
  !===================...ELSE CALCULATE TRAJECTORY====================
 
  call trajectory( lon_rad,lat_rad, u_ib,v_ib, new_u_ib,new_v_ib, &
		   lon_deg,lat_deg,old_lon,old_lat, dt/REAL(steps_per_FESOM_step))
  	   
  iceberg_elem=local_idx_of(iceberg_elem)  	!local
  call find_new_iceberg_elem(iceberg_elem, (/lon_deg, lat_deg/), left_mype)
  iceberg_elem=myList_elem2D(iceberg_elem)  	!global
  
  if(left_mype > 0.) then
   lon_rad = old_lon
   lat_rad = old_lat
   call parallel2coast(new_u_ib, new_v_ib, lon_rad,lat_rad, local_idx_of(iceberg_elem))
   call trajectory( lon_rad,lat_rad, new_u_ib,new_v_ib, new_u_ib,new_v_ib, &
		   lon_deg,lat_deg,old_lon,old_lat, dt/REAL(steps_per_FESOM_step))
   u_ib = new_u_ib
   v_ib = new_v_ib
		   
   iceberg_elem=local_idx_of(iceberg_elem)  	!local
   call find_new_iceberg_elem(iceberg_elem, (/lon_deg, lat_deg/), left_mype)
   iceberg_elem=myList_elem2D(iceberg_elem)  	!global
  end if		   
  !================END OF TRAJECTORY CALCULATION=====================
 end if ! iceberg stationary?
   
  !values for communication
  arr= (/ height_ib,length_ib,width_ib, u_ib,v_ib, lon_rad,lat_rad, &
          left_mype, old_lon,old_lat, frozen_in, dudt, dvdt, P_ib, conci_ib/) 
	
  !save in larger array	  
  arr_block((ib-1)*15+1 : ib*15)=arr
  elem_block(ib)=iceberg_elem
  	  
 end if !processor has element?
end if !... and first node belongs to processor?

 t1=MPI_Wtime()
 !if (mod(istep,logfile_outfreq)==0 .and. i_have_element .and. lastsubstep) write(*,*) 'dynamics  took', t1-t0
 ! =================== END OF ICEBERG CALCULATION ==================
 
 end subroutine iceberg_step1
 
 
subroutine iceberg_step2(arr, elem_from_block, ib, height_ib,length_ib,width_ib, lon_deg,lat_deg, &
			Co,Ca,Ci, Cdo_skin,Cda_skin, rho_icb, 		   &
			conc_sill,P_sill, rho_h2o,rho_air,rho_ice,	   & 
			u_ib,v_ib, iceberg_elem, find_iceberg_elem, 	   &
			lastsubstep, steps_per_FESOM_step, f_u_ib_old,	   &
			f_v_ib_old, l_semiimplicit, semiimplicit_coeff,    &
			AB_coeff)
			
 !============================= MODULES & DECLARATIONS =========================================!=
 												!=
 use o_param 		!for rad								!=
 use o_arrays		!for coriolis_param_elem2D						!=
 use o_mesh		!for nod2D, (cavities: for cavity_flag_nod2d)				!=
 use g_parsup		!for myDim_elem2D, myList_nod2D						!=
 use g_rotate_grid	!for subroutine g2r, logfile_outfreq					!=
 												!=
#ifdef use_cavity
 use iceberg_params, only: smallestvol_icb, buoy_props, bvl_mean, lvlv_mean, lvle_mean, lvlb_mean, ascii_out, l_geo_out, icb_outfreq, l_allowgrounding, draft_scale, reject_elem	!=
#else
 use iceberg_params, only: smallestvol_icb, buoy_props, bvl_mean, lvlv_mean, lvle_mean, lvlb_mean, ascii_out, l_geo_out, icb_outfreq, l_allowgrounding, draft_scale		!=
#endif
 												!=
 implicit none											!=
 
 real, 	  intent(in)	:: arr(15)
 integer, intent(in)	:: elem_from_block
 integer, intent(in)	:: ib
 real,    intent(inout)	:: height_ib,length_ib,width_ib
 real,    intent(inout)	:: lon_deg,lat_deg
 real, 	  intent(in)	:: Co,Ca,Ci, Cdo_skin,Cda_skin
 real,	  intent(in)	:: rho_icb, conc_sill,P_sill, rho_h2o,rho_air,rho_ice
 real,	  intent(inout)	:: u_ib,v_ib
 integer, intent(inout)	:: iceberg_elem !global
 logical, intent(inout)	:: find_iceberg_elem
 logical, intent(in)	:: lastsubstep 
 integer, intent(in)	:: steps_per_FESOM_step
 real,    intent(inout)	:: f_u_ib_old, f_v_ib_old
 logical, intent(in)	:: l_semiimplicit
 real,    intent(in)	:: semiimplicit_coeff
 real,    intent(in)	:: AB_coeff						
 
 
 integer, dimension(:), save, allocatable :: local_idx_of
 real      			:: depth_ib, volume_ib, mass_ib
 real				:: lon_rad, lat_rad, new_u_ib, new_v_ib
 real	   			:: old_lon,old_lat, frozen_in, P_ib, conci_ib
 real				:: lon_rad_out, lat_rad_out  !for unrotated output
 real				:: lon_deg_out, lat_deg_out  !for unrotated output
 real				:: u_ib_out, v_ib_out  	     !for unrotated output
 real				:: dudt_out, dvdt_out  	     !for unrotated output
 integer   			:: i, iceberg_node  
 real 				:: dudt, dvdt
 
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
 											!=
 !========================= MODULES & DECLARATIONS =====================================!=
 
 !all PEs enter here with identical array arr
  
 !**** check if iceberg melted in step 1 ****! 
 !mass_ib = arr(1) * arr(2) * arr(3) * rho_icb
 !if(mass_ib .le. 1.0e-6) then
 ! return
 !end if 
 
 !**** check if iceberg melted in step 1 ****!
 !call com_values(i_have_element, arr, iceberg_elem) 
 iceberg_elem= elem_from_block !update element as before in com_values
 old_element = elem_from_block !save if iceberg left model domain
 height_ib= arr(1)
 length_ib= arr(2)
 width_ib = arr(3)
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
 volume_ib = height_ib * length_ib * width_ib ! * rho_icb
 if(volume_ib .le. smallestvol_icb) then
  buoy_props(ib, :) = 0. ! for output: NaN or MissVal could be written here
  return
 end if 
 !**** check if iceberg melted in step 1 ****!

 t2=MPI_Wtime()
 
 if(left_mype > 0.) then
 call iceberg_elem4all(iceberg_elem, lon_deg, lat_deg) !Just PE changed?
 
  if(iceberg_elem == 0) then !IB left model domain
  
   !if (mod(istep,logfile_outfreq)==0 .and. mype==0 .and. lastsubstep) write(*,*) 'iceberg ',ib, ' left model domain'
   lon_rad = old_lon
   lat_rad = old_lat 
   lon_deg = lon_rad/rad
   lat_deg = lat_rad/rad
   iceberg_elem = old_element
   u_ib    = 0.
   v_ib    = 0.  
   !if(i_have_element) then
    !OCEAN VELOCITY uo_ib, voib is new velocity
    !call iceberg_avvelo(startu,startv,depth_ib,local_idx_of(iceberg_elem))
    !call FEM_3eval(u_ib,v_ib,lon_rad,lat_rad,startu,startv,local_idx_of(iceberg_elem))
   !end if
   
  else
   !if (mype==0) write(*,*) 'iceberg ',ib, ' changed PE or was very fast'
  end if
 end if
 
 t3=MPI_Wtime()
 

 !if(mype==0 .and. lastsubstep .and. mod(istep,icb_outfreq)==0 .and. ascii_out) then
 !if(mype==mod(ib,npes-1) .and. lastsubstep .and. mod(istep,icb_outfreq)==0 .and. ascii_out) then
 if(mype==0 .and. lastsubstep .and. mod(istep,icb_outfreq)==0) then

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

 if(ascii_out) then !use old ASCII output
  
 file_track='/work/ollie/lackerma/iceberg/iceberg_ICBref_'
 !convert ib integer to string
 write(ib_char,'(I10)') ib
 !left-adjust the string..
 ib_char = adjustl(ib_char)
 !.. and trim while concatenating:
 file_track	=  trim(file_track) // trim(ib_char) // '.dat'
   
   open(unit=42,file=file_track,position='append')
   write(42,'(I,12e15.7)') 	istep, lon_rad_out, lat_rad_out, lon_deg_out, lat_deg_out, &
   				u_ib_out, v_ib_out, frozen_in, P_sill, P_ib, conci_ib, dudt_out, dvdt_out
   close(42)

 else !write in array for faster netcdf output

  buoy_props(ib, 1) = lon_rad_out
  buoy_props(ib, 2) = lat_rad_out
  buoy_props(ib, 3) = lon_deg_out
  buoy_props(ib, 4) = lat_deg_out
  buoy_props(ib, 5) = frozen_in
  buoy_props(ib, 6) = dudt_out
  buoy_props(ib, 7) = dvdt_out
  buoy_props(ib, 8) = u_ib_out
  buoy_props(ib, 9) = v_ib_out
  buoy_props(ib,10) = height_ib
  buoy_props(ib,11) = length_ib 
  buoy_props(ib,12) = width_ib

 end if

 end if 
 
 t4=MPI_Wtime()
 
 !if (mod(istep,logfile_outfreq)==0 .and. mype==0 .and. lastsubstep) then 
  !write(*,*) 'comvalues took', t2-t1
  !write(*,*) 'left mype took', t3-t2
  !write(*,*) 'track out took', t4-t3
 !end if

end subroutine iceberg_step2


!****************************************************************************************************************************
!****************************************************************************************************************************

subroutine initialize_velo(i_have_element, ib, u_ib, v_ib, lon_rad, lat_rad, depth_ib, localelem)			
 
 use g_rotate_grid,  only: vector_g2r
 use iceberg_params, only: l_initial, l_iniuser, ini_u, ini_v
 implicit none
 
 logical, intent(in)	:: i_have_element
 integer, intent(in)	:: ib
 real, intent(inout) 	:: u_ib, v_ib
 real, intent(in)	:: lon_rad, lat_rad, depth_ib
 integer, intent(in)	:: localelem
 
 real, dimension(3)	:: startu, startv
 real 			:: ini_u_rot, ini_v_rot	
 
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
   		call iceberg_avvelo(startu,startv,depth_ib,localelem)
   		call FEM_3eval(u_ib,v_ib,lon_rad,lat_rad,startu,startv,localelem)
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

subroutine depth_bathy(Zdepth3, elem)  
  use o_mesh
  use o_param
  use i_therm_param
  use i_param
  use i_arrays
  use g_parsup
  
  use o_arrays         
  use g_clock
  use g_forcing_index
  use g_forcing_arrays
  use g_rotate_grid
  
  implicit none
  
  real, dimension(3), intent(OUT) 	:: Zdepth3 !depth in column below element
  integer, intent(IN)			:: elem	  !local element
  integer				:: m, n2, k, n_low
  
  
  Zdepth3=0.0
  
  ! loop over all nodes of the iceberg element
  do m=1, 3
   !for each 2D node of the iceberg element..
   n2=elem2D_nodes(m,elem)

   k=num_layers_below_nod2d(n2)+1
   n_low= nod3d_below_nod2d(k,   n2)	!deepest node below n2

   !..compute depth below this node: 
   Zdepth3(m) = abs(coord_nod3D(3, n_low))
  end do
  
end subroutine depth_bathy


!****************************************************************************************************************************
!****************************************************************************************************************************

subroutine parallel2coast(u, v, lon,lat, elem)
 use o_mesh		!for index_nod2D, (cavities: for cavity_flag_nod2d)
 use g_parsup		!for myDim_nod2D
#ifdef use_cavity
 use iceberg_params, only: coastal_nodes
#endif
 implicit none
 
 real, intent(inout) 	:: u, v 	!velocity
 real, intent(in)	:: lon, lat 	!radiant
 integer, intent(in)	:: elem
 
 integer, dimension(3) :: n
 integer :: node, m, i
 real, dimension(2) :: velocity, velocity1, velocity2
 real :: d1, d2
 
#ifdef use_cavity
  SELECT CASE ( coastal_nodes(elem) ) !num of "coastal" points
#else
  SELECT CASE ( sum( index_nod2D(elem2D_nodes(:,elem)) ) ) !num of coastal points
#endif
   CASE (0) !...coastal points: do nothing
    return
    
   CASE (1) !...coastal point
   n = 0
   i = 1
   velocity = (/ u, v /)
    do m = 1, 3
      node = elem2D_nodes(m,elem)
      !write(*,*) 'index ', m, ':', index_nod2D(node)
#ifdef use_cavity
      if( index_nod2D(node)==1 .OR. cavity_flag_nod2d(node)==1 ) then
#else
      if( index_nod2D(node)==1 ) then
#endif
       n(i) = node
       exit
      end if 
    end do 
    
   !write(*,*) 'one coastal node ', n(1)
   
   i = 2 
   if ( n(1) <= myDim_nod2D ) then	!all neighbours known
    do m = 1, nghbr_nod2D(n(1))%nmb
      node = nghbr_nod2D(n(1))%addresses(m) 
#ifdef use_cavity
      if ( (node /= n(1)) .and. ( (index_nod2D(node)==1) .OR. (cavity_flag_nod2d(node)==1) ) ) then   
#else
      if ( (node /= n(1)) .and. (index_nod2D(node)==1)) then
#endif
       n(i) = node
       i = i+1
       if(i==4) exit
      end if
    end do
    
    !write(*,*) 'nodes n(i) ', n
    
    d1 = sqrt( (lon - coord_nod2D(1, n(2)))**2 + (lat - coord_nod2D(2, n(2)))**2 )
    d2 = sqrt( (lon - coord_nod2D(1, n(3)))**2 + (lat - coord_nod2D(2, n(3)))**2 )
    !write(*,*) 'distances :' , d1, d2
    !write(*,*) 'velocity vor :' , velocity
    if (d1 < d2) then
      call projection(velocity, n(2), n(1))
    else
      call projection(velocity, n(3), n(1))
    end if
    !write(*,*) 'velocity nach:', velocity
    !call projection(velocity, n(3), n(2))
      
   else
    !if coastal point is not first node of element, the coastal point could be in eDim_nod2D,
    !so not all neighbours of this node are known to PE. WHAT SHOULD BE DONE?
   end if    
    
    
   CASE (2) !...coastal points
    n = 0
    i = 1
    velocity = (/ u, v /)
    do m = 1, 3
      node = elem2D_nodes(m,elem) 
#ifdef use_cavity
      if( (index_nod2D(node)==1) .OR. (cavity_flag_nod2d(node)==1)) then
#else
      if( index_nod2D(node)==1 ) then
#endif
       n(i) = node
       i = i+1
      end if
    end do   
    call projection(velocity, n(1), n(2))
    
   
   CASE DEFAULT 
    return  	!mesh element MUST NOT have 3 coastal points!

 END SELECT
 
 u = velocity(1)
 v = velocity(2)

end subroutine parallel2coast


!****************************************************************************************************************************
!****************************************************************************************************************************


subroutine projection(velocity, n1, n2)
 use o_mesh		!for coord_nod2D
 implicit none
 
 real, dimension(2), intent(inout) :: velocity
 integer, intent(in) :: n1, n2  
 
 real, dimension(2) :: direction
 real :: length, sp  
 
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


subroutine iceberg_restart
 use iceberg_params 
 use g_parsup		!for mype
 use g_config, only : ib_num

 implicit none
 integer :: icbID, ib
 LOGICAL :: file_exists
 INQUIRE(FILE=IcebergRestartPath, EXIST=file_exists) 
 icbID = mype+10
 
 call allocate_icb
 
 if(file_exists) then
  open(unit=icbID,file=IcebergRestartPath,status='old')
  
  do ib=1, ib_num 
  
   !read all parameters that icb_step needs:			
   read(icbID,'(18e15.7,I,L,3e15.7)')						&
   	height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),	&
	Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
	conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
	u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib),		&
	f_u_ib_old(ib), f_v_ib_old(ib), calving_day(ib) 

  end do
  close(icbID)

  if(mype==0) then
  write(*,*) 'read iceberg restart file'

  if(.NOT.ascii_out) call determine_save_count ! computed from existing records in netcdf file

  write(*,*) '*************************************************************'
  end if
 else

  if(mype==0) then
  write(*,*) 'no iceberg restart'

  if(.NOT.ascii_out) call init_buoy_output

  end if

  !call init_buoys ! all PEs read LON,LAT from files  
  !write(*,*) 'initialized positions from file'
  call init_icebergs ! all PEs read LON,LAT,LENGTH from files
  write(*,*) 'initialized positions and length/width from file'
  write(*,*) '*************************************************************'
 end if

end subroutine iceberg_restart


!****************************************************************************************************************************
!****************************************************************************************************************************


subroutine iceberg_out 
 use iceberg_params
 use g_parsup		!for mype
 use g_clock		!for dayold
 implicit none
 integer :: icbID, ib
 
 icbID = 42
 
 !calving_day has to be adjusted for restarts because calving_day gives the amount
 !of days (since the model FIRST has been started) after which icebergs are released
 !Criterion for calving is:
 !if( real(istep) > real(step_per_day)*calving_day(ib) -1 ) then !iceberg calved
 calving_day = calving_day - REAL(istep/step_per_day)
 where(calving_day <= 0.0)
 calving_day = 0.0	!to avoid negative calving_days
 end where
 
 if(mype==0) then
  open(unit=icbID,file=IcebergRestartPath,position='append', status='replace')
  
  do ib=1, ib_num 
  
   !write all parameters that icb_step needs:
   write(icbID,'(18e15.7,I,L,3e15.7)')						&
   	height_ib(ib),length_ib(ib),width_ib(ib), lon_deg(ib),lat_deg(ib),	&
	Co(ib),Ca(ib),Ci(ib), Cdo_skin(ib),Cda_skin(ib), rho_icb(ib), 		&
	conc_sill(ib),P_sill(ib), rho_h2o(ib),rho_air(ib),rho_ice(ib),	   	& 
	u_ib(ib),v_ib(ib), iceberg_elem(ib), find_iceberg_elem(ib),		&
	f_u_ib_old(ib), f_v_ib_old(ib), calving_day(ib)
	
  end do		
  close(icbID)
 end if
end subroutine iceberg_out

!========================================================================
! reads lon and lat values for buoys start position
! from files
! written by Madlen Kimmritz, 25.07.2015
!========================================================================
subroutine init_buoys
 use iceberg_params
 use g_config
 use g_parsup

 implicit none
 integer :: i
 integer :: io_error

!buoys_xlon_file > lon_deg
 open(unit=99, file=buoys_xlon_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file buoys_xlon_file'
 do i = 1, ib_num
    read(99,*) lon_deg(i)
 end do
 close (99)
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
 use iceberg_params
 use g_config
 use g_parsup

 implicit none
 integer :: i
 integer :: io_error

!buoys_xlon_file > lon_deg
 open(unit=99, file=buoys_xlon_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file buoys_xlon_file'
 do i = 1, ib_num
    read(99,*) lon_deg(i)
 end do
 close (99)
!buoys_ylat_file > lat_deg
 open(unit=98, file=buoys_ylat_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file buoys_ylat_file'
 do i = 1, ib_num
    read(98,*) lat_deg(i)
 end do
 close(98)
!length_icb_file > length_ib
 open(unit=99, file=length_icb_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file length_icb_file'
 do i = 1, ib_num
    read(99,*) length_ib(i)
 end do
 close (99)
!width_icb_file > lat_deg
 open(unit=98, file=width_icb_file,status='old',action='read',iostat=io_error)
 if ( io_error.ne.0) stop 'ERROR while reading file width_icb_file'
 do i = 1, ib_num
    read(98,*) width_ib(i)
 end do
 close(98)

end subroutine init_icebergs
!
!--------------------------------------------------------------------------------------------
!

subroutine determine_save_count
  ! computes save_count_buoys and prev_sec_in_year from records in existing netcdf file
  !-----------------------------------------------------------  
  use g_clock
  use iceberg_params, only : file_icb_netcdf, save_count_buoys, prev_sec_in_year
  !use iceberg_params, only : save_count_buoys, prev_sec_in_year
  use g_parsup
  implicit none

#include "netcdf.inc" 
  integer		    :: buoy_nrec
  integer                   :: status, ncid
  integer                   :: dimid_rec
  integer                   :: time_varid

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
subroutine init_buoy_output
  ! Initialize output file for buoys/icebergs
  ! written by Madlen Kimmritz, 24.07.2015
  ! reviewed by T. Rackow, 14.08.2015
  !-----------------------------------------------------------  
  use g_clock
  use g_config, only : ib_num
  use iceberg_params, only : file_icb_netcdf, save_count_buoys !ggf in namelist
  !use iceberg_params, only : ib_num, save_count_buoys !ggf in namelist
  use g_parsup
  implicit none

#include "netcdf.inc" 
  integer                   :: status,ncid
  integer                   :: dimid_ib, dimid_rec, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: lonrad_id, latrad_id, londeg_id, latdeg_id
  integer                   :: frozen_id, dudt_id, dvdt_id
  integer                   :: uib_id, vib_id
  integer                   :: height_id, length_id, width_id
  integer                   :: bvl_id, lvlv_id, lvle_id, lvlb_id
  character(100)            :: longname, att_text, description

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

subroutine write_buoy_props_netcdf
  ! write output file for buoys/icebergs for current time step
  ! written by Madlen Kimmritz, 25.07.2015
  ! reviewed by T. Rackow, 17.08.2015
  !-----------------------------------------------------------  

 use g_config
 use g_clock
 use g_parsup
 use iceberg_params, only : buoy_props, file_icb_netcdf, save_count_buoys, prev_sec_in_year, bvl_mean, lvlv_mean, lvle_mean, lvlb_mean
 !use iceberg_params, only : ib_num, buoy_props, save_count_buoys, prev_sec_in_year, bvl_mean, lvlv_mean, lvle_mean, lvlb_mean
  
  use o_arrays
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use i_arrays
  use g_forcing_param

 implicit none

#include "netcdf.inc" 

  integer                   :: status,ncid
  integer                   :: dimid_ib, dimid_rec, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: lonrad_id, latrad_id, londeg_id, latdeg_id
  integer                   :: frozen_id, dudt_id, dvdt_id
  integer                   :: uib_id, vib_id  
  integer                   :: height_id, length_id, width_id
  integer                   :: bvl_id, lvlv_id, lvle_id, lvlb_id
  integer                   :: start(2), count(2)
  real(kind=8)              :: sec_in_year

!   /gfs1/work/hbkkim15/output/slabt/buoys_track.nc
  
  if (mype==0) then 
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


     !close file
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     save_count_buoys=save_count_buoys+1
!==========================================================

  end if!mype==0


end subroutine write_buoy_props_netcdf
