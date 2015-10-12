module o_age_tracer_mod
  ! Ocean age tracer module
  ! Test version!
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------

  use o_param
  use o_arrays
  use o_mesh
  use g_config
  use g_clock
  use g_parsup
  implicit none

  integer, allocatable, dimension(:)   :: index_age_tracer  
  integer, allocatable, dimension(:,:) :: age_tracer_loc_index

contains


  subroutine age_tracer_init

    integer              :: i, j, k, n, n3, row, fileID
    integer              :: n_loc, num_nod
    integer, allocatable :: temp_arr2d(:), nodes_release(:)
    character(1)         :: cageind
    character(4)         :: tr_name
    character(100)       :: file_name

    !--------------------------------------------------------------  
    ! find index
    allocate(index_age_tracer(num_age_tracer))
    if(.not.use_passive_tracer) num_passive_tracer=0
    do j=1, num_age_tracer
       index_age_tracer(j)=j+num_passive_tracer
    end do

    !--------------------------------------------------------------  
    ! set initial values
    do j=1, num_age_tracer
       tracer(:,:,index_age_tracer(j))=0.0
    end do

    !-------------------------------------------------------------- 
    ! restore time scale at the release region
    if(zero_age_at_release) age_tracer_restore_time=dt

    !--------------------------------------------------------------    	
    ! set age tracer location index: 1 at release, 0 otherwise

    allocate(age_tracer_loc_index(myDim_nod2D+eDim_nod2D,num_age_tracer))
    age_tracer_loc_index=0

    allocate(temp_arr2d(nod2d))
    temp_arr2d=0
    do n=1, myDim_nod2D+eDim_nod2D
       temp_arr2d(myList_nod2D(n))=n
    end do

    do j=1, num_age_tracer
       if(age_tracer_global_surf) then
         age_tracer_loc_index=1
       else
          write(cageind,'(i1)') j
          tr_name='age'//cageind
          file_name=trim(meshpath)//'age_tracer_release_nodes_'//tr_name//'.out'
          fileID=160
          open(fileID, file=file_name)
          read(fileID,*) num_nod
          allocate(nodes_release(num_nod))
          read(fileID,*) nodes_release
          close(fileID)

          do n=1,num_nod
             n_loc=temp_arr2d(nodes_release(n))
             if(n_loc>0) then
                age_tracer_loc_index(n_loc,j)=1
             end if
          end do
          deallocate(nodes_release)
       end if
    end do

    deallocate(temp_arr2d)

    !--------------------------------------------------------------
    ! in case release in volume

    !if(age_release_in_volume) then
    !end if

    !backup
    !tracer0(:,index_age_tracer)=tracer(:,index_age_tracer)

  end subroutine age_tracer_init
  !
  !-------------------------------------------------------------------------
  !
  subroutine age_tracer_tendency

    integer      :: j, n, nz, column_ind(nl-1)
    real(kind=8) :: inside_val, outside_val

    outside_val=1.0/(86400.0*(365+fleapyear))
    
    do j=1, num_age_tracer
    do n=1, myDim_nod2D+eDim_nod2D
       if(age_release_in_volume) then
         column_ind=1
       else
         column_ind=0
         column_ind(1)=1
       end if
         column_ind=column_ind*age_tracer_loc_index(n,j) ! 1-inside, 0-outside   
       do nz=1,nlevels_nod2d(n)-1
          inside_val=-tracer(nz,n,index_age_tracer(j))/age_tracer_restore_time  
          !assume initial value 0, restore to it in the above formula
          tracer_rhs(nz,n,index_age_tracer(j))=tracer_rhs(nz,n,index_age_tracer(j)) &
               + inside_val*column_ind(nz)+outside_val*(1-column_ind(nz))
       end do      
    end do
    end do
    ! TODO: look carefully at what stage this update is done 
    ! (volume multiplicationmay not be required) 
  end subroutine age_tracer_tendency
  !
  !-------------------------------------------------------------------------
  !
  subroutine age_tracer_cutoff_restore

    integer   :: j, n, nz, column_ind(nl)

    do j=1, num_age_tracer
       do n=1,myDim_nod2D+eDim_nod2D
          column_ind=1          
          if(zero_age_at_release .and. age_tracer_loc_index(n,j)==1) then
            if(age_release_in_volume) then
              column_ind=0  
            else
             column_ind(1)=0
            end if 
          end if
          do nz=1, nlevels_nod2D(n)-1 
             tracer(nz,n,index_age_tracer(j))= &
               max(column_ind(nz)*tracer(nz,n,index_age_tracer(j)),0.0_WP)
          end do
	end do  
    end do

  end subroutine age_tracer_cutoff_restore

end module o_age_tracer_mod
