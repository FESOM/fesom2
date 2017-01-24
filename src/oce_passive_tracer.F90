module o_passive_tracer_mod
  ! Ocean passive tracer module
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------

  use o_param
  use o_arrays
  use o_mesh
  use g_config
  use g_forcing_arrays
  use g_clock
  use g_parsup
  implicit none

  integer, allocatable, dimension(:)         :: index_passive_tracer
  integer, allocatable, dimension(:,:)       :: passive_tracer_loc_index
  real(kind=WP), allocatable, dimension(:,:) :: ptr_sfc_force

contains

! =================================================================
  subroutine passive_tracer_init

    integer              :: i, j, k, n, n3, row, fileID
    integer              :: n_loc, num_nod
    integer, allocatable :: temp_arr2d(:), nodes_release(:)
    character(1)         :: cptrind
    character(4)         :: tr_name
    character(100)       :: file_name

    allocate(index_passive_tracer(num_passive_tracer))
    do j=1, num_passive_tracer
       index_passive_tracer(j)=j
    end do

    !--------------------------------------------------------------  
    ! initial values
    do j=1, num_passive_tracer
       tracer(:,:,index_passive_tracer(j))=ptr_background_value
    end do

    !--------------------------------------------------------------
    ! in case ptr is restored in a region
    if(passive_tracer_restore) then
       ! set passive tracer location index: 1 at release, 0 otherwise
       ! We label only horizontal locations, the rest can be done  
       ! in the restoring routine
       allocate(passive_tracer_loc_index(myDim_nod2d+eDim_nod2D,num_passive_tracer))
       passive_tracer_loc_index=0

       allocate(temp_arr2d(nod2d))
       temp_arr2d=0
       do n=1, myDim_nod2D+eDim_nod2D
          temp_arr2d(myList_nod2D(n))=n
       end do

       do j=1, num_passive_tracer
          write(cptrind,'(i1)') j
          tr_name='ptr'//cptrind
          file_name=trim(meshpath)//'passive_tracer_restore_nodes_'//tr_name//'.out'
          fileID=160
          open(fileID, file=file_name)
          read(fileID,*) num_nod
          allocate(nodes_release(num_nod))
          read(fileID,*) nodes_release
          close(fileID)
          do n=1,num_nod
             n_loc=temp_arr2d(nodes_release(n))
             if(n_loc>0) then
                passive_tracer_loc_index(n_loc,j)=1
                tracer(1,n_loc,index_passive_tracer(j))=ptr_restore_value
             end if
          end do
          deallocate(nodes_release)
       end do

       deallocate(temp_arr2d)

       !--------------------------------------------------------------
       ! in case restore volume    
       ! We should take care of restoring in the procedure 
       if(ptr_restore_in_volume) then
          do j=1, num_passive_tracer 
            do n=1,myDim_nod2D+eDim_nod2D
             if(passive_tracer_loc_index(n_loc,j)==1) then
               tracer(:,n,index_passive_tracer(j))=ptr_restore_value
             end if
            end do
          end do
       end if
    
    end if 
    !--------------------------------------------------------------
    ! in case that passive tracers enter through surface fluxes
    if(passive_tracer_flux) then
       allocate(ptr_sfc_force(myDim_nod2D+eDim_nod2D,num_passive_tracer))
    end if

    !backup 
    !tracer0(:,index_passive_tracer)=tracer(:,index_passive_tracer)

  end subroutine passive_tracer_init
  !
  !-------------------------------------------------------------------------
  !
  subroutine ptr_sfc_bc
    implicit none
    integer      :: n,j

    if(.not.passive_tracer_flux) return
       ptr_sfc_force=0.0
       do n=1,myDim_nod2d+eDim_nod2D             
          do j=1,num_passive_tracer
             ptr_sfc_force(n,j)=-area(1,n)*(tr_arr(1,n,2)+tracer(n,index_passive_tracer(j),1)) &
                  * runoff_landice(n)*landice_season(month)
          end do
       end do
  end subroutine ptr_sfc_bc
  !
  !-------------------------------------------------------------------------
  !
  subroutine ptr_cutoff_restore
    implicit none
    integer   :: j, n
    if(.not.passive_tracer_restore) return
    ! TODO: the information about 62 should be external, ideally 
    ! in  the same file where the restoring nodes are
    if(ptr_restore_in_volume) then               
    do j=1, num_passive_tracer
          do n=1,myDim_nod2d+eDim_nod2D
             if(passive_tracer_loc_index(n,j)==0) cycle
                tracer(:,n,index_passive_tracer)=ptr_background_value  !enforce others to be zero
                tracer(:,n,index_passive_tracer(j))=ptr_restore_value
                if(geo_coord_nod2D(2,n)<62.0*rad) then
	        tracer(:,n,index_passive_tracer(j))=ptr_background_value  !enforce each to be zero
	        endif	
         end do
     end do
     else
     do j=1, num_passive_tracer
          do n=1,myDim_nod2d+eDim_nod2D
             if(passive_tracer_loc_index(n,j)==0) cycle
                tracer(1,n,index_passive_tracer)=ptr_background_value  !enforce others to be zero
                tracer(1,n,index_passive_tracer(j))=ptr_restore_value
                if(geo_coord_nod2D(2,n)<62.0*rad) then
	        tracer(1,n,index_passive_tracer(j))=ptr_background_value  !enforce each to be zero
	        endif	
         end do
     end do
     end if
  end subroutine ptr_cutoff_restore


end module o_passive_tracer_mod
