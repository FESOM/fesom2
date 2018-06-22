!a set of auxuary routines for: 
!1. smoothing FESOM fields using mass matrix
!2. computing surface integrals of the FESOM fields
module g_support
  use o_mesh
  use g_parsup
  use g_comm_auto
  use o_ARRAYS
  implicit none

  private
  public :: smooth_nod, smooth_elem, integrate_nod
  real(kind=WP), dimension(:), allocatable  :: work_array
!
!--------------------------------------------------------------------------------------------
! generic interface for smoothing nodal fields
  INTERFACE smooth_nod
            MODULE PROCEDURE smooth_nod2D, smooth_nod3D
  END INTERFACE
!
!--------------------------------------------------------------------------------------------
! generic interface for smoothing fields on elements
  INTERFACE smooth_elem
            MODULE PROCEDURE smooth_elem2D, smooth_elem3D
  END INTERFACE
!
!--------------------------------------------------------------------------------------------
! computes 2D integral of a nodal field
  INTERFACE integrate_nod
            MODULE PROCEDURE integrate_nod_2D
  END INTERFACE
!
!--------------------------------------------------------------------------------------------
!
contains
!
!--------------------------------------------------------------------------------------------
!
subroutine smooth_nod2D(arr, N)
  IMPLICIT NONE
  integer, intent(in)                        :: N
  real(KIND=WP), dimension(:), intent(inout) :: arr
  integer                                    :: node, elem, j, q, elnodes(3)
  real(kind=WP)                              :: vol
  
  allocate(work_array(myDim_nod2D))
  DO q=1, N !apply mass matrix N times to smooth the field
     DO node=1, myDim_nod2D
        vol=0._WP
        work_array(node)=0._WP
        DO j=1, nod_in_elem2D_num(node)
           elem=nod_in_elem2D(j, node)
           elnodes=elem2D_nodes(:,elem)
           work_array(node)=work_array(node)+sum(arr(elnodes))/3._WP*elem_area(elem)
           vol=vol+elem_area(elem)
       END DO
       work_array(node)=work_array(node)/vol
    END DO
    DO node=1,myDim_nod2D
       arr(node)=work_array(node)
    ENDDO
    call exchange_nod(arr)
  END DO
  deallocate(work_array)
end subroutine smooth_nod2D

!--------------------------------------------------------------------------------------------
!
subroutine smooth_nod3D(arr, N_smooth)
  IMPLICIT NONE
  integer, intent(in)            :: N_smooth
  real(KIND=WP), intent(inout)   :: arr(:,:)

  integer                        :: n, el, nz, j, q, num_el, nlev, nl_loc
  real(kind=WP)                  :: vol(nl,myDim_nod2D)
  real(kind=WP), allocatable     :: work_array(:,:)

  nlev=ubound(arr,1)
  allocate(work_array(nlev,myDim_nod2D))
  
! Precompute area of patches on all levels (at the bottom, some neighbouring
! nodes may vanish in the bathymetry) in the first smoothing step
  DO n=1, myDim_nod2D
     
     vol(       1:min(nlev, nlevels_nod2d(n)),n) = 0._WP
     work_array(1:min(nlev, nlevels_nod2d(n)),n) = 0._WP

     DO j=1, nod_in_elem2D_num(n)
        el = nod_in_elem2D(j,n)
        nl_loc = min(nlev, minval(nlevels_nod2d(elem2D_nodes(1:3,el))))
        DO nz=1, nl_loc
           vol(nz,n) = vol(nz,n) + elem_area(el)
           work_array(nz,n) = work_array(nz,n) + elem_area(el) * (arr(nz, elem2D_nodes(1,el)) &
                                                                + arr(nz, elem2D_nodes(2,el)) &
                                                                + arr(nz, elem2D_nodes(3,el)))
        END DO
     ENDDO
     DO nz=1,nlevels_nod2d(n)
        vol(nz,n) = 1._WP / (3._WP * vol(nz,n))  ! Here, we need the inverse and scale by 1/3
     END DO
  END DO

  ! combined: scale by patch volume + copy back to original field
  DO n=1, myDim_nod2D
     DO nz=1, min(nlev, nlevels_nod2d(n))
        arr(nz, n) = work_array(nz, n) *vol(nz,n) 
     END DO
  end DO
  call exchange_nod(arr)

! And the remaining smoothing sweeps

  DO q=1,N_smooth-1
     DO n=1, myDim_nod2D

        work_array(1:min(nlev, nlevels_nod2d(n)),n) = 0._WP

        DO j=1,nod_in_elem2D_num(n)
           el = nod_in_elem2D(j,n)
           nl_loc = min(nlev, minval(nlevels_nod2d(elem2D_nodes(1:3,el))))
           DO nz=1, nl_loc
              work_array(nz,n) = work_array(nz,n) + elem_area(el) * (arr(nz, elem2D_nodes(1,el)) &
                                                                   + arr(nz, elem2D_nodes(2,el)) &
                                                                   + arr(nz, elem2D_nodes(3,el)))
           END DO
        ENDDO
     ENDDO
! combined: scale by patch volume + copy back to original field
     DO n=1, myDim_nod2D
        DO nz=1, min(nlev, nlevels_nod2d(n))
           arr(nz, n) = work_array(nz, n) *vol(nz,n) 
        END DO
     end DO
     call exchange_nod(arr)
  enddo

  deallocate(work_array)

end subroutine smooth_nod3D

!
!--------------------------------------------------------------------------------------------
!
subroutine smooth_elem2D(arr, N)
  IMPLICIT NONE
  integer, intent(in)                        :: N
  real(KIND=WP), dimension(:), intent(inout) :: arr
  integer                                    :: node, elem, j, q, elnodes(3)
  real(kind=WP)                              :: vol
  
  allocate(work_array(myDim_nod2D+eDim_nod2D))
  DO q=1, N !apply mass matrix N times to smooth the field
     DO node=1, myDim_nod2D+eDim_nod2D
        vol=0._WP
        work_array(node)=0._WP
        DO j=1, nod_in_elem2D_num(node)
           elem=nod_in_elem2D(j, node)
           elnodes=elem2D_nodes(:,elem)
           work_array(node)=work_array(node)+arr(elem)*elem_area(elem)
           vol=vol+elem_area(elem)
       END DO
       work_array(node)=work_array(node)/vol
    END DO

    DO elem=1, myDim_elem2D
       elnodes=elem2D_nodes(:, elem)
       arr(elem)=sum(work_array(elnodes))/3.0_WP  ! Here, we need the inverse and scale by 1/3
    ENDDO
    call exchange_elem(arr)
  END DO
  deallocate(work_array)
end subroutine smooth_elem2D
!
!--------------------------------------------------------------------------------------------
!
subroutine smooth_elem3D(arr, N)
  IMPLICIT NONE
  integer, intent(in)                          :: N
  real(KIND=WP), dimension(:,:), intent(inout) :: arr
  integer                                      :: node, elem, nl, nz, j, q, elnodes(3)
  real(kind=WP)                                :: vol
  
  allocate(work_array(myDim_nod2D+eDim_nod2D))
  nl=ubound(arr,1)
  DO q=1, N !apply mass matrix N times to smooth the field
     DO nz=1, nl
        DO node=1, myDim_nod2D+eDim_nod2D
           vol=0._WP
           work_array(node)=0._WP
           if (nlevels_nod2d(node) < nz) CYCLE
           DO j=1, nod_in_elem2D_num(node)
              elem=nod_in_elem2D(j, node)
              if (elem<=0) CYCLE
              elnodes=elem2D_nodes(:,elem)
              work_array(node)=work_array(node)+arr(nz, elem)*elem_area(elem)
              vol=vol+elem_area(elem)
          END DO
          work_array(node)=work_array(node)/vol
        END DO
     END DO

    DO elem=1, myDim_elem2D
       elnodes=elem2D_nodes(:, elem)
       arr(nz, elem)=sum(work_array(elnodes))/3.0_WP
    ENDDO
    call exchange_elem(arr)
  END DO
  deallocate(work_array)
end subroutine smooth_elem3D
!
!--------------------------------------------------------------------------------------------
!
subroutine integrate_nod_2D(data, int2D)
  use o_MESH
  use g_PARSUP
  use g_comm_auto

  IMPLICIT NONE
  real(kind=WP), intent(in)       :: data(:)
  real(kind=WP), intent(inout)    :: int2D

  integer       :: row
  real(kind=WP) :: lval

  lval=0.0
  do row=1, myDim_nod2D
     lval=lval+data(row)*area(1,row)
  end do

  int2D=0.0
  call MPI_AllREDUCE(lval, int2D, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
end subroutine integrate_nod_2D
end module g_support

