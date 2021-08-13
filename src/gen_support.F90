!a set of auxuary routines for: 
!1. smoothing FESOM fields using mass matrix
!2. computing surface integrals of the FESOM fields
module g_support
  USE MOD_MESH
  use g_parsup
  use g_comm_auto
  use o_ARRAYS
  use g_config, only: dummy
  implicit none

  private
  public :: smooth_nod, smooth_elem, integrate_nod, extrap_nod
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
            MODULE PROCEDURE integrate_nod_2D, integrate_nod_3D
  END INTERFACE
!
!--------------------------------------------------------------------------------------------
! fills dummy values in a scalar field
  INTERFACE extrap_nod
            MODULE PROCEDURE extrap_nod3D
  END INTERFACE
!
!--------------------------------------------------------------------------------------------
!
contains
!
!--------------------------------------------------------------------------------------------
!
subroutine smooth_nod2D(arr, N, mesh)
  IMPLICIT NONE
  type(t_mesh), intent(in)                  , target :: mesh
  integer, intent(in)                        :: N
  real(KIND=WP), dimension(:), intent(inout) :: arr
  integer                                    :: node, elem, j, q, elnodes(3)
  real(kind=WP)                              :: vol

#include "associate_mesh.h"
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
!
!--------------------------------------------------------------------------------------------
!
subroutine smooth_nod3D(arr, N_smooth, mesh)
  IMPLICIT NONE
  type(t_mesh), intent(in)      , target :: mesh
  integer, intent(in)            :: N_smooth
  real(KIND=WP), intent(inout)   :: arr(:,:)

  integer                        :: n, el, nz, j, q, num_el, nlev, nl_loc, nu_loc
  integer                        :: uln, nln, ule, nle
  real(kind=WP)                  :: vol(mesh%nl,myDim_nod2D)
  real(kind=WP), allocatable     :: work_array(:,:)

#include "associate_mesh.h"
  nlev=ubound(arr,1)
  allocate(work_array(nlev,myDim_nod2D))
  
! Precompute area of patches on all levels (at the bottom, some neighbouring
! nodes may vanish in the bathymetry) in the first smoothing step
  DO n=1, myDim_nod2D
     uln = ulevels_nod2d(n)
     nln = min(nlev,nlevels_nod2d(n))
     vol(       1:nln,n) = 0._WP
     work_array(1:nln,n) = 0._WP
     !!PS vol(       1:min(nlev, nlevels_nod2d(n)),n) = 0._WP
     !!PS work_array(1:min(nlev, nlevels_nod2d(n)),n) = 0._WP
     DO j=1, nod_in_elem2D_num(n)
        el = nod_in_elem2D(j,n)
!!PS         nl_loc = min(nlev, minval(nlevels_nod2d(elem2D_nodes(1:3,el))))
!!PS         nu_loc = maxval(ulevels_nod2D(elem2D_nodes(1:3,el)))
        ule = max( uln, ulevels(el) )
        nle = min( nln, min(nlev,nlevels(el)) )
        !!PS DO nz=1, nl_loc
        DO nz=ule, nle
           vol(nz,n) = vol(nz,n) + elem_area(el)
           work_array(nz,n) = work_array(nz,n) + elem_area(el) * (arr(nz, elem2D_nodes(1,el)) &
                                                                + arr(nz, elem2D_nodes(2,el)) &
                                                                + arr(nz, elem2D_nodes(3,el)))
        END DO
     ENDDO
     !!PS DO nz=1,nlevels_nod2d(n)
     DO nz=uln,nln
        vol(nz,n) = 1._WP / (3._WP * vol(nz,n))  ! Here, we need the inverse and scale by 1/3
     END DO
  END DO

  ! combined: scale by patch volume + copy back to original field
  DO n=1, myDim_nod2D
     !!PS DO nz=1, min(nlev, nlevels_nod2d(n))
     uln = ulevels_nod2d(n)
     nln = min(nlev,nlevels_nod2d(n))
     DO nz=uln,nln
        arr(nz, n) = work_array(nz, n) *vol(nz,n) 
!!PS         if (arr(nz,n)/=arr(nz,n)) then
!!PS             write(*,*) ' --> found NaN in smoothing'
!!PS             write(*,*) ' mype = ', mype
!!PS             write(*,*) ' n    = ', n
!!PS             write(*,*) ' nz,uln,nln      = ', nz,uln,nln
!!PS             write(*,*) ' arr(nz,n)       = ', arr(nz,n)
!!PS             write(*,*) ' work_array(nz,n)= ', work_array(nz,n)
!!PS             write(*,*) ' vol(nz,n)       = ', vol(nz,n)
!!PS         endif 
     END DO
  end DO
  
  call exchange_nod(arr)

! And the remaining smoothing sweeps

  DO q=1,N_smooth-1
     DO n=1, myDim_nod2D
        uln = ulevels_nod2d(n)
        nln = min(nlev,nlevels_nod2d(n))
        !!PS work_array(1:min(nlev, nlevels_nod2d(n)),n) = 0._WP
        work_array(1:nln,n) = 0._WP
        DO j=1,nod_in_elem2D_num(n)
           el = nod_in_elem2D(j,n)
           !!PS nl_loc = min(nlev, minval(nlevels_nod2d(elem2D_nodes(1:3,el))))
           !!PS nu_loc = maxval(ulevels_nod2D(elem2D_nodes(1:3,el)))
           !!PS DO nz=1, 
           ule = max( uln, ulevels(el) )
           nle = min( nln, min(nlev,nlevels(el)) )
           DO nz=ule,nle
              work_array(nz,n) = work_array(nz,n) + elem_area(el) * (arr(nz, elem2D_nodes(1,el)) &
                                                                   + arr(nz, elem2D_nodes(2,el)) &
                                                                   + arr(nz, elem2D_nodes(3,el)))
           END DO
        ENDDO
     ENDDO
! combined: scale by patch volume + copy back to original field
     DO n=1, myDim_nod2D
        !!PS DO nz=1, min(nlev, nlevels_nod2d(n))
        uln = ulevels_nod2d(n)
        nln = min(nlev,nlevels_nod2d(n))
        DO nz=uln,nln
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
subroutine smooth_elem2D(arr, N, mesh)
    IMPLICIT NONE
    type(t_mesh), intent(in)                  , target :: mesh
    integer, intent(in)                        :: N
    real(KIND=WP), dimension(:), intent(inout) :: arr
    integer                                    :: node, elem, j, q, elnodes(3)
    real(kind=WP)                              :: vol
#include "associate_mesh.h"    
    allocate(work_array(myDim_nod2D+eDim_nod2D))
    DO q=1, N !apply mass matrix N times to smooth the field
        DO node=1, myDim_nod2D
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
        call exchange_nod(work_array)
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
subroutine smooth_elem3D(arr, N, mesh)
    IMPLICIT NONE
    type(t_mesh), intent(in)                    , target :: mesh
    integer, intent(in)                          :: N
    real(KIND=WP), dimension(:,:), intent(inout) :: arr
    integer                                      :: node, elem, my_nl, nz, j, q, elnodes(3)
    real(kind=WP)                                :: vol
#include "associate_mesh.h"  

    allocate(work_array(myDim_nod2D+eDim_nod2D))
    
    my_nl=ubound(arr,1)
    DO q=1, N !apply mass matrix N times to smooth the field
        DO nz=1, my_nl
            work_array = 0.0_WP
            DO node=1, myDim_nod2D
                vol=0._WP
                if (nz > nlevels_nod2d(node)) CYCLE
                if (nz < ulevels_nod2D(node)    ) CYCLE
                DO j=1, nod_in_elem2D_num(node)
                    elem=nod_in_elem2D(j, node)
                    if (elem<=0) CYCLE
                    if (nz>nlevels(elem)  ) CYCLE
                    if (nz<ulevels(elem)) CYCLE
                    elnodes=elem2D_nodes(:,elem)
                    work_array(node)=work_array(node)+arr(nz, elem)*elem_area(elem)
                    vol=vol+elem_area(elem)
                END DO
                work_array(node)=work_array(node)/vol
            END DO
            call exchange_nod(work_array)
            DO elem=1, myDim_elem2D
                if (nz>nlevels(elem)  ) CYCLE
                if (nz<ulevels(elem)) CYCLE
                elnodes=elem2D_nodes(:, elem)
                arr(nz, elem)=sum(work_array(elnodes))/3.0_WP
            ENDDO
        END DO
        call exchange_elem(arr)
    END DO
    deallocate(work_array)

end subroutine smooth_elem3D
!
!--------------------------------------------------------------------------------------------
!
subroutine integrate_nod_2D(data, int2D, mesh)
  use o_MESH
  use g_PARSUP
  use g_comm_auto

  IMPLICIT NONE
  type(t_mesh),  intent(in)      , target :: mesh
  real(kind=WP), intent(in)       :: data(:)
  real(kind=WP), intent(inout)    :: int2D

  integer       :: row
  real(kind=WP) :: lval
#include "associate_mesh.h"
  lval=0.0_WP
  do row=1, myDim_nod2D
     !!PS lval=lval+data(row)*area(1,row)
     lval=lval+data(row)*areasvol(ulevels_nod2D(row),row)
  end do

  int2D=0.0_WP
  call MPI_AllREDUCE(lval, int2D, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
end subroutine integrate_nod_2D
!
!--------------------------------------------------------------------------------------------
!
subroutine integrate_nod_3D(data, int3D, mesh)
  use g_PARSUP
  use g_comm_auto

  IMPLICIT NONE
  type(t_mesh),  intent(in)      , target :: mesh
  real(kind=WP), intent(in)       :: data(:,:)
  real(kind=WP), intent(inout)    :: int3D

  integer       :: k, row
  real(kind=WP) :: lval
#include "associate_mesh.h"

  lval=0.0_WP
  do row=1, myDim_nod2D
     !!PS do k=1, nlevels_nod2D(row)-1
     do k=ulevels_nod2D(row), nlevels_nod2D(row)-1
        lval=lval+data(k, row)*areasvol(k,row)*hnode_new(k,row)  ! --> TEST_cavity
     end do
  end do
  int3D=0.0_WP
  call MPI_AllREDUCE(lval, int3D, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
end subroutine integrate_nod_3D
!
!--------------------------------------------------------------------------------------------
!
subroutine extrap_nod3D(arr, mesh)
    IMPLICIT NONE
    type(t_mesh),  intent(in)     , target :: mesh
    real(KIND=WP), intent(inout)   :: arr(:,:)
    integer                        :: n, nl1, nz, k, j, el, cnt, jj
    real(kind=WP), allocatable     :: work_array(:)
    real(kind=WP)                  :: val
    integer                        :: enodes(3)
    logical                        :: success
    real(kind=WP)                  :: loc_max, glob_max
    integer                        :: loc_sum, glob_sum, glob_sum_old

#include "associate_mesh.h"
    !___________________________________________________________________________
    allocate(work_array(myDim_nod2D+eDim_nod2D))
    call exchange_nod(arr)
    
    !___________________________________________________________________________
    loc_max=maxval(arr(1,:))
    glob_max=0._WP
    call MPI_AllREDUCE(loc_max, glob_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    glob_sum=-1

    !___________________________________________________________________________
    do while (glob_max>0.99_WP*dummy)  
        !_______________________________________________________________________
        ! extrapolate in horizontal direction
        do nz=1, nl-1
            work_array=arr(nz,:)
            success=.true.
            
            !___________________________________________________________________
            ! horizontal extrapolation 
            do while (success) ! --> do while runs as long as success==.true.
                success=.false.
                
                !_______________________________________________________________
                ! loop over local vertices n 
                do n=1, myDim_nod2D+eDim_nod2D
                    ! found node n that has to be extrapolated
                    if ( (work_array(n)>0.99_WP*dummy) .and.  (nlevels_nod2D(n)>nz)) then
                        cnt=0
                        val=0._WP
                        
                        !_______________________________________________________
                        ! loop over adjacent elements
                        do k=1, nod_in_elem2D_num(n)
                            el=nod_in_elem2D(k, n)
                            
                            if (nz>nlevels(el)) cycle
                            enodes=elem2D_nodes(:, el)
                            !___________________________________________________
                            ! loop over vertices of adjacent element 
                            do j=1, 3
                                if (enodes(j)==0) cycle
                                if ((work_array(enodes(j))<0.99_WP*dummy) .and. (nlevels_nod2D(enodes(j))>nz)) then
                                    val=val+work_array(enodes(j))
                                    cnt=cnt+1              
                                end if
                            end do ! --> do j=1, 3
                        end do ! --> do k=1, nod_in_elem2D_num(n)
                        
                        !_______________________________________________________
                        ! found valid neighbouring node that can be used for 
                        ! extrapooolation 
                        if (cnt>0) then
                            work_array(n)=val/real(cnt,WP)
                            success=.true. ! --> reason to stay in while loop
                        end if
                        
                    end if ! --> if ( (work_array(n)>0.99_WP*dummy) ....
                    
                end do ! --> do n=1, myDim_nod2D
                
            end do ! --> do while (success)
            arr(nz,:)=work_array
            
        end do ! --> do nz=1, nl-1
        
        !_______________________________________________________________________
        call exchange_nod(arr)
        
        !_______________________________________________________________________
        loc_max=maxval(arr(1,:))
        call MPI_AllREDUCE(loc_max, glob_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)   
        
    END DO ! -->  DO WHILE (glob_max>0.99_WP*dummy)  
    
    !___________________________________________________________________________
    ! extrapolate in vertical direction
    do n=1, myDim_nod2D
        nl1 = nlevels_nod2D(n)-1
        !!PS DO nz=2, nl1
        !!PS  do nz=ulevels_nod2D(n)+1, nl1
        do nz=2, nl1
            if (arr(nz,n)>0.99_WP*dummy) arr(nz,n)=arr(nz-1,n)
        end do
    end do
    call exchange_nod(arr)
    
    !___________________________________________________________________________
    deallocate(work_array)
    
end subroutine extrap_nod3D
!
!--------------------------------------------------------------------------------------------
!
end module g_support


