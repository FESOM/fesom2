!a set of auxuary routines for: 
!1. smoothing FESOM fields using mass matrix
!2. computing surface integrals of the FESOM fields
module g_support
  USE MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_comm_auto
  use o_ARRAYS
  use g_config, only: dummy
  implicit none

  private
  public :: smooth_nod, smooth_elem, integrate_nod, integrate_elem, extrap_nod, omp_min_max_sum1, omp_min_max_sum2
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
! computes 2D integral of a nodal field
  INTERFACE integrate_elem
            MODULE PROCEDURE integrate_elem_2D, integrate_elem_3D
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
subroutine smooth_nod2D(arr, N, partit, mesh)
  IMPLICIT NONE
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer, intent(in)                        :: N
  real(KIND=WP), dimension(:), intent(inout) :: arr
  integer                                    :: node, elem, j, q, elnodes(3)
  real(kind=WP)                              :: vol

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  allocate(work_array(myDim_nod2D))
  DO q=1, N !apply mass matrix N times to smooth the field
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, elem, j, q, elnodes, vol)
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
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    DO node=1,myDim_nod2D
       arr(node)=work_array(node)
    ENDDO
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP MASTER
    call exchange_nod(arr, partit)
!$OMP END MASTER
!$OMP BARRIER
  END DO
  deallocate(work_array)
end subroutine smooth_nod2D
!
!--------------------------------------------------------------------------------------------
!
subroutine smooth_nod3D(arr, N_smooth, partit, mesh)
  IMPLICIT NONE
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit

  integer, intent(in)            :: N_smooth
  real(KIND=WP), intent(inout)   :: arr(:,:)

  integer                        :: n, q, el, nz, j, nlev
  integer                        :: uln, nln, ule, nle
  real(kind=WP)                  :: vol(mesh%nl, partit%myDim_nod2D)
  real(kind=WP), allocatable     :: work_array(:,:)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  nlev=ubound(arr,1)
  allocate(work_array(nlev,myDim_nod2D))

! Precompute area of patches on all levels (at the bottom, some neighbouring
! nodes may vanish in the bathymetry) in the first smoothing step
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, q, el, nz, j, uln, nln, ule, nle)
!$OMP DO
  DO n=1, myDim_nod2D
     uln = ulevels_nod2d(n)
     nln = min(nlev,nlevels_nod2d(n))
     vol(       1:nln,n) = 0._WP
     work_array(1:nln,n) = 0._WP
     DO j=1, nod_in_elem2D_num(n)
        el = nod_in_elem2D(j,n)
        ule = max( uln, ulevels(el) )
        nle = min( nln, min(nlev,nlevels(el)) )
        DO nz=ule, nle
           vol(nz,n) = vol(nz,n) + elem_area(el)
           work_array(nz,n) = work_array(nz,n) + elem_area(el) * (arr(nz, elem2D_nodes(1,el)) &
                                                                + arr(nz, elem2D_nodes(2,el)) &
                                                                + arr(nz, elem2D_nodes(3,el)))
        END DO
     ENDDO
     DO nz=uln,nln
        vol(nz,n) = 1._WP / (3._WP * vol(nz,n))  ! Here, we need the inverse and scale by 1/3
     END DO
  END DO
!$OMP END DO
  ! combined: scale by patch volume + copy back to original field
!$OMP DO
  DO n=1, myDim_nod2D
     uln = ulevels_nod2d(n)
     nln = min(nlev,nlevels_nod2d(n))
     DO nz=uln,nln
        arr(nz, n) = work_array(nz, n) *vol(nz,n) 
     END DO
  END DO
!$OMP END DO  
!$OMP BARRIER
!$OMP MASTER
  call exchange_nod(arr, partit)
!$OMP END MASTER
!$OMP BARRIER
! And the remaining smoothing sweeps
  DO q=1,N_smooth-1
!$OMP DO
     DO n=1, myDim_nod2D
        uln = ulevels_nod2d(n)
        nln = min(nlev,nlevels_nod2d(n))
        work_array(1:nln,n) = 0._WP
        DO j=1,nod_in_elem2D_num(n)
           el = nod_in_elem2D(j,n)
           ule = max( uln, ulevels(el) )
           nle = min( nln, min(nlev,nlevels(el)) )
           DO nz=ule,nle
              work_array(nz,n) = work_array(nz,n) + elem_area(el) * (arr(nz, elem2D_nodes(1,el)) &
                                                                   + arr(nz, elem2D_nodes(2,el)) &
                                                                   + arr(nz, elem2D_nodes(3,el)))
           END DO
        ENDDO
     ENDDO
!$OMP END DO
! combined: scale by patch volume + copy back to original field
!$OMP DO
     DO n=1, myDim_nod2D
        !!PS DO nz=1, min(nlev, nlevels_nod2d(n))
        uln = ulevels_nod2d(n)
        nln = min(nlev,nlevels_nod2d(n))
        DO nz=uln,nln
           arr(nz, n) = work_array(nz, n) *vol(nz,n) 
        END DO
     END DO
!$OMP END DO
!$OMP BARRIER
!$OMP MASTER
     call exchange_nod(arr, partit)
!$OMP END MASTER
!$OMP BARRIER
  enddo
!$OMP END PARALLEL
  deallocate(work_array)
end subroutine smooth_nod3D

!
!--------------------------------------------------------------------------------------------
!
subroutine smooth_elem2D(arr, N, partit, mesh)
    IMPLICIT NONE
    type(t_mesh),   intent(in),    target         :: mesh
    type(t_partit), intent(inout), target         :: partit
    integer, intent(in)                        :: N
    real(KIND=WP), dimension(:), intent(inout) :: arr
    integer                                    :: node, elem, j, q, elnodes(3)
    real(kind=WP)                              :: vol
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    allocate(work_array(myDim_nod2D+eDim_nod2D))
    DO q=1, N !apply mass matrix N times to smooth the field
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, elem, j, q, elnodes, vol)
!$OMP DO
        DO node=1, myDim_nod2D
            vol=0._WP
            work_array(node)=0._WP
            DO j=1, nod_in_elem2D_num(node)
                elem=nod_in_elem2D(j, node)
                work_array(node)=work_array(node)+arr(elem)*elem_area(elem)
                vol=vol+elem_area(elem)
            END DO
            work_array(node)=work_array(node)/vol
        END DO
!$OMP END DO
!$OMP BARRIER
!$OMP MASTER
        call exchange_nod(work_array, partit)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
        DO elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:, elem)
            arr(elem)=sum(work_array(elnodes))/3.0_WP  ! Here, we need the inverse and scale by 1/3
        ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP MASTER
        call exchange_elem(arr, partit)
!$OMP END MASTER
!$OMP BARRIER
!$OMP END PARALLEL
    END DO
    deallocate(work_array)
end subroutine smooth_elem2D
!
!--------------------------------------------------------------------------------------------
!
subroutine smooth_elem3D(arr, N, partit, mesh)
    IMPLICIT NONE
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    integer, intent(in)                          :: N
    real(KIND=WP), dimension(:,:), intent(inout) :: arr
    integer                                      :: node, elem, my_nl, nz, j, q, elnodes(3)
    real(kind=WP)                                :: vol
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 

    allocate(work_array(myDim_nod2D+eDim_nod2D))
    
    my_nl=ubound(arr,1)
    DO q=1, N !apply mass matrix N times to smooth the field
        DO nz=1, my_nl
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, elem, j, q, elnodes, vol)
!$OMP DO
            DO node=1, myDim_nod2D+eDim_nod2D
               work_array(node) = 0.0_WP
            END DO
!$OMP END DO
!$OMP DO
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
!$OMP END DO
!$OMP BARRIER
!$OMP MASTER
            call exchange_nod(work_array, partit)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
            DO elem=1, myDim_elem2D
                if (nz>nlevels(elem)  ) CYCLE
                if (nz<ulevels(elem)) CYCLE
                elnodes=elem2D_nodes(:, elem)
                arr(nz, elem)=sum(work_array(elnodes))/3.0_WP
            ENDDO
!$OMP END DO
!$OMP END PARALLEL
        END DO
        call exchange_elem(arr, partit)
!$OMP BARRIER
    END DO
    deallocate(work_array)

end subroutine smooth_elem3D
!
!--------------------------------------------------------------------------------------------
!
subroutine integrate_nod_2D(data, int2D, partit, mesh)
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_comm_auto

  IMPLICIT NONE
  type(t_mesh),  intent(in),    target :: mesh
  type(t_partit),intent(inout), target :: partit
  real(kind=WP), intent(in)         :: data(:)
  real(kind=WP), intent(inout)      :: int2D

  integer       :: row
  real(kind=WP) :: lval
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 

lval=0.0_WP
#if !defined(__openmp_reproducible)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(row)
!$OMP DO REDUCTION (+: lval)
#endif
  do row=1, myDim_nod2D
     lval=lval+data(row)*areasvol(ulevels_nod2D(row),row)
  end do
#if !defined(__openmp_reproducible) 
!$OMP END DO
!$OMP END PARALLEL
#endif
  int2D=0.0_WP
  call MPI_AllREDUCE(lval, int2D, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
end subroutine integrate_nod_2D
!
!--------------------------------------------------------------------------------------------
!
subroutine integrate_nod_3D(data, int3D, partit, mesh)
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_comm_auto

  IMPLICIT NONE
  type(t_mesh),  intent(in), target :: mesh
  type(t_partit),intent(in), target :: partit
  real(kind=WP), intent(in)       :: data(:,:)
  real(kind=WP), intent(inout)    :: int3D

  integer       :: k, row
  real(kind=WP) :: lval
  real(kind=WP) :: lval_row


#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 

  lval=0.0_WP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row, k, lval_row) REDUCTION(+: lval)
  do row=1, myDim_nod2D
     lval_row = 0.
     do k=ulevels_nod2D(row), nlevels_nod2D(row)-1
        lval_row=lval_row+data(k, row)*areasvol(k,row)*hnode_new(k,row)  ! --> TEST_cavity
     end do
#if defined(__openmp_reproducible)
!$OMP ORDERED
#endif
     lval = lval + lval_row
#if defined(__openmp_reproducible)
!$OMP END ORDERED
#endif
  end do
!$OMP END PARALLEL DO

  int3D=0.0_WP
  call MPI_AllREDUCE(lval, int3D, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
end subroutine integrate_nod_3D
!
!--------------------------------------------------------------------------------------------
!
subroutine extrap_nod3D(arr, partit, mesh)
    IMPLICIT NONE
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(KIND=WP),  intent(inout)  :: arr(:,:)
    integer                        :: n, nl1, nz, k, j, el, cnt, jj
    real(kind=WP), allocatable     :: work_array(:)
    real(kind=WP)                  :: val
    integer                        :: enodes(3)
    logical                        :: success
    real(kind=WP)                  :: loc_max, glob_max
    integer                        :: loc_sum, glob_sum, glob_sum_old

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    !___________________________________________________________________________
    allocate(work_array(myDim_nod2D+eDim_nod2D))
    call exchange_nod(arr, partit)
    
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
        call exchange_nod(arr, partit)
        
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
    call exchange_nod(arr, partit)
    
    !___________________________________________________________________________
    deallocate(work_array)
    
end subroutine extrap_nod3D
!
!--------------------------------------------------------------------------------------------
! returns min/max/sum of a one dimentional array (same as minval) but with the support of OpenMP
FUNCTION omp_min_max_sum1(arr, pos1, pos2, what, partit, nan) 
  USE MOD_PARTIT
  implicit none
  real(kind=WP), intent(in)   :: arr(:)
  integer,       intent(in)   :: pos1, pos2
  character(3),  intent(in)   :: what
  real(kind=WP), optional     :: nan !to be implemented upon the need (for masked arrays)
  real(kind=WP)               :: omp_min_max_sum1
  real(kind=WP)               :: val
  integer                     :: n

  type(t_partit),intent(in), &
                       target :: partit

  SELECT CASE (trim(what))
    CASE ('sum')
       val=0.0_WP
#if !defined(__openmp_reproducible)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n)
!$OMP DO REDUCTION(+: val)
#endif
       do n=pos1, pos2
          val=val+arr(n)
       end do
#if !defined(__openmp_reproducible)
!$OMP END DO
!$OMP END PARALLEL
#endif

    CASE ('min')
       val=arr(1)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n)
!$OMP DO REDUCTION(min: val)
       do n=pos1, pos2
          val=min(val, arr(n))
       end do
!$OMP END DO
!$OMP END PARALLEL

    CASE ('max')
       val=arr(1)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n)
!$OMP DO REDUCTION(max: val)
       do n=pos1, pos2
          val=max(val, arr(n))
       end do
!$OMP END DO
!$OMP END PARALLEL

    CASE DEFAULT
       if (partit%mype==0) write(*,*) trim(what), ' is not implemented in omp_min_max_sum case!'
       call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
       STOP
  END SELECT

  omp_min_max_sum1=val
END FUNCTION
!
!--------------------------------------------------------------------------------------------
! returns min/max/sum of a two dimentional array (same as minval) but with the support of OpenMP
FUNCTION omp_min_max_sum2(arr, pos11, pos12, pos21, pos22, what, partit, nan) 
  implicit none
  real(kind=WP), intent(in)   :: arr(:,:)
  integer,       intent(in)   :: pos11, pos12, pos21, pos22
  character(3),  intent(in)   :: what
  real(kind=WP), optional     :: nan !to be implemented upon the need (for masked arrays)
  real(kind=WP)               :: omp_min_max_sum2
  real(kind=WP)               :: val, vmasked, val_part(pos11:pos12)
  integer                     :: i, j
  

  type(t_partit),intent(in), &
                       target :: partit
  
  IF (PRESENT(nan)) vmasked=nan
  
  SELECT CASE (trim(what))
    CASE ('min')
      if (.not. present(nan)) vmasked=huge(vmasked) !just some crazy number
      val=arr(1,1)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j)
!$OMP DO REDUCTION(min: val)
      do j=pos21, pos22
         do i=pos11, pos12
            if (arr(i,j)/=vmasked) val=min(val, arr(i,j))
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL

    CASE ('max')
      if (.not. present(nan)) vmasked=tiny(vmasked) !just some crazy number
      val=arr(1,1)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j)
!$OMP DO REDUCTION(max: val)
      do j=pos21, pos22
         do i=pos11, pos12
            if (arr(i,j)/=vmasked) val=max(val, arr(i,j))
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL

    CASE ('sum')
      if (.not. present(nan)) vmasked=huge(vmasked) !just some crazy number
      val=0.
#if  !defined(__openmp_reproducible)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j) REDUCTION(+: val)
      do j=pos21, pos22
         do i=pos11, pos12
            if (arr(i,j)/=vmasked) val=val+arr(i,j)
         end do
      end do
!$OMP END PARALLEL DO
#else
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) 
      do j=pos21, pos22
         val_part(j) = sum(arr(pos11:pos12,j), mask=(arr(pos11:pos12,j)/=vmasked))
      end do
!$OMP END PARALLEL DO 
      val = sum(val_part(pos21:pos22))
#endif

   CASE DEFAULT
      if (partit%mype==0) write(*,*) trim(what), ' is not implemented in omp_min_max_sum case!'
      call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
      STOP
   END SELECT

omp_min_max_sum2=val
END FUNCTION
!
!--------------------------------------------------------------------------------------------
!
subroutine integrate_elem_3D(data, int3D, partit, mesh)
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_comm_auto

  IMPLICIT NONE
  type(t_mesh),  intent(in), target :: mesh
  type(t_partit),intent(in), target :: partit
  real(kind=WP), intent(in)       :: data(:,:)
  real(kind=WP), intent(inout)    :: int3D

  integer       :: k, row
  real(kind=WP) :: lval
  real(kind=WP) :: lval_row

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 

  lval=0.0_WP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row, k, lval_row) REDUCTION(+: lval)
  do row=1, myDim_elem2D
     if(elem2D_nodes(1, row) > myDim_nod2D) cycle
     lval_row = 0.
     do k=ulevels(row), nlevels(row)-1
        lval_row=lval_row+data(k, row)*elem_area(row)*helem(k,row)
     end do
#if defined(__openmp_reproducible)
!$OMP ORDERED
#endif
     lval = lval + lval_row
#if defined(__openmp_reproducible)
!$OMP END ORDERED
#endif
  end do
!$OMP END PARALLEL DO

  int3D=0.0_WP
  call MPI_AllREDUCE(lval, int3D, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
end subroutine integrate_elem_3D
!
!--------------------------------------------------------------------------------------------
!
subroutine integrate_elem_2D(data, int2D, partit, mesh)
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_comm_auto

  IMPLICIT NONE
  type(t_mesh),  intent(in), target :: mesh
  type(t_partit),intent(in), target :: partit
  real(kind=WP), intent(in)         :: data(:)
  real(kind=WP), intent(inout)      :: int2D

  integer       :: row
  real(kind=WP) :: lval

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 

  lval=0.0_WP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) REDUCTION(+: lval)
  do row=1, myDim_elem2D
     if(elem2D_nodes(1, row) > myDim_nod2D) cycle
#if defined(__openmp_reproducible)
!$OMP ORDERED
#endif
     lval = lval + data(row)*elem_area(row)
#if defined(__openmp_reproducible)
!$OMP END ORDERED
#endif
  end do
!$OMP END PARALLEL DO

  int2D=0.0_WP
  call MPI_AllREDUCE(lval, int2D, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
end subroutine integrate_elem_2D
end module g_support


