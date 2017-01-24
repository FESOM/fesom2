      program main
      implicit none

#include "mpif.h"
#include "fparms.h"
#include "faux.h"      

      integer :: ierr, npro, myid
      integer :: mprocx, mprocy, mprocz, xnmesh, ynmesh
      integer :: nx, ny, nz, n, size, nloc, iov, iout
      integer :: nrow, start, pos, nnzrow
      integer, dimension(:), allocatable :: mapptr, maptmp, ja, ia
#if defined(DBL_CMPLX)      
      double complex, dimension(:), allocatable :: a
      double complex, dimension(100) :: stencil
#else
      double precision, dimension(:), allocatable :: a
      real*8, dimension(100) :: stencil
#endif      
      real*8 :: tpc, ttol, res0, res1, ratio
      integer :: nnz, i, j, its, len
      logical :: cnpro
      character(len=100) :: pcname, pciluname
      fprm :: prm
!!$variables related to pARMS
      parms_Map :: map
      parms_Mat :: mat
      parms_PC  :: pc
      parms_Timer :: t
      parms_Solver :: ksp
#if defined(DBL_CMPLX)      
      double complex, dimension(:), allocatable :: sol, rhs, y
#else
      double precision, dimension(:), allocatable :: sol, rhs, y
#endif      

      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
      call mpi_comm_size(MPI_COMM_WORLD, npro, ierr)

      call fread_param('inputs', prm, mprocx, mprocy, xnmesh, ynmesh)
      if (myid .eq. 0) then
            print *, 'mprocx = ', mprocx, 'mprocy = ', mprocy
            print *, 'xnmesh = ', xnmesh, 'ynmesh = ', ynmesh
      end if

      i = mprocx * mprocy 

      if (npro .ne. i) then
            if (myid .eq. 0) then
                  print *, 'ERROR: npro should be =', mprocx*mprocy
            end if
            call mpi_barrier(MPI_COMM_WORLD, ierr)
            call mpi_abort(MPI_COMM_WORLD, 98, ierr)
      end if

      nx = mprocx * xnmesh
      ny = mprocy * ynmesh
      nz = 1
      size = xnmesh*ynmesh*mprocx*mprocy


      if (.not. allocated(mapptr)) then
            allocate(mapptr(npro+1), stat=ierr)
      end if
      if (.not. allocated(maptmp)) then
            allocate(maptmp(size), stat=ierr)
      end if

      iov = 0
      mprocz = 1
   
      call part2(nx, ny, nz, mprocx, mprocy, mprocz, iov, maptmp, 
     & mapptr, iout)
      n = nx * ny * nz

!!$ Create map object and setup mapping from global labels to processors 
      call parms_mapcreatefromptr(map, n, maptmp, mapptr, 
     & MPI_COMM_WORLD, 1, NONINTERLACED, ierr)
! get local size
      call parms_mapgetlocalsize(map, nloc)


!! Allocate memory for distributed vectors
      allocate(rhs(nloc),    stat=ierr)
      allocate(sol(nloc),    stat=ierr)
      allocate(y(nloc),   stat=ierr)

      nnz = 10*xnmesh*ynmesh

      if (.not. allocated(a)) then
            allocate(a(nnz), stat=ierr)
      end if
      if (.not. allocated(ja)) then
            allocate(ja(nnz), stat=ierr)
      end if
      if (.not. allocated(ia)) then
            allocate(ia(nloc+1), stat=ierr)
      end if
      nrow  = mapptr(myid+2) - mapptr(myid+1)
      start = mapptr(myid+1)

! Generate the local matrix 
#if defined(DBL_CMPLX)
      call zgen5loc(nx, ny, nz, nloc, maptmp(start), a, ja, ia, stencil)  
#else
      call gen5loc(nx, ny, nz, nloc, maptmp(start), a, ja, ia, stencil)  
#endif 

! Create a matrix object based on map
      call parms_matcreate(mat, map, ierr) 

! Insert entries into the matrix
      call parms_matsetvalues(mat, nrow, maptmp(start), ia, ja, a, 
     & INSERT, ierr) 

! Free the matrix stored in CSR format (ia, ja, a)
      deallocate(ia)
      deallocate(ja)
      deallocate(a)

! Setup the matrix A
      call parms_matsetup(mat, ierr)
! Create a timer
      call parms_timercreate(t, ierr)

! Create a preconditioner object
      call parms_pccreate(pc, mat, ierr)
! Set parameters for pc
      call fset_pc_params(pc, prm)
! Reset the timer
      call parms_timerreset(t, ierr)
! Setup the preconditioner
      call parms_pcsetup(pc, ierr)

      call parms_timerget(t, tpc, ierr)

! Pause the timer
      call parms_timerpause(t, ierr)

! Create a krylov subspace object
      call parms_solvercreate(ksp, mat, pc, ierr)
! Set the type of solver ksp.
      call parms_solversettype(ksp, SOLFGMRES, ierr)
! Set parameters for the solver
      call fset_solver_params(ksp, prm)

      call fprm_free(prm)

! Set up the artifical right-hand-side vector
      do i=1,nloc
            sol(i) = 1.d0
      enddo
      call parms_matvec(mat, sol, rhs, ierr)

! Set initial guess 
      do i=1,nloc
            sol(i) = 0.d0
      enddo

! Get initial residual norm
      call parms_matvec(mat, sol, y, ierr)
      call parms_vecaxpy(y, rhs, (-1.0d0,0.d0), map, ierr)
      call parms_vecgetnorm2(y, res0, map, ierr)

      call parms_timerrestart(t, ierr)
      call parms_solverapply(ksp, rhs, sol, ierr)
      call parms_timerget(t, ttol, ierr)

! Get the residual error
      call parms_matmvpy(mat,(-1.0d0,0.d0),sol,(1.0d0,0.d0),rhs,y,ierr)
      call parms_vecgetnorm2(y, res1, map, ierr)
      call parms_pcgetratio(pc,ratio, ierr)
      if (myid .eq. 0) then
            call parms_solvergetits(ksp, its, ierr)
            call parms_pcgetname(pc, pcname, len, ierr)
            write(*, 60) pcname(1:len)
            call parms_pcilugetname(pc, pciluname, len, ierr)
            write(*, 66) pciluname(1:len)
            write(*, 78) ratio
            write(*, 98) npro
            write(*, 100) its
            write(*, 90)  tpc
            write(*, 200) ttol - tpc
            write(*, 300) ttol
            write(*, 400) res0
            write(*, 500) res1
   60       format('The preconditioner', 12x, A)
   66       format('The local preconditioner', 6x, A)
   78       format('The memory usage ', 11x, '=', f5.2)
   98       format('The number of iteration', 3x, '=', I4.1)
  100       format('The number of iteration', 5x, '=', I4.1)
   90       format('The time for creating pc ', 3x, '=', es9.2,'s')
  200       format('The solving time ', 11x, '=', es9.2,'s')
  300       format('The total time ', 13x, '=', es9.2,'s')
  400       format('The initial residual error  =', es9.2)
  500       format('The solution residual error =', es9.2)
      end if

! free memories for map, vec, mat, solver 
      deallocate(sol, stat=ierr)
      deallocate(rhs,stat=ierr)
      deallocate(y,stat=ierr)
      deallocate(mapptr, maptmp)
      call parms_solverfree(ksp, ierr)
      call parms_mapfree(map, ierr)
      call parms_matfree(mat, ierr)
      call parms_pcfree(pc, ierr)
      call parms_timerfree(t, ierr)

      call mpi_finalize(ierr)

      end program main

