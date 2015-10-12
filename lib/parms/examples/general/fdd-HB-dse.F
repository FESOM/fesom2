!!$ *                        Program dd-HB-dse
!!$ *--------------------------------------------------------------------
!!$ *
!!$ *  In this test program, each processor reads the whole matrix
!!$ *  from file. The matrix is assumed to be in Harwell-Boeing format.
!!$ *  Matrix graph is then  partitioned  using  DSE, a simple partitioning
!!$ *  routine, and scatters the local matrices to each processor. Once
!!$ *  these submatrices are received each processor solves the problem
!!$ *  using preconditioned FGMRES preconditioned with :
!!$ *                         BJ, RAS, SCHUR
!!$ *--------------------------------------------------------------------
      program main

      implicit none
#include "mpif.h"
#include "fparms.h"
#include "faux.h"

      integer, allocatable, dimension(:) :: ja, ia, iwk, riord, im
      integer, dimension(:), allocatable :: dom, idom, mask, jwk, link
      integer, dimension(:), allocatable :: jb, ib, jc, ic
#if defined(DBL_CMPLX)      
      double complex, dimension(:), allocatable :: a, b, mc, rhstmp
#else
      double precision, dimension(:), allocatable :: a, b, mc, rhstmp
#endif
      double precision :: tpc, ttol, res0, res1, ratio
      character(len=100) :: matrix
      character(len=2) :: guesol
      character(len=72) :: title
      character(len=8) :: key
      character(len=3) :: type
      integer :: myid, npro, nmax, nzmax, job, nrhs, nnz, usmy
      integer :: ierr, i, j, its, n, ncol, i1, len, nloc, matlen
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

!!$  Read matrix; either using user-defined function (routine uread) or 
!!$              SPARSKIT function for reading Harwell-Boeieng matrices

!!$  read matrix name from file input 
      call fread_param('inputs', prm, matrix, matlen)

!      print*,'len = ', len, 'matlen = ', matlen
!      print*,'matrix name = ',matrix
!      stop
!!$ variable "matrix" stores the name of the file in HB format 
!!$ Read a Harwell-Boeing matrix using wreadmtc from
!!$ SPARSKIT - call wreadmtc a first time to determine sizes
!!$ of arrays. read in values on the second call.
      job = 0
#if defined(DBL_CMPLX)      
!      call zwreadmtc(nmax, nzmax, job, matrix,matlen, a, ja, ia, 
!     & rhstmp, nrhs, guesol, n, ncol, nnz, title, key, type, ierr)
      call zreadmtc(nmax, nzmax, job, matrix, a, ja, ia, 
     & rhstmp, nrhs, guesol, n, ncol, nnz, title, key, type, ierr)     
#else
!      call wreadmtc(nmax, nzmax, job, matrix,matlen, a, ja, ia, 
!     & rhstmp, nrhs, guesol, n, ncol, nnz, title, key, type, ierr)
      call readmtc(nmax, nzmax, job, matrix, a, ja, ia, 
     & rhstmp, nrhs, guesol, n, ncol, nnz, title, key, type, ierr)     
#endif  
      allocate(ia(n+1),    stat=ierr)
      allocate(ja(nnz),    stat=ierr)
      allocate(a(nnz),     stat=ierr)
      allocate(rhstmp(n),  stat=ierr) 

!!$ Array sizes determined. Now call wreadmtc again for 
!!$ really reading
      nmax = n
      nzmax = nnz
      job = 3
#if defined(DBL_CMPLX)      
!      call zwreadmtc(n, nnz, job, matrix,matlen, a, ja, ia, rhstmp, 
!     &  nrhs, guesol, n, ncol, nnz, title, key, type, ierr)
      call zreadmtc(n, nnz, job, matrix, a, ja, ia, rhstmp, 
     &  nrhs, guesol, n, ncol, nnz, title, key, type, ierr)     
#else
!      call wreadmtc(n, nnz, job, matrix,matlen, a, ja, ia, rhstmp, 
!     &  nrhs, guesol, n, ncol, nnz, title, key, type, ierr)
      call readmtc(n, nnz, job, matrix, a, ja, ia, rhstmp, 
     &  nrhs, guesol, n, ncol, nnz, title, key, type, ierr)     
#endif

      if (ierr .ne. 0) then
           print *, 'ierr = ', ierr
           print *, 'cannot read matrix'
           call mpi_finalize(ierr)
      end if

      if (myid .eq. 0) then
            write(*, 70) key, type
   70       format('READ the matrix', 2X, A8, 1X, A3)
            write(*, 80) n, nnz
   80       format('Matrix dimension is ', I9, ' Number of nonzeros 
     &        is', I9)
      end if

!!$symmetric case or not 
      usmy = 0
      if (usmy .eq. 1) then   
           allocate(iwk(n),    stat=ierr)
           allocate(b(nnz),    stat=ierr)
           allocate(jb(nnz),   stat=ierr)
           allocate(ib(n+1),   stat=ierr)
           allocate(mc(2*nnz), stat=ierr)
           allocate(jc(2*nnz), stat=ierr)
           allocate(ic(n+1),   stat=ierr)
#if defined(DBL_CMPLX)           
           call zcsrcsc(n, 1, 1, a, ja, ia, b, jb, ib)
#else
           call csrcsc(n, 1, 1, a, ja, ia, b, jb, ib)
#endif
           do i = 1, nnz
              a(i) = 0.0d0
           end do
           i1 = 2*nnz;
!!$compute C = A + B 
#if defined(DBL_CMPLX)
           call zaplb(n, n, 1, b, jb, ib, a, ja, ia, mc, jc, ic, i1, 
     & iwk, ierr) 
#else
           call aplb(n, n, 1, b, jb, ib, a, ja, ia, mc, jc, ic, i1, 
     & iwk, ierr) 
#endif
           deallocate(iwk, stat=ierr)

           nnz = ic(n+1)-1
           deallocate(ja, stat=ierr)
           deallocate(a,  stat=ierr)
           allocate(ja(nnz), stat=ierr)
           allocate(a(nnz),  stat=ierr)
           do i = 1, n+1
              ia(i) = ic(i)
           end do

           do i = 1, nnz
              ja(i) = jc(i)
           end do
           do i = 1, nnz
              a(i) = mc(i)
           end do
           deallocate(b,  stat=ierr)
           deallocate(jb, stat=ierr)
           deallocate(ib, stat=ierr)
           deallocate(mc, stat=ierr)
           deallocate(jc, stat=ierr)
           deallocate(ic, stat=ierr)
      else
           allocate(b(nnz),    stat=ierr)
           allocate(jb(nnz),   stat=ierr)
           allocate(ib(n+1),   stat=ierr)
!!-------------- convert matrix from csc to csr format ------           
#if defined(DBL_CMPLX)           
           call zcsrcsc(n, 1, 1, a, ja, ia, b, jb, ib)
#else
           call csrcsc(n, 1, 1, a, ja, ia, b, jb, ib)
#endif

!!---------- copy csr matrix -------------------
           do i = 1, n+1
              ia(i) = ib(i)
           end do

           do i = 1, nnz
              ja(i) = jb(i)
           end do
           do i = 1, nnz
              a(i) = b(i)
           end do
           deallocate(b,  stat=ierr)
           deallocate(jb, stat=ierr)
           deallocate(ib, stat=ierr)    
      end if


!!$non-overlapping partitioning
      if (.not. allocated(idom)) then
         allocate(idom(npro+1), stat=ierr)
      end if

      if (.not. allocated(dom)) then
         allocate(dom(n), stat=ierr)
      end if

      if (npro .eq. 1) then
         do i = 1, n
         dom(i) = i
      end do
      idom(1) = 1
      idom(2) = n+1
      else
           if (.not. allocated(riord)) then
              allocate(riord(n), stat=ierr)
           end if
           if (.not. allocated(mask)) then
              allocate(mask(n), stat=ierr)
           end if
           if (.not. allocated(jwk)) then
              allocate(jwk(2*n), stat=ierr)
           end if
           if (.not. allocated(link)) then
              allocate(link(n), stat=ierr)
           end if
           call dse(n, ja, ia, npro, riord, dom, idom, mask, jwk, link)
           deallocate(riord, stat=ierr)
           deallocate(mask,  stat=ierr)
           deallocate(jwk,   stat=ierr)
           deallocate(link,  stat=ierr)
      end if
  
!!$create map object 
      call parms_mapcreatefromptr(map, n, dom, idom, MPI_COMM_WORLD, 
     & 1, NONINTERLACED, ierr);

!!$free dom and idom
      deallocate(dom)
      deallocate(idom)

!! get local size of matrix

      call parms_mapgetlocalsize(map, nloc)

!! Allocate memory for distributed vectors
      allocate(rhs(nloc),    stat=ierr)
      allocate(sol(nloc),    stat=ierr)
      allocate(y(nloc),   stat=ierr)
  

!!$ create a matrix based on map 
      call parms_matcreate(mat, map, ierr);

!!$insert entries into the matrix
      allocate(im(n), stat=ierr)
      do i = 1, n
            im(i) = i
      end do
      call parms_matsetvalues(mat, n, im, ia, ja, a, INSERT,ierr)
!      call parms_matsetelementmatrix(mat, n, im, ia, ja, a, INSERT,ierr)
!      call parms_matassembleelementmatrix(mat,ierr);
!!$free temporary arrays 
      deallocate(a, stat=ierr)
      deallocate(ja,stat=ierr)
      deallocate(ia,stat=ierr)

!!$assemble the matrix 
      call parms_matsetup(mat, ierr)

!!$create a timer
      call parms_timercreate(t, ierr)

!!$create a preconditioner
      call parms_pccreate(pc, mat, ierr)
!!$set parameters for the preconditioner
      call fset_pc_params(pc, prm)
      call parms_timerreset(t, ierr)
      call parms_pcsetup(pc, ierr)
      call parms_timerget(t, tpc, ierr)

!!$pause the timer
      call parms_timerpause(t, ierr)

!!$create a solver
      call parms_solvercreate(ksp, mat, pc, ierr)
!!$set the type of solver ksp
      call parms_solversettype(ksp, SOLFGMRES, ierr)
!!$set parameters for the solver
      call fset_solver_params(ksp, prm)
      call fprm_free(prm)

!!$set actual/ artificial right-hand-side vector
  
      if(nrhs .eq. 0) then
         do i=1,nloc
               sol(i) = 1.d0
         enddo
         call parms_matvec(mat, sol, rhs, ierr)
      else
            call parms_vecsetvalues(rhs, n, im, rhstmp, INSERT, map, 
     &      ierr)
      endif

      deallocate(rhstmp, stat=ierr)
      deallocate(im,stat=ierr)
!!$set initial guess 
      do i=1,nloc
            sol(i) = 0.d0
      enddo

!!$get initial residual error
      call parms_matvec(mat, sol, y, ierr)
      call parms_vecaxpy(y, rhs, (-1.0d0,0.d0), map, ierr)
      call parms_vecgetnorm2(y, res0, map, ierr)

!!$get the residual error
      call parms_timerrestart(t, ierr)
      call parms_solverapply(ksp, rhs, sol, ierr)
      call parms_timerget(t, ttol, ierr)
      call parms_matvec(mat, sol, y, ierr)
      call parms_vecaypx(y, rhs, (-1.0d0,0.d0), map, ierr)
      call parms_vecgetnorm2(y, res1, map, ierr)
      call parms_pcgetratio(pc,ratio, ierr)

      if (myid .eq. 0) then
            call parms_solvergetits(ksp, its, ierr)
            call parms_pcgetname(pc, pcname, len, ierr)
            write(*, 60) pcname(1:len)
            call parms_pcilugetname(pc, pciluname, len, ierr)
            write(*, 66) pciluname(1:len)
            write(*, 78) ratio

            write(*, 86) npro
            write(*, 90) its
            write(*, 100) tpc
            write(*, 200) ttol - tpc
            write(*, 300) ttol
            write(*, 400) res0
            write(*, 500) res1
   60       format('The preconditioner', 12x, A)
   66       format('The local preconditioner', 6x, A)
   78       format('The memory usage ', 10x, '=', f5.2)
   86       format('The number of processors ', 2x, '=', I3)
   90       format('The number of iteration', 4x, '=', I3)
  100       format('The time for creating pc ', 2x, '=', es9.2,'s')
  200       format('The solving time ', 10x, '=', es9.2, 's')
  300       format('The total time ', 12x, '=', es9.2, 's')
  400       format('The initial residual error =', es9.2)
  500       format('The solution residual error =', es9.2)
      end if

!!$ free memories for map, vec, mat, solver, etc.
      deallocate(sol, stat=ierr)
      deallocate(rhs,stat=ierr)
      deallocate(y,stat=ierr)
      call parms_solverfree(ksp, ierr)
      call parms_mapfree(map, ierr)
      call parms_matfree(mat, ierr)
      call parms_timerfree(t, ierr)
      call parms_pcfree(pc, ierr)

!!$exit MPI environment
      call mpi_finalize(ierr)
      end program main       
