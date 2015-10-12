      subroutine part1 (nx,ny,mpx,mpy,ovlp,lst,lstptr,iout) 
      integer  nx,ny,mpx,mpy,ovlp,iout,lst(*),lstptr(*)
c-----------------------------------------------------------------------
c	
c     does a trivial map of a square grid into a virtual mpx x mpy 
c     processor grid. Overlapping allowed: when assigning a block to a 
c     subdomain then `ovlp' extra lines  are added to each of the four 
c     sides of the subrectangle being considered, outward.
c 
c-----------------------------------------------------------------------
c on entry:
c--------- 
c nx, ny   = number of mesh points in x and y directions respectively
c mpx, mpy = number of processors in x and y directions respectively
c            with (mpx .le. nx)  and (mpy .le. ny)
c ovlp     = a nonnegative integer determing the amount of overlap 
c            (number of lines) in each direction. 
c
c on return:
c ---------
c lst      = node per processor list. The nodes are listed contiguously
c            from proc 1 to nproc = mpx*mpy. 
c lstptr   = pointer array for array lst. list for proc. i starts at 
c            lstptr(i) and ends at lstptr(i+1)-1 in array lst.
c-----------------------------------------------------------------------
c iout     = not used.
c-----------------------------------------------------------------------
      nproc = mpx*mpy 
      mx = nx/mpx
      my = ny/mpy 
      nod = 0
      ipr = 1 
      ko = 1
      lstptr(ipr) = 1 
c-----------------------------------------------------------------------
      do 1 jj = 1, mpy
         j = (jj-1)*my
         do 2 ii = 1, mpx
            i = (ii-1)*mx
c           write (iout,*) ' *** * proc = ',ipr, ' ****** '
            mymax = min(my+ovlp,ny-j)
            mxmax = min(mx+ovlp,nx-i) 
            mymin = - ovlp 
            if (mymin + j .lt. 0) mymin = 0
            mxmin = -ovlp 
            if (mxmin + i .lt. 0) mxmin = 0
c-----------------------------------------------------------------------
            do ky = mymin, mymax-1
               do kx = mxmin, mxmax-1                  
                  nod = (j+ky)*nx + i + kx + 1 
                  lst(ko) = nod
c                  write (iout,*) ' ko ', ko, ' nod ', nod 
                  ko = ko+1 
               enddo
            enddo
            ipr = ipr+1 
            lstptr(ipr) = ko
 2       continue
 1    continue 
c-----------------------------------------------------------------------
      return
c----------------------------------------------------------------------- 
      end

      subroutine part2 (nx,ny,nz,mpx,mpy,mpz,ovlp,lst,lstptr,iout) 
      integer  nx,ny,nz,mpx,mpy,mpz,ovlp,iout,lst(*),lstptr(*)
c-----------------------------------------------------------------------
c	
c     does a trivial map of a square grid into a virtual mpx x mpy 
c     processor grid. Overlapping allowed: when assigning a block to a 
c     subdomain then `ovlp' extra lines  are added to each of the four 
c     sides of the subrectangle being considered, outward.
c 
c-----------------------------------------------------------------------
c on entry:
c--------- 
c nx, ny,nz   = number of mesh points in x, y and z directions respectively
c mpx,mpy,mpz = number of processors in x, y and z directions respectively
c            with (mpx .le. nx), (mpy .le. ny) and (mpz .le. nz)
c ovlp     = a nonnegative integer determing the amount of overlap 
c            (number of lines) in each direction. 
c
c on return:
c ---------
c lst      = node per processor list. The nodes are listed contiguously
c            from proc 1 to nproc = mpx*mpy. 
c lstptr   = pointer array for array lst. list for proc. i starts at 
c            lstptr(i) and ends at lstptr(i+1)-1 in array lst.
c-----------------------------------------------------------------------
c iout     = not used.
c-----------------------------------------------------------------------
      nproc = mpx*mpy*mpz 
      mx = nx/mpx
      my = ny/mpy 
      mz = nz/mpz
      nod = 0
      ipr = 1 
      ko = 1
      lstptr(ipr) = 1 
c-----------------------------------------------------------------------
      do kk = 1, mpz
         k = (kk-1)*mz
         do 1 jj = 1, mpy
            j = (jj-1)*my
            do 2 ii = 1, mpx
               i = (ii-1)*mx
c           write (iout,*) ' *** * proc = ',ipr, ' ****** '
               mzmax = min(mz+ovlp,nz-k)
               mymax = min(my+ovlp,ny-j)
               mxmax = min(mx+ovlp,nx-i)
               mzmin = - ovlp
               if (mzmin + k .lt. 0) mzmin = 0
               mymin = - ovlp 
               if (mymin + j .lt. 0) mymin = 0
               mxmin = -ovlp 
               if (mxmin + i .lt. 0) mxmin = 0
c-----------------------------------------------------------------------
               do kz = mzmin, mzmax-1
                  do ky = mymin, mymax-1
                     do kx = mxmin, mxmax-1                  
                        nod = (k+kz)*nx*ny+ (j+ky)*nx + i + kx + 1 
                        lst(ko) = nod
c     write (iout,*) ' ko ', ko, ' nod ', nod 
                        ko = ko+1 
                     enddo
                  enddo
               end do
               ipr = ipr+1 
               lstptr(ipr) = ko
 2          continue
 1       continue 
      end do
c-----------------------------------------------------------------------
      return
c----------------------------------------------------------------------- 
      end
