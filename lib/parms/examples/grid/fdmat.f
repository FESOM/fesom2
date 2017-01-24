      subroutine gen5pt (nx,ny,nz,a,ja,ia,iau) 
      integer ja(*),ia(*),iau(*)
      real*8 a(*), h
c--------------------------------------------------------------------
c This subroutine computes the sparse matrix in compressed
c format for the elliptic operator
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + 
c       d delx ( u ) + e dely (u) + f delz( u ) + g u
c
c with Dirichlet Boundary conditions, on a rectangular 1-D, 
c 2-D or 3-D grid using centered difference schemes.
c 
c The functions a, b, ..., g are known through the
c subroutines  afun, bfun, ..., gfun.
c note that to obtain the correct matrix, any function that is not
c needed should be set to zero. For example for two-dimensional
c problems, nz should be set to 1 and the functions cfun and ffun
c should be zero functions. 
c
c uses natural ordering, first x direction, then y, then z
c mesh size h is uniform and determined by grid points 
c in the x-direction.
c--------------------------------------------------------------------
c parameters:
c
c nx      = number of points in x direction
c ny	  = number of points in y direction
c nz	  = number of points in z direction
c
c a, ja, ia =  resulting matrix in row-sparse format
c
c iau     = integer*n containing the poisition of the diagonal element
c           in the a, ja, ia structure
c
c stencil =  work array of size 7, used to store local stencils.
c 
c--------------------------------------------------------------------
c
c     stencil [1:7] has the following meaning:
c
c     center point = stencil(1)
c     west point = stencil(2)
c     east point = stencil(3)
c     south point = stencil(4)
c     north point = stencil(5)
c     front point = stencil(6) 
c     back point = stencil(7)
c
c
c                           st(5)
c                            |
c                            |  
c                            |
c                            |          .st(7)
c                            |     .
c                            | . 
c         st(2) ----------- st(1) ---------- st(3)
c                       .    |
c                   .        |
c               .            |
c            st(6)           |
c                            |
c                            |
c                           st(4)
c
c-------------------------------------------------------------------
      real*8 stencil(7) 
c
      h = 1.0/real(nx+1)	
      kx = 1
      ky = nx
      kz = nx*ny
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
               call getsten(nx,ny,nz,ix,iy,iz,stencil,h)
c     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
		  a(iedge) = stencil(2)
                  iedge=iedge + 1
               end if
c     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
		  a(iedge) = stencil(4)
                  iedge=iedge + 1
               end if
c     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
		  a(iedge) = stencil(6)
                  iedge=iedge + 1
               endif
c     center node
               ja(iedge) = node
               iau(node) = iedge
               a(iedge) = stencil(1)
               iedge = iedge + 1
c     -- upper part  
c     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
		  a(iedge) = stencil(3)
                  iedge=iedge + 1
               end if
c     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
		  a(iedge) = stencil(5)
                  iedge=iedge + 1
               end if
c     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
                  a(iedge) = stencil(7)
                  iedge=iedge + 1
               end if
c------next node -------------------------
               node=node+1
 80         continue
 90      continue
 100  continue
      ia(node)=iedge
      return
c----------end of gen5pt-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine gen5loc (nx,ny,nz,nloc,riord,a,ja,ia,stencil)
      integer nx, ny, nz, nloc, ja(*),ia(*),riord(*) 
      real*8 a(*), stencil(*), h
c--------------------------------------------------------------------
c Local version of gen5pt. This subroutine generates only the nloc
c rows of the matrix that are specified by the integer array riord
c rows  riord(1), ..., riord(nloc) of the matrix will be generated
c and put in the a, ja, ia  data structure. The column indices in 
c ja will still be given in the original labeling. 
c-----------------------------------------------------------------------
c see documentation for gen5pt for mode details.
c on input
c ------ 
c same arguments as for gen5pt. Plus:
c
c nloc   = number of rows to be generated
c riord  = integer array of length nloc. the code will generate rows
c          riord(1), ..., riord(nloc) of the matrix. 
c 
c on return
c----------
c a, ja, ia =  resulting matrix in CSR format.
c              can be  viewed as a rectangular matrix of size 
c              nloc x n, containing the rows riord(1),riord(2), ...
c              riord(nloc) of A (in this order)  in CSR format.
c-------------------------------------------------------------------
      h = 1.0/real(nx+1)	
      kx = 1
      ky = nx
      kz = nx*ny
      iedge = 1
      node = 1
      next = 1
      ia(next) = iedge
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               if (node .eq. riord(next)) then
                  call getsten(nx,ny,nz,ix,iy,iz,stencil,h)
c     west
                  if (ix.gt.1) then
                     ja(iedge)=node-kx
                     a(iedge) = stencil(2)
                     iedge=iedge + 1
                  end if
c     south
                  if (iy.gt.1) then
                     ja(iedge)=node-ky
                     a(iedge) = stencil(4)
                     iedge=iedge + 1
                  end if
c     front plane
                  if (iz.gt.1) then
                     ja(iedge)=node-kz
                     a(iedge) = stencil(6)
                     iedge=iedge + 1
                  endif
c     center node
                  ja(iedge) = node
                  a(iedge) = stencil(1)
                  iedge = iedge + 1
c     -- upper part  
c     east
                  if (ix.lt.nx) then
                     ja(iedge)=node+kx
                     a(iedge) = stencil(3)
                     iedge=iedge + 1
                  end if
c     north
                  if (iy.lt.ny) then
                     ja(iedge)=node+ky
                     a(iedge) = stencil(5)
                     iedge=iedge + 1
                  end if
c     back plane
                  if (iz.lt.nz) then
                     ja(iedge)=node+kz
                     a(iedge) = stencil(7)
                     iedge=iedge + 1
                  end if
c------next node -------------------------
                  next = next+1
                  ia(next) = iedge
               endif
               node = node + 1
               if (next .gt. nloc) return 
 80         continue 
 90      continue
 100  continue
      return
c----------end of gen5loc-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
       subroutine getsten (nx,ny,nz,kx,ky,kz,stencil,h)
c-----------------------------------------------------------------------
c     This subroutine calcultes the correct stencil values for
c     centered difference discretization of the elliptic operator
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + 
c	delx ( d u ) + dely (e u) + delz( f u ) + g u 
c
c   For 2-D problems the discretization formula that is used is:
c      
c h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
c	       a(i-1/2,j)*{u(i-1,j) - u(i,j)} + 
c              b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
c              b(i,j-1/2)*{u(i,j-1) - u(i,j)} + 
c              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + 
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + 
c               (h**2)*g(i,j)*u(i,j) 
c-----------------------------------------------------------------------
      real*8 stencil(*), h, hhalf,cntr, afun, bfun, cfun, dfun, 
     *      efun, ffun, gfun, x, y, z, coeff
c------------
      do 200 k=1,7
         stencil(k) = 0.0
 200  continue
c------------
      hhalf = h*0.5
      x = h*real(kx)
      y = h*real(ky)
      z = h*real(kz)
      cntr = 0.0
c     differentiation wrt x:
      coeff = afun(x+hhalf,y,z)
      stencil(3) = stencil(3) + coeff
      cntr = cntr + coeff
c     
      coeff = afun(x-hhalf,y,z)
      stencil(2) = stencil(2) + coeff
      cntr = cntr + coeff
c     
      coeff = dfun(x,y,z)*hhalf
      stencil(3) = stencil(3) + coeff
      stencil(2) = stencil(2) - coeff
      if (ny .le. 1) goto 99
c     
c     differentiation wrt y:
c     
      coeff = bfun(x,y+hhalf,z)
      stencil(5) = stencil(5) + coeff
      cntr = cntr + coeff
c     
      coeff = bfun(x,y-hhalf,z)
      stencil(4) = stencil(4) + coeff
      cntr = cntr + coeff
c     
      coeff = efun(x,y,z)*hhalf
      stencil(5) = stencil(5) + coeff
      stencil(4) = stencil(4) - coeff
      if (nz .le. 1) goto 99
c
c differentiation wrt z:
c 
      coeff = cfun(x,y,z+hhalf)
      stencil(7) = stencil(7) + coeff
      cntr = cntr + coeff
c
      coeff = cfun(x,y,z-hhalf)
      stencil(6) = stencil(6) + coeff
      cntr = cntr + coeff
c
      coeff = ffun(x,y,z)*hhalf
      stencil(7) = stencil(7) + coeff
      stencil(6) = stencil(6) - coeff
c
c discretization of  product by g:
c
 99   coeff = gfun(x,y,z)
      stencil(1) = h*h*coeff - cntr
c     
      return
c------end-of-getsten--------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
c-----------------------------------------------------------------------
c-------------- Routines to generate complex-valued problem ------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine zgen5pt (nx,ny,nz,a,ja,ia,iau) 
      integer ja(*),ia(*),iau(*)
      double complex a(*)
      real*8 h
c--------------------------------------------------------------------
c This subroutine computes the sparse matrix in compressed
c format for the elliptic operator
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + 
c       d delx ( u ) + e dely (u) + f delz( u ) + g u
c
c with Dirichlet Boundary conditions, on a rectangular 1-D, 
c 2-D or 3-D grid using centered difference schemes.
c 
c The functions a, b, ..., g are known through the
c subroutines  afun, bfun, ..., gfun.
c note that to obtain the correct matrix, any function that is not
c needed should be set to zero. For example for two-dimensional
c problems, nz should be set to 1 and the functions cfun and ffun
c should be zero functions. 
c
c uses natural ordering, first x direction, then y, then z
c mesh size h is uniform and determined by grid points 
c in the x-direction.
c--------------------------------------------------------------------
c parameters:
c
c nx      = number of points in x direction
c ny	  = number of points in y direction
c nz	  = number of points in z direction
c
c a, ja, ia =  resulting matrix in row-sparse format
c
c iau     = integer*n containing the poisition of the diagonal element
c           in the a, ja, ia structure
c
c stencil =  work array of size 7, used to store local stencils.
c 
c--------------------------------------------------------------------
c
c     stencil [1:7] has the following meaning:
c
c     center point = stencil(1)
c     west point = stencil(2)
c     east point = stencil(3)
c     south point = stencil(4)
c     north point = stencil(5)
c     front point = stencil(6) 
c     back point = stencil(7)
c
c
c                           st(5)
c                            |
c                            |  
c                            |
c                            |          .st(7)
c                            |     .
c                            | . 
c         st(2) ----------- st(1) ---------- st(3)
c                       .    |
c                   .        |
c               .            |
c            st(6)           |
c                            |
c                            |
c                           st(4)
c
c-------------------------------------------------------------------
      double complex stencil(7) 
c
      h = 1.0/real(nx+1)	
      kx = 1
      ky = nx
      kz = nx*ny
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
               call zgetsten(nx,ny,nz,ix,iy,iz,stencil,h)
c     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
		  a(iedge) = stencil(2)
                  iedge=iedge + 1
               end if
c     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
		  a(iedge) = stencil(4)
                  iedge=iedge + 1
               end if
c     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
		  a(iedge) = stencil(6)
                  iedge=iedge + 1
               endif
c     center node
               ja(iedge) = node
               iau(node) = iedge
               a(iedge) = stencil(1)
               iedge = iedge + 1
c     -- upper part  
c     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
		  a(iedge) = stencil(3)
                  iedge=iedge + 1
               end if
c     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
		  a(iedge) = stencil(5)
                  iedge=iedge + 1
               end if
c     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
                  a(iedge) = stencil(7)
                  iedge=iedge + 1
               end if
c------next node -------------------------
               node=node+1
 80         continue
 90      continue
 100  continue
      ia(node)=iedge
      return
c----------end of zgen5pt-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zgen5loc (nx,ny,nz,nloc,riord,a,ja,ia,stencil)
      integer nx, ny, nz, nloc, ja(*),ia(*),riord(*) 
      double complex a(*)
      real*8 h, stencil(*)
c--------------------------------------------------------------------
c Local version of gen5pt. This subroutine generates the complex version
c of the real counterpart. The imaginary part is set to zero. Can be 
c made to generate a fully complex system by defining stencil(*) as 
c 'double complex' and using the zgetsten function instead of getsten.
c This subroutine generates only the nloc
c rows of the matrix that are specified by the integer array riord
c rows  riord(1), ..., riord(nloc) of the matrix will be generated
c and put in the a, ja, ia  data structure. The column indices in 
c ja will still be given in the original labeling. 
c-----------------------------------------------------------------------
c see documentation for gen5pt for mode details.
c on input
c ------ 
c same arguments as for gen5pt. Plus:
c
c nloc   = number of rows to be generated
c riord  = integer array of length nloc. the code will generate rows
c          riord(1), ..., riord(nloc) of the matrix. 
c 
c on return
c----------
c a, ja, ia =  resulting matrix in CSR format.
c              can be  viewed as a rectangular matrix of size 
c              nloc x n, containing the rows riord(1),riord(2), ...
c              riord(nloc) of A (in this order)  in CSR format.
c-------------------------------------------------------------------
      h = 1.0/real(nx+1)	
      kx = 1
      ky = nx
      kz = nx*ny
      iedge = 1
      node = 1
      next = 1
      ia(next) = iedge
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               if (node .eq. riord(next)) then
                  call getsten(nx,ny,nz,ix,iy,iz,stencil,h)
c     west
                  if (ix.gt.1) then
                     ja(iedge)=node-kx
                     a(iedge) = stencil(2)
                     iedge=iedge + 1
                  end if
c     south
                  if (iy.gt.1) then
                     ja(iedge)=node-ky
                     a(iedge) = stencil(4)
                     iedge=iedge + 1
                  end if
c     front plane
                  if (iz.gt.1) then
                     ja(iedge)=node-kz
                     a(iedge) = stencil(6)
                     iedge=iedge + 1
                  endif
c     center node
                  ja(iedge) = node
                  a(iedge) = stencil(1)
                  iedge = iedge + 1
c     -- upper part  
c     east
                  if (ix.lt.nx) then
                     ja(iedge)=node+kx
                     a(iedge) = stencil(3)
                     iedge=iedge + 1
                  end if
c     north
                  if (iy.lt.ny) then
                     ja(iedge)=node+ky
                     a(iedge) = stencil(5)
                     iedge=iedge + 1
                  end if
c     back plane
                  if (iz.lt.nz) then
                     ja(iedge)=node+kz
                     a(iedge) = stencil(7)
                     iedge=iedge + 1
                  end if
c------next node -------------------------
                  next = next+1
                  ia(next) = iedge
               endif
               node = node + 1
               if (next .gt. nloc) return 
 80         continue 
 90      continue
 100  continue
      return
c----------end of zgen5loc-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
       subroutine zgetsten (nx,ny,nz,kx,ky,kz,stencil,h)
c-----------------------------------------------------------------------
c     This subroutine calcultes the correct stencil values for
c     centered difference discretization of the elliptic operator
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + 
c	delx ( d u ) + dely (e u) + delz( f u ) + g u 
c
c   For 2-D problems the discretization formula that is used is:
c      
c h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
c	       a(i-1/2,j)*{u(i-1,j) - u(i,j)} + 
c              b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
c              b(i,j-1/2)*{u(i,j-1) - u(i,j)} + 
c              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + 
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + 
c               (h**2)*g(i,j)*u(i,j) 
c-----------------------------------------------------------------------
      real*8 h, hhalf, x, y, z
      double complex stencil(*),cntr, zafun, zbfun, zcfun, zdfun, 
     *      zefun, zffun, zgfun, coeff
c------------
      do 200 k=1,7
         stencil(k) = 0.0
 200  continue
c------------
      hhalf = h*0.5
      x = h*real(kx)
      y = h*real(ky)
      z = h*real(kz)
      cntr = cmplx(0.0, 0.0)
c     differentiation wrt x:
      coeff = zafun(x+hhalf,y,z)
      stencil(3) = stencil(3) + coeff
      cntr = cntr + coeff
c     
      coeff = zafun(x-hhalf,y,z)
      stencil(2) = stencil(2) + coeff
      cntr = cntr + coeff
c     
      coeff = zdfun(x,y,z)*hhalf
      stencil(3) = stencil(3) + coeff
      stencil(2) = stencil(2) - coeff
      if (ny .le. 1) goto 99
c     
c     differentiation wrt y:
c     
      coeff = zbfun(x,y+hhalf,z)
      stencil(5) = stencil(5) + coeff
      cntr = cntr + coeff
c     
      coeff = zbfun(x,y-hhalf,z)
      stencil(4) = stencil(4) + coeff
      cntr = cntr + coeff
c     
      coeff = zefun(x,y,z)*hhalf
      stencil(5) = stencil(5) + coeff
      stencil(4) = stencil(4) - coeff
      if (nz .le. 1) goto 99
c
c differentiation wrt z:
c 
      coeff = zcfun(x,y,z+hhalf)
      stencil(7) = stencil(7) + coeff
      cntr = cntr + coeff
c
      coeff = zcfun(x,y,z-hhalf)
      stencil(6) = stencil(6) + coeff
      cntr = cntr + coeff
c
      coeff = zffun(x,y,z)*hhalf
      stencil(7) = stencil(7) + coeff
      stencil(6) = stencil(6) - coeff
c
c discretization of  product by g:
c
 99   coeff = zgfun(x,y,z)
      stencil(1) = h*h*coeff - cntr
c     
      return
c------end-of-zgetsten--------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
