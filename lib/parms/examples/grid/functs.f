c-----------------------------------------------------------------------
c contains the functions needed for defining the PDE problems. 
c 
c Extracted from SPARSEKIT
c problem is :
c
c format for the elliptic equation:
c       d    du    d    du    d    du      du     du     du
c L u = --(A --) + --(B --) + --(C --) + D -- + E -- + F -- + G u = .. 
c       dx   dx    dy   dy    dz   dz      dx     dy     dz
c
c Following 
c-----------------------------------------------------------------------
      function afun (x,y,z)
      real*8 afun, x,y, z 
      afun = -1.0
      return 
      end
c
      function bfun (x,y,z)
      real*8 bfun, x,y, z 
      bfun = -1.0
      return 
      end
c      
      function cfun (x,y,z)
      real*8 cfun, x,y, z 
      cfun = -1.0d0
      return 
      end
c
      function dfun (x,y,z)
      real*8 dfun, x,y, z 
      data gamma /100.0/ 
      dfun = gamma*dexp(x*y)
      return 
      end
c
      function efun (x,y,z)
      real*8 efun, x,y, z
      data gamma /100.0/ 
      efun = gamma*dexp(-x*y) 
      return 
      end
c      
      function ffun (x,y,z)
      real*8 ffun, x,y, z 
      ffun = 0.0
      return 
      end
c
      function gfun (x,y,z)
      real*8 gfun, x,y, z 
c
c     gfun negative may make the problem indefinite
c
      gfun = -10.0d0 
      return 
      end
c-----------------------------------------------------------------------
c functions for the block PDE's 
c-----------------------------------------------------------------------
      subroutine afunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = 0.0d0
 1       continue
         coeff((j-1)*nfree+j) = -1.0d0
 2    continue
      return 
      end
      subroutine bfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = 0.0d0
 1       continue
         coeff((j-1)*nfree+j) = -1.0d0
 2    continue
      return 
      end
        subroutine cfunbl (nfree,x,y,z,coeff)
        real*8 x, y, z, coeff(100) 
	do 2 j=1, nfree
	   do 1 i=1, nfree
              coeff((j-1)*nfree+i) = 0.0d0
 1         continue
           coeff((j-1)*nfree+j) = -1.0d0
 2	continue
        return 
	end
      subroutine dfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = 0.0d0
 1       continue
 2    continue
      return 
      end
c     
      subroutine efunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = 0.0d0
 1       continue
 2    continue
      return 
      end
c     
      subroutine ffunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
            coeff((j-1)*nfree+i) = 0.0d0
 1       continue
 2    continue
      return 
      end
c     
      subroutine gfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = 0.0d0
 1       continue
 2    continue
      return 
      end
c-----------------------------------------------------------------------
c     The material property function xyk for the 
c     finite element problem 
c-----------------------------------------------------------------------
      subroutine xyk(nel,xyke,x,y,ijk,node)
      implicit real*8 (a-h,o-z)
      dimension xyke(2,2), x(*), y(*), ijk(node,*)
c     
c     this is the identity matrix.
c     
      xyke(1,1) = 1.0d0
      xyke(2,2) = 1.0d0
      xyke(1,2) = 0.0d0
      xyke(2,1) = 0.0d0
c     
      return
c-----------------------------------------------------------------------
      end

      subroutine setinit(x,n)
      implicit none
      double precision x(*)
      integer n,j

      do j=1, n
         x(j) = sin(4.0*3.1415926*j/n) 
      enddo
      return 
      end



c----------------------------------------------------------------------- 
c-----------------------------------------------------------------------
c-------------- Routines to generate complex-valued problem ------------
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c contains the functions needed for defining the PDE problems. 
c 
c Extracted from SPARSEKIT
c problem is :
c
c format for the elliptic equation:
c       d    du    d    du    d    du      du     du     du
c L u = --(A --) + --(B --) + --(C --) + D -- + E -- + F -- + G u = .. 
c       dx   dx    dy   dy    dz   dz      dx     dy     dz
c
c Following 
c-----------------------------------------------------------------------
      function zafun (x,y,z)
      double complex zafun
      real*8 x,y, z 
      zafun = cmplx(-1.0,0.0d0)
      return 
      end
c
      function zbfun (x,y,z)
      real*8 x, y, z
      complex zbfun 
      zbfun = cmplx(-1.0, 0.0d0)
      return 
      end
c      
      function zcfun (x,y,z)
      complex zcfun
      real*8 x,y, z 
      zcfun = cmplx(-1.0d0, 0.0d0)
      return 
      end
c
      function zdfun (x,y,z)
      complex zdfun
      real*8 x,y, z 
      data gamma /100.0/ 
      zdfun = cmplx(gamma*dexp(x*y), 0.0)
      return 
      end
c
      function zefun (x,y,z)
      complex zefun
      real*8 x,y, z
      data gamma /100.0/ 
      zefun = cmplx(gamma*dexp(-x*y), 0.0) 
      return 
      end
c      
      function zffun (x,y,z)
      complex zffun
      real*8 x,y, z 
      zffun = cmplx(0.0d0, 0.0d0)
      return 
      end
c
      function zgfun (x,y,z)
      complex zgfun
      real*8 x,y, z 
c
c     gfun negative may make the problem indefinite
c
      zgfun = cmplx(-10.0d0, 0.0d0) 
      return 
      end
c-----------------------------------------------------------------------
c functions for the block PDE's 
c-----------------------------------------------------------------------
      subroutine zafunbl (nfree,x,y,z,coeff)
      real*8 x, y, z 
      double complex coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = cmplx(0.0d0, 0.0d0)
 1       continue
         coeff((j-1)*nfree+j) = cmplx(-1.0d0, 0.0d0)
 2    continue
      return 
      end
      subroutine zbfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z
      double complex coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = cmplx(0.0d0, 0.0d0)
 1       continue
         coeff((j-1)*nfree+j) = cmplx(-1.0d0, 0.0d0)
 2    continue
      return 
      end
        subroutine zcfunbl (nfree,x,y,z,coeff)
        real*8 x, y, z
        double complex coeff(100) 
	do 2 j=1, nfree
	   do 1 i=1, nfree
              coeff((j-1)*nfree+i) = cmplx(0.0d0, 0.0d0)
 1         continue
           coeff((j-1)*nfree+j) = cmplx(-1.0d0, 0.0d0)
 2	continue
        return 
	end
      subroutine zdfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z
      double complex coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = cmplx(0.0d0, 0.0d0)
 1       continue
 2    continue
      return 
      end
c     
      subroutine zefunbl (nfree,x,y,z,coeff)
      real*8 x, y, z
      double complex coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = cmplx(0.0d0, 0.0d0)
 1       continue
 2    continue
      return 
      end
c     
      subroutine zffunbl (nfree,x,y,z,coeff)
      real*8 x, y, z
      double complex coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
            coeff((j-1)*nfree+i) = cmplx(0.0d0, 0.0d0)
 1       continue
 2    continue
      return 
      end
c     
      subroutine zgfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z
      double complex coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
	    coeff((j-1)*nfree+i) = cmplx(0.0d0, 0.0d0)
 1       continue
 2    continue
      return 
      end
c-----------------------------------------------------------------------
c     The material property function xyk for the 
c     finite element problem 
c-----------------------------------------------------------------------
      subroutine zxyk(nel,xyke,x,y,ijk,node)
      implicit real*8 (a-h,o-z)
      dimension xyke(2,2), x(*), y(*), ijk(node,*)
c     
c     this is the identity matrix.
c     
      xyke(1,1) = 1.0d0
      xyke(2,2) = 1.0d0
      xyke(1,2) = 0.0d0
      xyke(2,1) = 0.0d0
c     
      return
c-----------------------------------------------------------------------
      end

      subroutine zsetinit(x,n)
      implicit none
      double precision x(*)
      integer n,j

      do j=1, n
         x(j) = sin(4.0*3.1415926*j/n) 
      enddo
      return 
      end

