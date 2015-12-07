! ========================================================================
! Halo exchange routines + broadcast routines that collect information
! on the entire field (needed for output)
! The routines here are very similar, difference is the data type and  
! exchange pattern.
! exchange_nod2D_i(arr(myDim_nod2D+eDim_nod2D))    INTEGER
! exchange_nod2D(arr(myDim_nod2D+eDim_nod2D))      WP
! exchange_nod3D(arr(nl-1,myDim_nod2D+eDim_nod2D)) WP
! exchange_nod3D_full(arr(nl,myDim_nod2D+eDim_nod2D)) WP
! exchange_edge2D(edge_array2D)     WP  not used currently  !!! no buffer!!!  
! exchange_edge3D(edge_array3D)     WP  not used currently  !!! no buffer!!!
! exchange_elem3D(elem_array3D)     WP
! exchange_elem3D_full(elem_array3D) WP (eDim+ eXDim are exchanged)
! exchange_elem2d_full
! exchange_elem2d_full_i
! ========================================================================
module g_comm
implicit none
contains

subroutine exchange_nod2D_i(nod_array2D)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 2D nodal fields
! of integer type
 
 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 integer  :: nod_array2D(:) 
 
  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum
  
  ! Put data to be communicated into send buffer 
 
  Do n=1,com_nod2D%sPEnum
   nini=com_nod2D%sptr(n)
   nend=com_nod2D%sptr(n+1)-1
   s_buff_nod2D_i(n)%array=nod_array2D(com_nod2D%slist(nini:nend))
  end do  
 
  DO n=1, sn
     dest=com_nod2D%sPE(n)
     nini=com_nod2D%sptr(n)
     offset=com_nod2D%sptr(n+1) - nini
     
     
     call MPI_ISEND(s_buff_nod2D_i(n)%array, offset, MPI_INTEGER, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1,rn    
     source=com_nod2D%rPE(n)
     nini=com_nod2D%rptr(n)
     offset=com_nod2D%rptr(n+1) - nini
     
      
     call MPI_IRECV(r_buff_nod2D_i(n)%array, offset, MPI_INTEGER, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
  Do n=1,com_nod2D%rPEnum
   nini=com_nod2D%rptr(n)
   nend=com_nod2D%rptr(n+1)-1
   nod_array2D(com_nod2D%rlist(nini:nend))=r_buff_nod2D_i(n)%array
  end do 
  
  
END SUBROUTINE exchange_nod2D_i

! ========================================================================
subroutine exchange_nod2D(nod_array2D)
USE o_MESH
USE g_PARSUP
USE o_PARAM 
IMPLICIT NONE

! General version of the communication routine for 2D nodal fields
 
 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer    :: n, sn, rn, dest, nini, nend, offset, source,tag
 !real*8 :: nod_array2D(myDim_nod2D+eDim_nod2D) 
 real*8 :: nod_array2D(:)
  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum
  ! Put data to be communicated into send buffer 
  !nod_array2D=0
  Do n=1,com_nod2D%sPEnum
   nini=com_nod2D%sptr(n)
   nend=com_nod2D%sptr(n+1)-1
   s_buff_nod2D(n)%array=nod_array2D(com_nod2D%slist(nini:nend))
  end do  
 
  DO n=1, sn
     dest=com_nod2D%sPE(n)
     nini=com_nod2D%sptr(n)
     offset=com_nod2D%sptr(n+1) - nini
     
     
     call MPI_ISEND(s_buff_nod2D(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
!call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  DO n=1,rn    
     source=com_nod2D%rPE(n)
     nini=com_nod2D%rptr(n)
     offset=com_nod2D%rptr(n+1) - nini
     
      
     call MPI_IRECV(r_buff_nod2D(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
  Do n=1,com_nod2D%rPEnum
   nini=com_nod2D%rptr(n)
   nend=com_nod2D%rptr(n+1)-1
   nod_array2D(com_nod2D%rlist(nini:nend))=r_buff_nod2D(n)%array
  end do 
   
END SUBROUTINE exchange_nod2D
! ========================================================================
subroutine exchange_nod3D(nod_array3D)
USE o_MESH
USE g_PARSUP
USE o_PARAM 
IMPLICIT NONE

! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
 
 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 integer  :: nz, nh, nc, nl1
 !real(kind=WP) :: nod_array3D(:,myDim_nod2D+eDim_nod2D) 
 real(kind=WP) :: nod_array3D(:,:) 
  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum
  ! Put data to be communicated into send buffer 
  nl1=ubound(nod_array3D,1)
  Do n=1,com_nod2D%sPEnum
   nini=com_nod2D%sptr(n)
   nend=com_nod2D%sptr(n+1)-1
   nc=0
   DO nh=nini, nend 
      DO nz=1,nl1
      nc=nc+1
      s_buff_nod3D(n)%array(nc)=nod_array3D(nz,com_nod2D%slist(nh))
      END DO
   END DO   
  end do  
 
  DO n=1, sn
     dest=com_nod2D%sPE(n)
     nini=com_nod2D%sptr(n)
     offset=(com_nod2D%sptr(n+1) - nini)*(nl1)
     
     
     call MPI_ISEND(s_buff_nod3D(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1,rn    
     source=com_nod2D%rPE(n)
     nini=com_nod2D%rptr(n)
     offset=(com_nod2D%rptr(n+1) - nini)*(nl1)
     
      
     call MPI_IRECV(r_buff_nod3D(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
  Do n=1,com_nod2D%rPEnum
   nini=com_nod2D%rptr(n)
   nend=com_nod2D%rptr(n+1)-1
   nc=0
   DO nh=nini, nend
      DO nz=1, nl1
      nc=nc+1
      nod_array3D(nz,com_nod2D%rlist(nh))=r_buff_nod3D(n)%array(nc)
      END DO
   END DO   
  end do 
   
END SUBROUTINE exchange_nod3D
! ========================================================================
subroutine exchange_edge3D(edge_array3D)
  use o_MESH
  use g_PARSUP 
  implicit none
  
 ! Communication of edge based data stored in (vertical, horizontal) format 

 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 integer  :: nz, nh, nc
 real(kind=WP) :: edge_array3D(nl-1,edge2D) 
  
  sn=com_edge2D%sPEnum
  rn=com_edge2D%rPEnum
  ! Put data to be communicated into send buffer 


  do n=1, sn
     nini=com_edge2D%sptr(n)
     nend=com_edge2D%sptr(n+1) - 1
      nc=0
      DO nh=nini, nend
         DO nz=1, nl-1
	 nc=nc+1
         s_buff_edge3D(n)%array(nc)=edge_array3D(nz,com_edge2D%slist(nh))
	 END DO
      END DO	 
  end do


  do n=1, sn
     dest=com_edge2D%sPE(n)
     nini=com_edge2D%sptr(n)
     offset=(com_edge2D%sptr(n+1) - nini)*(nl-1)

     call MPI_ISEND(s_buff_edge3D(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
          MPI_COMM_WORLD, sreq(n), MPIerr)
  end do
  do n=1, rn
     source=com_edge2D%rPE(n)
     nini=com_edge2D%rptr(n)
     offset=(com_edge2D%rptr(n+1) - nini)*(nl-1)

     call MPI_IRECV(r_buff_edge3D(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
          source, MPI_COMM_WORLD, rreq(n), MPIerr) 
  end do

  call MPI_WAITALL(sn,sreq,sstat, MPIerr)
  call MPI_WAITALL(rn,rreq,rstat, MPIerr)

  ! Put received data to their destination

  do n=1, rn
     nini=com_edge2D%rptr(n)
     nend=com_edge2D%rptr(n+1) - 1
      nc=0
      DO nh=nini, nend
         DO nz=1, nl-1
	 nc=nc+1
	 edge_array3D(nz,com_edge2D%rlist(nh))=r_buff_edge3D(n)%array(nc)
         END DO
      END DO	  
  end do

end subroutine exchange_edge3D
!==========================================================================

subroutine exchange_edge2D(edge_array2D)
  use o_MESH
  use g_PARSUP 
  implicit none
  INTEGER  :: sreq(maxPEnum)
  INTEGER  :: rreq(maxPEnum)
  INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
  INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
  integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
  real(kind=8) :: edge_array2D(myDim_edge2D+eDim_edge2D) 

  sn=com_edge2D%sPEnum
  rn=com_edge2D%rPEnum
  ! Put data to be communicated into send buffer 
 do n=1, sn
     nini=com_edge2D%sptr(n)
     nend=com_edge2D%sptr(n+1) - 1
     s_buff_edge2D(n)%array=edge_array2D(com_edge2D%slist(nini:nend))
  end do
  do n=1, sn
     dest=com_edge2D%sPE(n)
     nini=com_edge2D%sptr(n)
     offset=com_edge2D%sptr(n+1) - nini

     call MPI_ISEND(s_buff_edge2D(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
          MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1,rn 
     source=com_edge2D%rPE(n)
     nini=com_edge2D%rptr(n)
     offset=com_edge2D%rptr(n+1) - nini

     call MPI_IRECV(r_buff_edge2D(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
          source, MPI_COMM_WORLD, rreq(n), MPIerr) 
  end do
  call MPI_WAITALL(sn,sreq,sstat, MPIerr)
  call MPI_WAITALL(rn,rreq,rstat, MPIerr)

  ! Put received data to their destination

  do n=1, rn
     nini=com_edge2D%rptr(n)
     nend=com_edge2D%rptr(n+1) - 1
     edge_array2D(com_edge2D%rlist(nini:nend))=r_buff_edge2D(n)%array
  end do

end subroutine exchange_edge2D
!=============================================================================
subroutine exchange_elem3D(elem_array3D)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
 
 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 integer  :: nz, nh, nc, nl1
 real(kind=WP) :: elem_array3D(:,:)
 type(com_struct)   :: com_elem2D_temp
 type(com_array),dimension(:),allocatable    :: s_buff_temp,r_buff_temp
  nl1 = ubound(elem_array3D,1)
  if (ubound(elem_array3D,2)>myDim_elem2D+eDim_elem2D) then
        allocate(s_buff_temp(com_elem2D_full%sPEnum),  r_buff_temp(com_elem2D_full%rPEnum))
	com_elem2D_temp=com_elem2D_full
        s_buff_temp=s_buff_elem3D_full
        r_buff_temp=r_buff_elem3D_full
  else
        allocate(s_buff_temp(com_elem2D%sPEnum),  r_buff_temp(com_elem2D%rPEnum))
	com_elem2D_temp=com_elem2D
        s_buff_temp=s_buff_elem3D
        r_buff_temp=r_buff_elem3D
  endif
  sn=com_elem2D_temp%sPEnum
  rn=com_elem2D_temp%rPEnum
  ! Put data to be communicated into send buffer 
  Do n=1,com_elem2D_temp%sPEnum
   nini=com_elem2D_temp%sptr(n)
   nend=com_elem2D_temp%sptr(n+1)-1
   nc=0
   DO nh=nini, nend 
      DO nz=1,nl1
      nc=nc+1
      s_buff_temp(n)%array(nc)=elem_array3D(nz,com_elem2D_temp%slist(nh))
      END DO
   END DO   
  end do  
 
  DO n=1, sn
     dest=com_elem2D_temp%sPE(n)
     nini=com_elem2D_temp%sptr(n)
     offset=(com_elem2D_temp%sptr(n+1) - nini)*nl1
     
     
     call MPI_ISEND(s_buff_temp(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1,rn   
     source=com_elem2D_temp%rPE(n)
     nini=com_elem2D_temp%rptr(n)
     offset=(com_elem2D_temp%rptr(n+1) - nini)*nl1
     
      
     call MPI_IRECV(r_buff_temp(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
  Do n=1,com_elem2D_temp%rPEnum
   nini=com_elem2D_temp%rptr(n)
   nend=com_elem2D_temp%rptr(n+1)-1
   nc=0
   DO nh=nini, nend
      DO nz=1, nl1
      nc=nc+1
      elem_array3D(nz,com_elem2D_temp%rlist(nh))=r_buff_temp(n)%array(nc)
      END DO
   END DO   
  end do
  deallocate(s_buff_temp,r_buff_temp) 
  !deallocate(com_elem2D_temp)
END SUBROUTINE exchange_elem3D
!=============================================================================
subroutine exchange_elem3D_full(elem_array3D)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
 
 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 integer  :: nz, nh, nc
 real(kind=WP) :: elem_array3D(nl-1,elem2D) 
  sn=com_elem2D_full%sPEnum
  rn=com_elem2D_full%rPEnum
  ! Put data to be communicated into send buffer 
 
  Do n=1,com_elem2D_full%sPEnum
   nini=com_elem2D_full%sptr(n)
   nend=com_elem2D_full%sptr(n+1)-1
   nc=0
   DO nh=nini, nend 
      DO nz=1,nl-1
      nc=nc+1
      s_buff_elem3D_full(n)%array(nc)=elem_array3D(nz,com_elem2D_full%slist(nh))
      END DO
   END DO   
  end do  
 
  DO n=1, sn
     dest=com_elem2D_full%sPE(n)
     nini=com_elem2D_full%sptr(n)
     offset=(com_elem2D_full%sptr(n+1) - nini)*(nl-1)
     
     
     call MPI_ISEND(s_buff_elem3D_full(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1,rn   
     source=com_elem2D_full%rPE(n)
     nini=com_elem2D_full%rptr(n)
     offset=(com_elem2D_full%rptr(n+1) - nini)*(nl-1)
     
      
     call MPI_IRECV(r_buff_elem3D_full(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
  Do n=1,com_elem2D_full%rPEnum
   nini=com_elem2D_full%rptr(n)
   nend=com_elem2D_full%rptr(n+1)-1
   nc=0
   DO nh=nini, nend
      DO nz=1, nl-1
      nc=nc+1
      elem_array3D(nz,com_elem2D_full%rlist(nh))=r_buff_elem3D_full(n)%array(nc)
      END DO
   END DO   
  end do 
   
END SUBROUTINE exchange_elem3D_full
!========================================================================
subroutine exchange_elem2D(elem_array2D)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 3D elemental fields
! stored in (vertical, horizontal) format
 
 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 real(kind=WP) :: elem_array2D(:) 
 type(com_array),dimension(:),allocatable :: s_buff_temp,  r_buff_temp
 type(com_struct):: com_elem2D_temp
  if (ubound(elem_array2D,1)<=myDim_elem2D+eDim_elem2D) then
        com_elem2D_temp=com_elem2D
        allocate(s_buff_temp(com_elem2D%sPEnum),  r_buff_temp(com_elem2D%rPEnum))
        s_buff_temp=s_buff_elem2D
        r_buff_temp=r_buff_elem2D
  else
        com_elem2D_temp=com_elem2D_full
        allocate(s_buff_temp(com_elem2D_full%sPEnum),  r_buff_temp(com_elem2D_full%rPEnum))
        s_buff_temp=s_buff_elem2D_full
        r_buff_temp=r_buff_elem2D_full
  endif 
  sn=com_elem2D_temp%sPEnum
  rn=com_elem2D_temp%rPEnum
  ! Put data to be communicated into send buffer 
 
  Do n=1,com_elem2D_temp%sPEnum
   nini=com_elem2D_temp%sptr(n)
   nend=com_elem2D_temp%sptr(n+1)-1
   s_buff_temp(n)%array=elem_array2D(com_elem2D_temp%slist(nini:nend))
  end do  
 
  DO n=1, sn
     dest=com_elem2D_temp%sPE(n)
     nini=com_elem2D_temp%sptr(n)
     offset=(com_elem2D_temp%sptr(n+1) - nini)
     
     
     call MPI_ISEND(s_buff_temp(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1,rn   
     source=com_elem2D_temp%rPE(n)
     nini=com_elem2D_temp%rptr(n)
     offset=(com_elem2D_temp%rptr(n+1) - nini)
     
      
     call MPI_IRECV(r_buff_temp(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
  Do n=1,com_elem2D_temp%rPEnum
   nini=com_elem2D_temp%rptr(n)
   nend=com_elem2D_temp%rptr(n+1)-1
   elem_array2D(com_elem2D_temp%rlist(nini:nend))=r_buff_temp(n)%array
  end do 
  deallocate(s_buff_temp,r_buff_temp)
END SUBROUTINE exchange_elem2D
! ========================================================================
subroutine exchange_elem2D_i(elem_array2D)
!Exchange with ALL(!) the neighbours
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 2D elemental fields
! Integer arrays
 
 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 integer  :: elem_array2D(elem2D) 
 !type(com_array_i), allocatable :: s_buff_temp(:),  r_buff_temp(:)
 !type(com_struct), allocatable:: com_elem2D_temp


 ! if (ubound(elem_array2D,1)>myDim_elem2D+eDim_elem2D) then
 !       com_elem2D_temp=com_elem2D_temp
 !       s_buff_temp=s_buff_temp
 !       r_buff_temp=r_buff_temp
 ! else
 !       com_elem2D_temp=com_elem2D
 !       s_buff_temp=s_buff_elem2D
 !       r_buff_temp=r_buff_elem2D
 ! endif
 
  sn=com_elem2D_full%sPEnum
  rn=com_elem2D_full%rPEnum
  ! Put data to be communicated into send buffer 
 
  Do n=1,com_elem2D_full%sPEnum
   nini=com_elem2D_full%sptr(n)
   nend=com_elem2D_full%sptr(n+1)-1
   s_buff_elem2D_full_i(n)%array=elem_array2D(com_elem2D_full%slist(nini:nend))
  end do  
 
  DO n=1, sn
     dest=com_elem2D_full%sPE(n)
     nini=com_elem2D_full%sptr(n)
     offset=(com_elem2D_full%sptr(n+1) - nini)
     
     
     call MPI_ISEND(s_buff_elem2D_full_i(n)%array, offset, MPI_integer, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1,rn   
     source=com_elem2D_full%rPE(n)
     nini=com_elem2D_full%rptr(n)
     offset=(com_elem2D_full%rptr(n+1) - nini)
     
      
     call MPI_IRECV(r_buff_elem2D_full_i(n)%array, offset, MPI_integer, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
  Do n=1,com_elem2D_full%rPEnum
   nini=com_elem2D_full%rptr(n)
   nend=com_elem2D_full%rptr(n+1)-1
   elem_array2D(com_elem2D_full%rlist(nini:nend))=r_buff_elem2D_full_i(n)%array
  end do 
   
END SUBROUTINE exchange_elem2D_i
!=============================================================================

subroutine exchange_e2D(elem_array2D, ne)
USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

! General version of the communication routine for 2D elemental fields
! stored in (vertical, horizontal) format
 
 INTEGER  :: sreq(npes)
 INTEGER  :: rreq(npes)
 INTEGER  :: sstat(MPI_STATUS_SIZE,npes)
 INTEGER  :: rstat(MPI_STATUS_SIZE,npes)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 integer  :: nz, nh, nc, ne
 real(kind=8) :: elem_array2D(ne,elem2D) 
 type(com_array), allocatable :: s_buff(:), r_buff(:)
 
      allocate(s_buff(com_elem2D%sPEnum))
      allocate(r_buff(com_elem2D%rPEnum))
      do n=1, com_elem2D%sPEnum
        offset=com_elem2D%sptr(n+1) - com_elem2D%sptr(n)
        allocate(s_buff(n)%array(offset*ne))
      end do
      do n=1, com_elem2D%rPEnum
        offset=com_elem2D%rptr(n+1) - com_elem2D%rptr(n)
        allocate(r_buff(n)%array(offset*ne))
      end do
 
  sn=com_elem2D%sPEnum
  rn=com_elem2D%rPEnum
  ! Put data to be communicated into send buffer 
  
  do n=1, sn
     nini=com_elem2D%sptr(n)
     nend=com_elem2D%sptr(n+1) - 1
     nc=0
   DO nh=nini, nend 
      DO nz=1,ne
      nc=nc+1
      s_buff(n)%array(nc)=elem_array2D(nz,com_elem2D%slist(nh))
      END DO
   END DO   
  end do
  
  DO n=1, sn
     dest=com_elem2D%sPE(n)
     nini=com_elem2D%sptr(n)
     offset=(com_elem2D%sptr(n+1) - nini)
     
     
     call MPI_ISEND(s_buff(n)%array, ne*offset, MPI_DOUBLE_PRECISION, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1, rn   
     source=com_elem2D%rPE(n)
     nini=com_elem2D%rptr(n)
     offset=(com_elem2D%rptr(n+1) - nini)
     
      
     call MPI_IRECV(r_buff(n)%array, ne*offset, MPI_DOUBLE_PRECISION, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(sn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
   do n=1, rn
     nini=com_elem2D%rptr(n)
     nend=com_elem2D%rptr(n+1) - 1
     nc=0
      DO nh=nini, nend
       DO nz=1, ne
       nc=nc+1
       elem_array2D(nz,com_elem2D%rlist(nh))=r_buff(n)%array(nc)
       END DO
      END DO
   end do 
 ! ===============  
! Deallocate the buffers
! ===============
      do n=com_elem2D%rPEnum,1,-1
        deallocate(r_buff(n)%array)
      end do
      do n=com_elem2D%sPEnum,1,-1
        deallocate(s_buff(n)%array)
      end do
      deallocate(r_buff, s_buff)    
  
  
END SUBROUTINE exchange_e2D
! ========================================================================
! ========================================================================
subroutine exchange_n3D(nod_array3D,ne,vsize)
USE o_MESH
USE g_PARSUP
USE o_PARAM 
IMPLICIT NONE

! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
 
 INTEGER  :: sreq(maxPEnum)
 INTEGER  :: rreq(maxPEnum)
 INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
 INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
 integer  :: n, sn, rn, dest, nini, nend, offset, source,tag
 integer  :: nz, nh, nc, nn
 integer, intent(in)          :: ne,vsize 
 real(kind=WP), intent(inout) :: nod_array3D(ne,vsize,myDim_nod2D+eDim_nod2D)
 type(com_array), allocatable :: s_buff(:),  r_buff(:)

  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum
  
  ! ============
  ! Allocate buffers
  ! ============
     allocate(s_buff(com_nod2D%sPEnum),r_buff(com_nod2D%rPEnum))
     do n=1, com_nod2D%sPEnum
        offset=com_nod2D%sptr(n+1) - com_nod2D%sptr(n)
	allocate(s_buff(n)%array(offset*ne*vsize))
     end do
     do n=1, com_nod2D%rPEnum
        offset=com_nod2D%rptr(n+1) - com_nod2D%rptr(n)
	allocate(r_buff(n)%array(offset*ne*vsize))
     end do
  
  ! Put data to be communicated into send buffer 
 
  Do n=1,com_nod2D%sPEnum
   nini=com_nod2D%sptr(n)
   nend=com_nod2D%sptr(n+1)-1
   nc=0
   DO nh=nini, nend 
      DO nz=1,vsize
      DO nn=1,ne
      nc=nc+1
      s_buff(n)%array(nc)=nod_array3D(nn,nz,com_nod2D%slist(nh))
      END DO
      END DO
   END DO   
  end do  
 
  DO n=1, sn
     dest=com_nod2D%sPE(n)
     nini=com_nod2D%sptr(n)
     offset=(com_nod2D%sptr(n+1) - nini)*ne*vsize
     
     
     call MPI_ISEND(s_buff(n)%array, offset, MPI_DOUBLE_PRECISION, dest, mype, & 
               MPI_COMM_WORLD, sreq(n), MPIerr)
  END DO
  DO n=1,rn    
     source=com_nod2D%rPE(n)
     nini=com_nod2D%rptr(n)
     offset=(com_nod2D%rptr(n+1) - nini)*ne*vsize
     
      
     call MPI_IRECV(r_buff(n)%array, offset, MPI_DOUBLE_PRECISION, source, &
               source, MPI_COMM_WORLD, rreq(n), MPIerr) 
            
  END DO 
  
     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
  
  ! Put received data to their destination
  Do n=1,com_nod2D%rPEnum
   nini=com_nod2D%rptr(n)
   nend=com_nod2D%rptr(n+1)-1
   nc=0
   DO nh=nini, nend
      DO nz=1, vsize
       DO nn=1,ne
       nc=nc+1
       nod_array3D(nn,nz,com_nod2D%rlist(nh))=r_buff(n)%array(nc)
       END DO
      END DO
   END DO   
  end do 
! ===============  
! Deallocate the buffers
! ===============
      do n=com_nod2D%rPEnum,1,-1
        deallocate(r_buff(n)%array)
      end do
      do n=com_nod2D%sPEnum,1,-1
        deallocate(s_buff(n)%array)
      end do
      deallocate(r_buff, s_buff)   
END SUBROUTINE exchange_n3D
! ========================================================================





! ========================================================================
! Broadcast routines
! Many because of different sizes.
! ========================================================================
subroutine broadcast_nod3D(arr3D, arr3Dglobal)

! Make nodal information available to all PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER :: ireals, nz, counter,nl1
integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:) ::  isendbuf, irecvbuf

!real(kind=WP) ::  arr3D(nl-1,myDim_nod2D+eDim_nod2D)
!real(kind=WP) ::  arr3Dglobal(nl-1,nod2D)
real(kind=WP) ::  arr3D(:,:)
real(kind=WP) ::  arr3Dglobal(:,:)
real(kind=WP), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf

nl1=ubound(arr3D,1)
IF ( mype == 0 ) THEN
    if (npes>1) then
    arr3Dglobal(:,myList_nod2D(1:myDim_nod2D))=arr3D(:,1:myDim_nod2D)
    end if
    DO  n = 1, npes-1

       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, MPI_COMM_WORLD, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(recvbuf(1:nTS*nl1), irecvbuf(1:nTS) )
       CALL MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, MPI_COMM_WORLD, status, MPIerr )
       CALL MPI_RECV( recvbuf(1), nTS*nl1, MPI_DOUBLE_PRECISION, sender, &
                      2, MPI_COMM_WORLD, status, MPIerr )

       counter=0
       DO i = 1, nTS
          DO nz=1, nl1
	  counter=counter+1
          arr3Dglobal(nz,irecvbuf(i)) = recvbuf(counter)
	  ENDDO
       ENDDO
       DEALLOCATE( recvbuf, irecvbuf )

    ENDDO

ELSE

    ALLOCATE( sendbuf(1:myDim_nod2D*nl1), isendbuf(1:myDim_nod2D) )
    counter=0
    DO n = 1, myDim_nod2D
       isendbuf(n) = myList_nod2D(n)
       DO nz=1, nl1
       counter=counter+1
       sendbuf(counter)  = arr3D(nz,n)
       ENDDO
    ENDDO
    CALL MPI_SEND( myDim_nod2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( isendbuf(1), myDim_nod2D, MPI_INTEGER, 0, 1, &
                   MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( sendbuf(1), myDim_nod2D*nl1, MPI_DOUBLE_PRECISION, &
                   0, 2, MPI_COMM_WORLD, MPIerr )
    DEALLOCATE( sendbuf, isendbuf )

ENDIF
 DO n=1,nl1
 CALL MPI_BCAST( arr3Dglobal(n,:), nod2d, MPI_DOUBLE_PRECISION, 0, &
                 MPI_COMM_WORLD, MPIerr)
 END DO
end subroutine broadcast_nod3D
!
!============================================================================
subroutine broadcast_nod2D(arr2D, arr2Dglobal)
! A 2D version of the previous routine
use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER :: ireals
integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:) ::  isendbuf, irecvbuf

real(kind=WP) ::  arr2D(myDim_nod2D+eDim_nod2D)
real(kind=WP) ::  arr2Dglobal(nod2D)
real(kind=WP), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf

IF ( mype == 0 ) THEN
    if (npes>1) then
    arr2Dglobal(myList_nod2D(1:myDim_nod2D))=arr2D(1:myDim_nod2D)
    end if
    DO  n = 1, npes-1

       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, MPI_COMM_WORLD, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE( recvbuf(1:nTS), irecvbuf(1:nTS) )
       CALL MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, MPI_COMM_WORLD, status, MPIerr )
       CALL MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
                      2, MPI_COMM_WORLD, status, MPIerr )

       DO i = 1, nTS
          arr2Dglobal(irecvbuf(i)) = recvbuf(i)
       ENDDO
       DEALLOCATE( recvbuf, irecvbuf )

    ENDDO

ELSE

    ALLOCATE( sendbuf(1:myDim_nod2d), isendbuf(1:myDim_nod2D) )
    DO n = 1, myDim_nod2D
       isendbuf(n) = myList_nod2D(n)
       sendbuf(n)  = arr2D(n)
    ENDDO
    CALL MPI_SEND( myDim_nod2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( isendbuf(1), myDim_nod2D, MPI_INTEGER, 0, 1, &
                   MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( sendbuf(1), myDim_nod2D, MPI_DOUBLE_PRECISION, &
                   0, 2, MPI_COMM_WORLD, MPIerr )
    DEALLOCATE( sendbuf, isendbuf )

ENDIF

 CALL MPI_BCAST(arr2Dglobal, nod2d, MPI_DOUBLE_PRECISION, 0, &
                 MPI_COMM_WORLD, MPIerr)

end subroutine broadcast_nod2D
!===================================================================
subroutine broadcast_edge2D(arr2D, arr2Dglobal)
! A 2D version of the previous routine
use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER :: ireals
integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:) ::  isendbuf, irecvbuf

real(kind=WP) ::  arr2D(myDim_edge2D+eDim_edge2D)
real(kind=WP) ::  arr2Dglobal(nod2D)
real(kind=WP), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf

IF ( mype == 0 ) THEN
    if (npes>1) then
    arr2Dglobal(myList_edge2D(1:myDim_edge2D))=arr2D(1:myDim_edge2D)
    end if
    DO  n = 1, npes-1
    
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, MPI_COMM_WORLD, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE( recvbuf(1:nTS), irecvbuf(1:nTS) )
       CALL MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, MPI_COMM_WORLD, status, MPIerr )
       CALL MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
                      2, MPI_COMM_WORLD, status, MPIerr )

       DO i = 1, nTS
          arr2Dglobal(irecvbuf(i)) = recvbuf(i)
       ENDDO
       DEALLOCATE( recvbuf, irecvbuf )

    ENDDO

ELSE

    ALLOCATE( sendbuf(1:myDim_edge2d), isendbuf(1:myDim_edge2D) )
    DO n = 1, myDim_edge2D
       isendbuf(n) = myList_edge2D(n)
       sendbuf(n)  = arr2D(n)
    ENDDO
    CALL MPI_SEND( myDim_edge2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( isendbuf(1), myDim_edge2D, MPI_INTEGER, 0, 1, &
                   MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( sendbuf(1), myDim_edge2D, MPI_DOUBLE_PRECISION, &
                   0, 2, MPI_COMM_WORLD, MPIerr )
    DEALLOCATE( sendbuf, isendbuf )

ENDIF
    
 CALL MPI_BCAST(arr2Dglobal, edge2d, MPI_DOUBLE_PRECISION, 0, &
                 MPI_COMM_WORLD, MPIerr)

end subroutine broadcast_edge2D 
!===================================================================
subroutine broadcast_elem3D(arr3D, arr3Dglobal)
!Makes elemental information available to all PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER :: ireals, nz, counter
integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:) ::  isendbuf, irecvbuf

real(kind=8) ::  arr3D(nl-1,myDim_elem2D+eDim_elem2D)
real(kind=8) ::  arr3Dglobal(nl-1,elem2D)
real(kind=8), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf

IF ( mype == 0 ) THEN
    if (npes>1) then
    arr3Dglobal(:,myList_elem2D(1:myDim_elem2D))=arr3D(:,1:myDim_elem2D)
    end if
    DO  n = 1, npes-1
       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, MPI_COMM_WORLD, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(recvbuf(1:nTS*(nl-1)), irecvbuf(1:nTS) )
       CALL MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, MPI_COMM_WORLD, status, MPIerr )
       CALL MPI_RECV( recvbuf(1), nTS*(nl-1), MPI_DOUBLE_PRECISION, sender, &
                      2, MPI_COMM_WORLD, status, MPIerr )

       counter=0
       DO i = 1, nTS
          DO nz=1, nl-1
	  counter=counter+1
          arr3Dglobal(nz,irecvbuf(i)) = recvbuf(counter)
	  ENDDO
       ENDDO
       DEALLOCATE( recvbuf, irecvbuf )

    ENDDO

ELSE

    ALLOCATE( sendbuf(1:myDim_elem2D*(nl-1)), isendbuf(1:myDim_elem2D) )
    counter=0
    DO n = 1, myDim_elem2D
       isendbuf(n) = myList_elem2D(n)
       DO nz=1, nl-1
       counter=counter+1
       sendbuf(counter)  = arr3D(nz,n)
       ENDDO
    ENDDO
    CALL MPI_SEND( myDim_elem2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( isendbuf(1), myDim_elem2D, MPI_INTEGER, 0, 1, &
                   MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( sendbuf(1), myDim_elem2D*(nl-1), MPI_DOUBLE_PRECISION, &
                   0, 2, MPI_COMM_WORLD, MPIerr )
    DEALLOCATE( sendbuf, isendbuf )

ENDIF
 DO n=1,nl-1 
 CALL MPI_BCAST( arr3Dglobal(n,:), elem2d, MPI_DOUBLE_PRECISION, 0, &
                 MPI_COMM_WORLD, MPIerr)
 END DO 
end subroutine broadcast_elem3D
!
!===================================================================
subroutine broadcast_elem2D(arr2D, arr2D_global)
! PE0 collects all information. Only its copy of arrays will be
! correct
use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER :: ireals
integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:) ::  isendbuf, irecvbuf

real(kind=8) ::  arr2D(myDim_elem2D+eDim_elem2D)
real(kind=8) ::  arr2D_global(elem2D)
real(kind=8), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf

IF ( mype == 0 ) THEN
    arr2D_global(myList_elem2D(1:myDim_elem2D))=arr2D(1:myDim_elem2D)
    DO  n = 1, npes-1

       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     0, MPI_COMM_WORLD, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE( recvbuf(1:nTS), irecvbuf(1:nTS) )
       CALL MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
                      1, MPI_COMM_WORLD, status, MPIerr )
       CALL MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
                      2, MPI_COMM_WORLD, status, MPIerr )

       DO i = 1, nTS
          arr2D_global(irecvbuf(i)) = recvbuf(i)
       ENDDO
       DEALLOCATE( recvbuf, irecvbuf )

    ENDDO

ELSE

    ALLOCATE( sendbuf(1:myDim_elem2d), isendbuf(1:myDim_elem2D) )
    DO n = 1, myDim_elem2D
       isendbuf(n) = myList_elem2D(n)
       sendbuf(n)  = arr2D(n)
    ENDDO
    CALL MPI_SEND( myDim_elem2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( isendbuf(1), myDim_elem2D, MPI_INTEGER, 0, 1, &
                   MPI_COMM_WORLD, MPIerr )
    CALL MPI_SEND( sendbuf(1), myDim_elem2D, MPI_DOUBLE_PRECISION, &
                   0, 2, MPI_COMM_WORLD, MPIerr )
    DEALLOCATE( sendbuf, isendbuf )

ENDIF

 CALL MPI_BCAST( arr2D_global, elem2d, MPI_DOUBLE_PRECISION, 0, &
                 MPI_COMM_WORLD, MPIerr)

end subroutine broadcast_elem2D
!============================================================================
subroutine gather_nod3D(arr3D, arr3D_global)

! Make nodal information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(kind=WP)  ::  arr3D(:,:)
real(kind=WP)  ::  arr3D_global(:,:)
real(kind=WP), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D

CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,nod2D))

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_WORLD, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_nod2D*nl1, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, MPIerr )
   
ENDIF


end subroutine gather_nod3D
!==============================================
subroutine gather_nod2D(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(kind=WP)  ::  arr2D(:)
real(kind=WP)  ::  arr2D_global(:)
real(kind=WP)  :: recvbuf(nod2D)
integer        :: req(npes-1)
integer        :: start, n2D


CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then

      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_WORLD, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, MPIerr )
   
ENDIF


end subroutine gather_nod2D
!==============================================
!============================================================================
subroutine gather_elem3D(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(kind=WP)  ::  arr3D(:,:)
real(kind=WP)  ::  arr3D_global(:,:)
real(kind=WP), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D, ende, err_alloc
integer        :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( mype == 0 ) THEN
   
   if (npes>1) then
! 
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl1
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_WORLD, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   

      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                        = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_elem2D*nl1, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, MPIerr )
   
ENDIF


end subroutine gather_elem3D
!==============================================
subroutine gather_elem2D(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(kind=WP)  ::  arr2D(:)
real(kind=WP)  ::  arr2D_global(:)
real(kind=WP)  :: recvbuf(elem2D)
integer        :: req(npes-1)
integer        :: start, e2D


CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_WORLD, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_elem2D, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, MPIerr )
   
ENDIF


end subroutine gather_elem2D

!============================================================================
subroutine gather_real4_nod3D(arr3D, arr3D_global)

! Make nodal information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(kind=WP)  ::  arr3D(:,:)
real(kind=4)   ::  arr3D_global(:,:)
real(kind=4), allocatable :: recvbuf(:,:)
real(kind=4), allocatable :: sendbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D, ierr


CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,nod2D))
!      if (ierr/=0) print *,'allocate knall auf pe=0 in nod3D_real4',ierr

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_REAL, n, 2, MPI_COMM_WORLD, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                     = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   allocate(sendbuf(nl1,myDim_nod2D))
   sendbuf(1:nl1,1:myDim_nod2D) = arr3D(1:nl1,1:myDim_nod2D)
   
   call MPI_SEND(sendbuf, myDim_nod2D*nl1, MPI_REAL, 0, 2, MPI_COMM_WORLD, MPIerr )
   deallocate(sendbuf)
   
ENDIF


end subroutine gather_real4_nod3D
!==============================================
subroutine gather_real4_nod2D(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(kind=WP)  ::  arr2D(:)
real(kind=4)   ::  arr2D_global(:)
real(kind=4)   :: recvbuf(nod2D), sendbuf(myDim_nod2D)
integer        :: req(npes-1)
integer        :: start, n2D

! Consider MPI-datatypes to recv directly into arr2D_global!

CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
IF ( mype == 0 ) THEN
   
   if (npes>1) then

      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_REAL, n, 2, MPI_COMM_WORLD, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   sendbuf(1:myDim_nod2D) = real(arr2D(1:myDim_nod2D),4)

   call MPI_SEND(sendbuf, myDim_nod2D, MPI_REAL, 0, 2, MPI_COMM_WORLD, MPIerr )
   
ENDIF


end subroutine gather_real4_nod2D
!==============================================
!============================================================================
subroutine gather_real4_elem3D(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(kind=WP)  ::  arr3D(:,:)
real(kind=4)  ::  arr3D_global(:,:)
real(kind=4), allocatable :: recvbuf(:,:)
real(kind=4), allocatable :: sendbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D


CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_REAL, n, 2, MPI_COMM_WORLD, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                     = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   allocate(sendbuf(nl1,myDim_elem2D))
   sendbuf(1:nl1,1:myDim_elem2D) = arr3D(1:nl1,1:myDim_elem2D)
   
   call MPI_SEND(sendbuf, myDim_elem2D*nl1, MPI_REAL, 0, 2, MPI_COMM_WORLD, MPIerr )
   deallocate(sendbuf)
ENDIF


end subroutine gather_real4_elem3D
!==============================================
subroutine gather_real4_elem2D(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(kind=WP)  ::  arr2D(:)
real(kind=4)  ::  arr2D_global(:)
real(kind=4)  :: recvbuf(elem2D)
real(kind=4)  :: sendbuf(myDim_elem2D)
integer        :: req(npes-1)
integer        :: start, e2D



CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_REAL, n, 2, MPI_COMM_WORLD, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   sendbuf(1:myDim_elem2D) = real(arr2D(1:myDim_elem2D),4)
   call MPI_SEND(sendbuf, myDim_elem2D, MPI_REAL, 0, 2, MPI_COMM_WORLD, MPIerr )
   
ENDIF


end subroutine gather_real4_elem2D
!==============================================

end module g_comm



module g_comm_auto
use g_comm
implicit none
interface exchange_nod
      module procedure exchange_nod2D
      module procedure exchange_nod2D_i
      module procedure exchange_nod3D
!      module procedure exchange_nod3D_full
end interface exchange_nod

interface exchange_edge
      module procedure exchange_edge2D
      module procedure exchange_edge3D
end interface exchange_edge

interface exchange_elem
!TODO overhead in 2d exchange routine
      module procedure exchange_elem3D
!      module procedure exchange_elem3D_full
      module procedure exchange_elem2d
      module procedure exchange_elem2d_i
end interface exchange_elem

interface broadcast_nod
      module procedure broadcast_nod3D
      module procedure broadcast_nod2D
end interface broadcast_nod

interface broadcast_elem
      module procedure broadcast_elem3D
      module procedure broadcast_elem2D
end interface broadcast_elem

interface gather_nod
      module procedure gather_nod3D
      module procedure gather_nod2D
      module procedure gather_real4_nod3D
      module procedure gather_real4_nod2D
end interface gather_nod

interface gather_elem
      module procedure gather_elem3D
      module procedure gather_elem2D
      module procedure gather_real4_elem3D
      module procedure gather_real4_elem2D
end interface gather_elem


private  ! hides items not listed on public statement 
public :: exchange_nod,exchange_edge,exchange_elem,broadcast_nod,broadcast_elem, &
          gather_nod, gather_elem
end module g_comm_auto
