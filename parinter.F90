#define __MYFILE__ 'parinter.F90'
MODULE parinter

#if defined key_mpp_mpi
   USE mpi
#endif
   USE scripremap
   USE scrippar
   USE nctools

   IMPLICIT NONE

   ! Type to contains interpolation information
   ! (like what is in scripremaptype) and message
   ! passing information

   TYPE parinterinfo  
      ! Number of local links
      INTEGER :: num_links
      ! Destination side
      INTEGER, POINTER, DIMENSION(:) :: dst_address
      ! Source addresses and work array
      INTEGER, POINTER, DIMENSION(:) :: src_address
      ! Local remap matrix
      REAL(scripdp), POINTER, DIMENSION(:,:) :: remap_matrix
      ! Message passing information
      ! Array of local addresses for send buffer
      ! packing
      INTEGER, POINTER, DIMENSION(:) :: send_address
      ! Sending bookkeeping
      INTEGER :: nsendtot                
      INTEGER, POINTER, DIMENSION(:) :: nsend,nsdisp
      ! Receiving bookkeeping
      INTEGER :: nrecvtot
      INTEGER, POINTER, DIMENSION(:) :: nrecv,nrdisp
   END TYPE parinterinfo

CONTAINS

   SUBROUTINE parinter_init( mype, nproc, mpi_comm, &
      & nsrclocpoints, nsrcglopoints, srcmask, srcgloind,  &
      & ndstlocpoints, ndstglopoints, dstmask, dstgloind, &
      & remap, pinfo, lcommout, commoutprefix, iunit )

      ! Setup interpolation based on SCRIP format weights in
      ! remap and the source/destination grids information.

      ! Procedure:

      !   1) A global SCRIP remapping file is read on all processors.
      !   2) Find local destination points in the global grid.
      !   3) Find which processor needs source data and setup buffer
      !      information for sending data.
      !   4) Construct new src remapping for buffer received

      ! All information is stored in the TYPE(parinterinfo) output
      ! data type

      ! Input arguments.

      ! Message passing information
      INTEGER, INTENT(IN) :: mype, nproc, mpi_comm
      ! Source grid local and global number of grid points
      INTEGER, INTENT(IN) :: nsrclocpoints, nsrcglopoints
      ! Source integer mask (0/1) for SCRIP compliance
      INTEGER, INTENT(IN), DIMENSION(nsrclocpoints) :: srcmask
      ! Source global addresses of each local grid point
      INTEGER, INTENT(IN), DIMENSION(nsrclocpoints) :: srcgloind
      ! Destination grid local and global number of grid points
      INTEGER, INTENT(IN) :: ndstlocpoints, ndstglopoints
      ! Destination integer mask (0/1) for SCRIP compliance
      INTEGER, INTENT(IN), DIMENSION(ndstlocpoints) :: dstmask
      ! Destination global addresses of each local grid point
      INTEGER, INTENT(IN), DIMENSION(ndstlocpoints) :: dstgloind
      ! SCRIP remapping data
      TYPE(scripremaptype) :: remap
      ! Switch for output communication patterns
      LOGICAL :: lcommout
      CHARACTER(len=*) :: commoutprefix
      ! Unit to use for output
      INTEGER :: iunit

      ! Output arguments

      ! Interpolation and message passing information
      TYPE(parinterinfo), INTENT(OUT) :: pinfo

      ! Local variable

      ! Variable for glocal <-> local address/pe information
      INTEGER, DIMENSION(nsrcglopoints) :: ilsrcmppmap, ilsrclocind
      INTEGER, DIMENSION(nsrcglopoints) :: igsrcmppmap, igsrclocind
      INTEGER, DIMENSION(ndstglopoints) :: ildstmppmap, ildstlocind
      INTEGER, DIMENSION(ndstglopoints) :: igdstmppmap, igdstlocind
      INTEGER, DIMENSION(nsrcglopoints) :: isrcpe,isrcpetmp
      INTEGER, DIMENSION(nsrcglopoints) :: isrcaddtmp
      INTEGER, DIMENSION(0:nproc-1)     :: isrcoffset
      INTEGER, DIMENSION(nproc)         :: isrcno, isrcoff, isrccur
      INTEGER, DIMENSION(nproc)         :: ircvoff, ircvcur
      INTEGER, DIMENSION(:), ALLOCATABLE :: isrctot, ircvtot

      ! Misc variable
      INTEGER :: i,n,pe
      INTEGER :: istatus
      CHARACTER(len=256) :: cdfile

      ! Check that masks are consistent.

      ! Remark: More consistency tests between remapping information
      ! and input argument could be code, but for now we settle
      ! for checking the masks.

      ! Source grid

      DO i=1,nsrclocpoints
         IF (srcmask(i)/=remap%src%grid_imask(srcgloind(i))) THEN
            WRITE(iunit,*)'Source imask is inconsistent at '
            WRITE(iunit,*)'global index = ',srcgloind(i)
            WRITE(iunit,*)'Source mask  = ',srcmask(i)
            WRITE(iunit,*)'Remap  mask  = ',remap%src%grid_imask(srcgloind(i))
            WRITE(iunit,*)'Latitude     = ',remap%src%grid_center_lat(srcgloind(i))
            WRITE(iunit,*)'Longitude    = ',remap%src%grid_center_lon(srcgloind(i))
            CALL flush(iunit)
            CALL abort
         ENDIF
      ENDDO

      ! Destination grid
      
      DO i=1,ndstlocpoints
         IF (dstmask(i)/=remap%dst%grid_imask(dstgloind(i))) THEN
            WRITE(iunit,*)'Destination imask is inconsistent at '
            WRITE(iunit,*)'global index = ',dstgloind(i)
            WRITE(iunit,*)'Destin mask  = ',dstmask(i)
            WRITE(iunit,*)'Remap  mask  = ',remap%dst%grid_imask(dstgloind(i))
            WRITE(iunit,*)'Latitude     = ',remap%dst%grid_center_lat(dstgloind(i))
            WRITE(iunit,*)'Longitude    = ',remap%dst%grid_center_lon(dstgloind(i))
            CALL flush(iunit)
            CALL abort
         ENDIF
      ENDDO

      ! Setup global to local and vice versa mappings.

      ilsrcmppmap(:)=-1
      ilsrclocind(:)=0
      ildstmppmap(:)=-1
      ildstlocind(:)=0      
      
      DO i=1,nsrclocpoints
         ilsrcmppmap(srcgloind(i))=mype
         ilsrclocind(srcgloind(i))=i
      ENDDO

      DO i=1,ndstlocpoints
         ildstmppmap(dstgloind(i))=mype
         ildstlocind(dstgloind(i))=i
      ENDDO
      
#if defined key_mpp_mpi
      CALL mpi_allreduce(ilsrcmppmap,igsrcmppmap,nsrcglopoints, &
         &               mpi_integer,mpi_max,mpi_comm,istatus)
      CALL mpi_allreduce(ilsrclocind,igsrclocind,nsrcglopoints, &
         &               mpi_integer,mpi_max,mpi_comm,istatus)
      CALL mpi_allreduce(ildstmppmap,igdstmppmap,ndstglopoints, &
         &               mpi_integer,mpi_max,mpi_comm,istatus)
      CALL mpi_allreduce(ildstlocind,igdstlocind,ndstglopoints, &
         &               mpi_integer,mpi_max,mpi_comm,istatus)
#else
      igsrcmppmap(:)=ilsrcmppmap(:)
      igsrclocind(:)=ilsrclocind(:)
      igdstmppmap(:)=ildstmppmap(:)
      igdstlocind(:)=ildstlocind(:)
#endif

      ! Optionally construct an ascii file listing what src and
      ! dest points belongs to which task

      ! Since igsrcmppmap and igdstmppmap are global data only do
      ! this for mype==0.

      IF (lcommout.AND.(mype==0)) THEN
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_srcmppmap_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO i=1,nsrcglopoints
            WRITE(9,*)remap%src%grid_center_lat(i),&
               & remap%src%grid_center_lon(i), &
               & igsrcmppmap(i)+1,remap%src%grid_imask(i)
         ENDDO
         CLOSE(9)
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_dstmppmap_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO i=1,ndstglopoints
            WRITE(9,*)remap%dst%grid_center_lat(i),&
               & remap%dst%grid_center_lon(i), &
               & igdstmppmap(i)+1,remap%dst%grid_imask(i)
         ENDDO
         CLOSE(9)
      ENDIF

      !
      ! Standard interpolation in serial case is
      !
      ! DO n=1,remap%num_links
      !   zdst(remap%dst_address(n)) = zdst(remap%dst_address(n)) + &
      !      & remap%remap_matrix(1,n)*zsrc(remap%src_address(n))
      ! END DO
      !

      ! In parallel we need to first find local number of links
      
      pinfo%num_links=0
      DO i=1,remap%num_links
         IF (igdstmppmap(remap%dst_address(i))==mype) &
            & pinfo%num_links=pinfo%num_links+1
      ENDDO
      ALLOCATE(pinfo%dst_address(pinfo%num_links),&
         &     pinfo%src_address(pinfo%num_links),&
         &     pinfo%remap_matrix(1,pinfo%num_links))

      ! Get local destination addresses

      n=0
      DO i=1,remap%num_links
         IF (igdstmppmap(remap%dst_address(i))==mype) THEN
            n=n+1
            pinfo%dst_address(n)=&
               & igdstlocind(remap%dst_address(i))
            pinfo%remap_matrix(:,n)=&
               & remap%remap_matrix(:,i)
         ENDIF
      ENDDO

      ! Get sending processors maps.

      ! The same data point might need to be sent to many processors
      ! so first construct a map for processors needing the data

      isrcpe(:)=-1
      DO i=1,remap%num_links
         IF (igdstmppmap(remap%dst_address(i))==mype) THEN
            isrcpe(remap%src_address(i))=&
               & igsrcmppmap(remap%src_address(i))
         ENDIF
      ENDDO

      ! Optionally write a set if ascii file listing which tasks
      ! mype needs to send to communicate with

      IF (lcommout) THEN
         ! Destination processors
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_dsts_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO pe=0,nproc-1
            IF (pe==mype) THEN
               isrcpetmp(:)=isrcpe(:)
            ENDIF
#if defined key_mpp_mpi
            CALL mpi_bcast(isrcpetmp,nsrcglopoints,mpi_integer,pe,mpi_comm,istatus)
#endif
            DO i=1,nsrcglopoints
               IF (isrcpetmp(i)==mype) THEN
                  WRITE(9,*)remap%src%grid_center_lat(i),&
                     & remap%src%grid_center_lon(i), &
                     & pe+1,mype+1
               ENDIF
            ENDDO
         ENDDO
         CLOSE(9)
      ENDIF

      ! Get number of points to send to each processor

      ALLOCATE(pinfo%nsend(0:nproc-1))
      isrcno(:)=0
      DO i=1,nsrcglopoints
         IF (isrcpe(i)>=0) THEN
            isrcno(isrcpe(i)+1)=isrcno(isrcpe(i)+1)+1
         ENDIF
      ENDDO
#if defined key_mpp_mpi
      CALL mpi_alltoall(isrcno,1,mpi_integer, &
         &              pinfo%nsend(0:nproc-1),1,mpi_integer, &
         &              mpi_comm,istatus)
#else
      pinfo%nsend(0:nproc-1) = isrcno(1:nproc)
#endif
      pinfo%nsendtot=SUM(pinfo%nsend(0:nproc-1))

      ! Construct sending buffer mapping. Data is mapping in
      ! processor order.

      ALLOCATE(pinfo%send_address(pinfo%nsendtot))
      
      ! Temporary arrays for mpi all to all.

      ALLOCATE(isrctot(SUM(isrcno(1:nproc))))
      ALLOCATE(ircvtot(SUM(pinfo%nsend(0:nproc-1))))

      ! Offset for message parsing

      isrcoff(1)=0
      ircvoff(1)=0
      DO i=1,nproc-1
         isrcoff(i+1) = isrcoff(i) + isrcno(i)
         ircvoff(i+1) = pinfo%nsend(i-1) + ircvoff(i)
      ENDDO

      ! Pack indices i into a buffer

      isrccur(:)=0
      DO i=1,nsrcglopoints
         IF (isrcpe(i)>=0) THEN
            isrccur(isrcpe(i)+1)=isrccur(isrcpe(i)+1)+1
            isrctot(isrccur(isrcpe(i)+1)+isrcoff(isrcpe(i)+1)) = i
         ENDIF
      ENDDO
      
      ! Send the data

#if defined key_mpp_mpi
      CALL mpi_alltoallv(&
         & isrctot,isrccur,isrcoff,mpi_integer, &
         & ircvtot,pinfo%nsend(0:nproc-1),ircvoff,mpi_integer, &
         & mpi_comm,istatus)
#else
      ircvtot(:)=isrctot(:)
#endif

      ! Get the send address. ircvtot will at this point contain the 
      ! addresses in the global index needed for message passing

      DO i=1,pinfo%nsendtot
         pinfo%send_address(i)=igsrclocind(ircvtot(i))
      ENDDO

      ! Deallocate the mpi all to all arrays

      DEALLOCATE(ircvtot,isrctot)

      ! Get number of points to receive to each processor

      ALLOCATE(pinfo%nrecv(0:nproc-1))
      pinfo%nrecv(0:nproc-1)=0
      DO i=1,nsrcglopoints
        IF (isrcpe(i)>=0 .AND. isrcpe(i)<nproc) THEN
          pinfo%nrecv(isrcpe(i))=pinfo%nrecv(isrcpe(i))+1
        ENDIF
      ENDDO
      pinfo%nrecvtot=SUM(pinfo%nrecv(0:nproc-1))

      ! Find new src address mapping

      ! Setup local positions in the global array in a temporary array
      ! taking into accound the processor ordering

      isrcaddtmp(:)=-1
      isrcoffset(0)=0
      DO pe=1,nproc-1
        isrcoffset(pe)=isrcoffset(pe-1)+pinfo%nrecv(pe-1)
      ENDDO

      DO i=1,nsrcglopoints
        IF (isrcpe(i)>=0 .AND. isrcpe(i)<nproc) THEN
          isrcoffset(isrcpe(i))=isrcoffset(isrcpe(i))+1
          isrcaddtmp(i)=isrcoffset(isrcpe(i))
        ENDIF
      ENDDO 

      ! Find the local positions in the temporary array.

      n=0
      DO i=1,remap%num_links
         IF (igdstmppmap(remap%dst_address(i))==mype) THEN
            n=n+1
            pinfo%src_address(n)=isrcaddtmp(remap%src_address(i))
         ENDIF
      ENDDO

      ! MPI displacements for mpi_alltoallv

      ALLOCATE(pinfo%nsdisp(0:nproc-1),pinfo%nrdisp(0:nproc-1))
      pinfo%nsdisp(0)=0
      pinfo%nrdisp(0)=0
      DO pe=1,nproc-1
         pinfo%nsdisp(pe)=pinfo%nsdisp(pe-1)+pinfo%nsend(pe-1)
         pinfo%nrdisp(pe)=pinfo%nrdisp(pe-1)+pinfo%nrecv(pe-1)
      ENDDO

      ! Optionally construct an ascii file listing number of
      ! data to send and receive from which processor.

      IF (lcommout) THEN
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_nsend_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO pe=0,nproc-1
            WRITE(9,*)pe+1,pinfo%nsend(pe)
         ENDDO
         CLOSE(9)
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_nrecv_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO pe=0,nproc-1
            WRITE(9,*)pe+1,pinfo%nrecv(pe)
         ENDDO
         CLOSE(9)
      ENDIF

   END SUBROUTINE parinter_init

   SUBROUTINE parinter_fld( mype, nproc, mpi_comm, &
      & pinfo, nsrclocpoints, zsrc, ndstlocpoints, zdst )

      ! Perform a single interpolation from the zsrc field
      ! to zdst field based on the information in pinfk
      
      ! Input arguments

      ! Message passing information
      INTEGER, INTENT(IN) :: mype, nproc, mpi_comm
      ! Interpolation setup
      TYPE(parinterinfo), INTENT(IN) :: pinfo
      ! Source data/
      INTEGER, INTENT(IN) :: nsrclocpoints
      REAL, INTENT(IN), DIMENSION(nsrclocpoints) :: zsrc

      ! Output arguments

      ! Destination data
      INTEGER , INTENT(IN):: ndstlocpoints
      REAL, DIMENSION(ndstlocpoints) :: zdst

      ! Local variables

      ! MPI send/recv buffers
      REAL(scripdp) :: zsend(pinfo%nsendtot),zrecv(pinfo%nrecvtot)
      ! Misc variables
      INTEGER :: i,istatus

      ! Pack the sending buffer

      DO i=1,pinfo%nsendtot
         zsend(i)=zsrc(pinfo%send_address(i))
      ENDDO
      
      ! Do the message passing

#if defined key_mpp_mpi

      CALL mpi_alltoallv(&
         & zsend,pinfo%nsend(0:nproc-1),&
         & pinfo%nsdisp(0:nproc-1),mpi_double_precision, &
         & zrecv,pinfo%nrecv(0:nproc-1), &
         & pinfo%nrdisp(0:nproc-1),mpi_double_precision, &
         & mpi_comm,istatus)

#else

      zrecv(:)=zsend(:)

#endif

      ! Do the interpolation

      zdst(:)=0.0
      DO i=1,pinfo%num_links
         zdst(pinfo%dst_address(i)) = zdst(pinfo%dst_address(i)) + &
            & pinfo%remap_matrix(1,i)*zrecv(pinfo%src_address(i))
      END DO

   END SUBROUTINE parinter_fld

   SUBROUTINE parinter_write( mype, nproc, &
      & nsrcglopoints, ndstglopoints, &
      & pinfo, cdpath, cdprefix )

      ! Write pinfo information in a netCDF file in order to
      ! be able to read it rather than calling parinter_init

      ! Input arguments.

      ! Message passing information
      INTEGER, INTENT(IN) :: mype, nproc
      ! Source grid local global number of grid points
      INTEGER, INTENT(IN) :: nsrcglopoints
      ! Destination grid global number of grid points
      INTEGER, INTENT(IN) :: ndstglopoints
      ! Interpolation and message passing information
      TYPE(parinterinfo), INTENT(IN) :: pinfo
      ! Path and file prefix
      CHARACTER(len=*) :: cdpath, cdprefix

      ! Local variable

      ! Misc variable
      CHARACTER(len=1024) :: cdfile
      INTEGER :: ncid, dimnl, dimnw, dimnst, dimnrt, dimnpr
      INTEGER :: dims1(1), dims2(2)
      INTEGER :: idda, idsa, idrm, idns, idsaa, idnr, idnrp, idnsp

      WRITE(cdfile,'(A,2(I8.8,A),2(I4.4,A),A)') &
         & TRIM(cdpath)//'/'//TRIM(cdprefix)//'_', &
         & nsrcglopoints,'_',ndstglopoints,'_',mype,'_',nproc,'.nc'

      CALL nchdlerr(nf90_create(TRIM(cdfile),nf90_clobber,ncid), &
         &          __LINE__, __MYFILE__ )

      ! To avoid problems with multiple unlimited netCDF dimensions
      ! we don't write num_links, nsendtot and nrecvtot if
      ! they are 0. If they are not there we assume they are 0
      ! in parinter_read.

      IF (pinfo%num_links>0) THEN
         CALL nchdlerr(nf90_def_dim(ncid,'num_links',&
            &                       pinfo%num_links,dimnl),&
            &          __LINE__,__MYFILE__)
      ENDIF

      CALL nchdlerr(nf90_def_dim(ncid,'num_wgts',&
         &                       1,dimnw),&
         &          __LINE__,__MYFILE__)

      IF (pinfo%nsendtot>0) THEN
         CALL nchdlerr(nf90_def_dim(ncid,'nsendtot',&
            &                       pinfo%nsendtot,dimnst),&
            &          __LINE__,__MYFILE__)
      ENDIF
      
      IF (pinfo%nrecvtot>0) THEN
         CALL nchdlerr(nf90_def_dim(ncid,'nrecvtot',&
            &                       pinfo%nrecvtot,dimnrt),&
            &          __LINE__,__MYFILE__)
      ENDIF
      
      CALL nchdlerr(nf90_def_dim(ncid,'nproc',&
         &                       nproc,dimnpr),&
         &          __LINE__,__MYFILE__)

      IF (pinfo%num_links>0) THEN

         dims1(1)=dimnl
         CALL nchdlerr(nf90_def_var(ncid,'dst_address',&
            &                       nf90_int,dims1,idda),&
            &          __LINE__,__MYFILE__)
         
         dims1(1)=dimnl
         CALL nchdlerr(nf90_def_var(ncid,'src_address',&
            &                       nf90_int,dims1,idsa),&
            &          __LINE__,__MYFILE__)
         
         dims2(1)=dimnw
         dims2(2)=dimnl
         CALL nchdlerr(nf90_def_var(ncid,'remap_matrix',&
            &                       nf90_double,dims2,idrm),&
            &          __LINE__,__MYFILE__)

      ENDIF

      dims1(1)=dimnpr
      CALL nchdlerr(nf90_def_var(ncid,'nsend',&
         &                       nf90_int,dims1,idns),&
         &          __LINE__,__MYFILE__)

      IF (pinfo%nsendtot>0) THEN

         dims1(1)=dimnst
         CALL nchdlerr(nf90_def_var(ncid,'send_address',&
            &                       nf90_int,dims1,idsaa),&
            &          __LINE__,__MYFILE__)

      ENDIF

      dims1(1)=dimnpr
      CALL nchdlerr(nf90_def_var(ncid,'nrecv',&
         &                       nf90_int,dims1,idnr),&
         &          __LINE__,__MYFILE__)

      dims1(1)=dimnpr
      CALL nchdlerr(nf90_def_var(ncid,'nsdisp',&
         &                       nf90_int,dims1,idnsp),&
         &          __LINE__,__MYFILE__)

      dims1(1)=dimnpr
      CALL nchdlerr(nf90_def_var(ncid,'nrdisp',&
         &                       nf90_int,dims1,idnrp),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_enddef(ncid),__LINE__,__MYFILE__)


      IF (pinfo%num_links>0) THEN

         CALL nchdlerr(nf90_put_var(ncid,idda,pinfo%dst_address),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_put_var(ncid,idsa,pinfo%src_address),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_put_var(ncid,idrm,pinfo%remap_matrix),&
            &          __LINE__,__MYFILE__)

      ENDIF

      CALL nchdlerr(nf90_put_var(ncid,idns,pinfo%nsend(0:nproc-1)),&
         &          __LINE__,__MYFILE__)

      IF (pinfo%nsendtot>0) THEN

         CALL nchdlerr(nf90_put_var(ncid,idsaa,pinfo%send_address),&
            &          __LINE__,__MYFILE__)

      ENDIF

      CALL nchdlerr(nf90_put_var(ncid,idnr,pinfo%nrecv(0:nproc-1)),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,idnsp,pinfo%nsdisp(0:nproc-1)),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,idnrp,pinfo%nrdisp(0:nproc-1)),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_close(ncid),__LINE__, __MYFILE__ )

   END SUBROUTINE parinter_write

   SUBROUTINE parinter_read( mype, nproc, &
      & nsrcglopoints, ndstglopoints, &
      & pinfo, cdpath, cdprefix, lexists )

      ! Write pinfo information in a netCDF file in order to
      ! be able to read it rather than calling parinter_init

      ! Input arguments.

      ! Message passing information
      INTEGER, INTENT(IN) :: mype, nproc
      ! Source grid local global number of grid points
      INTEGER, INTENT(IN) :: nsrcglopoints
      ! Destination grid global number of grid points
      INTEGER, INTENT(IN) :: ndstglopoints
      ! Interpolation and message passing information
      TYPE(parinterinfo), INTENT(OUT) :: pinfo
      ! Does the information exists
      LOGICAL :: lexists
      ! Path and file prefix
      CHARACTER(len=*) :: cdpath, cdprefix

      ! Local variable

      ! Misc variable
      CHARACTER(len=1024) :: cdfile
      INTEGER :: ncid, dimid, varid, num_wgts

      WRITE(cdfile,'(A,2(I8.8,A),2(I4.4,A),A)') &
         & TRIM(cdpath)//'/'//TRIM(cdprefix)//'_', &
         & nsrcglopoints,'_',ndstglopoints,'_',mype,'_',nproc,'.nc'


      lexists=nf90_open(TRIM(cdfile),nf90_nowrite,ncid)==nf90_noerr

      IF (lexists) THEN

         ! If num_links is not present we assume it to be zero.
         
         IF (nf90_inq_dimid(ncid,'num_links',dimid)==nf90_noerr) THEN
            CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
               &                                 len=pinfo%num_links),&
               &          __LINE__,__MYFILE__)
         ELSE
            pinfo%num_links=0
         ENDIF

         CALL nchdlerr(nf90_inq_dimid(ncid,'num_wgts',dimid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
            &                                 len=num_wgts),&
            &          __LINE__,__MYFILE__)
         IF (num_wgts/=1) THEN
            WRITE(0,*)'parinter_read: num_wgts has to be 1 for now'
            CALL abort
         ENDIF

         ! If nsendtot is not present we assume it to be zero.

         IF (nf90_inq_dimid(ncid,'nsendtot',dimid)==nf90_noerr) THEN
            CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
               &                                 len=pinfo%nsendtot),&
               &          __LINE__,__MYFILE__)
         ELSE
            pinfo%nsendtot=0
         ENDIF

         IF(nf90_inq_dimid(ncid,'nrecvtot',dimid)==nf90_noerr) THEN
            CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
               &                                 len=pinfo%nrecvtot),&
               &          __LINE__,__MYFILE__)
         ELSE
            pinfo%nrecvtot=0
         ENDIF

         ALLOCATE(pinfo%dst_address(pinfo%num_links),&
            &     pinfo%src_address(pinfo%num_links),&
            &     pinfo%remap_matrix(num_wgts,pinfo%num_links),&
            &     pinfo%nsend(0:nproc-1),&
            &     pinfo%send_address(pinfo%nsendtot),&
            &     pinfo%nrecv(0:nproc-1),&
            &     pinfo%nsdisp(0:nproc-1),&
            &     pinfo%nrdisp(0:nproc-1))
 
         IF (pinfo%num_links>0) THEN
            CALL nchdlerr(nf90_inq_varid(ncid,'dst_address',varid),&
               &          __LINE__,__MYFILE__)
            CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%dst_address),&
               &          __LINE__,__MYFILE__)
            
            CALL nchdlerr(nf90_inq_varid(ncid,'src_address',varid),&
               &          __LINE__,__MYFILE__)
            CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%src_address),&
               &          __LINE__,__MYFILE__)
            
            CALL nchdlerr(nf90_inq_varid(ncid,'remap_matrix',varid),&
               &          __LINE__,__MYFILE__)
            CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%remap_matrix),&
               &          __LINE__,__MYFILE__)
         ENDIF

         CALL nchdlerr(nf90_inq_varid(ncid,'nsend',varid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%nsend(0:nproc-1)),&
            &          __LINE__,__MYFILE__)

         IF (pinfo%nsendtot>0) THEN

            CALL nchdlerr(nf90_inq_varid(ncid,'send_address',varid),&
               &          __LINE__,__MYFILE__)
            CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%send_address),&
               &          __LINE__,__MYFILE__)

         ENDIF

         CALL nchdlerr(nf90_inq_varid(ncid,'nrecv',varid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%nrecv(0:nproc-1)),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_inq_varid(ncid,'nsdisp',varid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%nsdisp(0:nproc-1)),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_inq_varid(ncid,'nrdisp',varid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%nrdisp(0:nproc-1)),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_close(ncid),__LINE__, __MYFILE__ )

      ENDIF

   END SUBROUTINE parinter_read
   
END MODULE parinter
