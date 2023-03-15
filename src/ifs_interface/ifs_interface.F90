!=====================================================
! IFS interface for calling FESOM2 as a subroutine.
!
! -Original code for NEMO by Kristian Mogensen, ECMWF.
! -Adapted to FESOM2 by Thomas Rackow, AWI, 2018.
!-----------------------------------------------------

MODULE nemogcmcoup_steps
   INTEGER :: substeps !per IFS timestep
END MODULE nemogcmcoup_steps

SUBROUTINE nemogcmcoup_init( mype, icomm, inidate, initime, itini, itend, zstp, &
   & lwaveonly, iatmunit, lwrite )

   ! Initialize the FESOM model for single executable coupling 

   USE par_kind !in ifs_modules.F90
   USE fesom_main_storage_module, only: fesom => f ! only: MPI_COMM_FESOM, mype (previously in g_parsup)
   USE fesom_module, ONLY : fesom_init
   USE g_config, only: dt
   USE g_clock, only: timenew, daynew, yearnew, month, day_in_month
   USE nemogcmcoup_steps, ONLY : substeps
   
   IMPLICIT NONE

   ! Input arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype ! was added to ifs/nemo/ininemo.F90 to allow diagnostics based on the first tasks only
   INTEGER, INTENT(IN) :: icomm
   ! Initial date (e.g. 20170906), time, initial timestep and final time step
   INTEGER, INTENT(OUT) ::  inidate, initime, itini, itend
   ! Length of the time step
   REAL(wpIFS), INTENT(OUT) :: zstp

   ! inherited from interface to NEMO, not used here:
   ! Coupling to waves only
   LOGICAL, INTENT(IN) :: lwaveonly
   ! Logfile unit (used if >=0)
   INTEGER :: iatmunit
   ! Write to this unit
   LOGICAL :: lwrite
   ! FESOM might perform substeps
   INTEGER :: itend_fesom
   INTEGER :: i
   NAMELIST/namfesomstep/substeps

   ! overwritten from value namelist
   substeps=2
   OPEN(9,file='namfesomstep.in')
   READ(9,namfesomstep)
   CLOSE(9)

   fesom%partit%MPI_COMM_FESOM=icomm

   itini = 1
   CALL fesom_init(itend_fesom) !also sets mype and npes 
   itend=itend_fesom/substeps
   if(fesom%mype==0) then
   WRITE(0,*)'!======================================'
   WRITE(0,*)'! FESOM is initialized from within IFS.'
   WRITE(0,*)'! get MPI_COMM_FESOM. ================='
   WRITE(0,*)'! main_initialize done. ==============='
   WRITE(0,*)'Thomas/Kristian parinter_mult version'
   endif

   ! Set more information for the caller
   
   ! initial date and time (time is not used)
   inidate = yearnew*10000 + month*100 + day_in_month ! e.g. 20170906
   initime = 0
   if(fesom%mype==0) then
   WRITE(0,*)'! FESOM initial date is ', inidate ,' ======'
   WRITE(0,*)'! FESOM substeps are ', substeps ,' ======'
   endif

   ! fesom timestep (as seen by IFS)
   zstp = REAL(substeps,wpIFS)*dt
   if(fesom%mype==0) then
   WRITE(0,*)'! FESOM timestep as seen by IFS is ', real(zstp,4), 'sec (',substeps,'xdt)'
   WRITE(0,*)'!======================================'
   endif

END SUBROUTINE nemogcmcoup_init


SUBROUTINE nemogcmcoup_coupinit( mypeIN, npesIN, icomm, &
   &                             npoints, nlocmsk, ngloind )

   ! FESOM modules
   USE fesom_main_storage_module, only: fesom => f ! only: mype, npes, myDim_nod2D, eDim_nod2D, myDim_elem2D, eDim_elem2D, eXDim_elem2D, &
							!  myDim_edge2D, eDim_edge2D, myList_nod2D, myList_elem2D

   ! Initialize single executable coupling 
   USE parinter
   USE scripremap
   USE interinfo
   IMPLICIT NONE

   ! Input arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mypeIN,npesIN,icomm
   ! Gaussian grid information   
   ! Number of points
   INTEGER, INTENT(IN) :: npoints
   ! Integer mask and global indices
   INTEGER, DIMENSION(npoints), INTENT(IN) :: nlocmsk, ngloind
   INTEGER :: iunit = 0

   ! Local variables
   !type(t_mesh), target :: mesh
   integer    , pointer :: nod2D 
   integer    , pointer :: elem2D   
   integer,     pointer :: myDim_nod2D, eDim_nod2D
   integer, dimension(:), pointer :: myList_nod2D
   integer, pointer               :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
   integer, dimension(:), pointer :: myList_elem2D

   ! Namelist containing the file names of the weights
   CHARACTER(len=256) :: cdfile_gauss_to_T, cdfile_gauss_to_UV, &
      &                  cdfile_T_to_gauss, cdfile_UV_to_gauss
   CHARACTER(len=256) :: cdpathdist
   LOGICAL :: lwritedist, lreaddist
   LOGICAL :: lcommout
   CHARACTER(len=128) :: commoutprefix
   NAMELIST/namfesomcoup/cdfile_gauss_to_T,&
      &                 cdfile_gauss_to_UV,&
      &                 cdfile_T_to_gauss,&
      &                 cdfile_UV_to_gauss,&
      &                 cdpathdist, &
      &                 lreaddist, &
      &                 lwritedist, &
      &                 lcommout, &
      &                 commoutprefix,&
      &                 lparbcast, &
      &                 lparinterp2p, &
      &                 lparintmultatm

   ! Global number of gaussian gridpoints
   INTEGER :: nglopoints
   ! Ocean grids accessed with NEMO modules
   INTEGER :: noglopoints,nopoints
   INTEGER, ALLOCATABLE, DIMENSION(:) :: omask,ogloind
   ! SCRIP remapping data structures.
   TYPE(scripremaptype) :: remap_gauss_to_T, remap_T_to_gauss, &
      & remap_gauss_to_UV, remap_UV_to_gauss
   ! Misc variables
   INTEGER :: i,j,k,ierr
   LOGICAL :: lexists

   ! associate the mesh, only what is needed here
   ! #include "associate_mesh.h"
   nod2D              => fesom%mesh%nod2D              
   elem2D             => fesom%mesh%elem2D  
   myDim_nod2D        => fesom%partit%myDim_nod2D
   eDim_nod2D         => fesom%partit%eDim_nod2D
   myList_nod2D(1:myDim_nod2D+eDim_nod2D) => fesom%partit%myList_nod2D
   myDim_elem2D       => fesom%partit%myDim_elem2D
   eDim_elem2D        => fesom%partit%eDim_elem2D
   eXDim_elem2D       => fesom%partit%eXDim_elem2D
   myList_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => fesom%partit%myList_elem2D

   ! here FESOM knows about the (total number of) MPI tasks

   if(fesom%mype==0) then
   write(*,*) 'MPI has been initialized in the atmospheric model'
   write(*, *) 'Running on ', fesom%npes, ' PEs'
   end if

   ! Read namelists
   
   cdfile_gauss_to_T = 'gausstoT.nc'
   cdfile_gauss_to_UV = 'gausstoUV.nc'
   cdfile_T_to_gauss = 'Ttogauss.nc'
   cdfile_UV_to_gauss = 'UVtogauss.nc'
   lcommout          = .FALSE.
   commoutprefix     = 'parinter_comm'
   cdpathdist        = './'
   lreaddist         = .FALSE.
   lwritedist        = .FALSE.
   lparintmultatm    = .TRUE.

   OPEN(9,file='namfesomcoup.in')
   READ(9,namfesomcoup)
   CLOSE(9)

   ! Global number of Gaussian gridpoints

   CALL mpi_allreduce( npoints, nglopoints, 1, &
      &                mpi_integer, mpi_sum, icomm, ierr)


   if(fesom%mype==0) then
   WRITE(0,*)'!======================================'
   WRITE(0,*)'! SCALARS ============================='

   WRITE(0,*)'Update FESOM global scalar points'
   endif

   noglopoints=nod2D
   nopoints=myDim_nod2d

   ! Ocean mask and global indicies
   
   ALLOCATE(omask(MAX(nopoints,1)),ogloind(MAX(nopoints,1)))
   omask(:)= 1			! all points are ocean points
   ogloind(1:myDim_nod2d)= myList_nod2D(1:myDim_nod2d)	! global index for local point number

   ! Could be helpful later:
   ! Replace global numbering with a local one
   ! tmp(1:nod2d)=0
   ! DO n=1, myDim_nod2D+eDim_nod2D
   ! tmp(myList_nod2D(n))=n

   ! Read the interpolation weights and setup the parallel interpolation
   ! from atmosphere Gaussian grid to ocean T-grid
   
   IF (lreaddist) THEN
      CALL parinter_read( fesom%mype, fesom%npes, nglopoints, noglopoints, gausstoT,  &
         & cdpathdist,'ifs_to_fesom_gridT',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_gauss_to_T,remap_gauss_to_T,&
            & fesom%mype,fesom%npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_gauss_to_T,remap_gauss_to_T)
      ENDIF
      CALL parinter_init( fesom%mype, fesom%npes, icomm, &
         & npoints, nglopoints, nlocmsk, ngloind, &
         & nopoints, noglopoints, omask, ogloind, & 
         & remap_gauss_to_T, gausstoT, lcommout, TRIM(commoutprefix)//'_gtoT', &
         & iunit )
      CALL scripremap_dealloc(remap_gauss_to_T)
      IF (lwritedist) THEN
         CALL parinter_write( fesom%mype, fesom%npes, nglopoints, noglopoints, gausstoT,  &
            & cdpathdist,'ifs_to_fesom_gridT')
      ENDIF
   ENDIF

   ! From ocean T-grid to atmosphere Gaussian grid

   IF (lreaddist) THEN
      CALL parinter_read( fesom%mype, fesom%npes, noglopoints, nglopoints, Ttogauss,  &
         & cdpathdist,'fesom_gridT_to_ifs',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_T_to_gauss,remap_T_to_gauss,&
            & fesom%mype,fesom%npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_T_to_gauss,remap_T_to_gauss)
      ENDIF

      CALL parinter_init( fesom%mype, fesom%npes, icomm, &
         & nopoints, noglopoints, omask, ogloind, & 
         & npoints, nglopoints, nlocmsk, ngloind, &
         & remap_T_to_gauss, Ttogauss, lcommout, TRIM(commoutprefix)//'_Ttog', &
         & iunit )
      CALL scripremap_dealloc(remap_T_to_gauss)
      IF (lwritedist) THEN
         CALL parinter_write( fesom%mype, fesom%npes, noglopoints, nglopoints, Ttogauss,  &
            & cdpathdist,'fesom_gridT_to_ifs')
      ENDIF
   ENDIF
   
   DEALLOCATE(omask,ogloind)


   if(fesom%mype==0) then
   WRITE(0,*)'!======================================'
   WRITE(0,*)'! VECTORS ============================='

   WRITE(0,*)'Update FESOM global vector points'
   endif
   noglopoints=elem2D
   nopoints=myDim_elem2D

   ! Ocean mask and global indicies
   
   ALLOCATE(omask(MAX(nopoints,1)),ogloind(MAX(nopoints,1)))

   omask(:)= 1			! all elements are in the ocean 	
   ogloind(1:myDim_elem2D) = myList_elem2D(1:myDim_elem2D) ! global index for local element number

   ! Read the interpolation weights and setup the parallel interpolation
   ! from atmosphere Gaussian grid to ocean UV-grid
   
   IF (lreaddist) THEN
      CALL parinter_read( fesom%mype, fesom%npes, nglopoints, noglopoints, gausstoUV,  &
         & cdpathdist,'ifs_to_fesom_gridUV',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_gauss_to_UV,remap_gauss_to_UV,&
            & fesom%mype,fesom%npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_gauss_to_UV,remap_gauss_to_UV)
      ENDIF
      CALL parinter_init( fesom%mype, fesom%npes, icomm, &
         & npoints, nglopoints, nlocmsk, ngloind, &
         & nopoints, noglopoints, omask, ogloind, & 
         & remap_gauss_to_UV, gausstoUV, lcommout, TRIM(commoutprefix)//'_gtoUV', &
         & iunit )
      CALL scripremap_dealloc(remap_gauss_to_UV)
      IF (lwritedist) THEN
         CALL parinter_write( fesom%mype, fesom%npes, nglopoints, noglopoints, gausstoUV,  &
            & cdpathdist,'ifs_to_fesom_gridUV')
      ENDIF
   ENDIF

   ! From ocean UV-grid to atmosphere Gaussian grid

   IF (lreaddist) THEN
      CALL parinter_read( fesom%mype, fesom%npes, noglopoints, nglopoints, UVtogauss,  &
         & cdpathdist,'fesom_gridUV_to_ifs',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_UV_to_gauss,remap_UV_to_gauss,&
            & fesom%mype,fesom%npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_UV_to_gauss,remap_UV_to_gauss)
      ENDIF

      CALL parinter_init( fesom%mype, fesom%npes, icomm, &
         & nopoints, noglopoints, omask, ogloind, & 
         & npoints, nglopoints, nlocmsk, ngloind, &
         & remap_UV_to_gauss, UVtogauss, lcommout, TRIM(commoutprefix)//'_UVtog', &
         & iunit )
      CALL scripremap_dealloc(remap_UV_to_gauss)
      IF (lwritedist) THEN
         CALL parinter_write( fesom%mype, fesom%npes, noglopoints, nglopoints, UVtogauss,  &
            & cdpathdist,'fesom_gridUV_to_ifs')
      ENDIF
   ENDIF
   
   DEALLOCATE(omask,ogloind)
         
END SUBROUTINE nemogcmcoup_coupinit


SUBROUTINE nemogcmcoup_lim2_get( mype, npes, icomm, &
   &                             nopoints, pgsst, pgist, pgalb, &
   &                             pgifr, pghic, pghsn, pgucur, pgvcur, &
   &                             pgistl, licelvls )

   ! Interpolate sst, ice; surf T; albedo; concentration; thickness,
   ! snow thickness and currents from the FESOM grid to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind ! in ifs_modules.F90
   USE fesom_main_storage_module, only: fesom => f
   !USE o_ARRAYS, ONLY : UV ! tr_arr is now tracers, UV in dynamics derived type
   !USE i_arrays, ONLY : m_ice, a_ice, m_snow
   !USE i_therm_param, ONLY : tmelt
   USE g_rotate_grid, only: vector_r2g
   USE parinter
   USE scripremap
   USE interinfo

   IMPLICIT NONE
   
   ! Arguments
   REAL(wpIFS), DIMENSION(nopoints) :: pgsst, pgist, pgalb, pgifr, pghic, pghsn, pgucur, pgvcur
   REAL(wpIFS), DIMENSION(nopoints,3) :: pgistl
   LOGICAL :: licelvls

   !type(t_mesh), target :: mesh
   real(kind=wpIFS), dimension(:,:), pointer :: coord_nod2D
   integer, dimension(:,:)         , pointer :: elem2D_nodes
   integer, pointer :: myDim_nod2D, eDim_nod2D
   integer, pointer :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
   real(kind=wpIFS), dimension(:), pointer :: a_ice, m_ice, m_snow, ice_temp, ice_alb
   real(kind=wpIFS)              , pointer :: tmelt
   
   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number Gaussian grid points
   INTEGER, INTENT(IN) :: nopoints

   ! Local variables
   INTEGER , PARAMETER :: maxnfield = 8
   INTEGER , PARAMETER :: maxnfielduv = 2
   INTEGER :: nfield = 0
   INTEGER :: nfielduv = 0
   REAL(wpIFS), DIMENSION(fesom%partit%myDim_nod2D,maxnfield)  :: zsendnf
   REAL(wpIFS), DIMENSION(fesom%partit%myDim_elem2D,maxnfielduv) :: zsendnfUV
   REAL(wpIFS), DIMENSION(nopoints,maxnfield)  :: zrecvnf
   REAL(wpIFS), DIMENSION(nopoints,maxnfielduv) :: zrecvnfUV
   INTEGER			     :: elnodes(3)
   REAL(wpIFS)			     :: rlon, rlat	

   ! Loop variables
   INTEGER :: n, elem, ierr, jf

   !#include "associate_mesh.h"
   ! associate what is needed only
   myDim_nod2D  => fesom%partit%myDim_nod2D
   eDim_nod2D   => fesom%partit%eDim_nod2D

   myDim_elem2D => fesom%partit%myDim_elem2D
   eDim_elem2D  => fesom%partit%eDim_elem2D
   eXDim_elem2D => fesom%partit%eXDim_elem2D

   coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)                  => fesom%mesh%coord_nod2D(:,:)   
   elem2D_nodes(1:3, 1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => fesom%mesh%elem2D_nodes(:,:)
   a_ice        => fesom%ice%data(1)%values(:)
   m_ice        => fesom%ice%data(2)%values(:)
   m_snow       => fesom%ice%data(3)%values(:)
   ice_temp     => fesom%ice%data(4)%values(:)
   ice_alb      => fesom%ice%atmcoupl%ice_alb(:)
   tmelt        => fesom%ice%thermo%tmelt ! scalar const.


   nfield = 0
   ! =================================================================== !
   ! Pack SST data and convert to K. 'pgsst' is on Gauss grid.
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=fesom%tracers%DATA(1)%values(1, n) +tmelt ! sea surface temperature [K], 
							 ! (1=surface, n=node, data(1/2)=T/S)
   ENDDO

   ! =================================================================== !
   ! Pack ice fraction data [0..1]
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=a_ice(n)
   ENDDO
   
   ! =================================================================== !
   ! Pack ice temperature data (already in K)
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=ice_temp(n)
   ENDDO

   ! =================================================================== !
   ! Pack ice albedo data and interpolate: 'pgalb' on Gaussian grid.
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=ice_alb(n)
   ENDDO

   ! =================================================================== !
   ! Pack ice thickness data and interpolate: 'pghic' on Gaussian grid.
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=m_ice(n)!/MAX(a_ice(n),0.01) ! ice thickness (mean over ice)
   ENDDO
   
   ! =================================================================== !
   ! Pack snow thickness data and interpolate: 'pghsn' on Gaussian grid.
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=m_snow(n)!/MAX(a_ice(n),0.01) ! snow thickness (mean over ice)
   ENDDO

   ! =================================================================== !
   ! Pack U surface currents; need to be rotated to geographical grid
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=fesom%dynamics%uvnode(1,1,n) ! (u/v,level,nod2D)
   ENDDO

   ! =================================================================== !
   ! Pack V surface currents; need to be rotated to geographical grid
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=fesom%dynamics%uvnode(2,1,n) ! (u/v,level,nod2D)
   ENDDO

   ! Rotate vectors (U,V) to geographical coordinates (r2g)
   DO n=1,myDim_nod2D
      rlon=coord_nod2D(1,n)
      rlat=coord_nod2D(2,n)
      CALL vector_r2g(zsendnf(n,nfield-1), zsendnf(n,nfield), rlon, rlat, 0) ! 0-flag for rot. coord
   ENDDO

   ! =================================================================== !
   ! Interpolate all fields
   IF (lparintmultatm) THEN
      CALL parinter_fld_mult( nfield, mype, npes, icomm, Ttogauss, &
         &                    myDim_nod2D, zsendnf, &
         &                    nopoints, zrecvnf )
   ELSE
      DO jf = 1, nfield
         CALL parinter_fld( mype, npes, icomm, Ttogauss, &
            &               myDim_nod2D, zsendnf(:,jf), &
            &               nopoints, zrecvnf(:,jf) )
      ENDDO
   ENDIF

   nfield = 0
   ! =================================================================== !
   ! Unpack 'pgsst' on Gauss.
   ! zsend(:)=a_ice(:)
   nfield = nfield + 1
   pgsst(:) = zrecvnf(:,nfield)
   !
   ! =================================================================== !
   ! Unpack 'pgifr' on Gauss.
   ! zsend(:)=a_ice(:)
   nfield = nfield + 1
   pgifr(:) = zrecvnf(:,nfield)
   !
   ! =================================================================== !
   ! Unpack ice temperature data (already in K)
   nfield = nfield + 1
   pgist(:) = zrecvnf(:,nfield)

   ! =================================================================== !
   ! Unpack ice albedo data pgalb on Gaussian grid.
   nfield = nfield + 1
   pgalb(:) = zrecvnf(:,nfield)

   ! =================================================================== !
   ! Unpack ice thickness data pghic on Gaussian grid.
   nfield = nfield + 1
   pghic(:) = zrecvnf(:,nfield)

   ! =================================================================== !
   ! Unpack snow thickness data pghsn on Gaussian grid.
   nfield = nfield + 1
   pghsn(:) = zrecvnf(:,nfield)

   ! =================================================================== !
   ! Unpack surface currents data pgucur/pgvcur on Gaussian grid.
   nfield = nfield + 1
   pgucur(:) = zrecvnf(:,nfield)
   nfield = nfield + 1
   pgvcur(:) = zrecvnf(:,nfield)

   ! Pack u(v) surface currents on elements
   !zsendnfUV(:,1)=fesom%dynamics%UV(1,1,1:myDim_elem2D)
   !zsendnfUV(:,2)=fesom%dynamics%UV(2,1,1:myDim_elem2D) !UV includes eDim, leave those away here
   !nfielduv = 2
   !
   !do elem=1, myDim_elem2D
   !
   !   ! compute element midpoints
   !   elnodes=elem2D_nodes(:,elem)
   !   rlon=sum(coord_nod2D(1,elnodes))/3.0_wpIFS
   !   rlat=sum(coord_nod2D(2,elnodes))/3.0_wpIFS
   !
   !   ! Rotate vectors to geographical coordinates (r2g)
   !   CALL vector_r2g(zsendnfUV(elem,1), zsendnfUV(elem,2), rlon, rlat, 0) ! 0-flag for rot. coord
   !
   !end do
   
#ifdef FESOM_TODO

   ! We need to sort out the non-unique global index before we
   ! can couple currents

   ! Interpolate: 'pgucur' and 'pgvcur' on Gaussian grid.
   IF (lparintmultatm) THEN
      CALL parinter_fld_mult( nfielduv, mype, npes, icomm, UVtogauss, &
         &                    myDim_nod2D, zsendnfUV, &
         &                    nopoints, zrecvnfUV )
   ELSE
      DO jf = 1, nfielduv
         CALL parinter_fld( mype, npes, icomm, UVtogauss, &
            &               myDim_nod2D, zsendnfUV(:,jf), &
            &               nopoints, zrecvnfUV(:,jf) )
      ENDDO
   ENDIF
   pgucur(:) = zrecvnfUV(:,1)
   pgvcur(:) = zrecvnfUV(:,2)
   
#else

   !pgucur(:) = 0.0
   !pgvcur(:) = 0.0

#endif

END SUBROUTINE nemogcmcoup_lim2_get


SUBROUTINE nemogcmcoup_exflds_get( mype, npes, icomm, &
   &                               nopoints, pgssh, pgmld, pg20d, pgsss, &
   &                               pgtem300, pgsal300 )

   ! Interpolate SSH, MLD, 20C isotherm, sea surface salinity, average T&S over upper 300m
   ! from the FESOM grid to IFS's Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind
   USE scripremap
   USE parinter
   USE interinfo
   USE fesom_main_storage_module, only: fesom => f
   USE o_ARRAYS, only : MLD1
   USE diagnostics, only : ldiag_extflds, zisotherm, saltzavg, tempzavg
   IMPLICIT NONE
   
   ! Arguments
   REAL(wpIFS), DIMENSION(nopoints) :: pgssh, pgmld, pg20d, pgsss, &
      & pgtem300, pgsal300
   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number Gaussian grid points
   INTEGER, INTENT(IN) :: nopoints

   ! Local variables
   INTEGER , PARAMETER :: maxnfield = 6
   INTEGER :: nfield = 0
   REAL(wpIFS), DIMENSION(fesom%partit%myDim_nod2D,maxnfield)  :: zsendnf
   REAL(wpIFS), DIMENSION(nopoints,maxnfield)  :: zrecvnf	
   real(kind=wpIFS), dimension(:,:), pointer :: coord_nod2D
   integer, pointer :: myDim_nod2D, eDim_nod2D

   ! Loop variables
   INTEGER :: n, elem, ierr, jf

   !#include "associate_mesh.h"
   ! associate what is needed only
   myDim_nod2D  => fesom%partit%myDim_nod2D
   eDim_nod2D   => fesom%partit%eDim_nod2D
   coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D) => fesom%mesh%coord_nod2D   
   

   nfield = 0
   ! =================================================================== !
   ! Pack SSH data 'pgssh' is on Gauss grid.
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=fesom%dynamics%eta_n(n)  ! in meters
   ENDDO

   ! =================================================================== !
   ! Pack MLD data
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=-MLD1(n) ! depth at which the density over depth differs 
        			 ! by 0.125 sigma units from the surface density (Griffies et al., 2009)
   ENDDO
   
   ! =================================================================== !
   ! Pack depth of 20C isotherm
   nfield = nfield + 1
   if (ldiag_extflds) then 
     DO n=1,myDim_nod2D
      zsendnf(n,nfield)=zisotherm(n) ! extra diagnostics active
     ENDDO
   else
     DO n=1,myDim_nod2D
      zsendnf(n,nfield)=-1. ! set to -1 as placeholder
     ENDDO
   end if

   ! =================================================================== !
   ! Pack sea surface salinity data: 'pgsss' on Gaussian grid.
   nfield = nfield + 1
   DO n=1,myDim_nod2D
      zsendnf(n,nfield)=fesom%tracers%data(2)%values(1,n) ! in psu
   ENDDO

   ! =================================================================== !
   ! Pack average temp over upper 300m: 'pgtem300' on Gaussian grid.
   nfield = nfield + 1
   if (ldiag_extflds) then
     DO n=1,myDim_nod2D
      zsendnf(n,nfield)=tempzavg(n) ! extra diagnostic active
     ENDDO
   else
     DO n=1,myDim_nod2D
      zsendnf(n,nfield)=-1. ! set to -1 as placeholder
     ENDDO
   end if

   ! =================================================================== !
   ! Pack average salinity over upper 300m: 'pgsal300' on Gaussian grid.
   nfield = nfield + 1
   if (ldiag_extflds) then
     DO n=1,myDim_nod2D
      zsendnf(n,nfield)=saltzavg(n) ! extra diagnostic active
     ENDDO
   else
     DO n=1,myDim_nod2D
      zsendnf(n,nfield)=-1. ! set to -1 as placeholder
     ENDDO
   end if

   ! =================================================================== !
   ! Interpolate all fields
   IF (lparintmultatm) THEN
      CALL parinter_fld_mult( nfield, mype, npes, icomm, Ttogauss, &
         &                    myDim_nod2D, zsendnf, &
         &                    nopoints, zrecvnf )
   ELSE
      DO jf = 1, nfield
         CALL parinter_fld( mype, npes, icomm, Ttogauss, &
            &               myDim_nod2D, zsendnf(:,jf), &
            &               nopoints, zrecvnf(:,jf) )
      ENDDO
   ENDIF

   nfield = 0
   ! =================================================================== !
   ! Unpack 'pgssh' on Gauss.
   nfield = nfield + 1
   pgssh(:) = zrecvnf(:,nfield)
   !
   ! =================================================================== !
   ! Unpack 'pgmld' on Gauss.
   nfield = nfield + 1
   pgmld(:) = zrecvnf(:,nfield)
   !
   ! =================================================================== !
   ! Unpack depth of 20C isotherm data
   nfield = nfield + 1
   pg20d(:) = zrecvnf(:,nfield)

   ! =================================================================== !
   ! Unpack sea surface salinity pgsss on Gaussian grid.
   nfield = nfield + 1
   pgsss(:) = zrecvnf(:,nfield)

   ! =================================================================== !
   ! Unpack average temp over upper 300m pgtem300 on Gaussian grid.
   nfield = nfield + 1
   pgtem300(:) = zrecvnf(:,nfield)

   ! =================================================================== !
   ! Unpack average salinity over upper 300m on Gaussian grid.
   nfield = nfield + 1
   pgsal300(:) = zrecvnf(:,nfield)


END SUBROUTINE nemogcmcoup_exflds_get


SUBROUTINE nemogcmcoup_lim2_update( mype, npes, icomm, &
   &                                npoints,  &
   &                                taux_oce, tauy_oce, taux_ice, tauy_ice, &
   &                                qs___oce, qs___ice, qns__oce, qns__ice, dqdt_ice, &
   &                                evap_tot, evap_ice, prcp_liq, prcp_sol, &
   &                                runoffIN, ocerunoff, tcc, lcc, tice_atm, &
   &                                kt, ldebug, loceicemix, lqnsicefilt )

   ! Update fluxes in nemogcmcoup_data by parallel
   ! interpolation of the input gaussian grid data
   
   USE par_kind !in ifs_modules.F90
   USE fesom_main_storage_module, only: fesom => f
   !USE o_PARAM, ONLY : WP, use wpIFS from par_kind (IFS)
   USE g_rotate_grid, 	only: vector_r2g, vector_g2r
   USE g_forcing_arrays, only: 	shortwave, prec_rain, prec_snow, runoff, & 
      			&	evap_no_ifrac, sublimation !'longwave' only stand-alone, 'evaporation' filled later
!    USE i_ARRAYS, 	only: stress_atmice_x, stress_atmice_y, oce_heat_flux, ice_heat_flux
!    USE i_ARRAYS, 	only: oce_heat_flux, ice_heat_flux 
   USE o_ARRAYS,        only: stress_atmoce_x, stress_atmoce_y
   USE g_comm_auto	! exchange_nod does the halo exchange
   
   ! all needed?
   USE parinter
   USE scripremap
   USE interinfo

   IMPLICIT NONE

   ! =================================================================== !
   ! Arguments ========================================================= !

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Fluxes on the Gaussian grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wpIFS), DIMENSION(npoints), INTENT(IN) :: &
      & taux_oce, tauy_oce, taux_ice, tauy_ice, &
      & qs___oce, qs___ice, qns__oce, qns__ice, &
      & dqdt_ice, evap_tot, evap_ice, prcp_liq, prcp_sol, &
      & runoffIN, ocerunoff, tcc, lcc, tice_atm

   ! Current time step
   INTEGER, INTENT(in) :: kt
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug
   ! QS/QNS mixed switch
   LOGICAL, INTENT(IN) :: loceicemix
   ! QNS ice filter switch (requires tice_atm to be sent)
   LOGICAL, INTENT(IN) :: lqnsicefilt

   !type(t_mesh), target :: mesh

   ! Local variables
   INTEGER		:: n, jf
   integer, pointer     :: myDim_nod2D, eDim_nod2D
   REAL(wpIFS), parameter 	:: rhofwt = 1000.        ! density of freshwater
   REAL(wpIFS), parameter  :: lfus = 333.7          ! latent heat of fusion [J/g]

   ! Packed send/receive buffers
   INTEGER , PARAMETER :: maxnfield = 11
   INTEGER :: nfield = 0
   REAL(wpIFS), DIMENSION(npoints,maxnfield) :: zsendnf
   REAL(wpIFS), DIMENSION(fesom%partit%myDim_nod2D,maxnfield) :: zrecvnf

   !#include "associate_mesh.h"
   ! associate only the necessary things
   real(kind=wpIFS), dimension(:,:), pointer :: coord_nod2D
   real(kind=wpIFS), dimension(:)  , pointer :: stress_atmice_x, stress_atmice_y
   real(kind=wpIFS), dimension(:)  , pointer :: oce_heat_flux, ice_heat_flux, a_ice
   real(kind=wpIFS), dimension(:)  , pointer :: enthalpyoffuse
   myDim_nod2D        => fesom%partit%myDim_nod2D
   eDim_nod2D         => fesom%partit%eDim_nod2D
   coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D) => fesom%mesh%coord_nod2D(:,:)  
   stress_atmice_x    => fesom%ice%stress_atmice_x
   stress_atmice_y    => fesom%ice%stress_atmice_y
   oce_heat_flux      => fesom%ice%atmcoupl%oce_flx_h(:)
   ice_heat_flux      => fesom%ice%atmcoupl%ice_flx_h(:)
   a_ice              => fesom%ice%data(1)%values(:)
   enthalpyoffuse     => fesom%ice%atmcoupl%enthalpyoffuse(:)
   ! =================================================================== !
   ! Sort out incoming arrays from the IFS and put them on the ocean grid

   ! TODO
   shortwave(:)=0.	! Done, updated below. What to do with shortwave over ice??
   !longwave(:)=0.	! Done. Only used in stand-alone mode.
   prec_rain(:)=0. 	! Done, updated below.
   prec_snow(:)=0.	! Done, updated below.
   evap_no_ifrac=0.	! Done, updated below. This is evap over ocean, does this correspond to evap_tot?
   sublimation=0.	! Done, updated below. 
   !
   ice_heat_flux=0. 	! Done. This is qns__ice currently. Is this the non-solar heat flux?    !   non solar heat fluxes below  !  (qns)
   oce_heat_flux=0. 	! Done. This is qns__oce currently. Is this the non-solar heat flux?
   !
   !runoff(:)=0.		! not used apparently. What is runoffIN, ocerunoff?
   !evaporation(:)=0.
   !ice_thermo_cpl.F90:  !---- total evaporation (needed in oce_salt_balance.F90)
   !ice_thermo_cpl.F90:  evaporation = evap_no_ifrac*(1.-a_ice) + sublimation*a_ice
   stress_atmice_x=0.   ! Done, taux_ice
   stress_atmice_y=0.   ! Done, tauy_ice
   stress_atmoce_x=0.   ! Done, taux_oce
   stress_atmoce_y=0.   ! Done, tauy_oce

   ! =================================================================== !
   ! Pack all arrays
   nfield = 0
   !1. Ocean solar radiation to T grid
   nfield = nfield + 1
   zsendnf(:,nfield) = qs___oce(:)

   ! =================================================================== !
   !2. Ice solar radiation to T grid
   ! DO NOTHING

   ! =================================================================== !
   !3. Ocean non-solar radiation to T grid (is this non-solar heat flux?)
   nfield = nfield + 1
   zsendnf(:,nfield) = qns__oce(:)

   ! =================================================================== !
   !4. Total flux over sea ice to T grid
   nfield = nfield + 1
   zsendnf(:,nfield) = qns__ice(:)+qs___ice(:)

   ! =================================================================== !
   !5. D(q)/dT to T grid
   ! DO NOTHING


   ! =================================================================== !
   !6. Total evaporation to T grid
   ! =================================================================== !
   !ice_thermo_cpl.F90:  total evaporation (needed in oce_salt_balance.F90)
   !ice_thermo_cpl.F90:  evaporation = evap_no_ifrac*(1.-a_ice) + sublimation*a_ice
   ! =================================================================== !
   nfield = nfield + 1
   zsendnf(:,nfield) = evap_tot(:)

   ! =================================================================== !
   !7. Sublimation (evaporation over ice) to T grid
   nfield = nfield + 1
   zsendnf(:,nfield) = evap_ice(:)

   ! =================================================================== !
   !8. Interpolate liquid precipitation to T grid
   nfield = nfield + 1
   zsendnf(:,nfield) = prcp_liq(:)

   ! =================================================================== !
   !9. Interpolate solid precipitation to T grid
   nfield = nfield + 1
   zsendnf(:,nfield) = prcp_sol(:)
   
   ! =================================================================== !
   !10. Interpolate runoff to T grid
   !
   !CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, runoff,  &
   !   &               myDim_nod2D, zrecv )
   !
   ! Unpack runoff, without halo
   !runoff(1:myDim_nod2D)=zrecv(1:myDim_nod2D) !conversion??
   !
   ! Do the halo exchange
   !call exchange_nod(runoff,fesom%partit)
   !
   !11. Interpolate ocean runoff to T grid
   !
   !CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, ocerunoff,  &
   !   &               myDim_nod2D, zrecv )
   !
   ! Unpack ocean runoff
   ! ??

   !12. Interpolate total cloud fractions to T grid (tcc)
   !
   !13. Interpolate low cloud fractions to T grid (lcc)


   ! =================================================================== !
   ! STRESSES

   ! OVER OCEAN:
   nfield = nfield + 1
   zsendnf(:,nfield) = taux_oce(:)

   nfield = nfield + 1
   zsendnf(:,nfield) = tauy_oce(:)

   ! =================================================================== !
   ! OVER ICE:
   nfield = nfield + 1
   zsendnf(:,nfield) = taux_ice(:)

   nfield = nfield + 1
   zsendnf(:,nfield) = tauy_ice(:)
   
   ! =================================================================== !
   ! Interpolate arrays
   IF (lparintmultatm) THEN
      CALL parinter_fld_mult( nfield, mype, npes, icomm, gausstoT, npoints, &
         &                    zsendnf, myDim_nod2D, &
         &                    zrecvnf )
   ELSE
      DO jf = 1, nfield
         CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, &
            &               zsendnf(:,jf), myDim_nod2D, &
            &               zrecvnf(:,jf) )
      ENDDO
   ENDIF
   
   ! =================================================================== !
   ! Unpack all arrays
   nfield = 0
   ! =================================================================== !
   !1. Unpack ocean solar radiation, without halo
   nfield = nfield + 1
   shortwave(1:myDim_nod2D)=zrecvnf(1:myDim_nod2D,nfield)

   ! Do the halo exchange
   call exchange_nod(shortwave,fesom%partit)


   ! =================================================================== !
   !2. Interpolate ice solar radiation to T grid
   ! DO NOTHING


   ! =================================================================== !
   !3. Unpack ocean non-solar, without halo
   nfield = nfield + 1
   oce_heat_flux(1:myDim_nod2D)=zrecvnf(1:myDim_nod2D,nfield)

   ! Do the halo exchange
   call exchange_nod(oce_heat_flux,fesom%partit)


   ! =================================================================== !
   !4. Unpack ice non-solar
   nfield = nfield + 1
   ice_heat_flux(1:myDim_nod2D)=zrecvnf(1:myDim_nod2D,nfield)
   where (a_ice<=1.e-12)
         ice_heat_flux=0.0
   end where
   ! Do the halo exchange
   call exchange_nod(ice_heat_flux,fesom%partit)


   ! =================================================================== !
   !5. D(q)/dT to T grid
   ! DO NOTHING


   ! =================================================================== !
   !6. Unpack total evaporation to T grid
   ! =================================================================== !
   !ice_thermo_cpl.F90:  total evaporation (needed in oce_salt_balance.F90)
   !ice_thermo_cpl.F90:  evaporation = evap_no_ifrac*(1.-a_ice) + sublimation*a_ice
   ! =================================================================== !
   ! Unpack total evaporation, without halo
   nfield = nfield + 1
   evap_no_ifrac(1:myDim_nod2D)=-zrecvnf(1:myDim_nod2D,nfield)/rhofwt	! kg m^(-2) s^(-1) -> m/s; change sign

   !7. Unpack sublimation (evaporation over ice), without halo
   nfield = nfield + 1
   sublimation(1:myDim_nod2D)=-zrecvnf(1:myDim_nod2D,nfield)/rhofwt	! kg m^(-2) s^(-1) -> m/s; change sign

   ! =================================================================== !
   sublimation(1:myDim_nod2D)  =sublimation(1:myDim_nod2D)*a_ice(1:myDim_nod2D)     ! sublimation   -> sublimation weighted with A
   evap_no_ifrac(1:myDim_nod2D)=evap_no_ifrac(1:myDim_nod2D)-sublimation(1:myDim_nod2D)            ! evap_no_ifrac -> evap weighted with (1-A)

   ! Do the halo exchange
   call exchange_nod(evap_no_ifrac,fesom%partit)
   call exchange_nod(sublimation,  fesom%partit)

   ! =================================================================== !
   !8. Unpack liquid precipitation, without halo
   nfield = nfield + 1
   prec_rain(1:myDim_nod2D)=zrecvnf(1:myDim_nod2D,nfield)/rhofwt	! kg m^(-2) s^(-1) -> m/s
   
   ! Do the halo exchange
   call exchange_nod(prec_rain,fesom%partit)


   ! =================================================================== !
   !9. Unpack solid precipitation, without halo
   nfield = nfield + 1
   prec_snow(1:myDim_nod2D)=zrecvnf(1:myDim_nod2D,nfield)/rhofwt	! kg m^(-2) s^(-1) -> m/s

   ! Do the halo exchange
   call exchange_nod(prec_snow,fesom%partit)


   ! =================================================================== !
   !10. Interpolate runoff to T grid
   !
   !CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, runoff,  &
   !   &               myDim_nod2D, zrecv )
   !
   ! Unpack runoff, without halo
   !runoff(1:myDim_nod2D)=zrecv(1:myDim_nod2D) !conversion??
   !
   ! Do the halo exchange
   !call exchange_nod(runoff,fesom%partit)
   !
   !11. Interpolate ocean runoff to T grid
   !
   !CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, ocerunoff,  &
   !   &               myDim_nod2D, zrecv )
   !
   ! Unpack ocean runoff
   ! ??

   !12. Interpolate total cloud fractions to T grid (tcc)
   !
   !13. Interpolate low cloud fractions to T grid (lcc)


   ! =================================================================== !
   ! STRESSES

   ! OVER OCEAN:
   nfield = nfield + 1
   ! Unpack x stress atm->oce, without halo; then do halo exchange
   stress_atmoce_x(1:myDim_nod2D)=zrecvnf(1:myDim_nod2D,nfield)
   call exchange_nod(stress_atmoce_x,fesom%partit)
   !
   ! Unpack y stress atm->oce, without halo; then do halo exchange
   nfield = nfield + 1
   stress_atmoce_y(1:myDim_nod2D)=zrecvnf(1:myDim_nod2D,nfield)
   call exchange_nod(stress_atmoce_y,fesom%partit)

   ! =================================================================== !
   ! OVER ICE:
   ! Unpack x stress atm->ice, without halo; then do halo exchange
   nfield = nfield + 1
   stress_atmice_x(1:myDim_nod2D)=zrecvnf(1:myDim_nod2D,nfield)
   call exchange_nod(stress_atmice_x,fesom%partit)
   !
   ! Unpack y stress atm->ice, without halo; then do halo exchange
   nfield = nfield + 1
   stress_atmice_y(1:myDim_nod2D)=zrecvnf(1:myDim_nod2D,nfield)
   call exchange_nod(stress_atmice_y,fesom%partit)


   ! =================================================================== !
   ! ROTATE VECTORS FROM GEOGRAPHIC TO FESOMS ROTATED GRID

   !if ((do_rotate_oce_wind .AND. do_rotate_ice_wind) .AND. rotated_grid) then
   do n=1, myDim_nod2D+eDim_nod2D
	   call vector_g2r(stress_atmoce_x(n), stress_atmoce_y(n), coord_nod2D(1, n), coord_nod2D(2, n), 0) !0-flag for rot. coord.
	   call vector_g2r(stress_atmice_x(n), stress_atmice_y(n), coord_nod2D(1, n), coord_nod2D(2, n), 0)
   end do
   !do_rotate_oce_wind=.false.
   !do_rotate_ice_wind=.false.
   !end if
   ! Enthalpy heat of fusion: take heat from the ocean in order to melt the snow that is falling into the ocean
   ! prec_snow*rho [kg/m2/s] * lfus [J/kg] = W/m2
   enthalpyoffuse(1:myDim_nod2D)= - rhofwt * prec_snow(1:myDim_nod2D) * lfus * 1000.
   call exchange_nod(enthalpyoffuse, fesom%partit)
END SUBROUTINE nemogcmcoup_lim2_update


SUBROUTINE nemogcmcoup_step( istp, icdate, ictime )

   USE g_clock, only: yearnew, month, day_in_month
   USE fesom_main_storage_module, only: fesom => f ! mype
   USE fesom_module, ONLY : fesom_runloop
   USE nemogcmcoup_steps, ONLY : substeps
   IMPLICIT NONE

   ! Arguments

   ! Time step
   INTEGER, INTENT(IN) :: istp

   ! Data and time from NEMO
   INTEGER, INTENT(OUT) :: icdate, ictime

   if(fesom%mype==0) then
   WRITE(0,*)'! IFS at timestep ', istp, '. Do ', substeps , 'FESOM timesteps...'
   endif
   CALL fesom_runloop(substeps)

   ! Compute date and time at the end of the time step

   icdate = yearnew*10000 + month*100 + day_in_month ! e.g. 20170906
   ictime = 0 ! (time is not used)

   if(fesom%mype==0) then
   WRITE(0,*)'! FESOM date at end of timestep is ', icdate ,' ======'
   endif

END SUBROUTINE nemogcmcoup_step


SUBROUTINE nemogcmcoup_final

   USE fesom_main_storage_module, only: fesom => f ! mype
   USE fesom_module, ONLY : fesom_finalize

   ! Finalize the FESOM model

   IMPLICIT NONE

   if(fesom%mype==0) then
   WRITE(*,*)'Finalization of FESOM from IFS.'
   endif
   CALL fesom_finalize

END SUBROUTINE nemogcmcoup_final
