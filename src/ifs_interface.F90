#if defined (__ifsinterface)
!=====================================================
! IFS interface for calling FESOM2 as a subroutine.
!
! -Original code for NEMO by Kristian Mogensen, ECMWF.
! -Adapted to FESOM2 by Thomas Rackow, AWI, 2018.
!-----------------------------------------------------

SUBROUTINE nemogcmcoup_init( icomm, inidate, initime, itini, itend, zstp, &
   & lwaveonly, iatmunit, lwrite )

   ! Initialize the FESOM model for single executable coupling 

   USE par_kind !in ifs_modules.F90
   USE g_PARSUP, only: MPI_COMM_FESOM
   USE g_config, only: dt
   USE g_clock, only: timenew, daynew, yearnew, month, day_in_month
   IMPLICIT NONE

   ! Input arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: icomm
   ! Initial date (e.g. 20170906), time, initial timestep and final time step
   INTEGER, INTENT(OUT) ::  inidate, initime, itini, itend
   ! Length of the time step
   REAL(wp), INTENT(OUT) :: zstp

   ! inherited from interface to NEMO, not used here:
   ! Coupling to waves only
   LOGICAL, INTENT(IN) :: lwaveonly
   ! Logfile unit (used if >=0)
   INTEGER :: iatmunit
   ! Write to this unit
   LOGICAL :: lwrite


   WRITE(0,*)'!======================================'
   WRITE(0,*)'! FESOM is initialized from within IFS.'

   WRITE(0,*)'! get MPI_COMM_FESOM. ================='
   MPI_COMM_FESOM=icomm

   itini = 1
   CALL main_initialize(itend) !also sets mype and npes
   WRITE(0,*)'! main_initialize done. ==============='

   ! Set more information for the caller
   
   ! initial date and time (time is not used)
   inidate = yearnew*10000 + month*100 + day_in_month ! e.g. 20170906
   initime = 0
   WRITE(0,*)'! FESOM initial date is ', inidate ,' ======'
   
   ! fesom timestep
   zstp = dt
   WRITE(0,*)'! FESOM timestep is ', real(zstp,4), 'sec'
   WRITE(0,*)'!======================================'

END SUBROUTINE nemogcmcoup_init


SUBROUTINE nemogcmcoup_coupinit( mypeIN, npesIN, icomm, &
   &                             npoints, nlocmsk, ngloind )

   ! FESOM modules
   USE g_PARSUP, only: mype, npes, myDim_nod2D, myDim_elem2D, myList_nod2D, myList_elem2D
   USE o_MESH,   only: nod2D, elem2D

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
      &                 lparbcast

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


   ! here FESOM knows about the (total number of) MPI tasks

   if(mype==0) then
   write(*,*) 'MPI has been initialized in the atmospheric model'
   write(*, *) 'Running on ', npes, ' PEs'
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

   OPEN(9,file='namfesomcoup.in')
   READ(9,namfesomcoup)
   CLOSE(9)

   ! Global number of Gaussian gridpoints

   CALL mpi_allreduce( npoints, nglopoints, 1, &
      &                mpi_integer, mpi_sum, icomm, ierr)


   WRITE(0,*)'!======================================'
   WRITE(0,*)'! SCALARS ============================='

   WRITE(0,*)'Update FESOM global scalar points'
   noglopoints=nod2D
   nopoints=myDim_nod2d

   ! Ocean mask and global indicies
   
   ALLOCATE(omask(MAX(nopoints,1)),ogloind(MAX(nopoints,1)))
   omask(:)= 1			! all points are ocean points
   ogloind = myList_nod2D	! global index for local point number

   ! Could be helpful later:
   ! Replace global numbering with a local one
   ! tmp(1:nod2d)=0
   ! DO n=1, myDim_nod2D+eDim_nod2D
   ! tmp(myList_nod2D(n))=n

   ! Read the interpolation weights and setup the parallel interpolation
   ! from atmosphere Gaussian grid to ocean T-grid
   
   IF (lreaddist) THEN
      CALL parinter_read( mype, npes, nglopoints, noglopoints, gausstoT,  &
         & cdpathdist,'ifs_to_fesom_gridT',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_gauss_to_T,remap_gauss_to_T,&
            & mype,npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_gauss_to_T,remap_gauss_to_T)
      ENDIF
      CALL parinter_init( mype, npes, icomm, &
         & npoints, nglopoints, nlocmsk, ngloind, &
         & nopoints, noglopoints, omask, ogloind, & 
         & remap_gauss_to_T, gausstoT, lcommout, TRIM(commoutprefix)//'_gtoT', &
         & iunit )
      CALL scripremap_dealloc(remap_gauss_to_T)
      IF (lwritedist) THEN
         CALL parinter_write( mype, npes, nglopoints, noglopoints, gausstoT,  &
            & cdpathdist,'ifs_to_fesom_gridT')
      ENDIF
   ENDIF

   ! From ocean T-grid to atmosphere Gaussian grid

   IF (lreaddist) THEN
      CALL parinter_read( mype, npes, noglopoints, nglopoints, Ttogauss,  &
         & cdpathdist,'fesom_gridT_to_ifs',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_T_to_gauss,remap_T_to_gauss,&
            & mype,npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_T_to_gauss,remap_T_to_gauss)
      ENDIF

      CALL parinter_init( mype, npes, icomm, &
         & nopoints, noglopoints, omask, ogloind, & 
         & npoints, nglopoints, nlocmsk, ngloind, &
         & remap_T_to_gauss, Ttogauss, lcommout, TRIM(commoutprefix)//'_Ttog', &
         & iunit )
      CALL scripremap_dealloc(remap_T_to_gauss)
      IF (lwritedist) THEN
         CALL parinter_write( mype, npes, noglopoints, nglopoints, Ttogauss,  &
            & cdpathdist,'fesom_gridT_to_ifs')
      ENDIF
   ENDIF
   
   DEALLOCATE(omask,ogloind)


   WRITE(0,*)'!======================================'
   WRITE(0,*)'! VECTORS ============================='

   WRITE(0,*)'Update FESOM global vector points'
   noglopoints=elem2D
   nopoints=myDim_elem2D

   ! Ocean mask and global indicies
   
   ALLOCATE(omask(MAX(nopoints,1)),ogloind(MAX(nopoints,1)))

   omask(:)= 1			! all elements are in the ocean 	
   ogloind = myList_elem2D	! global index for local element number

   ! Read the interpolation weights and setup the parallel interpolation
   ! from atmosphere Gaussian grid to ocean UV-grid
   
   IF (lreaddist) THEN
      CALL parinter_read( mype, npes, nglopoints, noglopoints, gausstoUV,  &
         & cdpathdist,'ifs_to_fesom_gridUV',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_gauss_to_UV,remap_gauss_to_UV,&
            & mype,npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_gauss_to_UV,remap_gauss_to_UV)
      ENDIF
      CALL parinter_init( mype, npes, icomm, &
         & npoints, nglopoints, nlocmsk, ngloind, &
         & nopoints, noglopoints, omask, ogloind, & 
         & remap_gauss_to_UV, gausstoUV, lcommout, TRIM(commoutprefix)//'_gtoUV', &
         & iunit )
      CALL scripremap_dealloc(remap_gauss_to_UV)
      IF (lwritedist) THEN
         CALL parinter_write( mype, npes, nglopoints, noglopoints, gausstoUV,  &
            & cdpathdist,'ifs_to_fesom_gridUV')
      ENDIF
   ENDIF

   ! From ocean UV-grid to atmosphere Gaussian grid

   IF (lreaddist) THEN
      CALL parinter_read( mype, npes, noglopoints, nglopoints, UVtogauss,  &
         & cdpathdist,'fesom_gridUV_to_ifs',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_UV_to_gauss,remap_UV_to_gauss,&
            & mype,npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_UV_to_gauss,remap_UV_to_gauss)
      ENDIF

      CALL parinter_init( mype, npes, icomm, &
         & nopoints, noglopoints, omask, ogloind, & 
         & npoints, nglopoints, nlocmsk, ngloind, &
         & remap_UV_to_gauss, UVtogauss, lcommout, TRIM(commoutprefix)//'_UVtog', &
         & iunit )
      CALL scripremap_dealloc(remap_UV_to_gauss)
      IF (lwritedist) THEN
         CALL parinter_write( mype, npes, noglopoints, nglopoints, UVtogauss,  &
            & cdpathdist,'fesom_gridUV_to_ifs')
      ENDIF
   ENDIF
   
   DEALLOCATE(omask,ogloind)
         
END SUBROUTINE nemogcmcoup_coupinit


SUBROUTINE nemogcmcoup_lim2_get( mype, npes, icomm, &
   &                             nopoints, pgsst, pgist, pgalb, &
   &                             pgifr, pghic, pghsn, pgucur, pgvcur, &
   &                             pgistl, licelvls )

   ! Interpolate sst, ice: surf T; albedo; concentration; thickness,
   ! snow thickness and currents from the FESOM grid to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind ! in ifs_modules.F90
   USE o_ARRAYS, ONLY : tr_arr, UV
   USE i_arrays, ONLY : m_ice, a_ice, m_snow
   USE i_therm_param, ONLY : tmelt
   USE g_PARSUP, only: myDim_nod2D, myDim_elem2D
   USE o_MESH, only: elem2D_nodes, coord_nod2D
   USE g_rotate_grid, only: vector_r2g
   USE parinter
   USE scripremap
   USE interinfo

   IMPLICIT NONE
   
   ! Arguments
   REAL(wp), DIMENSION(nopoints) :: pgsst, pgist, pgalb, pgifr, pghic, pghsn, pgucur, pgvcur
   REAL(wp), DIMENSION(nopoints,3) :: pgistl
   LOGICAL :: licelvls

   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number Gaussian grid points
   INTEGER, INTENT(IN) :: nopoints

   ! Local variables
   REAL(wp), DIMENSION(myDim_nod2D)  :: zsend
   REAL(wp), DIMENSION(myDim_elem2D) :: zsendU, zsendV
   INTEGER			     :: elnodes(3)
   REAL(wp)			     :: rlon, rlat	

   ! Loop variables
   INTEGER :: n, elem


   ! =================================================================== !
   ! Pack SST data and convert to K. 'pgsst' is on Gauss grid.
   do n=1,myDim_nod2D
      zsend(n)=tr_arr(1, n, 1)+tmelt	! sea surface temperature [K], 
					! (1=surface, n=node, 1/2=T/S)
   enddo

   ! Interpolate SST
   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               myDim_nod2D, zsend, &
      &               nopoints, pgsst )


   ! =================================================================== !
   ! Pack ice fraction data [0..1] and interpolate: 'pgifr' on Gauss.
   ! zsend(:)=a_ice(:)
   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               myDim_nod2D, a_ice, &
      &               nopoints, pgifr )


   ! =================================================================== !
   ! Pack ice temperature data (already in K)
   zsend(:)=273.15

   ! Interpolate ice surface temperature: 'pgist' on Gaussian grid. 
   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               myDim_nod2D, zsend, &
      &               nopoints, pgist )


   ! =================================================================== !
   ! Pack ice albedo data and interpolate: 'pgalb' on Gaussian grid.
   zsend(:)=0.7
   
   ! Interpolate ice albedo
   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               myDim_nod2D, zsend, &
      &               nopoints, pgalb )


   ! =================================================================== !
   ! Pack ice thickness data and interpolate: 'pghic' on Gaussian grid.
   zsend(:)=m_ice(:)/max(a_ice(:),0.01) ! ice thickness (mean over ice)

   ! Interpolation of average ice thickness
   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               myDim_nod2D, zsend, &
      &               nopoints, pghic ) 


   ! =================================================================== !
   ! Pack snow thickness data and interpolate: 'pghsn' on Gaussian grid.
   zsend(:)=m_snow(:)/max(a_ice(:),0.01) ! snow thickness (mean over ice)

   ! Interpolation of snow thickness
   CALL parinter_fld( mype, npes, icomm, Ttogauss, &
      &               myDim_nod2D, zsend, &
      &               nopoints, pghsn )


   ! =================================================================== !
   ! Surface currents need to be rotated to geographical grid

   ! Pack u(v) surface currents
   zsendU(:)=UV(1,1,1:myDim_elem2D)
   zsendV(:)=UV(2,1,1:myDim_elem2D) !UV includes eDim, leave those away here

   do elem=1, myDim_elem2D

      ! compute element midpoints
      elnodes=elem2D_nodes(:,elem)
      rlon=sum(coord_nod2D(1,elnodes))/3.0_WP
      rlat=sum(coord_nod2D(2,elnodes))/3.0_WP

      ! Rotate vectors to geographical coordinates (r2g)
      call vector_r2g(zsendU(elem), zsendV(elem), rlon, rlat, 0) ! 0-flag for rot. coord

   end do

   ! Interpolate: 'pgucur' and 'pgvcur' on Gaussian grid.
   CALL parinter_fld( mype, npes, icomm, UVtogauss, &
      &               myDim_elem2D, zsendU, &
      &               nopoints, pgucur )

   CALL parinter_fld( mype, npes, icomm, UVtogauss, &
      &               myDim_elem2D, zsendV, &
      &               nopoints, pgvcur )

#ifndef FESOM_TODO

   WRITE(0,*)'Everything implemented except ice level temperatures (licelvls).'

#else

   ! Ice level temperatures

   IF (licelvls) THEN

#if defined key_lim2

      DO jl = 1, 3
         
         ! Pack ice temperatures data at level jl(already in K)
         
         jk = 0 
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               jk = jk + 1
               zsend(jk) = tbif (ji,jj,jl)
            ENDDO
         ENDDO
         
         ! Interpolate ice temperature  at level jl
         
         CALL parinter_fld( mype, npes, icomm, Ttogauss, &
            &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zsend, &
            &               nopoints, pgistl(:,jl) )
         
      ENDDO

#else
      WRITE(0,*)'licelvls needs to be sorted for LIM3'
      CALL abort
#endif     

   ENDIF

   IF(nn_timing == 1) CALL timing_stop('nemogcmcoup_lim2_get')
   IF(lhook) CALL dr_hook('nemogcmcoup_lim2_get',1,zhook_handle)

#endif

END SUBROUTINE nemogcmcoup_lim2_get


SUBROUTINE nemogcmcoup_lim2_update( mype, npes, icomm, &
   &                                npoints,  &
   &                                taux_oce, tauy_oce, taux_ice, tauy_ice, &
   &                                qs___oce, qs___ice, qns__oce, qns__ice, dqdt_ice, &
   &                                evap_tot, evap_ice, prcp_liq, prcp_sol, &
   &                                runoff, ocerunoff, tcc, lcc, tice_atm, &
   &                                kt, ldebug, loceicemix, lqnsicefilt )

   ! Update fluxes in nemogcmcoup_data by parallel
   ! interpolation of the input gaussian grid data
   
   USE par_kind !in ifs_modules.F90
   USE g_PARSUP, only: myDim_nod2D, myDim_elem2D
   USE o_MESH, only: elem2D_nodes, coord_nod2D
   USE g_rotate_grid, only: vector_r2g

   IMPLICIT NONE

   ! =================================================================== !
   ! Arguments ========================================================= !

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Fluxes on the Gaussian grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wp), DIMENSION(npoints), INTENT(IN) :: &
      & taux_oce, tauy_oce, taux_ice, tauy_ice, &
      & qs___oce, qs___ice, qns__oce, qns__ice, &
      & dqdt_ice, evap_tot, evap_ice, prcp_liq, prcp_sol, &
      & runoff, ocerunoff, tcc, lcc, tice_atm

   ! Current time step
   INTEGER, INTENT(in) :: kt
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug
   ! QS/QNS mixed switch
   LOGICAL, INTENT(IN) :: loceicemix
   ! QNS ice filter switch (requires tice_atm to be sent)
   LOGICAL, INTENT(IN) :: lqnsicefilt

   ! Local variables

   ! Packed receive buffer
   REAL(wp), DIMENSION(myDim_nod2D) :: zrecv
   REAL(wp), DIMENSION(myDim_elem2D):: zrecvU, zrecvV


   ! =================================================================== !

#ifdef FESOM_TODO

   ! Packed receive buffer
   REAL(wp), DIMENSION((nlei-nldi+1)*(nlej-nldj+1)) :: zrecv
   ! Unpacked fields on ORCA grids
   REAL(wp), DIMENSION(jpi,jpj) :: zqs___oce, zqs___ice, zqns__oce, zqns__ice
   REAL(wp), DIMENSION(jpi,jpj) :: zdqdt_ice, zevap_tot, zevap_ice, zprcp_liq, zprcp_sol
   REAL(wp), DIMENSION(jpi,jpj) :: zrunoff, zocerunoff
   REAL(wp), DIMENSION(jpi,jpj) :: ztmp, zicefr
   ! Arrays for rotation
   REAL(wp), DIMENSION(jpi,jpj) :: zuu,zvu,zuv,zvv,zutau,zvtau 
   ! Lead fraction for both LIM2/LIM3
   REAL(wp), DIMENSION(jpi,jpj) :: zfrld
   ! Mask for masking for I grid
   REAL(wp) :: zmsksum
   ! For summing up LIM3 contributions to ice temperature
   REAL(wp) :: zval,zweig

   ! Loop variables
   INTEGER :: ji,jj,jk,jl
   ! netCDF debugging output variables
   CHARACTER(len=128) :: cdoutfile
   INTEGER :: inum
   REAL(wp) :: zhook_handle ! Dr Hook handle

   IF(lhook) CALL dr_hook('nemogcmcoup_lim2_update',0,zhook_handle)
   IF(nn_timing == 1) CALL timing_start('nemogcmcoup_lim2_update')
   
   ! Allocate the storage data

   IF (.NOT.lallociceflx) THEN
      ALLOCATE( &
         & zsqns_tot(jpi,jpj),   &
         & zsqns_ice(jpi,jpj),   &
         & zsqsr_tot(jpi,jpj),   &
         & zsqsr_ice(jpi,jpj),   &
         & zsemp_tot(jpi,jpj),   &
         & zsemp_ice(jpi,jpj),   &
	 & zsevap_ice(jpi,jpj),  &
         & zsdqdns_ice(jpi,jpj), &
         & zssprecip(jpi,jpj),   &
	 & zstprecip(jpi,jpj),   &
         & zstcc(jpi,jpj),       &
         & zslcc(jpi,jpj),       &
         & zsatmist(jpi,jpj),    &
         & zsqns_ice_add(jpi,jpj)&
         & )
      lallociceflx = .TRUE.
   ENDIF
   IF (.NOT.lallocstress) THEN
      ALLOCATE( &
         & zsutau(jpi,jpj),     &
         & zsvtau(jpi,jpj),     &
         & zsutau_ice(jpi,jpj), &
         & zsvtau_ice(jpi,jpj)  &
         & )
      lallocstress = .TRUE.
   ENDIF

   ! Sort out incoming arrays from the IFS and put them on the ocean grid
   
   !1. Interpolate ocean solar radiation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, qs___oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean solar radiation

   zqs___oce(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zqs___oce(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !2. Interpolate ice solar radiation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, qs___ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ice solar radiation

   zqs___ice(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zqs___ice(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !3. Interpolate ocean non-solar radiation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, qns__oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean non-solar radiation

   zqns__oce(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zqns__oce(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !4. Interpolate ice non-solar radiation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, qns__ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ice non-solar radiation

   zqns__ice(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zqns__ice(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !5. Interpolate  D(q)/dT to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, dqdt_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack D(q)/D(T) 

   zdqdt_ice(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zdqdt_ice(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !6. Interpolate total evaporation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, evap_tot,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack total evaporation

   zevap_tot(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zevap_tot(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !7. Interpolate evaporation over ice to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, evap_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack evaporation over ice

   zevap_ice(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zevap_ice(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
 
   !8. Interpolate liquid precipitation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, prcp_liq,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack liquid precipitation

   zprcp_liq(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zprcp_liq(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !9. Interpolate solid precipitation to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, prcp_sol,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack precipitation over ice

   zprcp_sol(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zprcp_sol(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   !10. Interpolate runoff to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, runoff,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack runoff

   zrunoff(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zrunoff(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !11. Interpolate ocean runoff to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, ocerunoff,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean runoff

   zocerunoff(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zocerunoff(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !12. Interpolate total cloud fractions to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, tcc,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean runoff

   zstcc(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zstcc(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   !13. Interpolate low cloud fractions to T grid

   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, lcc,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ocean runoff

   zslcc(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zslcc(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! get sea ice fraction and lead fraction

#if defined key_lim2
   zfrld(:,:) = frld(:,:)
   zicefr(:,:) = 1 - zfrld(:,:)
#else
   zicefr(:,:) = 0.0_wp
   DO jl = 1, jpl
      zicefr(:,:) = zicefr(:,:) + a_i(:,:,jl)
   ENDDO
   zfrld(:,:) = 1 - zicefr(:,:)
#endif
    
   zsemp_tot(:,:) = zevap_tot(:,:) - zprcp_liq(:,:) - zprcp_sol(:,:)
   zstprecip(:,:) = zprcp_liq(:,:) + zprcp_sol(:,:)
   ! More consistent with NEMO, but does changes the results, so
   ! we don't do it for now.
   ! zsemp_tot(:,:) = zevap_tot(:,:) - zstprecip(:,:)
   zsemp_ice(:,:) = zevap_ice(:,:) - zprcp_sol(:,:)
   zssprecip(:,:) = - zsemp_ice(:,:)
   zsemp_tot(:,:) = zsemp_tot(:,:) - zrunoff(:,:)
   zsemp_tot(:,:) = zsemp_tot(:,:) - zocerunoff(:,:)
   zsevap_ice(:,:) = zevap_ice(:,:)

   !   non solar heat fluxes   !  (qns)
   IF (loceicemix) THEN
      zsqns_tot(:,:) =  zqns__oce(:,:)
   ELSE
      zsqns_tot(:,:) =  zfrld(:,:) * zqns__oce(:,:) + zicefr(:,:) * zqns__ice(:,:)
   ENDIF
   zsqns_ice(:,:) =  zqns__ice(:,:)
   ztmp(:,:) = zfrld(:,:) * zprcp_sol(:,:) * lfus  ! add the latent heat of solid precip. melting

   zsqns_tot(:,:) = zsqns_tot(:,:) - ztmp(:,:)    ! over free ocean 
   !      solar heat fluxes    !   (qsr)
  
   IF (loceicemix) THEN
      zsqsr_tot(:,:) =  zqs___oce(:,:)
   ELSE
      zsqsr_tot(:,:) =  zfrld(:,:) * zqs___oce(:,:) + zicefr(:,:) * zqs___ice(:,:)
   ENDIF
   zsqsr_ice(:,:) =  zqs___ice(:,:)
   
   IF( ln_dm2dc ) THEN   ! modify qsr to include the diurnal cycle
      zsqsr_tot(:,:) = sbc_dcy( zsqsr_tot(:,:) )
      zsqsr_ice(:,:) = sbc_dcy( zsqsr_ice(:,:) )
   ENDIF
  
   zsdqdns_ice(:,:) = zdqdt_ice(:,:)

   ! Apply lateral boundary condition
   
   CALL lbc_lnk(zsqns_tot, 'T', 1.0)
   CALL lbc_lnk(zsqns_ice, 'T', 1.0)
   CALL lbc_lnk(zsqsr_tot, 'T', 1.0)
   CALL lbc_lnk(zsqsr_ice, 'T', 1.0)
   CALL lbc_lnk(zsemp_tot, 'T', 1.0)
   CALL lbc_lnk(zsemp_ice, 'T', 1.0)
   CALL lbc_lnk(zsdqdns_ice, 'T', 1.0)
   CALL lbc_lnk(zssprecip, 'T', 1.0)
   CALL lbc_lnk(zstprecip, 'T', 1.0)
   CALL lbc_lnk(zstcc, 'T', 1.0)
   CALL lbc_lnk(zslcc, 'T', 1.0)

   ! Interpolate  atmospheric ice temperature to T grid
      
   CALL parinter_fld( mype, npes, icomm, gausstoT, npoints, tice_atm,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack atmospheric ice temperature
      
   zsatmist(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zsatmist(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   CALL lbc_lnk(zsatmist, 'T', 1.0)
   
   zsqns_ice_add(:,:) = 0.0_wp

   ! Use the dqns_ice filter

   IF (lqnsicefilt) THEN

      ! Add filtr to qns_ice
      
#if defined key_lim2 
      ztmp(:,:) = tn_ice(:,:,1)
#else
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zval=0.0
            zweig=0.0
            DO jl = 1, jpl
               zval = zval + tn_ice(ji,jj,jl) * a_i(ji,jj,jl)
               zweig = zweig + a_i(ji,jj,jl)
            ENDDO
            IF ( zweig > 0.0 ) THEN
               ztmp(ji,jj) = zval /zweig
            ELSE
               ztmp(ji,jj) = rt0
            ENDIF
         ENDDO
      ENDDO
      CALL lbc_lnk(ztmp, 'T', 1.0)
#endif

      WHERE ( zicefr(:,:) > .001_wp )
         zsqns_ice_add(:,:) = zsdqdns_ice(:,:) * ( ztmp(:,:) - zsatmist(:,:) )
      END WHERE

      zsqns_ice(:,:) = zsqns_ice(:,:) + zsqns_ice_add(:,:)
      
   ENDIF

   ! Interpolate u-stress to U grid

   CALL parinter_fld( mype, npes, icomm, gausstoU, npoints,taux_oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack u stress on U grid

   zuu(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zuu(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate v-stress to U grid

   CALL parinter_fld( mype, npes, icomm, gausstoU, npoints, tauy_oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack v stress on U grid

   zvu(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zvu(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate u-stress to V grid

   CALL parinter_fld( mype, npes, icomm, gausstoV, npoints,taux_oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack u stress on V grid

   zuv(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zuv(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate v-stress to V grid

   CALL parinter_fld( mype, npes, icomm, gausstoV, npoints, tauy_oce,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack v stress on V grid

   zvv(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zvv(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   ! Rotate stresses from en to ij and put u,v stresses on U,V grids
   
   CALL repcmo( zuu, zvu, zuv, zvv, zsutau, zsvtau )

   ! Apply lateral boundary condition on u,v stresses on the U,V grids

   CALL lbc_lnk( zsutau, 'U', -1.0 )
   CALL lbc_lnk( zsvtau, 'V', -1.0 )

   ! Interpolate ice u-stress to U grid

   CALL parinter_fld( mype, npes, icomm, gausstoU, npoints,taux_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack ice u stress on U grid

   zuu(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zuu(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate ice v-stress to U grid

   CALL parinter_fld( mype, npes, icomm, gausstoU, npoints, tauy_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )
   
   ! Unpack ice v stress on U grid

   zvu(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zvu(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate ice u-stress to V grid

   CALL parinter_fld( mype, npes, icomm, gausstoV, npoints,taux_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack ice u stress on V grid

   zuv(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zuv(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO

   ! Interpolate ice v-stress to V grid

   CALL parinter_fld( mype, npes, icomm, gausstoV, npoints, tauy_ice,  &
      &               ( nlei - nldi + 1 ) * ( nlej - nldj + 1 ), zrecv )

   ! Unpack ice v stress on V grid

   zvv(:,:) = 0.0
   DO jj = nldj, nlej
      DO ji = nldi, nlei
         jk = ( jj - nldj ) * ( nlei - nldi + 1 ) + ( ji - nldi + 1 )
         zvv(ji,jj) = zrecv(jk)
      ENDDO
   ENDDO
   
   ! Rotate stresses from en to ij and put u,v stresses on U,V grids

   CALL repcmo( zuu, zvu, zuv, zvv, zutau, zvtau )

   ! Apply lateral boundary condition on u,v stresses on the U,V grids

   CALL lbc_lnk( zutau, 'U', -1.0 )
   CALL lbc_lnk( zvtau, 'V', -1.0 )

#if defined key_lim2_vp

   ! Convert to I grid for LIM2 for key_lim_vp
   DO jj = 2, jpjm1                                   ! (U,V) ==> I
      DO ji = 2, jpim1   ! NO vector opt.
         zmsksum = umask(ji-1,jj,1) + umask(ji-1,jj-1,1) 
         zsutau_ice(ji,jj) = ( umask(ji-1,jj,1) * zutau(ji-1,jj) + &
            &                  umask(ji-1,jj-1,1) * zutau(ji-1,jj-1) )
         IF ( zmsksum > 0.0 ) THEN
            zsutau_ice(ji,jj) = zsutau_ice(ji,jj) / zmsksum
         ENDIF
         zmsksum = vmask(ji,jj-1,1) + vmask(ji-1,jj-1,1) 
         zsvtau_ice(ji,jj) = ( vmask(ji,jj-1,1) * zvtau(ji,jj-1) + &
            &                  vmask(ji-1,jj-1,1) * zvtau(ji-1,jj-1) )
         IF ( zmsksum > 0.0 ) THEN
            zsvtau_ice(ji,jj) = zsvtau_ice(ji,jj) / zmsksum
         ENDIF
      END DO
   END DO

#else
   
   zsutau_ice(:,:) = zutau(:,:)
   zsvtau_ice(:,:) = zvtau(:,:)

#endif

   CALL lbc_lnk( zsutau_ice, 'I', -1.0 )
   CALL lbc_lnk( zsvtau_ice, 'I', -1.0 )
   
   ! Optionally write files write the data on the ORCA grid via IOM.

   IF (ldebug) THEN      
      WRITE(cdoutfile,'(A,I8.8)') 'zsutau_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsutau' , zsutau )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsvtau_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsvtau' , zsvtau )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsutau_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsutau_ice' , zsutau_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsvtau_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsvtau_ice' , zsvtau_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqns_tot_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqns_tot' , zsqns_tot )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqns_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqns_ice' , zsqns_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqsr_tot_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqsr_tot' , zsqsr_tot )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqsr_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqsr_ice' , zsqsr_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsemp_tot_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsemp_tot' , zsemp_tot )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsemp_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsemp_ice' , zsemp_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsdqdns_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsdqdns_ice' , zsdqdns_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zssprecip_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zssprecip' , zssprecip )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zstprecip_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zstprecip' , zstprecip )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsevap_ice_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsevap_ice' , zsevap_ice )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zstcc_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zstcc' , zstcc )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zslcc_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zslcc' , zslcc )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsatmist_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsatmist' , zsatmist )
      CALL iom_close( inum )
      WRITE(cdoutfile,'(A,I8.8)') 'zsqns_ice_add_',kt
      CALL iom_open( TRIM(cdoutfile), inum, ldwrt = .TRUE., kiolib = jprstlib)
      CALL iom_rstput( kt, kt, inum, 'zsqns_ice_add' , zsqns_ice_add )
      CALL iom_close( inum )
   ENDIF

   IF(nn_timing == 1) CALL timing_stop('nemogcmcoup_lim2_update')
   IF(lhook) CALL dr_hook('nemogcmcoup_lim2_update',1,zhook_handle)

#else

   WRITE(0,*)'nemogcmcoup_lim2_update not done for FESOM yet'
   CALL abort

#endif

END SUBROUTINE nemogcmcoup_lim2_update


SUBROUTINE nemogcmcoup_step( istp, icdate, ictime )

   IMPLICIT NONE

   ! Arguments

   ! Time step
   INTEGER, INTENT(IN) :: istp

   ! Data and time from NEMO
   INTEGER, INTENT(OUT) :: icdate, ictime

   ! Local variables
   
   ! Advance the FESOM model 1 time step

   WRITE(0,*)'Insert FESOM step here.'

   ! Compute date and time at the end of the time step.

#ifdef FESOM_TODO
   iye = ndastp / 10000
   imo = ndastp / 100 - iye * 100
   ida = MOD( ndastp, 100 )
   CALL greg2jul( 0, 0, 0, ida, imo, iye, zjul )
   zjul = zjul + ( nsec_day + 0.5_wp * rdttra(1) ) / 86400.0_wp
   CALL jul2greg( iss, imm, ihh, ida, imo, iye, zjul )
   icdate = iye * 10000 + imo * 100 + ida
   ictime = ihh * 10000 + imm * 100 + iss
#endif

END SUBROUTINE nemogcmcoup_step


SUBROUTINE nemogcmcoup_final

   ! Finalize the NEMO model

   IMPLICIT NONE

   WRITE(*,*)'Insert call to finalization of FESOM'
   CALL abort

END SUBROUTINE nemogcmcoup_final
#endif
