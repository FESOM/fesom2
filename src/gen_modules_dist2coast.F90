module g_dist2coast
    !___________________________________________________________________________
    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use g_config, only: flag_debug, cyclic_length
    use o_param           
    use o_arrays
    use g_comm_auto 
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    private
    public :: find_coast_pnts, build_coast_bins, nearest_coast_dist, compute_dist2coast
    
    interface allgatherv_coast
        module procedure allgatherv_coast_real4
        module procedure allgatherv_coast_real8
    end interface

    !___________________________________________________________________________
    contains
   
    !
    !
    !___________________________________________________________________________
    ! compute distance in [m] for each elemental surface points to the closest
    ! coastal point across partition domains
    subroutine compute_dist2coast(dist, mesh, partit, modein)
        !___INPUT/OUTPUT VARIABLES______________________________________________
        real(kind=WP)   , intent(inout)            :: dist(:)
        type(t_mesh)    , intent(in)   , target    :: mesh
        type(t_partit)  , intent(inout), target    :: partit
        character(len=*), intent(in), optional     :: modein
        !___LOCAL VARIABLES_____________________________________________________
        real(kind=WP) , allocatable                :: coast_lon(:), coast_lat(:)
        integer                                    :: n_coast
        real(kind=WP)                              :: bin_dlon=10.0*rad, bin_dlat=10.0*rad
        integer       , allocatable                :: bin_head(:,:), bin_next(:)
        integer                                    :: bin_nlon, bin_nlat
        integer                                    :: elem, node, elnodes(3) 
        real(kind=WP)                              :: lon_e, lat_e, lonr(3), xmin
        character(len=10)                          :: mode
        real(kind=WP)                              :: t0, t1, t2, t3
        
        !_______________________________________________________________________
        mode = 'elem'
        if (present(modein)) mode=trim(modein)
        if (partit%mype==0 .and. (trim(mode) .eq. 'elem')) then 
            write(*,*) '     └> compute dist2coast for elements'
        elseif (partit%mype==0 .and. (trim(mode) .eq. 'node')) then
            write(*,*) '     └> compute dist2coast for vertices'
        end if 
        
        !_______________________________________________________________________
        ! find all coastal points across all partitions and combine them in 
        ! lon/lat array that is known to all partitions
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'          --> call find_coast_pnts'//achar(27)//'[0m'
        t0=MPI_Wtime()
        call find_coast_pnts(coast_lon, coast_lat, n_coast, mesh, partit)
        
        !_______________________________________________________________________
        ! build coastal binning array with a linked list of all coastal points 
        ! within a list so that we do not need to loop everytime over all coastal
        ! points to compute the distance
        if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'          --> call build_coast_bins'//achar(27)//'[0m'
        t1=MPI_Wtime()
        call build_coast_bins(coast_lon, coast_lat  , &
                              bin_dlon , bin_dlat   , &
                              bin_nlon , bin_nlat   , & 
                              bin_head , bin_next)
        
        !_______________________________________________________________________
        ! compute distance for each local elem point to the coastline
        ! allocate(dist(partit%myDim_elem2D)) --> dist should be allocated and 
        ! initialised outside
        t2=MPI_Wtime()
        if (trim(mode) .eq. 'elem') then 
            !___________________________________________________________________
            if ( .not. (size(dist) == partit%myDim_elem2D+partit%eDim_elem2D .or. &
                        size(dist) == partit%myDim_elem2D)) then
                if (partit%mype==0) then         
                    print *, achar(27)//'[33m'
                    write(*,*) '____________________________________________________________________'
                    write(*,*) ' ERROR: size(dist)=', size(dist), ' does not agree to element computation'
                    write(*,*) '____________________________________________________________________'
                    print *, achar(27)//'[0m'
                    write(*,*)
                end if 
            end if 
            
            !___________________________________________________________________
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'          --> call nearest_coast_dist'//achar(27)//'[0m'
            do elem = 1, partit%myDim_elem2D
                elnodes = mesh%elem2d_nodes(:,elem)
                lat_e = sum(mesh%geo_coord_nod2D(2,elnodes))/3.0_WP
                
                ! do properly computing of elem centroid lon position across 
                ! periodic boundary
                lonr  = mesh%geo_coord_nod2D(1, elnodes)
                xmin  = minval(lonr)
                do node=1,3
                    if(lonr(node)-xmin>=cyclic_length/2.0_WP) lonr(node)=lonr(node)-cyclic_length
                    if(lonr(node)-xmin<-cyclic_length/2.0_WP) lonr(node)=lonr(node)+cyclic_length
                end do
                lon_e = sum(lonr)/3.0_WP
                
                ! compute elem centoid distance to nearest coast line point
                call nearest_coast_dist(lon_e    , lat_e        , &
                                        coast_lon, coast_lat    , &
                                        bin_dlon , bin_dlat     , &
                                        bin_head , bin_next     , &
                                        dist(elem))
            end do ! --> do elem = 1, partit%myDim_elem2D
            
        elseif (trim(mode) .eq. 'node') then
            !___________________________________________________________________
            if ( .not. (size(dist) == partit%myDim_nod2D+partit%eDim_nod2D .or. &
                        size(dist) == partit%myDim_nod2D)) then
                if (partit%mype==0) then         
                    print *, achar(27)//'[33m'
                    write(*,*) '____________________________________________________________________'
                    write(*,*) ' ERROR: size(dist)=', size(dist), ' does not agree to vertice computation'
                    write(*,*) '____________________________________________________________________'
                    print *, achar(27)//'[0m'
                    write(*,*)
                end if 
            end if 
            
            !___________________________________________________________________
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'          --> call nearest_coast_dist'//achar(27)//'[0m'
            do node = 1, partit%myDim_nod2D
                lat_e = mesh%geo_coord_nod2D(2,node) ! convert in degree
                
                ! do properly computing of elem centroid lon position across 
                ! periodic boundary
                lon_e = mesh%geo_coord_nod2D(1,node) ! convert in degree
                
                ! compute elem centoid distance to nearest coast line point
                call nearest_coast_dist(lon_e    , lat_e        , &
                                        coast_lon, coast_lat    , &
                                        bin_dlon , bin_dlat     , &
                                        bin_head , bin_next     , &
                                        dist(node))
            end do ! --> do node = 1, partit%myDim_nod2D
            
        else
            if (partit%mype==0) then 
                print *, achar(27)//'[33m'
                write(*,*) '____________________________________________________________________'
                write(*,*) ' ERROR: This mode=', trim(mode), ' is not supported when computing distance '
                write(*,*) '        to coast line.'
                write(*,*) '____________________________________________________________________'
                print *, achar(27)//'[0m'
                write(*,*)
            end if 
        end if ! --> if (trim(mode) .eq. 'elem') then 
        t3=MPI_Wtime()
        
        !_______________________________________________________________________
        if (partit%mype==0) then 
            write(*,*) '        └> time find coast pnts:', t1-t0
            write(*,*) '        └> time build bins     :', t2-t1
            write(*,*) '        └> time comp. dist     :', t3-t2
        end if 
        
        !_______________________________________________________________________
        deallocate(bin_head, bin_next)
        deallocate(coast_lon, coast_lat)
        
    end subroutine compute_dist2coast
    
    
    
    !
    !
    !___________________________________________________________________________
    ! create array with all coastal lon/lat coordinate points. This array has to be 
    ! visible to all partitions is need to compute smallest distance of every vertice
    ! with respect to the coast
    ! coast_lon/lat only needed for distance computation will be deallocated 
    ! immediately after
    subroutine find_coast_pnts(coast_lon, coast_lat, n_coast, mesh, partit)
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_mesh)  , intent(in)   , target      :: mesh
        type(t_partit), intent(inout), target      :: partit
        real(kind=WP) , intent(out)  , allocatable :: coast_lon(:), coast_lat(:)
        integer       , intent(out)                :: n_coast
        !___LOCAL VARIABLES_____________________________________________________
        integer                                    :: p, edge, ednodes(2), lcl_n_coast, cnt
        real(kind=WP), allocatable                 :: lcl_coast_lon(:), lcl_coast_lat(:)
        logical, allocatable                       :: lcl_coast_isT(:)
        integer, allocatable                       :: all_counts(:), all_displs(:)
        
        ! Identify local coastal vertices
        allocate(lcl_coast_isT(partit%myDim_nod2D))
        lcl_coast_isT = .false.
        do edge = 1, partit%myDim_edge2D
            if (partit%myList_edge2D(edge) <= mesh%edge2D_in) cycle  ! boundary edge
            ednodes = mesh%edges(:, edge)
            ! Check if this is a boundary edge
            if (ednodes(1)<=partit%myDim_nod2D) lcl_coast_isT(ednodes(1)) = .true.
            if (ednodes(2)<=partit%myDim_nod2D) lcl_coast_isT(ednodes(2)) = .true.
        end do
        
        ! How many coastal boundary point are in this partition
        lcl_n_coast = 0
        do p = 1, partit%myDim_nod2D
            if (lcl_coast_isT(p)) lcl_n_coast = lcl_n_coast+1
        end do    
        
        ! Determine coordiantes of local coastal points
        if (lcl_n_coast > 0) THEN
            allocate(lcl_coast_lon(lcl_n_coast), lcl_coast_lat(lcl_n_coast))
            cnt=1
            do p=1,partit%myDim_nod2D
                if (lcl_coast_isT(p)) then
                    lcl_coast_lon(cnt)=mesh%geo_coord_nod2D(1,p) ! is in radian
                    lcl_coast_lat(cnt)=mesh%geo_coord_nod2D(2,p) ! is in radian
                cnt=cnt+1
                end if
            end do
        else
            allocate(lcl_coast_lon(0), lcl_coast_lat(0))
        end if 
        
        ! Gather all counts from all processors into array
        allocate(all_counts(partit%npes))
        all_counts = 0
        call MPI_ALLGATHER(lcl_n_coast, 1, MPI_INTEGER, all_counts, 1, MPI_INTEGER, &
                           partit%MPI_COMM_FESOM, partit%MPIERR) 
        
        ! Calculate displacements across processors
        allocate(all_displs(partit%npes))
        all_displs=0
        do p = 2, partit%npes
            all_displs(p) = all_displs(p-1) + all_counts(p -1)
        end do
        
        ! Total number of all coastal points
        n_coast = sum(all_counts)
        
        ! Gather all coastal coordinates from all partitions
        allocate(coast_lon(n_coast), coast_lat(n_coast))
        call allgatherv_coast(lcl_coast_lon, lcl_n_coast            , &
                              coast_lon, all_counts, all_displs     , &
                              partit)
        call allgatherv_coast(lcl_coast_lat, lcl_n_coast            , &
                              coast_lat, all_counts, all_displs     , &
                              partit)
        
        !_______________________________________________________________________
        deallocate(lcl_coast_isT, lcl_coast_lon, lcl_coast_lat)
        deallocate(all_counts, all_displs)
        
    end subroutine find_coast_pnts

    
    
    !
    !
    !___________________________________________________________________________
    ! make sure MPI_ALLGATHERV for all coastal points is precision selectiv by 
    ! its self through interface for the case WP is changed
    subroutine allgatherv_coast_real8(send, sendcount, recv, counts, displs, partit)
        real(real64)  , intent(in)              :: send(:)
        real(real64)  , intent(out)             :: recv(:)
        integer       , intent(in)              :: sendcount, counts(:), displs(:)
        type(t_partit), intent(inout), target   :: partit
        call MPI_ALLGATHERV(send, sendcount, MPI_DOUBLE_PRECISION       , &
                            recv, counts, displs, MPI_DOUBLE_PRECISION  , &
                            partit%MPI_COMM_FESOM, partit%MPIERR)
    end subroutine allgatherv_coast_real8
    
    subroutine allgatherv_coast_real4(send, sendcount, recv, counts, displs, partit)
        real(real32)  , intent(in)              :: send(:)
        real(real32)  , intent(out)             :: recv(:)
        integer       , intent(in)              :: sendcount, counts(:), displs(:)
        type(t_partit), intent(inout), target   :: partit
        call MPI_ALLGATHERV(send, sendcount, MPI_REAL       , &
                            recv, counts, displs, MPI_REAL  , &
                            partit%MPI_COMM_FESOM, partit%MPIERR)
    end subroutine allgatherv_coast_real4
    
    
    
    !
    !
    !___________________________________________________________________________
    ! Build bin linked-lists for coastal points
    ! Inputs: coast_lon/coast_lat (degrees), dlon/dlat (degrees)
    ! Outputs:
    !   nlon,nlat, head(nlon,nlat), next(ncoast)
    subroutine build_coast_bins(coast_lon, coast_lat, dlon, dlat, nlon, nlat, head, next)
        !___INPUT/OUTPUT VARIABLES______________________________________________
        real(WP), intent(in)  :: coast_lon(:), coast_lat(:)
        real(WP), intent(in)  :: dlon, dlat
        integer,  intent(out) :: nlon, nlat
        integer,  intent(out), allocatable :: head(:,:), next(:)
        !___LOCAL VARIABLES_____________________________________________________
        integer :: ci, loni, lati, ncoast
        
        ! number of coastal points 
        ncoast = size(coast_lon)
        
        ! number of lon/lat points depending on binning resolution
        nlon = int(2.0_WP*pi/dlon)
        nlat = int(       pi/dlat)
        if (nlon < 1) nlon = 1
        if (nlat < 1) nlat = 1
        
        ! allocate fields for linked list
        allocate(head(nlon, nlat))
        allocate(next(ncoast))
        head = 0
        next = 0
        
        !  fill linked list head and next array with indices, How does it work: 
        !  next array has length number of coastal point, a next array indices position
        !  next(i) contains the index of the next point in coast_lon/coast_lat
        !  head(ib,jb) contains only the starting index from the list from which you 
        !  then proceed to next --> to next --> to next index until you reach next(i)=0
        !  which is the end of the list at this bin location 
        do ci = 1, ncoast
            loni             = lon_bin(coast_lon(ci), dlon, nlon)
            lati             = lat_bin(coast_lat(ci), dlat, nlat)
            next(ci)         = head(loni, lati)
            head(loni, lati) = ci
        end do
    end subroutine build_coast_bins
    
    
    
    !
    !
    !___________________________________________________________________________
    ! pure declares a Fortran procedure as side-effect free: it cannot modify global 
    ! data or perform I/O and depends only on its arguments. It guarantees that 
    ! the same inputs always produce the same output (its rekursive). This makes 
    ! the procedure thread-safe and safe for use in DO CONCURRENT and parallel regions.
    ! While it doesn’t make the function itself faster, it enables safer optimization 
    ! and vectorization.
    pure integer function lon_bin(lon, dlon, nlon) result(ib)
        real(WP), intent(in) :: lon, dlon
        integer,  intent(in) :: nlon
        real(WP) :: l
        l  = wrap_lon(lon)
        ib = int((l + pi) / dlon) + 1
        if (ib < 1) ib = 1
        if (ib > nlon) ib = nlon
    end function lon_bin
    
    
    
    !
    !
    !___________________________________________________________________________
    pure integer function lat_bin(lat, dlat, nlat) result(jb)
        real(WP), intent(in) :: lat, dlat
        integer,  intent(in) :: nlat
        real(WP) :: a
        a  = max(-0.5_WP*pi, min(0.5_WP*pi, lat))
        jb = int((a + 0.5_WP*pi) / dlat) + 1
        if (jb < 1) jb = 1
        if (jb > nlat) jb = nlat
    end function lat_bin
    
    
    
    !
    !
    !___________________________________________________________________________
    pure real(WP) function wrap_lon(lon) result(lw)
        real(WP), intent(in) :: lon  ! degrees, any range
        lw = modulo(lon + pi, 2.0_WP*pi) - pi
    end function wrap_lon
    
  

    !
    !
    !___________________________________________________________________________
    ! Find nearest coastal distance for one vertex by expanding bin rings
    ! If your bins are too fine, remote ocean expands r a bit until it finds coast.
    subroutine nearest_coast_dist(vlon, vlat, coast_lon, coast_lat, dlon, dlat, &
                                  head, next, dist)
        !___INPUT/OUTPUT VARIABLES______________________________________________
        real(WP), intent(in)   :: vlon, vlat
        real(WP), intent(in)   :: coast_lon(:), coast_lat(:)
        real(WP), intent(in)   :: dlon, dlat
        integer,  intent(in)   :: head(:,:), next(:)
        real(WP), intent(out)  :: dist
        !___LOCAL VARIABLES_____________________________________________________
        integer                :: nlon, nlat, lati0, loni0, r
        real(WP)               :: dist_best, minstep, cos_lat
        
        ! number of lon/lat bins 
        nlon = size(head, 1)
        nlat = size(head, 2)
        
        ! vlon, vlat are the coordinates of the reference points for which we want
        ! to compute the smallest distance to the coast line 
        loni0   = lon_bin(vlon, dlon, nlon)
        lati0   = lat_bin(vlat, dlat, nlat)
        cos_lat = max(0.0_WP, cos(vlat))
        minstep = min(dlat, dlon*cos_lat)
        
        ! returning huge number dist_best=1.7977E+308
        dist_best = huge(1.0_WP)
        
        ! Expand rings r = 0,1,2,... until we have a candidate and further rings
        ! cannot beat current dist_best (simple conservative stop rule).
        ! scan through concentric index rings around index location of reverence 
        ! point
        !   r=0    r=1        r=2
        !                  .-------.
        !         .---.    |       | 
        !    #    | # |    |   #   |  ...
        !         `---´    |       |
        !                  `-------´
        !
        ! do not rescan bins that already have been scanned scan only new outer 
        ! shell 
        do r = 0, max(nlon, nlat)
            
            ! top row
            call scan_row(lati0+r, -r, r)
            
            if (r > 0) then
                ! bottom row (avoid double count when r=0)
                call scan_row(lati0-r, -r       , r        )
                
                ! left/right columns excluding corners (already done by rows)
                call scan_col(-r     , lati0-r+1, lati0+r-1)
                call scan_col( r     , lati0-r+1, lati0+r-1)
            end if
            
            ! Conservative early stop (safe globally): lon/lat spacing
            ! do not accept center value r=0 as automatically being the best
            if (dist_best < huge(1.0_WP) .and. r>0) then
                if ( real(r,WP) * minstep * r_earth > dist_best ) exit
            end if
            
        end do
        dist = dist_best
        
        !_______________________________________________________________________
        contains
            !
            !___________________________________________________________________
            !  make sure lon indices properly wrap across the periodic boiundary
            pure integer function wrap_loni(i) result(iw)
                integer, intent(in) :: i
                iw = i
                if (iw < 1)    iw = iw + nlon
                if (iw > nlon) iw = iw - nlon
            end function wrap_loni
            !
            !___________________________________________________________________
            ! scan outer row of binning shell with radius r
            subroutine scan_row(lati, ii1, ii2)
                integer, intent(in) :: lati, ii1, ii2
                integer :: ii, loni, idx
                real(WP) :: d
                if (lati < 1 .or. lati > nlat) return
                
                do ii = ii1, ii2
                    loni = wrap_loni(loni0 + ii)
                    idx  = head(loni, lati)
                    do while (idx /= 0)
                        d = haversine_m(vlon, vlat, coast_lon(idx), coast_lat(idx))
                        if (d < dist_best) dist_best = d
                        idx = next(idx)
                    end do
                end do
            end subroutine scan_row
            !
            !___________________________________________________________________
            ! scan outer column of binning shell with radius r
            subroutine scan_col(ii, lat1, lat2)
                integer, intent(in) :: ii, lat1, lat2
                integer :: lati, loni, idx
                real(WP) :: d
                loni = wrap_loni(loni0 + ii)
                
                do lati = max(1, lat1), min(nlat, lat2)
                    idx = head(loni, lati)
                    do while (idx /= 0)
                        d = haversine_m(vlon, vlat, coast_lon(idx), coast_lat(idx))
                        if (d < dist_best) dist_best = d
                        idx = next(idx)
                    end do
                end do
            end subroutine scan_col    
        
    end subroutine nearest_coast_dist


    
    !
    !
    !___________________________________________________________________________
    ! computes the great-circle distance between two lon/lat points on a sphere (Earth), 
    ! using the haversine formula, and returns the distance in meters.
    pure real(WP) function haversine_m(lon1, lat1, lon2, lat2) result(d)
        real(WP), intent(in) :: lon1, lat1, lon2, lat2  ! degrees
        real(WP) :: dlam, dphi, a
        dphi = (lat2 - lat1)
        dlam = (wrap_lon(lon2 - lon1))
        
        ! The standard haversine expression:
        !       a=sin⁡2(dphi/2)+cos⁡(phi1)cos⁡(phi2)sin⁡2(dlam/2)
        ! where phi is latitude, lam is longitude.
        a  = sin(dphi*0.5_WP)**2 + cos(lat1)*cos(lat2)*sin(dlam*0.5_WP)**2
        
        ! central angle is: c=2arcsin⁡(a)
        ! then distance on a sphere is: d=R⋅c
        d  = 2.0_WP*r_earth*asin(min(1.0_WP, sqrt(a)))
        
    end function haversine_m

end module g_dist2coast