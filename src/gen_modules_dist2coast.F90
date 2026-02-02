module dist2coast
    !___________________________________________________________________________
    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use o_param           
    use o_arrays
    use g_comm_auto 
    implicit none
    private
    public :: find_coast_pnts, build_coast_bins, nearest_coast_dist, compute_dist2coast
    
    !___________________________________________________________________________
    contains

    !
    !
    !___________________________________________________________________________
    subroutine compute_dist2coast(dist, mesh, partit)
        !___INPUT/OUTPUT VARIABLES______________________________________________
        type(t_mesh)  , intent(inout), target      :: mesh
        type(t_partit), intent(inout), target      :: partit
        real(kind=WP) , intent(out  ), allocatable :: dist(:)
        !___LOCAL VARIABLES_____________________________________________________
        real(kind=WP) , allocatable                :: coast_lon(:), coast_lat(:)
        integer                                    :: n_coast
        real(kind=WP)                              :: bin_dlon=10.0, bin_dlat=10.0
        integer       , allocatable                :: bin_head(:,:), bin_next(:), visited(:,:)
        integer                                    :: bin_nlon, bin_nlat, vis_stamp
        integer                                    :: elem, elnodes(3) 
        real(kind=WP)                              :: lon_e_deg, lat_e_deg
        
        ! find all coastal points across all partitions and combine them in 
        ! lon/lat array that is known to all partitions
        call find_coast_pnts(coast_lon, coast_lat, n_coast, mesh, partit)
        
        ! build coastal binning array with a linked list of all coastal points 
        ! within a list so that we do not need to loop everytime over all coastal
        ! points to compute the distance
        call build_coast_bins(coast_lon, coast_lat  , &
                              bin_dlon, bin_dlat    , &
                              bin_nlon, bin_nlat    , & 
                              bin_head, bin_next)
        
        ! i do this here with visited and vis_stamp to skip the initialisation 
        ! that would be neccessary if i would use a logical array for visited
        ! each elem iteration i would needto reset it with false
        allocate(visited(bin_nlon, bin_nlat))
        visited   = 0
        vis_stamp = 0
        
        ! compute distance for each local elem point to the coastline
        allocate(dist(partit%myDim_elem2D))
        dist = 0.0_WP
        
        do elem = 1, partit%myDim_elem2D
            vis_stamp = vis_stamp + 1
            
            elnodes = mesh%elem2d_nodes(:,elem)
            lon_e_deg = sum(mesh%geo_coord_nod2D(1,elnodes))/3.0 / rad   
            lat_e_deg = sum(mesh%geo_coord_nod2D(2,elnodes))/3.0 / rad
                        
            call nearest_coast_dist(lon_e_deg, lat_e_deg,                   &
                                    coast_lon, coast_lat,                   &
                                    bin_dlon, bin_dlat, bin_head, bin_next, &
                                    visited, vis_stamp,                     &
                                    dist(elem))
        end do ! --> do elem = 1, partit%myDim_elem2D
        
        deallocate(bin_head, bin_next)
        deallocate(coast_lon, coast_lat)
        deallocate(visited)
        
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
        type(t_mesh)  , intent(inout), target      :: mesh
        type(t_partit), intent(inout), target      :: partit
        real(kind=WP) , intent(out)  , allocatable :: coast_lon(:), coast_lat(:)
        integer       , intent(out)                :: n_coast
        !___LOCAL VARIABLES_____________________________________________________
        integer                                    :: pi, edge, ednodes(2), lcl_n_coast
        real(kind=WP), allocatable                 :: lcl_coast_lon(:), lcl_coast_lat(:)
        logical, allocatable                       :: lcl_coast_isT(:)
        integer, allocatable                       :: all_counts(:), all_displs(:)
        
        ! Identify local coastal vertices
        allocate(lcl_coast_isT(partit%myDim_nod2D))
        lcl_coast_isT = .false.
        do edge = 1, partit%myDim_edge2D
            ednodes = mesh%edges(:, edge)
            ! Check if this is a boundary edge
            if (partit%myList_edge2D(edge) > mesh%edge2D_in) then  ! boundary edge
                lcl_coast_isT(ednodes(1)) = .true.
                lcl_coast_isT(ednodes(2)) = .true.
            end if
        end do
        
        ! How many coastal boundary point are in this partition
        lcl_n_coast = COUNT(lcl_coast_isT)
        
        ! Determine coordiantes of local coastal points
        if (lcl_n_coast > 0) THEN
            allocate(lcl_coast_lon(lcl_n_coast))
            allocate(lcl_coast_lat(lcl_n_coast))
            lcl_coast_lon = PACK(mesh%coord_nod2D(1,1:partit%myDim_nod2D), lcl_coast_isT)
            lcl_coast_lat = PACK(mesh%coord_nod2D(2,1:partit%myDim_nod2D), lcl_coast_isT)
        else
            allocate(lcl_coast_lon(0))
            allocate(lcl_coast_lat(0))
        end if 
        
        ! Gather all counts from all processors into array
        allocate(all_counts(partit%npes))
        call MPI_ALLGATHER(lcl_n_coast, 1, MPI_INTEGER, all_counts, 1, MPI_INTEGER, &
                           partit%MPI_COMM_FESOM, partit%MPIERR) 
        
        ! Calculate displacements across processors
        allocate(all_displs(partit%npes))
        all_displs(1) = 0
        DO pi = 2, partit%npes
            all_displs(pi) = all_displs(pi-1) + all_counts(pi-1)
        ENDDO
        
        ! Total number of all coastal points
        n_coast = sum(all_counts)
            
        ! Gather all coastal coordinates
        allocate(coast_lon(n_coast))
        allocate(coast_lat(n_coast))
        call MPI_ALLGATHERV(lcl_coast_lon, lcl_n_coast, MPI_DOUBLE_PRECISION, &
                            coast_lon, all_counts, all_displs, &
                            MPI_DOUBLE_PRECISION, partit%MPI_COMM_FESOM, partit%MPIERR)
        
        call MPI_ALLGATHERV(lcl_coast_lat, lcl_n_coast, MPI_DOUBLE_PRECISION, &
                            coast_lat, all_counts, all_displs, &
                            MPI_DOUBLE_PRECISION, partit%MPI_COMM_FESOM, partit%MPIERR)
        
        !_______________________________________________________________________
        deallocate(lcl_coast_isT, lcl_coast_lon, lcl_coast_lat)
        deallocate(all_counts, all_displs)
        
    end subroutine find_coast_pnts


    
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
        nlon = int(360.0_WP/dlon)
        nlat = int(180.0_WP/dlat)
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
            loni = lon_bin(coast_lon(ci), dlon, nlon)
            lati = lat_bin(coast_lat(ci), dlat, nlat)
            next(ci)          = head(loni, lati)
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
        ib = int((l + 180.0_WP) / dlon) + 1
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
        a  = max(-90.0_WP, min(90.0_WP, lat))
        jb = int((a + 90.0_WP) / dlat) + 1
        if (jb < 1) jb = 1
        if (jb > nlat) jb = nlat
    end function lat_bin
    
    
    
    !
    !
    !___________________________________________________________________________
    pure real(WP) function wrap_lon(lon) result(lw)
        real(WP), intent(in) :: lon  ! degrees, any range
        lw = modulo(lon + 180.0_WP, 360.0_WP) - 180.0_WP
    end function wrap_lon
    
  

    !
    !
    !___________________________________________________________________________
    ! Find nearest coastal distance for one vertex by expanding bin rings
    ! If your bins are too fine, remote ocean expands r a bit until it finds coast.
    subroutine nearest_coast_dist(vlon, vlat, coast_lon, coast_lat, dlon, dlat, &
                                  head, next, visited, vis_stamp, dist)
        !___INPUT/OUTPUT VARIABLES______________________________________________
        real(WP), intent(in)   :: vlon, vlat
        real(WP), intent(in)   :: coast_lon(:), coast_lat(:)
        real(WP), intent(in)   :: dlon, dlat
        integer,  intent(in)   :: head(:,:), next(:)
        integer,  intent(inout):: visited(:,:), vis_stamp
        real(WP), intent(out)  :: dist
        !___LOCAL VARIABLES_____________________________________________________
        integer                :: nlon, nlat, lati0, loni0, r, loni, lati, idx, ii
        real(WP)               :: d, best
        
        ! number of lon/lat bins 
        nlon = size(head, 1)
        nlat = size(head, 2)
        
        ! vlon, vlat are the coordinates of the reference points for which we want
        ! to compute the smallest distance to the coast line 
        loni0 = lon_bin(vlon, dlon, nlon)
        lati0 = lat_bin(vlat, dlat, nlat)
        
        ! returning huge number best=1.7977E+308
        best = huge(1.0_WP)
        
        ! Expand rings r = 0,1,2,... until we have a candidate and further rings
        ! cannot beat current best (simple conservative stop rule).
        ! scan through concentric index rings around index location of reverence 
        ! point
        do r = 0, max(nlon, nlat)
            ! Scan all bins in the (2r+1)x(2r+1) neighborhood.
            do lati = max(1, lati0-r), min(nlat, lati0+r)
                do ii = -r, r
                    ! wrap longitude bins cyclically
                    loni = loni0 + ii
                    if (loni < 1)    loni = loni + nlon
                    if (loni > nlon) loni = loni - nlon
                    
                    ! only scan the outer ring don't double the scanning of the inner rings
                    if (visited(loni,lati) == vis_stamp) cycle
                    visited(loni,lati) = vis_stamp
                    
                    idx = head(loni, lati)
                    
                    do while (idx /= 0)
                        d = haversine_m(vlon, vlat, coast_lon(idx), coast_lat(idx))
                        if (d < best) best = d
                        idx = next(idx)
                    end do
                end do !--> do ii = -r, r
            end do ! --> lati = max(1, lati0-r), min(nlat, lati0+r)
            
            if (best < huge(1.0_WP)) then
                ! Conservative early stop:
                ! Minimum possible distance to points beyond this ring is at least ~ r*min(dlon,dlat) degrees.
                ! If that exceeds best, stop.
                !if ( (real(r+1,WP) * min(dlon*cos(vlat*rad),dlat) * rad * r_earth) > best ) exit
                if ( (real(r+1,WP) * dlat * rad * r_earth) > best ) exit
            end if
        end do ! --> do r = 0, max(nlon, nlat)
        
        dist = best
    end subroutine nearest_coast_dist


    !
    !
    !___________________________________________________________________________
    ! computes the great-circle distance between two lon/lat points on a sphere (Earth), 
    ! using the haversine formula, and returns the distance in meters.
    pure real(WP) function haversine_m(lon1, lat1, lon2, lat2) result(d)
        real(WP), intent(in) :: lon1, lat1, lon2, lat2  ! degrees
        real(WP) :: phi1, phi2, dlam, dphi, a
        phi1 = lat1*rad
        phi2 = lat2*rad
        dphi = (lat2 - lat1)*rad
        dlam = (wrap_lon(lon2 - lon1))*rad
        
        ! The standard haversine expression:
        !       a=sin⁡2(dphi/2)+cos⁡(phi1)cos⁡(phi2)sin⁡2(dlam/2)
        ! where phi is latitude, lam is longitude.
        a  = sin(dphi*0.5_WP)**2 + cos(phi1)*cos(phi2)*sin(dlam*0.5_WP)**2
        
        ! central angle is: c=2arcsin⁡(a)
        ! then distance on a sphere is: d=R⋅c
        d  = 2.0_WP*r_earth*asin(min(1.0_WP, sqrt(a)))
        
  end function haversine_m


end module dist2coast