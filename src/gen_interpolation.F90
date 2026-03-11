module g_interp
    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use o_PARAM, only: WP
    implicit none
    
    !___________________________________________________________________________
    interface interp_e2n
        module procedure interp_e2n_3d
        module procedure interp_e2n_2d
        module procedure interp_e2n_1d
    end interface interp_e2n
    
    !___________________________________________________________________________
    contains

    


! routines doing 3D, 2D and 1D interpolation
subroutine interp_2d_field_v2(num_lon_reg, num_lat_reg, lon_reg, lat_reg, data_reg, missvalue, &
     num_mod, lon_mod, lat_mod, data_mod, partit)
  !-------------------------------------------------------------------------------------
  ! A second version of 2D interpolation.
  ! This routine does 2d interpolation from a regular grid to specified nodes
  ! on the surface grid. The regular grid is assumed to be global. 
  ! In this version the missing value will be checked. If a model grid point is outside
  ! the ocean, the nearest value will be assigned.
  ! The output is data_mod, others are input.
  ! Arguments:
  ! num_lon_reg				number of regular grid points in the longitude direction
  ! num_lat_reg				number of regular grid points in the latitude direction
  ! lon_reg(num_lon_reg)		longitude of the regular grid points
  ! lat_reg(num_lat_reg)        	latitude of the regular grid points
  ! data_reg(num_lon_reg, num_lat_reg) 	data on the regular grid
  ! missvalue                           missing value in the raw data
  ! num_mod      			number of interpolation nodes
  ! lon_mod(num_mod)			longitude of interpolation nodes
  ! lat_mod(num_mod)			latitude of interpolation nodes
  ! data_mod(num_mod)			output data on interpolation nodes
  ! Unit of lon_reg, lat_reg, lon_mod, lat_mod: degree
  ! Order of lon_reg: monotonically increasing in the range [0 360]  
  !                   (range [180,360] is associated with longitude [-180, 0])
  ! Order of lat_reg: monotonically increasing in the range [-90 90]
  ! Make sure that 'lon_mod' is also in the [0, 360] range.
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------------------------------
  implicit none
  integer             		:: n, i, ii, jj, k, nod_find
  integer			:: ind_lat_h, ind_lat_l, ind_lon_h, ind_lon_l
  integer,        intent(in)   	:: num_lon_reg, num_lat_reg, num_mod
  type(t_partit), intent(inout) :: partit
  real(kind=WP) 			:: x, y, diff, d, dmin
  real(kind=WP)			:: rt_lat1, rt_lat2, rt_lon1, rt_lon2
  real(kind=WP)                  :: data(2,2)
  real(kind=WP)                  :: data_lo, data_up
  real(kind=WP), intent(in)	:: lon_reg(num_lon_reg), lat_reg(num_lat_reg)
  real(kind=WP), intent(in)	:: data_reg(num_lon_reg, num_lat_reg), missvalue
  real(kind=WP), intent(in)	:: lon_mod(num_mod), lat_mod(num_mod)
  real(kind=WP), intent(out)  	:: data_mod(num_mod)
  !
  if(lon_reg(1)<0.0 .or. lon_reg(num_lon_reg)>360.) then
    if (partit%mype==0) then 
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) ' ERROR: in 2D interpolation found either lon_reg<0.0 or lon_reg>360'
        write(*,*) '        The regular grid is not in the proper longitude range.'
        write(*,*) '        The range should be from [0...360-dlon] '
        write(*,*) '        --> something like: np.arange(0,360,dlon) '
        print *, achar(27)//'[0m'
        write(*,*)
    end if 
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
  end if

  do n=1,num_mod
     x=lon_mod(n)
     y=lat_mod(n)
     ! find the surrounding rectangular box and get interpolation ratios
     ! 1) north-south direction
     if(y<lat_reg(1)) then
        ind_lat_h=2
        ind_lat_l=1
        y=lat_reg(1)
     elseif(y>lat_reg(num_lat_reg)) then
        ind_lat_h=num_lat_reg
        ind_lat_l=num_lat_reg-1
        y=lat_reg(num_lat_reg)
     else
        do i=2,num_lat_reg
           if(lat_reg(i)>=y) then
              ind_lat_h=i
              ind_lat_l=i-1
              exit
           end if
        end do
     end if
     diff=lat_reg(ind_lat_h)-lat_reg(ind_lat_l)
     rt_lat1=(lat_reg(ind_lat_h)-y)/diff
     rt_lat2=1.0_WP-rt_lat1
     ! 2) east_west direction
     if(x<lon_reg(1)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360._WP-lon_reg(ind_lon_l))
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0_WP-rt_lon1
     elseif(x>lon_reg(num_lon_reg)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360._WP-lon_reg(ind_lon_l))
        rt_lon2=(x-lon_reg(ind_lon_l))/diff
!!PS         rt_lon1=1.0_WP-rt_lon2
     else
        do i=2,num_lon_reg
           if(lon_reg(i)>=x) then
              ind_lon_h=i
              ind_lon_l=i-1
              exit
           end if
        end do
        diff=lon_reg(ind_lon_h)-lon_reg(ind_lon_l)
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0_WP-rt_lon1
     end if
     !
     data(1,1)=data_reg(ind_lon_l,ind_lat_l)
     data(1,2)=data_reg(ind_lon_l,ind_lat_h)
     data(2,1)=data_reg(ind_lon_h,ind_lat_l)
     data(2,2)=data_reg(ind_lon_h,ind_lat_h)
     !
     ! interpolate data
     if(any(data==missvalue)) then
        dmin=10000.0_WP
        nod_find=0
        do k=1,5
           do ii=max(1,ind_lon_l-k+1),min(num_lon_reg,ind_lon_l+k)
              do jj=max(1,ind_lat_l-k+1),min(num_lat_reg,ind_lat_l+k)
                 if(data_reg(ii,jj)==missvalue) cycle
                 d=(x-lon_reg(ii))**2 + (y-lat_reg(jj))**2
                 if(d<dmin) then
                    nod_find=1
                    data_mod(n)=data_reg(ii,jj)
                    dmin=d
                 end if
              end do
           end do
           if(nod_find==1) exit
        end do
     else
        data_mod(n)=(data(1,1)*rt_lon1 + data(2,1)*rt_lon2)*rt_lat1 + &
             (data(1,2)*rt_lon1 + data(2,2)*rt_lon2)*rt_lat2
     end if
  end do
end subroutine interp_2d_field_v2
!
!---------------------------------------------------------------------------
!
subroutine interp_2d_field(num_lon_reg, num_lat_reg, lon_reg, lat_reg, data_reg, &
     num_mod, lon_mod, lat_mod, data_mod, phase_flag, partit)
  !-------------------------------------------------------------------------------------
  ! This routine does 2d interpolation from a regular grid to specified nodes
  ! on the surface grid. The regular grid is assumed to be global.
  ! This version assumes that no dummy values exist in the raw data.
  ! The argument 'phase_flag' is used to set if a phase angle field is to be interpolated.
  ! The output is data_mod, others are input.
  ! Arguments:
  ! num_lon_reg				number of regular grid points in the longitude direction
  ! num_lat_reg				number of regular grid points in the latitude direction
  ! lon_reg(num_lon_reg)		longitude of the regular grid points
  ! lat_reg(num_lat_reg)        	latitude of the regular grid points
  ! data_reg(num_lon_reg, num_lat_reg) 	data on the regular grid
  ! num_mod      			number of interpolation nodes
  ! lon_mod(num_mod)			longitude of interpolation nodes
  ! lat_mod(num_mod)			latitude of interpolation nodes
  ! data_mod(num_mod)			output data on interpolation nodes
  ! phase_flag                          1: interpolate phase angle (0-360°); 0: otherwise
  ! Unit of lon_reg, lat_reg, lon_mod, lat_mod: degree
  ! Order of lon_reg: monotonically increasing in the range [0 360]  
  !                   (range [180,360] is associated with longitude [-180, 0])
  ! Order of lat_reg: monotonically increasing in the range [-90 90]
  ! Make sure that 'lon_mod' is also in the [0, 360] range.
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------------------------------
  implicit none
  integer             		:: n, i
  integer			:: ind_lat_h, ind_lat_l, ind_lon_h, ind_lon_l
  integer, intent(in)         	:: num_lon_reg, num_lat_reg, num_mod
  integer, intent(in)          	:: phase_flag
  type(t_partit), intent(inout) :: partit
  real(kind=WP) 			:: x, y, diff
  real(kind=WP)			:: rt_lat1, rt_lat2, rt_lon1, rt_lon2
  real(kind=WP)                 :: data_ll, data_lh, data_hl, data_hh
  real(kind=WP)                 :: data_lo, data_up
  real(kind=WP), intent(in)	:: lon_reg(num_lon_reg), lat_reg(num_lat_reg)
  real(kind=WP), intent(in)	:: data_reg(num_lon_reg, num_lat_reg)
  real(kind=WP), intent(in)	:: lon_mod(num_mod), lat_mod(num_mod)
  real(kind=WP), intent(out)  	:: data_mod(num_mod)
  !
  if(lon_reg(1)<0.0_WP .or. lon_reg(num_lon_reg)>360._WP) then
     write(*,*) 'Error in 2D interpolation!'
     write(*,*) 'The regular grid is not in the proper longitude range.'
     call par_ex(partit%MPI_COMM_FESOM, partit%mype)
     stop
  end if

  do n=1,num_mod
     x=lon_mod(n)
     y=lat_mod(n)
     ! find the surrounding rectangular box and get interpolation ratios
     ! 1) north-south direction
     if(y<lat_reg(1)) then
        ind_lat_h=2
        ind_lat_l=1
        y=lat_reg(1)
     elseif(y>lat_reg(num_lat_reg)) then
        ind_lat_h=num_lat_reg
        ind_lat_l=num_lat_reg-1
        y=lat_reg(num_lat_reg)
     else
        do i=2,num_lat_reg
           if(lat_reg(i)>=y) then
              ind_lat_h=i
              ind_lat_l=i-1
              exit
           end if
        end do
     end if
     diff=lat_reg(ind_lat_h)-lat_reg(ind_lat_l)
     rt_lat1=(lat_reg(ind_lat_h)-y)/diff
     rt_lat2=1.0_WP-rt_lat1
     ! 2) east_west direction
     if(x<lon_reg(1)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360._WP-lon_reg(ind_lon_l))
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0_WP-rt_lon1
     elseif(x>lon_reg(num_lon_reg)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360._WP-lon_reg(ind_lon_l))
        rt_lon2=(x-lon_reg(ind_lon_l))/diff
        rt_lon1=1.0_WP-rt_lon2
     else
        do i=2,num_lon_reg
           if(lon_reg(i)>=x) then
              ind_lon_h=i
              ind_lon_l=i-1
              exit
           end if
        end do
        diff=lon_reg(ind_lon_h)-lon_reg(ind_lon_l)
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0_WP-rt_lon1
     end if
     !
     data_ll=data_reg(ind_lon_l,ind_lat_l)
     data_lh=data_reg(ind_lon_l,ind_lat_h)
     data_hl=data_reg(ind_lon_h,ind_lat_l)
     data_hh=data_reg(ind_lon_h,ind_lat_h)
     !
     ! interpolate data
     if(phase_flag==1) then   ! interpolate phase value (0,360)
        if(abs(data_ll-data_hl)>180._WP) then
           if(data_ll<data_hl) then
              data_ll=data_ll+360._WP
           else
              data_hl=data_hl+360._WP
           end if
        end if
        if(abs(data_lh-data_hh)>180._WP) then
           if(data_lh<data_hh) then
              data_lh=data_lh+360._WP
           else
              data_hh=data_hh+360._WP
           end if
        end if
        data_lo = data_ll*rt_lon1 + data_hl*rt_lon2                
        data_up = data_lh*rt_lon1 + data_hh*rt_lon2
        if(abs(data_lo-data_up)>180._WP) then
           if(data_lo<data_up) then
              data_lo=data_lo+360._WP
           else
              data_up=data_up+360._WP
           end if
        end if
        data_mod(n)=data_lo*rt_lat1 + data_up*rt_lat2
        if(data_mod(n)>=360._WP) then
           data_mod(n)=mod(data_mod(n), 360._WP)
        end if
     else   ! other case
        data_mod(n)=(data_ll*rt_lon1 + data_hl*rt_lon2)*rt_lat1 + &
             (data_lh*rt_lon1 + data_hh*rt_lon2)*rt_lat2
     end if
  end do
end subroutine interp_2d_field
!
!---------------------------------------------------------------------------
!
subroutine interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
     lon_reg, lat_reg, lay_reg, data_reg, &
     num_mod_z, num_mod, lon_mod, lat_mod, lay_mod, data_mod, partit, mesh)
  !-------------------------------------------------------------------------------------
  ! This routine does 3d interpolation from a regular grid to specified nodes.
  ! The regular grid is assumed to be global.
  ! The output is data_mod, others are input.
  ! Arguments:
  ! num_lon_reg				number of regular grid points in the longitude direction
  ! num_lat_reg				number of regular grid points in the latitude direction
  ! num_lay_reg				number of regular grid points in the vertical direction
  ! lon_reg(num_lon_reg)		longitude of the regular grid points
  ! lat_reg(num_lat_reg)        	latitude of the regular grid points
  ! lay_reg(num_lay_reg)                depth of the regular grid points
  ! data_reg(:,:,:)  	                data on the regular grid
  ! num_mod      			number of interpolation nodes
  ! lon_mod(num_mod)			longitude of interpolation nodes
  ! lat_mod(num_mod)			latitude of interpolation nodes
  ! lay_mod(num_mod)                    depth of interpolation nodes
  ! data_mod(num_mod)			output data on interpolation nodes
  ! Unit of lon_reg, lat_reg, lon_mod, lat_mod: degree
  ! Unit of lay_reg, lay_mod: m
  ! Order of lon_reg: monotonically increasing in the range [0 360]  
  !                   (range [180,360] is associated with longitude [-180, 0])
  ! Order of lat_reg: monotonically increasing in the range [-90 90]
  ! Order of lay_reg: monotonically decreasing from surface to bottom
  ! Make sure that 'lon_mod' is also in the [0, 360] range.
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------------------------------
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  use o_param, only: WP
  implicit none
  integer             		:: n, i, flag,nz
  integer			:: ind_lat_h, ind_lat_l, ind_lon_h, ind_lon_l
  integer                       :: ind_lay_h, ind_lay_l
  integer, intent(in)         	:: num_lon_reg, num_lat_reg, num_lay_reg
  integer, intent(in)          	:: num_mod_z,num_mod
  real(kind=WP) 		:: x, y, zz, diff 
  real(kind=WP)			:: rt_lat1, rt_lat2, rt_lon1, rt_lon2
  real(kind=WP)                 :: rt_lay1, rt_lay2, v_dup, v_dlo
  real(kind=WP)                 :: data_ll, data_lh, data_hl, data_hh
  real(kind=WP)                 :: v_col(4), z_col(4), H, aux1, aux2
  real(kind=WP)                 :: dz, a, b, c, d
  real(kind=WP), intent(in)	:: lon_reg(num_lon_reg), lat_reg(num_lat_reg)
  real(kind=WP), intent(in)     :: lay_reg(num_lay_reg)
  real(kind=WP), intent(in)	:: data_reg(num_lon_reg, num_lat_reg, num_lay_reg)
  real(kind=WP), intent(in)	:: lon_mod(num_mod), lat_mod(num_mod), lay_mod(num_mod)
  real(kind=WP), intent(out)  	:: data_mod(num_mod_z,num_mod)
  type(t_mesh),   intent(in),    target :: mesh  
  type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  do n=1,num_mod
  !!PS do nz=1,nlevels_nod2D(n)-1
  do nz=ulevels_nod2D(n),nlevels_nod2D(n)-1
!  do nz=1,num_mod_z-1
     x=lon_mod(n)
     y=lat_mod(n)
     zz=lay_mod(nz)

     ! find the surrounding box and get interpolation ratios

     ! 1) north-south direction
     if(y<lat_reg(1)) then
        ind_lat_h=2
        ind_lat_l=1
        y=lat_reg(1)
     elseif(y>lat_reg(num_lat_reg)) then
        ind_lat_h=num_lat_reg
        ind_lat_l=num_lat_reg-1
        y=lat_reg(num_lat_reg)
     else
        do i=2,num_lat_reg
           if(lat_reg(i)>=y) then
              ind_lat_h=i
              ind_lat_l=i-1
              exit
           end if
        end do
     end if
     diff=lat_reg(ind_lat_h)-lat_reg(ind_lat_l)
     rt_lat1=(lat_reg(ind_lat_h)-y)/diff
     rt_lat2=1.0-rt_lat1

     ! 2) east_west direction
     if(x<lon_reg(1)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360._WP-lon_reg(ind_lon_l))
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0_WP-rt_lon1
     elseif(x>lon_reg(num_lon_reg)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360._WP-lon_reg(ind_lon_l))
        rt_lon2=(x-lon_reg(ind_lon_l))/diff
        rt_lon1=1.0_WP-rt_lon2
     else
        do i=2,num_lon_reg
           if(lon_reg(i)>=x) then
              ind_lon_h=i
              ind_lon_l=i-1
              exit
           end if
        end do
        diff=lon_reg(ind_lon_h)-lon_reg(ind_lon_l)
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0_WP-rt_lon1
     end if

     ! 3) up-down direction
     if(zz>lay_reg(1)) then
        ind_lay_h=1
        ind_lay_l=2
        zz=lay_reg(1)
     elseif(zz<lay_reg(num_lay_reg)) then
        ind_lay_h=num_lay_reg-1
        ind_lay_l=num_lay_reg
        zz=lay_reg(num_lay_reg)
     else
        do i=2,num_lay_reg
           if(lay_reg(i)<=zz) then
              ind_lay_h=i-1
              ind_lay_l=i
              exit
           end if
        end do
     end if
     diff=lay_reg(ind_lay_h)-lay_reg(ind_lay_l)
     rt_lay1=(zz-lay_reg(ind_lay_l))/diff
     rt_lay2=1.0_WP-rt_lay1
     data_ll=data_reg(ind_lon_l,ind_lat_l,ind_lay_h)
     data_lh=data_reg(ind_lon_l,ind_lat_h,ind_lay_h)
     data_hl=data_reg(ind_lon_h,ind_lat_l,ind_lay_h)
     data_hh=data_reg(ind_lon_h,ind_lat_h,ind_lay_h)    
     v_dup=(data_ll*rt_lon1 + data_hl*rt_lon2)*rt_lat1 + &
          (data_lh*rt_lon1 + data_hh*rt_lon2)*rt_lat2

     data_ll=data_reg(ind_lon_l,ind_lat_l,ind_lay_l)
     data_lh=data_reg(ind_lon_l,ind_lat_h,ind_lay_l)
     data_hl=data_reg(ind_lon_h,ind_lat_l,ind_lay_l)
     data_hh=data_reg(ind_lon_h,ind_lat_h,ind_lay_l)    
     v_dlo=(data_ll*rt_lon1 + data_hl*rt_lon2)*rt_lat1 + &
          (data_lh*rt_lon1 + data_hh*rt_lon2)*rt_lat2

     data_mod(nz,n)=v_dup*rt_lay1 + v_dlo*rt_lay2

  end do
  end do
end subroutine interp_3d_field
!
!
!_______________________________________________________________________________
! interpolate from node to elements 3d (2, nz, elem)
subroutine interp_e2n_3d(data_e, data_n, mesh, partit, do_overz_in)
    implicit none
    !___INPUT/OUTPUT VARIABLES______________________________________________
    real(kind=WP)   , intent(in   ), dimension(:,:,:) :: data_e
    real(kind=WP)   , intent(inout), dimension(:,:,:) :: data_n    
    type(t_mesh)    , intent(in   ), target           :: mesh
    type(t_partit)  , intent(inout), target           :: partit
    logical         , intent(in)   , optional         :: do_overz_in
    !___LOCAL VARIABLES_____________________________________________________
    integer                                           :: node, elem, nz, nl1, ul1, k
    real(kind=WP)                                     :: vx, vy
    logical                                           :: do_overz
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"     
    !
    !___________________________________________________________________________
    do_overz = .True.
    if (present(do_overz_in)) do_overz = do_overz_in
    
    !
    !___________________________________________________________________________
    ! Do interpolation over vertical depth dimension consider bottom topography 
    if (do_overz) then 
        do node = 1, myDim_nod2D
            nl1 = nlevels_nod2D(node)-1
            ul1 = ulevels_nod2D(node)
            do nz = ul1, nl1
                vx = 0.0_WP
                vy = 0.0_WP
                do k = 1, nod_in_elem2D_num(node)
                    elem = nod_in_elem2D(k, node)
                    if( nz.LE.(nlevels(elem)-1) .and. nz.GE.(ulevels(elem))) then
                        vx = vx + data_e(1, nz, elem) * elem_area(elem)
                        vy = vy + data_e(2, nz, elem) * elem_area(elem)
                    endif
                end do
                data_n(1, nz, node)=vx/3.0_WP/areasvol(nz, node)
                data_n(2, nz, node)=vy/3.0_WP/areasvol(nz, node)
            end do
        end do
    
    !
    !___________________________________________________________________________
    ! Do interpolation along spectral domain no bottom topography to consider
    else
        ul1=1
        nl1=size(data_e,2)
        do node = 1, myDim_nod2D
            do nz = ul1+1, nl1-1
                vx = 0.0_WP
                vy = 0.0_WP
                do k = 1, nod_in_elem2D_num(node)
                    elem = nod_in_elem2D(k, node)
                    vx = vx + data_e(1, nz, elem) * elem_area(elem)
                    vy = vy + data_e(2, nz, elem) * elem_area(elem)
                end do
                data_n(1, nz, node)=vx/3.0_WP/areasvol(1, node)
                data_n(2, nz, node)=vy/3.0_WP/areasvol(1, node)
            end do
        end do
    end if 
    
end subroutine interp_e2n_3d


!
!
!_______________________________________________________________________________
! interpolate from elem to nodes 2d (2, elem), (nz, elem), (nfbin, elem) --> (.., node)
subroutine interp_e2n_2d(data_e, data_n, mesh, partit, do_overz_in)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_param, only: WP
    implicit none
    !___INPUT/OUTPUT VARIABLES______________________________________________
    real(kind=WP)   , intent(in   ), dimension(:,:) :: data_e
    real(kind=WP)   , intent(inout), dimension(:,:) :: data_n    
    type(t_mesh)    , intent(in   ), target         :: mesh
    type(t_partit)  , intent(inout), target         :: partit
    logical         , intent(in)   , optional       :: do_overz_in
    !___LOCAL VARIABLES_____________________________________________________
    integer                                         :: node, elem, nz, nl1, ul1, k
    real(kind=WP)                                   :: vx, vy
    logical                                         :: do_overz
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"     
    !
    !___________________________________________________________________________
    do_overz = .True.
    if (present(do_overz_in)) do_overz = do_overz_in
    
    !___________________________________________________________________________
    ! Do interpolation for 2d vector data e.g 2d gradient (2, elem)
    if (size(data_e, 1) == 2) then
        do node = 1, myDim_nod2D
            vx = 0.0_WP
            vy = 0.0_WP
            do k = 1, nod_in_elem2D_num(node)
                elem = nod_in_elem2D(k, node)
                vx = vx + data_e(1, elem) * elem_area(elem)
                vy = vy + data_e(2, elem) * elem_area(elem)
            end do    
            data_n(1, node)=vx/3.0_WP/areasvol(1, node)
            data_n(2, node)=vy/3.0_WP/areasvol(1, node)
        end do
    
    !___________________________________________________________________________
    ! Do interpolation over vertical depth dimension consider bottom topography (nz,elem)
    else if (do_overz) then 
        do node = 1, myDim_nod2D
            nl1 = nlevels_nod2D(node)-1
            ul1 = ulevels_nod2D(node)
            do nz = ul1, nl1
                vx = 0.0_WP
                do k = 1, nod_in_elem2D_num(node)
                    elem = nod_in_elem2D(k, node)
                    if( nz.LE.(nlevels(elem)-1) .and. nz.GE.(ulevels(elem))) then
                        vx = vx + data_e(nz, elem) * elem_area(elem)
                    endif
                end do
                data_n(nz, node)=vx/3.0_WP/areasvol(nz, node) 
            end do
        end do
    
    !___________________________________________________________________________
    ! Do interpolation along spectral domain no bottom topography to consider (nfbin, elem)
    else
        ul1=1
        nl1=size(data_e,1)
        do node = 1, myDim_nod2D
            do nz = ul1+1, nl1-1
                vx = 0.0_WP
                do k = 1, nod_in_elem2D_num(node)
                    elem = nod_in_elem2D(k, node)
                    vx = vx + data_e(nz, elem) * elem_area(elem)
                end do
                data_n(nz, node)=vx/3.0_WP/areasvol(1, node)
            end do
        end do
    end if 
    
end subroutine interp_e2n_2d


!
!
!_______________________________________________________________________________
! interpolate from elem to nodes 1d (elem) --> (node)
subroutine interp_e2n_1d(data_e, data_n, mesh, partit)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_param, only: WP
    implicit none
    !___INPUT/OUTPUT VARIABLES______________________________________________
    real(kind=WP)   , intent(in   ), dimension(:)   :: data_e
    real(kind=WP)   , intent(inout), dimension(:)   :: data_n    
    type(t_mesh)    , intent(in   ), target         :: mesh
    type(t_partit)  , intent(inout), target         :: partit
    !___LOCAL VARIABLES_____________________________________________________
    integer                                         :: node, elem, nz, nl1, ul1, k
    real(kind=WP)                                   :: vx
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"     
    !
    !___________________________________________________________________________
    ! Do interpolation for 1d vector data (elem) --> (node)
    do node = 1, myDim_nod2D
        vx = 0.0_WP
        do k = 1, nod_in_elem2D_num(node)
            elem = nod_in_elem2D(k, node)
            vx = vx + data_e(elem) * elem_area(elem)
            data_n(node)=vx/3.0_WP/areasvol(1, node)
        end do
    end do    
end subroutine interp_e2n_1d

end module g_interp