!
!----------------------------------------------------------------------------
!
subroutine read_fesom(ncid, m, id, nxy, nz, res)
  use g_config
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  implicit none
  #include "netcdf.inc"

  integer, intent(in)                           :: ncid, m, nxy, nz
  character(15), intent(in)                     :: id
  real(kind=8), dimension(nz, nxy), intent(out) :: res

  integer                                       :: i, s, varid, istart(3), icount(3)
  INTEGER, DIMENSION(50)                        :: stat    !  status array
  write(*,*) 'reading netcdf input; month, xydim, yzdim=', m, nxy, nz
  s=1
  istart= (/1, 1, m/)
  icount=(/nz, nxy, 1/)
  stat(s)=nf_inq_varid(ncid, trim(id), varid); s=s+1
  stat(s)=nf_get_vara_double(ncid,varid,istart,icount,res); s=s+1
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in reading variable, i, stat=', i, stat(i)
  end do  
end subroutine read_fesom
!
!----------------------------------------------------------------------------
!
subroutine read_fesom_lev(ncid, m, id, nxy, lev, res)
  use g_config
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  implicit none
  #include "netcdf.inc"

  integer, intent(in)                           :: ncid, m, nxy, lev
  character(15), intent(in)                     :: id
  real(kind=8), dimension(nxy), intent(out)     :: res

  integer                                       :: i, s, varid, istart(3), icount(3)
  INTEGER, DIMENSION(50)                        :: stat    !  status array

  write(*,*) 'reading netcdf input; month, xydim=', m, nxy, ' at level ', lev
  s=1
  istart= (/lev, 1, m/)
  icount= (/1,  nxy, 1/)
  stat(s)=nf_inq_varid(ncid, trim(id), varid); s=s+1
  stat(s)=nf_get_vara_double(ncid,varid,istart,icount,res); s=s+1
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in reading variable, i, stat=', i, stat(i), &
        ' var=', trim(id), ', lev=', lev
  end do  
end subroutine read_fesom_lev
