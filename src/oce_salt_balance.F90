subroutine check_imb_freshwater
  use o_param
  use o_mesh 
  use o_elements
  use o_array
  use g_config
  use g_forcing_arrays
  use g_parfe
  implicit none

  integer      :: row
  real(kind=8) :: flux, corr

  flux=0.0
  do row=1,myDim_nod2D
#ifdef use_cavity
     if(cavity_flag_nod2d(row)==1) cycle   
#endif 
     flux=flux+(evaporation(row)+prec_rain(row)+ &
          prec_snow(row)+runoff(row))*cluster_area_2d(row)
  end do

  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  corr=0.0
  call MPI_AllREDUCE( flux, corr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  corr=corr/ocean_area
  water_flux=water_flux+corr  ! the + sign should be used here

end subroutine check_imb_freshwater
!
!-------------------------------------------------------------------------
!
subroutine check_imb_salt_flux
  use o_param
  use o_mesh 
  use o_elements
  use o_array
  use g_config
  use g_parfe
  implicit none

  integer      :: row
  real(kind=8) :: flux_rel, flux_vir, corr

  flux_rel=0.0
  flux_vir=0.0

  do row=1,myDim_nod2D
#ifdef use_cavity
     if(cavity_flag_nod2d(row)==1) cycle   
#endif 
     flux_rel=flux_rel+relax_salt(row)*cluster_area_2d(row)
#ifndef use_fullfreesurf
     flux_vir=flux_vir+virtual_salt(row)*cluster_area_2d(row)
#endif
  end do

  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  corr=0.0
  call MPI_AllREDUCE(flux_rel, corr, 1, MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  corr=corr/ocean_area
  relax_salt=relax_salt-corr

#ifndef use_fullfreesurf  
  corr=0.0
  call MPI_AllREDUCE(flux_vir, corr, 1, MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  corr=corr/ocean_area
  virtual_salt=virtual_salt-corr
#endif
  
end subroutine
