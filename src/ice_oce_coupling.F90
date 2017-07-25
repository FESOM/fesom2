!
!======================================================================================
!
subroutine ice2ocean
  ! transmits the relevant fields from the ice to the ocean model
  !
  use o_PARAM
  use o_ARRAYS
  use o_MESH
  use i_ARRAYS
  use g_PARSUP
  use i_PARAM
  USE g_CONFIG
  use g_comm_auto
  implicit none
  
  integer                 :: n, elem, elnodes(3),n1
  real(kind=WP)           :: aux, aux1
  ! ==================
  ! heat and freshwater
  ! ==================
     heat_flux_old = heat_flux !PS
     water_flux_old = water_flux !PS
     
     heat_flux   = -net_heat_flux 
     water_flux  = -fresh_wa_flux

     call exchange_nod_begin(heat_flux, water_flux)
  ! ==================
  ! momentum flux:
  ! ==================
  do n=1,myDim_nod2D+eDim_nod2D   
     if(a_ice(n)>0.001) then
       aux=sqrt((u_ice(n)-u_w(n))**2+(v_ice(n)-v_w(n))**2)*density_0*Cd_oce_ice
       stress_iceoce_x(n) = aux * (u_ice(n)-u_w(n))
       stress_iceoce_y(n) = aux * (v_ice(n)-v_w(n))
     else
       stress_iceoce_x(n)=0.0
       stress_iceoce_y(n)=0.0
     end if
  end do
  DO elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     stress_surf(1,elem)=sum(stress_iceoce_x(elnodes)*a_ice(elnodes) + &
                             stress_atmoce_x(elnodes)*(1.0_WP-a_ice(elnodes)))/3.0_WP
     stress_surf(2,elem)=sum(stress_iceoce_y(elnodes)*a_ice(elnodes) + &
                             stress_atmoce_y(elnodes)*(1.0_WP-a_ice(elnodes)))/3.0_WP
  END DO
  call exchange_nod_end()
  if (use_sw_pene) call cal_shortwave_rad
end subroutine ice2ocean
!
!======================================================================================
!
subroutine ocean2ice
  
  ! transmits the relevant fields from the ocean to the ice model

  use o_PARAM
  use o_ARRAYS
  use i_ARRAYS
  use o_MESH
  use g_PARSUP
  USE g_CONFIG
  use g_comm_auto
  implicit none

  integer :: n, elem, k
  real*8 :: uw,vw
  ! the arrays in the ice model are renamed
     
if (ice_update) then
  do n=1, myDim_nod2d+eDim_nod2d    
     T_oc_array(n)=tr_arr(1,n,1)
     S_oc_array(n)=tr_arr(1,n,2)  
  end do
  if ( .not. use_ALE ) then
     elevation(:)= eta_n(:)
  else
     elevation(:)= hbar(:)
  endif
else
  do n=1, myDim_nod2d+eDim_nod2d    
     T_oc_array(n)=(T_oc_array(n)*real(ice_steps_since_upd)+tr_arr(1,n,1))/real(ice_steps_since_upd+1)
     S_oc_array(n)=(S_oc_array(n)*real(ice_steps_since_upd)+tr_arr(1,n,2))/real(ice_steps_since_upd+1)
!NR !PS      elevation(n)=(elevation(n)*real(ice_steps_since_upd)+eta_n(n))/real(ice_steps_since_upd+1)
!NR     elevation(n)=(elevation(n)*real(ice_steps_since_upd)+hbar(n))/real(ice_steps_since_upd+1) !PS
  end do
  if ( .not. use_ALE ) then
     elevation(:)= (elevation(:)*real(ice_steps_since_upd)+eta_n(:))/real(ice_steps_since_upd+1)
  else
     elevation(:)= (elevation(:)*real(ice_steps_since_upd)+hbar(:))/real(ice_steps_since_upd+1)
  endif
endif
     do n=1, myDim_nod2d  
       uw=0.0_WP
       vw=0.0_WP
       DO k=1, nod_in_elem2D_num(n)
          elem=nod_in_elem2D(k,n)          
          uw = uw+ UV(1,1,elem)*elem_area(elem)
          vw = vw+ UV(2,1,elem)*elem_area(elem)
       END DO
       uw = uw/area(1,n)/3.0_WP	  
       vw = vw/area(1,n)/3.0_WP	 
     
if (ice_update) then
     u_w(n)=uw
     v_w(n)=vw
else
     u_w(n)=(u_w(n)*real(ice_steps_since_upd)+uw)/real(ice_steps_since_upd+1)
     v_w(n)=(v_w(n)*real(ice_steps_since_upd)+vw)/real(ice_steps_since_upd+1)
endif
     enddo
     call exchange_nod(u_w, v_w)
end subroutine ocean2ice
!
!======================================================================================
!
subroutine oce_fluxes
  use o_ARRAYS
  use i_ARRAYS
  use g_comm_auto
  use g_forcing_param, only: use_virt_salt
  use g_forcing_arrays
  use g_PARSUP,        only: myDim_nod2D, eDim_nod2D
  implicit none
  integer                   :: n, elem, elnodes(3),n1
  real(kind=WP)             :: rsss, net
  real(kind=8), allocatable :: flux(:)

  allocate(flux(myDim_nod2D+eDim_nod2D))
  ! ==================
  ! heat and freshwater
  ! ==================
  heat_flux_old  = heat_flux !PS
  water_flux_old = water_flux !PS
     
  heat_flux   = -net_heat_flux 
  water_flux  = -fresh_wa_flux

  call exchange_nod(heat_flux, water_flux) ! do we really need it?

  ! virtual salt flux
  if (use_virt_salt) then
     rsss=ref_sss
     do n=1, myDim_nod2D+eDim_nod2D
        if (ref_sss_local) rsss = tr_arr(1,n,2)
        virtual_salt(n)=rsss*water_flux(n) 
     end do
  end if

  ! SSS restoring to climatology
  do n=1, myDim_nod2D+eDim_nod2D
     relax_salt(n)=surf_relax_S*(Ssurf(n)-tr_arr(1,n,2))
  end do

  ! enforce the total freshwater/salt flux be zero

  ! 1. water flux ! if (.not. use_virt_salt) can be used!
  ! we conserve only the fluxes from the database plus evaporation.
  ! the rest (ocean/ice transformation etc. will follow from the conservation of volume)
  flux=evaporation+prec_rain+ prec_snow+runoff  
  call comp_net_imbalance(flux, net)
  water_flux=water_flux+net ! the + sign should be used here
    
  ! 2. virtual salt flux
  if (use_virt_salt) then ! virtual_salt array is not allocated otherwise !
     call comp_net_imbalance(virtual_salt, net)
     virtual_salt=virtual_salt-net
  end if

  ! 3. restoring to SSS climatology
  call comp_net_imbalance(relax_salt, net)
  relax_salt=relax_salt-net

  deallocate(flux)
end subroutine oce_fluxes
!
!======================================================================================
!
subroutine comp_net_imbalance(data, net)
  use o_MESH
  use g_PARSUP
  use g_comm_auto

  IMPLICIT NONE
  real(KIND=WP), dimension(:), intent(in) :: data
  real(kind=WP), intent(out)              :: net

  integer      :: row
  real(kind=8) :: flux, corr

  flux=0.0
  do row=1,myDim_nod2D
     flux=flux+data(row)
  end do

  corr=0.0
  call MPI_AllREDUCE(flux, corr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  net=corr/ocean_area
end subroutine comp_net_imbalance
!
!======================================================================================
!
