!========================================================================================
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
    
     heat_flux   = -net_heat_flux 
     water_flux  = -fresh_wa_flux
  ! ==================
  ! momentum flux:
  ! ==================
  if(ice_v_n) then 
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
  else               ! ice velocity, .... at elements 
  do n=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     aux1=sum(a_ice(elnodes))/3.0_WP
     if(aux1<0.001) then
        stress_iceoce_x(elem)=0.0_WP
	stress_iceoce_y(elem)=0.0_WP
     else	
        aux=sqrt((u_ice(elem)-u_w(elem))**2+(v_ice(elem)-v_w(elem))**2)* &
                density_0*Cd_oce_ice
        stress_iceoce_x(elem) = aux * (u_ice(elem)-u_w(elem))
        stress_iceoce_y(elem) = aux * (v_ice(elem)-v_w(elem))
     end if
     stress_surf(1,elem)=stress_iceoce_x(elem)*aux1 + &
                             sum(stress_atmoce_x(elnodes))*(1.0_WP-aux1)/3.0_WP
     stress_surf(2,elem)=stress_iceoce_y(elem)*aux1 + &
                             sum(stress_atmoce_y(elnodes))*(1.0_WP-aux1)/3.0_WP
  end do
  end if
  call exchange_nod(heat_flux)
  call exchange_nod(water_flux)
end subroutine ice2ocean
!=======================================================================================
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

  do n=1, myDim_nod2d+eDim_nod2d         
if (ice_update) then
     T_oc_array(n)=tr_arr(1,n,1)
     S_oc_array(n)=tr_arr(1,n,2)  
     elevation(n)= eta_n(n)
else
     T_oc_array(n)=(T_oc_array(n)*real(ice_steps_since_upd)+tr_arr(1,n,1))/real(ice_steps_since_upd+1)
     S_oc_array(n)=(S_oc_array(n)*real(ice_steps_since_upd)+tr_arr(1,n,2))/real(ice_steps_since_upd+1)
     elevation(n)=(elevation(n)*real(ice_steps_since_upd)+eta_n(n))/real(ice_steps_since_upd+1)
write(*,*) 'debugging mode, remove this stop statement in ocean2ice' !DS
call par_ex
stop
endif
  end do
  if(ice_v_n) then   ! ice velocity at nodes
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
write(*,*) 'debugging mode, remove this stop statement in ocean2ice' !DS
call par_ex
stop
     u_w(n)=(u_w(n)*real(ice_steps_since_upd)+uw)/real(ice_steps_since_upd+1)
     v_w(n)=(v_w(n)*real(ice_steps_since_upd)+vw)/real(ice_steps_since_upd+1)
endif
     enddo
     call exchange_nod(u_w)
     call exchange_nod(v_w)
  else  ! ice velocity at elements
!not really tested
if (ice_update) then
     do elem=1,myDim_elem2D
        u_w(elem)=UV(1,1,elem)
	v_w(elem)=UV(2,1,elem)
     end do
else
write(*,*) '??? u_w(n) shall be changed to u_w(elem), same for u_w(n) ???'
call par_ex
stop
     do elem=1,myDim_elem2D
      u_w(n)=(u_w(n)*real(ice_steps_since_upd)+UV(1,1,elem))/real(ice_steps_since_upd+1)
      v_w(n)=(v_w(n)*real(ice_steps_since_upd)+UV(2,1,elem))/real(ice_steps_since_upd+1)
     end do     
endif	
  end if

end subroutine ocean2ice
!=========================================================================================================
