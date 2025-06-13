!=======================================================================
!
! This submodule exchanges variables between icepack and fesom2
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!
!=======================================================================

      submodule (icedrv_main) icedrv_transfer

      contains

        module subroutine fesom_to_icepack(ice, mesh)

          ! fesom modules
          use g_forcing_arrays, only: Tair, shum, u_wind,   v_wind,       & ! Atmospheric forcing fields
               shortwave,  longwave, prec_rain, prec_snow,  press_air
          use g_sbf, only: l_mslp                     
          use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted, &
               icepack_query_parameters, icepack_sea_freezing_temperature
          use g_comm_auto, only: exchange_nod
          use icedrv_system, only: icedrv_system_abort
          use mod_mesh
          use mod_tracer
          use mod_ice
          use g_clock
          ! icepack modules
          use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted, &
	       icepack_query_parameters, icepack_sea_freezing_temperature
          use icedrv_system, only: icedrv_system_abort
 
          implicit none

          character(len=*), parameter :: subname='(fesom_to_icepack)'

          logical(kind=log_kind) :: calc_strair
          integer(kind=int_kind) :: i, n,  k,  elem
          real(kind=dbl_kind)    :: tx, ty, tvol, rhowat, aux, cprho

          type(t_mesh), target, intent(in)     :: mesh
          type(t_ice), target, intent(inout)   :: ice
          real(kind=WP), dimension(:), pointer :: u_ice, v_ice, S_oc_array, T_oc_array, & 
               u_w, v_w, stress_atmice_x, stress_atmice_y

          ! fesom parameter
          real (kind=dbl_kind), parameter :: &
               frcvdr = 0.28_dbl_kind,    & ! frac of incoming sw in vis direct band
               frcvdf = 0.24_dbl_kind,    & ! frac of incoming sw in vis diffuse band
               frcidr = 0.31_dbl_kind,    & ! frac of incoming sw in near IR direct band
               frcidf = 0.17_dbl_kind,    & ! frac of incoming sw in near IR diffuse band
               R_dry  = 287.05_dbl_kind,  & ! specific gas constant for dry air (J/K/kg)
               R_vap  = 461.495_dbl_kind, & ! specific gas constant for water vapo (J/K/kg)
               ex     = 0.286_dbl_kind
          ! PS rhowat = 1025.0_dbl_kind,  & ! Water density
          ! PS cc     = rhowat*4190.0_dbl_kind, & ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)

#include "associate_mesh.h"

          ! make sure we use the same water density to scale precip and snow in 
          ! fesom and icepack
          call icepack_query_parameters(rhow_out=rhowat) 

          u_ice           => ice%uice(:)
          v_ice           => ice%vice(:)
          S_oc_array      => ice%srfoce_salt(:)
          T_oc_array      => ice%srfoce_temp(:)
          u_w             => ice%srfoce_u(:)
          v_w             => ice%srfoce_v(:)
          stress_atmice_x => ice%stress_atmice_x(:)
          stress_atmice_y => ice%stress_atmice_y(:)
          ! Ice 
    
          uvel(:)   = u_ice(:)
          vvel(:)   = v_ice(:)
    
          ! Atmosphere 
    
          T_air(:)  = Tair(:) + 273.15_dbl_kind
          Qa(:)     = shum(:)
          uatm(:)   = u_wind(:)
          vatm(:)   = v_wind(:)
          fsw(:)    = shortwave(:)
          flw(:)    = longwave(:)
          frain(:)  = prec_rain(:) * rhowat
          fsnow(:)  = prec_snow(:) * rhowat
          !if (mype==7) write(nu_diag,*) 'icedrv.fesom_to_icepack: Max/Min air temperature = ', maxval(T_air), minval(T_air)
          !if (mype==7) write(nu_diag,*) 'icedrv.fesom_to_icepack: Max/Min humidity = ', maxval(Qa), minval(Qa)
          !if (mype==7) write(nu_diag,*) 'icedrv.fesom_to_icepack: Max/Min longwave = ', maxval(flw), minval(flw)

          wind(:)   = sqrt(uatm(:)**2 + vatm(:)**2)

          if ( l_mslp ) then
             potT(:) = T_air(:)*(press_air(:)/100000.0_dbl_kind)**ex             
             rhoa(:) = press_air(:) / (R_dry * T_air(:) * (c1 + ((R_vap/R_dry) * Qa) ))
          else 
             ! The option below is used in FESOM2
             potT(:) = T_air(:)
             rhoa(:) = 1.3_dbl_kind
          endif
    
          ! Ocean 
    
          sst(:)     = T_oc_array(:)
          sstdat(:)  = T_oc_array(:)
          sss(:)     = S_oc_array(:)
          uocn(:)    = u_w(:)
          vocn(:)    = v_w(:)

          ! divide shortwave into spectral bands
          swvdr      = fsw*frcvdr        ! visible direct
          swvdf      = fsw*frcvdf        ! visible diffuse
          swidr      = fsw*frcidr        ! near IR direct
          swidf      = fsw*frcidf        ! near IR diffuse

          call icepack_query_parameters(calc_strair_out=calc_strair, cprho_out=cprho)
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__,line= __LINE__)
    
          if (.not. calc_strair) then
             strax(:) = stress_atmice_x(:)
             stray(:) = stress_atmice_y(:)
          endif
    
          do i = 1, nx
             ! ocean - ice stress
             aux = sqrt((uvel(i)-uocn(i))**2+(vvel(i)-vocn(i))**2)*rhowat*ice%cd_oce_ice
             strocnxT(i) = aux*(uvel(i) - uocn(i))
             strocnyT(i) = aux*(vvel(i) - vocn(i))
             ! freezing - melting potential
             Tf(i)   = icepack_sea_freezing_temperature(sss(i))
             frzmlt(i) = min((Tf(i)-sst(i)) * cprho * hmix(i) / dt, 1000.0_dbl_kind)
             !if (geo_coord_nod2D(2, i) > 89.5 * 3.14159/180.) &
             !     write(*,*) 'pole index =', i, geo_coord_nod2D(1, i), geo_coord_nod2D(2, i), mype
          enddo

          ! Compute convergence and shear on the nodes

          do i = 1, nx_nh
             tvol = c0
             tx   = c0
             ty   = c0
             do k = 1, nod_in_elem2D_num(i)
                elem = nod_in_elem2D(k,i)
                tvol = tvol + elem_area(elem)
                tx = tx + rdg_conv_elem(elem)  * elem_area(elem)
                ty = ty + rdg_shear_elem(elem) * elem_area(elem)
             enddo
             rdg_conv(i)  = tx / tvol
             rdg_shear(i) = ty / tvol
          enddo
          
          call exchange_nod(rdg_conv, rdg_shear, p_partit)

          ! Clock variables
    
          days_per_year = ndpyr
          daymo         = num_day_in_month(fleapyear,:)
          if (fleapyear==1) then
             daycal     = daycal366
          else
             daycal     = daycal365
          end if
          istep1        = mstep
          time          = mstep*dt
          mday          = day_in_month
          month_i       = month
          nyr           = yearnew
          sec           = timenew
          yday          = real(ndpyr, kind=dbl_kind)
          dayyr         = real(days_per_year, kind=dbl_kind)
          secday        = real(sec, kind=dbl_kind)
          calendar_type = 'Gregorian'
          dt_dyn        = dt/real(ndtd,kind=dbl_kind) ! dynamics et al timestep
          
        end subroutine fesom_to_icepack

        !=======================================================================

        module subroutine fesom_to_icepack_para()

          ! fesom modules
          use g_forcing_param, only: ncar_bulk_z_wind, ncar_bulk_z_tair, ncar_bulk_z_shum
          use g_config, only: dt
          use o_param, only: mstep
          ! icepack modules
          use icedrv_system, only: icedrv_system_abort

          character(len=*), parameter :: subname='(fesom_to_icepack_para)'

          if ( ncar_bulk_z_tair .ne. ncar_bulk_z_shum ) then
             call icedrv_system_abort(string=subname, file=__FILE__,line= __LINE__)
          else
             zlvl_s    = ncar_bulk_z_tair ! should be equal to ncar_bulk_z_shum
             zlvl_v    = ncar_bulk_z_wind
          endif

        end subroutine fesom_to_icepack_para

        !=======================================================================
        
        module subroutine icepack_to_fesom( &
             nx_in, aice_out,  vice_out,  vsno_out,  &
             fhocn_tot_out, fresh_tot_out,    &
             strocnxT_out,  strocnyT_out,     &
             dhs_dt_out,    dhi_dt_out,       &
             fsalt_out,     evap_ocn_out,     &
             evap_out)
          
          implicit none
          
          integer (kind=int_kind), intent(in) :: &
               nx_in      ! block dimensions
          
          real (kind=dbl_kind), dimension(nx_in), intent(out), optional :: &
               aice_out, &  
               vice_out, &
               vsno_out, &
               fhocn_tot_out, &
               fresh_tot_out, &
               strocnxT_out,  &
               strocnyT_out,  &
               fsalt_out,     &
               dhs_dt_out,    &
               dhi_dt_out,    &
               evap_ocn_out,  &
               evap_out
          
          character(len=*),parameter :: subname='(icepack_to_fesom)'   
          
          if (present(aice_out)              ) aice_out         = aice
          if (present(vice_out)              ) vice_out         = vice
          if (present(vsno_out)              ) vsno_out         = vsno
          if (present(fresh_tot_out)         ) fresh_tot_out    = fresh_tot
          if (present(fhocn_tot_out)         ) fhocn_tot_out    = fhocn_tot
          if (present(strocnxT_out)          ) strocnxT_out     = strocnxT
          if (present(strocnyT_out)          ) strocnyT_out     = strocnyT
          if (present(dhi_dt_out)            ) dhi_dt_out       = dhi_dt
          if (present(dhs_dt_out)            ) dhs_dt_out       = dhs_dt
          if (present(fsalt_out)             ) fsalt_out        = fsalt
          if (present(evap_ocn_out)          ) evap_ocn_out     = evap_ocn
          if (present(evap_out)              ) evap_out         = evap
          
        end subroutine icepack_to_fesom

        !=======================================================================

        module subroutine icepack_to_fesom_single_point(nx_in, strength_out)
          
          implicit none

          integer (kind=int_kind), intent(in) :: nx_in      ! surface node or element
          real (kind=dbl_kind), intent(out), optional :: strength_out

          character(len=*),parameter :: subname='(icepack_to_fesom_single_point)'

          if (present(strength_out)          ) strength_out     = strength(nx_in)

        end subroutine icepack_to_fesom_single_point

        !=======================================================================

      end submodule icedrv_transfer
