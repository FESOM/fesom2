!=======================================================================
!
! This submodule exchanges variables between icepack and fesom2
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!
!=======================================================================

      submodule (icedrv_main) icedrv_transfer

      contains

      module subroutine fesom_to_icepack()

          implicit none

          character(len=*), parameter :: subname='(fesom_to_icepack)'
    
          real (kind=dbl_kind), parameter :: &
             frcvdr = 0.28_dbl_kind, & ! frac of incoming sw in vis direct band
             frcvdf = 0.24_dbl_kind, & ! frac of incoming sw in vis diffuse band
             frcidr = 0.31_dbl_kind, & ! frac of incoming sw in near IR direct band
             frcidf = 0.17_dbl_kind, & ! frac of incoming sw in near IR diffuse band
             ex     = 0.286_dbl_kind

          use g_forcing_arrays, only: Tair, shum, u_wind,   v_wind,       & ! Atmospheric forcing fields
                                      shortwave,  longwave, prec_rain,    &
                                      prec_snow,  press_air
          use g_forcing_param,  only: ncar_bulk_z_wind, ncar_bulk_z_tair, &
                                      ncar_bulk_z_shum
          use g_sbf,            only: l_mslp                     
          use i_arrays,         only: S_oc_array, T_oc_array,             & ! Ocean and sea ice field
                                      u_w,        v_w,                    &
                                      u_ice,      v_ice
    
          ! Ice 
    
          uvel(:)  = u_ice(:)
          vvel(:)  = v_ice(:)
    
          ! Atmosphere 
    
          T_air(:)  = Tair(:) + 273.15_dbl_kind
          Qa(:)     = shum(:)
          uatm(:)   = u_wind(:)
          vatm(:)   = v_wind(:)
          fsw(:)    = shortwave(:)
          flw(:)    = longwave(:)
          frain(:)  = prec_rain(:) * 1000.0_dbl_kind
          fsnow(:)  = prec_snow(:) * 1000.0_dbl_kind

          zlvl_t    = ncar_bulk_z_tair
          zlvl_q    = ncar_bulk_z_shum
          zlvl_v    = ncar_bulk_z_wind

          if ( l_mslp ) then
             potT(:) = T_air(:)*(press_air(:)/100000.0_dbl_kind)^ex             
          else 
             potT(:) = T_air(:)
          endif
    
          ! Ocean 
    
          sst(:)    = T_oc_array(:)
          sstdat(:) = T_oc_array(:)
          sss(:)    = S_oc_array(:)
          uocn(:)   = u_w(:)
          vocn(:)   = v_w(:)

          ! divide shortwave into spectral bands
          swvdr = fsw*frcvdr        ! visible direct
          swvdf = fsw*frcvdf        ! visible diffuse
          swidr = fsw*frcidr        ! near IR direct
          swidf = fsw*frcidf        ! near IR diffuse

      end subroutine fesom_to_icepack

!=======================================================================


      module subroutine icepack_to_fesom()

          implicit none


!=======================================================================

      end subroutine icepack_to_fesom

!=======================================================================

      end submodule icedrv_transfer 
