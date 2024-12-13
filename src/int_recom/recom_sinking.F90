module diff_ver_recom_expl_interface
  interface
    subroutine diff_ver_recom_expl(tr_num, tracer, partit, mesh) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        integer       , intent(in)   , target :: tr_num
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine
  end interface
end module
module ver_sinking_recom_interface
  interface
    subroutine ver_sinking_recom(tr_num, tracer, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        integer       , intent(in)   , target :: tr_num
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine
  end interface
end module
module ver_sinking_recom_benthos_interface
  interface
    subroutine ver_sinking_recom_benthos(tr_num, tracer, partit, mesh)
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        integer       , intent(in)   , target :: tr_num
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine
  end interface
end module
!===============================================================================
! YY: sinking of second detritus adapted from Ozgur's code
! but not using recom_det_tracer_id, since
! second detritus has a different sinking speed than the first
! define recom_det2_tracer_id to make it consistent???
!===============================================================================
subroutine ver_sinking_recom_benthos(tr_num, tracers, partit, mesh)

    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_TRACER

    use recom_declarations
    use recom_locvar
    use recom_glovar
    use recom_config
    use recom_ciso
    ! use ver_sinking_recom_benthos_interface

    use g_support
    use g_clock
    use o_PARAM
    use g_config
    use o_param           ! num_tracers
    use o_arrays
    use g_forcing_arrays  ! press_air
    use g_comm_auto
    implicit none

    integer       , intent(in)   , target :: tr_num
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh

    integer                   :: elem,k
    integer                   :: nl1,ul1,nz,n,nzmin, nzmax, net
    real(kind=WP)             :: Vben(mesh%nl),  aux(mesh%nl-1),  flux(mesh%nl), add_benthos_2d(partit%myDim_nod2D)
    integer                   :: nlevels_nod2D_minimum
    real(kind=WP)             :: tv
    real(kind=WP), dimension(:,:), pointer :: trarr

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    trarr=>tracers%data(tr_num)%values(:,:)

   do n=1, myDim_nod2D ! needs exchange_nod in the end
        nl1=nlevels_nod2D(n)-1
        ul1=ulevels_nod2D(n)

        aux=0._WP
        Vben=0._WP
        add_benthos_2d=0._WP

! Calculate sinking velociy for vertical sinking case
! ******************************************************
          if (any(recom_det_tracer_id == tracers%data(tr_num)%ID)) Vben = Vdet
          if (any(recom_phy_tracer_id == tracers%data(tr_num)%ID)) Vben = VPhy
          if (any(recom_dia_tracer_id == tracers%data(tr_num)%ID)) Vben = VDia
        if (allow_var_sinking) then
          Vben = Vdet_a * abs(zbar_3d_n(:,n)) + Vben
        end if

! Constant vertical sinking for the second detritus class
! *******************************************************

#if defined(__3Zoo2Det)
          if(tracers%data(tr_num)%ID==1025 .or. &  !idetz2n
             tracers%data(tr_num)%ID==1026 .or. &  !idetz2c
             tracers%data(tr_num)%ID==1027 .or. &  !idetz2si
             tracers%data(tr_num)%ID==1028 ) then  !idetz2calc
             Vben = VDet_zoo2
          endif
#endif

        Vben= Vben/SecondsPerDay ! conversion [m/d] --> [m/s] (vertical velocity, note that it is positive here)

        k=nod_in_elem2D_num(n)

        !! * Screening minimum depth in neigbouring nodes around node n*
        nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

        do nz=nlevels_nod2D_minimum, nl1
           tv = trarr(nz,n)*Vben(nz)
           aux(nz)= - tv*(area(nz,n)-area(nz+1,n))
        end do

        do nz=ul1,nl1
           str_bf(nz,n) = str_bf(nz,n) + (aux(nz))*dt/area(nz,n)/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
           add_benthos_2d(n) = add_benthos_2d(n) - (aux(nz))*dt    !!!!!!!!CHECK Maybe /area(nz,n) -> [mmol/m2]
        end do

        !! * Particulate Organic Nitrogen *
        if( tracers%data(tr_num)%ID==1004 .or. &  !iphyn
            tracers%data(tr_num)%ID==1007 .or. &  !idetn
            tracers%data(tr_num)%ID==1013 .or. &  !idian
            tracers%data(tr_num)%ID==1025 ) then  !idetz2n
            Benthos(n,1)= Benthos(n,1) +  add_benthos_2d(n) ![mmol]

            if (use_MEDUSA) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
               SinkFlx_tr(n,1,tr_num) = SinkFlx_tr(n,1,tr_num) + add_benthos_2d(n) / area(1,n)/dt ![mmol/m2]
        ! now SinkFlx hat the unit mmol/time step 
        ! but mmol/m2/time is needed for MEDUSA: thus /area
            endif
            if ((.not.use_MEDUSA).or.(sedflx_num.eq.0)) then  
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
               Benthos_tr(n,1,tr_num)= Benthos_tr(n,1,tr_num) +  add_benthos_2d(n) ![mmol]
            endif

        endif

        !! * Particulate Organic Carbon *
        if( tracers%data(tr_num)%ID==1005 .or. &  !iphyc
            tracers%data(tr_num)%ID==1008 .or. &  !idetc
            tracers%data(tr_num)%ID==1014 .or. &  !idiac
            tracers%data(tr_num)%ID==1026 ) then  !idetz2c
            Benthos(n,2)= Benthos(n,2) + add_benthos_2d(n)

            if (use_MEDUSA) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
               SinkFlx_tr(n,2,tr_num) = SinkFlx_tr(n,2,tr_num) + add_benthos_2d(n) / area(1,n)/dt
            endif
            if ((.not.use_MEDUSA).or.(sedflx_num.eq.0)) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
               Benthos_tr(n,2,tr_num)= Benthos_tr(n,2,tr_num) + add_benthos_2d(n)
            endif

        endif

        !! *Particulate Organic Silicon *
        if( tracers%data(tr_num)%ID==1016 .or. &  !idiasi
            tracers%data(tr_num)%ID==1017 .or. &  !idetsi
            tracers%data(tr_num)%ID==1027 ) then  !idetz2si
            Benthos(n,3)= Benthos(n,3) + add_benthos_2d(n)

            if (use_MEDUSA) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
               SinkFlx_tr(n,3,tr_num) = SinkFlx_tr(n,3,tr_num) + add_benthos_2d(n) / area(1,n)/dt
            endif
            if ((.not.use_MEDUSA).or.(sedflx_num.eq.0)) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
               Benthos_tr(n,3,tr_num)= Benthos_tr(n,3,tr_num) + add_benthos_2d(n)
            endif

        endif

        !! * Cal *
        if( tracers%data(tr_num)%ID==1020 .or. &  !iphycal
            tracers%data(tr_num)%ID==1021 .or. &  !idetcal
            tracers%data(tr_num)%ID==1028 ) then  !idetz2cal
            Benthos(n,4)= Benthos(n,4) + add_benthos_2d(n)

            if (use_MEDUSA) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
               SinkFlx_tr(n,4,tr_num) = SinkFlx_tr(n,4,tr_num) + add_benthos_2d(n) / area(1,n)/dt
            endif
            if ((.not.use_MEDUSA).or.(sedflx_num.eq.0)) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
               Benthos_tr(n,4,tr_num)= Benthos_tr(n,4,tr_num) + add_benthos_2d(n)
            endif

        endif

        ! flux of 13C into the sediment
        if (ciso) then             
            if( tracers%data(tr_num)%ID==1305 .or. & !iphyc_13
                tracers%data(tr_num)%ID==1308 .or. & !idetc_13
                tracers%data(tr_num)%ID==1314 ) then !idiac_14

                if (use_MEDUSA) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
                   SinkFlx_tr(n,5,tr_num) = SinkFlx_tr(n,5,tr_num) + add_benthos_2d(n) / area(1,n)/dt
                endif
                if ((.not.use_MEDUSA).or.(sedflx_num.eq.0)) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
                   Benthos_tr(n,5,tr_num)= Benthos_tr(n,5,tr_num) + add_benthos_2d(n)
                endif

            endif

           if( tracers%data(tr_num)%ID==1320 .or. &  !iphycal
               tracers%data(tr_num)%ID==1321 ) then  !idetcal

               if (use_MEDUSA) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
                  SinkFlx_tr(n,6,tr_num) = SinkFlx_tr(n,6,tr_num) + add_benthos_2d(n) / area(1,n)/dt
               endif
               if ((.not.use_MEDUSA).or.(sedflx_num.eq.0)) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
                  Benthos_tr(n,6,tr_num)= Benthos_tr(n,6,tr_num) + add_benthos_2d(n)
               endif

           endif

        endif
        
        ! flux of 14C into the sediment
        if (ciso .and. ciso_organic_14) then             
           if( tracers%data(tr_num)%ID==1405 .or. & !iphyc_13
               tracers%data(tr_num)%ID==1408 .or. & !idetc_13
               tracers%data(tr_num)%ID==1414 ) then !idiac_14

               if (use_MEDUSA) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
                  SinkFlx_tr(n,7,tr_num) = SinkFlx_tr(n,7,tr_num) + add_benthos_2d(n) / area(1,n)/dt
               endif
               if ((.not.use_MEDUSA).or.(sedflx_num.eq.0)) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
                  Benthos_tr(n,7,tr_num)= Benthos_tr(n,7,tr_num) + add_benthos_2d(n)
               endif

           endif

           if( tracers%data(tr_num)%ID==1420 .or. &  !iphycal
               tracers%data(tr_num)%ID==1421 ) then  !idetcal
               if (use_MEDUSA) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
                  SinkFlx_tr(n,8,tr_num) = SinkFlx_tr(n,8,tr_num) + add_benthos_2d(n) / area(1,n)/dt
               endif
               if ((.not.use_MEDUSA).or.(sedflx_num.eq.0)) then
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
                  Benthos_tr(n,8,tr_num)= Benthos_tr(n,8,tr_num) + add_benthos_2d(n)
               endif
           endif

        endif

   end do

   if(use_MEDUSA) then
        do n=1, bottflx_num
!           SinkFlx(:,n) = Sinkflx(:,n)/dt
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
           call exchange_nod(SinkFlx_tr(:,n,tr_num), partit)
        end do
   end if ! use_MEDUSA

   do n=1, benthos_num
! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
      call exchange_nod(Benthos_tr(:,n,tr_num), partit)

      call exchange_nod(Benthos(:,n), partit)
   end do

end subroutine ver_sinking_recom_benthos
!
!
!===============================================================================
subroutine diff_ver_recom_expl(tr_num, tracers, partit, mesh)
! Remineralization from benthos
! bottom_flux

    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_TRACER

    use recom_declarations
    use recom_locvar
    use recom_glovar
    use recom_config
    use recom_ciso
    ! use diff_ver_recom_expl_interface

    use g_clock
    use o_PARAM
    use g_config
    use o_param 
    use o_arrays
    use g_forcing_arrays
    use g_comm_auto

    IMPLICIT NONE

    integer       , intent(in)   , target  :: tr_num
    type(t_tracer), intent(inout), target  :: tracers
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh

    integer                                :: elem,k
    integer                                :: n2,nl1,nl2,nz,n,id,ul1
    real(kind=WP)                          :: vd_flux(mesh%nl)
    integer                                :: nlevels_nod2D_minimum
    real(kind=WP)                          :: bottom_flux(partit%myDim_nod2D+partit%eDim_nod2D)

    real(kind=WP), dimension(:,:), pointer :: trarr

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    trarr=>tracers%data(tr_num)%values(:,:)


    bottom_flux = 0._WP
    id = tracers%data(tr_num)%ID

#if defined(__recom)
if (use_MEDUSA .and. (sedflx_num .ne. 0)) then
   !CV update: the calculation later has been changed by Ozgur in such
   !a way  that now the  variable bottom_flux is in  (mol/time) units,
   !rather than  a flux in  (mol/time/area). I therefore  multiply the
   !Medusa fluxes by the area to get the same unit.

   SELECT CASE (id)
    CASE (1001)
      bottom_flux = GloSed(:,1) * area(1,:) ! DIN
    CASE (1002)
      bottom_flux = GloSed(:,2) * area(1,:) ! DIC
    CASE (1003)
      bottom_flux = GloSed(:,3) * area(1,:) ! Alk
    CASE (1018)
      bottom_flux = GloSed(:,4) * area(1,:) ! Si
    CASE (1019)
      bottom_flux = GloSed(:,1) * Fe2N_benthos * area(1,:)
    CASE (1022)
      bottom_flux = GloSed(:,5) * area(1,:) ! Oxy
    CASE (1302)
      if (ciso) then
        bottom_flux = GloSed(:,6) * area(1,:) ! DIC_13 and Calc: DIC_13
      end if
    CASE (1402)
      if (ciso) then
        bottom_flux = GloSed(:,7) * area(1,:) ! DIC_14 and Calc: DIC_14
      end if
    CASE DEFAULT
      if (partit%mype==0) then
        write(*,*) 'check specified in boundary conditions'
        write(*,*) 'the model will stop!'
      end if
      call par_ex(partit%MPI_COMM_FESOM, partit%mype)
      stop
  END SELECT
else
    SELECT CASE (id)
       CASE (1001)
          bottom_flux = GlodecayBenthos(:,1) !*** DIN [mmolN/m^2/s] ***
       CASE (1002)
          bottom_flux = GlodecayBenthos(:,2) + GlodecayBenthos(:,4) !*** DIC + calcification ***
       CASE (1003)
          bottom_flux = GlodecayBenthos(:,4) * 2.0_WP - 1.0625_WP * GlodecayBenthos(:,1) !*** Alk ***
       CASE (1018)
          bottom_flux = GlodecayBenthos(:,3) !*** Si ***
       CASE (1019)
          bottom_flux = GlodecayBenthos(:,1) * Fe2N_benthos !*** DFe ***
       CASE (1022)
          bottom_flux = -GlodecayBenthos(:,2) * redO2C !*** O2 ***
       CASE (1302)
         if (ciso) then
           bottom_flux = GlodecayBenthos(:,5) + GlodecayBenthos(:,6) !*** DIC_13 and Calc: DIC_13 ***
         end if
       CASE (1402)
         if (ciso) then
           bottom_flux = GlodecayBenthos(:,7) + GlodecayBenthos(:,8) !*** DIC_14 and Calc: DIC_14 ***
         end if
       CASE DEFAULT
          if (partit%mype==0) then
             write(*,*) 'check specified in boundary conditions'
             write(*,*) 'the model will stop!'
          end if
          call par_ex(partit%MPI_COMM_FESOM, partit%mype)
          stop
    END SELECT
endif ! (use_MEDUSA .and. (sedflux_num .gt. 0))  
#endif

    do n=1, myDim_nod2D

        nl1=nlevels_nod2D(n)-1
        ul1=ulevels_nod2D(n)

        vd_flux=0._WP

        k=nod_in_elem2D_num(n)
        ! Screening minimum depth in neigbouring nodes around node n
        nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

        !_______________________________________________________________________
        ! Bottom flux
        do nz=nlevels_nod2D_minimum, nl1
            vd_flux(nz)=(area(nz,n)-area(nz+1,n))* bottom_flux(n)/(area(1,n))  !!!!!!!!
        end do
        nz=nl1
        vd_flux(nz+1)= (area(nz+1,n))* bottom_flux(n)/(area(1,n))
        !_______________________________________________________________________
        ! writing flux into rhs
        do nz=ul1,nl1
            ! flux contribute only the cell through its bottom !!!
!            dtr_bf(nz,n) = dtr_bf(nz,n) + vd_flux(nz+1)*dt/area(nz,n)/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
            dtr_bf(nz,n) = dtr_bf(nz,n) + vd_flux(nz+1)*dt/areasvol(nz,n)/hnode_new(nz,n)
        end do
    end do
end subroutine diff_ver_recom_expl

subroutine ver_sinking_recom(tr_num, tracers, partit, mesh)
! Sinking in water column

    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_TRACER

    use REcoM_declarations
    use REcoM_LocVar
    use REcoM_GloVar
    use recom_config
    use REcoM_ciso
    ! use ver_sinking_recom_interface

    use g_clock
    use o_PARAM
    use g_config
    use o_param
    use o_arrays
    use g_forcing_arrays
    use g_comm_auto
    implicit none

    integer       , intent(in)   , target  :: tr_num
    type(t_tracer), intent(inout), target  :: tracers
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh

    integer                                :: node, nz,  id, nzmin, nzmax, n, k, nlevels_nod2D_minimum
    real(kind=WP)                          :: vd_flux(mesh%nl)
    real(kind=8)                           :: dz_trr(mesh%nl), aux
    real(kind=8)                           :: wLoc,wM,wPs
    real(kind=8)                           :: Rjp,Rj,Rjm

    real(kind=8)                           :: cfl, d0, d1, thetaP, thetaM, psiP, psiM
    real(kind=8)                           :: onesixth	= 1.d0/6.d0
    real(kind=8)                           :: dt_sink, c1, c2
    real(kind=8)                           :: Vsink, tv
    real(kind=8),dimension(mesh%nl)        :: Wvel_flux

    real(kind=WP), dimension(:,:), pointer :: trarr

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    trarr=>tracers%data(tr_num)%values(:,:)

!< calculate scaling factors
!< scaling_density1_3D, scaling_density2_3D
!< scaling_visc_3D
!< .OG. 04.11.2022

!< Constant sinking velocities (we prescribe them under namelist recom)
!< This hardcoded part is temporary
!< .OG. 07.07.2021

    Vsink=0.0_WP

    if (tracers%data(tr_num)%ID ==1007 .or.    &  !idetn
        tracers%data(tr_num)%ID ==1008 .or.    &  !idetc
        tracers%data(tr_num)%ID ==1017 .or.    &  !idetsi
        tracers%data(tr_num)%ID ==1021 ) then     !idetcal

            Vsink = VDet

    elseif(tracers%data(tr_num)%ID ==1004 .or. &  !iphyn
        tracers%data(tr_num)%ID ==1005 .or.    &  !iphyc
        tracers%data(tr_num)%ID==1006 ) then     !ipchl

            Vsink = VPhy

    elseif(tracers%data(tr_num)%ID==1013 .or. &  !idian
        tracers%data(tr_num)%ID==1014 .or.    &  !idiac
        tracers%data(tr_num)%ID==1016 .or.    &  !idiasi
        tracers%data(tr_num)%ID==1015 ) then     !idchl

            Vsink = VDia

#if defined (__coccos)
    elseif(tracers%data(tr_num)%ID == 1029 .or. &  !icocn
        tracers%data(tr_num)%ID == 1030 .or.    &  !icocc
        tracers%data(tr_num)%ID == 1031 ) then     !icchl

            Vsink = VCocco
#endif

    elseif(tracers%data(tr_num)%ID == 1020) then   !iphycal
       
#if defined (__coccos)
            Vsink = VCocco
#else
            Vsink = VPhy
#endif
            
#if defined (__3Zoo2Det)
    elseif(tracers%data(tr_num)%ID==1025 .or. &  !idetz2n
           tracers%data(tr_num)%ID==1026 .or. &  !idetz2c
           tracers%data(tr_num)%ID==1027 .or. &  !idetz2si
           tracers%data(tr_num)%ID==1028 ) then  !idetz2calc 
            
            Vsink = VDet_zoo2
#endif
    end if

!! ---- No sinking if Vsink < 0.1 m/day
if (Vsink .gt. 0.1) then 

   do n = 1,myDim_nod2D
      if (ulevels_nod2D(n)>1) cycle
      nzmin = ulevels_nod2D(n)
      nzmax = nlevels_nod2D(n)-1

      !! distance between tracer points, surface and bottom dz_trr is half
      !! the layer thickness
      dz_trr                = 0.0d0
      dz_trr(nzmin+1:nzmax) = abs(Z_3d_n(nzmin:nzmax-1,n)-Z_3d_n(nzmin+1:nzmax,n))
      dz_trr(nzmin)         = hnode(nzmin,n)/2.0d0
      dz_trr(nzmax+1)       = hnode(nzmax,n)/2.0d0

      Wvel_flux(nzmin:nzmax+1)= 0.d0  ! Vertical velocity for BCG tracers

      do nz=nzmin,nzmax+1

         Wvel_flux(nz) = -Vsink/SecondsPerDay ! allow_var_sinking = .false.

         if (allow_var_sinking) then
            Wvel_flux(nz) = -((Vdet_a * abs(zbar_3d_n(nz,n))/SecondsPerDay) + Vsink/SecondsPerDay)
            if (use_ballasting) then
                Wvel_flux(nz) = w_ref1 * scaling_density1_3D(nz,n) * scaling_visc_3D(nz,n)

                if (depth_scaling1.gt.0.0) Wvel_flux(nz) = Wvel_flux(nz) + (depth_scaling1 * abs(zbar_3d_n(nz,n)))

                if (abs(Wvel_flux(nz)) .gt. max_sinking_velocity) Wvel_flux(nz) = max_sinking_velocity

                !! * sinking velocity [m d-1] surface --> bottom (negative)* 
                Wvel_flux(nz) = -1.0d0 * Wvel_flux(nz)/SecondsPerDay ! now in [m s-1]
            endif
         end if

#if defined (__3Zoo2Det)

         !! ---- We assume *constant* sinking for second detritus
         if(tracers%data(tr_num)%ID ==1025 .or. &  !idetz2n
            tracers%data(tr_num)%ID ==1026 .or. &  !idetz2c
            tracers%data(tr_num)%ID ==1027 .or. &  !idetz2si
            tracers%data(tr_num)%ID ==1028 ) then  !idetz2calc
               Wvel_flux(nz) = -VDet_zoo2/SecondsPerDay

               if (use_ballasting) then

                  Wvel_flux(nz) = w_ref2*scaling_density2_3D(nz,n)*scaling_visc_3D(nz,n)

                  if (depth_scaling2.gt.0.0) Wvel_flux(nz) = Wvel_flux(nz) + (depth_scaling2 * abs(zbar_3d_n(nz,n)))

                  if (abs(Wvel_flux(nz)) .gt. max_sinking_velocity) Wvel_flux(nz) = max_sinking_velocity

                  !! * sinking velocity [m d-1] surface --> bottom (negative) *
                  Wvel_flux(nz) = -1.0d0 * Wvel_flux(nz)/SecondsPerDay ! now in [m s-1]
               end if

         endif
#endif
      end do

      dt_sink = dt
      vd_flux = 0.0d0

!FIXME: Having IF True and IF False is bad practice. Either throw away the old code, or make a namelist switch...
if (.TRUE.) then ! 3rd Order DST Sceheme with flux limiting. This code comes from old recom

      k=nod_in_elem2D_num(n)
      ! Screening minimum depth in neigbouring nodes around node n
      nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

      vd_flux(nzmin:nzmax+1)= 0.0_WP

      do nz=nzmax, nzmin+1,-1

         Rjp = trarr(nz,n)              - trarr(min(nz+1,nzmax),n)
         Rj  = trarr(max(nzmin,nz-1),n) - trarr(nz,n)
         Rjm = trarr(max(nzmin,nz-2),n) - trarr(max(nzmin,nz-1),n)

         cfl = abs(Wvel_flux(nz) * dt_sink / dz_trr(nz)) !(Z_n(nz-1)-Z_n(nz)))       ! [m/day] * [day] * [1/m]  ! NEW BALL changed dt to dt_sink

         wPs = Wvel_flux(nz) + abs(Wvel_flux(nz)) ! --> Positive vertical velocity
         wM  = Wvel_flux(nz) - abs(Wvel_flux(nz)) ! --> Negative vertical velocity

         d0 = (2.d0 - cfl)*(1.d0 - cfl)*onesixth
         d1 = (1.d0 - cfl*cfl)*onesixth

         thetaP = Rjm/(1.d-20+Rj)
         psiP = d0 + d1*thetaP
         psiP = max(0.d0, min(min(1.d0,psiP), &
            (1.d0-cfl)/(1.d-20+cfl)*thetaP))

         thetaM = Rjp/(1.d-20 + Rj)
         psiM = d0 + d1*thetaM
         psiM = max(0.d0, min(min(1.d0,psiM), &
            (1.d0-cfl)/(1.d-20-cfl)*thetaM))

         tv= (0.5 * wPs * (trarr(nz,n)              + psiM * Rj)+ &
	      0.5 * wM  * (trarr(max(nzmin,nz-1),n) + psiP * Rj))
         vd_flux(nz)= - tv*area(nz,n)
      end do
end if ! 3rd Order DST Sceheme with flux limiting

if (.FALSE.) then ! simple upwind

      ! Surface flux
      vd_flux(nzmin)= 0.0_WP

      ! Bottom flux
      vd_flux(nzmax+1)= 0.0_WP

      k=nod_in_elem2D_num(n)
      ! Screening minimum depth in neigbouring nodes around node n
      nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

      do nz=nzmin+1,nzmax !nlevels_nod2D_minimum-1
!         tv = trarr(nz,n)                                ! simple scheme       - test1
!         tv = 0.5_WP*(trarr(nz-1,n)+trarr(nz,n))        ! consider both layers - test2
!         tv = tv*Wvel_flux(nz) ! Wvel_flux is negative
         tv = - 0.5* & ! - test3
            (trarr(nz-1,n)*(Wvel_flux(nz)-abs(Wvel_flux(nz))) + &
             trarr(nz  ,n)*(Wvel_flux(nz)+abs(Wvel_flux(nz))))
         vd_flux(nz)= tv*area(nz,n)

      end do
end if ! simple upwind
      do nz=nzmin,nzmax
         vert_sink(nz,n) = vert_sink(nz,n) + (vd_flux(nz)-vd_flux(nz+1))*dt/areasvol(nz,n)/hnode_new(nz,n) !/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
      end do
   end do
end if ! Vsink .gt. 0.1

end subroutine ver_sinking_recom
!-------------------------------------------------------------------------------
! Subroutine calculate ballasting
!-------------------------------------------------------------------------------
subroutine ballast(tr_num, tracers, partit, mesh)

    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_TRACER

    use recom_config
    use recom_glovar

    USE o_PARAM
    USE o_ARRAYS
    USE g_CONFIG
    use g_forcing_arrays
    use g_comm_auto
    use g_clock
    use g_rotate_grid
    use mvars
    use mdepth2press
    use gsw_mod_toolbox, only: gsw_sa_from_sp,gsw_ct_from_pt,gsw_rho

    implicit none
    integer       , intent(in)   , target :: tr_num
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    integer                               :: row, k, nzmin, nzmax
    real(kind=8)                          :: depth_pos(1)
    real(kind=8)                          :: pres(1)
    real(kind=8)                          :: sa(1)
    real(kind=8)                          :: ct(1)
    real(kind=8)                          :: rho_seawater(1)
    real(kind=8)                          :: Lon_degree(1)
    real(kind=8)                          :: Lat_degree(1)

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

  ! For ballasting, calculate scaling factors here and pass them to FESOM, where sinking velocities are calculated
     ! -----
     ! If ballasting is used, sinking velocities are a function of a) particle composition (=density),
     ! b) sea water viscosity, c) depth (currently for small detritus only), and d) a constant reference sinking speed
     ! -----

     !___________________________________________________________________________
     ! loop over local nodes
     do row=1,myDim_nod2D
         ! max. number of levels at node n
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)
         !! lon
        Lon_degree(1)=geo_coord_nod2D(1,row)/rad !! convert from rad to degree
         !! lat
        Lat_degree(1)=geo_coord_nod2D(2,row)/rad !! convert from rad to degree

        ! get scaling vectors -> these need to be passed to FESOM to get sinking velocities
        ! get local seawater density
        do k=nzmin, nzmax

           !! level depth
           depth_pos(1) = abs(Z_3d_n(k,row))  ! take depth of tracers instead of levels abs(zbar_3d_n(k,row))

           call depth2press(depth_pos(1), Lat_degree(1), pres, 1)  ! pres is output of function,1=number of records
           sa           = gsw_sa_from_sp(tracers%data(2)%values(k,row), pres, Lon_degree(1), Lat_degree(1))
           ct           = gsw_ct_from_pt(sa,tracers%data(1)%values(k,row))
           rho_seawater = gsw_rho(sa, ct, pres)

           ! (i.e. no density scaling)
           scaling_density1_3D(k,row)=1.0
           scaling_density2_3D(k,row)=1.0

              if (use_density_scaling) then
                 if (tracers%data(tr_num)%ID ==1008)then !idetc
                    if (tracers%data(tr_num)%values(k,row)>0.001) then ! only apply ballasting above a certain biomass (OG Todo: remove) 
                       scaling_density1_3D(k,row) = (rho_particle1(k,row)-rho_seawater(1))/(rho_ref_part-rho_ref_water)
                    endif
                 endif
#if defined (__3Zoo2Det)

                    if (tracers%data(tr_num)%ID ==1026)then ! idetz2c
                       if (tracers%data(tr_num)%values(k,row)>0.001) then ! only apply ballasting above a certain biomass (OG Todo: remove) 
                          scaling_density2_3D(k,row) = (rho_particle2(k,row)-rho_seawater(1))/(rho_ref_part-rho_ref_water)
                       endif
                    endif
#endif
              endif

            scaling_visc_3D(k,row)=1.0

            if (use_viscosity_scaling) then
                if (seawater_visc_3D(k,row)==0) then
                    scaling_visc_3D(k,row)=1.0
                else
                    scaling_visc_3D(k,row)= visc_ref_water/seawater_visc_3D(k,row)
                endif
            endif

        end do
        rho_particle1(nzmax+1,row) = rho_particle1(nzmax,row)
        rho_particle2(nzmax+1,row) = rho_particle2(nzmax,row)
        scaling_visc_3D(nzmax+1,row) = scaling_visc_3D(nzmax,row)
     end do
    ! in the unlikely (if possible at all...) case that rho_particle(k)-rho_seawater(1)<0, prevent the scaling factor from being negative

    if (any(scaling_density1_3D(:,:) <= tiny)) scaling_density1_3D(:,:) = 1.0_WP      ! tiny = 2.23D-16
#if defined (__3Zoo2Det)
    if (any(scaling_density2_3D(:,:) <= tiny)) scaling_density2_3D(:,:) = 1.0_WP      ! tiny = 2.23D-16
#endif
end subroutine ballast
!-------------------------------------------------------------------------------
! Subroutine calculate density of particle
! depending on composition (detC, detOpal, detCaCO3) based on Cram et al. (2018)
!-------------------------------------------------------------------------------
subroutine get_particle_density(tracers, partit, mesh)

    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_TRACER

    use recom_config
    use recom_glovar
    USE o_PARAM
    USE o_ARRAYS
    USE g_CONFIG
    use g_forcing_arrays
    use g_comm_auto
    use g_clock
    use g_rotate_grid

    implicit none
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh

    integer                               :: row, k, nzmin, nzmax, tr_num, num_tracers

    real(kind=8)                          :: a1(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D) ! [n.d.] fraction of carbon in detritus class
    real(kind=8)                          :: a2(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D) ! [n.d.] fraction of nitrogen in detritus class
    real(kind=8)                          :: a3(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D) ! [n.d.] fraction of Opal in detritus class
    real(kind=8)                          :: a4(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D) ! [n.d.] fraction of CaCO3 in detritus class
    real(kind=8)                          :: b1(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=8)                          :: b2(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=8)                          :: b3(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=8)                          :: b4(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=8)                          :: aux(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    num_tracers=tracers%num_tracers

    rho_particle1 = 0.0
    b1 = 0.0
    b2 = 0.0
    b3 = 0.0
    b4 = 0.0
    aux = 0.0

    do tr_num=1,num_tracers
        if (tracers%data(tr_num)%ID==1008)  b1 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetc      ! [mmol m-3] detritus carbon
        if (tracers%data(tr_num)%ID==1007)  b2 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetn      ! [mmol m-3] detritus nitrogen
        if (tracers%data(tr_num)%ID==1017)  b3 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetsi     ! [mmol m-3] detritus Si
        if (tracers%data(tr_num)%ID==1021)  b4 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetcal    ! [mmol m-3] detritus CaCO3
    end do

    do row=1,myDim_nod2d
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)
        aux(nzmin:nzmax,row) = b1(nzmin:nzmax,row)+b2(nzmin:nzmax,row)+b3(nzmin:nzmax,row)+b4(nzmin:nzmax,row)
        a1(nzmin:nzmax,row)  = b1(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a2(nzmin:nzmax,row)  = b2(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a3(nzmin:nzmax,row)  = b3(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a4(nzmin:nzmax,row)  = b4(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        rho_particle1(nzmin:nzmax,row) = rho_CaCO3*a4(nzmin:nzmax,row) + rho_opal*a3(nzmin:nzmax,row) + rho_POC*a1(nzmin:nzmax,row) + rho_PON*a2(nzmin:nzmax,row)
    end do

#if defined (__3Zoo2Det)
    rho_particle2 = 0.0
    b1 = 0.0
    b2 = 0.0
    b3 = 0.0
    b4 = 0.0
    aux = 0.0
    do tr_num=1,num_tracers
        if (tracers%data(tr_num)%ID==1026)  b1 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetz2c
        if (tracers%data(tr_num)%ID==1025)  b2 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetz2n
        if (tracers%data(tr_num)%ID==1027)  b3 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetz2si
        if (tracers%data(tr_num)%ID==1028)  b4 = max(tiny,tracers%data(tr_num)%values(:,:)) !idetz2calc
    end do

    do row=1,myDim_nod2d+eDim_nod2D   ! myDim is sufficient
        !if (ulevels_nod2D(row)>1) cycle
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)
        aux(nzmin:nzmax,row) = b1(nzmin:nzmax,row)+b2(nzmin:nzmax,row)+b3(nzmin:nzmax,row)+b4(nzmin:nzmax,row)
        a1(nzmin:nzmax,row)  = b1(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a2(nzmin:nzmax,row)  = b2(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a3(nzmin:nzmax,row)  = b3(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        a4(nzmin:nzmax,row)  = b4(nzmin:nzmax,row)/aux(nzmin:nzmax,row)
        rho_particle2(nzmin:nzmax,row) = rho_CaCO3*a4(nzmin:nzmax,row) + rho_opal*a3(nzmin:nzmax,row) + rho_POC*a1(nzmin:nzmax,row) + rho_PON*a2(nzmin:nzmax,row)
    end do
#endif

end subroutine get_particle_density
!-------------------------------------------------------------------------------
! Subroutine to approximate seawater viscosity with current temperature
! based on Cram et al. (2018)
!-------------------------------------------------------------------------------

! neglecting salinity effects, which are much smaller than those of temperature
! https://bitbucket.org/ohnoplus/ballasted-sinking/src/master/tools/waterviscosity.m

subroutine get_seawater_viscosity(tr_num, tracers, partit, mesh)

    use MOD_MESH
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_TRACER

    use recom_config
    use recom_glovar
    USE o_PARAM
    USE o_ARRAYS
    USE g_CONFIG
    use g_forcing_arrays
    use g_comm_auto
    use g_clock
    use g_rotate_grid

    implicit none

!!  temp [degrees C] Ocean temperature
!!  salt [g/kg or n.d.] Ocean salinity
!!  seawater_visc_3D [kg m-1 s-1] Ocean viscosity

    real(kind=8),dimension(1)             :: A, B, mu_w
    integer                               :: row, k, nzmin, nzmax

    integer       , intent(in)   , target :: tr_num
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    seawater_visc_3D(:,:) = 0.0
    do row=1,myDim_nod2d
     !if (ulevels_nod2D(row)>1) cycle
! OG Do we need any limitation here?
! i.e., if (seawater_visc_3D(row)<=0.0_WP) cycle
        nzmin = ulevels_nod2D(row)
        nzmax = nlevels_nod2D(row)

        do k=nzmin, nzmax
     ! Eq from Sharaway 2010
     ! validity:
     !  0<temp<180 degC
     !  0<salt<0.15 kg/kg
     ! Note: because salinity is expected to be in kg/kg, use conversion factor 0.001 below!
            A(1) = 1.541 + 1.998*0.01*tracers%data(1)%values(k,row) - 9.52*1e-5*tracers%data(1)%values(k,row)*tracers%data(1)%values(k,row)
            B(1) = 7.974 - 7.561*0.01*tracers%data(1)%values(k,row) + 4.724*1e-4*tracers%data(1)%values(k,row)*tracers%data(1)%values(k,row)
            mu_w(1) = 4.2844*1.0e-5 + (1.0/(0.157*(tracers%data(1)%values(k,row)+64.993)*(tracers%data(1)%values(k,row)+64.993)-91.296))
            seawater_visc_3D(k,row) = mu_w(1) * (1.0 + A(1)*tracers%data(2)%values(k,row)*0.001 + B(1)*tracers%data(2)%values(k,row)*0.001*tracers%data(2)%values(k,row)*0.001)
        enddo
    end do

end subroutine get_seawater_viscosity
