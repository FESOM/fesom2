!
!
!===============================================================================
subroutine pressure_bv
! fill in the hydrostatic pressure and the Brunt-Vaisala frequency 
! in a single pass the using split form of the equation of state
! as proposed by NR
	USE o_PARAM
	USE o_MESH
	USE o_ARRAYS
	USE g_PARSUP
	USE g_config
	use i_arrays
	USE o_mixing_KPP_mod, only: dbsfc
	IMPLICIT NONE
	
	real(kind=WP)         :: dz_inv, bv,  a, rho_up, rho_dn, t, s
	integer               :: node, nz, nl1, nzmax
	real(kind=WP)         :: rhopot(nl), bulk_0(nl), bulk_pz(nl), bulk_pz2(nl), rho(nl), dbsfc1(nl), db_max
	real(kind=WP)         :: bulk_up, bulk_dn, smallvalue, buoyancy_crit, rho_surf
	real(kind=WP)         :: sigma_theta_crit=0.125   !kg/m3, Levitus threshold for computing MLD2
	logical               :: flag1, flag2, mixing_kpp

	smallvalue=1.0e-20
	buoyancy_crit=0.0003
	mixing_kpp = (trim(mix_scheme)=='KPP')  ! NR Evaluate string comparison outside the loop. It is expensive.
	!___________________________________________________________________________
	! Screen salinity
	a=0.0_WP
	do node=1, myDim_nod2D+eDim_nod2D
		do nz=1,nlevels_nod2d(node)-1
			a=min(a,tr_arr(nz,node,2))
		enddo
	enddo
	
	!___________________________________________________________________________
	if(a<0.) then
		write (*,*)' --> pressure_bv: s<0 happens!', a
		pe_status=1
		do node=1, myDim_nod2D+eDim_nod2D
			do nz=1, nlevels_nod2d(node)-1
				if (tr_arr(nz, node, 2) < 0) write (*,*) 'the model blows up at n=', mylist_nod2D(node), ' ; ', 'nz=', nz
			end do
		end do
	endif
	
	!___________________________________________________________________________
	if(use_ALE) then
		!_______________________________________________________________________
		do node=1, myDim_nod2D+eDim_nod2D
			nl1= nlevels_nod2d(node)-1
			! also compute the maximum buoyancy gradient between the surface and any depth
			! it will be used for computing MLD according to FESOM 1.4 implementation (after Large et al. 1997)
			db_max=0.0
			!___________________________________________________________________
			do nz=1, nl1
				t=tr_arr(nz, node,1)
				s=tr_arr(nz, node,2)
				call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz))
			enddo
			!NR split the loop here. The Intel compiler could not resolve that there is no dependency 
			!NR and did not vectorize the full loop. 
			do nz=1, nl1
				rho(nz)= bulk_0(nz)   + Z_3d_n(nz,node)*(bulk_pz(nz)   + Z_3d_n(nz,node)*bulk_pz2(nz))
				rho(nz)=rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*Z_3d_n(nz,node))-density_0
				density_m_rho0(nz,node) = rho(nz) 
				! squared buoyancy difference between the surface and the grid points blow (adopted from FESOM 1.4)
				rho_surf=bulk_0(1)   + Z_3d_n(nz,node)*(bulk_pz(1)   + Z_3d_n(nz,node)*bulk_pz2(1))
				rho_surf=rho_surf*rhopot(1)/(rho_surf+0.1_WP*Z_3d_n(nz,node))-density_0
				dbsfc1(nz) = -g * ( rho_surf - rho(nz) ) / (rho(nz)+density_0)      ! this is also required when KPP is ON
				db_max=max(dbsfc1(nz)/abs(Z_3d_n(1,node)-Z_3d_n(max(nz, 2),node)), db_max)
			end do
			dbsfc1(nl)=dbsfc1(nl1)
			if (mixing_kpp) then ! in case KPP is ON store the buoyancy difference with respect to the surface (m/s2)
				dbsfc(1:nl, node )=dbsfc1(1:nl)
			end if
			!___________________________________________________________________
			! Pressure
			hpressure(1, node)=-Z_3d_n(1,node)*rho(1)*g
			DO nz=2, nl1
				a=0.5_WP*g*(rho(nz-1)*(zbar_3d_n(nz-1,node)-zbar_3d_n(nz,node))+rho(nz)*(zbar_3d_n(nz,node)-zbar_3d_n(nz+1,node)))
				hpressure(nz, node)=hpressure(nz-1, node)+a
			END DO
			
			!___________________________________________________________________
			! BV frequency:  bvfreq(nl,:), squared value is stored   
			MLD1(node)=Z_3d_n(2,node)
			MLD2(node)=Z_3d_n(2,node)
			MLD1_ind(node)=2
			MLD2_ind(node)=2
			flag1=.true.
			flag2=.true.
			DO nz=2,nl1
				bulk_up = bulk_0(nz-1) + zbar_3d_n(nz,node)*(bulk_pz(nz-1) + zbar_3d_n(nz,node)*bulk_pz2(nz-1)) 
				bulk_dn = bulk_0(nz)   + zbar_3d_n(nz,node)*(bulk_pz(nz)   + zbar_3d_n(nz,node)*bulk_pz2(nz))
				rho_up = bulk_up*rhopot(nz-1) / (bulk_up + 0.1*zbar_3d_n(nz,node))  
				rho_dn = bulk_dn*rhopot(nz)   / (bulk_dn + 0.1*zbar_3d_n(nz,node))  
				dz_inv=1.0_WP/(Z_3d_n(nz-1,node)-Z_3d_n(nz,node))  
				bvfreq(nz,node)  = -g*dz_inv*(rho_up-rho_dn)/density_0
				! Define MLD following FESOM 1.4 implementation (after Large et al. 1997) 
				! MLD is the shallowest depth where the local buoyancy gradient matches the maximum buoyancy gradient 
				! between the surface and any discrete depth within the water column.
				if (bvfreq(nz, node) > db_max .and. flag1) then
					MLD1(node)    =Z_3d_n(nz, node)
					MLD1_ind(node)=nz
					flag1=.false.
				end if
				! another definition of MLD after Levitus
				if ((rhopot(nz)-rhopot(1) > sigma_theta_crit) .and. flag2) then
					MLD2(node)=MLD2(node)+(Z_3d_n(nz,node)-MLD2(node))/(rhopot(nz)-rhopot(nz-1)+1.e-20)*(rhopot(1)+sigma_theta_crit-rhopot(nz-1))
					MLD2_ind(node)=nz
					flag2=.false.
				elseif (flag2) then
					MLD2(node)=Z_3d_n(nz,node)
				end if
			END DO
			if (flag2) MLD2_ind(node)=nl1
			bvfreq(1,node)=bvfreq(2,node)
			bvfreq(nl1+1,node)=bvfreq(nl1,node) 
			!___________________________________________________________________
			! The mixed layer depth 
			! mixlay_depth    
			! bv_ref
		end do
		!_______________________________________________________________________
	else
		do node=1, myDim_nod2D+eDim_nod2D
			nl1= nlevels_nod2d(node)-1
			!___________________________________________________________________
			do nz=1, nl1
				t=tr_arr(nz, node,1)
				s=tr_arr(nz, node,2)
				call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz))
				rho(nz)= bulk_0(nz)   + Z(nz)*(bulk_pz(nz)   + Z(nz)*bulk_pz2(nz))
				rho(nz)=rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*Z(nz))-density_0
			end do
			
			!___________________________________________________________________
			! Pressure
			hpressure(1, node)=-Z(1)*rho(1)*g
			DO nz=2, nl1
				a=0.5_WP*g*(rho(nz-1)*(zbar(nz-1)-zbar(nz))+rho(nz)*(zbar(nz)-zbar(nz+1)))
				hpressure(nz, node)=hpressure(nz-1, node)+a
			END DO
			
			!___________________________________________________________________
			! BV frequency:  bvfreq(nl,:), squared value is stored   
			DO nz=2,nl1
				bulk_up = bulk_0(nz-1) + zbar(nz)*(bulk_pz(nz-1) + zbar(nz)*bulk_pz2(nz-1)) 
				bulk_dn = bulk_0(nz)   + zbar(nz)*(bulk_pz(nz)   + zbar(nz)*bulk_pz2(nz))
				rho_up = bulk_up*rhopot(nz-1) / (bulk_up + 0.1*zbar(nz))  
				rho_dn = bulk_dn*rhopot(nz)   / (bulk_dn + 0.1*zbar(nz))  
				dz_inv=1.0_WP/(Z(nz-1)-Z(nz))  
				bvfreq(nz,node)  = -g*dz_inv*(rho_up-rho_dn)/density_0
			END DO
			bvfreq(1,node)=bvfreq(2,node)
			bvfreq(nl1+1,node)=bvfreq(nl1,node) 
			!___________________________________________________________________
			! The mixed layer depth 
			! mixlay_depth    
			! bv_ref
		end do
	end if
	! BV is defined on full levels except for the first and the last ones.
end subroutine pressure_bv
!
!
!===============================================================================
! Calculate pressure gradient force (PGF) on elements by interpolating the nodal 
! density on nodal depthlayers on the depth of the elments using a cubic-spline 
! interpolation.
! First coded by Q. Wang for FESOM1.4, adapted by P. Scholz for FESOM2.0
! 26.04.2018
subroutine pressure_force
	use o_PARAM
	use o_MESH
	use o_ARRAYS
	use g_PARSUP
	use g_config
	implicit none
	
	integer 			:: elem, elnodes(3), nle, nln(3), nlz, flag
	integer 			:: ni, node, dens_ind,kk
	real(kind=WP)		:: ze
	integer				:: s_ind(4)
	real(kind=WP)		:: s_z(4), s_dens(4), s_H, aux1, aux2, aux(2), s_dup, s_dlo
	real(kind=WP)		:: a, b, c, d, dz, rho_n(3), rhograd_e(2), p_grad(2)
	
	if     (trim(which_ale)=='linfs' .and. (use_partial_cell==.false.) ) then
		!___________________________________________________________________________
		! loop over triangular elemments
		do elem=1, myDim_elem2D
			!_______________________________________________________________________
			! number of levels at elem
			nle=nlevels(elem)-1
			
			!_______________________________________________________________________
			! node indices of elem 
			elnodes = elem2D_nodes(:,elem)
			
			!_______________________________________________________________________
			! loop over depth levels at element elem
			do nlz=1,nle
				pgf_x(nlz,elem) = sum(gradient_sca(1:3,elem)*hpressure(nlz,elnodes)/density_0)
				pgf_y(nlz,elem) = sum(gradient_sca(4:6,elem)*hpressure(nlz,elnodes)/density_0)
			end do 
		end do !-->do elem=1, myDim_elem2D
	else 
		!___________________________________________________________________________
        ! loop over triangular elemments
        do elem=1, myDim_elem2D
            !_______________________________________________________________________
            ! calculate mid depth element level --> Z_e
            ! nle...number of mid-depth levels at elem
            nle          = nlevels(elem)-1
            zbar_n       = 0.0_WP
            Z_n          = 0.0_WP
            zbar_n(nle+1)= zbar_e_bot(elem)
            Z_n(nle)     = zbar_n(nle+1) + helem(nle,elem)/2.0_WP
            do nlz=nle,2,-1
                zbar_n(nlz) = zbar_n(nlz+1) + helem(nlz,elem)
                Z_n(nlz-1)  = zbar_n(nlz)   + helem(nlz-1,elem)/2.0_WP
            end do
            zbar_n(1) = zbar_n(2) + helem(1,elem)
            
            !_______________________________________________________________________
            ! node indices of elem 
            elnodes = elem2D_nodes(:,elem)
            
            !_______________________________________________________________________
            ! nln...number of mid depth levels at each node that corresponds to 
            ! element elem
            nln=nlevels_nod2D(elnodes)-1
            
            !_______________________________________________________________________
            ! loop over mid-depth levels at element elem
            p_grad=0.0_WP
            do nlz=1,nle
                
                !___________________________________________________________________
                rho_n = 0.0_WP
                ! loop over the three node points that span up element elem
                do ni=1,3
                    ! node...
                    node = elnodes(ni)
                    
                    !_______________________________________________________________
                    ! search for center index of cubic spline interpolation 
                    ! [ center-1, center, center+1, center +2 ] at each node to interpolate 
                    ! elemental mid-depth level value
                    dens_ind=nln(ni)
                    do kk=1,nln(ni) 
                        ! Z_n ... mid depth level at element
                        ! Z_3d_n ... mid depth levels at nodes
                        if (Z_3d_n(kk,node)<=Z_n(nlz)) then
                            dens_ind=kk-1
                            if (kk==1) dens_ind=1
                            exit
                        end if
                    end do
                    
                    !_______________________________________________________________
                    ! prepare cubic spline interpolation 
                    s_ind=(/dens_ind-1,dens_ind,dens_ind+1,dens_ind+2/)
                    
                    ! distinguish between bulg (flag=0), surface (flag=1) and
                    ! bottom (flag=-1) 
                    ! why s_ind(1)==0 ? so compiler doesn't complain of zero index
                    ! also s(1) is ot used when flag=1
                    flag=0
                    if (nlz==1 .or. s_ind(1)==0) then ! surface layer
                        flag = 1
                        s_ind(1)=nlz
                    elseif (nlz==nln(ni) ) then ! bottom layer
                        flag = -1
                        s_ind(4)=nlz+1
                    end if 
                    s_z    = Z_3d_n(s_ind,node)
                    s_dens = density_m_rho0(s_ind,node)
                    s_H    = s_z(3)-s_z(2)
                    aux1   = (s_dens(3)-s_dens(2))/s_H
                    
                    !_______________________________________________________________
                    ! calculate derivatives in a way to get monotonic profile
                    if     (flag==0) then ! subsurface/above bottom case
                        aux2=(s_dens(2)-s_dens(1))/(s_z(2)-s_z(1))
                        s_dup=0.0
                        if(aux1*aux2>0.)  s_dup=2.0*aux1*aux2/(aux1+aux2)
                        aux2=(s_dens(4)-s_dens(3))/(s_z(4)-s_z(3))
                        s_dlo=0.0
                        if(aux1*aux2>0.) s_dlo=2.0*aux1*aux2/(aux1+aux2)
                    elseif (flag==1) then ! surface case
                        aux2=(s_dens(4)-s_dens(3))/(s_z(4)-s_z(3))
                        s_dlo=0.0
                        if(aux1*aux2>0.) s_dlo=2.0*aux1*aux2/(aux1+aux2)
                        s_dup=1.5*aux1-0.5*s_dlo
                    elseif (flag==-1) then! bottom case
                        aux2=(s_dens(2)-s_dens(1))/(s_z(2)-s_z(1))
                        s_dup=0.0
                        if(aux1*aux2>0.)  s_dup=2.0*aux1*aux2/(aux1+aux2)
                        s_dlo=1.5*aux1-0.5*s_dup
                    end if
                    
                    !_______________________________________________________________
                    ! cubic polynomial coefficients
                    a=s_dens(2)
                    b=s_dup
                    c=-(2.0*s_dup+s_dlo)/s_H + 3.0*(s_dens(3)-s_dens(2))/s_H**2
                    d=(s_dup+s_dlo)/s_H**2 - 2.0*(s_dens(3)-s_dens(2))/s_H**3
                    
                    !_______________________________________________________________
                    ! interpolate
                    dz=Z_n(nlz)-s_z(2)
                    rho_n(ni)=a+b*dz+c*dz**2+d*dz**3
                    
                end do ! --> do ni=1,3
                
                !___________________________________________________________________
                ! calculate element wise density gradient
                rhograd_e(1) = sum(gradient_sca(1:3,elem)*rho_n)
                rhograd_e(2) = sum(gradient_sca(4:6,elem)*rho_n)
                
                !___________________________________________________________________
                ! calculate element wise pressure gradient force 
                ! helem ... here because of vertical integral 
                aux             = g*helem(nlz,elem)*rhograd_e/density_0
                
                ! *0.5 because pgf_xy is calculated at mid depth levels but at 
                ! this point p_grad is integrated pressure gradient force until
                ! full depth  layers of previouse depth layer
                pgf_x(nlz,elem) = p_grad(1) + aux(1)*0.5
                pgf_y(nlz,elem) = p_grad(2) + aux(2)*0.5
                
                ! integration to full depth levels
                p_grad          = p_grad    + aux
                
            end do ! --> do nlz=1,nle
        end do ! --> do elem=1, myDim_elem2D
    endif
end subroutine pressure_force
!
!
!
! ===========================================================================
SUBROUTINE densityJM_local(t, s, pz, rho_out)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP, only: par_ex,pe_status
IMPLICIT NONE

  !
  ! - calculates in-situ density as a function of potential temperature
  !   (relative to the surface)
  !   using the Jackett and McDougall equation of state
  !   (Copyright (c) 1992, CSIRO, Australia)
  ! - has been derived from the SPEM subroutine rhocal
  !
  !---------------------------------------------------------------------------

  real(kind=WP), intent(IN)  :: t,s,pz
  real(kind=WP), intent(OUT) :: rho_out                 
  real(kind=WP)              :: rhopot, bulk
  real(kind=WP)              :: bulk_0, bulk_pz, bulk_pz2
  !compute secant bulk modulus

  call densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot)

  bulk = bulk_0 + pz*(bulk_pz + pz*bulk_pz2) 

  rho_out = bulk*rhopot / (bulk + 0.1*pz) - density_0

end subroutine densityJM_local
		
! ===========================================================================
SUBROUTINE densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP, only: par_ex,pe_status
IMPLICIT NONE

  !
  ! - calculates in-situ density as a function of potential temperature
  !   (relative to the surface)
  !   using the Jackett and McDougall equation of state
  !   (Copyright (c) 1992, CSIRO, Australia)
  ! - has been derived from the SPEM subroutine rhocal
  !
  ! Ralph Timmermann, August 2005
  !---------------------------------------------------------------------------
  ! N. Rakowski 2014 the split form
  !---------------------------------------------------------------------------
  real(kind=WP), intent(IN)  :: t,s
  real(kind=WP), intent(OUT) :: bulk_0, bulk_pz, bulk_pz2, rhopot
  real(kind=WP)              :: s_sqrt

  real(kind=WP), parameter   :: a0    = 19092.56,     at   = 209.8925
  real(kind=WP), parameter   :: at2   = -3.041638,    at3  = -1.852732e-3
  real(kind=WP), parameter   :: at4   = -1.361629e-5
  real(kind=WP), parameter   :: as    = 104.4077,     ast  = -6.500517
  real(kind=WP), parameter   :: ast2  = .1553190,     ast3 = 2.326469e-4 
  real(kind=WP), parameter   :: ass   = -5.587545,    asst = 0.7390729 
  real(kind=WP), parameter   :: asst2 = -1.909078e-2
  real(kind=WP), parameter   :: ap    = -4.721788e-1, apt  = -1.028859e-2
  real(kind=WP), parameter   :: apt2  = 2.512549e-4,  apt3 = 5.939910e-7 
  real(kind=WP), parameter   :: aps   = 1.571896e-2,  apst = 2.598241e-4 
  real(kind=WP), parameter   :: apst2 = -7.267926e-6, apss = -2.042967e-3
  real(kind=WP), parameter   :: ap2   = 1.045941e-5,  ap2t = -5.782165e-10 
  real(kind=WP), parameter   :: ap2t2 = 1.296821e-7
  real(kind=WP), parameter   :: ap2s  = -2.595994e-7,ap2st = -1.248266e-9 
  real(kind=WP), parameter   :: ap2st2= -3.508914e-9
  
  real(kind=WP), parameter   :: b0 = 999.842594,    bt  = 6.793952e-2
  real(kind=WP), parameter   :: bt2 = -9.095290e-3, bt3 = 1.001685e-4
  real(kind=WP), parameter   :: bt4 = -1.120083e-6, bt5 = 6.536332e-9
  real(kind=WP), parameter   :: bs = 0.824493,      bst = -4.08990e-3
  real(kind=WP), parameter   :: bst2 = 7.64380e-5,  bst3 = -8.24670e-7		
  real(kind=WP), parameter   :: bst4 = 5.38750e-9	
  real(kind=WP), parameter   :: bss = -5.72466e-3,  bsst = 1.02270e-4
  real(kind=WP), parameter   :: bsst2 = -1.65460e-6,bss2 = 4.8314e-4

  !compute secant bulk modulus

  s_sqrt = sqrt(s)

  bulk_0 =  a0      + t*(at   + t*(at2  + t*(at3 + t*at4)))      &
          + s* (as  + t*(ast  + t*(ast2 + t*ast3))               &
               + s_sqrt*(ass  + t*(asst + t*asst2)))                

  bulk_pz =  ap  + t*(apt  + t*(apt2 + t*apt3))                  &
                  + s*(aps + t*(apst + t*apst2) + s_sqrt*apss)

  bulk_pz2 = ap2 + t*(ap2t + t*ap2t2)		                 &
                + s *(ap2s + t*(ap2st + t*ap2st2))

  rhopot =  b0 + t*(bt + t*(bt2 + t*(bt3  + t*(bt4  + t*bt5))))	 &
               + s*(bs + t*(bst + t*(bst2 + t*(bst3 + t*bst4)))  &
                  + s_sqrt*(bss + t*(bsst + t*bsst2))            &
                       + s* bss2)

end subroutine densityJM_components
! ===================================================================
function ptheta(s,t,p,pr)
  ! Compute local potential temperature at pr
  ! using bryden 1973 polynomial for adiabatic lapse rate
  ! and runge-kutta 4-th order integration algorithm.
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! fofonoff,n.,1977,deep-sea res.,24,489-491
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       reference prs   pr       decibars
  !       potential tmp.  theta    deg celsius
  ! checkvalue: theta= 36.89073 c,s=40 (ipss-78),t=40 deg c,
  ! p=10000 decibars,pr=0 decibars
  !
  ! Coded by ??
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_param, only: WP
  implicit none
  real(kind=WP) 			:: ptheta, s, t, p, pr
  real(kind=WP) 			:: h, xk, q
  real(kind=WP), external	        :: atg

  h = pr - p
  xk = h*atg(s,t,p)
  t = t + 0.5*xk
  q = xk
  p = p + 0.5*h
  xk = h*atg(s,t,p)
  t = t + 0.29289322*(xk-q)
  q = 0.58578644*xk + 0.121320344*q
  xk = h*atg(s,t,p)
  t = t + 1.707106781*(xk-q)
  q = 3.414213562*xk - 4.121320344*q
  p = p + 0.5*h
  xk = h*atg(s,t,p)
  ptheta = t + (xk-2.0*q)/6.0
  return
end function ptheta
!
!-----------------------------------------------------------------
!
function atg(s,t,p)
  ! adiabatic temperature gradient deg c per decibar
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       adiabatic       atg      deg. c/decibar
  ! checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
  ! t=40 deg c,p0=10000 decibars
  !
  ! Coded by ??
  ! Reviewed by ??
  !--------------------------------------------------------

  use o_param, only: WP
  implicit none
  real(kind=WP)  atg, s, t, p, ds

  ds = s - 35.0
  atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p   &
       +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t        &
       +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p             &
       +(-4.2393e-8*t+1.8932e-6)*ds                          &
       +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5

  return
end function atg
!
!----------------------------------------------------------------------------
!
subroutine sw_alpha_beta(TF1,SF1)
  ! DESCRIPTION:
  !   A function to calculate the thermal expansion coefficient
  !   and saline contraction coefficient. (elementwise)
  !
  ! INPUT:
  !   tracer(:,2) = salinity              [psu      (PSS-78)]
  !   tracer(:,1) = potential temperature [degree C (ITS-90)]
  !   z           = pressure (or -depth)  [db]
  !
  ! OUTPUT:
  !   sw_alpha = Thermal expansion coeff (alpha) [degree_C.^-1]
  !   sw_beta  = Saline contraction coeff (beta) [psu.^-1]
  !
  ! Qiang Wang, 25,11,2004
  !
  ! REFERENCE:
  !    McDougall, T.J. 1987.  Neutral Surfaces
  !    Journal of Physical Oceanography, vol 17, 1950-1964,
  !-----------------------------------------------------------------
  ! CHECK VALUE:
  !    sw_beta=0.72088e-3 psu^-1 @ S=40.0psu, ptmp=10.0C (ITS-90), p=4000db
  !    a_over_b=0.34765 psu*C^-1 @ S=40.0psu, ptmp=10.0C, p=4000db
  !-----------------------------------------------------------------
  use o_mesh
  use o_arrays
  use g_parsup
  use o_param
  implicit none
  !
  integer        :: n, nz
  real(kind=WP)  :: t1,t1_2,t1_3,t1_4,p1,p1_2,p1_3,s1,s35,s35_2 
  real(kind=WP)  :: a_over_b    
  real(kind=WP)  :: TF1(nl-1, myDim_nod2D+eDim_nod2D),SF1(nl-1, myDim_nod2D+eDim_nod2D)

  do n = 1,myDim_nod2d
     do nz=1, nlevels_nod2D(n)-1
     
     t1 = TF1(nz,n)*1.00024_WP
     s1 = SF1(nz,n)
     p1 = abs(Z(nz)) 
     
     t1_2 = t1*t1
     t1_3 = t1_2*t1
     t1_4 = t1_3*t1
     p1_2 = p1*p1
     p1_3 = p1_2*p1
     s35 = s1-35.0_WP
     s35_2 = s35*s35

     ! calculate beta
     sw_beta(nz,n) = 0.785567e-3 - 0.301985e-5*t1 &
          + 0.555579e-7*t1_2 - 0.415613e-9*t1_3 &
          + s35*(-0.356603e-6 + 0.788212e-8*t1 &
          + 0.408195e-10*p1 - 0.602281e-15*p1_2) &
          + s35_2*(0.515032e-8) & 
          + p1*(-0.121555e-7 + 0.192867e-9*t1 - 0.213127e-11*t1_2) &
          + p1_2*(0.176621e-12 - 0.175379e-14*t1) &
          + p1_3*(0.121551e-17)

     ! calculate the thermal expansion / saline contraction ratio
     a_over_b = 0.665157e-1 + 0.170907e-1*t1 &
          - 0.203814e-3*t1_2 + 0.298357e-5*t1_3 &
          - 0.255019e-7*t1_4 &
          + s35*(0.378110e-2 - 0.846960e-4*t1 &
          - 0.164759e-6*p1 - 0.251520e-11*p1_2) &
          + s35_2*(-0.678662e-5) &
          + p1*(0.380374e-4 - 0.933746e-6*t1 + 0.791325e-8*t1_2) &
          + p1_2*t1_2*(0.512857e-12) &
          - p1_3*(0.302285e-13)

     ! calculate alpha
     sw_alpha(nz,n) = a_over_b*sw_beta(nz,n)
   end do
 end do     
end subroutine sw_alpha_beta
!
!----------------------------------------------------------------------------
!
subroutine compute_sigma_xy(TF1,SF1)
  !--------------------------------------------------------------------
  ! DESCRIPTION:
  !   computes density gradient
  !
  ! INPUT:
  !   SF          = salinity              [psu      (PSS-78)]
  !   TF          = potential temperature [degree C (ITS-90)]
  !
  ! OUTPUT:
  ! based on thermal expansion and saline contraction coefficients
  ! computes density gradient sigma_xy
  !-------------------------------------------------------------------
  use o_mesh
  use o_param
  use o_arrays
  use g_parsup
  use g_comm_auto
  implicit none
  !
  real(kind=WP), intent(IN)   :: TF1(nl-1, myDim_nod2D+eDim_nod2D), SF1(nl-1, myDim_nod2D+eDim_nod2D)
  real(kind=WP)               :: tx(nl-1), ty(nl-1), sx(nl-1), sy(nl-1), vol(nl-1), testino(2)
  integer                     :: n, nz, elnodes(3),el, k, nl1
  
  !
  DO n=1, myDim_nod2D
        nl1 = nlevels_nod2D(n)-1
        vol(1:nl1) = 0.0_WP
	tx(1:nl1)  = 0.0_WP
	ty(1:nl1)  = 0.0_WP
	sx(1:nl1)  = 0.0_WP
	sy(1:nl1)  = 0.0_WP
        DO k=1, nod_in_elem2D_num(n)
           el=nod_in_elem2D(k, n)
           
           DO nz=1, nlevels(el)-1 

              vol(nz) = vol(nz)+elem_area(el)

              !NR  writing the sum over elem2D_nodes explicitly helps the compiler to vectorize the nz-loop

              tx(nz) = tx(nz)+(gradient_sca(1,el)*TF1(nz,elem2D_nodes(1,el)) &
                             + gradient_sca(2,el)*TF1(nz,elem2D_nodes(2,el)) &
                             + gradient_sca(3,el)*TF1(nz,elem2D_nodes(3,el)))*elem_area(el)

              ty(nz) = ty(nz)+(gradient_sca(4,el)*TF1(nz,elem2D_nodes(1,el)) &
                             + gradient_sca(5,el)*TF1(nz,elem2D_nodes(2,el)) &
                             + gradient_sca(6,el)*TF1(nz,elem2D_nodes(3,el)))*elem_area(el)

              sx(nz) = sx(nz)+(gradient_sca(1,el)*SF1(nz,elem2D_nodes(1,el)) &
                             + gradient_sca(2,el)*SF1(nz,elem2D_nodes(2,el)) &
                             + gradient_sca(3,el)*SF1(nz,elem2D_nodes(3,el)))*elem_area(el)

              sy(nz) = sy(nz)+(gradient_sca(4,el)*SF1(nz,elem2D_nodes(1,el)) &
                             + gradient_sca(5,el)*SF1(nz,elem2D_nodes(2,el)) &
                             + gradient_sca(6,el)*SF1(nz,elem2D_nodes(3,el)))*elem_area(el)
           END DO
        enddo
	  sigma_xy(1,1:nl1,n) = (-sw_alpha(1:nl1,n)*tx(1:nl1)+sw_beta(1:nl1,n)*sx(1:nl1))/vol(1:nl1)*density_0
	  sigma_xy(2,1:nl1,n) = (-sw_alpha(1:nl1,n)*ty(1:nl1)+sw_beta(1:nl1,n)*sy(1:nl1))/vol(1:nl1)*density_0
  END DO 

  call exchange_nod(sigma_xy)
end subroutine compute_sigma_xy
!===============================================================================
subroutine compute_neutral_slope
	use o_ARRAYS
	use g_PARSUP
	use o_MESH
	USE o_param
	use g_config
        use g_comm_auto
	IMPLICIT NONE
	real(kind=WP)   :: deltaX1,deltaY1,deltaX2,deltaY2
	integer         :: edge
	integer         :: n,nz,nl1,el(2),elnodes(3),enodes(2)
	real(kind=WP)   :: c, ro_z_inv,eps,S_cr,S_d

	!if sigma_xy is not computed
	eps=5.0e-6
	S_cr=1.0e-2
	S_d=1.0e-3
        slope_tapered=0.
	do n=1, myDim_nod2D
                nl1=nlevels_nod2d(n)-1
		do nz = 2,nl1
			ro_z_inv=2._WP*g/density_0/max(bvfreq(nz,n)+bvfreq(nz+1,n), eps**2) !without minus, because neutral slope S=-(nabla\rho)/(d\rho/dz)
			neutral_slope(1,nz,n)=sigma_xy(1,nz,n)*ro_z_inv
			neutral_slope(2,nz,n)=sigma_xy(2,nz,n)*ro_z_inv
			neutral_slope(3,nz,n)=sqrt(neutral_slope(1,nz,n)**2+neutral_slope(2,nz,n)**2)
			!tapering
                        c=1.0_WP
			c=0.5*(1.0_WP + tanh((S_cr - neutral_slope(3,nz,n))/S_d))
                        if ((bvfreq(nz,n) <= 0.0_WP) .or. (bvfreq(nz+1,n) <= 0.0_WP)) c=0.0_WP
			slope_tapered(:,nz,n)=neutral_slope(:,nz,n)*c
!                       slope_tapered(:,nl1-1:nl1,n)=0.
!                       slope_tapered(:,1:2,n)      =0.
		enddo
	enddo

        call exchange_nod(neutral_slope)
        call exchange_nod(slope_tapered)
end subroutine compute_neutral_slope
!===============================================================================
!converts insitu temperature to a potential one
!               tr_arr(:,:,1) will be modified!
subroutine insitu2pot
  use o_mesh
  use o_param
  use o_arrays
  use g_config
  use g_PARSUP
  implicit none
  real(kind=WP), external     :: ptheta
  real(kind=WP)               :: pp, pr, tt, ss
  integer                     :: n, nz
  ! Convert in situ temperature into potential temperature
  pr=0.0_WP
  do n=1,myDim_nod2d+eDim_nod2D
     do nz=1, nlevels_nod2D(n)-1    
        tt=tr_arr(nz,n,1)
        ss=tr_arr(nz,n,2)
        pp=abs(Z(nz))
        tr_arr(nz,n,1)=ptheta(ss, tt, pp, pr)
     end do	
  end do
end subroutine insitu2pot
