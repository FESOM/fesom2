!---------------------------------------------------------------------------
!Implementation of Gent & McWiliams parameterization after R. Ferrari et al., 2010
!Contains:
!  fer_solve_Gamma
!  fer_gamma2vel
!  fer_compute_C_K ! this subroutine shall be a subject of future tuning (with respect to fer_k)
!===========================================================================
subroutine fer_solve_Gamma
	USE o_MESH
	USE o_PARAM
	USE o_ARRAYS, ONLY: sigma_xy, fer_gamma, bvfreq, fer_c, fer_K, zbar_n, Z_n, hnode_new, zbar_n_bot
	USE g_PARSUP
	USE g_CONFIG
	use g_comm_auto
	IMPLICIT NONE
	
	integer                         :: nz, n, nzmax
	real(kind=WP)                          :: zinv1,zinv2, zinv, m, r
	real(kind=WP)                          :: a(nl), b(nl), c(nl)
	real(kind=WP)                          :: cp(nl), tp(2,nl)
	real(kind=WP), dimension(:,:), pointer :: tr
	
	DO n=1,myDim_nod2D
		tr=>fer_gamma(:,:,n)
! 		!_____________________________________________________________________
! 		! minimum number of levels below elements containing node n
! 		nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
		
		!_______________________________________________________________________
		! compute Z_n(:) and zbar_n(:) for the current ALE step
		nzmax=nlevels_nod2D(n)
		zbar_n=0.0_WP
		Z_n   =0.0_WP
		zbar_n(nzmax)=zbar_n_bot(n)                                 ! depth of the deepest level with respect to partial cell
		Z_n(nzmax-1) =zbar_n(nzmax) + hnode_new(nzmax-1,n)/2.0_WP ! depth of the deepest layer
		do nz=nzmax-1,2,-1
			zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)            ! depth of the nz level
			Z_n(nz-1)  = zbar_n(nz)   + hnode_new(nz-1,n)/2.0_WP   ! depth of the nz layer
		end do
		zbar_n(1) = zbar_n(2) + hnode_new(1,n)                    ! surface level height (/depth)
		
		!_____________________________________________________________________
		! minimum number of levels below elements containing node n
		nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
		
		! The first row
		c(1)=0.0_WP
		a(1)=0.0_WP
		b(1)=1.0_WP
		
		zinv2=1.0_WP/(zbar_n(1)-zbar_n(2))
		DO nz=2, nzmax-1
			zinv1=zinv2
			zinv2=1.0_WP/(zbar_n(nz)-zbar_n(nz+1))
			zinv =1.0_WP/(Z_n(nz-1)-Z_n(nz))
			a(nz)= fer_c(n)*zinv1*zinv
			c(nz)= fer_c(n)*zinv2*zinv
			b(nz)=-a(nz)-c(nz)-max(bvfreq(nz,n), 1.e-8)
		END DO
		! The last row
		nz=nzmax
		c(nz)=0.0_WP
		a(nz)=0.0_WP
		b(nz)=1.0_WP
		! ===========================================
		! The rhs:
		tr(:, 1)=0.
		tr(:, nzmax)=0.
          DO nz=2, nzmax-1
             r=g/density_0
             tr(1, nz)=r*0.5_WP*sum(sigma_xy(1,nz-1:nz,n))*fer_K(nz, n)
             tr(2, nz)=r*0.5_WP*sum(sigma_xy(2,nz-1:nz,n))*fer_K(nz, n)
          END DO
         ! =============================================
          ! The sweep algorithm
          ! initialize c-prime and s,t-prime
          cp(1) = c(1)/b(1)
          tp(:,1) = tr(:,1)/b(1)
! solve for vectors c-prime and t, s-prime
          DO nz = 2, nzmax
           m = b(nz)-cp(nz-1)*a(nz)
           cp(nz) = c(nz)/m
           tp(:,nz) = (tr(:,nz)-tp(:,nz-1)*a(nz))/m
          END DO
! initialize x
		tr(:,nzmax) = tp(:,nzmax)
		! solve for x from the vectors c-prime and d-prime
		do nz = nzmax-1, 1, -1
			tr(:,nz) = tp(:,nz)-cp(nz)*tr(:,nz+1)
		end do
	END DO   !!! cycle over nodes
	
	call exchange_nod(fer_gamma)
	
END subroutine fer_solve_Gamma
!====================================================================
subroutine fer_gamma2vel
  USE o_MESH
  USE o_PARAM
  USE o_ARRAYS, ONLY: fer_gamma, fer_uv, helem
  USE g_PARSUP
  USE g_CONFIG
  use g_comm_auto
  IMPLICIT NONE

   integer                         :: nz, nzmax, el, elnod(3)
   real(kind=WP)                          :: zinv
   real(kind=WP)                          :: onethird=1._WP/3._WP

   DO el=1, myDim_elem2D
      elnod=elem2D_nodes(:,el)
      ! max. number of levels at element el
      nzmax=nlevels(el)
      DO nz=1, nzmax-1
         zinv=onethird/helem(nz,el)
         fer_uv(1,nz,el)=sum(fer_gamma(1,nz,elnod)-fer_gamma(1,nz+1,elnod))*zinv
         fer_uv(2,nz,el)=sum(fer_gamma(2,nz,elnod)-fer_gamma(2,nz+1,elnod))*zinv
      END DO
   END DO
   call exchange_elem(fer_uv)
end subroutine fer_gamma2vel
!====================================================================
subroutine init_Redi_GM!fer_compute_C_K_Redi
  USE o_MESH
  USE o_PARAM
  USE o_ARRAYS, ONLY: fer_c, fer_k, Ki, bvfreq, MLD1_ind, neutral_slope, coriolis_node, hnode_new
  USE g_PARSUP
  USE g_CONFIG
  use g_comm_auto
  IMPLICIT NONE

   integer                         :: n, nz, nzmax
   real(kind=WP)                          :: reso, c1, rosb, scaling, rr_ratio
   real(kind=WP)                          :: x0=1.5, sigma=.15 ! Fermi function parameters to cut off GM where Rossby radius is resolved
   real(kind=WP)                          :: c_min=0.5, f_min=1.e-6, r_max=200000.
   real(kind=WP)                          :: zscaling(nl)
   real(kind=WP)                          :: bvref

! fill arrays for 3D Redi and GM coefficients: F1(xy)*F2(z)
!******************************* F1(x,y) ***********************************************************
   DO n=1, myDim_nod2D
      reso=mesh_resolution(n)
      if (Fer_GM) then
         c1=0._wp
         nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(n), n)), 1)
         DO nz=1, nzmax-1
            c1=c1+hnode_new(nz,n)*(sqrt(max(bvfreq(nz,n), 0._WP))+sqrt(max(bvfreq(nz+1,n), 0._WP)))/2.
         END DO
         c1=max(c_min, c1/pi) !ca. first baroclinic gravity wave speed limited from below by c_min
         scaling=1._WP
         !Cutoff K_GM depending on (Resolution/Rossby radius) ratio
         if (scaling_Rossby) then
            rosb=min(c1/max(abs(coriolis_node(n)), f_min), r_max)
            rr_ratio=min(reso/rosb, 5._WP)
            scaling=1._WP/(1._WP+exp(-(rr_ratio-x0)/sigma))
         end if
         !Scale K_GM with resolution (referenced to 100,000m)
         if (scaling_resolution) then
            scaling=scaling*(reso/100000._WP)**2 !put to repo
         end if
         if (reso < 40000.0_WP) then
            scaling=scaling*max((reso/10000.0_WP-3.0_WP), 0._WP) !no GM below 30km resolution
         end if
         fer_k(1,n)=min(K_GM*scaling, 2000.0_WP) !put to repo
         fer_k(1,n)=max(K_GM*scaling, 2.0_WP)    !put to repo
         fer_c(n)=c1*c1                          !put to repo
      end if
      !note, Redi can be used without GM and vise versa!
      ! if both are used it will be reset below
      if (Redi) then
         Ki(1,n)=K_hor*(reso/100000.0_WP)**2
      end if
   END DO
 
   !Like in FESOM 1.4 we make Redi equal GM
   if (Redi .and. Fer_GM) then
      Ki(1,:)=fer_k(1,:)
   end if
!******************************* F2(z) (e.g. Ferreira et al., 2005) *********************************
   DO n=1,myDim_nod2D
      nzmax=nlevels_nod2D(n)
      ! Allpy vertical scaling after Ferreira et al.(2005)
      if (scaling_Ferreira) then
         bvref=max(bvfreq(MLD1_ind(n)+1, n), 1.e-6_WP)
         DO nz=1, nzmax
            zscaling(nz)=max(bvfreq(nz, n)/bvref, 0.2_WP)
            zscaling(nz)=min(zscaling(nz), 1.0_WP)
         END DO
      else
         zscaling=1.0_WP
      end if
      ! Switch off GM and Redi within a BL in NH (a strategy following FESOM 1.4)
      if (scaling_FESOM14) then
         !zscaling(1:MLD1_ind(n)+1)=0.0_WP
         DO nz=1, nzmax
            if (neutral_slope(3, min(nz, nl-1), n) > 5.e-3) zscaling(nz)=0.0_WP
         END DO
      end if
      if (Fer_GM) then
         fer_k(:nzmax,n)=fer_k(1,n)*zscaling(1:nzmax)
      end if
      if (Redi) then
         Ki(1:nzmax-1,n)  =Ki(1,n)*0.5*(zscaling(1:nzmax-1)+zscaling(2:nzmax))
      end if
   END DO

   if (Fer_GM) call exchange_nod(fer_c)
   if (Fer_GM) call exchange_nod(fer_k)
   if (Redi)   call exchange_nod(Ki)
end subroutine init_Redi_GM
!====================================================================
