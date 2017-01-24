module g_config
  use o_param
  implicit none
  save

  ! *** Modelname ***
  character(5)             	:: runid='fesom'                ! a model/setup name
  integer                  	:: year_start=1948, year_end=1953, snap_per_year=12
  CHARACTER*100	:: datapath ='/home/h/hbkdsido/fvom/results/core2_sw_chl0.1/'
  CHARACTER*100	:: outpath  ='/home/h/hbkdsido/fvom/results/'
namelist /config/ runid, datapath, outpath, year_start, year_end, snap_per_year

  logical                  	:: do_UVnorm=.true.
  logical                  	:: do_MOC=.true.
  logical                  	:: do_UVcurl=.true.
namelist /todo/ do_UVcurl, do_UVnorm, do_MOC

  ! *** model geometry
  logical                  	:: cartesian=.false.
  real(kind=8)             	:: domain_length=360.    	![degree]
  !
  logical                  	:: rotated_grid=.false. !rotated mesh yes/no
  logical                  	:: rotated_rslt=.false. !rotated result yes/no (if the rotation was forced)
  real(kind=8)             	:: alphaEuler=50. 		![grad] Euler angles, convention:
  real(kind=8)             	:: betaEuler=15.  		![grad] first around z, then around new x,
  real(kind=8)			:: gammaEuler=-90.		![grad] then around new z.
  CHARACTER*100	                :: meshpath ='/gfs2/work/hbkdsido/input/mesh/mesh_CORE2_final/'
namelist /fesom_mesh/ meshpath, rotated_grid, rotated_rslt, alphaEuler, betaEuler, gammaEuler
end module g_config
