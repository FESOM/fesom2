&diag_list
ldiag_solver     =.false.
lcurt_stress_surf=.false.
ldiag_curl_vel3  =.false.
ldiag_energy     =.false.
ldiag_salt3D     =.false.
ldiag_dMOC       =.false.
ldiag_DVD        =.false.
ldiag_forc       =.false.
/

&nml_listsize
io_listsize=100 !number of streams to allocate. shallbe large or equal to the number of streams in &nml_list
/

! for sea ice related variables use_ice should be true, otherewise there will be no output
! for 'curl_surf' to work lcurt_stress_surf must be .true. otherwise no output
! for 'fer_C', 'bolus_u', 'bolus_v', 'bolus_w', 'fer_K' to work Fer_GM must be .true. otherwise no output
! 'otracers' - all other tracers if applicable
! for 'dMOC' to work ldiag_dMOC must be .true. otherwise no output
&nml_list
io_list =  'sst       ',1, 's', 4,
           'sss       ',1, 's', 4,
           'a_ice     ',1, 's', 4,
           'temp      ',1, 's', 4,
           'salt      ',1, 's', 8,
           'u         ',1, 's', 4,
/
&io_ctl
!lzarr_io=.true.
lzarr_io=.false.
/
