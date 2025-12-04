&diag_list
ldiag_solver     =.false.
lcurt_stress_surf=.false.
ldiag_curl_vel3  =.false.
ldiag_Ri         =.false.
ldiag_turbflux   =.false.
ldiag_salt3D     =.false.
ldiag_dMOC       =.false.
ldiag_DVD        =.false.
ldiag_forc       =.false.
ldiag_extflds    =.false.
/

&nml_general
io_listsize    =150 !number of streams to allocate. shallbe large or equal to the number of streams in &nml_list
vec_autorotate =.false.
/

! for sea ice related variables use_ice should be true, otherewise there will be no output
! for 'curl_surf' to work lcurt_stress_surf must be .true. otherwise no output
! for 'fer_C', 'bolus_u', 'bolus_v', 'bolus_w', 'fer_K' to work Fer_GM must be .true. otherwise no output
! 'otracers' - all other tracers if applicable
! for 'dMOC' to work ldiag_dMOC must be .true. otherwise no output
&nml_list
io_list =  'sst       ',1, 'm', 4,
           'sss       ',1, 'm', 4,
    	   'ssh       ',1, 'm', 4,
           'uice      ',1, 'm', 4,
           'vice      ',1, 'm', 4,
           'a_ice     ',1, 'm', 4,
           'm_ice     ',1, 'm', 4,
           'm_snow    ',1, 'm', 4,
           'MLD1      ',1, 'm', 4,
           'MLD2      ',1, 'm', 4,
           'MLD3      ',1, 'm', 4,
           'tx_sur    ',1, 'm', 4,
           'ty_sur    ',1, 'm', 4,
           'temp      ',1, 'm', 4,
           'salt      ',1, 'm', 8,
           'otracers  ',1, 'm', 4,
           'N2        ',1, 'y', 4,
           'Kv        ',1, 'y', 4,
           'u         ',1, 'm', 4,
           'v         ',1, 'm', 4,
           'unod      ',1, 'm', 4,
           'vnod      ',1, 'm', 4,
           'w         ',1, 'm', 4,
           'Av        ',1, 'y', 4,
           'bolus_u   ',1, 'y', 4,
           'bolus_v   ',1, 'y', 4,
           'bolus_w   ',1, 'y', 4,
           'dpCO2s    ',1, 'm', 4,
           'pCO2s     ',1, 'm', 4,
           'CO2f      ',1, 'm', 4,
           'Hp        ',1, 'm', 4,
           'aFe       ',1, 'm', 4,
           'aN        ',1, 'm', 4,
           'denb      ',1, 'm', 4,
           'benN      ',1, 'm', 4,
           'benC      ',1, 'm', 4,
           'benSi     ',1, 'm', 4,
           'benCalc   ',1, 'm', 4,
           'Chldegd   ',1, 'm', 4,
           'Chldegn   ',1, 'm', 4,
           'NNAd      ',1, 'm', 4,
           'NNAn      ',1, 'm', 4,
           'GPPd      ',1, 'm', 4,
           'GPPn      ',1, 'm', 4,
           'NPPd      ',1, 'm', 4,
           'NPPn      ',1, 'm', 4,
           'PAR       ',1, 'm', 4,
           'respmeso',1, 'm', 4,
           'respmacro',1, 'm', 4,
           'respmicro',1, 'm', 4,
           'calcdiss',1, 'm', 4,
           'calcif',1, 'm', 4,
           'aggn',1, 'm', 4,
           'aggd',1, 'm', 4,
           'docexn',1, 'm', 4,
           'docexd',1, 'm', 4,
           'respn',1, 'm', 4,
           'respd',1, 'm', 4,
/
