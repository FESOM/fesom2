&nml_listsize
io_listsize=100 !number of streams to allocate. shallbe large or equal to the number of streams in &nml_list
/
&nml_list
io_list =  'sst       ',1, 'm', 4, 
           'sss       ',1, 'm', 4,
           'vve_5     ',1, 'm', 4,
           'uice      ',1, 'm', 4, !use_ice must be .true. otherwise no output
           'vice      ',1, 'm', 4, !use_ice must be .true. otherwise no output
           'a_ice     ',1, 'm', 4, !use_ice must be .true. otherwise no output
           'm_ice     ',1, 'm', 4, !use_ice must be .true. otherwise no output
           'thdgr     ',1, 'm', 4, !use_ice must be .true. otherwise no output
           'thdgrsn   ',1, 'm', 4, !use_ice must be .true. otherwise no output
           'm_snow    ',1, 'm', 4, !use_ice must be .true. otherwise no output
           'MLD1      ',1, 'm', 4,
           'MLD2      ',1, 'm', 4,
           'fh        ',1, 'm', 4,
           'fw        ',1, 'm', 4,
           'atmice_x  ',1, 'm', 4,
           'atmice_y  ',1, 'm', 4,
           'atmoce_x  ',1, 'm', 4,
           'atmoce_y  ',1, 'm', 4,
           'alpha     ',1, 'm', 4,
           'beta      ',1, 'm', 4,
           'runoff    ',1, 'm', 4,
           'evap      ',1, 'm', 4,
           'prec      ',1, 'm', 4,
           'hbl       ',1, 'm', 4, !KPP must be .true. otherwise no output
           'Bo        ',1, 'm', 4, !KPP must be .true. otherwise no output
           'tx_sur    ',1, 'm', 4,
           'ty_sur    ',1, 'm', 4,
           'curl_surf ',1, 'm', 4, !lcurt_stress_surf must be .true. otherwise no output
           'fer_C     ',1, 'm', 4, !Fer_GM must be .true. otherwise no output
           'temp      ',1, 'y', 4,
           'salt      ',1, 'y', 4,
           'otracers  ',1, 'y', 4, !all other tracers if applicable
           'slope_x   ',1, 'y', 4,
           'slope_y   ',1, 'y', 4,
           'slope_z   ',1, 'y', 4,
           'N2        ',1, 'y', 4,
           'Kv        ',1, 'y', 4,
           'u         ',1, 'y', 4,
           'v         ',1, 'y', 4,
           'w         ',1, 'y', 4,
           'Av        ',1, 'y', 4,
           'bolus_u   ',1, 'y', 4, !Fer_GM must be .true. otherwise no output
           'bolus_v   ',1, 'y', 4, !Fer_GM must be .true. otherwise no output
           'bolus_w   ',1, 'y', 4, !Fer_GM must be .true. otherwise no output
           'fer_K     ',1, 'y', 4, !Fer_GM must be .true. otherwise no output
/
