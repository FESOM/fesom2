module info_module
! synopsis: query information from FESOM

  implicit none  
  public info
  private

  type :: info_type
  contains
    procedure, nopass :: print_definitions
  end type
  type(info_type) info

contains

  ! this is a list of preprocessor definitions from the FESOM Fortran source files
  ! it will probably become outdated at some point and should be reviewed
  ! the result will reflect the status of definitions as they are set when *this file* had been compiled
  subroutine print_definitions()
#ifdef __icepack
      print '(g0)', '__icepack is ON'
#else
      print '(g0)', '__icepack is OFF'
#endif  
#ifdef __oasis
      print '(g0)', '__oasis is ON'
#else
      print '(g0)', '__oasis is OFF'
#endif  
#ifdef __oifs
      print '(g0)', '__oifs is ON'
#else
      print '(g0)', '__oifs is OFF'
#endif  
#ifdef DEBUG
      print '(g0)', 'DEBUG is ON'
#else
      print '(g0)', 'DEBUG is OFF'
#endif  
#ifdef DISABLE_MULTITHREADING
      print '(g0)', 'DISABLE_MULTITHREADING is ON'
#else
      print '(g0)', 'DISABLE_MULTITHREADING is OFF'
#endif  
#ifdef false
      print '(g0)', 'false is ON'
#else
      print '(g0)', 'false is OFF'
#endif  
#ifdef FVOM_INIT
      print '(g0)', 'FVOM_INIT is ON'
#else
      print '(g0)', 'FVOM_INIT is OFF'
#endif  
#ifdef oifs
      print '(g0)', 'oifs is ON'
#else
      print '(g0)', 'oifs is OFF'
#endif  
#ifdef OMP_MAX_THREADS
      print '(g0)', 'OMP_MAX_THREADS is ON'
#else
      print '(g0)', 'OMP_MAX_THREADS is OFF'
#endif  
#ifdef PARMS
      print '(g0)', 'PARMS is ON'
#else
      print '(g0)', 'PARMS is OFF'
#endif  
#ifdef PETSC
      print '(g0)', 'PETSC is ON'
#else
      print '(g0)', 'PETSC is OFF'
#endif  
#ifdef use_cavity
      print '(g0)', 'use_cavity is ON'
#else
      print '(g0)', 'use_cavity is OFF'
#endif  
#ifdef use_fullfreesurf
      print '(g0)', 'use_fullfreesurf is ON'
#else
      print '(g0)', 'use_fullfreesurf is OFF'
#endif  
#ifdef VERBOSE
      print '(g0)', 'VERBOSE is ON'
#else
      print '(g0)', 'VERBOSE is OFF'
#endif  
  end subroutine

end module
