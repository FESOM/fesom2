file(REMOVE_RECURSE
  "libcvmix.a"
  "libcvmix.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/cvmix_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
