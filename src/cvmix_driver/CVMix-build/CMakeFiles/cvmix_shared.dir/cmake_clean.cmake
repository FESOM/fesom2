file(REMOVE_RECURSE
  "libcvmix.pdb"
  "libcvmix.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/cvmix_shared.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
