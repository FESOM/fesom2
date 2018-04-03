find_path(OASIS_Fortran_INCLUDE_DIRECTORIES mod_oasis.mod HINTS ${TOPLEVEL_DIR}/../oasis3-mct/ecconf/build/lib/psmile.MPI1)
find_library(OASIS_Fortran_LIBRARIES psmile HINTS ${TOPLEVEL_DIR}/../oasis3-mct/ecconf/lib)

find_path(MCT_Fortran_INCLUDE_DIRECTORIES mct_mod.mod HINTS ${TOPLEVEL_DIR}/../oasis3-mct/ecconf/build/lib/mct)
find_library(MCT_Fortran_LIBRARIES mct HINTS ${TOPLEVEL_DIR}/../oasis3-mct/ecconf/build/lib/mctdir/mct)

find_path(MPEU_Fortran_INCLUDE_DIRECTORIES m_mpout.mod HINTS ${TOPLEVEL_DIR}/../oasis3-mct//ecconf/build/lib/mct)
find_library(MPEU_Fortran_LIBRARIES mpeu HINTS ${TOPLEVEL_DIR}/../oasis3-mct/ecconf/build/lib/mctdir/mpeu)

find_path(SCRIP_Fortran_INCLUDE_DIRECTORIES remap_bicubic.mod HINTS ${TOPLEVEL_DIR}/../oasis3-mct/ecconf/build/lib/scrip)
find_library(SCRIP_Fortran_LIBRARIES scrip HINTS ${TOPLEVEL_DIR}/../oasis3-mct/ecconf/lib)
