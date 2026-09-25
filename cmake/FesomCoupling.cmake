# =============================================================================
# Coupling-interface selection
# =============================================================================
#
# FESOM uses exactly one interface to external components. The choice is
# compile-time: it decides which coupler library is linked and which API
# generation the driver code is written against. A single scalar variable
# makes "exactly one" structural rather than a separate validation step.
#
#   standalone  no external component
#   direct      IFS drives FESOM in-process; no coupler library
#   oasis28     OASIS3-MCT 2.8  (prism_* API generation)
#   oasis50     OASIS3-MCT 5.0  (oasis_* API generation)
#   yac         YAC
#
# The external components FESOM is coupled to are configured at run time,
# in config/namelist.cpl.
#
# This file is the single source of truth for the coupling configuration.
# Everything downstream reads the FESOM_COUPLING_* variables it exports, not
# the deprecated booleans below.

set(FESOM_COUPLING_VALUES standalone direct oasis28 oasis50 yac)
string(REPLACE ";" ", " FESOM_COUPLING_VALUES_TEXT "${FESOM_COUPLING_VALUES}")

# Empty default so that "the user asked for nothing" stays distinguishable
# from "the user asked for standalone"; resolved below.
set(FESOM_COUPLING_DOC
    "External-component interface: standalone|direct|oasis28|oasis50|yac")

set(FESOM_COUPLING "" CACHE STRING ${FESOM_COUPLING_DOC})
set_property(CACHE FESOM_COUPLING PROPERTY STRINGS ${FESOM_COUPLING_VALUES})

# -----------------------------------------------------------------------------
# Deprecated boolean switches
# -----------------------------------------------------------------------------
# Superseded by FESOM_COUPLING. Kept working for one release so that existing
# job scripts, machine configurations and downstream bundles do not break.

set(FESOM_COUPLED OFF CACHE BOOL
    "DEPRECATED: use FESOM_COUPLING=oasis28 or oasis50")
set(OIFS_COUPLED OFF CACHE BOOL
    "DEPRECATED: use FESOM_COUPLING=oasis50")
set(USE_YAC OFF CACHE BOOL
    "DEPRECATED: use FESOM_COUPLING=yac")

set(ENABLE_IFS_INTERFACE_DEFAULT OFF)
if(DEFINED BUILD_FESOM_AS_LIBRARY AND NOT DEFINED ENABLE_IFS_INTERFACE)
   message(DEPRECATION
           "BUILD_FESOM_AS_LIBRARY is deprecated, use FESOM_COUPLING=direct")
   set(ENABLE_IFS_INTERFACE_DEFAULT ${BUILD_FESOM_AS_LIBRARY})
endif()
option(ENABLE_IFS_INTERFACE "DEPRECATED: use FESOM_COUPLING=direct"
       ${ENABLE_IFS_INTERFACE_DEFAULT})
if(DEFINED FESOM_ENABLE_IFS_INTERFACE)
   # distinguishes the option in a nested cmake project (bundle)
   set(ENABLE_IFS_INTERFACE ${FESOM_ENABLE_IFS_INTERFACE})
endif()

# -----------------------------------------------------------------------------
# Resolve FESOM_COUPLING
# -----------------------------------------------------------------------------
if(FESOM_COUPLING STREQUAL "")
   set(_fesom_cpl_legacy "")
   set(_fesom_cpl_legacy_vars "")

   if(FESOM_COUPLED AND OIFS_COUPLED)
      set(_fesom_cpl_legacy oasis50)
      set(_fesom_cpl_legacy_vars "FESOM_COUPLED + OIFS_COUPLED")
   elseif(FESOM_COUPLED)
      set(_fesom_cpl_legacy oasis28)
      set(_fesom_cpl_legacy_vars "FESOM_COUPLED")
   elseif(OIFS_COUPLED)
      message(FATAL_ERROR
              "OIFS_COUPLED=ON requires FESOM_COUPLED=ON. "
              "Use -DFESOM_COUPLING=oasis50 instead.")
   endif()

   foreach(_pair "USE_YAC;yac" "ENABLE_IFS_INTERFACE;direct")
      list(GET _pair 0 _var)
      list(GET _pair 1 _val)
      if(${_var})
         if(_fesom_cpl_legacy)
            message(FATAL_ERROR
                    "${_fesom_cpl_legacy_vars} and ${_var} are mutually "
                    "exclusive: FESOM has exactly one coupling interface. "
                    "Use -DFESOM_COUPLING=<one of "
                    "${FESOM_COUPLING_VALUES_TEXT}>.")
         endif()
         set(_fesom_cpl_legacy ${_val})
         set(_fesom_cpl_legacy_vars ${_var})
      endif()
   endforeach()

   if(_fesom_cpl_legacy)
      message(DEPRECATION
              "${_fesom_cpl_legacy_vars} is deprecated; "
              "use -DFESOM_COUPLING=${_fesom_cpl_legacy} instead.")
      set(FESOM_COUPLING ${_fesom_cpl_legacy})
   else()
      set(FESOM_COUPLING standalone)
   endif()

   # Write the resolved value back so that ccmake/cmake-gui show it.
   set(FESOM_COUPLING ${FESOM_COUPLING} CACHE STRING
       ${FESOM_COUPLING_DOC} FORCE)
endif()

if(NOT FESOM_COUPLING IN_LIST FESOM_COUPLING_VALUES)
   message(FATAL_ERROR
           "FESOM_COUPLING must be one of: standalone, direct, oasis28, "
           "oasis50, yac (got '${FESOM_COUPLING}')")
endif()

# -----------------------------------------------------------------------------
# Derived properties of the selected interface
# -----------------------------------------------------------------------------
#   FESOM_COUPLING_ENABLED  not standalone
#   FESOM_COUPLING_COUPLER  a coupler library owns MPI init and the
#                           FESOM communicator (OASIS or YAC)
#   FESOM_COUPLING_IS_OASIS OASIS driver is compiled (either generation)
#   FESOM_CPL_MACROS        the compile definitions to hand the compiler

set(FESOM_COUPLING_ENABLED   OFF)
set(FESOM_COUPLING_COUPLER   OFF)
set(FESOM_COUPLING_IS_OASIS  OFF)
set(FESOM_COUPLING_IS_YAC    OFF)
set(FESOM_COUPLING_IS_DIRECT OFF)
set(FESOM_CPL_MACROS "")

if(FESOM_COUPLING STREQUAL "standalone")
   set(FESOM_CPL_MACROS __standalone)
elseif(FESOM_COUPLING STREQUAL "direct")
   set(FESOM_COUPLING_ENABLED   ON)
   set(FESOM_COUPLING_IS_DIRECT ON)
   set(FESOM_CPL_MACROS __cpl_direct __cpl_enabled)
elseif(FESOM_COUPLING STREQUAL "oasis28")
   set(FESOM_COUPLING_ENABLED  ON)
   set(FESOM_COUPLING_COUPLER  ON)
   set(FESOM_COUPLING_IS_OASIS ON)
   set(FESOM_CPL_MACROS __cpl_oasis28 __cpl_oasis __cpl_coupler __cpl_enabled)
elseif(FESOM_COUPLING STREQUAL "oasis50")
   set(FESOM_COUPLING_ENABLED  ON)
   set(FESOM_COUPLING_COUPLER  ON)
   set(FESOM_COUPLING_IS_OASIS ON)
   set(FESOM_CPL_MACROS __cpl_oasis50 __cpl_oasis __cpl_coupler __cpl_enabled)
elseif(FESOM_COUPLING STREQUAL "yac")
   set(FESOM_COUPLING_ENABLED ON)
   set(FESOM_COUPLING_COUPLER ON)
   set(FESOM_COUPLING_IS_YAC  ON)
   set(FESOM_CPL_MACROS __cpl_yac __cpl_coupler __cpl_enabled)
endif()

# Keep the deprecated booleans consistent with the resolved interface. These
# plain set()s shadow the cache entries for every subdirectory, so code that
# has not been migrated yet still sees the right value. Remove together with
# the last consumer.
set(FESOM_COUPLED ${FESOM_COUPLING_IS_OASIS})
set(USE_YAC ${FESOM_COUPLING_IS_YAC})
set(ENABLE_IFS_INTERFACE ${FESOM_COUPLING_IS_DIRECT})
if(FESOM_COUPLING STREQUAL "oasis50")
   set(OIFS_COUPLED ON)
else()
   set(OIFS_COUPLED OFF)
endif()

message(STATUS "FESOM_COUPLING: ${FESOM_COUPLING}")
