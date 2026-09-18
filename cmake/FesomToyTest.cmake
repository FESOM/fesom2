#===============================================================================
# FesomToyTest.cmake - CTest support for the idealized (toy_ocean) configurations
#===============================================================================
#
# Idealized configurations run with toy_ocean=.true. and supply their own
# forcing, so they read a mesh and nothing else: no atmospheric forcing and no
# climatology. Their meshes are in the repository under test/meshes/, so these
# tests point MeshPath straight at that tree instead of using a mesh name under
# tests/data/MESHES -- nothing is copied or linked.
#
# Requires FesomTesting.cmake, which tests/CMakeLists.txt already includes.

include_guard(GLOBAL)

# Set KEY=VALUE in namelist group GROUP of NAMELIST_FILE, inserting the key
# after the group header when it is not already present. VALUE_PATTERN is the
# regex for an existing value (a logical, a number, a quoted string).
#
# Insertion is needed because keys such as toy_ocean and state_equation are
# absent from the stock namelists; a plain REGEX REPLACE would silently do
# nothing. A missing group is an error for the same reason.
function(fesom_nml_set NAMELIST_FILE GROUP KEY VALUE VALUE_PATTERN)
    file(READ "${NAMELIST_FILE}" _content)

    if(NOT _content MATCHES "&${GROUP}")
        message(FATAL_ERROR "fesom_nml_set: no &${GROUP} group in ${NAMELIST_FILE}")
    endif()

    # Anchored on a non-word character and requiring '=' so a key cannot match
    # as the prefix of a longer one (use_ice vs use_icebergs).
    if(_content MATCHES "[^A-Za-z0-9_]${KEY}[ \t]*=")
        string(REGEX REPLACE "([^A-Za-z0-9_])${KEY}[ \t]*=[ \t]*${VALUE_PATTERN}"
               "\\1${KEY}=${VALUE}" _content "${_content}")
    else()
        string(REGEX REPLACE "(&${GROUP}[^\n]*\n)" "\\1${KEY}=${VALUE}\n"
               _content "${_content}")
    endif()

    file(WRITE "${NAMELIST_FILE}" "${_content}")
endfunction()

# add_fesom_toy_test(<name> WHICH_TOY <soufflet|neverworld2> MESH_DIR <path>
#                    CYCLIC_LENGTH <degrees> [NP <n>] [TIMEOUT <s>]
#                    [STEP_PER_DAY <n>])
#
# Namelist values follow setups/<config>/setup.yml, which is what CI runs. The
# config/namelist.*.toy_* templates are not used: they set a different
# step_per_day, run_length, which_ALE and cartesian, and there is no
# namelist.dyn.toy_soufflet at all.
function(add_fesom_toy_test TEST_NAME)
    set(oneValueArgs WHICH_TOY NP TIMEOUT MESH_DIR CYCLIC_LENGTH STEP_PER_DAY)
    cmake_parse_arguments(TOY "" "${oneValueArgs}" "" ${ARGN})

    if(NOT TOY_NP)
        set(TOY_NP 8)
    endif()
    if(NOT TOY_TIMEOUT)
        set(TOY_TIMEOUT 600)
    endif()
    if(NOT TOY_STEP_PER_DAY)
        set(TOY_STEP_PER_DAY 36)
    endif()

    # Skip before registering: a test that is already added cannot be removed,
    # and a missing mesh would otherwise surface as a Fortran read error deep in
    # startup.
    if(NOT EXISTS "${TOY_MESH_DIR}/nod2d.out" OR NOT EXISTS "${TOY_MESH_DIR}/dist_${TOY_NP}")
        message(STATUS "Toy test ${TEST_NAME} skipped: no mesh or no dist_${TOY_NP} in ${TOY_MESH_DIR}")
        return()
    endif()

    add_fesom_test_with_options(${TEST_NAME}
        "pi" "${TOY_STEP_PER_DAY}" "1" "d" "1" "d" "10" ".false." ".false."
        MPI_TEST
        NP ${TOY_NP}
        TIMEOUT ${TOY_TIMEOUT}
        LABEL toy
        MIX_SCHEME "PP"
    )

    # The namelists are written at configure time and read by the generated run
    # script at run time, so they can be adjusted here without touching the
    # shared test generator.
    set(_dir "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")
    set(_logical "\\.[a-zA-Z]+\\.")
    set(_number "[-0-9.eE+]+")

    # Use the mesh in place; the positional mesh name above is a placeholder.
    file(READ "${_dir}/namelist.config" _cfg)
    string(REGEX REPLACE "([^A-Za-z0-9_])MeshPath[ \t]*=[ \t]*'[^']*'"
           "\\1MeshPath='${TOY_MESH_DIR}/'" _cfg "${_cfg}")
    file(WRITE "${_dir}/namelist.config" "${_cfg}")

    fesom_nml_set("${_dir}/namelist.config" run_config use_ice       ".false."            "${_logical}")
    fesom_nml_set("${_dir}/namelist.config" run_config use_sw_pene   ".false."            "${_logical}")
    fesom_nml_set("${_dir}/namelist.config" run_config toy_ocean     ".true."             "${_logical}")
    fesom_nml_set("${_dir}/namelist.config" run_config which_toy     "'${TOY_WHICH_TOY}'" "'[^']*'")
    fesom_nml_set("${_dir}/namelist.config" geometry   cyclic_length "${TOY_CYCLIC_LENGTH}" "${_number}")
    fesom_nml_set("${_dir}/namelist.config" geometry   rotated_grid  ".false."            "${_logical}")

    fesom_nml_set("${_dir}/namelist.oce" oce_dyn state_equation "0"       "${_number}")
    fesom_nml_set("${_dir}/namelist.oce" oce_dyn Fer_GM         ".false." "${_logical}")
    fesom_nml_set("${_dir}/namelist.oce" oce_dyn Redi           ".false." "${_logical}")

    fesom_nml_set("${_dir}/namelist.tra" tracer_phys use_momix          ".false." "${_logical}")
    fesom_nml_set("${_dir}/namelist.tra" tracer_phys K_hor              "10."     "${_number}")
    fesom_nml_set("${_dir}/namelist.tra" tracer_phys surf_relax_S       "0.0"     "${_number}")
    fesom_nml_set("${_dir}/namelist.tra" tracer_phys balance_salt_water ".false." "${_logical}")

    fesom_nml_set("${_dir}/namelist.dyn" dynamics_general use_wsplit ".false." "${_logical}")
endfunction()
