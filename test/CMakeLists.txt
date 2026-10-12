# FESOM3 self-tests. add_fesom_test(<name> <nranks>) is defined in the top-level
# CMakeLists; each builds test/<name>.F90 against libfesom3 and registers a ctest.

# test_util.F90: scaffolding shared by the hybrid-mesh tests (mesh setup, check/finish,
# reductions, hval, set_depths); add_fesom_test links it into every test.
add_library(fesom3_test_util STATIC ${CMAKE_CURRENT_SOURCE_DIR}/test_util.F90)
apply_fesom_compile_flags(fesom3_test_util)
target_link_libraries(fesom3_test_util PUBLIC fesom3)

add_fesom_test(test_params 1)
add_fesom_test(test_types  1)
add_fesom_test(test_partit 1)
add_fesom_test(test_partit 2)
add_fesom_test(test_partit 8)
add_fesom_test(test_halo   1)
add_fesom_test(test_halo   2)
add_fesom_test(test_halo   8)
add_fesom_test(test_dump   1)
add_fesom_test(test_mesh   1)
add_fesom_test(test_vert_spline 1)
add_fesom_test(test_vinv 1)
add_fesom_test(test_vinv 2)
add_fesom_test(test_muscl_onthefly 1)
add_fesom_test(test_muscl_onthefly 2)
add_fesom_test(test_ivertvisc 1)
add_fesom_test(test_ivertvisc 2)
add_fesom_test(test_wsplit 1)
add_fesom_test(test_wsplit 2)
add_fesom_test(test_wimpl_tra 1)
add_fesom_test(test_wimpl_tra 2)
add_fesom_test(test_bottom 1)
add_fesom_test(test_bottom 2)
add_fesom_test(test_bottom 8)
add_fesom_test(test_io_decomp 1)
add_fesom_test(test_io_decomp 2)
add_fesom_test(test_io_decomp 8)
add_fesom_test(test_io_decomp_gather 1)
add_fesom_test(test_io_decomp_gather 2)
add_fesom_test(test_io_posix 1)
add_fesom_test(test_vector_rotate 1)
add_fesom_test(test_io_means 1)

# Hybrid meshes (docs/plans/completed/2026-10-08-hybrid-meshes.md). The writable meshes live
# under /work/ab0995/a270029/meshes. A test whose mesh (or IC file) is absent at configure
# time is registered DISABLED (ctest: "Not Run (Disabled)", reconfigure after adding a
# mesh); a present mesh with a missing dist_N/com file, or any failing file open, FAILS.
set(_hyb_meshes /work/ab0995/a270029/meshes)
set(_mesh_q   ${_hyb_meshes}/nemo025q_cut)
set(_mesh_mix ${_hyb_meshes}/nemo025mix_cut)
set(_mesh_piw ${_hyb_meshes}/pi_w)
set(_mesh_c2w ${_hyb_meshes}/core2_w)
set(_mesh_pi  /home/a/a270088/port2/fesom2/tests/data/MESHES/pi)
set(_icf  /pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc)
# require_files(<test> <file>...): disable <test> unless every file exists
function(require_files test)
    foreach(_f ${ARGN})
        if(NOT EXISTS ${_f})
            set_tests_properties(${test} PROPERTIES DISABLED TRUE)
            message(STATUS "${test}: ${_f} not found, test disabled")
            return()
        endif()
    endforeach()
endfunction()
# add_mesh_test(<test> <target> <mesh dir> <nranks> [env...]): <target> on <nranks> ranks
# with FESOM3_MESH_DIR=<mesh dir> and the extra environment entries
function(add_mesh_test test target mesh np)
    add_test(NAME ${test}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${np} $<TARGET_FILE:${target}>)
    set(_env "FESOM3_MESH_DIR=${mesh}" ${ARGN})
    set_tests_properties(${test} PROPERTIES ENVIRONMENT "${_env}")
    require_files(${test} ${mesh}/nod2d.out)
endfunction()
# add_default_mesh_test(<target> <mesh dir> <nranks>...): add_fesom_test at each rank
# count; <mesh dir> is the test's built-in default mesh
function(add_default_mesh_test target mesh)
    foreach(_np ${ARGN})
        add_fesom_test(${target} ${_np})
        require_files(${target}_np${_np} ${mesh}/nod2d.out)
    endforeach()
endfunction()

# Task 6: mesh + extended dist_N reading. test_mesh_read checks read_mesh against the mesh
# files (vertex counts, levels, double-halo nodes) and enforce_cw_orientation on flipped
# elements; default mesh = the NEMO-mixed cut-out. test_partit checks the appended dist_N
# sections (FESOM3_EXPECT_FORMAT=1: dist_N/format.out must be there).
add_default_mesh_test(test_mesh_read ${_mesh_mix} 1 2)
foreach(_case quad:${_mesh_q} mixed:${_mesh_mix} pi_w:${_mesh_piw})
    string(REPLACE ":" ";" _pair ${_case})
    list(GET _pair 0 _tag)
    list(GET _pair 1 _mesh)
    foreach(_np 1 2)
        add_mesh_test(test_mesh_read_${_tag}_np${_np} test_mesh_read ${_mesh} ${_np})
        add_mesh_test(test_partit_${_tag}_np${_np} test_partit ${_mesh} ${_np} FESOM3_EXPECT_FORMAT=1)
    endforeach()
endforeach()
# Task 7: mixed-mesh geometry identities (kite weights, Green-Gauss gradient, conventions,
# triangles == the previous P1 formulas); default mesh = the NEMO-mixed cut-out, the np 1
# run also checks the analytic mixed patch. pi_w at np 1/2 for the triangle comparison.
add_default_mesh_test(test_geom_mixed ${_mesh_mix} 1 2)
add_mesh_test(test_geom_mixed_pi_w_np1 test_geom_mixed ${_mesh_piw} 1 FESOM3_TEST_ROTATE=1)
add_mesh_test(test_geom_mixed_pi_w_np2 test_geom_mixed ${_mesh_piw} 2 FESOM3_TEST_ROTATE=1)
# Task 8: edge exchange (com_edge2D) and double node halo exchange (com_nod2D_full) round
# trips; default mesh = the NEMO quad cut-out (np 1 = the synthesis, empty lists).
add_default_mesh_test(test_halo_ext ${_mesh_q} 1)
foreach(_case quad:${_mesh_q} mixed:${_mesh_mix} pi_w:${_mesh_piw} core2_w:${_mesh_c2w})
    string(REPLACE ":" ";" _pair ${_case})
    list(GET _pair 0 _tag)
    list(GET _pair 1 _mesh)
    foreach(_np 2 4 8)
        add_mesh_test(test_halo_ext_${_tag}_np${_np} test_halo_ext ${_mesh} ${_np})
    endforeach()
endforeach()
# Task 9: SSH stiffness matrix on mixed meshes (sparsity with quad diagonals, row sums,
# symmetry, update, CG, checkerboard null mode, quadratic-field Laplacian on the
# parallelogram patch); default mesh = the NEMO quad cut-out, the np 1 run also checks the
# analytic mixed and parallelogram patches. Triangle meshes (rotated) compare against the
# previous assembly.
add_default_mesh_test(test_ssh_stiff ${_mesh_q} 1 2)
set_tests_properties(test_ssh_stiff_np1 PROPERTIES ENVIRONMENT "FESOM3_EXPECT_CHECKERBOARD=1")
foreach(_np 1 2 8)
    add_mesh_test(test_ssh_stiff_mixed_np${_np} test_ssh_stiff ${_mesh_mix} ${_np} FESOM3_TEST_ROTATE=0)
endforeach()
add_mesh_test(test_ssh_stiff_pi_w_np1    test_ssh_stiff ${_mesh_piw} 1 FESOM3_TEST_ROTATE=1)
add_mesh_test(test_ssh_stiff_pi_w_np2    test_ssh_stiff ${_mesh_piw} 2 FESOM3_TEST_ROTATE=1)
add_mesh_test(test_ssh_stiff_core2_w_np8 test_ssh_stiff ${_mesh_c2w} 8 FESOM3_TEST_ROTATE=1)
# Task 10: ocean dynamics kernels with the unified 4-slot formulas. test_pgf_mixed (PGF of a
# resting/linear state, node gather, SSH-gradient update; default mesh = the NEMO-mixed
# cut-out, the np 1 run also checks the analytic mixed patch); test_vinv / test_ivertvisc
# (exact uniform-flow / linear-shear / closed-form tests) on the analytic patches
# (FESOM3_TEST_PATCH=mixed|quad, np 1) and the unrotated NEMO cut-outs.
add_default_mesh_test(test_pgf_mixed ${_mesh_mix} 1 2)
add_mesh_test(test_pgf_mixed_quad_np2 test_pgf_mixed ${_mesh_q}   2 FESOM3_TEST_ROTATE=0)
add_mesh_test(test_pgf_mixed_pi_w_np1 test_pgf_mixed ${_mesh_piw} 1 FESOM3_TEST_ROTATE=1)
foreach(_t test_vinv test_ivertvisc)
    add_test(NAME ${_t}_patch_mixed_np1
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:${_t}>)
    set_tests_properties(${_t}_patch_mixed_np1 PROPERTIES ENVIRONMENT "FESOM3_TEST_PATCH=mixed")
endforeach()
add_test(NAME test_vinv_patch_quad_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:test_vinv>)
set_tests_properties(test_vinv_patch_quad_np1 PROPERTIES ENVIRONMENT "FESOM3_TEST_PATCH=quad")
add_mesh_test(test_vinv_mixed_np1       test_vinv      ${_mesh_mix} 1 FESOM3_TEST_ROTATE=0)
add_mesh_test(test_vinv_mixed_np2       test_vinv      ${_mesh_mix} 2 FESOM3_TEST_ROTATE=0)
add_mesh_test(test_vinv_quad_np2        test_vinv      ${_mesh_q}   2 FESOM3_TEST_ROTATE=0)
add_mesh_test(test_ivertvisc_mixed_np2  test_ivertvisc ${_mesh_mix} 2 FESOM3_TEST_ROTATE=0)
# Task 11: tracer kernels on mixed meshes. test_updown_mixed (MUSCL up/downwind element
# search, linear-tracer exactness of the reconstruction and of the advection; default mesh =
# the NEMO-mixed cut-out, the np 1 run also checks the analytic mixed and parallelogram
# patches); test_muscl_onthefly part R (kernel == the reference RK3 rule) on the unrotated
# NEMO cut-outs.
add_default_mesh_test(test_updown_mixed ${_mesh_mix} 1 2)
add_mesh_test(test_updown_mixed_quad_np2    test_updown_mixed   ${_mesh_q}   2 FESOM3_TEST_ROTATE=0)
add_mesh_test(test_updown_mixed_mixed_np8   test_updown_mixed   ${_mesh_mix} 8 FESOM3_TEST_ROTATE=0)
add_mesh_test(test_updown_mixed_pi_w_np2    test_updown_mixed   ${_mesh_piw} 2 FESOM3_TEST_ROTATE=1)
add_mesh_test(test_muscl_onthefly_mixed_np1 test_muscl_onthefly ${_mesh_mix} 1 FESOM3_TEST_ROTATE=0)
add_mesh_test(test_muscl_onthefly_mixed_np2 test_muscl_onthefly ${_mesh_mix} 2 FESOM3_TEST_ROTATE=0)
add_mesh_test(test_muscl_onthefly_quad_np2  test_muscl_onthefly ${_mesh_q}   2 FESOM3_TEST_ROTATE=0)
# Tasks 13, 15, 16, each with the default mesh = the NEMO-mixed cut-out at np 1/2 (the np 1
# run also checks the analytic patches), the NEMO quad cut-out at np 2 and the triangle
# meshes pi_w / core2_w rotated (triangles == the former formulas to roundoff):
#  - Task 13, test_coupling_mixed: forcing and ice-ocean coupling (kite-weighted
#    node->element surface stress of oce_fluxes_mom / forcing_stress_surf, momentum input
#    equal on nodes and elements, oce_fluxes balancing);
#  - Task 15, test_ice_dyn_mixed: sea-ice EVP / mEVP dynamics with the unified 4-slot
#    formulas (strain rates of uniform drift / rigid rotation / linear deformation,
#    stress2rhs = -adjoint of the strain rates, SSH-tilt node shares, rest state of uniform
#    ice, ocean2ice gather; analytic triangle, mixed and parallelogram patches);
#  - Task 16, test_ice_mass: sea-ice FCT mass matrix and lumped mass (user's specification:
#    matrix entries vs the element contributions, row sums == lumped mass == sum S/nv,
#    total == area, constant field stays constant, sum(lumped*f) conserved, FCT result
#    within the low-order bounds, triangles == the former FCT); a missing CSR column stops
#    the mass-matrix fill.
foreach(_t test_coupling_mixed test_ice_dyn_mixed test_ice_mass)
    add_default_mesh_test(${_t} ${_mesh_mix} 1 2)
    add_mesh_test(${_t}_quad_np2    ${_t} ${_mesh_q}   2 FESOM3_TEST_ROTATE=0)
    add_mesh_test(${_t}_pi_w_np1    ${_t} ${_mesh_piw} 1 FESOM3_TEST_ROTATE=1)
    add_mesh_test(${_t}_core2_w_np8 ${_t} ${_mesh_c2w} 8 FESOM3_TEST_ROTATE=1)
endforeach()
add_test(NAME test_ice_mass_missing_col_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:test_ice_mass>)
set_tests_properties(test_ice_mass_missing_col_np1 PROPERTIES
    ENVIRONMENT "FESOM3_TEST_MISSING_COL=1"
    PASS_REGULAR_EXPRESSION "not found in the ssh_stiff CSR row of node"
    FAIL_REGULAR_EXPRESSION "did NOT stop")
# Task 14: mesh diagnostics store in the hybrid layout (4-column face_nodes/face_edges/face_links/
# gradient_sca_x/y, a triangle's 4th = _FillValue -1 / -1 / -999), checked against elem2d.out and
# edge_nodes/edge_face_links, and partition-independent (np 1 == np N). Standard-library python.
foreach(_case mixed:${_mesh_mix}:1,2,8 quad:${_mesh_q}:1,2 pi_w:${_mesh_piw}:1,2)
    string(REPLACE ":" ";" _triple ${_case})
    list(GET _triple 0 _tag)
    list(GET _triple 1 _mesh)
    list(GET _triple 2 _nps)
    add_test(NAME meshdiag_${_tag}
             COMMAND python3 ${CMAKE_SOURCE_DIR}/tools/check_meshdiag.py ${_mesh}
                     ${CMAKE_CURRENT_BINARY_DIR}/meshdiag_${_tag} $<TARGET_FILE:fesom_meshdiagdump> ${_nps}
                     ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS})
    set_tests_properties(meshdiag_${_tag} PROPERTIES PASS_REGULAR_EXPRESSION "MESHDIAG CHECK PASS")
    require_files(meshdiag_${_tag} ${_mesh}/nod2d.out)
endforeach()
# without dist_N/format.out (the original pi dist_2) exchange_edge stops with a message
add_mesh_test(test_halo_ext_noformat_np2 test_halo_ext ${_mesh_pi} 2)
set_tests_properties(test_halo_ext_noformat_np2 PROPERTIES
    PASS_REGULAR_EXPRESSION "com_edge2D is not available: dist_N has no format.out"
    FAIL_REGULAR_EXPRESSION "did NOT stop")
# more than MAX_NEIGHBOR_PARTITIONS (32) receive/send PEs in a com list stop read_com
# (synthetic dist_2 in the build tree, no mesh needed)
foreach(_side recv:receive send:send)
    string(REPLACE ":" ";" _pair ${_side})
    list(GET _pair 0 _s)
    list(GET _pair 1 _word)
    add_test(NAME test_partit_many_neighbours_${_s}_np2
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 2 $<TARGET_FILE:test_partit>)
    set_tests_properties(test_partit_many_neighbours_${_s}_np2 PROPERTIES
        ENVIRONMENT "FESOM3_TEST_MANY_NEIGHBOURS=${_s}"
        PASS_REGULAR_EXPRESSION "com_nod2D_full: 33 ${_word} PEs, outside 0..MAX_NEIGHBOR_PARTITIONS=32"
        FAIL_REGULAR_EXPRESSION "did NOT stop")
endforeach()
# the old pi dist_2 (no format.out, FESOM2 layout only) is still read
add_mesh_test(test_mesh_read_pi_np2 test_mesh_read ${_mesh_pi} 2)
# a malformed elem2d.out stops read_mesh with a message (regex pair as test_minlayers);
# FESOM3_TEST_MALFORMED selects the defect
foreach(_case "1:element 3 has a wrong number of node ids: \"2 4\""
              "2:element 1 has neither 3 nor 4 node ids: \"1 2\""
              "3:element 2 is not a list of integers: \"1 3 x\""
              "4:element 2 has a node id outside 1..nod2D: \"1 3 9\""
              "5:element 3 is missing")
    string(FIND "${_case}" ":" _colon)
    string(SUBSTRING "${_case}" 0 ${_colon} _k)
    math(EXPR _colon "${_colon} + 1")
    string(SUBSTRING "${_case}" ${_colon} -1 _msg)
    set(_name test_mesh_read_malformed_np1)
    if(NOT _k EQUAL 1)
        set(_name test_mesh_read_malformed${_k}_np1)
    endif()
    add_test(NAME ${_name}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:test_mesh_read>)
    set_tests_properties(${_name} PROPERTIES
        ENVIRONMENT "FESOM3_TEST_MALFORMED=${_k}"
        PASS_REGULAR_EXPRESSION "the line of ${_msg}"
        FAIL_REGULAR_EXPRESSION "did NOT stop")
endforeach()
# fesom_ini (mesh_part/) on a copy of the NEMO quad cut-out in the build tree
# (tools/check_fesom_ini.sh): the FESOM2 level finder on quads (aux3d.out) and the error
# stops of use_depthfile="nemo", the >= 2-layer rule and the non-convex quad check.
foreach(_case "aux3d:check_fesom_ini: aux3d OK"
              "nodepth:use_depthfile=\"nemo\" needs the file"
              "onelem:set use_depthonelem=.false."
              "count1:depth.out layer counts outside 2..nl-1"
              "nonconvex:test_tri_ini: non-convex quad")
    string(FIND "${_case}" ":" _colon)
    string(SUBSTRING "${_case}" 0 ${_colon} _c)
    math(EXPR _colon "${_colon} + 1")
    string(SUBSTRING "${_case}" ${_colon} -1 _msg)
    add_test(NAME fesom_ini_${_c}
             COMMAND bash ${CMAKE_SOURCE_DIR}/tools/check_fesom_ini.sh $<TARGET_FILE:fesom_ini>
                     ${_mesh_q} ${CMAKE_CURRENT_BINARY_DIR}/fesom_ini_${_c} ${_c})
    set_tests_properties(fesom_ini_${_c} PROPERTIES
        PASS_REGULAR_EXPRESSION "${_msg}"
        FAIL_REGULAR_EXPRESSION "check_fesom_ini: dist_2 written;check_fesom_ini: aux3d FAILED")
    require_files(fesom_ini_${_c} ${_mesh_q}/nod2d.out ${_mesh_q}/depth.out)
endforeach()

# dump_diff.py self-check (pure python; no external data needed).
add_test(NAME dump_diff_selftest
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tools/dump_diff.py --selftest)

# End-to-end analytic driver (built from src/drivers/). Runs clean on 1 + 2 ranks.
add_test(NAME fesom_analytic_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:fesom_analytic>)
add_test(NAME fesom_analytic_np2
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 2 $<TARGET_FILE:fesom_analytic>)

# Conservation + no-leakage gate for the bottom-at-vertices change (docs/plans/completed/
# 20260910-fesom3-bottom-at-vertices.md). zstar conserves to round-off, so this is a hard
# gate; the driver itself error-stops when FESOM3_CONSERVE_TOL is exceeded. Needs the pi
# mesh and the pool IC file, so it is skipped (not failed) when those are absent.
foreach(_np 1 2)
    add_test(NAME fesom_conserve_zstar_np${_np}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np}
                     $<TARGET_FILE:fesom_conserve>)
    set_tests_properties(fesom_conserve_zstar_np${_np} PROPERTIES
        ENVIRONMENT "FESOM3_MESH_DIR=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi;FESOM3_IC_FILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc;FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;FESOM3_CONSERVE_TOL=1e-12"
        SKIP_REGULAR_EXPRESSION "cannot open;No such file")
endforeach()

# Redi variant: diff_ver_part_redi_expl is the only caller of the tr_xynodes node-average
# whose denominator this change touched, and it runs only when Redi is on.
add_test(NAME fesom_conserve_zstar_redi_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1
                 $<TARGET_FILE:fesom_conserve>)
set_tests_properties(fesom_conserve_zstar_redi_np1 PROPERTIES
    ENVIRONMENT "FESOM3_MESH_DIR=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi;FESOM3_IC_FILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc;FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;FESOM3_REDI=1;FESOM3_CONSERVE_TOL=1e-12"
    SKIP_REGULAR_EXPRESSION "cannot open;No such file")

# One case per vertical-mixing scheme: the area/areasvol consumers live in the vertical
# diffusion TDMA and the KPP non-local / shortwave terms, reachable only with the scheme on.
foreach(_mix KPP TKE)
    add_test(NAME fesom_conserve_zstar_${_mix}_np1
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1
                     $<TARGET_FILE:fesom_conserve>)
    set_tests_properties(fesom_conserve_zstar_${_mix}_np1 PROPERTIES
        ENVIRONMENT "FESOM3_MESH_DIR=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi;FESOM3_IC_FILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc;FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;FESOM3_MIX_${_mix}=1;FESOM3_CONSERVE_TOL=1e-12"
        SKIP_REGULAR_EXPRESSION "cannot open;No such file")
endforeach()

# Task 12 (hybrid meshes): the mixing / GM kernels (PP default, KPP, TKE, GM, Redi) on the
# NEMO-mixed cut-out (closed regional boundaries, triangles + quads), np 1 and 2. The
# element means of Av/viscAE/fer_uv are kite-weighted (elem_w); the driver's helem
# invariant checks the same elem_w mean on every element.
foreach(_np 1 2)
    foreach(_mixcase PP KPP:FESOM3_MIX_KPP TKE:FESOM3_MIX_TKE GM:FESOM3_FER_GM REDI:FESOM3_REDI)
        string(REPLACE ":" ";" _pair ${_mixcase})
        list(GET _pair 0 _tag)
        set(_knob "")
        list(LENGTH _pair _len)
        if(_len GREATER 1)
            list(GET _pair 1 _var)
            set(_knob "${_var}=1;")
        endif()
        add_test(NAME fesom_conserve_mixcut_${_tag}_np${_np}
                 COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np}
                         $<TARGET_FILE:fesom_conserve>)
        set_tests_properties(fesom_conserve_mixcut_${_tag}_np${_np} PROPERTIES
            ENVIRONMENT "FESOM3_MESH_DIR=${_mesh_mix};FESOM3_IC_FILE=${_icf};FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;${_knob}FESOM3_CONSERVE_TOL=1e-12")
        require_files(fesom_conserve_mixcut_${_tag}_np${_np} ${_mesh_mix}/nod2d.out ${_icf})
    endforeach()
endforeach()

# Task 17 (hybrid meshes): conservation with SEA ICE ON on the NEMO-quad and NEMO-mixed
# cut-outs, np 1 and 2 (EVP) + mEVP at np 2. The ocean physics is the production set the
# driver supports (TKE, GM + Redi, splines, vector-invariant momentum with rvo_upwind 1, w
# split) on zstar (linfs does not conserve); the ice is seeded and driven by a wind-on-ice
# stress, dynamics + FCT advection + ice-ocean drag every step, no thermodynamics
# (FESOM3_CONSERVE_ICE, fesom_conserve header). Heat/salt/volume and the three ice masses
# sum(fct_lumpedmass*f) within 1e-12; measured <= 1.2e-14 (max|uice| 0.35-0.39 m/s).
# _prod_ocean = prod_ocean_env of tools/prod_env.sh (keep the two lists equal).
set(_prod_ocean "FESOM3_MIX_TKE=1;FESOM3_REDI=1;FESOM3_FER_GM=1;FESOM3_SHEAR_SPLINES=1;FESOM3_N2_SPLINES=1;FESOM3_MOMADV_OPT=1;FESOM3_RVO_UPWIND=1.0;FESOM3_WSPLIT=1")
foreach(_cut q:${_mesh_q} mix:${_mesh_mix})
    string(REPLACE ":" ";" _pair ${_cut})
    list(GET _pair 0 _tag)
    list(GET _pair 1 _mesh)
    foreach(_case 1:0 2:0 2:1)
        string(REPLACE ":" ";" _c ${_case})
        list(GET _c 0 _np)
        list(GET _c 1 _evp)
        set(_name fesom_conserve_${_tag}cut_ice_np${_np})
        if(_evp EQUAL 1)
            set(_name fesom_conserve_${_tag}cut_ice_mevp_np${_np})
        endif()
        add_test(NAME ${_name}
                 COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np}
                         $<TARGET_FILE:fesom_conserve>)
        set_tests_properties(${_name} PROPERTIES
            ENVIRONMENT "FESOM3_MESH_DIR=${_mesh};FESOM3_IC_FILE=${_icf};FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;${_prod_ocean};FESOM3_CONSERVE_ICE=1;FESOM3_WHICHEVP=${_evp};FESOM3_CONSERVE_TOL=1e-12"
            PASS_REGULAR_EXPRESSION "sea ice ON.*CONSERVATION OK")
        require_files(${_name} ${_mesh}/nod2d.out ${_icf})
    endforeach()
endforeach()

# Smooth w = w_e + w_i split (use_wsplit, docs/plans/completed/2026-10-02-wsplit-smooth.md) with all
# its consumers live: momentum TDMA + advection, the tracer-diffusion TDMA and the FCT
# adv_tra_vert_impl. The 20-step cold start never reaches CFL_z ~ 1 (measured max
# 2.76e-2), so the cap is 0.005 (onset 0.0025 by the driver's 0.5*maxcfl default): 8632
# owned faces split, 500 with an implicit share >= 0.5. FESOM3_WSPLIT_EXPECT_SPLIT=1000
# error-stops unless >= 1000 faces split and at least one has f >= 0.5 (no vacuous pass;
# 8x below the measured count). The tuned triple is _wsplit_cfg below -- the same triple
# as WSPLIT_CFG in tools/run_conserve_pi.sh; re-tune the cap in both places.
# np 2 runs the next step's consumers of the halo w_e/w_i
# (the momentum TDMA and, in the vinv variant, the vector-invariant advection read them at
# the halo vertices of owned elements), but the conserve gate is BLIND to a wrong halo
# value (heat/salt/volume still telescope): the halo contract itself is pinned bitwise by
# test_wsplit S5 at np 2. The vinv variant runs the advective-form implicit w_i correction
# of the momentum TDMA.
set(_wsplit_cfg "FESOM3_WSPLIT=1;FESOM3_WSPLIT_MAXCFL=0.005;FESOM3_WSPLIT_EXPECT_SPLIT=1000")
foreach(_np 1 2)
    add_test(NAME fesom_conserve_zstar_wsplit_np${_np}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np}
                     $<TARGET_FILE:fesom_conserve>)
    set_tests_properties(fesom_conserve_zstar_wsplit_np${_np} PROPERTIES
        ENVIRONMENT "FESOM3_MESH_DIR=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi;FESOM3_IC_FILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc;FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;${_wsplit_cfg};FESOM3_CONSERVE_TOL=1e-12"
        SKIP_REGULAR_EXPRESSION "cannot open;No such file")
    add_test(NAME fesom_conserve_zstar_vinv_wsplit_np${_np}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np}
                     $<TARGET_FILE:fesom_conserve>)
    set_tests_properties(fesom_conserve_zstar_vinv_wsplit_np${_np} PROPERTIES
        ENVIRONMENT "FESOM3_MESH_DIR=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi;FESOM3_IC_FILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc;FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;FESOM3_MOMADV_OPT=1;${_wsplit_cfg};FESOM3_CONSERVE_TOL=1e-12"
        SKIP_REGULAR_EXPRESSION "cannot open;No such file")
endforeach()

# The teeth of the guards above, as ctests (a guard that was only ever tripped by hand is
# a guard that can die unnoticed):
#  * vacuous config: cap 0.5 on the cold start splits nothing, so EXPECT_SPLIT must print
#    "vacuous config" and stop BEFORE the drift report (FAIL_REGULAR_EXPRESSION). The pass
#    regex also pins the driver-side default mincfl = 0.5*maxcfl (0.2500 at the cap 0.5 --
#    the type's absolute 0.5 would print 0.5000; " *" so the f8.4 blanks are not counted);
#    CMake's "." matches newlines.
#  * bad params: (mincfl, maxcfl) = (0.9, 0.5) must be rejected before init.
#  * bad value: an unparsable FESOM3_WSPLIT_MAXCFL must be an error, not the default cap.
# Each is a PASS regex (the message: stopped for the RIGHT reason) plus a FAIL regex (a
# line the driver prints only if it went on). A PASS regex makes ctest ignore the exit
# code, and WILL_FAIL cannot restore it: CMake inverts the verdict AFTER the regexes, so
# WILL_FAIL + a FAIL regex would PASS a driver that prints the drift line and exits 0
# (checked with CMake 3.26) -- hence the regex pair, not WILL_FAIL.
# These three stop before the mesh is read (or never reach the drift), so they need no
# skip rule; the vacuous one keeps it because it does read the mesh.
add_test(NAME fesom_conserve_zstar_wsplit_vacuous_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1
                 $<TARGET_FILE:fesom_conserve>)
set_tests_properties(fesom_conserve_zstar_wsplit_vacuous_np1 PROPERTIES
    ENVIRONMENT "FESOM3_MESH_DIR=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi;FESOM3_IC_FILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc;FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;FESOM3_WSPLIT=1;FESOM3_WSPLIT_MAXCFL=0.5;FESOM3_WSPLIT_EXPECT_SPLIT=1;FESOM3_CONSERVE_TOL=1e-12"
    PASS_REGULAR_EXPRESSION "wsplit_mincfl= *0\\.2500, wsplit_maxcfl= *0\\.5000.*vacuous config"
    FAIL_REGULAR_EXPRESSION "relative drift over the run"
    SKIP_REGULAR_EXPRESSION "cannot open;No such file")
add_test(NAME fesom_conserve_wsplit_badparams_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1
                 $<TARGET_FILE:fesom_conserve>)
set_tests_properties(fesom_conserve_wsplit_badparams_np1 PROPERTIES
    ENVIRONMENT "FESOM3_WSPLIT=1;FESOM3_WSPLIT_MINCFL=0.9;FESOM3_WSPLIT_MAXCFL=0.5"
    PASS_REGULAR_EXPRESSION "need maxcfl > 0 and 0 <= mincfl <= maxcfl"
    FAIL_REGULAR_EXPRESSION "fesom_conserve: nod2D=")
add_test(NAME fesom_conserve_wsplit_badvalue_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1
                 $<TARGET_FILE:fesom_conserve>)
set_tests_properties(fesom_conserve_wsplit_badvalue_np1 PROPERTIES
    ENVIRONMENT "FESOM3_WSPLIT=1;FESOM3_WSPLIT_MAXCFL=0.5x"
    PASS_REGULAR_EXPRESSION "FESOM3_WSPLIT_MAXCFL: not a real"
    FAIL_REGULAR_EXPRESSION "fesom_conserve: nod2D=")

# The setup-time >= 2-layer column invariant (read_mesh's assert_min_layers, the
# precondition of every vertical kernel's row layout), tripped: test_minlayers shrinks one
# pi column to a single layer and re-walks derive_vertical_bounds, which must error stop
# with this message (same regex pair as the guards above).
add_fesom_test(test_minlayers 1)
set_tests_properties(test_minlayers_np1 PROPERTIES
    PASS_REGULAR_EXPRESSION "fewer than 2 layers"
    FAIL_REGULAR_EXPRESSION "did NOT stop")

# Shared-element sync (docs/plans/completed/20261012-shared-element-sync.md): test_shared_elem checks
# the list com_elem2D_shared against an independent oracle (S1), its exchange and the
# copy-spread probe (S2), and the mechanism through update_vel + compute_hbar_ale with the
# sync on and off (S3). pi has the FESOM2 dist layout of production dist_512; pi_w, core2_w
# and the NEMO cut-outs the fesom_ini layout.
# The triangle meshes are read rotated (FESOM3_TEST_ROTATE=1), the cut-outs unrotated.
# np 1 (the 1-rank synthesis): the list is empty.
add_fesom_test(test_shared_elem 1)
set_tests_properties(test_shared_elem_np1 PROPERTIES ENVIRONMENT "FESOM3_TEST_ROTATE=1")
foreach(_case pi:${_mesh_pi}:2,8:1 pi_w:${_mesh_piw}:2,4:1 mixed:${_mesh_mix}:2:0
              quad:${_mesh_q}:2:0 core2_w:${_mesh_c2w}:8:1)
    string(REPLACE ":" ";" _c ${_case})
    list(GET _c 0 _tag)
    list(GET _c 1 _mesh)
    list(GET _c 2 _nps)
    list(GET _c 3 _rot)
    string(REPLACE "," ";" _nps ${_nps})
    foreach(_np ${_nps})
        add_mesh_test(test_shared_elem_${_tag}_np${_np} test_shared_elem ${_mesh} ${_np}
                      FESOM3_TEST_ROTATE=${_rot})
    endforeach()
endforeach()
# a vertex without a known owner (com_nod2D receive blocks cleared) stops the builder
add_mesh_test(test_shared_elem_corrupt_np2 test_shared_elem ${_mesh_pi} 2
              FESOM3_TEST_SHARED_CORRUPT=1)
set_tests_properties(test_shared_elem_corrupt_np2 PROPERTIES
    PASS_REGULAR_EXPRESSION "has a vertex with no known owner"
    FAIL_REGULAR_EXPRESSION "did NOT stop")
# FESOM3_SHARED_SYNC: an invalid value stops the builder; 0 still builds the list and logs OFF
add_mesh_test(test_shared_elem_badknob_np2 test_shared_elem ${_mesh_pi} 2
              FESOM3_TEST_ROTATE=1 FESOM3_SHARED_SYNC=2)
set_tests_properties(test_shared_elem_badknob_np2 PROPERTIES
    PASS_REGULAR_EXPRESSION "FESOM3_SHARED_SYNC must be 0 or 1"
    FAIL_REGULAR_EXPRESSION "did NOT stop")
add_mesh_test(test_shared_elem_off_np2 test_shared_elem ${_mesh_pi} 2
              FESOM3_TEST_ROTATE=1 FESOM3_SHARED_SYNC=0)
set_tests_properties(test_shared_elem_off_np2 PROPERTIES
    PASS_REGULAR_EXPRESSION "shared-element sync: OFF, [0-9]+ shared elements.*ALL PASS"
    FAIL_REGULAR_EXPRESSION "FAIL")
# Shared-element gate in the full model (fesom_conserve, header SHARED ELEMENTS): linfs,
# vector-invariant momentum with rvo_upwind 1 (the CORE2 leak's trigger), sea ice, np 2.
# FESOM3_DIVSUM_TOL=1e-12: the divergence-sum ratio must stay below 1e-12 in every step and
# the uv / uv_rhsAB / sigma copy spreads must be 0 (tools/run_conserve_pi.sh runs np 8 and the
# scalar advection too).
add_test(NAME fesom_conserve_linfs_vinv_shsync_np2
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 2 $<TARGET_FILE:fesom_conserve>)
set_tests_properties(fesom_conserve_linfs_vinv_shsync_np2 PROPERTIES
    ENVIRONMENT "FESOM3_MESH_DIR=${_mesh_pi};FESOM3_IC_FILE=${_icf};FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=linfs;FESOM3_MOMADV_OPT=1;FESOM3_RVO_UPWIND=1.0;FESOM3_CONSERVE_ICE=1;FESOM3_DIVSUM_TOL=1e-12"
    PASS_REGULAR_EXPRESSION "shared elements \\(sync ON\\).* -- OK within tol"
    FAIL_REGULAR_EXPRESSION "VIOLATED"
    SKIP_REGULAR_EXPRESSION "cannot open;No such file")
