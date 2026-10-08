# FESOM3 self-tests. add_fesom_test(<name> <nranks>) is defined in the top-level
# CMakeLists; each builds test/<name>.F90 against libfesom3 and registers a ctest.

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

# Hybrid meshes (docs/plans/2026-10-08-hybrid-meshes.md, Task 6): mesh + extended dist_N
# reading. test_mesh_read checks read_mesh against the mesh files (vertex counts, levels,
# double-halo nodes); default mesh = the NEMO-mixed cut-out. test_partit checks the
# appended dist_N sections (dist_N/format.out). The writable meshes live under
# /work/ab0995/a270029/meshes; a missing mesh skips the test.
set(_hyb_meshes /work/ab0995/a270029/meshes)
add_fesom_test(test_mesh_read 1)
add_fesom_test(test_mesh_read 2)
foreach(_case quad:nemo025q_cut mixed:nemo025mix_cut pi_w:pi_w)
    string(REPLACE ":" ";" _pair ${_case})
    list(GET _pair 0 _tag)
    list(GET _pair 1 _mesh)
    foreach(_np 1 2)
        foreach(_t test_mesh_read test_partit)
            add_test(NAME ${_t}_${_tag}_np${_np}
                     COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np}
                             $<TARGET_FILE:${_t}>)
            set_tests_properties(${_t}_${_tag}_np${_np} PROPERTIES
                ENVIRONMENT "FESOM3_MESH_DIR=${_hyb_meshes}/${_mesh}"
                SKIP_REGULAR_EXPRESSION "cannot open;No such file")
        endforeach()
    endforeach()
endforeach()
# Task 7: mixed-mesh geometry identities (kite weights, Green-Gauss gradient, conventions,
# triangles == the previous P1 formulas); default mesh = the NEMO-mixed cut-out, the np 1
# run also checks the analytic mixed patch. pi_w at np 1/2 for the triangle comparison.
add_fesom_test(test_geom_mixed 1)
add_fesom_test(test_geom_mixed 2)
foreach(_np 1 2)
    add_test(NAME test_geom_mixed_pi_w_np${_np}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np} $<TARGET_FILE:test_geom_mixed>)
    set_tests_properties(test_geom_mixed_pi_w_np${_np} PROPERTIES
        ENVIRONMENT "FESOM3_MESH_DIR=${_hyb_meshes}/pi_w;FESOM3_TEST_ROTATE=1"
        SKIP_REGULAR_EXPRESSION "cannot open;No such file")
endforeach()
# Task 8: edge exchange (com_edge2D) and double node halo exchange (com_nod2D_full) round
# trips; default mesh = the NEMO quad cut-out (np 1 = the synthesis, empty lists).
add_fesom_test(test_halo_ext 1)
foreach(_case quad:nemo025q_cut:2,4,8 mixed:nemo025mix_cut:4 pi_w:pi_w:2,4,8 core2_w:core2_w:2,4,8)
    string(REPLACE ":" ";" _trip ${_case})
    list(GET _trip 0 _tag)
    list(GET _trip 1 _mesh)
    list(GET _trip 2 _nps)
    string(REPLACE "," ";" _nps ${_nps})
    foreach(_np ${_nps})
        add_test(NAME test_halo_ext_${_tag}_np${_np}
                 COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np} $<TARGET_FILE:test_halo_ext>)
        set_tests_properties(test_halo_ext_${_tag}_np${_np} PROPERTIES
            ENVIRONMENT "FESOM3_MESH_DIR=${_hyb_meshes}/${_mesh}"
            SKIP_REGULAR_EXPRESSION "cannot open;No such file")
    endforeach()
endforeach()
# Task 9: SSH stiffness matrix on mixed meshes (sparsity with quad diagonals, row sums,
# symmetry, update, CG, checkerboard null mode); default mesh = the NEMO quad cut-out, the
# np 1 run also checks the analytic mixed and parallelogram patches. Triangle meshes
# (rotated) compare against the previous assembly.
add_fesom_test(test_ssh_stiff 1)
add_fesom_test(test_ssh_stiff 2)
foreach(_case mixed:nemo025mix_cut:1,2,8:0 pi_w:pi_w:1,2:1 core2_w:core2_w:8:1)
    string(REPLACE ":" ";" _quad ${_case})
    list(GET _quad 0 _tag)
    list(GET _quad 1 _mesh)
    list(GET _quad 2 _nps)
    list(GET _quad 3 _rot)
    string(REPLACE "," ";" _nps ${_nps})
    foreach(_np ${_nps})
        add_test(NAME test_ssh_stiff_${_tag}_np${_np}
                 COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np} $<TARGET_FILE:test_ssh_stiff>)
        set_tests_properties(test_ssh_stiff_${_tag}_np${_np} PROPERTIES
            ENVIRONMENT "FESOM3_MESH_DIR=${_hyb_meshes}/${_mesh};FESOM3_TEST_ROTATE=${_rot}"
            SKIP_REGULAR_EXPRESSION "cannot open;No such file")
    endforeach()
endforeach()
# Task 10: ocean dynamics kernels with the unified 4-slot formulas. test_pgf_mixed (PGF of a
# resting/linear state, node gather, SSH-gradient update; default mesh = the NEMO-mixed
# cut-out, the np 1 run also checks the analytic mixed patch); test_vinv / test_ivertvisc
# (exact uniform-flow / linear-shear / closed-form tests) on the analytic patches
# (FESOM3_TEST_PATCH=mixed|quad, np 1) and the unrotated NEMO cut-outs.
add_fesom_test(test_pgf_mixed 1)
add_fesom_test(test_pgf_mixed 2)
foreach(_case quad:nemo025q_cut:2:0 pi_w:pi_w:1:1)
    string(REPLACE ":" ";" _quad ${_case})
    list(GET _quad 0 _tag)
    list(GET _quad 1 _mesh)
    list(GET _quad 2 _np)
    list(GET _quad 3 _rot)
    add_test(NAME test_pgf_mixed_${_tag}_np${_np}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np} $<TARGET_FILE:test_pgf_mixed>)
    set_tests_properties(test_pgf_mixed_${_tag}_np${_np} PROPERTIES
        ENVIRONMENT "FESOM3_MESH_DIR=${_hyb_meshes}/${_mesh};FESOM3_TEST_ROTATE=${_rot}"
        SKIP_REGULAR_EXPRESSION "cannot open;No such file")
endforeach()
foreach(_t test_vinv test_ivertvisc)
    add_test(NAME ${_t}_patch_mixed_np1
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:${_t}>)
    set_tests_properties(${_t}_patch_mixed_np1 PROPERTIES ENVIRONMENT "FESOM3_TEST_PATCH=mixed")
endforeach()
add_test(NAME test_vinv_patch_quad_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:test_vinv>)
set_tests_properties(test_vinv_patch_quad_np1 PROPERTIES ENVIRONMENT "FESOM3_TEST_PATCH=quad")
foreach(_case test_vinv:mixed:nemo025mix_cut:1 test_vinv:mixed:nemo025mix_cut:2
              test_vinv:quad:nemo025q_cut:2 test_ivertvisc:mixed:nemo025mix_cut:2)
    string(REPLACE ":" ";" _quad ${_case})
    list(GET _quad 0 _t)
    list(GET _quad 1 _tag)
    list(GET _quad 2 _mesh)
    list(GET _quad 3 _np)
    add_test(NAME ${_t}_${_tag}_np${_np}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np} $<TARGET_FILE:${_t}>)
    set_tests_properties(${_t}_${_tag}_np${_np} PROPERTIES
        ENVIRONMENT "FESOM3_MESH_DIR=${_hyb_meshes}/${_mesh};FESOM3_TEST_ROTATE=0"
        SKIP_REGULAR_EXPRESSION "cannot open;No such file")
endforeach()
# Task 11: tracer kernels on mixed meshes. test_updown_mixed (MUSCL up/downwind element
# search, linear-tracer exactness of the reconstruction and of the advection; default mesh =
# the NEMO-mixed cut-out, the np 1 run also checks the analytic mixed and parallelogram
# patches); test_muscl_onthefly part R (kernel == the reference RK3 rule) on the unrotated
# NEMO cut-outs.
add_fesom_test(test_updown_mixed 1)
add_fesom_test(test_updown_mixed 2)
foreach(_case test_updown_mixed:quad:nemo025q_cut:2:0 test_updown_mixed:mixed:nemo025mix_cut:8:0
              test_updown_mixed:pi_w:pi_w:2:1
              test_muscl_onthefly:mixed:nemo025mix_cut:1:0 test_muscl_onthefly:mixed:nemo025mix_cut:2:0
              test_muscl_onthefly:quad:nemo025q_cut:2:0)
    string(REPLACE ":" ";" _quint ${_case})
    list(GET _quint 0 _t)
    list(GET _quint 1 _tag)
    list(GET _quint 2 _mesh)
    list(GET _quint 3 _np)
    list(GET _quint 4 _rot)
    add_test(NAME ${_t}_${_tag}_np${_np}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np} $<TARGET_FILE:${_t}>)
    set_tests_properties(${_t}_${_tag}_np${_np} PROPERTIES
        ENVIRONMENT "FESOM3_MESH_DIR=${_hyb_meshes}/${_mesh};FESOM3_TEST_ROTATE=${_rot}"
        SKIP_REGULAR_EXPRESSION "cannot open;No such file")
endforeach()
# without dist_N/format.out (the original pi dist_2) exchange_edge stops with a message
add_test(NAME test_halo_ext_noformat_np2
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 2 $<TARGET_FILE:test_halo_ext>)
set_tests_properties(test_halo_ext_noformat_np2 PROPERTIES
    ENVIRONMENT "FESOM3_MESH_DIR=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi"
    PASS_REGULAR_EXPRESSION "com_edge2D is not available: dist_N has no format.out"
    FAIL_REGULAR_EXPRESSION "did NOT stop")
# the old pi dist_2 (no format.out, FESOM2 layout only) is still read
add_test(NAME test_mesh_read_pi_np2
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 2 $<TARGET_FILE:test_mesh_read>)
set_tests_properties(test_mesh_read_pi_np2 PROPERTIES
    ENVIRONMENT "FESOM3_MESH_DIR=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi"
    SKIP_REGULAR_EXPRESSION "cannot open;No such file")
# a malformed element line stops read_mesh with a message (regex pair as test_minlayers)
add_test(NAME test_mesh_read_malformed_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:test_mesh_read>)
set_tests_properties(test_mesh_read_malformed_np1 PROPERTIES
    ENVIRONMENT "FESOM3_TEST_MALFORMED=1"
    PASS_REGULAR_EXPRESSION "the line of element 3 has a wrong number of node ids: \"2 4\""
    FAIL_REGULAR_EXPRESSION "did NOT stop")

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
    foreach(_mix PP KPP:FESOM3_MIX_KPP TKE:FESOM3_MIX_TKE GM:FESOM3_FER_GM REDI:FESOM3_REDI)
        string(REPLACE ":" ";" _pair ${_mix})
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
            ENVIRONMENT "FESOM3_MESH_DIR=${_hyb_meshes}/nemo025mix_cut;FESOM3_IC_FILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc;FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;${_knob}FESOM3_CONSERVE_TOL=1e-12"
            SKIP_REGULAR_EXPRESSION "cannot open;No such file")
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
