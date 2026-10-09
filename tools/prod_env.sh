# The production physics of work/job_levante as FESOM3_* environment variables, sourced by
# tools/run_tri_regression.sh and tools/run_cutout_short.sh (one list for both). The ctests
# keep the ocean part as _prod_ocean in test/CMakeLists.txt (CMake cannot source this file).
PROD_POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
PROD_IC=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc

# prod_ocean_env: the ocean physics fesom_conserve also supports (TKE, GM + Redi, splines,
# vector-invariant momentum with rvo_upwind 1, w split)
prod_ocean_env() {
    export FESOM3_MIX_TKE=1 FESOM3_REDI=1 FESOM3_FER_GM=1 FESOM3_SHEAR_SPLINES=1 \
        FESOM3_N2_SPLINES=1 FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=1.0 FESOM3_WSPLIT=1
}

# prod_env <ssh_soltol>: the production driver set (fesom_lifecycle_native_mr): PHC winter IC,
# JRA55 from 1958-01-01, linfs, EVP sea ice, SW penetration, the ocean physics above, step
# diagnostics, the SSH solved to <ssh_soltol>
prod_env() {
    prod_ocean_env
    export FESOM3_IC_FILE=$PROD_IC FESOM3_FORCING_DIR=$PROD_POOL FESOM3_FORCING=JRA55 \
        FESOM3_RUNOFF_FILE=$PROD_POOL/CORE2_runoff.nc FESOM3_SSS_FILE=$PROD_POOL/PHC2_salx.nc \
        FESOM3_START_CLOCK="0 1 1958" FESOM3_WHICH_ALE=linfs FESOM3_WHICHEVP=0 \
        FESOM3_SW_PENE=1 FESOM3_CHL_SWEENEY=1 FESOM3_STEP_DIAG=1 \
        FESOM3_SSH_SOLTOL=$1 FESOM3_SSH_MAXITER=20000
}
