#!/bin/bash
#SBATCH --job-name=recom_ci
#SBATCH -p compute
#SBATCH --ntasks=2
#SBATCH --time=00:30:00
#SBATCH -o recom_ci_%j.out
#SBATCH -e recom_ci_%j.err
#SBATCH -A ba0989

set -e

REPO_DIR="/work/ab0246/a270092/model_codes/fesom-2.7"
SIF="${REPO_DIR}/fesom2-ci.sif"

module load singularity

# The SIF must be pre-pulled on the login node (SLURM compute nodes
# have a read-only TMPDIR that breaks singularity pull). To pull:
#   module load singularity
#   singularity pull "$SIF" docker://ghcr.io/fesom/fesom2_docker:fesom2_test_refactoring-master
if [ ! -f "$SIF" ]; then
    echo "ERROR: $SIF not found. Pull it on the login node first (see script comments)."
    exit 1
fi

echo "==> Running RECOM CI test inside container..."
singularity exec --bind "${REPO_DIR}:/fesom/fesom2" "$SIF" bash -l -c '
set -e
cd /fesom/fesom2

echo "==> Cleaning stale build cache..."
rm -rf build

echo "==> Compiling with REcoM..."
./configure.sh ubuntu -DRECOM_COUPLED=ON

echo "==> Creating test run directory..."
mkrun recom test_pi_recom -m docker

echo "==> Setting up REcoM namelists..."
cp config/bin_2p3z2d/namelist.recom work_recom/
cp config/bin_2p3z2d/namelist.tra   work_recom/

cd work_recom
sed -i "/^&nam_rsbc/,/^\//d" namelist.recom
sed -i "s|REcoM_restart.*=.*\.true\.|REcoM_restart         =.false.|" namelist.recom
sed -i "s|fe_pisces_opa_eq_init_3D_changed_name\.nc|fe5deg.nc|" namelist.tra
sed -i "s|woa18_all_o00_01_mmol_fesom2\.nc|oxy5deg.nc|" namelist.tra
sed -i "s|woa13_all_i00_01_fesom2\.nc|si5deg.nc|" namelist.tra
sed -i "s|GLODAPv2\.2016b\.TAlk_fesom2_mmol_fix_z_Fillvalue\.nc|talk5deg.nc|" namelist.tra
sed -i "s|GLODAPv2\.2016b\.TCO2_fesom2_mmol_fix_z_Fillvalue\.nc|tco2_5deg.nc|" namelist.tra
sed -i "s|woa13_all_n00_01_fesom2\.nc|din5deg.nc|" namelist.tra
sed -i "s|phc3\.0_winter\.nc|woa18_netcdf_5deg.nc|g" namelist.tra
sed -i "s|l_mslp *= *\.false\.|l_mslp=.true.|" namelist.forcing
cd ..

echo "==> Copying input data..."
cp test/input/recom/*.nc work_recom/
(cd test/input/recom && ln -sf DustClimMonthlyAlbani_pimesh.nc DustClimMonthlyAlbani.nc)
cp test/input/global/*.1948.nc work_recom/
cp test/input/global/runoff.nc work_recom/
cp test/input/global/PHC2_salx.nc work_recom/

echo "==> Running FESOM2 + REcoM (32 timesteps, 1 day)..."
cd work_recom
chmod +x job_docker_new
./job_docker_new

echo "==> Validating results..."
fcheck .

echo "==> RECOM CI test completed successfully!"
'
