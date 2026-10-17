#!/bin/bash
#SBATCH --job-name=rfc
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 32
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --output=logs/rfc/03_rfc_%j.out
#SBATCH --error=logs/rfc/03_rfc_%j.err
#SBATCH --mail-type=END,FAIL

# Random Forest analysis (notebooks/03_16S_ITS1_RFC.py), batch version.
#
# Runtime reference: on a 6-core laptop the importance step alone took ~15 h with
# everything else served from cache. It is the dominant cost -- for the Combined
# feature set that is 534 features x 10 repeats x 5 folds, about 27k model
# evaluations per target.
#
# Cores: parallelism lives in the loops, not in the trees -- N_JOBS_CV=-1 while each
# forest runs n_jobs=1, so the two never nest and oversubscribe. That matters here
# because the expensive step, permutation_importance, is itself called with
# n_jobs=N_JOBS_CV and so spreads across the whole allocation, which it could not do
# on the laptop run. 32 is a reasonable first guess; check
# `sinfo -p shared -o "%n %c %m"` for what the nodes actually offer before raising it.
#
# RFC_QUICK=1 gives a cheap pass that exercises every code path in minutes --
# run that first to confirm the environment before committing to the full job.

module load python/3.12.5

source ~/nobackup/venvs/rfc/bin/activate

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/rfc results/figures/random_forest results/tables/rf_cache

# Keep the numerical libraries inside the SLURM allocation. Without this they
# size their thread pools from the physical core count of the whole node and
# oversubscribe, which on this workload is slower than staying within the quota.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export LOKY_MAX_CPU_COUNT=${SLURM_CPUS_PER_TASK:-1}
export MPLBACKEND=Agg

# Without this the whole job runs on ONE core, silently.
#
# Hamilton sets OMP_PROC_BIND=spread and OMP_PLACES=cores in every job. When
# sklearn is imported it initialises libgomp, OpenMP then binds the master thread
# to a single place, and the process affinity mask collapses from 32 CPUs to 1.
# joblib reads that mask, so cpu_count() returns 1 and every n_jobs=-1 in the
# analysis quietly means "sequential". Measured by bisecting the imports:
#
#   numpy/pandas/scipy/seaborn/biom  affinity=32
#   sklearn                          affinity=1     <- here
#
# Measured effect on the dominant phase (64 permutations, farm_id/ITS1, warm
# worker pool): 185.1 s sequential vs 31.0 s parallel, a 6x difference.
export OMP_PROC_BIND=false

export RFC_QUICK=${RFC_QUICK:-0}

# Tag every output produced here, so cluster and laptop results can share a
# directory after an rsync instead of silently overwriting each other:
#   results/figures/random_forest/03_hpc_*.png
#   results/tables/03_hpc_*.csv
#   results/tables/rf_cache/v2_hpc_*.json
export RFC_TAG=${RFC_TAG:-_hpc}

echo "host        : $(hostname)"
echo "cores       : ${SLURM_CPUS_PER_TASK:-?}"
echo "RFC_QUICK   : ${RFC_QUICK}"
echo "python      : $(python -V 2>&1)"
echo "started     : $(date)"

# Affinity check, AFTER importing sklearn -- the import is what collapses it. This
# number must equal the core count above; if it reads 1, the job is about to run
# single-threaded and should be cancelled rather than left to burn a day.
python -c "
import os, sklearn
from joblib import cpu_count
print(f'usable cores: {len(os.sched_getaffinity(0))} (joblib sees {cpu_count()})')"
echo

cd notebooks            # relative paths in the script assume this directory
# -u: unbuffered, so the log shows progress while the job runs instead of staying
# empty until the process exits.
time python -u 03_16S_ITS1_RFC.py

echo
echo "finished    : $(date)"
python -c "import sklearn, shap, numpy, pandas, scipy, biom; print(
    'versions:', sklearn.__version__, shap.__version__, numpy.__version__,
    pandas.__version__, scipy.__version__, biom.__version__)"
