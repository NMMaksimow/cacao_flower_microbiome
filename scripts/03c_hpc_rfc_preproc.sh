#!/bin/bash
#SBATCH --job-name=rfc_preproc
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 32
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH --output=logs/rfc/03c_preproc_%j.out
#SBATCH --error=logs/rfc/03c_preproc_%j.err
#SBATCH --mail-type=END,FAIL

# Random Forest preprocessing check (notebooks/03c_16S_ITS1_RFC_preproc.py):
# feature representations, prevalence thresholds and a shuffled-label null, each
# scored by 5-fold CV x 5 repeats.
#
# Inputs besides the git checkout (rsync them first, results/ is not in git):
#   results/<version>/biom/          biom tables, as for notebook 03
#   results/<version>/tables/clr/    CLR matrices of scripts/08c_16S_ITS1_clr_export.R
# Outputs: results/<version>/tables/03c_hpc_{scores,null,diagnostics}_<set>.csv,
# written row by row; a resubmitted job skips the rows already there.
# Figures are drawn on the laptop by notebooks/03c_16S_ITS1_RFC_preproc_vis.py.
#
# Size of the full run, at 500 trees: 144 scored settings x 25 fits = 3 600 fits, plus
# the null, 2 markers x 2 targets x 100 shuffles x 25 fits = 10 000 smaller fits (63
# or 100 flowers). Laptop timing (~11 core-seconds per fit on 288 flowers) suggests
# well under 2 h on 32 cores.
# RFC_PREPROC_QUICK=1 runs every code path with 100 trees, 2 repeats, 5 shuffles.
#
# Submit, e.g.:
#   RFC_PREPROC_QUICK=1 sbatch scripts/03c_hpc_rfc_preproc.sh    # test pass first
#   sbatch scripts/03c_hpc_rfc_preproc.sh

module load python/3.12.5
source ~/nobackup/venvs/rfc/bin/activate
cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/rfc

# Same thread settings as 03_hpc_rfc.sh: libraries stay inside the allocation, and
# OMP_PROC_BIND=false stops the sklearn import from collapsing the job to one core.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export LOKY_MAX_CPU_COUNT=${SLURM_CPUS_PER_TASK:-1}
export OMP_PROC_BIND=false
export MPLBACKEND=Agg

export CFM_VERSION=${CFM_VERSION:-v2}
export CFM_SAMPLE_SET=${CFM_SAMPLE_SET:-main}
export RFC_PREPROC_QUICK=${RFC_PREPROC_QUICK:-0}
export RFC_TAG=${RFC_TAG:-_hpc}
echo "host: $(hostname) | cores: ${SLURM_CPUS_PER_TASK:-?} | data: $CFM_VERSION, $CFM_SAMPLE_SET | quick: $RFC_PREPROC_QUICK"

# Must equal the core count; 1 means the job would run single-threaded
python -c "
import os, sklearn
from joblib import cpu_count
print(f'usable cores: {len(os.sched_getaffinity(0))} (joblib sees {cpu_count()})')"

cd notebooks    # relative paths and the utils_rfc import assume this directory
time python -u 03c_16S_ITS1_RFC_preproc.py
