#!/bin/bash
#SBATCH --job-name=rfc_lofo
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 32
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --output=logs/rfc/03b_lofo_%j.out
#SBATCH --error=logs/rfc/03b_lofo_%j.err
#SBATCH --mail-type=END,FAIL

# Random Forest, leave-one-farm-out (notebooks/03b_16S_ITS1_RFC_lofo.py).
# Separate from the notebook 03 RFC: outputs go to results/<version>/rfc_lofo/ only.
# Quick laptop test (LOFO_QUICK=1): 2.5 min on 6 cores. Full run: 6 x (7 x 73 fits
# at 500 trees + 1000 permutations x 7 fits at 200 trees); expected well under 12 h.
#
# Submit, e.g.: CFM_SAMPLE_SET=main sbatch scripts/03b_hpc_rfc_lofo.sh

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
export LOFO_QUICK=${LOFO_QUICK:-0}
echo "host: $(hostname) | cores: ${SLURM_CPUS_PER_TASK:-?} | data: $CFM_VERSION, $CFM_SAMPLE_SET | quick: $LOFO_QUICK"

# Must equal the core count; 1 means the job would run single-threaded
python -c "
import os, sklearn
from joblib import cpu_count
print(f'usable cores: {len(os.sched_getaffinity(0))} (joblib sees {cpu_count()})')"

cd notebooks    # relative paths in the script assume this directory
time python -u 03b_16S_ITS1_RFC_lofo.py
