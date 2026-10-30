#!/bin/bash
# Spec 09 A3/A5: one array task per footprint chunk (row of the chunk list written by scripts/08c_make_chunks.R);
# rendered by scripts/slurm_submit.R --what=chunks. The job's CPUs are the workers of the chunk. A long list is split
# into several arrays (MaxArraySize): task i of a part runs chunk row ROW_OFFSET + i.
#SBATCH --job-name=pblh_ffp
#SBATCH --partition={{SLURM_PARTITION}}
#SBATCH --array=1-{{N_TASKS}}{{ARRAY_MAX}}
#SBATCH --cpus-per-task={{CHUNK_CPUS_PER_TASK}}
#SBATCH --mem-per-cpu={{CHUNK_MEM_PER_CPU}}
#SBATCH --time={{SLURM_TIME}}
#SBATCH --mail-user={{SLURM_MAIL_USER}}
#SBATCH --mail-type={{SLURM_MAIL_TYPE}}
#SBATCH --output={{LOG_DIR}}/%x_%A_%a.out

set -euo pipefail
{{R_SETUP}}
export R_CONFIG_ACTIVE={{R_CONFIG_ACTIVE}}
cd "{{PROJECT_DIR}}"
ROW=$(( SLURM_ARRAY_TASK_ID + {{ROW_OFFSET}} ))
echo "chunk row ${ROW} on $(hostname), ${SLURM_CPUS_PER_TASK} CPUs"
Rscript scripts/08b_footprint_run.R --chunk-list="{{LIST_FILE}}" --chunk-row="${ROW}" --workers="${SLURM_CPUS_PER_TASK}" {{EXTRA_ARGS}}
