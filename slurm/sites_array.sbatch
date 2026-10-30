#!/bin/bash
# Spec 09 A3: one array task per site (row of the site list); rendered by scripts/slurm_submit.R --what=sites.
# Placeholders (double braces) come from slurm/local.env (machine data, untracked) and from slurm_submit.R.
#SBATCH --job-name=pblh_sites
#SBATCH --partition={{SLURM_PARTITION}}
#SBATCH --array=1-{{N_TASKS}}{{ARRAY_MAX}}
#SBATCH --cpus-per-task={{SITES_CPUS_PER_TASK}}
#SBATCH --mem={{SITES_MEM}}
#SBATCH --time={{SLURM_TIME}}
#SBATCH --mail-user={{SLURM_MAIL_USER}}
#SBATCH --mail-type={{SLURM_MAIL_TYPE}}
#SBATCH --output={{LOG_DIR}}/%x_%A_%a.out

set -euo pipefail
{{R_SETUP}}
export R_CONFIG_ACTIVE={{R_CONFIG_ACTIVE}}
cd "{{PROJECT_DIR}}"
SITE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "{{LIST_FILE}}")
echo "task ${SLURM_ARRAY_TASK_ID}: ${SITE} on $(hostname)"
Rscript scripts/run_site.R "${SITE}" --steps={{STEPS}} {{EXTRA_ARGS}}
