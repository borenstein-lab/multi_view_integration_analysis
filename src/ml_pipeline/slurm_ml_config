#!/bin/bash

# To execute: 
# cd /specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis/src/ml_pipeline
# sbatch slurm_ml_config
# Check queue: squeue

########################
# - GENERAL SETTINGS - #
########################

#SBATCH --job-name=ml_efrat
#SBATCH --partition=cpu-elbo 			
#SBATCH --ntasks=1 									         
#SBATCH --cpus-per-task=5                   
#SBATCH --mail-type=ALL,TIME_LIMIT_80        # Notification policy
#SBATCH --time=3-00:00 # Max running time 
#SBATCH --output=slurm_logs/slurm.%A.%a.out
#SBATCH --error=slurm_logs/slurm.%A.%a.err
#SBATCH --array=0-10

########################
# ------- JOBS ------- #
########################

DATASETS=(
crc_s3_s4_yachida_2019
cirrhosis_qin_2014
esrd_wang_2020
uc_franzosa_2019
cd_franzosa_2019
crc_feng_2015
crc_yu_2015
metacardis_1_8
metacardis_3_8
sth_rubel_2020
uc_spain_nielsen_2014
)

DS_NAME=${DATASETS[$SLURM_ARRAY_TASK_ID]}
echo "SLURM DEBUG: now working on dataset ${DS_NAME}"

BASE_DIR=/specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis
cd ${BASE_DIR}

srun udocker run \
  --volume=/specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis \
  efrat_ubun_r Rscript \
  /specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis/src/ml_pipeline/ml_pipeline.R \
  ${DS_NAME} 

echo "SLURM DEBUG: finished working on dataset ${DS_NAME}"
