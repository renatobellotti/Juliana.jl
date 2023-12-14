#!/bin/bash
#SBATCH --cluster=gmerlin6
#SBATCH --account=gwendolen
#SBATCH --partition=gwendolen
#SBATCH --gpus=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1
#SBATCH --time=00:29:59
#SBATCH --mem-per-gpu=32G
#SBATCH --output=logs/juliana_%A_%a.out
#SBATCH --error=logs/juliana_%A_%a.err

source ~/.bashrc

module use unstable; module load anaconda; module load Java/12u2; module load gcc
conda activate /data/user/bellotti_r/pyftpp_reproduction_env

printf -v patientID "test_%02d" $SLURM_ARRAY_TASK_ID

planFile=$DATA_DIR/clinical_plans_300_iter/"$patientID"_0.json


papermill -p patient_ID $patientID \
          -p data_dir $DATA_DIR \
          -p output_dir output/"$patientID" \
          -p plan_file "$planFile" \
          optimise_julia.ipynb \
          executed_notebooks/optimise_julia_"$patientID"_"$BA"_BA.ipynb
