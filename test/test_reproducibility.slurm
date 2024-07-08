#!/bin/bash
#SBATCH --cluster=gmerlin6
#SBATCH --account=gwendolen
#SBATCH --partition=gwendolen
#SBATCH --gpus=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1
#SBATCH --time=00:29:59
#SBATCH --mem-per-gpu=32G
#SBATCH --output=../logs/test_reproducibility_%A.out
#SBATCH --error=../logs/test_reproducibility_%A.err

source ~/.bashrc

module use unstable; module load anaconda; module load Java/12u2; module load gcc
conda activate /data/user/bellotti_r/miniconda_env

patientID="test_00"
BA="clinical"

if [ "$BA" == "clinical" ]; then
    planFile=$DATA_DIR/clinical_plans_300_iter/"$patientID"_0.json
elif [ "$BA" == "predicted" ]; then
    planFile=/data/user/bellotti_r/feature_learning/gambas_output/plans_new/"$patientID".json
else
    echo "Unknown beam arrangement: $BA"
    exit 1;
fi


papermill -p patient_ID $patientID \
          -p data_dir $DATA_DIR \
          -p output_dir ../output/test_reproducibility1 \
          -p plan_file "$planFile" \
          -p dicom_output_dir ../output/test_reproducibility1_DICOM \
          -p early_stopping_delta 1 \
          -p early_stopping_patience 1 \
          ../optimise_julia.ipynb \
          ../executed_notebooks/test_reproducibility1.ipynb

papermill -p patient_ID $patientID \
          -p data_dir $DATA_DIR \
          -p output_dir /data/user/bellotti_r/differentiable-planning/output/test_reproducibility2 \
          -p plan_file "$planFile" \
          -p dicom_output_dir /data/user/bellotti_r/differentiable-planning/output/test_reproducibility2_DICOM \
          -p early_stopping_delta 1 \
          -p early_stopping_patience 1 \
          ../optimise_julia.ipynb \
          ../executed_notebooks/test_reproducibility2.ipynb

