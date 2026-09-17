#!/bin/bash
#PBS -N variant_filtering
#PBS -l select=1:ncpus=4:mem=150gb
#PBS -l walltime=15:00:00
#PBS -q cray
#PBS -o job_output.txt
#PBS -e job_error.txt

# Change to the directory where you submitted the job from
cd $PBS_O_WORKDIR

# Set numpy thread settings for reproducibility
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Set library path for Python
export LD_LIBRARY_PATH=/data/chengr/python-3.10.4/lib:$LD_LIBRARY_PATH

# Activate your virtualenv (FULL PATH)
source /home/chengr/env/bin/activate

# Print memory info at start
echo "=========================================="
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Memory allocated: 150GB"
echo "=========================================="

# Run Python script
python /home/chengr/Capstone/variant_effect_prediction/Variant-Effect-Prediction/phyloP_distribution.py
echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="