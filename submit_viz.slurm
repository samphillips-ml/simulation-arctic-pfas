#!/bin/bash
#SBATCH --job-name=pfas_viz
#SBATCH --output=output/slurm_viz_%j.out
#SBATCH --error=output/slurm_viz_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sphill86@charlotte.edu

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Partition   : $SLURM_JOB_PARTITION"
echo "======================================================"

set -euo pipefail
cd $SLURM_SUBMIT_DIR
source .pfas-sim/bin/activate
export PYTHONPATH=$SLURM_SUBMIT_DIR${PYTHONPATH:+:$PYTHONPATH}

echo "Python : $(which python)"
echo "ffmpeg : $(which ffmpeg || echo 'NOT FOUND')"

python -u presentation_video.py output/simulation_pfoa_89.nc output/presentation.mp4

echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"