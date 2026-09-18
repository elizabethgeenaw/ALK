#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=150G
#SBATCH --job-name=conv
#SBATCH --output=/project/gilad/geena/ALK_analysis/scripts/control/slurm_output/conv-%j.out
#SBATCH --error=/project/gilad/geena/ALK_analysis/scripts/control/slurm_output/conv-%j.err

module load python
python /project/gilad/geena/ALK_analysis/scripts/control/1_ebQTL_controls_splitoriginal.py