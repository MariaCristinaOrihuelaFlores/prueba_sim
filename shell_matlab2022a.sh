#!/bin/bash
#SBATCH --nodelist=worker9
#SBATCH --output="slurm-%j.out"
#SBATCH --time=200000
#SBATCH --partition=thinkstation-p360
#SBATCH --gres=gpu:1
srun /usr/local/MATLAB/R2023b/bin/matlab -nosplash -nodesktop -nodisplay -r "sim_attemp1; exit"