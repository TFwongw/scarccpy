#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=10
#SBATCH --mem=80g
#SBATCH --tmp=25g
#SBATCH --mail-type=NONE

module load python3
module load java/openjdk-13.0.2
module load gurobi/9.0.2
source ~/vscar/bin/activate
unset PYTHONPATH
unset PYTHONHOME
unset PYTHONSTARTUP

python generate_and_run.py