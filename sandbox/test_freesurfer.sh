#!/bin/bash
# Job name:
#SBATCH --job-name=kent
#
# Project:
#SBATCH --account=nn9279k
# Wall clock limit:
#SBATCH --time='01:00:00'
#
# Max memory usage per task:
#SBATCH --mem-per-cpu=3800M
#
# Number of tasks (cores):
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --partition=long
#SBATCH --output=output.$SCRATCH 


export FREESURFER_HOME=$HOME/freesurfer 
export PATH=$FREESURFER_HOME/bin:$PATH
export SUBJECTS_DIR=$FREESURFER_HOME/subjects

recon-all -subjid erika -i  $SUBJECTS_DIR/MRI_ERIKA  



