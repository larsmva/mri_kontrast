#!/bin/bash
# Job name:
#SBATCH --job-name=kent
#
# Project:
#SBATCH --account=nn9279k
# Wall clock limit:
#SBATCH --time='06:00:00'
#
# Max memory usage per task:
#SBATCH --mem-per-cpu=3800M
#
# Number of tasks (cores):
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --partition=long
##SBATCH --output=output.$SCRATCH 
source /cluster/bin/jobsetup
module load freesurfer
#chkfile subjects 

export FREESURFER_HOME=$HOME/freesurfer 
export SUBJECTS_DIR=$SCRATCH
#export PATH=$FREESURFER_HOME/bin:$PATH
echo $SCRATCH
echo `date`

cp -r $FREESURFER_HOME/subjects/MRI_ERIKA $SCRATCH
ls $SCRATCH 
ls $SCRATCH/MRI_ERIKA/*

recon-all -subjid erika -i   $SCRATCH/MRI_ERIKA/Dicom/IM_0594  -all 


echo `date`

