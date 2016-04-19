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
chkfile subjects erika 

#export PATH=$FREESURFER_HOME/bin:$PATH
echo $SCRATCH
echo `date`

cp -r $HOME/erika-T1 $SCRATCH
cp -r $HOME/erika-T2 $SCRATCH
ls $SCRATCH 

recon-all -subjid erika -sd $SCRATCH -i $SCRATCH/erika-T1/IM_0209 -T2 $SCRATCH/erika-T2/IM_0594 -T2pial -all 


echo `date`

