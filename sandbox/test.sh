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
##SBATCH --nodes=1 --ntasks=15
#SBATCH --ntasks=6
##SBATCH --hint=compute_bound
#SBATCH --cpus-per-task=1

#SBATCH --partition=long
#SBATCH --output=output.$SCRATCH 

## Set up job environment
source /cluster/bin/jobsetup

echo $1 $2

#module load gcc/4.9.2
#module load openmpi.gnu/1.8.4
#source ~oyvinev/intro/hashstack/fenics-1.5.0.abel.gnu.conf
source ~oyvinev/fenics1.6/fenics1.6

# Expand pythonpath with locally installed packages
export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python2.7/site-packages/

# Define what to do when job is finished (or crashes)
chkfile U_ref*

echo "SCRATCH is $SCRATCH"
# Copy necessary files to $SCRATCH
cp test.py pial_mesh.h5 pial_mesh.xdmf $SCRATCH
cd $SCRATCH
ls
echo $SCRATCH
mpirun --bind-to none python test.py No_refinements=$1 dt_val=$2 

cp $SCRATCH/*xdmf $PWD
cp $SCRATCH/*h5 $PWD
