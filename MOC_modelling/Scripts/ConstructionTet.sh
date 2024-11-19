#!/bin/bash
#SBATCH -J ConstructTet
#SBATCH -A NITSCHKE-SL3-CPU
#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=12:00:00

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

source /home/pcpt3/.bashrc
module purge
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module load rhel8/default-icl              # REQUIRED - loads the basic environment

# Python
conda init
conda activate env3-p11
export PYTHONPATH="/home/pcpt3/PythonCalculations/Tetrahedrons"

python --version

# Run python code
python ConstructionTet.py

