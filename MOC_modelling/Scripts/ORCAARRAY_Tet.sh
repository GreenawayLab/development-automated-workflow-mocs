#!/bin/bash
#SBATCH -J ORCAARRAY_Tet
#SBATCH -A NITSCHKE-SL3-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time 1:00:00


#SBATCH --array=0-7 

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
workdir="$SLURM_SUBMIT_DIR"

source /home/pcpt3/.bashrc
conda init
conda activate env3-p11
python --version

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment
module load python3
module load openmpi/4.1.1/gcc/mnop75he

export PATH=/rds/user/pcpt3/hpc-work/orca:$PATH
export LD_LIBRARY_PATH=/rds/user/pcpt3/hpc-work/orca:$LD_LIBRARY_PATH

export OPENMPI_HOME=/usr/local/software/spack/spack-views/rhel8-icelake-20211027_2/openmpi-4.1.1/gcc-11.2.0/mnop75hedc4hru2c22d6oqlu7opdc4jg 
export PATH=$OPENMPI_HOME/bin:$PATH
export LD_LIBRARY_PATH=$OPENMPI_HOME/lib:$LD_LIBRARY_PATH

Aldehydes_Br=("ald1" "ald2" "ald3" "ald4")
linker=("l0" "l1")

cage_names=()
t=0
#Define cages to iterate over by ORCA array
for a in "${Aldehydes_Br[@]}"; do
    for l in "${linker[@]}"; do
        cage_names[t]="Tetrahedron_${l}_${a}"
        (( t++ ))
    done
done

echo cage_names

#Create folder for ORCAjobs
export SLURM_SUBMIT_DIR
export SLURM_ARRAY_TASK_ID

i=$SLURM_ARRAY_TASK_ID

# Setup ORCAjob including loop variable
export CAGE_NAME="${cage_names[i]}"
echo "$CAGE_NAME"

dir_path="/rds/user/pcpt3/hpc-work/CageConstruction/GreenawayCollab/Tetrahedrons_ORCA_DFT/${cage_names[i]}_ORCA_DFT"
mkdir -p "$dir_path"
cd "$dir_path"
   
# Creating local scratch folder 
export scratchlocation=/rds/user/pcpt3/hpc-work/Orca-tmp
tdir=$(mktemp -d $scratchlocation/ORCAjob__$SLURM_JOB_ID-XXXX-${cage_names[i]}_ORCA_DFT)

# Copy only the necessary stuff in submit directory to scratch directory. Add more (e.g. .gbw) here if needed.
cp  $dir_path/*.inp $tdir/
cp  $dir_path/*.xyz $tdir/

# Creating nodefile in scratch
echo $SLURM_NODELIST > $tdir/$job.nodes

# cd to scratch
cd $tdir

rm -f ${dir_path}/${cage_names[i]}_ORCA_DFT.out
echo "Job execution start: $(date)" >>  ${dir_path}/${cage_names[i]}_ORCA_DFT.out
echo "Shared library path: $LD_LIBRARY_PATH" >>  ${dir_path}/${cage_names[i]}_ORCA_DFT.out
echo "Slurm Job ID is: ${SLURM_JOB_ID}" >>  ${dir_path}/${cage_names[i]}_ORCA_DFT.out
echo "Slurm Job name is: ${SLURM_JOB_NAME}" >>  ${dir_path}/${cage_names[i]}_ORCA_DFT.out
echo $SLURM_NODELIST >> ${dir_path}/${cage_names[i]}_ORCA_DFT.out
echo "Now in directory: $(pwd)"

# Execute the application
/rds/user/pcpt3/hpc-work/orca/orca ${cage_names[i]}_ORCA_DFT.inp > ${cage_names[i]}_ORCA_DFT.out
if [ $? -eq 0 ]; then
   echo "Single-point energy calculation with ORCA executed successfully for ${cage_names[i]}"
   # ORCA has finished here. Now copy important stuff back. Add more here if needed.
   cp $tdir/*.{xyz,gbw,out,nodes,final,system,new,info,mol} $dir_path

else
   echo "Single-point energy calculation with ORCA failed for ${cage_names[i]}, check ${cage_names[i]}_ORCA_DFT.out for details"
   cp $tdir/*.{xyz,gbw,out,nodes,final,system,new,info,mol} $dir_path
fi

rm -rf $tdir

