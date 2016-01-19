#!/bin/sh
#SBATCH --job-name=<JOB_NAME>
#SBATCH --time=12:00:00
#SBATCH --mail-user="corentin.cadiou@obspm.fr,bruno.pagani@obspm.fr" --mail-type=ALL
#SBATCH --partition=medium
#SBATCH --nodes=<NUM_NODES>

export OMP_NUM_THREADS=<OMP_NUM_THREADS>
export MPI_PROCS=<MPI_PROCS>
export PROCS_PER_NODE=<PROCS_PER_NODE>
export FLAG_MEMORY=<FLAG_MEMORY>
export FLAG_DIAG=<FLAG_DIAG>

# Get a fresh directory
RUNDIR=/mnt/lnec/travail/$USER/run.${SLURM_JOBID}
mkdir ${RUNDIR}
cd ${RUNDIR}
#GIT_SSL_NO_VERIFY=true git clone http://github.com/cphyc/n-body.git
cp /obs/$USER/n-body/*.f90 .
cp /obs/$USER/n-body/Makefile .

# Get informations for mpi
scontrol show hostname $SLURM_JOB_NODELIST > /tmp/hosts.$USER.${SLURM_JOBID}
export LD_LIBRARY_PATH=/opt/openmpi/1.4.3-icc/lib:/opt/intel/composerxe-2011.4.191/compiler/lib/intel64/

# Replace the
sed -i "s/flag_memory = .*/flag_memory = $FLAG_MEMORY/g" constants.f90
sed -i "s/flag_diag = .*/flag_diag = $FLAG_DIAG/g" constants.f90

# Compile
make clean
make gen
make

# Run initial conditions
mpirun --hostfile /tmp/hosts.$USER.${SLURM_JOBID} --n 1 --npernode 1 ./gen

# Print initial parameters
echo "${SLURM_NNODES} nodes, $PROCS_PER_NODE mpi/node ($MPI_PROCS total), $OMP_NUM_THREADS OpenMP/mpi ($PROCS_PER_NODE*$OMP_NUM_THREADS/node, $MPI_PROCS*$OMP_NUM_THREADS total)" > params.log

# Run the simulation
/usr/bin/time -f "Total: %es User: %Us System: %Ss CPU: %P" -o simul.log \
	      mpirun --hostfile /tmp/hosts.$USER.${SLURM_JOBID} --n $MPI_PROCS --npernode $PROCS_PER_NODE ./simul
