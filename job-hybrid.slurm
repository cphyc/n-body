#!/bin/sh
#SBATCH --job-name=<JOB_NAME>
#SBATCH --time=00:05:00
#SBATCH --mail-user=corentin.cadiou@obspm.fr --mail-type=ALL
#SBATCH --partition=medium
#SBATCH --nodes=<NUM_NODES>

export OMP_NUM_THREADS=<OMP_NUM_THREADS>
export MPI_PROCS=<MPI_PROCS>
export PROCS_PER_NODE=<PROCS_PER_NODE>
export FLAG_MPI=<FLAG_MPI>
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
cat constants.f90 |
    sed "s/flag_mpi = .*/flag_mpi = $FLAG_MPI/g" |
    sed "s/flag_diag = .*/flag_diag = $FLAG_DIAG/g" | tee constants.f90

# Compile
make clean
make gen
make

# Run initial conditions
mpirun --hostfile /tmp/hosts.$USER.${SLURM_JOBID} --n 1 --npernode 1 ./gen

# Run the simulation
/usr/bin/time -f "Total: %es User: %Us System: %Ss CPU: %P" -o simul.log \
	      mpirun --hostfile /tmp/hosts.$USER.${SLURM_JOBID} --n $MPI_PROCS --npernode $PROCS_PER_NODE ./simul
