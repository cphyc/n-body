#!/bin/sh
#SBATCH --job-name=test
#SBATCH --time=00:05:00
#SBATCH --mail-user=corentin.cadiou@obspm.fr --mail-type=ALL
#SBATCH --partition=medium
## Pour reserver 1 noeud
#SBATCH --nodes=2

export OMP_NUM_THREADS=12
RUNDIR=/mnt/lnec/travail/$USER/run.${SLURM_JOBID}
mkdir ${RUNDIR}
cd ${RUNDIR}
scontrol show hostname $SLURM_JOB_NODELIST > /tmp/hosts.ccadiou.${SLURM_JOBID}

cp -r /obs/$USER/n-body/* .

make clean
make gen
make

export LD_LIBRARY_PATH=/opt/openmpi/1.4.3-icc/lib:/opt/intel/composerxe-2011.4.191/compiler/lib/intel64/

mpirun --hostfile /tmp/hosts.ccadiou.${SLURM_JOBID} --n 1 --npernode 1 ./gen

/usr/bin/time -f "Total: %es User: %Us System: %Ss CPU: %P" -o simul.log \
	      mpirun --hostfile /tmp/hosts.ccadiou.${SLURM_JOBID} --n 2 --npernode 1 ./simul
