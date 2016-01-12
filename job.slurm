#!/bin/sh
#SBATCH --job-name=test
#SBATCH --time=00:05:00
#SBATCH --mail-user=corentin.cadiou@obspm.fr --mail-type=ALL
#SBATCH --partition=medium
## Pour reserver 1 noeud
#SBATCH --nodes=1

export OMP_NUM_THREADS=12
RUNDIR=/mnt/lnec/travail/$USER/run.${SLURM_JOBID}
mkdir ${RUNDIR}
cd ${RUNDIR}

cp -r /obs/$USER/n-body/* .

make clean
make gen
./gen
make
/usr/bin/time -o simul.log ./simul
