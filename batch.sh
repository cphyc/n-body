#!/bin/sh

mkdir -p tmp
for i in 1 2 3 4 5 6 7 8 9 10 11 12; do
    sed "s/OMP_NUM_THREADS=.*/OMP_NUM_THREADS=$i/g" job.slurm |
	sed "s/job-name=.*/job-name=batch-$i/g" > tmp/job$i.slurm
    sbatch tmp/job$i.slurm
done
