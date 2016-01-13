#!/bin/sh

mkdir -p tmp
for i in {1..12}; do
    sed "s/OMP_NUM_THREADS=.*/OMP_NUM_THREADS=$i/g" job.slurm |
	sed "s/job-name=.*/job-name=batch-$i/g" > tmp/job$i.slurm
    sbatch tmp/job$i.slurm
done
