#!/bin/sh

sed -i "s/time=12/time=01/g" jobs/i16j*
sed -i "s/partition=medium/partition=short/g" jobs/i16j*

for file in jobs/*; do
    sbatch $file
done
