#!/bin/sh

for i in $(seq 1 12); do
    for k in 1 2 3 4 6 12; do
        for j in $(seq 1 $k); do
            for l in $(seq 1 $((12/$k))); do
                echo "$i noeuds, $j mpi, $l OpenMP"
            done
        done
    done
done
