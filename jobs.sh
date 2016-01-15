#!/bin/sh

for i in 1 2 4 8 16; do
    for j in 1 2 4 8; do
        for k in $(seq 1 $((12/$j))); do
            echo "$i nodes, $j mpi/node ($(($i*$j)) total), $k OpenMP/mpi ($(($j*$k))/node, $(($i*$j*$k)) total)"
        done
    done
done
