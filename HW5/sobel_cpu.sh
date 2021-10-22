#!/bin/bash

for N in 1 2 4 8 12 16
do
    echo "Running OMP_NUM_THREADS=$N"
    export OMP_NUM_THREADS=$N
    export OMP_PLACES=threads
    export OMP_PROC_BIND=spread
    ./sobel_cpu > ../logs/cpu/sobel_cpu_$N.txt 2>&1
done
