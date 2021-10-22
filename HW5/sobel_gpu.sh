#!/bin/bash

for b in 1 4 16 64 256 1024 4096
do
    for t in 32 64 128 256 512 1024
    do
	export NUM_BLOCKS=$b
	export THREADS_PER_BLOCK=256
	echo "Running B=${b} T=${t}"
	nvprof ./sobel_gpu > ../sobel_gpu_${b}_${t}.txt 2>&1
    done
done
