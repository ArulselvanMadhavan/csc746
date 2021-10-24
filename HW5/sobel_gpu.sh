#!/bin/bash
export LIBOMPTARGET_INFO=4
for b in 1 4 16 64 256 1024 4096
do
    for t in 32 64 128 256 512 1024
    do
	export NUM_BLOCKS=$b
	export THREADS_PER_BLOCK=256
	echo "Running B=${b} T=${t}"
	nvprof -m sm_efficiency --csv --log-file ../logs/gpu/eff/sobel_gpu_eff_${b}_${t}.txt ./sobel_gpu
	nvprof --csv --log-file ../logs/gpu/exec/sobel_gpu_exec_${b}_${t}.txt ./sobel_gpu
    done
done
