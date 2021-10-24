#!/bin/bash -l
#SBATCH -N 1
#SBATCH -c 10
#SBATCH -G 1
#SBATCH -C gpu
#SBATCH -t 5:00
#SBATCH -J queue
#SBATCH --job-name=gpu-job
#SBATCH --output=gpu-job.o%j
#SBATCH --error=gpu-job.e%j

export LIBOMPTARGET_INFO=4
#
# note: you will need to modify below here to launch your specific program
# it is assumed your environment is set up properly for using the Cori GPUs
# prior to you launching this batch script
#
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

# nvprof -m sm_efficiency ./sobel_gpu

# nvprof -m sm_efficiency ./sobel_cpu_omp_offload


