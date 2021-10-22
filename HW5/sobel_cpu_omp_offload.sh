#!/bin/bash

export LIBOMPTARGET_INFO=4
./sobel_cpu_omp_offload > ../logs/offload/sobel_cpu_omp_offload.txt 2>&1
