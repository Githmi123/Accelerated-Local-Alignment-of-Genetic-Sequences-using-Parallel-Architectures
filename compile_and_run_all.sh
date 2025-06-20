#!/bin/bash

# Compile serial code
echo "Compiling serial code..."
gcc -o serial_code/serial_code serial_code/serial_code.c -lm

# Compile OpenMP code
echo "Compiling OpenMP code..."
gcc -fopenmp -o phase01/using_openmp phase01/using_openmp.c -lm

# Compile CUDA code
echo "Compiling CUDA code..."
nvcc -o phase02/cuda_app phase02/using_CUDA.cu -lm

# Run all
echo "Running serial code..."
./serial_code/serial_code

echo "Running OpenMP code..."
./phase01/using_openmp

echo "Running CUDA code..."
./phase02/cuda_app
