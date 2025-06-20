#!/bin/bash

# Run all
echo "Running serial code..."
./serial_code/serial_code

echo "Running OpenMP code..."
./phase01/using_openmp

echo "Running CUDA code..."
./phase02/cuda_app
