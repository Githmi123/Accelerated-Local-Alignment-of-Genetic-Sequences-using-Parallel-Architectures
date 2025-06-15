# Accelerated-Local-Alignment-of-Genetic-Sequences-using-Parallel-Architectures
## Description
Genetic sequence alignment is a fundamental bioinformatics problem that is used to compare genetic sequences. This project considers the Smith Waterman algorithm’s dynamic programming matrix computation part and finding it’s maximum value and will be done in three approaches:
1.	Phase1 – Shared Memory Approach : Implementation using OpenMP to parallelize dynamic programming matrix computation.
2.	Phase2 – Distributed Approach : Implementation using CUDA for GPU acceleration of the matrix computation.
3.	Phase3 - Combining OpenMP and CUDA to maximize performance. The alignment task is to be divided across CPU threads using OpenMP and each thread will launch CUDA kernels to process it’s assigned data on the GPU, achieving multi-level parallelism.

## Objective
To implement the same Genetic Sequence Alignment problem in 3 different programming models. First would be to use shared memory on a single machine and the second would be to use a distributed memory model with GPU memory, and the third would be a hybrid model with multi-level parallelism, aiming for a significantly better performance than the first two. The final goal is to compare their performance and show that the hybrid version outperforms the other two.

## Technologies Used
- C 
- OpenMP (for CPU parallelism)  
- CUDA (for GPU parallelization)  
- Git for version control 