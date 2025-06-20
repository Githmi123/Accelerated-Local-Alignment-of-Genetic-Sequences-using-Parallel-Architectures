#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../constants.h"
#include <sys/time.h>
#include <math.h>
#include <cuda_runtime.h>

char seq1_list[MAX_PAIRS][MAX_SEQ_LENGTH];
char seq2_list[MAX_PAIRS][MAX_SEQ_LENGTH];
int score_matrix[MAX_PAIRS];

int load_sequences(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if ( file == NULL ) {
        perror("File Open Error");
        printf("Error opening file %s: %d", filename, errno);
        exit(1);
    }

    int count = 0;

    while (fscanf( file, "%200[^,],%200s\n", seq1_list[count], seq2_list[count]) == 2) {
        count++;
    }

    fclose(file);
    return count;
}

__global__ void smith_waterman_kernel(char* d_seq1, char* d_seq2, int* d_offsets, int* d_scores)    
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int offset = d_offsets[idx];
    char* s1_pointer = &d_seq1[offset * MAX_SEQ_LENGTH];
    char* s2_pointer = &d_seq2[offset * MAX_SEQ_LENGTH];

    int len1 = MAX_SEQ_LENGTH;
    int len2 = MAX_SEQ_LENGTH;

    // High Scoring Local Alignment Matrix (H)
    int H[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH] = {0};

    int score_diagonal, score_up, score_left, max_score = 0;

    for (int i = 1; i <= len1 ; i++)
    {
        for (int j = 1; j <= len2 ; j++)
        {
            score_diagonal = H[i - 1][j - 1] + (s1_pointer[i - 1] == s2_pointer[j - 1] ? MATCHING_SCORE : MISMATCHING_SCORE);
            score_up = H[i - 1][j] + GAP_PENALTY;
            score_left = H[i][j - 1] + GAP_PENALTY;
            H[i][j] = max(0, max(score_diagonal, max(score_up, score_left)));

            if (H[i][j] > max_score) 
                max_score = H[i][j];
        }
    }

    d_scores[idx] = max_score;
}

void save_score_matrix(const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        perror("File Open Error");
        printf("Error opening file %s: %d\n", filename, errno);
        exit(1);
    }

    for(int i = 0; i < MAX_PAIRS; i++)
    {
        fprintf(file, "Index %d, Score: %d\n", i, score_matrix[i]);
    }

    fclose(file);
    printf("Score matrix saved to %s\n", filename);
}

int main()
{
    printf("\n"
           "========================================\n"
           "Smith-Waterman Algorithm - Using CUDA\n"
           "========================================\n");
    int n = load_sequences("data/DNASequences.txt");
    printf("Loaded %d pairs of sequences.\n", n);

    char *h_seq1 = (char*) malloc(n * MAX_SEQ_LENGTH * sizeof(char));
    char *h_seq2 = (char*) malloc(n * MAX_SEQ_LENGTH * sizeof(char));
    int *h_offsets = (int*) malloc(n * sizeof(int));
    int *h_scores = (int*) malloc(n * sizeof(int));

    char *d_seq1, *d_seq2;
    int *d_offsets, *d_scores;
    cudaMalloc(&d_seq1, n * MAX_SEQ_LENGTH * sizeof(char));
    cudaMalloc(&d_seq2, n * MAX_SEQ_LENGTH * sizeof(char));
    cudaMalloc(&d_offsets, n * sizeof(int));
    cudaMalloc(&d_scores, n * sizeof(int));

    for ( int i = 0; i < n; i++ )
    {
        memcpy(&h_seq1[i * MAX_SEQ_LENGTH], seq1_list[i], MAX_SEQ_LENGTH);
        memcpy(&h_seq2[i * MAX_SEQ_LENGTH], seq2_list[i], MAX_SEQ_LENGTH);
        h_offsets[i] = i;
    }

    cudaMemcpy(d_seq1, h_seq1, n * MAX_SEQ_LENGTH, cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, h_seq2, n * MAX_SEQ_LENGTH, cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, h_offsets, n * sizeof(int), cudaMemcpyHostToDevice);

    int threads_per_block = 256;
    int num_blocks = (n + threads_per_block - 1) / threads_per_block;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    smith_waterman_kernel<<<num_blocks, threads_per_block>>>(d_seq1, d_seq2, d_offsets, d_scores);

    gettimeofday(&end, NULL);

    printf("Completed alignment of %d sequence pairs.\n", n);
    
    long start_time = (start.tv_sec * 1000000 + start.tv_usec);
    long end_time = (end.tv_sec * 1000000 + end.tv_usec);
    long elapsed_time = end_time - start_time;
    printf("Total time taken: %0.6f seconds\n", (float)elapsed_time / 1000000);

    cudaMemcpy(h_scores, d_scores, n * sizeof(int), cudaMemcpyDeviceToHost);

    memcpy(score_matrix, h_scores, n * sizeof(int));

    cudaFree(d_seq1);
    cudaFree(d_seq2);
    cudaFree(d_offsets);
    cudaFree(d_scores);

    free(h_seq1);
    free(h_seq2);
    free(h_offsets);
    free(h_scores);

    save_score_matrix("output/phase02_code_output_max_scores.txt");
}