#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../constants.h"
#include <sys/time.h>
#include <math.h>
#include <omp.h>

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

void smith_waterman(char* seq1_list, char* seq2_list, int index)
{
    int len1 = strlen(seq1_list);
    int len2 = strlen(seq2_list);

    // High Scoring Local Alignment Matrix (H)
    int** H = malloc((len1 + 1) * sizeof(int*));

    for (int i = 0; i <= len1 ; i++) // Can be parallelized
    {
        H[i] = malloc((len2 + 1) * sizeof(int));
        memset(H[i], 0, (len2 + 1) * sizeof(int)); // Initialize all rows to zero
    }

    int score_diagonal, score_up, score_left, max_score = 0;

    for (int i = 1; i <= len1 ; i++)
    {
        for (int j = 1; j <= len2 ; j++)
        {
            score_diagonal = H[i - 1][j - 1] + (seq1_list[i - 1] == seq2_list[j - 1] ? MATCHING_SCORE : MISMATCHING_SCORE);
            score_up = H[i - 1][j] + GAP_PENALTY;
            score_left = H[i][j - 1] + GAP_PENALTY;
            H[i][j] = fmax(0, fmax(score_diagonal, fmax(score_up, score_left)));

            if (H[i][j] > max_score) 
                max_score = H[i][j];
        }
    }

    score_matrix[index] = max_score;
    // printf("Max score for pair %d: %d\n", index, max_score);

    for ( int i = 0; i <= len1; i++ ) // Can be parallelized
    {
        free(H[i]);
    }
    free(H);
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
           "Smith-Waterman Algorithm - Using OpenMP\n"
           "========================================\n");
    int n = load_sequences("data/DNASequences.txt");
    printf("Loaded %d pairs of sequences.\n", n);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++)
    {
        smith_waterman(seq1_list[i], seq2_list[i], i);
    }

    gettimeofday(&end, NULL);

    printf("Completed alignment of %d sequence pairs.\n", n);
    
    long start_time = (start.tv_sec * 1000000 + start.tv_usec);
    long end_time = (end.tv_sec * 1000000 + end.tv_usec);
    long elapsed_time = end_time - start_time;
    printf("Total time taken: %0.6f seconds\n", (float)elapsed_time / 1000000);

    save_score_matrix("output/phase01_code_output_max_scores.txt");
}