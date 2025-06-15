#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../constants.h"
#include <sys/time.h>

char seq1_list[MAX_PAIRS][MAX_SEQ_LENGTH];
char seq2_list[MAX_PAIRS][MAX_SEQ_LENGTH];

int load_sequences(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if ( file == NULL ) {
        perror("File Open Error");
        printf("Error opening file %s: %d", filename, errno);
        return -1;  
    }

    int count = 0;

    while (fscanf( file, "%200[^,],%200s\n", seq1_list[count], seq2_list[count]) == 2) {
        count++;
    }

    fclose(file);
    return count;
}

int main()
{
    int n = load_sequences("../data/DNASequences.txt");
    printf("Loaded %d pairs of sequences.\n", n);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (int i = 0; i < n; i++)
    {
        // Process
    }

    gettimeofday(&end, NULL);

    printf("Completed alignment of %d sequence pairs.\n", n);
    printf("Time taken in seconds: %ld seconds\n", end.tv_sec - start.tv_sec);
    printf("Time taken in microseconds: %ld microseconds\n", end.tv_usec - start.tv_usec);
    printf("Total time taken: %0.7f seconds\n", 
           (float)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000));
}
