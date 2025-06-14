#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../constants.h"

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
}
