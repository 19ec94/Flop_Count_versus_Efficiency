#include<stdlib.h>
#include <stdio.h>

double reallocate(double **input, int newsize) {
    *input = realloc( *input, newsize*sizeof(double) );
    if(*input!=NULL) {
        printf("successful\n");
    }
    else return(-1);
}

    