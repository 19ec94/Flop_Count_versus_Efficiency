#include<stdlib.h>

void initialize(double *locptr, int size){
       int i; 
    for (i=0; i<size; i++)
    locptr[i] = i;
    return ;
}