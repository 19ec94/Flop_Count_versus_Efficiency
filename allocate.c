#include<stdlib.h>

double* allocate(double *locptr,int size){
    int i; 
    for (i=0; i<size; i++)
    locptr[i] = i;
    return locptr;
 }