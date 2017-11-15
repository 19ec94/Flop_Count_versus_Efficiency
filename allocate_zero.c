#include<stdlib.h>
double* allocate_zero(double *locptr,int size){
    int i; 
    for (i=0; i<size; i++)
    locptr[i] = 0;
    return locptr;
 }