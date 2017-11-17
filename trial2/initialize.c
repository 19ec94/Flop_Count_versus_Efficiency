#include<stdlib.h>

double* initialize( int size){
      double *locptr;
      locptr = (double*)malloc(sizeof(double) * size);
      return locptr;
}
    