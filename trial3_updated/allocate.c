#include<stdlib.h>
#include <stdio.h>

double* allocate(int size){
    
    double* locptr = (double* )malloc(sizeof(double) * size);
      if (locptr == NULL){
           printf("Memory is not alloacted \n");
      }
      return locptr;
 }