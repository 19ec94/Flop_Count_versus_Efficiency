#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>
#include "sys/time.h"
#include "time.h"
double* initialize(int size);
double* allocate(double *locptr,int size);
double* allocate_zero(double *locptr,int size);

extern double dummy_indv_duration[64][30],dummy_no_operation[64][30],indv_gflops[64];


void multi(double* A,double* B,double *inter,double **tempptr, int m,int n,int k,int tt, int itr){ 
    
  *tempptr = (double *) realloc( *tempptr, m*n*sizeof(double) );
  if (*tempptr !=NULL)
       inter = *tempptr;
  else 
       printf("reallocation failed \n");
  struct timeval start,finish;

    gettimeofday(&start, NULL);    //start timer                          
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,k,B,n,0.0,inter,n);  
    gettimeofday(&finish, NULL);   //end timer
    
    dummy_indv_duration[tt][itr] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
    dummy_no_operation[tt][itr] = 2.0 * m *n*k;
    

    return ;
}




//indv_gflops[tt] = no_operation[tt]/indv_duration[tt]*1.0e-12;
