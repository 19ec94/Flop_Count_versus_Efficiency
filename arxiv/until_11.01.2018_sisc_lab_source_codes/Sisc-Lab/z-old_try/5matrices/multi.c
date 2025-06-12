#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>
#include "sys/time.h"
#include "time.h"
double* initialize(int size);
double* allocate(double *locptr,int size);
double* allocate_zero(double *locptr,int size);

extern double indv_duration[64],no_operation[64],indv_gflops[64];


void multi(double* A,double* B,double *inter1,double **tempptr, int m,int n,int k,int tt){ 
    *tempptr = (double *) realloc( *tempptr, m*n*sizeof(double) );
    if (*tempptr !=NULL)
        inter1 = *tempptr;
     else 
        printf("reallocation failed \n");
    struct timeval start,finish;

    gettimeofday(&start, NULL);  
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,k,B,n,0.0,inter1,n);  
    gettimeofday(&finish, NULL);
    
    indv_duration[tt] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
    no_operation[tt] = 2.0 * m *n*k;
    indv_gflops[tt] = no_operation[tt]/indv_duration[tt]*1.0e-12;

    return ;
}
