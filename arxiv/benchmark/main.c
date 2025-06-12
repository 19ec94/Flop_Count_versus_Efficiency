/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A C Program to calculate all the possible ways of multiplying                                %
% 5 matrices                                                                                   %
%                                                                                              %
% Input : matrix sizes                                                                         %
% output: number of combinations, time , number of operation by each combination, total gflops %
%                                                                                              %
%                                                                                              %
% Authors: Edilbert Chrsithuraj, Sadulla Aghayev                                               %
%                                                                                              %
%                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */



#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"
#include "cblas.h"
#include <string.h>
#define min(x,y) (((x) < (y)) ? (x) : (y))

double* allocate(int size);                                                                                           //alloacte memory for matrices                                                                      //initialize matrices to zero
void initialize( double* locptr, int size);                                                                            //Initialize matrices to random values                              
void printfunc (double* matrix, int row, int col);                                                                    //print resultant matrix 

int main(int argc, char* argv[])
{
 
  int i,j;
  double total_duration =0.0,total_gflops=0.0,total_no_operation=0.0;                                                   //stores values for each paths
  int  m = 200,n= 200,k =200;
  int  d_m = 10000,d_n= 10000,d_k =10000;                                                              
  double *dummy_a=NULL,*dummy_b=NULL,*dummy_c=NULL,*A=NULL,*B=NULL,*C=NULL;
  struct timeval start,finish;
  
  dummy_a=allocate(m*n);
  dummy_b=allocate(n*k);
  dummy_c=allocate(k*n);
  A =allocate(d_m*d_n) ;
  B =allocate(d_n*d_k) ;
  C =allocate(d_m*d_k) ;
  srand((unsigned)time(NULL));
  initialize(dummy_a,m*n);
  initialize(dummy_b,n*k);
  initialize(dummy_c,m*k);
  initialize(A,d_m*d_n);
  initialize(B,d_n*d_k);
  initialize(C,d_m*d_k);
  
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,d_m,d_n,d_k,1.0,A,d_k,B,d_n,0.0,C,d_n); 
  
  gettimeofday(&start,NULL);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n); 
  gettimeofday(&finish,NULL);
   
  total_no_operation = 2*m*n*k;
  total_duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  total_gflops = total_no_operation/total_duration;

  printf("%lf \t %f \t %f \n",total_no_operation,total_duration,total_gflops);

  printf("Everything is successful\n");
  
  free( dummy_a);
  free (dummy_b);
  free (dummy_c);
  free (A);
  free (B);
  free (C);


return 0;

}   //end of main

double* allocate(int size){
  double* locptr = (double* )malloc(sizeof(double) * size);
  if (locptr == NULL){
    printf("Memory is not alloacted \n");
  }
  return locptr;
}

void initialize(double *locptr, int size){
  int i; 
  for (i=0; i<size; i++)
    locptr[i] = 1;
  return ;
}
void printfunc (double* matrix, int row, int col){
  int i, j;
  printf("____________________________________________________________\n");
  for (i=0; i<min(row,6); i++) {
      for (j=0; j<min(col,6); j++) {
        printf ("%.2f ", matrix[j+i*col ]);
      }
    printf ("\n");
  }
  printf("____________________________________________________________\n");
}


