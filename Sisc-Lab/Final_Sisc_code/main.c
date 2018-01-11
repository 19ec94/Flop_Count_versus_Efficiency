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
#define OUTPUT_in_FILE 1    //change it to 1, if output needs be written in a file
#define OUTPUT_on_SCREEN 1  //change it to Zero (0) , if you don't want to print it on screen
#define MATRIX 0            //change it to 1, if you want to check for the final matrix of each path
#define global_iteration 5
#define INTERMEDIATE_TIMINGS 0

double* allocate(int size);                                                                                           //alloacte memory for matrices
void initialize_zero(double *locptr,int size);                                                                        //initialize matrices to zero
void initialize( double* locptr, int size);                  //Initialize matrices to random values                              
//void multi(double* A,double*B, double*inter1, double **tempptr,int m, int n, int k,int tt,int itr );                          //matrix multiplication  
double* reallocate(double** tempptr, int m, int n);
void printfunc (double* matrix, int row, int col);                                                                    //print resultant matrix 
void sorting (double* total);

static int mat_size[7];                                                                                                // stores the matrix sizes
double dummy_indv_duration[64][30],indv_gflops[64],dummy_no_operation[64],indv_duration[64],no_operation[64];


int main(int argc, char* argv[])
{
  int i,j,outer_itr, my_iteration=1;
  double total_duration[14],total_gflops[14],total_no_operation[14];                                                   //stores values for each paths
  char *all_paths[] ={"((((A*B)*C)*D)*E)","((A*(B*(C *D)))*E)","(A*(B*(C*(D*E))))","(A*(((B*C)*D)*E))",
                 "(A*((B*C)*(D*E)))","(A*(B*((C*D)*E)))","(A*((B*(C*D))*E))","((A*B)(C*(D*E)))",
                 "((A*B)((C*D)*E))","(((A*B)*C)(D*E))","((A*(B*C))*(D*E))","(((A*(B*C))*D)*E)","((A*((B*C)*D))*E)","(((A*B)*(C*D))*E)"};
  char time_path[20],flop_path[20];
  initialize_zero(total_duration,14);
  initialize_zero(total_gflops,14);
 
 
 int outer_time_index[global_iteration];
 int outer_flop_index[global_iteration];
 double outer_min_time[global_iteration];
 double outer__min_flops[global_iteration];
 double outer_deviation[global_iteration];
 double total_duration_1[60],final_duration[60],total_gflops_1[60],no_operation[60];
 int duration_index[20],gflops_index[20];  //dummy variables for sorting; 


 
 printf("test!\n");
  if(argc<7){
    printf("Input Error\n");
    return 1;
  }
 
   for (i=0; i<7; i++)
    mat_size[i]=0;

  mat_size[0] = atoi(argv[1]);                                                                //row of matrix A and row of result matrix =m
  mat_size[1] = atoi(argv[2]);                                                                // column of matrix A and row of matrix B = k
  mat_size[2] = atoi(argv[3]);                                                                // column of matrix B and row of matrix C =n
  mat_size[3] = atoi(argv[4]);                                                                // column of matrix C and row of matrix D
  mat_size[4] = atoi(argv[5]);                                                                // column of matrix D and row of matrix E
  mat_size[5] = atoi(argv[6]);                                                                // column of matrix E and column of result matrix 
  int  m = 2000,n= 2000,k =2000;                                                               //for dummy matrices 
 
  
  double *A,*B,*C,*D,*E,*inter1=NULL,*inter2=NULL,*inter3=NULL,*inter4=NULL,*dummy_a,*dummy_b,*dummy_c;
   struct timeval start,finish;

  A=allocate(mat_size[0]*mat_size[1]);
  B=allocate(mat_size[1]*mat_size[2]);
  C=allocate(mat_size[2]*mat_size[3]);
  D=allocate(mat_size[3]*mat_size[4]);
  E=allocate(mat_size[4]*mat_size[5]);
  dummy_a=allocate(m*n);
  dummy_b=allocate(n*k);
  dummy_c=allocate(k*n);
  
                                                                                                //initialize matrices with random values
  srand((unsigned)time(NULL));
  initialize(A,mat_size[0]*mat_size[1]);
  initialize(B,mat_size[1]*mat_size[2]);
  initialize(C,mat_size[2]*mat_size[3]);
  initialize(D,mat_size[3]*mat_size[4]);
  initialize(E,mat_size[4]*mat_size[5]);
  initialize(dummy_a,m*n);
  initialize(dummy_b,n*k);
  initialize(dummy_c,k*m);


double inter_min_time ;


for (outer_itr=0; outer_itr<global_iteration; outer_itr++ )
{
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n); 
  //Path 0 - ((((A*B)*C)*D)*E)
  inter1=reallocate(&inter1,mat_size[0],mat_size[2]);
  inter2=reallocate(&inter2,mat_size[0],mat_size[3]);
  inter3=reallocate(&inter3,mat_size[0],mat_size[4]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);
   
  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n); 
 
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[2],mat_size[1],1.0,A,mat_size[1],B,mat_size[2],0.0,inter1,mat_size[2]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[3],mat_size[2],1.0,inter1,mat_size[2],C,mat_size[3],0.0,inter2,mat_size[3]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[4],mat_size[3],1.0,inter2,mat_size[3],D,mat_size[4],0.0,inter3,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[4],1.0,inter3,mat_size[4],E,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 

  dummy_indv_duration[0][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[0] = 2.0 *( (mat_size[0]*mat_size[2]*mat_size[1]) + (mat_size[0]*mat_size[3]*mat_size[2]) );
  dummy_no_operation[0] = dummy_no_operation[0]+ 2.0 *( ( mat_size[0]*mat_size[4]*mat_size[3] )+ (mat_size[0]*mat_size[5]*mat_size[4]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);
  
  //path 1 - ((A*(B*(C *D)))*E)
  inter1=reallocate(&inter1,mat_size[2],mat_size[4]);
  inter2=reallocate(&inter2,mat_size[1],mat_size[4]);
  inter3=reallocate(&inter3,mat_size[0],mat_size[4]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
 
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[2],mat_size[4],mat_size[3],1.0,C,mat_size[3],D,mat_size[4],0.0,inter1,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[4],mat_size[2],1.0,B,mat_size[2],inter1,mat_size[4],0.0,inter2,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[4],mat_size[1],1.0,A,mat_size[1],inter2,mat_size[4],0.0,inter3,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[4],1.0,inter3,mat_size[4],E,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[1][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[1] = 2.0 *( (mat_size[2]*mat_size[4]*mat_size[3]) + (mat_size[1]*mat_size[4]*mat_size[2]) );
  dummy_no_operation[1]=dummy_no_operation[1]+ 2.0* ( ( mat_size[0]*mat_size[4]*mat_size[1] )+ (mat_size[0]*mat_size[5]*mat_size[4]) );
  if(MATRIX!=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  //path 2 - (A*(B*(C*(D*E))))
  inter1=reallocate(&inter1,mat_size[3],mat_size[5]);
  inter2=reallocate(&inter2,mat_size[2],mat_size[5]);
  inter3=reallocate(&inter3,mat_size[1],mat_size[5]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[3],mat_size[5],mat_size[4],1.0,D,mat_size[4],E,mat_size[5],0.0,inter1,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[2],mat_size[5],mat_size[3],1.0,C,mat_size[3],inter1,mat_size[5],0.0,inter2,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[5],mat_size[2],1.0,B,mat_size[2],inter2,mat_size[5],0.0,inter3,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[1],1.0,A,mat_size[1],inter3,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[2][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[2] = 2.0 *( (mat_size[3]*mat_size[5]*mat_size[4]) + (mat_size[2]*mat_size[5]*mat_size[3]) );
  dummy_no_operation[2] = dummy_no_operation[2]+ 2.0 * ( ( mat_size[1]*mat_size[5]*mat_size[2] )+ (mat_size[0]*mat_size[5]*mat_size[1]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


  //path 3 - (A*(((B*C)*D)*E))
  inter1=reallocate(&inter1,mat_size[1],mat_size[3]);
  inter2=reallocate(&inter2,mat_size[1],mat_size[4]);
  inter3=reallocate(&inter3,mat_size[1],mat_size[5]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);
 

  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[3],mat_size[2],1.0,B,mat_size[2],C,mat_size[3],0.0,inter1,mat_size[3]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[4],mat_size[3],1.0,inter1,mat_size[3],D,mat_size[4],0.0,inter2,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[5],mat_size[4],1.0,inter2,mat_size[4],E,mat_size[5],0.0,inter3,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[1],1.0,A,mat_size[1],inter3,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[3][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[3]  = 2.0 * (mat_size[1]*mat_size[3]*mat_size[2]);
  dummy_no_operation[3] += 2.0 * (mat_size[1]*mat_size[4]*mat_size[3]);
  dummy_no_operation[3] += 2.0 * (mat_size[1]*mat_size[5]*mat_size[4]);
  dummy_no_operation[3] += 2.0 * (mat_size[0]*mat_size[5]*mat_size[1]);
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  
  //path 4 - (A*((B*C)*(D*E)))
  inter1=reallocate(&inter1,mat_size[1],mat_size[3]);
  inter2=reallocate(&inter2,mat_size[3],mat_size[5]);
  inter3=reallocate(&inter3,mat_size[1],mat_size[5]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[3],mat_size[2],1.0,B,mat_size[2],C,mat_size[3],0.0,inter1,mat_size[3]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[3],mat_size[5],mat_size[4],1.0,D,mat_size[4],E,mat_size[5],0.0,inter2,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[5],mat_size[3],1.0,inter1,mat_size[3],inter2,mat_size[5],0.0,inter3,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[1],1.0,A,mat_size[1],inter3,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[4][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[4] = 2.0 *( (mat_size[1]*mat_size[3]*mat_size[2]) + (mat_size[3]*mat_size[5]*mat_size[4]) );
  dummy_no_operation[4] += 2.0*( ( mat_size[1]*mat_size[5]*mat_size[3] )+ (mat_size[0]*mat_size[5]*mat_size[1]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  //path 5 - (A*(B*((C*D)*E)))  
  inter1=reallocate(&inter1,mat_size[2],mat_size[4]);
  inter2=reallocate(&inter2,mat_size[2],mat_size[5]);
  inter3=reallocate(&inter3,mat_size[1],mat_size[5]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[2],mat_size[4],mat_size[3],1.0,C,mat_size[3],D,mat_size[4],0.0,inter1,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[2],mat_size[5],mat_size[4],1.0,inter1,mat_size[4],E,mat_size[5],0.0,inter2,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[5],mat_size[2],1.0,B,mat_size[2],inter2,mat_size[5],0.0,inter3,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[1],1.0,A,mat_size[1],inter3,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[5][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[5] = 2.0 *( (mat_size[2]*mat_size[4]*mat_size[3]) + (mat_size[2]*mat_size[5]*mat_size[4]) );
  dummy_no_operation[5] += 2.0 *( ( mat_size[1]*mat_size[5]*mat_size[2] )+ (mat_size[0]*mat_size[5]*mat_size[1]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


  //path 6 - (A*((B*(C*D))*E))
  inter1=reallocate(&inter1,mat_size[2],mat_size[4]);
  inter2=reallocate(&inter2,mat_size[1],mat_size[4]);
  inter3=reallocate(&inter3,mat_size[1],mat_size[5]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[2],mat_size[4],mat_size[3],1.0,C,mat_size[3],D,mat_size[4],0.0,inter1,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[4],mat_size[2],1.0,B,mat_size[2],inter1,mat_size[4],0.0,inter2,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[5],mat_size[4],1.0,inter2,mat_size[4],E,mat_size[5],0.0,inter3,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[1],1.0,A,mat_size[1],inter3,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[6][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[6]  = 2.0 * (mat_size[2]*mat_size[4]*mat_size[3]);
  dummy_no_operation[6] += 2.0 * (mat_size[1]*mat_size[4]*mat_size[2]);
  dummy_no_operation[6] += 2.0 * (mat_size[1]*mat_size[5]*mat_size[4]);
  dummy_no_operation[6] += 2.0 * (mat_size[0]*mat_size[5]*mat_size[1]);
  if (MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  //path 7 - ((A*B)(C*(D*E)))
  inter1=reallocate(&inter1,mat_size[0],mat_size[2]);
  inter2=reallocate(&inter2,mat_size[3],mat_size[5]);
  inter3=reallocate(&inter3,mat_size[2],mat_size[5]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
 
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[2],mat_size[1],1.0,A,mat_size[1],B,mat_size[2],0.0,inter1,mat_size[2]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[3],mat_size[5],mat_size[4],1.0,D,mat_size[4],E,mat_size[5],0.0,inter2,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[2],mat_size[5],mat_size[3],1.0,C,mat_size[3],inter2,mat_size[5],0.0,inter3,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[2],1.0,inter1,mat_size[2],inter3,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[7][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[7] = 2.0 *( (mat_size[0]*mat_size[2]*mat_size[1]) + (mat_size[3]*mat_size[5]*mat_size[4]) );
  dummy_no_operation[7] += 2.0 *( ( mat_size[2]*mat_size[5]*mat_size[3] )+ (mat_size[0]*mat_size[5]*mat_size[2]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


  //path 8 - ((A*B)((C*D)*E))
  inter1=reallocate(&inter1,mat_size[0],mat_size[2]);
  inter2=reallocate(&inter2,mat_size[2],mat_size[4]);
  inter3=reallocate(&inter3,mat_size[2],mat_size[5]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[2],mat_size[1],1.0,A,mat_size[1],B,mat_size[2],0.0,inter1,mat_size[2]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[2],mat_size[4],mat_size[3],1.0,C,mat_size[3],D,mat_size[4],0.0,inter2,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[2],mat_size[5],mat_size[4],1.0,inter2,mat_size[4],E,mat_size[5],0.0,inter3,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[2],1.0,inter1,mat_size[2],inter3,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[8][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[8] = 2.0 *( (mat_size[0]*mat_size[2]*mat_size[1]) + (mat_size[2]*mat_size[4]*mat_size[3]) );
  dummy_no_operation[8] += 2.0* ( ( mat_size[2]*mat_size[5]*mat_size[4] )+ (mat_size[0]*mat_size[5]*mat_size[2]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  //path 9 - (((A*B)*C)(D*E))  
  inter1=reallocate(&inter1,mat_size[0],mat_size[2]);
  inter2=reallocate(&inter2,mat_size[3],mat_size[5]);
  inter3=reallocate(&inter3,mat_size[0],mat_size[3]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[2],mat_size[1],1.0,A,mat_size[1],B,mat_size[2],0.0,inter1,mat_size[2]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[3],mat_size[5],mat_size[4],1.0,D,mat_size[4],E,mat_size[5],0.0,inter2,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[3],mat_size[2],1.0,inter1,mat_size[2],C,mat_size[3],0.0,inter3,mat_size[3]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[3],1.0,inter3,mat_size[3],inter2,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[9][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[9] = 2.0 *( (mat_size[0]*mat_size[2]*mat_size[1]) + (mat_size[3]*mat_size[5]*mat_size[4]) );
  dummy_no_operation[9] += 2.0 * ( ( mat_size[0]*mat_size[3]*mat_size[2] )+ (mat_size[0]*mat_size[5]*mat_size[3]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  //path 10 - ((A*(B*C))*(D*E))
  inter1=reallocate(&inter1,mat_size[1],mat_size[3]);
  inter2=reallocate(&inter2,mat_size[3],mat_size[5]);
  inter3=reallocate(&inter3,mat_size[0],mat_size[3]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[3],mat_size[2],1.0,B,mat_size[2],C,mat_size[3],0.0,inter1,mat_size[3]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[3],mat_size[5],mat_size[4],1.0,D,mat_size[4],E,mat_size[5],0.0,inter2,mat_size[5]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[3],mat_size[1],1.0,A,mat_size[1],inter1,mat_size[3],0.0,inter3,mat_size[3]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[3],1.0,inter3,mat_size[3],inter2,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[10][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[10] = 2.0 *( (mat_size[1]*mat_size[3]*mat_size[2]) + (mat_size[3]*mat_size[5]*mat_size[4]) );
  dummy_no_operation[10] += 2.0 *( ( mat_size[0]*mat_size[3]*mat_size[1] )+ (mat_size[0]*mat_size[5]*mat_size[3]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  //path 11 - (((A*(B*C))*D)*E)
  inter1=reallocate(&inter1,mat_size[1],mat_size[3]);
  inter2=reallocate(&inter2,mat_size[0],mat_size[3]);
  inter3=reallocate(&inter3,mat_size[0],mat_size[4]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[3],mat_size[2],1.0,B,mat_size[2],C,mat_size[3],0.0,inter1,mat_size[3]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[3],mat_size[1],1.0,A,mat_size[1],inter1,mat_size[3],0.0,inter2,mat_size[3]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[4],mat_size[3],1.0,inter2,mat_size[3],D,mat_size[4],0.0,inter3,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[4],1.0,inter3,mat_size[4],E,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[11][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[11] = 2.0 *( (mat_size[1]*mat_size[3]*mat_size[2]) + (mat_size[0]*mat_size[3]*mat_size[1]) );
  dummy_no_operation[11] += 2.0 *( ( mat_size[0]*mat_size[4]*mat_size[3] )+ (mat_size[0]*mat_size[5]*mat_size[4]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  //path 12 - ((A*((B*C)*D))*E)
  inter1=reallocate(&inter1,mat_size[1],mat_size[3]);
  inter2=reallocate(&inter2,mat_size[1],mat_size[4]);
  inter3=reallocate(&inter3,mat_size[0],mat_size[4]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
 
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[3],mat_size[2],1.0,B,mat_size[2],C,mat_size[3],0.0,inter1,mat_size[3]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[1],mat_size[4],mat_size[3],1.0,inter1,mat_size[3],D,mat_size[4],0.0,inter2,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[4],mat_size[1],1.0,A,mat_size[1],inter2,mat_size[4],0.0,inter3,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[4],1.0,inter3,mat_size[4],E,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[12][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[12] = 2.0 *( (mat_size[1]*mat_size[3]*mat_size[2]) + (mat_size[1]*mat_size[4]*mat_size[3]) );
  dummy_no_operation[12] += 2.0 * ( ( mat_size[0]*mat_size[4]*mat_size[1] )+ (mat_size[0]*mat_size[5]*mat_size[4]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  //path 13 - (((A*B)*(C*D))*E)
  inter1=reallocate(&inter1,mat_size[0],mat_size[2]);
  inter2=reallocate(&inter2,mat_size[2],mat_size[4]);
  inter3=reallocate(&inter3,mat_size[0],mat_size[4]);
  inter4=reallocate(&inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
 
  gettimeofday(&start, NULL); 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[2],mat_size[1],1.0,A,mat_size[1],B,mat_size[2],0.0,inter1,mat_size[2]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[2],mat_size[4],mat_size[3],1.0,C,mat_size[3],D,mat_size[4],0.0,inter2,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[4],mat_size[2],1.0,inter1,mat_size[2],inter2,mat_size[4],0.0,inter3,mat_size[4]);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,mat_size[0],mat_size[5],mat_size[4],1.0,inter3,mat_size[4],E,mat_size[5],0.0,inter4,mat_size[5]);
  gettimeofday(&finish, NULL); 
  dummy_indv_duration[13][i] = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  }
  dummy_no_operation[13] = 2.0 *( (mat_size[0]*mat_size[2]*mat_size[1]) + (mat_size[2]*mat_size[4]*mat_size[3]) );
  dummy_no_operation[13] += 2.0 *( ( mat_size[0]*mat_size[4]*mat_size[2] )+ (mat_size[0]*mat_size[5]*mat_size[4]) );
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);                                                      //Actual matrix multiplication ends
  
  int time_index =0, flop_index =0;                                //stores index of min_time, min_flops among all paths
  double min_time=0,  path_minTime=0;
  double min_flops=0, path_minFlops=0;
  double deviation=0;


 for(i=0; i<14; i++) {
  inter_min_time =dummy_indv_duration[i][0];
  for (j=0; j<my_iteration; j++){
    if(dummy_indv_duration[i][j] <=inter_min_time){
      inter_min_time =dummy_indv_duration[i][j];
    }
    indv_duration[i] =inter_min_time;
  }
 }

 if (INTERMEDIATE_TIMINGS != 0)
 {
  for(i=0; i<14; i++){
    for(j=0; j<my_iteration; j++){  
    printf("%f \t",dummy_indv_duration[i][j]);
    }
  
    printf("%f %d \n",indv_duration [i],i);
  }

  for(i=0; i<14; i++)
  printf("%f \n",dummy_no_operation[i]);
 }



 for(i=0; i<14; i++)
 {                                                             //calculates total time, flops for each path          
   total_duration[i] = indv_duration[i];
   total_no_operation[i] = dummy_no_operation[i];
   total_gflops[i] = total_no_operation[i]/total_duration[i];
 }
  







  min_time =total_duration[0];                      //finds min_time                              
 for (i=1; i<14; i++){
   if(total_duration[i] <=min_time){
      min_time =total_duration[i];
      time_index =i;
   }
 }
 min_flops =total_no_operation[0];                  //finds min_flops
 for (i=1; i<14; i++){
  if(total_no_operation[i] <=min_flops){
    min_flops =total_no_operation[i];
    flop_index=i;
  }
 }

  //deviation = (min_time - timeof min_flop path)/min_time
  deviation = ( total_duration[flop_index] -total_duration[time_index] )/total_duration[time_index];


  for(i=0; i<14; i++)
  {
   total_duration_1[i] = total_duration[i];
   final_duration[i] = total_duration [i]; 
 }


  for(i=0; i<14; i++)
  {
     total_gflops_1[i] = total_no_operation[i];
     no_operation[i] = total_no_operation[i];
  }



 sorting(final_duration);
 sorting(no_operation);

  for (i=0; i<14; i++) 
  {
       for (j=0; j<14; j++) 
        {
         if (final_duration[i] == total_duration_1[j])
         {
            duration_index[i] = j;
         }
         else continue ;
        } 
  } 



  for (i=0; i<14; i++) 
  {
       for (j=0; j<14; j++) 
        {
          if (no_operation[i] == total_gflops_1[j])
          {
             gflops_index[i] = j;
          }
          else continue ;
        } 
  } 






if (OUTPUT_in_FILE !=0){
  FILE *fp;                                                           //writes output to a file 
  fp = fopen("result.txt", "a");
  for(i=0; i<14; i++)                                                                     
   fprintf(fp,"path[%d] \t\t%lf (s)\t\t%lf \t\t\t%lf TFLOPS\t\t\t\t %lf (s)\t[%d] \t\t%lf\t[%d]\n",i,total_duration[i],total_no_operation[i],total_gflops[i],final_duration[i],duration_index[i],no_operation[i],gflops_index[i]);
   //fprintf(fp,"path[%d] \t%lf s\t%lf \t%lf TFLOPS\n",i,total_duration[i],total_no_operation[i],total_gflops[i]);
  
  fprintf(fp,"\n");
  fprintf(fp,"\t\t  min_flops path[%d] takes %f s \n",flop_index,total_duration[flop_index]);
  fprintf(fp,"\t\t  min_time  path[%d] takes %f s \n",time_index,min_time );
  fprintf(fp,"\t\t  deviation is =%f \n",deviation);   
  fprintf(fp,"____******_______*******________________________________________*****************_____________________\n");
  fclose(fp);
}

if (OUTPUT_on_SCREEN !=0){
  for(i=0; i<14; i++)                                                                     //prints output on screen
   printf("path[%d] \t\t%lf (s)\t\t%lf \t\t\t%lf TFLOPS\t\t\t\t %lf (s)[%d] \t\t%lf\t[%d]\n",i,total_duration[i],total_no_operation[i],total_gflops[i],final_duration[i],duration_index[i],no_operation[i],gflops_index[i]);
   //printf("path[%d] \t%lf (s)\t%lf \t%lf TFLOPS\t\t\t\t %lf (s)[%d]\t %lf[%d]\n  ",i,total_duration[i],total_no_operation[i],total_gflops[i],final_duration[i],duration_index[i],no_operation[i],gflops_index[i]);
  printf("\n");
  
  printf("\t\t  min_flops path[%d] takes %f s \n",flop_index,total_duration[flop_index]);
  printf("\t\t  min_time  path[%d] takes %f s \n",time_index,min_time );
  printf("\t\t  deviation is =%f \n",deviation);
  printf("____******_______*******________________________________________*****************_____________________\n");
 }


 outer_time_index[outer_itr] = time_index;
 outer_flop_index[outer_itr] = flop_index ;
 outer_min_time[outer_itr]   = total_duration[time_index];
 outer__min_flops[outer_itr] = total_duration[flop_index];
 outer_deviation[outer_itr]  = deviation;
} //end of global iteration


printf("\n\n");
if (OUTPUT_on_SCREEN != 0)
{
for(i=0; i<global_iteration; i++)
printf( " iteration[%d] PATH[%d]-min_time=%f \t PATH[%d]-min_flops=%f \t deviation=%f \n",i,outer_time_index[i],outer_min_time[i],outer_flop_index[i],outer__min_flops[i],outer_deviation[i]);
}

if (OUTPUT_in_FILE !=0){
  FILE *fp;                                                           //writes output to a file 
  fp = fopen("result.txt", "a");
  for(i=0; i<global_iteration; i++)                                                                     
  fprintf( fp," iteration[%d] PATH[%d]-min_time=%f \t PATH[%d]-min_flops=%f \t deviation=%f \n",i,outer_time_index[i],outer_min_time[i],outer_flop_index[i],outer__min_flops[i],outer_deviation[i]);  
  fprintf(fp,"____******_______*******________________________________________*****************_____________________\n");
  fclose(fp);
}



  free(A);
  free(B);
  free(C);
  free(D);
  free(E);
  free(inter1);
  free(inter2);
  free(inter3);
  free(inter4);
  free(dummy_a);
  free(dummy_b);
  free(dummy_c);
  printf("Everything is successful\n");
  return 0;
}   //end of main

double* allocate(int size){
  double* locptr = (double* )malloc(sizeof(double) * size);
  if (locptr == NULL){
    printf("Memory is not alloacted \n");
  }
  return locptr;
}
double* reallocate(double** tempptr, int m, int n){
  *tempptr = (double *) realloc( *tempptr, m*n*sizeof(double) );
    if (*tempptr ==NULL)
        printf("reallocation failed \n");
   return *tempptr;
}



void initialize(double *locptr, int size){
  int i; 
  for (i=0; i<size; i++)
    locptr[i] = (rand() %2 );
  return ;
}
void initialize_zero(double *locptr,int size){
  int i; 
  for (i=0; i<size; i++)
    locptr[i] = 0.0;
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

void sorting (double* total)
{
  int i,j;
  double temp;
  for (i = 0; i < 14; i++)
    {
        for (j = 0; j < (14 - i - 1); j++)
        {
            if (total[j] > total[j + 1])

            {

                temp = total[j];

                total[j] = total[j + 1];

                total[j + 1] = temp;  
            }
        }
    }
}


