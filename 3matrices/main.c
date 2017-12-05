/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A C Program to calculate all the possible ways of multiplying                                %
% 3 matrices                                                                                   %
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
#define OUTPUT_in_FILE 0    //change it to 1, if output needs be written in a file
#define OUTPUT_on_SCREEN 1  //change it to Zero (0) , if you don't want to print it on screen
#define MATRIX 0             //change it to 1, if you want to check for the final matrix of each path


double* allocate(int size);                                                                   // Allocate memory for matrix
void initialize_zero(double *locptr,int size);                                                //Initialize matrix to zero
void initialize( double* locptr, int size);                                                   //Fill in matrices with random numbers
void multi(double* A,double*B, double*inter1, double **tempptr,int m, int n, int k,int tt );  //matrix multiplication
void printfunc (double* matrix, int row, int col);                                            //print resultant matrices

static int mat_size[5];                                                                       // stores the matrix sizes                                             
double indv_duration[5],indv_gflops[5],no_operation[5];                                       //stores values for each path


int main(int argc, char* argv[])
{
  int i,j;
  double total_duration[3],total_gflops[3],total_no_operation[3];
  double min_time=0,  path_minTime=0;
  double min_flops=0, path_minFlops=0;
  char *real_paths[]={"((A*B)*C)","(A*(B*C))"};
  char time_path[20],flop_path[20];
  initialize_zero(total_duration,3);
  initialize_zero(total_gflops,3);
  
  printf("test!\n");
  
  if(argc<5){
    printf("Input Error\n");
    return 1;
  }
 
   for (i=0; i<5; i++)
    mat_size[i]=0;

  mat_size[0] = atoi(argv[1]);                                                         //row of matrix A and row of result matrix =m
  mat_size[1] = atoi(argv[2]);                                                         // column of matrix A and row of matrix B = k
  mat_size[2] = atoi(argv[3]);                                                         // column of matrix B and row of matrix C =n
  mat_size[3] = atoi(argv[4]);                                                         // column of matrix C and Column of resultant matrix 
  
  
  double *A,*B,*C,*inter1=NULL,*inter2=NULL;                                           //Three main matrices two intermediate matrices
 
                                                                                       //Allocate memeory for main matrices 
  A=allocate(mat_size[0]*mat_size[1]);
  B=allocate(mat_size[1]*mat_size[2]);
  C=allocate(mat_size[2]*mat_size[3]);
  srand((unsigned)time(NULL));
                                                                                      //initialize main matrices with random values                                                                                    
  initialize(A,mat_size[0]*mat_size[1]);
  initialize(B,mat_size[1]*mat_size[2]);
  initialize(C,mat_size[2]*mat_size[3]);
   
  //path 0 - ((A*B)*C)                                                                //actual multiplication starts
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],0);
  multi(inter1,C,inter2,&inter2,mat_size[0],mat_size[3],mat_size[2],1);
  if(MATRIX !=0)
  printfunc(inter2,mat_size[0],mat_size[3]);
 
  //path 1 - (A*(B*C))
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],2);
  multi(A,inter1,inter2,&inter2,mat_size[0],mat_size[3],mat_size[1],3);
  if (MATRIX !=0)
  printfunc(inter2,mat_size[0],mat_size[3]);                                          //actual multiplication ends


 
  for(i=0; i<2; i++){                                                                 //calculation for total time,gflops taken by each path
    total_duration[i] = indv_duration[2*i]+indv_duration[2*i+1];
    total_no_operation[i] = no_operation[2*i]+ no_operation[2*i+1];
    total_gflops[i] = total_no_operation[i]/total_duration[i];
  }
  
  int time_index =0,flop_index =0;                                                    //stores index of path that takes min time, min flop operation
  
  min_time =total_duration[0];                                                        //finds minimum time path
  for (i=1; i<2; i++){
    if(total_duration[i] <=min_time){
      min_time =total_duration[i];
      time_index =i;
    }
  }
 
  min_flops =total_no_operation[0];                                                    //finds minimum flop path
 for (i=1; i<2; i++){
   if(total_no_operation[i] <=min_flops){
      min_flops =total_no_operation[i];
      flop_index=i;
    }
  }

  double deviation=0;                                                                   //calculates the deviation between min_time & time of min_flop
  deviation = ( total_duration[flop_index] -total_duration[time_index] )/total_duration[time_index];
  //strcpy(time_path,real_paths[time_index]);
  //strcpy(flop_path,real_paths[flop_index]);

if (OUTPUT_in_FILE !=0){
  FILE *fp;                                                                            //write output to a file "result.txt"
  fp = fopen("result.txt", "a");
  for(i=0; i<2; i++)
      fprintf(fp,"path[%d] \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
      
  fprintf(fp,"\n");
  fprintf(fp,"min_time =%f  \n",min_time );
  fprintf(fp,"min_flops =%f  and takes %f s \n",min_flops,total_duration[flop_index]);
  fprintf(fp,"deviation is =%f \n",deviation);
  fclose(fp);
}

if (OUTPUT_on_SCREEN !=0){
for(i=0; i<2; i++)                                                                     //print on screen
    printf("path[%d] \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
    
  printf("\n");
  printf("path[%d] min_time =%f  \n",time_index,min_time );
  printf("path[%d] min_flops =%f  and takes %f s \n",flop_index,min_flops,total_duration[flop_index]);
  printf("deviation is =%f \n",deviation);
}



  free(A);
  free(B);
  free(C);
  free(inter1);
  free(inter2);
  printf("Everything is successful\n");
  
return 0;
}

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
    locptr[i] =1.0;
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
