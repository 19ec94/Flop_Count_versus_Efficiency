#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"
#include "cblas.h"
#include <string.h>
#define min(x,y) (((x) < (y)) ? (x) : (y))

double* allocate(int size);
void initialize_zero(double *locptr,int size);
void initialize( double* locptr, int size);
void multi(double* A,double*B, double*inter1, double **tempptr,int m, int n, int k,int tt );
void printfunc (double* matrix, int row, int col);

static int mat_size[5]; // stores the matrix sizes
double indv_duration[5],indv_gflops[5],no_operation[5];


int main(int argc, char* argv[])
{
  int i,j;
  double total_duration[3],total_gflops[3],total_no_operation[3];
  double min_time=0,  path_minTime=0;
  double min_flops=0, path_minFlops=0;
  //char *real_paths[3]={"((A*B)*C)","(A*(B*C))"};
  //char time_path[3],flop_path[3];
  initialize_zero(total_duration,14);
  initialize_zero(total_gflops,14);
  printf("test!\n");
  if(argc<5){
    printf("Input Error\n");
    return 1;
  }
 
   for (i=0; i<5; i++)
    mat_size[i]=0;

  mat_size[0] = atoi(argv[1]);//row of matrix A and row of result matrix =m
  mat_size[1] = atoi(argv[2]);// column of matrix A and row of matrix B = k
  mat_size[2] = atoi(argv[3]);// column of matrix B and row of matrix C =n
  mat_size[3] = atoi(argv[4]);
  //struct timeval start,finish;
  
  double *A,*B,*C,*inter1=NULL,*inter2=NULL;
 
  //Allocate memeory for main matrices 
  A=allocate(mat_size[0]*mat_size[1]);
  B=allocate(mat_size[1]*mat_size[2]);
  C=allocate(mat_size[2]*mat_size[3]);
    //initialize matrices with random values
  initialize(A,mat_size[0]*mat_size[1]);
  initialize(B,mat_size[1]*mat_size[2]);
  initialize(C,mat_size[2]*mat_size[3]);
   
  //path 0 - ((A*B)*C)
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],0);
  multi(inter1,C,inter2,&inter2,mat_size[0],mat_size[3],mat_size[2],1);
  printfunc(inter2,mat_size[0],mat_size[3]);
 //path 1 - (A*(B*C))
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],2);
  multi(A,inter1,inter2,&inter2,mat_size[0],mat_size[3],mat_size[1],3);
  printfunc(inter2,mat_size[0],mat_size[3]); 


 
for(i=0; i<2; i++){
  total_duration[i] = indv_duration[2*i]+indv_duration[2*i+1];
  total_no_operation[i] = no_operation[2*i]+ no_operation[2*i+1];
  total_gflops[i] = total_no_operation[i]/total_duration[i];
}
  int time_index =0;
  int flop_index =0;
 min_time =total_duration[0];
  for (i=1; i<2; i++){
    if(total_duration[i] <=min_time){
      min_time =total_duration[i];
      time_index =i;
    }
  }
 min_flops =total_no_operation[0];
for (i=1; i<2; i++){
  if(total_no_operation[i] <=min_flops){
    min_flops =total_no_operation[i];
    flop_index=i;
  }
}
printf("%d \n", flop_index);
printf("%d \n", time_index);

double deviation=0;
deviation = ( total_duration[flop_index] -total_duration[time_index] )/total_duration[time_index];
//strcpy(time_path,real_paths[time_index]);
//strcpy(flop_path,real_paths[flop_index]);

/*
FILE *fp;
  fp = fopen("result.txt", "a");
  for(i=0; i<2; i++)
  fprintf(fp,"path[%d] \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
  fprintf(fp,"\n");
  fprintf(fp,"min_time =%f  \n",min_time );
  fprintf(fp,"min_flops =%f  and takes %f s \n",min_flops,total_duration[flop_index]);
  fprintf(fp,"deviation is =%f \n",deviation);
  fclose(fp);
 */

for(i=0; i<2; i++)
 printf("path[%d] \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
 printf("\n");
 printf("path[%d] min_time =%f  \n",time_index,min_time );
 printf("path[%d] min_flops =%f  and takes %f s \n",flop_index,min_flops,total_duration[flop_index]);
 printf("deviation is =%f \n",deviation);




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
