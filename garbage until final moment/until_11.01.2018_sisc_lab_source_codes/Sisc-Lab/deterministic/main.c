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
#define OUTPUT_in_FILE 0    //change it to 1, if output needs be written in a file
#define OUTPUT_on_SCREEN 1  //change it to Zero (0) , if you don't want to print it on screen
#define MATRIX 0            //change it to 1, if you want to check for the final matrix of each path

double* allocate(int size);                                                                                           //alloacte memory for matrices
void initialize_zero(double *locptr,int size);                                                                        //initialize matrices to zero
void initialize( double* locptr, int size);                  //Initialize matrices to random values                              
void multi(double* A,double*B, double*inter, double **tempptr,int m, int n, int k,int tt,int itr );                          //matrix multiplication  
void printfunc (double* matrix, int row, int col);                                                                    //print resultant matrix 


static int mat_size[7];                                                                                                // stores the matrix sizes
double dummy_indv_duration[64][30],indv_gflops[64],dummy_no_operation[64][30],indv_duration[64],no_operation[64];


int main(int argc, char* argv[])
{
  int i,j;
  double total_duration[14],total_gflops[14],total_no_operation[14];                                                   //stores values for each paths
  char *all_paths[] ={"((((A*B)*C)*D)*E)","((A*(B*(C *D)))*E)","(A*(B*(C*(D*E))))","(A*(((B*C)*D)*E))",
                 "(A*((B*C)*(D*E)))","(A*(B*((C*D)*E)))","(A*((B*(C*D))*E))","((A*B)(C*(D*E)))",
                 "((A*B)((C*D)*E))","(((A*B)*C)(D*E))","((A*(B*C))*(D*E))","(((A*(B*C))*D)*E)","((A*((B*C)*D))*E)","(((A*B)*(C*D))*E)"};
  char time_path[20],flop_path[20];
  initialize_zero(total_duration,14);
  initialize_zero(total_gflops,14);
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

int my_iteration = 5;
double inter_min_time ;

  for(i=0; i<my_iteration;  i++){
  //Path 0 - ((((A*B)*C)*D)*E)             
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n); 
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],0,i);
  multi(inter1,C,inter2,&inter2,mat_size[0],mat_size[3],mat_size[2],1,i);
  multi(inter2,D,inter3,&inter3,mat_size[0],mat_size[4],mat_size[3],2,i);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],3,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);
  
  for(i=0; i<my_iteration;  i++){
  //path 1 - ((A*(B*(C *D)))*E)
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],4,i);
  multi(B,inter1,inter2,&inter2,mat_size[1],mat_size[4],mat_size[2],5,i);
  multi(A,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[1],6,i);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],7,i);
  }
  if(MATRIX!=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  //path 2 - (A*(B*(C*(D*E))))
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(D,E,inter1,&inter1,mat_size[3],mat_size[5],mat_size[4],8,i);
  multi(C,inter1,inter2,&inter2,mat_size[2],mat_size[5],mat_size[3],9,i);
  multi(B,inter2,inter3,&inter3,mat_size[1],mat_size[5],mat_size[2],10,i);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],11,i);
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  //path 3 - (A*(((B*C)*D)*E))
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],12,i);
  multi(inter1,D,inter2,&inter2,mat_size[1],mat_size[4],mat_size[3],13,i);
  multi(inter2,E,inter3,&inter3,mat_size[1],mat_size[5],mat_size[4],14,i);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],15,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  //path 4 - (A*((B*C)*(D*E)))
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],16,i);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],17,i);
  multi(inter1,inter2,inter3,&inter3,mat_size[1],mat_size[5],mat_size[3],18,i);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],19,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  //path 5 - (A*(B*((C*D)*E)))  
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],20,i);  
  multi(inter1,E,inter2,&inter2,mat_size[2],mat_size[5],mat_size[4],21,i);
  multi(B,inter2,inter3,&inter3,mat_size[1],mat_size[5],mat_size[2],22,i);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],23,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  //path 6 - (A*((B*(C*D))*E))
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],24,i);
  multi(B,inter1,inter2,&inter2,mat_size[1],mat_size[4],mat_size[2],25,i);
  multi(inter2,E,inter3,&inter3,mat_size[1],mat_size[5],mat_size[4],26,i);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],27,i); 
  }
  if (MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  //path 7 - ((A*B)(C*(D*E)))
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],28,i);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],29,i);
  multi(C,inter2,inter3,&inter3,mat_size[2],mat_size[5],mat_size[3],30,i);
  multi(inter1,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[2],31,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  //path 8 - ((A*B)((C*D)*E))
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],32,i);
  multi(C,D,inter2,&inter2,mat_size[2],mat_size[4],mat_size[3],33,i);
  multi(inter2,E,inter3,&inter3,mat_size[2],mat_size[5],mat_size[4],34,i); 
  multi(inter1,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[2],35,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  //path 9 - (((A*B)*C)(D*E))
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],36,i);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],37,i);
  multi(inter1,C,inter3,&inter3,mat_size[0],mat_size[3],mat_size[2],38,i);
  multi(inter3,inter2,inter4,&inter4,mat_size[0],mat_size[5],mat_size[3],39,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

  for(i=0; i<my_iteration;  i++){
  //path 10 - ((A*(B*C))*(D*E))
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],40,i);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],41,i);
  multi(A,inter1,inter3,&inter3,mat_size[0],mat_size[3],mat_size[1],42,i);
  multi(inter3,inter2,inter4,&inter4,mat_size[0],mat_size[5],mat_size[3],43,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);
 
  for(i=0; i<my_iteration;  i++){
  //path 11 - (((A*(B*C))*D)*E)
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],44,i);
  multi(A,inter1,inter2,&inter2,mat_size[0],mat_size[3],mat_size[1],45,i);
  multi(inter2,D,inter3,&inter3,mat_size[0],mat_size[4],mat_size[3],46,i);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],47,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  //path 12 - ((A*((B*C)*D))*E)
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],48,i);
  multi(inter1,D,inter2,&inter2,mat_size[1],mat_size[4],mat_size[3],49,i);
  multi(A,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[1],50,i);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],51,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


  for(i=0; i<my_iteration;  i++){
  //path 13 - (((A*B)*(C*D))*E)
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],52,i);
  multi(C,D,inter2,&inter2,mat_size[2],mat_size[4],mat_size[3],53,i);
  multi(inter1,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[2],54,i); 
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],55,i); 
  }
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);                                                      //Actual matrix multiplication ends
  
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);


  for(i=0; i<56; i++) {
  inter_min_time =dummy_indv_duration[i][0];
  for (j=0; j<my_iteration; j++){
    if(dummy_indv_duration[i][j] <=inter_min_time){
      inter_min_time =dummy_indv_duration[i][j];
    }
    indv_duration[i] =inter_min_time;
  }
 }


 for(i=0; i<56; i++){
  for(j=0; j<my_iteration; j++){  
   printf("%f \t",dummy_indv_duration[i][j]);
  }
  printf("%f Multi(%d) \n",indv_duration [i],i);
 }

 /*
 for(i=0; i<56; i++){
  for(j=0; j<my_iteration; j++){  
   printf("%f \t",dummy_no_operation[i][j]);
  }
  printf("%d \n",i);
 }

 */


  for(i=0; i<14; i++){                                                             //calculates total time, flops for each path          
   total_duration[i] = indv_duration[4*i]+indv_duration[4*i+1]+indv_duration[4*i+2]+indv_duration[4*i+3];
   total_no_operation[i] = dummy_no_operation[4*i][0]+ dummy_no_operation[4*i+1][0]+ dummy_no_operation[4*i+2][0]+ dummy_no_operation[4*i+3][0];
   total_gflops[i] = total_no_operation[i]/total_duration[i];
  }
  int time_index =0, flop_index =0;                                //stores index of min_time, min_flops among all paths
  double min_time=0,  path_minTime=0;
  double min_flops=0, path_minFlops=0;
  
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


  double deviation=0;                            //deviation = (min_time - timeof min_flop path)/min_time
  deviation = ( total_duration[flop_index] -total_duration[time_index] )/total_duration[time_index];
  strcpy(time_path,all_paths[time_index]);
  strcpy(flop_path,all_paths[flop_index]);

 
if (OUTPUT_in_FILE !=0){
  FILE *fp;                                                           //writes output to a file 
  fp = fopen("result.txt", "a");
  for(i=0; i<14; i++)                                                                    
   fprintf(fp,"path[%d] \t%lf s\t%lf \t%lf TFLOPS\n",i,total_duration[i],total_no_operation[i],total_gflops[i]);
  fprintf(fp,"\n");
  fprintf(fp,"\t\t  min_flops path[%d] takes %f s \n",flop_index,total_duration[flop_index]);
  fprintf(fp,"\t\t  min_time  path[%d] takes %f s \n",time_index,min_time );
  fprintf(fp,"\t\t  deviation is =%f \n",deviation);   
  fclose(fp);
}

 if (OUTPUT_on_SCREEN !=0){
  for(i=0; i<14; i++)                                                                     //prints output on screen
   printf("path[%d] \t%lf s\t%lf \t%lf TFLOPS\n",i,total_duration[i],total_no_operation[i],total_gflops[i]);
  printf("\n");
  
  printf("\t\t  min_flops path[%d] takes %f s \n",flop_index,total_duration[flop_index]);
  printf("\t\t  min_time  path[%d] takes %f s \n",time_index,min_time );
  printf("\t\t  deviation is =%f \n",deviation);
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
    locptr[i] = 1.0;
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
