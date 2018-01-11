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
#define MATRIX 1            //change it to 1, if you want to check for the final matrix of each path

double* allocate(int size);                                                                                           //alloacte memory for matrices
void initialize_zero(double *locptr,int size);                                                                        //initialize matrices to zero
void initialize( double* locptr, int size);                                                  
void multi(double* A,double*B, double*inter, double **tempptr,int m, int n, int k,int tt,int dum_itr );                          
void printfunc (double* matrix, int row, int col);                                                                    //print resultant matrix 


static int mat_size[7];                                                                                                // stores the matrix sizes
double dummy_indv_duration[64][8],indv_gflops[64],no_operation[64],indv_duration[64];            
int main(int argc, char* argv[])
{
  int i,j;
  double inter_min_time;
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
  int  m =4000,n=4000,k=4000;

                                                         
  double *A,*B,*C,*D,*E,*inter1=NULL,*inter2=NULL,*inter3=NULL,*inter4=NULL,*dummy_a,*dummy_b,*dummy_c; 
                                                 
  A=allocate(mat_size[0]*mat_size[1]);
  B=allocate(mat_size[1]*mat_size[2]);
  C=allocate(mat_size[2]*mat_size[3]);
  D=allocate(mat_size[3]*mat_size[4]);
  E=allocate(mat_size[4]*mat_size[5]);
  dummy_a=allocate(m*n);
  dummy_b=allocate(n*k);
  dummy_c=allocate(k*n);
  
  srand((unsigned)time(NULL));
  initialize(A,mat_size[0]*mat_size[1]);
  initialize(B,mat_size[1]*mat_size[2]);
  initialize(C,mat_size[2]*mat_size[3]);
  initialize(D,mat_size[3]*mat_size[4]);
  initialize(E,mat_size[4]*mat_size[5]);
  initialize(dummy_a,m*n);
  initialize(dummy_b,n*k);
  initialize(dummy_c,k*m);

      int my_iteration =5;
      int kiwi;

  //Path 0 - ((((A*B)*C)*D)*E)                                                                            
for(kiwi=0; kiwi<my_iteration; kiwi++) {
cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],0,kiwi);
  multi(inter1,C,inter2,&inter2,mat_size[0],mat_size[3],mat_size[2],1,kiwi);
  multi(inter2,D,inter3,&inter3,mat_size[0],mat_size[4],mat_size[3],2,kiwi);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],3,kiwi); 
 }
if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);
  
 

//path 1 - ((A*(B*(C *D)))*E)
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],4,kiwi);
  multi(B,inter1,inter2,&inter2,mat_size[1],mat_size[4],mat_size[2],5,kiwi);
  multi(A,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[1],6,kiwi);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],7,kiwi);
  }
if(MATRIX!=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

//path 2 - (A*(B*(C*(D*E))))
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(D,E,inter1,&inter1,mat_size[3],mat_size[5],mat_size[4],8,kiwi);
  multi(C,inter1,inter2,&inter2,mat_size[2],mat_size[5],mat_size[3],9,kiwi);
  multi(B,inter2,inter3,&inter3,mat_size[1],mat_size[5],mat_size[2],10,kiwi);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],11,kiwi);
}
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


 //path 3 - (A*(((B*C)*D)*E))
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],12,kiwi);
  multi(inter1,D,inter2,&inter2,mat_size[1],mat_size[4],mat_size[3],13,kiwi);
  multi(inter2,E,inter3,&inter3,mat_size[1],mat_size[5],mat_size[4],14,kiwi);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],15,kiwi); 
} 
 if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


 //path 4 - (A*((B*C)*(D*E)))
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n); 
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],16,kiwi);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],17,kiwi);
  multi(inter1,inter2,inter3,&inter3,mat_size[1],mat_size[5],mat_size[3],18,kiwi);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],19,kiwi); 
} 
 if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

 //path 5 - (A*(B*((C*D)*E)))
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],20,kiwi);  
  multi(inter1,E,inter2,&inter2,mat_size[2],mat_size[5],mat_size[4],21,kiwi);
  multi(B,inter2,inter3,&inter3,mat_size[1],mat_size[5],mat_size[2],22,kiwi);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],23,kiwi); 
}  
if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);



//path 6 - (A*((B*(C*D))*E))
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n); 
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],24,kiwi);
  multi(B,inter1,inter2,&inter2,mat_size[1],mat_size[4],mat_size[2],25,kiwi);
  multi(inter2,E,inter3,&inter3,mat_size[1],mat_size[5],mat_size[4],26,kiwi);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],27,kiwi); 
} 
 if (MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


 //path 7 - ((A*B)(C*(D*E)))
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],28,kiwi);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],29,kiwi);
  multi(C,inter2,inter3,&inter3,mat_size[2],mat_size[5],mat_size[3],30,kiwi);
  multi(inter1,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[2],31,kiwi); 
} 
 if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);

 //path 8 - ((A*B)((C*D)*E))
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],32,kiwi);
  multi(C,D,inter2,&inter2,mat_size[2],mat_size[4],mat_size[3],33,kiwi);
  multi(inter2,E,inter3,&inter3,mat_size[2],mat_size[5],mat_size[4],34,kiwi); 
  multi(inter1,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[2],35,kiwi); 
}
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


//path 9 - (((A*B)*C)(D*E))
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],36,kiwi);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],37,kiwi);
  multi(inter1,C,inter3,&inter3,mat_size[0],mat_size[3],mat_size[2],38,kiwi);
  multi(inter3,inter2,inter4,&inter4,mat_size[0],mat_size[5],mat_size[3],39,kiwi); 
} 
 if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);


 //path 10 - ((A*(B*C))*(D*E))
for(kiwi=0; kiwi<my_iteration; kiwi++) { 
 cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],40,kiwi);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],41,kiwi);
  multi(A,inter1,inter3,&inter3,mat_size[0],mat_size[3],mat_size[1],42,kiwi);
  multi(inter3,inter2,inter4,&inter4,mat_size[0],mat_size[5],mat_size[3],43,kiwi); 
}
  if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);
 

 //path 11 - (((A*(B*C))*D)*E)
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],44,kiwi);
  multi(A,inter1,inter2,&inter2,mat_size[0],mat_size[3],mat_size[1],45,kiwi);
  multi(inter2,D,inter3,&inter3,mat_size[0],mat_size[4],mat_size[3],46,kiwi);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],47,kiwi); 
} 
 if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);



//path 12 - ((A*((B*C)*D))*E)
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],48,kiwi);
  multi(inter1,D,inter2,&inter2,mat_size[1],mat_size[4],mat_size[3],49,kiwi);
  multi(A,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[1],50,kiwi);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],51,kiwi); 
} 
 if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);



  //path 13 - (((A*B)*(C*D))*E)
for(kiwi=0; kiwi<my_iteration; kiwi++) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],52,kiwi);
  multi(C,D,inter2,&inter2,mat_size[2],mat_size[4],mat_size[3],53,kiwi);
  multi(inter1,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[2],54,kiwi); 
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],55,kiwi); 
} 
 if(MATRIX !=0)
  printfunc(inter4,mat_size[0],mat_size[5]);   

printf("___________________________________\n");
                                                                                                                            
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
printf("%f \t",dummy_indv_duration [i][j]);
}
printf("%f %d \n",indv_duration [i],i);
}



  for(i=0; i<14; i++){                                                                                                           
   total_duration[i] = indv_duration[4*i]+indv_duration[4*i+1]+indv_duration[4*i+2]+indv_duration[4*i+3];
   total_no_operation[i] = no_operation[4*i]+ no_operation[4*i+1]+ no_operation[4*i+2]+ no_operation[4*i+3];
   total_gflops[i] = total_no_operation[i]/total_duration[i];
  }
  int time_index =0, flop_index =0;                                                        
  double min_time=0,  path_minTime=0;
  double min_flops=0, path_minFlops=0;
  
  min_time =total_duration[0];                     
  for (i=1; i<14; i++){
    if(total_duration[i] <=min_time){
      min_time =total_duration[i];
      time_index =i;
    }
  }
 min_flops =total_no_operation[0];                                                                                               
  for (i=1; i<14; i++){
   if(total_no_operation[i] <=min_flops){
    min_flops =total_no_operation[i];
    flop_index=i;
   }
  }

  double deviation=0;                    
  deviation = ( total_duration[flop_index] -total_duration[time_index] )/total_duration[time_index];
  //strcpy(time_path,all_paths[time_index]);
  //strcpy(flop_path,all_paths[flop_index]);

if (OUTPUT_in_FILE !=0){
  FILE *fp;          
  fp = fopen("result.txt", "a");
   for(i=0; i<14; i++)
    fprintf(fp,"path[%d] \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
  fprintf(fp,"\n");
  fprintf(fp," min_time =%f  \n",min_time );
  fprintf(fp,"min_flops =%f  and takes %f s \n",min_flops,total_duration[flop_index]);
  fprintf(fp,"deviation is =%f \n",deviation);
  fclose(fp);
}

if (OUTPUT_on_SCREEN !=0){
  for(i=0; i<14; i++)           
   printf("path[%d] \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
  printf("\n");
  printf(" path[%d] min_time =%f  \n",time_index,min_time );
  printf("path[%d] min_flops =%f  and takes %f s \n",flop_index,min_flops,total_duration[flop_index]);
  printf("deviation is =%f \n",deviation);
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
