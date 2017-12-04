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


static int mat_size[7]; // stores the matrix sizes
 double indv_duration[15],indv_gflops[15],no_operation[15];


int main(int argc, char* argv[])
{
  int i,j;
  double total_duration[6],total_gflops[6],total_no_operation[6];
  double min_time=0,  path_minTime=0;
  double min_flops=0, path_minFlops=0;
  char *real_paths[] ={"((((A*B)*C)*D)*E)","((A*(B*(C *D)))*E)","(A*(B*(C*(D*E))))","(A*(((B*C)*D)*E))",
                 "(A*((B*C)*(D*E)))","(A*(B*((C*D)*E)))","(A*((B*(C*D))*E))","((A*B)(C*(D*E)))",
                 "((A*B)((C*D)*E))","(((A*B)*C)(D*E))","((A*(B*C))*(D*E))","(((A*(B*C))*D)*E)","((A*((B*C)*D))*E)","(((A*B)*(C*D))*E)"};
  char time_path[20],flop_path[20];
  initialize_zero(total_duration,6);
  initialize_zero(total_gflops,6);
  printf("test!\n");
  if(argc<6){
    printf("Input Error\n");
    return 1;
  }
 
   for (i=0; i<7; i++)
    mat_size[i]=0;

  mat_size[0] = atoi(argv[1]);//row of matrix A and row of result matrix =m
  mat_size[1] = atoi(argv[2]);// column of matrix A and row of matrix B = k
  mat_size[2] = atoi(argv[3]);// column of matrix B and row of matrix C =n
  mat_size[3] = atoi(argv[4]);// column of matrix C and row of matrix D
  mat_size[4] = atoi(argv[5]);// column of matrix D and row of matrix E
 
  
  double *A,*B,*C,*D,*E,*inter1=NULL,*inter2=NULL,*inter3=NULL,*dummy_a,*dummy_b,*dummy_c;
 
  //Allocate memeory for main matrices 
  A=allocate(mat_size[0]*mat_size[1]);
  B=allocate(mat_size[1]*mat_size[2]);
  C=allocate(mat_size[2]*mat_size[3]);
  D=allocate(mat_size[3]*mat_size[4]);
  dummy_a=allocate(2000*2000);
  dummy_b=allocate(2000*2000);
  dummy_c=allocate(2000*2000);
  
//initialize matrices with random values
  srand((unsigned)time(NULL));
  initialize(A,mat_size[0]*mat_size[1]);
  initialize(B,mat_size[1]*mat_size[2]);
  initialize(C,mat_size[2]*mat_size[3]);
  initialize(D,mat_size[3]*mat_size[4]);
 
  initialize(dummy_a,2000*2000);
  initialize(dummy_b,2000*2000);
  initialize(dummy_c,2000*2000);
 int  m =2000;
 int  n= 2000;
 int  k =2000;
 cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);  

  //Path 0 - (((A*B)*C)*D)
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],0);
  multi(inter1,C,inter2,&inter2,mat_size[0],mat_size[3],mat_size[2],1);
  multi(inter2,D,inter3,&inter3,mat_size[0],mat_size[4],mat_size[3],2);
printfunc(inter3,mat_size[0],mat_size[4]);
cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);

 //path 1 - ((A*(B*(C *D)))
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],3);
  multi(B,inter1,inter2,&inter2,mat_size[1],mat_size[4],mat_size[2],4);
  multi(A,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[1],5);
printfunc(inter3,mat_size[0],mat_size[4]);

cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  //path 2 - (A*((B*C)*D))
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],6);
  multi(inter1,D,inter2,&inter2,mat_size[1],mat_size[4],mat_size[3],7);
  multi(A,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[1],8);
printfunc(inter3,mat_size[0],mat_size[4]);


cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
 //path 3 - (A*((B*C)*D))
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],9);
  multi(A,inter1,inter2,&inter2,mat_size[0],mat_size[3],mat_size[1],10);
  multi(inter2,D,inter3,&inter3,mat_size[0],mat_size[4],mat_size[3],11);
printfunc(inter3,mat_size[0],mat_size[4]);


cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);
  //path 4 -((A*B)*(C*D))
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],12);
  multi(C,D,inter2,&inter2,mat_size[2],mat_size[4],mat_size[3],13);
  multi(inter1,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[2],14);
printfunc(inter3,mat_size[0],mat_size[4]);


 cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,dummy_a,k,dummy_b,n,0.0,dummy_c,n);


for(i=0; i<5; i++){
  total_duration[i] = indv_duration[3*i]+indv_duration[3*i+1]+indv_duration[3*i+2];
  total_no_operation[i] = no_operation[3*i]+ no_operation[3*i+1]+ no_operation[3*i+2];
  total_gflops[i] = total_no_operation[i]/total_duration[i];
}
  int time_index =0;
  int flop_index =0;
 min_time =total_duration[0];
  for (i=1; i<5; i++){
    if(total_duration[i] <=min_time){
      min_time =total_duration[i];
      time_index =i;
    }
  }
 min_flops =total_no_operation[0];
for (i=1; i<5; i++){
  if(total_no_operation[i] <=min_flops){
    min_flops =total_no_operation[i];
    flop_index=i;
  }
}
//flop_index =5;
//time_index =6;


double deviation=0;
deviation = ( total_duration[flop_index] -total_duration[time_index] )/total_duration[time_index];
strcpy(time_path,real_paths[time_index]);
strcpy(flop_path,real_paths[flop_index]);

FILE *fp;
  fp = fopen("result.txt", "a");
  for(i=0; i<5; i++)
 fprintf(fp,"path[%d] \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
  fprintf(fp,"\n");
 fprintf(fp," min_time =%f  \n",min_time );
  fprintf(fp,"min_flops =%f  and takes %f s \n",min_flops,total_duration[flop_index]);
  fprintf(fp,"deviation is =%f \n",deviation);
  fclose(fp);
 

for(i=0; i<5; i++)
 printf("path[%d] \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
  printf("\n");
 printf(" path[%d] min_time =%f  \n",time_index,min_time );
  printf("path[%d] min_flops =%f  and takes %f s \n",flop_index,min_flops,total_duration[flop_index]);
  printf("deviation is =%f \n",deviation);




  free(A);
  free(B);
  free(C);
  free(D);
  free(inter1);
  free(inter2);
  free(inter3);
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
    locptr[i] =1.0;
    return ;
}
void initialize_zero(double *locptr,int size){
    int i; 
    for (i=0; i<size; i++)
    locptr[i] = 0;
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
