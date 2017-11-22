#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"
#include "cblas.h"
#define min(x,y) (((x) < (y)) ? (x) : (y))

double* allocate(int size);
void initialize_zero(double *locptr,int size);
void initialize( double* locptr, int size);
void multi(double* A,double*B, double*inter1, double **tempptr,int m, int n, int k,int tt );
void printfunc (double* matrix, int row, int col);

static int mat_size[7]; // stores the matrix sizes
 double indv_duration[64],indv_gflops[64],no_operation[64];

int main(int argc, char* argv[])
{
  int i,j;
  double total_duration[14],total_gflops[14],total_no_operation[14];
  initialize_zero(total_duration,14);
  //initialize_zero(total_gflops,14);
  printf("test!\n");
  if(argc<7){
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
  mat_size[5] = atoi(argv[6]);// column of matrix E and column of result matrix 
  
  //struct timeval start,finish;
  
  double *A,*B,*C,*D,*E,*inter1=NULL,*inter2=NULL,*inter3=NULL,*inter4=NULL;
 
  //Allocate memeory for main matrices 
  A=allocate(mat_size[0]*mat_size[1]);
  B=allocate(mat_size[1]*mat_size[2]);
  C=allocate(mat_size[2]*mat_size[3]);
  D=allocate(mat_size[3]*mat_size[4]);
  E=allocate(mat_size[4]*mat_size[5]);
  //initialize matrices with random values
  srand((unsigned)time(NULL));
  initialize(A,mat_size[0]*mat_size[1]);
  initialize(B,mat_size[1]*mat_size[2]);
  initialize(C,mat_size[2]*mat_size[3]);
  initialize(D,mat_size[3]*mat_size[4]);
  initialize(E,mat_size[4]*mat_size[5]);
  
  //Tree 1 - ((((A*B)*C)*D)*E)
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],0);
  multi(inter1,C,inter2,&inter2,mat_size[0],mat_size[3],mat_size[2],1);
  multi(inter2,D,inter3,&inter3,mat_size[0],mat_size[4],mat_size[3],2);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],3); 
//printfunc(inter4,mat_size[0],mat_size[5]);
 //Tree 2 - ((A*(B*(C *D)))*E)
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],4);
  multi(B,inter1,inter2,&inter2,mat_size[1],mat_size[4],mat_size[2],5);
  multi(A,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[1],6);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],7);
//printfunc(inter4,mat_size[0],mat_size[5]);
  //Tree 3 - (A*(B*(C*(D*E))))
  multi(D,E,inter1,&inter1,mat_size[3],mat_size[5],mat_size[4],8);
  multi(C,inter1,inter2,&inter2,mat_size[2],mat_size[5],mat_size[3],9);
  multi(B,inter2,inter3,&inter3,mat_size[1],mat_size[5],mat_size[2],10);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],11);
//printfunc(inter4,mat_size[0],mat_size[5]);
 //Tree 4 - (A*(((B*C)*D)*E))
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],12);
  multi(inter1,D,inter2,&inter2,mat_size[1],mat_size[4],mat_size[3],13);
  multi(inter2,E,inter3,&inter3,mat_size[1],mat_size[5],mat_size[4],14);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],15); 
//printfunc(inter4,mat_size[0],mat_size[5]);
  //Tree 5 - (A*((B*C)*(D*E)))
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],16);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],17);
  multi(inter1,inter2,inter3,&inter3,mat_size[1],mat_size[5],mat_size[3],18);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],19); 
 //printfunc(inter4,mat_size[0],mat_size[5]);
  //Tree 6 - (A*(B*((C*D)*E)))
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],20);
  multi(inter1,E,inter2,&inter2,mat_size[2],mat_size[5],mat_size[4],21);
  multi(B,inter2,inter3,&inter3,mat_size[1],mat_size[5],mat_size[2],22);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],23); 
//printfunc(inter4,mat_size[0],mat_size[5]);
  //Tree 7 - (A*((B*(C*D))*E))
  multi(C,D,inter1,&inter1,mat_size[2],mat_size[4],mat_size[3],24);
  multi(B,inter1,inter2,&inter2,mat_size[1],mat_size[4],mat_size[2],25);
  multi(inter2,E,inter3,&inter3,mat_size[1],mat_size[5],mat_size[4],26);
  multi(A,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[1],27); 
 //printfunc(inter4,mat_size[0],mat_size[5]);
 //Tree 8 - ((A*B)(C*(D*E)))
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],28);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],29);
  multi(C,inter2,inter3,&inter3,mat_size[2],mat_size[5],mat_size[3],30);
  multi(inter1,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[2],31); 
   
//printfunc(inter4,mat_size[0],mat_size[5]);
 //Tree 9 - ((A*B)((C*D)*E))
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],32);
  multi(C,D,inter2,&inter2,mat_size[2],mat_size[4],mat_size[3],33);
  multi(inter2,E,inter3,&inter3,mat_size[2],mat_size[5],mat_size[4],34); 
  multi(inter1,inter3,inter4,&inter4,mat_size[0],mat_size[5],mat_size[2],35); 
  
 //printfunc(inter4,mat_size[0],mat_size[5]);

  //Tree 10 - (((A*B)*C)(D*E))
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],36);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],37);
  multi(inter1,C,inter3,&inter3,mat_size[0],mat_size[3],mat_size[2],38);
  multi(inter3,inter2,inter4,&inter4,mat_size[0],mat_size[5],mat_size[3],39); 
 //printfunc(inter4,mat_size[0],mat_size[5]);
 //Tree 11 - ((A*(B*C))*(D*E))
 multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],40);
  multi(D,E,inter2,&inter2,mat_size[3],mat_size[5],mat_size[4],41);
  multi(A,inter1,inter3,&inter3,mat_size[0],mat_size[3],mat_size[1],42);
  multi(inter3,inter2,inter4,&inter4,mat_size[0],mat_size[5],mat_size[3],43); 
 // printfunc(inter4,mat_size[0],mat_size[5]);
 
//Tree 12 - (((A*(B*C))*D)*E)
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],44);
  multi(A,inter1,inter2,&inter2,mat_size[0],mat_size[3],mat_size[1],45);
  multi(inter2,D,inter3,&inter3,mat_size[0],mat_size[4],mat_size[3],46);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],47); 
//printfunc(inter4,mat_size[0],mat_size[5]);
//Tree 13 - ((A*((B*C)*D))*E)
  multi(B,C,inter1,&inter1,mat_size[1],mat_size[3],mat_size[2],48);
  multi(inter1,D,inter2,&inter2,mat_size[1],mat_size[4],mat_size[3],49);
  multi(A,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[1],50);
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],51); 
 //printfunc(inter4,mat_size[0],mat_size[5]);
 //Tree 14 - (((A*B)*(C*D))*E)
  multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1],52);
  multi(C,D,inter2,&inter2,mat_size[2],mat_size[4],mat_size[3],53);
  multi(inter1,inter2,inter3,&inter3,mat_size[0],mat_size[4],mat_size[2],54); 
  multi(inter3,E,inter4,&inter4,mat_size[0],mat_size[5],mat_size[4],55); 
 
//printfunc(inter4,mat_size[0],mat_size[5]);
 
for(i=0; i<14; i++){
  total_duration[i] = indv_duration[4*i]+indv_duration[4*i+1]+indv_duration[4*i+2]+indv_duration[4*i+3];
  total_no_operation[i] = no_operation[4*i]+ no_operation[4*i+1]+ no_operation[4*i+2]+ no_operation[4*i+3];
  total_gflops[i] = total_no_operation[i]/total_duration[i];
//total_gflops[i] =  (no_operation[4*i]+ no_operation[4*i+1]+ no_operation[4*i+2]+ no_operation[4*i+3])/( indv_duration[4*i]+indv_duration[4*i+1]+indv_duration[4*i+2]+indv_duration[4*i+3]);
  
}
/*
  FILE *fp;
  fp = fopen("timeDGEMM.txt", "a");
  for(i=0; i<14; i++)
  fprintf(fp, "path[%]d \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
  fprintf("\n");
  fclose(fp);
 FILE *fp1;
  fp1 = fopen("timeDGEMM_f.txt", "a");
  for(i=0; i<56; i++)
  fprintf(fp1, "path[%d] \t%lf s \t%lf  \t%lf TFLOPS\n",i,indv_duration[i],no_operation[i],indv_gflops[i]);
  fprintf("\n");
  fclose(fp1);
 */

 for(i=0; i<14; i++)
  printf("path[%d] \t%lf s\t%lf TFLOPS \t%lf\n",i,total_duration[i],total_gflops[i],total_no_operation[i]);
  printf("\n");
 
/*  for(i=0; i<56; i++)
  printf("path[%d] \t%lf s \t%lf  \t%lf MFLOPS\n",i,indv_duration[i],no_operation[i],indv_gflops[i]);
  printf("\n");
*/
  free(A);
  free(B);
  free(C);
  free(D);
  free(E);
  free(inter1);
  free(inter2);
  free(inter3);
  free(inter4);
  printf("Everything is successful\n");
  return 0;
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
