#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"
#include "cblas.h"
#define min(x,y) (((x) < (y)) ? (x) : (y))

double* allocate(int size);
void initialize_zero(double *locptr,int size);
void initialize( double* locptr, int size);
void multi(double* A,double*B, double*inter1, double **tempptr,int m, int n, int k );

int mat_size[7]; // stores the matrix sizes

int main(int argc, char* argv[])
{
  int i,j;
  printf("test!\n");
  if(argc<4){
    printf("Input Error\n");
    return 1;
  }
 
   for (i=0; i<7; i++)
    mat_size[i]=0;

  mat_size[0] = atoi(argv[1]);//row of first matrix and row of result matrix
  mat_size[1] = atoi(argv[2]);// column of first matrix and row of second matrix matrix
  mat_size[2] = atoi(argv[3]);// column of second matrix and column of result matrix
 
  struct timeval start,finish;
  double duration;
  double *A,*B,*inter1=NULL;
 
  //Allocate memeory for main matrices 
  A=allocate(mat_size[0]*mat_size[1]);
  B=allocate(mat_size[1]*mat_size[2]);
  
  //initialize matrices with random values
  srand((unsigned)time(NULL));
  initialize(A,mat_size[0]*mat_size[1]);
  initialize(B,mat_size[1]*mat_size[2]);
  

 multi(A,B,inter1,&inter1,mat_size[0],mat_size[2],mat_size[1]);

  
 /* gettimeofday(&start, NULL);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,k,B,n,0.0,C,n);
  gettimeofday(&finish, NULL);

  duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  double gflops = 2.0 * m *n*k;
  gflops = gflops/duration*1.0e-6;

  FILE *fp;
  fp = fopen("timeDGEMM.txt", "a");
  fprintf(fp, "%dx%dx%d\t%lf s\t%lf MFLOPS\n", m, k, n, duration, gflops);
  fclose(fp);
*/

  for (i=0; i<min(mat_size[0],6); i++) {
      for (j=0; j<min(mat_size[1],6); j++) {
        printf ("%12.0f", A[j+i*mat_size[1]]);
      }
      printf ("\n");
    }
  printf("____________________________________________________________\n");
   for (i=0; i<min(mat_size[1],6); i++) {
      for (j=0; j<min(mat_size[2],6); j++) {
        printf ("%12.0f", B[j+i*mat_size[2]]);
      }
      printf ("\n");
    }

  printf("____________________________________________________________\n");
 for (i=0; i<min(mat_size[0],6); i++) {
      for (j=0; j<min(mat_size[2],6); j++) {
        printf ("%12.5G", inter1[j+i*mat_size[2] ]);
      }
      printf ("\n");
    }




  free(A);
  free(B);
  free(inter1);
  return 0;
}



