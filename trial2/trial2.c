#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"
#include "cblas.h"

double* initialize(int size);
double* allocate(double *locptr,int size);
double* allocate_zero(double *locptr,int size);

//extern void dgemm_(char*, char*, int*, int*,int*, double*, double*, int*, double*, int*, double*, double*, int*);

int main(int argc, char* argv[])
{
  int i,j;
  printf("test!\n");
  if(argc<4){
    printf("Input Error\n");
    return 1;
  }

  int m = atoi(argv[1]);
  int k = atoi(argv[2]);
  int n = atoi(argv[3]);
  int sizeofa = m * k;
  int sizeofb = k * n;
  int sizeofc = m * n;
  char ta = 'N';
  char tb = 'N';
  double alpha = 1.0;
  double beta = 0.0;

  struct timeval start,finish;
  double duration;
  double *A,*B,*C;
  
  A = initialize(sizeofa);
  B = initialize(sizeofb);
  C = initialize(sizeofc);
  srand((unsigned)time(NULL));
  A = allocate(A,sizeofa);
  B = allocate(B,sizeofb);
  C = allocate_zero(C,sizeofc);


  

  gettimeofday(&start, NULL);
  //dgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,alpha,A,k,B,n,beta,C,n);
  gettimeofday(&finish, NULL);

  duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  double gflops = 2.0 * m *n*k;
  gflops = gflops/duration*1.0e-6;

for (i=0; i<m; i++){
    for(j=0; j<k; j++){
      printf("%.2f \t",A[j+i*k] );
    }
    printf("\n");
  }
  printf("____________________________________________________________\n");
  for (i=0; i<k; i++){
    for(j=0; j<n; j++){
      printf("%.2f \t",B[j+i*n] );
    }
    printf("\n");
  }
  printf("____________________________________________________________\n");
  for (i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf("%.2f \t",C[j+i*n] );
    }
    printf("\n");
  }


  FILE *fp;
  fp = fopen("timeDGEMM.txt", "a");
  fprintf(fp, "%dx%dx%d\t%lf s\t%lf MFLOPS\n", m, n, k, duration, gflops);
  fclose(fp);

  free(A);
  free(B);
  free(C);
  return 0;
}