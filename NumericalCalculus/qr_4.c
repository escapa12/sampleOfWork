#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"qrres.h"
#define N 3000
/*
 * Per compilar:
        gcc -o qr_4 -g -Wall qr_4.c qrres.c -lm
 */

int main(){
  int i,m=N,n=N;  
  double *a,*b,*x,dr[n-1];

  if((a=(double *)malloc((m*n-1)*sizeof(double)))==NULL) {printf("Error de memoria");}

  b=NULL;
  srand(time(0));
  for(i=0; i<m*n; i++)
     a[i]=((double)rand())/RAND_MAX;
 
  time_t t0, t1;  /* Mesurem des d’aquí ... */
  t0=clock();
  
  qrres(m,n,a,dr,b,x);

   t1=clock(); /* ... fins a aquí*/
   printf("n:%d temps: %.3f s\n", n,((double)(t1-t0))/CLOCKS_PER_SEC);
  return 0;
}


