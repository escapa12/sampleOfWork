#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"qrres.h"

#define N 300
/*
 * Per compilar:
        gcc -o qr_3 -g -Wall qr_3.c qrres.c -lm
 */

int main(){
  int i,j,m=N,n=N;  
  double a[m*n-1],x[n-1],dr[n-1];
  double b[m-1];
  
  srand(time(0));
  for(i=0; i<m*n; i++)
     a[i]=((double)rand())/RAND_MAX;
  for(i=0; i<m; i++){ //elecio de b tal que x=1...1
     b[i]=0;
     for(j=0; j<n; j++){
        b[i]+= a[i+j*m];
     }
  }
   qrres(m,n,a,dr,b,x);
 double dif,max=0;
   for(i=0;i<n;i++){
     dif=fabs(x[i]-1);
    if(dif>max) max=dif; 
   }
  printf("error max: %.24G\n",max);
  
  return 0;
}