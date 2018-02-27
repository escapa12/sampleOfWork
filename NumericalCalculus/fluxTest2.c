#include <stdio.h>
#include "flux.h"
#include "camp.h"
/*
 * Per compilar:
        gcc -o fluxTest2 -g -Wall fluxTest2.c oscilharm.c flux.c rk78.c -lm
 */
#define N 2
#define PI 3.1415926535897932384626
  int main(){
  double pasmin, pasmax, npasmx,tol, t, x[N], h, T;
   pasmin= 1E-3;   pasmax= 4;   npasmx=1E4;   tol=1E-12; h=0.1;  
   t=0;  T=2*PI; 
   x[0]=0;   x[1]=1;
   //solucio del camp: y(t)=(sin(t),cos(t)) y(0)=(0,1)=y(2*PI)
   
   flux(&t, x, &h, T, pasmin, pasmax, tol, npasmx,N,oscilharm,NULL);
   printf("t0+T=%.16G x=%.16G,  %.16G\n",t,x[0],x[1]);
   printf("Tol=%.16G Error= %.16G,   %.16G\n",tol, x[0],1-x[1]); 
   return 0;
}