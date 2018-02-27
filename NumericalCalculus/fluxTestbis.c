#include <stdio.h>
#include "flux.h"
#include "camp.h"
/*
 * Per compilar:
        gcc -o fluxTestbis -g -Wall fluxTestbis.c camp.c flux.c rk78.c  -lm
 */
#define N 1

  int main(){
  double pasmin, pasmax, npasmx,tol, t, t0,x[N], h, T;
   pasmin= 1E-10;   pasmax= 4;   npasmx=1E4;   tol=1E-12; 
   t0=1; t=t0; T=4;
   x[0]=1;   h=pasmin;  //cal que h tingui el mateix signe que T
   
   flux(&t, x, &h, T, pasmin, pasmax, tol, npasmx,N,camp1,NULL);   
   printf("t0+T=%.16G x=%.16G\n",t,x[0]);
   printf("Tol=%.16G Error= %.16G\n",tol, (t0+T)*(t0+T)-x[0]); 
   printf("\n");
   
   x[0]=1;   h=-pasmin; t=t0; T=-0.5;
   
   flux(&t, x, &h, T, pasmin, pasmax, tol, npasmx,N,camp1,NULL);   
   printf("t0+T=%.16G x=%.16G\n",t,x[0]);
   printf("Tol=%.16G Error= %.16G\n",tol, (t0+T)*(t0+T)-x[0]); 
   printf("\n");
   
   t0=-1;
   
   x[0]=1;   h=pasmin; t=t0; T=0.5;
   flux(&t, x, &h, T, pasmin, pasmax, tol, npasmx,N,camp1,NULL);   
   printf("t0+T=%.16G x=%.16G\n",t,x[0]);
   printf("Tol=%.16G Error= %.16G\n",tol, (t0+T)*(t0+T)-x[0]); 
   printf("\n");
   
   x[0]=1;   h=-pasmin; t=t0; T=-4;
   flux(&t, x, &h, T, pasmin, pasmax, tol, npasmx,N,camp1,NULL);   
   printf("t0+T=%.16G x=%.16G\n",t,x[0]);
   printf("Tol=%.16G Error= %.16G\n",tol, (t0+T)*(t0+T)-x[0]); 
   printf("\n");
   
   return 0;
}
