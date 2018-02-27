#include <stdio.h>
#include "pendol.h"
#include "cmani.h"

#define maxit 5
#define npasmx 1E4
#define pasmax 4
#define pasmin 1E-4
#define tolfl 1E-14
#define tol 1E-13

/*
 * Per compilar:
        gcc -o cmani_pendol -g -Wall cmani_pendol.c cmani.c pendol.c qrres.c  flux.c rk78.c -lm
 */

int main(){
  double x0[2],xf[2],dv[2],pas0=0.1, dt=1.57079633; 
  x0[0]=1 ;x0[1]=0 ;xf[0]=0 ; xf[1]= -0.95885108 ;  dv[0]=0;dv[1]=0;
  cmani(1, x0,xf,dt,dv,tol, maxit,pas0, pasmin,pasmax, tolfl,npasmx, pendol ,NULL);
 
  return 0;
}
