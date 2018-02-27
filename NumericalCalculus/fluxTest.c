#include <stdio.h>
#include "rk78.h"
#include "camp.h"
/*
 * Per compilar:
        gcc -o fluxTest -g -Wall fluxTest.c camp.c rk78.c -lm
 */
#define N 1

int flux (double *t, double x[], double *h, double T,double pasmin, double pasmax, double tol,int npasmx,int n,int (*camp)(int n, double t, double x[], double f[], void *prm),void *prm);
 
int main(){
  double pasmin, pasmax, npasmx,tol, t, x[N], h, T;
   pasmin= 1E-3;   pasmax= 4;   npasmx=1E4;   tol=1E-12; t=1;  T=4;
   x[0]=1;   h=0.5;  
   
   flux(&t, x, &h, T, pasmin, pasmax, tol, npasmx,N,camp1,NULL);
   
   printf("t0+T=%.16G x=%.16G\n",t,x[0]);
   printf("Tol=%.16G Error= %.16G\n",tol, 25-x[0]); 
   return 0;
}

int flux (double *t, double x[], double *h, double T,double pasmin, double pasmax, double tol,int npasmx,int n,int (*camp)(int n, double t, double x[], double f[], void *prm),void *prm){
  int i=0;
  double tf=*t+T;  //temps final
  if( (*h)*T<=0) {
    printf("flux falla: h i T tenen signe diferents o alguna és 0\n");
    return -1;
  }
   while(i<npasmx){     
       if((*t+*h<tf && T>0) || (*t+*h>tf && T<0)){ 
	 if( rk78(t,x,h,pasmin,pasmax, tol, n, camp, prm)!=0) printf("ha fallat l'rk78\n"); 
	 printf("t:%.16G x:%.16G h:%.16G\n", *t, x[0], *h);
       }
       else break;
       i++;
  }
  while(i<npasmx){ printf("el pas previst fa superar el tf\n");
     *h=tf-*t;
     if( rk78(t,x,h,pasmin,pasmax, tol, n, camp, prm)!=0) printf("ha fallat l'rk78\n"); 
     printf("t:%.16G x:%.16G h:%.16G\n", *t, x[0], *h);
     if(*t==tf) return 0;
     i++;
  }
  printf("flux falla: s'ha superat el numero maxim d'iteracions\n");
  return -1;
}