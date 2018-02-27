#include<stdio.h>
#include "rk78.h"
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
       }
       else break;
       i++;
  }
  while(i<npasmx){
     *h=tf-*t;
     if( rk78(t,x,h,pasmin,pasmax, tol, n, camp, prm)!=0) printf("ha fallat l'rk78\n"); 
     if(*t==tf) return 0;
     i++;
  }
  printf("flux falla: s'ha superat el numero maxim d'iteracions\n");
  return -1;
}