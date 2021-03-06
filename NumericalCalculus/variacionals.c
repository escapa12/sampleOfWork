#include <stdio.h>
#include "flux.h"
#include "camp.h"

 /*Per compilar:
        gcc -o variacionals -g -Wall variacionals.c campvar.c flux.c rk78.c -lm
   */
#define tol 1E-12
#define pasmax 10
#define pasmin 1E-4
#define npasmx 1E5
#define T 0.5
int main(){
    int N=6;
    double t=0, h=pasmin, x0[2], x[N], alpha=0.5;
    
    x0[0]=0.5;    x0[1]=0.5;
    x[0]=x0[0]; x[1]=x0[1]; x[2]=1; x[3]=0;  x[4]=0; x[5]=1;
    
    flux(&t,x,&h,T,pasmin,pasmax,tol,npasmx,N,campvar,&alpha);//  integraem el PVI i les variacionals per trobar Dflux
    printf("DxVar=%.16G %.16G %.16G %.16G\n",x[2],x[3],x[4],x[5]);

    N=2;
    int j;
    double y[N],Dfl[4], sgm=1E-6;
    for(j=0; j<N; j++){//aproximem Dflux amb la def de derivada
	y[0]=x0[0]; y[1]=x0[1]; t=0; h=pasmin; 	 
	y[j]+=sgm;
	flux(&t,y,&h,T,pasmin,pasmax,tol,npasmx,N,campvar,&alpha);
	Dfl[N*j]=y[0]; Dfl[1+N*j]=y[1];
	y[0]=x0[0]; y[1]=x0[1]; t=0; h=pasmin; 
	y[j]-=sgm;
	flux(&t,y,&h,T,pasmin,pasmax,tol,npasmx,N,campvar,&alpha);
	Dfl[N*j]= (Dfl[N*j]-y[0])/(2*sgm) ; Dfl[1+N*j]= (Dfl[1+N*j]-y[1])/(2*sgm) ;
    }
    
    printf("DFdif: %.16G %.16G %.16G %.16G\n ", Dfl[0], Dfl[1], Dfl[2], Dfl[3]);
    
    printf("Err:");
    for(j=0;j<4; j++)
         printf(" %.16G ", x[2+j]-Dfl[j]);
     printf("\n");
    return 0;
}