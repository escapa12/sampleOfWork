#include <stdio.h>
#include "rtbps.h"
#include "cmani.h"

#define N 6
#define npasmx 1E4
#define pasmax 4
#define pasmin 1E-4
#define tolfl 1E-13

/* Per compilar:
        gcc -o cmani_rtbp -g -Wall cmani_rtbp.c cmani.c rtbps.c qrres.c  flux.c rk78.c -lm
        
   Per executar:
        ./cmani_rtbp 1.215058560962404e-2 1e-12=tolnwt 10=maxitnwt < cmani_rtbp.inp >cmani_rtbp.out
 */

int main(int argc, char *argv[]) {
  
  double mu,tolnwt,maxitnwt,x0[6],xf[6],dv[6],pas0=pasmin, dt; 
  
     if (argc!=4
         || sscanf(argv[1],"%lf",&mu)!=1
         || sscanf(argv[2],"%lf",&tolnwt)!=1
         || sscanf(argv[3],"%lf",&maxitnwt)!=1
      ) {
      fprintf(stderr, "./cmani_rtbp mu tolnwt maxitnwt\n");
      return -1;
    }
    while (scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  &dt,&x0[0],&x0[1],&x0[2],&x0[3],&x0[4],&x0[5],
	          &xf[0],&xf[1],&xf[2],&xf[3],&xf[4],&xf[5])==13) {   
     int i;
     for(i=0; i<N;i++)dv[i]=0;
     cmani(3, x0,xf,dt,dv,tolnwt, maxitnwt,pas0, pasmin,pasmax, tolfl,npasmx, rtbps ,&mu);
     for(i=0; i<N;i++) printf("%.16G ",dv[i]); printf("\n"); 
    }
  return 0;
}
