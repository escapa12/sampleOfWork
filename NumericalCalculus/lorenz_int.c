#include <stdio.h>
#include "flux.h"
#include "camp.h"

/*
 * Per compilar:
        gcc -o lorenz_int -g -Wall lorenz_int.c camp_lorenz.c flux.c rk78.c -lm
   Per executar:
        echo "-0.4164607449115608 -0.9089362634520914 0.01438311162938695" | ./lorenz_int >lorenz.txt 3 26.5 1.0 150 10000
      
   Gnuplot:
     set term x11 size 700,700
     set view equal xyz
     set ticslevel 0.1
     unset key
     splot 'lorenz.txt' u 2:3:4 w l
     set output 'lorenz.eps'
     set term post eps color solid
     replot
     set term x11
     set output
 */

#define N 3
int main (int argc, char *argv[]) {
 double sgm, rho, bet, tf, t,h,x[3],prm[3];
 int nt;
 if(argc!=6
         || sscanf(argv[1],"%lf",&sgm)!=1
         || sscanf(argv[2],"%lf",&rho)!=1
         || sscanf(argv[3],"%lf",&bet)!=1
         || sscanf(argv[4],"%lf",&tf)!=1
         || sscanf(argv[5],"%d",&nt)!=1
      ) {
      fprintf(stderr, "./lorenz_int sgm rho bet tf nt\n");
      return -1;
   }
  prm[0]=sgm;prm[1]=rho;prm[2]=bet;
 
   while (scanf("%lf %lf %lf",&x[0],&x[1],&x[2])==3) { 
     int i;
     t=0; h=1E-4;
     printf("0 %.16G %.16G %.16G\n", x[0], x[1], x[2]); //primer punt 
     for(i=1; i<=nt ; i++){
       flux(&t, x, &h,tf/nt, 1E-6, 10, 1E-12, 1E4,N,lorenz,prm);
       printf("%.16G %.16G %.16G %.16G\n",t, x[0], x[1], x[2]);
     }
     printf("\n");
   }
   return 0;
}