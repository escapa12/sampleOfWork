#include <stdio.h>
#include "rk78.h"
#include "pendol.h"

/*
 * Per compilar:
        gcc -o rf_pendol -g -Wall rf_pendol.c pendol.c rk78.c -lm
   Per executar(exercici 2)
       ./rf_pendol 0.01 3 1E-12 <rf_inic.txt >rfasependol.txt
   Comandes gnuplot:
       set xrange [-10:10]; set yrange [-3:3]
       plot'rfasependol.txt' u 2:3 w l
       set output ’lorenz.eps’
       set term post eps color solid
       set out 'rf_pendol.eps'
       replot
       set out
 */

#define N 2

int main (int argc, char *argv[]) {
   double h0, hmin, hmax, tol, t, x[N], h;
   int np, i;
   if (argc!=4
         || sscanf(argv[1],"%lf",&hmin)!=1
         || sscanf(argv[2],"%lf",&hmax)!=1
         || sscanf(argv[3],"%lf",&tol)!=1
      ) {
      fprintf(stderr, "./rf_pendol hmin hmax tol\n");
      return -1;
   }
   while (scanf("%lf %lf %lf %d",&x[0],&x[1],&h0,&np)==4) {
      t=0; h=h0;
      printf("%.16G %.16G %.16G %G\n", t, x[0], x[1], h);
      for (i=0; i<np; i++) {
         rk78(&t,x,&h,hmin,hmax,tol,N,pendol,NULL);
         printf("%.16G %.16G %.16G %G\n", t, x[0], x[1], h);
      }
      printf("\n");
   }
   return 0;
}
