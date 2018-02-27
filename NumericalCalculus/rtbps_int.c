#include <stdio.h>
#include "flux.h"
#include "rtbps.h"
/*
Per compilar:
        gcc -o rtbps_int -g -Wall rtbps_int.c rtbps.c flux.c rk78.c -lm
   Per executar:
        ./rtbps_int 1.215058560962404e-2 < halos.inp > halos.txt
   Gnuplot:
       set term x11 size 1400,1400
       set view equal xyz
       set ticslevel 0.1
       set ytics .1
       set xlabel 'x' ; set ylabel 'y' ; set zlabel 'z'
       splot 'halos.txt' u 2:3:4 w l,'L1.txt' w p,'trr.txt' w p,'lln.txt' w p
       set output 'halos.eps'
       set term post eps color solid
       replot
       set term x11
       set output
       
   */
#define N 6

int main (int argc, char *argv[]) {

     double mu,t,h,x[N],tf;
     int nt;
     if(argc!=2  || sscanf(argv[1],"%lf",&mu)!=1 ) {
        fprintf(stderr, "./rtbsps_int mu\n");
        return -1;
     }  
     while (scanf("%lf %lf %lf %lf %lf %lf %lf %d ",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&tf ,&nt)==8 ){ 
     int i;
     t=0; h=1E-6;
     printf("0 %.16G %.16G %.16G\n", x[0], x[1], x[2]); //primer punt 
     for(i=1; i<=nt ; i++){
         flux(&t, x, &h, tf/nt, 1E-6, 10, 1E-12, 1E4,N,rtbps,&mu);
         printf("%.16G %.16G %.16G %.16G\n",t, x[0], x[1], x[2]);
     }
     printf("\n");
   }
   return 0;
}
  
   