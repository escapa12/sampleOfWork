#include <math.h>

#define X1 x[0]
#define X2 x[1]
#define X1D f[0]
#define X2D f[1]

#define AD(i,j) f[2+(j)*2+(i)]
#define A(i,j) x[2+(j)*2+(i)]

int pendol (int n, double t, double *x, double *f, void *prm) {

   X1D=X2;
   X2D=-sin(X1);
   
   if(n>2){
    
   double ydx=-cos(X1);
   AD(0,0)= A(1,0);
   AD(1,0)=A(0,0)*ydx;
   AD(0,1)= A(1,1);
   AD(1,1)=A(0,1)*ydx ;
   }
   return 0;
}

#undef X1
#undef X2
#undef X1D
#undef X2D
#undef AD
#undef A

