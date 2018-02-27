#include <stdio.h>
#include <math.h>
 
#define X x[0]
#define Y x[1]
#define XD f[0]
#define YD f[1]
#define AD(i,j) f[2+(j)*2+(i)]
#define A(i,j) x[2+(j)*2+(i)]
#define r2 X*X-Y*Y

int campvar (int n, double t, double *x, double *f, void *prm){
double alpha;
alpha=*((double *)prm);

XD=alpha*(1-r2)*X-Y;
YD=alpha*(1-r2)*Y+X;

if (n>2){      //variacionals
   double dxdx,dxdy,dydx,dydy;
   dxdx=(1-3*X*X-Y*Y)*alpha;     dxdy=-2*X*Y*alpha-1;
   dydx= 1-2*X*Y*alpha;          dydy= (1-X*X-3*Y*Y)*alpha;
   AD(0,0)=A(0,0)*dxdx + A(1,0)*dxdy;
   AD(1,0)=A(0,0)*dydx + A(1,0)*dydy;
   AD(0,1)=A(0,1)*dxdx + A(1,1)*dxdy;
   AD(1,1)=A(0,1)*dydx + A(1,1)*dydy;
}
 return 0;
}
#undef X 
#undef Y 
#undef XD
#undef YD
#undef AD
#undef A
#undef r2