#include <math.h>
 
#define X1 x[0]
#define X2 x[1]
#define X3 x[2]
#define X1D f[0]
#define X2D f[1]
#define X3D f[2]
 
#define sgm par[0]
#define rho par[1]
#define bet par[2]

int lorenz (int n, double t, double *x, double *f, void *prm) {
  double *par = prm;
   X1D=sgm*(X2-X1);
   X2D=-X1*X3+rho*X1-X2;
   X3D=X1*X2-bet*X3;
return 0;
}
 
#undef X1
#undef X2
#undef X3
#undef X1D
#undef X2D
#undef X3D
#undef sgm 
#undef rho 
#undef bet 