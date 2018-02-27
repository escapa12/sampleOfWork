#include <math.h>

#define X1 x[0]
#define X2 x[1]
#define X1D f[0]
#define X2D f[1]

int oscilharm (int n, double t, double *x, double *f, void *prm) {

  X1D = X2;
  X2D = -X1;

   return 0;
}

#undef X1
#undef X2
#undef X1D
#undef X2D



