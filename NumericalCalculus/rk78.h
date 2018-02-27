int rk78 (double *t, double x[], double *h,
      double hmin, double hmax, double tol,
      int n, int (*camp)(int n, double t, double x[], double f[], void *prm),
      void *prm);