
int cmani_gdg (int m, double x0[], double xf[], double dt, double dv[],double g[],double dg[],double pas0, 	double pasmin, double pasmax, double tolfl, int npasmx,
	      int (*camp)(int n, double t, double x[], double f[], void *prm),void *prm);



int cmani (int m, double x0[], double xf[], double dt, double dv[],double tol, int maxit,double pas0, double pasmin, double pasmax, double tolfl, int npasmx,
	   
	   int (*camp)(int n, double t, double x[], double f[], void *prm),void *prm);

void y_0(double *x0,double *y,int m);