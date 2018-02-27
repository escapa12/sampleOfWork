#include <math.h>

/*!
 * \file rtbps.c
 * 
 * \brief RTBP en sinòdiques (posicions + velocitats)
 * 
 * NO té suport per derivades respecte de mu
 */

#define N 6
#define MU mu
#define M1 (1-MU)
#define M2 MU
#define U x[3]
#define V x[4]
#define W x[5]

int rtbps (int n, double t, double *x, double *f, void *prm) {
   double mu, xmmu, xmmup1, y, z, y2, z2, td1, td2, rh12, rh22, rh1, rh2,
	  rh13, rh23, pot3, potx3;
   mu=*((double *)prm);
   xmmu=x[0]-MU; xmmup1=xmmu+1;
   y=x[1]; y2=y*y; z=x[2]; z2=z*z;
   td1=y2+z2;
   rh12=xmmu*xmmu+td1; rh22=xmmup1*xmmup1+td1;
   rh1=sqrt(rh12); rh2=sqrt(rh22);
   rh13=rh12*rh1; rh23=rh22*rh2;
   td1=-M1/rh13; td2=-M2/rh23;
   pot3=td1+td2;
   potx3=td1*xmmu+td2*xmmup1;
/* Equacions */
   f[0]=U; f[1]=V; f[2]=W;
   f[3]=2*V+x[0]+potx3;
   f[4]=-2*U+y*(1+pot3);
   f[5]=z*pot3;
/* Variacionals primeres */
   if (n>N) {
#define VARX(i,j) x[(1+(j))*N+(i)]
#define VARF(i,j) f[(1+(j))*N+(i)]
      double pot5, potx5, omgxx, omgxy, omgxz, omgyy, omgyz, omgzz;
      int j;
      td1=M1/(rh13*rh12); td2=M2/(rh23*rh22);
      pot5=td1+td2;
      td1*=xmmu; td2*=xmmup1;
      potx5=td1+td2;
      omgxx=1+pot3+3*(td1*xmmu+td2*xmmup1);
      omgxy=3*y*potx5;
      omgxz=3*z*potx5;
      omgyy=1+pot3+3*y2*pot5;
      omgyz=3*y*z*pot5;
      omgzz=pot3+3*z2*pot5;
      for (j=0; j<N; j++) {
	 VARF(0,j)=VARX(3,j); VARF(1,j)=VARX(4,j); VARF(2,j)=VARX(5,j);
	 VARF(3,j)=omgxx*VARX(0,j)+omgxy*VARX(1,j)+omgxz*VARX(2,j)
		     +2*VARX(4,j);
	 VARF(4,j)=omgxy*VARX(0,j)+omgyy*VARX(1,j)+omgyz*VARX(2,j)
		     -2*VARX(3,j);
	 VARF(5,j)=omgxz*VARX(0,j)+omgyz*VARX(1,j)+omgzz*VARX(2,j);
      }
#undef VARX
#undef VARF
   }
   return 0;
}

#undef W
#undef V
#undef U
