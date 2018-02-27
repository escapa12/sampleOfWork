#include <stdio.h>
#include <math.h>
#include "flux.h"
#include "qrres.h"

#define N 2*m
#define M N*(1+N)
#define A(i,j) y[N+(j)*N+(i)]
#define DG(i,j) dg[(j)*N+(i)]
#define DGaux(i,j) DGaux[(j)*N+(i)]

  void y_0(double *g,double *y,int m){    //escriu les condicions inicials de y per fer Dflux
	int i,j;
	for(i=0; i<N; i++) y[i]=g[i];    //y_1...y_2m=x_0
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){      //punt inicial de la part de les variacionals-->identitat
			if(i==j) A(i,j)=1;
			else A(i,j)=0;
		}
	}
}

int cmani_gdg (int m, double x0[], double xf[], double dt, double dv[],double g[],
              double dg[],double pas0, 	double pasmin, double pasmax, double tolfl, int npasmx,
	          int (*camp)(int n, double t, double x[], double f[], void *prm),void *prm){
//avaluacio de G i DG simultania
//calcular Dflux requereix vector de 2m+2m*2m:=M per contenir les variacionals
  int i,j,k;
  double t=0, y[M] , DGaux[2*m*m];
  
  for(i=0; i<m; i++) g[i]=x0[i];  
  for(i=0; i<m; i++) g[m+i]=x0[m+i]+dv[i];  //x_0.5
  y_0(g,y,m);          //prepara les condicions inicials d y per calcular DG: y_1...y_2m=x_0.5 y_M=Id
  if(flux(&t,y,&pas0,dt/2,pasmin,pasmax,tolfl,npasmx,M,camp,prm)!=0) return-1;  //necessitarem les m últimes columnes per trobar les m primeres columes de DG
  for (i=0; i<N ; i++){             
	  for (j=m; j<N ; j++){     //guardem les m ultimes columnes de y a DGaux
		DGaux(i,j-m)= A(i,j);  
	  }
  }
  if(flux(&t,g,&pas0,dt/2,pasmin,pasmax,tolfl,npasmx,N,camp,prm)!=0) return -1;  //x_1
  for(i=0; i<m; i++) g[m+i]+=dv[m+i]; //x_1.5 
  y_0(g,y,m);//prepara les cond inicials de y per calcular DG
  if(flux(&t,y,&pas0,dt/2,pasmin,pasmax,tolfl,npasmx,M,camp,prm)!=0) return -1; //les m últimes columnes de y són les m últimes columes de DG
  for (i=0; i<N ; i++){ 
	for (j=m; j<N ; j++){
		DG(i,j)= A(i,j);
	}
  }
  if(flux(&t,g,&pas0,dt/2,pasmin,pasmax,tolfl,npasmx,N,camp,prm)!=0) return -1;  //x_2
  for(i=0; i<N; i++) g[i]-=xf[i];
  for (i=0; i<N ; i++){   //primeres m columnes de DG= Dflux(x_1.5)*DGaux (producte matricial)NxN*Nxm-->Nxm
	  for (j=0; j<m ; j++){
		DG(i,j)=0;
		for(k=0; k<N; k++)
			DG(i,j)+= A(i,k)*DGaux(k,j);
	}
  }
	return 0;
  }

  int cmani (int m, double x0[], double xf[], double dt, double dv[],double tol, int maxit,double pas0, double pasmin, 
	double pasmax, double tolfl, int npasmx,int (*camp)(int n, double t, double x[], double f[], void *prm),void *prm){
  int j,i=0;
  double dx[N],dr[N],dg[4*m*m],g[N],ng;
  while (i<maxit){ //els comentaris que contenen printf són per els prints utilitzats per fer el test cmani_pendol
      /*met Newton;
      *resoldre dg(dv)*dx= g(dv)
      *dv-=dx;*/
      if (cmani_gdg(m,x0,xf,dt,dv,g,dg,pas0,pasmin,pasmax,tolfl,npasmx,camp,prm)!=0 ) return -1; //  G(dv), DG(dv)
     ng=0;
     for(j=0; j<N; j++) ng+=g[j]*g[j];
     ng= sqrt(ng);
//     printf("it %d: ng:%.16G ",i,ng); per a cmani_pendol
     if(ng<tol) { //printf("\n%.16G %.16G\n",dv[0],dv[1]); 
       return 0;}
     qrres(N, N, dg, dr, g, dx); //resol: dg(dv)*dx= g(dv) i guarda la solucio a les dx
     ng=0;
     for(j=0; j<N; j++) ng+=dx[j]*dx[j];
     ng= sqrt(ng);
//     printf("nc:%.16G \n",ng); per a cmani_pendol
     if(ng<tol) { //printf("\n%.16G %.16G\n",dv[0],dv[1]);
       return 0;}
     for(j=0;j<N;j++) dv[j]-=dx[j];
	  i++;
  } 
      printf("Error: S'ha assolit el nombre màxim d'iteracions del metode de Newton \n");
	  return -1;
}