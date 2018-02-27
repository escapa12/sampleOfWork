#include<stdio.h>
#include<math.h>

#define SGN(x) ((x)>=0?1:-1)

void qrres(int m,int n,double *a, double *dr, double *b, double*x){
  if(m<n){printf("Error: m<n!\n");/*return -1;*/}
  int i,j,k=0;
  double s,alfa,beta;
  while(k<n){
    s=0;
     for(i=k;i<m;i++){
        s+=a[i+k*m]*a[i+k*m];
     }
    s=SGN(a[(m+1)*k])*sqrt(s);
    beta=1/(s*s+s*a[(m+1)*k]);
    a[(m+1)*k]+=s; 
    dr[k]=-s;
    for(j=k+1; j<n; j++){
       alfa=0;
    
       for(i=k;i<m; i++) {
	 alfa+=a[i+k*m]*a[i+j*m];
      }
       alfa*=beta;
       for(i=k; i<m ; i++){
	 a[i+j*m]-=alfa*a[i+k*m];
       }
    }
   if(b!=NULL){
   alfa=0;
    for(i=k;i<m;i++)
      alfa+=a[i+k*m]*b[i];
    alfa*=beta;
      for(i=k; i<m; i++)
         b[i]-=alfa*a[i+k*m];
    }
    k++;
  }
  if(b!=NULL){
    for(i=n-1; i>=0; i--){
       alfa=0;
       for(j=i+1;j<n;j++){
          alfa+=a[i+j*m]*x[j];
       }
       x[i]=(b[i]-alfa)/dr[i] ;
     }
  }
}