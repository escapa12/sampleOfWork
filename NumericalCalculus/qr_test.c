#include<stdio.h>
#include<math.h>

/*
 * Per compilar:
        gcc -o qr_test -g -Wall qr_test.c -lm
 */

int qrres (int m, int n, double *a, double *dr, double *b, double *x);
void printm(double*a,int m,int n);

int main(){
  int m=3,n=2;
  double a[m*n-1],x[n-1],dr[n-1];
  double b[m-1];
  /*double *b;
  b=NULL;*/
  a[0]=0;  a[1]=0;  a[2]=5;  a[3]=-4;  a[4]=0;  a[5]=-2;  
  b[0]=1;  b[1]=3;  b[2]=2;

  qrres(m,n,a,dr,b,x);
  
  double r[m*n-1];
  int i,j;
  for (i=0; i<m;i++){         //escriu la matriu R del QR
    for(j=0; j <n; j++){
      if(i==j) r[i+j*m]=dr[i];
      if(i<j) r[i+j*m]=a[i+j*m];
      if(j<i) r[i+j*m]=0;
    }
  }
  printm(a,m,n);
  printm(r,m,n);
  printf("b=%.16G, %.16G, %.16G\n",b[0],b[1],b[2]);
  printf("x=%.16G, %.16G\n",x[0],x[1]);
  return 0;
}

int qrres(int m,int n,double *a, double *dr, double *b, double*x){   //matriu a es guarda per columnes!! -->a[i+j*m]= a[i][j]
  if(m<n){printf("Error QR: m<n!\n");}
  int sig,i,j,k=0;
  double s,alfa,beta;
  while(k<n){
    s=0;
    sig=1;
    if(a[(m+1)*k]<0) sig=-1;
     for(i=k;i<m;i++){
        s+=a[i+k*m]*a[i+k*m];
     }
    s=sig*sqrt(s);
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
 printm(a,m,n);


  return 0;
}

 void printm(double*a,int m,int n){
   int i, j;
  for (i=0; i<m;i++){
    for(j=0; j <n; j++){
      printf("%.16G ",a[i+j*m]);
    }
    printf("\n");
  }
    printf("\n");
 }
