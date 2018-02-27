import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

### Strategy one - Solving linear system using spase matrices###
def str1(G,m=0.15): 
    n=G.shape[0]
    c=np.asarray(G.sum(axis=0)).squeeze()
    D=sparse.diags(np.where(c==0,0,1/c) )
    e=np.ones(n);
    x=spsolve(sparse.eye(n)-(1-m)*G.dot(D),e)
    return x/sum(x)

### Strategy two - Using power method###
def get_Dandz(n,G,m):
     ##Observe that the column $j$ of A is all zeros  if and only  that the column $j$ of G is all zeros. So there is no need to compute $A$.   
    z=np.empty(n)
    c=np.array(G.sum(axis=0))[0].astype(float)
    for i in range (n):
        if c[i] !=0 :
            c[i]= 1/c[i]
            z[i]=m/n
        else: ##if sum c[i]==0 means that the column i is 0.
            z[i]=1/n
    return c,z

def get_Dandz2(n,G,m):
    c=np.asarray(G.sum(axis=0)).squeeze()
    #we read c twice but 2 np.where is more efficient that  1 loop in Python
    z=np.where(c==0,1/n,m/n) 
    c=np.where(c==0,0,1/c) 
    return c,z

def str2(G,m=0.15,tol=1e-16): 
    n=G.shape[0]
    x_k1=np.empty(n); x_k1.fill(1/n) ##x0: equal rank for all pages
    x_k=np.zeros(n)
    D,z=get_Dandz2(n,G,m)
    while np.linalg.norm(x_k1-x_k,np.inf) > tol :
        x_k=x_k1
        x_k1= G.dot((1-m)*D*x_k)+np.dot(z,x_k) ##no need to construct e
    return(x_k1/sum(x_k1))

### Strategy three - Using power method without storing M ###
def get_L_c(n,G):
    indG=G.nonzero()
    L=[]
    c=np.empty(n,dtype=int)
    for j in range(n):
        inx=np.where(indG[1]==j)[0]
        L.append(indG[0][inx])
        c[j]=len(indG[0][inx])
    return(L,c)

def Mx(m,n,c,L,x):
    xc=x
    x=np.zeros(n)
    for j in range (n):
        if(c[j]==0):
            x=x+xc[j]/n
        else:
            for i in L[j]:
                x[i]+= xc[j]/c[j]
    x=(1-m)*x+m/n
    return(x,xc)

def str3(G,m=0.15,tol=1e-16): 
    n=G.shape[0]
    L,c=get_L_c(n,G)
    x_k1=np.empty(n); x_k1.fill(1/n) ##x0: equal rank for all pages
    x_k=np.zeros(n)
    while (np.linalg.norm(x_k1-x_k,ord=np.inf) > tol):
        x_k1,x_k= Mx(m,n,c,L,x_k1)
    return( x_k1/sum(x_k1))


##BONUS:  improvements trying to perform better loops in Mx

import itertools as it

def Mx2(m,n,c,L,x,inx0):
    xc=x
    x=np.zeros(n)
    x+=np.sum(xc[inx0]/n)   
    for i,j in L:
        x[i]+= xc[j]/c[j]
    x=(1-m)*x+m/n
    return(x,xc)

def str3bis(G,m=0.15,tol=1e-16): 
    n=G.shape[0]
    r,col=G.nonzero()
    L = list(it.zip_longest(*[r,col])) #we compute L with a fancy function
    c=np.asarray(G.sum(axis=0)).squeeze()                         
    inx0=np.flatnonzero(c==0)
    x_k1=np.empty(n); x_k1.fill(1/n) ##x0: equal rank for all pages
    x_k=np.zeros(n)
    while (np.linalg.norm(x_k1-x_k,ord=np.inf) > tol):
        x_k1,x_k= Mx2(m,n,c,L,x_k1,inx0)
    return( x_k1/sum(x_k1))
        

