import numpy as np
from scipy import sparse
import scipy.io as sio
from lib import str1,str2,get_Dandz2 

G=sio.mmread('Data/p2p-Gnutella30.mtx')
G=sparse.csr_matrix(G)

sol1=str1(G,m=0.15) #we take as exact solution the one we get with str1 and str2 with tol=1e-16
sol2=str2(G,tol=1e-16)

n=G.shape[0] ; m=0.15
x_k1=np.empty(n); x_k1.fill(1/n) ##x0: equal rank for all pages
x_k=np.zeros(n)
D,z=get_Dandz2(n,G,m)
e1=[np.linalg.norm(x_k1-sol1)]
e2=[np.linalg.norm(x_k1-sol2)]
while np.linalg.norm(x_k1-x_k,np.inf) > 1e-10 :
    x_k=x_k1
    x_k1= G.dot((1-m)*D*x_k)+np.dot(z,x_k)
    x_k1=x_k1/sum(x_k1)
    e1.append(np.linalg.norm(x_k1-sol1))
    e2.append(np.linalg.norm(x_k1-sol2))
l=len(e1)
tail=10
q1=0; q2=0
##if we just take the last approximation of q we may be unlucky and have a bad approximation
for i in range(tail):
    q1+=e1[l-tail+i]/e1[l-tail+i-1]
    q2+=e2[l-tail+i]/e2[l-tail+i-1]

q1=q1/tail; q2=q2/tail    
print("Rate of convergence with str1 as exact solution:", q1)
print("Rate of convergence with str2(tol=1e-16) as exact solution:", q2)   
    

