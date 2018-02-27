import numpy as np
import scipy.io as sio
from scipy import sparse
from lib import str2,str3bis
import time

G=sio.mmread('Data/web-Google.mtx')
G=sparse.csr_matrix(G)

def str1bis(G,m=0.15,tol=1e-8): 
    n=G.shape[0]
    c=np.asarray(G.sum(axis=0)).squeeze()
    D=sparse.diags(np.where(c==0,0,1/c) )
    e=np.empty(n); e.fill(1.)
    x,success=sparse.linalg.bicgstab(sparse.eye(n)-(1-m)*G.dot(D),e)
    if(success==0):
        return x/sum(x)
    else:
        print("str1 did not converge")

t0=time.clock()
sol1=str1bis(G)
tf=time.clock()
print("\nStrategy bicgstab:\n    Time:",tf-t0,"sec")
t0=time.clock()
sol2=str2(G,tol=1e-8)
tf=time.clock()
print("\nStrategy 2:\nTol: 1e-8\n    Time:",tf-t0,"sec")

t0=time.clock()
sol3=str3bis(G,tol=1e-8)
tf=time.clock()
print("\nStrategy 3 :\nTol: 1e-8\n    Time:",tf-t0,"sec")

