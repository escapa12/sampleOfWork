import numpy as np
import scipy.linalg as spla
import scipy.sparse as spsp
import pandas as pd
from sparsesvd import sparsesvd
import time

df=pd.read_csv("RCsGoff.csv")
df.drop("gene",1, inplace=True)
Sample=df.columns
X=df.values
m,n=X.shape
t0=time.clock()
for i in range(m): ##mean of variables =0
    X[i,:]=X[i,:]-np.mean(X[i,:])

Y= spsp.csc_matrix(X.T/np.sqrt(n-1))

U,S,Vt= sparsesvd(Y, m) ###svd for sparse matrices

Var_component= np.array(S**2/sum(S**2))
NewData=np.matrix(np.dot(Vt,X)).T

with open('PCA-ex2-output.txt','wb') as f:
    i=0
    for line in NewData:
        f.write((Sample[i]+',').encode())
        np.savetxt(f, line, fmt='%.6f',delimiter=',',newline=',')
        f.write(str(Var_component[i]).encode())
        f.write('\n'.encode())
        i+=1
print("PCA done. Results at 'PCA-ex2-output.txt'.  Computacional time:",time.clock()-t0)