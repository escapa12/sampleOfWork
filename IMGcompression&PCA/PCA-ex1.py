import numpy as np
import scipy.linalg as spla

TOL=1e-14
X = np.loadtxt("exemple.dat")
m,n=X.shape
#n observations
#m variables
for i in range(m): ##for each variable:  mean of variable is 0
    X[i,:]=X[i,:]-np.mean(X[i,:])

### Using covariance ###
print("Using covariance matrix: \n")
Y=X.T/np.sqrt(n-1)
U,S,Vt=spla.svd(Y,full_matrices=False)
S=S[S>TOL]
print("   Variance for component:\n",S**2/sum(S**2))
print("   SD for component:\n",S)
print("   Dataset in the new PCA coordinates:\n",np.dot(Vt[0:len(S),:],X))

### test ###
C=np.dot(X,X.T)/(n-1)
lam,v=np.linalg.eig(C)
L=np.flip(np.sort(lam[:len(S)].real),0)

if (np.all(L/sum(L)-(S**2/sum(S**2)) <TOL)):
    print("\nResult tested \n")
else:
    print("\nTest failed !!!\n")
    print(L/sum(L))

### Using correlation ###
print("Using correlation matrix: \n")
for i in range(m):
    X[i,:]=X[i,:]/np.std(X[i,:])

Y=X.T/np.sqrt(n-1)
U,S,Vt=spla.svd(Y,full_matrices=False)
S=S[S>TOL]
print("   Variance for component:\n",S**2/sum(S**2))
print("   SD for component:\n",S)
print("   Dataset in the new PCA coordinates:\n",np.dot(Vt[0:len(S),:],X))