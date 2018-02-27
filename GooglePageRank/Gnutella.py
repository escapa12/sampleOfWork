import numpy as np
import scipy.io as sio
from scipy import sparse
from lib import str1,str2,str3,str3bis
import time

G=sio.mmread('Data/p2p-Gnutella30.mtx')
G=sparse.csr_matrix(G)
    
t0=time.clock()
sol1=str1(G)
tf=time.clock()
print("Computacional time of strategy 1:",tf-t0)
sol2=str2(G,tol=1e-16)
sort1=np.argsort(sol1)
sort2=np.argsort(sol2)

### Strategy 2: for different TOLS##
file=open('Results/nutella_str2.txt','w')
tol=1e-4
j=4
for i in range(13):
    t0=time.clock()
    PR=str2(G,tol=tol)
    t= time.clock()-t0
    err1=np.linalg.norm(sol1-PR)
    err2=np.linalg.norm(sol2-PR)
    ranking=np.argsort(PR)
    acc1=(sum(sort1==ranking))/36682
    acc2=(sum(sort2==ranking))/36682
    file.write('{} {} {} {} {} {} {}\n'.format(j, tol, t,err1,err2, acc1, acc2 ))
    tol=tol/10
    j+=1
file.close()

### Strategy 3: for different TOLS###

file=open('Results/nutella_str3.txt','w')
tol=1e-4
j=4
for i in range(13):
    t0=time.clock()
    PR=str3bis(G,tol=tol)
    t= time.clock()-t0
    err1=np.linalg.norm(sol1-PR)
    err2=np.linalg.norm(sol2-PR)
    ranking=np.argsort(PR)
    acc1=(sum(sort1==ranking))/36682
    acc2=(sum(sort2==ranking))/36682
    file.write('{} {} {} {} {} {} {}\n'.format(j, tol, t,err1,err2, acc1, acc2 ))
    tol=tol/10
    j+=1
file.close()

### Bonus: comparing computacional time with improvements###
t0=time.clock()
s3=str3(G,tol=1e-16)
tf=time.clock()
print("Time str3 with given code:", tf-t0, "Norm solution error:", np.linalg.norm(s3-sol2) )

t0=time.clock()
s3bis=str3bis(G,tol=1e-16)
tf=time.clock()
print("Time str3 with own code:", tf-t0, "Norm solution error:", np.linalg.norm(s3bis-sol2) )