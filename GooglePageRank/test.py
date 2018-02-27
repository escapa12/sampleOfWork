import numpy as np
from scipy import sparse
from lib import str1,str2,str3

G1 = sparse.csr_matrix([[0.,0.,1.,1.],[1.,0.,0.,0.],[1.,1.,0.,1.],[1.,1.,0.,0.]])
G2 = sparse.csr_matrix([[0, 0, 1,1,0], [1, 0, 0,0,0], [1, 1, 0,1,1],[1,1,0,0,0],[0,0,1,0,0]])
G3 = sparse.csr_matrix([[0, 0, 0,1], [1, 0, 0,0], [1, 1, 0,1],[1,1,0,0]])

### Strategy 1###
sol11=str1(G1)
sol12=str1(G2)
sol13=str1(G3)
### Strategy 2###
sol21=str2(G1)
sol22=str2(G2)
sol23=str2(G3)
### Strategy 3###
sol31=str3(G1)
sol32=str3(G2)
sol33=str3(G3)
### Matrix G1 with m=0###
sol11m=str1(G1,m=0)
sol21m=str2(G1,m=0)
sol31m=str3(G1,m=0)

A1 = np.array([[0, 0, 1,1/2], [1/3, 0, 0,0], [1/3, 1/2, 0,1/2],[1/3,1/2,0,0]])
A2 = np.array([[0, 0, 1/2,1/2,0], [1/3, 0, 0,0,0], [1/3, 1/2, 0,1/2,1],[1/3,1/2,0,0,0],[0,0,1/2,0,0]])
A3 = np.array([[0, 0, 0,1/2], [1/3, 0, 0,0], [1/3, 1/2, 0,1/2],[1/3,1/2,0,0]])

vaps,veps = np.linalg.eig(A1)
sol1=veps[:,0]/sum(veps[:,0])

n11=np.linalg.norm(sol1-sol11); n21=np.linalg.norm(sol1-sol21);n31=np.linalg.norm(sol1-sol31);  
n11m=np.linalg.norm(sol1-sol11m); n21m=np.linalg.norm(sol1-sol21m);n31m=np.linalg.norm(sol1-sol31m);

vaps,veps = np.linalg.eig(A2)
sol2=veps[:,0]/sum(veps[:,0])
n12=np.linalg.norm(sol2-sol12); n22=np.linalg.norm(sol2-sol22); n32=np.linalg.norm(sol2-sol32);  
#n12m=np.linalg.norm(sol2-sol12m); n22m=np.linalg.norm(sol2-sol22m); n32m=np.linalg.norm(sol2-sol32m);

vaps,veps = np.linalg.eig(A3)
sol3=veps[:,3]/sum(veps[:,3])
n13=np.linalg.norm(sol3-sol13); n23=np.linalg.norm(sol3-sol23);n33=np.linalg.norm(sol3-sol33);  

print("\n\nOriginal matrix of exercise 1:")
print("\nExact solution(np.linalg.eig):",sol1)
print("\nStrategy 1:",sol11, "\nNorm error:", n11)
print("\nStrategy 2:",sol21, "\nNorm error:", n21)
print("\nStrategy 3:",sol31, "\nNorm error:",n31)

print("\n\nSetting m=0:")
print("\nStrategy 1:",sol11m,  "\nNorm error:", n11m)
print("\nStrategy 2:",sol21m,  "\nNorm error:", n21m)
print("\nStrategy 3:",sol31m,  "\nNorm error:", n31m)

print("\n\nExtended matrix of exercise 1:")
print("\nExact solution(np.linalg.eig):",sol2)
print("\nStrategy 1:",sol12,  "\nNorm error:", n12)
print("\nStrategy 2:",sol22,  "\nNorm error:", n22)
print("\nStrategy 3:",sol32,  "\nNorm error:", n32)

print("\n\nMatrix of exercise 2 (dangling nodes):")
print("\nParron eigenvalue(np.linalg.eig):",sol3)
print("\nStrategy 1:",sol13, "\nNorm error respect Parron:", n13)
print("\nStrategy 2:",sol23, "\nNorm error respect Parron:", n23)
print("\nStrategy 3:",sol33, "\nNorm error respect Parron:", n33)