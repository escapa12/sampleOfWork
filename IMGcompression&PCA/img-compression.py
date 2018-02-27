import scipy.linalg as spla
import scipy.ndimage as spimg
import scipy.misc as spmi
import time
import numpy as np

file_path1="./IMG/gauss.jpg"
file_path2="./IMG/pi.jpg"

K=[100,50,20,10,3]

def image_compresion_BW(file_path,K):
    t0=time.clock()
    IMG=spimg.imread(file_path)
    shape_img=np.shape(IMG)
    if len(shape_img)==3: ### hence, color image
        BW=np.empty((shape_img[0],shape_img[1]), dtype=int)
        for i in range(shape_img[0]):
            for j in range(shape_img[1]):
                BW[i,j]=IMG[i,j,0]
    else:                 #### black and white image
        BW=IMG
    U,S,Vt=spla.svd(BW,full_matrices=False)
    Frob_norm=np.sqrt(sum(S**2))
    for k in K:
        Ak=np.dot(U[:,0:k],np.einsum('ij,i->ij',Vt[0:k,:],S[:k])) ###no explicit construction of S. 
        prop=np.sqrt(sum(S[:k]**2))/Frob_norm
        spmi.imsave(file_path[:-4]+'BW'+str(prop)+'.jpg',Ak)
    print("\nBlack and white image compresion of "+file_path)
    print("Size original img:",shape_img)
    print("Compressed img: ", Ak.shape)
    print("Computacional time:",time.clock()-t0)
    
def image_compresion_color(file_path,K):
    t0=time.clock()
    IMG=spimg.imread(file_path)
    shape_img=np.shape(IMG)

    R=np.empty((shape_img[0],shape_img[1]), dtype=int)
    G=np.empty((shape_img[0],shape_img[1]), dtype=int)
    B=np.empty((shape_img[0],shape_img[1]), dtype=int)
    
    for i in range(shape_img[0]):
        for j in range(shape_img[1]):
            R[i,j]=IMG[i,j,0]
            G[i,j]=IMG[i,j,1]
            B[i,j]=IMG[i,j,2]

    RU,RS,RVt=spla.svd(R,full_matrices=False)
    GU,GS,GVt=spla.svd(G,full_matrices=False)
    BU,BS,BVt=spla.svd(B,full_matrices=False)

    Frob_norm=np.sqrt(sum(RS**2))+np.sqrt(sum(GS**2))+np.sqrt(sum(BS**2))
    K=[100,50,20,10,3]
    for k in K:
        RAk=np.dot(RU[:,0:k],np.einsum('ij,i->ij',RVt[0:k,:],RS[:k])) ###no explicit construction of RS. 
        GAk=np.dot(GU[:,0:k],np.einsum('ij,i->ij',GVt[0:k,:],GS[:k])) ###no explicit construction of GS.
        BAk=np.dot(BU[:,0:k],np.einsum('ij,i->ij',BVt[0:k,:],BS[:k])) ###no explicit construction of BS.

        prop=(np.sqrt(sum(RS[:k]**2))+np.sqrt(sum(GS[:k]**2))+np.sqrt(sum(BS[:k]**2)))/Frob_norm

        IMG_compr = np.empty((np.shape(RAk)[0],np.shape(RAk)[1], 3))
        IMG_compr[:,:,0]=RAk; IMG_compr[:,:,1]=GAk; IMG_compr[:,:,2]=BAk; 

        spmi.imsave(file_path[:-4]+str(prop)+'.jpg',IMG_compr)

    print("\nColor image compresion of "+file_path)
    print("Size img:",shape_img)
    print("Compressed img: ", IMG_compr.shape[:2])
    print("Computacional time:",time.clock()-t0)

image_compresion_BW(file_path1,K)
image_compresion_BW(file_path2,K)
image_compresion_BW('./IMG/galois.jpg',K)

image_compresion_color(file_path1,K)
image_compresion_color('./IMG/uni.jpg',K=[400,100,50,5])

