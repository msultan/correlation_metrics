def kdistANN(N,ke,data):
    ##calculate the distance to the kth neighbor for every pt in the dataset 
    #using euclidean distances
    from scipy import spatial
    from scipy.spatial import KDTree
    kdist=numpy.zeros(N)
    kdTree = spatial.KDTree(data)
    for i in range(0,N):
        pt=np.array([data[i][0],data[i][1]])
        temp=kdTree.query(pt,k=ke+1,eps=0,p=2.0)
        kdist[i]=temp[0][k]

    return kdist

def count_dist_x(data,kdist,refPt,N):
    #This routine count how many points in the dataset data have 
    #distance lesser then dist to the point x.
    countX=0
    countY=0
    for i in range(0,N):
        dTempX=abs(data[i][0]-refPt[0])
        dTempY=abs(data[i][1]-refPt[1])
        if(dTempX < kdist):
            countX +=1
        if(dTempY < kdist):
            countY += 1
    return countX-1,countY-1


def MIstrict(N,k,kdist,data):
    import scipy.special

    out=np.zeros(N)
    #for every pt find its neighbors in x and y regime based on the distance to the
    #kth neighbor
    d=1
    cdy=cdx=numpy.pi**(d/2)/scipy.special.gamma(1+d/2)/2**d
    d=2
    
    cdxy=numpy.pi**(d/2)/scipy.special.gamma(1+d/2)/2**d
    for i in range(0,N):
        nx,ny=count_dist_x(data,kdist[i],data[i],N)
        out[i]=-scipy.special.psi(nx+1)-scipy.special.psi(ny+1)+log(cdx)+log(cdy)-log(cdxy)
    return out
    
    
def mutual_nearest_neighbors(N,k,data):
    #make the tree and calculate the distance to the kth neighbor
    kdist=kdistANN(N,k,data)
    #populate the nx+ny
    out=MIstrict(N,k,kdist,data)
    return np.mean(out)+scipy.special.psi(N)+scipy.special.psi(k)