#!/usr/bin/env python
from __future__ import division
from scipy import spatial
from scipy.spatial import KDTree
import numpy as np 
import scipy.special
import numpy as np 

def kdistANN(N,kN,data,dim):
    ##calculate the distance to the kth neighbor for every pt in the dataset 
    #using euclidean distances
    pt=np.zeros((1,2*dim))
    kdist=np.zeros(N)
    kdTree = spatial.KDTree(data)
    for i in range(0,N):
        pt[0,:dim] = data[i][:dim] 
        pt[0,dim:] = data[i][dim:]  
        temp=kdTree.query(pt.flatten(),kN+1,'inf')
        kdist[i]=temp[0][kN]

    return kdist
 
def count_dist_x(data,kdist,refPt,N,dim):
    from scipy.spatial.distance import euclidean as euc
    #This routine count how many points in the dataset data have 
    #distance lesser then dist to the point x.
    countX=0
    countY=0
    for i in range(0,N):
        dTempX=euc(data[i][:dim],refPt[:dim])
        dTempY=euc(data[i][dim:],refPt[dim:])
        if(dTempX < kdist):
            countX +=1
        if(dTempY < kdist):
            countY += 1
    #-1 because otherwise we would count the actual pt as well. 
    return countX-1,countY-1


def MIstrict(N,k,kdist,data,dim):


    out=np.zeros(N)
    #for every pt find its neighbors in x and y regime based on the distance to the
    #kth neighbor
    d=dim
    cdy=cdx=np.pi**(d/2.0)/scipy.special.gamma(1+d/2.0)/2.0**d
    d=2*d
    cdxy=np.pi**(d/2.0)/scipy.special.gamma(1+d/2.0)/2.0**d
    correctionFactor=np.log(cdx)+np.log(cdy)-np.log(cdxy)
    
    for i in range(0,N):
        nx,ny=count_dist_x(data,kdist[i],data[i],N,dim)
        out[i]=-scipy.special.psi(nx+1)-scipy.special.psi(ny+1)+correctionFactor
    return out
    

def mutual_nearest_neighbors(N,k,i,j,data):
    dim = int(data.shape[1]/2)
    kdist=kdistANN(N,k,data,dim)
    out=MIstrict(N,k,kdist,data,dim)
    return (i,j,np.mean(out)+scipy.special.psi(N)+scipy.special.psi(k))