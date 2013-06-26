from __future__ import division
import posMutualCode
from posMutualCode import * 
import numpy as np

def printMessage():
	print"Testing the script for a set of gaussians"
	return

def testPositionMutualInformationCode():
	import matplotlib.pyplot as mlp

	tN=5
	n_frm=N=300
	k=5
	mean=[1,1]
	cov=[[1,0.9],[0.9,1]]
	
	printMessage()
	tD=np.zeros(tN)
	print "Testing Position MutualInformation Code for %d iterations" %tN
	dim=1
	for t in range(tN):
		#get random data	
		data=np.zeros((N,2))
		coords=np.random.multivariate_normal(mean,cov,size=n_frm)
		#set it to data
		data[:,0]=coords[:,0]
		data[:,1]=coords[:,1]
		#subtract the mean from the data
		#data=data-[np.mean(coords[:,0]),0,0,0,np.mean(coords[:,1]),0]
		#data[:,1]=data[:,3]=0
		#data=np.reshape(data.flatten(),(n_frm,6))
		mutual=mutual_nearest_neighbors(N,k,0,0,data)
		print t,mutual[-1]
		tD[t]=mutual[-1]


	detxy=np.linalg.det(np.nan_to_num(cov))
	print "Analytical Solution:",-0.5*np.log(detxy)
	print 'Average of %d Trial for %d datapoints is %f with standard deviation \
	%f'%(tN,N, np.mean(tD),np.std(tD))
	#mlp.plot(tD)
	#mlp.axhline(-0.5*np.log(det))
	#mlp.show()

	return
