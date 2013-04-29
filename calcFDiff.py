import numpy as np

from genTools import dataManip.aveUnqTime as aveUnqTime

def calcForDiff(time,X):
	"""Calulates the finite differance 
	derivative of variables wrt time
	using forward differenc (last point uses
	backward difference).  If multiple points
	for each time point then the fd is calculated  	
	using the diffence to the average value at the 
	nest time point.
	time	vector of time points 
	X	n obs by m vars data matrix
	"""
	# lets do the average first
	timeUnq,XUnq = aveUnqTime(time,X)
	# lets no calculate the fd for each time point 
	n = len(time)
	dXdt = np.zeros(X.shape)
	for i in range(len(time)):
		t = time[i]
		ind = np.argmin(np.abs(t-timeUnq))
		if ind < len(timeUnq)-1:
			tDel = timeUnq[ind+1] - t
			xDel = XUnq[ind+1] - X[i]
		else:
			tDel = t - timeUnq[ind-1]
			xDel = X[i] - XUnq[ind-1]
		
		dXdt[i] = xDel/tDel

	return dXdt

def do4Grp(time,grp,X):
	dXdt = np.zeros(X.shape)
	unqGrp = np.array(list(set(grp)))
	for i in range(len(unqGrp)):
		g = unqGrp[i]
		ind = grp==g
		dXdt[ind] = calcForDiff(time[ind],X[ind])

	return dXdt
		
	
		
		
	
