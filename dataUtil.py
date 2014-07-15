import numpy as np 
from operator import itemgetter

def aveUnqTime(time,X,method='mean'):
	"""Look for points at the same time 
	and average the all variables at that
	point.
	time	vector of time points 
	X	n obs by m vars data matrix
	method	are we doing the defulat mean or 'median'
	return	vector of unquie time points in order
	return	matrix corrisponding average varriables 
	"""
	timeUnq = np.array(list(set(time)))
	n = len(timeUnq)
	
	XUnq = np.zeros(X.shape)
	err = np.zeros(X.shape)
	for i in range(n):
		if method == 'mean': 
			XUnq[i] = np.mean(X[time==timeUnq[i]],0)
			err[i] = np.sqrt(np.var(X[time==timeUnq[i]],0))
		elif method == 'median': 
			XUnq[i] = np.median(X[time==timeUnq[i]],0)
			err[i] = np.sqrt(np.var(X[time==timeUnq[i]],0))
		else: print 'wrong method type'

	# sort them 
	ind = np.argsort(timeUnq)
		
	return timeUnq[ind], XUnq[ind], err[ind]

def interpMissing(t1,t2,x1,x2):
	"""find the missing observations, 
	ie time points, between 2 sets of data
	and fill in the points by interpolation,
	such that we have one set of observations 
	(time points) consistent between 2 sets of data.
	return the time and data vectors for full 
	data set.
	assumes the t vars are the ordered, non 
	repeateing time points for each set (assending)
	xs corrisponding data, vectors or matrix, must be
	same dimmensions.
	Noe extrapolation, those points will be removed
	"""
	# lets get a complete list of unique time points
	tNew = np.sort(np.array(list(set(np.append(t1,t2)))))
	# no extrap so find smallest t possible
	if t1[0]!=t2[0]:
		tMin = np.max([t1[0],t2[0]])
		# has to be in list, just find it 
		# could be very very close but not exact 
		# do to numerical errors
		ind = np.argmin(np.abs(tNew - tMin))
		# keep this and verything larger
		tNew = tNew[ind:]
		
	# no extrap so find largest t possible
	if t1[-1]!=t2[-1]:
		tMax = np.min([t1[-1],t2[-1]])
		# has to be in list, just find it 
		# could be very very close but not exact 
		# do to numerical errors
		ind = np.argmin(np.abs(tNew - tMax))
		# keep this and everything smaller
		tNew = tNew[:ind+1]

	n = len(tNew)
	# now we have our set of t values set up x's
	# cannot think of a better way to deal with vectors or matrices
	m = x1.shape
	if len(m)>1:
		x1New = np.zeros((n,m[1]))
		x2New = np.zeros((n,m[1]))
	else:
		x1New = np.zeros(n)
		x2New = np.zeros(n)

	for i in range(n):
		t = tNew[i]
		# lets find closest t
		x1New[i] = getOrInterp(t1,x1,t)
		x2New[i] = getOrInterp(t2,x2,t)

	return tNew, x1New, x2New
	
def getOrInterp(t1,x1,t):
	"""Sub method for interMissing
	done here because done twice above,
	May be useful elsewhere??
	given a set of time course data t1 and x1
	and wanting to get a value for t
	will return corrisponding values from x1
	or, if not avalibe, will interpolate.
	assumes that t1 is ordered not repeating 
	and t dose not lie outside the range of t1
	t	scalar
	t1	vector of avalible time points
	x1	corrisponding data for t1,
		times along rows and vars along cols
	"""
	
	ind = np.argmin(np.abs(t1-t))
	val = t1[ind]-t
	if np.abs(val)<1E-21:
		# same so just get value
		x = x1[ind]
	else:
		# diffrent so we need to do some interpolation
		# need indices for that:
		if val>0:
			# t1 is larger so we need this a before
			indLow = ind-1
			indHi = ind
		else:
			indLow = ind
			indHi = ind+1

		# now lets do the interpolation
		x = interp((t1[indLow],x1[indLow]),(t1[indHi],x1[indHi]),t)
	return x
 
	

def interp(pt1,pt2,x):
	"""Given 2 sets of points (x,y) and 
	a given point between x, find corrisponding y.
	The ponts are tuples for corrisponding x,y, and 
	can handel a vector of y values for a given point x
	and will return all y values via indipendent interploation 
	for disiered x.
	May handel vectors for other varriables, but not tested.
	
	
	"""
	x1,y1 = pt1
	x2,y2 = pt2
	if x<x1 or x>x2:
		raise ValueError('x point for interpolation must be between pt1[1] and pt2[1]')
	slope = (y2-y1)/(x2-x1)
	step = x-x1
	return y1+step*slope

def trapAUC(t,x):
	""" use trapizoid rule to calculate area under
	the curve x(t).
	t	vector of time points
	x	values corisponding to times t
		difrent times along the rows
		varriables along the columns
	NOTE multiple x values and t will be averaged
	"""
	t,x = aveUnqTime(t,x)
	tDelta = t[1:] - t[:-1]
	intAUC = (.5*tDelta*(x[1:] + x[:-1]).T).T
	return np.sum(intAUC,0)

def rank2(x):
	"""very inefficent method for 
	ranking based on 2 columns.  col
	one is the primary and col 2 is for ties
	"""
	n = len(x)
	stuff = []
	for i in range(n):
		stuff.append((i,x[i,0],x[i,1]))

	sortedStuff = sorted(stuff,key=itemgetter(1,2))
	
	ranks = np.zeros(n)
	for i in range(n):
		index,_,_ = sortedStuff[i]
		ranks[index]=i
	return ranks

def FD(x,t,method='back'):
	"""Preforms finite difference (back only) 
	to estimate the discrete derivatives 
	corrisponding to t.  To do this properly
	the values need to be in order and 
	no repeting values of t.  To this end
	a preprocessing step is implimented 
	to order and combine (mean) the x values.
	t	vector of time points
	x	values corisponding to times t
	returns 
	dxdt	derivatives of x
	tNew	t after the preprocessing, 
		corrisponds to dxdt
	"""
	t,x = aveUnqTime(t,x)
	if method=='back':
		dxdt = np.zeros(len(x))
		dx = x[1:]-x[:-1]
		dt = t[1:]-t[:-1]
		dxdt[1:] = dx/dt
		dxdt[0] = dxdt[1]
	else:
		raise ValueError('method not implemented')
	
	return(dxdt,t)

def rankAgg(X,w=[]):
	"""A simple rank aggrigation method.
	This is similar to Broda and modified to
	have geometrically decreasing weights.
	Additional weight terms can be included.
	X has n features which will be ranked over m tests
	We note that high positive numbers rank high.
	"""
	n,m = X.shape
	if len(w) ==0:
		w = np.ones(m)

	ranks = np.argsort(np.argsort(-1*X,0),0)
	scores = np.sum((1/(ranks+1.))*w,1)
	return(scores)

		
