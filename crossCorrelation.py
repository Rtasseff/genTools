import numpy as np
import genTools.dataUtil as du

def crossCorrUniform(x1,x2,tstop=0):
	"""Calculate the cross correlation between
	two dynamic trajectories
	x_i,t^hat = x_i,t-mean(x_i)
	cc(\tau) = 1/N \sum_t (x_1,t^hat* x_2,t+\tau^hat) 
	Assumes the two trajectories corrispond to the same
	EVENLY spaced time seris.
	x's	trajectory, 1-d numpy array
	return 	cc, cross-corrilation as a function
		of \tau, 0:del t:stop t
		del t is the time step in the 
		trajectory and stop t is 1/2 of the trajectory.
	return	ccSdErr, an approximated standard error of cc
	""" 
	# get size
	n = len(x1)
	if n!=len(x2):
		raise ValueError('len(x1) must equal len(x2)')

	# pre process
	x1Hat = x1-np.mean(x1)
	x2Hat = x2-np.mean(x2)

	if tstop==0:
		tstop = int(.75*n)

	cc = np.zeros(tstop)
	ccSdErr = np.zeros(tstop)

	for t in range(tstop):
		summ = 0;
		summSq = 0;
		count = 0;
		for ind in range(n-t):
			summ = summ + x1Hat[ind]*x2Hat[ind+t]
			summSq = summSq + (x1Hat[ind]*x2Hat[ind+t])**2
			count = count+1
		cc[t] = summ/float(count)
		var = summSq/float(count) - cc[t]**2
		ccSdErr[t] = np.sqrt(var/float(count))

	return cc,ccSdErr
		
def crossCorr(traj1,traj2,tstep=0,taumax=0):
	"""Performs corss correlation between two 
	dynamic trajectories with non uniform time spacing,
	and possibally diffrent time courses.
	will only use the overlaping time interval.
	will use a step size for tau equal to 
	the smallest step size in the data
	and will apply linear interpolation to 
	compute a uniform time course.
	A preprocessing step is used to 
	get a unique time seris for each 
	trajectory, meaidan value is used.
	"""
	# get unique time courses
	t1,x1_tmp = du.aveUnqTime(traj1[0],traj1[1],'median')
	t2,x2_tmp = du.aveUnqTime(traj2[0],traj2[1],'median')
	# find common range
	tmin = np.max(np.array([t1[0],t2[0]]))
	tmax = np.min(np.array([t1[-1],t2[-1]]))
	# we could use any step, but no reason to go lower
	# then smallest interval
	if tstep==0:
		tstep = np.min(np.array([np.min(t1[1:]-t1[:-1]),np.min(t2[1:]-t2[:-1])]))
	
	# the range to do linear interp
	t = np.arange(tmin,tmax,tstep)
	x1 = np.zeros(len(t))
	x2 = np.zeros(len(t))
	# do linear interp
	for i in range(len(t)):
		x1[i] = du.getOrInterp(t1,x1_tmp,t[i])
		x2[i] = du.getOrInterp(t2,x2_tmp,t[i])

	# calculate the number of tau steps if a max is specified
	if taumax>0:
		# how many steps to get to max?
		tausteps = int(taumax/tstep)+1

	else:
		tausteps = 0
	
	# now we have an evenly spaced time seris for two varriables
	# run the corss corelation
	cc,_ = crossCorrUniform(x1,x2,tstop=tausteps)
	tau = np.arange(float(len(cc)))
	tau = tau*tstep
	return(cc,tau)

	

			
		

	


