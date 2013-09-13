#!/usr/bin/env python
#
# 
#     Copyright (C) 2003-2012 Institute for Systems Biology
#                             Seattle, Washington, USA.
# 
#     This library is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 2.1 of the License, or (at your option) any later version.
# 
#     This library is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public
#     License along with this library; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
# 
# 20120904 RAT
import numpy as np
import scipy.stats as stats
#import rpy2.robjects as robjects
import scipy.linalg
import operator as op
import itertools

def iqr(x):
	"""returns the inter quantile distance for
	the varriables in x.  
	if x is a matrix we assum the columns are
	difrent distributions and the iqr returned
	will be a vector corrisponding to the columns of x.
	"""
	n = len(x)
	if len(x.shape)==1:
		xSort = np.sort(x)
		q1 = xSort[int(.25*n)]
		q3 = xSort[int(.75*n)]
	elif len(x.shape)==2:
		n,m = x.shape
		q1 = np.zeros(m)
		q3 = np.zeros(m)
		for i in range(m):
			xSort = np.sort(x[:,i])
			q1[i] = xSort[int(.25*n)]
			q3[i] = xSort[int(.75*n)]

	return(q3-q1)

def fdr_bh(p,alpha=.05):
	"""Performs the Benjamini & Hochberg 1995
	multiple test correction for controlling
	the false discovery rate in familywise 
	analysis.  Tests must be independent or 
	positivly corrilated.
	p	original pvalues, np 1d array
	alpha 	threshold FDR, scalar float
	returns	h, regect or accept, 1d np bool array
	returus	p_adj, adjusted pvalues, 1d np array
	returns	pCrit, the critial p-value cut off
	"""
	m = len(p)
	sortInd = np.argsort(p)
	pSort = p[sortInd]
	unsortInd = np.argsort(sortInd)
	pAdj = np.zeros(m)*np.nan
	gamma = (np.arange(m)+1)*(alpha/m)
	pTmp = m*pSort/(np.arange(m)+1)
	for i in range(m):
		pAdj[i] = np.min(pTmp[i:])

	pAdjUnsort = pAdj[unsortInd]
	rejSort = pSort<=gamma
	# find the largest value still under threshold
	# note they are sorted
	maxInd = np.sum(rejSort)-1
	if maxInd<0:
		pCrit = 0
	else:
		pCrit = pSort[maxInd]

	h = p<=pCrit
	
	return h,pAdjUnsort,pCrit

def fisherComb(p):
	"""Apply fisher method to combine 
	the pvaules in the np array p, to
	get a single pvalue and return"""
	n = len(p)
	c = -2*np.sum(np.log(p))
	pComb = 1-stats.chi.cdf(np.sqrt(c),2*n)
	return pComb


def mi(bins,calcP=False,v=False):
	"""calculates mutual information 
	\sum_i\sum_j p(i,j)*log_2[p(i,j)/(p(i)*p(j))]
	for 2 variable histogram defined by bins
	bins	int 2d np array
	returns MI
	"""
	tot = float(np.sum(bins))
	n,m = bins.shape
	p_x = np.sum(bins,1)/tot
	p_y = np.sum(bins,0)/tot
	null_chi = 0
	mi = 0
	for i in range(n):
		for j in range(m):
			p_xy = bins[i,j]/tot
			frac = p_xy/(p_x[i]*p_y[j])
			if frac>0:
				mi = mi + p_xy*np.log2(frac)
			elif v:
				print 'warning some bins are empty, skipping them'
			if calcP:
				E_xy = tot*(p_x[i]*p_y[j])
				if E_xy == 0 and v: print 'WARNING 0 expectation in pvalue calc'
				if E_xy<1 and v: print 'WARNING expectation < 1 in pvalue calc, not good'
				if E_xy<5 and v: print 'WARNING expectation < 5 in pvalue calc, Cochran criterion may be violated'
				null_chi = null_chi + (bins[i,j]-E_xy)**2/E_xy

	if calcP:
		v = (n)*(m)
		p = 1-stats.chi.cdf(np.sqrt(null_chi),v)
		return mi,p
	else: return mi


def sampleIndexWR(n):
	"""Given an array size, n, returns
	the indicies for random uniform sampling 
	of the array, with replacment.
	"""
	indFloat = sampleWR(range(n),n)
	ind = map(int,indFloat)
	return(ind)	

def sampleWR(pop,k=0):
	"""given a population, pop,
	return a new sample, size=k,
	by sampling pop with replacment
	"""
	if k<1:k=len(pop)
	n = len(pop)
	sel = np.zeros(k)
	for i in range(k):
		index = np.random.randint(0,n)
		sel[i] = pop[index]

	return sel

def cmds(D):
	"""Preform clasical multidimensional scaling.
	using D as the pair-wsie distances Matrix
	returns Y the column matrix defining the
	new corrdinates that retain maximum conserved 
	distance and E the eigan values that describe
	the contribution of each dimension (Y col).
	"""
	[n,m] = D.shape
	eps = 1E-21
	if n!=m:
		raise ValueError('Wrong size matrix')

	# Construct an n x n centering matrix
	# The form is P = I - (1/n) U where U is a matrix of all ones
	P = np.eye(n) - (1/float(n) * np.ones((n,n)))
	B = np.dot(np.dot(P,-.5*D**2),P)
	# Calculate the eigenvalues/vectors
	[E, V] = scipy.linalg.eig((B+B.T)/2) # may help with round off error??
	E = np.real(E) # these come out as complex but imaginary part should be ~eps
	V = np.real(V) # same argument here
	# sort that mo fo
	ind = np.argsort(E)[::-1]
	E = E[ind]
	V = V[:,ind]
	# lets now create our return matrix 
	if np.sum(E>eps)==0:
		Y = 0
	else:
		Y = V[:,E>eps]*np.sqrt(E[E>eps])
	return(Y,E)


def makeDesingMat(y):
	"""transfomr a class label vector into the designe / indicator matrix Y
	y	class vector corrisponding to observations
	return	matrix with observations on rows and
		classes on cols, with 1 indicating obs in class 
	"""
	if np.sum(np.unique(y) != np.array(range(len(np.unique(y)))))>0:
		raise ValueError('Class vector y must be a numeric vector, with values as 0,1,2,...')
	return(np.eye(len(np.unique(y)))[y,:])


def cat2int(y):
	"""change a set of str catigory lables with n
	unique values into an int vector with unique 
	values of numbers from 0 to n-1
	returns the new int vector
		a list of catigories corrisponding to int value
	"""
	unique = list(set(y))
	yNew = np.array(np.zeros(len(y)),dtype=int)
	for i in range(len(unique)):
		yNew[y==unique[i]]=i

	return(yNew,unique)

def softThresh(X,thresh):
	"""Take X, and do soft 
	threshold on it and return that.
	"""
	tmp = np.abs(X) - thresh
	tmp[tmp<0] = 0
	return(tmp*np.sign(X))

def enrich(n_pos,n_draw,total_pos,total_all):
	"""Standard enrichment test usign hypergeometric 
	distribution.
	ie I got n_pos red balls out of n_draw draws from
	a bag with total_pos red balls and total all balls of 
	any color, calculate the probability of drawing n_pos
	red balls at random.
	"""
	
	p = stats.hypergeom.sf(n_pos-1,total_all,total_pos,n_draw)
	return(p)


def rankSum(x,y,forceExact=False):
	"""This is a two-sided Wilcoxon rank sum test,
	robust and nonparametric way to estimate if 
	two samples x and y have diffrent medians.
	Test assumes that two samples x,y are independent but 
	from distributions with the same median.  
	This is tested agains the alternate hyp that 
	the two samples are indipendent but come from 
	distributions with diffrent medians.
	Identical to Mann-Whitney U test.
	By defult the normal approximation is
	made for the test statistic IF the 
	sum of the samples has 15 or more obs;
	ortherwise, we enumerate the combinations
	to find the exact p-value.
	If forceExact==True then the exact p-value
	is calculated regardless of sample size.
	returns 
	p-value - sacalr
	z, approximate test statistic - scalar

	This function is similar to matlab ranksum.
	"""
	# find smallest 
	xLen = len(x)
	yLen = len(y)
	if xLen<=yLen:
		n = xLen
		sampSm = x
		sampLrg = y
	else:
		n = xLen
		sampSm = x
		sampLrg = y
	
	
	sampComb = np.append(sampSm,sampLrg)
	# check for no varriation
	if np.var(sampComb) == 0: return(np.nan,np.nan)
	# get ranks, possibly tied
	ranks,tieAdj = tiedRank(sampComb)
	ranks = ranks +1
	# get the stat
	xRank = ranks[:n]
	w = np.sum(xRank)
	# calculate the z statistic
	cor = 2*tieAdj/( (xLen+yLen)*(xLen+yLen-1) )
	wVar = xLen*yLen*((xLen+yLen+1) - cor)/12.
	wMean = n*(xLen+yLen+1)/2.
	wC = w-wMean
	z = (wC - .5 * np.sign(wC))/np.sqrt(wVar)

	# more complex methods exist, but I am using this 
	# simple way to deal with small samples
	if forceExact or (xLen+yLen)<15:
		allComb = chooseAllComb(ranks,n)
		allCombSum = np.sum(allComb,1)
		allCombLen = len(allCombSum)
		plo = np.sum(allCombSum<=w)/float(allCombLen)
		phi = np.sum(allCombSum>=w)/float(allCombLen)
		p = np.min([plo,phi]) 
		p = np.min([2*p,1])
	else:
		p = 2 * stats.norm.cdf(-np.abs(z))

	
		
	return(p,z)

def rankSum1S(x,y,forceExact=False):
	"""This is a one-sided Wilcoxon rank sum test,
	robust and nonparametric way to estimate if 
	the median of sample x > y.
	Test assumes that x,y are independent.  
	This is tested aginst the alternate hyp that 
	the two samples are indipendent but come from 
	distributions with the same medians.
	Identical to a one-sided Mann-Whitney U test.
	By defult the normal approximation is
	made for the test statistic IF the 
	sum of the samples has 15 or more obs;
	ortherwise, we enumerate the combinations
	to find the exact p-value.
	If forceExact==True then the exact p-value
	is calculated regardless of sample size.
	returns 
	p-value - sacalr
	z, approximate test statistic - scalar

	This function is similar to matlab ranksum.
	"""
	# find smallest 
	xLen = len(x)
	yLen = len(y)
	n = xLen
	
	sampComb = np.append(x,y)
	# check for no varriation
	if np.var(sampComb) == 0: return(np.nan,np.nan)
	# get ranks, possibly tied
	ranks,tieAdj = tiedRank(sampComb)
	ranks = ranks +1
	# get the stat
	xRank = ranks[:n]
	w = np.sum(xRank)
	# calculate the z statistic
	cor = 2*tieAdj/( (xLen+yLen)*(xLen+yLen-1) )
	wVar = xLen*yLen*((xLen+yLen+1) - cor)/12.
	wMean = n*(xLen+yLen+1)/2.
	wC = w-wMean
	z = (wC - .5 * np.sign(wC))/np.sqrt(wVar)

	# more complex methods exist, but I am using this 
	# simple way to deal with small samples
	if forceExact or (xLen+yLen)<15:
		allComb = chooseAllComb(ranks,n)
		allCombSum = np.sum(allComb,1)
		allCombLen = len(allCombSum)
		p = np.sum(allCombSum>=w)/float(allCombLen)
	else:
		p = stats.norm.cdf(-z)

	
		
	return(p,z)


def tiedRank(x):
	"""Calculates ranks for vector x while considering 
	ties.  If values are the same then the rank assigned
	is the average of the ranks they would span.
	Also returns the tie adjustment used in some tests.
	"""
	n = len(x)
	# sor the values 
	ind = np.argsort(x)
	xSort = x[ind]
	# get the nan values, at the end
	nNan = np.sum(np.isnan(x))
	xLen = n-nNan
	
		
	# find all ties
	tie = xSort[:xLen-1] == xSort[1:xLen]	
	tieInd = np.arange(xLen-1)[tie]
	tieInd = np.append(tieInd,xLen+2)
	nTies = len(tieInd)

	rankTmp = np.arange(float(n)) # the nan does not seem to work with ints 
	rankTmp[xLen:] = np.nan
	tieAdj = 0
	count = 0
	while count<nTies-1:
		tieStart = tieInd[count]
		nTied = 2
		# count number of ties, multiple ties will have consecuitive numbers
		while tieInd[count+1] == tieInd[count]+1:
			count = count + 1
			nTied = nTied + 1
			
		tieAdj = tieAdj + nTied*(nTied-1)*(nTied+1)/2.
		# compute average 
		rankTmp[tieStart:tieStart+nTied] = np.sum(rankTmp[tieStart:tieStart+nTied])/nTied
		count = count+1

	# put in order
	rank = np.zeros(n)
	rank[ind] = rankTmp
	return(rank,tieAdj)

	

def nck(n, k):
	"""Calculate n choose k
	N!/K!(N-K)!.  No warning for 
	large numbers that will take 
	forever!"""
	k = min(k, n-k)
	if k == 0: return 1
	numer = reduce(op.mul, xrange(n, n-k, -1))
	denom = reduce(op.mul, xrange(1, k+1))
	return numer//denom

def chooseAllComb(V,k):
	"""choose all possible combinations of 
	k items from the vector V
	blows up fast so be careful
	"""
	return(np.array(list(itertools.combinations(V,k))))
	
