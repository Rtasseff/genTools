# python 2.7 script to test the pairwise function in statsUtil

import statsUtil
import numpy as np

print '****** test B B *******'
x = np.zeros(100,dtype=int)
y = np.zeros(100,dtype=int)
tmp = np.random.randn(100)
x[:] = 1
x[tmp>0]=2
tmp = np.random.randn(100)
y[:] = -1
y[tmp>0]=-2
r,p,_ = statsUtil.runPairwise(x,y,'B','B')
print 'r should be low'
print r
print 'p should be high'
print p
tmp = np.random.randn(100)
tmp[x==1] = 0
y[tmp>.1] = -1
print 'r should be moderate'
print r
print 'p should be moderate'
print p
tmp = np.random.randn(100)
tmp[x==2] = 0
y[tmp>.1] = -2
r,p,_ = statsUtil.runPairwise(x,y,'B','B')
print 'r should be high'
print r
print 'p should be sig'
print p
print '****** test B N *******'
y = np.random.rand(100)
r,p,_ = statsUtil.runPairwise(x,y,'B','N')
print 'r should be low'
print r
print 'p should be high'
print p
y[x==1] = y[x==1]+.1
r,p,_ = statsUtil.runPairwise(x,y,'B','N')
print 'r should be moderate'
print r
print 'p should be moderate'
print p
y[x==1] = y[x==1]+.4
r,p,_ = statsUtil.runPairwise(x,y,'B','N')
print 'r should be high'
print r
print 'p should be sig'
print p
y[x==1] = y[x==1]+10
r,p,_ = statsUtil.runPairwise(x,y,'B','N')
print 'r should be higher'
print r
print 'p should be very sig'
print p
print '****** test C N *******'
tmp = np.random.randn(100)
x[tmp>.666] = 3
y = np.random.rand(100)
r,p,_ = statsUtil.runPairwise(x,y,'C','N')
print 'r should be low'
print r
print 'p should be high'
print p
y[x==1] += .1
y[x==2] -= .1
r,p,_ = statsUtil.runPairwise(x,y,'C','N')
print 'r should be moderate'
print r
print 'p should be moderate'
print p
y[x==2] -= .4
y[x==1] += .4
r,p,_ = statsUtil.runPairwise(x,y,'C','N')
print 'r should be high'
print r
print 'p should be sig'
print p
print '****** test C C *******'
y = np.zeros(100,dtype=int)
tmp = np.random.rand(100)
y[tmp>.25] += 1
y[tmp>.5] += 1
y[tmp>.75] += 1
r,p,_ = statsUtil.runPairwise(x,y,'C','C')
print 'r should be low'
print r
print 'p should be high'
print p
tmp = np.random.rand(100)
tmp[x==1] = 0
y[tmp>.2] = 2
r,p,_ = statsUtil.runPairwise(x,y,'C','C')
print 'r should be high'
print r
print 'p should be sig'
print p
tmp = np.random.rand(100)
tmp[x==2] = 0
y[tmp>.2] = 1
r,p,_ = statsUtil.runPairwise(x,y,'C','C')
print 'r should be higher'
print r
print 'p should be very sig'
print p

