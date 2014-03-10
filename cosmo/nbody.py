import numpy as np
from matplotlib import pyplot as plt
from random import random

# function to calculate acceleration 
def a(r):
    G = 1E-4
    epsilon = 1E-4
    acc = np.zeros((3,3))
    
    for i in range(3):
        for j in range (3):
            if i == j: 
                continue
            else:
                acc[i] += G*(r[j]-r[i])/(np.linalg.norm(r[j]-r[i])+epsilon)**3
    return acc

h=.01 #stepsize
tmax=500 #max steps
steps=int(tmax/h)+2
rcurr= np.array([[1.500,1.00,1.00],[.750,1.133,1.001],[.751,.567,.799]]) #array with starting positions
vcurr= np.zeros((3,3)) #array with starting velocities (3 objects, 3 dimensions)

for i in range(2):
    for j in range(3):
        vcurr[i][j] = random()*1E-4
        vcurr[2][j] -= vcurr[i][j]

tcurr=0

acurr = a(rcurr)

rlist=np.zeros((3,steps,3))
vlist=np.zeros((3,steps,3))
tlist=np.zeros(steps)
rlist[:,0]=rcurr
vlist[:,0]=vcurr
tlist[0]=tcurr


n=1

while tcurr < tmax:
    
    rnext = rcurr + vcurr*h + .5*acurr*h**2
    
    anext = a(rnext)
    vnext = vcurr + .5*(acurr+anext)*h
    tnext = tcurr + h
    
       
    rlist[:,n]=rnext
    vlist[:,n]=vnext
    tlist[n]=tnext
    
    rcurr,vcurr,acurr,tcurr = rnext,vnext,anext,tnext
    
    n=n+1
    
    if n % 1000/h == 0: print 'finished t = %0.f' % tcurr
    
plt.clf()
for i in range(3):
    plt.plot(rlist[i].T[0],rlist[i].T[1])
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.axes().set_aspect('equal')
plt.show()
