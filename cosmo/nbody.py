import numpy as np
from matplotlib import pyplot as plt
from random import random

# function to calculate acceleration 
def a(r):
    G = 1E-4
    epsilon = 1E-3
    acc = np.zeros((3,3))
    
    for i in range(3):
        for j in range (i+1,3):
            acc[i] += G*(r[j]-r[i])/(np.linalg.norm(r[j]-r[i])+epsilon)**3
            acc[j] -= acc[i]
    return acc

h=.01 #stepsize
tmax=1000 #max time to prevent infinite loop if particles don't exit box
#steps=int(tmax/h)+2
rcurr= np.array([[1.500,1.00,1.00],[.750,1.133,1.001],[.751,.567,.799]])-1 #array with starting positions
vcurr= np.zeros((3,3)) #array with starting velocities (3 objects, 3 dimensions)

#for now no random velocities for reproducability
#for i in range(2):
#    for j in range(3):
#        vcurr[i][j] = random()*1E-5
#        vcurr[2][j] -= vcurr[i][j]

tcurr=0
acurr = a(rcurr)

rlist=np.array([])
vlist=np.array([])
tlist=np.array([tcurr])

rlist=np.append(rlist,rcurr,1)
print rlist
vlist=vcurr
tlist=np.append(tlist,tcurr)

vnext=np.zeros((3,3))+1

print vlist,'\n\n',vnext

exit(0)

n=1

while np.amax(abs(rcurr)) < 1 and tcurr < tmax:
    
    rnext = rcurr + vcurr*h + .5*acurr*h**2
    
    anext = a(rnext)
    vnext = vcurr + .5*(acurr+anext)*h
    tnext = tcurr + h
    
       
    rlist[:,n]=rnext
    vlist[:,n]=vnext
    tlist[n]=tnext
    
    rcurr,vcurr,acurr,tcurr = rnext,vnext,anext,tnext
    
    n += 1
    
    if n % 1000/h == 0: print 'finished t = %0.f' % tcurr
    
    
plt.clf()
for i in range(3):
    plt.plot(rlist[i].T[0],rlist[i].T[1])
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.axes().set_aspect('equal')
plt.show()
