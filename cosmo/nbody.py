import numpy as np
from matplotlib import pyplot as plt

# function to calculate acceleration 
def a(r):
    G = 1E-4
    acc = np.zeros((3,3))
    
    for i in range(3):
        for j in range(3):
            if i == j: 
                continue
            else:
                acc[i] += -G*r[i]/np.linalg.norm(r[j]-r[i])**3
    return acc

h=.1 #stepsize
tmax=100 #max steps
steps=int(tmax/h)+1
rcurr= np.array([[1.500,1.00,1.00],[.750,1.133,1.001],[.751,.567,.799]]) #array with starting positions
vcurr= np.zeros((3,3)) #array with starting velocities (3 objects, 3 dimensions)
tcurr=0

acurr=np.array(3)
for i in range(3):
    acurr = a(rcurr)


rlist=np.zeros((3,steps,3))
vlist=np.zeros((3,steps,3))
tlist=np.zeros(steps)
rlist[0]=rcurr
vlist[0]=vcurr
tlist[0]=tcurr

n=1

while tcurr < tmax:
    
    rnext = rcurr + vcurr*h + .5*acurr*h**2
    
    anext = a(sun,rnext,msun)
    vnext = vcurr + .5*(acurr+anext)*h
    tnext = tcurr + h
    
       
    rlist[n]=rnext
    vlist[n]=vnext
    tlist[n]=tnext
    
    rcurr,vcurr,acurr,tcurr = rnext,vnext,anext,tnext
    
    n=n+1
    
    if n % 100 == 0: print 'finished t = %0.f' % tcurr

plt.clf()
plt.plot(rlist.T[0],rlist.T[1],color='blue')
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.axes().set_aspect('equal')
plt.show()
