#use leapfrog:
#xii = xi + vi dt + .5 ai dt^2
#vii = vi + .5 (ai + aii) dt
#ai = - G / xi^2
# note that in order to calculate vii, both ai and aii need to be known. aii depends on xii, so xii has to be calculated first.

import numpy as np
from matplotlib import pyplot as plt


# test with earth around sun
au=1.5E11 #m
msun=1.9E30 #kg
G = 6.67384E-11
sun=np.array([0.,0.,0.]) #xyz
earth=np.array([au,0.,0.])

   
def a(loc1,loc2,m1):
    G = 6.67384E-11
    acc=-G*m1*loc2/(np.linalg.norm(loc2-loc1))**3
    return acc

h=3600 # seconds
tmax=8760*3600 #1 year
steps=int(tmax/h)+1
rcurr=earth
vcurr=np.array([0,np.sqrt(G*msun/au),0])
tcurr=0
acurr=a(sun,rcurr,msun)
rlist=np.zeros((steps,3))
vlist=np.zeros((steps,3))
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
plt.xlim(-1.1*au,1.1*au)
plt.ylim(-1.1*au,1.1*au)
plt.axes().set_aspect('equal')
plt.show()
