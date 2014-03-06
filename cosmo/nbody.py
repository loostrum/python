#use leapfrog:
#xii = xi + vi dt + .5 ai dt^2
#vii = vi + .5 (ai + aii) dt
#ai = - G / xi^2
# note that in order to calculate vii, both ai and aii need to be known. aii depends on xii, so xii has to be calculated first.

import numpy as np
from matplotlib import pyplot as plt


# test with earth around sun
au=1.5E11
msun=1.9E30 #kg
sun=np.array([0,0,0]) #xyz
earth=np.array([au,0,0])

def a(loc1,loc2,m1):
    G = 6.67E-11
    acc=np.array([])
    for i in range(3):
        if loc1[i] == loc2[i]:
            acc=np.append(acc,0)
        else:    
            acc=np.append(acc,-G*m1/(loc2[i]-loc1[i])**2)
    return acc
    


h=3600*24 #1 day
tmax=365*3600*24 #1 year
steps=366
rcurr=earth
vcurr=np.array([0,29780,0])
tcurr=0
acurr=a(sun,rcurr,msun)
print acurr
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
    print anext
    rlist[n]=rnext
    vlist[n]=vnext
    tlist[n]=tnext
    
    rcurr,vcurr,acurr,tcurr = rnext,vnext,anext,tnext
    
    n=n+1

plt.clf()
plt.plot(rlist.T[0],rlist.T[1],color='blue')
plt.xlim(-1.1*au,1.1*au)
plt.ylim(-1.1*au,1.1*au)
plt.show()
