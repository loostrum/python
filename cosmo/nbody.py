import numpy as np
from matplotlib import pyplot as plt
from random import random
from mpl_toolkits.mplot3d.axes3d import Axes3D

# function to calculate acceleration 
def a(r):
    G = 1E-4
    epsilon = 1E-3
    acc = np.zeros((3,3))
    
    for i in range(3):
        for j in range (i+1,3):
            dist=np.linalg.norm(r[i]-r[j]) #distance between objects
            rhat = (r[i]-r[j])/dist #unit vector
            acc[i] += -G*rhat/(dist+epsilon)**2
            acc[j] -= acc[i] #particles have the same mass, so acc(i->j) = -acc(j->i)
    return acc

h=.002 #stepsize
tmax=1500 #max time to prevent infinite loop if particles don't exit box
steps=int(tmax/h)+2

rcurr= np.array([[.500,0.000,0.000],[-.250,.133,.001],[-.249,-.433,-.201]]) #array with starting positions. Shifted so box is [-1,1] instead of [0,2], which makes the boundary check easier.
vcurr= np.zeros((3,3)) #array with starting velocities (3 objects, 3 dimensions)

acurr = a(rcurr)
tcurr=0

rlist=np.zeros((3,steps,3))
rlist[:,0]=rcurr
vlist=np.zeros((3,steps,3))
vlist[:,0]=vcurr
tlist=np.zeros(steps)
tlist[0]=tcurr

#for now no random velocities for reproducability
#for i in range(2):
#    for j in range(3):
#        vcurr[i][j] = (random()-.5)*1E-4
#        vcurr[2][j] -= vcurr[i][j]


def epot(rcurr):
    epot=0
    G=1E-4
    for i in range(3):
        for j in range(3):
            if i==j:
                continue
            else:
                epot += -G/np.linalg.norm(rcurr[i]-rcurr[j])
    return epot
    
def ekin(vcurr):
    ekin=0
    for i in range(3):
        ekin += .5*np.linalg.norm(vcurr[i])**2
    return ekin
    
def etot(rcurr,vcurr):
    return epot(rcurr)+ekin(vcurr)


print 'kinetic energy in system: %0.2e \npotential energy in system: %0.2e\ntotal energy in system %0.2e' % (ekin(vcurr),epot(rcurr),etot(rcurr,vcurr))

n=1
while np.max(abs(rcurr)) <= 1 and tcurr < tmax:
    
    rnext = rcurr + vcurr*h + .5*acurr*h**2
    
    anext = a(rnext)
    vnext = vcurr + .5*(acurr+anext)*h
    tnext = tcurr + h
    
    rcurr,vcurr,acurr,tcurr = rnext,vnext,anext,tnext
    
    rlist[:,n]=rnext
    vlist[:,n]=vnext
    tlist[n]=tnext
    
    n += 1
    
    if n % 1E3 == 0: 
        print 'finished t = %0.f' % tcurr


print 'kinetic energy in system: %0.2e \npotential energy in system: %0.2e\ntotal energy in system %0.2e' % (ekin(vcurr),epot(rcurr),etot(rcurr,vcurr))



        
#xy plot
fig = plt.figure(1)
ax = fig.add_subplot(111)
for i in range(3):
    ax.plot(rlist[i].T[0][:n],rlist[i].T[1][:n])
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('XY plane')
ax.set_aspect('equal')
plt.show()


##xz plot
#fig = plt.figure(2)
#ax = fig.add_subplot(111)
#for i in range(3):
    #ax.plot(rlist[i].T[0][:n],rlist[i].T[2][:n])
#ax.set_xlim(-1,1)
#ax.set_ylim(-1,1)
#ax.set_xlabel('x')
#ax.set_ylabel('z')
#ax.set_title('XZ plane')
#ax.set_aspect('equal')
#plt.show()

##3d plot
#fig = plt.figure(3)
#ax = fig.add_subplot(111,projection='3d')
#for i in range(3):
    #ax.plot(rlist[i].T[0][:n],rlist[i].T[1][:n],rlist[i].T[2][:n])
#ax.set_xlim(-1,1)
#ax.set_ylim(-1,1)
#ax.set_zlim(-1,1)
#ax.set_aspect('equal')
#plt.show()
