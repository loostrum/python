from random import uniform
import math
from matplotlib import pyplot as plt
from scipy.stats import itemfreq

def rand_mu():
	return uniform(-1,1)

def rand_theta():
	return math.cos(uniform(0,1)*2*math.pi)

n=int(1E5)

mlist=[]
thlist=[]
for i in range(n):
    mlist.append(rand_mu())
    thlist.append(rand_theta())

mlist = itemfreq([ round(item,1) for item in mlist ])
thlist = itemfreq([ round(item,1) for item in thlist ])


plt.plot(mlist[:,0],mlist[:,1])
plt.plot(thlist[:,0],thlist[:,1])
plt.show()


