import math
import numpy as np
from random import uniform
from matplotlib import pyplot as plt
from scipy.stats import itemfreq
from time import time

def frange(limit1, limit2 = None, increment = 1.):
  """
  Range function that accepts floats (and integers).

  Usage:
  frange(-2, 2, 0.1)
  frange(10)
  frange(10, increment = 0.5)

  The returned value is an iterator.  Use list(frange) for a list.
  """


  if limit2 is None:
    limit2, limit1 = limit1, 0.
  else:
    limit1 = float(limit1)

  count = int(math.ceil(limit2 - limit1)/increment)
  return (limit1 + n*increment for n in range(count))



#generate random angle in lab frame
#needs to be random mu, but why?
def rand_mu():
    return uniform(-1,1)

#angle aberration for cos(theta) (to electron frame)
def mu_el(mu_l,beta):
    return (mu_l-beta)/(1-beta*mu_l)

#angle aberration for cos(theta) (to lab frame)
def mu_la(mu_e,beta):
    return (mu_e+beta)/(1+beta*mu_e)
    
#lorentz transformation for energy (relative to e_l)
def e_el(mu_l,beta,gamma):
    return gamma*(1-beta*mu_l)

#generate scattered photon energy in lab frame for given theta1_e (=scatter angle in e frame), relative to e_l
def e1_la(mu1_e,e_e,beta,gamma):
    return e_e*gamma*(1+beta*mu1_e)


def main():
    #constants (CGS)
    m=1E-27
    c=3E10
    E=1.6E-6
    gamma=E/(m*c**2)
    beta=math.sqrt(1-gamma**-2)

    print 'gamma, beta, emax/e: ',gamma,beta,(gamma*(1+beta))**2
    #niter=int(1E6)
    niter=10**int(raw_input('Give Log10(iterations): '))
    disc=0
    tstart=time()
    #make empty list for final scattered energies
    e1_list=[]
    for i in range(niter):
        #generate uniform photon direction
        mu_l=rand_mu()
        #calc incident energy in e frame
        e_e=e_el(mu_l,beta,gamma)
        #calc scattered energy in lab frame, assuming no recoil (e_e=e_1 in electron frame)
        #scattered angle is uniform in electron rest frame. 
        mu1_e=math.cos(2*math.pi*uniform(0,1))
        mu1_l=mu_la(mu1_e,beta)
        e_1=e1_la(mu1_e,e_e,beta,gamma)
        #discard if impossible 
        if e_1 < 1/(gamma**2*(1+beta)*(1-beta*mu1_l)) or e_1 > 1/(gamma*(1-beta)*(1-beta*mu1_l)):
            #print 'discarded:',e_1
            pass
        else:
            #print 'approved'
	    disc+=1
            e1_list.append(e_1)

        if 100*float(i+1)/niter%5==0:
		print 100*(i+1)/niter, '%'
    
    tend=time()
    print 'Time used for generating photons: ', round(tend-tstart,1), 's'
    
 

    ebins=list(frange(0,10,.1))
    #ebins=list(frange(-3,1,.1))
    #ebins = [ 10**item for item in ebins ]


    finallist=[]
    for i in range(len(ebins)):
        finallist.append([])
 
    tstart=time()
   
    #add items to bins
    for item in e1_list:
        for i in range(len(ebins)):
            if item > ebins[i]:
                continue
            else:
                finallist[i].append(item)
                break
              
    tend=time()
    print 'Time used for binning: ', round(tend-tstart,1), 's'
    print 'discarded: ', 100*float(disc)/niter, '%'
    ylist = [ len(elem) for elem in finallist ]
    #ylist=[ float(len(elem))/niter for elem in finallist ]
    #plt.loglog(ebins,ylist)
    plt.plot(ebins,ylist)
    plt.show()

if __name__=="__main__":
    main()
