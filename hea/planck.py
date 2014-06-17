import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from random import uniform
from bisect import bisect_left

rc('text', usetex=True) # use latex for greek letters with a nicer font

def planck(nu):
    T=2.936E7
    h=6.626E-27
    k=1.381E-16
    c=3E10
    #acc disk radius  = 500 r_g = 1.5E6 cm
    R=1.5E6
    #units: erg/s/cm^2/sr/Hz
    #multiply by 2 pi sr (assuming isotropic distribution, but only half is pointing upwards)
    #multiply by source area (=pi R**2) 
    #divide by h*nu to get # photons instead of their energy
    #final units: #/s/Hz
    power = 2*(math.pi*R*nu/c)**2/math.expm1(h*nu/(k*T))
    return power
    
def get_planck_photon():
    h=6.626E-27
    #generate a random photon
    nu=10**uniform(10,20)
    #get power for given nu and
    w = planck(nu)
    return nu,w
    
def main():
    nu_bins=np.logspace(10,20,100)
    weights=np.zeros(len(nu_bins))
    
    niter = int(1E5)
    
    for i in range(niter):
        #get random photon and weight
        nu,w=get_planck_photon()
        #find correct bin
        b = bisect_left(nu_bins,nu)
        #add weight to corresponding bin
        weights[b]+=w
        
        #progress tracker
        if 100*float(i+1)/niter%5==0:
            print 100*(i+1)/niter, '%'
        
    #normalise weights
    weights = [w/niter for w in weights]
    
  
    data=[nu_bins,weights]
    
    #np.savetxt('planck.txt',data)
    
    plt.loglog(nu_bins,weights)
    plt.show()
    
    
if __name__=='__main__':
    main()
    
    
#divide by h*nu to get #photons instead of energy
