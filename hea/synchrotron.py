import math
import numpy as np
from matplotlib import pyplot as plt
from random import uniform
from scipy.special import kv
from scipy.integrate import quad
from bisect import bisect_left

def pow_spec_rand():
    #input electron power spectrum
    g_min=1.3
    g_max=8.6
    g_bins=np.linspace(1.2,8.7,100)
    p=2
    norm=(1-p)/(g_max**(1-p)-g_min**(1-p))
    #make binned CDF (CDF = integral(gmin,g) of g*(-p)
    pow_bins=[norm*(g**(1-p)-g_min**(1-p))/(1-p) for g in g_bins]
    rand_pow = uniform(0,1)
    b = bisect_left(pow_bins,rand_pow)
    return g_bins[b]    
    
    
    
def sync_power(nu,gamma):
    #power spectrum from single electron
    #relativistic case (eq 6.18)
    B = 1E6
    e = 5E-10
    m = 1E-27
    c = 3E10
    sina = .5*math.sqrt(2)
    nu_c= 3*gamma**2*e*B*sina/(4*math.pi*m*c)
    x = nu/nu_c
    F = x*quad(lambda x: kv(5./3,x),x,np.inf)[0]
    #plotting F returns fig on p. 179, despite the convergence warning
    power = math.sqrt(3)*e**3*B*sina/(m*c**2)*F
    return power
    
def get_sync_photon():
    #generate a random photon
    nu=10**uniform(10,20)
    #get random gamma 
    gamma=pow_spec_rand()
    #get power(=weight) for given nu and gamma
    w = sync_power(nu,gamma)
    return nu,w
    
    
    
def main():
    nu_bins=np.logspace(10,20,100)
    weights=np.zeros(len(nu_bins))
    
    niter=int(1E5)
    
    for i in range(niter):
        #get random photon and weight
        nu,w=get_sync_photon()
        #find correct bin
        b = bisect_left(nu_bins,nu)
        #add weight to corresponding bin
        weights[b]+=w
        
        #progress tracker
        if 100*float(i+1)/niter%5==0:
            print 100*(i+1)/niter, '%'
        
    #normalise weights
    weights = [w/niter for w in weights]
    
    #ignore bins <1E-6*max contribution (negligible contribution to spectrum)
    wmax=max(weights)
    for i in range(len(weights)):
        if weights[i] < wmax*1E-6:
            weights[i] = 0.
    
    data=[nu_bins,weights]
    
    #np.savetxt('sync.txt',data)
    
    
    
    #plot
    plt.loglog(nu_bins,weights)
    plt.show()
    

if __name__=='__main__':
    main()
