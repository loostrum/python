import math
import numpy as np
from random import uniform
from matplotlib import pyplot as plt
import scipy.stats as stats
from scipy.integrate import quad
from time import time
from bisect import bisect_left
#import argparse

# arguments that can be given to the program
#argparser=argparse.ArgumentParser(description='Calculates an inverse Compton spectrum')
#argparser.add_argument('-s','--save', dest='save',action='store_true', help='Save figures instead of displaying them')
#argparser.add_argument('-v','--verbose', dest='v',action='store_true', help='Switch on verbosity')
#argparser.add_argument('N',type=int ,help='Number of iterations (Log)')
#args=argparser.parse_args()



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
    

def gen_photon(beta,gamma):
    #generate random photon direction
    mu_l=rand_mu()
    #calc incident energy in e frame
    e_e=e_el(mu_l,beta,gamma)
    #calc scattered energy in lab frame, assuming no recoil (e_e=e_1 in electron frame)
    mu1_e=rand_mu()
    mu1_l=mu_la(mu1_e,beta)
    e_1=e1_la(mu1_e,e_e,beta,gamma)
    #discard if impossible 
    if e_1 < 1/(gamma**2*(1+beta)*(1-beta*mu1_l)) or e_1 > 1/(gamma*(1-beta)*(1-beta*mu1_l)):
        #print 'discarded:',e_1
        return 0
    else:
        #print 'approved'
        return e_1
        
def MJ(gamma):
    #returns the probability of finding energy gamma in a MJ distribution
    #T = 10^9 K
    theta=.168637 # =kt/mc^2
    const = 3238.1 # = 1/(theta * K_2(1/theta))
    return const*math.sqrt(1-gamma**-2)*gamma**2*math.exp(-gamma/theta)

def MJCDF(gamma):
    # returns CDF of the MJ PDF
    return quad(MJ,1,gamma)[0]
    
        
def get_electron(electronbins,mjdist):
    #steps to select gamma_electron:
    #select random number between 0 and 1
    #find in which bin this belongs using mjdist (e.g. index 8)
    #take that element from electronbins (e.g. electronbins[[8]).
    #the result is the gammafactor of the elecron.
    b = bisect_left(mjdist,uniform(0,1))
    if b == len(mjdist):
        b-=1
    gamma = electronbins[b]
    beta = math.sqrt(1-gamma**-2)
    return beta,gamma

def create_MJ():
    #create binned MJ CDF
    #electronbins = binned gamma
    electronbins= list(frange(1,5,1E-3))
    mjdist=[]
    for i in electronbins:
        mjdist.append(MJCDF(i))
        
    return electronbins,mjdist
    
def planck(nu):
    #T = 10^9 K, kT = 1.3806E-7 erg
    #create planck spectrum, converted to a PDF
    #which is B_nu(T)*4pi(isotropic)/c(to energy density)/h nu (to number density) / n_tot (to PDF) 
    k=1.3806E-16
    T=1E9
    h=6.626E-27 # erg*s
    c=2.998E10 #cm/s
    #Zeta(3)
    z3=1.20206
    #ntot is integral over modified B_nu
    ntot=16*math.pi*z3*(k*T/(h*c))**3
    return 8*math.pi*nu**2/(ntot*c**3*math.expm1(h*nu/(k*T)))
    #uses expm1 == exp(x)-1, which avoids loss of precision
    
def planckCDF(nu):
    return quad(planck,0,nu)[0]
    #max 1E22 Hz, to avoid overflow error
    
    
def create_planck():
    photonbins= [10**i for i in frange(1,20,.1)]
    planckdist=[]
    for i in photonbins:
        planckdist.append(planckCDF(i))
        
    return photonbins,planckdist
    
def get_seed_photon():
    #generate thermal photon from planck's law
    pass


def main():

    electronbins,mjdist = create_MJ()     
    photonbins,planckdist = create_planck()
 
    
    n=int(raw_input('Number of iterations (Log): '))
    niter=10**n


    print 'Producing 1E'+str(n)+' photons'

            
    #initialize stuff
    disc=0
    e1_list=[]
    ebins=np.asarray(list(frange(0,10,.1)))
    finallist=np.zeros(len(ebins))

    
    tstart=time()
    
    #generate photons, only e/e0 for a given gamma,beta
    for i in range(niter):
        #choose a random electron
        beta,gamma=get_electron(electronbins,mjdist)
        #tmp is e/e0 for the photon. Need to multiply by random energy from synchrotron/thermal distribution to get e_final.
        tmp=gen_photon(beta,gamma)
        #gen_photon returns 0 if the angle/energy combination was not allowed
        if tmp == 0:
            disc+=1
        else:
            #e1_list.append(tmp)
            b = bisect_left(ebins,tmp)
            #fix if it belongs in the last bin (gives max index + 1)
            if b == len(ebins):
                b-=1
            finallist[b]+=1
        if 100*float(i+1)/niter%5==0:
            print 100*(i+1)/niter, '%'
    
    tend=time()
    
    print 'Time used for generating photons: ', round(tend-tstart,1), 's'
    
 
    print 'discarded: ', 100*float(disc)/niter, '%'
    

    
   
    plt.clf()
    plt.plot(ebins,finallist)
    plt.xlabel('e1/e0')
    plt.ylabel('dn/de')
    plt.show()
    
if __name__=="__main__":
    main()
