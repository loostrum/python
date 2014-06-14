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
    

def scatter(beta,gamma):
    while True:
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
            continue
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
    #take that element from electronbins (e.g. electronbins[8]).
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
    electronbins= np.linspace(1,10,1E3)
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
    #ntot=16*math.pi*z3*(k*T/(h*c))**3
    ntot=2.0284E28
    return 8*math.pi*nu**2/(ntot*c**3*math.expm1(h*nu/(k*T)))
    #uses expm1 == exp(x)-1, which avoids loss of precision
    
def planckCDF(nu):
    return quad(planck,0,nu)[0]
    #max 1E22 Hz, to avoid overflow error
    
    
def create_planck():
    photonbins= np.logspace(0,21,1000)
    planckdist=[]
    for i in photonbins:
        planckdist.append(planckCDF(i))
        
    return photonbins,planckdist
    
def get_seed_photon(photonbins,planckdist):
    #h=6.626E-27 #to convert to energy
    #generate thermal photon from planck's law
    while True:
        b = bisect_left(planckdist,uniform(0,1))
        if b == len(planckdist):
            #photon falls outside bin range
            continue
        else:
            nu = photonbins[b]
            #return photon freq
            return nu
    
def escape():
    rand=uniform(0,1)
    tau=.3
    escape_chance=math.exp(-tau)
    if rand <= escape_chance:
        return True
    else:
        return False
        
def main():

    #create binned distributions
    electronbins,mjdist = create_MJ()     
    photonbins,planckdist = create_planck()
    
    
    n=int(raw_input('Number of iterations (Log): '))
    niter=10**n
    
    #init list for final photon frequencies
    fbins=np.logspace(15,25,100)

    print 'Producing 1E'+str(n)+' photons'
    tstart=time()



    #initialize stuff
    disc=0 #=amount of discarded photons
    finallist=np.zeros(len(fbins)) #list that will contain escaped photon frequencies.
    
    #generate niter photons and scatter them 
    for i in range(niter):
        #get photon
        photon_freq=get_seed_photon(photonbins,planckdist)
        #check if it escapes
        while not escape():
            #choose a random electron
            beta,gamma=get_electron(electronbins,mjdist)
            #tmp is e/e0=f/f0 for the photon. Need to multiply by random energy from synchrotron/thermal distribution to get e_final.
            # for now only thermal photon.
            #multiply photon freq with scattering
            photon_freq *= scatter(beta,gamma)
            #now repeat until the photon escapes
            
            
        #the photon has escaped, so it needs to be binned        
        #find bin to put photon in.
        b = bisect_left(fbins,photon_freq)
        #fix if freq is higher than allowed in bins
        if b == len(fbins):
            print 'f larger than maximum bin: ',photon_freq
            disc+=1
        else:
            #add photon to correct bin
            finallist[b]+=1
            
        if 100*float(i+1)/niter%5==0:
            print 100*(i+1)/niter, '%'
            
    print 'discarded: ', 100*float(disc)/niter, '%'
    

    #normalise results
    nphotons=float(niter-disc)
    finallist = [item/nphotons for item in finallist]
    tend=time()
    print 'Time used for generating photons: ', round(tend-tstart,1), 's'
    
    
    
    
   
    plt.clf()
    plt.plot(fbins,finallist)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('f (Hz)')
    plt.ylabel('dn/df')
    #plt.xlim(1E18,1E22)
    #plt.ylim(1E-4,1E0)
    
    plt.show()
    
if __name__=="__main__":
    main()
