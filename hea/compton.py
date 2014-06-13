import math
import numpy as np
from random import uniform
from matplotlib import pyplot as plt
import scipy.stats as stats
from scipy.integrate import quad
from bisect import bisect_left


#in definitions:
#_l,_e = lab,electron frame (variables)
#_la,_el = lab,electron frame (functions)
#generate random angle in lab frame
def rand_mu():
    return uniform(-1,1)

#angle aberration for mu (to electron frame)
def mu_el(mu_l,beta):
    return (mu_l-beta)/(1-beta*mu_l)

#angle aberration for mu (to lab frame)
def mu_la(mu_e,beta):
    return (mu_e+beta)/(1+beta*mu_e)
    
#lorentz transformation for energy (relative to e_l)
def e_el(mu_l,beta,gamma):
    return gamma*(1-beta*mu_l)

#generate scattered photon energy in lab frame for given mu1_e (=scatter angle in e frame), relative to e_l
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
        #this will restart the loop
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
    
        
def create_MJ():
    #create binned MJ CDF
    #electronbins = binned gamma
    electronbins= np.linspace(1,10,1E3)
    mjdist=[]
    for i in electronbins:
        mjdist.append(MJCDF(i))
        
    return electronbins,mjdist

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

 
def planck(nu):
    #T=10^9
    #h=6.626E-27
    #k=1.3806E-16
    #hkt = h / k T
    hkt=4.79924E-20
    #normalise to [0,1]
    #norm = integral over nu^2/(exp(hv/kt)-1)
    #norm=2*z3*(hkt)**-3
    lognorm=58.3374
    
    #calc the log10, to avoid overflows
    # -log(norm) + 2 log(nu) - log(expm1(h nu / k T))
    
    #when above nu ~1E22, expm1 overflows (exponent ~1000)
    #these values are not realistic anyway (P<<<1). 
    #
    if nu > 1E22:
        raise ValueError('nu is larger than 1E22')
    else:
        return -lognorm + 2*math.log10(nu) - math.log10(math.expm1(hkt*nu))
          

    
def get_seed_photon():
    #h=6.626E-27 #to convert to energy
    #generate thermal photon from planck's law
    lognu_min=15 #chance of E<1E15: 4.8E-10
    lognu_max=21 #chance of E>1E21: 0.
    nu=10**uniform(lognu_min,lognu_max)
    #weight is value of distribution at chosen nu
    #still need to divide by nu?
    w=10**planck(nu)
    
    return nu,w
    
def absorbed(nu):
    h= 4.135668E-15 #ev/Hertz
    e=h*nu #photon energy in eV
    #absorption chance depends on E and has two cases
    if e > 1:
        p_abs = .5*e**-3
        
    else:
        p_abs = .5
    
    rand=uniform(0,1)
    #if rand < p_abs, the photon will be absorbed
    if rand < p_abs:
        return True
    else:
        return False

        
def main():

    tau=.3
    exptau=math.exp(-tau)
    n=3
    niter=int(10**n)
    
    photons=[]
    nu_bins=np.logspace(10,55,100)
    weights=np.zeros(len(nu_bins))
    electronbins,mjdist=create_MJ()
    n_absorbed=0
    n_scattered=0
    n_photons=0
    for i in range(niter):
        #get seed photon
        nu,w=get_seed_photon()
        
        #keep scattering until the photon is absorbed or w < 1E-35
        while w > 1E-35:
            print w*exptau
            #escape chance is w*exp(-tau)
            w_esc=w*exptau
            #add weight of escaped photon fraction to correct freq bin
            b=bisect_left(nu_bins,nu)
            
         
            if b == len(nu_bins):
                #photon falls outside bin range
                print 'Photon energy falls outside bin range: ',nu
            else:
                #bin the photon
                weights[b]+=w_esc
                n_photons += 1
            
                
            #check if the photon is absorbed
            if absorbed(nu):
                #the photon is absorbed, so start over with a new photon
                n_absorbed += 1
                break
            else:
                #the photon will be scattered
                #the weight of the remainig photon is w(1-exptau)
                w *= (1-exptau)
                #get an electron scattering partner
                beta,gamma=get_electron(electronbins,mjdist)
                #scatter the photon (needs to be refined with 1+mu**2 and separate e and photon direction)
                nu *= scatter(beta,gamma)       
                #change the weight
                w *= exptau 
                n_scattered += 1
                #now restart from adding the escaped fraction to the spectrum
        
        #progress tracker        
        if 100*float(i+1)/niter%5==0:
            print 100*(i+1)/niter, '%'
            
  
    #normalise the weights
    total_weights=sum(weights)
    weights = [ item/total_weights for item in weights ]

    #every bin with a weight < 1E-6 has negligible influence on the spectrum
    #to improve the plots, these bins are set to zero
    for i in range(len(weights)):
        if weights[i] < 1E-6:
            weights[i]=0
    
    
    print 'absorbed photons: ',n_absorbed
    print 'scatterings: ',n_scattered
    print 'total number of photons in spectrum: ',n_photons
    
    
    
    
    #create input planck spectrum
    planck_bins = np.logspace(10,21,100)
    planck_values = []
    for nu in planck_bins:
        planck_values.append(10**planck(nu))
        
    #normalise
    nplanck = float(sum(planck_values))
    planck_values = [item/nplanck for item in planck_values]
    
    #remove items with value < 1E-6
    for i in range(len(planck_values)):
        if planck_values[i] < 1E-6:
            planck_values[i]=0
    
    
    
    
    
    #plot the planck input spectrum  
    plt.loglog(planck_bins,planck_values,color='red')
        
    plt.loglog(nu_bins,weights,color='blue')
    plt.xlim(1E14,1E25)
    plt.ylim(1E-4,1)
    plt.xlabel('nu (Hz)')
    plt.ylabel('dn/dnu')
    plt.show()
    
if __name__=="__main__":
    main()
