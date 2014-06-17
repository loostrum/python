import math
import numpy as np
from random import uniform,randint
from matplotlib import pyplot as plt
import scipy.stats as stats
from scipy.integrate import quad
from bisect import bisect_left
from time import time


def scatter(beta,gamma):
    while True:
        rho_min = (1.0-beta)**2/(2*beta)
        rho_max = (1.0+beta)**2/(2*beta)
        rho = uniform(rho_min,rho_max)
        muth = (1-math.sqrt(2*beta*rho))/beta
        muth_p = (muth-beta)/(1-beta*muth)    
        
        #photon energy in electron frame
        e1_p = gamma*(1.0-beta*muth)
        
        a = uniform(-1,1)
        b = uniform(0,2)
        while b >= 1.0+a**2:
            a = uniform(-1,1) #a is distributed like a cosine
            b = uniform(0,2)
        mua_p = a
        
        muphi_p = math.cos(uniform(0,1)*2.0*math.pi)		# Uniform distribution for muphi_p
        sin_theta_p = math.sqrt(1.0-muth_p**2)		# outgoing theta prime
        sin_alpha_p = math.sqrt(1.0-mua_p**2)		# outgoing alpha prime
        muth1_p = muth_p*mua_p-sin_alpha_p*sin_theta_p*muphi_p		# scattering angle in electron rest frame: muth1_p = cos_theta1_p
    
        e1 = e1_p*gamma*(1.0+beta*muth1_p)		# Final photon energy in lab frame
    
        return e1

        
        
def MJCDF(gamma):
    #returns the cumulative probability of finding energy gamma in a MJ distribution
    #T = 10^9 K
    theta=.168637 # =kt/mc^2
    const = 3238.1 # = 1/(theta * K_2(1/theta))
    return quad(lambda gamma: const*math.sqrt(1-gamma**-2)*gamma**2*math.exp(-gamma/theta),1,gamma)[0]

        
def create_MJ():
    #create binned MJ CDF
    #electronbins = binned gamma
    electronbins= np.linspace(1,15,1E3)
    mjdist=[]
    for i in electronbins:
        mjdist.append(MJCDF(i))
        
    return electronbins,mjdist

def get_sync_photon():
    #minimum freq: 1.3E13, maximum freq: 8E19
    #use 1E13 and 1E20
    lognu_min = 13
    lognu_max = 20
    #The multi-wavelength polarization of Cygnus X-1
    #Russel & Shahbaz, 2013
    #norm follows from fig 2.
    norm = 5E-19
    p = .69
    nu = 10**uniform(lognu_min,lognu_max)
    #-p-1 to get # photons instead of ergs
    # final units: #/cm^2/s/Hz
    w = norm*nu**(-p-1)
    
    return nu,w

def planck(nu):
    T=2.936E7
    h=6.626E-27
    k=1.381E-16
    norm=1E-82
    #units: erg/s/cm^2/sr/Hz
    #divide by h*nu to get # photons instead of their energy
    #final units: #/cm^2/s/Hz (sr are fixed in norm)
    power = nu**2/math.expm1(h*nu/(k*T))
    #planck peak is about 2 oom higher than sync peak
    power *= norm
    return power
    
def get_planck_photon():
    #generate a random photon in sensibile range
    nu=10**uniform(10,20)
    #get power for given nu
    w = planck(nu)
    return nu,w


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

    
   
    
def absorbed(nu):
    h= 4.135668E-15 #ev/Hertz
    e=h*nu #photon energy in eV
    #absorption chance depends on E and has two cases
    if e > 1:
        p_abs = .5*e**-3
        
    else:
        p_abs = .5
    
    #if rand < p_abs, the photon will be absorbed
    if uniform(0,1) < p_abs:
        return True
    else:
        return False

 
def main():
    
    tstart=time()
    
    print 'Initializing'
    
           
    h= 4.135668E-18 #keV / Hz
    tau=.3
    exptau=math.exp(-tau)
    #txt files made with = 8, takes .5 hrs
    n=6
    niter=int(10**n)
    
    
    electronbins,mjdist=create_MJ()
    
    n_absorbed=0
    n_photons=0
    
    nu_bins=np.logspace(10,30,200)
    n_planck=np.zeros(len(nu_bins))
    n_sync=np.zeros(len(nu_bins))
    n_non=np.zeros(len(nu_bins))
    n_comp=np.zeros(len(nu_bins))
    
    
    
    print 'Generating photons'
    while n_photons < niter:
        #get seed photon
        #50/50 sync or planck
        if uniform(0,1) < .5:
            #get planck photon
            nu,w = get_planck_photon()
            w0 = w
            #bin the planck photon for reference
            b = bisect_left(nu_bins,nu)
            n_planck[b] += w

        else:
            #get sync photon
            nu,w = get_sync_photon()
            w0 = w
            #bin the sync photon for reference
            b = bisect_left(nu_bins,nu)
            n_sync[b] += w
            
            
                        
        #escape chance is w*exp(-tau)
        w_esc=w*exptau
        #add weight of escaped photon fraction to correct freq bin
        b=bisect_left(nu_bins,nu)
        #bin the photon to the non-scattered list
        n_non[b] += w*exptau
        n_photons += 1
                                
                
        #progress tracker        
        if 100*float(n_photons+1)/niter%5==0:
            print 100*(n_photons+1)/niter, '%'
                    
                    
            
        #keep scattering until the photon is absorbed or w/w0 < 1E-6
        while w/w0 > 1E-6:
            
                
            #check if the photon is absorbed
            if absorbed(nu):
                #the photon is absorbed, so start over with a new photon
                n_absorbed += 1
                break
            else:
                #the photon will be scattered
                #the weight needs to be changed, as only 1-exptau part of the photon still exists
                w *= (1-exptau)
                #get an electron scattering partner
                beta,gamma=get_electron(electronbins,mjdist)
                #scatter the photon
                nu *= scatter(beta,gamma)       
                
                #add the escaped fraction to the compt list
                #escape chance is w*exp(-tau)
                w_esc=w*exptau
                #add weight of escaped photon fraction to correct freq bin
                b=bisect_left(nu_bins,nu)
                
                #bin the photon to the scattered list
                n_comp[b] += w_esc
                n_photons += 1
                    
                    
            #progress tracker        
            if 100*float(n_photons+1)/niter%5==0:
                print 100*(n_photons+1)/niter, '%'
                

    
    
        
    
    
    print 'absorbed photons: ',n_absorbed
    print 'total number of photons in spectrum: ',n_photons
    
    print 'Manipulating data'
    
    #transform weights to nu * fnu and nu_bins to energy
    #also normalise to niter
    e_bins = nu_bins
    for i in range(len(nu_bins)):
        n_planck[i] *= nu_bins[i]/niter
        n_sync[i] *= nu_bins[i]/niter
        n_non[i] *= nu_bins[i]/niter
        n_comp[i] *= nu_bins[i]/niter
        e_bins[i] *= h
    
    n_input = np.add(n_planck,n_sync)
    n_output = np.add(n_non,n_comp)
    
    #save data
    #d='data'
    #np.savetxt(d+'/bins.txt',e_bins)
    #np.savetxt(d+'/sync-input.txt',n_sync)
    #np.savetxt(d+'/planck-input.txt',n_planck)
    #np.savetxt(d+'/non-scattered.txt',n_non)
    #np.savetxt(d+'/scattered.txt',n_comp)
    
    tend=time()
    print 'Done!'
    print 'Total time used: ',tend-tstart
        
    #plot spectra (input, comp, full input, full)
    plt.loglog(nu_bins,n_input,color='green',label='Sync + planck input')
    plt.loglog(nu_bins,n_non,color='red',ls='--',label='Non-scattered output')
    plt.loglog(nu_bins,n_comp,color='blue',ls='--',label='Comptonized output')
    plt.loglog(nu_bins,n_output,color='black',label='Full output')
    plt.xlim(1E-2,1E4)
    plt.ylim(1E-35,1E-30)
    plt.xlabel(r'E (keV)')
    plt.ylabel(r'$\nu \, F_\nu$' )
    plt.legend()
    plt.show()
    
    
if __name__=="__main__":
    main()
