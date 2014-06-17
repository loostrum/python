import math
import numpy as np
from random import uniform,randint
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
    # log nu F_nu = -9 @ log nu = 20
    # log F_nu = -29 = lognorm - p*lognu
    #lognorm = log F_nu + p*lognu = -29 + .69*20 = -15.2 = 6.3E-16 
    #wrong, 2 OOM
    norm = 5E-19
    p = .69
    nu = 10**uniform(lognu_min,lognu_max)
    #-p-1 to get # photons instead of ergs
    # final units: #/cm^2/s/Hz (?)
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
    
    print 'Initializing'
    
           
    h= 4.135668E-18 #keV / Hz
    tau=.3
    exptau=math.exp(-tau)
    n=7
    niter=int(10**n)
    
    
    electronbins,mjdist=create_MJ()
    
    n_absorbed=0
    n_scattered=0
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
                     
        if b == len(nu_bins):
            #photon falls outside bin range
            print nu
            raise ValueError('Photon energy falls outside bin range: ')
                
        else:
            #bin the photon to the non-scattered list
            n_non[b] += w_esc
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
                #scatter the photon (needs to be refined with 1+mu**2 and separate e and photon direction)
                nu *= scatter(beta,gamma)       
                n_scattered += 1
                
                #add the escaped fraction to the compt list
                #escape chance is w*exp(-tau)
                w_esc=w*exptau
                #add weight of escaped photon fraction to correct freq bin
                b=bisect_left(nu_bins,nu)
                     
                if b == len(nu_bins):
                    #photon falls outside bin range
                    print nu
                    raise ValueError('Photon energy falls outside bin range: ')
                
                else:
                    #bin the photon to the scattered list
                    n_comp[b] += w_esc
                    n_photons += 1
                    
                    
            #progress tracker        
            if 100*float(n_photons+1)/niter%5==0:
                print 100*(n_photons+1)/niter, '%'
                

    
    
        
    
    
    print 'absorbed photons: ',n_absorbed
    print 'scatterings: ',n_scattered
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
    np.savetxt('bins.txt',e_bins)
    np.savetxt('sync-input.txt',n_sync)
    np.savetxt('planck-input.txt',n_planck)
    np.savetxt('non-scattered.txt',n_non)
    np.savetxt('scattered.txt',n_comp)
    
    
    print 'Done!'
        
    #plot spectra (input, comp, full input, full)
    #plt.loglog(nu_bins,n_input,color='green',label='sync + planck input')
    plt.loglog(nu_bins,n_non,color='red',ls='--',label='non-scattered output')
    plt.loglog(nu_bins,n_comp,color='blue',ls='--',label='comptonized output')
    plt.loglog(nu_bins,n_output,color='black',label='Full output')
    plt.xlim(1E-2,1E4)
    plt.ylim(1E-35,1E-30)
    plt.xlabel(r'E (keV)')
    plt.ylabel(r'$\nu \, F_\nu$' )
    plt.legend()
    plt.show()
    
    
if __name__=="__main__":
    main()
