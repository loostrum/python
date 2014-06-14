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
    electronbins= np.linspace(1,15,1E3)
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
    #bnu = 2 h nu^3 / c^2    1/expm1(hnu/kT) (erg / sr / s / cm^2 / hz)
    #need to convert to # / cm^2 / s / hz
    #b = 8 pi nu^2 / c^2 1/expm1(hnu/kt) 
    #normalize to [0,1]
    #T=2.936E7
    #h=6.626E-27
    #k=1.3806E-16
    #hkt = h / k T
    hkt=1.63462E-18
    #norm = 1/int(nu**2/expm1(hnu/kt)) = 1/(2 zeta[3] (kt/h)^3)
    z3=1.20206
    norm=hkt**3/(2*z3)
    
    #expm1 = exp(x)-1, without losing precision
    #when above nu ~ 1E20 , expm1 overflows (exponent ~700)
    #could get around this using logs, but
    #these values are not realistic anyway (P<<<1). 
    if nu > 1E20:
        raise ValueError('nu is larger than 1E20')
    else:
        return norm*nu**2/(math.expm1(nu*hkt))
    
    
def get_seed_photon_planck():
    #generate thermal photon from planck's law
    lognu_min=12 #chance of E<1E15: 4.8E-10
    lognu_max=20 #chance of E>1E20: 0.
    nu=10**uniform(lognu_min,lognu_max)
    #weight is value of distribution at chosen nu
    #still need to divide by nu (=large particles approach)
    w=planck(nu)/nu
    return nu,w
    
    
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

#synchrotron power
def sync_power(nu):
    #from Malzac et al., 2006 
    #g_min=1.3
    #g_max=8.6
    #p=2
    #C=1/int(gamma**-p,{gamma,g_min,g_max})
    #inc angle = 45 deg
    #B = 10^6 G (del Santo et al., 2012)
    #see R&L 6.36
    #nu = omega/2pi
    #P_tot(omega) = sqrt(3) e^3 C B sin(a) / 2pi m c^2 (p+1)    gamma(p/4 + 19/12) gamma(p/4 -1/12) (mc omega / 3 e b sin(a))^-(p-1)/2
    # = const*nu^(-(p-1)/2)
    # divide by h*nu to get units of # /s / cm^3 / Hz
    #h=6.626E-27
    #B=1E6
    #C=1.53151
    #sina=.5*math.sqrt(2)
    #e=5E-10
    #p=2.
    #c=3E10
    #m=1E-27
    #const=math.sqrt(3)*e**3*C*B*sina*math.gamma(p/4+19./12)*math.gamma(p/4-1./12)*(2*math.pi*m*c/(3*e*B*sina))**(-.5)/(2*math.pi*m*c**2*3)/h
    
    #p ~ nu^-((p-1)/2) = nu^-.5
    #need to divide by hnu, to get # / s /cm^3 /Hz
    #normalize to [0,1]
    norm=5E5
    return norm*nu**-1.5
        

def get_seed_photon_sync():
    # nu_ssa=5.6E11
    # stay above this as SSA is not modelled
    lognu_min = 12
    lognu_max = 20
    nu=10**uniform(lognu_min,lognu_max)
    #weight is value of distribution at chosen nu
    #still need to divide by nu (=large particles approach)
    w=sync_power(nu)/nu
    return nu,w
        
def fix(x_list,y_list,width):
    h= 1#4.135668E-18 #keV / Hz
    #normalise the binned values
    total_values=sum(y_list)
    values = [ item/total_values for item in y_list ]
    
    
    #every bin with a value < 1E-6 has negligible influence on the spectrum
    #to improve the plots, these bins are set to zero
    #in the E F(E) plot, there now is a sudden break
    #when plotting everything, the uncertainty in dn/de can be seen (large wiggles)
    for i in range(len(values)):
        if values[i] < 1E-6:
            values[i]=0
    
    
    #divide the values by the nu_bin width (log) to get dn/de
    for i in range(len(values)):
        values[i] /= width
        
    #convert xlist to energy instead of freq:
    x_list = [ h*i for i in x_list ]

    #multiply each weight by h*nu to get E F(E)
    #needs to be done after normalising!
    for i in range(len(values)):
        values[i] *= x_list[i]
        
    return x_list,values
    

def main():
    
    print 'Initializing'
    
    h= 4.135668E-18 #keV / Hz
    tau=.3
    exptau=math.exp(-tau)
    n=5
    niter=int(10**n)
    
    photons=[]
    
    step=100
    mi=10
    ma=35
    width=float(ma-mi)/(step-1) #width of a bin
    nu_bins=np.logspace(mi,ma,step)
    weights=np.zeros(len(nu_bins))
    electronbins,mjdist=create_MJ()
    n_absorbed=0
    n_scattered=0
    n_photons=0
    
    pbins=nu_bins
    planckphotons=np.zeros(len(pbins))
    
    print 'Generating photons'
    while n_photons < niter:
        #get seed photon
        #50/50 synch or planck
        #if uniform(0,1) < .5:
        #get planck photon
        nu,w=get_seed_photon_planck()
        w0=w
        #add photon to planck spectrum
        b=bisect_left(pbins,nu)
        if b==len(pbins):
            print nu
            raise ValueError('Photon falls outside bin range: ')
        else:
            planckphotons[b]+=w
        #else:
            #get sync photon
            
            
        #keep scattering until the photon is absorbed or w/w0 < 1E-6
        while w/w0 > 1E-6:
            
            #escape chance is w*exp(-tau)
            w_esc=w*exptau
            #add weight of escaped photon fraction to correct freq bin
            b=bisect_left(nu_bins,nu)
                     
            if b == len(nu_bins):
                #photon falls outside bin range
                print nu
                raise ValueError('Photon energy falls outside bin range: ')
                
            else:
                #bin the photon
                weights[b]+=w_esc
                n_photons += 1
                
                #progress tracker        
                if 100*float(n_photons+1)/niter%5==0:
                    print 100*(n_photons+1)/niter, '%'
            
                
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
                #now goto adding the escaped fraction to the spectrum

    
    
        
    
    
    print 'absorbed photons: ',n_absorbed
    print 'scatterings: ',n_scattered
    print 'total number of photons in spectrum: ',n_photons
    
    print 'Manipulating data'
    
    e_bins,weights = fix(nu_bins,weights,width)
    pbins,planckphotons = fix(pbins,planckphotons,width)
    
    
       
    #plt.loglog(nu_bins,weights,color='blue')
    plt.loglog(pbins,planckphotons,color='red')
    plt.loglog(e_bins,weights)
    #plt.xlim(1E14,1E22)
    #plt.ylim(1E14,1E18)
    plt.xlabel('E (keV)')
    plt.ylabel('E dn/dE')
    plt.show()
    
if __name__=="__main__":
    main()
