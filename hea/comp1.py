import math
import numpy as np
from random import random

#generate random angle in lab frame
def rand_mu():
    return math.cos(2*math.pi*random())

#angle aberration for cos(theta)
def mu_el(mu_l,beta):
    return (mu_l-beta)/(1-beta*mu_l)
    
#lorentz transformation for energy (relative to e_l)
def e_el(mu_l,beta,gamma):
    return gamma*(1-beta*mu_l)

#generate scattered photon energy in lab frame for given theta1_e (=scatter angle in e frame), relative to e_l
def e1_l(theta1_e,e_e,beta,gamma):
    return e_e*gamma*(1+beta*math.cos(theta1_e))


def main():
    #constants (CGS)
    m=1E-27
    c=3E10
    E=1.6E-6
    gamma=E/(m*c**2)
    beta=math.sqrt(1-gamma**-2)
    niter=100
    
   
    #create scattered angle list
    thlist=np.arange(0,2*math.pi,1E-3)
    #creat temp list to save scattered energies
    tmp=np.zeros(len(thlist))
    
    #make empty list for final scattered energies
    e1_list=[]
    
    for i in range(niter):
        #generate random photon direction
        mu_l=rand_mu()
        #calc incident energy in e frame
        e_e=e_el(mu_l,beta,gamma)
        #generate possible scatter energies for different theta
        for j in range(len(thlist)):
            tmp[j]=e1_l(thlist[j],e_e,beta,gamma)
        e1_list.append(tmp)
    #flatten e1_list
    e1_list=[item for sublist in e1_list for item in sublist]
    print len(e1_list)
   
    
    
    
if __name__=="__main__":
    main()
