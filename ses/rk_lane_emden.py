from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
from math import pi
import argparse

argparser=argparse.ArgumentParser(description='Solves the Lane-Emden equation')
argparser.add_argument('-s','--save', dest='save',action='store_true', help='Save figures instead of displaying them')
argparser.add_argument('-v','--verbose', dest='v',action='store_true', help='Switch on verbosity')
args=argparser.parse_args()

rc('text', usetex=True) #use latex for greek letters with a nicer font

# define Cash-Karp parameters
a2,   a3,  a4,  a5,  a6      =        1/5.,    3/10.,       3/5.,            1.,        7/8.
b21, b31, b32, b41, b42, b43 =        1/5.,    3/40.,      9/40.,         3/10.,      -9/10., 6/5.
b51, b52, b53, b54           =     -11/54.,     5/2.,    -70/27.,        35/27.
b61, b62, b63, b64, b65      = 1631/55296., 175/512., 575/13824., 44275/110592.,   253/4096. 
c1,   c2,  c3,  c4,  c5, c6  =     37/378.,       0.,   250/621.,      125/594.,          0.,  512/1771.
c1star, c2star, c3star, c4star, c5star, c6star = 2825/27648., 0.,  18575/48384.,13525/55296., 277/14336., 1/4.

def stepper(derivx,derivy, n, t, x, y, h, tol): # we have functions where x' is not a function of x, but of y and t so we have y'(x,t) and x'(y,t)

   k1x = h*derivx(n,t,y)
   k1y = h*derivy(n,t,x)
   k2x = h*derivx(n,t+a2*h,y+b21*k1y)
   k2y = h*derivy(n,t+a2*h,x+b21*k1x)
   k3x = h*derivx(n,t+a3*h,y+b31*k1y+b32*k2y)
   k3y = h*derivy(n,t+a3*h,x+b31*k1x+b32*k2x)
   k4x = h*derivx(n,t+a4*h,y+b41*k1y+b42*k2y+b43*k3y)
   k4y = h*derivy(n,t+a4*h,x+b41*k1x+b42*k2x+b43*k3x)
   k5x = h*derivx(n,t+a5*h,y+b51*k1y+b52*k2y+b53*k3y+b54*k4y)
   k5y = h*derivy(n,t+a5*h,x+b51*k1x+b52*k2x+b53*k3x+b54*k4x)
   k6x = h*derivx(n,t+a6*h,y+b61*k1y+b62*k2y+b63*k3y+b64*k4y+b65*k5y)
   k6y = h*derivy(n,t+a6*h,x+b61*k1x+b62*k2x+b63*k3x+b64*k4x+b65*k5x)
   x_n_plus_1      = x +     c1*k1x +     c2*k2x +     c3*k3x +     c4*k4x +     c5*k5x +     c6*k6x
   y_n_plus_1      = y +     c1*k1y +     c2*k2y +     c3*k3y +     c4*k4y +     c5*k5y +     c6*k6y
   x_n_plus_1_star = x + c1star*k1x + c2star*k2x + c3star*k3x + c4star*k4x + c5star*k5x + c6star*k6x
   y_n_plus_1_star = y + c1star*k1y + c2star*k2y + c3star*k3y + c4star*k4y + c5star*k5y + c6star*k6y
   DELTAx          = x_n_plus_1 - x_n_plus_1_star
   DELTAy          = y_n_plus_1 - y_n_plus_1_star           
   try:
       hx = h*abs(tol/DELTAx)**0.2    # Finds step size required to meet given tolerance
   except ZeroDivisionError:
       hx = h                        # When you are very close to ideal step, DELTA can be zero
       
   try:
       hy = h*abs(tol/DELTAy)**0.2
   except ZeroDivisionError:
       hy = h
   
   h1 = min(hx,hy)
   return x_n_plus_1,y_n_plus_1, h1


def control(derivx,derivy,n,t0,x0,y0,t_max,h,tol):

        #set initial solutions lists
        t = np.array(t0)
        x = np.array(x0)
        y = np.array(y0)
       
        t_curr,x_curr,y_curr = t0,x0,y0
        count_iter,count_total = 0,0
        
        adaptive = True #enables change in stepsize if the ideal h changes by more than 10%
     
    
        while x_curr > 1E-15: # stop if theta is close enough to zero
            if t_curr > t_max:  #failsafe if theta doesn't go to zero
                break
            x_next,y_next,h1 = stepper(derivx,derivy,n,t_curr,x_curr,y_curr,h,tol)
            t_next=t_curr+h
            
            if h1 < 0.9*h and adaptive:
                if args.v: print 'Decreased stepsize from %0.2e to %0.2e at theta = %0.5f' % (h,h1,x_curr) 
                h = h1
            elif h1 > 1.1*h and adaptive:
                if args.v: print 'Increased stepsize from %0.2e to %0.2e at theta = %0.5f' % (h,h1,x_curr)
                h = h1
            elif x_next < 0:    # if h is ok for precision, check if it doesn't get theta below zero
                if args.v: print 'Theta is getting below zero, adapting stepsize. Theta: %0.4e' % x_curr
                h=.01*h
                adaptive=False #disables optimal h check as this could cause an endless loop
            else:            
                t = np.append(t,t_next)
                x = np.append(x,x_next)
                y = np.append(y,y_next)
                t_curr,x_curr,y_curr = t_next,x_next,y_next
                count_iter += 1
            
                
            count_total += 1
        
        if args.v: print 'Calculation completed with %0.1f iterations in %0.1f steps for n= %0.2f' % (count_iter,count_total,n)
        return t,x,y


def dthetadxi(n,xi,phi):
    try: 
        return phi/xi**2.
    except ZeroDivisionError:
        # the boundary conditions state theta'(0)=0
        return 0.
        
def dphidxi(n,xi,theta):
    try:
        return -xi**2.*theta**n
    except ValueError:
        # fails when theta<0 and n!=int. Then the surface of the star is reached anyway, so we can safely return zero. This should not be reached as the integration stops at positive theta
        return 0.
        
def analyt(n,xi):
    if n==0:
        return 1.-(xi**2.)/6.
    elif n==1:
        if xi==0.:
            return 1.
        else:    
            return np.sin(xi)/xi
        
    elif n==5:
        return (1.+(xi**2.)/3.)**-.5
    else:
        raise ValueError('This value of n has no analytical solution: '+str(n))    
        
vanalyt=np.vectorize(analyt) #vectorize so an numpy array can be used as input. The if statement in n=1 doesn't work in an array.

def plotter(n,x,y,ylabel):
    plt.clf()
    plt.plot(x,y,color='blue',label='Numerical')
    plt.title('n = %0.2f' %(n))
    plt.xlabel(r'$\xi$')
    plt.ylabel(ylabel)
    plt.ylim(0,1)
    plt.legend()
    plt.show()

def plotall(nlist,x,y,ylabel,xlim,ymin=0,ymax=1):
    plt.clf()
    for i in range(len(nlist)):
        plt.plot(x[i],y[i],label='n = '+str(nlist[i]))
    plt.xlim(0,xlim)
    plt.ylim(ymin,ymax)
    plt.title('All numerical solutions for '+ylabel)
    plt.xlabel(r'$\xi$')
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()
        
def solver(nlist):
    xi0=0.
    ximax=20.
    phi0=0.
    theta0=1.
    h=1E-4
    tol=1E-15
    
    xi,theta,phi,rho = [],[],[],[]
    
    for n in nlist:
        x,y,z = control(dthetadxi,dphidxi,n,xi0,theta0,phi0,ximax,h,tol)
        xi.append(x)
        theta.append(y)
        phi.append(z)
        rho.append(y**n)
        
        
    return xi,theta,phi,rho

# equation for rho_c as given in eq. 8.28
def D_n(n,xi_0,phi_0):
    return 1/(-3*dthetadxi(n,xi_0,phi_0)/xi_0)

def rho_core(n,xi_0,phi_0,M,R):
    average=3*M/(4*pi*R**3)
    central=average*D_n(n,xi_0,phi_0)
    return central

# solve neutron star rho_c
def solve_ns(nlist,zeros,phi_at_zero):    
    mass = 1.4*1.988435E33 #gram
    radius = 1E6 #centimeter
    rho_nuc = 4E14 #gram/cm3
    n = 1
    xi_0=zeros[nlist.index(n)]
    phi_0 = phi_at_zero[nlist.index(n)]
    rho = rho_core(n,xi_0,phi_0,mass,radius)
    return rho/rho_nuc
    
# solve high rho WD mass  
def solve_wd(nlist,zeros,phi_at_zero):
    n = 3
    xi_0 = zeros[nlist.index(n)]
    phi_0 = phi_at_zero[nlist.index(n)]
    
    #mass-radius relation for n=3 gives GM/m_3 = 16K**3/pi*G where m_3 = -xi_0**2 dthetadxi = - phi_at_zero
    K=4.88735E14 #cgs
    G=6.67384E-8 #cgs
    solarmass=1.988435E33
    M=-phi_0*(16*K**3/(pi*G**3))**.5
    return M/solarmass
    
def print_table(nlist,zeros,phi_at_zero):
    print '\nn      xi            |phi|'
    for i in range(len(nlist)):
        if nlist[i]==4:
            print '%0.1f    %0.8f   %0.8f' % (nlist[i],zeros[i],-phi_at_zero[i])
        else:
            print '%0.1f    %0.8f    %0.8f' % (nlist[i],zeros[i],-phi_at_zero[i])



if __name__=='__main__':
    
    nlist = [0,1,1.5,2,3,4] #n for which to solve the l-e equation
    xi,theta,phi,rho = solver(nlist) #solutions for n in nlist
    
    #useful to have lists with the values for phi and xi at theta=0:    
    zeros=[]
    phi_at_zero=[]
    for i in range(len(nlist)):
        zeros.append(xi[i][-1])
        phi_at_zero.append(phi[i][-1])
        
    #print a table with calculated values    
    print_table(nlist,zeros,phi_at_zero)

    #solve neutron star central density
    ns = solve_ns(nlist,zeros,phi_at_zero)
    print '\nrho_ns/rho_nuc: %0.4f' % (ns)

    #solve high-rho WD mass
    wd = solve_wd(nlist,zeros,phi_at_zero)
    print '\nM_WD/M_sun: %0.4f' % (wd)
    
    reldiff,diff = [],[]
    for nn in [1]:
        n=nlist.index(nn)
        for i in range(len(xi[n])):
            minus=(theta[n][i]-vanalyt(n,xi[n][i]))
            diff.append(minus)
            reldiff.append(minus/theta[n][i])
        plt.clf()
        plt.plot(xi[n],diff,label='Absolute difference',color='blue')
        plt.plot(xi[n],reldiff,label='Relative difference',color='red')
        plt.ylim(-1E-13,1E-13)
        plt.legend()
        plt.show()

    
    #plot individual solutions for theta
    #for i in range(len(nlist)):
    #    plotter(nlist[i],xi[i],theta[i],r'$\theta$')
        
    #plot all solutions for theta
    plotall(nlist,xi,theta,r'$\theta$',16)
    
    #plot all solutions for rho/rho_c
    #plotall(nlist,xi,rho,r'$\rho/\rho_c$',6)
    
    #plot all solutions for phi
    #plotall(nlist,xi,phi,r'$\phi$',.16,-1,0)
    
    #K = const.hbar*const.c*(3*pi**2/(4**4*const.m_p**4))**.3333333333333333 #.333... is used instead of 1/3. to fix unit
