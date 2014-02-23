from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
#no astropy @ masterroom
#from astropy import constants as const
#from astropy import units as u
from math import pi
from scipy.interpolate import UnivariateSpline #interpolation
from scipy.optimize import brentq #rootfinding 

rc('text', usetex=True) #use latex for greek letters with a nicer font

# define Cash-Karp parameters
a2,   a3,  a4,  a5,  a6      =        1/5.,    3/10.,       3/5.,            1.,        7/8.
b21, b31, b32, b41, b42, b43 =        1/5.,    3/40.,      9/40.,         3/10.,      -9/10., 6/5.
b51, b52, b53, b54           =     -11/54.,     5/2.,    -70/27.,        35/27.
b61, b62, b63, b64, b65      = 1631/55296., 175/512., 575/13824., 44275/110592.,   253/4096. 
c1,   c2,  c3,  c4,  c5, c6  =     37/378.,       0.,   250/621.,      125/594.,          0.,  512/1771.
c1star, c2star, c3star, c4star, c5star, c6star = 2825/27648., 0.,  18575/48384.,13525/55296., 277/14336., 1/4.

def stepper(derivx, n, t, x, y, h, tol): # we have functions where x' is not a function of x, but of y and t so we have y'(x,t) and x'(y,t)

   k1 = h*derivx(n,t,y)
   k2 = h*derivx(n,t+a2*h,y+b21*k1)
   k3 = h*derivx(n,t+a3*h,y+b31*k1+b32*k2)
   k4 = h*derivx(n,t+a4*h,y+b41*k1+b42*k2+b43*k3)
   k5 = h*derivx(n,t+a5*h,y+b51*k1+b52*k2+b53*k3+b54*k4)
   k6 = h*derivx(n,t+a6*h,y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5)
   x_n_plus_1      = x +     c1*k1 +     c2*k2 +     c3*k3 +     c4*k4 +     c5*k5 +     c6*k6
   x_n_plus_1_star = x + c1star*k1 + c2star*k2 + c3star*k3 + c4star*k4 + c5star*k5 + c6star*k6
   DELTA           = x_n_plus_1 - x_n_plus_1_star
   try:
       h1 = h*abs(tol/DELTA)**0.2    # Finds step size required to meet given tolerance
   except ZeroDivisionError:
       h1 = h                        # When you are very close to ideal step, DELTA can be zero
   return x_n_plus_1, h1


def control(derivx,derivy,n,t0,x0,y0,t_max,h,tol,hmax):

        #set initial solutions lists
        t = np.array(t0)
        x = np.array(x0)
        y = np.array(y0)
       
        t_curr,x_curr,y_curr = t0,x0,y0
        count_iter,count_total = 0,0
        
        while x_curr >= -.1: #to -.1 to be able to interpolate the zero
            if t_curr > t_max:  #failsafe if theta doesn't go to zero
                break
            x_next,hx = stepper(derivx,n,t_curr,x_curr,y_curr,h,tol)
            y_next,hy = stepper(derivy,n,t_curr,y_curr,x_curr,h,tol)
            h1 = min(hx,hy)
            t_next=t_curr+h
            
            if h1 < 0.9*h:
                print 'Changed stepsize from %0.2e to %0.2e' % (h,h1) 
                h = h1
            elif h1 > 1.1*h:
                if h1 > hmax:
                    print 'Stepsize is now hmax'
                    h = hmax
                    
                    #still need to continue to prevent endless loop
                    t = np.append(t,t_next)
                    x = np.append(x,x_next)
                    y = np.append(y,y_next)
                    t_curr,x_curr,y_curr = t_next,x_next,y_next
                    count_iter += 1
                    
                else:
                    print 'Changed stepsize from %0.2e to %0.2e' % (h,h1)
                    h = h1
            else:            
                t = np.append(t,t_next)
                x = np.append(x,x_next)
                y = np.append(y,y_next)
                t_curr,x_curr,y_curr = t_next,x_next,y_next
                count_iter += 1
                
            count_total += 1
        
        print 'Calculation completed with %0.1f iterations in %0.1f steps' % (count_iter,count_total)
        return t,x,y


def dthetadxi(n,xi,phi):
    try: 
        return phi/xi**2
    except ZeroDivisionError:
        # the boundary conditions state theta'(0)=0
        return 0.
        
def dphidxi(n,xi,theta):
    try:
        return -xi**2*theta**n
    except ValueError:
        # fails when theta<0 and n!=int. Then the surface of the star is reached anyway, so we can safely return zero
        return 0.
        
def analyt(n,xi):
    if n==0:
        return 1-(xi**2)/6.
    elif n==1:
        try:    # xi can be zero
            return np.sin(xi)/xi
        except ZeroDivisionError:
            return 1.
    elif n==5:
        return (1+(xi**2)/3.)**-.5
    else:
        raise ValueError('This value of n has no analytical solution: '+str(n))        

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
    plt.title('All numerical solutions')
    plt.xlabel(r'$\xi$')
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()
        
def solver(nlist):
    xi0=0.
    ximax=20
    phi0=0.
    theta0=1.
    h=1E-5
    hmax=5E-4
    tol=1E-11
    
    xi,theta,phi,rho = [],[],[],[]
    
    for n in nlist:
        x,y,z = control(dthetadxi,dphidxi,n,xi0,theta0,phi0,ximax,h,tol,hmax)
        xi.append(x)
        theta.append(y)
        phi.append(z)
        if n == 3/2.: #exception because negative number ** non-int is imaginary. Last value is removed and 0 added to keep same length.
            y = np.delete(y,-1)
            y = np.append(y,0.)
        rho.append(y**n)
        
        
    return xi,theta,phi,rho

def find_zeros(x,y,z):
    zeros = np.array([])
    phi = np.array([])
    for i in range(len(x)):
        intpy = UnivariateSpline(x[i],y[i]) #interpolation of theta values to find where theta=0
        intpz = UnivariateSpline(x[i],z[i]) #interpolation of phi value to evaluate phi(xi=xi_0)
        zero = brentq(intpy,x[i][0],x[i][-1]) #calculate where theta=0
        zeros = np.append(zeros,zero)
        z_value = intpz(zero) #calculate phi at xi_0
        phi = np.append(phi,z_value)
        
    return zeros,phi
    
# equation for rho_c as given in eq. 8.28
def D_n(xi_0,phi_0):
    return 1/(-3*dthetadxi(n,xi_0,phi_0)/xi_0)

def rho_core(n,xi_0,phi_0,M,R):
    average=3*M/(4*pi*R**3)
    central=average*D_n(xi_0,phi_0)
    return central

if __name__=='__main__':
    
    #solve numerical
    nlist=[0,1,2,3,4]
    xi,theta,phi,rho = solver(nlist)
    
    #plot individual solutions for theta
    #for i in range(len(nlist)):
    #    plotter(nlist[i],xi[i],theta[i],r'$\theta$')
        
    #plot all solutions for theta
    #plotall(nlist,xi,theta,r'$\theta$',16)
    
    #plot all solutions for rho/rho_c
    #plotall(nlist,xi,rho,r'$\rho/\rho_c$',6)
    
    #plot all solutions for phi
    #plotall(nlist,xi,phi,r'$\phi$',16,-10,0)
    
    #K = const.hbar*const.c*(3*pi**2/(4**4*const.m_p**4))**.3333333333333333 #.333... is used instead of 1/3. to fix unit
    
    zeros,phi_at_zero = find_zeros(xi,theta,phi) #where theta=0 and corresponding phi(=xi**2*dthetadxi) values
    
    print 'n    xi1    phi'
    for i in range(len(nlist)):
        print '%0.2f %0.2f %0.2f ' % (nlist[i],zeros[i],phi_at_zero[i])
    
    #solve neutron star central density
    mass = 1.4*1.988435E33 #gram
    radius = 1E6 #centimeter
    rho_nuc = 4E14 #gram/cm3
    n = 1
    xi_0=zeros[1]
    phi_0 = phi_at_zero[1]
    rho = rho_core(n,xi_0,phi_0,mass,radius)
    #print 'rho_ns/rho_nuc: %0.2f' % (rho/rho_nuc)

    #solve high-rho WD mass
    n = 3
    xi_0 = zeros[nlist.index(n)]
    phi_0 = phi_at_zero[nlist.index(n)]
    #print 'xi_0: %0.2e' % xi_0
    #print 'phi_0: %0.2e' % phi_0
    
    #mass-radius relation for n=3 gives GM/m_3 = 16K**3/pi*G where m_3 = -xi_0**2 dthetadxi = - phi_at_zero
    K=4.88735E14 #cgs
    G=6.67384E-8 #cgs
    M=-phi_0*(16*K**3/(pi*G**3))**.5
    solarmass=1.988435E33
    print '\nM_WD/M_sun: %0.2f' % (M/solarmass)
