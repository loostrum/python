from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np

rc('text', usetex=True) #use latex for greek letters with a nicer font

# define Cash-Karp parameters
a2,   a3,  a4,  a5,  a6      =        1/5.,    3/10.,       3/5.,            1.,        7/8.
b21, b31, b32, b41, b42, b43 =        1/5.,    3/40.,      9/40.,         3/10.,      -9/10., 6/5.
b51, b52, b53, b54           =     -11/54.,     5/2.,    -70/27.,        35/27.
b61, b62, b63, b64, b65      = 1631/55296., 175/512., 575/13824., 44275/110592.,   253/4096. 
c1,   c2,  c3,  c4,  c5, c6  =     37/378.,       0.,   250/621.,      125/594.,          0.,  512/1771.
c1star, c2star, c3star, c4star, c5star, c6star = 2825/27648., 0.,  18575/48384.,13525/55296., 277/14336., 1/4.

def stepper(derivx, n, t, x, y, h, tol): # we have functions where x' is not a function of x, but of y and t so we have y'(x,t) and x'(y,t)
   '''
   This function is called by the control function to take
   a single step forward. The inputs are the derivative function,
   the previous time and function value, the step size in time (h),
   and the tolerance for error between 5th order Runge-Kutta and 4th
   order Runge-Kutta.
   '''

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


def control(derivx,derivy,n,t0,x0,y0,t_max,h,tol):

        #set initial solutions lists
        t = np.array(t0)
        x = np.array(x0)
        y = np.array(y0)
       
        t_curr,x_curr,y_curr = t0,x0,y0

        while x_curr >= 0:
            if t_curr > t_max:  #failsafe if theta doesn't go to zero
                break
            x_next,hx = stepper(derivx,n,t_curr,x_curr,y_curr,h,tol)
            y_next,hy = stepper(derivy,n,t_curr,y_curr,x_curr,h,tol)
            h1 = min(hx,hy)
            t_next=t_curr+h
            
            if h1 < 0.9*h: 
                h = h1
            elif h1 > 1.1*h:
                h = h1
            else:            
                t = np.append(t,t_next)
                x = np.append(x,x_next)
                y = np.append(y,y_next)
                t_curr,x_curr,y_curr = t_next,x_next,y_next
            
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
    plt.xlim(0,x[-1])
    plt.ylim(0,1)
    plt.legend()
    plt.show()

def plotall(nlist,x,y,ylabel,xlim):
    plt.clf()
    for i in range(len(nlist)):
        plt.plot(x[i],y[i],label='n = '+str(nlist[i]))
    plt.xlim(0,xlim)
    plt.ylim(0,1)
    plt.title('All numerical solutions')
    plt.xlabel(r'$\xi$')
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()
        
def solver():
    nlist=[0,1,1.5,2,3,4]
    xi0=0.
    ximax=20
    phi0=0.
    theta0=1.
    h=1E-4
    tol=1E-10
    
    xi,theta,phi,rho = [],[],[],[]
    
    for n in nlist:
        x,y,z = control(dthetadxi,dphidxi,n,xi0,theta0,phi0,ximax,h,tol)
        xi.append(x)
        theta.append(y)
        phi.append(z)
        rho.append(y**n)
        
    return nlist,xi,theta,phi,rho

def main():
    
    #solve numerical
    nlist,xi,theta,phi,rho = solver()
    
    #plot individual solutions for theta
    #for i in range(len(nlist)):
    #    plotter(nlist[i],xi[i],theta[i],r'$\theta$')
        
    #plot all solutions for theta
    #plotall(nlist,xi,theta,r'$\theta$',16)
    
    #plot all solutions for rho/rho_c
    #plotall(nlist,xi,rho,r'$\rho/\rho_c$',6)
    
if __name__=='__main__':
    main()
