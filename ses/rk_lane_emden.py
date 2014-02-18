from matplotlib import pyplot as plt
from matplotlib import rc #to be able to print greek letters
from time import sleep
import numpy as np

rc('text', usetex=True) #use latex for greek letters
plt.close('all') #close previous plots if the exist


### Cash-Karp Parameters - From literature

a2,   a3,  a4,  a5,  a6      =        1/5.,    3/10.,       3/5.,            1.,        7/8.
b21, b31, b32, b41, b42, b43 =        1/5.,    3/40.,      9/40.,         3/10.,      -9/10., 6/5.
b51, b52, b53, b54           =     -11/54.,     5/2.,    -70/27.,        35/27.
b61, b62, b63, b64, b65      = 1631/55296., 175/512., 575/13824., 44275/110592.,   253/4096. 
c1,   c2,  c3,  c4,  c5, c6  =     37/378.,       0.,   250/621.,      125/594.,          0.,  512/1771.
c1star, c2star, c3star, c4star, c5star, c6star = 2825/27648., 0.,  18575/48384.,13525/55296., 277/14336., 1/4.

def control(derivy, derivz, tol, y0, z0, t0, h, tmax, v=False):
    '''
    This function takes in a python function that returns the derivative,
    a tolerance for error between RK5 and RK4, initial conditions on y and z (dependant) 
    and t (independant) as well as an initial step size.
    
    Keyword arguments:
    v - Verbose
    '''
    if v==True: print "Solving with initial condition (%0.2f, %0.2f), step size of %0.2f from t=0...%0.2f" % (y0, t0, h, tmax); sleep(1)

    y = [y0]      # Set up the initial conditions on y
    z = [z0]
    t = [t0]      # and t while creating the output lists
   
    t_curr, y_curr, z_curr, count, ncount = t0, y0, z0, 0, 0 # Setup counters and trackers

    while t_curr < tmax:
        t_next, y_next, z_next, h, h1 = stepper(derivy, derivz, t_curr, y_curr, z_curr, h, tol)
        if h1 < 0.9*h: 
            if v==True: print "Reduced step size from %0.4e to %0.4e at t = %0.2f" % (h, h1, t_curr)
            h = h1
        elif h1 > 1.1*h:
            if v==True: print "Increased step size from %0.4e to %0.4e at t = %0.2f" % (h, h1, t_curr)
            h = h1
        else:
            y.append(y_next)
            z.append(z_next)
            t.append(t_next)
            y_curr, z_curr, t_curr = y_next, z_next, t_next
            ncount += 1
        count += 1
    if v==True: print "Done. %i iterations, %i points" % (count, ncount )
    return y, z, t

def stepper(derivy, derivz, t, y, z, h, tol):
   '''
   This function is called by the control function to take
   a single step forward. The inputs are the derivative function,
   the previous time and function value, the step size in time (h),
   and the tolerance for error between 5th order Runge-Kutta and 4th
   order Runge-Kutta.
   '''
   # watch out! derivy has z as argument and vice versa
   k1 = h*derivy(t,z)
   k2 = h*derivy(t+a2*h,z+b21*k1)
   k3 = h*derivy(t+a3*h,z+b31*k1+b32*k2)
   k4 = h*derivy(t+a4*h,z+b41*k1+b42*k2+b43*k3)
   k5 = h*derivy(t+a5*h,z+b51*k1+b52*k2+b53*k3+b54*k4)
   k6 = h*derivy(t+a6*h,z+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5)
   y_n_plus_1      = y +     c1*k1 +     c2*k2 +     c3*k3 +     c4*k4 +     c5*k5 +     c6*k6
   y_n_plus_1_star = y + c1star*k1 + c2star*k2 + c3star*k3 + c4star*k4 + c5star*k5 + c6star*k6
   DELTAy           = y_n_plus_1 - y_n_plus_1_star
   
   try:
       hy = h*abs(tol/DELTAy)**0.2    # Finds step size required to meet given tolerance
   except ZeroDivisionError:
       hy = h                        # When you are very close to ideal step, DELTA can be zero
       
   k1 = h*derivz(t,y)
   k2 = h*derivz(t+a2*h,y+b21*k1)
   k3 = h*derivz(t+a3*h,y+b31*k1+b32*k2)
   k4 = h*derivz(t+a4*h,y+b41*k1+b42*k2+b43*k3)
   k5 = h*derivz(t+a5*h,y+b51*k1+b52*k2+b53*k3+b54*k4)
   k6 = h*derivz(t+a6*h,y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5)
   z_n_plus_1      = z +     c1*k1 +     c2*k2 +     c3*k3 +     c4*k4 +     c5*k5 +     c6*k6
   z_n_plus_1_star = z + c1star*k1 + c2star*k2 + c3star*k3 + c4star*k4 + c5star*k5 + c6star*k6
   DELTAz           = z_n_plus_1 - z_n_plus_1_star
       
   try:
       hz = h*abs(tol/DELTAz)**0.2
   except ZeroDivisionError:
       hz = h
       
   h1 = min([hy,hz])
       
   return t+h, y_n_plus_1, z_n_plus_1, h, h1


# define the lane emden equation as 2 ODEs

def dthetadxi(xi, phi):
    try:
        return phi/xi**2
    except ZeroDivisionError:
        return 1.   #when xi=0, phi=0 and the b.c. states dtheta/dxi = 1.
        
def dphidxi(xi, theta):
    n=1
    return -xi**2*theta**n

def analyt(xi):
	if n==0:
		return 1-(xi**2)/6.
	elif n==1:
		try: 	# xi can be zero
			return np.sin(xi)/xi
		except ZeroDivisionError:
			return 1.
	elif n==5:
		return (1+(xi**2)/3.)**-.5
	else:
		raise ValueError('This value of n has no analytical solution: '+str(n))

def plotter(a,b,c):
    plt.plot(c,a, color='red', label='Numerical')
    plt.plot(c,np.sin(c)/c,color='blue', label='Analytical')
    plt.xlim(0,5)
    plt.ylim(0,1)
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$\theta (\xi) $')
    plt.legend()
    plt.show()
    
def main():
    xi0=1E-10 #have to watch out for divide by zero
    phi0=0.
    theta0=1.
    h=.001
    ximax=5
    npts=int(ximax/h)
    tol=1E-13

    a,b,c = control(dthetadxi,dphidxi,tol,theta0,phi0,xi0,h,ximax,v=True) # control returns theta,phi,xi
    
    plotter(a,b,c)
    
if __name__=='__main__':
    main()
