#Use Runge-Kutta integrator to 4th order
#3 levels:
#1: driver (starts/stops the integration)
#2: stepper (oversees incrementation of variable and calls algorithm) and output (container in which the stepper writes its output)
#3: algorithm (implements formulas)

#DthetaDxi=phi/xi**2
#DphiDxi=-xi**2*theta**n

#k1=hf(xn,yn)
#k2=hf(xn+h/2,yn+k1/2)
#k3=hf(xn+h/2,yn+k2/2)
#k4=hf(xn+h,yn+k3)
#yn+1=yn+k1/6+k2/3+k3/3+k4/6

from scipy.integrate import odeint
import numpy as np
from pylab import *

def analyt(xi):
	if n==0:
		return 1-(xi**2)/6.
	elif n==1:
		return np.sin(xi)/xi
	elif n==5:
		return (1+(xi**2)/3.)**-.5
	else:
		return 0

def deriv(z,xi):
	return np.array([1,z[2]/z[0]**2,-z[0]**2*z[1]**n])  #z[0]=xi, z[1]=theta, z[2]=phi

n=1
xirange=np.linspace(0.0,10.0,1000)
zinit=np.array([1E-10,1.,0.])
z=odeint(deriv,zinit,xirange)
	
x=np.linspace(0.0001,10.,500)
plot(x,analyt(x), color='red')

plot(xirange,z[:,1],color='blue')
xlabel('xi')
ylabel('theta')
show()
