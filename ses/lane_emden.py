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
import math

def analyt(xi):
	if n==0:
		ans=1-(xi**2)/6.
	elif n==1:
		try: 	# xi can be zero
			ans=np.sin(xi)/xi
		except ZeroDivisionError:
			ans=1.
		
	elif n==5:
		ans=(1+(xi**2)/3.)**-.5
	else:
		ans=0
	return ans
		
def deriv(z,xi):
    try:
        th=z[2]/z[0]**2
    except ZeroDivisionError:
        th=0.
	return np.array([1,th,-z[0]**2*z[1]**n])  #z[0]=xi, z[1]=theta, z[2]=phi

#define vectorized versions of funcs (essentially for-loop to apply the function to the elements of an array)
vanalyt=np.vectorize(analyt) #n is just a number an should be input as a single value
#vderiv=np.vectorize(deriv, excluded=['z']) #z should be input as an array
n=1 #value of polytrope index

steps=10000 #number of points at which the ODE is solved
xirange=np.linspace(0.0,10.0,steps) #discrete value for which de ode will be solved
zinit=np.array([1E-10,1.,0.]) # initial conditions. First one should be zero, but dividing by zero goes horribly wrong, so an extremely small value is chosen.
z=odeint(deriv,zinit,xirange) #acutally solve the ode for the given xi values. z is an array with xi, theta and phi. Phi isn't needed in the end, only for solving

real=analyt(xirange) #calculate analytical theta for given xi
diff=real-z[:,1] #absolute difference between analytical and numeric values
reldiff=diff/xirange #relative difference between analytical and numeric values


print z[1:100]
#figure(1)
#plot(xirange,analyt(xirange), color='red') #plots analytical solution for given values of xi
#plot(xirange,z[:,1],color='blue') #plots numeric solution for given values of xi
#plot(xirange,reldiff, color='red') #plots relative diff between analytical and numeric solutions
#xlabel('xi')
#ylabel('theta(analytical)/theta(numeric)')
#show()

#we see that the solutions differ for small xi, tot about 10^-4. Change the yrange to see the diff.
#figure(1)
#plot(xirange,reldiff, color='red') #plots relative diff between analytical and numeric solutions
#xlabel('xi')
#ylabel('theta(analytical)/theta(numeric)')
#ylim(-1E-7,1E-7)
#show()
#to get physical values, we use r=a xi and rho=rho_c theta**n
nlist=[0,1,3/2.,2,3,4]
sols=[]
rho=[] # will be rho/rho_c = theta**n
for i in nlist:
	n=i
	item=odeint(deriv,zinit,xirange)
	sols.append(item[:,1])
	rho.append(item[:,1]**n)

#we now have the numerical solution for 6 different n in the array sols.
#subplot(111) #row,col,figure
#produce files for theta vs xi and rho/rho_c (=theta**n) vs xi
for i in range(len(nlist)):
	plot(xirange,sols[i]) #theta vs xi
	xlim(0,10)
	ylim(0,1)
	xlabel('xi')
	ylabel('theta')
	savefig('theta_xi_'+str(i+1), bbox_inches='tight')
	clf() #clears the figure
	plot(xirange,rho[i])
	xlim(0,10)
	ylim(0,1)
	xlabel('xi')
	ylabel('rho/rho_c')
	savefig('rho_xi_'+str(i+1), bbox_inches='tight')
	clf()
	
#we see the density sometimes becomes negative, which is of course nonsense. The first the rho reaches zero is the outside of the star.

#NS are modelled by n=1. Numerical solution for n=1 is sols[1] or rho[1]
#show this one for clarity
clf()
plot(xirange,rho[1]) 
xlabel('xi')
ylabel('rho/rho_c')
xlim(0,4)
ylim(0,1)
show()

#find where rho goes to zero
temp=min(enumerate(rho[1][:.4*steps]), key=lambda x: abs(x[1])) #gives index, value. 4000 is chosen because the zero point lies before that and there are more zero points
minimum=[xirange[temp[0]],temp[1]] #replaces index by xi value of that index
print minimum
#the value for xi is remarkably close to pi
#this is correct, as the analytical solution is sin(xi)/xi, which is zero at xi=pi
