#import math
#print dir(math) #lists functions within math
#print help(math) #lists functions within math with explanation

#numpy is a package for multi-dimensional arrays designed for scientific computation
#numpy is ususally faster than standard python functions. 

import numpy as np

a=[1,2,5,5,6.888]
print a
print type(a) # a is a list!
print '\n'


b = np.array(a)
print b
print type(b) # b is now a numpy.ndarray
print b.shape # gives dimensions
print '\n'


c=np.matrix([[1,2],[3,4]])
print type(c)
print c.shape
#get values from matrix with np.get()


print np.arange(3) #faster than normal range()
print np.linspace(1,20,20) #last one is number of points in between
print np.logspace(0,2,10)
print '\n'

y = lambda x: x**2 #fast way of creating a function instead of using def
def y(x):
    return x**2

a=np.linspace(1,10,10)
print y(a)

#np.loadtxt() and np.savetxt()

#clever loading of txt and csv:
#np.genfromtxt()
#np.recfromcsv()

#start interactive session with pylab (from terminal):
#ipython --pylab

#from matplotlib import pyplot
from pylab import *

t=arange(0.0, 2.0, 0.01)
s=sin(2*pi*t)
xlabel('XLABEL')
ylabel('YLABEL')
title('TITLE')
grid(True)
#savefig('test.png')
plot(t,s)
#show()

#gnuplot can be imported in python

#scipy contains various scientific functions like fft, interpolation, numerical integration, signal processing, ...
from scipy.integrate import quad #most basic integration method. (1/2 (a+b)*n)

y= lambda x: x**2
print quad(y,0,4) #is answer,precision

#romberg integration:
from scipy.integrate import romb
