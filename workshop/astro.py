from astropy import constants as const
from astropy import units as u

print const.G,'\n'

print u.m
s=10*u.m #s = 10 meters
print s
print s.cgs #gives quantity in cgs units
print s + s.cgs #chooses first unit
print 'type(s)= ',type(s)
print 's.value= ',s.value
print 's.unit= ',s.unit

t=5*u.s
v=s/t
print 'v= ',v,'\n'

cms=u.cm/u.s
from astropy.units import imperial
mph=imperial.mile/u.h

from astropy.cosmology import Planck13 #imports 2013 results of Planck
#from astropy import cosmology #imporst cosmology
#cosmology.core.set_current(cosmology.core.Planck13) #sets current cosmology model
from math import pi

print Planck13.H(0) #gives H at certain z
print Planck13.age(0)

from astropy.cosmology import FlatLambdaCDM
cosmo=FlatLambdaCDM(H0=70, Om0=.3)
print cosmo.H(0)


# read/write files with astropy
import numpy as np
from astropy.table import table

x=np.asarray([1,6,5,3,29])
y_measurements=np.asarray([10,20,30,40,50])

measurements=table.Table([x,y_measurements],names=('x','y'))
print measurements #saving this as a file will organize it as a table

#there is an FITS module in astropy.
#HDF5 is a file format for storing of scientific data. Readable by multiple languages (c++, fortran, python, ...) e.g. AMUSE uses it
#astropy can read/write latex tables

measurements.write("table.tex",format='latex') #creates a latex table :D


#AMUSE www.amusecode.org



