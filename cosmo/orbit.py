## This can be run from a python shell by typing
##    import orbit
## or from unix command line with
##    python orbit.py
from visual import *

## By default, "visual" creates a 3-D object called scene
scene.autoscale=0
scene.range=2

## create three objects, set their initial
##   position, radius, color, and other
##   properties: mass, momentum("p")
giant = sphere()
giant.pos = vector(.5,0,0)
giant.radius = 0.05 ; giant.color = color.cyan
giant.mass = 1
giant.p = vector(0, 0, 0)

dwarf = sphere()
dwarf.pos = vector(-.25,.133,.001)
dwarf.radius = 0.05 ; dwarf.color = color.green
dwarf.mass = 1
dwarf.p = vector(0,0,0)

moon = sphere()
moon.pos = vector(-.249,-.433,-.201)
moon.radius = 0.05 ; moon.color = color.red
moon.mass = 1
moon.p = vector(0,0,0)

## create 'curve' objects showing where we've been
for a in [giant, dwarf, moon]:
  a.orbit = curve(color=a.color, radius = 0.01)


def pstep( giant, dwarf ): 
  dist = dwarf.pos - giant.pos
  force = G * giant.mass * dwarf.mass * dist / ((mag(dist))*(mag(dist)+1E-3)**2)
  giant.p = giant.p + force*dt
  dwarf.p = dwarf.p - force*dt
  dist = dwarf.pos - giant.pos

dt = 0.005
G = 1E-4
while 1:
  ## set the picture update rate (100 times per second)
  rate(10000)

  pstep( giant, dwarf )
  pstep( giant, moon )
  pstep( moon, dwarf )

  for a in [giant, dwarf, moon]:
    a.pos = a.pos + a.p/a.mass * dt
    a.orbit.append(pos=a.pos)

## For an intro to visual python see
## http://wiki.aims.ac.za/mediawiki/index.php/Vpython:Getting_Started
