from swampy.TurtleWorld import *
from math import pi

world=TurtleWorld()
bob=Turtle()


def polygon(t, n, r):
	angle=360./n
	length=2*pi*r/n
	for i in range(n):
		fd(t, length)
		lt(t, angle)

r=100
input=raw_input('How many sides do you want?\n')
n=int(input)

print 'Generating an',n,'sided polygon using turtle',bob,'...'
polygon(bob,n,r)
print 'Done!'
wait_for_user()
