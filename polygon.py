from swampy.TurtleWorld import *

world=TurtleWorld()
bob=Turtle()
bob.delay=.05


def polygon(t, n, l):
	angle=360./n
	for i in range(n):
		fd(t, l)
		rt(t, angle)

r=100
input=raw_input('How many sides do you want?\n')
n=int(input)

print 'Generating an '+str(n)+'-sided polygon using turtle '+str(bob)+'...'
polygon(bob,n,50)
print 'Done!'
wait_for_user()
