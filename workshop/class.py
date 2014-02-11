class physics(object):
    def __init__(self,m,s,p):
        self.ma=m
        self.spe=s
        self.positi=p
     
    def mass(self):
        return self.ma
        
    def speed(self):
        return self.spe
        
    def position(self):
        return self.positi
        
    def ekin(self):
        return .5*self.ma*self.spe**2
    
    def epot(self):
        G=6.67E-11
        if (self.positi == 0):
            return 'Position is zero, potential energy undefined'
        return G*self.ma/self.positi
        
star = physics(10E30,10,100)

opts=[star.mass(),star.speed(),star.position(),star.ekin(),star.epot()]

while True:
    input=int(raw_input("What do you want to know?\n1) mass \n2) speed\n3) position\n4) ekin\n5) epot\n"))-1

    if input > 5:
        print 'index out of bounds'
    for i in range(5):
        if (i==input):
            print opts[i],'\n'            

