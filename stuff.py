a=0
b=10

while (a != b):  #repeats until value in () is false
    a += 1
    print a

for i in range(1, 10): # repeats replacing i by values in range(), can be any list
    print i
    
    
mylist=["Apple","Pear"]
for item in mylist:
    print item

def p():
    print mylist
    
anotherlist=["Carrots","Spinach"]
p()
mylist.extend(anotherlist) #adds anotherlist to mylist (mylist.append also works)
p()
mylist.remove("Spinach") #remove first occurence of Spinage
p()
mylist.index("Carrots") #gives index of first occurence of Bananap
p()
mylist.insert(3, "Food") # insert Food at index 3 (previous index 3 item is now the 4th)
p()

#copying is creating another pointer to the same space in memory!! e.g.:
mylist=["a","b"]
mylist_copy=mylist #mylist and mylist_copy are the same thing, so modifying one also changes the other
print mylist
print mylist_copy
mylist.pop()
print mylist
print mylist_copy

#it does work in this way:
mylist_copy=[]
mylist_copy.extend(mylist)
#file=open("file",'r') #open file read-only
#a=file.readline() #a is first line of file, then 'cursor' jumps to next line. Use something else if you don't want this
#a=file.readlines() #gives list of all lines. newline is noted by \n
#b=a[0].split() #a[0] is first line, split removes spaces and makes separate items
#a=file.read() #puts whole file in one string
#file.close() #closes file

import numpy as np
#np.loadtxt(file, skiprows=1, unpack=True) #opens file, skips first row and assigns each column to a string in the list
mydict={} #creates a dictionary, which is an unordened list
mydict["word"]="meaning"
mydict["age"]=21
print mydict

#classes
class my_class(object):
    def __init__(self, name, age): #self is no variable, but refers to created object
        self.name=name
        self.age=age
        
#a method is a function defined within in a class

person=my_class("Leon",21)
print person.name, person.age

person_2=my_class("Dude",20)
print person_2.name, person_2.age

class cube(object):
    def __init__(self,side):
        self.side=side
        
    def area(self):
        self.area=6*self.side**2
        return self.area
        
    def volume(self):
        self.volume=self.side**3
        return self.volume
        
import math

class sphere(object):
    def __init__(self,radius):
        self.radius=radius
        
    def area(self):
        self.area=4*math.pi*self.radius**2
        return self.area
    
    def volume(self):
        self.volume=4./3.*math.pi*self.radius**3
        return self.volume

print "\n"
cube_length=float(raw_input("Choose cube size\n"))
first_cube=cube(cube_length)
print "Area is",first_cube.area(),". Volume is",first_cube.volume(),"."

print "\n"
sphere_radius=float(raw_input("Choose sphere size\n"))
first_sphere=sphere(sphere_radius)
print "Area is",first_sphere.area(),". Volume is",first_sphere.volume(),"."


