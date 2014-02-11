import numpy as np

file_id='students.txt'
file_grades='grades.txt'

ids=np.loadtxt(file_id, dtype='string')
grades=np.loadtxt(file_grades, dtype='string')

length=len(ids)

students={}
gradelist={}
for i in range(length):
    students[ids[i][0]]=ids[i][1]   #ids[i][0] = name, ids[i][1] = studentnumber
    gradelist[grades[i][0]]=grades[i][1] # studentnumber, grade
    
#print students['Jan'] # gives studentnumber of Jan

print 'Name Grade'
for i in range(length):
    print ids[i][0],gradelist[students[ids[i][0]]]
