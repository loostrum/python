import file_with_functions

from sys import argv

def main(filename="grades"):

    if len(argv)==2:
        scipt,filename=argv
        
    grades= file_with_functions.Read_Grades_data(filename)
    
if __name__ == '__main__':
    main()
    

#correcter way is using classes:
#class student(object):
#    def __init__(self,name):
#        self.name=name_
#        self.grade_list={}
