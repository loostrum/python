def swap(l,a):
    ''' 
    list -> list
    swaps l[a] with l[a-1]
    
    >>> swap([1, 2, 3, 4],1)
    [2, 1, 3, 4]
    
    >>> swap([1, 2, 3, 4],2)
    [1, 3, 2, 4]
    '''
    
    item=l[a]
    l[a]=l[a-1]
    l[a-1]=item
    return l
    
def need_swap(l,a):
    '''
    True if l[a] < l[a-1]
    
    >>> need_swap([1, 2, 3, 4],1)
    False
    >>> need_swap([2, 1, 3, 4],1)
    True
    '''
    
    if (l[a]<l[a-1]):
        return True
    else:
        return False
    

if __name__ == '__main__':
    import random
    
    length=int(raw_input("Give length of list\n"))
    unsorted_list=[]
    sorted_list=[]
    
    for i in range(length):
        unsorted_list.append(random.randint(0,10))
    
    print unsorted_list
            
    sorted_list.append(unsorted_list[0]) #sort first element
    unsorted_list.pop(0)        
            
    for i in range(1,length): #first element is already sorted so skip index 0
        sorted_list.append(unsorted_list[0])
        unsorted_list.pop(0)
            
        #now apply bubble sort to sorted_list
        count=None
                
        while (count!=0):
            count=0
            for n in range(1,i):
                if need_swap(sorted_list,n):
                    swap(sorted_list,n)
                    count+=1
        
            
    print sorted_list
        
        
        
