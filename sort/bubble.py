def swap(l,a):
    ''' 
    list -> list
    swaps l[a] with l[a+1]
    
    >>> swap([1, 2, 3, 4],1)
    [1, 3, 2, 4]
    
    >>> swap([1, 2, 3, 4],2)
    [1, 2, 4, 3]
    '''
    
    item=l[a]
    l[a]=l[a+1]
    l[a+1]=item
    return l

def need_swap(l,a):
    '''
    True if l[a] > l[a+1]
    
    >>> need_swap([1, 2, 3, 4],1)
    False
    >>> need_swap([1, 3, 2, 4],1)
    True
    '''
    
    if (l[a]>l[a+1]):
        return True
    else:
        return False


if __name__ == '__main__':
    import random
    import doctest
    doctest.testmod()
    
    
    length=int(raw_input("Give length of list\n"))
    to_sort=[]

    for i in range(length):
        to_sort.append(random.randint(0,10))

    count=None
    print "Random list of length",length,":",to_sort, "\n"
    
    while (count!=0):
        count=0
        for n in range(length-1):
            if need_swap(to_sort,n):
                swap(to_sort,n)
                count+=1

    print "Sorted list :",to_sort
