if __name__ == '__main__':
    import random
    
    length=int(raw_input("Give length of list\n"))
    unsorted_list=[]
    sorted_list=[]
    
    for i in range(length):
        unsorted_list.append(random.randint(0,10))
        
    print unsorted_list
    
    for i in range(length):
        item=min(unsorted_list)
        unsorted_list.remove(item)
        sorted_list.append(item)
        
    print sorted_list
