def count_vowels(s):
    ''' (str) -> int
    Returns the number of vowels in str
    '''
    vowels=["a","e","i","o","u","A","E","I","O","U"]
    found=0
    for letter in vowels:
        found += s.count(letter)
    print "There are "+str(found)+" vowels in "+s+".\n"
    
while True:
    input=raw_input("Give a word ('stop' to exit)\n")
    if (input == 'stop'):
        break
    count_vowels(input)
