# check if some string is a palindrome (str = reverse(str))
#basic way: n steps to take reverse of string of length n

def reverse(s):
    ''' str -> str
    Returns the reverse of a string
    
    >>> reverse('hello')        
    'olleh'
    '''
    
    # the docstring example can be tested by doctest with doctest.testmod (first import doctest)
    
    rev=''
    for char in s:
        rev = char + rev
    return rev
    
def is_palindrome(s):
    ''' str -> bool
    Returns true is the string is a palindrome
    '''
    return s == reverse(s)
    
    
#second way: reverse only second part of string and compare to first part. This is twice as fast. 
#Use integer division to determine what is first and second half
    
def is_palindrome_2(s):
    n = len(s)
    return s[:(n//2)] == reverse(s[(n-(n//2)):]) #// is integer division. e.g. racecar: n//2 is 3

#third way: compare last char to first one, n-1 with 2nd one etc.

#don't run the code if imported in another file, but only if main code:
if __name__ == '__main__':
    print 'this is the main code'
