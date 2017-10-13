#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python


# Class for combinatorial looping


class CombinatorialLoop(object):
    def __init__(self, DictOfLists):
    	'Constructor for the iterator, takes a dictionary of lists as input'
    	# Create a ListOfLists which stores the possible values for each
    	# dimension. The iterator loops over the outer 
    	# product of all the dimensions.
    	# Heterogeneous list types are possible, no type
    	# checking is done.
    	#
    	# The .self value at each iteration is a list
    	# of elemements, equal in length to the size of 
    	# the ListOfLists. Each iteration produces a
    	# different element from the outer product.
    	#
    	
    	self.vars = DictOfLists.keys()
    	ListOfLists = []
    	for cur_var in self.vars:
    		ListOfLists.append(DictOfLists[cur_var])
    	
    	self.data   = ListOfLists
    	self.NLists = len( ListOfLists )
    	self.first  = True
    	
    	# How many items in each list?
    	self.ListLengths = []
    	for i in range(0,self.NLists):
    		self.ListLengths.append( len(ListOfLists[i]) )
    	
    	# Initialize counter for each list to 0	
        self.CurrentIndex = [0]*self.NLists
        
        # Initialize current state
        self.current = []
        for i in range(0,self.NLists):
        	self.current.append( ListOfLists[i][0] )
        
        # Find the final element of each list 
        self.high = []
        for i in range(0,self.NLists):
        	self.high.append( ListOfLists[i][ self.ListLengths[i]-1 ] )
        	
        	
        	
    
    def __iter__(self):
        'Returns itself as an iterator object'
        return self
        
        
        
        
	# The next function finds the next suitable n-tuple from the
	# outer product. 
	# The return value is an array (not a tuple) containing 
	# an element from the outer product.
	# The increment algorithm is a couunter, with the first list
	# representing the least-significant digit, the next list
	# represeting the next digit etc.
	# The 'counter' algorithm proceeds until all list elements are 
	# incremented to max values.

    def next(self):
        'Returns the next value until no more items remain'
        
        # Is this the very first use? If so, return
        # initial state.
        # Incrementing happens with next use of 'next'
        if self.first:
        	self.first = False
        	ret_dict = {}
        	for index, cur_var in enumerate(self.vars):
        		ret_dict[cur_var] = self.current[index]
        	return ret_dict
        
        # Increment lowest index; if already at max, reset 
        # and incement next index
        incremented = False
        
        for digit in range(0,self.NLists):
        	# Try to increment this index
        	self.CurrentIndex[digit] += 1
        	                
        	if self.CurrentIndex[digit] == self.ListLengths[digit]:
        		# Carry
        		# The current digit is not incremented
        		# Instead, it's reset to the zero index
        		# 'incremented' is not set to true, allowing
        		# the next digit to be incremented in the 
        		# next iteration of the loop
        		self.CurrentIndex[digit] = 0
        		self.current[digit] = \
        			self.data[digit][self.CurrentIndex[digit]]
        		continue
        	else:
        		# No carry, actually increment this digit
        		# Done with 'next', no more looping over digits
        		# Declare 'incremented' True and break out of loop
        		self.current[digit] = self.data[digit][self.CurrentIndex[digit]]
        		incremented  = True
        		break
        		
        # Termination condition: no digits could be incremented,
        # every digit is at max index
        # If termination condition is not met, return the current
        # state instead of raising exception	 
        if not incremented:
            raise StopIteration
        else:
            # Build a dictionary which contains the current
            # values of each variable
            ret_dict = {}
            for index, cur_var in enumerate(self.vars):
        		ret_dict[cur_var] = self.current[index]
            return ret_dict
            
            
    # The reset function is used for restarting the iterator
    # This is useful if one wishes to re-use the iterator, rather
    # than creating a new one.              
            
    def reset(self):
    	# Initialize counter for each list to 0	
        self.CurrentIndex = [0]*self.NLists
        
        # Initialize current state
        self.current = []
        for i in range(0,self.NLists):
        	self.current.append( ListOfLists[i][0] )
       	
       	# Restart count 	
       	self.first  = True
              
      
# Testsuite
      
if __name__ == "__main__":      
     # Do some testcases       
	counter = CombinatorialLoop( { 'alpha1' : ['a','b','c'],
	                               'int1'   : [1,2],
	                               'alpha2' : ['x','y','z'] } )             
		        
	for curVal in counter:
		print 'curVal=', curVal
