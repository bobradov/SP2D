# Class for looping over 'split tables'
# The split table is presented as a collecttion of
# iterables
# The iterables are expected to produce a dictionary at every iteration



class SplitTable:
	def __init__(self, iterables, default_dict={}):
		self.iterables    = iterables
		self.default_dict = default_dict
		
	def __iter__(self):
		return self
		
	# Get complete split dictionary
	# No checks are made whether the data for the
	# split is available from the current iterable
	# This is checked in 'next', possibly switching
	# to the next iterable	
	def build_split(self):
	
		# Get new split values
		curVal =  self.iterables[-1].next()
		cur_split = {}
		
		# Accumulate default values and split values
		# into a single return dictionary
		for def_var in self.default_dict:
			cur_split[def_var] = self.default_dict[def_var]
		for cur_var in curVal:
			cur_split[cur_var] = curVal[cur_var]
		return cur_split	
		
		
		
	# Standard iterator 'next' function
	# Switches between iterables as the data
	# from each is exhausted	
	def next(self):
		try:
			# If this iterable is not empty, it will be the next split
			return self.build_split()
		
		except(StopIteration):
			# Done with this iterable, look for next one
			# and return first split of the new iterable
			self.iterables.pop()
			
			if len(self.iterables) > 0:
				return self.build_split()
			else:
				raise StopIteration
			
			
