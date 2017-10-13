import pandas as pnd
import numpy  as np
import glob
import os
import os.path



# Read a set of data files, collate into a single Pandas DataFrame.
# The data files are all expected to have the same set of variables.
# Return value: Pandas DataFrame.

def FOMtoPandas(root='Split_', separator='_', fom_name='FOM.dat'):
		
	# Helper function to determine if a string represents a number		
	def is_number(s):
		try:
		    float(s)
		    return True
		except ValueError:
		    return False	
				

	# 1. Find all directories named 'Split_*'
	split_dirs = glob.glob(root+'*')
	split_dirs.sort()
	
	
	# Save all information into a dictionary, convert to Pandas DF
	table_dict = {}

	for cur_split in split_dirs:
		# Get split number
		(dummy, split_num) = cur_split.split(separator)
	
		# Now parse the FOM file
		fom_file = cur_split + '/' + fom_name
		if not os.path.isfile(fom_file):
			continue 
		fom = open(fom_file, 'r')
		for cur_line in fom:
			(name, val) = cur_line.split()
			# Add to table
			# Is this the first time we see this variable?
			if not name in table_dict:
				table_dict[name] = []
			# Is it an expected to be non-numerical?
			if is_number(val):
				val = float(val)
			table_dict[name].append(val)
		
		fom.close()

	# Create a pandas table
	df = pnd.DataFrame(table_dict)
	return df
