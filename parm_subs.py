#+-------------------------------------------------------------
#|
#| Function for parameter substitution
#| Takes a template file name (infile), output file 
#| name (outfile), and a dictionary of name, value pairs
#| (val_dict).
#| Expects the template to use the following substitution 
#| format:
#| text @keyword@ text ...
#| with the result:
#| text keyword_value text ...
#|
#+-------------------------------------------------------------- 



def parm_subs(infile, outfile, val_dict):

	inf  = open(infile, 'r')
	outf = open(outfile, 'w')
	
	for cur_line in inf:
		#print 'Got: ', cur_line
		
		# Does the line contain substitution targets?
		if '@' not in cur_line:
			outf.write(cur_line)
		else:
			# There may be substitution targets
			for cur_word in val_dict.keys():
				if '@'+cur_word+'@' in cur_line:
					#print 'Substitution for ', cur_word
					cur_line = cur_line.replace('@'+cur_word+'@', 
							    str(val_dict[cur_word]))
			outf.write(cur_line)
											
		
	outf.close()
	inf.close()
