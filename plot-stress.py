#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python 

import pandas as pnd
import numpy  as np
import glob
import os
import matplotlib.pyplot as plt
from   FOMtoPandas import FOMtoPandas
import plot_fun as pf

# Create a new table, add a row for each split
# Use data from FOM.dat in each split

plot_types = ['pgaa', 'fin']

df = FOMtoPandas()
print df
#exit()
df = df[df['converged']=='True']
#df = df[df['vsp'] !=3.0 ]

# Average charge in channel slice
df['average_charge'] = 1e21*df['tot_charge']/(df['L']*df['W'])

# Charge per unit height
df['charge_height'] = 1e14*df['tot_charge']/(df['W']+ 2.0*df['vsp'])

# Create sidewall charge
df['hsp']    = 2.49
#df['fin_ns'] = df['tot_charge']*1e14/(2.0*df['W'] + df['L'] +
#										2.0*df['hsp'] + 2.0*df['vsp'])
										
df['fin_ns'] = df['tot_charge']*1e14/(2.0*df['W'] + df['L'])
										
df['sw_ns']  = df['tot_charge']*1e14/(2.0*(df['W']+2.0*df['vsp']))

df['gaa_ns'] = df['tot_charge']*1e14/(2.0*(df['W']+df['L'] + 
										   2.0*df['vsp'] + 2.0*df['hsp']))
										 
df['gate_ns'] = 0.0
df.ix[(df['dev_type']=='fin'), 'gate_ns'] = df.ix[ (df['dev_type']=='fin'), 'fin_ns']
df.ix[(df['dev_type']=='pgaa'),'gate_ns'] = df.ix[ (df['dev_type']=='pgaa'),'sw_ns'] 
df.ix[(df['dev_type']=='nw'), 'gate_ns'] = df.ix[ (df['dev_type']=='nw'), 'gaa_ns']  
		
df['vel'] = df['vel']*1.05								   
						 

df.to_csv(path_or_buf='done.csv')
print df.head()


# Make some evil plots
lines   = ['-','--','-.',':','-o', '-D']
colors  = ['r','k','b','m']
symbols = ['o', '']
markersize = 7
lw  = 3

# Plotting
pf.set_plt_defaults()

cat_var = 'd2_shift'

print df['d2_shift'].unique()
#exit()

# Create a list of data frames, one for each plot type
dev_df = [ df[df['dev_type']==x] for x in plot_types ]

for dev_index, cur_dev_df in enumerate(dev_df):
	cur_col = colors[dev_index]
	
	for d2_index, d2_val in enumerate(sorted(cur_dev_df['d2_shift'].unique())):
		print ' >>> Working on d2 = ', d2_val
		cat_df  = cur_dev_df[cur_dev_df['d2_shift']==d2_val]
		cur_symbol = symbols[d2_index]
		if cur_col == 'k':
			cur_symbol = ''
	
		for index, cur_height in enumerate(np.sort(cat_df['W'].unique())):
			sdf = (cat_df[cat_df['W']==cur_height]).sort(columns=['Vg'], ascending=True)
			cur_line = lines[index] + cur_symbol + cur_col	
		
			pf.plot_fun(sdf, cur_line, lw)
	
plt.show()
	
	
	
	


		
	
