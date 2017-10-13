#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python 

import pandas as pnd
import numpy  as np
import glob
import os
import matplotlib.pyplot as plt
import scipy.interpolate as sp_interp
import scipy.optimize    as optimize
from   FOMtoPandas import FOMtoPandas

# Create a new table, add a row for each split
# Use data from FOM.dat in each split

df = FOMtoPandas()
print df.head()
df = df[df['converged']=='True']
df = df[df['vsp'] ==2.4923]
df.to_csv(path_or_buf='done.csv')
print df.head()

# Make some plots
# Make some evil plots
def set_plt_defaults():
	plt.rcParams['figure.facecolor'] = 'w'
	plt.rc('xtick', labelsize=18) 
	plt.rc('ytick', labelsize=18) 
	font = {'family' : 'Bitstream Vera Sans',
		    'weight' : 'bold',
		    'size'   : 16}

	plt.rc('font', **font)
	plt.rcParams['toolbar'] = 'None'

	plt.tick_params(which='major', length=15, width=2)
	plt.tick_params(which='minor', length=5, width=2)

	plt.rcParams['axes.linewidth'] = 2 

lines   = ['-','--','-.',':','-o', '-D']
colors  = ['r','b','k','m']
symbols = ['-o','D','s','^','D']
markersize = 7
lw  = 3

# Plotting
set_plt_defaults()

cat_var = 'vsp'


for cat_index, cat_val in enumerate(np.sort(df[cat_var].unique())):
	cat_df = df[df[cat_var]==cat_val]
	cur_col = colors[cat_index]

	for index, cur_height in enumerate(np.sort(df['W'].unique())):
		sdf = (cat_df[cat_df['W']==cur_height]).sort(columns=['Vg'], ascending=True)
		cur_line = lines[index] + cur_col	
		#cur_line = symbols[index] + 'b'
		fignum = 1	
		
		# Charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['sw_charge']
		plt.semilogy(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('ns')
		plt.grid(True)
		
		
	
		# Gate Capacitance
		plt.figure(fignum)
		fignum += 1
		vvals = sdf['Vg'].values
		qvals = sdf['sw_charge']
			# Interpolate onto sparse mesh for smoothing
		q_interpolator = sp_interp.interp1d(vvals, qvals, kind='cubic')
		v_i = np.linspace(vvals[0], vvals[-1], 10)
		q_i = q_interpolator(v_i)
		dV  = v_i[1]-v_i[0]
		cvals = 1e6*1.6e-19*np.gradient(q_i, dV)
			# Now interpolate back to original grid
		c_interpolator = sp_interp.interp1d(v_i, cvals, kind='cubic')
		cvals_final = c_interpolator(vvals)
		plt.plot(vvals, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('Gate Capacitance')
		plt.grid(True)
		
		
		
		# Gate capacitance with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg'].values
		dV    = xvals[1]-xvals[0]
		qvals = sdf['sw_charge']
		cvals = 1e6*1.6e-19*np.gradient(qvals, dV)
			# Interpolant for capacitance
		c_interp = sp_interp.interp1d(xvals, cvals)
		rhs = lambda x : c_interp(x)-1.0
			# Now find shift require to put cval=1.5 to Vg=0.0
		v0 = 0.0
		vf= optimize.fsolve(rhs, v0, xtol=1e-5)
		print 'Got vf=', vf
		yvals = cvals_final
		plt.plot(vvals-vf, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('Gate Capacitance')
		plt.grid(True)
		
		# Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['sw_charge']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('ns')
		plt.grid(True)
		
		
		# Delta2 Fraction vs Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['d2']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('d2_frac')
		plt.grid(True)
		
		# Delta2 Fraction vs sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['sw_charge']
		yvals = sdf['d2']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('d2_frac')
		plt.grid(True)
		plt.xlim(xmin=1e11, xmax=1e14)
	
	
		# Velocity vs. Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['vel']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('velocity')
		plt.grid(True)
	
		# Velocity vs. sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['sw_charge']
		yvals = sdf['vel']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('ns')
		plt.ylabel('velocity')
		plt.grid(True)
		plt.xlim(xmin=1e11, xmax=1e14)
	
plt.show()
	
	
	
	


		
	
