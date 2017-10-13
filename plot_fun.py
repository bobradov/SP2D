import matplotlib.pyplot as plt
from   matplotlib.ticker import FormatStrFormatter
import numpy             as np
import scipy.interpolate as sp_interp
import scipy.optimize    as optimize


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



def plot_fun(sdf, cur_line, lw):



		fignum = 1

		# Charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['gate_ns']
		plt.semilogy(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg [V]')
		plt.ylabel('Gate ns [1/cm^2]')
		plt.grid(True)



		# Gate Capacitance
		plt.figure(fignum)
		fignum += 1
		vvals = sdf['Vg'].values
		qvals = sdf['gate_ns']
			# Interpolate onto sparse mesh for smoothing
		print('interpolating: vvals:', vvals.shape, ' qvals:', qvals.shape)
		#print sdf
		q_interpolator = sp_interp.interp1d(vvals, qvals, kind='cubic')
		v_i = np.linspace(vvals[0], vvals[-1], 10)
		q_i = q_interpolator(v_i)
		dV  = v_i[1]-v_i[0]
		cvals = 1e6*1.6e-19*np.gradient(q_i, dV)
			# Now interpolate back to original grid
		c_interpolator = sp_interp.interp1d(v_i, cvals, kind='cubic')
		cvals_final = c_interpolator(vvals)
		plt.plot(vvals, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg [V]')
		plt.ylabel('Gate Capacitance [uF/cm^2]')
		plt.grid(True)
		plt.ylim(ymin=0.0, ymax=3.5)

		# Gate Capacitance for total charge
		plt.figure(fignum)
		fignum += 1
		vvals = sdf['Vg'].values
		qvals_tot = sdf['tot_charge']*1e14
			# Interpolate onto sparse mesh for smoothing

		q_interpolator_tot = sp_interp.interp1d(vvals, qvals_tot, kind='cubic')
		v_i_tot = np.linspace(vvals[0], vvals[-1], 10)
		q_i_tot = q_interpolator_tot(v_i_tot)
		dV_tot  = v_i[1]-v_i[0]
		cvals_tot = 1e6*1.6e-19*np.gradient(q_i_tot, dV_tot)
			# Now interpolate back to original grid
		c_interpolator_tot = sp_interp.interp1d(v_i_tot, cvals_tot, kind='cubic')
		cvals_final_tot = c_interpolator_tot(vvals)
		if sdf['dev_type'].values[0] != 'fin':
			plt.plot(vvals, cvals_final_tot, cur_line, linewidth=lw)
		plt.xlabel('Vg [V]')
		plt.ylabel('Gate Capacitance [uF/cm]')
		plt.grid(True)
		plt.ylim(ymin=0.0, ymax=100)



		# Gate capacitance with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg'].values
		dV    = xvals[1]-xvals[0]
		qvals = sdf['gate_ns']
		cvals = 1e6*1.6e-19*np.gradient(qvals, dV)
			# Interpolant for capacitance
		#c_interp = sp_interp.interp1d(xvals, cvals)
		rhs = lambda x : c_interpolator(x)-0.15
			# Now find shift require to put cval=1.5 to Vg=0.0
		v0 = -0.05
		vf= optimize.fsolve(rhs, v0, xtol=1e-6)
		
		yvals = cvals_final
		plt.plot(vvals-vf, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Gate Capacitance [uF/cm^2]')
		plt.grid(True)
		plt.ylim(ymin=0.0, ymax=3.5)


		# Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['gate_ns']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Gate ns [1/cm^2]')
		plt.grid(True)

		# Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['charge_height']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Charge / height')
		plt.grid(True)

		# Total Charge with Vt shift
		# Don't plot for FinFET
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = 1e7*sdf['tot_charge']
		if sdf['dev_type'].values[0] != 'fin':
			plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Total Charge [1/nm]')
		plt.grid(True)

		# Average Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['average_charge']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Average Charge Conc. [1/cm^3]')
		plt.grid(True)


		# Delta2 Fraction vs Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['d2']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt, [V]')
		plt.ylabel('Delta 2 Fraction')
		plt.grid(True)

		# Delta2 Fraction vs sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['gate_ns']
		yvals = sdf['d2']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Gate ns [1/cm^2]')
		plt.ylabel('Delta 2 Fraction')
		plt.grid(True)
		plt.xlim(xmin=2e11, xmax=2e13)


		# Velocity vs. Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['vel']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Injection Velocity [cm/s]')
		plt.grid(True)

		# Velocity vs. sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['gate_ns']
		yvals = sdf['vel']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Gate ns [1/cm^2]')
		plt.ylabel('Injection Velocity [cm/s]')
		plt.grid(True)
		plt.xlim(xmin=2e11, xmax=2e13)

		# Velocity vs. average charge conc.



		plt.figure(fignum)
		plt.tick_params(axis='x', which='minor')
		ax = plt.gca()
		#ax.set_yscale('log')
		plt.tick_params(axis='x', which='minor')
		ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
		fignum += 1
		xvals = sdf['average_charge']
		yvals = sdf['vel']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Average Charge Conc. [1/cm^3]')
		plt.ylabel('Injection Velocity [cm/s]')
		plt.grid(True)
		plt.xlim(xmin=2e18, xmax=2e20)
