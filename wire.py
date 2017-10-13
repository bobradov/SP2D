#! /opt/local/bin/python2.7

import valley as valley
import numpy  as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import phys
import matplotlib.pyplot as plt



# Generate DOS for a single valley
n_points = 400
Delta2 = valley.valley(4.0, 20.0, 0.92, 0.19, 2)
#Delta2 = valley(4.0, 40.0, 0.05, 0.05, 1)
energy = np.linspace(0, 0.5, n_points)
DOS    = Delta2.computeDOS(energy, 5, 50)

# Generate charge density for the valley
Ef = 0.1
occDOS = DOS * Delta2.occupancy( energy, Ef )
DOS_interp = interp.interp1d( energy, DOS )


fermi_points = 50
fermi_vals  = np.linspace(0, 0.2, fermi_points)
charge_vals = np.zeros(fermi_points)
for i in range(0, fermi_points):
    if i % 5 == 0:
        print 'Working on Fermi level:', fermi_vals[i]
    (charge, abserr) = integral.quad( lambda x :
                        phys.q * DOS_interp(x) * Delta2.occupancy(x, fermi_vals[i]),
                        0.0, fermi_vals[i] + 5.0 * phys.kT,
                        limit=200,
                        epsrel=1e-3)
    #print 'charge: ', charge, ' with relerr=', abserr/charge
    charge_vals[i] = 1e14 * charge / (2.0*Delta2.L + 2.0*Delta2.W)

#print 'fermi:', fermi_vals
#print 'charge:', charge_vals


plt.figure(1)
plt.plot(energy, DOS, '-b', linewidth=3)
plt.plot(energy, occDOS, '-r', linewidth=3)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS')
plt.grid(True)

plt.figure(2)
plt.plot(fermi_vals, charge_vals, '-k', linewidth=3)
plt.xlabel('Fermi Level [eV]')
plt.ylabel('Electron Density [1/nm]')
plt.grid(True)

plt.show()
