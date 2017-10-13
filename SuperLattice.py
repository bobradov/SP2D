#+---------------------------------------------------------
#|
#| Function for defining a superlattice CB offset
#| This needs to be added to the CB derived from
#| the electrostatic potential.
#|
#+---------------------------------------------------------

import matplotlib.pyplot as plt
import math as math

class SuperLattice:

	def __init__(self, period, thickness, energy, offset):
		self.period    = period
		self.thickness = thickness
		self.energy    = energy
		self.offset    = offset


	def evaluate(self, pos):
		# Determine the position of the evaluation point
		# within the supercell
		frac = ((pos-self.offset)/self.period - 
		               math.floor((pos-self.offset)/self.period))
		if frac < (self.period-self.thickness)/self.period:
			return 0.0
		else:
			return self.energy
		
		
if __name__ == "__main__":
	n_points  = 100
	period    = 5.0
	thickness = 1.0
	energy    = 1.0
	offset    = 3.0
	sl = SuperLattice(period, thickness, energy, offset)
	xvals = [ i*(10.0/(n_points-1)) for i in range(0, n_points) ]
	yvals = [ sl.evaluate(x) for x in xvals ]


	plt.figure(1)
	plt.plot(xvals, yvals, '-b')
	plt.ylim(ymin=-0.1, ymax = 1.5)
	plt.show()
	
	
