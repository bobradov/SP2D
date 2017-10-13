import numpy as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import matplotlib.pyplot as plt
import phys as phys


#+------------------------------------------------------------------
#|
#| Base valley class
#|
#+------------------------------------------------------------------

class valley(object):
    """ 
        Base valley class
        This is an abstract class which contains generic 
        valley functionality, but nothing specific to the
        valley type. 
        
        The functions which need to be defined by sub-classes:
        1. __init__       : parameters for the valley types are different
        2. singleStateDOS : this depends on the shape of the valley
        3. singleStateVel : this depends on the shape of the valley
        
        Other functions, such as the computation of total charge etc.
        are provided in the generic functionality
    """

    def __init__(self, name, L, W, mstar, mult, shift):
        self.name = name
        self.L = L
        self.W = W
        self.mstar    = mstar
        self.mstarDOS = self.mstar[2,2]
        self.mult     = mult
        self.shift    = shift

    
    def estimateEnergyLevel(self, n, m):
        mxx    = self.mstar[0,0]
        myy    = self.mstar[1,1]
        prefac = phys.hbar * phys.hbar / (2.0 * phys.m0)
        kx     = phys.pi / (self.L*1e-9)
        ky     = phys.pi / (self.W*1e-9)
        Ex     = n * n * prefac * kx * kx / (mxx * phys.q)
        Ey     = m * m * prefac * ky * ky / (myy * phys.q)
        return Ex + Ey

    
    def singleStateDOS(self, E, Ec):
        """ 
        Abstract method 
        Compute DOS for a single 1-D sub-band
        NOTE: singularity at E==Ec, and invalid for
        energies below that.
        Must make sure only energies > Ec are computed
        Result in #/(J-nm)
        """
        raise NotImplementedError

    # Compute DOS, based on a spectrum of energy levels
    # Use the energy argument (array) as the x-axis grid
    def computeDOS(self, energy, energy_levels):
        n_points = energy.shape
        DOS    = np.zeros(n_points)
        for cur_energy in energy_levels:
            DOS[energy > cur_energy] += self.singleStateDOS(
                            energy[energy > cur_energy], cur_energy)
        return DOS * self.mult

    # Fermi-Dirac occupancy of a single energy level
    def occupancy(self, energy, Ef):
        return 1.0/(1+np.exp((energy-Ef)/phys.kT))

    # Compute charge in a single 1-D state
    # This is done by integrating the DOS function
    # Integrand is singula and requires one-sided integration
    # returns charge value
    # Units are #/nm
    def singleStateCharge(self, Ef, Ec, limit=300, epsrel=1e-5):
        integrand = lambda x : self.singleStateDOS(x, Ec) * \
                               self.occupancy(x, Ef)

        (charge, abserr) = integral.quad( integrand,
                                            Ec, np.inf,
                                            limit=limit,
                                            epsrel=epsrel)
        return charge
        
    
    
    
    def singleStateVel(self, Ef, Ec, limit=300, epsrel=1e-5):
       """ Abstract function """
       raise NotImplementedError    



    # Compute QM charge
    # Ef is Fermi level
    # evals is an array of eigen-energies
    # n_energy_points is the size of the energy grid
    # on which to evalue the DOS (basis for charge calc.)
    def computeSidewallCharge(self, Ef, evals,
                        limit=300,
                        epsrel=1e-5):

        # Integrate charge across all bound states,
        # accumulate charge
        charge = 0.0
        for cur_state in evals:
            charge += self.singleStateCharge(Ef, cur_state,
                                                limit=limit,
                                                epsrel=epsrel)

        # Divide total charge by sidewall length
        return  1e14 * charge / (2.0*self.W)
        #return  1e14 * charge / (2.0*self.L)
        
    
    def computeQMVelocity(self, Ef, evals,
                        limit=300,
                        epsrel=1e-5):

        # Integrate velocity*charge product across all bound states,
        # accumulate charge
        vel_sum    = 0.0
        
        for cur_state in evals:
            vel = self.singleStateVel(Ef, cur_state,
                                                limit=limit,
                                                epsrel=epsrel)
            vel_sum    += vel
            
       
        return  vel_sum
