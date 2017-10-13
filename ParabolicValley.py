import numpy as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import phys as phys
from   valley import *

class ParabolicValley(valley):
   
    def __init__(self, name, L, W, mstar, mult, shift=0.0):
        """
        Constructor calls base class for base things, but
        adds custom construction explicity
        """
        super(ParabolicValley, self).__init__(name, L, W, 
                                              mstar, mult, shift)
        
        
        
    def singleStateDOS(self, E, Ec):
        """
        Compute single-state DOS at energy E
        Usese bound-state energy Ec
        """
        prefac = 2.0 / (phys.pi * phys.hbar)
        root   = np.sqrt(phys.m0 * self.mstarDOS/(2.0 * phys.q * (E-Ec)))
        return 1e-9 * prefac * root * self.mult * phys.q
        
    
    
    def singleStateVel(self, Ef, Ec, limit=300, epsrel=1e-5):
        """
        Computes the velocity integral for the single bound state
        Returns both the velocity integral and the associated charge
        Overall vel is assembled from calls to singleStateVel:
        
        vel = sum(vel_i)/sum(q_i)
        
        where vel_i is the velocity integral from the single state, and
        q_i is the charge integral from the single state.
        
        NOTE: vel_i actually has units of velocity * charge
        """
        v_integrand = lambda x : \
                      1e-7 * 1e-4 * 1e21 * self.singleStateDOS(x, Ec) * \
                      self.occupancy(x, Ef) * \
                      np.sqrt((2.0*phys.kT*phys.q /
                              (phys.pi*self.mstar[2, 2]*phys.m0)) * (x-Ec))
                               

        (velnum, abserr) = integral.quad( v_integrand,
                                            Ec, np.inf,
                                            limit=limit,
                                            epsrel=epsrel )
                                            
                                            
                                            
        return velnum
