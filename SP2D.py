

import valley as valley
import FDMesh as fdmesh
import Schroedinger as sch

import numpy  as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import scipy.optimize    as optimize
import phys



class SP2D:

    def __init__(self, valley_array, fdmesh):
        self.valley_array = valley_array
        self.fdmesh       = fdmesh
        self.Sch_array    = []
        self.nvalleys     = 0
        self.eigen_solution = []
        self.min_energy   = 0.0
        self.max_energy   = 0.0

        # Build Schroedinger solvers, one for each valley
        for cur_valley in self.valley_array:
            self.Sch_array.append( sch.Schroedinger(cur_valley, self.fdmesh) )
            self.nvalleys += 1
            # Each valley stores the "eigen_solution"
            # The eigen_solution is a dictionary
            # with the format: { 'evals' : [], 'wfns' : (X,Y,Z) }

    def getValley(self, index):
        return self.valley_array[index]

    def solveWithCB(self, CB_function, n_eigen=12):
        # Clear solution
        self.eigen_solution = []

    # Helper function for determining the top of the CB
        def maxCB(self, CB_function):
            L = self.fdmesh.L
            W = self.fdmesh.W
            return CB_function(L / 2.0, W / 2.0)

        # Solve Schroedinger eigenproblem for each valley
        for cur_schroed in self.Sch_array:
            cb_max = maxCB(self, CB_function)
            (evals, (X,Y,Z)) = cur_schroed.solve(n_eigen,
                                                CB=CB_function,
                                                sigma=cb_max)

            self.eigen_solution.append(
                    { 'evals' : evals, 'wfns' : (X,Y,Z) }
            )
            #print 'Finished with numerical eigensolve ...'
            # Find energy range across all valleys
            all_evals = self.eigen_solution[0]['evals']
            for cur_eig in self.eigen_solution:
                all_evals = np.concatenate((all_evals, cur_eig['evals']))
            self.min_energy = np.min(all_evals)
            self.max_energy = np.max(all_evals)

    def getSolution(self, valley_index=-1):
        if valley_index == -1:
            return self.eigen_solution
        else:
            return self.eigen_solution[valley_index]

    def computeDOS(self, n_energy_points=300, e_min='none', e_max='none'):
        # Build energy axis
        if e_min == 'none':
            e_min = self.min_energy-0.05

        if e_max == 'none':
            e_max = self.max_energy+0.1

        energy = np.linspace(e_min, e_max, n_energy_points)
        ret_val = []
        for index, cur_valley in enumerate(self.valley_array):
            DOS = cur_valley.computeDOS(energy,
                                self.eigen_solution[index]['evals'])
            ret_val.append(DOS)

        return( (energy, ret_val) )



    def computeChargePlot(self, fmax, fermi_points=50, n_energy_points=300):
        # energy grid for integration
        energy = np.linspace(self.min_energy-0.05,
                             self.max_energy+0.1, n_energy_points)

        # Fermi energy grid
        fermi_vals  = np.linspace(energy[0]-3.0*phys.kT,
                                  fmax, fermi_points)

        charge_vals = np.zeros((self.nvalleys, fermi_points))
        tot_charge  = np.zeros(fermi_points)

        for i in range(0, fermi_points):
            if i % 5 == 0:
                print('Working on Fermi level:', fermi_vals[i])
            for valley_index, cur_valley in enumerate(self.valley_array):
                # Compute valley charge
                charge_vals[valley_index, i] = cur_valley.computeQMCharge(
                                fermi_vals[i],
                                self.eigen_solution[valley_index]['evals'])
                # Compute total charge
                tot_charge[i] += charge_vals[valley_index, i]

        return (fermi_vals, tot_charge, charge_vals)

    def computeSidewallCharge(self, Fermi):
        # energy grid for integration
        charge_vals = np.zeros(self.nvalleys)
        tot_charge = 0.0

        for valley_index, cur_valley in enumerate(self.valley_array):
            # Compute valley charge
            charge_vals[valley_index] = cur_valley.computeSidewallCharge(
                Fermi,
                self.eigen_solution[valley_index]['evals'])
            # Compute total charge
            tot_charge += charge_vals[valley_index]

        return (tot_charge, charge_vals)

    def computeChargeAndVelocity(self, Fermi):
        # Compute total charge
        tot_charge = 0.0
        valley_charge = []

        for valley_index, cur_valley in enumerate(self.valley_array):
            val_dict = self.getSolution(valley_index)
            valley_charge.append(0.0)

            print('Working on valley ',
                  valley_index, ' = ', valley_charge[valley_index])

            # Energies
            evals = val_dict['evals']

            # Loop over all states in valley
            # Break from loop when state energy is Ef+5*kT
            for index, cur_Ec in enumerate(evals):

                # Have we included enough states?
                # Don't break until some charge has been added
                # This prevents the problem of getting absolute zero
                # charge under very low bias conditions

                # compute the charge in each state
                state_charge = cur_valley.singleStateCharge(Fermi, cur_Ec)
                tot_charge += state_charge
                valley_charge[valley_index] += state_charge
                print('state charge=', state_charge, ' valley_charge=',
                      valley_charge[valley_index], ' tot_charge=', tot_charge)

        # Compute velocity integral, divide by total charge
        velocity = 1e-7 * self.computeVelocity(Fermi) / tot_charge

        return (velocity, tot_charge, valley_charge)

    def computeChargeDensity(self, Ef):
        ZCharge = np.zeros((self.fdmesh.ny, self.fdmesh.nx))

        for valley_index, cur_valley in enumerate(self.valley_array):
            val_dict = self.getSolution(valley_index)
            # Wavefunctions - not yet normalized
            (X, Y, Z) = val_dict['wfns']

            # Energies
            evals = val_dict['evals']

            # Loop over all states in valley
            # Break from loop when state energy is Ef+5*kT
            for index, cur_Ec in enumerate(evals):
                if cur_Ec > Ef + 5 * phys.kT:
                   break

                    # compute the charge in each state
                    # Normalize wavefunctions to box
                norm_sqr =  self.fdmesh.dx * self.fdmesh.dy * np.sum(np.sum(
                                Z[:,:,index] * Z[:,:,index] ) )
                state_charge = 1e21*cur_valley.singleStateCharge(Ef, cur_Ec)

                ZCharge += state_charge * Z[:,:,index] * Z[:,:,index] / norm_sqr
                #print 'valley: ', valley_index, ' including state ', cur_Ec, \
                #          ' charge=', state_charge, ' norm_sqr=', norm_sqr
        return ZCharge


    # Compute the total velocity integral across all modes and valleys
    # Must be normalized by total Q to get velocity
    def computeVelocity(self, Fermi):
        # energy grid for integration

        vel_sum = 0.0

        for valley_index, cur_valley in enumerate(self.valley_array):
            # Compute valley velocity for each valley
            vel_sum += cur_valley.computeQMVelocity(
                            Fermi,
                            self.eigen_solution[valley_index]['evals'])

        return vel_sum
