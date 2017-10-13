
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as la_sparse

import phys as phys
import FVMMesh as mesh


class FVMPoisson:
    def __init__(self, p_mesh):
        self.pmesh = p_mesh


        self.Lap_mat = self.buildLaplaceMat()
        # Make sure Lap-matrix is Hermitian
        assert (z in self.Lap_mat for z in np.transpose(self.Lap_mat))



    def buildRHS(self, bc_dict, charge_func='none', charge_grid='none'):
        rows = self.pmesh.nVert
        rhs  = np.zeros(rows)

        # Add Dirichlet BCs
        for cur_cont in self.pmesh.contacts.keys():
            val = bc_dict[cur_cont]

            # Loop over Dirichlet nodes
            for cur_node_index in self.pmesh.contacts[cur_cont]:
                rhs[cur_node_index] = val

        #
        # Add charge to all nodes
        # Two options:
        #             1. Use function (x, y)-->charge_dens (charge_func)
        #             2. Map directly from a grid of semiconductor nodes
        #
        # Can be both, so charge is accumulated from function to grid.
        # If both options are 'none', Laplace solution is obtained
        #

        # Add charge using function
        if charge_func != 'none':
            charge_dens_func = charge_func
        elif charge_grid !='none':
            charge_dens_func = lambda x, y: charge_grid.nn_interp(x, y)

        # Loop over all elements
        # If a semiconductor element, add charge to all nodes of element
        for cur_elem in self.pmesh.getElementIter():
            if cur_elem.material != 'Silicon':
                continue

            # Shortcut references for local vertices
            vIter = cur_elem.getVertexIter(self.pmesh)
            (gb0, v0) = next(vIter)
            (gb1, v1) = next(vIter)
            (gb2, v2) = next(vIter)
            (gb3, v3) = next(vIter)

            # Get spacings for the element
            dx = v1.x-v0.x
            dy = v3.y-v0.y

            # Loop over all nodes of this element
            # Only a quarter of the nodal area of each node is
            # in the element
            for gb_index, vertex in cur_elem.getVertexIter(self.pmesh):
                #vertex = self.pmesh.vertex_list[global_index]
                x    = vertex.x
                y    = vertex.y
                dA   = 0.25 * dx * dy # Node gets one quarter of element
                qval = -1.0*100*charge_dens_func(x, y)*dA*1e-14*phys.q/phys.eps0
                rhs[gb_index] += qval


        return rhs



    def buildLaplaceMat(self):
        rows  = self.pmesh.nVert
        Lap_mat_dok = sparse.dok_matrix((rows, rows))

        # Loop over all elements
        for cur_elem in self.pmesh.getElementIter():
              # All vertex indices in element
              #print 'Working on element: ', cur_elem


            # Shortcut references for local vertices
            vIter = cur_elem.getVertexIter(self.pmesh)
            (gb0, v0) = next(vIter)
            (gb1, v1) = next(vIter)
            (gb2, v2) = next(vIter)
            (gb3, v3) = next(vIter)

              # Get spacings for the element
            dx = v1.x-v0.x
            dy = v3.y-v0.y

              # Coupling coefficients
            eps = self.pmesh.prop_dict[cur_elem.material]['eps']
            tx = eps * 0.5 * dy / dx
            ty = eps * 0.5 * dx / dy


            for loc_index, (global_index, cur_vert) in enumerate(
                                        cur_elem.getVertexIter(self.pmesh)):
                  # Sequence is:
                  # 3 -- 2
                # |    |
                # 0 -- 1
                #
                # Construct partial control volume for each node

                #--------------------------------------------------
                # Global row is equal to the global vertex index
                #--------------------------------------------------

                # Is it a Dirichlet node?
                # If so, simply assign a 1.0 into the matrix
                if cur_vert.boundary:
                    Lap_mat_dok[global_index, global_index] = 1.0

                # Not a Dirchlet node
                else:
                    if loc_index == 0:
                        # Couple right
                        row_r = v1.index
                        Lap_mat_dok[global_index, global_index] += tx
                        Lap_mat_dok[global_index, row_r]        -= tx

                        # Couple up
                        row_u = v3.index
                        Lap_mat_dok[global_index, global_index] += ty
                        Lap_mat_dok[global_index, row_u]        -= ty

                    if loc_index == 1:
                        # Couple up
                        row_u = v2.index
                        Lap_mat_dok[global_index, global_index] += ty
                        Lap_mat_dok[global_index, row_u]        -= ty

                        # Couple left
                        row_l = v0.index
                        Lap_mat_dok[global_index, global_index] += tx
                        Lap_mat_dok[global_index, row_l]        -= tx

                    if loc_index == 2:
                        # Couple down
                        row_d = v1.index
                        Lap_mat_dok[global_index, global_index] += ty
                        Lap_mat_dok[global_index, row_d]        -= ty

                        # Couple left
                        row_l = v3.index
                        Lap_mat_dok[global_index, global_index] += tx
                        Lap_mat_dok[global_index, row_l]        -= tx

                    if loc_index == 3:
                        # Couple down
                        row_d = v0.index
                        Lap_mat_dok[global_index, global_index] += ty
                        Lap_mat_dok[global_index, row_d]        -= ty

                        # Couple right
                        row_r = v2.index
                        Lap_mat_dok[global_index, global_index] += tx
                        Lap_mat_dok[global_index, row_r]        -= tx




        return Lap_mat_dok.tocsr()
