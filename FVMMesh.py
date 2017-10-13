import numpy as np

class FVMVertex:
    def __init__(self, index, x, y, ind_i, ind_j, boundary=False):
        self.index = index
        self.x     = float(x)
        self.y     = float(y)
        self.ind_i = ind_i
        self.ind_j = ind_j
        self.boundary = boundary
        self.contact  = 'none'

    def __str__(self):
        tstr = '((' + str(self.index) +  ') '
        tstr += '(' + str(self.x) + ', ' + str(self.y) + ') '
        tstr += str(self.boundary)
        if self.boundary:
            tstr += ' ' + self.contact
        return tstr



class FVMElement:
    def __init__(self, Index, Material, VertexList):
        self.index      = Index
        self.material   = Material
        self.vertexList = VertexList

    def __str__(self):
        tstr = 'Element: ' + str(self.index) + ' '
        tstr += str(self.material) + ' '
        tstr += str(self.vertexList)
        return tstr

    def getVertexIter(self, mesh):
        return ElemVertexIter(self, mesh)

#+-----------------------------------------------------------------------+
#|
#| Finite Volume Mesh Class
#|
#| Supports cartesian meshes with arbitrary line distributions
#| Element properties are mapped onto the mesh via a zoner function
#| A dictionary linking material to properties is included in prop_dict
#|
#+-----------------------------------------------------------------------+


class FVMMesh:
    def __init__(self, xlines, ylines, zoner_func, prop_dict):
        # Check to make sure there are no duplicate lines
        xdpl = self.duplicates(xlines)
        if xdpl !=-1:
            print('Duplicated lines found in xlines:')
            print('Duplicated line:', xlines[xdpl])
            raise ValueError('Duplicated xline')

        ydpl = self.duplicates(ylines)
        if ydpl !=-1:
            print('Duplicated lines found in ylines:')
            print('Duplicated line:', ylines[ydpl])
            raise ValueError('Duplicated yline')

        self.xlines = xlines   # x-coordinates of lines
        self.ylines = ylines   # y-coordinates of lines
        self.nxl    = len(xlines)
        self.nyl    = len(ylines)
        self.prop_dict = prop_dict
        self.nVert     = self.nxl * self.nyl
        self.nElem     = (self.nxl-1)*(self.nyl-1)
        self.zoner     = zoner_func

        # Allocate mesh arrays
        self.vertex_list = []
        self.elem_list   = []
        self.contacts    = {}

        self.buildElements()

    # Determine of any of the lines are duplicated
    # If so, return the index of the duplicated line
    # Its duplicate is at index+1
    # Assumes the input array is sorted
    def duplicates(self, lines):
        for index in range(0, len(lines)-1):
            if lines[index]==lines[index+1]:
                return index
        return -1

    def getElementIter(self):
        return ElementIter(self)

    def mapSolutionToGrid(self, z):
        X, Y = np.meshgrid(self.xlines, self.ylines)
        Z    = z.reshape(len(self.ylines), len(self.xlines))
        return (X, Y, Z)

    def mapFunctionToGrid(self, func):
        X, Y = np.meshgrid(self.xlines, self.ylines)
        Z    = func(X, Y)
        return (X, Y, Z)

    def mapper(self, i, j):
        return i*self.nxl + j


    def loadAccumulator(self, solution):
        self.solution = solution

    # Nearest-neighbor interpolation for x, y
    def nn_interp(self, x, y):
        # find the xline which most closely matches x
        index_x      = (np.abs(self.xlines-x)).argmin()
        index_y      = (np.abs(self.ylines-y)).argmin()
        global_index = self.mapper(index_y, index_x)

        return self.solution[global_index]



    def __str__(self):
        # Header
        tstr = 'FVM Mesh:\n'
        tstr += 'nxl: ' + str(self.nxl) + '\n'
        tstr += 'nyl: ' + str(self.nyl) + '\n'
        tstr += 'nElem: ' + str(self.nElem) + '\n'

        # All elements
        for cur_elem in self.elem_list:
            tstr += cur_elem.__str__() + '    '
            # Get vertex coordinates
            for curv_index in cur_elem.vertexList:
                tstr += self.vertex_list[curv_index].__str__() + ' '
            tstr += '\n'

        # Contacts
        for cur_con in self.contacts.keys():
            tstr += cur_con + ' : '
            for cur_ind in self.contacts[cur_con]:
                tstr += str(cur_ind) + ' '
            tstr += '\n'

        return tstr

    def buildElements(self):
        # Loop over lines, create vertices, edges, elements
        # Vertical lines are left edge of element
        #
        v_index = 0
        for ind_i, yi in enumerate(self.ylines):
            for ind_j, xi in enumerate(self.xlines):
                # Add a new vertex
                self.vertex_list.append(FVMVertex(v_index,
                                                  xi, yi,
                                                  ind_i, ind_j))
                v_index += 1


        # Add elements
        # Loop over vertices, add element which has the vertex for
        # its top-left corner
        v_index    = 0
        elem_index = 0
        for yi in range(0, self.nyl-1):
            for xi in range(0, self.nxl-1):
                # Create nodelist
                # Local node numbering: counter-clockwise, start with bot left
                # 3 -- 2
                # |    |
                # 0 -- 1
                #
                node_list = [v_index,
                         v_index+1,
                         v_index+self.nxl+1,
                         v_index+self.nxl]

                # Determine material based on zoner function
                # Sample at the center of the element

                #print 'doing v_index=', v_index, ' nxl=', self.nxl, ' nyl=', self.nyl

                mat = self.zoner( 0.5*(self.vertex_list[v_index].x +
                                       self.vertex_list[v_index+1].x),
                                  0.5*(self.vertex_list[v_index].y +
                                         self.vertex_list[v_index+self.nxl].y ) )

                self.elem_list.append( FVMElement( elem_index,
                                                   mat,
                                                   node_list ) )
                v_index    += 1
                elem_index += 1
            # Completed row
            v_index += 1 # Move to next row

    def assignContact(self, cont_func, name):
        # Find all the nodes that belong to this contact
        self.contacts[name] = []

        for cur_vert in self.vertex_list:
            if cont_func(cur_vert.x, cur_vert.y):
                cur_vert.boundary = True
                cur_vert.contact  = name
                self.contacts[name].append(cur_vert.index)

#+---------------------------------------------------------------
#|
#| Iterators for FVMMesh
#|
#+---------------------------------------------------------------

# Return elements of mesh
# Returns mesh object, not index
class ElementIter:
    def __init__(self, mesh):
        self.mesh = mesh
        self.i    = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.i<self.mesh.nElem:
            self.i += 1
            return self.mesh.elem_list[self.i-1]
        else:
            raise StopIteration()

# Return vertices of an element
# Returns a tuple:
# (global_index_of_vertex, vertex_object)
class ElemVertexIter:
    def __init__(self, elem, mesh):
        self.mesh = mesh
        self.elem = elem
        self.i    = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.i < 4:
            self.i += 1
            global_node = self.elem.vertexList[self.i-1]
            return (global_node, self.mesh.vertex_list[global_node])
        else:
            raise StopIteration()
