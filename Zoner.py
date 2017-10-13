import FDMesh as fdmesh
import numpy as np

class WireZoner:
    def __init__(self, L, W, fdmesh, vsp, hsp, nh, nv, no, IL=0.5, Eps_HiK=22.2):
        self.L         = L
        self.W         = W
        self.fdmesh = fdmesh
        self.vsp     = vsp
        self.hsp     = hsp
        self.nh     = nh
        self.nv     = nv
        self.no     = no
        self.IL     = IL

        # Property dictionary defines material properties
        self.prop_dict = { 'Silicon' : { 'eps' : 11.9 },
                           'Oxide'   : { 'eps' : 3.9  },
                           'HiK'     : { 'eps' : Eps_HiK },
                           'SiGe'    : { 'eps' : 16.0 } }

        self.buildMesh()


    def buildMesh(self):
        #------ Horizontal grid
        self.dx  = (self.hsp-self.IL) / (self.nh-1)
        self.dxo = self.IL / (self.no-1)

        # Additional vertical grids on left side
        xl = np.linspace(-self.hsp, -self.IL-self.dx, self.nh)


        # Additional vertical grids in the oxide
        xlo = np.linspace(-self.IL, -self.dxo, self.no)

        # Additional vertical grids on right side IL
        xro = np.linspace(self.fdmesh.L+self.dxo,
                         self.fdmesh.L+self.IL-self.dxo, self.no)
        #print 'xro: ', xro

        xr  = np.linspace(self.L+self.IL, self.L+self.hsp, self.nh)
        #print 'xr: ', xr

        # FDMesh horizontal grid
        xc = np.linspace(0.0, self.fdmesh.L, self.fdmesh.nx)


        self.xmesh = np.concatenate((xl, xlo, xc, xro, xr))
        #print 'xmesh: ', self.xmesh

        #------ Vertical grid
        self.dy  = (self.vsp-self.IL) / (self.nv-1)
        self.dyo = self.IL / (self.no-1)

        # Additional vertical grids on left side
        yl = np.linspace(-self.vsp, -self.IL-self.dy, self.nv)

        # Additional vertical grids in the oxide
        ylo = np.linspace(-self.IL, -self.dyo, self.no)

        # Additional vertical grids on right side IL
        yro = np.linspace(self.fdmesh.W+self.dyo,
                         self.fdmesh.W+self.IL-self.dyo, self.no)

        yr  = np.linspace(self.W+self.IL, self.W+self.vsp, self.nv)

        # FDMesh vertical grid
        yc = np.linspace(0.0, self.fdmesh.W, self.fdmesh.ny)

        self.ymesh = np.concatenate((yl, ylo, yc, yro, yr))
        #print 'ymesh: ', self.ymesh


    def getXMesh(self):
        return self.xmesh

    def getYMesh(self):
        return self.ymesh

    def zoner_func(self, x, y):
        print('Calling abstract function zoner_func.')
        raise RuntimeException

    def contact_func(self, x, y):
        print('Calling abstract function contact_func.')
        raise RuntimeException

    def plotting_rectangle_tuple(self):
        print('Calling abstract function plotting_rectangle_list')
        raise RuntimeException


class pGAAZoner(WireZoner):

    # The zoner function has to return a material type string for each x, y
    # coordinate pair
    def zoner_func(self, x, y):
        if (     x>=0.0
             and x<=self.L
             and y>=0.0
             and y<=self.W ):
            return 'Silicon'

        elif (   x>=0.0-self.IL
             and x<=self.L+self.IL
             and y>=0.0-self.IL
             and y<=self.W+self.IL ):
            return 'Oxide'
        else:
            return 'HiK'

    # The contact function returns True if a given coordinate is contained
    # within the gate electrode, False otherwise
    def contact_func(self, x, y):
        if (x<-self.hsp+1e-8) or (x>self.L+self.hsp-1e-8):
            return True
        else:
            return False


    def plotting_rectangle_tuple(self):
        return ( ((0.0,           0.0), self.L,             self.W),
                 ((-self.IL, -self.IL), self.L+2.0*self.IL, self.W+2.0*self.IL)
                )


class finFETZoner(WireZoner):

    # The zoner function has to return a material type string for each x, y
    # coordinate pair
    def zoner_func(self, x, y):
        if (     x>=0.0
             and x<=self.L
             and y>=0.0
             and y<=self.W ):
            return 'Silicon'

        elif (   x>=0.0-self.IL
             and x<=self.L+self.IL
             and y>=0.0-self.IL
             and y<=self.W+self.IL
             or  y<=0.0 ):
            return 'Oxide'
        else:
            return 'HiK'

    # The contact function returns True if a given coordinate is contained
    # within the gate electrode, False otherwise
    def contact_func(self, x, y):
        if (   (x<-self.hsp+1e-8         and y > 0.0)
            or (x>self.L+self.hsp-1e-8   and y > 0.0)
            or (y>self.W+self.vsp-1e-8)  and y > 0.0):
            return True
        else:
            return False

    def plotting_rectangle_tuple(self):
        return ( ( (0.0,             0.0), self.L,              self.W ),
                 ( (-self.IL,        0.0), self.L+2.0*self.IL,  self.W+self.IL ),
                 ( (-self.hsp, -self.vsp), self.L+2.0*self.hsp, self.vsp )
                )

class NWZoner(WireZoner):

    # The zoner function has to return a material type string for each x, y
    # coordinate pair
    def zoner_func(self, x, y):
        if (     x>=0.0
             and x<=self.L
             and y>=0.0
             and y<=self.W ):
            return 'Silicon'

        elif (   x>=0.0-self.IL
             and x<=self.L+self.IL
             and y>=0.0-self.IL
             and y<=self.W+self.IL ):
            return 'Oxide'
        else:
            return 'HiK'

    # The contact function returns True if a given coordinate is contained
    # within the gate electrode, False otherwise
    def contact_func(self, x, y):
        if (   (x<-self.hsp+1e-8     )
            or (x>self.L+self.hsp-1e-8)
            or (y>self.W+self.vsp-1e-8)
            or (y<-self.vsp+1e-8)     ):
            return True
        else:
            return False

    def plotting_rectangle_tuple(self):
        return ( ((0.0,           0.0), self.L,             self.W),
                 ((-self.IL, -self.IL), self.L+2.0*self.IL, self.W+2.0*self.IL)
                )


class MLFETZoner(WireZoner):

    # The zoner function has to return a material type string for each x, y
    # coordinate pair
    def zoner_func(self, x, y):
        if (     x>=0.0
             and x<=self.L
             and y>=0.0
             and y<=self.W ):
            return 'Silicon'

        elif (   x>=0.0
             and x<=self.L
             and (y<0.0 or y>self.W)):
            return 'SiGe'

        elif (   x>=0.0-self.IL
             and x<=self.L+self.IL):
            return 'Oxide'
        else:
            return 'HiK'

    # The contact function returns True if a given coordinate is contained
    # within the gate electrode, False otherwise
    def contact_func(self, x, y):
        if (x<-self.hsp+1e-8) or (x>self.L+self.hsp-1e-8):
            return True
        else:
            return False

    def plotting_rectangle_tuple(self):
        return ( ((0.0,           0.0), self.L,             self.W),
                 ((0.0,        self.W), self.L,           self.vsp),
                 ((0.0,     -self.vsp), self.L,           self.vsp),
                 ((-self.IL, -self.vsp),self.IL, self.W+2.0*self.vsp),
                 ((self.L,   -self.vsp),self.IL, self.W+2.0*self.vsp)
                )
