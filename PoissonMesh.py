# Mesh for Poisson
# Superset of the Schroedinger mesh

import FDMesh as fdmesh
import numpy  as np

class PoissonMesh:
	
	def __init__(self, schr_mesh, vsp, hsp, eps_ox=3.9, 
							eps_semi=11.9, nh=3, nv=3):
		self.schr_mesh    = schr_mesh
		self.eps_ox       = eps_ox
		self.eps_semi     = eps_semi
		self.vsp          = vsp
		self.hsp          = hsp
		self.nh           = nh
		self.nv           = nv
		if self.nh <2:
			self.nh = 2
		if self.nv <2:
			self.nv = 2
			
		self.nx = self.schr_mesh.nx + 2*(self.nh-1)
		self.ny = self.schr_mesh.ny + 2*(self.nv-1) 
		
		self.buildMesh()
		
	def getXMesh(self):
		return self.xmesh	
		
	def getYMesh(self):
		return self.ymesh
		
	def getXInterfaces(self):
		return (self.nh-1, self.nh-1+self.schr_mesh.nx-1)
		
	def getYInterfaces(self):
		return (self.nv-1, self.nv-1+self.schr_mesh.ny-1)
		
	def getMeshSpacing(self):
		return (self.dx, self.dy)
		
	
		
	# Build the arrays of x and y coordinates
	# This mesh is non-uniform
	# The "center" part of the mesh comes from th Schroedinger Mesh
	# The outer portions are insulator regions, used only for Poisson	
	def buildMesh(self):
	
		#------ Horizontal grid
		self.dx = self.hsp / (self.nh-1)
		
		# Additional vertical grids on left side
		xl = np.linspace(-self.hsp, 0-self.dx, self.nh-1)
		
		# Additional vertical grids on right side
		xr = np.linspace(self.schr_mesh.L+self.dx, 
						 self.schr_mesh.L+self.hsp, self.nh-1)
		
		# FDMesh horizontal grid
		xc = np.linspace(0.0, self.schr_mesh.L, self.schr_mesh.nx)
		
		self.xmesh = np.concatenate((xl, xc, xr))
		print 'xmesh: ', self.xmesh 
		
		#------ Vertical grid
		self.dy = self.vsp / (self.nv-1)
		
		# Additional vertical grids on bottom
		yb = np.linspace(-self.vsp, 0-self.dy, self.nv-1)
		
		# Additional vertical grids on top
		yt = np.linspace(self.schr_mesh.W+self.dy, 
						 self.schr_mesh.W+self.vsp, self.nv-1)
		
		# FDMesh vertical grid
		yc = np.linspace(0.0, self.schr_mesh.W, self.schr_mesh.ny)
		
		self.ymesh = np.concatenate((yb, yc, yt))
		print 'ymesh: ', self.ymesh
	
