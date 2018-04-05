import numpy as np
from pymlpg import Circle
from scipy.spatial import cKDTree

class Nodes:
	methods = ['Irregular', 'Regular']
	def __init__(self, domain, N, dc = None, method = 'Regular'):
		"""
		Arguments:
			domain - instance of Domain
			N - number of nodes

		"""
		self.domain = domain
		
		if N:
			# Estimated: Number of Nodes
			self.N = N
			# Estimated: Average Nodal Spacing
			self.dc = np.sqrt(self.domain.area)/(np.sqrt(N)-1)
		elif dc:
			self.dc = dc
			self.N = self.domain.area/dc**2 + 2*self.domain.area/dc + 1
		else:
			raise ZeroDivisionError

		self.coor = np.zeros([self.N,2])

		if method == 'Regular':
			self.distribute_nodes_reg()
		else:
			self.distribute_nodes_irreg()

		self.tree = cKDTree(self.coor)
	
	def __getitem__(self,i):
		return self.coor[i]

	def __len__ (self):
		return len(self.coor)


	def distribute_nodes_reg(self):
		bounds = self.domain.bounds()
		
		lx = bounds[1,0] - bounds[0,0]
		ly = bounds[1,1] - bounds[0,1]
		
		bounds_area = lx*ly
		
		area = self.domain.area
		N = self.N*bounds_area/area
		
		nx = np.sqrt(lx/ly*N)
		ny = np.sqrt(ly/lx*N)

		x = np.linspace(bounds[0,0], bounds[1,0], nx)
		y = np.linspace(bounds[0,1], bounds[1,1], ny)
		xv, yv = np.meshgrid(x, y)
		self.coor = np.array(zip(xv.ravel(), yv.ravel()))

		mask = self.domain.contains(self.coor)
		for sub in self.domain.subdomains:
			mask = mask*(np.logical_not(sub.contains(self.coor)))
		self.coor = self.coor[mask]

	
	def distribute_nodes_irreg(self):
		if isinstance(self.domain, Circle):
			self.distribute_nodes_circle()

	
	def distribute_nodes_circle(self):
		'''
		This function Computes Fibonacci grid points inside a disk.

		Reference:
			Swinbank, Richard, and R. James Purser. "Fibonacci grids: 
			A novel approach to global modelling." Quarterly Journal 
			of the Royal Meteorological Society 132.619 (2006): 1769-1793.

		'''
		r0 = self.domain.radius/np.sqrt (float(self.N) - 0.5)
		c = self.domain.center
		
		phi = (1.0 + np.sqrt (5.0))/2.0
		self.coor = np.zeros ((self.N,2))
		idx = np.arange(self.N,dtype = np.float)
		gr = r0*np.sqrt(idx+0.5)
		gt = 2.0*np.pi*(idx+1)/phi
		self.coor[:,0] = c[0] + gr * np.cos(gt)
		self.coor[:,1] = c[1] + gr * np.sin(gt)


	def plot(self):
		from matplotlib import pylab
		self.domain.plot(True)
		pylab.scatter(self.coor[:,0], 
					self.coor[:,1],
					s = 4,
					color = 'k')
		pylab.show()
		