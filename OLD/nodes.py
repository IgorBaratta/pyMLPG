from mshr import *
from dolfin import *
import numpy as np
from matplotlib import pyplot as plt
from functools import reduce


class nodes:
	def __init__(self, coord):
		dmed = 0.1
		self.coord = coord
		self.xmax = max(coord[:,0])
		self.ymax = max(coord[:,1])
		self.xmin = min(coord[:,0])
		self.ymin = min(coord[:,1])
		self.size = coord.shape[0]
		self.type = np.zeros(self.size)
		self.r = np.ones(self.size)*dmed
		self.left = []
		self.right = []
		self.bot = []
		self.top =[]
		self.set_radius()
		self.on_boundary()
		self.mark(self.left,1)
		self.mark(self.right,2)
		self.mark(self.bot,3)
		self.mark(self.top,4)
		self.boundary = reduce(np.union1d,(self.left[0],self.right[0],self.bot[0],self.top[0]))

	def mark(self,ind,type):
		self.type[ind] = type
	
	def on_boundary(self):
		self.left.append(np.where(self.coord[:,0]==self.xmin)[0])
		self.right.append(np.where(self.coord[:,0]==self.xmax)[0])
		self.bot.append(np.where(self.coord[:,1]==self.ymin)[0])
		self.top.append(np.where(self.coord[:,1]==self.ymax)[0])
	
	def plot(self):
		colors = []
		plt.scatter(self.coord[:,0],self.coord[:,1])
		plt.scatter(self.coord[self.left,0],self.coord[self.left,1],color='red',s=50)
		plt.scatter(self.coord[self.right,0],self.coord[self.right,1],color='red',s=50)
		plt.scatter(self.coord[self.top,0],self.coord[self.top,1],color='red',s=50)
		plt.scatter(self.coord[self.bot,0],self.coord[self.bot,1],color='red',s=50)
		fig = plt.gcf()
		ax = fig.gca()
		for idx, node in enumerate(self.coord):
			if idx not in self.boundary:
				circle = plt.Circle((node[0], node[1]), self.r[idx], color='b',alpha=0.1)
				ax.add_artist(circle)
		plt.axis('equal')
	
	def set_neighborhood(self,neighborhood):
		self.neighbor = neighborhood
	
	def get_support(self,idx):
		support = self.neighbor[idx]
		j = 0
		while support.size<6:
			neighbor = self.neighbor[support[j]]
			types = self.type[neighbor]
			ind = (types == 0)
			support = np.append(support,neighbor[ind])
			_, index = np.unique(support, return_index=True)
			support = support[np.sort(index)]
			support = support[np.where(support != idx)[0]]	
			j = j + 1
		support[5] = idx
		support = support[0:6]
		return support
	

	def set_radius(self):
		xmax = self.xmax
		ymax = self.ymax
		xmin = self.xmin
		ymin = self.ymin
		for i in range(self.size):
			if (self.coord[i,0]+self.r[i])>xmax:
				self.r[i] = xmax - self.coord[i,0]
			if (self.coord[i,0]-self.r[i])<xmin:
				self.r[i] = self.coord[i,0] -xmin
			if (self.coord[i,1]+self.r[i])>ymax:
				self.r[i] = ymax - self.coord[i,1]
			if (self.coord[i,1]-self.r[i])<ymin:
				self.r[i] = self.coord[i,1] -ymin


def connect(mesh):
		# Init vertex-edge connectivity
	mesh.init(0,1)
	nei = []
	for v in vertices(mesh):
	    idx = v.index()
	    neighborhood = [Edge(mesh, i).entities(0) for i in v.entities(1)]
	    neighborhood = np.array(neighborhood).flatten()
	    neighborhood = neighborhood[np.where(neighborhood != idx)[0]]
	    nei.append(neighborhood)
	return nei
			