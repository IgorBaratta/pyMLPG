from mshr import *
from dolfin import *
import matplotlib.tri as tri
import numpy as np
import random
from matplotlib import pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse import dok_matrix
from nodes import *
from PIM import *
from mesh import *
from matplotlib import pyplot as plt
import matplotlib.tri as tri
from functools import reduce

parameters['reorder_dofs_serial'] = False

def mesh2triang(mesh):
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())


def pairwise_distance(Xn):
	n = Xn.size/2
	r = np.zeros([n,n])
	for (x,i) in zip(Xn,range(n)):
		for (y,j) in zip(Xn,range(n)):
			r[i,j] = norm(x-y)
	return r


def submeshes(mesh,nd,d_x, ovr, Verbose=True):
	limites = np.linspace(0, d_x, num=nd + 1)
	domains = MeshFunction("size_t", mesh, mesh.topology().dim())
	domains.set_all(50)
	submeshes = []
	r = lambda x: random.randint(0,x)
	for i in range(nd):
		d1 = Domain(limites[i]-ovr,limites[i+1]+ovr)
		d1.mark(domains,i)
		submesh = SubMesh(mesh, domains, i)
		if Verbose == True:
			color = '#%02X%02X%02X' % (r(10),r(100),r(100))
			plt.gca().set_aspect('equal')
			plt.triplot(mesh2triang(submesh), color = color, alpha	= 0.5)
		submeshes.append(submesh)
	if Verbose == True:
		plt.show()
	return submeshes


class Domain(SubDomain):
    def __init__(self,xmin,xmax):
        super(Domain, self).__init__()
        self.xmin = xmin
        self.xmax = xmax
    def inside(self, x, on_boundary):
        return x[0]<=self.xmax and x[0]>=self.xmin


def indices(list_mesh, mesh):
	r = []
	nd = len(list_mesh)
	coord = mesh.coordinates()
	for i in range(nd):
		subi = list_mesh[i]
		nv = subi.num_vertices() 
		ri = np.zeros(nv,dtype=np.int)
		for ind, point in enumerate(subi.coordinates()):
			in_x = point[0]==coord[:,0]
			in_y = point[1]==coord[:,1]
			ri[ind] = np.where(in_x* in_y)[0][0]
		r.append(ri)
	rd = []
	for i in range(nd-1):
		rd.append(np.intersect1d(r[i],r[i+1]))
	rd = reduce(np.union1d,rd)
	return r, rd


# Restriction Operator
def restriction(m1,mt,ind,ind_overlap):
	row =[]
	col =[]
	cold = []
	data = []
	datad = []
	rowd = []
	for i in xrange(m1):
		row.append(i)
		col.append(ind[i])
		data.append(1.0)
		if ind[i] in ind_overlap:
			rowd.append(i)
			cold.append(i)
			datad.append(0.5)
		else:
			rowd.append(i)
			cold.append(i)
			datad.append(1.0)
	R1 = csr_matrix((data,(row,col)),shape=(m1,mt))
	D1 = csr_matrix((datad,(rowd,cold)),shape=(m1,m1))
	
	return R1, D1

def interface(x, bcoord):
	eps = 0.01
	for y in bcoord:
		if all(x-y<eps):
			return True
	return False

# def assemble(Xn,k):
# 	f = np.zeros(Xn.size,dtype=np.complex)
# 	K =  dok_matrix((Xn.size, Xn.size), dtype=np.complex)

# 	for I in range(Xn.size):
# 		support = Xn.coord[Xn.get_support(I),:]
# 		Phi, dPhix, dPhiy, dPhix2, dPhiy2  = shape_function(Xn.coord[I],support)
# 		x = Xn.coord[I,0]
# 		y = Xn.coord[I,1]
# 		if Xn.type[I] >= 1: # Boundary Nodes
# 			if x == 0 or x == 5:
# 				K[I,Xn.get_support(I)] = ((-1)**x)*dPhix + 1j*k* Phi
# 				f[I] = ((-1)**x)*1j*k*np.exp(-1j*k*x) - 1j*k*np.exp(-1j*k*x)
# 			else:
# 				K[I,Xn.get_support(I)] = ((-1)**y)*dPhiy + 1j*k*Phi
# 				f[I] = -1j*k*np.exp(-1j*k*x)		
# 		else: # Internal Nodes
# 			K[I,Xn.get_support(I)] = dPhix2 + dPhiy2 +(k**2)*Phi
# 	K = K.tocsr()
# 	return K, f


class Counter_Iter(object):
    def __init__(self):
        self.niter = 0
        self.rk=[]
    def __call__(self, rk=None, x=None):
        self.niter += 1
        self.rk.append(rk)



def assemble(Xn,k):
	f = np.zeros(Xn.size,dtype=np.complex)
	K =  dok_matrix((Xn.size, Xn.size), dtype=np.complex)

	for I in range(Xn.size):
		support = Xn.coord[Xn.get_support(I),:]
		Phi, dPhix, dPhiy, dPhix2, dPhiy2  = shape_function(Xn.coord[I],support)
		x = Xn.coord[I,0]
		y = Xn.coord[I,1]
		if Xn.type[I] >= 1: # Boundary Nodes
			if Xn.type[I]==1:
				n = 1
				K[I,Xn.get_support(I)] = n*dPhix + 1j*k* Phi
				f[I] = n*1j*k*np.exp(-1j*k*x) - 1j*k*np.exp(-1j*k*x)
			elif Xn.type[I]==2:
				n = -1
				K[I,Xn.get_support(I)] = n*dPhix + 1j*k* Phi
				f[I] = n*1j*k*np.exp(-1j*k*x) - 1j*k*np.exp(-1j*k*x)
			elif Xn.type[I]==3:
				n = 1
				K[I,Xn.get_support(I)] = n*dPhiy + 1j*k*Phi
				f[I] = -1j*k*np.exp(-1j*k*x)
			else:
				n = -1
				K[I,Xn.get_support(I)] = n*dPhiy + 1j*k*Phi
				f[I] = -1j*k*np.exp(-1j*k*x)		
		else: # Internal Nodes
			K[I,Xn.get_support(I)] = dPhix2 + dPhiy2 +(k**2)*Phi
	K = K.tocsr()
	return K, f