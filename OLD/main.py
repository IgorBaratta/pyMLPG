from mshr import *
from dolfin import *
from nodes import *
from scipy.sparse import dok_matrix
from ddm import *
from matplotlib import pyplot as plt
import matplotlib.tri as tri
from scipy.sparse.linalg import spsolve, gmres, splu, LinearOperator, inv
from scipy.linalg import expm_cond

parameters['reorder_dofs_serial'] = False
geo = 'geo12'
mesh = Mesh('Geometria/geo12.xml');
mt = mesh.num_vertices()
coord = mesh.coordinates()
nei = connect(mesh)
dmed = (mesh.hmax() + mesh.hmin())/2
Xn = nodes(coord)
Xn.set_neighborhood(nei)

# Number of subdomains
nd = 10 
K =  dok_matrix((Xn.size, Xn.size), dtype=np.complex)
k = 2*np.pi
f = np.zeros(Xn.size,dtype=np.complex)

A, f = assemble(Xn,k)
u_d = spsolve(A,f)

triang = mesh2triang(mesh)
plt.gca().set_aspect('equal')
plt.tripcolor(triang, np.real(u_d),shading='gouraud')
plt.show()

dx = nd
ovl = 0.1
list_mesh = submeshes(mesh,nd,dx,ovl,Verbose=False)
r,rd = indices(list_mesh, mesh)

R = [];D = [];K = []; Kinv = []
for j in range(nd):
	Ri, Di = restriction(len(r[j]),mt,r[j],rd)
	submesh = list_mesh[j]
	nei = connect(submesh)
	Xj = nodes(submesh.coordinates())
	Xj.set_neighborhood(nei)
	Kj =  dok_matrix((Xj.size, Xj.size), dtype=np.complex)
	
	# Assemble submatrix
	for I in range(Xj.size):
		support = Xj.coord[Xj.get_support(I),:]
		Phi, dPhix, dPhiy, dPhix2, dPhiy2  = shape_function(Xj.coord[I],support)
		if Xj.type[I] == 1: 	# Left Boundary
			n = 1
			if j >= 1:
				Kj[I,Xj.get_support(I)] = n*dPhix + 1j*k*Phi - 1j/(2*k)*dPhiy2
			else:
				Kj[I,Xj.get_support(I)] = n*dPhix + 1j*k* Phi
		elif Xj.type[I] == 2: 	# Right Boundary
			n = -1
			if j < nd:
				Kj[I,Xj.get_support(I)] = n*dPhix + 1j*k* Phi - 1j/(2*k)*dPhiy2
			else:
				Kj[I,Xj.get_support(I)] = n*dPhix + 1j*k* Phi
		elif Xj.type[I] == 3: 	# Bottom Boundary
			n = 1
			Kj[I,Xj.get_support(I)] = n*dPhiy + 1j*k*Phi
		elif Xj.type[I] == 4:	# Top Boundary
			n = -1
			Kj[I,Xj.get_support(I)] = n*dPhiy + 1j*k*Phi
		else: 	# Internal Nodes
			Kj[I,Xj.get_support(I)] = dPhix2 + dPhiy2 +(k**2)*Phi
	
	R.append(Ri)
	D.append(Di)
	K.append(Kj)
	Kinv.append(splu(Kj.tocsc()))

maxiter = 1000
tol = 1e-5
b = rhs


def ddm_operator(r1):
	u = 1j*np.zeros(A.shape[0])
	residual = r1 - A*u
	for i in range(nd):
		v1 = Kinv[i].solve(R[i]*residual)
		u = u + R[i].transpose()*(D[i]*v1)
	return u

M_x = ddm_operator
M1 = LinearOperator(A.shape, M_x)

counter_ddm = Counter_Iter()
M_oras = sum(R[i].transpose()*D[i]*inv(K[i])*R[i] for i in range(nd))

usol_gmres, info = gmres(A,f,M = M_oras,restart = 2000, maxiter=maxiter, \
							callback=counter_ddm, tol=tol)

rn = np.array(counter_ddm.rk)
np.save(geo, rn)
plt.semilogy(rn/max(rn))


counter_ddm = Counter_Iter()

usol_gmres, info = gmres(A,f,restart = 100, maxiter=maxiter, \
							callback=counter_ddm, tol=tol)
