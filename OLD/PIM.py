import numpy as np
from numpy.linalg import inv, norm, matrix_rank
import itertools


p = lambda x, y: np.array([1.0, x, y, x*y, x**2, y**2])

dpx  = lambda x, y: np.array([0, 1, 0, y, 2*x, 0])
dpx2 = lambda x, y: np.array([0, 0, 0, 0, 2, 0])
dpy = lambda x, y: np.array([0, 0, 1, x, 0, 2*y])
dpy2 = lambda x, y: np.array([0, 0, 0, 0, 0, 2])

def moment_matrix(nodes):
	N = nodes.size/2
	Pq = np.zeros([N,N])
	for (point,i) in zip(nodes,range(N)):
		Pq[i,:] = p(point[0],point[1])
	return Pq

def shape_function(x,support):
	Pq = moment_matrix(support)
	p_x = p(x[0],x[1])
	dp_x = dpx(x[0],x[1])
	dp_y = dpy(x[0],x[1])
	dp_x2 = dpx2(x[0],x[1])
	dp_y2 = dpy2(x[0],x[1])
	Pq_inv =  inv(Pq)
	Phi = np.dot(p_x,Pq_inv)
	dPhix = np.dot(dp_x,Pq_inv)
	dPhiy = np.dot(dp_y,Pq_inv)
	dPhix2 = np.dot(dp_x2,Pq_inv)
	dPhiy2 = np.dot(dp_y2,Pq_inv)

	return Phi, dPhix, dPhiy, dPhix2, dPhiy2