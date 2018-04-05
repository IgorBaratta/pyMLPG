import quadpy
import numpy
from PIM import *

def kahan_sum(a, axis=0):
    '''Kahan summation of the numpy array a along axis k.
    '''
    # See <https://en.wikipedia.org/wiki/Kahan_summation_algorithm> for
    # details.
    k = axis % len(a.shape)
    s = numpy.zeros(a.shape[:axis] + a.shape[k+1:])
    c = numpy.zeros(s.shape)
    for i in range(a.shape[axis]):
        # http://stackoverflow.com/a/42817610/353337
        y = a[(slice(None),) * k + (i,)] - c
        t = s + y
        c = (t - s) - y
        s = t.copy()
    return s

# def integrate(f, center, radius, rule, sumfun=kahan_sum):
#     center = numpy.array(center)
#     rr = numpy.multiply.outer(radius, rule.points)
#     rr = numpy.swapaxes(rr, 0, -2)
#     ff = numpy.array(f((rr + center).T))
#     a = numpy.outer(ff,rule.weights)
#     out = sumfun(a, axis=-1)
#     return numpy.array(radius)**2 * out



def integrate(center, radius, support):
	# int (u dOmega)
	rule = quadpy.disk.Lether(3)
	center = numpy.array(center)
	rr = numpy.multiply.outer(radius, rule.points)
	rr = numpy.swapaxes(rr, 0, -2)
	n = support.shape[0]
	a = numpy.zeros([rr.shape[0],n])
	for i , r in enumerate(rr):
		x = numpy.array(r + center)
		phi, _ , _= shape_function(x,support)
		a[i,:] = phi*rule.weights[i]
	Omega = numpy.array(radius)**2 * kahan_sum(a)

	# int (u*n dGamma)
	rule = quadpy.circle.Equidistant(5)
	rr = numpy.multiply.outer(radius, rule.points)
	dphix = numpy.zeros([rr.shape[0],n])
	dphiy = numpy.zeros([rr.shape[0],n])
	for i , r in enumerate(rr):
		x = numpy.array(r + center)
		phi, dphix , dphiy = shape_function(x,support)
		dphidn = -dphix*(r[0]/radius) - dphiy*(r[1]/radius)
		a[i,:] = dphidn*rule.weights[i]
	Gamma = radius*kahan_sum(a)
	
	return Omega, Gamma
