import numpy as np

# TODO - Use Numba to speed up 
class PolynomialBasis:
    '''
    Thies class creates a complete polynomial basis function 
    of monomials of arbitrary order 'm' in 2 dimensions
    
    Usage Example:
        P = PolynomialBasis(3)
        x = np.array([2,1])
        p = P(x)
        array([ 1.,  2.,  1.,  4.,  2.,  1.,  8.,  4.,  2.,  1.])

        d1pdx0 = P.diff(x)
        array([  0.,   1.,   0.,   4.,   1.,   0.,  12.,   4.,   1.,   0.])

    '''
    def __init__(self, m):
        self.m = m
        self.l = (m+1)*(m+2)/2

        #Initialize exponent array
        self.e = self.exponent()

    def exponent(self):
        exponent = np.zeros([self.l,2])
        monome = 0
        for k in range(self.m + 1):
            for i in range(k + 1):
                exponent[monome,:] = [k-i,i]
                monome += 1
        return exponent

    def __repr__(self):
        return "P%d(x,y)" %self.m

    def __str__(self):
        return "Polynomial Basis of order " + str(self.m)

    def __len__(self):
        '''
        Return the number of monomes in the polynomial basis P
        '''
        return self.l

    def __call__(self,x):
        '''
        Return the comple polynomial basis [P] evaluated at the point
        [x].
        attr:
            x - np.array(float) of shape(2,)
        output:
            P - np.array(float) of shape(len(self),)
        '''
        P = np.prod(np.power(x,self.e),axis=1)
        return  P
    


    def diff(self, x, axis = 0, order = 1):
        '''
        Derivative of 'order' on the direction 'axis' at
        the point 'x'

        attr:
            x - np.array(float) of shape(2,)
            axis - 0 for x[0] or (x), 1 for x[1] or (y)
            order - order of the derivative

        output:
            dP_dx -  np.array(float) of shape(len(self),)
        '''
        from scipy.misc import factorial
        eps= 9.7656e-10

        den = factorial(self.e[:,axis]-order)
        den[den<eps] = np.infty
        alpha = factorial(self.e[:,axis])/den
        
        order_array = np.zeros(2)
        order_array[axis] = order

        dP_dx = alpha*np.prod(np.power(x,self.e-order_array[axis]),axis=1)

        return dP_dx