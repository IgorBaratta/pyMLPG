from pymlpg import *
import numpy as np

# # Square 2x2
# p1 = np.array([0,0])
# p2 = np.array([2,2])
# domain = Rectangle(p1,p2)

# x = np.linspace(0.5, 1.5, 3)
# y = np.linspace(0.5, 1.5, 3)
# xv, yv = np.meshgrid(x, y)
# xv = xv.flatten()
# yv = yv.flatten()
# centers = np.vstack((xv,yv))

# for center in centers.transpose():
# 	circle = Circle(center,radius=0.1)
# 	domain.addMaterial(circle,'Empty')

# nodes = Nodes(domain,N=1000)
# nodes.plot()



domain = Circle()
nodes = Nodes(domain, N=1000, method = 'Irregular')
nodes.plot()