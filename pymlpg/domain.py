import numpy as np
from matplotlib.patches import Polygon

class Domain:
    def __init__(self, vertices, subdomains = []):
        '''
        This class describes a 2D polygonal domain 
        defined by n vertices following a closed
        path:
        vertices = np.array([[x0,y0],
                             [x1,y1],
                             [x2,y2],
                              ...
                             [xn,yn]])

        '''
        self.vertices = vertices
        self.polygon = Polygon(vertices)
        self.path = self.polygon.get_path()
        self.subdomains = subdomains

    def __add__ (self,other):
        if isinstance(other, Domain):
            return Domain(self.vertices,self.subdomains.append(other))
        else:
            return NotImplemented

    def __sub__ (self, other):
        return None

    def __len__(self):
        return len(self.subdomains) + 1

    @property
    def area(self):
        ''' 
        Implementation of Shoelace formula using Numpy
        https://en.wikipedia.org/wiki/Shoelace_formula
        '''
        x = self.vertices[:,0]
        y = self.vertices[:,1]
        return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

    def contains(self, points):
        """
        Return True for points inside Domain
        """
        return np.array(self.path.contains_points(points))

    def contains_all(self,points):
        """
        Return True if all points are inside Domain
        """
        return self.contains(points).all()

    def bounds(self):
        '''
        Returns a (minx, miny, maxx, maxy) tuple (float values) that bounds the object.
        '''
        return np.array([np.min(self.vertices,axis=0), np.max(self.vertices,axis=0)])

    def addMaterial(self, Subdomain,Value):
        if self.contains_all(Subdomain.vertices):
            Subdomain.Value = Value
            self.subdomains.append(Subdomain)
        else:
            raise TypeError

    def plot(self, show=True):
        from matplotlib import pylab
        from matplotlib.lines import Line2D
        color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        fig, ax = pylab.subplots()
        self.polygon.set_alpha(0.3)
        x, y = zip(*self.polygon.xy)
        
        line = Line2D(x, y, linestyle= '-',
                            marker='o', 
                            markerfacecolor='k')
        ax.add_line(line)
        ax.set_title('Computational Domain')
        ax.set_xlim((min(x)-0.1*max(x), 1.1*max(x)))
        ax.set_ylim((min(y)-0.1*max(y), 1.1*max(y)))
        line.set_color('k')
        ax.add_patch(self.polygon)
        for i, sub in enumerate(self.subdomains):
            x, y = zip(*sub.polygon.xy)
            line = Line2D(x, y, linestyle= '--',
                            color = 'k')
            ax.add_line(line)
            sub.polygon.set_color(color[1])
            sub.polygon.set_alpha(0.3)
            ax.add_patch(sub.polygon)       
        if show:
            pylab.show()


class Rectangle(Domain):
    def __init__(self, Point1=[0, 0], Point2=[1, 1]):
        
        '''
        This class describes a 2D rectangle domain 
        defined by two opposite corners points:
        lower left  -    Point1 = np.array(x0, y0) 
        upper right -    Point2 = np.array(x1, y1).
        By the default, this class creates a unit
        square domain.
        
        '''
        vertices = self.rectangle_vertices(Point1, Point2)
        Domain.__init__(self, vertices)

    def rectangle_vertices(self, Point1, Point2):
        vertices = np.zeros((4,2))
        vertices[0] = Point1
        vertices[1] = [Point2[0],Point1[1]]
        vertices[2] = Point2
        vertices[3] = [Point1[0],Point2[1]]
        return vertices


class Circle(Domain):
    def __init__(self, Center=np.array([0,0]), radius=1.0, resolution=100):
        '''
        This class describes a 2D circular domain 
        defined by center Point and the radius:
        Center Point  -   Center = np.array(x, y) 
        radius        -   float(radius)

        The circle is approximated by a regular polygon 
        with resolution sides.

        By default a circle centered at the origin with unitary
        radius is created
        '''
        self.radius = radius
        self.center = Center
        vertices = [[radius*np.sin(x)+Center[0],radius*np.cos(x)+Center[1]] 
                    for x in np.linspace(0,2*np.pi,resolution)[:-1]]
        Domain.__init__(self, vertices)

    @property
    def area(self):
        self.area = np.pi*(self.radius**2)
        return self.area

