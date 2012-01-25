from dolfin import * 
from math import exp, sqrt, pi

import sw_lib

params=sw_lib.parameters({
    'depth' : 2.,
    'g' : 9.81,
    'f' : 0.0,
    'dump_period' : 1,
    'eta0' : 2 # Wave height
    })

# Basin radius.
basin_x=3000 # The length of the basin
basin_y=1000 # The width of the basin
# Long wave celerity.
c=sqrt(params["g"]*params["depth"])

params["finish_time"]=100
params["dt"]=params["finish_time"]/4000.
params["k"]=pi/basin_x 

class InitialConditions(Expression):
    def __init__(self):
        pass
    def eval(self, values, X):
        values[0]=params['eta0']*sqrt(params['g']*params['depth'])*cos(params["k"]*X[0])
        values[1]=0.
        values[2]=params['eta0']*cos(params["k"]*X[0])
    def value_shape(self):
        return (3,)


def generate_mesh(nx=20, ny=3):
  ''' Generates a rectangular mesh for the divett test
      nx = Number of cells in x direction
      ny = Number of cells in y direction  '''
  mesh = Rectangle(0, 0, basin_x, basin_y, nx, ny)
  mesh.order()
  mesh.init()
  return mesh

class Left(SubDomain):
      def inside(self, x, on_boundary):
           return on_boundary and near(x[0], 0.0)

class Right(SubDomain):
      def inside(self, x, on_boundary):
           return on_boundary and near(x[0], basin_x)

class Sides(SubDomain):
      def inside(self, x, on_boundary):
           return on_boundary and (near(x[1], 0.0) or near(x[1], basin_y))

# Initialize sub-domain instances
left = Left()
right = Right()
sides = Sides()
mesh = generate_mesh()

# Initialize mesh function for boundary domains
boundaries = FacetFunction("uint", mesh)
boundaries.set_all(0)
left.mark(boundaries, 1)
right.mark(boundaries, 2)
sides.mark(boundaries, 3)
ds = Measure("ds")[boundaries]