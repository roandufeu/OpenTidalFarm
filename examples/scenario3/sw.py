import configuration 
import numpy
import IPOptUtils
from reduced_functional import ReducedFunctional
from dolfin import *
from scipy.optimize import fmin_slsqp
set_log_level(ERROR)

# We set the perturbation_direction with a constant seed, so that it is consistent in a parallel environment.
numpy.random.seed(21) 

# Some domain information extracted from the geo file
basin_x = 1200
basin_y = 1000
land_x = 600
land_y = 300
land_site_delta = 100
site_x = 150
site_y = 100
site_x_start = basin_x - land_x
site_y_start = land_y + land_site_delta 
config = configuration.ScenarioConfiguration("mesh.xml", inflow_direction = [0,1])
config.set_site_dimensions(site_x_start, site_x_start + site_x, site_y_start, site_y_start + site_y)

# Place some turbines 
IPOptUtils.deploy_turbines(config, nx = 3, ny = 3) 

model = ReducedFunctional(config, scaling_factor = -10**1, plot = True)
m0 = model.initial_control()

g = lambda m: []
dg = lambda m: []

# Get the upper and lower bounds for the turbine positions
lb, ub = IPOptUtils.position_constraints(config) 
bounds = [(lb[i], ub[i]) for i in range(len(lb))]

f_ieqcons, fprime_ieqcons = IPOptUtils.get_minimum_distance_constraint_func(config)

fmin_slsqp(model.j, m0, fprime = model.dj, bounds = bounds, f_ieqcons = f_ieqcons, fprime_ieqcons = fprime_ieqcons, iprint = 2, full_output = True)
