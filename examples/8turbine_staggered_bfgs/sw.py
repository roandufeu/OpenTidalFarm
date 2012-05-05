''' This example optimises the position of three turbines using the hallow water model. '''

import sys
import configuration 
import numpy
import IPOptUtils
from dirichlet_bc import DirichletBCSet
from helpers import test_gradient_array
from animated_plot import *
from reduced_functional import ReducedFunctional
from dolfin import *
from scipy.optimize import fmin_l_bfgs_b
set_log_level(ERROR)

# An animated plot to visualise the development of the functional value
plot = AnimatedPlot(xlabel='Iteration', ylabel='Functional value')

# We set the perturbation_direction with a constant seed, so that it is consistent in a parallel environment.
numpy.random.seed(21) 
config = configuration.ConstantInflowPeriodicSidesPaperConfiguration(nx=100, ny=33)

# The turbine position is the control variable 
turbine_pos = [[30, 38], [50, 28], [70, 38], [90, 28], [110, 38], [130, 28], [150, 38], [170, 28]] 

config.set_turbine_pos(turbine_pos)
info_blue("Deployed " + str(len(turbine_pos)) + " turbines.")

model = ReducedFunctional(config, scaling_factor = -10**-6, plot = True)
m0 = model.initial_control()

g = lambda m: []
dg = lambda m: []

# Get the upper and lower bounds for the turbine positions
lb, ub = IPOptUtils.position_constraints(config)
bounds = [(lb[i], ub[i]) for i in range(len(lb))]
fmin_l_bfgs_b(model.j, m0, fprime = model.dj, bounds = bounds, iprint = 2)