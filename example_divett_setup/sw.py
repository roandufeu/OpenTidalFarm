import configuration 
import numpy
from dirichlet_bc import DirichletBCSet
from reduced_functional import ReducedFunctional
from dolfin import *
set_log_level(PROGRESS)

config = configuration.ConstantInflowPeriodicSidesPaperConfiguration(nx = 100, ny = 33, basin_x = 200, basin_y = 66)

# The turbine position is the control variable 
#config.params["turbine_pos"] = [[100, 33]] 
border_x = 80.
border_y = 20.
for x_r in numpy.linspace(0.+border_x, config.params["basin_x"]-border_x, 2):
    for y_r in numpy.linspace(0.+border_y, config.params["basin_y"]-border_y, 2):
      config.params["turbine_pos"].append((float(x_r), float(y_r)))
      

info_blue("Deployed " + str(len(config.params["turbine_pos"])) + " turbines.")
# Choosing a friction coefficient of > 0.25 ensures that overlapping turbines will lead to
# less power output.
config.params["turbine_friction"] = numpy.ones(len(config.params["turbine_pos"]))

model = ReducedFunctional(config)

m0 = model.initial_control()
j0 = model.j(m0, forward_only = True)
info("Power outcome: %f" % (j0, ))
info("Timing summary:")
timer = Timer("NULL")
timer.stop()

list_timings()

