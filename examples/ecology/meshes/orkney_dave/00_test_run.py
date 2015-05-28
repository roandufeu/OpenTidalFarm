from opentidalfarm import *

domain = FileDomain("mesh/earth_orkney_converted_coarse.xml")

site_x_start = 10307240. 
site_x = 1100
site_y_start = 6522134.31356 
site_y = 700

# Specify boundary conditions.
bcs = BoundaryConditionSet()
bcs.add_bc("u", Constant((2, 0)), facet_id=1)
bcs.add_bc("eta", Constant(0), facet_id=2)
# The free-slip boundary conditions.
bcs.add_bc("u", Constant((0, 0)), facet_id=3, bctype="weak_dirichlet")

# Set the shallow water parameters
prob_params = SteadySWProblem.default_parameters()
prob_params.domain = domain
prob_params.bcs = bcs
prob_params.viscosity = Constant(60) 
prob_params.depth = Constant(50)
prob_params.friction = Constant(0.0025)

# Set up the turbines
turbine = BumpTurbine(diameter=20.0, friction=12.0)

# And the farm
farm = RectangularFarm(domain, 
                       site_x_start=site_x_start, 
                       site_x_end=site_x_start+site_x,
                       site_y_start=site_y_start, 
                       site_y_end=site_y_start+site_y, 
                       turbine=turbine)
farm.add_regular_turbine_layout(num_x=8, num_y=4)
prob_params.tidal_farm = farm

# Set up the problem
problem = SteadySWProblem(prob_params)

# Set up the solver parameters
sol_params = CoupledSWSolver.default_parameters()
sol_params.dump_period = 1
solver = CoupledSWSolver(problem, sol_params)

# Create the functional and reduced functional
functional = PowerFunctional(problem)
control = TurbineFarmControl(farm)
rf_params = ReducedFunctional.default_parameters()
rf_params.automatic_scaling = 5
rf = ReducedFunctional(functional, control, solver, rf_params)

print rf_params

lb, ub = farm.site_boundary_constraints()
f_opt = maximize(rf, 
                 bounds=[lb, ub], 
                 method="SLSQP", 
                 options={'maxiter': 100})



