from opentidalfarm import *
set_log_level(INFO)

# build the domain and the location of the turbine site - info should be in the
# geo file
domain = FileDomain('mesh/mesh.xml')

site_x_start = 160
site_y_start = 80
site_x_end = 480
site_y_end = 240
basin_y = 320

# determine the eco domain
eco_x_start = 500
eco_y_start = 180
eco_x_end = 640
eco_y_end = 320

# set the boundary conditions
bcs = BoundaryConditionSet()
bcs.add_bc("u", Constant((2, 0)), facet_id=1)
bcs.add_bc("eta", Constant(0), facet_id=2)
bcs.add_bc("u", Constant((0, 0)), facet_id=3, bctype="weak_dirichlet")

# set parameters
prob_params = SteadySWProblem.default_parameters()
prob_params.domain = domain
prob_params.bcs = bcs
# prob_params.viscosity = Constant(60)
# prob_params.friction = Constant(0.0025)
prob_params.eco_penalty_multiplier = 1
prob_params.sed_penalty_multiplier = 200
prob_params.ecology = True
prob_params.power = True
prob_params.sediment = True

# Create a farm with no effect in order to determine the ambient flow
turbine_amb = BumpTurbine(diameter=20.0, friction=0.0)
farm_amb = RectangularFarm(domain, site_x_start, site_x_end, site_y_start,
                           site_y_end, turbine=turbine_amb)
farm_amb.add_regular_turbine_layout(num_x=1, num_y=1)
prob_params.tidal_farm = farm_amb

# create amb problem & solver
prob_params.tidal_farm.update()
problem_amb = SteadySWProblem(prob_params)
sol_params = CoupledSWSolver.default_parameters()
sol_params.dump_period = 1
# sol_params.output_turbine_power = True
solver_amb = CoupledSWSolver(problem_amb, sol_params)

# solve for ambient flow
if prob_params.ecology or prob_params.sediment:
    import os.path
    if not os.path.exists('ambient.pvd'):
        # Initialise the solver
        solver = solver_amb.solve()
        state = solver.next()

        # Solve the current state of the problem
        state = solver.next()

        # save current state of problem
        u_ambient_state = state['state']

        File('ambient.xml') << u_ambient_state
        File('ambient.pvd') << u_ambient_state

        import opentidalfarm as otf
        V, H = otf.finite_elements.p2p1(domain.mesh)
        func_space = MixedFunctionSpace([V, H])
        u_ambient = Function(func_space, u_ambient_state)
    else:
        import opentidalfarm as otf
        V, H = otf.finite_elements.p2p1(domain.mesh)
        func_space = MixedFunctionSpace([V, H])
        u_ambient = Function(func_space, 'ambient.xml')

# create eco variables
if prob_params.ecology or prob_params.sediment:
    # create eco domain
    temp = FunctionSpace(domain.mesh, 'CG', 1)
    eco_fn = Eco_function_square(eco_x_start, eco_y_start, eco_x_end, eco_y_end)
    unit_function = interpolate(eco_fn, temp)
    prob_params.eco_domain = unit_function
    # create an eco penatly multiplier that takes into account the area and flow
    # speed over the eco domain with the aim of balancing the equations and
    # getting a multiplier with greater meaning
    ua_x = u_ambient[0]
    ua_y = u_ambient[1]
    Ua = (ua_x**2 + ua_y**2)**0.5
    flow = assemble(Ua*unit_function*prob_params.domain.dx)
    eco_mult = 55000000/flow
    prob_params.multiplier = eco_mult

# build the turbine farm
turbine = BumpTurbine(diameter=20.0, friction=21.0)
farm = RectangularFarm(domain, site_x_start, site_x_end, site_y_start,
                       site_y_end, turbine=turbine)
farm.add_regular_turbine_layout(num_x=8, num_y=4)
prob_params.tidal_farm = farm
prob_params.tidal_farm.update()

# set up problem and solver
problem = SteadySWProblem(prob_params)
solver = CoupledSWSolver(problem, sol_params)

# build functional
functional = PowerFunctional(problem) \
    + SedimentFunctional(problem, u_ambient)

# build reduced functional
control = TurbineFarmControl(farm)
rf_params = ReducedFunctional.default_parameters()
rf_params.automatic_scaling = 5
rf = ReducedFunctional(functional, control, solver, rf_params)

print rf_params

# optimise
lb, ub = farm.site_boundary_constraints()
f_opt = maximize(rf, bounds=[lb, ub], method="L-BFGS-B",
                 options={'maxiter': 150})
