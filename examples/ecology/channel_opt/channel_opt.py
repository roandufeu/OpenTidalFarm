from opentidalfarm import *
set_log_level(INFO)

# build the domain and the location of the turbine site - info should be in the
# geo file
domain = FileDomain('mesh/mesh.xml')

site_x_start = 160
site_y_start = 80
site_x_end = 480
site_y_end = 240


class FarmDomain(SubDomain):
    def inside(self, x, on_boundary):
        return (site_x_start <= x[0] <= site_x_end and
                site_y_start <= x[1] <= site_y_end)
farm_domain = FarmDomain()

# determine the eco domain
eco_x_start = 500
eco_y_start = 180
eco_x_end = 640
eco_y_end = 320


class EcoDomain(SubDomain):
    def inside(self, x, on_boundary):
        return (eco_x_start <= x[0] <= eco_x_end and
                eco_y_start <= x[1] <= eco_y_end)
eco_domain = EcoDomain()

# set the boundary conditions
bcs = BoundaryConditionSet()
bcs.add_bc("u", Constant((2, 0)), facet_id=1)
bcs.add_bc("eta", Constant(0), facet_id=2)
bcs.add_bc("u", Constant((0, 0)), facet_id=3, bctype="weak_dirichlet")

# set parameters
prob_params = SteadySWProblem.default_parameters()
prob_params.domain = domain
prob_params.bcs = bcs
prob_params.eco_penalty_multiplier = 1
prob_params.sed_penalty_multiplier = 200

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
solver_amb = CoupledSWSolver(problem_amb, sol_params)

# solve for ambient flow
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
eco_mult = 1000
prob_params.multiplier = eco_mult

# build the turbine farm
turbine = BumpTurbine(diameter=20.0, friction=21.0)
V = FunctionSpace(domain.mesh, "DG", 0)
farm = Farm(domain, turbine, function_space=V)
farm._regular_turbine_layout(8, 4, site_x_start, site_x_end, site_y_start,
                             site_y_end)

# set up domains
domains = MeshFunction("size_t", domain.mesh, domain.mesh.topology().dim())
domains.set_all(0)
farm_domain.mark(domains, 1)
eco_domain.mark(domains, 2)
domain_dx = Measure("dx")[domains]
farm.site_dx = domain_dx(1)
prob_params.eco_dx = domain_dx(2)
prob_params.sed_dx = domain_dx(2)

# update tidal farm
prob_params.tidal_farm = farm
prob_params.tidal_farm.update()

# set up problem and solver
problem = SteadySWProblem(prob_params)
solver = CoupledSWSolver(problem, sol_params)

# build functional
functional = PowerFunctional(problem) \
    - SedimentFunctional(problem, u_ambient) \
    - EcologyFunctional(problem, u_ambient)

# build reduced functional
control = TurbineFarmControl(farm)
rf_params = ReducedFunctional.default_parameters()
rf_params.automatic_scaling = 5
rf = ReducedFunctional(functional, control, solver, rf_params)

print rf_params

# optimise
n_turbines = len(farm.turbine_positions)
radius = farm._turbine_specification.radius
lb = n_turbines*[Constant(site_x_start + radius),
                 Constant(site_y_start + radius)]
ub = n_turbines*[Constant(site_x_end - radius),
                 Constant(site_y_end - radius)]

f_opt = maximize(rf, bounds=[lb, ub], method="L-BFGS-B",
                 options={'maxiter': 150})
