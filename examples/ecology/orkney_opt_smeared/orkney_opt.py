from opentidalfarm import *
from model_turbine import ModelTurbine
set_log_level(INFO)

model_turbine = ModelTurbine()
print model_turbine

# build the domain
domain = FileDomain('mesh/earth_orkney_converted_coarse.xml')

prob_params = SteadySWProblem.default_parameters()
prob_params.domain = domain

# determine the location of the farm
site_x_start = 10307240.
site_y_start = 6522134.31356
site_x_end = site_x_start + 1100
site_y_end = site_y_start + 700


class FarmDomain(SubDomain):
    def inside(self, x, on_boundary):
        return (site_x_start <= x[0] <= site_x_end and
                site_y_start <= x[1] <= site_y_end)
farm_domain = FarmDomain()

# determine the location of the eco domain
eco_x_start = site_x_end + 200
eco_y_start = site_y_start + 300
eco_x_end = eco_x_start + 300
eco_y_end = eco_y_start + 300


class EcoDomain(SubDomain):
    def inside(self, x, on_boundary):
        return (eco_x_start <= x[0] <= eco_x_end and
                eco_y_start <= x[1] <= eco_y_end)
eco_domain = EcoDomain()

# create the farm and the domains
turbine = SmearedTurbine()
V = FunctionSpace(domain.mesh, "DG", 0)
farm = Farm(domain, turbine, function_space=V)

domains = MeshFunction("size_t", domain.mesh, domain.mesh.topology().dim())
domains.set_all(0)
farm_domain.mark(domains, 1)
eco_domain.mark(domains, 2)
site_dx = Measure("dx")[domains]
farm.site_dx = site_dx(1)
prob_params.eco_dx = site_dx(2)
prob_params.sed_dx = site_dx(2)

prob_params.tidal_farm = farm

# set the boundary conditions
bcs = BoundaryConditionSet()
bcs.add_bc("u", Constant((2, 0)), facet_id=1)
bcs.add_bc("eta", Constant(0), facet_id=2)
bcs.add_bc("u", Constant((0, 0)), facet_id=3, bctype="weak_dirichlet")

# set parameters
prob_params.bcs = bcs
prob_params.viscosity = Constant(60)
prob_params.friction = Constant(0.0025)
prob_params.eco_penalty_multiplier = 4
prob_params.sed_penalty_multiplier = 0.0001
prob_params.initial_condition = Expression(("1e-7", "0", "0"))
prob_params.multiplier = 4000

# set up ambient problem
problem_amb = SteadySWProblem(prob_params)
sol_params = CoupledSWSolver.default_parameters()
sol_params.dump_period = 1
sol_params.cache_forward_state = False
solver_amb = CoupledSWSolver(problem_amb, sol_params)

# solve for ambient flow
import os.path
if not os.path.exists('ambient.pvd'):
    # Initialise the solver
    solver2 = solver_amb.solve()
    state = solver2.next()

    # Solve the current state of the problem
    state = solver2.next()

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

# set up initial farm friction


class InitialFriction(Expression):
    def eval(self, values, x):
        if x[0] > site_x_start and x[1] > site_y_start and \
                x[0] < site_x_end and x[1] < site_y_end:
            values[:] = Constant(model_turbine.maximum_smeared_friction)
        else:
            values[:] = Constant(0)

initial_friction = InitialFriction()
prob_params.tidal_farm.friction_function.assign(initial_friction)

# build final problem and solver
problem = SteadySWProblem(prob_params)
solver = CoupledSWSolver(problem, sol_params)

# build functional
functional = PowerFunctional(problem) - EcologyFunctional(problem, u_ambient)
#    - SedimentFunctional(problem, u_ambient)

# build reduced functional
control = TurbineFarmControl(farm)
rf_params = ReducedFunctional.default_parameters()
rf_params.automatic_scaling = None
rf = ReducedFunctional(functional, control, solver, rf_params)

print rf_params

# optimise
f_opt = maximize(rf, bounds=[0, model_turbine.maximum_smeared_friction],
                 method="L-BFGS-B", options={'maxiter': 30})
