from opentidalfarm import *
from model_turbine import ModelTurbine
set_log_level(INFO)


model_turbine = ModelTurbine()
print model_turbine

# build the domain and the turbine site
domain = FileDomain('mesh/earth_orkney_converted_coarse.xml')

prob_params = SteadySWProblem.default_parameters()
prob_params.domain = domain

turbine = SmearedTurbine()
V = FunctionSpace(domain.mesh, "DG", 0)
farm = Farm(domain, turbine, function_space=V)

site_x_start = 10307240.
site_y_start = 6522134.31356
site_x_end = site_x_start + 1100
site_y_end = site_y_start + 700


class FarmDomain(SubDomain):
    def inside(self, x, on_boundary):
        return (site_x_start <= x[0] <= site_x_end and
                site_y_start <= x[1] <= site_y_end)

farm_domain = FarmDomain()

domains = MeshFunction("size_t", domain.mesh, domain.mesh.topology().dim())
domains.set_all(0)
farm_domain.mark(domains, 1)
site_dx = Measure("dx")[domains]
farm.site_dx = site_dx(1)

prob_params.tidal_farm = farm

# set the boundary conditions
bcs = BoundaryConditionSet()
bcs.add_bc("u", Constant((2, 0)), facet_id=1)
bcs.add_bc("eta", Constant(0), facet_id=2)
bcs.add_bc("u", Constant((0, 0)), facet_id=3, bctype="weak_dirichlet")

# set parameters and build problem
prob_params.bcs = bcs
prob_params.viscosity = Constant(60)
prob_params.friction = Constant(0.0025)
prob_params.initial_condition = Expression(("1e-7", "0", "0"))
problem = SteadySWProblem(prob_params)

# set up problem and solver
sol_params = CoupledSWSolver.default_parameters()
sol_params.dump_period = 1
sol_params.cache_forward_state = False
solver = CoupledSWSolver(problem, sol_params)

# build functional
functional = PowerFunctional(problem)

# build reduced functional
control = TurbineFarmControl(farm)
rf_params = ReducedFunctional.default_parameters()
rf_params.automatic_scaling = False
rf = ReducedFunctional(functional, control, solver, rf_params)

print rf_params

# optimise
f_opt = maximize(rf, bounds=[0, model_turbine.maximum_smeared_friction],
                 method="L-BFGS-B", options={'maxiter': 150})
