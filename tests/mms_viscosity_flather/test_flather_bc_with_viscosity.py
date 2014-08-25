import math
from opentidalfarm import *
from opentidalfarm.initial_conditions import SinusoidalInitialCondition
from dolfin_adjoint import adj_reset
from dolfin import log, INFO, ERROR


class TestFlatherBoundaryConditionsWithViscosity(object):

    def error(self, problem, config, eta0, k):
        state = Function(config.function_space)
        ic_expr = SinusoidalInitialCondition(eta0, k,
                                             problem.parameters.depth, 
                                             problem.parameters.start_time)

        ic = project(ic_expr, state.function_space())
        state.assign(ic, annotate=False)
        # The analytical veclocity of the shallow water equations has been
        # multiplied by depth to account for the change of variable (\tilde u =
        # depth u) in this code.
        u_exact = "eta0*sqrt(g/depth) * cos(k*x[0]-sqrt(g*depth)*k*t)"
        ddu_exact = "(viscosity * eta0*sqrt(g/depth) * \
                     cos(k*x[0]-sqrt(g*depth)*k*t) * k*k)"
        eta_exact = "eta0*cos(k*x[0]-sqrt(g*depth)*k*t)"

        # The source term
        source = Expression((ddu_exact,
                            "0.0"),
                            eta0=eta0, g=problem.parameters.g,
                            depth=problem.parameters.depth,
                            t=problem.parameters.current_time,
                            k=k, viscosity=problem.parameters.viscosity)

        adj_reset()
        parameters = ShallowWaterSolver.default_parameters()
        parameters.dump_period = -1
        solver = ShallowWaterSolver(problem, parameters, config)
        solver.solve(state, annotate=False,
                                     u_source=source)

        analytic_sol = Expression((u_exact,
                                  "0",
                                  eta_exact),
                                  eta0=eta0, g=problem.parameters.g,
                                  depth=problem.parameters.depth,
                                  t=problem.parameters.current_time,
                                  k=k)
        return errornorm(analytic_sol, state)


    def compute_spatial_error(self, linear_problem_params, refinement_level):
        nx = 2 * 2**refinement_level
        ny = 1

        config = configuration.DefaultConfiguration(nx=nx, ny=ny)
        domain = domains.RectangularDomain(3000, 1000, nx, ny)
        config.set_domain(domain)

        eta0 = 2.0
        k = pi/config.domain.basin_x
        linear_problem_params.start_time = 0.0
        linear_problem_params.finish_time = (pi/(sqrt(linear_problem_params.g *
                                         linear_problem_params.depth) * k) / 1000)
        linear_problem_params.dt = linear_problem_params.finish_time / 2
        linear_problem_params.include_viscosity = True
        linear_problem_params.viscosity = 10.0
        linear_problem_params.flather_bc_expr = Expression(
            ("2*eta0*sqrt(g/depth)*cos(-sqrt(g*depth)*k*t)", "0"),
            eta0=eta0,
            g=linear_problem_params.g,
            depth=linear_problem_params.depth,
            t=linear_problem_params.current_time,
            k=k
        )

        problem = ShallowWaterProblem(linear_problem_params)

        return self.error(problem, config, eta0, k)

    def test_spatial_convergence_is_two(self, sw_linear_problem_parameters):
        errors = []
        tests = 4
        for refinement_level in range(tests):
            errors.append(self.compute_spatial_error(sw_linear_problem_parameters,
                                                     refinement_level))
        # Compute the order of convergence
        conv = []
        for i in range(len(errors)-1):
            conv.append(abs(math.log(errors[i+1]/errors[i], 2)))

        log(INFO, "Errors: %s" % errors)
        log(INFO, "Spatial order of convergence (expecting 2.0): %s" % conv)
        assert min(conv) > 1.8