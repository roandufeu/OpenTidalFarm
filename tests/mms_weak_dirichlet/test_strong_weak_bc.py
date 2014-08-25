''' Tests the spatial order of convergence with weakly imposed
    Dirichlet conditions. '''

import sys
import math
from opentidalfarm import *
from opentidalfarm.initial_conditions import SinusoidalInitialCondition
from dolfin_adjoint import adj_reset
from dolfin import log, INFO, ERROR


class TestWeakDirichletBoundaryConditions(object):
    
    def error(self, problem, config, eta0, k):
        state = Function(config.function_space)
        ic_expr = SinusoidalInitialCondition(eta0, k, 
                                             problem.parameters.depth,
                                             problem.parameters.start_time)
        ic = project(ic_expr, state.function_space())
        state.assign(ic, annotate=False)

        adj_reset()
        parameters = ShallowWaterSolver.default_parameters()
        parameters.dump_period = -1
        solver = ShallowWaterSolver(problem, parameters, config)
        solver.solve(state, annotate=False)

        analytic_sol = Expression(
            ("eta0*sqrt(g/depth)*cos(k*x[0]-sqrt(g*depth)*k*t)", "0",
             "eta0*cos(k*x[0]-sqrt(g*depth)*k*t)"),
            eta0=eta0, g=problem.parameters.g,
            depth=problem.parameters.depth,
            t=problem.parameters.current_time, k=k)
        return errornorm(analytic_sol, state)


    def compute_spatial_error(self, problem_params, refinement_level):
        nx = 2**refinement_level
        ny = 1
        config = configuration.DefaultConfiguration(
            nx, ny, finite_element=finite_elements.p1dgp2) 
        domain = domains.RectangularDomain(3000, 1000, nx, ny)
        config.set_domain(domain)
        eta0 = 2.0
        k = pi/config.domain.basin_x
        problem_params.finish_time = pi / (sqrt(problem_params.g *
                                        problem_params.depth) * k) / 20
        problem_params.dt = problem_params.finish_time / 10
        problem_params.bctype = "dirichlet"
        problem_params.u_weak_dirichlet_bc_expr = Expression(
            ("eta0*sqrt(g/depth)*cos(k*x[0]-sqrt(g*depth)*k*t)", "0"), 
            eta0=eta0,
            g=problem_params.g,
            depth=problem_params.depth,
            t=problem_params.current_time,
            k=k)

        problem = ShallowWaterProblem(problem_params)

        return self.error(problem, config, eta0, k)

    def compute_temporal_error(self, problem_params, refinement_level):
        nx = 2**3
        ny = 1
        config = configuration.DefaultConfiguration(
            nx, ny, finite_element=finite_elements.p1dgp2) 
        domain = domains.RectangularDomain(3000, 1000, nx, ny)
        config.set_domain(domain)
        eta0 = 2.0
        k = pi/config.domain.basin_x
        problem_params.finish_time = 2 * pi / (sqrt(problem_params.g * 
                                        problem_params.depth) * k)
        problem_params.dt = problem_params.finish_time/(4*2**refinement_level)
        problem_params.theta = 0.5
        problem_params.bctype = "dirichlet"
        problem_params.u_weak_dirichlet_bc_expr = Expression(
            ("eta0*sqrt(g/depth)*cos(k*x[0]-sqrt(g*depth)*k*t)", "0"), 
            eta0=eta0,
            g=problem_params.g,
            depth=problem_params.depth,
            t=problem_params.current_time,
            k=k)

        problem = ShallowWaterProblem(problem_params)

        return self.error(problem, config, eta0, k)

    def test_spatial_convergence_is_two(self, sw_linear_problem_parameters):
        errors = []
        tests = 4
        for refinement_level in range(tests):
            error = self.compute_spatial_error(sw_linear_problem_parameters,
                                               refinement_level)
            errors.append(error)
        # Compute the order of convergence 
        conv = [] 
        for i in range(len(errors)-1):
            conv.append(abs(math.log(errors[i+1] / errors[i], 2)))

        log(INFO, "Spatial order of convergence (expecting 2.0): %s" % str(conv))
        assert min(conv) > 1.8

    def test_temporal_convergence_is_two(self, sw_linear_problem_parameters):
        errors = []
        tests = 4
        for refinement_level in range(tests):
            error = self.compute_temporal_error(sw_linear_problem_parameters,
                                                refinement_level)
            errors.append(error)
        # Compute the order of convergence 
        conv = [] 
        for i in range(len(errors)-1):
            conv.append(abs(math.log(errors[i+1] / errors[i], 2)))

        log(INFO, "Temporal order of convergence (expecting 2.0): %s" % str(conv))
        assert min(conv) > 1.8