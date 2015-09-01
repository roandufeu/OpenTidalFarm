from dolfin import dot, Constant, dx
from dolfin import *
from prototype_functional import PrototypeFunctional


class EcologyFunctional(PrototypeFunctional):

    def __init__(self, problem, u_ambient, alpha=0.8, beta=3):
        self.problem = problem
        self.u_ambient = u_ambient
        self.alpha = alpha
        self.beta = beta
        self.farm = problem.parameters.tidal_farm
        self.rho = problem.parameters.rho
        self.farm.update()

    def ecology(self, state, times):

        n = len(times)
        #m  = self.solver.optimisation_iteration * 100
        p = n % 100 # n-m

#        if self.solver.optimisation_iteration == 0:
#            state >> File(os.path.join("ambient/", p, "_u.xml")
#            u_ambient = state
#        else:
#            Vamb = VectorFunctionSpace(self.problem.parameters.domain.mesh,
#                                      'CG', 1, dim=2)

      #  from IPython import embed; embed()
        ua_x = self.u_ambient[p][0]
        ua_y = self.u_ambient[p][1]
        Ua = (ua_x**2 + ua_y**2)**0.5

        # condition the current flow
        u_x = state[0]
        u_y = state[1]
        U = (u_x**2 + u_y**2)**0.5

        # import parameters
        penalty_multiplier = self.problem.parameters.eco_penalty_multiplier


        penalty_function = penalty_multiplier*((0.5*(dolfin.tanh(
            self.beta*(((U-Ua)**2)**0.5 - self.alpha))+1))**2)
        return penalty_function

    def Jt(self, state, turbine_field, times):
        eco_value = assemble(
            self.ecology(state, times)*self.problem.parameters.eco_dx)
        print "ecology = ", eco_value

        # the objective functional
        return self.ecology(state, times)*self.problem.parameters.eco_dx

    def Jt_individual(self, state, i):
        """ Computes the power output of the i'th turbine.

        :param state: Current solution state
        :type state: UFL
        :param i: refers to the i'th turbine
        :type i: Integer

        """
        print 'Jt_individual not yet implemented for ecology'
        return None
