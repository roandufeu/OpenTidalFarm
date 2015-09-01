from dolfin import dot, Constant, dx
from dolfin import *
from prototype_functional import PrototypeFunctional


class FishingFunctional(PrototypeFunctional):

    def __init__(self, problem):
        self.problem = problem
        self.farm = problem.parameters.tidal_farm
        self.rho = problem.parameters.rho
        self.farm.update()

    def fishing(self, state, times):

        n = len(times)
        p = n % 100 # n-m

	friction = self.problem.parameters.tidal_farm.turbine_frictions	

	from IPython import embed; embed()     

 	# condition the current flow                                             
        u_x = state[0]                                                           
        u_y = state[1]                                                           
        U = (u_x**2 + u_y**2)**0.5   


        return 1

    def Jt(self, state, turbine_field, times):
        eco_value = assemble(
            self.fishing(state, times)*self.problem.parameters.eco_dx)
        print "ecology = ", eco_value

        # the objective functional
        return self.fishing(state, times)*self.problem.parameters.eco_dx

    def Jt_individual(self, state, i):
        """ Computes the power output of the i'th turbine.

        :param state: Current solution state
        :type state: UFL
        :param i: refers to the i'th turbine
        :type i: Integer

        """
        print 'Jt_individual not yet implemented for ecology'
        return None
