from dolfin import dot, dx
from dolfin import *
from prototype_functional import PrototypeFunctional


class SedimentFunctional(PrototypeFunctional):

    def __init__(self, problem, u_ambient, beta=3):
        self.problem = problem
        self.u_ambient = u_ambient
        self.beta = beta
        self.farm = problem.parameters.tidal_farm
        self.rho = problem.parameters.rho
        self.farm.update()

    def sediment(self, state):

        # condition the ambient flow
        ua_x = self.u_ambient[0]
        ua_y = self.u_ambient[1]
        Ua = (ua_x**2 + ua_y**2)**0.5

        # condition the current flow
        u_x = state[0]
        u_y = state[1]
        U = (u_x**2 + u_y**2)**0.5

        # divide domain into sediment regions
        boulders = conditional(ge(Ua, 3.19), 1, 2)
        coarse_gravel_upper = conditional(lt(Ua, 3.19), 1, 2)
        coarse_gravel_lower = conditional(ge(Ua, 2.16), 1, 2)
        med_gravel_upper = conditional(lt(Ua, 2.16), 1, 2)
        med_gravel_lower = conditional(ge(Ua, 1.49), 1, 2)
        fine_gravel_upper = conditional(lt(Ua, 1.49), 1, 2)
        fine_gravel_lower = conditional(ge(Ua, 0.701), 1, 2)
        coarse_sand_upper = conditional(lt(Ua, 0.701), 1, 2)
        coarse_sand_lower = conditional(ge(Ua, 0.325), 1, 2)
        med_sand_upper = conditional(lt(Ua, 0.325), 1, 2)
        med_sand_lower = conditional(ge(Ua, 0.275), 1, 2)
        fine_sand_upper = conditional(lt(Ua, 0.275), 1, 2)
        #        from IPython import embed; embed()

        if boulders == 1:
            upper_limit = 10
            lower_limit = 3.19
        elif (coarse_gravel_upper == 1) and (coarse_gravel_lower == 1):
            upper_limit = 3.19
            lower_limit = 2.16
        elif (med_gravel_upper == 1) and (med_gravel_lower == 1):
            upper_limit = 2.16
            lower_limit = 1.49
        elif (fine_gravel_upper == 1) and (fine_gravel_lower == 1):
            upper_limit = 1.49
            lower_limit = 0.701
        elif (coarse_sand_upper == 1) and (coarse_sand_lower == 1):
            upper_limit = 0.701
            lower_limit = 0.325
        elif (med_sand_upper == 1) and (med_sand_lower == 1):
            upper_limit = 0.325
            lower_limit = 0.275
        elif (fine_sand_upper == 1):
            upper_limit = 0.275
            lower_limit = 0.0
        else:
            upper_limit = 10
            lower_limit = 0
#            raise ValueError(
#                "cannot split into sediment domains - invalid flow velocity?")

        # import parameters
        penalty_multiplier = (self.problem.parameters.sed_penalty_multiplier *
            self.problem.parameters.multiplier)
        eco_domain = self.problem.parameters.eco_domain

        penalty_function = eco_domain * (
            penalty_multiplier * 0.5 * (
                dolfin.tanh(self.beta * ((U) - upper_limit-0.2)) +
                dolfin.tanh(self.beta * (lower_limit - 0.2 - (U))) + 2)**2)
        return penalty_function

    def Jt(self, state, turbine_field):

        sed_value = assemble(self.sediment(state) *
                             self.problem.parameters.domain.dx)

        print "sediment = ", sed_value

        # the objective functional
        return -self.sediment(state)*self.problem.parameters.domain.dx

    def Jt_individual(self, state, i):
        """ Computes the power output of the i'th turbine
        :param state: Current solution state
        :type state: UFL
        :param i: refers to the i'th turbine
        :type i: Integer

        """
        print 'Jt_individual not yet implemented for sediment'
        return None
