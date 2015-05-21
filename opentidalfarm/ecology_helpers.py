from dolfin import CellFunction, MeshFunction
from opentidalfarm import *


class Eco_function_square(Expression):
    def __init__(self, eco_x_start, eco_y_start, eco_x_end, eco_y_end):
        self.eco_x_start = eco_x_start
        self.eco_y_start = eco_y_start
        self.eco_x_end = eco_x_end
        self.eco_y_end = eco_y_end

    def eval(self, values, x):
        if x[0] > self.eco_x_start and x[1] > self.eco_y_start and \
           x[0] < self.eco_x_end and x[1] < self.eco_y_end:
            values[:] = Constant(1)
        else:
            values[:] = Constant(0)


class Eco_function_corner(Expression):
    def __init__(self, eco_x_start, eco_y_start, site_x_end, site_y_end):
        self.eco_x_start = eco_x_start
        self.eco_y_start = eco_y_start
        self.site_x_end = site_x_end
        self.site_y_end = site_y_end

    def eval(self, values, x):
        if x[0] > self.eco_x_start and x[1] > self.site_y_end or \
           x[0] > self.site_x_end and x[1] > self.eco_y_start:
            values[:] = Constant(1)
        else:
            values[:] = Constant(0)


class Bathymetry(Expression):
    def __init__(self, prob_params, max_depth, min_depth, channel_width):
        self.prob_params = prob_params
        self.step = (channel_width)/3
        self.max_depth = max_depth
        self.min_depth = min_depth
        self.m = (self.max_depth - self.min_depth)/self.step
        self.m2 = (self.min_depth - self.max_depth)/self.step

    def eval(self, values, x):
        if x[1] > 2 * self.step:
            values[:] = self.m2 * (x[1] - 2 * self.step) + self.max_depth
        elif x[1] > self.step and x[1] <= 2 * self.step:
            values[:] = self.max_depth
        else:
            values[:] = (self.m * x[1] + self.min_depth)


class Eco_domain(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] > 500 and x[0] < 640 and x[1] > 160 \
            and x[1] < 320 else False
