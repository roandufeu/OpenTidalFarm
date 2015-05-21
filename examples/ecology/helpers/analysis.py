from dolfin import *
from opentidalfarm import *
import opentidalfarm as otf
import numpy as np
import scipy.interpolate

domain = FileDomain('mesh/earth_orkney_converted_coarse.xml')
V, H = otf.finite_elements.p2p1(domain.mesh)
func_space = MixedFunctionSpace([V, H])
u_ambient = Function(func_space, 'ambient.xml')

u_res_sed = Function(func_space, 'iter_100/p2p1_ures.xml')
u_res = Function(func_space, '../channel_opt_orkney_control/iter_100/p2p1_ures.xml')

u_amb_np = u_ambient.vector().array()
u_res_np = u_res.vector().array()
u_res_sed_np = u_res_sed.vector().array()

mesh = domain.mesh

# attempting to isolate the scalar flow field...
Output_V = VectorFunctionSpace(mesh, 'CG', 1, dim=2)
u_out = TrialFunction(Output_V)
v_out = TestFunction(Output_V)

M_out_res = assemble(inner(v_out, u_out)*dx)
out_state_res = Function(Output_V)
M_out_amb = assemble(inner(v_out, u_out)*dx)
out_state_amb = Function(Output_V)
M_out_res_sed = assemble(inner(v_out, u_out)*dx)
out_state_res_sed = Function(Output_V)

rhs_amb = assemble(inner(v_out, u_ambient.split()[0]) * dx)
rhs_res = assemble(inner(v_out, u_res.split()[0]) * dx)
rhs_res_sed = assemble(inner(v_out, u_res_sed.split()[0]) * dx)

solve(M_out_amb, out_state_amb.vector(), rhs_amb, "cg", "sor", annotate=False)
out_state_amb_np = out_state_amb.vector().array()

solve(M_out_res, out_state_res.vector(), rhs_res, "cg", "sor", annotate=False)
out_state_res_np = out_state_res.vector().array()

solve(M_out_res_sed, out_state_res_sed.vector(), rhs_res_sed, "cg", "sor", annotate=False)
out_state_res_sed_np = out_state_res_sed.vector().array()

coords = mesh.coordinates()
x = coords[:, 0]
y = coords[:, 1]

xi, yi = numpy.linspace(x.min(), x.max(), 3000), numpy.linspace(y.min(), y.max(), 3000)
xi, yi = numpy.meshgrid(xi, yi)

# recreate eco domain

eco_x_start = 10307240. + 1100 +600
eco_y_start = 6522134.31356 + 700 + 600
eco_x_end = eco_x_start + 400
eco_y_end = eco_y_start + 500

temp = FunctionSpace(domain.mesh, 'CG', 1)
eco_fn = Eco_function_square(eco_x_start, eco_y_start, eco_x_end, eco_y_end)
eco_domain = interpolate(eco_fn, temp)

# from IPython import embed; embed()

# Ua_scalar = scipy.interpolate.griddata((x, y), u_amb_np, (xi, yi), method='linear')
# U_scalar = scipy.interpolate.griddata((x, y), u_res_np, (xi, yi), method='linear')

# from IPython import embed; embed()

# condition the ambient flow
ua_x = u_ambient[0]
ua_y = u_ambient[1]
Ua = (ua_x**2 + ua_y**2)**0.5
# condition the current flow
u_x = u_res[0]
u_y = u_res[1]
U = (u_x**2 + u_y**2)**0.5

us_x = u_res_sed[0]
us_y = u_res_sed[1]
Us = (us_x**2 + us_y**2)**0.5

scoured_amb = conditional(ge(Ua, 3.19), 1, 0)
coarse_gravel_amb = conditional(ge(Ua, 2.16), 1, 0)*conditional(lt(Ua, 3.19), 1, 0)
med_gravel_amb = conditional(lt(Ua, 2.16), 1, 0)*conditional(ge(Ua, 1.49), 1, 0)
fine_gravel_amb = conditional(lt(Ua, 1.49), 1, 0)*conditional(ge(Ua, 0.701), 1, 0)
coarse_sand_amb = conditional(lt(Ua, 0.701), 1, 0)*conditional(ge(Ua, 0.325), 1, 0)
med_sand_amb = conditional(lt(Ua, 0.325), 1, 0)*conditional(ge(Ua, 0.275), 1, 0)
fine_sand_amb = conditional(lt(Ua, 0.275), 1, 0)

scoured_res = conditional(ge(U, 3.19), 1, 0)
coarse_gravel_res = conditional(ge(U, 2.16), 1, 0)*conditional(lt(U, 3.19), 1, 0)
med_gravel_res = conditional(lt(U, 2.16), 1, 0)*conditional(ge(U, 1.49), 1, 0)
fine_gravel_res = conditional(lt(U, 1.49), 1, 0)*conditional(ge(U, 0.701), 1, 0)
coarse_sand_res = conditional(lt(U, 0.701), 1, 0)*conditional(ge(U, 0.325), 1, 0)
med_sand_res = conditional(lt(U, 0.325), 1, 0)*conditional(ge(U, 0.275), 1, 0)
fine_sand_res = conditional(lt(U, 0.275), 1, 0)


scoured_res_sed = conditional(ge(Us, 3.19), 1, 0)
coarse_gravel_res_sed = conditional(ge(Us, 2.16), 1, 0)*conditional(lt(Us, 3.19), 1, 0)
med_gravel_res_sed = conditional(lt(Us, 2.16), 1, 0)*conditional(ge(Us, 1.49), 1, 0)
fine_gravel_res_sed = conditional(lt(Us, 1.49), 1, 0)*conditional(ge(Us, 0.701), 1, 0)
coarse_sand_res_sed = conditional(lt(Us, 0.701), 1, 0)*conditional(ge(Us, 0.325), 1, 0)
med_sand_res_sed = conditional(lt(Us, 0.325), 1, 0)*conditional(ge(Us, 0.275), 1, 0)
fine_sand_res_sed = conditional(lt(Us, 0.275), 1, 0)

# from IPython import embed; embed()

Ua_dom = Ua*eco_domain
U_dom = U*eco_domain
Us_dom = Us*eco_domain

# Ua_dom_np =  Ua_dom.vector().array()
# U_dom_np = U_dom.vector().array()
# Us_dom_np = Us_dom.vector().array()

from IPython import embed; embed()

# class Sediment(Expression):
#    def __init__(self, U)
#        self.U = U
#    def eval(self, values, x):
