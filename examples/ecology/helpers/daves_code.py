    def find_ambient_flow_field(self):
        """ Solve a dummy of the real problem to yield the ambient flow (i.e.
        the flow over the domain in the abscence of turbines)
        """
        ambient_solve = self.solver.solve()
        ambient_state = ambient_solve.next()
        ambient_state = ambient_solve.next()
        self.ambient_state = ambient_state['state']


        state = self.ambient_state

        mesh = self.problem_params.domain.mesh
        coords = mesh.coordinates()

        Output_V = dolfin.VectorFunctionSpace(mesh, 'CG', 1, dim=2)
        u_out = dolfin.TrialFunction(Output_V)
        v_out = dolfin.TestFunction(Output_V)
        M_out = dolfin.assemble(dolfin.inner(v_out, u_out) * dolfin.dx)

        out_state = dolfin.Function(Output_V)
        rhs = dolfin.assemble(dolfin.inner(v_out, state.split()[0]) * dolfin.dx)
        dolfin.solve(M_out, out_state.vector(), rhs, 'cg', 'sor')

        out_state = out_state.vector().array()
        velocity = out_state.reshape(len(out_state)/2,2).transpose()
        velocity = numpy.sqrt(velocity[0]**2 + velocity[1]**2)

        x = coords[:,0]
        y = coords[:,1]



