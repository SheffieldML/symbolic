import numpy as np
import sympy as sym
from GPy.util.symbolic import differfln
from sympy import Function, ImmutableMatrix, Piecewise, And, Or, Eq, Gt, Lt

from symbolic import Symbolic


class Ode1_lfm(Symbolic):
    """
    A symbolic covariance based on a first order differential equation being driven by a latent force. 

    """

    def __init__(self, k_uu, k_fu=None, k_ff=None, output_dim=1, parameters=None, name='Ode1_lfm', func_modules=[], cse=True):

        x_0, x_1, z_0, z_1 = sym.symbols('x_:2, z_:2', positive=True)

#       does not work yet, because of a bug in sympy.cse
#        k = Sel(ImmutableMatrix([[k_ff, k_fx], [k_fx, k_xx]]), x_2, z_1)
        

        

        if k_fu is None:
            x_0, z_0, decay = sym.symbols('x_0, z_0, decay', positive=True)
            scale = sym.symbols('scale')
            k_fu = scale*sym.exp(-decay*x_0)*sym.integrate(sym.exp(decay*x_0) * k_uu, (x_0, 0, x_0))
        if k_ff is None:
            x_0, z_0, decay = sym.symbols('x_0, z_0, decay', positive=True)
            scale = sym.symbols('scale')
            k_ff = scale**2*sym.exp(-decay*(x_0+z_0))*sym.integrate(sym.exp(decay*(x_0+z_0)) * k_uu, (x_0, 0, x_0), (z_0, 0, z_0))

        k_uf = k_fu.subs( ((x_0,z_0), (z_0,x_0)),simultaneous=True)

        k = Piecewise(
                (k_ff, And(Eq(x_1,1), Eq(z_1,1))),
                (k_fu, And(Eq(x_1,0), Eq(z_1,1))),
                (k_uf, And(Eq(x_1,1), Eq(z_1,0))),
                (k_uu, True)

            )
        super(Ode1_lfm, self).__init__(input_dim=2, k=k, output_dim=output_dim, parameters=parameters, name=name, func_modules=func_modules, cse=cse)


