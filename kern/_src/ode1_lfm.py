import numpy as np
import sympy as sym
from GPy.util.symbolic import differfln
from sympy import Function, ImmutableMatrix, Piecewise, And, Or, Eq

from symbolic import Symbolic


class Ode1_lfm(Symbolic):
    """
    A symbolic covariance based on a first order differential equation being driven by a latent force. 

    """

    def __init__(self, k_uu, k_fu, k_ff, output_dim=1, parameters=None, name='Ode1_lfm', func_modules=[]):

        x_1, z_1 = sym.symbols('x_1, z_1', positive=True)

#       does not work yet, because of a bug in sympy.cse
#        k = Sel(ImmutableMatrix([[k_ff, k_fx], [k_fx, k_xx]]), x_1, z_1)
        
        k = Piecewise(

                (k_ff, And(Eq(x_1,1), Eq(z_1,1))),
                (k_fu, Or(And(Eq(x_1,0), Eq(z_1,1), And(Eq(x_1,1), Eq(z_1,0))))), #(x_1 == 0 && z_1 == 1) or (x_1 == 1 && z_1 == 0)
                (k_uu, True)

            )
        super(Ode1_lfm, self).__init__(input_dim=2, k=k, output_dim=output_dim, parameters=parameters, name=name, func_modules=func_modules)


