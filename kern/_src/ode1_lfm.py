import numpy as np
import sympy as sym
from GPy.util.symbolic import differfln
from sympy import Function, ImmutableMatrix, Piecewise, And, Or, Eq, Gt, Lt

from symbolic import Symbolic


class Ode1_lfm(Symbolic):
    """
    A symbolic covariance based on a first order differential equation being driven by a latent force. 

    """

    def __init__(self, k_uu, k_fu, k_ff, output_dim=1, parameters=None, name='Ode1_lfm', func_modules=[],cse=True):

        x_0, x_1, x_2, z_0, z_1, z_2 = sym.symbols('x_:3, z_:3', positive=True)
#        x_1, x_2, z_1, z_2 = sym.symbols('x_1, x_2, z_1, z_2', positive=True)

#       does not work yet, because of a bug in sympy.cse
#        k = Sel(ImmutableMatrix([[k_ff, k_fx], [k_fx, k_xx]]), x_2, z_1)
        
        k = Piecewise(

                (1, And(Gt(x_1,0.5), Gt(z_1,0.5))),
                (x_0 + z_0, Or(And(Lt(x_1,0.5), Gt(z_1,0.5), And(Gt(x_1,0.5), Lt(z_1,0.5))))), #(x_2 == 0 && z_1 == 1) or (x_2 == 1 && z_1 == 0)
                (3, True)

            )
        super(Ode1_lfm, self).__init__(input_dim=2, k=k, output_dim=output_dim, parameters=parameters, name=name, func_modules=func_modules, cse=cse)


