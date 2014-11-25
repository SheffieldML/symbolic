import numpy as np
import sympy as sym
import GPy 
from sympy import sqrt, exp, erf

from ode1_lfm import Ode1_lfm
from ...util.symbolic import *

class Ode1_eq_lfm(Ode1_lfm):
    """
    A symbolic covariance based on a first order differential equation being driven by a latent force that is an exponentiated quadratic. 

    """
    def __init__(self, output_dim=1, decay=1, variance=1, lengthscale=1, name='Ode1_eq_lfm', cse=True):

        parameters = {'decay' : decay, 'variance' : variance, 'lengthscale' : lengthscale}

        input_dim = 2

        x_0, x_1, x_2, z_0, z_1, z_2, decay, variance, lengthscale = sym.symbols('x_:3, z_:3, decay, variance, lengthscale', positive=True)
#        x_1, x_2, z_1, z_2, decay, variance, lengthscale = sym.symbols('x_1:3, z_1:3, decay, variance, lengthscale', positive=True)
        scale = sym.symbols('scale')
        
        k_uu = variance*sym.exp( -(x_1-z_1)**2/lengthscale**2 )
        k_fu = ( sym.sqrt(sym.pi)*lengthscale*scale/2 *
                 sym.exp((decay*lengthscale/2)**2
                         -decay*(z_1-x_1) +
                         +differfln( (z_1-x_1)/lengthscale - decay*lengthscale/2,
                                         -x_1 /lengthscale - decay*lengthscale/2
                                   )
                         )
               )
        k_ff = scale**2*(    h(x_1, z_1, decay, decay, lengthscale)
                           + h(z_1, x_1, decay, decay, lengthscale)
                        )


        super(Ode1_eq_lfm, self).__init__(k_uu, k_fu, k_ff, output_dim=output_dim, parameters=parameters, name=name, func_modules=[{'differfln':GPy.util.functions.differfln}], cse=cse)
