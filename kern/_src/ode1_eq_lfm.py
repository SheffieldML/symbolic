import numpy as np
import sympy as sym
import GPy 

from ode1_lfm import Ode1_lfm
from ...util.symbolic import *

class Ode1_eq_lfm(Ode1_lfm):
    """
    A symbolic covariance based on a first order differential equation being driven by a latent force that is an exponentiated quadratic. 

    """
    def __init__(self, output_dim=1, decay=1, variance=1, lengthscale=1, name='Ode1_eq_lfm'):

        parameters = {'decay' : decay, 'variance' : variance, 'lengthscale' : lengthscale}

        input_dim = 2

        x_0, x_1, z_0, z_1, decay, variance, lengthscale = sym.symbols('x_:2, z_:2, decay, variance, lengthscale', positive=True)
        scale = sym.symbols('scale')
        
        k_ff = variance*sym.exp( -(x_0-z_0)**2/lengthscale**2 )
        k_fx = ( sym.sqrt(sym.pi)*scale/2 *
                 sym.exp((decay*lengthscale/2)**2) *
                 sym.exp(-decay*(z_0-x_0)) *( sym.erf((z_0-x_0)/lengthscale - decay*lengthscale/2) +
                   sym.erf(x_0/lengthscale + decay*lengthscale/2)
                 )
               )
        k_xx = scale**2*(    h(x_0, z_0, decay, decay, lengthscale)
                           + h(z_0, x_0, decay, decay, lengthscale)
                        )

        super(Ode1_eq_lfm, self).__init__(k_ff=k_ff, k_fx=k_fx, k_xx=k_xx, output_dim=output_dim, parameters=parameters, name=name, func_modules=[{'differfln':GPy.util.functions.differfln}])
