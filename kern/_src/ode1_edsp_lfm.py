import numpy as np
import sympy as sym
import GPy 
from sympy import sqrt, exp, erf, Le

from ode1_lfm import Ode1_lfm
from ...util.symbolic import *





class Ode1_edsp_lfm(Ode1_lfm):
    """
    A symbolic covariance based on a first order differential equation being driven by a latent force that smoothly switches between two random levels. 

    """
    def __init__(self, output_dim=1, variance=1, lengthscale=1, scale=1, decay=1, switchingpoint=0, name='Ode1_edsp_lfm', cse=True):

        parameters = {'variance' : variance, 'lengthscale' : lengthscale, 'scale' : scale, 'decay' : decay, 'switchingpoint' : switchingpoint}

        input_dim = 2

        x_0, z_0, variance, lengthscale, sp = sym.symbols('x_0, z_0, variance, lengthscale, sp', positive=True)
        x, scale = sym.symbols('x, scale')


        transition = Piecewise( (1, Le(x,sp)), (sym.exp(-x) / lengthscale, True) )
        
        k_uu = variance * (transition.subs(x,x_0) * transition.subs(x,z_0) + (1-transition.subs(x,x_0)) * (1-transition.subs(x,z_0)))



        super(Ode1_edsp_lfm, self).__init__(k_uu, output_dim=output_dim, parameters=parameters, name=name, func_modules=[{'differfln':GPy.util.functions.differfln}], cse=cse)
