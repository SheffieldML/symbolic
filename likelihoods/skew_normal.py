# Copyright (c) 2014 The GPy authors (see AUTHORS.txt)
# Licensed under the BSD 3-clause license (see LICENSE.txt)


import sympy as sym
from sympy.utilities.lambdify import lambdify

import numpy as np
from scipy import stats

from GPy.util.univariate_Gaussian import std_norm_pdf, std_norm_cdf
from GPy.util.functions import clip_exp
from GPy.likelihoods import link_functions

from ..util.symbolic import normcdfln, normcdf
from symbolic import Symbolic

class Skew_normal(Symbolic):
    """
    Skew Normal distribution.

    .. math::

    .. Note::
    Y takes real values.
    link function is identity

    .. See also::
    symbolic.py, for the parent class
    """
    def __init__(self, gp_link=None, shape=1.0, scale=1.0):
        parameters = {'scale': scale, 'shape':shape}
        if gp_link is None:
            gp_link = link_functions.Identity()

        # # this likelihood has severe problems with likelihoods saturating exponentials, so clip_exp is used in place of the true exp as a solution for dealing with the numerics.
        # func_modules = [{'exp':clip_exp}]
        func_modules = []

        scale = sym.Symbol('scale', positive=True, real=True)
        shape = sym.Symbol('shape', real=True)
        y_0 = sym.Symbol('y_0', real=True)
        f_0 = sym.Symbol('f_0', real=True) 
        log_pdf=-sym.log(scale)-1./2*sym.log(2*sym.pi)-1./2*((y_0-f_0)/scale)**2 + sym.log(2) + normcdfln(shape*(y_0-f_0)/scale) 
        super(Skew_normal, self).__init__(log_pdf=log_pdf, parameters=parameters, gp_link=gp_link, name='Skew_normal', func_modules=func_modules)

        self.log_concave = True
        


            
