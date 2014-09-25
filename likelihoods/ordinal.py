# Copyright (c) 2014 The GPy authors (see AUTHORS.txt)
# Licensed under the BSD 3-clause license (see LICENSE.txt)


import sympy as sym
from ..util import gammaln, normcdfln, normcdf, create_selector
import numpy as np
from GPy.likelihoods import link_functions
from symbolic import Symbolic
from scipy import stats

class Ordinal(Symbolic):
    """
    Ordinal

    .. math::
        p(y_{i}|\pi(f_{i})) = \left(\frac{r}{r+f_i}\right)^r \frac{\Gamma(r+y_i)}{y!\Gamma(r)}\left(\frac{f_i}{r+f_i}\right)^{y_i}

    .. Note::
        Y takes non zero integer values..
        link function should have a positive domain, e.g. log (default).

    .. See also::
        symbolic.py, for the parent class
    """
    def __init__(self, categories=3, widths=None, gp_link=None):
        if gp_link is None:
            gp_link = link_functions.Identity()
        
        if widths is None:
            widths = np.ones(categories)*1./categories
        parameters = {'w':widths}

        
        dispersion = sym.Symbol('width', positive=True, real=True)
        y_0 = sym.Symbol('y_0', nonnegative=True, integer=True)
        f_0 = sym.Symbol('f_0', positive=True, real=True) 
        if categories>2:
            w = sym.symbols('w_:' + str(categories-2), positive=True, real=True)
        
        log_pdf = []
        log_pdf.append(normcdfln(-f_0))
        from sympy.parsing.sympy_parser import parse_expr
        if categories>2:
            old_w_sum=parse_expr('0')
            for i in range(1, categories-1):
                w_sum = old_w_sum+parse_expr('w_%i' %i) 
                log_pdf.append(sym.log(normcdf(old_w_sum + f_0) - normcdf(w_sum-f_0) ))
                old_w_sum = w_sum
            log_pdf.append(normcdfln(w_sum + f_0))
        else:
            log_pdf.append(normcdfln(f_0))
        log_pdf = create_selector(log_pdf, y_0)
        super(Ordinal, self).__init__(log_pdf=log_pdf, gp_link=gp_link, parameters=parameters, name='Ordinal')

        # TODO: Check this.
        self.log_concave = True

