try:
    import sympy as sym
    sympy_available=True
except ImportError:
    sympy_available=False

import numpy as np
from symbolic import Symbolic

class Ode1(Symbolic):
    """
    A symbolic covariance based on ODE1 equation with an exponentiated quadratic covariance.

    """
    def __init__(self, input_dim=1, rank=1, output_dim=1, param=None, name='ode1'):
        #for LFM first and second order input dimension is 1
        x = sym.symbols('x')
        z = sym.symbols('z')
        decay_i, decay_j = sym.var('decay_i decay_j', positive=True)
        Si = sym.var('S_i:' + str(rank)) #Sensitivities for output i
        Sj = sym.var('S_j:' + str(rank)) #Sensititvities for ouput j
        shared_lengthscale = sym.var('shared_lengthscale_:' + str(rank), positive=True)
        tmp = sym.var('tmp') #temp variable that allows interchange between length scale
        decay_t = sym.var('decay_t')
        #code to define h functions
        h_0 = sym.exp((tmp*decay_j)**2/4.)*sym.exp(-decay_j*z)*(sym.exp(decay_j*x)*
              (sym.erf((z-x)/tmp - tmp*decay_j/2) + sym.erf(x/tmp + tmp*decay_j/2))
              - sym.exp(-decay_i*x)*(sym.erf(z/tmp - tmp*decay_j/2) + sym.erf(tmp*decay_j/2)))
              
        from sympy.parsing.sympy_parser import parse_expr
        #Initialization of kernel value
        kff = 0
        for q in range(rank): #for each latent function
            #Selection of lengthscale regarding to k        
            ls = parse_expr('shared_lengthscale_' + str(q))
            #Product of sensitivities            
            Sprod = parse_expr('S_i' + str(q))*parse_expr('S_j' + str(q))
            #Symbolic h calculation            
            h0_new = h_0.subs({tmp: ls})
            h1_new = h0_new.subs({ decay_i: decay_t, decay_j: decay_i})
            h1_new = h1_new.subs({decay_t: decay_j})

            dist = h0_new + h1_new
            # this is the covariance function    
            kff += Sprod*sym.sqrt(sym.pi)*ls/(2.*(decay_i + decay_j))*dist
        
        # extra input dim is to signify the output dimension.
        super(Ode1, self).__init__(input_dim=input_dim+1, k=kff, output_dim=output_dim, name=name)