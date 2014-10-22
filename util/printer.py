from sympy.printing.lambdarepr import LambdaPrinter
from sympy.utilities import default_sort_key

class VectorizedPrinter(LambdaPrinter):


#TODO!
#    def _print_Piecewise(self, expr):
#        from sympy.sets.sets import Interval
#        result = []
#        i = 0
#        for arg in expr.args:
#            e = arg.expr
#            c = arg.cond
#            result.append('((')
#            result.append(self._print(e))
#            result.append(') if (')
#            result.append(self._print(c))
#            result.append(') else (')
#            i += 1
#        result = result[:-1]
#        result.append(') else None)')
#        result.append(')'*(2*i - 2))
#        return ''.join(result)

    def _print_And(self, expr):
        result = ['(']
        for arg in sorted(expr.args, key=default_sort_key):
            result.extend(['(', self._print(arg), ')'])
            result.append(' & ')
        result = result[:-1]
        result.append(')')
        return ''.join(result)

    def _print_Or(self, expr):
        result = ['(']
        for arg in sorted(expr.args, key=default_sort_key):
            result.extend(['(', self._print(arg), ')'])
            result.append(' | ')
        result = result[:-1]
        result.append(')')
        return ''.join(result)

    def _print_Not(self, expr):
        result = ['(', '~ (', self._print(expr.args[0]), '))']
        return ''.join(result)

