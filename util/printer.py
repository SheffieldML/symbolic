from sympy.printing.lambdarepr import LambdaPrinter
from sympy.printing.str import StrPrinter
from sympy.printing.precedence import precedence
from sympy.utilities import default_sort_key

class VectorizedPrinter(StrPrinter):

    def _print_Piecewise(self, expr):
        result = [expr.func.__name__, '(']
        result.append( ','.join([ '(%s,%s)' % (self._print(e), self._print(c)) for e,c in expr.args]) )
        result.append(')')
        return ''.join(result)

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


    def _print_BooleanTrue(self, expr):
        return "True"

    def _print_BooleanFalse(self, expr):
        return "False"

