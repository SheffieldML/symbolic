
from numpy import *

def piecewise(*pairs):
    expressions, conditions = zip(*pairs)
    bc = broadcast(*conditions)
    first_true = empty(bc.shape, dtype=int)
    first_true.flat = [argmax(conds) for conds in bc]
    return choose(first_true, expressions)
