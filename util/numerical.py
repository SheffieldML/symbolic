
import numpy as np

def piecewise(*pairs):
    expressions, conditions = zip(*pairs)
    bc = np.broadcast(*conditions)
    first_true = np.empty(bc.shape, dtype=int)
    first_true.flat = [np.argmax(conds) for conds in bc]
    return np.choose(first_true, expressions)
