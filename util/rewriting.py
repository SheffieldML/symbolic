from collections import defaultdict
import sympy as sym
from itertools import product


def split_on_condition(to_split, condition):
    a,b = [], []
    for x in to_split:
        (a if condition(x) else b).append(x)
    return a,b

def is_product_with_erf(expr):  
    return expr.func == sym.erf or (expr.func == sym.Mul and any(x.func == sym.erf for x in expr.args))

def is_sum_with_erf(expr):
    return is_product_with_erf(expr) or (expr.func == sym.Add and any(is_product_with_erf(x) for x in expr.args))

def parse_sum_with_erf(expr):
    w = sym.Wild("w")
    c = sym.collect(expr, sym.erf(w), evaluate=False)
    return (c.pop(1) if 1 in c else 1, {key.args[0]: value for key, value in c.items()})


def rewrite_erfs(expr):
    c, erfdict = parse_sum_with_erf(expr)
    return c+differfln_terms(erfdict)


#splits a product into an integer part and the rest
def split_off_int(expr):
    if expr.func == sym.Mul:
        return tuple( sym.Mul(*x) for x in split_on_condition(expr.args, lambda x:x.is_integer) )
    else:
        return (expr,1) if expr.is_integer else (1,expr)
  
def pairterms(counts):
    positives, negatives = split_on_condition(counts.items(), lambda x:x[1]>0)
    print positives
    print negatives
    negatives = [ (y,-n) for y,n in negatives]
#    positives = sorted(positives, key=lambda x:tuple(reversed(x)), reverse=True)
#    negatives = sorted(positives, key=lambda x:tuple(reversed(x)), reverse=True)
    while positives:
        x,p =  positives[0]
        y,n = negatives[0]
        if n > p:
            negatives[0] = (y, n - p)
            positives.pop()
            c = p
        elif p > n:
            positives[0] = (x, p - n)
            negatives.pop()
            c = n
        else:
            negatives.pop()
            positives.pop()
            c = p
        yield c*differfln(x,y)




def reorder(erfdict):
    invdict = defaultdict(dict)
    for key, value in erfdict.items():
        intpart, exprpart = split_off_int(value)        
        invdict[exprpart][key] = intpart
    return invdict

def differfln_terms(erfdict):
    invdict = reorder(erfdict)
    #return sum( factor*sym.Mul(*pairterms(counts)) for factor, counts in invdict )
    return sum( factor * sum(pairterms(flip_erfs(counts))) for factor, counts in invdict.items() )

def test_hack(kernel):
    return [ (x, is_sum_with_erf(x)) for x in kernel.args]

def rewrite_hack(kernel):
    return sym.Mul( *(rewrite_erfs(expr) if is_sum_with_erf(expr) else expr for expr in kernel.args))

def flip_erfs(counts):
    args, weights = zip(*counts.items())
    for signs in product([1,-1], repeat=len(counts)):
        if sum( a*b for a,b in zip(signs, weights) ) == 0:
            return {sign*arg: sign*weight for arg, weight, sign in zip(args, weights, signs)}


    

