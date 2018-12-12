#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/04 16:16
# @Author   : Iydon
# @File     : course4.1.py


from Poly import *
import numpy as np

def n_point_differentiation(xs, fxs):
    """
    Midpoint differentiation.
    """
    mid = len(xs)//2 + 1
    def Lagrange(i):
        def prod(lst, result=1):
            return prod(lst[:-1], result*lst[-1]) if lst else result
        num = prod([Poly([1,-x]) for x in xs[:i]+xs[i+1:]])
        den = prod([(xs[i]-x) for x in xs[:i]+xs[i+1:]])
        return num / den
    result = 0
    for i in range(2*mid-1):
        result += fxs[i] * Lagrange(i)
    return result.diff()

"""
xs = [1.8, 1.9, 2.0, 2.1, 2.2]
d  = n_point_differentiation(xs, [np.exp(x)*x for x in xs])
print(d)
"""

def three_point(xs:list, fxs:list, pos:int):
    """
    Args:
        pos: -1(end), 0(mid), 1(first).
    """
    h = (xs[-1]-xs[0]) / (len(xs)-1)
    d = {-1: -3*fxs[0]+4*fxs[1]-fxs[2],
          0: fxs[2]-fxs[0],
          1: fxs[0]-4*fxs[1]+3*fxs[2],}
    return d[pos] / 2 / h

from sympy import *
x = symbols("x")
f = 3*x*exp(x)-cos(x)
g = diff(diff(f))
print(g.subs(x, 1.3))

h = 0.1
d = lambda x: 1/h/h*(x[0]-2*x[1]+x[2])
print(d([11.59006,14.04276,16.86187])-g.subs(x, 1.3))

h = 0.01
d = lambda x: 1/h/h*(x[0]-2*x[1]+x[2])
print(d([13.78176,14.04276,14.30741])-g.subs(x, 1.3))

