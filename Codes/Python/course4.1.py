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

xs = [1.8, 1.9, 2.0, 2.1, 2.2]
d  = n_point_differentiation(xs, [np.exp(x)*x for x in xs])
print(d.lambdify()(2))
