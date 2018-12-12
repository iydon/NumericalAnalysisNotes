#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/09/23 09:23
# @Author   : Iydon
# @File     : course3.4.py

import numpy as np
from Poly import *


def Hermite(xs:list, fxs:list, dfxs:list):
    """
    Hermite polynomial
    ----------------------------
    H_{2n+1}(x) = \sum_{j=0}^nf(x_j)H_{n,j}(x)+\sum_{j=0}^nf'(x_j)\hat{H}_{n,j}(x)
    """
    def H(j):
        L  = Lagrange(xs, j)
        dL = L.diff().lambdify()
        H1 = (1+Poly([-2,2*xs[j]])*dL(xs[j])) * L**2
        H2 = Poly([1, -xs[j]]) * L**2
        return H1, H2

    def Lagrange(xs, j):
        def prod(lst, result=1):
            return prod(lst[:-1], result*lst[-1]) if lst else result
        num = prod([Poly([1,-x]) for x in xs[:j]+xs[j+1:]])
        den = prod([(xs[j]-x) for x in xs[:j]+xs[j+1:]])
        return num / den

    result = Poly([])
    for j in range(len(xs)):
        HH = H(j)
        result += fxs[j]*HH[0] + dfxs[j]*HH[1]
    return result



xs   = [1.3, 1.6, 1.9]
fxs  = [0.6200860, 0.4554022, 0.2818186]
dfxs = [-0.5220232, -0.5698959, -0.5811571]
h  = Hermite(xs, fxs, dfxs)
dh = h.diff()
print(h.lambdify()(1.5))
print(dh.lambdify()(1.5))