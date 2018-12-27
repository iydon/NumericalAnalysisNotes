#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/12/27 09:54
# @Author   : Iydon
# @File     : course8.2.py

import Poly
import sympy as sp



def least_square_poly(x, y, a, b, n, sympy, Poly, display=False):
    A  = sympy.zeros(n+1)
    b_ = sympy.zeros(n+1,1)
    for i in range(n+1):
        for j in range(n+1):
            A[i,j] = sympy.integrate(x**(i+j), (x,a,b))
        b_[i] = sympy.integrate(x**i*y, (x,a,b))
    x = A**-1 * b_
    f = Poly.Poly(x.T.tolist()[0][::-1])
    if display:
        print(A)
        print(b_.T)
        print(x.T)
    return f


x = sp.symbols("x")
y = sp.exp(x)
result = least_square_poly(x, y, 0, 2, 1, sp, Poly, display=True)
print(result)
