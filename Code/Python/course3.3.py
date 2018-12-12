#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/09 04:05
# @Author   : Iydon
# @File     : course3.3.py


import numpy as np
from Poly import *


def divided_difference(xs:list, fxs:list):
    if len(xs)==1:
        return fxs[0]
    elif len(xs)>1:
        f = divided_difference
        return (f(xs[1:],fxs[1:])-f(xs[:-1],fxs[:-1]))/(xs[-1]-xs[0])


def newton_divided_difference(xs:list, fxs:list, forward:bool=True):
    """
    Newton's Divided-Difference Formula.
    """
    length = len(xs)
    result = [i for i in fxs]
    sign   = -1 if forward else 1
    for i in range(length-1):
        for j in range(length-i-1, 0, -1):
            idx = -length+i+j if forward else length-i-j-1
            result[idx]  = result[idx+sign]-result[idx]
            result[idx] /= xs[idx+sign*(i+1)]-xs[idx]
        # print(result)
    return result[::-sign]


def newton_forward(xs:list, fxs:list, degree:int, x0:float):
    def prod(lst, result=1):
        return prod(lst[:-1], result*lst[-1]) if lst else result
    def s_k(s, k):
        return prod([s-i for i in range(k)])
    # body
    h = xs[1] - xs[0]
    s = (x0-xs[0]) / h
    result = fxs[0]
    for k in range(1, degree+1):
        # print(s_k(s,k), divided_difference(xs[:k+1],fxs[:k+1]), h**k)
        result += s_k(s,k) * divided_difference(xs[:k+1],fxs[:k+1]) * h**k
    return result

def newton_backward(xs:list, fxs:list, degree:int, x0:float):
    def prod(lst, result=1):
        return prod(lst[:-1], result*lst[-1]) if lst else result
    def s_k(s, k):
        return prod([s-i for i in range(k)])
    # body
    h = xs[1] - xs[0]
    s = (x0-xs[-1]) / h
    result = fxs[-1]
    print(result)
    for k in range(1, degree+1):
        print((-1)**k, s_k(-s,k), divided_difference(xs[-k-1:],fxs[-k-1:]), h**k)
        result += (-1)**k * s_k(-s,k) * divided_difference(xs[-k-1:],fxs[-k-1:]) * h**k
    return result

"""
xs  = [1.0, 1.3, 1.6, 1.9, 2.2]
fxs = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
n   = newton_divided_difference(xs, fxs, forward=True)
print(n)
"""

xs  = [-0.75, -0.5, -0.25, 0]
fxs = [-0.07181250, -0.02475000, 0.33493750, 1.10100000]
n_f = newton_forward(xs, fxs, 3, -1/3)
print(n_f)

print("*"*10)

xs  = [1.0, 1.3, 1.6, 1.9, 2.2]
fxs = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
n_b = newton_backward(xs, fxs, 4, 2)
print(n_b)

print("*"*10)

# Newton
def newton_divided_difference(xs:list, fxs:list, degree:int=0):
    """
    Newton's Divided-Difference Formula.
    """
    def prod(lst, result=1):
        return prod(lst[:-1], result*lst[-1]) if lst else result
    result = fxs[0]
    for i in range(1, degree+1):
        print(divided_difference(xs[:i+1],fxs[:i+1]), [Poly([1,-x]) for x in xs[:i]])
        result += divided_difference(xs[:i+1],fxs[:i+1]) * prod([Poly([1,-x]) for x in xs[:i]])
    return result

xs  = [-.1, .0, .2, .3, 0.35]
fxs = [5.30, 2.00, 3.19, 1.00, 0.9726]
n = newton_divided_difference(xs, fxs, 4)
print(n)

print("*"*10)

from sympy import *

x = symbols("x")
a,b,c,d,e = symbols("a,b,c,d,e")
a,b,c,d = 1,-5,8,0
f = a*x**4 + b*x**3 + c*x**2 + d*x + e
d = lambda g: g.subs(x, x+1) - g
dd = lambda g: d(d(g))
ddd = lambda g: d(dd(g))
dddd = lambda g: dd(dd(g))
print(expand(dd(f)).subs(x, 0))
