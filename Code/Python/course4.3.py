#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/15 21:09
# @Author   : Iydon
# @File     : course4.3.py


import numpy as np
from Poly import Poly


def numerical_quadrature(xs:list, fxs:list, ab:list=[]):
    # Lagrange element
    def Lagrange(i):
        def prod(lst, result=1):
            return prod(lst[:-1], result*lst[-1]) if lst else result
        num = prod([Poly([1,-x]) for x in xs[:i]+xs[i+1:]])
        den = prod([(xs[i]-x) for x in xs[:i]+xs[i+1:]])
        return num / den
    # function main
    if not ab:
        ab = [xs[0], xs[-1]]
    result = 0
    for i in range(len(xs)):
        F = (Lagrange(i)*Poly([1])).integrate().lambdify()
        result += (F(ab[-1])-F(ab[0])) * fxs[i]
    return result

def trapezoidal_rule(xs:list, fxs:list) -> float:
    """
    Trapezoidal Rule:
    /int_{x_0}^{x_1}f(x)dx = (h/2)[f(x_0)+f(x_1)] - E(f)
    E(f) = (h^3/12)f''(/xi)
    """
    return (xs[-1]-xs[0]) / 2 * sum(fxs)

def simpson_rule(xs:list, fxs:list) -> float:
    """
    Simpson's Rule:
    /int_{x_0}^{x_2}f(x)dx = (h/3)[f(x_0)+4f(x_1)+f(x_2)] - E(f)
    E(f) = (h^5/90)f^{(4)}(/xi)
    """
    return (xs[-1]-xs[0]) / 6 * (fxs[0]+4*fxs[1]+fxs[2])

def closed_newton_cotes(xs, fxs):
    """
    len(xs)==len(fxs)
    """
    return numerical_quadrature(xs, fxs, [xs[0],xs[-1]])

def open_newton_cotes(xs, fxs):
    """
    len(xs)==len(fxs)+2
    """
    return numerical_quadrature(xs[1:-1], fxs, [xs[0],xs[-1]])



# numberical quadrature
r = numerical_quadrature([1,2,3], [1, 8, 27])
print(r)

# trapezoidal and simpson
f = lambda x: x*x
xs = [0, 2]
print(trapezoidal_rule(xs, [f(x) for x in xs]))
xs = [0, 1, 2]
print(simpson_rule(xs, [f(x) for x in xs]))

# closed newton-cotes formula
pi = np.pi
f  = np.sin
xs = [0, pi/4]
print(closed_newton_cotes(xs, [f(x) for x in xs]))
xs = [0, pi/8, pi/4]
print(closed_newton_cotes(xs, [f(x) for x in xs]))
xs = [0, pi/12, pi/6, pi/4]
print(closed_newton_cotes(xs, [f(x) for x in xs]))
xs = [0, pi/16, pi/8, 3*pi/16, pi/4]
print(closed_newton_cotes(xs, [f(x) for x in xs]))

# open newton-cotes formula
pi = np.pi
f  = np.sin
xs = [0, pi/8, pi/4]
print(open_newton_cotes(xs, [f(x) for x in xs[1:-1]]))
xs = [0, pi/12, pi/6, pi/4]
print(open_newton_cotes(xs, [f(x) for x in xs[1:-1]]))
xs = [0, pi/16, pi/8, 3*pi/16, pi/4]
print(open_newton_cotes(xs, [f(x) for x in xs[1:-1]]))
xs = [0, pi/20, pi/10, 3*pi/20, pi/5, pi/4]
print(open_newton_cotes(xs, [f(x) for x in xs[1:-1]]))

# exercise
f  = lambda x: x**4
xs = [0.5, 1]
print(trapezoidal_rule(xs, [f(x) for x in xs]))
xs = [0.5, 0.75, 1]
print(simpson_rule(xs, [f(x) for x in xs]))

pi = np.pi
f  = lambda x: x*np.sin(x)
xs = [0, pi/4]
print(trapezoidal_rule(xs, [f(x) for x in xs]))
xs = [0, pi/8, pi/4]
print(simpson_rule(xs, [f(x) for x in xs]))

xs  = [0, 2]
fxs = [2, 2]
print(trapezoidal_rule(xs, fxs))
xs  = [0, 1, 2]
fxs = [2, 1/2, 2]
print(simpson_rule(xs, fxs))
