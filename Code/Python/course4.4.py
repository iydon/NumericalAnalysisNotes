    #!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/16 19:48
# @Author   : Iydon
# @File     : course4.4.py


import numpy as np


# Composite Numerical Integration
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

def midpoint_rule(xs:list, fxs:list) -> float:
    return (xs[-1]-xs[0]) * fxs[1]

def composite_trapezoidal_rule(xs:list, fxs:list) -> float:
    result = 0
    for i in range(len(xs)-1):
        result += trapezoidal_rule(xs[i:i+2], fxs[i:i+2])
    return result

def composite_simpson_rule(xs:list, fxs:list) -> float:
    result = 0
    for i in range(len(xs)//2):
        result += simpson_rule(xs[2*i:2*i+3], fxs[2*i:2*i+3])
    return result

def composite_midpoint_rule(xs:list, fxs:list) -> float:
    result = 0
    for i in range(len(xs)//2):
        result += midpoint_rule(xs[2*i:2*i+3], fxs[2*i:2*i+3])
    return result




f  = np.exp
xs = [i/2 for i in range(9)]
print(composite_simpson_rule(xs, [f(x) for x in xs]))
print(composite_trapezoidal_rule(xs, [f(x) for x in xs]))
print(composite_midpoint_rule(xs, [f(x) for x in xs]))

pi = np.pi
f  = np.sin
xs = [i*pi/359 for i in range(360)]
print(composite_midpoint_rule(xs, [f(x) for x in xs]))

pi = np.pi
f  = np.sin
xs = [i*pi/20 for i in range(21)]
print(composite_simpson_rule(xs, [f(x) for x in xs]))
print(composite_trapezoidal_rule(xs, [f(x) for x in xs]))
print(composite_midpoint_rule(xs, [f(x) for x in xs]))

f  = lambda x: x**3*np.exp(x)
xs = [-2, -1, 0, 1, 2]
print(composite_trapezoidal_rule(xs, [f(x) for x in xs]))
xs = [i/2 for i in range(-4,5)]
print(composite_simpson_rule(xs, [f(x) for x in xs]))

f  = lambda x: np.exp(2*x)*np.sin(3*x)
xs = [2*i/2169 for i in range(2170)]
print(composite_trapezoidal_rule(xs, [f(x) for x in xs]))
xs = [2*i/54 for i in range(55)]
print(composite_simpson_rule(xs, [f(x) for x in xs]))
 
print(-14.213977129862521)
