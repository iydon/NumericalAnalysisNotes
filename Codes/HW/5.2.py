#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/23 14:02
# @Author   : Iydon
# @File     : 5.2.py


import matplotlib.pyplot as plt
import numpy as np


def linspace(start, end, num=None, step=None):
    """
    Args:
        num:  numpy.linspace(start, end, num)
        step: start:step:end
    """
    if num:
        sepera = (end-start) / (num-1)
        return [start+i*sepera for i in range(num)]
    elif step:
        length = int((end-start)/step)
        return [start+i*step for i in range(length+1)]
    else:
        return [start, end]


def euler_method(dy, dy0, start, end, step):
    """
    Approximate the solution of the initial-value
    problem.
    | y' = f(t,y), a<=t<=b
    | y(a) = /alpha
    ----------------------------
    Args:
        dy:    y'.
        dy0:   /alpha.
        start: a.
        end:   b.
        step:  h.
    """
    xs = linspace(start, end, step=step)
    ys = [dy0 for i in xs]
    for i in range(len(xs)-1):
        ys[i+1] = ys[i] + step*(dy(xs[i],ys[i]))
    return ys

dy = lambda t,y: t*np.exp(3*t)-2*y
result = euler_method(dy, 0, 0, 1, 0.5)
print(result)

dy = lambda t,y: 1+(t-y)**2
result = euler_method(dy, 1, 2, 3, 0.5)
print(result)

dy = lambda t,y: np.exp(t-y)
result = euler_method(dy, 1, 0, 1, 0.5)
print(result)

dy = lambda t,y: (1+t) / (1+y)
result = euler_method(dy, 2, 1, 2, 0.5)
print(result)

dy = lambda t,y: -10*y
result = euler_method(dy, 1, 0, 2, 0.1)
print(result)

print("*")
ans = 1-np.exp(-1)
dy = lambda t,y: -y+1
h = 0.1
result = euler_method(dy, 0, 0, 1, h)
print(result[-1]-ans)
h = 0.01
result = euler_method(dy, 0, 0, 1, h)
print(result[-1]-ans)

print(ans)
