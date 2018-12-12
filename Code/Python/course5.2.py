#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/18 14:54
# @Author   : Iydon
# @File     : course5.2.py



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

def taylor_method(dy0, start, end, step, **args):
    """
    Approximate the solution of the initial-value
    problem.
    ----------------------------
    Args:
        args: "d%dy", like "d1y", "d2y".
    """
    xs = linspace(start, end, step=step)
    ys = [dy0 for i in xs]
    for i in range(len(xs)-1):
        tmp,pro = 0,1
        for j in range(len(args)):
            tmp += pro * args["d%dy"%(j+1)](xs[i],ys[i])
            pro *= step / (j+2)
        ys[i+1] = ys[i] + step*tmp
    return ys

def modified_euler_method(dy, dy0, start, end, step):
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
        tmp = ys[i] + step*dy(xs[i],ys[i])
        ys[i+1] = ys[i] + step*(dy(xs[i],ys[i])+dy(xs[i+1],tmp))/2
    return ys



dy = lambda t,y: y - t*t + 1
result = euler_method(dy, 0.5, 0, 2, 0.2)
print(result)

dy = lambda t,y: y - t*t + 1
result = modified_euler_method(dy, 0.5, 0, 2, 0.2)
print(result)

d1y = lambda t,y: y - t**2 + 1
d2y = lambda t,y: y - t**2 + 1 - 2*t
d3y = lambda t,y: y - t**2 - 1 - 2*t
d4y = lambda t,y: y - t**2 - 1 - 2*t
result = taylor_method(0.5, 0, 2, 0.2, d1y=d1y, d2y=d2y, d3y=d3y, d4y=d4y)
print(result)
