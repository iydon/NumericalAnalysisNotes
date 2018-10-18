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
    xs = linspace(start,end,step=step)
    ys = [dy0 for i in xs]
    for i in range(len(xs)-1):
        ys[i+1] = ys[i] + step*(dy(xs[i],ys[i]))
    return ys

dy = lambda t,y: y - t*t + 1
result = euler_method(dy, 0.5, 0, 2, 0.2)
print(result)
