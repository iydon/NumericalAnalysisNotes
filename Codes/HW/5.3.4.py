#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/23 14:02
# @Author   : Iydon
# @File     : 5.3.4.py

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

def runge_kutta(dy, dy0, start, end, step, order:int=2):
    """
    Runge-Kutta Methods of Order Two.
    """
    if order not in [1,2,3,4]:
        raise ValueError
    xs = linspace(start, end, step=step)
    ys = [dy0 for i in xs]
    for i in range(len(xs)-1):
        k1 = step * dy(xs[i],ys[i])
        if order<2:
            ys[i+1] = ys[i] + k1
        elif order<3:
            k2 = step * dy(xs[i]+step/2,ys[i]+k1/2)
            ys[i+1] = ys[i] + k2
        elif order<4:
            k2 = step * dy(xs[i]+step/2,ys[i]+k1/2)
            k3 = step * dy(xs[i]+step,ys[i]-k1+2*k2)
            ys[i+1] = ys[i] + (k1+4*k2+k3)/6
        elif order<5:
            k2 = step * dy(xs[i]+step/2,ys[i]+k1/2)
            k3 = step * dy(xs[i]+step/2,ys[i]+k2/2)
            k4 = step * dy(xs[i]+step,ys[i]+k3)
            ys[i+1] = ys[i] + (k1+2*k2+2*k3+k4)/6
    return ys



d1y = lambda t,y: t*np.exp(3*t)-2*y
d2y = lambda t,y: (t+1)*np.exp(3*t)-4*y
result = taylor_method(0, 0, 1, 0.1, d1y=d1y, d2y=d2y)
print(result)

print("*"*ord("*"))

y = lambda t: t*np.exp(3*t)/5 - np.exp(3*t)/25 + np.exp(-2*t)/25
true_val = [y(i) for i in [0, 0.5, 1]]
print(true_val)
true_val = np.array(true_val)

dy = lambda t,y: t*np.exp(3*t) - 2*y
result = modified_euler_method(dy, 0, 0, 1, 0.5)
print(result, sum(abs(true_val-result)))

result = runge_kutta(dy, 0, 0, 1, 0.5, 2)
print(result, sum(abs(true_val-result)))

result = runge_kutta(dy, 0, 0, 1, 0.5, 4)
print(result, sum(abs(true_val-result)))
