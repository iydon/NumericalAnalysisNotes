#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/19 15:14
# @Author   : Iydon
# @File     : course5.4.py


# See `course5.2.py'.

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


dy = lambda t,y: y-2*t/y
result = runge_kutta(dy, 1, 0, 1, 0.2, order=4)
print(result)