#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/31 17:03
# @Author   : Iydon
# @File     : 5.9.py


import numpy as np

def linspace(start, end, step=None):
    """
    Args:
        step: start:step:end
    """
    length = (end-start)/step
    return [start+i*step for i in range(int(length[0,0])+1)]

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



t1,t2 = lambda t: -(2*t*t+1)*np.exp(2*t), lambda t: (t*t+2*t-4)*np.exp(2*t)
dy    = lambda t,y: np.matrix([[t1(t[0,0])],[t2(t[1,0])]]) + np.matrix([[3,2],[4,1]])*y
dy0   = np.matrix([[1],[1]])
start = np.matrix([[0],[0]])
end   = np.matrix([[1],[1]])
step  = 0.2
result = runge_kutta(dy, dy0, start, end, step, order=4)
print(result)
u1 = lambda t: np.exp(5*t)/3-np.exp(-t)/3+np.exp(2*t)
u2 = lambda t: np.exp(5*t)/3+2*np.exp(-t)/3+t*t*np.exp(2*t)
t  = [0, 0.2, 0.4, 0.6, 0.8, 1]
print([u1(i) for i in t])
print([u2(i) for i in t])


t1,t2 = lambda t: np.cos(t)+4*np.sin(t), lambda t: -3*np.sin(t)
dy    = lambda t,y: np.matrix([[t1(t[0,0])],[t2(t[1,0])]]) + np.matrix([[-4,-2],[3,1]])*y
dy0   = np.matrix([[0],[-1]])
start = np.matrix([[0],[0]])
end   = np.matrix([[2],[2]])
step  = 0.1
result = runge_kutta(dy, dy0, start, end, step, order=4)
print(result)
u1 = lambda t: 2*np.exp(-t)-2*np.exp(-2*t)+np.sin(t)
u2 = lambda t: -3*np.exp(-t)+2*np.exp(-2*t)
t  = [i/10 for i in range(21)]
print([u1(i) for i in t])
print([u2(i) for i in t])
