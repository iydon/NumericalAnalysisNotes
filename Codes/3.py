#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/09/13 14:25
# @Author   : Iydon
# @File     : 3.py
import numpy as np
import time, random, threading


def tic(): globals()['TIC_TIME'] = time.time()
toc = lambda :time.time()-globals()['TIC_TIME']


def Newton_method(f, df, start:float=0.0, max_step:int=32, sign_dig:int=6)->float:
    """
    Newton method.
    ---------------------------
    Args:
        None.

    Returns:
        None.

    Raises:
        None.
    """
    fun = lambda x: x - f(x)/df(x)
    return fixed_point(fun, start, max_step, sign_dig)


def fixed_point(fun, start:float, max_step:int, sign_dig:int)->float:
    """
    Fixed-point iteration to find the zero.
    ----------------------------
    Args:
        fun: Function.
        start: Float, the first iteration point.
        max_step: Integer, max number of iteration.
        sign_dig: Integer, significant digits.

    Returns:
        Float zero.
        zero | fun(zero)\\sim 0.

    Raises:
        None.
    """
    fl = lambda x: round(x, 100)
    eps = 10**(-sign_dig)
    new_val = fun(start)
    for i in range(0, max_step):
        old_val = fl(new_val)
        new_val = fl(fun(old_val))
        if abs(old_val-new_val)<=2*eps:
            return (i, new_val)
    return "Max_step..."


def secant_method(f, start:list, max_step:int=32, eps:float=1e-6) -> float:
    """
    Secant method.
    ----------------------------
    Args:
        f: Function.
        start: List of float, the first iteration point.
        max_step: Integer, max number of iteration.

    Returns:
        Float zero.
        zero | fun(zero)\\sim 0.

    Raises:
        None.
    """
    p = [i    for i in start]
    q = [f(i) for i in start]
    for i in range(max_step):
    	_p = p[-1] - q[-1]*(p[-1]-p[0])/(q[-1]-q[0])
    	if abs(_p-p[-1])<eps:
    		return (i, _p)
    	p = [p[-1], _p]
    	q = [q[-1], f(_p)]
    return False


def false_position(f, start:list, max_step:int=32, eps:float=1e-6) -> float:
    """
    False position.
    ----------------------------
    Args:
        f: Function.
        start: List of float, the first iteration point.
        max_step: Integer, max number of iteration.

    Returns:
        Float zero.
        zero | fun(zero)\\sim 0.

    Raises:
        None.
    """
    p = [i    for i in start]
    q = [f(i) for i in start]
    for i in range(max_step):
    	_p = p[-1] - q[-1]*(p[-1]-p[0])/(q[-1]-q[0])
    	if abs(_p-p[-1]) < eps:
    		return (i, _p)
    	_q = f(_p)
    	if _q*q[-1] < 0:
    		p[0] = p[-1]
    		q[0] = q[-1]
    	p[-1] = _p
    	q[-1] = _q
    return False


def diff(f, eps:float=1e-6, left_flag=True):
    """
    \\frac{f(x+eps)-f(x)}{eps}.
    """
    if isinstance(f, type(lambda:0)):
        if left_flag:
            return lambda x: (f(x+eps)-f(x))/eps
        else:
            return lambda x: (f(x)-f(x-eps))/eps





f = lambda x: np.cos(x)
df = lambda x: -np.sin(x)
result = Newton_method(f, df, 9)
print(result)

df = diff(f, eps=1e-3)
result = Newton_method(f, df, 9)
print(result)

result = secant_method(f, [0,1])
print(result)

result = false_position(f, [0,1])
print(result)
