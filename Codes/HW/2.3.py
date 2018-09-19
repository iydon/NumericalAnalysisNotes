#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/09/19 15:10
# @Author   : Iydon
# @File     : 2.3.py
# @Abstract : Homework auxiliary code.
import numpy as np


def Newton_method(f, df, start:float=0.0, step:int=1)->float:
	"""
	From 2.2.py.
	"""
	fun = lambda x: x-f(x)/df(x)
	if step:
		val = fun(start)
		for i in range(1,step):
			val = fun(val)
		return val
	else:
		return start

def secant_method(f, start:list, step:int=1) -> float:
    """
    From 3.py.
    """
    p = [i    for i in start]
    q = [f(i) for i in start]
    for i in range(step):
    	_p = p[-1] - q[-1]*(p[-1]-p[0])/(q[-1]-q[0])
    	p = [p[-1], _p]
    	q = [q[-1], f(_p)]
    return _p

def false_position(f, start:list, step:int=1) -> float:
    """
    From 3.py.
    """
    p = [i    for i in start]
    q = [f(i) for i in start]
    for i in range(step):
    	_p = p[-1] - q[-1]*(p[-1]-p[0])/(q[-1]-q[0])
    	_q = f(_p)
    	if _q*q[-1] < 0:
    		p[0] = p[-1]
    		q[0] = q[-1]
    	p[-1] = _p
    	q[-1] = _q
    return False

# ------- HOMEWORK -------
f  = lambda x: x**2-6
df = lambda x: 2*x
result = Newton_method(f, df, 1, 2)
print(result)

f  = lambda x: -x**3-np.cos(x)
df = lambda x: -3*x**2+np.sin(x)
result = Newton_method(f, df, -1, 2)
print(result)

f = lambda x: -x**3-np.cos(x)
result = secant_method(f, [-1,0], 2)
print(result)

