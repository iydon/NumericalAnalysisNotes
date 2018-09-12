#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/09/12 09:01
# @Author   : Iydon
# @File     : 2.2.py
# @Abstract : Homework auxiliary code.
import numpy as np

def error(p:float, p_star:float) -> tuple:
    """
    From 1.2.py.
    """
    absolute_error = abs(p - p_star)
    relative_error = absolute_error / p
    return (absolute_error, relative_error)

def fixed_point(fun, start:float=0.0, max_step:int=128, sign_dig:int=6)->float:
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
    fl = lambda x: round(x, sign_dig)
    eps = 10**(-sign_dig)
    new_val = fun(start)
    for i in range(0, max_step):
        old_val = fl(new_val)
        new_val = fl(fun(old_val))
        if abs(old_val-new_val)<=eps:
            return (i,new_val)
    return "Max_step..."

title = """
11. For each of the following equations, determine an interval [a,b] on which fixed-point iteration will
    converge. Estimate the number of iterations necessary to obtain approximations accurate to within
    1e-5, and perform the calculations.
"""
print(title)
f = lambda x: (2-np.exp(x)+x**2)/3
df = lambda x: np.exp(x)/3+2*x/3

f = lambda x: 5/x**2+2
df = lambda x: -10*x**(-3)

f = lambda x: 5**(-x)
df = lambda x: -5**(-x)*np.log(5)

title = """
19. Use Theorem 2.4 to show that the sequence defined by
    a) x_n = \\frac{1}{2}x_{n-1}+\\frac{1}{x_{n-1}},\\quad for\\, n\\geq 1,
       converges to \\sqrt(x) whenever x_0>\\sqrt(2).
    b) Use the fact that 0<(x_0-sqrt(2))^2 whenever x_0\\neq\\sqrt(2) to show
       that if 0<x_0<\\sqrt(2)
"""
print(title)

title = """
24. Let g\\in C^1[a,b] and p be in (a,b) with g(p)=p and |g'(p)|>1. Show that there exists a \\delta >0 such
    that if 0<|p_0-p|<\\delta, then |p_0-p|<|p_1-p|. Thus, no matter how close the initial approximation
    p_0 is to p, the next iterate p_1 is farther away, so the fixed-point iteration does not converge if p_0\\neq p.
"""
print(title)