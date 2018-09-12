#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/09/12 09:01
# @Author   : Iydon
# @File     : 2.1.py
# @Abstract : Homework auxiliary code.
import numpy as np


def error(p:float, p_star:float) -> tuple:
    """
    From 1.2.py.
    """
    absolute_error = abs(p - p_star)
    relative_error = absolute_error / p
    return (absolute_error, relative_error)


def Bisection(fun, a:float, b:float, max_step:int=128, eps:float=1e-6) -> float:
    """
    Bisection method to find the zero.
    ----------------------------
    Args:
        fun: Function.
        a,b: Float, upper and lower bound.

    Returns:
        x | fun(x)\\sim 0.

    Raises:
        None.
    """
    mid_last = a
    if fun(a)*fun(b) < 0:
        for i in range(0, max_step):
            mid = (a+b) / 2
            if abs(mid-mid_last)<eps or abs(fun(mid))<eps:
                _error = error(fun(mid),0)
                print("Step: %d\nZero: %f.\nAbsolute error: %f\n%f"\
                    %(i, mid, _error[0], _error[1]))
                return mid
            else:
                if fun(mid)*fun(a)<0:
                    b = mid
                else:
                    a = mid
            mid_last = mid
        print('Bisection cannot be convergent within the pre-set steps.')



title = """
1. Use the Bisection method to find p_3 for f(x)=sqrt(x)-cos(x) on [0, 1].
"""
print(title)
f = lambda x: np.sqrt(x)-np.cos(x)
print(Bisection(f, 0, 1))