#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/31 22:30
# @Author   : Iydon
# @File     : first_order_linear_differential_equation.py

from sympy import *

def solver(P, Q, x, y, x0, y0):
	intP  = integrate(P, x)
	right = integrate(Q*exp(intP), x)
	left  = exp(-intP)
	# Calculate `C'.
	C = y0/left - right
	return left * (right + C.subs(x, x0))


x,y = symbols("x,y")
P = 2
Q = x * exp(3*x)
result = solver(P, Q, x, y, 0, 0)
print(result)
