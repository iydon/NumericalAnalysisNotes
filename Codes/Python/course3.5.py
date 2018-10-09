#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/04 14:23
# @Author   : Iydon
# @File     : course3.5.py

import numpy as np
from Poly import *
import matplotlib.pyplot as plt

def natural_cubic_spline(xs:list, fxs:list, display:bool=False):
	"""
	Cubic spline interpolation.
	"""
	n  = len(xs)
	hs = [xs[i+1]-xs[i] for i in range(n-1)]
	A  = np.diag([0]+hs[1:],1) + np.diag(hs[:-1]+[0],-1)
	A += np.diag([1]+[2*(hs[i+1]+hs[i]) for i in range(n-2)]+[1])
	bs = [0]+[3/hs[i+1]*(fxs[i+2]-fxs[i+1])-3/hs[i]*(fxs[i+1]-fxs[i]) for i in range(n-2)]+[0]
	# a, b, c, d: end with 'x'.
	cx = [i[0] for i in (np.linalg.inv(A) * np.matrix(bs).transpose()).tolist()]
	bx = [1/hs[i]*(fxs[i+1]-fxs[i])-hs[i]/3*(2*cx[i]+cx[i+1]) for i in range(n-1)]
	dx = [1/3/hs[i]*(cx[i+1]-cx[i]) for i in range(n-1)]
	# S_i(x)
	Ss = [fxs[i]+bx[i]*Poly([1,-xs[i]])+cx[i]*Poly([1,-xs[i]])**2+dx[i]*Poly([1,-xs[i]])**3 for i in range(n-1)]
	if display: print(fxs, bx, cx, dx, sep="\n\n\n")
	return Ss

def clamped_cubic_spline(xs:list, fxs:list, boundray:list=[0,0]):
	"""
	Cubic spline interpolation.
	"""
	n  = len(xs)
	hs = [xs[i+1]-xs[i] for i in range(n-1)]
	A  = np.diag(hs,1) + np.diag(hs,-1)
	A += np.diag([2*hs[0]]+[2*(hs[i+1]+hs[i]) for i in range(n-2)]+[2*hs[-1]])
	head = [3/hs[0]*(fxs[1]-fxs[0]) - 3*boundray[0]]
	tail = [3*boundray[-1] - 3/hs[-1]*(fxs[-1]-fxs[-2])]
	bs = head+[3/hs[i+1]*(fxs[i+2]-fxs[i+1])-3/hs[i]*(fxs[i+1]-fxs[i]) for i in range(n-2)]+tail
	# a, b, c, d: end with 'x'.
	cx = [i[0] for i in (np.linalg.inv(A) * np.matrix(bs).transpose()).tolist()]
	bx = [1/hs[i]*(fxs[i+1]-fxs[i])-hs[i]/3*(2*cx[i]+cx[i+1]) for i in range(n-1)]
	dx = [1/3/hs[i]*(cx[i+1]-cx[i]) for i in range(n-1)]
	# S_i(x)
	Ss = [fxs[i]+bx[i]*Poly([1,-xs[i]])+cx[i]*Poly([1,-xs[i]])**2+dx[i]*Poly([1,-xs[i]])**3 for i in range(n-1)]
	return Ss

def cubic_spline_lambdify(S:str, xs:list):
	"""
	Lambdify the cubic spline function.
	"""
	f = ["%s[%d].lambdify()(x)*(%s<=x<%s)"%(S, i, xs[i], xs[i+1]) for i in range(len(xs)-1)]
	return eval("lambda x: %s"%"+".join(f))

xs  = [0.9,1.3,1.9,2.1,2.6,3.0,3.9,4.4,4.7,5.0,6.0,7.0,8.0,9.2,10.5,11.3,11.6,12.0,12.6,13.0,13.3]
fxs = [1.3,1.5,1.85,2.1,2.6,2.7,2.4,2.15,2.05,2.1,2.25,2.3,2.25,1.95,1.4,0.9,0.7,0.6,0.5,0.4,0.25]

S = natural_cubic_spline(xs, fxs)

f = cubic_spline_lambdify("S", xs)
plt.plot(xs, fxs, marker="*", color="orange")

x = np.linspace(0.9, 13.29, 100)
y = [f(x) for x in x]
plt.plot(x, y, color="blue")
plt.axis("equal")
plt.grid()
plt.show()
