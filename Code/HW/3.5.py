#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/04 14:23
# @Author   : Iydon
# @File     : course3.5.py

import numpy as np
from Poly import *
import matplotlib.pyplot as plt

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



xs1  = [1, 2, 5, 6, 7, 8, 10, 13, 17]
fxs1 = [3.0, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5]
bou1 = [1, -0.67]
S1   = clamped_cubic_spline(xs1, fxs1, bou1)
f1   = cubic_spline_lambdify("S1", xs1)

xs2  = [17, 20, 23, 24, 25, 27, 27.7]
fxs2 = [4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1]
bou2 = [3.0, -4.0]
S2   = clamped_cubic_spline(xs2, fxs2, bou2)
f2   = cubic_spline_lambdify("S2", xs2)

xs3  = [27.7, 28, 29, 30]
fxs3 = [4.1, 4.3, 4.1, 3.0]
bou3 = [0.33, -1.5]
S3   = clamped_cubic_spline(xs3, fxs3, bou3)
f3   = cubic_spline_lambdify("S3", xs3)


# Origin
plt.scatter(xs1, fxs1, marker="*", s=7, color="red")
plt.scatter(xs2, fxs2, marker="*", s=7, color="green")
plt.scatter(xs3, fxs3, marker="*", s=7, color="blue")

# Interpolation
eps = 1e-6

x1  = np.linspace(min(xs1), max(xs1)-eps, 100)
y1  = [f1(x) for x in x1]
plt.plot(x1, y1, color="red")

x2  = np.linspace(min(xs2), max(xs2)-eps, 100)
y2  = [f2(x) for x in x2]
plt.plot(x2, y2, color="green")

x3  = np.linspace(min(xs3), max(xs3)-eps, 100)
y3  = [f3(x) for x in x3]
plt.plot(x3, y3, color="blue")


plt.axis("equal")
plt.grid(linestyle='-.')
plt.show()


xs3  = [0, 1, 2]
fxs3 = [0, 1, 2]
bou3 = [0.33, -1.5]
S3   = clamped_cubic_spline(xs3, fxs3, bou3)
f3   = cubic_spline_lambdify("S3", xs3)