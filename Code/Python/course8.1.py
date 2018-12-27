#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/12/27 09:54
# @Author   : Iydon
# @File     : course8.1.py

import numpy as np
import Poly



def least_square_poly(xs, ys, n, numpy, Poly, display=False, error=False):
    m = len(xs)
    A = numpy.zeros((n+1,n+1))
    b = numpy.zeros((n+1,1))
    for i in range(n+1):
        for j in range(n+1):
            A[i,j] = sum(map(lambda x:x**(i+j), xs))
        b[i] = sum([y*x**i for x,y in zip(xs,ys)])
    x = numpy.matmul(numpy.linalg.inv(A), b)
    f = Poly.Poly(x.T.tolist()[0][::-1])
    if display:
        print(A)
        print(b.T)
        print(x.T)
    if error:
        xs_ = [f.lambdify()(x) for x in xs]
        print(sum([(x-y)**2 for x,y in zip(xs_,ys)]))
    return f


xs = [0, 0.25, 0.50, 0.75, 1.00]
ys = [1.0000, 1.2840, 1.6487, 2.1170, 2.7183]
result = least_square_poly(xs, ys, 2, np, Poly, display=True, error=True)
print(result)

print()

xs = [1.0, 1.1, 1.3, 1.5, 1.9, 2.1]
ys = [1.84, 1.96, 2.21, 2.45, 2.94, 3.18]
for i in [1,2,3]:
    result = least_square_poly(xs, ys, i, np, Poly, display=True, error=True)
    print(result)

print()

E  = 5.3
ls = [7, 9.4, 12.3]
Fs = [2, 4, 6]
ks = [Fs[i]/(ls[i]-E) for i in range(len(ls))]
result = least_square_poly(ls, ks, 0, np, Poly, display=True, error=True)
print(result)

print()

ls += [8.3, 11.3, 14.4, 15.9]
Fs += [3, 5, 8, 10]
ks = [Fs[i]/(ls[i]-E) for i in range(len(ls))]
result = least_square_poly(ls, ks, 0, np, Poly, display=True, error=True)
print(result)
