#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/28 11:41
# @Author   : Iydon
# @File     : Runge-Kutta.py

from sympy import *


def diff_n(f, x, y, n):
    dnf = [f for i in range(N)]
    for i in range(N-1):
        dnf[i+1] = diff(dnf[i],x) + diff(dnf[i],y)*f
        dnf[i+1] = simplify(dnf[i+1])
    return dnf

# variable
x,y,h = symbols("x,y,h")
# constant
N = 3
c = symbols("c%d:%d"%(1,N+1))
λ = symbols("λ%d:%d"%(1,N+1))
μ = [symbols("μ(%d;1:%d)"%(i,i)) for i in range(2,N+1)]
# function
f = Function("f")
K = [f(x,y) for i in range(N)]
for i in range(N-1):
    tmp = 0
    for j in range(i+1):
        tmp += μ[i][j]*K[j]
    K[i+1] = f(x+λ[i+1]*h, y+h*tmp)

for k in K:
    print(k)
