#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/28 11:41
# @Author   : Iydon
# @File     : Runge-Kutta-N.py

from sympy import *

# variable
x,h = symbols("x,h")
# constant
N = 3
c = symbols("c%d:%d"%(1,N+1))
λ = symbols("λ%d:%d"%(1,N+1))
μ = [symbols("μ(%d;1:%d)"%(i,i)) for i in range(2,N+1)]
# function
y = Function("y") # y(x)
f = diff(y(x), x) # f(x,y)
K = [f for i in range(N)]
for i in range(N-1):
    result = [μ[i][j]*K[j] for j in range(i+1)]
    K[i+1] = f.subs({x:x+λ[i+1]*h, y(x):y(x)+h*sum(result)})
    # K[i+1] = Subs(f, (x,y(x)), (x+λ[i+1]*h,y(x)+h*sum(result)))
    # K[i+1] = f.subs(x,x+λ[i+1]*h).subs(y(x),y(x)+h*sum(result))

result = 0
factor = 1
for i in range(1,N+1):
    factor *= i
    result += h**i/factor*diff(y(x),x,i)
    result -= h*c[i-1]*K[i-1]


for k in K:
    print(pretty(k), "\n"*3)

__import__("os").system("pause")