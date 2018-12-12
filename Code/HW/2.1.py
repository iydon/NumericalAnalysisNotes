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


def Bisection(fun, a:float, b:float, max_step:int=32, sign_dig:int=6) -> float:
    """
    Bisection method to find the zero.
    ----------------------------
    Args:
        fun: Function.
        a,b: Float, upper and lower bound.

    Returns:
        Tuple: (Step, zero).
        zero | fun(zero)\\sim 0.

    Raises:
        None.
    """
    fl = lambda x: round(x, sign_dig)
    eps = 10**(-sign_dig)
    mid_last = fl(a)
    if fun(a)*fun(b) < 0:
        for i in range(0, max_step):
            mid = fl((a+b) / 2)
            if abs(mid-mid_last)<eps or abs(fun(mid))<eps:
                # print("Step: %d\nZero: %f."%(i, mid))
                return (i,mid)
            else:
                if fun(mid)*fun(a)<0:
                    b = fl(mid)
                else:
                    a = fl(mid)
            mid_last = mid
        return (False, "Bisection cannot be convergent within the pre-set steps.")
    return (False, False)

def Theorem_2_1():
    """
    Theorem 2.1.
    """
    string = "Suppose that $f\\in C[a,b]$ and $f(a)\\dot{}f(b)<0$. The Bisection method " +\
             "generates a sequence $\\{p_n\\}_{n=1}^{\\infty}$ approximating a zero $p$ " +\
             "of $f$ with" +\
             "\\[\\left|p_n-p\\right|\\leq\\frac{b-a}{2^n},\\quad when\\quad n\\geq 1.\\]"
    return string

title = """
1. Use the Bisection method to find p_3 for f(x)=sqrt(x)-cos(x) on [0, 1].
"""
print(title)
f = lambda x: np.sqrt(x)-np.cos(x)
step,zero = Bisection(f, 0, 1)
print("    Step: %d\n    Zero: %f"%(step,zero))

title = """
10. Let f(x)=(x+2)(x+1)^2x(x-1)^3(x-2). To which zero of f does the Bisection method converge
    when applied on the following intervals?
"""
print(title)
f = lambda x: (x+2)*(x+1)**2*x*(x-1)**3*(x-2)
for interval in [[-1.5,2.5],[-0.5,2.4],[-0.5,3],[-3,-0.5]]:
    step,zero = Bisection(f, interval[0],interval[1])
    if step and isinstance(zero,float):
        print("    Interval: %s"%(interval))
        print("        Steps: %d\n        Zero: %f."%(step, zero))

title = """
11. Let f(x)=(x+2)(x+1)x(x-1)^3(x-2). To which zero of f does the Bisection method converge
    when applied on the following intervals?
"""
print(title)
f = lambda x: (x+2)*(x+1)*x*(x-1)**3*(x-2)
for interval in [[-3,2.5],[-2.5,3],[-1.75,1.5],[-1.5,1.75]]:
    step,zero = Bisection(f, interval[0],interval[1])
    if step and isinstance(zero,float):
        print("    Interval: %s"%(interval))
        print("        Steps: %d\n        Zero: %f."%(step, zero))

title = """
14. Use Theorem 2.1 to find a bound for the number of iterations needed to achieve an approximation
    with accuracy 1e-3 to the solution of x^3+x-4=0 lying in the interval [1,4]. Find an approximation
    to the root with this degree of accuracy.
"""
print(title)
f = lambda x: x**3+x-4
interval = [1, 4]
step,zero = Bisection(f, interval[0],interval[1], sign_dig=3)
if step and isinstance(zero,float):
    print("    Interval: %s"%(interval))
    print("        Steps: %d\n        Zero: %.3f."%(step, zero))

title = """
15. Use Theorem 2.1 to find a bound for the number of iterations needed to achieve an approximation
    with accuracy 1e-4 to the solution of x^3-x-1=0 lying in the interval [1,2]. Find an approximation
    to the root with this degree of accuracy.
"""
print(title)
f = lambda x: x**3-x-1
interval = [1, 2]
step,zero = Bisection(f, interval[0],interval[1], sign_dig=4)
if step and isinstance(zero,float):
    print("    Interval: %s"%(interval))
    print("        Steps: %d\n        Zero: %.4f."%(step, zero))

title = """
16. Let f(x)=(x-1)^{10}, and p_n=1+1/n. Show that |f(p_n)|<1e-3 whenever n>1 but that
    |p-p_n|<1e-3 requires that n>1000.
"""
print(title)
proof = """
\\[f(p_n)=n^{-10}<10^{-3}\\Rightarrow n>\\sqrt[10]{100}=1.585\\]
\\[\\left|1-(1+\\frac{1}{n})\\right|<10^{-3}\\Rightarrow n>10^3\\]
""".replace("\n", "\n    ")
print(proof)

title = """
17. Let {p_n} be the sequence defined by p_n=\\sum_{k=1}^{n}\\frac{1}{k}. Show that {p_n} diverges even though \\lim_{n\\rightarrow\\infty}(p_n-p_{n-1})=0.
"""
print(title)
proof = """
\\begin{align*}
& 1+\\frac{1}{2}+\\underbrace{\\frac{1}{3}+\\frac{1}{4}}+\\ldots + \\underbrace{\\frac{1}{2^k+1}+\\ldots+\\frac{1}{2^{k+1}}}+\\ldots \\\\
>& 1+\\frac{1}{2}+\\underbrace{\\frac{1}{4}+\\frac{1}{4}}+\\ldots + \\underbrace{\\frac{1}{2^{k+1}}+\\ldots+\\frac{1}{2^{k+1}}}+\\ldots \\\\
=& 1+\\frac{1}{2}+\\frac{1}{2}+\\ldots+\\frac{1}{2}+\\ldots
\\end{align*}
""".replace("\n", "\n    ")
print(proof)