#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/12/06 12:41
# @Author   : Iydon
# @File     : 9.3.py

import numpy as np
import scipy as sp


def gersgorin_circle(A, numpy):
    r"""
    |z-a_{ii}|<=\sum_{j!=i}|a_{ij}|
    """
    shape = A.shape
    if shape[0] != shape[-1]:
        raise Exception("Dimensions do not match.")
    result = [[] for i in range(shape[0])]
    for i in range(shape[0]):
        summation = numpy.sum(numpy.abs(A[i,])) - numpy.abs(A[i,i])
        result[i] = [A[i,i], summation]
    return result


def power_method(A, numpy, x0, judge=None, eps=1e-6, max_loop=64, disply=False, p=0):
    r"""
    Numpy
    B = A - pI.
    """
    if judge==None: judge = lambda x,y: numpy.linalg.norm(x-y, float("inf"))
    A_    = A.copy() - p*numpy.eye(A.__len__())
    count = 0
    mu    = numpy.max(numpy.abs(x0))
    x0    = x0 / mu
    last  = x0
    while count<max_loop:
        count += 1
        uv = numpy.matmul(A_, x0)
        mu = numpy.max(numpy.abs(uv)) 
        x0 = uv / mu
        if disply:
            print("%%%dd: %%s"%(len(str(max_loop)))%(count, mu))
        if judge(last,x0)<eps:
            break
        last = x0
    return mu + p


def inverse_power_method(A, scipy, x0, judge=None, eps=1e-6, max_loop=64, disply=False, p=0):
    r"""
    Scipy
    pass
    """
    from scipy import linalg
    if judge==None:
        norm  = lambda x: scipy.linalg.norm(x, float("inf"))
        judge = lambda x,y: abs(norm(x-y))
    A_    = A.copy() - p*scipy.eye(len(A))
    LU    = linalg.lu_factor(A_)
    count = 0
    mu    = max(abs(x0))
    x0    = x0 / mu
    last  = x0
    while count<max_loop:
        count += 1
        uv = linalg.lu_solve(LU, x0)
        mu = max(abs(x0))[0]
        x0 = uv / mu
        if disply:
            print("%%%dd: %%s"%(len(str(max_loop)))%(count, p+1/mu))
        if judge(last,x0)<eps:
            break
        last = x0
    return p + 1/mu


    

A  = np.matrix([[1,1,1],[1,1,0],[1,0,1]])
x0 = np.matrix([-1,0,1]).T
result = power_method(A, np, x0=x0, judge=None, eps=1e-6, max_loop=64, disply=True, p=0)
print(result)

result = inverse_power_method(A, np, x0=x0, judge=None, eps=1e-6, max_loop=64, disply=True, p=0)
print(result)

print(np.max(np.linalg.eigvals(A)))
print(np.min(np.linalg.eigvals(A)))



A  = np.matrix([[2,1,1],[1,2,1],[1,1,2]])
x0 = np.matrix([1,-1,2]).T
result = power_method(A, np, x0=x0, judge=None, eps=1e-6, max_loop=64, disply=True, p=0)
print(result)
result = inverse_power_method(A, np, x0=x0, judge=None, eps=1e-6, max_loop=64, disply=True, p=0)
print(result)

print(np.max(np.linalg.eigvals(A)))
print(np.min(np.linalg.eigvals(A)))

print("\n"*10)

A = np.matrix([[0,0,2,4],[1/2,0,0,0],[0,1/4,0,0],[0,0,1/8,0]])
x0 = np.matrix([[1],[1],[1],[1]])
result = power_method(A, np, x0=x0, judge=None, eps=1e-6, max_loop=128, disply=True, p=0)
print(result)
print(np.max(np.linalg.eigvals(A)))

print("\n"*10)

M = 11
a = 3/4
A = np.diag([(1+2*a) for i in range(M-1)]) - np.diag([a for i in range(M-2)],1) - np.diag([a for i in range(M-2)],-1)
x0 = np.ones((M-1,1))
result = inverse_power_method(A, np, x0=x0, judge=None, eps=1e-6, max_loop=128, disply=True, p=0)

print(result)
print(np.min(np.linalg.eigvals(A)))

print("\n"*10)

M = 11
a = 3/4
A = np.diag([(1+a) for i in range(M-1)]) - np.diag([-a/2 for i in range(M-2)],1) - np.diag([-a/2 for i in range(M-2)],-1)
B = np.diag([(1+a) for i in range(M-1)]) - np.diag([a/2 for i in range(M-2)],1) - np.diag([a/2 for i in range(M-2)],-1)
C = np.matmul( np.linalg.inv(A), B )
x0 = np.ones((M-1,1))
result = power_method(C, np, x0=x0, judge=None, eps=1e-6, max_loop=1024, disply=True, p=0)

print(result)
print(np.max(np.linalg.eigvals(C)))
