#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/12/04 15:05
# @Author   : Iydon
# @File     : course9.py

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


def power_method(A, numpy, x0=None, judge=None, eps=1e-6, max_loop=64, disply=False, p=0):
    r"""
    Numpy
    B = A - pI.
    """
    if x0==None:    x0 = np.ones((A.shape[0], 1))
    if judge==None: judge = lambda x,y: numpy.linalg.norm(x-y, float("inf"))
    A_    = A.copy() - p*numpy.eye(A.__len__())
    count = 0
    mu    = numpy.max(x0)
    x0    = x0 / mu
    last  = x0
    while count<max_loop:
        count += 1
        uv = numpy.matmul(A_, x0)
        mu = numpy.max(uv)
        x0 = uv / mu
        if disply:
            print("%%%dd: %%s"%(len(str(max_loop)))%(count, mu))
        if judge(last,x0)<eps:
            break
        last = x0
    return mu + p


def inverse_power_method(A, scipy, x0=None, judge=None, eps=1e-6, max_loop=64, disply=False, p=0):
    r"""
    Scipy
    pass
    """
    from scipy import linalg
    if x0==None:    x0 = np.ones((A.shape[0], 1))
    if judge==None:
        norm  = lambda x: scipy.linalg.norm(x, float("inf"))
        judge = lambda x,y: abs(norm(x-y))
    A_    = A.copy() - p*scipy.eye(len(A))
    LU    = linalg.lu_factor(A_)
    count = 0
    mu    = max(x0)
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


    


A = np.matrix([[0.25,1,0.5],[1,0.25,0.25],[0.5,0.25,1.25]])
print(A)
result = gersgorin_circle(A, np)
print(result)

result = power_method(A, np, disply=True, p=0)
print(result.T)
print(np.max(np.linalg.eigvals(A)))

A = sp.matrix([[1,2,4],[1,3,9],[1,5,1]])
result = inverse_power_method(A, sp, disply=True, eps=1e-8, max_loop=64, p=-7)
print(result.T)
print(np.min(np.linalg.eigvals(A)))
