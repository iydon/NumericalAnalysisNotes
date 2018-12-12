#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/11/25 14:27
# @Author   : Iydon
# @File     : 7.3.py

import numpy as np


def jacobi_method(A, b, x0=None):
    """
    A = L(ow) + D(iag) + U(pper)
    |ρ(inv(D)*(L+U))| < 1
    => convergent.
    """
    A_,b_ = A.copy(),b.copy()
    diag  = np.diag(A_)
    shap  = A.shape
    A_ = np.diag(diag) - A_
    if x0 is not None:
        x_ = x0.copy()
    else:
        x_ = np.zeros(b.shape)
    for i in range(shap[0]):
        A_[i,:] /= diag[i]
        b_[i,:] /= diag[i]
    while True:
        yield x_
        x_ = A_*x_ + b_


def gauss_seidel(A, b, x0=None):
    """
    A = L(ow) + D(iag) + U(pper)
    |ρ(inv(D-L)*U)| < 1
    => convergent.
    """
    A_,b_ = A.copy(),b.copy()
    DL_   = np.tril(A_)
    U_    = np.triu(A_, 1)
    if x0 is not None:
        x_ = x0.copy()
    else:
        x_ = np.zeros(b.shape)
    invDL_ = np.linalg.inv(DL_)
    A_ = np.matmul( invDL_, U_ )
    b_ = np.matmul( invDL_, b_ )
    while True:
        yield x_
        x_ = np.matmul(-A_, x_) + b_


A  = np.matrix([[3.,-1,1],[3,6,2],[3,3,7]])
b  = np.matrix([[1.],[0],[4]])
x0 = np.matrix([[0.],[0],[0]])
result = jacobi_method(A, b, x0)
t = 0
for r in result:
    print(t, r.T)
    t += 1
    if t > 10: break

result = gauss_seidel(A, b, x0)
t = 0
for r in result:
    print(t, r.T)
    t += 1
    if t > 10: break
