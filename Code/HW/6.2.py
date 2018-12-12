#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/11/09 19:02
# @Author   : Iydon
# @File     : 6.1.py


import numpy as np
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)


def chopping(num, dig):
    eps = 10**-dig / 2
    return np.round(num-eps, dig)


def gaussian_elimination(A, b, num_dig):
    """
    Chopping.
    """
    np_result = A**-1 * b
    print(np_result.T)
    n = A.shape[0]
    x = np.zeros((n,1))
    tmp = 0
    for k in range(n-1):
        for i in range(k+1,n):
            m = A[i,k] / A[k,k]
            A[i,:] = chopping( A[i,:] - m*A[k,:], num_dig )
            b[i,0] = chopping( b[i,0] - m*b[k,0], num_dig )
        print(np.hstack((A,b)))
    x[-1,0] = chopping( b[-1,0] / A[-1,-1], num_dig )
    for i in range(n-2,-1,-1):
        for j in range(i+1,n):
            tmp += chopping( A[i,j] * x[j,0], num_dig )
        x[i,0] = chopping( (b[i,0] - tmp) / A[i,i], num_dig )
        tmp = 0
    return x


def gaussian_elimination(A, b, num_dig):
    """
    Rounding.
    """
    np_result = A**-1 * b
    print(np_result.T)
    n = A.shape[0]
    x = np.zeros((n,1))
    tmp = 0
    for k in range(n-1):
        for i in range(k+1,n):
            m = A[i,k] / A[k,k]
            A[i,:] = np.round( A[i,:] - m*A[k,:], num_dig )
            b[i,0] = np.round( b[i,0] - m*b[k,0], num_dig )
        print(np.hstack((A,b)))
    x[-1,0] = np.round( b[-1,0] / A[-1,-1], num_dig )
    for i in range(n-2,-1,-1):
        for j in range(i+1,n):
            tmp += np.round( A[i,j] * x[j,0], num_dig )
        x[i,0] = np.round( (b[i,0] - tmp) / A[i,i], num_dig )
        tmp = 0
    return x


def gaussian_elimination_partial_pivoting(A, b, num_dig):
    """
    Rounding.
    """
    np_result = A**-1 * b
    print(np_result.T)
    n = A.shape[0]
    x = np.zeros((n,1))
    tmp = 0
    for k in range(n-1):
        M = k
        for m in range(k+1,n):
            if A[m,k] > A[M,k]: M = m
        A[[k,M]] = A[[M,k]]
        b[[k,M]] = b[[M,k]]
        for i in range(k+1,n):
            m = A[i,k] / A[k,k]
            A[i,:] = np.round( A[i,:] - m*A[k,:], num_dig )
            b[i,0] = np.round( b[i,0] - m*b[k,0], num_dig )
        print(np.hstack((A,b)))
    x[-1,0] = np.round( b[-1,0] / A[-1,-1], num_dig )
    for i in range(n-2,-1,-1):
        for j in range(i+1,n):
            tmp += np.round( A[i,j] * x[j,0], num_dig )
        x[i,0] = np.round( (b[i,0] - tmp) / A[i,i], num_dig )
        tmp = 0
    return x



A = np.matrix([[0.03, 58.9],[5.31, -6.10]])
b = np.matrix([[59.2], [47.0]])
result = gaussian_elimination_partial_pivoting(A, b, num_dig=3)
print(result.T)

print(np.linalg.norm(result-np.matrix([[10],[1]])))
