#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/11/07 20:33
# @Author   : Iydon
# @File     : course6.1.py


import numpy as np



def test():
    def matrix_is(A):
        # 0:no, 1:vector, 2:matrix
        c = len(A[0])
        for a in A:
            if len(a) is not c: return False
        return True

    def matrix_shape(A):
        return len(A), len(A[0])

    def matrix_init(size):
        return [[0 for j in range(size[-1])] for i in range(size[0])]

    def matrix_add(A, B):
        tA,tB = matrix_shape(A),matrix_shape(B)
        if tA != tB: return False
        C = matrix_init(tA)
        for i in range(tA[0]):
            for j in range(tA[-1]):
                C[i][j] = A[i][j] + B[i][j]
        return C

    def matrix_prod(A, B):
        tA,tB = matrix_shape(A),matrix_shape(B)
        if tA[-1] != tB[0]: return False
        C = matrix_init([tA[0],tB[-1]])
        for i in range(tA[0]):
            for j in range(tB[-1]):
                for k in range(tA[-1]):
                    C[i][j] += A[i][k] * B[k][j]
        return C

    A = matrix_init([3,2])
    B = [[1,2], [3,4], [5,6]]
    print(matrix_add(A,B))
    A = [[1,1,1], [1,1,1]]
    print(matrix_prod(A, B))

# ----------------------------

def gaussian_elimination(A, b):
    np_result = A**-1 * b
    print(np_result.T)
    n = A.shape[0]
    x = np.zeros((n,1))
    tmp = 0
    for k in range(n-1):
        for i in range(k+1,n):
            m = A[i,k] / A[k,k]
            for j in range(k,n):
                A[i,j] = A[i,j] - m*A[k,j]
            b[i,0] = b[i,0] - m*b[k,0]
    x[-1,0] = b[-1,0] / A[-1,-1]
    for i in range(n-2,-1,-1):
        for j in range(i+1,n):
            tmp += A[i,j] * x[j,0]
        x[i,0] = (b[i,0] - tmp) / A[i,i]
        tmp = 0
    return x

def gaussian_elimination_partial_pivoting(A, b):
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
            for j in range(k,n):
                A[i,j] = A[i,j] - m*A[k,j]
            b[i,0] = b[i,0] - m*b[k,0]
    x[-1,0] = b[-1,0] / A[-1,-1]
    for i in range(n-2,-1,-1):
        for j in range(i+1,n):
            tmp += A[i,j] * x[j,0]
        x[i,0] = (b[i,0] - tmp) / A[i,i]
        tmp = 0
    return x

def gaussian_elimination_scaled_partial_pivoting(A, b):
    np_result = A**-1 * b
    print(np_result.T)
    n = A.shape[0]
    x = np.zeros((n,1))
    tmp = 0
    for k in range(n):
        M = np.max(A[k,:])
        A[k,:] /= M
        b[k,0] /= M
    for k in range(n-1):
        for i in range(k+1,n):
            m = A[i,k] / A[k,k]
            for j in range(k,n):
                A[i,j] = A[i,j] - m*A[k,j]
            b[i,0] = b[i,0] - m*b[k,0]
    x[-1,0] = b[-1,0] / A[-1,-1]
    for i in range(n-2,-1,-1):
        for j in range(i+1,n):
            tmp += A[i,j] * x[j,0]
        x[i,0] = (b[i,0] - tmp) / A[i,i]
        tmp = 0
    return x

A = np.matrix([[1,2,3], [2,5,7],[3,7,11]])
b = np.matrix([[1.],[2],[3]])
result = gaussian_elimination(A, b)
print(result.T)
print("*")
A = np.matrix([[1.,2,3], [2,5,7],[3,7,11]])
b = np.matrix([[1.],[2],[3]])
result = gaussian_elimination_partial_pivoting(A, b)
print(result.T)
print("*")
A = np.matrix([[1.,2,3], [2,5,7],[3,7,11]])
b = np.matrix([[1.],[2],[3]])
result = gaussian_elimination_scaled_partial_pivoting(A, b)
print(result.T)
