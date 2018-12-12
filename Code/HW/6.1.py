#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/11/09 19:02
# @Author   : Iydon
# @File     : 6.1.py


import numpy as np


def gauss_jordan_method(A, b, num_dig:int):
    np_result = A**-1 * b
    print(np_result.T)
    n = A.shape[0]
    x = np.zeros((n,1))
    for k in range(n):
        for i in range(n):
            if k == i:
                continue
            m = A[i,k] / A[k,k]
            A[i,:] = np.round( A[i,:] - m*A[k,:], num_dig )
            b[i,0] = np.round( b[i,0] - m*b[k,0], num_dig )
            print(np.hstack((A,b)))
            print()
    for i in range(n):
    	x[i,0] = np.round( b[i,0]/A[i,i], num_dig )
    return x


A = np.matrix([[4.,-1,1],[2,5,2],[1,2,4]])
b = np.matrix([[8.],[3],[11]])
result = gauss_jordan_method(A, b, num_dig=2)
print(result.T)
