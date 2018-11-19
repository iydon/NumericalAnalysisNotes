#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/11/19 15:30
# @Author   : Iydon
# @File     : course7.3.py

import numpy as np


def jacobi_method(A, b):
	A_,b_ = A.copy(),b.copy()
	diag  = np.diag(A)
	shap  = A.shape
	A_ = np.diag(diag) - A_
	x_ = np.zeros(b.shape)
	for i in range(shap[0]):
		A_[i,:] /= diag[i]
		b_[i,:] /= diag[i]
	while True:
		yield x_
		x_ = A_*x_ + b_


A = np.matrix([[10.,-1,2,0],[-1,11,-1,3],[2,-1,10,-1],[0,3,-1,8]])
b = np.matrix([[6.],[25],[-11],[15]])
result = jacobi_method(A, b)
t = 0
for r in result:
	print(t, r.T, "\n")
	t += 1
	if t > 10: break
