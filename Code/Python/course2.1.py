# -*- coding=utf-8 -*-
import math


def Bisection(fun, a:float, b:float, max_step:int=128, eps:float=1e-6)->float:
    """
    pass
    """
    mid_last = a
    if fun(a)*fun(b) < 0:
        for i in range(0, max_step):
            mid = (a+b) / 2
            if abs(mid-mid_last)<eps or abs(fun(mid))<eps:
                print("Step: %d\nZero: %fc"%(i, mid))
                return mid
            else:
                if fun(mid)*fun(a)<0:
                    b = mid
                else:
                    a = mid
            mid_last = mid
        print('Bisection cannot be convergent within the pre-set steps.')


def fixed_point(fun, start:float=0, max_step:int=128, eps:float=1e-6)->float:
    """
    pass
    """
    new_val = fun(start)
    for i in range(0, max_step):
        old_val = new_val
        new_val = fun(old_val)
        if -eps<old_val-new_val<eps:
            print(i)
            return new_val


f = lambda x: x**2-3*x+1
print(Bisection(f, 0, 1))

f = lambda x: math.sqrt(10-x**3)/2
f = lambda x: math.sqrt(10/(x+4))
f = lambda x: x-(x**3+4*x**2-10)/(3*x**2+8*x)
result = fixed_point(f, 1)
print(result)