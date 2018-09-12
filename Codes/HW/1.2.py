#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/09/10 19:11
# @Author   : Iydon
# @File     : 1.2.py
# @Abstract : Homework auxiliary code.
import numpy as np


def error(p:float, p_star:float) -> tuple:
    """
    Compute the absolute error and
    relative error in approximations
    of p by p^\\star.
    ----------------------------

    Args:
        p, p_star: Float.

    Returns:
        (absolute error, relative error).

    Raises:
        None.
    """
    absolute_error = abs(p - p_star)
    relative_error = absolute_error / p
    return (absolute_error, relative_error)

def largest_interval_with_relative_error(p:float, error:float=1e-5, eps=1e-6) -> tuple:
    """
    Find the largest interval in which p^\\star
    must lie to approximate p with relative error
    at most `error' for each value of p.
    ----------------------------

    Args:
        p: Float, precise value.
        error: Float, like `epsilon'.
        eps: Float, epsilon.

    Returns:
        Interval: (a, b).

    Raises:
        None.
    """
    δ = error*p
    chopping = lambda x: eps * int(x/eps)
    interval_left  = chopping(p-δ+eps)
    interval_right = chopping(p+δ)
    return (interval_left, interval_right)

def chopping(p:float, sign_dig:int=6) -> float:
    """
    Chopping.
    ----------------------------
    Args:
        p: Float, `precise' value.
        sign_dig, Integer, significant digit(s).

    Returns:
        Float.

    Raises:
        None.
    """
    result = int(p*10**sign_dig) * 10**(-sign_dig)
    return round(result, sign_dig)

def rounding(p:float, sign_dig:int=6) -> float:
    """
    Rounding.
    ----------------------------
    Args:
        p: Float, `precise' value.
        sign_dig, Integer, significant digit(s).

    Returns:
        Float.

    Raises:
        None.
    """
    # `round(p, sign_dig).'
    result = p + 5*10**(-sign_dig-1)
    return chopping(result,sign_dig)

def bit2float(s:str="0"*64) -> float:
    """
    Use the 64-bit long real format to find the decimal
    equivalent of the following floating-point machine
    numbers.
    ----------------------------
    Args:
        s: String, only `0' and `1'.

    Returns:
        Float.

    Raises:
        None.
    """
    if len(s)==64:
        sign      = s[0:1]
        expontial = s[1:12]
        mantissa  = s[12:64]
        return (-1)**int(sign) * 2**(int(expontial,base=2)-1023) * (1+int(mantissa,base=2))

def solve_linear_system(_A:list, _b:list, sign_dig:int=6) -> list:
    """
    Solve the 2-by-2 linear systems using
    `sign_dig' rounding arithmetic.
    Solve `_A*x=_b'.
    ----------------------------
    Args:
        _A: List, 2-by-2.
        _b: List, 1-by-2.
        sign_dig, Integer, significant digit(s).

    Returns:
        List, 1-by-2, solutions.

    Raises:
        None.
    """
    if len(_A)==2 and len(_A[0])==2 and len(_A[1])==2 and len(_b)==2:
        fl = lambda x: round(x, sign_dig)
        a,b = [fl(i) for i in _A[0]]
        c,d = [fl(i) for i in _A[1]]
        e,f = [fl(i) for i in _b]
        m=fl(c/a); d1=fl(d-fl(m*b)); f1=fl(f-fl(m*e));
        y=fl(f1/d1); x=fl(fl(e-fl(b*y))/a)
        return (x,y)

class Rational:
    """
    `Rational' of Julia.
    Only consider real number.
    """
    # ------------- CONSTRUCTOR -------------
    def __init__(self, num:int, den:int):
        """
        Initialize the class `Rational'.
        """
        if den==0:
            string = "Denominator cannot be zero."
            raise ZeroDivisionError(string)
        elif den<0:
            string = "Denominator cannot be negative."
            raise ArithmeticError(string)
        div = self.__gcd(abs(num), den)
        self.num = num // div
        self.den = den // div
        self.__precise = 1e6

    def __new__(self, num:int, den:int):
        """
        New the object of `Rational'.
        """
        obj = super(Rational, self).__new__(self)
        obj.__init__(num, den)
        return obj

    def __str__(self):
        """
        str(self)
        """
        return "\\frac{%s}{%s}"%(self.num, self.den)

    # ------------- OPERATOR -------------
    def __abs__(self):
        """
        abs(self)
        """
        self.num = abs(self.num)

    def __add__(self, value:object):
        """
        Return self+value.
        """
        if isinstance(value, Rational):
            num = self.num*value.den + self.den*value.num
            den = self.den * value.den
            return self.__new__(Rational, num, den)
        elif isinstance(value, float):
            num,den = self.__float2rational(value)
            return self.__add__(Rational(num, den))
        elif isinstance(value, int):
            num = self.den*value + self.num
            return self.__new__(Rational, num, self.den)

    def __bool__(self):
        """
        self != 0.
        """
        return self.num != 0

    def __divmod__(self, value:int):
        """
        Return divmod(self, value).
        """
        div = int(float(self) / value)
        mod = self - div*value
        return (div, float(mod))
        # return divmod(float(self), value)

    def __eq__(self, value:object):
        """
        Return self==value.
        """
        if isinstance(value, Rational):
            return self.num==value.num and self.den==value.den
        elif isinstance(value, float):
            return float(self)==value
        elif isinstance(value, int):
            return float(self)==float(value)

    def __float__(self):
        """
        float(self).
        """
        return self.num / self.den

    def __floordiv__(self, value:int):
        """
        Return self//value.
        """
        return self.__new__(Rational, self.num, self.den*value)

    def __ge__(self, value:object):
        """
        Return self>=value.
        """
        return float(self) >= float(value)

    def __gt__(self, value:object):
        """
        Return self>value.
        """
        return float(self) > float(value)

    def __hash__(self):
        """
        Return hash(self).
        """
        return hash(float(self))

    def __int__(self):
        """
        int(self).
        """
        return self.num // self.den

    def __le__(self, value:object):
        """
        Return self<=value.
        """
        return float(self) <= float(value)

    def __lt__(self, value:object):
        """
        Return self<value.
        """
        return float(self) < float(value)

    def __mod__(self, value:int):
        """
        Return self%value.
        """
        div = int(float(self) / value)
        num = self.num - div*value*self.den
        return self.__new__(Rational, num, self.den)

    def __mul__(self, value:object):
        """
        Return self*value.
        """
        if isinstance(value, Rational):
            num = self.num * value.num
            den = self.den * value.den
            return self.__new__(Rational, num, den)
        elif isinstance(value, float):
            num,den = self.__float2rational(value)
            return self.__mul__(Rational(num, den))
        elif isinstance(value, int):
            return self.__new__(Rational, self.num*value, self.den)

    def __ne__(self, value:object):
        """
        Return self!=value.
        """
        return float(self) != float(value)

    def __neg__(self):
        """
        -self.
        """
        return self.__new__(Rational, -self.num, self.den)

    def __pos__(self):
        """
        +self.
        """
        return self
        # return self.__new__(Rational, self.num, self.den)

    def __pow__(self, value:object, mod=None):
        """
        Return pow(self, value, mod).
        """
        if isinstance(value, int):
            num = self.num ** value
            den = self.den ** value
            if mod:
                return self.__new__(Rational, num, den) % mod
            else:
                return self.__new__(Rational, num, den)
        elif isinstance(value, (Rational, float)):
            return float(self) ** float(value)

    def __radd__(self, value:object):
        """
        Return value+self.
        """
        return self + value

    def __rdivmod__(self, value:int):
        """
        Return divmod(value, self)
        """
        return divmod(value, float(self))

    def __rfloordiv__(self, value:int):
        """
        Return value//self.
        """
        num = self.den
        den = self.num
        return self.__new__(Rational, value*num, den)

    def __rmod__(self, value:object):
        """
        Return value%self.
        """
        return value % float(self)

    def __rmul__(self, value:object):
        """
        Return value*self.
        """
        return self * value

    def __rpow__(self, value:object, mod=None):
        """
        Return pow(value, self, mod).
        """
        return pow(value, float(self), mod)

    def __rsub__(self, value:object):
        """
        Return value-self.
        """
        return -self + value

    def __rtruediv__(self, value:object):
        """
        Return value/self.
        """
        num,den = self.den,self.num
        return Rational(num,den) * value

    def __sub__(self, value:object):
        """
        Return self-value.
        """
        return self + (-value)

    def __truediv__(self, value:object):
        """
        Return self/value.
        """
        if isinstance(value, Rational):
            num = self.num * value.den
            den = self.den * value.num
            return self.__new__(Rational, num, den)
        elif isinstance(value, int):
            den = self.den * value
            return self.__new__(Rational, self.num, den)
        elif isinstance(value, float):
            num,den = self.__float2rational(value)
            return self / Rational(num, den)

    # ------------- PRIVATE -------------
    def __gcd(self, num1:int, num2:int, flag=1) -> int:
        """
        Calculate the `gcd' between
        num1 and num2.
        ----------------------------
        Args:
            num1, num2: Integer.

        Returns:
            GCD between numbers.

        Raises:
            None.
        """
        if flag:
            return self.__gcd_recursive(num1, num2)
        else:
            result = 1
            minimu = min(num1, num2)
            for i in range(2, minimu+1):
                if num1%i==0 and num2%i==0:
                    result = i
            return result

    def __gcd_recursive(self, num1:int, num2:int) -> int:
        """
        Semiliar to `__gcd'.
        """
        if num1 % num2 == 0:
            return num2
        else:
            return self.__gcd_recursive(num2, num1%num2)

    def __float2rational(self, value:float) -> tuple:
        """
        Convert the float value to `Rational'.
        ----------------------------
        Args:
            value: Float.

        Returns:
            Tuple, (num, den).

        Raises:
            None.
        """
        num = round(value*self.__precise)
        den = int(self.__precise)
        return (num, den)


# ------------- Homework -------------
title = """
1. Compute the absolute error and relative error in approximations of p by p^\\star.
"""
print(title)
    # a
result = error(np.pi, 22/7)
print("    a) \\pi and 22/7:    ",result)
    # c
result = error(np.e, 2.718)
print("    c) e and 2.718:     ", result)
    # e
result = error(np.exp(10), 22000)
print("    e) e^{10} and 22000:", result)
    # f
result = error(8*7*6*5*4*3*2*1, 39900)
print("    g) 8! and 39900:    ", result)

title = """
2. Find the largest interval in which p^\\star must lie to approximate p with relative error at most 1e-4 for
   each value of p.
"""
print(title)
    # a
result = largest_interval_with_relative_error(np.pi, 1e-4)
print("    a) The largest interval of \\pi:     ", result)
    # c
result = largest_interval_with_relative_error(np.sqrt(2), 1e-4)
print("    c) The largest interval of \\sqrt(2):", result)

title = """
4. Perform the following computations (i) exactly, (ii) using three-digit chopping arithmetic, and (iii)
   using three-digit rounding arithmetic. (iv) Compute the relative errors in parts (ii) and (iii).
"""
print(title)
    # a
r = Rational(4,5) + Rational(1,3)
chopping_val = chopping(chopping(4/5,3) + chopping(1/3), 3)
rounding_val = rounding(rounding(4/5,3) + rounding(1/3), 3)
print("    a) Exactly: %s, 3-digit chopping: %.3f, 3-digit rounding: %.3f, \n\
       relative error: %f and %f." \
       %(r, chopping_val, rounding_val, error(float(r),chopping_val)[-1], error(float(r),rounding_val)[-1]))
    # c
r = (Rational(1,3) - Rational(3,11)) + Rational(3,20)
chopping_val = chopping(chopping(chopping(1/3,3)-chopping(3/11,3),3) + chopping(3/20,3), 3)
rounding_val = rounding(rounding(rounding(1/3,3)-rounding(3/11,3),3) + rounding(3/20,3), 3)
print("    c) Exactly: %s, 3-digit chopping: %.3f, 3-digit rounding: %.3f, \n\
       relative error: %f and %f." \
       %(r, chopping_val, rounding_val, error(float(r),chopping_val)[-1], error(float(r),rounding_val)[-1]))


title = """
15. Use the 64-bit long real format to find the decimal equivalent of the following floating-point machine
    numbers.
"""
print(title)
    # a
bit64 = "0 10000001010 1001001100000000000000000000000000000000000000000000".replace(" ","")
print("    a)", bit64, "=", bit2float(bit64))
    # c
bit64 = "0 01111111111 0101001100000000000000000000000000000000000000000000".replace(" ","")
print("    a)", bit64, "=", bit2float(bit64))

title = """
19. Solve the following linear systems using four-digit rounding arithmetic.
"""
print(title)
    # a
A = [[1.130,-6.990],[1.013,-6.099]]
b = [14.20, 14.22]
solution = solve_linear_system(A, b, sign_dig=4)
print("    a) Solution of linear system:", solution)

title = """
24. Suppose that fl(y) is a k-digit rounding approximation to y. Show that
    \\left|\\frac{y-fl(y)}{y}\\right| \\leq 0.5\\times 10^{-k+1}.
"""
print(title)
proof = """
If $d_{k+1}<5$, then $fl(y)=0.d_1d_2\\ldots d_k\\times 10^n$.
\\begin{align*}
\\left|\\frac{y-fl(y)}{y}\\right| &= \\left|\\frac{d_{k+1}.d_{k+2}\\ldots\\times 10^{n-k+1}}{0.d_1d_2\\ldots\\times 10^n}\\right| \\\\
&= \\frac{d_{k+1}.d_{k+2}\\ldots}{0.d_1d_2\\ldots}\\times 10^{-k+1} \\\\
&\\leq 0.5\\times 10^{-k+1}
\\end{align*}
If $d_{k+1}\\geq 5$, then $fl(y)=0.d_1d_2\\ldots d_k\\times 10^n+10^{n-k}$.
\\begin{align*}
\\left|\\frac{y-fl(y)}{y}\\right| &= \\left|\\frac{1-0.d_{k+1}d_{k+2}\\ldots\\times 10^{n-k}}{0.d_1d_2\\ldots\\times 10^n}\\right| \\\\
&= \\frac{1-0.d_{k+1}d_{k+2}\\ldots}{0.d_1d_2\\ldots}\\times 10^{-k} \\\\
&\\leq 0.5\\times 10^{-k+1}
\\end{align*}
"""
print(proof)