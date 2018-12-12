#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/09/27 01:00
# @Author   : Iydon
# @File     : Poly.py
import re
# str2poly



class Poly:
    """
    Polynomial.
    """
    # ------------- CONSTRUCTOR -------------
    def __init__(self, poly:list):
        """
        Initialize the class `Rational'.
        """
        while poly and poly[0]==0:
            poly.remove(0)
        if not poly:
            poly = [0]
        self.poly   = poly
        self.degree = len(poly)-1

    def __new__(self, poly:list):
        """
        New the object of `Rational'.
        """
        obj = super(Poly, self).__new__(self)
        obj.__init__(poly)
        return obj

    def __str__(self):
        """
        str(self)
        """
        parts = ["%s*X^%d"%(self.poly[-i-1],i) for i in range(self.degree+1)]
        return " + ".join(parts[::-1]).replace("*X^0", "")

    def __repr__(self):
        """
        repr(self)
        """
        return self.__str__()

    # ------------- PUBLIC -------------
    def eval(self, x):
        """
        Horner nested polynomial calculation.

        Args:
            poly: List, store the coefficient of the polynomial.
            x: Float, specify the variable in the polynomial.

        Returns:
            Float, result.

        Raises:
            If `poly' is empty, raise IndexError.
            If type(args) does not correspond, raise TypeError.
        """
        result = self.poly[0]
        for i in range(len(self)):
            result = x*result + self.poly[i+1]
        return result

    def lambdify(self):
        """
        Lambdify `self'.
        """
        return lambda x: self.eval(x)

    def diff(self):
        """
        Return the diff of the polynomial.
        """
        ranger = range(self.degree, 0, -1)
        result = [self.poly[-i-1]*i for i in ranger]
        return Poly(result)

    def integrate(self):
        """
        Return the integrate of the polynomial.
        """
        ranger = range(self.degree+1, 0, -1)
        result = [self.poly[-i]/i for i in ranger] + [0]
        return Poly(result)

    # ------------- OPERATOR -------------
    def __add__(self, value:object):
        """
        Return self+value.
        """
        if isinstance(value, Poly):
            if self.degree+1 >= len(value.poly):
                result = [i for i in self.poly]
                for i in range(len(value.poly)):
                    result[-i-1] += value.poly[-i-1]
            else:
                result = [i for i in value.poly]
                for i in range(self.degree+1):
                    result[-i-1] += self.poly[-i-1]
            return Poly(result)
        elif isinstance(value, list):
            return self.__add__(Poly(value))
        else:
            self.poly[-1] += value
            return Poly(self.poly)

    def __bool__(self):
        """
        self != 0.
        """
        return self.degree >= 0

    def __eq__(self, value:object):
        """
        Return self==value.
        """
        if isinstance(value, Poly):
            if len(self) == len(value):
                for e1,e2 in zip(self.poly,value.poly):
                    if e1 != e2:
                        return False
                return True
            else:
                return False
        elif isinstance(value, list):
            return self.__eq__(Poly(value))
        else:
            return self.__eq__(Poly([value]))

    def __hash__(self):
        """
        Return hash(self).
        """
        return hash(str(self))

    def __len__(self):
        """
        Return len(self).
        """
        return self.degree

    def __mul__(self, value:object):
        """
        Return self*value.
        """
        if isinstance(value, Poly):
            result = self.__convolution(self.poly, value.poly)
            return Poly(result)
        elif isinstance(value, list):
            return self.__mul__(Poly(value))
        else:
            return Poly([value*i for i in self.poly])

    def __ne__(self, value:object):
        """
        Return self!=value.
        """
        return not self==value

    def __neg__(self):
        """
        -self.
        """
        return Poly([-i for i in self.poly])

    def __pos__(self):
        """
        +self.
        """
        return Poly(self.poly)

    def __pow__(self, value:int, mod=None):
        """
        Return pow(self, value, mod).
        """
        if isinstance(value, int):
            if value==0:
                return Poly([1])
            elif value>0:
                result = self.poly
                for i in range(value-1):
                    result = self.__convolution(result, self.poly)
                return Poly(result)

    def __radd__(self, value:object):
        """
        Return value+self.
        """
        return self + value

    def __rmul__(self, value:object):
        """
        Return value*self.
        """
        return self * value

    def __rsub__(self, value:object):
        """
        Return value-self.
        """
        return -self + value

    def __truediv__(self, value):
        """
        Return self/value.
        """
        if isinstance(value, (Poly,list)):
            raise ArithmeticError
        else:
            return self * (1/value)

    def __sub__(self, value:object):
        """
        Return self-value.
        """
        return self + (-value)

    # ------------- PRIVATE -------------
    def __convolution(self, list1:list, list2:list) -> list:
        """
        Calculate the `convolution' between
        num1 and num2.
        ----------------------------
        Args:
            list1, list2: List.

        Returns:
            Convolution between lists.

        Raises:
            None.
        """
        n,m = len(list1)-1,len(list2)-1
        result = [0 for i in range(m+n+1)]
        for k in range(m+n+1):
            s = 0
            for i in range(max(0,k-m), min(k,n)+1):
                s += list1[-i-1]*list2[i-k-1]
            result[-k-1] = s
        return result



class Poly_Frac:
    """
    Fraction of Polynomial.
    """
    # ------------- CONSTRUCTOR -------------
    def __init__(self, poly:list, frac:list=[]):
        """
        Initialize.
        """
        self.poly = poly if isinstance(poly,Poly) else Poly(poly)
        for i in range(len(frac)):
            if frac[i][0]!=[0]:
                if isinstance(frac[i][0], Poly):
                    num,den = self.__fraction_simplify([frac[i][0],frac[i][1]])
                else:
                    num,den = self.__fraction_simplify([Poly(frac[i][0]),Poly(frac[i][1])])
                self.poly += num
                frac[i] = den if den[0]!=[0] else None
            else:
                frac[i] = None
        for i in range(len(frac)-1):
            if not frac[i] or frac[i][0]==[0]:
                frac[i] = None
                continue
            for j in range(i+1, len(frac)):
                if frac[j]:
                    if frac[i][-1]==frac[j][-1]:
                        frac[i][0] += frac[j][0]
                        frac[j] = None
        while frac and None in frac: frac.remove(None)
        self.frac = frac

    def __str__(self):
        """
        str(self)
        """
        decimal = ["(%s)/(%s)"%(f[0],f[-1]) for f in self.frac]
        return "%s%s%s"%(self.poly.__str__(), " + " if decimal else "", " + ".join(decimal))

    def __repr__(self):
        """
        repr(self)
        """
        return self.__str__()

    # ------------- OPERATOR -------------
    def __add__(self, value):
        """
        Return self+value.
        """
        if isinstance(value, Poly_Frac):
            pass
        else:
            return Poly_Frac(self.poly+value, self.frac)

    def __bool__(self):
        """
        self != 0.
        """
        pass

    def __eq__(self, value:object):
        """
        Return self==value.
        """
        pass

    def __hash__(self):
        """
        Return hash(self).
        """
        return hash(str(self))

    def __len__(self):
        """
        Return len(self).
        """
        pass

    def __mul__(self, value:object):
        """
        Return self*value.
        """
        pass

    def __ne__(self, value:object):
        """
        Return self!=value.
        """
        return not self==value

    def __neg__(self):
        """
        -self.
        """
        pass

    def __pos__(self):
        """
        +self.
        """
        pass

    def __pow__(self, value:int, mod=None):
        """
        Return pow(self, value, mod).
        """
        pass

    def __radd__(self, value:object):
        """
        Return value+self.
        """
        pass

    def __rmul__(self, value:object):
        """
        Return value*self.
        """
        pass

    def __rsub__(self, value:object):
        """
        Return value-self.
        """
        return -self + value

    def __truediv__(self, value):
        """
        Return self/value.
        """
        pass

    def __sub__(self, value:object):
        """
        Return self-value.
        """
        return self + (-value)

    # ------------- PRIVATE -------------
    def __fraction_simplify(self, frac):
        """
        Simplify the fraction of polynomial.
        ----------------------------
        Args:
            [Poly, Poly].
        """
        length = max(0, frac[0].degree-frac[-1].degree+1)
        result = [0 for i in range(length)]
        for i in range(length):
            result[i] = frac[0].poly[0]/frac[-1].poly[0]
            frac[0]  -= result[i] * Poly(frac[-1].poly+[0]*(length-i-1))
        return (result, [frac[0]/frac[-1].poly[0],frac[-1]/frac[-1].poly[0]])






'''

def str2poly(string:str):
    """
    1*x^3 + 0*x^2 + -2*x^1 + 1
    """
    return [float(re.findall("^[\d\.\-]+",part)[0]) for part in string.split(" + ")]

poly = [1, -2, 1]
frac = [[[1,3,3,2],[1,1]], [[1],[0,1]], [[1],[1,1]]]
p = Poly_Frac(poly, frac)
print(p)

from math import *

def interpolation(f, xs):
    """
    Interpolation.
    """
    def prod(lst, result=1):
        return prod(lst[:-1], result*lst[-1]) if lst else result
    result = Poly([0])
    for i in range(len(xs)):
        num = prod([Poly([1,-x]) for x in xs[:i]+xs[i+1:]])
        den = prod([(xs[i]-x) for x in xs[:i]+xs[i+1:]])
        result += f(xs[i]) * num / den
    return result

f = lambda x: cos(x)
res = interpolation(f, [0,0.9])
res = res.lambdify()
print(res(0.45))

'''
