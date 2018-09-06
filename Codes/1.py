# -*- coding=utf-8 -*-

def nested(poly:list=[1], x:float=0.0)->float:
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
    result = poly[0]
    for i in range(1, poly.__len__()):
        result = x*result + poly[i]
    return result


poly = [1,-2,1]
print(nested(poly, []))
print(poly)