# -*- coding=utf-8 -*-

def chopping(number, digits):
    pass

def rounding(number, digits):
    pass

def operator(num1, num2, operator):
    pass


class machine_numbers(object):
	"""
	pass
	"""

	def __init__(self):
		pass

	def __str__(self):
		pass

	def __add__(self):
		pass

	def __mul__(self):
		pass

	def __len__(self):
		return 64

	def to_bits(self):
		pass

	def chopping(self, dig_num):
		pass

	def rounding(self, dig_num):
		pass


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