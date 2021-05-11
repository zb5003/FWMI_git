"""
Test out some of Python's integration methods.
"""
import scipy.integrate as spi
import numpy as np

def test_func(x, y, func=None):
    if func is not None:
        x = func(x)
    return x * y

def test_2(x):
    return x**3

integral = spi.nquad(test_func, [[0, 1], [0, 1]], [test_2])
print(integral)