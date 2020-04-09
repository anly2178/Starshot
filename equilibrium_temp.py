import numpy as np

def iter_factorial(n):
    result = 1
    value = n
    while value > 0:
        result = value * result
        value -= 1
    return result

print(iter_factorial(2))
print(iter_factorial(3))
print(iter_factorial(4))
print(iter_factorial(5))

"""
Attempt at iterative solver for equilibirum temperature
assuming absorption depends on temperature.
Need to know:
    • Absorption at some base temperature
    • How absorption varies with Temperature
    • Power of laser
    • Beam width
    • Area of sail
    • Emissivity
"""

def find_T_equilibrium():
