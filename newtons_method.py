import numpy as np
from numpy import tan, cos, sin

def sec(x):
    return 1 / cos(x) 

def f(x):
    return np.tan(x) - 1/x

def df(x):
    return sec(x)**2 + 1/x**2

def newton_method(function: function, derivative: function, tolerance, max_iterations):
    x = 1
    while isinstance(function, function) == KeyError:
        x +=1
    for _ in range(max_iterations):
        x_new = x - f(x) / df(x)
        if abs(x_new - x) < tolerance:
            return round(x_new, 4)
        x = x_new
    return None