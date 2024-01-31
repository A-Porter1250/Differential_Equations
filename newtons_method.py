import numpy as np
from numpy import tan, cos, sin

def sec(x):
    return 1 / cos(x) 

def f(x):
    return 1/x

def df(x):
    return 1/x**2

def newton_method_answer(function, derivative, tolerance, max_iterations):
    x = 1
    while f(x) == ZeroDivisionError:
        x += 1
    for _ in range(max_iterations):
        x_new = x - function(x) / derivative(x)
        if abs(x_new - x) < tolerance:
            return round(x_new, 4)
        x = x_new
    return x