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

def bissection_method(tolerance, max_iterations, range, func):
    #check convergence condition
    x = (range[0] + range[1]) / 2
    if abs(((func(x+0.0001)-func(x-0.0001))/(0.0002))) > 1:
        print('local derivative causes divergence')
        return False


    #range in form [upper bound, lower bound]
    error = np.inf
    iterations = 0
    while(iterations < max_iterations):
        x = (range[0] + range[1]) / 2
        print(x)
        val = func(x)
        if abs(val) < tolerance:
            return x
        elif val > 0:
            range[0] = x
        else:
            range[1] = x
        iterations += 1
    print('unable to converge in max iterations')
    return False