import matplotlib.pyplot as plt
import numpy as np

def bisection_step(lower_bound, upper_bound, func):
    midpoint = (lower_bound + upper_bound) / 2
    f_midpoint = func(midpoint)
    return midpoint, f_midpoint

def bisection(lower_bound, upper_bound, func, max_iterations, tolerance, show_plot=True):
    a = lower_bound
    b = upper_bound
    fa = func(a)
    as_values = np.zeros(max_iterations)
    bs_values = np.zeros(max_iterations)
    midpoints = np.zeros(max_iterations)
    f_midpoints = np.zeros(max_iterations)

    for k in range(max_iterations):
        midpoint, f_midpoint = bisection_step(a, b, func)
        as_values[k] = a
        bs_values[k] = b
        midpoints[k] = midpoint
        f_midpoints[k] = f_midpoint

        if abs(f_midpoint) < tolerance:
            return midpoint

        if fa * f_midpoint > 0:
            a = midpoint
            fa = f_midpoint
        else:
            b = midpoint

    return None  # No root found within the maximum iterations

def nr_step(prev_x, func, dfunc):
    return prev_x - func(prev_x) / dfunc(prev_x)

def newton_raphson(x_guess, func, dfunc, max_iteration, tolerance):
    p0 = x_guess
    for _ in range(max_iteration):
        p1 = nr_step(p0, func, dfunc)

        if abs(p1 - p0) < tolerance:
            return p1

        p0 = p1

    return None  # No root found within the maximum iterations

# Bisection where f(x) = x^2 - 4
a = -1  # value of a should be so that f(a) is negative
b = 100  # value of b should be so that f(b) is positive
TOL = 0.1
MaxIteration = 100
f = lambda x: x**2 - 4
zero = bisection(a, b, f, MaxIteration, TOL)
print(f'Bisection method solution: {zero:.4f}')
