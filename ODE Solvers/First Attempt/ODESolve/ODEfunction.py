#ODE solving functions

# Import libraries
import sympy as sp

# Checks if a given ode is linear
def is_linear_ode(ode, dep_var, ind_var):
    # Set condition marker as 'True'
    is_linear = True

    # Condition 1: Check if all coefficients are functions of independent variables ONLY
    ode_lhs = ode.lhs
    ode_rhs = ode.rhs
    
    coefficients_dict_lhs = ode_lhs.as_coefficients_dict() 
    coefficients_dict_rhs = ode_rhs.as_coefficients_dict()
    
    coefficients_list_lhs = [coeff for term, coeff in coefficients_dict_lhs.items()]
    coefficients_list_rhs = [coeff for term, coeff in coefficients_dict_rhs.items()]
    
    term_list_lhs = [term for term, coeff in coefficients_dict_lhs.items()]
    term_list_rhs = [term for term, coeff in coefficients_dict_rhs.items()]


    for item in coefficients_list_lhs:
        if item.has(dep_var):
            is_linear = False
    for item in coefficients_list_rhs:
        if item.has(dep_var):
            is_linear = False

    # Condition 2: Check if all derivative terms are of degree one or zero
    for item in term_list_lhs:
        if item.args[1][0] > 1:
            is_linear = False

    # Condition 3: Check if there are no products of derivative terms
    if poly.degree() > 1:
        is_linear = False

    # Condition 4: Check if coefficient a_n is non-zero for any value of independent 
    #               variable in the solution range
    coeff_an = poly.coeff_monomial(dep_var.diff(ind_var, poly.degree()))
    if coeff_an == 0:
        is_linear = False

    # Condition 5: Check if the last term on the left-hand side is the dependent variable 
    #                   and not a function of the dependent variable
    last_term = poly.as_expr().as_coefficients_dict()[dep_var.diff(ind_var, poly.degree())]
    if last_term.has(dep_var):
        is_linear = False

    return is_linear

# Check if the ODE is homogeneous
def is_homogeneous(ode, dep_var, ind_var):
    """
        If of the form: 
        a.n(x)*(d^n)/dx^n + a.(n-1)(x)*a^(n-1)y/dx^(n-1) + ⋯ + a.1(x)*dy/dx + a.0(x)y = F(x)
        Homogeneous ODE if F(x) = 0
        Elif dy/dx = F(x,y) = F(hy,hx)
        Else forced
    """
    
    # Intitialize check variable
    ode_force_check = False
    
    # Check condition
    if ode.sp.rhs == 0:
        ode_force_check = True
    
    return ode_force_check

# Solve ODE as seperable
def solve_seperable(ode, y, x):
        
    seperated_eq_1 = sp.eq(ode.lhs/sp.diff(y),ode.rhs/sp.diff(y))
    seperated_eq_2 = sp.eq()
        
    # Separate the variables and integrate
    separated_eq = sp.integrate(1 / y, y) - sp.integrate(x, x)

    # Solve for y(x)
    solution = sp.simplify(separated_eq)
    
    solution = sp.dsolve(ode, y(x))
    
    return solution

def solve_homogeneous_sub(ode, y, x):

    # Perform the homogeneous substitution: y = x * v
    v = symbols('v')(x)
    homogeneous_ode = ode.subs(y, x*v)

    # Solve the homogeneous ODE
    homogeneous_solution = solve_seperable(homogeneous_ode)

    # Substitute back the original variable
    solution = [y.subs(v, sol) for sol in homogeneous_solution]
    
    return solution

def check_exact(ode):
    
    # Make nominator and denominator variables then multiply
    
    #nominator = #
    #denominator = #
    
    #ode1 = Eq(ode.sp.lhs*denominator, ode.sp.rhd*denominator)
    
    # Extract the dependent variable and the functions from the ODE
    x = sp.Symbol('x')
    y = sp.Function('y')(x)
    M = ode.coeff(sp.diff(y, x), 1)
    N = ode.rhs

    # Check if the partial derivatives satisfy the condition for exactness
    dM_dy = sp.diff(M, y)
    dN_dx = sp.diff(N, x)
    
    return dM_dy == dN_dx

def solve_exact_ODE(ode, y, x):
    # Check if the ODE is exact
    if not ode.check_exact():
        print("The ODE is not exact. Cannot solve using this method.")
        return

    # Determine M(x, y) and N(x, y)
    M = ode.coeff(sp.diff(y, x), 1)
    N = ode.rhs

    # Integrate M with respect to x
    U_M = sp.integrate(M, x)

    # Integrate N with respect to y
    U_N = sp.integrate(N, y)

    # Integrating Constant
    c = sp.Symbol('c')

    # Solve for the function U(x, y)
    U = U_M + U_N

    # Solve for the solution
    solution = Eq(c, U)

    return solution

def solve_exact_int_factor(ode, x, y):

    # Make function based on this procedure:
    """
        ODE made exact using Integration Factor
        ∂M/∂y = ∂/∂y(∂U/∂x) <> ∂N/∂x = ∂/∂y(∂U/∂y)
        u(x,y) factor multiplies through ODE to make it exact
        ∂(uM)/∂y = ∂(uN)/∂x
        Four methods to find integration factor
        1. Factor M&N into functions of x,y (M->f(x)*f(y) N->g(x)*g(y)) 
            seperatly and multiply ODE by 1/(f(x)*g(y)) or 1/(f(y)*g(x)) to make equatioin seperable
            - Directly integrate
        2. Find if (1/N)(∂M/∂y-∂N/∂x)=f(x) only, then u(x)=e^∫f(x)dx
        3. Find if (1/M)(∂N/∂x-∂M/∂y)=f(y) only, then u(y)=e^∫f(y)dy
        4. Assume form u(x,y) = (x^a)*(y^b), multiply through ODE and solve for a, b to make equation exact
    """
    
    # Make sure ode is in form N(x) * (dy/dx) = M(x)  

    M = sp.Symbol('M')(x, y)
    N = sp.Symbol('N')(x, y)
    
    M = ode.sp.coeff(diff(y, x), 1)
    N = ode.sp.rhs
    
    check_1 = sp.Eq((1/M)*(sp.diff(N, x) - sp.diff(M, y)), 1)
    check_2 = sp.Eq((1/N)*(sp.diff(M, y) - sp.diff(N, x)), 1)
    
    #if check_1 #is function of y ONLY
    #if check_2 #is function of x ONLY
    
    #If
    
    
    return

def solve_linear_ODE(ode, is_linear, y, x):

    p = ode.sp.coeff(diff(y, x))
    g = ode.sp.rhs
    c = sp.Symbol('c')
    
    # Solve the linear ODE using the first-order method
    # [e^int(p(x)dx)]*[int(g(x)*e^int(p(x)dx))dx+c]
    solution = sp.Eq(y, sp.E**sp.integrate(p,x) * sp.itegrate(g*E**sp.integrate(p,x),x))
        
    return solution

def linearize_ODE(ode, is_linear, y, x):
    # Linearize the non-linear ODE using the alternative methods
        
    # Option 1: Invert dependent/independent variables
    linear_ode_option1 = ode.subs(y, 1/x)
    solution_option1 = solve(linear_ode_option1, y)
        
    # Option 2: Use a function of the dependent variable
    linear_ode_option2 = ode.subs(y, Function('f')(x))
    linear_ode_option2 = simplify(linear_ode_option2.diff(x)/linear_ode_option2)
    solution_option2 = solve(linear_ode_option2, Function('f')(x))
        
    # Option 3: Use a special variable substitution (Bernoulli ODE)
    z = symbols('z')
    linear_ode_option3 = ode.subs(y, z/x)
    linear_ode_option3 = simplify(linear_ode_option3.diff(x)/linear_ode_option3)
    solution_option3 = solve(linear_ode_option3, z)
        
    # Print the linearized solutions
    print("Linearized Solution (Option 1):", solution_option1)
    print("Linearized Solution (Option 2):", solution_option2)
    print("Linearized Solution (Option 3):", solution_option3)

# Missing term (homogeneous ODE)
def RoO_homogeneous_ODE(ode, y, x):
        
    # Assume a second linearly independent solution v(x)
    v = Function('v')(x)
        
    # Integrate above functions to solve ode with substituted v(x)
        
    # Return the solutions
    return solution

# Case 2: Factor operators
def RoO_Factor_Operators(ode, y, x):
    
    # Assume a second linearly independent solution v(x)
    v = Function('v')(x)
    
    # sub v = (D+b)y, solve for v, solve for y
    # Use above functions to accomplish this

    return solution

# Case 3: D'Alembert method (non-homogeneous ODE)
def RoO_DAlemberts_method(ode, y, x): 
    # Assume a second solution of the form v(x) * y1(x), where y1(x) is the first solution
    v = Function('v')(x)
    
    # Sub in ode and get ode in v(x), Solve for v(x) then solution is Y(x) = y.1(x) + y.2(x)
    
    return solution